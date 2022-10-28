use crate::png::lzss::{lzss, LzssElement};

use super::lzss::LzssBuffer;

/// Writer to allow appending values to a u8 vec that may be smaller than one byte.
struct DeflateWriter {
    /// The bytes to write to a file.
    data: Vec<u8>,

    /// A temporary buffer of bits not yet added to the data vector.
    byte_to_write: [bool; 8],

    /// The number of stored bits in the byte_to_write array, when that array
    /// fills up it will be written to data.
    stored_bits: usize,
}

impl DeflateWriter {
    /// Create a new empty deflate writer
    fn new() -> Self {
        Self {
            data: vec![],
            byte_to_write: [false; 8],
            stored_bits: 0,
        }
    }

    /// Get the number of bytes of data in the writer
    fn len(&self) -> usize {
        let len = self.data.len();
        if self.stored_bits == 0 {
            len
        } else {
            len + 1
        }
    }

    /// Finish writing the data and return the data stream
    fn finish(mut self) -> Vec<u8> {
        self.finish_byte();
        self.data
    }

    /// fill the remaining space in the temporary buffer with zeros, so the next
    /// bit written is in a new byte.
    fn finish_byte(&mut self) {
        while self.stored_bits != 0 {
            self.bit(false);
        }
    }

    /// Write a single bit to the vector
    fn bit(&mut self, bit: bool) {
        self.byte_to_write[self.stored_bits] = bit;
        self.stored_bits += 1;

        // temporary buffer is full, write the byte out
        if self.stored_bits >= 8 {
            // Bits should be written with the first one on the right in the final
            // stream, the reverse of the order they are added to the temporary buffer.
            self.byte_to_write.reverse();

            let mut result = 0;

            for (idx, bit) in self.byte_to_write.into_iter().enumerate() {
                result |= (bit as u8) << idx;
            }

            self.data.push(result.reverse_bits());
            self.stored_bits = 0;
        }
    }

    /// Write the low order `bit_count` bits of data in little endian order (least
    /// significant bit first).
    fn write_little_endian(&mut self, mut data: u8, bit_count: usize) {
        for _ in 0..bit_count {
            self.bit((data & 1) != 0);
            data >>= 1;
        }
    }

    /// Write the low order `bit_count` bits of data in big endian order (most
    /// significant bit first).
    fn write_big_endian(&mut self, mut data: u8, bit_count: usize) {
        for _ in 0..bit_count {
            self.bit(((data >> 7) & 1) != 0);
            data <<= 1;
        }
    }

    /// Write lzss encoded data to the file using the given huffman code
    fn write_huffman_lzss(&mut self, data: &[LzssElement], code: &HuffmanCode) {}

    /// Zero fill the remaining bits of the temporary buffer and push the given
    /// bytes directly to the output.
    fn write_bytes(&mut self, data: &[u8]) {
        self.finish_byte();
        self.data.extend(data);
    }
}

/// The two Huffman codes used to encode an lzss buffer, represented as
/// (data, bit length of data) pairs where only the least significant `bit length`
/// bits of the data are written to an output buffer.
struct HuffmanCode {
    literal: Vec<(u16, u8)>,
    distance: Vec<(u16, u8)>,
}

impl HuffmanCode {
    /// construct the fixed huffman code as stated in rfc 1951:
    /// - 0-143 => bit length 8,
    /// - 144-255 => bit length 9,
    /// - 256-279 => bit length 7,
    /// - 280-287 => bit length 8,
    /// - all distance codes => bit length 5,
    fn fixed() -> Self {
        // get the bit length lists
        let rep = |val, count| std::iter::once(val).cycle().take(count);

        let literal: Vec<_> = rep(8, (0u16..=143).len())
            .chain(rep(9, (144u16..=255).len()))
            .chain(rep(7, (256u16..=279).len()))
            .chain(rep(8, (280u16..=287).len()))
            .collect();
        let distance: Vec<_> = rep(5, 32).collect();

        let literal_code = canonical_huffman(&literal);
        let distance_code: Vec<_> = (0..=31).collect();

        Self {
            literal: literal_code.into_iter().zip(literal).collect(),
            distance: distance_code.into_iter().zip(distance).collect(),
        }
    }

    /// Convert the frequency data created when writing an lzss buffer into a
    /// huffman code
    fn from_lzss(data: &LzssBuffer) -> Self {
        let literal = package_merge(&data.literal_length);
        let distance = package_merge(&data.distance);
        let literal_code = canonical_huffman(&literal);
        let distance_code = canonical_huffman(&distance);

        Self {
            literal: literal_code.into_iter().zip(literal).collect(),
            distance: distance_code.into_iter().zip(distance).collect(),
        }
    }

    /// Get the number of distance codes that participated in the huffman code
    fn distance_code_count(&self) -> u32 {
        // no literal codes are removed, even the ones with 0 probability are
        // included in the huffman tree (todo?)
        self.distance.len() as _
    }

    /// Get the number of distance codes that participated in the huffman code
    fn literal_code_count(&self) -> u32 {
        // no literal codes are removed, even the ones with 0 probability are
        // included in the huffman tree (todo?)
        self.literal.len() as _
    }
}

/// Implement the DEFLATE compression algorithm.
/// Assumes that only one DEFLATE block will be written, either with no compression,
/// fixed huffman code as given in rfc 1951 or a dynamic huffman code.  No attempt
/// at encoding different blocks is made.
pub fn deflate(data: &[u8]) -> Vec<u8> {
    // create output stream
    let mut output_buffer = DeflateWriter::new();

    // convert the input into lzss encoded data with a 32k block size
    let lzss = lzss(7, data);

    // create the huffman codes for the two alphabets if using dynamic coding
    let huffman = HuffmanCode::from_lzss(&lzss);

    // encode the dynamic header and the data encoded using the dynamic coding.
    {
        let header_alphabet = HeaderAlphabet::new(&huffman);
        let header_code = header_alphabet.order_bit_lengths();

        // below subtractions only valid while not removing unused codes from
        // either the literal/length or the distance alphabets. (todo)
        let h_lit = huffman.literal_code_count() - 257;
        let h_dist = huffman.distance_code_count() - 1;

        // get the number of zeros at the end of the header code
        let trailing_zeros = header_code.iter().rev().take_while(|i| **i == 0).count();
        // the header requires at least 4 header code elements remaining out of 19
        let trailing_zeros = trailing_zeros.max(19 - 4) as u32;
        // convert trailing zeros to number of remaining elements - 4 for the header
        let h_c_len = (19 - trailing_zeros) - 4;

        //     // always encoding only one block (todo?) so always record as final block.
        //     output.push_bits(&[true]);

        //     //
        //     output.push_bits(data);
    }

    // if the dynamic header is > 10% of the encoded block, try encoding the block
    // using the fixed huffman coding to see if it produces a shorter encoding and
    // if so use the fixed one instead.
    let fixed_huffman = HuffmanCode::fixed();
    let mut fixed_buffer = DeflateWriter::new();
    fixed_buffer.bit(true); // this is the final block
    fixed_buffer.write_little_endian(0b01, 2); // fixed compression block
    fixed_buffer.write_huffman_lzss(&lzss.data, &fixed_huffman);

    // if the compressed data length >= original data length then discard compression
    // and output the data as an uncompressed block.

    // uncompressed blocks are limited to 65535 bytes each, so split data into
    // multiple blocks.
    let chunk_max_size = 65535;

    // 5 bytes of header for each block + the data its self
    let uncompressed_data_len = (data.len() / chunk_max_size + 1) * 5 + data.len();

    if true || fixed_buffer.len() > uncompressed_data_len {
        let mut uncompressed = DeflateWriter::new();
        let last_block_idx = data.len() / chunk_max_size;
        for (block_idx, data) in data.chunks(chunk_max_size).enumerate() {
            uncompressed.bit(block_idx == last_block_idx); // is this the final block
            uncompressed.write_little_endian(0b00, 2); // non compression block

            uncompressed.write_bytes(&(data.len() as u16).to_le_bytes()); // length of data
            uncompressed.write_bytes(&(!(data.len() as u16)).to_le_bytes()); // one's complement
            uncompressed.write_bytes(data);
        }

        output_buffer = uncompressed;
    }

    output_buffer.finish()
}

/// Implementation of the package merge algorithm for calculating the number of
/// bits to use in a length-limited huffman encoding. Input data is the number of
/// occurrences of each symbol and output will have the bit length written to it,
/// so must be the same length as the input data.  If filter is true, then any data
/// points listed as having 0 occurrence or probability will not participate in
/// the construction of the huffman code.
fn package_merge(data: &[usize]) -> Vec<u8> {
    // the individual symbols in use
    let mut symbols: Vec<_> = data
        .iter()
        .enumerate()
        .map(|(idx, prob)| (*prob, vec![idx]))
        .collect();

    if symbols.is_empty() {
        return vec![0; data.len()];
    }

    let minimum_packages = 2 * symbols.len() - 2;

    while symbols.len() < minimum_packages {
        // sort in reverse order, lowest probability first
        symbols.sort_unstable_by_key(|a| a.0);

        // package adjacent symbols, if length is odd, the largest probability
        // symbol will be dropped.
        symbols = symbols
            .chunks_exact(2)
            .map(|chunk| {
                let a = &chunk[0];
                let b = &chunk[1];
                let nums = a.1.iter().copied().chain(b.1.iter().copied()).collect();

                (a.0 + b.0, nums)
            })
            .collect();

        // merge original symbols into the list
        symbols.extend(
            data.iter()
                .enumerate()
                .map(|(idx, prob)| (*prob, vec![idx])),
        );
    }

    // count occurrences of each original symbol and write to output
    let mut output = vec![0; data.len()];
    for (_, symbol) in symbols {
        for num in symbol {
            output[num] += 1;
        }
    }

    output
}

/// Converts a list of bit lengths of huffman codes into the bits of data to be
/// written or read using a canonical huffman encoding
fn canonical_huffman(bit_lengths: &[u8]) -> Vec<u16> {
    struct Data {
        idx: usize,
        bit_length: u8,
        output: u16,
    }

    let mut symbols: Vec<_> = bit_lengths
        .iter()
        .enumerate()
        .map(|(idx, &bit_length)| Data {
            idx,
            bit_length,
            output: 0,
        })
        .filter(|e| e.bit_length > 0)
        .collect();

    // sort by bit length, then by by symbol index
    symbols.sort_unstable_by(|a, b| a.bit_length.cmp(&b.bit_length).then(a.idx.cmp(&b.idx)));

    let mut prev_len = 0;
    let mut prev_code = None;
    for symbol in &mut symbols {
        let mut new_code = prev_code.map(|c: u16| c + 1).unwrap_or(0);

        while prev_len < symbol.bit_length {
            new_code <<= 1;
            prev_len += 1;
        }
        prev_code = Some(new_code);
        symbol.output = new_code;
    }

    // copy the output data to a new array, effectively doing the reverse of the
    // map and sort done at the start
    let mut output = vec![0; bit_lengths.len()];
    for symbol in symbols {
        output[symbol.idx] = symbol.output;
    }

    output
}

/// Data returned when encoding using the header alphabet
struct HeaderAlphabet {
    /// The data encoded using the header alphabet
    data: Vec<HeaderData>,

    /// The bit lengths of the header alphabet's huffman code in order 0..=18
    bit_lengths: [u8; 19],

    /// The huffman code values for each element of the header alphabet
    encoding: [u8; 19],
}

/// Simple encoding of the header alphabet as an enum, rather than minimal bit lengths
#[derive(Debug)]
enum HeaderData {
    /// Code lengths from 0 to 15 as literals
    CodeLength(u8),

    /// Repeat the previous code between 3 and 6 times
    Repeat(u8),

    /// Repeat the digit zero between 3 and 10 times
    RepeatZero(u8),

    /// Repeat the digit zero between 11 and 138 times
    RepeatZeroLong(u8),
}

impl HeaderAlphabet {
    /// Encode the huffman code lengths for the literal and distance alphabets using
    /// the code length alphabet as defined in rfc 1951.
    fn new(input: &HuffmanCode) -> HeaderAlphabet {
        // the data to encode is actually 1 stream containing both alphabets lengths,
        // so repetitions can cross between the two.
        let mut header_data = Vec::with_capacity(input.literal.len() + input.distance.len());
        header_data.extend(input.literal.iter().map(|x| x.1));
        header_data.extend(input.distance.iter().map(|x| x.1));

        // Encode the data using the header alphabet: because we are not removing
        // anything from either the literal or distance alphabets, even if it is not
        // used (todo?) it is very unlikely to ever need the RepeatZero or RepeatLongZero
        // data points, but regardless, try to encode them anyway.

        let mut data = vec![];

        let mut idx = 0;
        while idx < header_data.len() {
            match header_data[idx] {
                0 => {
                    // find as many zeros as possible
                    let zero_count = header_data[idx..]
                        .iter()
                        .take(138) // there isn't any encoding for 139 zeros...
                        .take_while(|e| **e == 0)
                        .count();
                    idx += zero_count;
                    match zero_count {
                        1 | 2 => {
                            for _ in 0..zero_count {
                                data.push(HeaderData::CodeLength(0))
                            }
                        }
                        3..=10 => data.push(HeaderData::RepeatZero(zero_count as _)),
                        11..=138 => data.push(HeaderData::RepeatZeroLong(zero_count as _)),
                        _ => panic!(),
                    }
                }
                item @ 1..=15 => {
                    let count = header_data[idx..]
                        .iter()
                        .take(7) // first item = 1, repeat = 3..=6
                        .take_while(|e| **e == item)
                        .count();
                    idx += count;
                    match count {
                        1 | 2 | 3 => {
                            for _ in 0..count {
                                data.push(HeaderData::CodeLength(item))
                            }
                        }
                        4..=7 => {
                            data.push(HeaderData::CodeLength(item));
                            data.push(HeaderData::Repeat(count as u8 - 1))
                        }
                        _ => panic!(),
                    }
                }
                _ => panic!(),
            }
        }

        // compute the encoding frequencies
        let mut frequencies = [0usize; 19];
        for elem in &data {
            match elem {
                HeaderData::CodeLength(u) => frequencies[*u as usize] += 1,
                HeaderData::Repeat(_) => frequencies[16] += 1,
                HeaderData::RepeatZero(_) => frequencies[17] += 1,
                HeaderData::RepeatZeroLong(_) => frequencies[18] += 1,
            }
        }

        // remove the unused symbols from the alphabet
        let mut reduced_alphabet = [(0, 0); 19];
        let mut idx = 0;
        for (symbol_idx, &freq) in frequencies.iter().enumerate() {
            if freq != 0 {
                reduced_alphabet[idx] = (symbol_idx, freq);
                idx += 1;
            }
        }

        // convert the reduced alphabet into a huffman code
        let reduced_frequencies = reduced_alphabet.map(|(_, freq)| freq);
        let reduced_frequencies = &reduced_frequencies[0..idx];
        let reduced_bit_lengths = package_merge(reduced_frequencies);
        let reduced_encoding = canonical_huffman(&reduced_bit_lengths);

        // expand the reduced alphabets huffman code adding zeros for unused symbols.
        let mut bit_lengths = [0; 19];
        let mut encoding = [0; 19];
        for ((&(symbol, _freq), bit_length), bit_encoding) in reduced_alphabet
            .iter()
            .zip(reduced_bit_lengths)
            .zip(reduced_encoding)
        {
            bit_lengths[symbol] = bit_length;
            encoding[symbol] = bit_encoding as _;
        }

        HeaderAlphabet {
            data,
            bit_lengths,
            encoding,
        }
    }

    /// order the symbol code lengths as specified for the dynamic huffman code
    /// header in rfc 1951
    fn order_bit_lengths(&self) -> [u8; 19] {
        let idx = [
            16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
        ];

        idx.map(|idx| self.bit_lengths[idx])
    }
}

struct ZlibHeader {
    compression_info: u8,
    compression_level: u8,
}

pub fn zlib(data: &[u8]) -> Vec<u8> {
    let mut output = vec![];
    output.extend_from_slice(
        &ZlibHeader {
            compression_info: 7,
            compression_level: 0,
        }
        .as_bytes(),
    );

    output.extend(deflate(&data));
    output.extend(adler32(data.as_ref()).to_be_bytes());

    output
}

fn adler32(data: &[u8]) -> u32 {
    let mut s1 = 1;
    let mut s2 = 0;

    for ch in data {
        s1 = (s1 + *ch as u32) % 65521;
        s2 = (s2 + s1) % 65521;
    }

    (s2 << 16) + s1
}

impl ZlibHeader {
    fn as_bytes(&self) -> [u8; 2] {
        let method = 8 & 0b1111; // method 8 = deflate
        let info = self.compression_info & 0b1111; // lz77 window size
        let cmf = info << 4 | method;

        let level = self.compression_level & 0b11;
        let dict = 0; // no preset dictionary
        let flg = level << 6 | dict << 5;

        let rem = ((cmf as u16) * 256 + flg as u16) % 31;
        let flg = flg | (31 - rem as u8);

        assert_eq!(0, ((cmf as u16) * 256 + flg as u16) % 31);
        [cmf, flg]
    }
}
