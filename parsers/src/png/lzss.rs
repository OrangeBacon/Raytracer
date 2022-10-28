use std::collections::VecDeque;

/// A buffer containing data written by the LZSS implementation.
pub struct LzssBuffer {
    /// The actual encoded output data, including the end of stream marker.
    pub data: Vec<LzssElement>,

    /// Number of occurrences of each value in the combined literal/length alphabet
    /// in the output data stream.
    pub literal_length: [usize; 286],

    /// Number of occurrences of each value in the distance alphabet.
    pub distance: [usize; 32],
}

/// Either a single byte or a back reference to previous bytes.
pub enum LzssElement {
    /// A single output bytes.
    Byte(u8),

    /// A marker signifying the end of the data stream.
    EndOfStream,

    /// A reference (distance, length) pair representing a distance behind and
    /// number of bytes to copy.
    Reference(usize, u16),
}

/// Compress the source text using LZSS encoding.  The block size is the size
/// of the window to examine when encoding, represented as log2(length)-8, so
/// a block length of 32k should be input as log2(32k)-8 = 15-8 = 7.
pub fn lzss(block_size: u8, data: &[u8]) -> LzssBuffer {
    let block_size = 1 << (block_size + 8);
    let _window: VecDeque<u8> = VecDeque::with_capacity(block_size);

    // currently don't compress, just forward the input data
    let mut buffer = LzssBuffer::new();
    for &byte in data {
        buffer.byte(byte);
    }

    buffer.end();

    buffer
}

impl LzssBuffer {
    /// Create a new empty buffer
    fn new() -> Self {
        Self {
            data: vec![],
            literal_length: [0; 286],
            distance: [0; 32],
        }
    }

    /// Add a single byte to the buffer
    fn byte(&mut self, byte: u8) {
        self.data.push(LzssElement::Byte(byte));
        self.literal_length[byte as usize] += 1;
    }

    fn end(&mut self) {
        self.data.push(LzssElement::EndOfStream);
        self.literal_length[256] += 1;
    }
}
