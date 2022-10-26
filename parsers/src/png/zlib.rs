struct DeflateHeader {
    is_final: bool,
    block_type: DeflateBlockType,
}

#[repr(u8)]
enum DeflateBlockType {
    None = 0b0,
    HuffmanFixed = 0b01,
    HuffmanDynamic = 0b10,
}

/// Implement the DEFLATE compression algorithm
pub fn deflate(data: impl AsRef<[u8]>) -> Vec<u8> {
    let data = data.as_ref();

    // simple non-compressing output
    let mut output = vec![];
    let chunk_count = data.len() / u16::MAX as usize + 1;
    for (idx, chunk) in data.chunks(u16::MAX as _).enumerate() {
        output.push(
            DeflateHeader {
                is_final: idx == chunk_count,
                block_type: DeflateBlockType::None,
            }
            .as_byte(),
        );

        // push length
        let len = (chunk.len() as u16).to_le_bytes();
        output.push(len[0]);
        output.push(len[1]);

        // push ones complement of length
        let len = (!(chunk.len() as u16)).to_le_bytes();
        output.push(len[0]);
        output.push(len[1]);

        output.extend_from_slice(chunk);
    }

    output
}

impl DeflateHeader {
    fn as_byte(&self) -> u8 {
        (self.is_final as u8) << 7 | self.block_type as u8
    }
}

struct ZlibHeader {
    compression_info: u8,
    compression_level: u8,
}

pub fn zlib(data: impl AsRef<[u8]>) -> Vec<u8> {
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
        let dict = 0 & 0b1; // no preset dictionary
        let flg = level << 6 | dict << 5;

        let rem = ((cmf as u16) * 256 + flg as u16) % 31;
        let flg = flg | (31 - rem as u8);

        assert_eq!(0, ((cmf as u16) * 256 + flg as u16) % 31);
        [cmf, flg]
    }
}
