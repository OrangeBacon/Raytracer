use std::fmt::Display;

use once_cell::sync::Lazy;

mod zlib;
mod lzss;

/// A single png image
pub struct Png<T: AsRef<[u8]>> {
    width: u32,
    height: u32,
    bit_depth: u8,
    colour_type: ColourType,
    image_data: T,
}

/// Information to create a new png image.
pub struct PngCreateInfo<T: AsRef<[u8]>> {
    /// The raw image data
    pub image_data: T,

    /// The number of pixels wide the supplied image data is
    pub width: usize,

    /// The number of pixels tall the image is
    pub height: usize,

    /// The kind of data stored in the image data
    pub colour_type: InputColourType,

    /// The image data should be interpreted as u16 components, rather than u8.
    /// Assumes that the data is written to the buffer in big endian order.
    pub u16: bool,
}

/// Possible colour types describing the input data
pub enum InputColourType {
    /// One component per pixel
    Greyscale = 1,

    /// Two components per pixel
    GreyscaleAlpha = 2,

    /// Three components per pixel
    RGB = 3,

    /// Four components per pixel
    RGBA = 4,
}

#[derive(Debug)]
pub enum PngError {
    ImageDataLengthWrong(usize, usize),
    WidthTooBig(usize),
    HeightTooBig(usize),
}

impl std::error::Error for PngError {}

impl Display for PngError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PngError::ImageDataLengthWrong(expected, actual) => write!(
                f,
                "Expected {expected} bytes of image data, received {actual} bytes"
            ),
            PngError::WidthTooBig(width) => write!(
                f,
                "The image width, {width} is too large to be represented in a png image"
            ),
            PngError::HeightTooBig(height) => write!(
                f,
                "The image height, {height} is too large to be represented in a png image"
            ),
        }
    }
}

#[repr(u8)]
enum ColourType {
    Greyscale = 0,
    TrueColour = 2,
    GreyscaleAlpha = 4,
    TrueColourAlpha = 6,
}

impl<T: AsRef<[u8]>> Png<T> {
    /// Create a new png image
    pub fn new(info: PngCreateInfo<T>) -> Result<Self, PngError> {
        let data = info.image_data.as_ref();

        let component_size = if info.u16 { 2 } else { 1 };
        let pixel_size = component_size * info.colour_type as usize;
        let image_size = info.width * info.height * pixel_size;

        if image_size != data.len() {
            return Err(PngError::ImageDataLengthWrong(image_size, data.len()));
        }

        if info.width > u32::MAX as usize {
            return Err(PngError::WidthTooBig(info.width));
        }

        if info.height > u32::MAX as usize {
            return Err(PngError::HeightTooBig(info.height));
        }

        let colour_type = match info.colour_type {
            InputColourType::Greyscale => ColourType::Greyscale,
            InputColourType::GreyscaleAlpha => ColourType::GreyscaleAlpha,
            InputColourType::RGB => ColourType::TrueColour,
            InputColourType::RGBA => ColourType::TrueColourAlpha,
        };

        Ok(Self {
            width: info.width as _,
            height: info.height as _,
            bit_depth: if info.u16 { 16 } else { 8 },
            colour_type,
            image_data: info.image_data,
        })
    }

    /// Write the image to a data stream
    pub fn write(&self) -> Vec<u8> {
        let data = self.image_data.as_ref();

        let pixel_components = match self.colour_type {
            ColourType::Greyscale => 1,
            ColourType::TrueColour => 3,
            ColourType::GreyscaleAlpha => 2,
            ColourType::TrueColourAlpha => 4,
        };
        let scanline_width = pixel_components * (self.bit_depth as usize / 8) * self.width as usize;

        let mut filtered_data = Vec::with_capacity(data.len() + self.height as usize);

        // simple filtering: always use filter type 0 = no filter
        for scanline in data.chunks_exact(scanline_width) {
            filtered_data.push(0);
            filtered_data.extend(scanline);
        }
        debug_assert_eq!(filtered_data.len(), data.len() + self.height as usize);

        let compressed_data = zlib::zlib(&filtered_data);

        let mut output = vec![137, 80, 78, 71, 13, 10, 26, 10];
        write_chunk(b"IHDR", &mut output, |vec| {
            vec.extend(self.width.to_be_bytes());
            vec.extend(self.height.to_be_bytes());
            vec.push(self.bit_depth);
            vec.push(self.colour_type as _);
            vec.push(0); // compression method 0 == zlib/deflate
            vec.push(0); // filter method 0 == adaptive filtering with 5 types
            vec.push(0); // interlace 0 == no interlacing
        });
        write_chunk(b"IDAT", &mut output, |vec| vec.extend(compressed_data));

        // IEND chunk is always empty
        write_chunk(b"IEND", &mut output, |_| {});

        println!("{:?}", zlib::deflate(b"a"));

        output
    }
}

static CRC_TABLE: Lazy<[u32; 256]> = Lazy::new(|| {
    let mut result = [0; 256];

    for (idx, result) in result.iter_mut().enumerate() {
        let mut c = idx as u32;
        for _ in 0..8 {
            if c & 1 != 0 {
                c = 0xedb88320 ^ (c >> 1);
            } else {
                c >>= 1;
            }
        }
        *result = c;
    }

    result
});

fn crc(data: &[u8]) -> u32 {
    let mut crc = 0xffffffff;
    for elem in data {
        crc = CRC_TABLE[((crc ^ *elem as u32) & 0xff) as usize] ^ (crc >> 8);
    }

    crc ^ 0xffffffff
}

fn write_chunk<T>(name: &[u8; 4], data: &mut Vec<u8>, write: impl FnOnce(&mut Vec<u8>) -> T) -> T {
    let start_idx = data.len();
    data.extend(0u32.to_be_bytes()); // write 0 length to update later

    let name_idx = data.len();
    data.extend(name);

    let data_start_len = data.len();
    let result = write(data);
    let length = ((data.len() - data_start_len) as u32).to_be_bytes();

    // replace length field with actual length
    for (len, write) in length.iter().zip(&mut data[start_idx..]) {
        *write = *len;
    }

    let crc = crc(&data[name_idx..]).to_be_bytes();

    data.extend(crc);

    result
}
