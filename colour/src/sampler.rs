use geometry::{ConstZero, Number, Point2, Point2i};

/// A n object to get locations to send rays within a pixel
pub trait Sampler<T: Number> {
    /// Start sampling at a given pixel location.  Must call the free function
    /// start pixel in the implementation of this method.
    fn start_pixel(&mut self, point: Point2i) {
        start_pixel(self, point)
    }

    /// Get the next sample value for the next dimension of the sample vector
    fn one(&mut self) -> T;

    /// Get the next two sample values for the next dimensions of the sample vector
    fn two(&mut self) -> Point2<T>;

    /// Initialise a camera sample
    fn camera_sample(&mut self, _point: Point2i) {
        todo!()
    }

    /// Adjust the number of samples to one more favourable to be passed to the
    /// request_1D and 2D_array functions
    fn round_count(&mut self, count: usize) -> usize {
        count
    }

    /// Start a new sample from the sampler, returns false if the number of samples
    /// requested when the sampler was created has already been created. Implementations
    /// must call the free function start_next_sample.
    fn start_next_sample(&mut self) -> bool {
        start_next_sample(self)
    }

    /// Set the index of the next sample to generate.  Implementations must call
    /// the free function set_sample_number.
    fn set_sample_number(&mut self, sample_num: usize) -> bool {
        set_sample_number(self, sample_num)
    }

    /// Clone this sampler, but set the random seed to the given seed
    fn clone_seed(&mut self, seed: u64) -> Self;

    /// Get the generic sampler data
    fn sample_data(&mut self) -> &mut SamplerData<T>;
}

#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SamplerData<T: Number> {
    /// How many samples should be taken per pixel
    pub samples_per_pixel: usize,

    /// Location of the current pixel
    pub(super) current_pixel: Point2i,

    /// Number of samples generated in the current pixel
    pub(super) current_pixel_sample_index: usize,

    /// Sizes of 1D sample array to generate
    pub(super) samples_one_array_sizes: Vec<usize>,

    /// Sizes of 2D sample array to generate
    pub(super) samples_two_array_sizes: Vec<usize>,

    /// The 1D samples
    pub(super) sample_array_one: Vec<Vec<T>>,

    /// The 2D samples
    pub(super) sample_array_two: Vec<Vec<Point2<T>>>,

    /// offset into the 1D sample array to use
    array_one_offset: usize,

    /// offset into the 2D sample array to use
    array_two_offset: usize,
}

impl<T: Number> SamplerData<T> {
    /// Create new sampler data
    pub fn new(samples_per_pixel: usize) -> Self {
        Self {
            samples_per_pixel,
            ..Default::default()
        }
    }

    /// Specify that a given 1D array size should be used.
    pub fn request_one_array(&mut self, count: usize) {
        self.samples_one_array_sizes.push(count);
        let vec = vec![T::ZERO; count * self.samples_per_pixel];
        self.sample_array_one.push(vec);
    }

    /// Specify that a given 2D array size should be used.
    pub fn request_two_array(&mut self, count: usize) {
        self.samples_two_array_sizes.push(count);
        let vec = vec![Point2::ZERO; count * self.samples_per_pixel];
        self.sample_array_two.push(vec);
    }

    /// Get an array of 1D samples of length specified from request_one_array
    pub fn one_array(&mut self, count: usize) -> &[T] {
        if self.array_one_offset == self.sample_array_one.len() {
            return &[];
        }

        self.array_one_offset += 1;
        let res = &self.sample_array_one[self.array_one_offset - 1]
            [self.current_pixel_sample_index * count..];
        res
    }

    /// Get an array of 2D samples of length specified from request_two_array
    pub fn two_array(&mut self, count: usize) -> &[Point2<T>] {
        if self.array_two_offset == self.sample_array_two.len() {
            return &[];
        }

        self.array_two_offset += 1;
        let res = &self.sample_array_two[self.array_two_offset - 1]
            [self.current_pixel_sample_index * count..];
        res
    }
}

/// Start sampling at a given pixel location.
pub fn start_pixel<T: Sampler<F> + ?Sized, F: Number>(this: &mut T, point: Point2i) {
    let this = this.sample_data();
    this.current_pixel = point;
    this.current_pixel_sample_index = 0;
    this.array_one_offset = 0;
    this.array_two_offset = 0;
}

/// Start a new sample from the sampler, returns false if the number of samples
/// requested when the sampler was created has already been created.
pub fn start_next_sample<T: Sampler<F> + ?Sized, F: Number>(this: &mut T) -> bool {
    let this = this.sample_data();
    this.array_one_offset = 0;
    this.array_two_offset = 0;
    this.current_pixel_sample_index += 1;
    this.current_pixel_sample_index < this.samples_per_pixel
}

/// Set the index of the next sample to generate.  
pub fn set_sample_number<T: Sampler<F> + ?Sized, F: Number>(
    this: &mut T,
    sample_num: usize,
) -> bool {
    let this = this.sample_data();
    this.array_one_offset = 0;
    this.array_two_offset = 0;
    this.current_pixel_sample_index = sample_num;
    this.current_pixel_sample_index < this.samples_per_pixel
}
