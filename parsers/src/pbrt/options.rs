use geometry::Float;

#[derive(Debug, Clone)]
pub struct RenderOptions {
    start_time: Float,
    end_time: Float,
}

impl RenderOptions {
    pub fn new() -> Self {
        Self {
            start_time: 0.0,
            end_time: 1.0,
        }
    }
}

impl Default for RenderOptions {
    fn default() -> Self {
        Self::new()
    }
}
