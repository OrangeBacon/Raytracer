use std::{fmt::Debug, sync::Arc};

use glam::{DVec2, DVec3};

pub trait Texture: Send + Sync + Debug {
    fn value(&self, uv: DVec2, point: DVec3) -> DVec3;
}

#[derive(Debug)]
pub struct SolidColour {
    pub colour: DVec3,
}

impl Texture for SolidColour {
    fn value(&self, _: DVec2, _: DVec3) -> DVec3 {
        self.colour
    }
}

#[derive(Debug)]
pub struct Checker {
    pub odd: Arc<dyn Texture>,
    pub even: Arc<dyn Texture>,
}

impl Texture for Checker {
    fn value(&self, uv: DVec2, point: DVec3) -> DVec3 {
        let sines = point * 10.0;
        let sines = sines.x.sin() * sines.y.sin() * sines.z.sin();

        if sines < 0.0 {
            self.odd.value(uv, point)
        } else {
            self.even.value(uv, point)
        }
    }
}
