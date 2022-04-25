use glam::Vec2;
use crate::utils::*;

#[derive(Debug, Clone, Copy, Default)]
pub struct Particle {
    pub x: Vec2,
    pub v: Vec2,
    pub a: Vec2,
    pub t: i32,
    pub z: i32,
    pub m: f32,
    pub d: f32,
    pub p: f32,
    pub c: f32,
    pub u: f32,
}

impl Particle {
    #[must_use]
    pub fn new(x: f32, y: f32, t: i32) -> Self {
        Self {
            x: Vec2::new(x, y) / 1000.0,
            a: Vec2::ZERO,
            t: t,
            m: if t == 1 {MASS * 1.138} else {MASS},
            d: if t == 1 {1138.0} else {1000.0},
            u: 0.001,
            ..Particle::default()
        }
    }

    pub fn assign_v(&mut self, v: Vec2) {
        self.v = v;
    }

    #[must_use]
    pub fn position(&self) -> Vec2 {
        self.x * 1000.0
    }

    pub fn get_cell(&self) -> (i32,i32) {
        let pos = self.x / 1000.0;
        ((pos.x / (2.0 * H)) as i32, (pos.y / (2.0 * H)) as i32)
    }
}