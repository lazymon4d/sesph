use std::f32::consts::PI;
use glam::Vec2;
use crate::utils::*;

pub fn W(rij: Vec2) -> f32 {
    let q: f32 = rij.length() / (H/1000.0);
    if q > 2.0 {
        0f32
    } else {
        7f32 / (4f32 * PI * HSQ/(1000.0*1000.0)) * (2.0 * q + 1.0) * (1.0 - q / 2.0).powi(4)
    }
}

pub fn dW(rij: Vec2) -> Vec2 {
    let q: f32 = rij.length() / (H/1000.0);
    if q > 2.0 {
        Vec2::ZERO
    } else {
        -35.0 / (4.0 * PI * H3/(1000000000.0)) * q * (1.0 - q / 2.0).powi(3) * rij.normalize()
    }
}