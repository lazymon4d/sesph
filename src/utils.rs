pub const REST_DENS1: f32 = 1000.0;
pub const REST_DENS2: f32 = 1138.0;
pub const GAS_CONST: f32 = 320572000.0;

pub const H: f32 = 10.0;
pub const HSQ: f32 = H * H;
pub const H3: f32 = H * H * H;
pub const MASS: f32 = 0.1;
pub const VISC: f32 = 0.001;

pub const EPS: f32 = H;
pub const BOUND_DAMPING: f32 = -0.5;
pub const CONC: f32 = 0.7;
pub const MXCI: f32 = 0.0896;
pub const CAP: usize = 5000;
pub const DT: f32 = 0.07;

pub const WINDOW_WIDTH: u32 = 600;
pub const WINDOW_HEIGHT: u32 = 600;
pub const VIEW_WIDTH: f32 = 3.0 * WINDOW_WIDTH as f32;
pub const VIEW_HEIGHT: f32 = 3.0 * WINDOW_HEIGHT as f32;

pub const G: glam::Vec2 = glam::const_vec2!([0.0, -9.81]);