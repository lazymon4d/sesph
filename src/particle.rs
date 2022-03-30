use crate::kernels::{
    Wpoly6,
    dWspiky,
    d2Wvisc,
    h
};

pub const d0: f64 = 1000.0f64;
pub const mass: f64 = (8.0f64/27.0f64) * d0 * h * h * h;
pub const kd: f64 = 0.7321f64;    // correct the value of state equation constant later
pub const cf: f64 = 0.7320f64;    // correct the value of coefficient of viscosity later
pub const ag: (f64,f64) = (0.0f64,-9.8130f64);    // acceleration due to gravity, balance si units everywhere else

fn dis(r1: &(f64,f64), r2: &(f64,f64)) -> f64 {
    let arg = (r1.0-r2.0)*(r1.0-r2.0) + (r1.1-r2.1)*(r1.1-r2.1);
    f64::sqrt(arg)
}

#[derive(Copy,Clone)]
pub enum Ptype {
    Fluid,
    Bound,
}

#[derive(Copy,Clone)]
pub struct Part {
    id    : u32,
    ptype : Ptype,
    f     : (f64,f64),
    r     : (f64,f64),
    v     : (f64,f64),
    d     : f64,
    p     : f64,
    m     : f64,              // mass is constant for all particles, usually taken as either d0*(h)^3 or d0*(2h/3)^3
}

impl Part {
    pub fn new(id: u32) -> Self {
        Part {
            id: id,
            ptype: Ptype::Fluid,
            r: (0f64,0f64),
            v: (0f64,0f64),
            f: (0.0f64,0.0f64),
            d: 0.0f64,
            p: 0.0f64,
            m: mass,
        }
    }

    pub fn part_type(&mut self, t: Ptype) -> &mut Self {
        self.ptype = t;
        self
    } 

    pub fn position(&mut self, r: (f64,f64)) -> &mut Self {
        self.r = r;
        self
    }

    pub fn velocity(&mut self, v: (f64,f64)) -> &mut Self {
        self.v = v;
        self
    }

    pub fn force(&mut self, f: (f64,f64)) -> &mut Self {
        self.f = f;
        self
    }

    pub fn density(&mut self, d: f64) -> &mut Self {
        self.d = d;
        self
    }

    pub fn pressure(&mut self, p: f64) -> &mut Self {
        self.p = p;
        self
    }

    pub fn mass(&mut self, m: f64) -> &mut Self {
        self.m = m;
        self
    }
}

impl Part {
    // interaction for pressure and density calculation
    pub fn interact1(&mut self, partj: &Part) {
        if self.id == partj.id {
            return;
        }
        self.d = self.d + partj.m * Wpoly6(dis(&(self.r), &(partj.r))); // could replace partj.m with const mass 
    }

    pub fn calc_p(&mut self) {
        self.p = kd * (self.d / d0 - 1.000f64);
    }

    // interaction for force calculation
    pub fn interact2(&mut self, partj: &Part) {
        let dr = dis(&(self.r),&(partj.r));
        let cx = (self.r.0-partj.r.0) / dr;
        let cy = (self.r.1-partj.r.1) / dr;

        self.f.0 = self.f.0 - cx * self.m * partj.m * (self.p / self.d.powi(2)) + (partj.p / partj.d.powi(2)) * dWspiky(dr);
        self.f.1 = self.f.1 - cy * self.m * partj.m * (self.p / self.d.powi(2)) + (partj.p / partj.d.powi(2)) * dWspiky(dr);

        self.f.0 = self.f.0 + cf * partj.m * ((partj.v.0 - self.v.0) / partj.d) * d2Wvisc(dr);
        self.f.1 = self.f.1 + cf * partj.m * ((partj.v.1 - self.v.1) / partj.d) * d2Wvisc(dr);
    }

    // position(r) and velocity(v) integration
    pub fn integrate(&mut self, step: &f64) {
        // include external forces (ag.0,ag.1)
        self.f.0 = self.f.0 + self.d * ag.0;
        self.f.1 = self.f.1 + self.d * ag.1;

        // velocity integration
        self.v.0 = self.v.0 + (self.f.0 / self.m) * (*step);
        self.v.1 = self.v.1 + (self.f.1 / self.m) * (*step);

        // position integration
        self.r.0 = self.r.0 + self.v.0 * (*step);
        self.r.1 = self.r.1 + self.v.1 * (*step); 
    }
}