mod kernels;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}

use crate::kernels::{Wpoly6,dWspiky,d2Wvisc,h};
use rayon::prelude::*;

const d0: f64 = 1000.0f64;
const mass: f64 = (8.0f64/27.0f64) * d0 * h * h * h;
const kd: f64 = 0.7321f64;    // correct the value of state equation constant later
const cf: f64 = 0.7320f64;    // correct the value of coefficient of viscosity later
const ag: (f64,f64) = (0.0f64,-9.8130f64);    // acceleration due to gravity, balance si units everywhere else

#[derive(Copy,Clone)]
enum Ptype {
    Fluid,
    Bound,
}

#[derive(Copy,Clone)]
struct Part {
    id    : u32,
    ptype : Ptype,
    r     : (f64,f64),
    v     : (f64,f64),
    f     : (f64,f64),
    d     : f64,
    p     : f64,
    m     : f64,              // mass is constant for all particles, usually taken as either d0*(h)^3 or d0*(2h/3)^3
}

impl Part {
    // interaction for pressure and density calculation
    fn interact1(&mut self, partj: &Part) {
        if self.id == partj.id {
            return;
        }
        self.d = self.d + partj.m * Wpoly6(dis(&(self.r), &(partj.r))); // could replace partj.m with const mass 
    }

    // interaction for force calculation
    fn interact2(&mut self, partj: &Part) {
        let dr = dis(&(self.r),&(partj.r));
        let cx = (self.r.0-partj.r.0) / dr;
        let cy = (self.r.1-partj.r.1) / dr;

        self.f.0 = self.f.0 - cx * self.m * partj.m * (self.p / self.d.powi(2)) + (partj.p / partj.d.powi(2)) * dWspiky(dr);
        self.f.1 = self.f.1 - cy * self.m * partj.m * (self.p / self.d.powi(2)) + (partj.p / partj.d.powi(2)) * dWspiky(dr);

        self.f.0 = self.f.0 + cf * partj.m * ((partj.v.0 - self.v.0) / partj.d) * d2Wvisc(dr);
        self.f.1 = self.f.1 + cf * partj.m * ((partj.v.1 - self.v.1) / partj.d) * d2Wvisc(dr);
    }

    // position(r) and velocity(v) integration
    fn integrate(&mut self, step: &f64) {
        // include external forces (ag.0,ag.1)
        self.f.0 = self.f.0 + self.d * ag.0;
        self.f.1 = self.f.1 + self.d * ag.1;

        // velocity integration
        self.v.0 = self.v.0 + (self.f.0 / self.m) * (*step);
        self.v.1 = self.v.1 + (self.f.1 / self.m) * (*step);

        // // position integration
        self.r.0 = self.r.0 + self.v.0 * (*step);
        self.r.1 = self.r.1 + self.v.1 * (*step);        
    }
}

fn dis(r1: &(f64,f64), r2: &(f64,f64)) -> f64 {
    let arg = (r1.0-r2.0)*(r1.0-r2.0) + (r1.1-r2.1)*(r1.1-r2.1);
    f64::sqrt(arg)
}

fn get_case(vpart: &mut Vec<Part>) {
    let n: usize = 1007;
    let rk: f64 = 170000f64;
    let vk: f64 = 110000f64;
    let mk: f64 = 200000f64;
    for i in 1..=n {
        vpart.push(Part {
            id: i as u32,
            ptype: Ptype::Fluid,
            r: ({(i as f64-(n as f64/2 as f64)) / rk},{((n as f64/2 as f64)-i as f64) / rk}),
            v: ({(i as f64-(n as f64/2 as f64)) / vk},{((n as f64/2 as f64)-i as f64) / vk}),
            f: (0.0f64,0.0f64),
            d: 0.0f64,
            p: 0.0f64,
            m: mass,
        });
    }
}

fn sph(vpart: &mut Vec<Part>, timesteps: &Vec<f64>) {
    let n: usize = vpart.len(); // balance si units everywhere    

    let vi: Vec<usize> = vec![0;n];
    let mut vpart2: Vec<Part> = vpart.to_vec();

    for step in timesteps.iter() {
        vpart.par_iter_mut().for_each(|i| *i = Part {
            f: (0.0f64,0.0f64),
            d: 0.0f64,
            p: 0.0f64,
            ..*i
        });

        vpart.par_iter_mut().for_each(|i| {
            for j in vpart2.iter() {
                i.interact1(j);
            }
            i.p = kd * (i.d / d0 - 1.000f64);
            for j in vpart2.iter() {
                i.interact2(j);
            }
            i.integrate(step);
        });

        vpart2 = vpart.to_vec();
    }
}

fn main() {
    let timesteps: Vec<f64> = vec![0.1;100];
    let mut vpart: Vec<Part> = Vec::<Part>::new(); 
    get_case(&mut vpart);
    sph(&mut vpart, &timesteps);
    println!("done!!");
}