mod grid;

use std::f32::consts::PI;
use crate::grid::*;
use arrayvec::ArrayVec;
use glam::Vec2;
use lazy_static::lazy_static;
use log::info;
use rand::random;
use rayon::prelude::*;

// all units in SI
const REST_DENS: f32 = 1000.0;
const GAS_CONST: f32 = 320572000.0;
const H: f32 = 10.0;
const HSQ: f32 = H * H;
const H3: f32 = H * H * H;
const MASS: f32 = 0.1;
const VISC: f32 = 0.001;
const DT: f32 = 0.07;
const EPS: f32 = H;
const BOUND_DAMPING: f32 = -0.5;
const CONC: f32 = 0.7;
const MXCI: f32 = 0.0896;
const CAP: usize = 50000;

pub const WINDOW_WIDTH: u32 = 600;
pub const WINDOW_HEIGHT: u32 = 600;
pub const VIEW_WIDTH: f32 = 3.0 * WINDOW_WIDTH as f32;
pub const VIEW_HEIGHT: f32 = 3.0 * WINDOW_HEIGHT as f32;
pub const G: Vec2 = glam::const_vec2!([0.0, -9.81]);

fn W(rij: Vec2) -> f32 {
    let q: f32 = rij.length() / H;
    if q > 2.0 {
        0f32
    } else {
        7f32 / (4f32 * PI * HSQ) * (2.0 * q + 1.0) * (1.0 - q / 2.0).powi(4)
    }
}

fn dW(rij: Vec2) -> Vec2 {
    let q: f32 = rij.length() / H;
    if q > 2.0 {
        Vec2::ZERO
    } else {
        -35.0 / (4.0 * PI * H3) * q * (1.0 - q / 2.0).powi(3) * rij.normalize()
    }
}

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
            d: 1000.0,
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

#[derive(Debug, Default)]
pub struct State<const M: usize> {
    // Track both an initial and final state buffer so we don't need to copy in order
    // to access the current state during a mutable parallel iteration. Simulation
    // methods read from initial and modify final. At the end of the simulation step,
    // copy the final state into the initial (to be used in the next step).
    pub i: Box<ArrayVec<Particle, M>>,
    pub f: Box<ArrayVec<Particle, M>>,
    pub g: Box<Grid>,
}

impl<const M: usize> State<M> {
    #[must_use]
    pub fn new() -> Self {
        State {
            ..Default::default()
        }
    }

    pub fn mix_inflow(&mut self) {
        let pi = &mut self.i;
        let pf = &mut self.f;

        if pi.len() > CAP || pf.len() > CAP {
            return;
        }

        for i in 876..925 {
            for j in 1500..1550 {
                let jitter = random::<f32>();
                let mut part = Particle::new(i as f32 + jitter, j as f32, 1);
                part.assign_v(Vec2::new(0f32,-1.33f32));
                pi.push(part);pf.push(part);
            }
        }
    }

    pub fn init_dam_break(&mut self, dam_max_particles: usize) {
        let particles = &mut self.i;
        let mut y = EPS;
        'outer: while y < (VIEW_HEIGHT * 0.67 - EPS * 2.0) {
            y += H;
            let mut x = EPS;
            while x <= VIEW_WIDTH - EPS * 2.0 {
                x += H;
                if particles.len() < 10000usize {
                    let jitter = random::<f32>();
                    particles.push(Particle::new(x + jitter, y, 0));
                } else {
                    break 'outer;
                }
            }
        }
        self.f.clone_from(particles);
        info!("Initialized dam break with {} particles", particles.len());
    }

    pub fn init_block(&mut self, max_block_particles: usize) {
        todo!();
    }

    pub fn integrate(&mut self) {
        self.f.par_iter_mut().for_each(|p| {
            // p.v += DT * p.a;
            p.x += DT * p.v;

            let mut pos = p.x;
            pos *= 1000.0; 
            // enforce boundary conditions
            if pos.x - EPS < 0.0 {
                p.v.x *= BOUND_DAMPING;
                pos.x = EPS;
            }
            if pos.x + EPS > VIEW_WIDTH {
                p.v.x *= BOUND_DAMPING;
                pos.x = VIEW_WIDTH - EPS;
            }
            if pos.y - EPS < 0.0 {
                p.v.y *= BOUND_DAMPING;
                pos.y = EPS;
            }
            if pos.y + EPS > VIEW_HEIGHT {
                p.v.y *= BOUND_DAMPING;
                pos.y = VIEW_HEIGHT - EPS;
            }
            p.x = pos / 1000.0;
        });
    }

    pub fn compute_density_pressure(&mut self) {
        let si = &self.i;
        let sf = &mut self.f;
        let L = si.len();
        assert_eq!(L,sf.len());

        for i in 0usize..L {
            let mut d = 0.0;
            for j in self.g.get_nbrs(si[i].get_cell()) {
                if i == j {
                    continue;
                }
                let rij = si[i].x-si[j].x;
                d += si[j].m * W(rij);
            }
            sf[i].d = d;
            sf[i].p = GAS_CONST * ((si[i].d / REST_DENS).powf(7.0) - 1.0);
        } 

        // si.par_iter().zip(sf.par_iter_mut()).for_each(|(pi, pf)| {
        //     let mut rho = 0.0f32;
        //     pi
        //     for pj in si.iter() {
        //         let rij = pj.x - pi.x;
        //         rho += pj.m * W(rij);
        //     }
        //     pf.d = rho;
        //     pf.p = GAS_CONST * ((pi.d / REST_DENS).powf(7.0) - 1.0);
        // });
    }

    pub fn compute_concentration_viscosity(&mut self) {
        let si = &self.i;
        let sf = &mut self.f;
        let L = si.len();
        assert_eq!(L,sf.len());

        for i in 0usize..L {
            let (mut N,mut D) = (0.0,0.0);
            for j in self.g.get_nbrs(si[i].get_cell()) {
                if i == j {
                    continue;
                }
                let rij = si[i].x-si[j].x;
                if si[j].t == si[i].t {
                    N += si[j].m/si[j].d;
                }
                D += si[j].m/si[j].d;
            }
            sf[i].c = MXCI * N/D;;
            sf[i].u = VISC * (1.0 - MXCI * N / (D * CONC)).powf(- 2.5 * CONC);
        }

        // si.par_iter().zip(sf.par_iter_mut()).for_each(|(pi, pf)| {
        //     let (mut N, mut D) = (0.0,0.0);
        //     for pj in si.iter() {
        //         let rij = pj.x - pi.x;
        //         if rij.length_squared() > 4.0 * HSQ {
        //             continue;
        //         } else if pj.t == pi.t {
        //             N += pj.m/pj.d;
        //         }
        //         D += pj.m/pj.d;
        //     }
        //     pf.c = MXCI * N/D;
        //     pf.u = VISC * (1.0 - pf.c / CONC).powf(- 2.5 * CONC);
        // });
    }

    pub fn compute_accl(&mut self) {
        let si = &self.i;
        let sf = &mut self.f;

        let L = si.len();
        assert_eq!(L,sf.len());

        for i in 0usize..L {
            let mut accl = G;
            for j in self.g.get_nbrs(si[i].get_cell()) {
                if i == j {
                    continue;
                }                
                let rij = si[i].x-si[j].x;
                accl += -si[j].m * (si[j].p/si[j].d.powi(2) + si[i].p/si[i].d.powi(2)) * dW(rij);
                accl += si[j].m * 4.0 * si[i].u * si[j].u * rij.dot(dW(rij)) * (si[i].v-si[j].v) / (si[j].d * si[i].d * (si[i].u + si[j].u) * rij.length_squared());
            }
        }

        // si.par_iter().zip(sf.par_iter_mut()).enumerate().for_each(|(i, (pi, pf))| {
        //     let mut accl = G;
        //     for (j, pj) in si.iter().enumerate() {
        //         if i == j {
        //             continue;
        //         }
        //         let rij = pj.x - pi.x;
        //         let r = rij.length();
        //         accl += -pj.m * (pj.p/pj.d.powi(2) + pi.p/pi.d.powi(2)) * dW(rij);
        //         accl += pj.m * 4.0 * pi.u * pj.u * rij.dot(dW(rij)) * (pi.v-pj.v) / (pj.d * pi.d * (pi.u + pj.u) * rij.length_squared());
        //     }
        //     pf.a = accl;  
        // });
    }

    pub fn update_z(&mut self) {
        let L = self.i.len();
        for i in 0usize..L {
            self.g.update_access(i,self.i[i].z,self.i[i].get_cell());
        }
    }

    pub fn update(&mut self) {
        self.mix_inflow();
        self.compute_density_pressure();
        self.compute_concentration_viscosity();
        self.compute_accl();
        self.integrate();
        self.i.clone_from(&self.f);
    }
}
