mod grid;
mod kernel;
mod utils;
mod particle;

use crate::grid::*;
use crate::kernel::*;
use crate::particle::*;
pub use crate::utils::*;

use arrayvec::ArrayVec;
use glam::Vec2;
use glam::Vec4;
use log::info;
use rand::random;
use rayon::prelude::*;

pub const DT: f32 = 0.0007;
pub const C1: Vec4 = glam::const_vec4!([0.149,0.141,0.912,1.000]);
pub const C2: Vec4 = glam::const_vec4!([1.000,0.833,0.224,1.000]);

#[derive(Debug, Default)]
pub struct State<const M: usize> {
    pub i: Box<ArrayVec<Particle, M>>,
    pub f: Box<ArrayVec<Particle, M>>,
    pub g: Box<Grid>,
}

impl<const M: usize> State<M> {
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

        for i in 896..906 {
            for j in 595..606 {
                let jitter = random::<f32>();
                let mut part = Particle::new(i as f32 + jitter, j as f32, 1);
                part.clr = C2;
                part.assign_v(Vec2::new(0f32,-1.33f32));
                pi.push(part);pf.push(part);
            }
        }
    }
    pub fn init_fluids(&mut self, dam_max_particles: usize) {
        let particles = &mut self.i;
        let mut y = EPS;
        'outer: while y < (VIEW_HEIGHT * 0.67 - EPS * 2.0) {
            y += H;
            let mut x = EPS;
            while x <= VIEW_WIDTH - EPS * 2.0 {
                x += H;
                if particles.len() < 2000usize {
                    let jitter = random::<f32>();
                    let mut part = Particle::new(x + jitter, y, 0);
                    part.clr = C1;
                    particles.push(part);
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
            p.v += DT * p.a;
            p.x += DT * p.v;
            let mut pos = p.x;
            pos *= 1000.0; 

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
            let d0 = if sf[i].t == 1 {REST_DENS2} else {REST_DENS1};
            sf[i].p = GAS_CONST * ((si[i].d / d0).powf(7.0) - 1.0);
        }
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
            sf[i].c = MXCI * N/D;
            sf[i].u = VISC * (1.0 - MXCI * N / (D * CONC)).powf(- 2.5 * CONC);
        }
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
                // accl += - si[j].m * (si[j].p/si[j].d.powi(2) + si[i].p/si[i].d.powi(2)) * dW(rij);
                // accl += si[j].m * 4.0 * si[i].u * si[j].u * rij.dot(dW(rij)) * (si[i].v-si[j].v) / (si[j].d * si[i].d * (si[i].u + si[j].u) * rij.length_squared());
            }
            sf[i].a = accl;
        }
    }
    pub fn update_z(&mut self) {
        let L = self.i.len();
        for i in 0usize..L {
            self.g.update_access(i,self.i[i].z,self.i[i].get_cell());
        }
    }
    pub fn update(&mut self) {
        self.mix_inflow();
        self.update_z();
        self.compute_density_pressure();
        self.compute_concentration_viscosity();
        self.compute_accl();
        self.integrate();
        self.i.clone_from(&self.f);
    }
}