use arrayvec::ArrayVec;
use std::collections::{HashMap,HashSet};
use crate::*;

#[derive(Debug, Default)]
pub struct Grid {
	zcells	: HashMap::<(i32,i32),i32>,
	zaccess	: HashMap::<i32,HashSet<usize>>, 
	idomain : (i32,i32),
}

impl Grid {
	pub fn z_index(&mut self, x: (i32,i32)) -> i32 {
		match self.zcells.get(&x) {
			Some(v) => *v,
			None => {
				self.idomain = ((VIEW_WIDTH/(2.0*H)) as i32, (VIEW_HEIGHT/(2.0*H)) as i32);
				let v = self.z_hash(x);
				self.zcells.insert(x,*&v);
				v
			},
		}
	}

	pub fn get_nbrs(&mut self, c: (i32,i32)) -> ArrayVec::<usize,100> {
		let mut nbrs = ArrayVec::<usize,100>::new();
		for i in self.nbr_zi(c).iter() {
			match self.zaccess.get(i) {
				Some(set) => {
					for j in set.iter() {
						if nbrs.len() < 95 {
							nbrs.push(*j);
						}
					}
				},
				None => {}
			}
		}
		nbrs
	}

	pub fn update_access(&mut self, idx: usize, z: i32, c: (i32,i32)) {
		match self.zaccess.get_mut(&z) {
			Some(set) => {
				set.remove(&idx);
			},
			None => {}
		}
		let z = self.z_hash(c);
		match self.zaccess.get_mut(&z) {
			Some(set) => {
				set.insert(idx);
			},
			None => {
				let mut set = self.zaccess.insert(z,HashSet::<usize>::new()).unwrap();
				set.insert(idx);
			},
		}
	}
	
	fn z_hash(&self, C: (i32,i32)) -> i32 {
		/*
			Z index:

			the grid is filled with a space filling z-curve, and indices are numbered accordingly,
			here, this serves no special purpose, but this technique, in more sophisticated implementations, 
			can be used to speed up the computation of neighbours.
		*/

		let mut M = self.idomain;
		let mut C = C;
		let mut z = 0;
		while M.0 > 0 && M.1 > 0 {
		    z <<= 2;
		    z += 2 * (C.1/M.1) + (C.0/M.0);
		    C.1 %= M.1; C.0 %= M.0;
		    M.1 >>= 1; M.0 >>= 1;
		}
		z
	}

	fn nbr_zi(&mut self, c: (i32,i32)) -> ArrayVec::<i32,9> {
        let affinity = [
            (-1,-1),(-1,0),(-1,1),
            (0,-1),(0,0),(0,1),
            (1,-1),(1,0),(1,1)
        ];
        let mut nz = ArrayVec::<i32,9>::new();

        for (dx, dy) in affinity {
            let nzi = self.z_index((c.0+dx, c.1+dy));
            nz.push(nzi);
        }
        return nz;
    }
}