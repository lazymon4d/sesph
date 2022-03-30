mod kernels;
mod particle;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}

use rayon::prelude::*;
use crate::particle::*;

fn get_case(vpart: &mut Vec<Part>) {
    let n: usize = 1007;
    let rk: f64 = 170000f64;
    let vk: f64 = 110000f64;
    let mk: f64 = 200000f64;
    for i in 1..=n {
        vpart.push(
            *Part::new(i as u32)
            .part_type(Ptype::Fluid)
            .position(({(i as f64-(n as f64/2 as f64)) / rk},{((n as f64/2 as f64)-i as f64) / rk}))
            .velocity(({(i as f64-(n as f64/2 as f64)) / vk},{((n as f64/2 as f64)-i as f64) / vk}))
        );
    }
}

fn sph(vpart: &mut Vec<Part>, timesteps: &Vec<f64>) {
    let n: usize = vpart.len(); // balance si units everywhere    

    let vi: Vec<usize> = vec![0;n];
    let mut vpart2: Vec<Part> = vpart.to_vec();

    for step in timesteps.iter() {
        vpart.par_iter_mut().for_each(|i| {
            i.force((0f64,0f64))
            .pressure(0f64)
            .density(0f64);
        }); 

        vpart.par_iter_mut().for_each(|i| {
            for j in vpart2.iter() {
                i.interact1(j);
            }
            i.calc_p();
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