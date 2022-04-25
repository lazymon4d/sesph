use mueller_sph_rs::{State,DT};

const DAM_PARTICLES: usize = 235000;
const BLOCK_PARTICLES: usize = 1000;
const MAX_PARTICLES: usize = DAM_PARTICLES + 25 * BLOCK_PARTICLES;
const M: usize = 1000000usize;

fn main() -> Result<(), String> {
    let mut simulation = State::<M>::new();
    simulation.init_fluids(DAM_PARTICLES);

    let mut time = 0.0;
    let num_iter = 50;
    let mut v = Vec::<f32>::new();

    for i in 0..num_iter {
        time += DT;
        let mut dia = 0.0;
        simulation.update(&mut dia);
        v.push(dia);
    }

    println!("{:?}",v);

    Ok(())
}