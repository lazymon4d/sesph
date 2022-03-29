#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}

#[derive(Copy,Clone)]
struct Part {
    r: (f64,f64),
    v: (f64,f64),
    m: f64
}

fn Wpoly6(r: f64) -> f64 {
    let h = 0.0000083189f64;
    r*h
}

fn dWspiky(r: f64) -> f64 {
    let h = 0.0000083189f64;
    r*h
}

fn d2Wvisc(r: f64) -> f64 {
    let h = 0.0000083189f64;
    r*h
}

fn dis(i: &(f64,f64), j: &(f64,f64)) -> f64 {
    let arg = (i.0-j.0)*(i.0-j.0) + (i.1-j.1)*(i.1-j.1);
    f64::sqrt(arg)
}

fn get_case(n: u32) -> Vec<Part> {
    let mut vpart: Vec<Part> = Vec::new();
    let rk: f64 = 170000f64;
    let vk: f64 = 110000f64;
    let mk: f64 = 200000f64;
    for i in 1..=n {
        vpart.push(Part {
            r: ({(i as f64-(n as f64/2 as f64)) / rk},{((n as f64/2 as f64)-i as f64) / rk}),
            v: ({(i as f64-(n as f64/2 as f64)) / vk},{((n as f64/2 as f64)-i as f64) / vk}),
            m: {i as f64 / mk},
        });
    }
    vpart
}

fn sph(vpart: Vec<Part>, timesteps: Vec<f64>) {
    let n: usize = vpart.len();
    let d0: f64 = 1000.0f64;    // balance si units everywhere    
    let kd: f64 = 0.7321f64;    // correct the value of state equation constant later
    let cf: f64 = 0.7320f64;    // correct the value of coefficient of viscosity later
    let ag: f64 = 9.8130f64;    // acceleration due to gravity, balance si units everywhere else

    let mut rx: Vec<f64> = vec![0.0f64;n];
    let mut ry: Vec<f64> = vec![0.0f64;n];
    let mut vx: Vec<f64> = vec![0.0f64;n];
    let mut vy: Vec<f64> = vec![0.0f64;n];
    let mut mi: Vec<f64> = vec![0.0f64;n];

    for i in 0..n {
        rx[i] = vpart[i].r.0;
        ry[i] = vpart[i].r.1;
        vx[i] = vpart[i].v.0;
        vy[i] = vpart[i].v.1;
        mi[i] = vpart[i].m;
    }

    let mut pi: Vec<f64> = vec![0.0f64;n];
    let mut di: Vec<f64> = vec![0.0f64;n];
    let mut fx: Vec<f64> = vec![0.0f64;n];
    let mut fy: Vec<f64> = vec![0.0f64;n];

    for step in timesteps.iter() {
        pi = vec![0.0f64;n];
        di = vec![0.0f64;n];
        fx = vec![0.0f64;n];
        fy = vec![0.0f64;n];
        
        for i in 0..n {
            for j in 0..n {
                if j==i {
                    continue;
                }
                let xij = (rx[i],rx[j]);
                let yij = (ry[i],ry[j]);
                di[i] = di[i] + mi[j]*Wpoly6(dis(&xij,&yij));
            }
            pi[i] = kd * (di[i] / d0 - 1.000f64);

            for j in 0..n {
                if j==i {
                    continue;
                }
                let xij = (rx[i],rx[j]);
                let yij = (ry[i],ry[j]);

                let cx = (rx[i]-rx[j])/dis(&xij,&yij);
                let cy = (ry[i]-ry[j])/dis(&xij,&yij);

                // pressure force calculation
                fx[i] = fx[i] - cx*mi[i]*mi[j]*(pi[i]/(di[i]*di[i]) + pi[j]/(di[j]*di[j]))*dWspiky(dis(&xij,&yij));
                fy[i] = fy[i] - cy*mi[i]*mi[j]*(pi[i]/(di[i]*di[i]) + pi[j]/(di[j]*di[j]))*dWspiky(dis(&xij,&yij));

                // visc force calculation
                fx[i] = fx[i] + cf*mi[j]*((vx[j]-vx[i])/di[j])*d2Wvisc(dis(&xij,&yij));
                fy[i] = fy[i] + cf*mi[j]*((vy[j]-vy[i])/di[j])*d2Wvisc(dis(&xij,&yij));
            }

            // ext force calculation
            fy[i] = fy[i] - di[i]*ag;

            // velocity integration
            vx[i] = vx[i] + fx[i]/mi[i] * (*step);
            vy[i] = vy[i] + fy[i]/mi[i] * (*step);

            // position integration
            rx[i] = rx[i] + vx[i] * (*step);
            ry[i] = ry[i] + vy[i] * (*step);
        }
    }
}

fn main() {
    let timesteps: Vec<f64> = vec![0.1;100];
    sph(get_case(200), timesteps);
    println!("done!!");
}
