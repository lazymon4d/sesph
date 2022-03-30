// Kernels

pub const h : f64 = 0.0000083189f64;
const h2: f64 = h*h;
const h3: f64 = h*h*h;
const pi: f64 = 3.1415926536f64;

const kpoly6: f64 = 315.0f64 / (64.0f64 * pi * (h3*h3*h3));
const kspiky: f64 = -45.0f64 / (pi * (h3*h3));
const kviscs: f64 = 45.0f64 / (pi * (h3*h3));


// Poly6 Kernel -> to be used in density calculation
pub fn Wpoly6(r: f64) -> f64 {
	if r > h {
		0.0f64
	} else {
		(kpoly6) * (h2 - r.powi(2)).powi(3) 
	}
}

// Del(Spiky Kernel) (is a vector) -> returns magnitude of vector
// Used to calculate pressure forces
pub fn dWspiky(r: f64) -> f64 {
	if r > h {
		0.0f64
	} else {
		(kspiky) * (h - r).powi(2)
	}
}

// Del2(Viscosity Kernel)
// Used to calculate Viscous forces
pub fn d2Wvisc(r: f64) -> f64 {
	if r > h {
		0.0f64
	} else {
		(kviscs) * (h - r)
	}
}