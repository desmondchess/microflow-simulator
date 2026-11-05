// Plane Poiseuille flow: u(y) = 6U (y/H)(1 - y/H),  dp/dx = -12 mu U / H^2
const TPF = require('../solver.js');
const H = 0.2e-3, L = 1.0e-3, U = 0.01, mu = 1e-3, h = 10e-6;
const cfg = TPF.defaults();
Object.assign(cfg, {
  Lx: L, Ly: H, h,
  ports: [
    { kind: 'inlet', side: 'left', a: 0, b: H, U, phase: 'fluid2' },
    { kind: 'outlet', side: 'right', a: 0, b: H, p: 0 },
  ],
  fluids: { rho1: 1000, mu1: mu, rho2: 1000, mu2: mu, sigma: 0 },
  init: { type: 'front', fill: 'fluid2', front: -1 },
});
cfg.numerics = Object.assign({}, cfg.numerics, { slipRatio: 0, ramp: 0 });
const sim = new TPF.Simulation(cfg);
console.log('grid', sim.g.Nx, 'x', sim.g.Ny, 'dt', sim.dt.toExponential(2), 'limited by', sim.dtLimitedBy);
const tEnd = 0.08; const t0 = Date.now();
while (sim.t < tEnd) sim.step();
console.log('steps', sim.nstep, 'wall', (Date.now() - t0) / 1000, 's, maxDiv', sim.maxDiv().toExponential(2));
const f = sim.fields(), { Nx, Ny } = sim.g;
const i = Math.round(0.6 * Nx);
let err = 0, umax = 0;
for (let j = 0; j < Ny; j++) {
  const y = (j + 0.5) * h, ue = 6 * U * (y / H) * (1 - y / H);
  err = Math.max(err, Math.abs(f.u[i + j * Nx] - ue)); umax = Math.max(umax, ue);
}
const jm = Ny >> 1, x1 = (0.3 * Nx | 0), x2 = (0.8 * Nx | 0);
const dpdx = (f.p[x2 + jm * Nx] - f.p[x1 + jm * Nx]) / ((x2 - x1) * h);
const dpdxE = -12 * mu * U / (H * H);
console.log('max |u - u_exact| / u_max =', (err / umax * 100).toFixed(3), '%');
console.log('dp/dx =', dpdx.toFixed(2), 'Pa/m ; exact', dpdxE.toFixed(2), '; error', ((dpdx / dpdxE - 1) * 100).toFixed(3), '%');
