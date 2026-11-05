// (1) Laplace pressure of a 2D water droplet in cyclohexane: dp = sigma / R
// (2) Sessile water droplet on a wall with prescribed contact angle (through water)
const TPF = require('../solver.js');
const which = process.argv[2] || 'both';

function base(L, Hh, h) {
  const cfg = TPF.defaults();
  Object.assign(cfg, { Lx: L, Ly: Hh, h,
    ports: [{ kind: 'outlet', side: 'top', a: 0.45 * L, b: 0.55 * L, p: 0 }],
    contact: { bottom: [90, 90], top: [90, 90], left: [90, 90], right: [90, 90] } });
  cfg.numerics = Object.assign({}, cfg.numerics, { ramp: 0 });
  return cfg;
}

if (which !== 'sessile') {
  const L = 0.4e-3, h = 8e-6, R = 0.1e-3;
  const cfg = base(L, L, h);
  cfg.init = { type: 'droplet', inside: 'fluid1', cx: L / 2, cy: L / 2, R };
  const sim = new TPF.Simulation(cfg);
  const t0 = Date.now();
  while (sim.t < 4e-3) sim.step();
  const f = sim.fields(), { Nx, Ny } = sim.g;
  const Req = Math.sqrt(sim.fluid1Area() / Math.PI);
  // pressure inside (avg over r < 0.5R) and outside (avg over r > 1.6R)
  let pi = 0, ni = 0, po = 0, no = 0, umax = 0;
  for (let j = 0; j < Ny; j++) for (let i = 0; i < Nx; i++) {
    const k = i + j * Nx, r = Math.hypot((i + 0.5) * h - L / 2, (j + 0.5) * h - L / 2);
    if (r < 0.5 * Req) { pi += f.p[k]; ni++; } else if (r > 1.6 * Req) { po += f.p[k]; no++; }
    umax = Math.max(umax, f.speed[k]);
  }
  const dp = pi / ni - po / no, dpE = cfg.fluids.sigma / Req;
  console.log(`[Laplace] steps ${sim.nstep}, wall ${(Date.now() - t0) / 1000}s, R_eq ${(Req * 1e6).toFixed(1)} um`);
  console.log(`[Laplace] dp = ${dp.toFixed(1)} Pa, sigma/R = ${dpE.toFixed(1)} Pa, error ${((dp / dpE - 1) * 100).toFixed(2)} %, spurious |u|max = ${(umax * 1e3).toFixed(3)} mm/s`);
}

if (which !== 'laplace') {
  const L = 0.6e-3, Hh = 0.3e-3, h = 7.5e-6, R = 0.12e-3;
  for (const th of (process.argv[3] ? [Number(process.argv[3])] : [60, 120])) {
    const cfg = base(L, Hh, h);
    cfg.contact.bottom = [th, th];
    cfg.init = { type: 'droplet', inside: 'fluid1', cx: L / 2, cy: 0, R };
    const sim = new TPF.Simulation(cfg);
    const t0 = Date.now(); const tEnd = Number(process.argv[4] || 0.03);
    while (sim.t < tEnd) sim.step();
    const f = sim.fields(), { Nx, Ny } = sim.g;
    // base width from the phi=0 crossing in the first row; height from the phi=0 crossing on the centre column
    const row = (j) => { const xs = []; for (let i = 0; i < Nx - 1; i++) { const a = f.phi[i + j * Nx], b = f.phi[i + 1 + j * Nx]; if (a * b < 0) xs.push((i + 0.5 + a / (a - b)) * h); } return xs; };
    const xs = row(0); const Wb = xs[xs.length - 1] - xs[0];
    const ic = Nx >> 1; let Hd = 0;
    for (let j = 0; j < Ny - 1; j++) { const a = f.phi[ic + j * Nx], b = f.phi[ic + (j + 1) * Nx]; if (a * b < 0) { Hd = (j + 0.5 + a / (a - b)) * h; break; } }
    // cap geometry: tan(theta/2) = 2H/W ; base row is h/2 above the wall -> extrapolate width using circle fit
    const thetaRaw = 2 * Math.atan(2 * Hd / Wb) * 180 / Math.PI;
    // circle through the two base points (at y=h/2) and the apex (0,Hd) -> intersect with y=0
    const a0 = Wb / 2, y0 = h / 2;
    const yc = (Hd * Hd - a0 * a0 - y0 * y0) / (2 * (Hd - y0)), Rc = Hd - yc;
    const theta = Math.acos(-yc / Rc) * 180 / Math.PI;
    console.log(`[Sessile ${th}°] steps ${sim.nstep}, wall ${(Date.now() - t0) / 1000}s, t=${sim.t.toFixed(4)}s: measured theta = ${theta.toFixed(1)}° (raw cap estimate ${thetaRaw.toFixed(1)}°)`);
  }
}
