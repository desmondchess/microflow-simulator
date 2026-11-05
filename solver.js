/*
 * TPF — 2D two-phase laminar flow solver
 * Incompressible Navier–Stokes (variable rho, mu) coupled to a Cahn–Hilliard
 * phase field, written in the same form COMSOL's "Two-Phase Flow, Phase Field"
 * interface uses, so results can be compared term by term.
 *
 *   phi = -1 : Fluid 1 (water)        phi = +1 : Fluid 2 (cyclohexane)
 *
 *   dphi/dt + u.grad(phi) = div( (gamma*lambda/eps^2) grad(psi) )
 *   psi = -eps^2 lap(phi) + (phi^2 - 1) phi
 *   lambda = 3 eps sigma / sqrt(8),   gamma = chi eps^2
 *   rho du/dt + rho (u.grad)u = -grad p + div[mu (grad u + grad u^T)] + G grad(phi) + rho g
 *   G = lambda psi / eps^2,   div u = 0
 *
 * Discretisation: uniform staggered (MAC) grid; explicit Euler; MUSCL/van Leer
 * phase advection; central momentum advection; constant-coefficient pressure
 * projection (Dodd & Ferrante 2014) solved with a banded Cholesky factorisation
 * computed once. Geometry is a cell mask, with inlet/outlet "ports" tagged on the
 * boundary, so imported shapes plug into the same solver.
 *
 * Works in the browser (window.TPF) and in Node (require('./solver.js')).
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.TPF = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const SOLID = 0, FLUID = 1, PORT = 2;
  const WALL = 0, INLET = 1, OUTLET = 2;
  const SIDE = { left: 0, right: 1, bottom: 2, top: 3 };
  const SIDES = ['left', 'right', 'bottom', 'top'];
  const DEG = Math.PI / 180;

  function defaults() {
    return {
      Lx: 2.0e-3, Ly: 0.5e-3, h: 12.5e-6,
      ports: [
        { kind: 'inlet', side: 'left', a: 0, b: 0.5e-3, U: 0.01, phase: 'fluid1', waterFrac: 0.5 },
        { kind: 'outlet', side: 'right', a: 0.3e-3, b: 0.5e-3, p: 0 },
        { kind: 'outlet', side: 'right', a: 0.0, b: 0.2e-3, p: 0 },
      ],
      fluids: { name1: 'Water', rho1: 998.2, mu1: 1.002e-3, name2: 'Cyclohexane', rho2: 778.5, mu2: 0.98e-3, sigma: 0.050 },
      // contact angle measured through Fluid 1 (water), degrees; [start, end] varies linearly along the side
      contact: { bottom: [45, 45], top: [135, 135], left: [90, 90], right: [90, 90], solid: 90 },
      init: { type: 'front', fill: 'fluid2', front: 0.15e-3 },
      numerics: { epsRatio: 1.0, chi: 1.0, slipRatio: 1.0, safety: 0.5, ramp: 2e-3, depth: 0, gx: 0, gy: 0 },
      tEnd: 0.25,
    };
  }

  // ---------------------------------------------------------------- geometry
  function buildGeometry(cfg) {
    const h = cfg.h;
    const Nx = Math.max(4, Math.round(cfg.Lx / h));
    const Ny = Math.max(4, Math.round(cfg.Ly / h));
    const W = Nx + 2, H = Ny + 2, N = W * H;
    const type = new Uint8Array(N);          // SOLID everywhere by default
    const port = new Int16Array(N).fill(-1);
    const side = new Int8Array(N).fill(-1);  // which outer side a ghost cell belongs to
    const s = new Float64Array(N);           // coordinate along that side
    for (let j = 1; j <= Ny; j++) for (let i = 1; i <= Nx; i++) {
      const c = i + j * W;
      type[c] = FLUID;
      if (cfg.mask && !cfg.mask((i - 0.5) * h, (j - 0.5) * h)) type[c] = SOLID;
    }
    for (let j = 1; j <= Ny; j++) {
      side[j * W] = 0; s[j * W] = (j - 0.5) * h;
      side[Nx + 1 + j * W] = 1; s[Nx + 1 + j * W] = (j - 0.5) * h;
    }
    for (let i = 1; i <= Nx; i++) {
      side[i] = 2; s[i] = (i - 0.5) * h;
      side[i + (Ny + 1) * W] = 3; s[i + (Ny + 1) * W] = (i - 0.5) * h;
    }
    const ports = (cfg.ports || []).map((p, id) => Object.assign({}, p, { id, cells: [] }));
    for (const p of ports) {
      const sd = SIDE[p.side];
      if (sd === undefined) throw new Error('Unknown side "' + p.side + '"');
      for (let c = 0; c < N; c++) {
        if (side[c] !== sd || type[c] === PORT) continue;
        if (s[c] > p.a && s[c] < p.b) { type[c] = PORT; port[c] = p.id; p.cells.push(c); }
      }
      if (!p.cells.length) throw new Error('Port ' + (p.id + 1) + ' (' + p.kind + ', ' + p.side + ') covers no boundary cells — check its extent.');
    }
    return { Nx, Ny, W, H, N, h, Lx: Nx * h, Ly: Ny * h, type, port, side, s, ports };
  }

  const vanLeer = (a, b) => (a * b > 0 ? (2 * a * b) / (a + b) : 0);

  // -------------------------------------------------------------- simulation
  class Simulation {
    constructor(cfg) {
      this.cfg = cfg;
      const g = (this.g = buildGeometry(cfg));
      const { N, h } = g;
      const fl = cfg.fluids, nm = cfg.numerics;
      this.r1 = fl.rho1; this.r2 = fl.rho2; this.m1 = fl.mu1; this.m2 = fl.mu2; this.sigma = fl.sigma;
      this.eps = nm.epsRatio * h;
      this.lam = (3 * this.eps * this.sigma) / (2 * Math.SQRT2); // mixing energy density lambda
      this.Mc = nm.chi * this.lam;                                // gamma*lambda/eps^2 = chi*lambda
      this.Gfac = this.lam / (this.eps * this.eps);
      this.rho0 = Math.min(this.r1, this.r2);
      const br = nm.slipRatio;                                    // slip length / h
      this.slipFac = (br - 0.5) / (br + 0.5);
      this.depth = nm.depth || 0;
      this.gx = nm.gx || 0; this.gy = nm.gy || 0;
      this.ramp = nm.ramp || 0;
      const F = () => new Float64Array(N);
      Object.assign(this, {
        phi: F(), phiN: F(), psi: F(), G: F(), rho: F(), mu: F(), p: F(), pOld: F(), ph: F(),
        u: F(), v: F(), us: F(), vs: F(), txy: F(), lap: F(), rate: F(), acc: F(),
      });
      this._classify();
      this._buildPoisson();
      this._computeDt(nm.safety);
      this.t = 0; this.nstep = 0;
      this.portStats = g.ports.map(() => ({ Q: 0, Qw: 0, V: 0, Vw: 0 }));
      this.initPhi(cfg.init || { type: 'front', fill: 'fluid2', front: 0 });
    }

    // classify every staggered face: active / boundary (wall, inlet, outlet) / ghost
    _classify() {
      const g = this.g, { W, N, h, type, port, side, s, ports } = g;
      const ct = this.cfg.contact || {};
      const fluid = [];
      for (let c = 0; c < N; c++) if (type[c] === FLUID) fluid.push(c);
      this.fc = Int32Array.from(fluid);
      const sideLen = [g.Ly, g.Ly, g.Lx, g.Lx];
      const thetaAt = (oc) => {
        const sd = side[oc];
        if (sd < 0) return ct.solid !== undefined ? ct.solid : 90;
        const pr = ct[SIDES[sd]] || [90, 90];
        return pr[0] + (pr[1] - pr[0]) * (s[oc] / sideLen[sd]);
      };
      const inletInfo = {};
      for (const p of ports) if (p.kind === 'inlet') {
        let lo = Infinity, hi = -Infinity;
        for (const c of p.cells) { lo = Math.min(lo, s[c] - h / 2); hi = Math.max(hi, s[c] + h / 2); }
        const prof = new Map(); let sum = 0;
        for (const c of p.cells) { const z = (2 * (s[c] - lo)) / (hi - lo) - 1; const val = 1 - z * z; prof.set(c, val); sum += val; }
        const scale = (p.U * p.cells.length) / sum;              // discrete mean = U exactly
        for (const [c, val] of prof) prof.set(c, val * scale);
        inletInfo[p.id] = { lo, hi, prof };
      }
      const eps = this.eps;
      const phiIn = (p, c) => {
        if (p.phase === 'fluid1') return -1;
        if (p.phase === 'fluid2') return 1;
        const { lo, hi } = inletInfo[p.id];
        return Math.tanh((s[c] - lo - p.waterFrac * (hi - lo)) / (Math.SQRT2 * eps));
      };
      const uA = [], vA = [];
      const B = { f: [], axis: [], cell: [], dir: [], kind: [], pid: [], cos: [], vin: [], phin: [] };
      const uK = new Int8Array(N).fill(-1), vK = new Int8Array(N).fill(-1);
      const add = (f, a, b, axis) => {
        const ta = type[a], tb = type[b];
        const K = axis === 0 ? uK : vK;
        if (ta === FLUID && tb === FLUID) { (axis === 0 ? uA : vA).push(f); K[f] = 1; return; }
        if (ta !== FLUID && tb !== FLUID) return;
        const cell = ta === FLUID ? a : b, oc = ta === FLUID ? b : a;
        const dir = ta === FLUID ? 1 : -1;     // outward normal of the fluid cell, along +axis
        let kind = WALL, pid = -1, cs = 0, vin = 0, phin = 0;
        if (type[oc] === PORT) {
          pid = port[oc]; const p = ports[pid];
          if (p.kind === 'inlet') { kind = INLET; vin = -dir * inletInfo[pid].prof.get(oc); phin = phiIn(p, oc); }
          else kind = OUTLET;
        } else cs = Math.cos(thetaAt(oc) * DEG);
        B.f.push(f); B.axis.push(axis); B.cell.push(cell); B.dir.push(dir); B.kind.push(kind);
        B.pid.push(pid); B.cos.push(cs); B.vin.push(vin); B.phin.push(phin);
        K[f] = 2;
      };
      for (let c = 0; c < N; c++) {
        if (c % W < W - 1) add(c, c, c + 1, 0);
        if (c + W < N) add(c, c, c + W, 1);
      }
      this.uA = Int32Array.from(uA); this.vA = Int32Array.from(vA);
      this.B = {
        f: Int32Array.from(B.f), axis: Uint8Array.from(B.axis), cell: Int32Array.from(B.cell),
        dir: Int8Array.from(B.dir), kind: Uint8Array.from(B.kind), pid: Int16Array.from(B.pid),
        cos: Float64Array.from(B.cos), vin: Float64Array.from(B.vin), phin: Float64Array.from(B.phin),
        n: B.f.length,
      };
      // ghost (tangential) velocities: u_ghost = fac * u_interior
      //   wall: Navier slip  fac = (b/h - 1/2)/(b/h + 1/2)   (no-slip -> -1)
      //   inlet: zero tangential velocity (-1); outlet: zero normal gradient (+1)
      const facOf = (a, b) => {
        const pa = type[a] === PORT ? a : type[b] === PORT ? b : -1;
        if (pa >= 0) return ports[port[pa]].kind === 'inlet' ? -1 : 1;
        return this.slipFac;
      };
      const gu = { d: [], s: [], f: [] }, gv = { d: [], s: [], f: [] };
      for (let f = 0; f < N; f++) {
        if (f % W < W - 1 && uK[f] === -1) {
          let src = -1;
          if (f + W < N && uK[f + W] === 1) src = f + W; else if (f - W >= 0 && uK[f - W] === 1) src = f - W;
          if (src >= 0) { gu.d.push(f); gu.s.push(src); gu.f.push(facOf(f, f + 1)); }
        }
        if (f + W < N && vK[f] === -1) {
          let src = -1;
          if (f % W < W - 1 && vK[f + 1] === 1) src = f + 1; else if (f % W > 0 && vK[f - 1] === 1) src = f - 1;
          if (src >= 0) { gv.d.push(f); gv.s.push(src); gv.f.push(facOf(f, f + W)); }
        }
      }
      const pack = (o) => ({ d: Int32Array.from(o.d), s: Int32Array.from(o.s), f: Float64Array.from(o.f) });
      this.gu = pack(gu); this.gv = pack(gv);
    }

    // constant-coefficient Laplacian on the fluid mask, factorised once (banded Cholesky)
    _buildPoisson() {
      const g = this.g, { Nx, Ny, W, N, type, ports } = g;
      const idx = new Int32Array(N).fill(-1);
      const cells = [];
      if (Nx >= Ny) { for (let i = 1; i <= Nx; i++) for (let j = 1; j <= Ny; j++) { const c = i + j * W; if (type[c] === FLUID) { idx[c] = cells.length; cells.push(c); } } }
      else { for (let j = 1; j <= Ny; j++) for (let i = 1; i <= Nx; i++) { const c = i + j * W; if (type[c] === FLUID) { idx[c] = cells.length; cells.push(c); } } }
      const n = cells.length;
      let bw = 1;
      for (const f of this.uA) bw = Math.max(bw, Math.abs(idx[f + 1] - idx[f]));
      for (const f of this.vA) bw = Math.max(bw, Math.abs(idx[f + W] - idx[f]));
      const R = bw + 1, L = new Float64Array(n * R);
      const bOut = new Float64Array(n);
      const link = (c1, c2) => {
        const k1 = idx[c1], k2 = idx[c2], k = Math.max(k1, k2), m = Math.min(k1, k2);
        L[k1 * R + bw] += 1; L[k2 * R + bw] += 1; L[k * R + m - k + bw] = -1;
      };
      for (const f of this.uA) link(f, f + 1);
      for (const f of this.vA) link(f, f + W);
      const B = this.B;
      for (let q = 0; q < B.n; q++) if (B.kind[q] === OUTLET) {
        const k = idx[B.cell[q]]; L[k * R + bw] += 2; bOut[k] += 2 * ports[B.pid[q]].p;
      }
      for (let k = 0; k < n; k++) {
        const rk = k * R, m0 = Math.max(0, k - bw);
        for (let m = m0; m <= k; m++) {
          const rm = m * R;
          let sum = L[rk + m - k + bw];
          for (let q = m0; q < m; q++) sum -= L[rk + q - k + bw] * L[rm + q - m + bw];
          if (m === k) {
            if (!(sum > 1e-12)) throw new Error('Pressure matrix is singular: every connected fluid region needs at least one outlet.');
            L[rk + bw] = Math.sqrt(sum);
          } else L[rk + m - k + bw] = sum / L[rm + bw];
        }
      }
      Object.assign(this, { pIdx: idx, pCells: Int32Array.from(cells), nP: n, bw, Lch: L, bOut, bvec: new Float64Array(n), xvec: new Float64Array(n) });
    }

    _solvePoisson() {
      const n = this.nP, bw = this.bw, R = bw + 1, L = this.Lch, b = this.bvec, x = this.xvec;
      for (let k = 0; k < n; k++) {
        const rk = k * R; let sum = b[k];
        for (let m = Math.max(0, k - bw); m < k; m++) sum -= L[rk + m - k + bw] * x[m];
        x[k] = sum / L[rk + bw];
      }
      for (let k = n - 1; k >= 0; k--) {
        let sum = x[k]; const mEnd = Math.min(n - 1, k + bw);
        for (let m = k + 1; m <= mEnd; m++) sum -= L[m * R + k - m + bw] * x[m];
        x[k] = sum / L[k * R + bw];
      }
    }

    _computeDt(safety) {
      const h = this.g.h, r1 = this.r1, r2 = this.r2;
      const rA = 0.5 * (r1 + r2), rMin = Math.min(r1, r2), muMax = Math.max(this.m1, this.m2);
      let Umax = 0;
      for (const p of this.g.ports) if (p.kind === 'inlet') Umax = Math.max(Umax, 1.5 * Math.abs(p.U));
      const lim = {
        capillary: this.sigma > 0 ? Math.sqrt((rA * h * h * h) / (2 * Math.PI * this.sigma)) : Infinity,
        viscous: (0.25 * h * h * rMin) / muMax,
        advective: Umax > 0 ? h / (2 * Umax) : Infinity,
        cahnHilliard: this.Mc > 0 ? 2 / (this.Mc * ((64 * this.eps * this.eps) / h ** 4 + 16 / (h * h))) : Infinity,
      };
      if (this.depth > 0) lim.heleShaw = (rMin * this.depth * this.depth) / (12 * muMax);
      this.dtLimits = lim;
      this.dt = (safety || 0.5) * Math.min(...Object.values(lim));
      let which = ''; let best = Infinity;
      for (const k in lim) if (lim[k] < best) { best = lim[k]; which = k; }
      this.dtLimitedBy = which;
    }

    initPhi(init) {
      const { W, h } = this.g, e = Math.SQRT2 * this.eps, phi = this.phi;
      const sg = (f) => (f === 'fluid1' ? -1 : 1);
      for (const c of this.fc) {
        const i = c % W, j = (c - i) / W, x = (i - 0.5) * h, y = (j - 0.5) * h;
        let val;
        if (init.fn) val = init.fn(x, y);
        else if (init.type === 'droplet') {
          const r = Math.hypot(x - init.cx, y - init.cy);
          val = sg(init.inside) * Math.tanh((init.R - r) / e);
        } else if (init.type === 'layers') {
          val = Math.tanh((y - init.y0) / e);           // Fluid 1 below y0
        } else {
          val = sg(init.fill) * Math.tanh((x - (init.front || 0)) / e); // other fluid for x < front
        }
        phi[c] = Math.max(-1, Math.min(1, val));
      }
      this._props();
    }

    _props() {
      const { phi, rho, mu, r1, r2, m1, m2 } = this;
      for (const c of this.fc) {
        const vf2 = 0.5 * (1 + Math.max(-1, Math.min(1, phi[c])));
        rho[c] = r1 + (r2 - r1) * vf2; mu[c] = m1 + (m2 - m1) * vf2;
      }
    }

    _fillGhosts() {
      const { u, v, gu, gv } = this;
      for (let n = 0; n < gu.d.length; n++) u[gu.d[n]] = gu.f[n] * u[gu.s[n]];
      for (let n = 0; n < gv.d.length; n++) v[gv.d[n]] = gv.f[n] * v[gv.s[n]];
    }

    step() {
      const g = this.g, W = g.W, N = g.N, h = g.h, type = g.type, ports = g.ports, dt = this.dt;
      const { phi, phiN, psi, G, rho, mu, u, v, us, vs, txy, lap, rate, p, pOld, ph, acc, fc, uA, vA, B } = this;
      const ih = 1 / h, ih2 = ih * ih, eps2 = this.eps * this.eps, Mc = this.Mc, Gfac = this.Gfac;
      const iew = 1 / (Math.SQRT2 * this.eps);
      const ramp = this.ramp > 0 && this.t < this.ramp ? 0.5 * (1 - Math.cos((Math.PI * this.t) / this.ramp)) : 1;

      // 1) mixture properties from phi^n
      this._props();

      // 2) chemical potential psi and surface-tension coefficient G
      for (let n = 0; n < fc.length; n++) lap[fc[n]] = 0;
      for (let n = 0; n < uA.length; n++) { const f = uA[n], d = phi[f + 1] - phi[f]; lap[f] += d; lap[f + 1] -= d; }
      for (let n = 0; n < vA.length; n++) { const f = vA[n], d = phi[f + W] - phi[f]; lap[f] += d; lap[f + W] -= d; }
      for (let q = 0; q < B.n; q++) {
        const c = B.cell[q], k = B.kind[q];
        if (k === WALL) { const pc = phi[c]; lap[c] += h * (-B.cos[q] * Math.max(0, 1 - pc * pc) * iew); } // wetted wall
        else if (k === INLET) lap[c] += 2 * (B.phin[q] - phi[c]);                                       // Dirichlet phi
      }
      for (let n = 0; n < fc.length; n++) {
        const c = fc[n], pc = phi[c];
        psi[c] = -eps2 * lap[c] * ih2 + pc * (pc * pc - 1);
        G[c] = Gfac * psi[c];
        rate[c] = 0;
      }

      // 3) phase transport with u^n (conservative fluxes)
      for (let n = 0; n < uA.length; n++) {
        const a = uA[n], b = a + 1, un = u[a];
        let pf;
        if (un >= 0) { const sl = type[a - 1] === FLUID ? vanLeer(phi[a] - phi[a - 1], phi[b] - phi[a]) : 0; pf = phi[a] + 0.5 * sl; }
        else { const sl = type[b + 1] === FLUID ? vanLeer(phi[b] - phi[b + 1], phi[a] - phi[b]) : 0; pf = phi[b] + 0.5 * sl; }
        const F = un * pf - Mc * (psi[b] - psi[a]) * ih;
        rate[a] -= F; rate[b] += F;
      }
      for (let n = 0; n < vA.length; n++) {
        const a = vA[n], b = a + W, vn = v[a];
        let pf;
        if (vn >= 0) { const sl = type[a - W] === FLUID ? vanLeer(phi[a] - phi[a - W], phi[b] - phi[a]) : 0; pf = phi[a] + 0.5 * sl; }
        else { const sl = type[b + W] === FLUID ? vanLeer(phi[b] - phi[b + W], phi[a] - phi[b]) : 0; pf = phi[b] + 0.5 * sl; }
        const F = vn * pf - Mc * (psi[b] - psi[a]) * ih;
        rate[a] -= F; rate[b] += F;
      }
      const ps = this.portStats;
      for (const st of ps) { st.Q = 0; st.Qw = 0; }
      for (let q = 0; q < B.n; q++) {
        const k = B.kind[q]; if (k === WALL) continue;
        const f = B.f[q], c = B.cell[q], vel = B.axis[q] === 0 ? u[f] : v[f], out = B.dir[q] * vel;
        const pf = k === INLET ? B.phin[q] : phi[c];
        rate[c] -= out * pf;
        const st = ps[B.pid[q]], Q = out * h;
        st.Q += Q; st.Qw += Q * 0.5 * (1 - pf);
      }
      for (const st of ps) { st.V += st.Q * dt; st.Vw += st.Qw * dt; }
      for (let n = 0; n < fc.length; n++) { const c = fc[n]; phiN[c] = phi[c] + dt * rate[c] * ih; }

      // 4) momentum predictor (pressure-free)
      this._fillGhosts();
      for (let c = 0; c < N - W - 1; c++) {
        let ms = 0, cnt = 0;
        if (type[c] === FLUID) { ms += mu[c]; cnt++; }
        if (type[c + 1] === FLUID) { ms += mu[c + 1]; cnt++; }
        if (type[c + W] === FLUID) { ms += mu[c + W]; cnt++; }
        if (type[c + W + 1] === FLUID) { ms += mu[c + W + 1]; cnt++; }
        txy[c] = cnt ? (ms / cnt) * (u[c + W] - u[c] + v[c + 1] - v[c]) * ih : 0;
      }
      const hs = this.depth > 0 ? 12 / (this.depth * this.depth) : 0, gx = this.gx, gy = this.gy;
      for (let n = 0; n < uA.length; n++) {
        const a = uA[n], b = a + 1, uc = u[a];
        const adv = uc * (u[a + 1] - u[a - 1]) * 0.5 * ih + 0.25 * (v[a] + v[b] + v[a - W] + v[b - W]) * (u[a + W] - u[a - W]) * 0.5 * ih;
        const visc = (2 * mu[b] * (u[a + 1] - uc) - 2 * mu[a] * (uc - u[a - 1])) * ih2 + (txy[a] - txy[a - W]) * ih;
        const rf = 0.5 * (rho[a] + rho[b]);
        const fst = 0.5 * (G[a] + G[b]) * (phi[b] - phi[a]) * ih;
        const fr = -hs * 0.5 * (mu[a] + mu[b]) * uc;
        us[a] = uc + dt * (-adv + (visc + fst + fr) / rf + gx);
      }
      for (let n = 0; n < vA.length; n++) {
        const a = vA[n], b = a + W, vc = v[a];
        const adv = vc * (v[a + W] - v[a - W]) * 0.5 * ih + 0.25 * (u[a] + u[a - 1] + u[b] + u[b - 1]) * (v[a + 1] - v[a - 1]) * 0.5 * ih;
        const visc = (2 * mu[b] * (v[a + W] - vc) - 2 * mu[a] * (vc - v[a - W])) * ih2 + (txy[a] - txy[a - 1]) * ih;
        const rf = 0.5 * (rho[a] + rho[b]);
        const fst = 0.5 * (G[a] + G[b]) * (phi[b] - phi[a]) * ih;
        const fr = -hs * 0.5 * (mu[a] + mu[b]) * vc;
        vs[a] = vc + dt * (-adv + (visc + fst + fr) / rf + gy);
      }
      for (let q = 0; q < B.n; q++) {
        const k = B.kind[q]; if (k === OUTLET) continue;
        (B.axis[q] === 0 ? us : vs)[B.f[q]] = k === INLET ? ramp * B.vin[q] : 0;
      }
      for (let q = 0; q < B.n; q++) if (B.kind[q] === OUTLET) {           // zero normal gradient
        const f = B.f[q], ax = B.axis[q], arr = ax === 0 ? us : vs;
        arr[f] = arr[f - B.dir[q] * (ax === 0 ? 1 : W)];
      }

      // 5) pressure Poisson (Dodd–Ferrante constant-coefficient splitting)
      const r0 = this.rho0, ir0 = 1 / r0, cdt = r0 / dt;
      for (let n = 0; n < fc.length; n++) { const c = fc[n]; ph[c] = 2 * p[c] - pOld[c]; acc[c] = 0; }
      for (let n = 0; n < uA.length; n++) {
        const a = uA[n], b = a + 1, rf = 0.5 * (rho[a] + rho[b]);
        const val = cdt * us[a] + (1 - r0 / rf) * (ph[b] - ph[a]) * ih;
        acc[a] += val; acc[b] -= val;
      }
      for (let n = 0; n < vA.length; n++) {
        const a = vA[n], b = a + W, rf = 0.5 * (rho[a] + rho[b]);
        const val = cdt * vs[a] + (1 - r0 / rf) * (ph[b] - ph[a]) * ih;
        acc[a] += val; acc[b] -= val;
      }
      for (let q = 0; q < B.n; q++) {
        const f = B.f[q], c = B.cell[q], out = B.dir[q] * (B.axis[q] === 0 ? us[f] : vs[f]);
        acc[c] += cdt * out;
        if (B.kind[q] === OUTLET) acc[c] += (1 - r0 / rho[c]) * (ports[B.pid[q]].p - ph[c]) * 2 * ih;
      }
      const pc = this.pCells, bv = this.bvec, bo = this.bOut;
      for (let k = 0; k < this.nP; k++) bv[k] = bo[k] - h * acc[pc[k]];
      this._solvePoisson();
      const x = this.xvec;
      for (let k = 0; k < this.nP; k++) { const c = pc[k]; pOld[c] = p[c]; p[c] = x[k]; }

      // 6) projection / velocity correction
      for (let n = 0; n < uA.length; n++) {
        const a = uA[n], b = a + 1, rf = 0.5 * (rho[a] + rho[b]);
        u[a] = us[a] - dt * ((p[b] - p[a]) * ih * ir0 + (1 / rf - ir0) * (ph[b] - ph[a]) * ih);
      }
      for (let n = 0; n < vA.length; n++) {
        const a = vA[n], b = a + W, rf = 0.5 * (rho[a] + rho[b]);
        v[a] = vs[a] - dt * ((p[b] - p[a]) * ih * ir0 + (1 / rf - ir0) * (ph[b] - ph[a]) * ih);
      }
      for (let q = 0; q < B.n; q++) {
        const f = B.f[q], ax = B.axis[q], src = ax === 0 ? us : vs, dst = ax === 0 ? u : v;
        if (B.kind[q] === OUTLET) {
          const c = B.cell[q], po = ports[B.pid[q]].p;
          const gn = (po - p[c]) * 2 * ih, gh = (po - ph[c]) * 2 * ih;
          dst[f] = src[f] - dt * B.dir[q] * (gn * ir0 + (1 / rho[c] - ir0) * gh);
        } else dst[f] = src[f];
      }

      // 7) advance phase field
      for (let n = 0; n < fc.length; n++) { const c = fc[n]; phi[c] = phiN[c]; }
      this.t += dt; this.nstep++;
      if (!Number.isFinite(phi[fc[0]]) || !Number.isFinite(u[uA[uA.length >> 1]])) throw new Error('Simulation diverged at t = ' + this.t.toExponential(3) + ' s');
    }

    run(n) { for (let k = 0; k < n; k++) this.step(); }

    // ---------------------------------------------------------- diagnostics
    maxDiv() {
      const { W, h } = this.g, { u, v } = this; let m = 0;
      for (const c of this.fc) m = Math.max(m, Math.abs(u[c] - u[c - 1] + v[c] - v[c - W]) / h);
      return m;
    }
    fluid1Area() { let A = 0; const h2 = this.g.h ** 2; for (const c of this.fc) A += 0.5 * (1 - this.phi[c]) * h2; return A; }

    // cell-centred fields on the Nx x Ny grid (row j=0 is the bottom); NaN outside the fluid
    fields() {
      const { Nx, Ny, W, type } = this.g, { u, v, p, phi } = this;
      const n = Nx * Ny, out = { u: new Float32Array(n), v: new Float32Array(n), speed: new Float32Array(n), p: new Float32Array(n), phi: new Float32Array(n) };
      for (let j = 1; j <= Ny; j++) for (let i = 1; i <= Nx; i++) {
        const c = i + j * W, k = i - 1 + (j - 1) * Nx;
        if (type[c] !== FLUID) { out.u[k] = out.v[k] = out.speed[k] = out.p[k] = out.phi[k] = NaN; continue; }
        const uc = 0.5 * (u[c - 1] + u[c]), vc = 0.5 * (v[c - W] + v[c]);
        out.u[k] = uc; out.v[k] = vc; out.speed[k] = Math.hypot(uc, vc); out.p[k] = p[c]; out.phi[k] = phi[c];
      }
      return out;
    }
  }

  return { Simulation, buildGeometry, defaults, SIDES, FLUID, SOLID, PORT };
});
