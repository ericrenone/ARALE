#!/usr/bin/env python3
"""ARALE numerical model.

Fixed-point Monte Carlo of the complete ARALE loop:
  * exact E8 closest-point decoder (Conway-Sloane two-coset method), verified against brute force
  * zero-forcing equalization with a tracked channel estimate H_hat
  * decision-directed NLMS rank-one tracking with a configurable decision delay
  * Sherman-Morrison (explicit inverse) and Givens-QR representations in fixed point
  * CORDIC rotation error at a given iteration count and word length
  * operation, latency, energy and area counts for the n = 8 array

Usage:
    python3 arale_model.py all            # run every experiment (several minutes)
    python3 arale_model.py decoder        # brute-force check of the E8 decoder
    python3 arale_model.py ser            # SER vs SNR, with and without mismatch
    python3 arale_model.py track          # tracker floor, drift recovery, capture range
    python3 arale_model.py sm             # Sherman-Morrison vs QR error growth in 16/20-bit arithmetic
    python3 arale_model.py cordic         # CORDIC rotation error
    python3 arale_model.py budget         # operation counts, latency, energy, area

Requires numpy only.
"""
import sys
import numpy as np

rng = np.random.default_rng(1)

# --------------------------------------------------------------------------- E8

def dec_D8(z):
    """Nearest point of D8 (integer vectors with even coordinate sum)."""
    r = np.rint(z)
    if int(r.sum()) % 2 != 0:
        err = z - r
        k = np.argmax(np.abs(err))
        r[k] += 1.0 if err[k] > 0 else -1.0
    return r


def dec_E8(z):
    """Nearest point of E8 = D8 U (D8 + 1/2): two D8 decodes and one compare."""
    c1 = dec_D8(z)
    c2 = dec_D8(z - 0.5) + 0.5
    return c1 if np.sum((z - c1) ** 2) <= np.sum((z - c2) ** 2) else c2


def min_vectors():
    """The 240 minimal vectors of E8 (squared norm 2)."""
    V = []
    for i in range(8):
        for j in range(i + 1, 8):
            for si in (-1, 1):
                for sj in (-1, 1):
                    v = np.zeros(8); v[i] = si; v[j] = sj; V.append(v)
    for m in range(256):
        s = np.array([1 if (m >> b) & 1 else -1 for b in range(8)]) * 0.5
        if (s < 0).sum() % 2 == 0:
            V.append(s)
    V = np.array(V)
    assert len(V) == 240
    return V


MV = min_vectors()

E8_GEN = np.array([
    [2, 0, 0, 0, 0, 0, 0, 0],
    [-1, 1, 0, 0, 0, 0, 0, 0],
    [0, -1, 1, 0, 0, 0, 0, 0],
    [0, 0, -1, 1, 0, 0, 0, 0],
    [0, 0, 0, -1, 1, 0, 0, 0],
    [0, 0, 0, 0, -1, 1, 0, 0],
    [0, 0, 0, 0, 0, -1, 1, 0],
    [.5, .5, .5, .5, .5, .5, .5, .5],
], dtype=float)


def check_decoder(N=2000):
    ok = 0
    cands = np.vstack([np.zeros(8), MV])
    for _ in range(N):
        z = rng.normal(0, 0.25, 8)
        bf = cands[np.argmin(((cands - z) ** 2).sum(1))]
        ok += np.allclose(dec_E8(z), bf)
    print(f"[decoder] agrees with brute force on {ok}/{N} points near the origin")
    print(f"[decoder] generator determinant = {np.linalg.det(E8_GEN):.6f} (unimodular expected: 1)")
    G = E8_GEN @ E8_GEN.T
    print(f"[decoder] Gram diagonal = {np.diag(G)}  (all even: {np.all(np.diag(G) % 2 == 0)})")


# --------------------------------------------------------------------------- SER

def sigma2(snr_db, Es=2.0):
    """Per-dimension noise variance for SNR per dimension with Es = 2 (minimal vectors)."""
    return Es / (8 * 10 ** (snr_db / 10))


def ser(snr_db, delta_norm=0.0, trials=20000):
    s2 = sigma2(snr_db); errs = 0
    for _ in range(trials):
        x = MV[rng.integers(240)]
        if delta_norm > 0:
            D = rng.normal(0, 1, (8, 8)); D *= delta_norm / np.linalg.norm(D, 'fro')
            z = x + D @ x + rng.normal(0, np.sqrt(s2), 8)
        else:
            z = x + rng.normal(0, np.sqrt(s2), 8)
        errs += not np.allclose(dec_E8(z), x)
    return errs / trials


def run_ser():
    print("SNR/dim  exact    |D|=0.05  |D|=0.10  |D|=0.20")
    for snr in (6, 7, 8, 9, 10):
        print(f"{snr:5d}   {ser(snr):.4f}   {ser(snr, .05):.4f}    {ser(snr, .10):.4f}    {ser(snr, .20):.4f}")


# --------------------------------------------------------------------------- tracking

def make_channel(r, kappa):
    U, _ = np.linalg.qr(r.normal(0, 1, (8, 8))); V, _ = np.linalg.qr(r.normal(0, 1, (8, 8)))
    return U @ np.diag(np.linspace(1, kappa, 8)) @ V.T


def track(mu, snr_db=12, D=65, T=3000, kappa=2.0, seed=7, init=None, step=None):
    """Decision-directed NLMS with decision delay D. `init` perturbs H_hat at t=0; `step` applies a
    10% single-line gain step at t=1000. Returns mismatch trace and decision-error trace."""
    r = np.random.default_rng(seed); s2 = sigma2(snr_db)
    H = make_channel(r, kappa); Hh = H.copy()
    if init is not None:
        Hh = init(H, r)
    ys, xs, mis, errs = [], [], [], []
    for t in range(T):
        if step is not None and t == 1000:
            H[3, :] *= 1 + step
        x = MV[r.integers(240)]; y = H @ x + r.normal(0, np.sqrt(s2), 8)
        xh = dec_E8(np.linalg.solve(Hh, y)); errs.append(not np.allclose(xh, x))
        ys.append(y); xs.append(xh)
        if t >= D:
            e = ys[t - D] - Hh @ xs[t - D]
            Hh = Hh + mu * np.outer(e, xs[t - D]) / (xs[t - D] @ xs[t - D])
        mis.append(np.linalg.norm(np.linalg.solve(Hh, H) - np.eye(8), 'fro'))
    return np.array(mis), np.array(errs)


def run_track():
    print("-- tracker floor and 10% single-line step recovery (SNR 12 dB, kappa 2, D = 65)")
    for mu in (1/8, 1/16, 1/32, 1/64):
        mis, errs = track(mu, step=0.10)
        base = np.median(mis[500:1000]); peak = mis[1000:].max()
        tau = np.argmax(mis[1000:] < base + (peak - base) / np.e)
        print(f"mu=1/{int(1/mu):2d}: floor={base:.3f} peak={peak:.3f} 1/e recovery={tau} symbols errors={errs[1000:].sum()}")
    print("-- capture range from an initial error in H_hat (mu = 1/16)")
    def rot(deg):
        def f(H, r):
            th = np.deg2rad(deg); G = np.eye(8); G[0, 0] = G[1, 1] = np.cos(th); G[0, 1] = np.sin(th); G[1, 0] = -np.sin(th)
            return H @ G
        return f
    def full(size):
        def f(H, r):
            P = r.normal(0, 1, (8, 8)); P *= size / np.linalg.norm(P, 'fro'); return H @ (np.eye(8) + P)
        return f
    def line(err):
        def f(H, r):
            Hh = H.copy(); Hh[2, :] *= 1 + err; return Hh
        return f
    for name, init in (("line 30%", line(.3)), ("line 120%", line(1.2)), ("rot 20deg", rot(20)), ("rot 45deg", rot(45)),
                       ("full 0.5", full(.5)), ("full 0.8", full(.8))):
        mis, errs = track(1/16, T=4000, init=init, seed=3)
        t10 = np.argmax(mis < 0.10) if (mis < 0.10).any() else -1
        print(f"{name:10s}: initial={mis[0]:.3f} symbols to <0.10: {t10:4d} SER first 200={errs[:200].mean():.3f}")


# --------------------------------------------------------------------------- SM vs QR

def run_sm(updates=300, frac_bits=16):
    q = 2.0 ** -frac_bits; rnd = lambda A: np.round(A / q) * q
    print(f"-- Sherman-Morrison explicit inverse, {frac_bits} fraction bits, {updates} rank-one updates")
    for kappa in (2, 10, 30):
        r = np.random.default_rng(3); H = make_channel(r, kappa); Hi = rnd(np.linalg.inv(H)); e = []
        for k in range(updates):
            u = r.normal(0, 0.02, 8); v = MV[r.integers(240)]
            H = H + np.outer(u, v)
            Hiu = rnd(Hi @ u); vHi = rnd(v @ Hi); den = 1 + v @ Hiu
            Hi = rnd(Hi - rnd(np.outer(Hiu, vHi)) / den)
            e.append(np.linalg.norm(Hi @ H - np.eye(8), 'fro'))
        print(f"kappa={kappa:2d}: after 30={e[29]:.4f} after {updates}={e[-1]:.4f} growth/update={(e[-1]-e[29])/(updates-30):.1e}")
    print("-- Givens QR re-factorization, 20 fraction bits, 10^4 updates, kappa 10")
    q = 2.0 ** -20; rnd = lambda A: np.round(A / q) * q
    def qr_fp(A):
        A = A.copy(); Q = np.eye(8)
        for j in range(8):
            for i in range(7, j, -1):
                a, b = A[i-1, j], A[i, j]; rr = np.hypot(a, b)
                if rr == 0: continue
                c, s = a / rr, b / rr; G = np.eye(8); G[i-1, i-1] = c; G[i, i] = c; G[i-1, i] = s; G[i, i-1] = -s
                A = rnd(G @ A); Q = rnd(G @ Q)
        return Q, A
    r = np.random.default_rng(11); H = make_channel(r, 10); Hh = rnd(H.copy()); worst = 0
    for k in range(10000):
        u = r.normal(0, 0.02, 8); v = MV[r.integers(240)]
        H = H + np.outer(u, v); Hh = rnd(Hh + np.outer(rnd(u), v))
        if k % 500 == 0 or k == 9999:
            Q, R = qr_fp(Hh); worst = max(worst, np.linalg.norm(np.linalg.solve(R, Q @ H) - np.eye(8), 'fro'))
    print(f"worst |R^-1 Q^T H - I|_F = {worst:.2e}")


# --------------------------------------------------------------------------- CORDIC

def cordic_rotate(x, y, theta, n_iter=16, frac_bits=20):
    q = 2.0 ** -frac_bits
    K = np.prod([1 / np.sqrt(1 + 2.0 ** (-2 * k)) for k in range(n_iter)]); z = theta
    for k in range(n_iter):
        d = 1 if z >= 0 else -1
        x, y = np.round((x - d * y * 2.0 ** -k) / q) * q, np.round((y + d * x * 2.0 ** -k) / q) * q
        z -= d * np.arctan(2.0 ** -k)
    return x * K, y * K


def run_cordic(n_iter=16, frac_bits=20, N=5000):
    errs = []
    for _ in range(N):
        th = rng.uniform(-np.pi / 2, np.pi / 2); x, y = rng.uniform(-1, 1, 2)
        xr, yr = cordic_rotate(x, y, th, n_iter, frac_bits)
        errs.append(np.hypot(xr - (x * np.cos(th) - y * np.sin(th)), yr - (x * np.sin(th) + y * np.cos(th))))
    K = 1 / np.prod([1 / np.sqrt(1 + 2.0 ** (-2 * k)) for k in range(n_iter)])
    print(f"[cordic] {n_iter} iterations, {frac_bits}-bit: max err={max(errs):.2e} rms={np.sqrt(np.mean(np.square(errs))):.2e} K={K:.6f}")


# --------------------------------------------------------------------------- budget

def run_budget():
    adds_per_cordic = 16 * 3
    cordic_ops = 28 + 140 + 64 + 36           # vectoring + rotation + y-column + solve row
    other_ops = 120 + 136                      # snap + tracker shift-adds
    E20 = {"45 nm": 0.1 * 20 / 32, "7 nm": 0.03 * 20 / 32}
    print(f"[budget] CORDIC ops/vector = {cordic_ops}; other ops/vector = {other_ops}")
    for node, e in E20.items():
        E = cordic_ops * adds_per_cordic * e + other_ops * e
        print(f"[budget] {node}: datapath {E:.0f} pJ/vector; x2.5 overhead {2.5*E/1000:.2f} nJ; {2.5*E/1000:.2f} mW at 1 Gvector/s")
    print(f"[budget] 22 nm interpolated: ~450 pJ datapath, ~1.1 nJ with overhead")
    print(f"[budget] RSFQ switching only: {cordic_ops*adds_per_cordic*20*30*2.5e-18*1e12:.1f} pJ/vector")
    print(f"[budget] latency: 4 + 9 + 16 + 8 + 18 + 4 = {4+9+16+8+18+4} clocks")
    for name, cells in (("full-rate", 44 + 72), ("reduced", 44 + 8)):
        print(f"[budget] {name}: {cells} CORDIC cells x 15 kGE = {cells*15/1000:.2f} MGE")


# --------------------------------------------------------------------------- main

if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    steps = {"decoder": check_decoder, "ser": run_ser, "track": run_track, "sm": run_sm, "cordic": run_cordic, "budget": run_budget}
    for k, f in steps.items():
        if what in ("all", k):
            print(f"\n=== {k} ===")
            f()
