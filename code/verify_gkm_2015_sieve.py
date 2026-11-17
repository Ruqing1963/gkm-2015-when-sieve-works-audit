"""Verification script for the audit note on
A. Granville, D. Koukoulopoulos and K. Matomaki, "When the sieve works",
Duke Math. J. 164 (2015), no. 10, 1935-1969 (arXiv:1205.0413v2; page numbers of the arXiv version).

Run:  python verify_gkm_2015_sieve.py
Writes ../data/verify_gkm_2015_sieve.txt and ../figures/fig_gkm_2015_audit.{pdf,png}
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT_TXT = os.path.join(ROOT, "data", "verify_gkm_2015_sieve.txt")
OUT_PDF = os.path.join(ROOT, "figures", "fig_gkm_2015_audit.pdf")
OUT_PNG = os.path.join(ROOT, "figures", "fig_gkm_2015_audit.png")
_lines = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    _lines.append(s)


def check(name, ok):
    say(f"[{'PASS' if ok else 'FAIL'}] {name}")
    if not ok:
        check.failures += 1


check.failures = 0

# ===========================================================================
say("=" * 76)
say("1. Proposition 2.6 (pp. 10-11): exponents in (2.7) and (2.8)")
say("=" * 76)


def expo(L, c, scale=1.0, half=True):
    """-log((v log v)^{-v/2}) / log x  (or with full exponent v if half=False),
    v = c * L' / log L' where L' = scale * L = log of the point where Prop. 2.5 is applied;
    the result is normalised by L = log x."""
    Lp = scale * L
    v = c * Lp / math.log(Lp)
    p = v / 2 if half else v
    return p * math.log(v * math.log(v)) / L


Ls = [1e2, 1e3, 1e4, 1e6, 1e10, 1e20, 1e40, 1e80]
say(f"   {'log x':>8} {'2nd claim v^-v':>15} {'(2.7) v^-v/2':>13} {'(2.8)':>8} {'(2.7) opt':>10} {'(2.8) opt':>10}")
for L in Ls:
    say(f"   {L:>8.0e} {expo(L, 0.4, half=False):>15.4f} {expo(L, 0.4):>13.4f} {expo(L, 0.4, 0.5):>8.4f}"
        f" {expo(L, 0.499):>10.4f} {expo(L, 0.499, 0.5):>10.4f}")
say("   claimed: 1/3 = 0.3333 (2nd claim and (2.7)), 1/6 = 0.1667 ((2.8)); limits: 2/5, 1/5, 1/10, 1/4, 1/8")
check("second claim: exponent -> 2/5 > 1/3 (x^-1/3 is valid there for large x)",
      expo(1e40, 0.4, half=False) > 1 / 3 and abs(expo(1e80, 0.4, half=False) - 0.4) < 0.01)
check("(2.7): exponent < 1/5 < 1/3 for all tested log x, and -> 1/5",
      all(expo(L, 0.4) < 0.2 for L in Ls) and abs(expo(1e80, 0.4) - 0.2) < 0.005)
check("(2.8): exponent < 1/10 < 1/6 for all tested log x, and -> 1/10",
      all(expo(L, 0.4, 0.5) < 0.1 for L in Ls) and abs(expo(1e80, 0.4, 0.5) - 0.1) < 0.003)
check("with the largest admissible v = (1/2-eps)log x/loglog x the limits are 1/4 and 1/8 (< 1/3, < 1/6)",
      abs(expo(1e80, 0.499) - 0.2495) < 0.005 and abs(expo(1e80, 0.499, 0.5) - 0.12475) < 0.003)
# exact asymptotics: (v/2) log(v log v) = (log x)/5 - (log x)(log 2.5)/(5 loglog x) + O(log x logloglog x/(loglog x)^2)
L = 1e40
v = 0.4 * L / math.log(L)
lhs = (v / 2) * math.log(v * math.log(v))
approx = L / 5 * (1 - math.log(2.5) / math.log(L))
check("(v/2)log(v log v) = (log x/5)(1 - log 2.5/loglog x + o(1/loglog x)) at log x = 1e40",
      abs(lhs / approx - 1) < 1e-2)
say("   hence e^{O(v)}(v log v)^{-v/2} = x^{-1/5 + O(1/loglog x)}: not << x^{-1/3},"
    " and not << (2u)^{-u/2} when u > v (that term is <= (2v)^{-v/2}).")
uu = 2 * v
check("for u = 2v: (2u)^{-u/2} is smaller than (v log v)^{-v/2} by many orders (log scale)",
      (uu / 2) * math.log(2 * uu) > (v / 2) * math.log(v * math.log(v)) + 10)

# ===========================================================================
say("=" * 76)
say("2. Lemma 8.5 (p. 25): counterexamples and the corrected statement")
say("=" * 76)


def ell_sum_interval(T, ell):
    return (ell * T[0], ell * T[1])


def window_hits(T, win, ellmax):
    """ell for which ell*T (open) meets the window (a, b]"""
    out = []
    for ell in range(1, ellmax + 1):
        lo, hi = ell_sum_interval(T, ell)
        if hi > win[0] and lo < win[1]:
            out.append(ell)
    return out


u, v_, e = 1.0, 4.0, math.e
for name, T, w in (("(a)", (0.55, 0.9), 0.5), ("(b)", (0.6, 0.7), 1.0)):
    lo_ok = T[0] > 1 / (e * v_) and T[1] <= 1 / u
    ellmax = int(e * v_ * w)
    hits = window_hits(T, (w - 1 / v_, w), ellmax)
    say(f"   {name}: u=1, v=4, w={w}, T={T} in (1/(4e),1] = ({1/(e*v_):.4f},1]: {lo_ok};"
        f" ell <= evw = {e*v_*w:.2f}; ell with ell*T meeting (w-1/v,w]: {hits}")
    check(f"{name}: left-hand side of Lemma 8.5 vanishes for every ell <= evw (and every ell >= 1)",
          lo_ok and hits == [] and window_hits(T, (w - 1 / v_, w), 50) == [])
check("(b): with the window (w-1/u, w] = (0,1], ell = 1 meets it",
      window_hits((0.6, 0.7), (0.0, 1.0), 10) == [1])


def lemma_check(T_ints, lo_t, w, width, grid=20000):
    """Numerically compute beta_l = mu^l(s_l in (w-width, w]) / tau^l for the measure dt/t on T,
    l = 1..K = floor(w/lo_t), via FFT convolution on a grid; returns (max beta_l, 1/K, 1/(w/lo_t))."""
    h = w / grid
    xs = (np.arange(grid + 1) + 0.5) * h
    dens = np.zeros(grid + 1)
    for a, b in T_ints:
        m = (xs > a) & (xs < b)
        dens[m] = 1 / xs[m] * h
    tau = dens.sum()
    K = int(math.floor(w / lo_t + 1e-12))
    p = dens / tau
    cur = p.copy()
    betas = []
    n = len(p)
    for ell in range(1, K + 1):
        pos = (np.arange(len(cur)) + 0.5 * ell) * h
        betas.append(float(cur[(pos > w - width) & (pos <= w)].sum()))
        nfft = 1 << int(np.ceil(np.log2(len(cur) + n)))
        cur = np.fft.irfft(np.fft.rfft(cur, nfft) * np.fft.rfft(p, nfft), nfft)[: len(cur) + n - 1]
        cur = np.clip(cur, 0, None)
        cur = cur[: int((w + 1) / h)]
    return max(betas), betas


rng = np.random.default_rng(11)
ok_corr = True
worst = 1e9
for trial in range(60):
    u_ = float(rng.uniform(1, 3))
    v2 = float(rng.uniform(u_, 6))
    lo_t, hi_t = 1 / (math.e * v2), 1 / u_
    cuts = np.sort(rng.uniform(lo_t, hi_t, size=2 * int(rng.integers(1, 4))))
    T_ints = [(cuts[2 * i], cuts[2 * i + 1]) for i in range(len(cuts) // 2)]
    w = float(rng.uniform(1 / u_, 3.0))
    mb, _ = lemma_check(T_ints, lo_t, w, 1 / u_, grid=6000)
    target = 1 / (math.e * v2 * w)
    worst = min(worst, mb / target)
    ok_corr &= mb >= target * 0.999
check("corrected Lemma 8.5' (window (w-1/u, w], w >= 1/u): max_l beta_l >= 1/(evw) in 60 random cases", ok_corr)
say(f"   smallest ratio max_l beta_l / (1/(evw)) over the trials: {worst:.3f}")

# ===========================================================================
say("=" * 76)
say("3. Proposition 8.1 (p. 24): w log w versus sigma^{-1} log(1/sigma)")
say("=" * 76)
lam1 = 1.0
say(f"   lambda_1 = {lam1}: w = (1+lambda_1)/sigma")
say(f"   {'sigma':>8} {'w log w':>10} {'s^-1 log(1/s)':>14} {'s^-1(1+log 1/s)':>16}")
rows3 = []
for s in (0.01, 0.1, 0.3, 0.6, 0.9, 0.99, 0.999):
    wv = (1 + lam1) / s
    a, b, c = wv * math.log(wv), math.log(1 / s) / s, (1 + math.log(1 / s)) / s
    rows3.append((s, a, b, c))
    say(f"   {s:>8} {a:>10.4f} {b:>14.4f} {c:>16.4f}")
check("w log w / (sigma^-1 log(1/sigma)) is unbounded as sigma -> 1",
      rows3[-1][1] / rows3[-1][2] > 1000)
check("w log w <= C sigma^-1 (1 + log(1/sigma)) with C = 2(1+lambda_1) on (0, 1]",
      all(r[1] <= 2 * (1 + lam1) * r[3] for r in rows3))
say("   For sigma >= 1/2, w <= 2(1+lambda_1) and w^{O(w)} = O(1); for sigma < 1/2, 1 <= log(1/sigma)/log 2.")
say("   So w^{O(w)} << sigma^{-O(1/sigma)} up to a constant, and Proposition 8.1 holds as stated.")

# ===========================================================================
say("=" * 76)
say("4. Hypothesis P => Hypothesis A (p. 15): endpoints of P")
say("=" * 76)
ok4 = True
for N, u4, v4 in ((1000, 3.0, 4.0), (5000, 1.7, 3.0), (20000, 2.5, 5.0)):
    ev = math.e * v4
    A = [a for a in range(int(N / ev) + 1, int(N / u4) + 1)]
    x_log = math.floor(N) + 1
    up_bad = [a for a in A if a + 1 > x_log / u4]
    lo_bad = [a for a in A if a < x_log / ev]
    loss = sum(1 / a for a in up_bad + lo_bad)
    say(f"   N={N}, u={u4}, v={v4}: a with e^(a+1) > x^(1/u): {up_bad};  a with e^a < x^(1/ev): {lo_bad};"
        f"  sum 1/a removed = {loss:.2e}  (<= 2ev/N = {2*ev/N:.2e})")
    ok4 &= len(up_bad) <= 1 and len(lo_bad) <= 1 and loss <= 2 * ev / N
check("at most one a at each end violates P in (x^(1/ev), x^(1/u)]; removing them costs <= 2ev/N", ok4)
check("at u = 1 the upper endpoint is never violated", all(
    not [a for a in range(int(N / (math.e * 3)) + 1, N + 1) if a + 1 > N + 1] for N in (100, 1000)))

say("=" * 76)
say(f"TOTAL FAILURES: {check.failures}")
with open(OUT_TXT, "w", encoding="utf-8") as fh:
    fh.write("\n".join(_lines) + "\n")

# ===========================================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 7.5,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42, "axes.linewidth": 0.7})
C1, C2, C3, C4 = "#1f4e79", "#c0392b", "#2e8b57", "#7f7f7f"
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.4))
fig.subplots_adjust(left=0.08, right=0.985, bottom=0.14, top=0.74, wspace=0.28)

Lg = np.logspace(2, 60, 300)
h1, = ax1.plot(Lg, [expo(L, 0.4) for L in Lg], color=C1, lw=1.6, label=r"(2.7) as proved: $\to 1/5$")
h2, = ax1.plot(Lg, [expo(L, 0.4, 0.5) for L in Lg], color=C2, lw=1.6, label=r"(2.8) as proved: $\to 1/10$")
h3, = ax1.plot(Lg, [expo(L, 0.499) for L in Lg], color=C1, lw=1.0, ls="--",
               label=r"optimal $v$: $\to 1/4$, $1/8$")
ax1.plot(Lg, [expo(L, 0.499, 0.5) for L in Lg], color=C2, lw=1.0, ls="--")
h4 = ax1.axhline(1 / 3, color=C1, lw=1.0, ls=":", label=r"claimed $1/3$ (2.7) and $1/6$ (2.8)")
ax1.axhline(1 / 6, color=C2, lw=1.0, ls=":")
ax1.set_xscale("log")
ax1.set_xlim(1e2, 1e60)
ax1.set_ylim(0, 0.36)
ax1.set_xlabel(r"$\log x$")
ax1.set_ylabel(r"exponent $\theta$ in $x^{-\theta}$")
ax1.text(0.97, 0.04, "(a)", transform=ax1.transAxes, ha="right", va="bottom", fontsize=10, fontweight="bold")

T = (0.6, 0.7)
for ell in range(1, 4):
    lo, hi = ell * T[0], ell * T[1]
    ax2.plot([lo, hi], [ell, ell], color=C1, lw=4, solid_capstyle="butt")
h5, = ax2.plot([], [], color=C1, lw=4, label=r"$\ell T=(0.6\ell,0.7\ell)$")
h6 = ax2.axvspan(0.75, 1.0, color=C2, alpha=0.25, lw=0, label=r"window $(w-1/v,w]=(0.75,1]$")
h7 = ax2.axvspan(0.0, 1.0, ymin=0.0, ymax=0.08, color=C3, alpha=0.5, lw=0,
                 label=r"corrected window $(w-1/u,w]=(0,1]$")
ax2.set_xlim(0, 2.3)
ax2.set_ylim(0.5, 3.6)
ax2.set_yticks([1, 2, 3])
ax2.set_xlabel(r"$t_1+\dots+t_\ell$")
ax2.set_ylabel(r"$\ell$")
ax2.text(0.97, 0.95, "(b)", transform=ax2.transAxes, ha="right", va="top", fontsize=10, fontweight="bold")
ax2.text(2.25, 0.62, r"$u=1,\ v=4,\ w=1$", ha="right", va="bottom", fontsize=7.5)

fig.legend(handles=[h1, h2, h3, h4, h5, h6, h7], loc="upper center", ncol=2, frameon=False,
           bbox_to_anchor=(0.5, 1.0), columnspacing=1.2, handlelength=2.2, borderaxespad=0.2)
fig.savefig(OUT_PDF)
fig.savefig(OUT_PNG, dpi=200)
sys.exit(1 if check.failures else 0)
