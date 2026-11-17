# A note on the smooth-number exponents and the additive windows in "When the sieve works"

Ruqing Chen (2026). doi:[10.5281/zenodo.24297466](https://doi.org/10.5281/zenodo.24297466)

A. Granville, D. Koukoulopoulos and K. Matomäki, *When the sieve works*, Duke Math. J. 164 (2015), no. 10,
1935–1969, doi:[10.1215/00127094-3120891](https://doi.org/10.1215/00127094-3120891); arXiv:1205.0413.
Page numbers refer to arXiv:1205.0413v2 (26 pages).

* **Proposition 2.6 (pp. 10–11).** With `v = log x/(2.5 loglog x)`, the factor `e^{O(v)}(v log v)^{-v/2}` in
  the proof is `x^{-1/5+o(1)}`, not `<< x^{-1/3}` as in (2.7); the rescaling in (2.8) gives `x^{-1/10+o(1)}`,
  not `x^{-1/6}`. With the largest admissible `v` the method gives `x^{-1/4}` and `x^{-1/8}`; the corrected
  Proposition 2.6 holds with `x^{-1/8+eps}`. Only a power saving is used later.
* **Lemma 8.5 (p. 25, stated without proof)** is false: for `u = 1, v = 4, w = 1, T = (0.6, 0.7)` the sums
  `l*T` skip the window `(0.75, 1]` for every `l`. A corrected lemma (window width = maximal step) is proved;
  its version for `T in (1/ev, 1/v]` is what Proposition 8.4 uses, so Proposition 8.4 is unaffected.
* **Hypothesis P => Hypothesis A (p. 15).** The set of primes can leave `(x^{1/ev}, x^{1/u}]` at both ends
  (e.g. `N = 1000, u = 3, v = 4`: `a = 333` and `a = 92`); removing at most two elements of `A` costs
  `<= 2ev/N`, absorbed by the existing `O(v^2/N)`.
* **Proposition 8.1 (p. 24).** `w log w << sigma^{-1} log(1/sigma)` fails as `sigma -> 1`, but there
  `w^{O(w)} = O(1)`; the proposition holds as stated.

No main result of the paper is affected.

| path | content |
|---|---|
| `gkm_when_sieve_works_audit.tex`, `.pdf` | the note (4 pages) |
| `code/verify_gkm_2015_sieve.py` | exponent asymptotics, the Lemma 8.5 counterexamples, 60 numerical tests of the corrected lemma, Prop. 8.1 and P => A checks, and the figure |
| `data/verify_gkm_2015_sieve.txt` | output of the script |
| `figures/fig_gkm_2015_audit.pdf`, `.png` | Figure 1 |

Requirements: Python 3 with numpy, matplotlib. Run `python code/verify_gkm_2015_sieve.py`.
