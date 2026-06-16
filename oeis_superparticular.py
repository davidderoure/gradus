"""
Generate and verify the superparticular companion sequence:
a(n) = T(n, n-1) = n + Omega*(n-1) for n = 2, 3, 4, ...

Superparticular (or epimoric) ratios n/(n-1) are the consecutive-harmonic
intervals: octave (2/1), fifth (3/2), fourth (4/3), major third (5/4), ...
gcd(n, n-1) = 1 always, so no coprimality check needed.
"""

from math import gcd
from sympy import factorint

def omega_star(n):
    if n <= 1: return 0
    return sum(e * (p - 1) for p, e in factorint(n).items())

def a(n):
    """a(n) = n + Omega*(n-1) for n >= 2."""
    assert n >= 2
    return n + omega_star(n - 1)

# ── Generate sequence ────────────────────────────────────────────────────────

terms = [(n, a(n)) for n in range(2, 302)]
values = [v for _, v in terms]

print("=" * 65)
print("Superparticular sequence a(n) = n + Omega*(n-1), n >= 2")
print("=" * 65)
print()
print("First 20 terms with interval names:")
names = {2:"Octave",3:"Fifth",4:"Fourth",5:"Maj 3rd",6:"Min 3rd",
         7:"Sep min 3rd",8:"Maj 2nd",9:"Min whole",10:"Maj whole",
         11:"Undecimal",12:"Min 2nd (just)"}
for n, v in terms[:15]:
    name = names.get(n, f"{n}/{n-1}")
    print(f"  n={n:3d}: a({n}) = {n} + Omega*({n-1}) = {n} + {omega_star(n-1)} = {v:3d}  [{name}]")

print()
print("First 80 terms (for OEIS sequence field):")
print(", ".join(map(str, values[:80])))
print()

# ── b-file ───────────────────────────────────────────────────────────────────

with open("/Users/davidderoure/gradus/b397105.txt", "w") as f:
    f.write("# A397105 a(n) = n + Omega*(n-1), n >= 2. Two-stage dissonance of superparticular interval n/(n-1).\n")
    for n, v in terms:
        f.write(f"{n} {v}\n")

print(f"b-file written: {len(terms)} terms (n=2..301)")
print()

# ── Spot checks ──────────────────────────────────────────────────────────────

checks = [(2,2),(3,4),(4,6),(5,7),(6,10),(7,10),(8,14),(9,12),(12,22),(16,22)]
print("Spot-checks:")
all_ok = True
for n, expected in checks:
    got = a(n)
    ok = "✓" if got == expected else f"✗ GOT {got}"
    if got != expected: all_ok = False
    print(f"  a({n:2d}) = {n} + Omega*({n-1}) = {n} + {omega_star(n-1)} = {got:3d}  {ok}")
print(f"All correct: {all_ok}")
print()

# ── Properties ───────────────────────────────────────────────────────────────

print("=" * 65)
print("PROPERTIES")
print("=" * 65)

# Is it monotone?
non_mono = [(n, values[i-1], v) for i,(n,v) in enumerate(terms)
            if i > 0 and v < values[i-1]][:10]
print(f"\nNon-monotone steps (a(n) < a(n-1)), first 10:")
for n, prev, cur in non_mono:
    print(f"  a({n})={cur} < a({n-1})={prev}  [{n}/{n-1} vs {n-1}/{n-2}]")

# Relation to Omega*(n-1) alone
from scipy.stats import spearmanr
import numpy as np
ns = np.array([n for n,_ in terms[:50]])
av = np.array([v for _,v in terms[:50]])
on = np.array([omega_star(n-1) for n,_ in terms[:50]])
rho_full, _ = spearmanr(ns, av)
print(f"\nSpearman rho of a(n) vs n: {rho_full:.4f}  (n dominates)")

# ── OEIS draft ───────────────────────────────────────────────────────────────

print()
print("=" * 65)
print("OEIS SUBMISSION DRAFT")
print("=" * 65)

draft = f"""
NAME
  a(n) = n + A275314(n-1) - 1 for n >= 2; two-stage dissonance of the
  superparticular musical interval n/(n-1).

SEQUENCE (first 40 terms, offset 2)
  {", ".join(map(str, values[:40]))}

OFFSET
  2, 1

COMMENTS
  The superparticular (epimoric) ratio n/(n-1) is the musical interval
  between consecutive partials n and n-1 of a harmonic series: the octave
  (2/1), perfect fifth (3/2), perfect fourth (4/3), major third (5/4),
  minor third (6/5), and so on.

  a(n) is the value of the two-stage dissonance triangle T(n,k) (A397104)
  along its rightmost populated diagonal, k = n-1.

  a(n) = n + Omega*(n-1), where Omega*(m) = Sum_{{p^e || m}} e*(p-1)
  is Euler's weighted prime-omega function (A275314(m) - 1).

  Since gcd(n, n-1) = 1 for all n, no coprimality condition is needed.

  Two-stage perceptual interpretation: to perceive the interval n/(n-1),
  the auditory system must (Stage 1) extrapolate from the lower note
  (partial n-1 of an implied fundamental) down to that fundamental,
  at cost Omega*(n-1); and (Stage 2) confirm the upper note as the
  n-th partial of the fundamental, at cost n. Total: a(n) = n + Omega*(n-1).
  The brain does not literally factorise integers; it performs subharmonic
  template matching or autocorrelation across candidate periods (cf. Parncutt
  1988, Terhardt 1979). Omega*(n-1) acts as an analytical proxy for the
  complexity of that template search.

  The sequence is not monotone: a(8) = 14 > a(9) = 12 because
  Omega*(8) = Omega*(2^3) = 3, making 9/8 (major second) cheaper than
  8/7 (septimal major second) despite having a larger n. In general,
  when n-1 is a highly composite number (many small prime factors),
  a(n) dips below neighboring values.

  Comparison with Euler's Gradus for superparticular ratios:
  Euler's G(n/(n-1)) = 1 + Omega*(n) + Omega*(n-1); the present formula
  replaces Omega*(n) with n, consistently giving larger values and a wider
  dynamic range (Tritone 45/32 gives a(45)=50 vs Gradus=14).

FORMULA
  a(n) = n + Omega*(n-1) = n + A001414(n-1) - A001222(n-1).
  a(n) = T(n, n-1) where T is the two-stage dissonance triangle, A397104.
  a(n) = n + A275314(n-1) - 1.

EXAMPLE
  a(2) = 2: octave 2/1. Upper = P2 (cost 2), Omega*(1) = 0. Total 2.
  a(3) = 4: fifth 3/2. Upper = P3 (cost 3), Omega*(2) = 1. Total 4.
  a(4) = 6: fourth 4/3. Upper = P4 (cost 4), Omega*(3) = 2. Total 6.
  a(5) = 7: major third 5/4. Upper = P5 (cost 5), Omega*(4) = 2. Total 7.
  a(6) = 10: minor third 6/5. Upper = P6 (cost 6), Omega*(5) = 4. Total 10.
  a(9) = 12: major second 9/8. Upper = P9 (cost 9), Omega*(8) = 3. Total 12.
    Note a(9) < a(8) = 14 because 8 = 2^3 has low Omega* despite being large.

MATHEMATICA
  OmegaStar[1] = 0;
  OmegaStar[n_Integer /; n > 1] :=
    Total[#[[2]] (#[[1]] - 1) & /@ FactorInteger[n]];
  a[n_] := n + OmegaStar[n - 1];
  Table[a[n], {{n, 2, 50}}]

PARI/GP
  OmegaStar(n) = {{
    if(n <= 1, return(0));
    my(f = factor(n));
    sum(i = 1, #f~, f[i, 2] * (f[i, 1] - 1))
  }};
  a(n) = n + OmegaStar(n - 1);
  vector(50, n, a(n + 1))

PYTHON
  from sympy import factorint
  def omega_star(n):
      if n <= 1: return 0
      return sum(e*(p-1) for p, e in factorint(n).items())
  def a(n): return n + omega_star(n - 1)
  print([a(n) for n in range(2, 51)])

CROSSREFS
  A397104: rightmost diagonal of the two-stage dissonance triangle T(n,k) (companion entry).
  A275314: Euler's Gradus Suavitatis; a(n) = n + A275314(n-1) - 1.
  A001414: sopfr(n); a(n) = n + A001414(n-1) - A001222(n-1).
  A001222: bigomega(n).
  A007188: Tenney height of superparticular ratio n/(n-1) = n*(n-1).
  Cf. A002265 (superparticular ratios), A054519.

LINKS
  David De Roure, An Asymmetric Formula for Interval Consonance and its
  Relation to Harmonic Coincidence, arXiv:2606.16412.

KEYWORDS
  nonn, hear

AUTHOR
  David De Roure
"""
print(draft)
