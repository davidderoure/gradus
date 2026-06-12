"""
Explore whether harmonic overlap predicts Euler's Gradus Suavitatis.

For two notes at frequency ratio p:q (integers, coprime), harmonics of note 1
are at multiples of 1, harmonics of note 2 are at multiples of p/q. An exact
coincidence occurs when n1 = n2 * p/q, i.e. n1*q = n2*p. Since gcd(p,q)=1,
this requires p | n1 and q | n2, so coincidences fall at harmonics p, 2p, 3p…
of note 1. Within N harmonics of note 1, there are floor(N/p) such coincidences.

Euler's Gradus for a ratio expressed as n = p*q (the product of numerator and
denominator in lowest terms) is: G(n) = 1 + sum_i e_i*(prime_i - 1)

OEIS A275314: https://oeis.org/A275314
"""

from math import gcd
from fractions import Fraction
from sympy import factorint
import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Euler's Gradus Suavitatis
# ---------------------------------------------------------------------------

def gradus(p: int, q: int) -> int:
    """Euler's Gradus for frequency ratio p:q (need not be coprime)."""
    g = gcd(p, q)
    p, q = p // g, q // g
    n = p * q  # characteristic number
    result = 1
    for prime, exp in factorint(n).items():
        result += exp * (prime - 1)
    return result


# ---------------------------------------------------------------------------
# Harmonic overlap metrics
# ---------------------------------------------------------------------------

def exact_coincidences(p: int, q: int, N: int) -> int:
    """
    Count harmonics of note 1 (up to N) that exactly coincide with a harmonic
    of note 2, for ratio p:q in lowest terms.

    Coincidences occur at n1 = p, 2p, 3p, ... so count = floor(N/p).
    By symmetry, if we count from note 2's side it's floor(N/q).
    """
    g = gcd(p, q)
    p, q = p // g, q // g
    return N // p  # harmonics of note 1 that match a harmonic of note 2


def near_coincidences(p: int, q: int, N: int, tolerance: float = 0.01) -> int:
    """
    Count near-coincidences using a tolerance in semitones (as fraction of
    semitone width). Useful to model finite-bandwidth ear perception.
    """
    g = gcd(p, q)
    p, q = p // g, q // g
    # Harmonics of note 1: 1, 2, ..., N  (at relative freq k)
    # Harmonics of note 2: p/q, 2p/q, ..., N*p/q
    harmonics1 = np.arange(1, N + 1, dtype=float)
    harmonics2 = np.arange(1, N + 1, dtype=float) * p / q

    count = 0
    for h2 in harmonics2:
        if h2 > N:
            break
        # nearest harmonic of note 1
        nearest = round(h2)
        if nearest >= 1 and abs(h2 - nearest) / nearest <= tolerance:
            count += 1
    return count


def overlap_score(p: int, q: int, N: int) -> float:
    """
    Roughness-inspired score: sum over all harmonic pairs of a proximity
    function. Closer pairs contribute more weight. This mimics Helmholtz/Plomp
    roughness models where nearby partials create beating dissonance.
    """
    g = gcd(p, q)
    p, q = p // g, q // g
    harmonics1 = np.arange(1, N + 1, dtype=float)
    harmonics2 = np.arange(1, N + 1, dtype=float) * p / q

    score = 0.0
    for h2 in harmonics2:
        if h2 > N * 1.1:
            break
        # Gaussian proximity to each harmonic of note 1
        diffs = np.abs(harmonics1 - h2) / harmonics1  # relative distance
        score += np.sum(np.exp(-diffs**2 / (2 * 0.02**2)))
    return score


# ---------------------------------------------------------------------------
# Survey of musical intervals
# ---------------------------------------------------------------------------

INTERVALS = [
    ("Unison",           1, 1),
    ("Octave",           2, 1),
    ("Perfect fifth",    3, 2),
    ("Perfect fourth",   4, 3),
    ("Major sixth",      5, 3),
    ("Major third",      5, 4),
    ("Minor third",      6, 5),
    ("Minor sixth",      8, 5),
    ("Major second",     9, 8),
    ("Minor seventh",    9, 5),
    ("Major seventh",   15, 8),
    ("Minor second",    16, 15),
    ("Tritone",         45, 32),
    ("Aug second",      75, 64),
]

N_HARMONICS = 64

print(f"{'Interval':<20} {'Ratio':>8}  {'Gradus':>6}  {'ExactN':>7}  {'NearN':>7}  {'Overlap':>9}")
print("-" * 65)

for name, p, q in INTERVALS:
    g_val   = gradus(p, q)
    exact   = exact_coincidences(p, q, N_HARMONICS)
    near    = near_coincidences(p, q, N_HARMONICS, tolerance=0.005)
    olap    = overlap_score(p, q, N_HARMONICS)
    print(f"{name:<20} {p}/{q:>5}  {g_val:>6}  {exact:>7}  {near:>7}  {olap:>9.2f}")


# ---------------------------------------------------------------------------
# Scatter plot: Gradus vs overlap metrics
# ---------------------------------------------------------------------------

gradus_vals  = [gradus(p, q)                         for _, p, q in INTERVALS]
exact_vals   = [exact_coincidences(p, q, N_HARMONICS) for _, p, q in INTERVALS]
overlap_vals = [overlap_score(p, q, N_HARMONICS)      for _, p, q in INTERVALS]
names        = [n for n, _, _ in INTERVALS]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle("Euler's Gradus vs Harmonic Overlap metrics\n"
             f"(first {N_HARMONICS} harmonics)", fontsize=13)

ax = axes[0]
ax.scatter(gradus_vals, exact_vals, zorder=3)
for i, name in enumerate(names):
    ax.annotate(name, (gradus_vals[i], exact_vals[i]),
                textcoords="offset points", xytext=(4, 2), fontsize=7)
ax.set_xlabel("Euler's Gradus (higher = more dissonant)")
ax.set_ylabel(f"Exact harmonic coincidences (of {N_HARMONICS})")
ax.set_title("Exact coincidences")
ax.grid(True, alpha=0.3)

ax = axes[1]
ax.scatter(gradus_vals, overlap_vals, color="darkorange", zorder=3)
for i, name in enumerate(names):
    ax.annotate(name, (gradus_vals[i], overlap_vals[i]),
                textcoords="offset points", xytext=(4, 2), fontsize=7)
ax.set_xlabel("Euler's Gradus (higher = more dissonant)")
ax.set_ylabel("Gaussian overlap score")
ax.set_title("Gaussian proximity score")
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("harmonic_overlap.png", dpi=150)
plt.show()
print("\nPlot saved to harmonic_overlap.png")


# ---------------------------------------------------------------------------
# Exact relationship for exact coincidences
# ---------------------------------------------------------------------------

print("\nExact formula analysis:")
print("For ratio p:q (coprime), exact coincidences in N harmonics = floor(N/p)")
print("Compare: Gradus = 1 + sum_i e_i*(prime_i-1) over factorisation of p*q")
print()
print("Ratios ranked by Gradus vs by coincidence count (most to least consonant):")
ranked_g = sorted(INTERVALS, key=lambda x: gradus(x[1], x[2]))
ranked_c = sorted(INTERVALS, key=lambda x: -exact_coincidences(x[1], x[2], N_HARMONICS))
print(f"{'By Gradus':<25}  {'By coincidences':<25}")
for (ng, _, _), (nc, _, _) in zip(ranked_g, ranked_c):
    print(f"  {ng:<23}    {nc:<23}")
