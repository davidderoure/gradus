"""
Mathematical analysis of the relationship between Euler's Gradus Suavitatis
and harmonic coincidence counts.

Key question: are they equivalent, or merely correlated?

For ratio p:q (coprime):
  - Coincidences of note 1's harmonics with note 2's: floor(N/p)
  - Coincidences of note 2's harmonics with note 1's: floor(N/q)
  - Shared harmonics (both notes): floor(N/(p*q))  [since lcm(p,q)=p*q when coprime]

  Gradus = 1 + Omega*(p*q)  where Omega*(n) = sum_i e_i*(prime_i - 1)
  Since Omega* is completely additive:  Gradus = 1 + Omega*(p) + Omega*(q)
"""

from math import gcd, log
from sympy import factorint
from itertools import product
import numpy as np
from scipy.stats import spearmanr, kendalltau
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def omega_star(n: int) -> int:
    """Weighted prime-omega: sum of e_i*(p_i - 1) over prime factorisation."""
    if n == 1:
        return 0
    return sum(exp * (prime - 1) for prime, exp in factorint(n).items())

def gradus(p: int, q: int) -> int:
    g = gcd(p, q); p, q = p//g, q//g
    return 1 + omega_star(p) + omega_star(q)

def coincidence_rate(p: int, q: int) -> float:
    """Asymptotic rate of coincidences per harmonic of note 1: 1/p."""
    g = gcd(p, q); p, q = p//g, q//g
    return 1.0 / p

def symmetric_rate(p: int, q: int) -> float:
    """Mean coincidence rate from both notes' perspectives: (1/p + 1/q)/2."""
    g = gcd(p, q); p, q = p//g, q//g
    return (1/p + 1/q) / 2

def joint_rate(p: int, q: int) -> float:
    """Rate of shared harmonics (both notes): 1/(p*q)."""
    g = gcd(p, q); p, q = p//g, q//g
    return 1.0 / (p * q)


# ---------------------------------------------------------------------------
# Enumerate all ratios p:q with p,q <= MAX_INT, p >= q (above unison)
# ---------------------------------------------------------------------------

MAX_INT = 16  # covers common just intonation intervals

ratios = []
for p, q in product(range(1, MAX_INT+1), range(1, MAX_INT+1)):
    if p >= q and gcd(p, q) == 1:
        ratios.append((p, q))

ratios.sort(key=lambda x: x[0]/x[1])  # sort by pitch ratio

# Compute metrics
data = []
for p, q in ratios:
    g_val   = gradus(p, q)
    c_rate  = coincidence_rate(p, q)
    s_rate  = symmetric_rate(p, q)
    j_rate  = joint_rate(p, q)
    char_n  = p * q  # Euler's characteristic number
    data.append({
        'p': p, 'q': q, 'ratio': p/q,
        'gradus': g_val,
        'coincidence_rate': c_rate,    # 1/p
        'symmetric_rate': s_rate,      # (1/p + 1/q)/2
        'joint_rate': j_rate,          # 1/(p*q)
        'log_char': log(p*q),          # log of characteristic number
        'omega_p': omega_star(p),
        'omega_q': omega_star(q),
    })

g_arr  = np.array([d['gradus']           for d in data])
c_arr  = np.array([d['coincidence_rate'] for d in data])
s_arr  = np.array([d['symmetric_rate']   for d in data])
j_arr  = np.array([d['joint_rate']       for d in data])
lc_arr = np.array([d['log_char']         for d in data])

print("=" * 70)
print("RANK CORRELATIONS with Euler's Gradus")
print("(Spearman rho and Kendall tau; negative = higher overlap => lower Gradus)")
print("=" * 70)

metrics = [
    ("1/p  (coincidence rate, note 1)",    -c_arr),
    ("(1/p+1/q)/2  (symmetric rate)",      -s_arr),
    ("1/(p*q)  (joint/shared rate)",       -j_arr),
    ("log(p*q)  (log characteristic n)",    lc_arr),
]
for label, arr in metrics:
    rho, p_rho = spearmanr(g_arr, arr)
    tau, p_tau = kendalltau(g_arr, arr)
    print(f"  {label}")
    print(f"    Spearman rho = {rho:.4f}  (p={p_rho:.2e})")
    print(f"    Kendall tau  = {tau:.4f}  (p={p_tau:.2e})")
    print()


# ---------------------------------------------------------------------------
# Where do they DISAGREE? Same Gradus, different coincidence rank (and vice versa)
# ---------------------------------------------------------------------------

print("=" * 70)
print("CASES WHERE GRADUS AND COINCIDENCE RATE DISAGREE")
print("(same Gradus, different coincidence rate — sorted by Gradus)")
print("=" * 70)

from itertools import groupby
sorted_data = sorted(data, key=lambda d: d['gradus'])
disagreements = []
for g_val, group in groupby(sorted_data, key=lambda d: d['gradus']):
    group = list(group)
    rates = [d['coincidence_rate'] for d in group]
    if len(set(rates)) > 1:  # same Gradus but different rates
        for d in group:
            disagreements.append(d)

if disagreements:
    print(f"{'Ratio':<10} {'Gradus':>6}  {'1/p rate':>9}  {'1/(pq) rate':>12}  {'Omega*(p)':>9}  {'Omega*(q)':>9}")
    prev_g = None
    for d in sorted(disagreements, key=lambda x: (x['gradus'], -x['coincidence_rate'])):
        if d['gradus'] != prev_g:
            print()
            prev_g = d['gradus']
        print(f"  {d['p']}/{d['q']:<7}  {d['gradus']:>6}  {d['coincidence_rate']:>9.4f}  "
              f"{d['joint_rate']:>12.5f}  {d['omega_p']:>9}  {d['omega_q']:>9}")
else:
    print("  None found in this range!")


# ---------------------------------------------------------------------------
# The algebra: why they're so similar but not identical
# ---------------------------------------------------------------------------

print()
print("=" * 70)
print("ALGEBRAIC RELATIONSHIP")
print("=" * 70)
print("""
Gradus(p/q) = 1 + Omega*(p) + Omega*(q)    [Omega* additive over p,q]

Coincidence rate from note 1: 1/p
  -> inversely related to p alone, ignores q's structure beyond its value

Gradus grows with Omega*(p) which is roughly log(p) but weights
large primes more heavily than repeated small primes:
  e.g.  p=8  = 2^3  -> Omega*(8) = 3*(2-1) = 3
        p=12 = 2^2*3 -> Omega*(12) = 2*(2-1)+1*(3-1) = 4  (larger!)
  But 1/8 < 1/12, so coincidence rate gives the OPPOSITE ordering here.

The metrics agree when Omega*(p) is monotone in p (which it roughly is),
but disagree at cases like these.
""")

# Show the 8 vs 12 case explicitly
for p, q in [(8,7),(8,5),(8,3),(12,11),(12,7),(12,5)]:
    if gcd(p,q)==1:
        print(f"  {p}/{q}: Gradus={gradus(p,q)}, 1/p={1/p:.4f}, "
              f"Omega*(p)={omega_star(p)}, Omega*(q)={omega_star(q)}")


# ---------------------------------------------------------------------------
# Plot: comprehensive scatter grid
# ---------------------------------------------------------------------------

fig = plt.figure(figsize=(14, 10))
fig.suptitle(f"Gradus vs Harmonic Overlap — all coprime ratios p/q, p,q ≤ {MAX_INT}\n"
             f"({len(data)} ratios)", fontsize=13)
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

plot_specs = [
    (gs[0,0], -c_arr, "−1/p  (coincidence rate, note 1)"),
    (gs[0,1], -s_arr, "−(1/p+1/q)/2  (symmetric rate)"),
    (gs[1,0], -j_arr, "−1/(p·q)  (joint rate)"),
    (gs[1,1], lc_arr, "log(p·q)  (log characteristic number)"),
]

for cell, y_arr, ylabel in plot_specs:
    ax = fig.add_subplot(cell)
    # jitter Gradus slightly to reveal overlapping points
    jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(g_arr))
    ax.scatter(g_arr + jitter, y_arr, alpha=0.4, s=15)
    rho, _ = spearmanr(g_arr, y_arr)
    ax.set_xlabel("Euler's Gradus")
    ax.set_ylabel(ylabel)
    ax.set_title(f"Spearman ρ = {rho:.3f}")
    ax.grid(True, alpha=0.25)

plt.savefig("/Users/davidderoure/gradus/gradus_analysis.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"\nPlot saved to gradus_analysis.png")


# ---------------------------------------------------------------------------
# Table: low-Gradus ratios with all metrics
# ---------------------------------------------------------------------------

print()
print("=" * 70)
print("LOW-GRADUS RATIOS — full metric table")
print("=" * 70)
print(f"{'Ratio':<8} {'Gradus':>6}  {'Ω*(p)':>6}  {'Ω*(q)':>6}  {'1/p':>7}  {'1/q':>7}  "
      f"{'sym':>8}  {'1/pq':>9}  {'log pq':>8}")
print("-" * 75)
for d in sorted(data, key=lambda x: (x['gradus'], -x['coincidence_rate'])):
    if d['gradus'] > 10:
        break
    print(f"  {d['p']}/{d['q']:<5} {d['gradus']:>6}  {d['omega_p']:>6}  {d['omega_q']:>6}  "
          f"{d['coincidence_rate']:>7.4f}  {1/d['q']:>7.4f}  "
          f"{d['symmetric_rate']:>8.4f}  {d['joint_rate']:>9.5f}  {d['log_char']:>8.3f}")
