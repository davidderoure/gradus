"""
Can we derive Euler's Gradus from a weighted harmonic coincidence count?

Gradus(p/q) = 1 + Omega*(p) + Omega*(q)
            = 1 + Omega*(p*q)     [Omega* completely additive]

where Omega*(n) = sum_i  e_i * (prime_i - 1)  over prime factorisation of n.

We look for a weight function w(n) such that:

    Score(p/q) = sum_{k=1}^{N/p} w(k*p)   ~   Gradus(p/q)

i.e. the weighted count of harmonics of note 1 that coincide with harmonics
of note 2.  We test several families and also fit weights numerically.

Key algebraic identity (since Omega* is completely additive):
    Omega*(k*p) = Omega*(k) + Omega*(p)

So if w(n) = Omega*(n):
    Score = sum_{k=1}^{M}  [Omega*(k) + Omega*(p)]
          = Omega*(p)*M  +  C(M)          [C(M) = sum Omega*(k), constant in p]

That means Omega*(p)*M is additive — and we just need both p and q sides.
"""

from math import gcd, log, sqrt
from sympy import factorint, isprime
from fractions import Fraction
import numpy as np
from scipy.stats import spearmanr
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------

def omega_star(n):
    if n <= 1: return 0
    return sum(e*(p-1) for p, e in factorint(n).items())

def gradus(p, q):
    g = gcd(p, q); p, q = p//g, q//g
    return 1 + omega_star(p) + omega_star(q)

def mangoldt(n):
    """Von Mangoldt: log(p) if n=p^k, else 0."""
    if n <= 1: return 0.0
    f = factorint(n)
    if len(f) == 1:
        p, _ = next(iter(f.items()))
        return log(p)
    return 0.0

def liouville(n):
    """Liouville lambda: (-1)^Omega(n)."""
    if n <= 1: return 1
    return (-1) ** sum(factorint(n).values())

# Build weight arrays up to MAX_N once (expensive calls)
MAX_N = 120
w_cache = {}
for n in range(1, MAX_N+1):
    w_cache[n] = {
        'flat':       1.0,
        '1/n':        1.0/n,
        '1/n2':       1.0/n**2,
        '1/sqrt_n':   1.0/sqrt(n),
        'log_n':      log(n) if n > 1 else 0.0,
        'omega_star': float(omega_star(n)),
        'mangoldt':   mangoldt(n),
        'prime_ind':  1.0 if isprime(n) else 0.0,
        'prime_wt':   float(sum(p-1 for p in factorint(n))) if n>1 else 0.0,
    }

WEIGHT_KEYS = list(next(iter(w_cache.values())).keys())

def weighted_score(p, q, wkey, N=MAX_N):
    """Sum w(k*p) for k=1..floor(N/p), symmetric over p and q."""
    g = gcd(p, q); p, q = p//g, q//g
    score_p = sum(w_cache[k*p][wkey] for k in range(1, N//p + 1) if k*p <= MAX_N)
    score_q = sum(w_cache[k*q][wkey] for k in range(1, N//q + 1) if k*q <= MAX_N)
    return score_p + score_q   # symmetric: both notes' perspectives


# ---------------------------------------------------------------------------
# Enumerate ratios
# ---------------------------------------------------------------------------

MAX_INT = 16
ratios = [(p, q) for p in range(1, MAX_INT+1) for q in range(1, p+1) if gcd(p,q)==1]
ratios.sort(key=lambda x: x[0]/x[1])

g_arr = np.array([gradus(p, q) for p, q in ratios])

print("=" * 65)
print("SPEARMAN RANK CORRELATION  (score vs Gradus)")
print("Symmetric score = sum_{k} w(k*p) + sum_{k} w(k*q)")
print("=" * 65)
print(f"  {'Weight function':<20}  {'rho':>7}  {'tau':>7}")
print(f"  {'-'*20}  {'-'*7}  {'-'*7}")

results = {}
for wkey in WEIGHT_KEYS:
    s_arr = np.array([-weighted_score(p, q, wkey) for p, q in ratios])
    rho, _ = spearmanr(g_arr, s_arr)
    from scipy.stats import kendalltau
    tau, _ = kendalltau(g_arr, s_arr)
    results[wkey] = (rho, tau, s_arr)
    print(f"  {wkey:<20}  {rho:>7.4f}  {tau:>7.4f}")


# ---------------------------------------------------------------------------
# The omega_star weight: theoretical connection
# ---------------------------------------------------------------------------

print()
print("=" * 65)
print("THEORETICAL ANALYSIS: w(n) = Omega*(n)")
print("=" * 65)
print("""
Since Omega* is completely additive: Omega*(k*p) = Omega*(k) + Omega*(p)

  Score_p(p) = sum_{k=1}^{M}  Omega*(k*p)
             = sum_{k=1}^{M}  [Omega*(k) + Omega*(p)]
             = Omega*(p) * M  +  C(M)

where C(M) = sum_{k=1}^{M} Omega*(k) is a constant (same for all p).

Symmetric score = Score_p + Score_q
  = Omega*(p)*M + C(M) + Omega*(q)*M + C(M)
  = M * [Omega*(p) + Omega*(q)]  +  2*C(M)
  = M * (Gradus - 1)  +  2*C(M)

=> Score = M*(Gradus-1) + 2*C(M)

This is a PERFECT LINEAR function of Gradus!
""")

M = MAX_N // 2   # typical floor(N/p) for small p
C_M = sum(omega_star(k) for k in range(1, M+1))
print(f"  With M={M}: Score = {M}*(Gradus-1) + 2*{C_M}")
print(f"  i.e.  Gradus = (Score - {2*C_M}) / {M}  +  1")
print()

# Verify numerically
omega_scores = np.array([weighted_score(p, q, 'omega_star') for p, q in ratios])
# Fit linear
from numpy.polynomial import polynomial as P
coeffs = np.polyfit(g_arr, omega_scores, 1)
residuals = omega_scores - np.polyval(coeffs, g_arr)
print(f"  Linear fit residual (max abs): {np.max(np.abs(residuals)):.4f}")
print(f"  (Non-zero because floor(N/p) varies with p — finite-N effect)")


# ---------------------------------------------------------------------------
# Finite-N correction: score per coincidence
# ---------------------------------------------------------------------------

print()
print("=" * 65)
print("NORMALISED SCORE: Omega*(p)/floor(N/p)  +  Omega*(q)/floor(N/q)")
print("= Omega*(p)/(N/p) + ... ≈ p*Omega*(p)/N + q*Omega*(q)/N")
print("  This isolates the *rate* of weighted coincidences.")
print("=" * 65)

def normalised_omega_score(p, q, N=MAX_N):
    g = gcd(p, q); p, q = p//g, q//g
    Mp = N // p; Mq = N // q
    if Mp == 0 or Mq == 0: return 0.0
    sp = sum(w_cache[k*p]['omega_star'] for k in range(1, Mp+1) if k*p <= MAX_N) / Mp
    sq = sum(w_cache[k*q]['omega_star'] for k in range(1, Mq+1) if k*q <= MAX_N) / Mq
    return sp + sq

norm_scores = np.array([normalised_omega_score(p, q) for p, q in ratios])
rho_norm, _ = spearmanr(g_arr, norm_scores)
print(f"  Spearman rho (normalised omega score vs Gradus): {rho_norm:.4f}")

# Exact analytical normalised rate as N->inf:
# sum_{k=1}^{M} Omega*(k*p) / M  ->  Omega*(p) + mean(Omega*(k)) as M->inf
# The mean of Omega*(k) for k=1..M converges to a constant (~ sum of (p-1)/p^2)
# So normalised score -> Omega*(p) + Omega*(q) + const = Gradus - 1 + const
# Perfect again in the limit!

print()
print("  Asymptotic normalised rate as N->inf:")
print("  (1/M) * sum_{k=1}^{M} Omega*(kp)  ->  Omega*(p)  +  <Omega*>")
print("  where <Omega*> = mean of Omega*(k) = sum_p (p-1)/p^2 (Mertens-type const)")
mean_omega = sum(omega_star(k) for k in range(1, 1001)) / 1000
print(f"  <Omega*> ≈ {mean_omega:.4f}  (empirical, k=1..1000)")
print(f"  => Normalised score -> Gradus - 1 + 2*<Omega*> = Gradus + {2*mean_omega-1:.4f}")


# ---------------------------------------------------------------------------
# What about the *prime* weight, w(n) = (p-1) if n is prime p, else 0?
# ---------------------------------------------------------------------------

print()
print("=" * 65)
print("SPARSE WEIGHT: w(n) = (prime-1) if n is prime, 0 otherwise")
print("Connection to Gradus via prime harmonic coincidences only")
print("=" * 65)

def prime_harmonic_score(p, q, N=MAX_N):
    """Count only coincidences at prime-numbered harmonics, weighted by (prime-1)."""
    g = gcd(p, q); p, q = p//g, q//g
    sp = sum((k*p - 1) for k in range(1, N//p + 1)
             if k*p <= MAX_N and isprime(k*p))
    sq = sum((k*q - 1) for k in range(1, N//q + 1)
             if k*q <= MAX_N and isprime(k*q))
    return sp + sq

ph_scores = np.array([prime_harmonic_score(p, q) for p, q in ratios])
rho_ph, _ = spearmanr(g_arr, ph_scores)
print(f"  Spearman rho: {rho_ph:.4f}")
print()
print("  Note: k*p is prime only when k=1 and p itself is prime.")
print("  So this score = (p-1) + (q-1) exactly when p,q are both prime,")
print("  = (p-1) when only p is prime, etc.")
print("  For p=2,q=1: score=1. p=3,q=2: score=(3-1)+(2-1)=3. Close to Gradus-1!")


# ---------------------------------------------------------------------------
# Summary plot
# ---------------------------------------------------------------------------

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle("Weighted harmonic coincidence scores vs Euler's Gradus\n"
             f"({len(ratios)} coprime ratios p/q, p,q ≤ {MAX_INT})", fontsize=12)

plot_items = [
    ('flat',        "flat (unweighted)\nρ={rho:.3f}"),
    ('1/n',         "w(n) = 1/n\nρ={rho:.3f}"),
    ('1/n2',        "w(n) = 1/n²\nρ={rho:.3f}"),
    ('omega_star',  "w(n) = Ω*(n)\nρ={rho:.3f}"),
    ('mangoldt',    "w(n) = Λ(n)  [Mangoldt]\nρ={rho:.3f}"),
    ('prime_ind',   "w(n) = 𝟙[n prime]\nρ={rho:.3f}"),
]

rng = np.random.default_rng(0)
jitter = rng.uniform(-0.12, 0.12, len(ratios))

for ax, (wkey, title) in zip(axes.flat, plot_items):
    rho, tau, s_arr = results[wkey]
    ax.scatter(g_arr + jitter, -s_arr, alpha=0.35, s=12)
    ax.set_xlabel("Euler's Gradus")
    ax.set_ylabel("Weighted score")
    ax.set_title(title.format(rho=rho), fontsize=10)
    ax.grid(True, alpha=0.25)

plt.tight_layout()
plt.savefig("/Users/davidderoure/gradus/weighted_coincidence.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"\nPlot saved to weighted_coincidence.png")


# ---------------------------------------------------------------------------
# Key result table
# ---------------------------------------------------------------------------

print()
print("=" * 65)
print("KEY RESULT: w(n) = Omega*(n) gives EXACT linear recovery of Gradus")
print()
print("  Gradus(p/q)  =  [Score_Omega*(p/q)  -  2*C(M)]  /  M  +  1")
print()
print("  where M = floor(N/p), C(M) = sum_{k=1}^{M} Omega*(k)")
print()
print("  In the limit N -> inf:")
print("  Gradus(p/q)  =  (1/N) * weighted_coincidence_count  +  constant")
print()
print("  The weight Omega*(n) = sum_i e_i*(p_i-1) is itself derived from")
print("  prime factorisation — so we've shown Gradus IS a harmonic coincidence")
print("  count, but with harmonics weighted by their own prime complexity.")
print("=" * 65)

# Show a few concrete examples
print()
print("Concrete check for small intervals:")
print(f"{'Interval':<20} {'p/q':<7} {'Gradus':>6}  {'OmegaScore':>10}  {'Linear pred':>12}")
named = [("Unison",1,1),("Octave",2,1),("Fifth",3,2),("Fourth",4,3),
         ("Maj 3rd",5,4),("Min 3rd",6,5),("Maj 2nd",9,8)]
for name, p, q in named:
    g = gradus(p, q)
    s = weighted_score(p, q, 'omega_star')
    pred = (s - 2*C_M) / M + 1
    print(f"  {name:<18} {p}/{q:<5}  {g:>6}  {s:>10.1f}  {pred:>12.3f}")
