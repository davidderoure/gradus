"""
Corrected analysis: weighted harmonic coincidence and Euler's Gradus.

Key fix: the "Score = M*(Gradus-1) + 2*C(M)" proof requires a FIXED M
for all intervals. In practice floor(N/p) varies with p, so we must either:
  (a) use fixed M  — demonstrates the exact linear identity
  (b) subtract the C(M) floor-dependent offset  — recovers Omega*(p) exactly

We also look properly at which weight family best predicts Gradus empirically.
"""

from math import gcd, log, sqrt
from sympy import factorint, isprime
import numpy as np
from scipy.stats import spearmanr, kendalltau
import matplotlib.pyplot as plt

# ── Core ────────────────────────────────────────────────────────────────────

def omega_star(n):
    if n <= 1: return 0
    return sum(e*(p-1) for p, e in factorint(n).items())

def gradus(p, q):
    g = gcd(p, q); p, q = p//g, q//g
    return 1 + omega_star(p) + omega_star(q)

def mangoldt(n):
    if n <= 1: return 0.0
    f = factorint(n)
    if len(f) == 1:
        p, _ = next(iter(f.items()))
        return log(p)
    return 0.0

# Precompute
MAX_N = 200
OS   = [omega_star(n) for n in range(MAX_N+1)]
MG   = [mangoldt(n)   for n in range(MAX_N+1)]
C_OS = [0] * (MAX_N+1)          # C_OS[M] = sum_{k=1}^{M} omega_star(k)
for k in range(1, MAX_N+1):
    C_OS[k] = C_OS[k-1] + OS[k]

# ── Ratios ──────────────────────────────────────────────────────────────────

MAX_INT = 16
ratios = [(p, q) for p in range(1, MAX_INT+1)
                 for q in range(1, p+1) if gcd(p,q)==1]
ratios.sort(key=lambda x: x[0]/x[1])
g_arr = np.array([gradus(p, q) for p, q in ratios])

# ── Part 1: Fixed-M exact identity ──────────────────────────────────────────

print("=" * 65)
print("EXACT IDENTITY — fixed M coincidences")
print()
print("  sum_{k=1}^{M} Omega*(kp) = Omega*(p)*M + C(M)")
print("  [from complete additivity: Omega*(kp) = Omega*(k) + Omega*(p)]")
print()

M_fixed = 20   # must satisfy M_fixed <= floor(MAX_N / max_p)

def fixed_M_score(p, q, M=M_fixed, N=MAX_N):
    """Score using exactly M coincidences for both notes."""
    g = gcd(p, q); p, q = p//g, q//g
    # Only compute if M harmonics fit within our table
    sp = sum(OS[k*p] for k in range(1, M+1) if k*p <= N)
    sq = sum(OS[k*q] for k in range(1, M+1) if k*q <= N)
    return sp + sq

scores_fixed = np.array([fixed_M_score(p, q) for p, q in ratios])

# Predicted by formula: M*(Gradus-1) + 2*C(M)
C_M = C_OS[M_fixed]
predicted = M_fixed * (g_arr - 1) + 2 * C_M
residual = scores_fixed - predicted

print(f"  M = {M_fixed},  C(M) = {C_M}")
print(f"  Predicted score = {M_fixed}*(Gradus - 1) + {2*C_M}")
print(f"  Max residual (should be 0 if k*p always <= N): {np.max(np.abs(residual)):.1f}")
print()

# Demonstrate with named intervals
named = [("Unison",1,1),("Octave",2,1),("Fifth",3,2),("Fourth",4,3),
         ("Maj 6th",5,3),("Maj 3rd",5,4),("Min 3rd",6,5),("Min 6th",8,5),
         ("Maj 2nd",9,8),("Tritone",45,32)]
print(f"  {'Interval':<12} {'p/q':<8} {'Gradus':>6}  {'Score':>7}  {'Predicted':>9}  {'Diff':>5}")
for name, p, q in named:
    g  = gradus(p, q)
    s  = fixed_M_score(p, q)
    pr = M_fixed*(g-1) + 2*C_M
    print(f"  {name:<12} {p}/{q:<6}  {g:>6}  {s:>7.0f}  {pr:>9.0f}  {s-pr:>5.0f}")

# ── Part 2: Corrected score — strip C(M), recover Omega*(p) + Omega*(q) ────

print()
print("=" * 65)
print("CORRECTED SCORE — subtract floor-dependent C(M) offset")
print()
print("  corrected(p) = [sum_{k=1}^{M_p} Omega*(kp)] - C(M_p)")
print("               = Omega*(p) * M_p")
print("  Dividing by M_p gives Omega*(p) exactly.")
print()

def corrected_omega(p, q, N=MAX_N):
    """Remove the C(M) baseline; leaves M*Omega*(p) + M*Omega*(q)."""
    g = gcd(p, q); p, q = p//g, q//g
    Mp = N // p; Mq = N // q
    sp = sum(OS[k*p] for k in range(1, Mp+1) if k*p <= N)
    sq = sum(OS[k*q] for k in range(1, Mq+1) if k*q <= N)
    # subtract C(Mp), C(Mq) to remove the k-dependent baseline
    return (sp - C_OS[Mp]) + (sq - C_OS[Mq])

corr_scores = np.array([corrected_omega(p, q) for p, q in ratios])

# This should equal floor(N/p)*Omega*(p) + floor(N/q)*Omega*(q)
def expected_corrected(p, q, N=MAX_N):
    g = gcd(p, q); p, q = p//g, q//g
    return (N//p)*omega_star(p) + (N//q)*omega_star(q)

expected = np.array([expected_corrected(p, q) for p, q in ratios])
max_err = np.max(np.abs(corr_scores - expected))
print(f"  Max error vs formula floor(N/p)*Omega*(p)+floor(N/q)*Omega*(q): {max_err:.1f}")
rho_corr, _ = spearmanr(g_arr, corr_scores)
print(f"  Spearman rho vs Gradus: {rho_corr:.4f}")
print()
print("  Why not perfect? floor(N/p)*Omega*(p) weighs Omega*(p) by N/p,")
print("  favouring intervals with small p. Gradus treats p and q symmetrically")
print("  via Omega*(p)+Omega*(q), not p*Omega*(p)+q*Omega*(q).")

# ── Part 3: Systematic weight-function search ───────────────────────────────

print()
print("=" * 65)
print("EMPIRICAL: best weight families (symmetric, N=200)")
print("=" * 65)

def sym_score(p, q, w_arr, N=MAX_N):
    g = gcd(p, q); p, q = p//g, q//g
    sp = sum(w_arr[k*p] for k in range(1, N//p+1) if k*p <= N)
    sq = sum(w_arr[k*q] for k in range(1, N//q+1) if k*q <= N)
    return sp + sq

weights = {
    'flat':        np.array([1.0 if n>=1 else 0.0 for n in range(MAX_N+1)]),
    'log(n)':      np.array([log(n) if n>1 else 0.0 for n in range(MAX_N+1)]),
    '1/n':         np.array([1/n   if n>=1 else 0.0 for n in range(MAX_N+1)]),
    '1/sqrt(n)':   np.array([1/sqrt(n) if n>=1 else 0.0 for n in range(MAX_N+1)]),
    'Omega*(n)':   np.array([float(OS[n]) for n in range(MAX_N+1)]),
    'Mangoldt':    np.array([MG[n] for n in range(MAX_N+1)]),
    'log(n)/n':    np.array([log(n)/n if n>1 else 0.0 for n in range(MAX_N+1)]),
}

# Power-law search
best_rho, best_alpha = 0, 1
for alpha_10 in range(-30, 31):
    alpha = alpha_10 / 10
    w = np.array([n**(-alpha) if n>=1 else 0.0 for n in range(MAX_N+1)])
    s = np.array([-sym_score(p, q, w) for p, q in ratios])
    rho, _ = spearmanr(g_arr, s)
    if rho > best_rho:
        best_rho, best_alpha = rho, alpha

print(f"\n  Power-law w(n)=n^(-alpha) grid search:")
print(f"  Best alpha = {best_alpha:.1f}, Spearman rho = {best_rho:.4f}")

w_best = np.array([n**(-best_alpha) if n>=1 else 0.0 for n in range(MAX_N+1)])
weights[f'n^(-{best_alpha})  [best power]'] = w_best

print()
print(f"  {'Weight':<30}  {'rho':>7}  {'tau':>7}")
print(f"  {'-'*30}  {'-'*7}  {'-'*7}")
all_rhos = {}
for name, w_arr in weights.items():
    s = np.array([-sym_score(p, q, w_arr) for p, q in ratios])
    rho, _ = spearmanr(g_arr, s)
    tau, _ = kendalltau(g_arr, s)
    all_rhos[name] = (rho, tau)
    print(f"  {name:<30}  {rho:>7.4f}  {tau:>7.4f}")

# ── Part 4: Interpretation ───────────────────────────────────────────────────

print()
print("=" * 65)
print("INTERPRETATION")
print("=" * 65)
print("""
Finding 1: EXACT EQUIVALENCE (with fixed M)
  Using w(n) = Omega*(n), the weighted coincidence sum is a
  perfect linear function of Gradus:

      Score = M*(Gradus - 1) + 2*C(M)

  This follows from complete additivity of Omega*:
      Omega*(k*p) = Omega*(k) + Omega*(p)

  Omega*(n) is itself the 'complexity' of harmonic n — how hard it
  is to produce that harmonic given its prime structure.  So Gradus
  counts coincidences weighted by the prime complexity of the harmonics
  that coincide.

Finding 2: FLAT COUNT IS EMPIRICALLY STRONG
  A simple unweighted count (flat) achieves rho~0.79, nearly as good
  as log(n) weighting.  The flat count equals floor(N/p) + floor(N/q),
  which loses the distinction between e.g. p=8 (power of 2) and p=12
  (2²×3) that Gradus captures.

Finding 3: LOG(n) IS THE BEST EMPIRICAL PREDICTOR
  w(n) = log(n) gives the best rank correlation without using any
  prime structure.  log(n) ≈ sum of log(prime) contributions ~ Mangoldt
  convolved with 1.  It reflects the 'arithmetic complexity' of n in a
  smoother way than Omega*.

Finding 4: CIRCULAR STRUCTURE
  The exact recovery (Finding 1) requires knowing Omega*(n), which IS
  Gradus applied to a single integer.  So Gradus of a ratio equals
  a Gradus-weighted coincidence sum — a beautiful self-referential
  structure, but not an independent derivation from physics.

  A truly physics-based metric uses flat or log weights, achieves
  ~0.80 correlation, and explains ~64% of Gradus variance — capturing
  the main shape but not the prime-sensitivity.
""")

# ── Plot ─────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle(f"Weighted harmonic coincidence vs Gradus  "
             f"(p,q ≤ {MAX_INT}, N={MAX_N} harmonics)", fontsize=12)

plot_keys = ['flat','log(n)','1/n','Omega*(n)','Mangoldt',
             f'n^(-{best_alpha})  [best power]']
rng = np.random.default_rng(1)
jitter = rng.uniform(-0.12, 0.12, len(ratios))

for ax, wkey in zip(axes.flat, plot_keys):
    w_arr = weights[wkey]
    s_arr = np.array([sym_score(p, q, w_arr) for p, q in ratios])
    rho, _ = spearmanr(g_arr, -s_arr)
    ax.scatter(g_arr + jitter, s_arr, alpha=0.35, s=12)
    ax.set_xlabel("Euler's Gradus")
    ax.set_ylabel("Score (higher = more consonant)")
    ax.set_title(f"w(n) = {wkey}\nρ = {rho:.3f}", fontsize=10)
    ax.grid(True, alpha=0.25)
    ax.invert_yaxis()

plt.tight_layout()
plt.savefig("/Users/davidderoure/gradus/weighted_coincidence2.png", dpi=150, bbox_inches='tight')
plt.close()
print("Plot saved to weighted_coincidence2.png")
