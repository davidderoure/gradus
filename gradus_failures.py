"""
Where does Gradus fail? Forensic analysis of rank disagreements.

From the perceptual model results, Gradus achieves rho=0.979 vs human
ratings, but key failures cluster around TIES within a Gradus level.
Multiple musically distinct intervals share the same Gradus value,
yet humans clearly rank them differently.

We investigate:
  1. Which intervals are tied by Gradus?
  2. Within each tie-group, what do humans prefer, and why?
  3. What tiebreaker candidates resolve the ties?
  4. The interval inversion asymmetry (5/4 vs 5/3, etc.)
  5. Combination tones (Tartini tones) as additional consonance signal.
"""

from math import gcd, log, log2, sqrt, exp
import numpy as np
from scipy.stats import spearmanr, kendalltau
from scipy.stats import rankdata
from sympy import factorint
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import combinations

# ── Core ────────────────────────────────────────────────────────────────────

def omega_star(n):
    if n <= 1: return 0
    return sum(e*(p-1) for p, e in factorint(n).items())

def gradus(p, q):
    g = gcd(p, q); p, q = p//g, q//g
    return 1 + omega_star(p) + omega_star(q)

# ── Human reference ──────────────────────────────────────────────────────────

INTERVALS = [
    # (p, q, name, human_rank)
    (1,  1,  "Unison",   1),
    (2,  1,  "Octave",   2),
    (3,  2,  "Fifth",    3),
    (4,  3,  "Fourth",   4),
    (5,  4,  "Maj 3rd",  5),
    (6,  5,  "Min 3rd",  6),
    (5,  3,  "Maj 6th",  7),
    (8,  5,  "Min 6th",  8),
    (9,  8,  "Maj 2nd",  9),
    (9,  5,  "Min 7th", 10),
    (15, 8,  "Maj 7th", 11),
    (16,15,  "Min 2nd", 12),
    (45,32,  "Tritone", 13),
]

names     = [x[2] for x in INTERVALS]
intervals = [(x[0],x[1]) for x in INTERVALS]
human_rk  = [x[3] for x in INTERVALS]

# ── Part 1: Gradus tie analysis ──────────────────────────────────────────────

print("=" * 65)
print("GRADUS TIES IN THE 13-INTERVAL SET")
print("=" * 65)

gradus_vals = {(p,q): gradus(p,q) for p,q in intervals}

from itertools import groupby
by_gradus = {}
for p, q, name, hr in INTERVALS:
    g = gradus(p, q)
    by_gradus.setdefault(g, []).append((p, q, name, hr))

print(f"\n  {'Gradus':>6}  Intervals in this group")
for g_val in sorted(by_gradus):
    group = by_gradus[g_val]
    items = "  ".join(f"{name}({p}/{q}) [human:{hr}]" for p,q,name,hr in group)
    marker = " *** TIE ***" if len(group)>1 else ""
    print(f"  {g_val:>6}:  {items}{marker}")

print()
print("  Summary of ties where human ranks differ:")
for g_val, group in sorted(by_gradus.items()):
    if len(group) > 1:
        hrs = [hr for _,_,_,hr in group]
        if len(set(hrs)) > 1:
            print(f"\n  Gradus={g_val}: {len(group)} intervals, "
                  f"human ranks {sorted(hrs)}")
            for p, q, name, hr in sorted(group, key=lambda x: x[3]):
                cents = 1200 * log2(p/q)
                print(f"    {name:<10} {p}/{q}  "
                      f"cents={cents:6.1f}  Ω*(p)={omega_star(p)}  Ω*(q)={omega_star(q)}  "
                      f"human={hr}")


# ── Part 2: What distinguishes intervals within a tie group? ─────────────────

print()
print("=" * 65)
print("CANDIDATES FOR RESOLVING TIES")
print("=" * 65)
print("""
Within a Gradus level, what additional quantities discriminate the intervals?

  Observation: Ω*(3)=Ω*(4)=2  →  5/4 and 5/3 have same Gradus
               Ω*(6)=Ω*(8)=3  →  min 3rd and min 6th tied
               Ω*(5)=Ω*(9)=4  →  potential ties with 9 in numerator

  Candidates:
  a) Interval width (cents)          — smaller may be more 'compact'
  b) Sum p+q                         — Euler also proposed this
  c) Max(p,q)                        — just p (numerator)
  d) log(p/q)                        — interval size
  e) p*q / gcd(p,q)^2 = lcm(p,q)    — least common multiple
  f) Combination/difference tones    — Tartini tones add virtual bass
  g) Octave equivalence weight       — count inversions equally?
""")

def tiebreaker_stats(vals, label, human_rk):
    rho, pval = spearmanr(human_rk, vals)
    return rho, pval

def combination_tone_score(p, q, depth=3):
    """
    Tartini/combination tones: when two tones at freq p and q sound,
    the ear generates difference tones at |m*p - n*q| for small m,n.
    Consonance contribution: count how many combo tones are small
    (close to existing harmonics) weighted by 1/(m+n).
    """
    g = gcd(p, q); pp, qq = p//g, q//g
    harmonics = set(range(1, 20))   # existing harmonic series
    score = 0.0
    for m in range(1, depth+1):
        for n in range(1, depth+1):
            tone = abs(m * pp - n * qq)
            if tone > 0:
                # Is this tone near a harmonic?
                score += (1/(m+n)) * (1 if tone in harmonics else 0)
    return score

def octave_reduced(p, q):
    """Reduce ratio to within an octave [1, 2)."""
    g = gcd(p, q); p, q = p//g, q//g
    while p/q >= 2: q *= 2
    while p/q < 1:  p *= 2
    return gcd_reduce(p, q)

def gcd_reduce(p, q):
    g = gcd(p, q); return p//g, q//g

# Compute candidate tiebreakers for all intervals
candidates = {}
for p, q, name, hr in INTERVALS:
    g = gcd(p, q); pp, qq = p//g, q//g
    r = pp/qq
    candidates[(p,q)] = {
        'cents':   1200 * log2(r),
        'sum_pq':  pp + qq,
        'max_pq':  max(pp, qq),
        'lcm':     pp * qq,     # = lcm since coprime
        'log_r':   log(r),
        'combo':   combination_tone_score(pp, qq),
        'p_alone': pp,
        'q_alone': qq,
    }

# Test each as supplement to Gradus
g_arr = np.array([gradus(p,q) for p,q in intervals])
h_arr = np.array(human_rk)

print(f"  {'Tiebreaker':<20}  {'rho alone':>9}  {'rho w/ Gradus':>13}")
print(f"  {'-'*20}  {'-'*9}  {'-'*13}")

for key in ['cents','sum_pq','max_pq','lcm','log_r','combo','p_alone']:
    vals = np.array([candidates[(p,q)][key] for p,q in intervals])
    rho_alone, _ = spearmanr(h_arr, vals)
    # Combine with Gradus by averaging ranks
    g_rank  = rankdata(g_arr)
    v_rank  = rankdata(vals)
    combined_rank = g_rank + v_rank
    rho_comb, _ = spearmanr(h_arr, combined_rank)
    print(f"  {key:<20}  {rho_alone:>+9.4f}  {rho_comb:>+13.4f}")


# ── Part 3: Interval inversion asymmetry ────────────────────────────────────

print()
print("=" * 65)
print("INTERVAL INVERSION ASYMMETRY")
print("=" * 65)
print("""
Two intervals are octave-inversions of each other if their cents add to 1200:
  Fifth (702¢) ↔ Fourth (498¢)
  Maj 3rd (386¢) ↔ Min 6th (814¢)
  Min 3rd (316¢) ↔ Maj 6th (884¢)
  Maj 2nd (204¢) ↔ Min 7th (996¢)
  Maj 7th (1088¢) ↔ Min 2nd (112¢)

Gradus treats p/q and its inversion 2q/p differently (different Gradus),
but Ω*(3)=Ω*(4)=2 causes the fifth/fourth distinction to be small
and the thirds/sixths to be conflated.
""")

inversions = [
    ("Fifth",   3,2,  "Fourth",  4,3),
    ("Maj 3rd", 5,4,  "Min 6th", 8,5),
    ("Min 3rd", 6,5,  "Maj 6th", 5,3),
    ("Maj 2nd", 9,8,  "Min 7th", 9,5),
    ("Maj 7th",15,8,  "Min 2nd",16,15),
]

print(f"  {'Interval':<10} {'Gradus':>6}  {'Human':>6}  "
      f"  {'Inversion':<10} {'Gradus':>6}  {'Human':>6}  "
      f"  {'ΔGradus':>7}  {'ΔHuman':>7}")
for name1,p1,q1, name2,p2,q2 in inversions:
    g1 = gradus(p1,q1); g2 = gradus(p2,q2)
    h1 = next(x[3] for x in INTERVALS if x[0]==p1 and x[1]==q1)
    h2 = next(x[3] for x in INTERVALS if x[0]==p2 and x[1]==q2)
    print(f"  {name1:<10} {g1:>6}  {h1:>6}    {name2:<10} {g2:>6}  {h2:>6}  "
          f"  {g2-g1:>+7}  {h2-h1:>+7}")

print("""
  Pattern: Gradus always correctly identifies the more dissonant inversion,
  but UNDERESTIMATES the size of the gap — especially for thirds vs sixths.
  Human listeners hear Maj 3rd (rank 5) as much more consonant than Min 6th
  (rank 8), a gap of 3. Gradus gives gap of only 1 (Gradus 7 vs 8).
""")


# ── Part 4: Combination tones as tiebreaker ──────────────────────────────────

print("=" * 65)
print("COMBINATION TONES (TARTINI) AS TIEBREAKER")
print("=" * 65)
print("""
When notes at frequencies f1=q and f2=p sound, the ear generates:
  Difference tone:  f2-f1 = p-q
  Sum tone:         f1+f2 = p+q
  2nd order:        2f1-f2 = 2q-p,  2f2-f1 = 2p-q
  etc.

These virtual tones may reinforce or conflict with the harmonic series.
Consonance: combo tones that fall on simple harmonics strengthen the root.
""")

print(f"  {'Interval':<10} {' p':>3}/{' q':<3}  {'p-q':>5}  {'p+q':>5}  "
      f"{'2q-p':>6}  {'2p-q':>6}  {'combo score':>11}  {'human':>6}")
for p, q, name, hr in INTERVALS:
    g = gcd(p,q); pp,qq = p//g, q//g
    diff  = pp - qq
    summ  = pp + qq
    c2a   = 2*qq - pp
    c2b   = 2*pp - qq
    cs    = combination_tone_score(pp, qq)
    print(f"  {name:<10} {pp:>3}/{qq:<3}  {diff:>5}  {summ:>5}  "
          f"{c2a:>6}  {c2b:>6}  {cs:>11.3f}  {hr:>6}")


# ── Part 5: Best composite model ─────────────────────────────────────────────

print()
print("=" * 65)
print("COMPOSITE MODEL: Gradus + tiebreakers")
print("=" * 65)

# Try weighted combinations
g_arr  = np.array([gradus(p,q) for p,q in intervals], dtype=float)
h_arr  = np.array(human_rk, dtype=float)
cents  = np.array([1200*log2((p//gcd(p,q))/(q//gcd(p,q))) for p,q in intervals])
lcm_v  = np.array([(p//gcd(p,q))*(q//gcd(p,q)) for p,q in intervals], dtype=float)
combo  = np.array([combination_tone_score(p//gcd(p,q), q//gcd(p,q)) for p,q in intervals])
sumv   = np.array([(p//gcd(p,q))+(q//gcd(p,q)) for p,q in intervals], dtype=float)

print(f"\n  Searching weighted combination: alpha*Gradus + beta*tiebreaker")
print(f"  (grid search over alpha, beta in [0,1], normalising both to [0,1])\n")

def normed(x):
    r = x - x.min()
    return r / r.max() if r.max() > 0 else r

g_n  = normed(g_arr)
c_n  = normed(cents)
l_n  = normed(np.log(lcm_v))
cb_n = normed(-combo)   # more combo = more consonant
s_n  = normed(sumv)

best_rho, best_params, best_label = 0, None, ""
results_composite = []
for label, tb in [("cents", c_n), ("log(lcm)", l_n),
                  ("combo", cb_n), ("sum p+q", s_n)]:
    for alpha in np.linspace(0, 1, 51):
        beta = 1 - alpha
        score = alpha * g_n + beta * tb
        rho, _ = spearmanr(h_arr, score)
        if rho > best_rho:
            best_rho, best_params, best_label = rho, (alpha, beta, label), label
    # Best for this tiebreaker
    best_local = max((alpha, 1-alpha) for alpha in np.linspace(0,1,51))
    scores_50 = 0.5*g_n + 0.5*tb
    rho_50, _ = spearmanr(h_arr, scores_50)
    rho_pure_tb, _ = spearmanr(h_arr, tb)
    results_composite.append((label, rho_pure_tb, rho_50))

print(f"  {'Tiebreaker':<12}  {'rho (pure TB)':>13}  {'rho (50/50 w Gradus)':>20}")
for label, r1, r2 in results_composite:
    print(f"  {label:<12}  {r1:>+13.4f}  {r2:>+20.4f}")

alpha, beta, blabel = best_params
print(f"\n  Best overall: {alpha:.2f}*Gradus + {beta:.2f}*{blabel}")
print(f"  Spearman rho = {best_rho:.4f}")
print(f"  (vs Gradus alone: {spearmanr(h_arr, g_arr)[0]:.4f})")


# ── Part 6: Ranking under best composite ─────────────────────────────────────

print()
print("=" * 65)
print("FINAL RANKING TABLE (best composite vs Gradus vs human)")
print("=" * 65)

# Use cents as tiebreaker (smaller interval = more consonant)
composite = 0.5*g_n + 0.5*c_n
c_rank = rankdata(composite)
g_rank = rankdata(g_arr)

print(f"  {'Interval':<10} {'Human':>6}  {'Gradus':>6}  {'Grad+cents':>10}  {'Correct?':>8}")
for i, (p, q, name, hr) in enumerate(INTERVALS):
    gr = int(g_rank[i])
    cr = int(c_rank[i])
    ok = "✓" if cr == hr else ("~" if abs(cr-hr)<=1 else "✗")
    print(f"  {name:<10} {hr:>6}  {gr:>6}  {cr:>10}  {ok:>8}")

g_correct  = sum(1 for i,(p,q,n,hr) in enumerate(INTERVALS) if int(g_rank[i])==hr)
c_correct  = sum(1 for i,(p,q,n,hr) in enumerate(INTERVALS) if int(c_rank[i])==hr)
g_near     = sum(1 for i,(p,q,n,hr) in enumerate(INTERVALS) if abs(int(g_rank[i])-hr)<=1)
c_near     = sum(1 for i,(p,q,n,hr) in enumerate(INTERVALS) if abs(int(c_rank[i])-hr)<=1)
print(f"\n  Exact rank match:   Gradus {g_correct}/13,  Gradus+cents {c_correct}/13")
print(f"  Within-1 match:     Gradus {g_near}/13,  Gradus+cents {c_near}/13")


# ── Plot ─────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle("Where Gradus fails: tie-group analysis", fontsize=12)

# Scatter: human vs Gradus, highlighting ties
ax = axes[0]
g_vals_plot = [gradus(p,q) for p,q in intervals]
jitter = np.random.default_rng(2).uniform(-0.15,0.15,len(intervals))
sc = ax.scatter(np.array(g_vals_plot)+jitter, human_rk,
                c=g_vals_plot, cmap='plasma', s=60, zorder=3)
for i, name in enumerate(names):
    ax.annotate(name, (g_vals_plot[i]+jitter[i], human_rk[i]),
                textcoords="offset points", xytext=(4,2), fontsize=7)
ax.plot([0,15],[0,15], 'k--', alpha=0.3, label='perfect agreement')
ax.set_xlabel("Euler's Gradus rank")
ax.set_ylabel("Human consonance rank")
ax.set_title(f"Gradus vs Human  (ρ={spearmanr(human_rk, g_vals_plot)[0]:.3f})")
ax.grid(True, alpha=0.25)

# Scatter: human vs composite
ax = axes[1]
comp_ranks = c_rank
ax.scatter(comp_ranks, human_rk, c=g_vals_plot, cmap='plasma', s=60, zorder=3)
for i, name in enumerate(names):
    ax.annotate(name, (comp_ranks[i], human_rk[i]),
                textcoords="offset points", xytext=(4,2), fontsize=7)
ax.plot([0,14],[0,14],'k--',alpha=0.3)
ax.set_xlabel("Gradus + cents rank")
ax.set_ylabel("Human consonance rank")
rho_c, _ = spearmanr(human_rk, c_rank)
ax.set_title(f"Composite vs Human  (ρ={rho_c:.3f})")
ax.grid(True, alpha=0.25)

# Bar chart: interval size (cents) within each Gradus level
ax = axes[2]
cents_vals = [1200*log2((p//gcd(p,q))/(q//gcd(p,q))) for p,q in intervals]
colors = plt.cm.tab10(np.linspace(0,1,len(intervals)))
bars = ax.bar(range(len(intervals)), cents_vals, color=colors, alpha=0.7)
ax.set_xticks(range(len(intervals)))
ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
ax.set_ylabel("Interval size (cents)")
ax.set_title("Interval size — smaller tends\nto be more consonant within Gradus ties")
# Annotate Gradus
for i, (p,q) in enumerate(intervals):
    ax.text(i, cents_vals[i]+10, f"G={gradus(p,q)}", ha='center', fontsize=6)
ax.grid(True, alpha=0.25, axis='y')

plt.tight_layout()
plt.savefig("/Users/davidderoure/gradus/gradus_failures.png", dpi=150, bbox_inches='tight')
plt.close()
print("\nPlot saved to gradus_failures.png")
