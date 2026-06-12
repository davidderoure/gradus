"""
Search for a simple dissonance formula that outperforms Gradus on human ratings.

Strategy:
  1. Parameterised Gradus: optimise prime weights w(2), w(3), w(5)
     (allowing asymmetry between numerator and denominator)
  2. Hybrid: max(p,q) + correction terms
  3. Interpret the best formula musically

Human ranking (gold standard):
  Unison=1, Octave=2, Fifth=3, Fourth=4, Maj3rd=5, Min3rd=6,
  Maj6th=7, Min6th=8, Maj2nd=9, Min7th=10, Maj7th=11, Min2nd=12, Tritone=13
"""

from math import gcd, log2, log
from sympy import factorint
import numpy as np
from scipy.stats import spearmanr, kendalltau
from scipy.optimize import minimize, differential_evolution
from scipy.stats import rankdata
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Data ────────────────────────────────────────────────────────────────────

INTERVALS = [
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

intervals  = [(x[0],x[1]) for x in INTERVALS]
names      = [x[2] for x in INTERVALS]
human_rk   = np.array([x[3] for x in INTERVALS], dtype=float)

def r(p,q): g=gcd(p,q); return p//g, q//g

def omega_star(n):
    if n<=1: return 0
    return sum(e*(p-1) for p,e in factorint(n).items())

def gradus(p,q):
    pp,qq=r(p,q); return 1+omega_star(pp)+omega_star(qq)

# ── Prime decomposition helpers ──────────────────────────────────────────────

def prime_exponents(n, primes=(2,3,5,7)):
    f = factorint(n)
    return tuple(f.get(p,0) for p in primes)

# Precompute exponents for all p, q
PRIMES = (2, 3, 5)
pq_exponents = []
for p,q in intervals:
    pp,qq = r(p,q)
    ep = prime_exponents(pp, PRIMES)
    eq = prime_exponents(qq, PRIMES)
    pq_exponents.append((ep, eq))


# ── 1: Optimise prime weights ────────────────────────────────────────────────
# Gradus: score = 1 + sum_primes e_p*(p-1) for both numerator and denominator
# Here: score = sum_primes [ep_num * wp_num[p] + ep_den * wp_den[p]]
# Parameters: wp_num[2], wp_num[3], wp_num[5], wp_den[2], wp_den[3], wp_den[5]
# Gradus uses wp_num = wp_den = (1, 2, 4) and adds 1.

def custom_score(params, exponents=pq_exponents):
    # params: [w_num_2, w_num_3, w_num_5, w_den_2, w_den_3, w_den_5]
    scores = []
    for ep, eq in exponents:
        s = sum(ep[i]*params[i] for i in range(3)) + \
            sum(eq[i]*params[i+3] for i in range(3))
        scores.append(s)
    return np.array(scores)

def neg_spearman(params):
    scores = custom_score(params)
    rho, _ = spearmanr(human_rk, scores)
    return -rho

# Grid search: try many weight combinations
print("=" * 65)
print("1. OPTIMISED PRIME WEIGHTS (asymmetric numerator/denominator)")
print("=" * 65)

# Gradus weights: w_num = w_den = (1, 2, 4)
gradus_params = [1, 2, 4, 1, 2, 4]
gradus_scores = custom_score(gradus_params)
rho_gradus, _ = spearmanr(human_rk, gradus_scores)
print(f"\n  Gradus (w2=1,w3=2,w5=4 symmetric): rho={rho_gradus:.4f}")

# Optimise using differential evolution (global)
bounds = [(0, 10)]*6
result = differential_evolution(neg_spearman, bounds, seed=42, maxiter=500,
                                tol=1e-8, popsize=20)
opt_params = result.x
opt_scores = custom_score(opt_params)
rho_opt, _ = spearmanr(human_rk, opt_scores)

print(f"  Optimised (asymmetric):             rho={rho_opt:.4f}")
print(f"  Optimal weights:")
print(f"    Numerator:   w(2)={opt_params[0]:.3f}  w(3)={opt_params[1]:.3f}  w(5)={opt_params[2]:.3f}")
print(f"    Denominator: w(2)={opt_params[3]:.3f}  w(3)={opt_params[4]:.3f}  w(5)={opt_params[5]:.3f}")

# Try symmetric optimisation
def neg_spearman_sym(params):
    params6 = list(params)*2
    scores = custom_score(params6)
    rho, _ = spearmanr(human_rk, scores)
    return -rho

result_sym = differential_evolution(neg_spearman_sym, [(0,10)]*3, seed=42,
                                    maxiter=500, tol=1e-8, popsize=20)
sym_params = list(result_sym.x)*2
sym_scores = custom_score(sym_params)
rho_sym, _ = spearmanr(human_rk, sym_scores)
print(f"\n  Optimised (symmetric):              rho={rho_sym:.4f}")
print(f"  Symmetric weights: w(2)={result_sym.x[0]:.3f}  w(3)={result_sym.x[1]:.3f}  w(5)={result_sym.x[2]:.3f}")
print(f"  Ratio to Gradus:   w(2)/w(3)={result_sym.x[0]/result_sym.x[1]:.3f}  "
      f"w(3)/w(5)={result_sym.x[1]/result_sym.x[2]:.3f}")


# ── 2: Try simple candidate formulas ────────────────────────────────────────

print()
print("=" * 65)
print("2. SIMPLE CANDIDATE FORMULAS — exhaustive survey")
print("=" * 65)

def v2(n):
    """2-adic valuation."""
    if n==0: return float('inf')
    k=0
    while n%2==0: n//=2; k+=1
    return k

def odd_part(n):
    while n%2==0: n//=2
    return n

candidates = {}
for p,q in intervals:
    pp,qq = r(p,q)
    e_p = prime_exponents(pp); e_q = prime_exponents(qq)
    candidates[(p,q)] = {
        'Gradus':               1+omega_star(pp)+omega_star(qq),
        'max(p,q)':             max(pp,qq),
        'p+q':                  pp+qq,
        'p*q':                  pp*qq,
        'log2(p*q)':            log2(pp*qq),
        # Asymmetric
        'p + Ω*(q)':            pp + omega_star(qq),
        'Ω*(p) + q':            omega_star(pp) + qq,
        'p + 2*Ω*(q)':          pp + 2*omega_star(qq),
        # Hybrids
        'max + (p-q)':          max(pp,qq) + (pp-qq),
        'max + Ω*(min)':        max(pp,qq) + omega_star(min(pp,qq)),
        'max + Ω*(q)':          max(pp,qq) + omega_star(qq),
        # Two-component
        'p*odd_part(q)':        pp * odd_part(qq),
        'odd_p*odd_q':          odd_part(pp) * odd_part(qq),
        'odd_p + odd_q':        odd_part(pp) + odd_part(qq),
        # v2-based
        'p + v2(q)':            pp + v2(qq),
        'max - v2(q)':          max(pp,qq) - v2(qq),
        # Log-based
        'p + log2(q)':          pp + log2(qq) if qq>0 else pp,
        # Stern-Brocot
        'SB depth':             sum(factorint(pp*qq).values()),
        # Denominaor-weighted
        'p + q/2':              pp + qq/2,
        'p + q/3':              pp + qq/3,
    }

print(f"\n  {'Formula':<22}  {'rho':>7}  {'tau':>7}")
print(f"  {'-'*22}  {'-'*7}  {'-'*7}")
formula_rhos = {}
for key in candidates[intervals[0]]:
    vals = np.array([candidates[(p,q)][key] for p,q in intervals])
    rho,_ = spearmanr(human_rk, vals)
    tau,_ = kendalltau(human_rk, vals)
    formula_rhos[key] = rho
    print(f"  {key:<22}  {rho:>7.4f}  {tau:>7.4f}")


# ── 3: Focus on p + Ω*(q) — the asymmetric Gradus ──────────────────────────

print()
print("=" * 65)
print("3. BEST SIMPLE FORMULA: f(p,q) = p + Ω*(q)")
print("   (Galileo's coincidence measure + Euler's denominator complexity)")
print("=" * 65)

def f_asym(p,q):
    pp,qq=r(p,q); return pp + omega_star(qq)

asym_vals = np.array([f_asym(p,q) for p,q in intervals])
rho_asym,_=spearmanr(human_rk, asym_vals)
tau_asym,_=kendalltau(human_rk, asym_vals)

print(f"\n  Spearman rho: {rho_asym:.4f}")
print(f"  Kendall tau:  {tau_asym:.4f}")
print()
print(f"  {'Interval':<10} {'p/q':>7}  {'p':>3}  {'Ω*(q)':>6}  {'f=p+Ω*(q)':>10}  "
      f"{'Gradus':>7}  {'Human':>6}  {'Diff':>5}")
asym_ranks = rankdata(asym_vals)
g_ranks = rankdata([gradus(p,q) for p,q in intervals])
for i,(p,q,name,hr) in enumerate(INTERVALS):
    pp,qq=r(p,q)
    fval=f_asym(p,q)
    g=gradus(p,q)
    diff_from_human = asym_ranks[i]-hr
    print(f"  {name:<10} {pp}/{qq:<5}  {pp:>3}  {omega_star(qq):>6}  {fval:>10}  "
          f"{g:>7}  {hr:>6}  {diff_from_human:>+5.1f}")


# ── 4: Can we fix the remaining failures? ───────────────────────────────────

print()
print("=" * 65)
print("4. REMAINING FAILURES AND REFINED FORMULA")
print("=" * 65)

# Where does p + Ω*(q) fail?
print("\n  Failures of f(p,q) = p + Ω*(q):")
for i,(p,q,name,hr) in enumerate(INTERVALS):
    diff = abs(asym_ranks[i]-hr)
    if diff >= 1.5:
        pp,qq=r(p,q)
        print(f"    {name} ({pp}/{qq}): formula rank {asym_ranks[i]:.1f}, human {hr}, diff={diff:.1f}")

# Try parameterised version: f(p,q) = p + alpha*Ω*(q)
print("\n  Optimising: f(p,q) = p + α*Ω*(q)")
best_rho, best_alpha = 0, 1
for alpha_100 in range(0, 400):
    alpha = alpha_100 / 100
    vals = np.array([r(p,q)[0] + alpha*omega_star(r(p,q)[1]) for p,q in intervals])
    rho,_ = spearmanr(human_rk, vals)
    if rho > best_rho:
        best_rho, best_alpha = rho, alpha

print(f"  Best α = {best_alpha:.2f}, rho = {best_rho:.4f}")

# Show the α sweep
alphas = np.linspace(0, 4, 400)
rhos_alpha = []
for alpha in alphas:
    vals = np.array([r(p,q)[0] + alpha*omega_star(r(p,q)[1]) for p,q in intervals])
    rho,_ = spearmanr(human_rk, vals)
    rhos_alpha.append(rho)

# Also try: f(p,q) = p + Ω*(q) + β*v2(q)  — add octave correction for denominator
print("\n  Also try: f(p,q) = p + Ω*(q) + β*v2(q)  [restore 2's contribution]")
best_rho2, best_b = 0, 0
for b_100 in range(-100, 300):
    b = b_100 / 100
    vals = np.array([r(p,q)[0] + omega_star(r(p,q)[1]) + b*v2(r(p,q)[1])
                     for p,q in intervals])
    rho,_ = spearmanr(human_rk, vals)
    if rho > best_rho2:
        best_rho2, best_b = rho, b

print(f"  Best β = {best_b:.2f}, rho = {best_rho2:.4f}")

def f_refined(p,q, b=best_b):
    pp,qq=r(p,q); return pp + omega_star(qq) + b*v2(qq)

ref_vals = np.array([f_refined(p,q) for p,q in intervals])
rho_ref,_=spearmanr(human_rk, ref_vals)
print()
print(f"  f(p,q) = p + Ω*(q) + {best_b:.2f}*v₂(q)  [rho={rho_ref:.4f}]")
print()
ref_ranks = rankdata(ref_vals)
print(f"  {'Interval':<10} {'p/q':>7}  {'f value':>8}  {'rank':>6}  {'human':>6}  {'ok?':>5}")
for i,(p,q,name,hr) in enumerate(INTERVALS):
    pp,qq=r(p,q)
    fv=f_refined(p,q)
    rk=ref_ranks[i]
    ok = "✓" if rk==hr else ("~" if abs(rk-hr)<=1 else "✗")
    print(f"  {name:<10} {pp}/{qq:<5}  {fv:>8.2f}  {rk:>6.1f}  {hr:>6}  {ok:>5}")


# ── 5: The "double Gradus" — what if Gradus applied separately? ─────────────

print()
print("=" * 65)
print("5. INTERPRETING f(p,q) = p + Ω*(q)")
print("=" * 65)
print(f"""
  Euler's Gradus:  G(p/q) = 1 + Ω*(p) + Ω*(q)
  New formula:     f(p/q) =     p      + Ω*(q)

  These differ only in how the NUMERATOR is weighted:
    Gradus uses Ω*(p): prime-complexity of the upper partial
    New formula uses p: the partial number itself (Galileo's measure)

  The denominator treatment is the same: Ω*(q) in both.

  Musically: the upper note is judged by its ordinal position in the
  harmonic series (how high up it sits); the lower note is judged by
  the prime complexity of which partial of the fundamental it is.

  Why this asymmetry?
  — The upper note sets the 'reach' of the interval: higher partial =
    less likely to coincide, harder to tune, harder to hear as unified.
    Partial number p is the direct measure of this reach.
  — The lower note (the bass) sets harmonic context. A bass note that
    is a 'pure octave' (q=2^k, Ω*(q)=k) provides a simpler context
    than a bass note involving a large odd prime (Ω*(q) high).
    Ω*(q) measures this contextual complexity.

  The formula is: reach of upper note + complexity of bass context.
""")

# Show what changes between Gradus and new formula
print(f"  {'Interval':<10}  {'Ω*(p)':>6}  {'p':>6}  {'diff':>5}  {'impact'}")
for p,q,name,_ in INTERVALS:
    pp,qq=r(p,q)
    op=omega_star(pp)
    diff=pp-op
    impact=""
    if diff>0: impact=f"new formula rates MORE dissonant ({pp}>{op})"
    elif diff<0: impact=f"new formula rates LESS dissonant ({pp}<{op})"
    else: impact="same"
    print(f"  {name:<10}  {op:>6}  {pp:>6}  {diff:>+5}  {impact}")


# ── 6: Final comparison of top formulas ─────────────────────────────────────

print()
print("=" * 65)
print("6. FINAL RANKING: all top formulas vs human")
print("=" * 65)

formulas = {
    'Human':          human_rk,
    'Gradus':         np.array([gradus(p,q)        for p,q in intervals]),
    'max(p,q)':       np.array([max(*r(p,q))       for p,q in intervals]),
    'p + Ω*(q)':      np.array([f_asym(p,q)        for p,q in intervals]),
    f'p+Ω*(q)+{best_b:.1f}v₂(q)': np.array([f_refined(p,q) for p,q in intervals]),
    'p+q (sum)':      np.array([sum(r(p,q))        for p,q in intervals]),
}

print(f"\n  {'Interval':<10} " +
      "".join(f"  {k[:12]:>12}" for k in formulas))
print(f"  {'-'*10} " + "  ------------" * len(formulas))

ranked = {k: rankdata(v) for k,v in formulas.items()}
for i,(p,q,name,hr) in enumerate(INTERVALS):
    row = f"  {name:<10} "
    for k in formulas:
        row += f"  {ranked[k][i]:>12.1f}"
    print(row)

print()
print(f"  {'Formula':<25}  {'rho':>7}  {'tau':>7}  {'exact matches':>13}")
for k,vals in formulas.items():
    if k=='Human': continue
    rho,_=spearmanr(human_rk, vals)
    tau,_=kendalltau(human_rk, vals)
    rk = rankdata(vals)
    exact = sum(rk[i]==human_rk[i] for i in range(len(intervals)))
    print(f"  {k:<25}  {rho:>7.4f}  {tau:>7.4f}  {exact:>13}")


# ── Plot ─────────────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(16, 10))
fig.suptitle("Formula search: deriving a simple dissonance measure", fontsize=12)
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

# 1. α sweep for p + α*Ω*(q)
ax = fig.add_subplot(gs[0,0])
ax.plot(alphas, rhos_alpha, color='steelblue', lw=2)
ax.axvline(x=best_alpha, color='red', ls='--', alpha=0.7, label=f'best α={best_alpha:.2f}')
ax.axhline(y=rho_gradus, color='green', ls=':', alpha=0.7, label=f'Gradus ρ={rho_gradus:.3f}')
ax.axhline(y=formula_rhos['max(p,q)'], color='orange', ls=':', alpha=0.7,
           label=f'max(p,q) ρ={formula_rhos["max(p,q)"]:.3f}')
ax.set_xlabel("α")
ax.set_ylabel("Spearman ρ")
ax.set_title("f(p,q) = p + α·Ω*(q)", fontsize=10)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.25)

# 2. Ranking comparison scatter
ax = fig.add_subplot(gs[0,1])
jitter = np.random.default_rng(3).uniform(-0.15,0.15,len(intervals))
for formula, color, label in [
    ('Gradus', 'darkorange', 'Gradus'),
    ('p + Ω*(q)', 'steelblue', 'p+Ω*(q)')]:
    rk = rankdata(formulas[formula])
    ax.scatter(human_rk+jitter, rk, alpha=0.6, s=40, label=f'{label} ρ={spearmanr(human_rk,rk)[0]:.3f}')
ax.plot([1,13],[1,13],'k--',alpha=0.3)
ax.set_xlabel("Human rank")
ax.set_ylabel("Formula rank")
ax.set_title("Gradus vs p+Ω*(q) rankings", fontsize=10)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.25)

# 3. Formula values visualised as number line
ax = fig.add_subplot(gs[0,2])
for i,(p,q,name,hr) in enumerate(INTERVALS):
    pp,qq=r(p,q)
    g_val=gradus(p,q)
    a_val=f_asym(p,q)
    ax.scatter(g_val, i, c='darkorange', s=50, zorder=3)
    ax.scatter(a_val, i, c='steelblue', s=50, marker='D', zorder=3)
    ax.plot([g_val,a_val],[i,i],'gray',alpha=0.4,lw=1)
ax.set_yticks(range(len(INTERVALS)))
ax.set_yticklabels([f"{n} (h={hr})" for _,_,n,hr in INTERVALS], fontsize=8)
ax.set_xlabel("Dissonance value")
ax.set_title("Gradus (●) vs p+Ω*(q) (◆)", fontsize=10)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([],[],marker='o',color='darkorange',label='Gradus',ls=''),
                   Line2D([],[],marker='D',color='steelblue',label='p+Ω*(q)',ls='')],
          fontsize=8)
ax.grid(True, axis='x', alpha=0.25)

# 4. Numerator contribution: Ω*(p) vs p
ax = fig.add_subplot(gs[1,0])
pp_vals = [r(p,q)[0] for p,q in intervals]
op_vals = [omega_star(r(p,q)[0]) for p,q in intervals]
ax.scatter(pp_vals, op_vals, s=60, color='steelblue')
for i,(_,_,name,_) in enumerate(INTERVALS):
    ax.annotate(name, (pp_vals[i], op_vals[i]),
                textcoords="offset points", xytext=(4,2), fontsize=7)
# add y=x line (if they were equal)
m=max(pp_vals)
ax.plot([0,m],[0,m],'k--',alpha=0.2, label='Ω*(p)=p (never)')
ax.set_xlabel("p  (partial number of upper note)")
ax.set_ylabel("Ω*(p)  (Gradus complexity of upper note)")
ax.set_title("Numerator: p vs Ω*(p)\nGradus uses Ω*(p); new formula uses p", fontsize=10)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.25)

# 5. Denominator: Ω*(q) for each interval
ax = fig.add_subplot(gs[1,1])
qq_vals = [r(p,q)[1] for p,q in intervals]
oq_vals = [omega_star(r(p,q)[1]) for p,q in intervals]
ax.scatter(qq_vals, oq_vals, s=60, color='forestgreen')
for i,(_,_,name,_) in enumerate(INTERVALS):
    ax.annotate(name, (qq_vals[i], oq_vals[i]),
                textcoords="offset points", xytext=(4,2), fontsize=7)
ax.set_xlabel("q  (partial number of lower note)")
ax.set_ylabel("Ω*(q)  (Euler complexity of lower note)")
ax.set_title("Denominator: q vs Ω*(q)\nBoth formulas use Ω*(q) here", fontsize=10)
ax.grid(True, alpha=0.25)

# 6. Summary bar: all formula rhos
ax = fig.add_subplot(gs[1,2])
fkeys = ['Gradus','max(p,q)','p + Ω*(q)',f'p+Ω*(q)+{best_b:.1f}v₂(q)',
         'p+q (sum)','p*q']
frhos = [spearmanr(human_rk, formulas.get(k, np.array([candidates[(p,q)].get(k,0)
         for p,q in intervals])))[0] for k in fkeys]
# Recompute properly
frhos = []
for k in fkeys:
    if k in formulas:
        v=formulas[k]
    else:
        v=np.array([candidates[(p,q)].get(k,0) for p,q in intervals])
    frhos.append(spearmanr(human_rk,v)[0])
colors_b=['#FF9800','#2196F3','#4CAF50','#9C27B0','#795548','#F44336']
bars=ax.barh(fkeys, frhos, color=colors_b, alpha=0.8)
ax.set_xlim(0.95, 1.01)
ax.set_xlabel("Spearman ρ vs human ratings")
ax.set_title("Formula comparison", fontsize=10)
for bar,rho in zip(bars,frhos):
    ax.text(rho+0.0005, bar.get_y()+bar.get_height()/2,
            f'{rho:.4f}', va='center', fontsize=8)
ax.grid(True, axis='x', alpha=0.25)

plt.savefig("/Users/davidderoure/gradus/formula_search.png", dpi=150, bbox_inches='tight')
plt.close()
print("\nPlot saved to formula_search.png")
