"""
Perceptual account of f(p,q) = p + Ω*(q) as a two-stage model.

Stage 1 — Bass establishes harmonic context: Ω*(q)
  The lower note is the q-th partial of an implied fundamental.
  The brain must "extrapolate" down to that fundamental.
  Ω*(q) measures the prime-complexity of that descent.

Stage 2 — Upper note fits into context: p
  Once the fundamental is implied, the upper note must be recognised
  as its p-th partial.  Higher partials are quieter, more narrowly
  tuned, and harder to confirm.  p is the direct cost.

Total perceptual effort: p + Ω*(q)
"""

from math import gcd, log2
from sympy import factorint
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.gridspec import GridSpec

def omega_star(n):
    if n <= 1: return 0
    return sum(e*(p-1) for p, e in factorint(n).items())

def gradus(p, q):
    g = gcd(p,q); p,q = p//g, q//g
    return 1 + omega_star(p) + omega_star(q)

def f(p, q):
    g = gcd(p,q); p,q = p//g, q//g
    return p + omega_star(q)

INTERVALS = [
    (1,  1,  "Unison"),
    (2,  1,  "Octave"),
    (3,  2,  "Fifth"),
    (4,  3,  "Fourth"),
    (5,  4,  "Maj 3rd"),
    (6,  5,  "Min 3rd"),
    (5,  3,  "Maj 6th"),
    (8,  5,  "Min 6th"),
    (9,  8,  "Maj 2nd"),
    (9,  5,  "Min 7th"),
    (15, 8,  "Maj 7th"),
    (16,15,  "Min 2nd"),
    (45,32,  "Tritone"),
]

def reduce(p, q):
    g = gcd(p, q); return p//g, q//g

# ── Figure 1: The two-stage model diagram ────────────────────────────────────

fig = plt.figure(figsize=(16, 13))
fig.patch.set_facecolor('#fafafa')
gs = GridSpec(3, 3, figure=fig, hspace=0.42, wspace=0.35,
              top=0.93, bottom=0.06, left=0.07, right=0.97)

fig.suptitle("A two-stage perceptual account of dissonance", fontsize=14,
             fontweight='bold', y=0.97)

# ── Panel A: conceptual diagram of the two stages ───────────────────────────

ax = fig.add_subplot(gs[0, :])
ax.set_xlim(0, 10); ax.set_ylim(0, 3); ax.axis('off')
ax.set_title("The two-stage model  —  illustrated for the Major third (5/4)",
             fontsize=11, fontweight='bold', loc='left')

# Harmonic series ladder on left (for Fifth 3/2 and Maj 3rd 5/4)
p_ex, q_ex, name_ex = 5, 4, "Major third"
pp, qq = reduce(p_ex, q_ex)

# Draw the harmonic series as rungs
ladder_x = 1.5
rung_ys = {k: 0.35 + (k-1)*0.38 for k in range(1, 9)}
ax.text(ladder_x, 2.75, "Harmonic series\n(implied fundamental at bottom)",
        ha='center', va='top', fontsize=9, color='#444')
for k in range(1, 7):
    y = rung_ys[k]
    col = '#dddddd'
    lw = 0.8
    if k == qq:
        col = '#2196F3'; lw = 2.5
    elif k == pp:
        col = '#E91E63'; lw = 2.5
    ax.plot([ladder_x-0.35, ladder_x+0.35], [y, y], color=col, lw=lw, solid_capstyle='round')
    ax.text(ladder_x-0.55, y, f"P{k}", ha='right', va='center', fontsize=8,
            color=col if col != '#dddddd' else '#999')

# Fundamental
ax.plot([ladder_x-0.35, ladder_x+0.35], [0.1, 0.1], color='#4CAF50', lw=2.5)
ax.text(ladder_x-0.55, 0.1, "F", ha='right', va='center', fontsize=9,
        color='#4CAF50', fontweight='bold')

# Stage 1 arrow: F → P_q
ax.annotate("", xy=(ladder_x+0.5, rung_ys[qq]),
            xytext=(ladder_x+0.5, 0.1),
            arrowprops=dict(arrowstyle='->', color='#2196F3', lw=2))
ax.text(ladder_x+0.7, (rung_ys[qq]+0.1)/2,
        f"Stage 1\nBass = P{qq}\nΩ*({qq}) = {omega_star(qq)}\n(octave of F)",
        ha='left', va='center', fontsize=9, color='#2196F3',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#E3F2FD', edgecolor='#2196F3', alpha=0.8))

# Stage 2 arrow: P_q → P_p
ax.annotate("", xy=(ladder_x-0.5, rung_ys[pp]),
            xytext=(ladder_x-0.5, rung_ys[qq]),
            arrowprops=dict(arrowstyle='->', color='#E91E63', lw=2))
ax.text(ladder_x-0.7, (rung_ys[pp]+rung_ys[qq])/2,
        f"Stage 2\nUpper = P{pp}\nCost = {pp}",
        ha='right', va='center', fontsize=9, color='#E91E63',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#FCE4EC', edgecolor='#E91E63', alpha=0.8))

# Formula display
ax.text(3.8, 1.8,
        f"f({pp}/{qq})  =  p  +  Ω*(q)",
        fontsize=14, ha='left', va='center', fontweight='bold')
ax.text(3.8, 1.35,
        f"            =  {pp}  +  Ω*({qq})",
        fontsize=13, ha='left', va='center', color='#555')
ax.text(3.8, 0.9,
        f"            =  {pp}  +  {omega_star(qq)}",
        fontsize=13, ha='left', va='center', color='#555')
ax.text(3.8, 0.45,
        f"            =  {f(pp,qq)}",
        fontsize=14, ha='left', va='center', fontweight='bold', color='#333')

# Stage descriptions
ax.text(6.5, 2.6, "Stage 1: Bass establishes context", fontsize=10,
        fontweight='bold', color='#2196F3')
ax.text(6.5, 2.25,
        "The lower note is partial q of an implied\n"
        "fundamental F.  The brain extrapolates down\n"
        "to F by 'dividing out' prime factors of q.\n"
        "Cost = Ω*(q) = Σ eᵢ(pᵢ−1): the total prime\n"
        "complexity of q.  Powers of 2 (octaves) cost\n"
        "1 each; factor of 3 costs 2; factor of 5 costs 4.",
        fontsize=9, va='top', color='#333',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#E3F2FD', alpha=0.5))

ax.text(6.5, 1.1, "Stage 2: Upper note reaches into series", fontsize=10,
        fontweight='bold', color='#E91E63')
ax.text(6.5, 0.75,
        "The upper note must be recognised as partial\n"
        "p of F.  Higher partials are quieter (~1/p),\n"
        "harder to tune, and less salient.  Cost = p\n"
        "(the partial number itself).",
        fontsize=9, va='top', color='#333',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#FCE4EC', alpha=0.5))


# ── Panel B: Stage 1 — Ω*(q) for each interval ──────────────────────────────

ax2 = fig.add_subplot(gs[1, 0])
ax2.set_title("Stage 1: Ω*(q)\n(bass context complexity)", fontsize=10,
              fontweight='bold')

names_short = [x[2] for x in INTERVALS]
oq_vals = [omega_star(reduce(p,q)[1]) for p,q,_ in INTERVALS]
colors_b = plt.cm.YlOrRd(np.linspace(0.15, 0.85, len(INTERVALS)))

bars = ax2.barh(names_short, oq_vals, color=colors_b, edgecolor='white', lw=0.5)
ax2.set_xlabel("Ω*(q)  — prime complexity of lower note's partial", fontsize=9)
ax2.invert_yaxis()

# Annotate with q value
for i, (p,q,_) in enumerate(INTERVALS):
    pp, qq = reduce(p,q)
    ax2.text(oq_vals[i]+0.08, i, f"q={qq}", va='center', fontsize=8, color='#555')

ax2.set_xlim(0, 8)
ax2.grid(True, axis='x', alpha=0.3)
ax2.text(0.98, 0.02, "Powers of 2 → low cost\nOdd primes → high cost",
         transform=ax2.transAxes, ha='right', va='bottom', fontsize=8, color='#777')


# ── Panel C: Stage 2 — p for each interval ──────────────────────────────────

ax3 = fig.add_subplot(gs[1, 1])
ax3.set_title("Stage 2: p\n(upper note's partial number)", fontsize=10,
              fontweight='bold')

p_vals = [reduce(p,q)[0] for p,q,_ in INTERVALS]
# log scale because tritone has p=45
bars2 = ax3.barh(names_short, p_vals, color=colors_b, edgecolor='white', lw=0.5)
ax3.set_xlabel("p  — partial number of upper note", fontsize=9)
ax3.invert_yaxis()
ax3.set_xscale('log')
ax3.set_xticks([1,2,3,5,10,20,45])
ax3.get_xaxis().set_major_formatter(plt.ScalarFormatter())

for i, pv in enumerate(p_vals):
    ax3.text(pv*1.05, i, str(pv), va='center', fontsize=8, color='#555')

ax3.grid(True, axis='x', alpha=0.3)
ax3.text(0.98, 0.02, "Low p → salient, easy to tune\nHigh p → faint, hard to confirm",
         transform=ax3.transAxes, ha='right', va='bottom', fontsize=8, color='#777')


# ── Panel D: Combined f = p + Ω*(q) vs Human rank ───────────────────────────

ax4 = fig.add_subplot(gs[1, 2])
ax4.set_title("f = p + Ω*(q)  vs  human rank", fontsize=10, fontweight='bold')

human_rks = [1,2,3,4,5,6,7,8,9,10,11,12,13]
f_vals    = [f(p,q) for p,q,_ in INTERVALS]

scatter_colors = plt.cm.RdYlGn_r(np.linspace(0.1, 0.9, 13))
ax4.scatter(human_rks, f_vals, c=scatter_colors, s=60, zorder=3)
for i, (_,_,name) in enumerate(INTERVALS):
    ax4.annotate(name, (human_rks[i], f_vals[i]),
                 textcoords="offset points", xytext=(5, 0), fontsize=7)

from scipy.stats import spearmanr
rho, _ = spearmanr(human_rks, f_vals)
ax4.set_xlabel("Human consonance rank  (1=most consonant)", fontsize=9)
ax4.set_ylabel("f(p,q) = p + Ω*(q)", fontsize=9)
ax4.set_title(f"f = p + Ω*(q)  vs  human rank\n(ρ = {rho:.3f})", fontsize=10,
              fontweight='bold')
ax4.grid(True, alpha=0.25)


# ── Panel E: Harmonic series diagram for all intervals ──────────────────────

ax5 = fig.add_subplot(gs[2, :])
ax5.set_title(
    "Every interval as two partials of a common implied fundamental  "
    "— cost shown as reach (p, red) and depth (Ω*(q), blue)",
    fontsize=10, fontweight='bold', loc='left')
ax5.set_xlim(-0.5, len(INTERVALS)-0.5)
ax5.set_ylim(-1.2, 6.5)
ax5.set_xticks(range(len(INTERVALS)))
ax5.set_xticklabels([x[2] for x in INTERVALS], rotation=30, ha='right', fontsize=9)
ax5.set_ylabel("Partial number (log₂ scale)", fontsize=9)
ax5.set_yticks([0,1,2,3,4,5,6])
ax5.set_yticklabels(['F','P2','P3','P4','P5','P6','P7+'], fontsize=8)
ax5.grid(True, axis='y', alpha=0.2)

for i, (p, q, name) in enumerate(INTERVALS):
    pp, qq = reduce(p, q)
    # y positions in log₂ scale relative to fundamental
    y_fund = 0
    y_q = log2(qq)   # bass note
    y_p = log2(pp)   # upper note

    # Draw vertical stem from fundamental to bass (Stage 1)
    ax5.plot([i, i], [y_fund, y_q], color='#2196F3', lw=2.5, alpha=0.7, solid_capstyle='round')
    # Draw from bass to upper (Stage 2)
    ax5.plot([i, i], [y_q, y_p], color='#E91E63', lw=2.5, alpha=0.7, solid_capstyle='round')

    # Fundamental dot
    ax5.scatter(i, y_fund, color='#4CAF50', s=60, zorder=4)
    # Bass dot
    ax5.scatter(i, y_q, color='#2196F3', s=80, zorder=5)
    # Upper dot
    ax5.scatter(i, y_p, color='#E91E63', s=80, marker='D', zorder=5)

    # Annotate cost
    mid_stage1 = (y_fund + y_q) / 2
    mid_stage2 = (y_q + y_p) / 2
    oq = omega_star(qq)
    if oq > 0:
        ax5.text(i+0.12, mid_stage1, str(oq), fontsize=7, color='#2196F3',
                 va='center', fontweight='bold')
    ax5.text(i-0.18, mid_stage2, str(pp), fontsize=7, color='#E91E63',
             va='center', fontweight='bold')

    # Total
    fval = f(p, q)
    ax5.text(i, y_p + 0.25, str(fval), ha='center', fontsize=8,
             fontweight='bold', color='#333')

# Legend
leg = [mpatches.Patch(color='#4CAF50', label='Implied fundamental'),
       mpatches.Patch(color='#2196F3', label='Bass note (lower): cost = Ω*(q)'),
       mpatches.Patch(color='#E91E63', label='Upper note: cost = p'),
       mpatches.Patch(color='white', label='Number above = f = p + Ω*(q)')]
ax5.legend(handles=leg, loc='upper left', fontsize=8, framealpha=0.85)

plt.savefig("/Users/davidderoure/gradus/perceptual_account.png",
            dpi=150, bbox_inches='tight')
plt.close()
print("Saved to perceptual_account.png")
