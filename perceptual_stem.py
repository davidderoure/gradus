"""
Stem diagram: every interval as two partials of a common implied fundamental.
Blue stem = Stage 1 (bass establishes context, cost = Ω*(q)).
Red stem  = Stage 2 (upper note reaches into series, cost = p).
"""

from math import gcd, log2
from sympy import factorint
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import spearmanr

def omega_star(n):
    if n <= 1: return 0
    return sum(e*(pr-1) for pr, e in factorint(n).items())

def f(p, q):
    g = gcd(p, q); p, q = p//g, q//g
    return p + omega_star(q)

def reduce(p, q):
    g = gcd(p, q); return p//g, q//g

INTERVALS = [
    (1,  1,  "Unison",        1),
    (2,  1,  "Octave",        2),
    (3,  2,  "Fifth",         3),
    (4,  3,  "Fourth",        4),
    (5,  4,  "Major third",   5),
    (6,  5,  "Minor third",   6),
    (5,  3,  "Major sixth",   7),
    (8,  5,  "Minor sixth",   8),
    (9,  8,  "Major second",  9),
    (9,  5,  "Minor seventh",10),
    (15, 8,  "Major seventh",11),
    (16,15,  "Minor second", 12),
    (45,32,  "Tritone",      13),
]

fig, ax = plt.subplots(figsize=(15, 7))
fig.patch.set_facecolor('#fafafa')

ax.set_xlim(-0.7, len(INTERVALS) - 0.3)
ax.set_ylim(-0.7, 6.2)
ax.set_xticks(range(len(INTERVALS)))
ax.set_xticklabels(
    [f"{name}\n{p}/{q}\nf={f(p,q)}" for p,q,name,_ in INTERVALS],
    fontsize=9.5
)
ax.set_ylabel("Partial number of implied fundamental  (log₂ scale)", fontsize=10)
ax.set_yticks([0, 1, log2(3), 2, log2(5), log2(6), log2(8), log2(9),
               log2(15), log2(16), log2(32), log2(45)])
ax.set_yticklabels(['F (1)', 'P2', 'P3', 'P4', 'P5', 'P6', 'P8',
                    'P9', 'P15', 'P16', 'P32', 'P45'], fontsize=8)
ax.grid(True, axis='y', alpha=0.18, color='#aaa')
ax.set_facecolor('#fafafa')

# Colour gradient: consonant green → dissonant red
cmap = plt.cm.RdYlGn_r
n = len(INTERVALS)

for i, (p, q, name, human) in enumerate(INTERVALS):
    pp, qq = reduce(p, q)
    col = cmap((human - 1) / 12)

    y_fund = 0
    y_q    = log2(qq)
    y_p    = log2(pp)

    # Stage 1: fundamental → bass  (blue tint)
    ax.plot([i, i], [y_fund, y_q],
            color='#1565C0', lw=4, alpha=0.55, solid_capstyle='round', zorder=2)

    # Stage 2: bass → upper  (red tint)
    ax.plot([i, i], [y_q, y_p],
            color='#B71C1C', lw=4, alpha=0.55, solid_capstyle='round', zorder=2)

    # Nodes
    ax.scatter(i, y_fund, color='#2E7D32', s=90,  zorder=5)   # fundamental
    ax.scatter(i, y_q,    color='#1565C0', s=110, zorder=5)   # bass
    ax.scatter(i, y_p,    color='#B71C1C', s=110,
               marker='D', zorder=5)                           # upper

    oq = omega_star(qq)
    # Stage 1 cost label
    if oq > 0:
        mid1 = (y_fund + y_q) / 2
        ax.text(i + 0.14, mid1, f"Ω*={oq}", fontsize=8,
                color='#1565C0', va='center', fontweight='bold')

    # Stage 2 cost label
    mid2 = (y_q + y_p) / 2 if y_p != y_q else y_p + 0.15
    ax.text(i - 0.16, mid2, f"p={pp}", fontsize=8,
            color='#B71C1C', va='center', ha='right', fontweight='bold')

    # Human rank in a small box at top
    ax.text(i, y_p + 0.28, f"h={human}", ha='center', fontsize=8,
            color='#333',
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                      edgecolor='#ccc', alpha=0.8))

# Legend
leg = [
    mpatches.Patch(color='#2E7D32', label='Implied fundamental F'),
    mpatches.Patch(color='#1565C0',
                   label='Stage 1 (blue): bass = partial q of F,  cost = Ω*(q)'),
    mpatches.Patch(color='#B71C1C',
                   label='Stage 2 (red): upper = partial p of F,  cost = p'),
    mpatches.Patch(color='white', edgecolor='#ccc',
                   label='h = human consonance rank  (1 = most consonant)'),
]
ax.legend(handles=leg, loc='upper left', fontsize=9.5, framealpha=0.9,
          edgecolor='#ccc')

ax.set_title(
    "Dissonance as two-stage perceptual effort:  f(p/q) = p + Ω*(q)\n"
    "Every interval shown as two partials of a common implied fundamental",
    fontsize=12, fontweight='bold', pad=12)

plt.tight_layout()
plt.savefig("/Users/davidderoure/gradus/perceptual_stem.png",
            dpi=150, bbox_inches='tight')
plt.close()
print("Saved to perceptual_stem.png")
