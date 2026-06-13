"""
Consonance wells: how G, f, and max behave near a just interval.

Scatter plot of formula values vs cents offset for all rationals p/q
with denominator <= 32 within ±40 cents of the just ratio.
Shows the deep trough at the just ratio and the steep walls on either
side, with f rising most sharply.
"""

from math import gcd, log2
from sympy import factorint
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# ── Core functions ────────────────────────────────────────────────────────────

def omega_star(n):
    if n <= 1: return 0
    return sum(e * (p - 1) for p, e in factorint(n).items())

def gradus(p, q):
    g = gcd(p, q); p, q = p // g, q // g
    return 1 + omega_star(p) + omega_star(q)

def f(p, q):
    g = gcd(p, q); p, q = p // g, q // g
    return p + omega_star(q)

def maxpq(p, q):
    g = gcd(p, q); p, q = p // g, q // g
    return max(p, q)

def neighbours(p0, q0, max_denom=32, window=40):
    seen = set()
    result = []
    ratio0 = p0 / q0
    lo = ratio0 * 2 ** (-window / 1200)
    hi = ratio0 * 2 ** ( window / 1200)
    for q in range(1, max_denom + 1):
        for p in range(max(1, int(lo * q)), int(hi * q) + 2):
            if p <= 0:
                continue
            if abs(1200 * log2((p / q) / ratio0)) > window:
                continue
            g = gcd(p, q)
            pr, qr = p // g, q // g
            if (pr, qr) in seen:
                continue
            seen.add((pr, qr))
            c = 1200 * log2((pr / qr) / ratio0)
            result.append((c, pr, qr,
                           gradus(pr, qr),
                           f(pr, qr),
                           maxpq(pr, qr)))
    result.sort(key=lambda x: x[0])
    return result

# ── Figure ────────────────────────────────────────────────────────────────────

SHOWCASES = [
    (3, 2, "Fifth  (3/2, 702¢)"),
    (5, 4, "Major third  (5/4, 386¢)"),
]

FORMULAS = [
    ('G',   'G(p/q)',         '#1565C0', 'o'),
    ('f',   'f(p/q)',         '#B71C1C', 's'),
    ('max', 'max(p,q)',       '#2E7D32', '^'),
]
IDX = {'G': 3, 'f': 4, 'max': 5}

CLIP   = 75
WINDOW = 40
MDENOM = 32

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor('#fafafa')

for ax, (p0, q0, title) in zip(axes, SHOWCASES):
    pts = neighbours(p0, q0, max_denom=MDENOM, window=WINDOW)
    just_pt  = [(c,pr,qr,G,fv,m) for c,pr,qr,G,fv,m in pts
                if pr == p0 and qr == q0]
    other_pts = [(c,pr,qr,G,fv,m) for c,pr,qr,G,fv,m in pts
                 if not (pr == p0 and qr == q0)]

    ax.set_facecolor('#fafafa')
    ax.axvspan(-15, 15, color='gold', alpha=0.18, zorder=0)
    ax.axvline(0, color='#aaa', lw=0.9, ls='--', zorder=1)
    ax.axhline(0, color='#ccc', lw=0.5, zorder=1)

    # ── Neighbouring rationals ────────────────────────────────────────────────
    for key, label, col, marker in FORMULAS:
        idx = IDX[key]
        xs = [r[0] for r in other_pts]
        ys = [min(r[idx], CLIP) for r in other_pts]
        clipped = [r[idx] > CLIP for r in other_pts]

        ax.scatter(xs, ys, color=col, marker=marker,
                   s=38, alpha=0.55, zorder=3,
                   edgecolors='none')
        # Mark clipped points with an upward arrow
        for x, y, clip in zip(xs, ys, clipped):
            if clip:
                ax.annotate('', xy=(x, CLIP), xytext=(x, CLIP - 6),
                            arrowprops=dict(arrowstyle='->', color=col,
                                            lw=1.2))

        # Draw faint vertical lines from 0 to each neighbour value
        for x, y in zip(xs, ys):
            ax.plot([x, x], [0, y], color=col,
                    alpha=0.12, lw=1.0, zorder=2)

    # ── Just ratio — prominent markers ───────────────────────────────────────
    if just_pt:
        c0, pr0, qr0, G0, f0, m0 = just_pt[0]
        vals = [('G', G0, '#1565C0'), ('f', f0, '#B71C1C'), ('max', m0, '#2E7D32')]
        # Horizontal guide lines from just ratio to y-axis
        for key, val, col in vals:
            ax.plot([c0, c0], [0, val], color=col,
                    lw=3.5, alpha=0.85, zorder=5, solid_capstyle='round')
            ax.scatter([c0], [val], color=col,
                       s=160, marker='D', zorder=7, edgecolors='white', lw=1)
            ax.text(c0 + 1.5, val + 0.8, str(val),
                    fontsize=10, fontweight='bold',
                    color=col, va='bottom', zorder=8)

    # ── Axis labels ───────────────────────────────────────────────────────────
    ax.set_xlim(-WINDOW - 2, WINDOW + 2)
    ax.set_ylim(-3, CLIP + 8)
    ax.set_xlabel("Cents offset from just ratio", fontsize=10)
    ax.set_ylabel("Formula value", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    ax.set_xticks(range(-WINDOW, WINDOW + 1, 10))
    ax.grid(True, axis='y', alpha=0.18, color='#bbb')

    # ── Annotate nearest rivals ───────────────────────────────────────────────
    if other_pts:
        nearest = sorted(other_pts, key=lambda r: abs(r[0]))[:3]
        for c, pr, qr, G, fv, m in nearest:
            best_y = max(G, fv, m)
            ax.text(c, min(best_y, CLIP) + 3,
                    f"${pr}/{qr}$",
                    ha='center', fontsize=7.5, color='#444',
                    va='bottom')

    # Info box
    ax.text(0.98, 0.97,
            f"{len(other_pts)} neighbours\n(denom ≤ {MDENOM}, ±{WINDOW}¢)",
            transform=ax.transAxes, ha='right', va='top',
            fontsize=8, color='#555',
            bbox=dict(boxstyle='round,pad=0.3', fc='white',
                      ec='#ccc', alpha=0.8))

# ── Legend ────────────────────────────────────────────────────────────────────
handles = [
    mlines.Line2D([0],[0], color='#1565C0', marker='o', ms=7, ls='none',
                  label=r"Euler's Gradus  $G = 1+\Omega^*(p)+\Omega^*(q)$"),
    mlines.Line2D([0],[0], color='#B71C1C', marker='s', ms=7, ls='none',
                  label=r"New formula  $f = p+\Omega^*(q)$"),
    mlines.Line2D([0],[0], color='#2E7D32', marker='^', ms=7, ls='none',
                  label=r"Galileo / Farey  $\max(p,q)$"),
    mpatches.Patch(color='gold', alpha=0.5,
                   label='±15¢ perceptual tolerance zone'),
    mlines.Line2D([0],[0], color='#1565C0', marker='D', ms=8, ls='none',
                  markeredgecolor='white', markeredgewidth=1,
                  label='Just ratio (diamond)'),
]
fig.legend(handles=handles, loc='lower center', ncol=3,
           fontsize=9, framealpha=0.95, edgecolor='#ccc',
           bbox_to_anchor=(0.5, 0.0))

fig.suptitle(
    "Consonance wells: formula values near two just ratios\n"
    f"Rationals with denominator ≤ {MDENOM} within ±{WINDOW} cents; "
    "arrows indicate values above 75",
    fontsize=11, fontweight='bold'
)
plt.tight_layout(rect=[0, 0.10, 1, 0.96])
plt.savefig("well_plot.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved well_plot.png")

# ── Sensitivity table ─────────────────────────────────────────────────────────
print()
print("=" * 68)
print("Sensitivity: just ratio vs three nearest rivals")
print("=" * 68)
for p0, q0, title in SHOWCASES:
    pts = neighbours(p0, q0, max_denom=MDENOM, window=WINDOW)
    just_pt  = [(c,pr,qr,G,fv,m) for c,pr,qr,G,fv,m in pts
                if pr==p0 and qr==q0][0]
    rivals   = sorted([(c,pr,qr,G,fv,m) for c,pr,qr,G,fv,m in pts
                       if not (pr==p0 and qr==q0)],
                      key=lambda r: abs(r[0]))
    print(f"\n{title}")
    print(f"  Just {p0}/{q0}:  G={just_pt[3]},  f={just_pt[4]},  max={just_pt[5]}")
    print(f"  {'Rival':<10}  {'Offset':>8}  {'G':>5}  {'f':>5}  {'max':>5}")
    for c,pr,qr,G,fv,m in rivals[:5]:
        print(f"  {pr}/{qr:<8}  {c:>+7.1f}¢  {G:>5}  {fv:>5}  {m:>5}")
