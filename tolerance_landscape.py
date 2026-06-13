"""
Consonance landscape across the octave.

For each cent value c in [0, 1200], find the simplest rational approximation
to the frequency ratio 2^(c/1200) within a tolerance window, then evaluate
G(p/q), f(p/q) = p + Omega*(q), and max(p,q).

Shows how the three functions vary as a continuous function of pitch interval,
and how 12-TET and just intervals sit within the consonance landscape.
"""

from math import gcd, log2, ceil
from fractions import Fraction
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

def best_rational(cents, max_denom=64, tol_cents=15):
    """
    Find the rational p/q with q <= max_denom that is closest to
    2^(cents/1200) and within tol_cents of it.
    Returns (p, q, error_cents) or None if nothing is within tolerance.
    """
    ratio = 2 ** (cents / 1200)
    best = None
    best_err = tol_cents + 1
    for q in range(1, max_denom + 1):
        p = round(ratio * q)
        if p <= 0:
            continue
        r = p / q
        err = abs(1200 * log2(r / ratio))
        if err < best_err:
            best_err = err
            best = (p, q, err)
    return best  # may still be outside tolerance if nothing qualifies

# ── Named intervals ───────────────────────────────────────────────────────────

JUST = [
    (0,    "Unison",        1,  1),
    (204,  "Maj 2nd",       9,  8),
    (316,  "Min 3rd",       6,  5),
    (386,  "Maj 3rd",       5,  4),
    (498,  "Fourth",        4,  3),
    (590,  "Tritone",      45, 32),
    (702,  "Fifth",         3,  2),
    (814,  "Min 6th",       8,  5),
    (884,  "Maj 6th",       5,  3),
    (996,  "Min 7th",       9,  5),
    (1088, "Maj 7th",      15,  8),
    (1200, "Octave",        2,  1),
]

# 12-TET semitone positions and names
TET12 = [
    (0,    "C"),
    (100,  "C#"),
    (200,  "D"),
    (300,  "Eb"),
    (400,  "E"),
    (500,  "F"),
    (600,  "F#"),
    (700,  "G"),
    (800,  "Ab"),
    (900,  "A"),
    (1000, "Bb"),
    (1100, "B"),
    (1200, "C'"),
]

# ── Sweep ─────────────────────────────────────────────────────────────────────

CENTS = np.arange(0, 1201, 1)
TOL   = 15    # cents tolerance for finding a rational approximation
MDENOM = 64   # max denominator to consider

G_vals   = []
F_vals   = []
M_vals   = []
found    = []   # (cents, p, q, err) for each point where a match was found

CLIP_G = 20    # clip display at these values to keep scale readable
CLIP_F = 55
CLIP_M = 50

for c in CENTS:
    hit = best_rational(c, max_denom=MDENOM, tol_cents=TOL)
    if hit and hit[2] <= TOL:
        p, q, err = hit
        G_vals.append(min(gradus(p, q), CLIP_G))
        F_vals.append(min(f(p, q),      CLIP_F))
        M_vals.append(min(maxpq(p, q),  CLIP_M))
        found.append((c, p, q, err))
    else:
        G_vals.append(np.nan)
        F_vals.append(np.nan)
        M_vals.append(np.nan)
        found.append(None)

G_vals = np.array(G_vals, dtype=float)
F_vals = np.array(F_vals, dtype=float)
M_vals = np.array(M_vals, dtype=float)

# ── Plot ──────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(3, 1, figsize=(16, 11), sharex=True)
fig.patch.set_facecolor('#fafafa')

colours = {'G': '#1565C0', 'f': '#B71C1C', 'max': '#2E7D32'}
labels  = {
    'G':   r"Euler's Gradus  $G(p/q) = 1 + \Omega^*(p) + \Omega^*(q)$",
    'f':   r"New formula  $f(p/q) = p + \Omega^*(q)$",
    'max': r"Galileo / Farey  $\max(p,q)$",
}
data  = [('G', G_vals, CLIP_G), ('f', F_vals, CLIP_F), ('max', M_vals, CLIP_M)]

for ax, (key, vals, clip) in zip(axes, data):
    ax.set_facecolor('#fafafa')

    # Shade ±TOL cent windows around each just interval
    for cents_j, name, pj, qj in JUST:
        ax.axvspan(cents_j - TOL, cents_j + TOL,
                   color='gold', alpha=0.18, zorder=0)

    # Plot the function as vertical stems from top (dissonance downward)
    # Actually plot as stems upward from 0 — lower = more consonant
    for i, c in enumerate(CENTS):
        if not np.isnan(vals[i]):
            ax.plot([c, c], [0, vals[i]], color=colours[key],
                    alpha=0.6, lw=1.0, zorder=2)

    # Mark 12-TET positions
    for c_tet, name_tet in TET12:
        ax.axvline(c_tet, color='#888', lw=0.6, ls=':', zorder=1)

    # Mark just interval positions with a label
    for cents_j, name_j, pj, qj in JUST:
        v = f(pj, qj) if key == 'f' else (gradus(pj, qj) if key == 'G'
                                           else maxpq(pj, qj))
        ax.scatter(cents_j, v, color=colours[key], s=50, zorder=5)
        ax.text(cents_j, v + clip * 0.03, name_j,
                ha='center', va='bottom', fontsize=6.5,
                color='#333', rotation=45)

    ax.set_ylabel(f"Value (clipped at {clip})", fontsize=9)
    ax.set_ylim(0, clip * 1.18)
    ax.set_yticks(range(0, clip + 1, 5))
    ax.set_title(labels[key], fontsize=10, fontweight='bold', pad=6)
    ax.grid(True, axis='y', alpha=0.2)

axes[-1].set_xlabel("Interval (cents)", fontsize=10)
axes[-1].set_xticks(range(0, 1201, 100))
axes[-1].set_xticklabels(
    [f"{c}\n{name}" for c, name in TET12] + [],
    fontsize=8
)
# Overlay TET names on x-axis
for c_tet, name_tet in TET12:
    axes[-1].text(c_tet, -clip * 0.13, name_tet,
                  ha='center', va='top', fontsize=7.5, color='#555')

# Legend patches
gold_patch  = mpatches.Patch(color='gold',  alpha=0.5,
                              label=f'±{TOL} cent window around just interval')
grey_line   = mlines.Line2D([], [], color='#888', ls=':', lw=0.8,
                             label='12-TET semitone')
fig.legend(handles=[gold_patch, grey_line],
           loc='lower center', ncol=2, fontsize=9,
           framealpha=0.9, edgecolor='#ccc', bbox_to_anchor=(0.5, 0.01))

fig.suptitle(
    f"Consonance landscape across the octave\n"
    f"Rational approx. with denominator ≤ {MDENOM}, tolerance ±{TOL} cents",
    fontsize=12, fontweight='bold', y=0.99
)

plt.tight_layout(rect=[0, 0.05, 1, 0.99])
plt.savefig("tolerance_landscape.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved tolerance_landscape.png")

# ── Print table: just vs 12-TET comparison ────────────────────────────────────

print()
print("=" * 75)
print("12-TET intervals: nearest just ratio and formula values")
print(f"(denominator ≤ {MDENOM}, tolerance ±{TOL} cents)")
print("=" * 75)
print(f"{'TET':>5}  {'Cents':>5}  {'Nearest just':>12}  {'Err(¢)':>7}  "
      f"{'Gradus':>6}  {'f(p/q)':>6}  {'max':>4}  {'Just name'}")
print("-" * 75)

for c_tet, name_tet in TET12:
    hit = best_rational(c_tet, max_denom=MDENOM, tol_cents=TOL)
    if hit:
        p, q, err = hit
        g = gcd(p, q); pr, qr = p // g, q // g
        # find just name
        just_name = ""
        for cj, nj, pj, qj in JUST:
            if pj == pr and qj == qr:
                just_name = nj
                break
        print(f"{name_tet:>5}  {c_tet:>5}  {pr}/{qr:>10}  {err:>7.2f}  "
              f"{gradus(pr,qr):>6}  {f(pr,qr):>6}  {maxpq(pr,qr):>4}  {just_name}")
    else:
        print(f"{name_tet:>5}  {c_tet:>5}  {'(none)':>12}  {'':>7}  "
              f"{'—':>6}  {'—':>6}  {'—':>4}")

print()
print("=" * 75)
print("Variation within the tolerance window: just interval ± 10 cents")
print("=" * 75)
for cents_j, name_j, pj, qj in JUST[1:-1]:   # skip unison and octave
    v_G  = gradus(pj, qj)
    v_f  = f(pj, qj)
    v_m  = maxpq(pj, qj)
    # nearest competitor within ±50 cents with higher complexity
    rivals = []
    for dc in range(1, 51):
        for sign in [-1, 1]:
            c2 = cents_j + sign * dc
            if c2 < 0 or c2 > 1200:
                continue
            hit = best_rational(c2, max_denom=MDENOM, tol_cents=5)
            if hit:
                p2, q2, err2 = hit
                g2 = gcd(p2, q2); p2r, q2r = p2 // g2, q2 // g2
                if (p2r, q2r) != (pj, qj):
                    rivals.append((abs(dc), p2r, q2r,
                                   gradus(p2r, q2r),
                                   f(p2r, q2r),
                                   maxpq(p2r, q2r)))
    if rivals:
        rivals.sort()
        rd, rp, rq, rG, rf_, rm = rivals[0]
        print(f"\n  {name_j:12s} ({cents_j}¢):  just {pj}/{qj} → "
              f"G={v_G}, f={v_f}, max={v_m}")
        print(f"    Nearest rival: {rp}/{rq} at ±{rd}¢ → "
              f"G={rG}, f={rf_}, max={rm}")
        print(f"    f contrast: {v_f} → {rf_}  (+{rf_-v_f})")
        print(f"    G contrast: {v_G} → {rG}  (+{rG-v_G})")
