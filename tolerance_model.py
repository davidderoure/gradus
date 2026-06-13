"""
Partial-beating tolerance model.

For a just interval p/q, a deviation of Δ cents produces a beat rate at
the first coincident partial of approximately:

    beat_rate ≈ p · f₀ · ln(2)/1200 · |Δ|   Hz

Setting a perceptual threshold τ Hz gives a tolerance half-width:

    Δ_tol(p/q, f₀) ≈ 1730 · τ / (p · f₀)   cents

Tolerance ∝ 1/p — the same p that appears in f(p/q) = p + Ω*(q).

Three figures:
  1. Tolerance wells by interval, showing ET deviations vs tolerance bands
     for three threshold values (4, 8, 15 Hz).
  2. Tolerance as a function of register (f₀) for each interval, with ET
     deviation marked — shows where each interval crosses into roughness.
  3. Summary: tolerance half-width vs p, annotated with intervals.
"""

from math import gcd, log, log2
from sympy import factorint
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# ── Constants ─────────────────────────────────────────────────────────────────

LN2_OVER_1200 = log(2) / 1200   # converts cents deviation to fractional ratio

def delta_tol(p, f0, tau):
    """Tolerance half-width in cents for partial-beating threshold tau Hz."""
    return 1730 * tau / (p * f0)

def beat_rate(p, f0, delta_cents):
    """Beat rate in Hz at first coincident partial."""
    return p * f0 * LN2_OVER_1200 * abs(delta_cents)

# ── Intervals ─────────────────────────────────────────────────────────────────

INTERVALS = [
    # (name, p, q, ET_deviation_cents)
    ("Unison",       1,  1,   0.0),
    ("Octave",       2,  1,   0.0),
    ("Fifth",        3,  2,   1.96),   # 700¢ vs 702¢
    ("Fourth",       4,  3,   1.96),   # 500¢ vs 498¢
    ("Maj 3rd",      5,  4,  13.69),   # 400¢ vs 386¢
    ("Min 3rd",      6,  5,  15.64),   # 300¢ vs 316¢
    ("Maj 6th",      5,  3,  15.64),   # 900¢ vs 884¢  (same p as maj 3rd)
    ("Min 6th",      8,  5,  13.69),   # 800¢ vs 814¢
    ("Maj 2nd",      9,  8,   3.91),   # 200¢ vs 204¢
    ("Min 7th",      9,  5,   3.91),   # 1000¢ vs 996¢
    ("Maj 7th",     15,  8,  11.73),   # 1100¢ vs 1088¢
    ("Min 2nd",     16, 15,  11.73),   # 100¢ vs 112¢
]

# Remove unison and octave from plots (trivial)
INTERVALS_PLOT = INTERVALS[2:]

THRESHOLDS = [4, 8, 15]          # Hz — tight, moderate, generous
TAU_COLS   = ['#2E7D32', '#E65100', '#7B1FA2']  # green, orange, purple
F0_REF     = 220.0               # A3 — reference register

# ── Figure 1: Tolerance wells per interval ────────────────────────────────────

fig1, ax1 = plt.subplots(figsize=(13, 6))
fig1.patch.set_facecolor('#fafafa')
ax1.set_facecolor('#fafafa')

n = len(INTERVALS_PLOT)
xs = np.arange(n)
width = 0.22

for ti, (tau, col) in enumerate(zip(THRESHOLDS, TAU_COLS)):
    tols = [delta_tol(iv[1], F0_REF, tau) for iv in INTERVALS_PLOT]
    offset = (ti - 1) * width
    bars = ax1.bar(xs + offset, tols, width,
                   color=col, alpha=0.55, label=f'τ = {tau} Hz tolerance',
                   zorder=3)

# ET deviations as horizontal marks
for i, iv in enumerate(INTERVALS_PLOT):
    name, p, q, et_dev = iv
    # Draw ET deviation as a thick horizontal line
    ax1.plot([i - 0.38, i + 0.38], [et_dev, et_dev],
             color='#B71C1C', lw=2.5, zorder=5, solid_capstyle='round')
    ax1.scatter(i, et_dev, color='#B71C1C', s=60, zorder=6)

# Shade "safe" (ET within most generous tolerance) vs "stretched"
for i, iv in enumerate(INTERVALS_PLOT):
    name, p, q, et_dev = iv
    tol_generous = delta_tol(p, F0_REF, THRESHOLDS[-1])
    if et_dev > tol_generous:
        ax1.axvspan(i - 0.5, i + 0.5, color='#FFCDD2', alpha=0.35, zorder=0)
    elif et_dev > delta_tol(p, F0_REF, THRESHOLDS[1]):
        ax1.axvspan(i - 0.5, i + 0.5, color='#FFF9C4', alpha=0.45, zorder=0)

ax1.set_xticks(xs)
ax1.set_xticklabels(
    [f"{iv[0]}\n{iv[1]}/{iv[2]}\nET: {iv[3]:.1f}¢" for iv in INTERVALS_PLOT],
    fontsize=8.5
)
ax1.set_ylabel("Tolerance half-width  (cents)", fontsize=10)
ax1.set_ylim(0, 42)
ax1.set_title(
    f"Partial-beating tolerance  Δ_tol = 1730·τ / (p·f₀),  reference f₀ = {int(F0_REF)} Hz (A3)\n"
    "Coloured bars = tolerance for three thresholds   ·   Red line = 12-TET deviation from just",
    fontsize=10, fontweight='bold', pad=8
)
ax1.grid(True, axis='y', alpha=0.2)

# Legend
handles = [mpatches.Patch(color=c, alpha=0.7, label=f'τ = {t} Hz')
           for t, c in zip(THRESHOLDS, TAU_COLS)]
handles += [
    mlines.Line2D([0],[0], color='#B71C1C', lw=2.5,
                  label='12-TET deviation from just'),
    mpatches.Patch(color='#FFCDD2', alpha=0.6,
                   label='ET outside all tolerances'),
    mpatches.Patch(color='#FFF9C4', alpha=0.7,
                   label='ET outside strict tolerance only'),
]
ax1.legend(handles=handles, fontsize=8.5, loc='upper right',
           framealpha=0.92, edgecolor='#ccc')

plt.tight_layout()
plt.savefig("tolerance_model_1.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved tolerance_model_1.png")

# ── Figure 2: Tolerance vs register ──────────────────────────────────────────

# For each interval, the tolerance in cents = 1730*tau/(p*f0)
# and the ET deviation is fixed. The crossing point f0* = 1730*tau/(p*delta_ET)
# is the register above which the ET deviation exceeds tolerance.

F0_RANGE = np.linspace(55, 880, 500)   # A1 to A5

fig2, axes2 = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
fig2.patch.set_facecolor('#fafafa')
axes2 = axes2.flatten()

SHOW = [("Fifth",   3, 2, 1.96),
        ("Maj 3rd", 5, 4, 13.69),
        ("Min 3rd", 6, 5, 15.64),
        ("Maj 2nd", 9, 8, 3.91)]

octave_ticks = [55, 110, 220, 440, 880]
octave_names = ["A1", "A2", "A3\n(220 Hz)", "A4\n(440 Hz)", "A5"]

for ax, (name, p, q, et_dev) in zip(axes2, SHOW):
    ax.set_facecolor('#fafafa')

    for tau, col in zip(THRESHOLDS, TAU_COLS):
        tols = delta_tol(p, F0_RANGE, tau)
        ax.plot(F0_RANGE, tols, color=col, lw=2,
                label=f'τ = {tau} Hz', alpha=0.85)
        # Crossing point
        f_cross = 1730 * tau / (p * et_dev)
        if 55 < f_cross < 880:
            ax.axvline(f_cross, color=col, lw=1, ls=':', alpha=0.7)
            ax.text(f_cross + 8, delta_tol(p, f_cross, tau) + 0.4,
                    f'{f_cross:.0f} Hz', fontsize=7, color=col)

    # ET deviation as horizontal line
    ax.axhline(et_dev, color='#B71C1C', lw=2.2, ls='-',
               label=f'ET deviation ({et_dev:.1f}¢)', zorder=5)

    # Shade where ET exceeds generous tolerance
    tol_gen = delta_tol(p, F0_RANGE, THRESHOLDS[-1])
    ax.fill_between(F0_RANGE, et_dev, tol_gen,
                    where=(tol_gen < et_dev),
                    color='#FFCDD2', alpha=0.4, zorder=0)

    ax.set_xlim(55, 880)
    ax.set_ylim(0, min(45, et_dev * 4))
    ax.set_xscale('log')
    ax.set_xticks(octave_ticks)
    ax.set_xticklabels(octave_names, fontsize=8)
    ax.set_ylabel("Tolerance half-width (¢)", fontsize=9)
    ax.set_title(f"{name}  ({p}/{q}),  ET deviation = {et_dev:.1f}¢",
                 fontsize=10, fontweight='bold')
    ax.grid(True, alpha=0.18)
    ax.legend(fontsize=8, loc='upper right', framealpha=0.9, edgecolor='#ccc')

fig2.suptitle(
    "Tolerance vs register: at what pitch does the 12-TET deviation exceed partial-beating threshold?\n"
    "Red line = ET deviation (fixed).  Dotted verticals = crossing frequency for each threshold.",
    fontsize=10, fontweight='bold'
)
plt.tight_layout()
plt.savefig("tolerance_model_2.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved tolerance_model_2.png")

# ── Figure 3: Tolerance vs p ──────────────────────────────────────────────────

fig3, ax3 = plt.subplots(figsize=(10, 6))
fig3.patch.set_facecolor('#fafafa')
ax3.set_facecolor('#fafafa')

p_range = np.linspace(1, 18, 300)
for tau, col in zip(THRESHOLDS, TAU_COLS):
    tols = delta_tol(p_range, F0_REF, tau)
    ax3.plot(p_range, tols, color=col, lw=2.2,
             label=f'τ = {tau} Hz', alpha=0.9)

# Mark each interval
for name, p, q, et_dev in INTERVALS_PLOT:
    for tau, col in zip(THRESHOLDS, TAU_COLS):
        tol = delta_tol(p, F0_REF, tau)
    # Use middle threshold for label placement
    tol_mid = delta_tol(p, F0_REF, THRESHOLDS[1])
    # ET deviation dot
    colour = '#B71C1C' if et_dev > delta_tol(p, F0_REF, THRESHOLDS[-1]) \
             else ('#E65100' if et_dev > delta_tol(p, F0_REF, THRESHOLDS[1])
                  else '#2E7D32')
    ax3.scatter(p, et_dev, color=colour, s=80, zorder=6,
                marker='D' if colour == '#B71C1C' else 'o')
    ax3.text(p + 0.15, et_dev + 0.3, name,
             fontsize=8, color='#333', va='bottom')

ax3.set_xlabel("Numerator  p  of just ratio  p/q", fontsize=10)
ax3.set_ylabel("Cents", fontsize=10)
ax3.set_xlim(1.5, 18)
ax3.set_ylim(0, 34)
ax3.set_xticks(range(2, 17))
ax3.set_title(
    f"Tolerance Δ_tol = 1730·τ/(p·f₀) vs partial number p,  f₀ = {int(F0_REF)} Hz (A3)\n"
    "Coloured curves = tolerance.  Points = 12-TET deviation "
    "(green ✓ = within all, orange = marginal, red ✗ = outside all)",
    fontsize=10, fontweight='bold', pad=8
)
ax3.grid(True, alpha=0.18)

handles = [mlines.Line2D([0],[0], color=c, lw=2.2, label=f'τ = {t} Hz')
           for t, c in zip(THRESHOLDS, TAU_COLS)]
handles += [
    mlines.Line2D([0],[0], color='#2E7D32', marker='o', ms=8, ls='none',
                  label='ET within all tolerances'),
    mlines.Line2D([0],[0], color='#E65100', marker='o', ms=8, ls='none',
                  label='ET exceeds strict tolerance'),
    mlines.Line2D([0],[0], color='#B71C1C', marker='D', ms=8, ls='none',
                  label='ET outside all tolerances'),
]
ax3.legend(handles=handles, fontsize=9, loc='upper right',
           framealpha=0.92, edgecolor='#ccc')

plt.tight_layout()
plt.savefig("tolerance_model_3.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved tolerance_model_3.png")

# ── Print summary table ───────────────────────────────────────────────────────

print()
print("=" * 78)
print(f"Partial-beating tolerance at f₀ = {int(F0_REF)} Hz (A3)")
print(f"{'Interval':<12} {'p':>3} {'ET dev':>7} "
      + "".join(f"  tol@{t}Hz" for t in THRESHOLDS)
      + "  Status")
print("-" * 78)
for name, p, q, et_dev in INTERVALS_PLOT:
    tols = [delta_tol(p, F0_REF, t) for t in THRESHOLDS]
    if et_dev <= tols[0]:
        status = "✓ within all"
    elif et_dev <= tols[1]:
        status = "~ within 8 Hz"
    elif et_dev <= tols[2]:
        status = "~ within 15 Hz"
    else:
        status = "✗ outside all"
    row = f"{name:<12} {p:>3} {et_dev:>6.1f}¢"
    row += "".join(f"  {t:>8.1f}¢" for t in tols)
    row += f"  {status}"
    print(row)

print()
print("=" * 78)
print("Register at which ET deviation crosses each threshold (Hz)")
print(f"{'Interval':<12} {'p':>3} {'ET dev':>7} "
      + "".join(f"  cross@{t}Hz" for t in THRESHOLDS))
print("-" * 78)
for name, p, q, et_dev in INTERVALS_PLOT:
    if et_dev == 0:
        continue
    crossings = [1730 * tau / (p * et_dev) for tau in THRESHOLDS]
    row = f"{name:<12} {p:>3} {et_dev:>6.1f}¢"
    for f_cross in crossings:
        if f_cross < 55:
            row += "  < A1 (always)"
        elif f_cross > 3520:
            row += "  > C8 (never) "
        else:
            row += f"  {f_cross:>8.0f} Hz"
    print(row)
