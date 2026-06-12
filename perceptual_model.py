"""
Perceptual consonance model based on Plomp-Levelt / Helmholtz roughness.

The ear hears dissonance not just from exact harmonic coincidences but from
beating between near-coinciding harmonics.  Two partials at frequencies f1, f2
beat at rate |f1-f2| Hz.  Roughness is maximal when the beat rate is ~25% of
the critical bandwidth CB(f) — the frequency resolution of the basilar membrane.

Model (after Sethares 2005, Plomp & Levelt 1965):
  For each pair of harmonics (n1 from note 1, n2 from note 2):
    df   = |f1 - f2|       frequency difference in Hz
    f_cb = min(f1, f2)     reference for CB
    CB(f) ≈ 1.72 * f^0.65  critical bandwidth (Hz)
    x    = df / (0.25*CB)   normalised beating rate
    r(x) = x * exp(1-x)    roughness shape  [peaks at x=1]

  Roughness contribution = a(n1) * a(n2) * r(x)
  where a(n) = 1/n  (harmonic amplitude rolloff)

Total roughness = sum over all harmonic pairs within hearing range.
Consonance = negative roughness (less rough = more consonant).

We compare:
  1. Euler's Gradus
  2. Exact flat coincidence count
  3. Roughness model (Helmholtz/Plomp-Levelt)
  4. Human consonance ratings (literature)
"""

from math import gcd, log, exp
import numpy as np
from scipy.stats import spearmanr, kendalltau
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sympy import factorint

# ── Core ────────────────────────────────────────────────────────────────────

def omega_star(n):
    if n <= 1: return 0
    return sum(e*(p-1) for p, e in factorint(n).items())

def gradus(p, q):
    g = gcd(p, q); p, q = p//g, q//g
    return 1 + omega_star(p) + omega_star(q)

def critical_bandwidth(f_hz):
    """Bark-scale CB approximation (Zwicker & Terhardt)."""
    return 1.72 * f_hz**0.65

def roughness_shape(x):
    """Plomp-Levelt roughness curve: peaks at x=1 (df = 0.25*CB)."""
    if x <= 0: return 0.0
    return x * exp(1.0 - x)

def roughness(p, q, f0_hz=110.0, N=32, amp_rolloff=1.0):
    """
    Total roughness between two notes at ratio p:q.
    f0_hz: fundamental of note 1 in Hz
    N: number of harmonics per note
    amp_rolloff: amplitude of harmonic n = 1/n^amp_rolloff
    """
    g = gcd(p, q); p, q = p//g, q//g
    freq_ratio = p / q

    total = 0.0
    for n1 in range(1, N+1):
        f1 = n1 * f0_hz
        a1 = n1 ** (-amp_rolloff)
        for n2 in range(1, N+1):
            f2 = n2 * f0_hz * freq_ratio
            a2 = n2 ** (-amp_rolloff)
            df = abs(f1 - f2)
            if df == 0:
                continue   # exact coincidence: no roughness
            f_cb = min(f1, f2)
            cb = critical_bandwidth(f_cb)
            x = df / (0.25 * cb)
            total += a1 * a2 * roughness_shape(x)
    return total


# ── Human consonance reference rankings ─────────────────────────────────────
# From Krumhansl (1990) and Huron (1994) listener studies.
# Lower rank = more consonant.

HUMAN_REFERENCE = {
    # (p, q) : (name, human_rank)   [1=most consonant]
    (1, 1):  ("Unison",       1),
    (2, 1):  ("Octave",       2),
    (3, 2):  ("Fifth",        3),
    (4, 3):  ("Fourth",       4),
    (5, 4):  ("Maj 3rd",      5),
    (6, 5):  ("Min 3rd",      6),
    (5, 3):  ("Maj 6th",      7),
    (8, 5):  ("Min 6th",      8),
    (9, 8):  ("Maj 2nd",      9),
    (9, 5):  ("Min 7th",     10),
    (15, 8): ("Maj 7th",     11),
    (16,15): ("Min 2nd",     12),
    (45,32): ("Tritone",     13),
}

intervals = list(HUMAN_REFERENCE.keys())
names     = [HUMAN_REFERENCE[k][0]  for k in intervals]
human_rk  = [HUMAN_REFERENCE[k][1]  for k in intervals]

# ── Compute all metrics ──────────────────────────────────────────────────────

N_HARM  = 32
F0      = 110.0   # Hz, A2

gradus_vals    = [gradus(p, q)                                    for p, q in intervals]
rough_vals     = [roughness(p, q, f0_hz=F0, N=N_HARM)            for p, q in intervals]
exact_vals     = [N_HARM // (p // gcd(p,q))                      for p, q in intervals]

# Sensitivity analysis: different f0 values
rough_55  = [roughness(p, q, f0_hz=55,  N=N_HARM) for p, q in intervals]
rough_220 = [roughness(p, q, f0_hz=220, N=N_HARM) for p, q in intervals]
rough_440 = [roughness(p, q, f0_hz=440, N=N_HARM) for p, q in intervals]

# Sensitivity: amplitude rolloff exponent
rough_flat = [roughness(p, q, f0_hz=F0, N=N_HARM, amp_rolloff=0) for p, q in intervals]
rough_1    = [roughness(p, q, f0_hz=F0, N=N_HARM, amp_rolloff=1) for p, q in intervals]
rough_2    = [roughness(p, q, f0_hz=F0, N=N_HARM, amp_rolloff=2) for p, q in intervals]

# ── Rank comparison ──────────────────────────────────────────────────────────

h = np.array(human_rk)
g = np.array(gradus_vals)
r = np.array(rough_vals)
e = np.array(exact_vals)

print("=" * 65)
print("RANK CORRELATIONS WITH HUMAN CONSONANCE RATINGS")
print(f"(Spearman rho; n={len(intervals)} intervals)")
print("=" * 65)

metrics = [
    ("Euler Gradus",                     g,                True),
    ("Exact coincidence count",          e,                False),
    ("Roughness f0=55Hz  (low bass)",    np.array(rough_55),  True),
    ("Roughness f0=110Hz (A2)",          r,                True),
    ("Roughness f0=220Hz (middle)",      np.array(rough_220), True),
    ("Roughness f0=440Hz (high)",        np.array(rough_440), True),
    ("Roughness amp~flat",               np.array(rough_flat),True),
    ("Roughness amp~1/n",                np.array(rough_1),   True),
    ("Roughness amp~1/n²",               np.array(rough_2),   True),
]

for label, arr, higher_is_worse in metrics:
    sign = 1 if higher_is_worse else -1
    rho, pval = spearmanr(h, sign * arr)
    print(f"  {label:<35}  rho = {rho:+.4f}  (p={pval:.3f})")

# ── Detailed ranking table ───────────────────────────────────────────────────

print()
print("=" * 65)
print("RANKING TABLE (rank 1 = most consonant)")
print("=" * 65)
print(f"  {'Interval':<10} {'Human':>6}  {'Gradus':>6}  {'Exact':>6}  {'Rough':>6}")
print(f"  {'-'*10} {'-'*6}  {'-'*6}  {'-'*6}  {'-'*6}")

# Convert to ranks
from scipy.stats import rankdata
g_rank = rankdata(g)
r_rank = rankdata(r)
e_rank = rankdata(-np.array(exact_vals))

for i, name in enumerate(names):
    print(f"  {name:<10} {human_rk[i]:>6}  {g_rank[i]:>6.0f}  "
          f"{e_rank[i]:>6.0f}  {r_rank[i]:>6.0f}")

# ── Where roughness and Gradus disagree ──────────────────────────────────────

print()
print("=" * 65)
print("DISAGREEMENTS: roughness rank vs Gradus rank")
print("=" * 65)
print(f"  {'Interval':<10} {'Gradus rk':>9}  {'Rough rk':>8}  {'Human rk':>8}  {'Diff':>5}")
diffs = []
for i, name in enumerate(names):
    diff = r_rank[i] - g_rank[i]
    diffs.append((abs(diff), name, g_rank[i], r_rank[i], human_rk[i], diff))

for _, name, gr, rr, hr, diff in sorted(diffs, reverse=True):
    print(f"  {name:<10} {gr:>9.0f}  {rr:>8.0f}  {hr:>8}  {diff:>+5.0f}")

# ── Effect of fundamental frequency ─────────────────────────────────────────

print()
print("=" * 65)
print("ROUGHNESS VS HUMAN: effect of fundamental frequency")
print("(rho with human ratings for different registers)")
print("=" * 65)
for f0, rv in [(55, rough_55), (110, rough_1), (220, rough_220), (440, rough_440)]:
    rho, _ = spearmanr(h, rv)
    print(f"  f0 = {f0:4d} Hz: rho = {rho:+.4f}")

# ── Plot ─────────────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(15, 10))
fig.suptitle("Consonance models vs human ratings\n"
             f"(f0={F0}Hz, {N_HARM} harmonics, amp~1/n)", fontsize=12)
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

# 1. Roughness vs human rank
ax = fig.add_subplot(gs[0, 0])
ax.scatter(human_rk, rough_vals, color='steelblue')
for i, name in enumerate(names):
    ax.annotate(name, (human_rk[i], rough_vals[i]),
                textcoords="offset points", xytext=(3,2), fontsize=7)
rho, _ = spearmanr(human_rk, rough_vals)
ax.set_xlabel("Human rank (1=consonant)")
ax.set_ylabel("Roughness (higher=more dissonant)")
ax.set_title(f"Roughness model\nρ={rho:+.3f}")
ax.grid(True, alpha=0.25)

# 2. Gradus vs human rank
ax = fig.add_subplot(gs[0, 1])
jitter = np.random.default_rng(7).uniform(-0.1, 0.1, len(names))
ax.scatter(human_rk, gradus_vals, color='darkorange')
for i, name in enumerate(names):
    ax.annotate(name, (human_rk[i], gradus_vals[i]),
                textcoords="offset points", xytext=(3,2), fontsize=7)
rho, _ = spearmanr(human_rk, gradus_vals)
ax.set_xlabel("Human rank (1=consonant)")
ax.set_ylabel("Gradus")
ax.set_title(f"Euler's Gradus\nρ={rho:+.3f}")
ax.grid(True, alpha=0.25)

# 3. Roughness vs Gradus
ax = fig.add_subplot(gs[0, 2])
ax.scatter(gradus_vals, rough_vals, color='forestgreen')
for i, name in enumerate(names):
    ax.annotate(name, (gradus_vals[i], rough_vals[i]),
                textcoords="offset points", xytext=(3,2), fontsize=7)
rho, _ = spearmanr(gradus_vals, rough_vals)
ax.set_xlabel("Gradus")
ax.set_ylabel("Roughness")
ax.set_title(f"Roughness vs Gradus\nρ={rho:+.3f}")
ax.grid(True, alpha=0.25)

# 4. Roughness at different registers
ax = fig.add_subplot(gs[1, 0])
f0s = [55, 110, 220, 440]
rhos_f0 = []
for f0 in f0s:
    rv = [roughness(p, q, f0_hz=f0, N=N_HARM) for p, q in intervals]
    rho, _ = spearmanr(human_rk, rv)
    rhos_f0.append(rho)
ax.plot(f0s, rhos_f0, 'o-', color='steelblue')
ax.axhline(y=spearmanr(human_rk, gradus_vals)[0], color='darkorange',
           linestyle='--', label="Gradus")
ax.set_xscale('log')
ax.set_xlabel("Fundamental frequency (Hz)")
ax.set_ylabel("Spearman ρ vs human")
ax.set_title("Roughness ρ vs register")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.25)

# 5. Effect of amplitude rolloff
ax = fig.add_subplot(gs[1, 1])
alphas = [0, 0.5, 1.0, 1.5, 2.0, 3.0]
rhos_amp = []
for alpha in alphas:
    rv = [roughness(p, q, f0_hz=F0, N=N_HARM, amp_rolloff=alpha) for p, q in intervals]
    rho, _ = spearmanr(human_rk, rv)
    rhos_amp.append(rho)
ax.plot(alphas, rhos_amp, 'o-', color='forestgreen')
ax.axhline(y=spearmanr(human_rk, gradus_vals)[0], color='darkorange',
           linestyle='--', label="Gradus")
ax.set_xlabel("Amplitude rolloff exponent α  (amp = 1/n^α)")
ax.set_ylabel("Spearman ρ vs human")
ax.set_title("Roughness ρ vs harmonic rolloff")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.25)

# 6. N harmonics effect
ax = fig.add_subplot(gs[1, 2])
Ns = [4, 8, 16, 32, 64]
rhos_N = []
for Nv in Ns:
    rv = [roughness(p, q, f0_hz=F0, N=Nv) for p, q in intervals]
    rho, _ = spearmanr(human_rk, rv)
    rhos_N.append(rho)
ax.plot(Ns, rhos_N, 'o-', color='purple')
ax.axhline(y=spearmanr(human_rk, gradus_vals)[0], color='darkorange',
           linestyle='--', label="Gradus")
ax.set_xlabel("Number of harmonics N")
ax.set_ylabel("Spearman ρ vs human")
ax.set_title("Roughness ρ vs harmonic count")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.25)

plt.savefig("/Users/davidderoure/gradus/perceptual_model.png", dpi=150, bbox_inches='tight')
plt.close()
print("\nPlot saved to perceptual_model.png")

# ── Summary ──────────────────────────────────────────────────────────────────

print()
print("=" * 65)
print("SUMMARY")
print("=" * 65)

g_rho, _  = spearmanr(human_rk, gradus_vals)
r_rho, _  = spearmanr(human_rk, rough_vals)
e_rho, _  = spearmanr(human_rk, [-x for x in exact_vals])
gr_rho, _ = spearmanr(gradus_vals, rough_vals)

print(f"""
  Gradus vs human ratings:      ρ = {g_rho:+.3f}
  Roughness vs human ratings:   ρ = {r_rho:+.3f}
  Exact count vs human ratings: ρ = {e_rho:+.3f}
  Roughness vs Gradus:          ρ = {gr_rho:+.3f}

  Key cases where roughness and Gradus diverge:
  - Tritone (45/32): very high Gradus (14), but roughness may not be
    extreme because its harmonics are widely spaced (few beats)
  - Maj 2nd (9/8) and Min 3rd (6/5): similar Gradus, but roughness
    captures the 'grinding' of the second vs the relative smoothness
    of the third
  - Register matters: higher f0 raises roughness of all intervals
    because CB grows slower than frequency, so more partials conflict
""")
