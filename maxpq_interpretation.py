"""
Musical interpretations of max(p,q) as a consonance metric.

For ratio p:q (p > q, coprime), max(p,q) = p achieves rho=0.989 vs human
ratings — slightly better than Gradus.  What does it mean musically?

We explore five perspectives:
  1. First coincidence: how far up the harmonic series until first meeting
  2. Implied fundamental: p is the partial number of the upper note
  3. Farey/Stern-Brocot: the 'depth' of the ratio in the harmonic lattice
  4. Galileo/Mersenne pulse coincidence: period of reinforcement in time
  5. Comparison with Tenney height, odd limit, Benedetti height
"""

from math import gcd, log, log2, ceil
from sympy import factorint, isprime
from fractions import Fraction
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

# ── Core ────────────────────────────────────────────────────────────────────

def omega_star(n):
    if n <= 1: return 0
    return sum(e*(p-1) for p, e in factorint(n).items())

def gradus(p, q):
    g = gcd(p, q); p, q = p//g, q//g
    return 1 + omega_star(p) + omega_star(q)

def odd_part(n):
    while n % 2 == 0: n //= 2
    return n

def odd_limit(p, q):
    """Max odd number appearing in the reduced ratio."""
    g = gcd(p, q); p, q = p//g, q//g
    return max(odd_part(p), odd_part(q))

def tenney_height(p, q):
    """Tenney harmonic distance: log2(p*q)."""
    g = gcd(p, q); p, q = p//g, q//g
    return log2(p * q)

def benedetti_height(p, q):
    """Benedetti's 'beating number': p*q."""
    g = gcd(p, q); p, q = p//g, q//g
    return p * q

def stern_brocot_depth(p, q):
    """
    Depth of p/q in the Stern-Brocot tree (number of steps from root 1/1).
    This equals the sum of partial quotients in the continued fraction of p/q.
    """
    g = gcd(p, q); p, q = p//g, q//g
    # Continued fraction of p/q
    depth = 0
    a, b = p, q
    while b:
        depth += a // b
        a, b = b, a % b
    return depth - 1   # subtract 1 so root 1/1 has depth 0

def first_coincidence_harmonic(p, q):
    """Harmonic number of lower note at which first coincidence occurs."""
    g = gcd(p, q); p, q = p//g, q//g
    return p   # = max(p,q) since p > q

def coincidence_period_ms(p, q, f0_hz=110.0):
    """Time between harmonic coincidences in milliseconds."""
    g = gcd(p, q); p, q = p//g, q//g
    # First coincidence at harmonic p of lower note, frequency p*f0
    # Period of that shared harmonic = 1/(p*f0)
    # Coincidences recur every p*(1/f0) seconds = p/f0 s
    return 1000.0 * p / f0_hz   # ms between coincidences

# ── Human reference ──────────────────────────────────────────────────────────

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

intervals = [(x[0],x[1]) for x in INTERVALS]
names     = [x[2] for x in INTERVALS]
human_rk  = np.array([x[3] for x in INTERVALS])

# ── Compute all metrics ──────────────────────────────────────────────────────

def r(p,q): g=gcd(p,q); return p//g, q//g

metrics = {}
for p,q,name,hr in INTERVALS:
    pp,qq = r(p,q)
    metrics[(p,q)] = {
        'Gradus':         gradus(p,q),
        'max(p,q)':       max(pp,qq),
        'min(p,q)':       min(pp,qq),
        'sum p+q':        pp+qq,
        'p*q (Benedetti)':pp*qq,
        'log2(p*q) (Tenney)': log2(pp*qq),
        'odd limit':      odd_limit(p,q),
        'SB depth':       stern_brocot_depth(p,q),
        'coincide harmonic': first_coincidence_harmonic(p,q),
        'coincide period ms': coincidence_period_ms(p,q),
    }

# ── Part 1: Metric comparison table ─────────────────────────────────────────

print("=" * 80)
print("METRIC VALUES FOR ALL 13 INTERVALS")
print("=" * 80)
mkeys = ['Gradus','max(p,q)','odd limit','sum p+q','p*q (Benedetti)','log2(p*q) (Tenney)','SB depth']
print(f"  {'Interval':<10} {'Human':>5} " +
      "".join(f"  {k[:10]:>10}" for k in mkeys))
print(f"  {'-'*10} {'-'*5} " + "  ----------" * len(mkeys))
for p,q,name,hr in INTERVALS:
    vals = metrics[(p,q)]
    row = f"  {name:<10} {hr:>5} "
    for k in mkeys:
        v = vals[k]
        row += f"  {v:>10.2f}" if isinstance(v,float) else f"  {v:>10}"
    print(row)

# ── Part 2: Correlation with human ratings ───────────────────────────────────

print()
print("=" * 80)
print("RANK CORRELATIONS WITH HUMAN RATINGS")
print("=" * 80)
print(f"  {'Metric':<25}  {'Spearman rho':>12}")
print(f"  {'-'*25}  {'-'*12}")
for k in mkeys + ['coincide harmonic','coincide period ms']:
    vals = np.array([metrics[(p,q)][k] for p,q,_,_ in INTERVALS])
    rho, _ = spearmanr(human_rk, vals)
    print(f"  {k:<25}  {rho:>+12.4f}")

# ── Part 3: The first-coincidence / Galileo interpretation ──────────────────

print()
print("=" * 80)
print("INTERPRETATION 1: GALILEO/MERSENNE PULSE COINCIDENCE (1638)")
print("=" * 80)
print("""
Galileo Galilei (Dialogues, 1638) described consonance in terms of the
coincidence of vibrations.  Two strings vibrating at ratio p:q (coprime)
produce coinciding pulses every p vibrations of the lower string.

  "The more frequent the coincidences, the sweeter the consonance."
  — Galileo

For ratio p:q, coincidences occur at times: p/f, 2p/f, 3p/f, ...
(where f is the lower frequency).  max(p,q) = p is the number of
vibrations of the lower string between successive coincidences.

This is IDENTICAL to the first-coincidence harmonic metric.
Galileo's pulse model is equivalent to max(p,q).
""")

print(f"  {'Interval':<10} {'p:q':>7}  {'Pulses between coincidences':>28}  "
      f"{'Period at A2=110Hz':>20}")
for p,q,name,hr in INTERVALS:
    pp,qq = r(p,q)
    ms = coincidence_period_ms(p,q, 110)
    bar = '█' * min(pp, 20) + ('…' if pp>20 else '')
    print(f"  {name:<10} {pp}:{qq:<5}  {bar:<28}  {ms:>16.1f} ms")

# ── Part 4: Harmonic series / implied fundamental ────────────────────────────

print()
print("=" * 80)
print("INTERPRETATION 2: POSITION IN HARMONIC SERIES")
print("=" * 80)
print("""
For ratio p:q (coprime), the two notes are the p-th and q-th partials
of a common implied fundamental at 1/p below the lower note.

  Lower note  = partial q  of the fundamental
  Upper note  = partial p  of the fundamental

max(p,q) = p = the partial number of the UPPER note.
The higher a partial, the less strongly the ear hears it as part of
a single fused tone — higher partials are progressively more 'foreign.'

This matches Rameau's theory of the corps sonore (1722): consonance
arises from tones that appear low in a single harmonic series.
""")

print(f"  {'Interval':<10} {'p:q':>7}  lower=partial  upper=partial  "
      f"implied fund. (re lower, cents)")
for p,q,name,hr in INTERVALS:
    pp,qq = r(p,q)
    # Implied fundamental is qq steps below lower note (lower note = qq-th partial)
    # Fundamental frequency = f_lower / qq
    # Fundamental in cents below lower note:
    fund_cents = -1200 * log2(qq)  # negative = below lower note
    print(f"  {name:<10} {pp}:{qq:<5}  lower=P{qq:<6}  upper=P{pp:<6}  "
          f"{fund_cents:>+10.1f}¢ (implied fund)")

# ── Part 5: Stern-Brocot / Farey depth ──────────────────────────────────────

print()
print("=" * 80)
print("INTERPRETATION 3: STERN-BROCOT TREE DEPTH")
print("=" * 80)
print("""
The Stern-Brocot tree arranges all positive rationals in a binary tree.
The depth of p/q in this tree equals the sum of continued-fraction
partial quotients — a measure of arithmetic 'complexity.'

The Farey sequence F_n contains all fractions a/b with b ≤ n (in [0,1]).
The fraction q/p (our interval ratio inverted to be < 1) first appears
in F_p.  So p = max(p,q) is the 'Farey index' of the interval.

This gives max(p,q) a precise lattice-theoretic meaning:
  p is the smallest n such that the interval appears in the Farey sequence F_n.
""")

print(f"  {'Interval':<10} {'p/q':>7}  {'max(p,q)':>8}  {'SB depth':>9}  "
      f"{'Farey index':>11}  {'human rk':>8}")
for p,q,name,hr in INTERVALS:
    pp,qq = r(p,q)
    sb = stern_brocot_depth(p,q)
    print(f"  {name:<10} {pp}/{qq:<5}  {pp:>8}  {sb:>9}  {pp:>11}  {hr:>8}")

print()
sb_vals = np.array([metrics[(p,q)]['SB depth'] for p,q,_,_ in INTERVALS])
max_vals = np.array([metrics[(p,q)]['max(p,q)'] for p,q,_,_ in INTERVALS])
rho_sb, _  = spearmanr(human_rk, sb_vals)
rho_mx, _  = spearmanr(human_rk, max_vals)
print(f"  SB depth correlation with human: {rho_sb:+.4f}")
print(f"  max(p,q) correlation with human: {rho_mx:+.4f}")
print(f"  (They differ for 9/8 vs 9/5: same max but different SB depths)")
# Show the CF expansions
for p,q,name,hr in INTERVALS:
    pp,qq = r(p,q)
    # CF of pp/qq
    cf = []
    a,b = pp,qq
    while b:
        cf.append(a//b)
        a,b = b,a%b
    sb = sum(cf)-1
    if name in ["Maj 2nd","Min 7th"]:
        print(f"    {name} ({pp}/{qq}): CF={cf}, SB depth={sb}")

# ── Part 6: Where max(p,q) still fails ──────────────────────────────────────

print()
print("=" * 80)
print("REMAINING FAILURES OF max(p,q)")
print("=" * 80)

from scipy.stats import rankdata
max_vals = np.array([max(*r(p,q)) for p,q,_,_ in INTERVALS])
max_ranks = rankdata(max_vals)

print(f"\n  {'Interval':<10} {'Human':>6}  {'max rk':>7}  {'Diff':>5}  note")
for i,(p,q,name,hr) in enumerate(INTERVALS):
    pp,qq = r(p,q)
    mr = max_ranks[i]
    diff = mr - hr
    note = ""
    if abs(diff) >= 1.5:
        note = "← still tied" if max_vals[i] == max_vals[i-1] or \
               (i<len(INTERVALS)-1 and max_vals[i]==max_vals[i+1]) else "← residual error"
    print(f"  {name:<10} {hr:>6}  {mr:>7.1f}  {diff:>+5.1f}  {note}")

print("""
  Remaining issue: max(p,q)=5 for BOTH Maj 3rd (5/4) and Maj 6th (5/3).
  The metric cannot distinguish them because the upper note is the 5th
  partial in both cases — same distance from the fundamental.

  The ONLY difference between 5/4 and 5/3 is the lower note:
    5/4: lower = 4th partial  (2 octaves above fundamental)
    5/3: lower = 3rd partial  (octave + fifth above fundamental)

  Humans prefer 5/4. Why might the 4th partial feel more stable than the 3rd?
  - 4 = 2² is a pure octave-stack, fully subsumed by octave equivalence
  - 3 is the first 'new' prime — the fifth introduces a new dimension
  Alternatively: 5/4 is a narrower interval (386¢ vs 884¢), and
  narrower intervals within the same partial class may feel more 'direct.'
""")

# ── Part 7: A refined metric — max(odd_part(p), odd_part(q)) ────────────────

print("=" * 80)
print("REFINED METRIC: max(odd_part(p), odd_part(q))  — 'odd numerator'")
print("=" * 80)
print("""
Since powers of 2 represent octave equivalences (trivially consonant),
stripping factors of 2 and taking the max odd part may better reflect
'genuine' harmonic complexity:

  odd_part(4) = 1  (just octaves)
  odd_part(3) = 3  (a fifth above the octave-equivalent)

  5/4: max(odd_part(5), odd_part(4)) = max(5, 1) = 5
  5/3: max(odd_part(5), odd_part(3)) = max(5, 3) = 5  ← still tied!

  But: odd_part distinguishes 8/5 from 6/5:
  8/5: max(odd_part(8), odd_part(5)) = max(1, 5) = 5
  6/5: max(odd_part(6), odd_part(5)) = max(3, 5) = 5  ← also tied

  This IS the odd-limit metric from just intonation theory (Harry Partch).
  It groups intervals by the highest odd number in their ratio.
""")

odd_lim_vals = np.array([odd_limit(p,q) for p,q,_,_ in INTERVALS])
rho_ol, _ = spearmanr(human_rk, odd_lim_vals)
print(f"  Odd limit correlation with human: {rho_ol:+.4f}")
print()
print(f"  {'Interval':<10} {'Human':>6}  {'max(p,q)':>9}  {'odd limit':>9}  {'Gradus':>7}")
for p,q,name,hr in INTERVALS:
    pp,qq = r(p,q)
    print(f"  {name:<10} {hr:>6}  {pp:>9}  {odd_limit(p,q):>9}  {gradus(p,q):>7}")

print(f"""
  Odd limit groups: 1-limit: {{Unison,Octave}}
                    3-limit: {{Fifth, Fourth}}
                    5-limit: {{Maj3, Min3, Maj6, Min6}}
                    9-limit: {{Maj2, Min7}}
                   15-limit: {{Maj7, Min2}}
                   45-limit: {{Tritone}}

  Within the 5-limit group, humans rank: Maj3 > Min3 > Maj6 > Min6
  Within the 9-limit group: Maj2 > Min7
  Odd limit (like max(p,q)) cannot resolve these within-group orderings.
  Gradus resolves some (Min3 vs Min6, Maj2 vs Min7) but not all (Maj3 vs Maj6).
""")


# ── Plot ─────────────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(16, 11))
fig.suptitle("Musical interpretations of max(p,q) as consonance metric", fontsize=13)
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.35)

# 1. Galileo pulse diagram — show coincidences for selected intervals
ax = fig.add_subplot(gs[0, :2])
ax.set_title("Galileo pulse coincidence (first 45 vibrations of lower note)", fontsize=10)
selected = [("Fifth 3:2",3,2,"#2196F3"),("Maj 3rd 5:4",5,4,"#4CAF50"),
            ("Maj 6th 5:3",5,3,"#FF9800"),("Min 3rd 6:5",6,5,"#9C27B0"),
            ("Maj 2nd 9:8",9,8,"#F44336")]
ypos = {name: i for i, (name,_,_,_) in enumerate(selected)}
ax.set_xlim(0, 46); ax.set_ylim(-0.7, len(selected)-0.3)
ax.set_xlabel("Vibration number of lower note")
ax.set_yticks(range(len(selected)))
ax.set_yticklabels([name for name,_,_,_ in selected], fontsize=9)
ax.grid(True, axis='x', alpha=0.3)
for name, p, q, col in selected:
    y = ypos[name]
    # Lower note pulses: at every integer
    for n in range(1, 46):
        ax.plot(n, y-0.15, '|', color='gray', alpha=0.3, markersize=8)
    # Upper note pulses: at every p/q
    n2 = 1
    while n2*p/q <= 45:
        ax.plot(n2*p/q, y+0.15, '|', color=col, alpha=0.6, markersize=8)
        n2 += 1
    # Mark coincidences
    for k in range(1, 46//p + 1):
        ax.plot(k*p, y, 'o', color=col, markersize=9, zorder=5)
    ax.axvline(p, color=col, alpha=0.4, linestyle='--', linewidth=0.8)
ax.text(0.98, 0.02, "● = coincidence  | = pulse\nFirst coincidence marked with dashed line",
        transform=ax.transAxes, ha='right', va='bottom', fontsize=8, color='gray')

# 2. Metric correlations bar chart
ax = fig.add_subplot(gs[0, 2])
metric_names = ['max(p,q)\n=Galileo','odd limit\n(Partch)','Gradus\n(Euler)',
                'Tenney\nlog₂(pq)','SB depth','Benedetti\np·q']
metric_vals_list = [
    [max(*r(p,q)) for p,q,_,_ in INTERVALS],
    [odd_limit(p,q) for p,q,_,_ in INTERVALS],
    [gradus(p,q) for p,q,_,_ in INTERVALS],
    [tenney_height(p,q) for p,q,_,_ in INTERVALS],
    [stern_brocot_depth(p,q) for p,q,_,_ in INTERVALS],
    [benedetti_height(p,q) for p,q,_,_ in INTERVALS],
]
rhos = [spearmanr(human_rk, v)[0] for v in metric_vals_list]
colors_bar = ['#2196F3','#FF9800','#4CAF50','#9C27B0','#795548','#F44336']
bars = ax.barh(metric_names, rhos, color=colors_bar, alpha=0.8)
ax.set_xlim(0.5, 1.0)
ax.set_xlabel("Spearman ρ vs human ratings")
ax.set_title("Consonance metric comparison", fontsize=10)
ax.axvline(rhos[2], color='#4CAF50', linestyle='--', alpha=0.5)  # Gradus line
for bar, rho in zip(bars, rhos):
    ax.text(rho+0.003, bar.get_y()+bar.get_height()/2,
            f'{rho:.3f}', va='center', fontsize=8)
ax.grid(True, axis='x', alpha=0.25)

# 3. Harmonic series diagram
ax = fig.add_subplot(gs[1, :2])
ax.set_title("Harmonic series position: every interval as two partials of a common fundamental",
             fontsize=10)
fundamentals = {}
max_partial = 10  # only show up to partial 10 for clarity
selected2 = [("Octave",2,1),("Fifth",3,2),("Fourth",4,3),
             ("Maj 3rd",5,4),("Maj 6th",5,3),("Min 3rd",6,5)]

# Draw harmonic series
partials = range(1, max_partial+1)
ax.set_xlim(0.5, max_partial+0.5)
ax.set_ylim(-0.5, len(selected2)-0.5)
ax.set_yticks(range(len(selected2)))
ax.set_yticklabels([f"{n} ({p}/{q})" for n,p,q in selected2], fontsize=9)
ax.set_xlabel("Partial number")
ax.set_xticks(partials)
ax.grid(True, axis='x', alpha=0.2)

colors2 = ['#2196F3','#4CAF50','#FF9800','#9C27B0','#F44336','#795548']
for i, (name, p, q) in enumerate(selected2):
    pp, qq = r(p,q)
    # Draw all partials lightly
    for n in partials:
        ax.plot(n, i, 'o', color='lightgray', markersize=12, zorder=1)
        ax.text(n, i, str(n), ha='center', va='center', fontsize=7, color='gray')
    # Highlight the two notes
    ax.plot(qq, i, 'o', color=colors2[i], markersize=18, zorder=3, alpha=0.8)
    ax.text(qq, i, str(qq), ha='center', va='center', fontsize=8,
            color='white', fontweight='bold')
    ax.plot(pp, i, 's', color=colors2[i], markersize=18, zorder=3, alpha=0.8)
    ax.text(pp, i, str(pp), ha='center', va='center', fontsize=8,
            color='white', fontweight='bold')
    ax.annotate("", xy=(pp,i), xytext=(qq,i),
                arrowprops=dict(arrowstyle='->', color=colors2[i], lw=1.5, alpha=0.5))

circle = mpatches.Patch(color='gray', label='● lower note (partial q)')
square = mpatches.Patch(color='gray', label='■ upper note (partial p) = max(p,q)')
ax.legend(handles=[circle,square], loc='lower right', fontsize=8)

# 4. Summary table
ax = fig.add_subplot(gs[1, 2])
ax.axis('off')
rows = []
for p,q,name,hr in INTERVALS:
    pp,qq = r(p,q)
    g = gradus(p,q)
    mx = max(pp,qq)
    ol = odd_limit(p,q)
    rows.append([name, f"{pp}:{qq}", str(hr), str(mx), str(ol), str(g)])
table = ax.table(cellText=rows,
                 colLabels=['Interval','ratio','human','max','odd lim','Gradus'],
                 loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(8)
table.scale(1, 1.3)
ax.set_title("Metric values", fontsize=10, pad=15)

plt.savefig("/Users/davidderoure/gradus/maxpq_interpretation.png", dpi=150,
            bbox_inches='tight')
plt.close()
print("Plot saved to maxpq_interpretation.png")

# ── Final summary ────────────────────────────────────────────────────────────

print()
print("=" * 80)
print("SUMMARY: THE FOUR FACES OF max(p,q)")
print("=" * 80)
print(f"""
  max(p,q) = p  (for ratio p:q in lowest terms, p > q) is simultaneously:

  1. GALILEO (1638):  the number of vibrations of the lower string between
     successive pulse coincidences.  "More frequent coincidences = sweeter."

  2. HARMONIC SERIES (Rameau 1722, Riemann 1874):  the partial number of
     the upper note within the implied harmonic series.  Notes low in the
     series fuse naturally; high partials sound foreign.

  3. FAREY SEQUENCE (Cauchy 1816):  the smallest n such that the interval
     q/p appears in the Farey sequence F_n.  An index of arithmetic depth.

  4. FIRST COINCIDENCE HARMONIC:  the position in the harmonic series of
     the lower note where the two notes first share a frequency.

  All four are mathematically identical.  max(p,q) is the most natural
  single-number summary of how 'deep' an interval sits in the harmonic lattice.

  HOW IT DIFFERS FROM GRADUS:
  Gradus = 1 + Ω*(p) + Ω*(q) weights prime factors: large primes penalised
  more than their size would suggest (e.g. 3 costs 2, but 4=2² also costs 2).
  max(p,q) = p treats all primes equally — only the magnitude matters.

  REMAINING FAILURE (both metrics):
  5/4 (Maj 3rd) and 5/3 (Maj 6th) both have max=5 and Gradus=7.
  The only difference is whether the lower note is the 4th (=2²) or 3rd partial.
  Humans prefer 5/4.  Resolution may require:
    — octave-equivalence weighting (4 collapses to 1 under octave reduction)
    — a two-dimensional metric: (max upper partial, max lower partial)
    — cultural/statistical familiarity with thirds over sixths in tonal music
""")

# Final rho summary
print(f"  Metric correlations vs human ratings:")
for name, vals in [("max(p,q)",  [max(*r(p,q)) for p,q,_,_ in INTERVALS]),
                   ("odd limit", [odd_limit(p,q) for p,q,_,_ in INTERVALS]),
                   ("Gradus",    [gradus(p,q)    for p,q,_,_ in INTERVALS]),
                   ("Tenney",    [tenney_height(p,q) for p,q,_,_ in INTERVALS]),
                   ("SB depth",  [stern_brocot_depth(p,q) for p,q,_,_ in INTERVALS])]:
    rho, _ = spearmanr(human_rk, vals)
    print(f"    {name:<20}  ρ = {rho:+.4f}")
