from math import gcd
from sympy import factorint
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

def omega_star(n):
    if n <= 1: return 0
    return sum(e*(p-1) for p, e in factorint(n).items())

def gradus(p, q):
    g = gcd(p, q); p, q = p//g, q//g
    return 1 + omega_star(p) + omega_star(q)

def f(p, q):
    g = gcd(p, q); p, q = p//g, q//g
    return p + omega_star(q)

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

rows = []
for p, q, name, human in INTERVALS:
    g = gcd(p,q); pp, qq = p//g, q//g
    rows.append([name, f"{pp}/{qq}", gradus(p,q), f(p,q), human])

col_labels = ["Interval", "Ratio", "Gradus", "p + Ω*(q)", "Human rank"]
col_widths  = [0.22,       0.10,    0.14,      0.14,        0.14]

fig, ax = plt.subplots(figsize=(9, 5.2))
ax.axis('off')

# Colour scale: consonant (green) → dissonant (red)
max_human = 13
def row_colour(human_rank, alpha=0.18):
    t = (human_rank - 1) / (max_human - 1)   # 0 = consonant, 1 = dissonant
    r = t
    g = 1 - t
    return (r, g*0.7, 0, alpha)

# Header
header_y = 1.0
header_x = 0.02
for j, (label, w) in enumerate(zip(col_labels, col_widths)):
    x = header_x + sum(col_widths[:j])
    ax.text(x + w/2, header_y, label,
            ha='center', va='center', fontsize=11, fontweight='bold',
            transform=ax.transAxes)

# Divider line under header
ax.plot([0.02, 0.98], [header_y - 0.04, header_y - 0.04],
        color='#333333', linewidth=1.2, transform=ax.transAxes)

row_height = 0.072
for i, (row, (p,q,name,human)) in enumerate(zip(rows, INTERVALS)):
    y = header_y - 0.085 - i * row_height
    bg = row_colour(human)
    # Background rectangle
    rect = plt.Rectangle((0.02, y - row_height*0.48), 0.96, row_height*0.92,
                          transform=ax.transAxes, color=bg, zorder=0)
    ax.add_patch(rect)
    for j, (val, w) in enumerate(zip(row, col_widths)):
        x = header_x + sum(col_widths[:j])
        align = 'left' if j == 0 else 'center'
        xpos  = x + 0.01 if j == 0 else x + w/2
        ax.text(xpos, y, str(val),
                ha=align, va='center', fontsize=10,
                transform=ax.transAxes)

# Bottom divider
y_last = header_y - 0.085 - (len(rows)-1)*row_height
ax.plot([0.02, 0.98], [y_last - row_height*0.5, y_last - row_height*0.5],
        color='#aaaaaa', linewidth=0.8, transform=ax.transAxes)

# Colour legend
ax.text(0.02, 0.03,
        "Colour: green = consonant → red = dissonant (human ranking)",
        fontsize=8, color='#555555', transform=ax.transAxes)

ax.set_title("Dissonance measures for small whole-number ratios",
             fontsize=13, fontweight='bold', pad=12)

plt.tight_layout()
plt.savefig("/Users/davidderoure/gradus/table_plot.png", dpi=150, bbox_inches='tight')
plt.show()
print("Saved to table_plot.png")
