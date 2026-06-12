from math import gcd
from sympy import factorint

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

print(f"{'Interval':<16}  {'Ratio':>7}  {'Gradus':>6}  {'p+Ω*(q)':>8}  {'Human rank':>10}")
print(f"{'-'*16}  {'-'*7}  {'-'*6}  {'-'*8}  {'-'*10}")
for p, q, name, human in INTERVALS:
    g = gcd(p,q); pp, qq = p//g, q//g
    print(f"{name:<16}  {pp}/{qq:<5}  {gradus(p,q):>6}  {f(p,q):>8}  {human:>10}")
