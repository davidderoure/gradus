"""
Generate and verify the superparticular companion sequence A397106:
a(n) = n + A275314(n) for n >= 1, where n is the lower partial of the
superparticular interval (n+1)/n.

Equivalently a(n) = (n+1) + Omega*(n), the two-stage dissonance of (n+1)/n.
gcd(n+1, n) = 1 always, so no coprimality check needed.
"""

from sympy import factorint

def omega_star(n):
    if n <= 1: return 0
    return sum(e * (p - 1) for p, e in factorint(n).items())

def a(n):
    """a(n) = n + 1 + Omega*(n) = n + A275314(n), for n >= 1."""
    assert n >= 1
    return n + 1 + omega_star(n)

# ── Generate sequence ────────────────────────────────────────────────────────

terms = [(n, a(n)) for n in range(1, 10001)]
values = [v for _, v in terms]

print("=" * 65)
print("Superparticular sequence a(n) = n + A275314(n), n >= 1")
print("=" * 65)
print()
print("First 20 terms with interval names:")
names = {1:"Octave",2:"Fifth",3:"Fourth",4:"Maj 3rd",5:"Min 3rd",
         6:"Sep min 3rd",7:"Maj 2nd",8:"Min whole",9:"Maj whole",
         10:"Undecimal",11:"Min 2nd (just)"}
for n, v in terms[:15]:
    name = names.get(n, f"{n+1}/{n}")
    print(f"  n={n:3d}: a({n}) = {n}+1 + Omega*({n}) = {n+1} + {omega_star(n)} = {v:3d}  [{name}]")

print()
print("First 80 terms (for OEIS sequence field):")
print(", ".join(map(str, values[:80])))
print()

# ── b-file ───────────────────────────────────────────────────────────────────

with open("/Users/davidderoure/gradus/b397106.txt", "w") as f:
    f.write("# A397106 a(n) = n + A275314(n), n >= 1. Two-stage dissonance of superparticular interval (n+1)/n.\n")
    for n, v in terms:
        f.write(f"{n} {v}\n")

print(f"b-file written: {len(terms)} terms (n=1..10000)")
print()

# ── Spot checks ──────────────────────────────────────────────────────────────

checks = [(1,2),(2,4),(3,6),(4,7),(5,10),(6,10),(7,14),(8,12),(11,22),(15,22)]
print("Spot-checks:")
all_ok = True
for n, expected in checks:
    got = a(n)
    ok = "✓" if got == expected else f"✗ GOT {got}"
    if got != expected: all_ok = False
    print(f"  a({n:2d}) = {n}+1 + Omega*({n}) = {n+1} + {omega_star(n)} = {got:3d}  {ok}")
print(f"All correct: {all_ok}")
print()

# ── Properties ───────────────────────────────────────────────────────────────

print("=" * 65)
print("PROPERTIES")
print("=" * 65)

non_mono = [(n, values[i-1], v) for i,(n,v) in enumerate(terms)
            if i > 0 and v < values[i-1]][:10]
print(f"\nNon-monotone steps (a(n) < a(n-1)), first 10:")
for n, prev, cur in non_mono:
    print(f"  a({n})={cur} < a({n-1})={prev}  [{n+1}/{n} vs {n}/{n-1}]")
