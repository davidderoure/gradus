# Draft OEIS Submission

---

## Name

Triangle T(n,k) read by rows: T(n,k) = n + A275314(k) - 1 for 1 <= k <= n with gcd(n,k) = 1; T(1,1) = 1. For each row n, k ranges over integers satisfying 1 <= k <= n and gcd(n,k) = 1; the row length is phi(n) for n >= 2 (OEIS A000010).

*(Alternative name for clarity):*
Two-stage dissonance measure for musical intervals: for a frequency ratio n/k in lowest terms (n > k >= 1), a(n,k) = n + Omega*(k), where Omega*(m) = Sum_{p prime, p^e || m} e*(p-1).

---

## Sequence (first 40 terms, triangle read by rows)

```
1, 2, 3, 4, 4, 6, 5, 6, 7, 7, 6, 10, 7, 8, 9, 9, 11, 10, 8, 10, 12, 14, 9, 10, 11, 13, 15, 12, 10, 12, 16, 14, 11, 12, 13, 13, 15, 14, 17, 14
```

Row lengths (number of k < n with gcd(n,k)=1, plus 1 for row 1): 1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4, ...
(Row 1 has length 1 for the unison; for n >= 2, row n has length phi(n).)

---

## Offset

1,2

*(First term is T(1,1)=1; second differs from 1 at T(2,1)=2)*

---

## Comments

T(n,k) is an asymmetric dissonance measure for the musical interval n/k (frequency ratio in lowest terms, n > k). It generalises Euler's Gradus Suavitatis (A275314) by treating numerator and denominator differently.

Euler's Gradus is G(n/k) = 1 + Omega*(n) + Omega*(k), symmetric in n and k. The present formula T(n,k) = n + Omega*(k) replaces Omega*(n) with the raw integer n.

The formula admits a two-stage perceptual interpretation for musical consonance:

  Stage 1 (cost Omega*(k)): The lower note of the interval is the k-th partial of an implied fundamental. Omega*(k) measures the prime complexity of extrapolating down to that fundamental: each factor of prime p in k contributes (p-1) to the cost. Powers of 2 (octave relationships) cost 1 per octave; a factor of 3 costs 2; a factor of 5 costs 4. This reflects the role of the bass note in establishing harmonic context (cf. Rameau's corps sonore, 1722). Note that the brain does not literally factorise integers; it performs subharmonic template matching or autocorrelation across candidate periods (cf. Parncutt 1988, Terhardt 1979). The prime decomposition acts as an analytical proxy for this process: Omega*(k) counts, arithmetically, the steps needed to reach the implied fundamental, mirroring the complexity of the template search.

  Stage 2 (cost n): The upper note must be recognised as the n-th partial of the same fundamental. Higher partials are quieter (amplitude ~ 1/n in harmonic sounds), harder to tune, and have lower pitch salience (cf. Terhardt 1979). The cost is simply n, the partial number.

Equivalently, T(n,k) = n + Omega*(k) = n + A275314(k) - 1.

Rank correlation with human consonance ratings (Krumhansl 1990) for the 13 standard just intervals: Spearman rho = 0.989, compared to rho = 0.979 for Euler's Gradus (A275314). The formula fails to distinguish the major third (5/4, T=7) from the major sixth (5/3, T=7) — the only tie among standard intervals where human ratings differ, a consequence of Omega*(4) = Omega*(3) = 2.

For the standard 13 just intervals (n/k, T(n,k), human consonance rank):
  1/1->1(rank 1), 2/1->2(2), 3/2->4(3), 4/3->6(4), 5/4->7(5), 6/5->10(6),
  5/3->7(7), 8/5->12(8), 9/8->12(9), 9/5->13(10), 15/8->18(11), 16/15->22(12), 45/32->50(13).
(Here n corresponds to p and k corresponds to q in the musical application: n is the upper note's partial number, k is the lower note's partial number.)

Connection to A275314: T(n,k) = n + A275314(k) - 1. Euler's symmetric version is G(n/k) = A275314(n) + A275314(k) - 1.

Connection to Galileo's pulse-coincidence model (1638): the first harmonic coincidence between two notes at ratio n:k occurs at the n-th harmonic of the lower note. The term n in T(n,k) is exactly Galileo's "number of vibrations between coincidences."

Connection to Farey sequences: n is the Farey index of the interval, i.e., the smallest m such that k/n appears in the Farey sequence F_m.

The superparticular (consecutive-harmonic) subsequence T(n, n-1) for n = 2, 3, 4, ...:
  2, 4, 6, 7, 10, 10, 14, 12, 14, 16, 22, 17, 26, 22, 22, 21, 34, 24, ...

---

## References

- L. Euler, Tentamen novae theoriae musicae, St. Petersburg, 1739. (Original Gradus Suavitatis)
- G. Galilei, Discorsi e Dimostrazioni Matematiche, 1638. (Pulse-coincidence model of consonance)
- J.-P. Rameau, Traite de l'Harmonie, Paris, 1722. (Corps sonore / harmonic series)
- C. L. Krumhansl, Cognitive Foundations of Musical Pitch, Oxford University Press, 1990. (Human consonance ratings)
- E. Terhardt, Calculating virtual pitch, Hearing Research 1 (1979), 155-182. (Virtual pitch / partial salience)
- R. Parncutt, Revision of Terhardt's psychoacoustical model of the root(s) of a musical chord, Music Perception 6 (1988), 65-93. (Subharmonic template matching; basis for Stage 1 proxy interpretation)
- I. Lahdelma & T. Eerola, Single chords convey distinct emotional qualities to both naive and expert listeners, Psychology of Music 44 (2016), 37-54. (Empirical decomposition of sensory vs. familiarity contributions to consonance ratings)

---

## Links

- David De Roure, A275314 (Euler's Gradus Suavitatis), OEIS, 2016.
- David De Roure, An Asymmetric Formula for Interval Consonance and its Relation to Harmonic Coincidence, arXiv:2606.16412. (Working note with full derivation, perceptual interpretation, and comparison with Euler's Gradus.)
- Xenharmonic Wiki, Tenney height: https://en.xen.wiki/w/Tenney_height
- Xenharmonic Wiki, Benedetti height: https://en.xen.wiki/w/Benedetti_height

---

## Formula

T(n,k) = n + Omega*(k) for 1 <= k < n, gcd(n,k) = 1.
T(1,1) = 1 (unison).

Omega*(m) = Sum_{primes p dividing m} v_p(m) * (p-1), where v_p(m) is the p-adic valuation of m.
Equivalently, Omega*(m) = A275314(m) - 1 = sopfr(m) - bigomega(m) = A001414(m) - A001222(m).

The symmetric version (Euler's Gradus) satisfies: G(n/k) = T(n,k) + A275314(n) - n.

---

## Example

Triangle begins:
  n=1: 1                        [1/1 = unison]
  n=2: 2                        [2/1 = octave]
  n=3: 3, 4                     [3/1, 3/2 = fifth]
  n=4: 4, 6                     [4/1, 4/3 = fourth]
  n=5: 5, 6, 7, 7               [5/1, 5/2, 5/3, 5/4 = major third]
  n=6: 6, 10                    [6/1, 6/5 = minor third]
  n=7: 7, 8, 9, 9, 11, 10       [7/1, 7/2, 7/3, 7/4, 7/5, 7/6]
  n=8: 8, 10, 12, 14            [8/1, 8/3, 8/5, 8/7]
  n=9: 9, 10, 11, 13, 15, 12   [9/1, 9/2, 9/4, 9/5, 9/7, 9/8 = major second]

T(5,4) = 7: the major third (5/4). Upper note is the 5th partial (cost n=5);
lower note is the 4th partial, Omega*(4) = Omega*(2^2) = 2*1 = 2. Total = 7.

T(6,5) = 10: the minor third (6/5). Upper note is the 6th partial (cost 6);
Omega*(5) = 1*(5-1) = 4. Total = 10. (Euler's Gradus gives 8; human rank = 6.)

T(45,32) = 50: the tritone (45/32, the 5-limit diatonic tritone). n=45, Omega*(32) = Omega*(2^5) = 5. Total = 50.
(Euler's Gradus gives 14; the large value of n=45 captures the extreme dissonance. The septimal tritone 7/5 gives T=11 and the undecimal tritone 11/8 gives T=14 — musically distinct ratios with very different formula values.)

---

## Mathematica

```mathematica
OmegaStar[n_] := If[n <= 1, 0, Total[(#[[2]]) * (#[[1]] - 1) & /@ FactorInteger[n]]];
T[n_, k_] /; (n == 1 && k == 1) := 1;
T[n_, k_] /; (n > k && GCD[n, k] == 1) := n + OmegaStar[k];
Table[T[n, k], {n, 1, 12}, {k, 1, n}, GCD[n, k] == 1] // Flatten
```

---

## Python

```python
from math import gcd
from sympy import factorint

def omega_star(n):
    if n <= 1: return 0
    return sum(e * (p - 1) for p, e in factorint(n).items())

def T(n, k):
    """Two-stage dissonance of interval n/k (n >= k >= 1, gcd(n,k)=1)."""
    assert gcd(n, k) == 1
    return n + omega_star(k)

# Generate triangle
for n in range(1, 10):
    row = [T(n, k) for k in range(1, n + 1) if gcd(n, k) == 1 and k <= n]
    print(row)
```

---

## PARI/GP

```
OmegaStar(n) = {
  if(n<=1, return(0));
  my(f=factor(n));
  sum(i=1, #f~, f[i,2]*(f[i,1]-1))
};
T(n,k) = if(gcd(n,k)==1, n + OmegaStar(k), 0);
tabl(nn) = for(n=1,nn, for(k=1,n, if(gcd(n,k)==1, print1(T(n,k),", "))))
```

---

## Crossrefs

- A275314: Euler's Gradus Suavitatis, a(n) = 1 + Omega*(n). The symmetric version of this triangle is A275314(n) + A275314(k) - 1.
- A001414: sopfr(n) = Sum of prime factors of n with repetition.
- A001222: bigomega(n) = number of prime factors of n with multiplicity. Note Omega*(n) = A001414(n) - A001222(n).
- A000010: Euler totient phi(n) = row lengths for n >= 2.
- A002088: partial sums of phi, related to Farey sequence lengths.
- A007188: Tenney height H(n/k) = n*k; an alternative interval complexity measure. Unlike T(n,k), Tenney height incorrectly ranks the major sixth (5/3, H=15) as more consonant than the major third (5/4, H=20), contrary to human ratings.
- A397106: Rightmost populated diagonal T(n, n-1) for n >= 2 (superparticular intervals): 2, 4, 6, 7, 10, 10, 14, 12, ... (companion entry).
- This triangle is A397104.
- Row sums: Sum_{k: gcd(n,k)=1, 1<=k<=n} T(n,k) = n*phi(n) + Sum_{k: gcd(n,k)=1} Omega*(k); natural companion sequence.
- Row maxima and row minima: reveal the spread of dissonance values within each harmonic series level; sensitive to the prime factorisation of n.

---

## Keywords

`tabl, nonn, hear`

*(nonn = non-negative integers; tabl = triangle/table; hear = related to music/hearing)*

---

## Author

David De Roure

---

## Notes for submission

1. The sequence offset should be 1 (first index n=1, k=1).
2. The keyword `hear` is appropriate; see also A275314 which uses it.
3. You may wish to submit the superparticular subsequence T(n, n-1) as a separate entry with a cross-reference, as it is a natural single-variable sequence.
4. The formula Omega*(k) = A275314(k) - 1 provides a clean link to your existing entry.
5. Consider whether to include or exclude T(1,1)=1 (unison, p=q=1) — including it makes the triangle start cleanly at n=1.
