# Verification of the Werner projective-locality certificate

**Verdict:** the claimed lower bound survives independent reconstruction. The uploaded certificate supports

\[
v_c^{\mathrm{proj}}>\frac{6914213000}{10002162779}
=0.6912717931882401\ldots>0.6912.
\]

This is a lower bound on the critical visibility, not a determination of its exact value. No claim that the actual threshold equals 0.6912 is established.

Audit date: 24 October 2026 (India). All unpublished inputs and calculations were kept private. Only the already-public universal Grothendieck bound was checked against public literature; no candidate coefficients, manuscript text, numerical research claim, or distinctive construction was submitted to a search service.

## What was independently checked

| Component | Result |
|---|---|
| Exact candidate | 21 strictly positive rational weights sum exactly to one; all 42 shapes have the specified 55 degree-nine monomials. |
| Exact scalar response | Reconstructed the degree-1101 source polynomial, its exact rational Chebyshev coefficients, the omitted tail, and every stored degree-111 response coefficient. All identities pass. |
| Probability range | Recomputed the full 1025-by-2048 shape grid and its Bernstein extension. All 42 shapes are globally below 0.9, so every response probability is valid. |
| Uploaded trigonometric table | Independently regenerated with 110-digit Decimal/Taylor arithmetic. Every uploaded sine/cosine enclosure checked contains the corresponding exact value. |
| Polar nodes and CC weights | Independently regenerated, without treating the uploaded table as an input to this calculation. |
| Associated Legendre functions | Recomputed using 110-digit arithmetic and propagated error balls, independently of the original double-precision recurrence envelope. |
| Harmonic projections | Recomputed from direct degree-nine polynomial evaluation and direct 2048-term Fourier sums. The independent calculation does not use the original 32-point shape reconstruction or FFT. |
| All 76 archived harmonic intervals | Each contains the entire corresponding independently computed conservative interval. |
| All 42 archived tail-energy bounds | Each exceeds the independently computed conservative upper bound. Smallest energy slack: approximately 1.84836e-5. |
| Dipole normalization | A separate Cartesian-moment calculation agrees with the spherical-harmonic dipole to 1.11e-16. |
| Response evaluation | 3001 values compared with 100-digit independent Clenshaw evaluation; maximum observed discrepancy 3.541e-15. This is an additional check, not the proof of the uniform error bound. |
| Final inequality | Exact integer/rational arithmetic reproduces the manuscript fraction. |

## Independent numerical enclosures

The rebuilt conservative certificate gives:

| Quantity | Reconstructed value |
|---|---:|
| Dipole center c1 | -0.6914215317532924 |
| Dipole radius | 5.3919823283155586e-8 |
| Lower bound on v0 | 0.6914214778334691 |
| Resolved residual upper bound, ell = 3,...,151 | 0.00006647894353998339 |
| Unresolved residual upper bound, ell > 151 | 0.00005185680015812362 |
| Total residual upper bound D | 0.00011833574369810702 |
| Corrected visibility lower bound v0/(1+1.783 D) | 0.691275623770854 |

For conservative readable inequalities, these imply

\[
v_0>0.6914214778,\qquad D<0.00011833575,
\qquad v>0.69127562.
\]

The uploaded rational theorem deliberately uses weaker endpoints:

\[
v_0>\frac{6914213}{10^7},\qquad
D<\frac{1213}{10^7}.
\]

The independent calculation clears these by approximately 1.77833e-7 and 2.96426e-6, respectively. Consequently,

\[
\frac{6914213/10^7}{1+(1783/1000)(1213/10^7)}
=\frac{6914213000}{10002162779}.
\]

There is no need to promote the slightly stronger audit value into a new manuscript claim. The published-looking shortened fraction already has a margin.

## Why the integration rules do not miss an aliasing or quadrature tail

The composed response is the restriction of a Cartesian polynomial of degree at most 9 times 111 = 999. Its azimuthal Fourier degree is therefore at most 999.

For a projection through ell = 151, multiplying by the harmonic introduces azimuthal frequencies at most 1150. For total response energy, squaring introduces frequencies at most 1998. Both are strictly below 2048, so the 2048-point azimuthal rule integrates the necessary products exactly in exact arithmetic.

After selecting the matching Fourier mode m, a degree-nine-template composition has the form

\[
f_m(z)=(1-z^2)^{m/2}q_m(z),\qquad \deg q_m\le999-m.
\]

The associated Legendre function has the same factor, multiplied by a polynomial of degree ell-m. Their product is a polynomial in z of degree at most 999+ell, hence at most 1150. Likewise, azimuthally averaged squared response has polar degree at most 1998. The Clenshaw-Curtis rule of degree 2048 therefore integrates these polynomials exactly. This is an algebraic exactness argument, not a convergence assumption.

Only harmonics through 151 are individually resolved. The remainder is bounded by total energy minus a certified lower bound on resolved energy, followed by Cauchy-Schwarz. No neglected high harmonic is assumed to vanish.

## Numerical error accounting used by the independent route

The scripts contain the specific allowances and formulas. Their main ingredients are:

1. **Tables:** Machin's formula for pi and Taylor series for sine/cosine at 110 significant digits; series terms stop below 1e-115. An absolute Decimal radius of 1e-95 dominates series remainders and accumulated arithmetic error. Binary64 conversions are widened by nextafter. CC weights are recomputed from their rational cosine formula at the same precision.
2. **Normalized Legendre functions:** a local 1e-90 error allowance is propagated through the high-precision recurrence. Exact functions obey |Pbar_lm| < 13 at the degrees used; intermediate expressions are bounded, and the local allowance exceeds Decimal rounding and node error by many orders of magnitude. Largest propagated Decimal radius was 8.281e-33. A further absolute 2e-15 covers conversion to binary64.
3. **Shape evaluation:** coordinates have absolute error below 5u, with u=2^-53. Degree-nine monomial, coefficient-conversion, and length-55 dot-product errors are covered by 200u times the largest coefficient row l1 norm. This gives shape error below 1.230e-11.
4. **Response evaluation:** the exact derivative bound sum(k^2 |alpha_k|)/0.9 is below 630.606. Stability of the odd Chebyshev recurrence follows from |U_k(x)| <= k+1 on [-1,1]. Together with a 2e-11 arithmetic allowance, the derived uniform response error is below 7.773e-9; the certificate uses the larger 2e-8.
5. **Direct Fourier sums and projections:** perturb response samples by the uniform bound and add explicit length-2048 and length-2049 dot-product bounds. The Fourier arithmetic allowance is 1e-9, comfortably larger than the usual gamma_n absolute-sum estimate. Projection errors also include table, weight, Legendre-function, normalization, and summation errors.
6. **Energy subtraction and mixture:** additional 1e-9 guards for energy totals/resolved sums and 1e-11 guards for each orbit overlap dominate rounding in the remaining positive sums and scalar arithmetic. These guards prevent subtraction of nearly equal energies from silently producing an underestimated tail.

The original Legendre recurrence was also compared directly with the independent high-precision enclosure. Its largest discrepancy was 9.345e-12, and it used at most 0.003074 of the original stated arithmetic allowance. No underestimated recurrence error was found for this dataset.

These are conventional computer-assisted numerical enclosures under standard IEEE-754 arithmetic assumptions, not a proof-assistant formalization or an exact-rational computation of all spherical integrals. The final shortened visibility inequality itself is exact rational arithmetic.

## Mathematical correction step

For real orthonormal spherical harmonics, the Haar-twirled coefficient is

\[
e_\ell=\frac1{4\pi}\sum_r w_r\langle a_\ell^{(r)},b_\ell^{(r)}\rangle.
\]

Odd response functions make both marginals vanish and remove all even sectors. Thus the primary local kernel is -v0 P1 plus the unwanted odd multipoles.

The addition theorem realizes P_ell as the inner product of unit vectors in dimension 2ell+1. The universal real Grothendieck upper bound supplies a local correlation plus or minus P_ell/K for K=1.783. This is the universal constant, not an assumed improved value of K_G(3), so the correction argument is not circular. The classical bound pi/[2 log(1+sqrt(2))] is approximately 1.782213978 and is below 1.783. Public reference: Braverman, Makarychev, Makarychev and Naor, https://web.math.princeton.edu/~naor/homepage%20files/GroKri.pdf .

Mixing the primary branch with relative correction masses K|d_ell| cancels each unwanted mode. The normalization is 1+KD, yielding v=v0/(1+KD). The shared-sign symmetrization enforces zero marginals for correcting branches. Finite-setting locality extends to the full input sphere by compactness, as in the manuscript. No missing factor 2ell+1 or 4pi was found.

For an independent dipole normalization check, define the unnormalized Cartesian moments I_A=integral n A(n) dOmega and I_B analogously. Then

\[
e_1=\frac3{16\pi^2}\sum_r w_r I_A^{(r)}\cdot I_B^{(r)}.
\]

This gives -0.6914215317532925, matching the harmonic result.

## Weaknesses in the uploaded proof software

The concern about verification was justified: passing the short checker alone is not enough.

- `verify_harmonic_certificate.py` checks downstream arithmetic conditional on the supplied intervals and tail energies. It does not establish those intervals from the candidate.
- `verify_all.py` runs that conditional check; it does not run harmonic regeneration.
- The filenames and directory layout in the upload differ from those assumed by the scripts. They were restored in a separate audit directory without modifying the uploaded originals.
- The summary JSON says the raw generator is absent. That note is stale for this upload, which includes its component Python scripts. The named shell launchers and checksum manifest were not attached, but they were not needed for the verification performed here.
- The raw generator uses several generously padded numerical envelopes without fully deriving them in code comments. Its half-angle table routine also does not explicitly direct every intermediate operation outward. These are presentation/proof-software weaknesses. Independent high-precision table checks found no actual failed enclosure, and the independent direct route removes these dependencies.
- The original raw computation requires numba, which was absent here. A NumPy implementation of exactly the same odd Chebyshev recurrence was used to reproduce that pipeline. The independent route needs only Python and NumPy.

None of these weaknesses invalidated the uploaded numerical data in this audit. The archived 76 coefficient intervals and 42 tail-energy bounds were independently substantiated, not merely accepted because a script printed PASS.

## Reproduction

Extract this archive and run:

```bash
python run_audit.py
```

A clean extraction of the packaged audit was rerun successfully, with every stage passing and exactly the same final reconstructed bound. The audit environment was Python 3.12.14, NumPy 2.3.5, x86-64.

Dependencies: Python 3.11 or later and NumPy. No network access is needed. Intermediate arrays require approximately 0.7 GB and are regenerated under `work/`; they are excluded from this archive. The run writes detailed logs under `logs/`. Representative logs and small certificate JSON outputs from this audit are included under `recorded_results/`.

The original candidate and certificate are under `data/`; the source inputs and original checker scripts are unchanged except for restoring expected filenames. Independent verification scripts begin with `independent_`. `compute_harmonic_data_numpy.py` is the compatibility version of the supplied raw computation, and is not needed by the main independent run.

A SHA-256 manifest records every packaged input, script, report and recorded result. It is checked before reproduction starts.
