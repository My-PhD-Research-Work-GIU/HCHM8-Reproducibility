# Plan: functional graphs, in-degree distributions and shear-family robustness of the discretizations (2026-11-15)

Written before the registered runs, after one smoke test on profile 99 (not a registered profile) that only checked that the
script runs and measured its speed. This plan cannot be pre-registered in the checkable sense used for the 56-bit test: the
current instruction forbids commits and pushes. Instead, the SHA-256 of this file and of the script, and the time, are written
to `DATA/run_log.txt` before the first registered run.

## Questions

Q1. Which part of the discretization makes the rounded emulator F_Q behave like a random mapping: the nonlinear shears, their
    strength, the rounding rule, or the place where quantization is applied?
Q2. Is the in-degree distribution of F_Q close to that of a random mapping (Poisson(1)), and how large are its collisions?
Q3. Do the finite-state results of the paper hold for other smooth periodic shear families with the same linear core?

## Maps (lattice L_Q = (Z/2^Q Z)^8, N = 2^(8Q) states; the linear core A is that of the paper)

Shear families g_i (each g_i is smooth, 1-periodic in every argument, does not depend on z_i, and |g_i| <= 1):

- family 0, none: g_i = 0 (control);
- family 1, HCHM-8: the eight shear functions of the paper (`f_scalar` of the emulator, imported unchanged);
- family 2, single sine: g_i = sin(2 pi z_a(i) + phi_i), a = (2, 5, 7, 0, 6, 1, 3, 4), the first argument of the HCHM-8 shear;
- family 3, sine of a difference: g_i = sin(2 pi (z_a(i) - z_b(i)) + phi_i), b = (6, 3, 1, 6, 7, 0, 5, 2).

Amplitudes eps_i = s * eps_i(profile), with eps_i(profile) in [0.07, 0.15] from the profile of the paper and strength s.
Phases phi_i = 2 pi u, u from SHAKE256("R01-B3-PHASES-2026-11-15" | family | profile | i).
Profiles 0 to 9 of the finite-state section (profile 0 is the reference profile).

Discretizations:

- F_Q (end quantization; the emulator of the paper for family 1, s = 1, nearest rounding): lattice point / 2^Q -> Float64
  T = H^-1 A H (shears, core modulo 1, inverse shears) -> quantize once: nearest = floor(2^Q y + 1/2) mod 2^Q, or
  truncation = floor(2^Q y) mod 2^Q;
- G_Q (grid-exact conjugate of the paper): H_Q^-1 A_Q H_Q with each shear increment rounded to a lattice step,
  z_i <- z_i + floor(2^Q eps_i g_i + 1/2) mod 2^Q, the inverse subtracting the same increments in reverse order;
- S_Q (stepwise quantization): every shear increment quantized separately, forward z_i <- z_i + R(2^Q eps_i g_i(z)), core
  A_Q, backward y_i <- y_i + R(-2^Q eps_i g_i(y)); with nearest rounding S_Q equals G_Q except at exact ties; with truncation
  every step is still a bijection of the lattice, so S_Q is a permutation, but not one conjugate to A_Q by H_Q;
- A_Q: the core modulo 2^Q.

## Registered runs

Part 1, every lattice state, Q = 2 (N = 2^16) and Q = 3 (N = 2^24), profiles 0 to 9:

| Condition | Family | s | Maps |
|---|---|---|---|
| C0 | 0 | - | F (nearest), F (truncation), checked equal to A on every state |
| C1 | 1 | 1 | F nearest, F truncation, G, S nearest, S truncation |
| C2 | 2 | 1 | F nearest, G, S truncation |
| C3 | 3 | 1 | F nearest, G, S truncation |
| C4 | 1 | 0.25, 0.5, 2 | F nearest, G |

For every map: successor of every state (all checked to lie in L_Q), in-degree histogram, image fraction, states with two or
more predecessors, maximum in-degree, number of components (= cycles), cycle-length distribution, tail-length histogram, mean
rho (tail + cycle) over all starting states, longest tail, basin size of every cycle and the largest basin fraction.
For every G: H_Q^-1(H_Q(x)) = x and G_Q(x) = H_Q^-1(A_Q(H_Q(x))) on every state, the cycle type compared with that of A_Q, and
the period of every state compared with the A_Q period of H_Q(x).

Part 2, sampled orbits (not exhaustive), Q = 4 (32 state bits) and Q = 5 (40 state bits), profiles 0 to 9, starting states
0 to 7 of the uncapped measurement (`initial(Q, run)` of the emulator), Brent's algorithm without cap: F nearest for families
1, 2, 3 at s = 1; F truncation for family 1 at s = 1; F nearest for family 1 at s = 0.5 and 2.

Independent checker: for Q = 2, family 1, s = 1, profiles 0 and 1, a separately written numpy implementation (vectorized
evaluation, its own cycle search) recomputes the successors of F, G and A and their in-degree histograms and cycle types.

## Predictions and criteria (fixed here)

- H1 (implementation check, every G of Part 1): H_Q is a bijection with the stated inverse, G_Q = H_Q^-1 A_Q H_Q on every
  state, G_Q is a permutation with the cycle type of A_Q, and every state has the A_Q period of its image under H_Q.
  A failure is an implementation error.
- H2 (in-degree, C1, F nearest, Q = 3, means over the ten profiles): the frequencies of in-degree 0 to 4 lie within 0.01 of
  Poisson(1) (0.3679, 0.3679, 0.1839, 0.0613, 0.0153), and the image fraction lies in [0.60, 0.66].
- H3 (shear families, C2 and C3 at Q = 3, F nearest): image fraction in [0.60, 0.66] and mean rho / sqrt(pi N / 2) in
  [0.7, 1.3] (means over profiles); for families 1, 2 and 3, the mean over profiles of the slope of log2(mean rho) on
  log2 N over Q = 2 to 5 lies in [0.45, 0.55].
- H4 (rounding rule, C1, Q = 3): the image fractions of F with nearest rounding and with truncation differ by less than 0.01
  (means over profiles), and both rho ratios lie in [0.7, 1.3].
- Exploratory, no criterion: the strength series C4 (with C0 as s = 0), the stepwise discretization S (agreement of S nearest
  with G; for S truncation the number of cycles against the random-permutation value H_N = ln N + 0.5772 and the mean cycle
  length of a random state against (N + 1) / 2), and every Q = 2 result.

Replication unit: the profile. Starting states within a profile are not independent units. Intervals: percentile bootstrap
over the ten profiles (10,000 resamples). No profile or condition is excluded; a crash is reported and the condition is run
once more with the same inputs. Every outcome is reported, whether it meets its criterion or not.

## Budget

About one to two hours on the laptop (16 threads, numba). Q >= 4 is not enumerated exhaustively (2^32 states and more).

## Not covered

Integer-only evaluation of the shear functions (all variants evaluate g_i in Float64); state sizes above 40 bits for the new
families; any cryptographic use of these maps.
