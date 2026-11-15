# Plan: a grid-exact conjugate of HCHM-8 against the rounded emulator (written 2026-11-10, before the registered runs)

## Question

The finite-state section of the paper measures the rounded emulator F_Q (lattice to Float64, T = H^-1 A H, rounding to the
nearest lattice point) and finds orbit lengths growing like sqrt(N), as for a random mapping, while the linear core alone
repeats after 3 * 2^(Q-1) iterations. Which property of the discretization produces the random-mapping behaviour: the
nonlinear shears, or the loss of bijectivity by rounding?

## Construction

G_Q = H_Q^-1 o A_Q o H_Q on the lattice L_Q = (Z / 2^Q Z)^8, where H_Q applies the eight shears of the paper in the same
order with each increment rounded to a lattice step, z_i <- z_i + round(2^Q eps_i g_i(z / 2^Q)) mod 2^Q, H_Q^-1 subtracts the
same increments in reverse order, and A_Q is the integer core modulo 2^Q. Because g_i does not depend on z_i, each lattice
shear is a bijection of L_Q, so H_Q is a bijection and G_Q is conjugate to A_Q. Consequently (proved, not tested):

- G_Q is a permutation of L_Q (no tails), and its cycle type equals that of A_Q;
- the period of the G_Q orbit of x equals the period of the A_Q orbit of H_Q(x), and divides 3 * 2^(Q-1) for Q >= 3
  (12 for Q = 2), the order of A modulo 2^Q established in the paper.

## Runs (fixed here)

Script: `CODE/grid_conjugacy.py <out> --profiles 0:10 --exhaustive-q 2,3 --orbit-q 4,5,6,7`. F_Q is `step()` of
`08_FUTURE_WORK/2026-10-11_FINITE_STATE_CYCLE_SCALING/CODE/cycle_scaling_profiles.py`, imported unchanged, with its
`profile_params()` and `initial()`.

- Part 1, every lattice state, Q = 2 (N = 2^16) and Q = 3 (N = 2^24), profiles 0 to 9: for F_Q, G_Q and A_Q the image
  fraction, the number of cyclic states and of cycles, the mean over all N starting states of tail plus cycle (rho), the
  longest tail, the cycle type of G_Q and A_Q, and the fraction of states on which F_Q and G_Q agree.
- Part 2, Q = 4, 5, 6, 7, profiles 0 to 9, the starting states of the uncapped measurement (32, 32, 13 and 2 per profile): the
  period of the G_Q orbit and of the A_Q orbit of H_Q(x0), with a cap of the core period plus one.

## Predictions and criteria

1. Theorem check (implementation test). In every case of Parts 1 and 2: G_Q is a permutation, its cycle type equals that of
   A_Q, and every period equals the corresponding A_Q period and divides the core period. A failure would be an
   implementation error and is reported as such.
2. Random-mapping signature of F_Q (empirical). For a random mapping of N points the expected image fraction is
   1 - 1/e = 0.632, and the expected number of cyclic points and the expected rho of a random start are both about
   sqrt(pi N / 2). Criterion, fixed here, at Q = 3: the mean over the ten profiles of the image fraction lies in [0.60, 0.66],
   and the mean over the profiles of mean rho / sqrt(pi N / 2) lies in [0.7, 1.3] (the band of the 56-bit test). Q = 2 is
   reported descriptively (the lattice has four points per coordinate, so the shear increments are mostly rounded to zero).
3. Descriptive: agreement of F_Q and G_Q, the number of cyclic states and cycles of F_Q, and the cycle types.

Whatever the outcome, it is reported. Both the plan and the script are committed and pushed to GitHub before the registered
runs start; the commit is opened as a pull request so that the hosting service records the time of the first workflow run.

## Smoke test before this plan

The script was run once to check that it executes, on profile 99 (not among the registered profiles 0 to 9) with Part 1 at
Q = 2 and Part 2 at Q = 4 and 7. Result, disclosed here because it was seen before the plan was written: G_Q a permutation
with the cycle type of A_Q; all periods equal and dividing the core period; for F_Q an image fraction of 0.608, rho ratio
0.70, cyclic-state ratio 0.50. A separate check with all shear amplitudes set to zero gave F_Q = G_Q = A_Q on every state.

## Not covered

An integer-only evaluation of the shear functions (G_Q still evaluates g_i in Float64, so G_Q itself may differ between
computing environments; its cycle type cannot, being that of A_Q), Q >= 4 exhaustively, and any cryptographic use of G_Q.
