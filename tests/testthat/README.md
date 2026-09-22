# Binary uptake robustness tests

From D:/prince_BART2026, run:

```r
testthat::test_local(".", filter = "binary-robustness")
```

This loads the source checkout, including the test helpers. Running the test
file with source() alone is not equivalent. The existing tests/testthat.R also
includes these tests during package checks.

There are seven test blocks:

1. A two-chain fit: array dimensions, probability bounds, and principal-stratum
   memberships constrained by the observed binary instrument and uptake.
2. A one-chain fit: the mixed effect retains its chain dimension and equals the
   complier-probability-weighted contrast for each draw.
3. Invalid binary values, mismatched input lengths, and propensity boundaries.
4. Missing Y, Z, W, or propensity: informative validation errors.
5. Constant outcomes (all zero or all one).
6. Perfect observed compliance (W = Z), which creates an empty initial
   noncomplier subset.
7. Bayes-rule calculations at interior and boundary outcome probabilities.

Each fit uses 96 deterministic observations, five trees, three warmup iterations,
four retained draws per chain, and a supplied propensity of 0.5. All fits run
sequentially with a fixed seed. The helper restores the previous future plan and
R random seed. Binary sampling currently has internal thinning of 20, so four
retained draws do not mean only four underlying tree updates.

These are execution and arithmetic checks, not tests of MCMC convergence,
credible-interval coverage, or effect recovery. No external dataset is needed.
They exercise binary uptake only.

The tests assert desired behavior and can fail on existing bugs. Source
inspection suggests failures in the one-chain reduction, missing-value
validation, perfect-compliance initialization, and boundary Bayes calculation.
Do not invert these assertions or skip the cases merely to obtain a green run.

The test suite was added without executing it. No model implementation was
changed as part of this addition.
