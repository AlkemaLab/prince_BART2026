# princeBART robustness and implementation audit

Date: 2026-09-22

## Scope and limitations

The configured checkout, C:/Users/bioRabbit/OneDrive - University of Massachusetts/Leontine Lab/prince_BART-main, was unavailable to this session. This review therefore covers the installed princeBART 0.2.0 snapshot, built 2026-09-14 05:36:58 UTC, in C:/Users/bioRabbit/AppData/Local/R/win-library/4.5/princeBART. Installed dbarts is 0.9.34.

The installed functions retain source references and source text, including roxygen comments. Locations below refer to the original R source files recovered from that snapshot, not a claim that the current IDE files were inspected. Newer fork changes may already address some findings.

I inspected input preparation, both samplers, posterior estimands, tree prediction, generalization/transportability, survey bootstrap, sensitivity optimization, segmentation, S3 methods, dependencies, and installed examples. I ran two isolated diagnostic scripts and codetools::checkUsagePackage(). I also loaded the bundled ordinal fit and exercised coef(). I could not inspect the checkout's tests, run its R CMD check, or verify all paper-to-code equivalences. This is not a formal statistical validation of the package.

No installed package or project files were modified. The accompanying reproduce.R and followup.R are diagnostic probes, not a production test suite: expected failures are caught and printed so later probes run. Both completed with process exit code 0. The sensitivity probe changes only a local copy's function environment to use an available solver; it does not modify the package namespace.

Evidence labels:
- Reproduced: an executable counterexample demonstrated the behavior.
- Source: directly established by inspected code, but not necessarily exercised end-to-end.
- Design: an estimand/model contract needs a decision or derivation; not automatically a programming error.

## Executive assessment

There are correctness problems, not just opportunities for cleaner code. Highest priorities are ordinal grid alignment and latent-pair probabilities, prediction-scale consistency, survey aggregation, and posterior-array dimensions. These can return plausible-looking numerical results for the wrong calculation.

Do not use affected ordinal fits or target-population/sensitivity results as validated scientific results until the relevant fixes and regression tests pass. This does not establish that every previously obtained estimate is wrong; impact depends on the path and data. Refactoring binary and ordinal into one sampler should wait until their model contracts are explicit.

## 1. Ordinal sampler and estimand

### 1.1 Outcome predictions are indexed by person instead of expanded grid row

Priority: critical. Evidence: reproduced and source.
Location: R/fit_psbart_ordinal.R:477, :478.

The outcome samplers predict on X_w_grid, which has multiple candidate (w0,w1) rows per person. build_strata_posterior_dt() assigns mean_y0[id] and mean_y1[id], although these vectors are already grid-row-aligned. The id vector indexes people, not prediction rows.

A four-row grid with id = c(1,1,2,2) and predictions c(.1,.2,.7,.8) becomes c(.1,.1,.2,.2). Candidate strata for a person share an outcome likelihood; it cancels from their relative posterior probabilities. Later people can receive predictions from another person's grid rows. The sampled stored outcome means inherit the same error.

Recommendation: assign predictions directly in grid order; explicitly distinguish unit_id from grid_row_id; assert prediction length and ordering; test that changing Y changes stratum odds when candidate strata have different outcome means.

### 1.2 Pair-probability calculation does not match latent-state sampling

Priority: critical. Evidence: reproduced and source.
Locations: R/fit_psbart_ordinal.R:435, :590, :599.

sample_latent_z_states() uses bounds derived from a_j() and bc(). p_wpair() uses a separate, inconsistent construction:

- For uptake zero, latent sampling uses (-Inf,1) at lambda=1, but p_wpair() uses an upper bound of zero.
- Coordinates are ordered as (w0,w1), but means/variances are paired as (mu_z1,mu_z0).
- p_wpair() has no lambda or w_max parameters, so it cannot match transformed bounds or the open-ended highest category.

At independent standard-normal means, P(W0=0,W1=0) is .25 in p_wpair(), versus pnorm(1)^2 = .707861 under the sampler's bounds. With unequal means, a nonzero-bin example gives .007304857 instead of .046390499.

Recommendation: use one shared bounds function and one unambiguous coordinate order. Test zero, interior, and top categories; unequal means/variances; transformed bounds; Monte Carlo agreement; and total mass over the full grid. If monotonicity truncates support, test normalization over permitted pairs separately.

### 1.3 Nonzero rho is not propagated through stratum membership

Priority: high. Evidence: source.
Locations: R/fit_psbart_ordinal.R:431, :467, :482, :554.

rho is used for correlated latent-state draws, but build_strata_posterior_dt() hardcodes rho=0 in p_wpair(). A nonzero-rho run consequently mixes dependent latent-state sampling with independent pair probabilities. There is also no clear entry validation requiring a finite scalar strictly between -1 and 1; conditional standard deviations involve sqrt(1-rho^2).

Recommendation: propagate the realized rho consistently; derive/check every affected conditional, including latent-response tree updates, before advertising sensitivity to rho. If only independence has been validated, reject nonzero rho explicitly until support is correct.

### 1.4 Ordinal covariates are scaled twice, but prediction reuses only the first recipe

Priority: high. Evidence: reproduced and source.
Locations: R/mc_psbart.R:352, :409; R/fit_psbart_ordinal.R:302; R/general_bart.R:507.

The public function scales X and appends qnorm(propensity). prepare_ordinal_data() scales the entire matrix again, including the propensity feature. Stored scaling reflects only the first operation; external prediction applies that first recipe and appends an unstandardized probit score.

Second scaling changes e even without trimming. Trimming also changes the means/variances used for other predictors. A perfectly legitimate known constant propensity of .5 yields e=0, then NaN after the ordinal second scaling.

Recommendation: fit exactly one preprocessing recipe; store and reuse it verbatim. Explicitly record probability versus probit scale. Handle constant columns without division by zero. Test training-data predictions through the external prediction path.

### 1.5 Effect orientation is inconsistent with the additional-child interpretation

Priority: high interpretation risk. Evidence: reproduced, source, and installed ordinal example.
Locations: R/estimands.R:285; R/fit_psbart_ordinal.R:396; installed extra/ordinal.html:386.

Outcome component y0 corresponds to Z=0 and y1 to Z=1. The ordinal affected group is W0-W1=1, so increasing Z reduces uptake by one. The implemented m_y1-m_y0 is the effect of that reduction. It is the negative of the effect of one additional unit of uptake under the documented interpretation.

Synthetic W0=2, W1=1, m_y0=.8, m_y1=.2 produces -.6. The corresponding extra-unit contrast is +.6.

Recommendation: define whether the API reports an instrument-induced reduction or an additional-unit treatment effect. Either can be a legitimate estimand, but labels, signs, plots, segmentation, and transportability must agree. Explicitly name conditioning events and distinguish potential outcomes indexed by instrument from those indexed by uptake.

## 2. Binary sampler and posterior summaries

### 2.1 Empty latent strata can leave stale likelihood data in active samplers

Priority: high. Evidence: reproduced and source.
Locations: R/fit_psbart_binary.R:264, :332.

update_samplers() skips setData() for empty classes; subsequent run_all_samplers() still advances those samplers using their previous datasets. That conditions on outdated memberships, rather than implementing a no-observation update.

The condition sum(co)>0 also gates both the all-person co model and the noncomplier atnoco model. An all-noncomplier probe updated only y0nt and y1at: neither co nor atnoco was refreshed, despite their current responses being well-defined. Perfect-compliance initialization (W=Z, 20 persons) failed with 'length of y must be greater than 0'.

Recommendation: separate masks for each component, always update the full-person class regression, and implement explicit prior/no-likelihood transitions for empty outcome components. Do not keep old likelihoods or silently add artificial observations. Diagnose missing instrument arms and weak/absent stratum information before fitting.

### 2.2 Bayes-rule arithmetic can return NaN for a valid observed-outcome likelihood

Priority: high. Evidence: reproduced.
Location: R/fit_psbart_binary.R:327.

compute_posterior_class_prob() forms both Y=1 and Y=0 fractions, then combines them using Y and 1-Y. For Y=1, equal class priors, and both outcome probabilities equal to 1, the correct posterior is .5; the code returns NaN because the irrelevant branch is 0/0 and zero times NaN remains NaN.

Recommendation: select the observed likelihood before normalization, preferably on the log scale. Distinguish numerical underflow from genuinely impossible observations. Add tests at probabilities 0, 1, and near each boundary.

### 2.3 get_mix_tau() silently treats units as chains when there is one chain

Priority: high. Evidence: reproduced.
Location: R/estimands.R:118.

Later default-dropping array slices undo the protection of earlier drop=FALSE slices. With 4 iterations, 1 chain, and 3 persons, the returned outcome summaries have shape (4,3,5), not (4,1,5). The three person-specific effects .2,.4,.6 appear as three chains; the proper aggregated effect is .4. The same test with two chains gives the expected dimensions.

Recommendation: centralize array extraction/reduction helpers with explicit dimension contracts. Test iterations=1, chains=1, units=1, and one-person subsets independently. Review similar dropping operations in imputation and segment extraction. Zero posterior mass for a stratum also produces nonfinite mixed estimands; return a documented unavailable result and diagnostic instead of unqualified numeric output.

### 2.4 Sample estimands mix never-takers and always-takers and ignore a subset for proportions

Priority: high. Evidence: reproduced and source.
Locations: R/estimands.R:174, :189; R/impute.R:25.

get_sample_tau() computes stratum proportions before the treated subset is applied. A subset containing only a never-taker and an always-taker should have proportions (0,.5,.5); output remained (1/3,1/3,1/3).

For the never/always outcome summaries it removes compliers, but does not restrict each summary to the intended remaining class. A synthetic never-taker outcome of 0 and always-taker outcome of 1 produces .5 for both labels.

The supplied imputation helper uses complier outcome models and holds observed outcomes by Z, appropriate for the complier calculation. That is not a general never/always potential-outcome imputer.

Recommendation: apply the same subset before all calculations, use class-specific masks and appropriate observed-uptake consistency, and omit unsupported summaries until correctly defined. Consider storing imputed draws: summary(type='sample') and coef(type='sample') currently can re-impute independently.

## 3. Prediction and generalization

### 3.1 Binary propensity predictions receive pnorm() twice

Priority: high. Evidence: reproduced with real dbarts and source.
Locations: R/general_bart.R:223, :499, :778, :1194.

For the installed dbarts, default predict() for a binary bart fit already returns probabilities. The package wraps this output in pnorm(). A test verified predict(default) equals pnorm(predict(type='bart')) exactly. Probabilities spanning approximately .000037 to .999996 were transformed into .500015 to .841344.

This changes target propensity features and overlap estimates. It is not a harmless rescaling because the source outcome trees were trained on another feature scale. Some ordinal overlap paths similarly apply pnorm() to continuous-regression predictions of probabilities (around lines 845 and 1254).

Recommendation: specify the prediction scale explicitly; apply inverse links exactly once. Test source/target feature parity and overlap calibration.

### 3.2 Extracted continuous tree sums are not on the latent response scale

Priority: high. Evidence: reproduced with real dbarts and source.
Locations: R/general_bart.R:540, :625.

predict_one_sample_raw() simply sums exported terminal-node values. For continuous-response dbarts fits those values are on its internal scaled response scale, whereas sampler$run()$test is on the response scale. Ordinal prediction uses the raw sums directly as mu_z0 and mu_z1.

In the follow-up test, raw sums ranged from -.4561 to .3859 while native test predictions ranged from 55.3615 to 121.5077; maximum discrepancy was 121.1218.

Recommendation: retain and apply each component's response transformation, potentially draw-specific when data/responses change, or use a supported native prediction mechanism that preserves sampler state. Require saved-tree versus native prediction equivalence tests for binary and continuous components.

### 3.3 Ordinal transport prediction substitutes one pair instead of integrating over the targeted stratum

Priority: high statistical-design concern. Evidence: source.
Location: R/general_bart.R:529.

The external ordinal path floors/transforms predicted latent means to a single w0_hat,w1_hat and predicts the outcome difference at that pair. It does not integrate residual latent uncertainty or explicitly condition that pair on W0-W1=1. This is not in general the same quantity as an affected-unit average used in the fitted ordinal summaries.

Recommendation: specify the target estimand algebraically, then integrate outcome contrasts against the appropriate posterior pair distribution. Retain the parameters required for that distribution. Test against direct enumeration in a tiny ordinal example. This issue is separate from, and remains after fixing, the raw-tree response scale.

### 3.4 Overlap tree prediction can collapse different chains into one

Priority: high. Evidence: reproduced and source.
Locations: R/general_bart.R:803, :1213.

These prediction loops select iteration without selecting chain. Tree IDs repeat across chains; the single-tree traversal then uses the first matching tree's structure/leaves. A two-chain probe with latent leaf means -1 and +1 returned .1586553 instead of the mean probability .5.

Recommendation: key trees by component, chain, iteration, and tree. Average probabilities only after distinct draws are predicted.

### 3.5 Prediction interfaces disagree about tree schema and fail for one target person

Priority: medium/high. Evidence: reproduced and source.
Locations: R/predict.R:20; R/general_bart.R:613; R/fit_psbart_binary.R:141; R/fit_psbart_ordinal.R:247.

Exported predict_trees() expects a sample column, while fitted tree assembly removes that column and supplies chain and iteration. Passing that fitted-style schema yields 'Non-numeric argument to mathematical function'. A model component must also be selected; the full tree table contains different regressions.

predict_one_sample() uses sapply() followed by rowSums(); one target row simplifies to a vector and errors. Tree traversal uses positional split variables, so predictor order and encoded levels must be validated, not merely column presence.

Recommendation: provide predict.prince_bart() with a documented component/type contract, reconstruct the saved design matrix in training order, and keep a fixed matrix result for singleton inputs. Reject unknown levels and unavailable trees clearly.

### 3.6 The fitted propensity model is not retained as part of the prediction recipe

Priority: model consistency concern. Evidence: source.
Location: R/general_bart.R propensity-fitting paths; R/mc_psbart.R:369.

Generalization refits propensity rather than predicting with the source fit used to construct its training feature; some paths use only predictors shared with the target. Thus the generated feature can change meaning in addition to the double-link error.

Recommendation: retain the source propensity fit and feature schema, or require a compatible externally supplied prediction function. State whether propensity uncertainty is treated as plug-in or propagated.

## 4. Survey bootstrap and sensitivity analysis

### 4.1 Survey weights are lost within PSUs, and subpopulation denominators count excluded people

Priority: high. Evidence: reproduced.
Location: R/general_bart.R:901.

survey_pate() computes an unweighted mean outcome within each PSU, then weights it using mean(survey weight) times the number of all PSU members. This only recovers the desired individual-weighted calculation under restrictive within-PSU weight conditions.

One PSU, effects (0,1), weights (1,9): actual .5, correct weighted mean .9.

Subpopulation masking sets excluded outcomes to NA but leaves their counts/weight mass in PSU aggregation. With equal eligible weight in two PSUs, an excluded nine-person addition to one PSU moved the bootstrap mean to .8217455 rather than approximately .5.

For posterior draw d and one shared bootstrap multiplier b[q,d] per PSU, a suitable individual-weighted domain mean is:

    sum_q b[q,d] * sum_{i in q AND domain} w[i] * tau[i,d]
    ---------------------------------------------------
          sum_q b[q,d] * sum_{i in q AND domain} w[i]

The actual bootstrap construction must also reflect the intended sampling design; pooling all PSUs is not automatically correct for every stratified complex survey.

Recommendation: implement the weighted numerator and eligible denominator directly, validate finite nonnegative weights with positive eligible mass, PSU IDs and lengths, and explicitly document supported designs. Uniform-spacing bayesian_bootstrap() itself is a valid way to generate Dirichlet(1,...,1) weights for positive integer n; the confirmed error is their subsequent use.

### 4.2 Sensitivity weighting is not the same functional as ordinary aggregation

Priority: high. Evidence: reproduced.
Locations: R/general_bart.R:980, :1015, :1035.

find_shift_weights() returns final weights w_i*r_i. compute_shift_pate() multiplies outcomes by those weights, takes a PSU mean, and applies the original PSU weights again. It does not calculate a normalized mean under the optimized final weights. The optimized objective and reported aggregation also use different weighting structures.

A local solver-enabled diagnostic at gamma=1, the no-shift boundary, gave baseline .25 and sensitivity .65 for outcomes (.1,.4), weights (1,3), one PSU. The correct individual-weighted value is .325. The public API currently restricts gamma>1; using 1 here is an internal identity test, not a claim that the public API accepts it.

Recommendation: distinguish multipliers r_i from final weights, use one correctly normalized aggregation functional for baseline and bounds, align the optimization objective with the quantity reported, and enforce the no-shift and common-weight-rescaling invariants. Avoid rounding optimized multipliers before checking feasibility.

### 4.3 An optional optimization path assumes an unavailable solver

Priority: medium. Evidence: reproduced.
Location: R/general_bart.R:1027.

CVXR was installed with CLARABEL, SCS, OSQP, and HIGHS available, but ECOS was unavailable. The hardcoded solver='ECOS' errors before the status fallback. CVXR also emitted deprecation warnings for the value-extraction interface. A fallback to uniform weights is not a valid successfully computed sensitivity bound.

Recommendation: detect supported installed solvers or accept an explicit solver argument; catch solver errors; return an explicit unavailable/failed result rather than misleading bounds. Test the supported dependency versions. Consider making CVXR a Suggests dependency if sensitivity is genuinely optional.

## 5. Input, object, and segmentation contracts

### 5.1 Separate model frames can silently pair different people's variables

Priority: high. Evidence: reproduced.
Location: R/utils.R:26.

parse_psbart_formula() constructs separate model frames with independent missing-row omission. Length checks cannot detect different omitted rows. A six-person example with different missing rows in Y, X, Z, W produced length 5 for every component, yet X rows were 1,3,4,5,6 and Y rows were 2,3,4,5,6.

Recommendation: build one all-variable frame with na.fail or one joint complete-case mask. Carry stable row IDs through encoding, propensity handling, and trimming. Never align observations by vector length alone.

### 5.2 Validation accepts nonfinite data and inconsistently coerces factors

Priority: medium/high. Evidence: reproduced and source.
Locations: R/fit_psbart_binary.R:179, :190, :200; R/mc_psbart.R:343.

Confirmed: validate_binary() accepts NA; validate_propensity() does not reject NA; constant X produces NaN during scaling; factor(c(0,1,0)) is converted to internal codes and rejected; infinite ordinal uptake gives a generic missing TRUE/FALSE error.

Additional source-level gaps deserving tests: scalar/integer constraints on iteration and worker counts; rho bounds; ordered finite overlap bounds; all rows trimmed or one instrument arm remaining; unnamed/duplicate feature names; collision with generated e, w0, w1 and reserved data columns; direct-interface categorical encoding versus formula encoding.

Recommendation: common validators should reject or explicitly handle missingness, nonfinite values, invalid scalar lengths and ranges. Accept binary factor levels intentionally or reject factors with a clear conversion instruction. Preserve a reusable factor-level/design-matrix schema. Treat zero-variance columns deliberately.

### 5.3 Segment display text is used as a data identifier

Priority: medium/high. Evidence: reproduced at helper level and source.
Locations: R/segment_heterogeneity.R:340, :415.

assign_segments() uses factor(nodes, labels=pretty_rules). Different rules can round to the same display string; duplicate factor labels merge levels. The helper maps x<.11 and x<.12 to the same x<.1 label. This demonstrates the collision mechanism, not a claim that every fitted tree contains colliding full paths.

The formatter also changes actual inequalities: age<=40 becomes age<40.0, and age>40 becomes age>=40.0.

Recommendation: retain immutable terminal-node IDs, use labels only for display, and preserve original inequalities/precision. Reconstruct categorical rules from saved encoding metadata instead of name-prefix heuristics.

### 5.4 Segment ranking and uncertainty need clearer definitions

Priority: medium; partly design.
Locations: R/segment_heterogeneity.R:184, :200, :241, :383.

Highest/lowest segment selection uses an unweighted mean cate, whereas reported segment effects use posterior group weights. Those can rank segments differently. Empty-draw summaries return ci90 but the caller expects ci, risking malformed output.

For ordinal segmentation, marginal posterior mean contrasts and marginal affected probabilities are not generally interchangeable with the conditional affected-unit contrast. Post hoc learned segments also require an explicit statement about whether partition-selection uncertainty is included.

Recommendation: use the same estimand for tree construction/ranking/reporting; unify empty-result schemas; label conditional-on-partition uncertainty accurately.

### 5.5 Public arguments and stored state do not fully describe the computation

Priority: medium. Evidence: reproduced and source.
Locations: R/princebart-class.R:231, :294; chain assembly in R/mc_psbart.R.

Ordinal summary/coef accept type and treated_only without applying them. In a synthetic fit, coef(type='nonsense', treated_only=TRUE) equaled default coef(). The summary message claims these options are interpreted, which is misleading.

Binary and ordinal thinning defaults differ; ordinal trees are saved on an internal ten-iteration schedule. Short ordinal chains may save none. combine_chain_trees(list(list(trees=NULL)), TRUE) returned a one-column chain object instead of NULL or an empty valid tree table.

Realized rho, transform settings, latent residual/response-scaling information, and the exact ordinal support are not all preserved in a stable model configuration. A stored call can contain unevaluated symbols and is not a substitute.

Recommendation: validate or reject unsupported arguments; expose/document thinning and tree-save schedules; validate tree availability before downstream use; use a versioned fit-object schema with realized configuration, preprocessing, row IDs, draw schema, and prediction metadata.

## 6. Naming and documentation

Naming differences are not automatically bugs. They matter when one name hides different scales, objects, or conditioning sets.

| Current convention | Problem | Suggested contract |
| --- | --- | --- |
| princeBART / prince_BART() / prince_bart | Package, function, and class use different conventions | Document the distinction; keep compatibility aliases if standardizing |
| co/c, nt/n, at/a | Different class abbreviations across samplers, arrays, and estimands | Use one canonical class vocabulary and named conversion helpers |
| Z versus z0/z1 | Observed instrument versus continuous latent uptake variables | Distinguish Z, W0/W1, and latent_W0/latent_W1 |
| y0co/y1co versus y0/y1; m_y0c versus m_y0 | Outcome component labels and conditioning differ by pipeline | Define component metadata, not substring inference |
| e / propensity | Probability in some places, probit transform elsewhere | propensity_prob and propensity_probit |
| X_model / X_raw / X | Stored X_model is unscaled encoded X plus probit e, not the exact fitted matrix | Explicit raw, encoded, transformed matrices or a saved recipe |
| sample / iteration / chain | Fitted trees and predictor disagree on draw identifier | Fixed chain + iteration keys and optional unique draw_id |
| weights / shift_wts / r | Survey weights, final weights, and multipliers conflated | survey_weight, weight_multiplier, final_weight |
| k / n_trees / n_thin / n_samples | Shrinkage, tree count, thinning, retained draws must remain distinct | Define units and defaults for both pipelines |
| treated_only and printed Z=1 labels | Instrument assignment is not necessarily actual treatment uptake | Name the exact selected event; reserve treated for W-based events |

Specific documentation defects:

1. R/mc_psbart.R:280 uses data=mydata without defining or bundling that object. There is no exported dataset generator in the installed namespace. extdata contains fitted/result RDS files, not a documented mydata dataset. Installed HTML examples contain simulation code, which could be factored into a small supported generator.
2. R/princebart-class.R examples similarly refer to df; generalization examples refer to source_data/survey_data. dontrun prevents routine example execution but does not make examples reproducible.
3. R/mc_psbart.R:267 promises data$X as a compatibility alias; constructed and bundled objects contain X_model, X_raw, Y, Z, W, e, with no X. The description of X_model as the matrix used for fitting also conflicts with storing unscaled encoded covariates.
4. The introduction computes Y0=plogis(lin0), Y1=plogis(lin0+tau), but labels mean(tau[complier]) as true_ate_c (installed extra/introduction.html:480). tau is a log-odds shift. For the probability-scale mixed effect, the reference truth should be mean((Y1-Y0)[complier]). A realized sample estimand needs realized potential outcomes, not just probabilities.
5. The ordinal one-additional-child interpretation conflicts with the current contrast sign; ordinal summary claims type/treated_only are interpreted when they are ignored.
6. Ordinal estimand comments describe a fixed grouping up to level 5 although the implementation uses adaptive/manual thresholds. Formulas must clearly distinguish an adjacent-level contrast from a fixed pair with ambiguous w0/w1 labels.
7. Describe Y~X|Z|W as this package's formula convention, not interchangeable ivreg/2SLS syntax.
8. Document actual parallel-plan precedence: workers=1 does not necessarily override an existing nonsequential future plan. Unused n_cores arguments in helpers should be implemented or deprecated.
9. Some roxygen strings have doubled backslash markup; render and inspect the generated help rather than checking source comments alone. Regenerate reference pages and cached outputs from executable examples after numerical fixes.

## 7. Maintainability and recommended implementation order

A small wrapper around dbarts can be justified, but the existing dbarts_binary() reaches into private internals and has unused retrieved helpers. codetools reported unused validateArgumentsInEnvironment, redirectCall, addCallArgument, parsePriors, and setDefaultsFromFormals in R/utils.R:145. Other usage warnings in segmentation are largely consistent with nonstandard evaluation, not proof of runtime errors.

Prefer supported public constructors where they meet requirements; otherwise isolate one compatibility adapter, document why it exists, constrain/test dependency versions, and compare its behavior with native dbarts. Avoid distributing private API assumptions across samplers and prediction.

Suggested order:

1. Fix ordinal grid alignment, shared pair bounds, rho consistency, and binary stale-data/likelihood arithmetic.
2. Define and test sign, stratum, and subset contracts; repair posterior-array dimensions and sample estimands.
3. Save one preprocessing/prediction recipe and validate native versus exported-tree predictions.
4. Correct survey/domain aggregation, then rederive sensitivity calculations using the same functional.
5. Add input/schema validators and fix segmentation IDs, unsupported arguments, and empty-tree handling.
6. Make examples executable, correct simulation truth, regenerate documentation, and run the real checkout's test suite/R CMD check.

Do not first collapse the two samplers solely because binary uptake has two categories. Shared preparation, draw storage, validation, and prediction interfaces are good candidates for reuse, but the current binary model directly models principal-class probabilities whereas the ordinal model discretizes latent uptake regressions. Equal support does not establish equal priors, likelihood parameterization, or posterior inference.

## 8. Minimum regression-test matrix

| Area | Essential cases/invariants |
| --- | --- |
| Inputs | NA in different formula terms; Inf; binary factors; constant columns; one arm; empty data; generated-name collisions |
| Arrays | iterations/chains/units each equal 1 separately; one-person subsets; stable dimnames; zero stratum mass |
| Binary sampler | all/none compliers; empty never/always class; stale-data prevention; Bayes-rule probability boundaries |
| Ordinal sampler | grid-row versus person alignment; zero/top bins; unequal means/SDs; rho 0/nonzero; mass and simulation agreement |
| Prediction | native versus extracted trees; probability/latent scales; same training recipe; reordered columns; unseen levels; multiple chains |
| Estimands | known synthetic contrasts and signs; subset-specific proportions; separate never/always means |
| Survey | one PSU with unequal weights; constant effects; common-weight-rescaling invariance; excluded-domain invariance |
| Sensitivity | gamma=1 identity; nested bounds; optimizer feasibility; objective/report agreement; unavailable solver |
| Segments | distinct node IDs despite identical display labels; exact boundaries; weighted ranking; empty draws |
| Public API/docs | unsupported arguments error; short keep_trees runs; runnable data generator; correct simulated truth |

The two accompanying scripts provide starting counterexamples for many of these tests. They should be converted to assertions against the corrected behavior, then run against the actual source checkout and supported dependency versions.
