# Read-only diagnostics against the installed package. No namespace mutation.
library(princeBART)
ns <- asNamespace("princeBART")
f <- function(name) get(name, ns)
check <- function(label, expr) {
  cat("\n", label, "\n", sep = "")
  tryCatch(print(force(expr)), error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}
cat("princeBART", as.character(packageVersion("princeBART")), "dbarts",
    as.character(packageVersion("dbarts")), "\n")

make_probs <- function(chains = 2L, iterations = 4L, n = 3L) {
  a <- array(0, c(iterations, chains, 6, n),
    dimnames = list(iteration = NULL, chain = NULL,
      variable = c("p_a", "p_n", "m_y0c", "m_y1c", "m_y0n", "m_y1a"), unit = NULL))
  a[, , 1:2, ] <- .2
  a[, , 3, ] <- .1
  a[, , 5:6, ] <- .4
  for (i in seq_len(n)) a[, , 4, i] <- .1 + .2 * i
  a
}
check("mixed array two chains: expected outcomes 4 x 2 x 5", {
  r <- f("get_mix_tau")(make_probs()); list(dims = lapply(r, dim), tau = r[[2]][, , 5])
})
check("mixed array one chain: expected outcomes 4 x 1 x 5, tau=.4", {
  r <- f("get_mix_tau")(make_probs(1)); list(dims = lapply(r, dim), tau = r[[2]][, , 5])
})
check("mixed zero-mass stratum", {
  p <- make_probs(); p[, , 1:2, ] <- .5
  r <- f("get_mix_tau")(p); sum(!is.finite(r[[2]]))
})
check("NA response accepted", f("validate_binary")(c(0, 1, NA), "Y"))
check("NA propensity accepted", f("validate_propensity")(c(.2, NA), 2))
check("constant X scaled to nonfinite", f("validate_and_prepare_X")(cbind(x = 1:4, constant = 1)))
check("factor binary coercion", f("validate_binary")(factor(c(0, 1, 0)), "Y"))
check("infinite ordinal uptake", f("resolve_and_validate_uptake")(c(0, Inf), "ordinal"))
check("formula missingness and row alignment", {
  d <- data.frame(Y=c(0,1,0,1,0), x=c(NA,20,30,40,50), Z=c(0,NA,0,1,0), W=c(0,1,0,1,0))
  r <- f("parse_psbart_formula")(Y ~ x | Z | W, d)
  list(lengths=lapply(r, NROW), x=r$X, y=r$Y, z=r$Z, w=r$W)
})
check("Bayes rule: observed Y=1, both outcome means=1, expected .5", {
  f("compute_posterior_class_prob")(1, .5, .5, 1, 1)
})
check("ordinal grid prediction indexing: expected .1,.2,.7,.8", {
  grid <- cbind(x=c(10,10,20,20),w0=c(1,2,1,2),w1=c(0,0,0,0))
  d <- f("build_strata_posterior_dt")(grid,c(1,1,2,2),c(1,0),c(1,1),
    c(.1,.2,.7,.8),c(.3,.4,.5,.6),c(0,0),c(0,0),1,1)
  as.data.frame(d[, c("id","w0","w1","my0","my1","ppw"), with=FALSE])
})
check("ordinal p_wpair(0,0) vs latent bin (-Inf,1) squared", {
  c(actual=f("p_wpair")(0,0,0,0,1,1), expected=pnorm(1)^2)
})
check("ordinal p_wpair axes with unequal latent means", {
  expected <- (pnorm(3,2,1)-pnorm(2,2,1))*(pnorm(2,0,1)-pnorm(1,0,1))
  c(actual=f("p_wpair")(2,1,0,2,1,1), expected=expected)
})
check("ordinal second scaling changes e", {
  d <- data.frame(Y=c(0,1,0,1),Z=c(0,1,0,1),W=0:3,
    x=as.vector(scale(1:4)),e=qnorm(c(.1,.2,.6,.9)))
  cbind(before=d$e,after=f("prepare_ordinal_data")(d)$X[,"e"])
})
check("ordinal constant known propensity becomes NaN", {
  d <- data.frame(Y=c(0,1,0,1),Z=c(0,1,0,1),W=0:3,x=1:4,e=0)
  f("prepare_ordinal_data")(d)$X
})
check("sample stratum proportions ignore selected units", {
  g <- array(0,c(4,2,2,3),dimnames=list(NULL,NULL,c("nt","at"),NULL))
  g[, , "nt", 2] <- 1; g[, , "at", 3] <- 1
  o <- array(.5,c(4,2,4,3),dimnames=list(NULL,NULL,c("y0","y1","cy0","cy1"),NULL))
  r <- f("get_sample_tau")(g,o,treated=c(FALSE,TRUE,TRUE),include_corr=FALSE)
  r[[1]][1,1,]
})
check("predict one target row", {
  trees <- data.frame(tree=1:2,var=-1,value=c(.2,.3))
  f("predict_one_sample")(trees,matrix(0,1,1))
})
check("overlap merges chains: expected average of Phi(-1),Phi(1)=.5", {
  trees <- data.frame(tree=c(1,1),var=-1,value=c(-1,1),chain=1:2)
  f("predict_one_sample")(trees,matrix(0,2,1))
})
check("raw predict_trees missing sample field", {
  trees <- data.frame(tree=1:2,var=-1,value=c(.2,.3),chain=1,iteration=1,m="y0co")
  f("predict_trees")(trees,matrix(0,2,1))
})
check("survey within-PSU varying weights: expected .9", {
  f("survey_pate")(matrix(c(0,1),2,2),c(TRUE,TRUE),c(1,1),c(1,9),seed=1)
})
check("survey subpopulation counts: equal selected mass expected BB mean near .5", {
  tau <- matrix(c(1,rep(0,9),0),11,5000)
  r <- f("survey_pate")(tau,c(TRUE,rep(FALSE,9),TRUE),c(rep(1,10),2),rep(1,11),seed=1)
  r$estimate
})
check("sensitivity solver available", CVXR::installed_solvers())
check("sensitivity routine ECOS dependency", f("find_shift_weights")(TRUE,1.1,c(.1,.4),c(1,3)))
check("sensitivity gamma=1 identity using installed alternative solver in a local copy", {
  env <- new.env(parent=ns)
  env$solve <- function(a, solver, ...) base::solve(a,solver="CLARABEL",...)
  shifted <- f("find_shift_weights"); environment(shifted) <- env
  env$find_shift_weights <- shifted
  compute <- f("compute_shift_pate"); environment(compute) <- env
  tau <- matrix(c(.1,.4),2,2)
  list(shift_output=shifted(TRUE,1,c(.1,.4),c(1,3)),
    baseline=f("survey_pate")(tau,c(TRUE,TRUE),c(1,1),c(1,3),seed=1)$estimate,
    sensitivity=compute(tau,c(TRUE,TRUE),c(1,1),c(1,3),1,TRUE,seed=1)$estimate)
})
check("pretty labels change boundary inequalities", f("simplify_and_pretty_rules")(c("age <= 40", "age > 40")))
check("duplicate pretty labels merge distinct segment IDs", {
  rules <- f("simplify_and_pretty_rules")(c("x < 0.11", "x < 0.12"))
  list(rules=rules,segments=factor(c(2,3),labels=rules))
})

set.seed(2026)
check("dbarts binary predict already returns probabilities", {
  x <- matrix(rnorm(80),ncol=1)
  m <- dbarts::bart2(x,as.numeric(x[,1]>0),n.trees=5L,n.chains=1L,n.threads=1L,
    n.burn=10L,n.samples=5L,keepTrees=TRUE,verbose=FALSE)
  raw <- predict(m,x,type="bart"); prob <- predict(m,x)
  list(max_error_probability_vs_pnorm_raw=max(abs(prob-pnorm(raw))),
    true_propensity_range=range(prob),double_pnorm_range=range(pnorm(prob)))
})
check("dbarts continuous tree values agree with sampler test predictions", {
  x <- matrix(as.numeric(seq_len(40)),ncol=1); y <- 50+2*x[,1]+rnorm(40)
  m <- dbarts::dbarts(x,y,test=x,control=dbarts::dbartsControl(n.trees=5L,n.chains=1L,
    n.threads=1L,n.burn=0L,n.samples=1L,keepTrees=TRUE,verbose=FALSE))
  s <- m$run(); t <- m$getTrees()
  list(max_abs_difference=max(abs(f("predict_one_sample_raw")(t,x)-s$test)),
    tree_prediction_range=range(f("predict_one_sample_raw")(t,x)),sampler_prediction_range=range(s$test))
})
