library(princeBART)
ns <- asNamespace('princeBART')
f <- function(name) get(name, ns)
check <- function(label, expr) {
  cat('\n',label,'\n',sep='')
  tryCatch(print(force(expr)),error=function(e) cat('ERROR:',conditionMessage(e),'\n'))
}
check('sensitivity gamma=1 identity; local copy uses available solver', {
  env <- new.env(parent=ns)
  env$solve <- function(a, solver, ...) base::solve(a,solver='CLARABEL',...)
  shifted <- f('find_shift_weights'); environment(shifted) <- env
  env$find_shift_weights <- shifted
  compute <- f('compute_shift_pate'); environment(compute) <- env
  tau <- matrix(c(.1,.4),2,2)
  list(shift_output=shifted(TRUE,1,c(.1,.4),c(1,3)),
    baseline=f('survey_pate')(tau,c(TRUE,TRUE),c(1,1),c(1,3),seed=1)$estimate,
    sensitivity=compute(tau,c(TRUE,TRUE),c(1,1),c(1,3),1,TRUE,seed=1)$estimate)
})
check('continuous tree extraction units vs dbarts test predictions', {
  set.seed(8)
  x <- matrix(as.numeric(seq_len(40)),ncol=1); y <- 50+2*x[,1]+rnorm(40)
  m <- dbarts::dbarts(x,y,test=x,control=dbarts::dbartsControl(n.trees=5L,n.chains=1L,
    n.threads=1L,n.burn=0L,n.samples=1L,keepTrees=TRUE,verbose=FALSE))
  s <- m$run(); t <- m$getTrees()
  list(max_abs_difference=max(abs(f('predict_one_sample_raw')(t,x)-s$test)),
    tree_prediction_range=range(f('predict_one_sample_raw')(t,x)),sampler_prediction_range=range(s$test))
})
check('equal lengths but misaligned formula rows', {
  d <- data.frame(Y=c(NA,1,0,1,0,1),x=c(10,NA,30,40,50,60),
    Z=c(0,1,NA,1,0,1),W=c(0,1,0,NA,0,1))
  r <- f('parse_psbart_formula')(Y ~ x | Z | W,d)
  list(lengths=sapply(r,NROW),X_rows=rownames(r$X),Y_rows=names(r$Y),Z=r$Z,W=r$W)
})
check('binary update all noncompliers: missing co and atnoco updates', {
  audit <- new.env(); audit$calls <- character()
  mk <- function(nm) list(setData=function(d) { audit$calls <- c(audit$calls,paste(nm,length(d@y))); invisible(NULL) })
  samplers <- setNames(lapply(c('co','atnoco','y0nt','y1at','y0co','y1co'),mk),c('co','atnoco','y0nt','y1at','y0co','y1co'))
  f('update_samplers')(samplers,matrix(as.numeric(1:4),ncol=1),c(0,1,0,1),c(0,1,0,1),
    rep(0,4),c(1,0,1,0),c(0,1,0,1),1,0,FALSE)
  audit$calls
})
check('binary perfect compliance initialization', {
  set.seed(9)
  x <- matrix(rnorm(40),20,2)
  f('.fit_psbart_binary')(x,rep(c(0,1),10),rep(c(0,1),10),rep(c(0,1),10),
    n_warmup=0L,n_samples=1L,n_trees=2L)
})
check('ordinal short saved chain tree combination', {
  f('combine_chain_trees')(list(list(trees=NULL)),TRUE)
})
check('ordinal sign and ignored public arguments', {
  p <- array(0,c(4,2,2,3),dimnames=list(iteration=NULL,chain=NULL,variable=c('m_y0','m_y1'),unit=NULL))
  p[, , 1, ] <- .8; p[, , 2, ] <- .2
  imp <- array(0,c(4,2,2,3),dimnames=list(iteration=NULL,chain=NULL,variable=c('w0','w1'),unit=NULL))
  imp[, , 1, ] <- 2; imp[, , 2, ] <- 1
  fit <- structure(list(probs=p,imp=imp,data=list(Z=c(0,1,1))),class=c('prince_bart','prince_bart_ordinal'))
  list(default=coef(fit),unsupported=coef(fit,type='nonsense',treated_only=TRUE))
})
check('sample outcome summaries mix never/always takers', {
  g <- array(0,c(4,2,2,3),dimnames=list(NULL,NULL,c('nt','at'),NULL))
  g[, , 'nt', 2] <- 1; g[, , 'at', 3] <- 1
  o <- array(0,c(4,2,4,3),dimnames=list(NULL,NULL,c('y0','y1','cy0','cy1'),NULL))
  o[,,,3] <- 1
  r <- f('get_sample_tau')(g,o,include_corr=FALSE)
  as.numeric(r[[2]][1,1,])
})
