zscore=function(x) qnorm((rank(x, na.last = "keep") - 0.5)/sum(!is.na(x)))

c.z.hglm <- function(kin){
  relmat <- kin
  relmat[upper.tri(relmat)] <- t(relmat)[upper.tri(relmat)]#上三角变成下三角
  svd <- svd(relmat)
  Z <- svd$u %*% diag(sqrt(svd$d)) #左边的矩阵乘奇异值
  return(Z)
}

get_go_fun=function(gene_name){
  GO_result2 = GO_result[gene_name,"GO.description"]
  GO_list <- unlist(lapply(GO_result2, function(x) {
    matches <- str_extract_all(x, "\\(([^)]+)\\)")
    if (length(matches) > 0) {
      extracted <- unlist(str_split(matches[[1]], ","))
      extracted <- str_replace_all(extracted, "\\(|\\)", "")
      extracted <- extracted[str_detect(extracted, "^GO")]
    } else {
      NULL
    }
  }), use.names = FALSE)
  return(GO_list)
}

mixed <- function(x,y,kk,method="REML",eigen=FALSE){
  loglike <- function(theta) {
    lambda <- exp(theta)
    logdt <- sum(log(lambda * delta + 1))
    h <- 1 / (lambda * delta + 1)
    yy <- sum(yu * h * yu)
    yx <- matrix(0, s, 1)
    xx <- matrix(0, s, s)
    for (i in 1:s) {
      yx[i] <- sum(yu * h * xu[, i])
      for (j in 1:s) {
        xx[i, j] <- sum(xu[, i] * h * xu[, j])
      }
    }
    #xx
    if (method == "REML") {
      loglike <- -0.5 * logdt - 0.5 * (n - s) * log(yy - t(yx) %*% solve(xx) %*% yx) - 0.5 * log(det(xx))
    } else {
      loglike <- -0.5 * logdt - 0.5 * n * log(yy - t(yx) %*% solve(xx) %*% yx)
    }
    return(-loglike)
  }
  
  fixed <- function(lambda) {
    h <- 1 / (lambda * delta + 1)
    yy <- sum(yu * h * yu)
    yx <- matrix(0, s, 1)
    xx <- matrix(0, s, s)
    for (i in 1:s) {
      yx[i] <- sum(yu * h * xu[, i])
      for (j in 1:s) {
        xx[i, j] <- sum(xu[, i] * h * xu[, j])
      }
    }
    beta <- solve(xx, yx)
    if (method == "REML") {
      sigma2 <- (yy - t(yx) %*% solve(xx) %*% yx) / (n - s)
    } else {
      sigma2 <- (yy - t(yx) %*% solve(xx) %*% yx) / n
    }
    var <- diag(solve(xx) * drop(sigma2))
    stderr <- sqrt(var)
    return(c(beta, stderr, sigma2))
  }
  
  n <- length(y)
  qq <- eigen(kk, symmetric = TRUE)
  delta <- qq[[1]]
  if(any(delta< 0))
    delta[delta <0] <- 1e-8
  uu <- qq[[2]]
  s <- ncol(x)
  yu <- t(uu) %*% y
  xu <- t(uu) %*% x
  theta <- 0
  parm <- optim(par = theta, fn = loglike, NULL, hessian = TRUE, method = "Brent", lower = -10, upper = 10)
  lambda <- exp(parm$par)
  conv <- parm$convergence
  fn1 <- parm$value
  fn0 <- loglike(-Inf)
  lrt <- 2 * (fn0 - fn1)
  hess <- parm$hessian
  parmfix <- fixed(lambda)
  beta <- parmfix[1:s]
  stderr <- parmfix[(s + 1):(2 * s)]
  ve <- parmfix[2 * s + 1]
  lod <- lrt / 4.61
  p_value <- 1 - pchisq(lrt, 1)
  va <- lambda * ve
  h2 <- va / (va + ve)
  par <- data.frame(method, beta, stderr, va, ve, lambda, h2, conv, fn1, fn0, lrt, lod, p_value)
  if(eigen){
    return(list(par,qq))
  } else {
    return(list(par))
  }
}


mixed.solve2 <- function (y, Z = NULL, K = NULL, X = NULL, method = "REML", bounds = c(1e-09, 
                                                                                       1e+09), SE = FALSE, return.Hinv = FALSE) 
{
  pi <- 3.14159
  n <- length(y)
  y <- matrix(y, n, 1)
  not.NA <- which(!is.na(y))
  if (is.null(X)) {
    p <- 1
    X <- matrix(rep(1, n), n, 1)
  }
  p <- ncol(X)
  if (is.null(p)) {
    p <- 1
    X <- matrix(X, length(X), 1)
  }
  if (is.null(Z)) {
    Z <- diag(n)
  }
  m <- ncol(Z)
  if (is.null(m)) {
    m <- 1
    Z <- matrix(Z, length(Z), 1)
  }
  stopifnot(nrow(Z) == n)
  stopifnot(nrow(X) == n)
  if (!is.null(K)) {
    stopifnot(nrow(K) == m)
    stopifnot(ncol(K) == m)
  }
  Z <- as.matrix(Z[not.NA, ])
  X <- as.matrix(X[not.NA, ])
  n <- length(not.NA)
  y <- matrix(y[not.NA], n, 1)
  XtX <- crossprod(X, X)
  rank.X <- qr(XtX)$rank
  if(!require(Matrix))
    require(Matrix)
  rank.X2 <- rankMatrix(XtX)
  if (rank.X2 < p) {
    stop("X not full rank")
  }
  XtXinv <- solve(XtX)
  S <- diag(n) - tcrossprod(X %*% XtXinv, X)
  if (n <= m + p) {
    spectral.method <- "eigen"
  }
  else {
    spectral.method <- "cholesky"
    if (!is.null(K)) {
      diag(K) <- diag(K) + 1e-06
      B <- try(chol(K), silent = TRUE)
      if (inherits(B, what = "try-error")) {
        stop("K not positive semi-definite.")
      }
    }
  }
  if (spectral.method == "cholesky") {
    if (is.null(K)) {
      ZBt <- Z
    }
    else {
      ZBt <- tcrossprod(Z, B)
    }
    svd.ZBt <- svd(ZBt, nu = n)
    U <- svd.ZBt$u
    phi <- c(svd.ZBt$d^2, rep(0, n - m))
    SZBt <- S %*% ZBt
    svd.SZBt <- try(svd(SZBt), silent = TRUE)
    if (inherits(svd.SZBt, what = "try-error")) {
      svd.SZBt <- svd(SZBt + matrix(1e-10, nrow = nrow(SZBt), 
                                    ncol = ncol(SZBt)))
    }
    QR <- qr(cbind(X, svd.SZBt$u))
    Q <- qr.Q(QR, complete = TRUE)[, (p + 1):n]
    R <- qr.R(QR)[p + 1:m, p + 1:m]
    ans <- try(solve(t(R^2), svd.SZBt$d^2), silent = TRUE)
    if (inherits(ans, what = "try-error")) {
      spectral.method <- "eigen"
    }
    else {
      theta <- c(ans, rep(0, n - p - m))
    }
  }
  if (spectral.method == "eigen") {
    offset <- sqrt(n)
    if (is.null(K)) {
      Hb <- tcrossprod(Z, Z) + offset * diag(n)
    }
    else {
      Hb <- tcrossprod(Z %*% K, Z) + offset * diag(n)
    }
    Hb.system <- eigen(Hb, symmetric = TRUE)
    phi <- Hb.system$values - offset
    if (min(phi) < -1e-06) {
      stop("K not positive semi-definite.")
    }
    U <- Hb.system$vectors
    SHbS <- S %*% Hb %*% S
    SHbS.system <- eigen(SHbS, symmetric = TRUE)
    theta <- SHbS.system$values[1:(n - p)] - offset
    Q <- SHbS.system$vectors[, 1:(n - p)]
  }
  omega <- crossprod(Q, y)
  omega.sq <- omega^2
  if (method == "ML") {
    f.ML <- function(lambda, n, theta, omega.sq, phi) {
      n * log(sum(omega.sq/(theta + lambda))) + sum(log(phi + 
                                                          lambda))
    }
    soln <- optimize(f.ML, interval = bounds, n, theta, omega.sq, 
                     phi)
    lambda.opt <- soln$minimum
    df <- n
  }
  else {
    f.REML <- function(lambda, n.p, theta, omega.sq) {
      n.p * log(sum(omega.sq/(theta + lambda))) + sum(log(theta + 
                                                            lambda))
    }
    soln <- optimize(f.REML, interval = bounds, n - p, theta, 
                     omega.sq)
    lambda.opt <- soln$minimum
    df <- n - p
  }
  Vu.opt <- sum(omega.sq/(theta + lambda.opt))/df
  Ve.opt <- lambda.opt * Vu.opt
  Hinv <- U %*% (t(U)/(phi + lambda.opt))
  W <- crossprod(X, Hinv %*% X)
  beta <- array(solve(W, crossprod(X, Hinv %*% y)))
  rownames(beta) <- colnames(X)
  if (is.null(K)) {
    KZt <- t(Z)
  }
  else {
    KZt <- tcrossprod(K, Z)
  }
  KZt.Hinv <- KZt %*% Hinv
  u <- array(KZt.Hinv %*% (y - X %*% beta))
  if (is.null(K)) {
    rownames(u) <- colnames(Z)
  }
  else {
    rownames(u) <- rownames(K)
  }
  LL = -0.5 * (soln$objective + df + df * log(2 * pi/df))
  if (!SE) {
    if (return.Hinv) {
      return(list(Vu = Vu.opt, Ve = Ve.opt, beta = beta, 
                  u = u, LL = LL, Hinv = Hinv))
    }
    else {
      return(list(Vu = Vu.opt, Ve = Ve.opt, beta = beta, 
                  u = u, LL = LL))
    }
  }
  else {
    Winv <- solve(W)
    beta.SE <- array(sqrt(Vu.opt * diag(Winv)))
    rownames(beta.SE) <- rownames(beta)
    WW <- tcrossprod(KZt.Hinv, KZt)
    WWW <- KZt.Hinv %*% X
    if (is.null(K)) {
      u.SE <- array(sqrt(Vu.opt * (rep(1, m) - diag(WW) + 
                                     diag(tcrossprod(WWW %*% Winv, WWW)))))
    }
    else {
      u.SE <- array(sqrt(Vu.opt * (diag(K) - diag(WW) + 
                                     diag(tcrossprod(WWW %*% Winv, WWW)))))
    }
    rownames(u.SE) <- rownames(u)
    if (return.Hinv) {
      return(list(Vu = Vu.opt, Ve = Ve.opt, beta = beta, 
                  beta.SE = beta.SE, u = u, u.SE = u.SE, LL = LL, 
                  Hinv = Hinv))
    }
    else {
      return(list(Vu = Vu.opt, Ve = Ve.opt, beta = beta, 
                  beta.SE = beta.SE, u = u, u.SE = u.SE, LL = LL))
    }
  }
}