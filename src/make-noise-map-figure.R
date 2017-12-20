#!/usr/bin/Rscript

EndemicEquilSIR <- function(beta=(R0 * (mu + gamma)), eta=17/5e4,
                            gamma=365/22, mu=1/50, p=0,  R0=17, verbose=FALSE) {
  ## Computes the endemic equilibrium of an SIR model with immigration, Eq. 6.
  ##
  ## Args:
  ##   beta: numeric. The transmission rate.
  ##   eta: numeric. The rate of infection from outside.
  ##   gamma: numeric. The recovery rate.
  ##   mu: numeric. The birth rate.
  ##   p: numeric. The vaccination uptake.
  ##   R0: numeric. The basic reproduction number.
  ##
  ## Returns:
  ##   A list with numeric elements S, I, and R, coresponding to the
  ##   equilibrium fractions of the population in the
  ##   susceptible, infected, and removed states.
  stopifnot(c(beta, eta, gamma, p, R0) >= 0)
  stopifnot(p <= 1)
  a <- - beta * (gamma + mu)
  b <- beta * mu * (1 - p) - (gamma + mu) * (eta + mu)
  c <- mu * (1 - p) * eta
  eq <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
  i.star <- ifelse(p == 1, 0, eq)
  s.star <- ifelse(p == 1, 0, mu * (1 - p)/ (beta * i.star + eta + mu))
  if (verbose) {
    ds.star <- mu *(1 - p) - beta * s.star * i.star - eta * s.star - mu * s.star
    di.star <- beta * s.star * i.star + eta * s.star - (gamma + mu) * i.star
    cat('dS = ', ds.star, '\n')
    cat('dI = ', di.star, '\n')
  }
  return(list(S=s.star, I=i.star, R=1 - i.star - s.star))
}

DiffusionMatrixSIR <- function(beta=R0 * (gamma + mu), eta, gamma, p, R0,
                               eval.pars, mu){
  S <- eval.pars['S']
  I <- eval.pars['I']
  ## D from Table 9
  rate.trans <- (beta * I + eta ) * S
  ret <- rbind(c(rate.trans + mu * (S + (1 - p)), -rate.trans),
               c(-rate.trans, rate.trans + (gamma + mu) * I))
  dimnames(ret) <- list(c('S', 'I'), c('S', 'I'))
  ret
}

DriftMatrixSIR <- function(beta=R0 * (gamma + mu), eta, gamma, p, R0,
                           eval.pars, mu){
  S <- eval.pars['S']
  I <- eval.pars['I']
  ## a from Table 9
  ret <- rbind(c(-beta * I - eta - mu, -beta * S),
               c(beta * I + eta, beta * S - (gamma + mu)))
  dimnames(ret) <- list(c('S', 'I'), c('S', 'I'))
  ret
}


MakeSIR <- function(def.pars=list(beta=282.3855, eta=(1 / 5e3), gamma=(365 /22),
                        p=0, R0=17, mu=(1 / 50), p.s=0, p.0=0, imports=0,
                        pop.size=1e5, init.vars=c('X'=1, 'Y'=1),
                        init.vars.determ=c(S=0.1, I=0.1),
                        init.vars.ensemble=c(S=0.1, I=0.1, MSS=0, MII=0, MSI=0),
                        eval.pars=c(S=0.1, I=0.1), n.intervals=1e3,
                        start.time=0, stop.time=100,
                        time.steps=seq(0, 500, len=5001), verbose=FALSE)) {
  ## Creates an object representing the SIR model given by Eq. 7 with default
  ## parameters supplied by the arguments.
  ##
  ## Args:
  ##   def.pars: a list with default values for all parameters for all methods
  ##             you will use on the object
  ##
  ## Returns:
  ##   An object of class 'sir', a list containing the list def.pars

  ret <- list(def.pars=def.pars)
  attr(ret, 'class') <- 'sir'
  return(ret)
}

GetPars <- function(def, over, need){
  stopifnot(names(need) %in% names(def))
  pars <- def
  if (length(over) > 0){
    stopifnot(names(over) %in% names(def))
    pars[names(over)] <- over
  }
  pars[names(need)]
}

EvalPars.sis <- function(x, ...){
  over <- list(...)
  stopifnot(names(over) %in% names(x$def.pars))
  x$def.pars[names(over)] <- over
  x$def.pars$eval.pars
}


EvalPars.sir <- function(x, ...){
  EvalPars.sis(x, ...)
}

AnalyticVar.sir <- function(x, vars=c('S', 'I'), ...){
  ## check that we are not at trivial equilibrium
  EvalI <- EvalPars(x, ...)['I']
  if (is.na(EvalI)) browser()
  if (EvalI <= 0){
    ret <- NA
  } else {
    KroneckerSum <- function(M, N){
      Im <- diag(nrow=nrow(M))
      In <- diag(nrow=nrow(N))
      kronecker(M, In) + kronecker(Im, N)
    }
    A <- DiffusionMatrix(x, ...)
    B <- DriftMatrix(x, ...)
    nvar <- nrow(A)
    foo <- try(solve(KroneckerSum(-B, -B)) %*% as.vector(A))
    if (inherits(foo, "try-error")){
      ret <- NA
    } else {
      cov <- matrix(foo, nvar, nvar, dimnames=dimnames(A))
      ret <- cov[vars, vars]
    }
  }
  ret
}

DiffusionMatrix.sir <- function(x, ...){
  def <- x$def.pars
  over <- list(...)
  need <- formals(DiffusionMatrixSIR)
  pars <- GetPars(def, over, need)
  do.call(DiffusionMatrixSIR, pars)
}

DriftMatrix.sir <- function(x, ...){
  def <- x$def.pars
  over <- list(...)
  need <- formals(DriftMatrixSIR)
  pars <- GetPars(def, over, need)
  do.call(DriftMatrixSIR, pars)
}

EndemicEquil.sir <- function(x, ...){
  def <- x$def.pars
  over <- list(...)
  need <- formals(EndemicEquilSIR)
  pars <- GetPars(def, over, need)
  do.call(EndemicEquilSIR, pars)
}

EvalPars <- function(x, ...) UseMethod("EvalPars")
EndemicEquil <- function(x, ...) UseMethod("EndemicEquil")
DiffusionMatrix <- function(x, ...) UseMethod("DiffusionMatrix")
DriftMatrix <- function(x, ...) UseMethod("DriftMatrix")
AnalyticVar <- function(x, ...) UseMethod("AnalyticVar")

## pre-flight checks

m <- MakeSIR()
m$def.pars$p <- 0
m$def.pars$beta <- 42
eval.pars <- unlist(EndemicEquil(m))

Sigma <- AnalyticVar(m, eval.pars=eval.pars)
D <- DiffusionMatrix(m, eval.pars=eval.pars)
F <- DriftMatrix(m, eval.pars=eval.pars)

T <- eigen(F)$vectors
Lambdas <- eigen(F)$values
Dt <- solve(T) %*% D %*% solve(t(T))

Qt11 <- 0
Qt22 <- 0
Qt12 <- (Lambdas[1] - Lambdas[2]) / sum(Lambdas) * Dt[1,2]
Qt21 <- (Lambdas[2] - Lambdas[1]) / sum(Lambdas) * Dt[2,1]
Qt <- matrix(c(Qt11, Qt21, Qt12, Qt22), ncol=2)
Sigmat <- solve(T) %*% Sigma %*% solve(t(T))

get_sigmat_alt <- function(D, T, eigs){
    n <- ncol(D)
    Dt <- solve(T) %*% D %*% solve(t(T))
    ret <- matrix(NA, ncol=n, nrow=n)
    for (i in 1:n){
        for (j in 1:n){
            ret[i, j] <- -Dt[i, j] * (1 + (eigs[i] - eigs[j]) / (eigs[i] + eigs[j])) / (2 * eigs[i])
        }
    }
    ret
}
all.equal(Sigmat, get_sigmat_alt(D, T, Lambdas))

check_abcd_formula <- function(Sigma, F){
    a <- F[1,1]
    b <- F[1,2]
    c <- F[2,1]
    d <- F[2,2]
    l <- eigen(F)$values
    e1 <- c(b, l[1] - a)
    e2 <- c(b, l[2] - a)
    T <- cbind(e1, e2)
    St <- solve(T) %*% Sigma %*% solve(t(T))
    ret <- matrix(nrow=2, ncol=2)
    ret[1,1] <- b*b*sum(St)
    ret[2,1] <- b * ((l[1] - a) * sum(St[1, ]) + (l[2] - a) * sum(St[2,]))
    ret[1,2] <- ret[2,1]
    ret[2,2] <- (l[1] - a) * ((l[1] - a) * St[1, 1] + (l[2] - a) * St[1,2]) + (l[2] - a) * ((l[1] - a) * St[2, 1] + (l[2] - a) * St[2,2])
    ret
}
all.equal(Sigma, check_abcd_formula(Sigma, F), check.attributes=FALSE)

check_autocovar_abcd_formula <- function(Sigma, F, tau){
    a <- F[1,1]
    b <- F[1,2]
    c <- F[2,1]
    d <- F[2,2]
    stopifnot(b != 0)
    l <- eigen(F)$values
    e1 <- c(b, l[1] - a)
    e2 <- c(b, l[2] - a)
    T <- cbind(e1, e2)
    St <- solve(T) %*% Sigma %*% solve(t(T))
    ret <- matrix(nrow=2, ncol=2)
    ret[1,1] <- b * b * sum(exp(l * tau) * St)
    ret[1,2] <- b * ((l[1] - a) * (exp(l[1] * tau) * St[1, 1] + exp(l[2] * tau) * St[1, 2])
                   + (l[2] - a) * (exp(l[1] * tau) * St[2, 1] + exp(l[2] * tau) * St[2, 2]))
    ret[2,1] <- b * ((l[1] - a) * (exp(l[1] * tau) * St[1, 1] + exp(l[1] * tau) * St[1, 2])
                   + (l[2] - a) * (exp(l[2] * tau) * St[2, 1] + exp(l[2] * tau) * St[2, 2]))
    ret[2,2] <- ((l[1] - a) * exp(l[1] * tau) * ((l[1] - a) * St[1, 1] + (l[2] - a) * St[1,2])
                 + (l[2] - a) * exp(l[2] * tau) * ((l[1] - a) * St[2, 1] + (l[2] - a) * St[2,2]))
    ret
}

all.equal(expm::expm(F * 1.0) %*% Sigma,
          check_autocovar_abcd_formula(Sigma, F, tau=1.0), check.attributes=FALSE)

all.equal(expm::expm(F * 1.5) %*% Sigma,
          check_autocovar_abcd_formula(Sigma, F, tau=1.5), check.attributes=FALSE)

Sigma2 <- AnalyticVar(m, eval.pars=eval.pars, mu=.2)
F2 <- DriftMatrix(m, eval.pars=eval.pars, mu=.2)

all.equal(expm::expm(F2 * 3.5) %*% Sigma2,
          check_autocovar_abcd_formula(Sigma2, F2, tau=3.5), check.attributes=FALSE)

check_w_formula <- function(Sigma, F){
    a <- F[1,1]
    b <- F[1,2]
    c <- F[2,1]
    d <- F[2,2]
    l <- eigen(F)$values
    e1 <- c(b, l[1] - a)
    e1 <- e1 / sqrt(sum(e1^2))
    e2 <- c(b, l[2] - a)
    e2 <- e2 / sqrt(sum(e2^2))
    W <- cbind(e1, e2)
    St <- solve(W) %*% Sigma %*% solve(t(W))
    ret <- matrix(nrow=2, ncol=2)
    ret[1,1] <- W[1,1]^2 * St[1, 1] + 2 * W[1, 1] * W[1, 2] * St[1, 2] + W[1, 2]^2 * St[2, 2]
    ret[2,1] <- W[1, 1] * W[2, 1] * St[1, 1] + (W[1, 1] * W[2, 2] + W[1, 2] * W[2, 1]) * St[1, 2] + W[1, 2] * W[2, 2] * St[2, 2]
    ret[1,2] <- ret[2,1]
    ret[2,2] <- W[2,1]^2 * St[1, 1] + 2 * W[2, 1] * W[2, 2] * St[1, 2] + W[2, 2]^2 * St[2, 2]
    ret
}
all.equal(Sigma, check_w_formula(Sigma, F), check.attributes=FALSE)

check_autocovar_w_formula <- function(Sigma, F, tau){
    a <- F[1,1]
    b <- F[1,2]
    c <- F[2,1]
    d <- F[2,2]
    l <- eigen(F)$values
    e1 <- c(b, l[1] - a)
    e1 <- e1 / sqrt(sum(e1^2))
    e2 <- c(b, l[2] - a)
    e2 <- e2 / sqrt(sum(e2^2))
    W <- cbind(e1, e2)
    St <- solve(W) %*% Sigma %*% solve(t(W))
    ret <- matrix(nrow=2, ncol=2)
    ret[1,1] <- exp(l[1] * tau) * W[1,1]^2 * St[1, 1] + (exp(l[1] * tau) + exp(l[2] * tau)) * W[1, 1] * W[1, 2] * St[1, 2] + exp(l[2] * tau) * W[1, 2]^2 * St[2, 2]
    ret[1,2] <- exp(l[1] * tau) * W[1, 1] * W[2, 1] * St[1, 1] + (exp(l[1] * tau) * W[1, 1] * W[2, 2] + exp(l[2] * tau) * W[1, 2] * W[2, 1]) * St[1, 2] + exp(l[2] * tau) * W[1, 2] * W[2, 2] * St[2, 2]
    ret[2,1] <- exp(l[1] * tau) * W[1, 1] * W[2, 1] * St[1, 1] + (exp(l[2] * tau) * W[1, 1] * W[2, 2] + exp(l[1] * tau) * W[1, 2] * W[2, 1]) * St[1, 2] + exp(l[2] * tau) * W[1, 2] * W[2, 2] * St[2, 2]
    ret[2,2] <- exp(l[1] * tau) * W[2,1]^2 * St[1, 1] + (exp(l[1] * tau) + exp(l[2] * tau)) * W[2, 1] * W[2, 2] * St[1, 2] + exp(l[2] * tau) * W[2, 2]^2 * St[2, 2]
    ret
}
all.equal(expm::expm(F * 1.5) %*% Sigma,
          check_autocovar_w_formula(Sigma, F, tau=1.5), check.attributes=FALSE)

## define functions for surface plots

get_boundary_coefs2 <- function(F, tau=1, eps=.1, focal_var="first"){
    a <- F[1, 1]
    b <- F[1, 2]
    c <- F[2, 1]
    d <- F[2, 2]
    l <- eigen(F)$values
    e1 <- c(b, l[1] - a)
    e1 <- e1 / sqrt(sum(e1^2))
    e2 <- c(b, l[2] - a)
    e2 <- e2 / sqrt(sum(e2^2))
    W <- unname(cbind(e1, e2))

    if (focal_var == "first") {
      n1 <- eps * c(W[1,1]^2 * (-1) / (2 * l[1]),
                   2 * W[1,1] * W[1,2] * (-1) * (1 + (l[1]  - l[2]) / sum(l)) / (2 * l[1]),
                   W[1,2]^2 * (-1) / (2 * l[2]))
      n2 <- c(0,
              (exp(l[2] * tau) - exp(l[1] * tau)) * W[1,1] * W[1,2] * (-1) * (1 + (l[1]  - l[2]) / sum(l)) / (2 * l[1]),
              (exp(l[2] * tau) - exp(l[1] * tau)) * W[1,2]^2 * (-1) / (2 * l[2]))
    } else {
      n1 <- eps * c(W[2,1]^2 * (-1) / (2 * l[1]),
                   2 * W[2,1] * W[2,2] * (-1) * (1 + (l[1]  - l[2]) / sum(l)) / (2 * l[1]),
                   W[2,2]^2 * (-1) / (2 * l[2]))
      n2 <- c(0,
             (exp(l[2] * tau) - exp(l[1] * tau)) * W[2,1] * W[2,2] * (-1) * (1 + (l[1]  - l[2]) / sum(l)) / (2 * l[1]),
             (exp(l[2] * tau) - exp(l[1] * tau)) * W[2,2]^2 * (-1) / (2 * l[2]))
    }

    A <- matrix(nrow=3,ncol=3)
    A[1,1] <- W[1,1]^2
    A[1,2] <- W[1,1] * W[1,2] * 2
    A[1,3] <- W[1,2]^2
    A[2,1] <- W[2,1] * W[1,1]
    A[2,2] <- W[2,2] * W[1,1] + W[2,1] * W[1,2]
    A[2,3] <- W[2,2] * W[1,2]
    A[3,1] <- W[2,1]^2
    A[3,2] <- W[2,1]*W[2,2]*2
    A[3,3] <- W[2,2]^2

    alphat <- n2[2] - n1[2]
    betat <- n2[3] - n1[3]
    gammat <- n2[1] - n1[1]

    nt <- c(c=gammat, a=alphat, b=betat)
    n <- drop(t(solve(A)) %*% nt)
    names(n) <- names(nt)
    n <- n / sqrt(sum(n^2))
    nt <- nt / sqrt(sum(nt^2))

    list(n1=n1, n2=n2, nt=nt, n=n, A=A)
}

calc_error_FD <- function(F, D, tau=1){
    a <- F[1,1]
    b <- F[1,2]
    c <- F[2,1]
    d <- F[2,2]
    l <- eigen(F)$values
    e1 <- c(b, l[1] - a)
    e1 <- e1 / sqrt(sum(e1^2))
    e2 <- c(b, l[2] - a)
    e2 <- e2 / sqrt(sum(e2^2))
    W <- cbind(e1, e2)
    St <- get_sigmat_alt(D=D, T=W, eigs=l)
    Sigma <- W %*% St %*% t(W)
    acov <- expm::expm(F * tau) %*% Sigma
    l <- eigen(F)$values
    acor <- c(acov[1,1] / Sigma[1,1], acov[2,2] / Sigma[2,2])
    list(acor - exp(l[1] * tau), acor - exp(l[2] * tau), W=W)
}

repnum_to_mats <- function(repnum, m){
  m$def.pars$p <- 0
  m$def.pars$beta <- repnum * (m$def.pars$mu + m$def.pars$gamma)
  eval.pars <- unlist(EndemicEquil(m))
  ret <- list()
  ret$F <- DriftMatrix(m, eval.pars=eval.pars)
  ret$D <- DiffusionMatrix(m, eval.pars=eval.pars)
  ret$pars <- c(m$def.pars[c("beta", "mu", "gamma", "eta", "mu")], eval.pars)
  ret
}

## make surface plot data

repnum_levels <- c(0.1, 0.5, 0.9)
m <- MakeSIR()

mats <- lapply(repnum_levels, repnum_to_mats, m=m)
Flist <- lapply(mats, "[[", "F")
Dlist <- lapply(mats, "[[", "D")

cat("Parameters for panels:\n")
print(sapply(mats, "[[", "pars"))

stopifnot(!any(sapply(Flist, function(x) is.complex(eigen(x)$values))))

## plot polygon containing surfaces

make_prism <- function(eps=c(-.01, 0.01), d11points=c(0, 0.1),
                       d22points=c(0, 0.01), F, D, focal_var="first",
                       topmar_text=NULL, leftmar_text=NULL){

  l <- eigen(F)$values
  tau <- (log(abs(l[2])) - log(abs(l[1]))) / (l[1] - l[2])
  get_d12_point <- function(F, epsilon, d11, d22){
    foo <- get_boundary_coefs2(F=F, tau=tau, eps=epsilon, focal_var=focal_var)
    (-foo$n["b"] * d22 - foo$n["c"] * d11)/ foo$n["a"]
  }

  get_angle <- function(eps1, eps2, F){
    foo1 <- get_boundary_coefs2(F=F, tau=tau, eps=eps1)
    foo2 <- get_boundary_coefs2(F=F, tau=tau, eps=eps2)
    acos(sum(foo1$n * foo2$n))
  }

  theta <- get_angle(eps[1], eps[2], F)
  stopifnot(theta < pi/2)

  vgrid <- expand.grid(eps=eps, d11=d11points, d22=d22points)
  vgrid$d12 <- with(vgrid, mapply(get_d12_point, d11=d11, d22=d22, epsilon=eps,
                                  MoreArgs=list(F=F)))

  reord <- function(x, cols){
      ord <- order(x[, cols[1]])
      x <- x[ord, ]
      if (x[2, cols[2]] != x[3, cols[2]]){
          x <- x[c(1,2,4,3), ]
      }
      x
  }

  topbot <- split(vgrid, vgrid$eps)
  topbot <- lapply(topbot, reord, cols=c(2,3))

  xl <- lapply(topbot, "[[", "d11")
  yl <- lapply(topbot, "[[", "d22")
  zl <- lapply(topbot, "[[", "d12")

  ew <- split(vgrid, vgrid$d22)
  ew <- lapply(ew, reord, cols=c(2,4))

  xl <- c(xl, lapply(ew, "[[", "d11"))
  yl <- c(yl, lapply(ew, "[[", "d22"))
  zl <- c(zl, lapply(ew, "[[", "d12"))

  ns <- split(vgrid, vgrid$d11)
  ns <- lapply(ns, reord, cols=c(3,4))

  xl <- c(xl, lapply(ns, "[[", "d11"))
  yl <- c(yl, lapply(ns, "[[", "d22"))
  zl <- c(zl, lapply(ns, "[[", "d12"))

  nabind <- function(x, y) c(x, NA, y)
  x <- Reduce(nabind, xl)
  y <- Reduce(nabind, yl)
  z <- Reduce(nabind, zl)

  plot3D::polygon3D(x, y, z, border="black", lwd=2, alpha=0.25, bty="b",
                    xlab="$d_{11}$", ylab="$d_{22}$", zlab="$d_{12}$",
                    ticktype="detailed", nticks=3)
  plot3D::scatter3D(z=D[1,2], x=D[1,1], y=D[2,2], add=TRUE)
  mtext(topmar_text, side=3)
  mtext(leftmar_text, side=2, line=2)

  z1 <- get_d12_point(F=F, eps[1], d11=D[1,1], d22=D[2,2])
  z1b <- get_d12_point(F=F, eps[1], d11=D[1,1], d22=D[2,2])
  if (abs (D[1,2] - z1b) < abs(D[1,2] - z1)){
    z1 <- z1b
  }
  plot3D::segments3D(x0=D[1,1], y0=D[2,2], z0=D[1,2], z1=z1, col="red", add=TRUE)
}

#pdf("noise-map.pdf", width=6.5, height=4.5)

tikzDevice::tikz("noise-map.tex", width = 6.5, height = 4.5, standAlone = TRUE)

par(mfrow=c(2, 3), mar=c(0.1 ,5.1, 2.1,0))
eps <- c(-0.1, 0.1)
stopifnot(repnum_levels == c(0.1, 0.5, 0.9))
make_prism(F=Flist[[1]], D=Dlist[[1]], focal_var="first", eps=eps,
           topmar_text="$R_0 = 0.1$",
           leftmar_text="Low-error region for\nsusceptibles")
make_prism(F=Flist[[2]], D=Dlist[[2]], focal_var="first", eps=eps,
           topmar_text="$R_0 = 0.5$")
make_prism(F=Flist[[3]], D=Dlist[[3]], focal_var="first", eps=eps,
           topmar_text="$R_0 = 0.9$")
make_prism(F=Flist[[1]], D=Dlist[[1]], focal_var="second", eps=eps,
           leftmar_text="Low-error region for\ninfecteds")
make_prism(F=Flist[[2]], D=Dlist[[2]], focal_var="second", eps=eps)
make_prism(F=Flist[[3]], D=Dlist[[3]], focal_var="second", eps=eps)
dev.off()

#plot3Drgl::plotrgl(new=FALSE)

q("no")

## check that choice of lag gets maximum error

get_Dmat <- function(F, epsilon, d11, d22, focal_var){
  l <- eigen(F)$values
  tau <- (log(abs(l[2])) - log(abs(l[1]))) / (l[1] - l[2])
  foo <- get_boundary_coefs2(F=F, tau=tau, eps=epsilon, focal_var=focal_var)
  d12 <- (-foo$n["b"] * d22 - foo$n["c"] * d11)/ foo$n["a"]
  list(D=matrix(c(d11, d12, d12, d22), ncol=2), taumax=tau)
}

foo <- get_Dmat(Flist[[1]], eps=.1, d11=.004, d22=.02, focal_var="first")
tmpf <- function(x) calc_error_FD(Flist[[1]], foo$D, tau=x)[[1]][1]

taus <- seq(0, 4, by=0.01)
errs <- sapply(taus, tmpf)

plot(taus, errs, type='l')
abline(v=foo$taumax)
abline(h=0.1)

plot(taus, errs, type='l', xlim=foo$taumax*c(.5,2), ylim=c(.09, .11))
abline(v=foo$taumax)
abline(h=0.1)

foo <- get_Dmat(Flist[[2]], eps=.01, d11=.01, d22=.03, focal_var="second")
tmpf <- function(x) calc_error_FD(Flist[[2]], foo$D, tau=x)[[1]][2]

taus <- seq(0, 4, by=0.01)
errs <- sapply(taus, tmpf)

plot(taus, errs, type='l')
abline(v=foo$taumax)
abline(h=0.1)

plot(taus, errs, type='l', xlim=foo$taumax*c(.5,2), ylim=c(.009, .011))
abline(v=foo$taumax)
abline(h=0.01)
