library(Matrix)

L <- 1
Tmax <- 1
Nz <- 30
Nt <- 200
hz <- L / Nz
ht <- Tmax / Nt

eps <- 1.0
mu  <- 1.0
lambda2 <- 1.0
g_val <- 0.0

z <- seq(0, L, length.out = Nz + 1)
t <- seq(0, Tmax, length.out = Nt + 1)

f2_t <- function(tt) sin(2 * pi * tt)
f_data <- function(tt) 0.5 * sin(2 * pi * tt)

p_init <- rep(0.5, Nz + 1)

build_tridiag <- function(p_vec) {
  n <- Nz + 1
  main <- numeric(n)
  upper <- numeric(n - 1)
  lower <- numeric(n - 1)
  for (i in 2:Nz) {
    a <- 1 / (mu * hz^2)
    b <- 1 / (mu * hz^2)
    c0 <- a + b + lambda2 / mu + eps / ht^2 + p_vec[i] / ht
    main[i] <- c0
    lower[i - 1] <- -a
    upper[i] <- -b
  }
  main[1] <- 1
  main[n] <- 1
  sparseMatrix(
    i = c(1:n, 2:n, 1:(n - 1)),
    j = c(1:n, 1:(n - 1), 2:n),
    x = c(main, lower, upper),
    symmetric = FALSE
  )
}

solve_direct <- function(p_vec) {
  v <- matrix(0, nrow = Nz + 1, ncol = Nt + 1)
  A <- build_tridiag(p_vec)
  v_prev <- rep(0, Nz + 1)
  v_curr <- rep(0, Nz + 1)
  for (n in 1:Nt) {
    rhs <- numeric(Nz + 1)
    rhs[1] <- 0
    rhs[Nz + 1] <- 0
    for (i in 2:Nz) {
      term_tt <- 2 * v_curr[i] - v_prev[i]
      rhs[i] <- eps * term_tt / ht^2 + g_val
    }
    rhs[2] <- rhs[2] + f2_t(t[n + 1]) / (hz)
    v_next <- solve(A, rhs)
    v_prev <- v_curr
    v_curr <- v_next
    v[, n + 1] <- v_curr
  }
  v
}

build_tridiag_phi <- function(p_vec) {
  n <- Nz + 1
  main <- numeric(n)
  upper <- numeric(n - 1)
  lower <- numeric(n - 1)
  for (i in 2:Nz) {
    a <- 1 / (mu * hz^2)
    b <- 1 / (mu * hz^2)
    c0 <- a + b + lambda2 / mu + eps / ht^2 + p_vec[i] / ht
    main[i] <- c0
    lower[i - 1] <- -a
    upper[i] <- -b
  }
  main[1] <- 1
  main[n] <- 1
  sparseMatrix(
    i = c(1:n, 2:n, 1:(n - 1)),
    j = c(1:n, 1:(n - 1), 2:n),
    x = c(main, lower, upper),
    symmetric = FALSE
  )
}

solve_adjoint <- function(p_vec, v) {
  phi <- matrix(0, nrow = Nz + 1, ncol = Nt + 1)
  A <- build_tridiag_phi(p_vec)
  phi_next <- rep(0, Nz + 1)
  phi_curr <- rep(0, Nz + 1)
  for (n in Nt:1) {
    rhs <- numeric(Nz + 1)
    rhs[Nz + 1] <- 0
    for (i in 2:Nz) {
      term_tt <- 2 * phi_curr[i] - phi_next[i]
      rhs[i] <- eps * term_tt / ht^2
    }
    rhs[1] <- 2 * (v[1, n] - f_data(t[n]))
    phi_prev <- solve(A, rhs)
    phi_next <- phi_curr
    phi_curr <- phi_prev
    phi[, n] <- phi_curr
  }
  phi
}

compute_J <- function(v) {
  sum((v[1, ] - f_data(t))^2) * ht
}

compute_grad <- function(p_vec, v, phi) {
  grad <- numeric(Nz + 1)
  for (i in 1:(Nz + 1)) {
    acc <- 0
    for (n in 2:(Nt + 1)) {
      vt <- (v[i, n] - v[i, n - 1]) / ht
      acc <- acc + vt * phi[i, n] * ht
    }
    grad[i] <- acc
  }
  grad
}

max_iter <- 10
alpha <- 0.01

p <- p_init
for (k in 1:max_iter) {
  v <- solve_direct(p)
  phi <- solve_adjoint(p, v)
  J_val <- compute_J(v)
  grad <- compute_grad(p, v, phi)
  p <- p - alpha * grad
  cat("iter =", k, "J =", J_val, "\n")
}

plot(z, p, type = "l")
