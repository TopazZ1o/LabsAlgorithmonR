set.seed(123)

# Формирование матрицы и вектора
make_system <- function(n = 10, p = 2, q = 2) {
  A <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        # Диагональные элементы
        A[i, j] <- 5 * i^(p / 2)
      } else {
        # Внедиагональные: случайные знаки
        sign1 <- sample(c(-1, 1), 1)
        sign2 <- sample(c(-1, 1), 1)
        A[i, j] <- sign1 * 0.01 * (i^(p / 2) + sign2 * j^(q / 2))
      }
    }
  }
  
  # Правая часть
  b <- 4.5 * (1:n)^(p / 2)
  return(list(A = A, b = b))
}

# Итерационный метод Якоби
jacobi <- function(A, b, tol = 1e-8, maxiter = 10000) {
  n <- nrow(A)
  x <- rep(0, n)
  D <- diag(A)
  R <- A - diag(D)
  
  for (k in 1:maxiter) {
    x_new <- (b - R %*% x) / D
    if (sqrt(sum((x_new - x)^2)) < tol) {
      return(list(x = as.numeric(x_new), iter = k, conv = TRUE))
    }
    x <- x_new
  }
  return(list(x = as.numeric(x), iter = maxiter, conv = FALSE))
}

solve_system <- function(n = 10, p = 2, q = 3) {
  cat("==============================================\n")
  cat("n =", n, ", p =", p, ", q =", q, "\n")
  
  sys <- make_system(n, p, q)
  A <- sys$A
  b <- sys$b
  
  # Прямой метод
  t1 <- proc.time()
  x_direct <- solve(A, b)
  time_direct <- (proc.time() - t1)[3]
  res_direct <- norm(A %*% x_direct - b, type = "2")
  
  # Метод Якоби
  t2 <- proc.time()
  jac <- jacobi(A, b, tol = 1e-9, maxiter = 20000)
  time_iter <- (proc.time() - t2)[3]
  res_iter <- norm(A %*% jac$x - b, type = "2")
  
  # Вывод
  cat("Прямой метод (solve):\n")
  cat("  ||Ax - b|| =", format(res_direct, scientific = TRUE),
      "   Время =", round(time_direct, 6), "сек\n\n")
  
  cat("Итерационный метод (Якоби):\n")
  cat("  ||Ax - b|| =", format(res_iter, scientific = TRUE),
      "   Итераций =", jac$iter,
      "   Время =", round(time_iter, 6), "сек\n")
  cat("  Сходимость =", jac$conv, "\n")
}

solve_system(n = 10, p = 2, q = 3)