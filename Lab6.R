a <- 0.1
b <- 1.0

f <- function(x) {
  arg <- x + sqrt(x * (1 + x))
  arg <- pmin(pmax(arg, -1), 1)  # ограничение в [-1, 1]
  sqrt(x) * asin(arg)
}

# Метод прямоугольников (средних)
rectangles_midpoint <- function(f, a, b, m) {
  h <- (b - a) / m
  x_mid <- a + ((0:(m-1)) + 0.5) * h
  y_mid <- f(x_mid)
  I <- h * sum(y_mid)
  return(list(I = I))
}

# Квадратурная формула Гаусса–Лежандра
gauss_legendre <- function(f, a, b, n) {
  suppressMessages({
    if (!require(statmod, quietly = TRUE)) {
      install.packages("statmod", repos = "https://cloud.r-project.org", quiet = TRUE)
      library(statmod, quietly = TRUE)
    }
  })
  res <- gauss.quad(n, kind = "legendre")
  t <- res$nodes
  A <- res$weights
  c <- (b - a) / 2.0
  d <- (a + b) / 2.0
  vals <- f(c * t + d)
  I <- c * sum(A * vals)
  return(list(I = I))
}

# "Точное" значение интеграла
I_true <- integrate(f, a, b)$value

cat(sprintf("Численно (integrate): I ≈ %.10f\n\n", I_true))

# Метод прямоугольников
cat("Метод прямоугольников (средних):\n")
for (m in c(10, 20, 40, 80, 100)) {
  I_rect <- rectangles_midpoint(f, a, b, m)$I
  err <- abs(I_rect - I_true)
  cat(sprintf("  m=%3d  I≈%.10f   |ошибка|=%.3e\n", m, I_rect, err))
}
cat("\n")

# Квадратурная формула Гаусса–Лежандра
cat("Квадратурная формула Гаусса–Лежандра:\n")
for (n in 5:11) {
  I_gauss <- gauss_legendre(f, a, b, n)$I
  err <- abs(I_gauss - I_true)
  cat(sprintf("  c_m=%2d  I≈%.10f   |ошибка|=%.3e\n", n, I_gauss, err))
}
cat("\n(Замечание: t — стандартные узлы на [-1,1], A — соответствующие веса.)\n")