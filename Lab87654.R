phi <- function(x) {
  1.4 * x * (1 - x)
}

a <- -1
b <- 1
eps <- 1e-3
c <- (3 - sqrt(5)) / 2  # коэффициент золотого сечения

# Итерационный процесс
k <- 0
while (abs(b - a) > eps) {
  x1 <- a + c * (b - a)
  x2 <- a + (1 - c) * (b - a)
  
  if (phi(x1) > phi(x2)) {
    a <- x1
  } else {
    b <- x2
  }
  
  k <- k + 1
}

xmin <- (a + b) / 2
phi_min <- phi(xmin)

cat("Результаты метода золотого сечения:\n")
cat("Минимум функции в точке x =", round(xmin, 6), "\n")
cat("Значение функции Φ(x) =", round(phi_min, 6), "\n")
cat("Количество итераций =", k, "\n")
