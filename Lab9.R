# Функция степенного метода
power_method <- function(A, eps = 1e-6, k_max = 1000) {
  n <- nrow(A)
  x <- runif(n)                # Случайный начальный вектор
  x <- x / sqrt(sum(x^2))      # Нормализация
  
  lambda_old <- 0
  for (k in 1:k_max) {
    y <- A %*% x               # Умножаем матрицу на вектор
    lambda_new <- sum(y * x)   # Оценка собственного значения
    x <- y / sqrt(sum(y^2))    # Нормализация нового вектора
    
    # Проверка сходимости
    if (abs(lambda_new - lambda_old) < eps) {
      cat("Сошлось за", k, "итераций\n")
      return(list(lambda = lambda_new, eigenvector = x))
    }
    
    lambda_old <- lambda_new
  }
  
  cat("Не сошлось за", k_max, "итераций\n")
  return(list(lambda = lambda_old, eigenvector = x))
}

n <- 5
a <- 1.4
b <- 0.1

# Задаём матрицу A
A <- matrix(0, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    if (i == j) {
      A[i, j] <- a * i^3
    } else {
      A[i, j] <- b * abs(i - j)^(-1/3)
    }
  }
}

# Запуск метода
result <- power_method(A, eps = 1e-6, k_max = 1000)

cat("Максимальное собственное значение λ =", result$lambda, "\n")
cat("Собственный вектор:\n")
print(result$eigenvector)