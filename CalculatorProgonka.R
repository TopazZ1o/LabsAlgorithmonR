# Калькулятор решения СЛАУ методом прогонки
solve_tridiag <- function(a, b, c, d, outfile = "Answer.txt") {
  n <- length(b)
  
  # Проверка размеров
  if (length(a) != n || length(c) != n || length(d) != n) {
    stop("Все векторы должны иметь одинаковую длину")
  }
  
  # Прегоночные коэффициенты
  A <- numeric(n)
  B <- numeric(n)
  x <- numeric(n)
  
  # Прямой ход
  A[1] <- c[1] / b[1]
  B[1] <- d[1] / b[1]
  
  for (i in 2:(n - 1)) {
    denom <- b[i] - a[i - 1] * A[i - 1]
    A[i] <- c[i] / denom
    B[i] <- (d[i] - a[i - 1] * B[i - 1]) / denom
  }
  
  # Последний элемент
  x[n] <- B[n]
  
  # Обратный ход
  for (i in (n - 1):1) {
    x[i] <- B[i] - A[i] * x[i + 1]
  }
  
  content <- paste("Решение СЛАУ методом прогонки:\n",
                   paste("x[", 1:n, "] = ", x, collapse = "\n", sep = ""),
                   sep = "")
  
  writeLines(content, outfile)
  
  message("Решение сохранено в файл: ", outfile)
  
  return(x)
}

n <- as.integer(readline("Введите размерность n: "))

cat("Введите нижнюю диагональ a (", n-1, " значений):\n")
a <- scan(what = numeric(), n = n - 1)
a <- c(a, 0)   # дополняем до длины n

cat("Введите главную диагональ b (", n, " значений):\n")
b <- scan(what = numeric(), n = n)

cat("Введите верхнюю диагональ c (", n-1, " значений):\n")
c <- scan(what = numeric(), n = n - 1)
c <- c(0, c)   # смещаем, чтобы длина была n

cat("Введите правую часть d (", n, " значений):\n")
d <- scan(what = numeric(), n = n)

x <- solve_tridiag(a, b, c, d, outfile = "Answer.txt")

cat("Решение:\n")
print(x)

