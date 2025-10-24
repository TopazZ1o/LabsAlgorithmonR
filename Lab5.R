a <- 0
b <- 1
n <- 50 # число делений (узлов n+1)
h <- (b - a) / n
x <- seq(a, b, by = h)

y_fun <- function(x) { cosh(x^2) }
dy_exact <- function(x) { 2*x*sinh(x^2) } # первая производная y' = x
d2y_exact <- function(x) { 2*sinh(x^2) + 4*x^2*cosh(x^2) } # вторая производная y'' = 1

# значения функции в узлах
y <- y_fun(x)

# 1. Численные производные
#    - Для внутренних узлов используем центральную разность (2-ой порядок)
#    - На границах — односторонние разности (порядок 1 или 2)

# Первая производная: инициализируем вектор
dy_num <- numeric(length(x))

# центральные разности для внутренних узлов j=2..n
for (j in 2:(length(x)-1)) {
  dy_num[j] <- (y[j+1] - y[j-1]) / (2*h)   # центральная разность
}
# Односторонние на концах
dy_num[1] <- (y[2] - y[1]) / h             # прямой вперёд (1-й порядок)
dy_num[length(x)] <- (y[length(x)] - y[length(x)-1]) / h  # назад

# Вторая производная: центральная формула
d2y_num <- numeric(length(x))
for (j in 2:(length(x)-1)) {
  d2y_num[j] <- (y[j+1] - 2*y[j] + y[j-1]) / (h^2)
}
# Границы: используем "экстраполяцию" простым подходом, копируем ближайшие
d2y_num[1] <- d2y_num[2]
d2y_num[length(x)] <- d2y_num[length(x)-1]

# 2. Ошибки относительно аналитического решения
dy_exact_vals <- dy_exact(x)
d2y_exact_vals <- d2y_exact(x)

abs_err_dy <- abs(dy_exact_vals - dy_num)
abs_err_d2y <- abs(d2y_exact_vals - d2y_num)

cat(sprintf("Первая производная: max abs err = %.3e ; RMSE = %.3e\n",
            max(abs_err_dy), sqrt(mean((dy_exact_vals - dy_num)^2))))
cat(sprintf("Вторая производная: max abs err = %.3e ; RMSE = %.3e\n",
            max(abs_err_d2y), sqrt(mean((d2y_exact_vals - d2y_num)^2))))

# 3. Графики
# 3.1 График исходной функции y = cosh(x^2)
plot(x, y, type = "l", lwd = 2, col = "blue",
     main = "Исходная функция y = cosh(x^2)",
     xlab = "x", ylab = "y")

# 3.2 первая производная: точная vs приближённая
plot(x, dy_exact_vals, type = "l", lwd = 2, main = "Первая производная: точная vs прибл.",
     xlab = "x", ylab = "y'")
lines(x, dy_num, lty = 2, lwd = 2)
legend("topleft", legend = c("Точная y'", "Численная y'"), lty = c(1,2), bty = "n")

# 3.3 вторая производная
plot(x, d2y_exact_vals, type = "l", lwd = 2, main = "Вторая производная: точная vs прибл.",
     xlab = "x", ylab = "y''")
lines(x, d2y_num, lty = 2, lwd = 2)
legend("topleft", legend = c("Точная y''", "Численная y''"), lty = c(1,2), bty = "n")

# 3.4 график ошибок первой производной
plot(x, abs_err_dy, type = "l", lwd = 2, main = "Абсолютная ошибка y'",
     xlab = "x", ylab = "|ошибка|")

# 4. таблица значений в узлах (первые 10 строк)
res_df <- data.frame(x = x, y = y, y_prime_exact = dy_exact_vals, y_prime_num = dy_num,
                     abs_err_dy = abs_err_dy, y_d2_exact = d2y_exact_vals, y_d2_num = d2y_num,
                     abs_err_d2y = abs_err_d2y)
print(head(res_df, 10))
#write.csv(res_df, "lab5_diff_results.csv", row.names = FALSE)