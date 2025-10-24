a <- 0      # левый конец отрезка
b <- 1      # правый конец отрезка
n <- 30     # число делений (узлов будет n+1)
deg <- 2    # степень аппроксимирующего полинома (здесь квадратичный)

y_fun <- function(x) {
  x*(1 - x) + x^(1/2.5)
}

# 1. Построение сетки узлов
h <- (b - a) / n           # шаг сетки
x <- seq(a, b, by = h)     # узлы x_j, длина n+1
# seq гарантирует, что длина = n+1 (при корректном h)

# 2. Значения функции в узлах
y <- y_fun(x) # значения y_j = y(x_j)

# 3. Составляем матрицу для решения задачи наименьших квадратов
# ищем коэффициенты c0 + c1*x + c2*x^2 (deg=2)
# Матрица A: строки — узлы, столбцы — базисные функции 1, x, x^2
A <- cbind(1, x, x^2)

# Решение нормальной системы (A^T A) c = A^T y
fit <- lm(y ~ x + I(x^2))  # модель y = c0 + c1*x + c2*x^2
coeffs <- coef(fit)        # извлекаем коэффициенты c0, c1, c2

cat("Коэффициенты аппроксимации (c0,c1,c2):\n")
print(coeffs)

# 4. Значения аппроксимации в узлах и на плотной сетке для графика
y_approx_nodes <- predict(fit, newdata = data.frame(x = x))
# плотная сетка для гладкого графика
xx <- seq(a, b, length.out = 500)
yy_approx <- coeffs[1] + coeffs[2]*xx + coeffs[3]*xx^2
yy_exact  <- y_fun(xx)

# 5. Оценки погрешности
abs_err_nodes <- abs(y - y_approx_nodes)
max_err_nodes <- max(abs_err_nodes)
rmse_nodes <- sqrt(mean((y - y_approx_nodes)^2))

cat(sprintf("Макс. абсолютнаdя погрешность на узлах: %g\n", max_err_nodes))
cat(sprintf("RMSE на узлах: %g\n", rmse_nodes))

# оценим ошибку на плотной сетке (для графика)
abs_err_dense <- abs(yy_exact - yy_approx)
cat(sprintf("Макс. абсолютная погрешность на плотной сетке: %g\n",
            max(abs_err_dense)))

# 6. Графики
# График: точные значения и аппроксимация
plot(xx, yy_exact, type = "l", lwd = 2, main = "Аппроксимация: точная функция и аппрокс.",
     xlab = "x", ylab = "y")
lines(xx, yy_approx, lwd = 2, lty = 2)
points(x, y, pch = 16, cex = 0.8)         # узловые точки
legend("topleft", legend = c("Точная", "Аппроксимация", "Узлы"),
       lty = c(1,2,NA), pch = c(NA, NA, 16), bty = "n")

# График погрешности
plot(xx, abs_err_dense, type = "l", lwd = 2, main = "Абсолютная погрешность аппроксимации",
     xlab = "x", ylab = "|ошибка|")
abline(h = max_err_nodes, col = "red", lty = 3)
mtext(sprintf("max_err = %.3e ; RMSE = %.3e", max_err_nodes, rmse_nodes), side = 3)

# 7. Таблица значений (узлы, y, y_approx, abs_err)
result_table <- data.frame(x = x, y = y, y_approx = y_approx_nodes, abs_err = abs_err_nodes)
print(head(result_table, 10))
#write.csv(result_table, "lab4_results.csv", row.names = FALSE)
