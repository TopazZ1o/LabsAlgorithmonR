a <- 0
b <- 1
# Тестовая задача из задания: I = 4 * ∫_0^1 sqrt(1 - x^2) dx = π
# Но обычно интеграл просто ∫0^1 f(x) dx; проверим оба варианта.
f <- function(x) { sqrt(1 - x^2) }   # подынтегральная для π/4
I_true_scaled <- pi / 4              # ∫_0^1 sqrt(1-x^2) dx = π/4

# Параметр для составного метода прямоугольников
n_rect <- 200   # увеличение n -> уменьшение погрешности

# Параметр для метода Гаусса: возьмём m = 5 и m = 11
m_list <- c(5, 11)

# 1. Составной метод прямоугольников (midpoint)
composite_midpoint <- function(f, a, b, n) {
  h <- (b - a) / n
  mid_points <- seq(a + h/2, b - h/2, length.out = n)
  return(h * sum(f(mid_points)))
}

I_mid <- composite_midpoint(f, a, b, n_rect)
cat(sprintf("Composite midpoint (n=%d): I = %.12f ; точное = %.12f ; abs err = %.3e\n",
            n_rect, I_mid, I_true_scaled, abs(I_mid - I_true_scaled)))

# 2. Метод Гаусса: таблица узлов и весов для m=5 и m=11 (стандартные узлы на [-1,1]).
#    Для интеграла по [a,b] нужна афинная замена переменной.

gauss_table <- function(m) {
  if (m == 5) {
    # узлы и веса для 5-точечной формулы Гаусса (по [-1,1])
    t <- c(-0.9061798459, -0.5384693101, 0.0, 0.5384693101, 0.9061798459)
    A <- c(0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851)
  } else if (m == 11) {
    # узлы и веса для 11-точечной формулы (значения округлены)
    t <- c(-0.9782286581, -0.8870625998, -0.7301520056, -0.5190961292, -0.2695431559,
           0.0,
           0.2695431559, 0.5190961292, 0.7301520056, 0.8870625998, 0.9782286581)
    A <- c(0.0556685671, 0.1255803695, 0.1862902109, 0.2331937646, 0.2628045442,
           0.2729250868,
           0.2628045442, 0.2331937646, 0.1862902109, 0.1255803695, 0.0556685671)
  } else {
    stop("Таблица весов/узлов реализована только для m=5 и m=11")
  }
  return(list(t = t, A = A))
}

gauss_integral <- function(f, a, b, m) {
  tab <- gauss_table(m)
  t <- tab$t; A <- tab$A
  # преобразование отрезка [-1,1] -> [a,b]: x = (b+a)/2 + (b-a)/2 * t
  x_nodes <- (a + b)/2 + (b - a)/2 * t
  integral <- sum(A * f(x_nodes)) * (b - a)/2
  return(integral)
}

# Посчитаем для каждого m
for (m in m_list) {
  I_gauss <- gauss_integral(f, a, b, m)
  cat(sprintf("Gauss m=%d: I = %.12f ; abs err = %.3e\n", m, I_gauss, abs(I_gauss - I_true_scaled)))
}

# 3. Вычисление I_scaled = 4 * integral и сравнение с π
I_mid_full <- 4 * I_mid
I_gauss_full <- sapply(m_list, function(m) 4 * gauss_integral(f, a, b, m))
cat(sprintf("Composite midpoint scaled (4*I): %.12f ; error vs pi = %.3e\n", I_mid_full, abs(I_mid_full - pi)))
for (i in seq_along(m_list)) {
  cat(sprintf("Gauss m=%d scaled (4*I): %.12f ; error vs pi = %.3e\n",
              m_list[i], I_gauss_full[i], abs(I_gauss_full[i] - pi)))
}

# 4. Графики: функция и точки Гаусса
xx <- seq(a, b, length.out = 500)
plot(xx, f(xx), type = "l", lwd = 2, main = "f(x) = sqrt(1 - x^2) и узлы Гаусса",
     xlab = "x", ylab = "f(x)")
# отметить узлы Гаусса для m=11
tab11 <- gauss_table(11)
xnodes11 <- (a + b)/2 + (b - a)/2 * tab11$t
points(xnodes11, f(xnodes11), pch = 16)
legend("topright", legend = c("f(x)", "Узлы Gauss m=11"), pch = c(NA, 16), lty = c(1,NA), bty = "n")

# 5. Табличный вывод
res_integration <- data.frame(
  method = c("midpoint", paste0("gauss_m", m_list)),
  I_value = c(I_mid, sapply(m_list, function(m) gauss_integral(f, a, b, m))),
  abs_err = c(abs(I_mid - I_true_scaled), sapply(m_list, function(m) abs(gauss_integral(f, a, b, m) - I_true_scaled)))
)
print(res_integration)
