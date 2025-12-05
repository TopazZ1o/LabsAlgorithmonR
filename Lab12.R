lx <- 1.0         # Длина области по x
ly <- 1.0         # Длина области по y (квадратная область)
N <- 40           # Число шагов по x
M <- 40           # Число шагов по y
hx <- lx / N      # Шаг сетки по x
hy <- ly / M      # Шаг сетки по y
eps <- 1e-5       # Точность сходимости (условие А2 из методички)
max_iter <- 10000 # Максимальное число итераций

cat("Параметры области")
cat("Длина по x: lx =", lx)
cat("Длина по y: ly =", ly)
cat("Область: G = [0,", lx, "] × [0,", ly, "]")
cat("Параметры сетки")
cat("Число шагов: N =", N, ", M =", M)
cat("Шаги сетки: hx =", round(hx, 4), ", hy =", round(hy, 4))
cat("Всего узлов:", (N+1) * (M+1))
cat("Внутренних узлов:", (N-1) * (M-1))

# Правая часть уравнения Пуассона: f(x,y)
# Для тестовой задачи: f(x,y) = -2π² sin(πx) sin(πy)
# Точное решение: u(x,y) = sin(πx) sin(πy)
f_func <- function(x, y) {
  return(-2 * pi^2 * sin(pi * x) * sin(pi * y))
}

# (условия Дирихле)
# Левая граница: u(0, y) при x = 0
boundary_left <- function(y) {
  return(0.0)
}

# Правая граница: u(lx, y) при x = lx
boundary_right <- function(y) {
  return(0.0)
}

# Нижняя граница: u(x, 0) при y = 0
boundary_bottom <- function(x) {
  return(0.0)
}

# Верхняя граница: u(x, ly) при y = ly
boundary_top <- function(x) {
  return(0.0)
}

# Точное решение для проверки (для тестовой задачи)
exact_solution <- function(x, y) {
  return(sin(pi * x) * sin(pi * y))
}

cat("Тестовая задача")
cat("f(x,y) = -2π² sin(πx) sin(πy)")
cat("Граничные условия: u = 0 на всех границах Г")
cat("Точное решение: u(x,y) = sin(πx) sin(πy)")

# Узлы сетки
x <- seq(0, lx, length.out = N + 1)
y <- seq(0, ly, length.out = M + 1)

# Матрица решения u[i,j]
# i - индекс по y (строки), j - индекс по x (столбцы)
u <- matrix(0, nrow = M + 1, ncol = N + 1)

# Граничные условия
for (j in 1:(N+1)) {
  u[1, j] <- boundary_bottom(x[j])      # y = 0
  u[M+1, j] <- boundary_top(x[j])       # y = ly
}

for (i in 1:(M+1)) {
  u[i, 1] <- boundary_left(y[i])        # x = 0
  u[i, N+1] <- boundary_right(y[i])     # x = lx
}

# Начальное приближение
for (i in 2:M) {
  for (j in 2:N) {
    u[i, j] <- 0.0
  }
}

cat("Начальные и граничные условия установлены")

# Для квадратной сетки N = M вычисляем:
# λ_max = cos(π/N) × cos(π/M)
# ω_opt = 2 / (1 + √(1 - λ_max²))

lambda_max <- cos(pi / N) * cos(pi / M)
omega_opt <- 2 / (1 + sqrt(1 - lambda_max^2))

cat("Параметры метода SOR")
cat("λ_max =", sprintf("%.8f", lambda_max))
cat("ω_opt =", sprintf("%.8f", omega_opt))

# Расчетная формула
# u[j,i]^(k+1) = ω × u_temp + (1-ω) × u[j,i]^(k)

# Коэффициенты разностной схемы на пятиточечном шаблоне
coef_x <- 1 / hx^2
coef_y <- 1 / hy^2
coef_center <- 2 * (coef_x + coef_y)

# Переменные для контроля сходимости
iter <- 0
diff_max <- Inf
convergence <- c()  # История сходимости

while (diff_max > eps && iter < max_iter) {
  diff_max <- 0.0
  iter <- iter + 1
  
  # Проход по всем внутренним узлам сетки
  for (i in 2:M) {
    for (j in 2:N) {
      # Сохраняем текущее значение
      u_old <- u[i, j]
      
      # Вычисление по методу Гаусса-Зейделя (предварительное значение)
      # Используем уже обновленные значения слева и снизу
      numerator <- coef_x * (u[i, j+1] + u[i, j-1]) +
        coef_y * (u[i+1, j] + u[i-1, j]) -
        f_func(x[j], y[i])
      
      u_temp <- numerator / coef_center
      
      # Применение релаксации (формула SOR)
      u[i, j] <- omega_opt * u_temp + (1 - omega_opt) * u_old
      
      # Вычисление максимального изменения
      diff <- abs(u[i, j] - u_old)
      if (diff > diff_max) {
        diff_max <- diff
      }
    }
  }
  
  # Сохранение истории для графика
  convergence <- c(convergence, diff_max)
  
  # Вывод промежуточных результатов каждые 50 итераций
  if (iter %% 50 == 0) {
    cat("Итерация", sprintf("%4d", iter), 
        ": max_diff =", sprintf("%.6e", diff_max), "\n")
  }
}

cat("Количество итераций:", iter)
cat("Достигнутая точность:", sprintf("%.6e", diff_max))

if (iter >= max_iter) {
  cat("Лимит итераций")
} else {
  cat("Сходимость достигнута (условие А2 выполнено)")
}

# Вычисление точного решения и погрешности
u_exact <- matrix(0, nrow = M + 1, ncol = N + 1)
error_abs <- matrix(0, nrow = M + 1, ncol = N + 1)

for (i in 1:(M+1)) {
  for (j in 1:(N+1)) {
    u_exact[i, j] <- exact_solution(x[j], y[i])
    error_abs[i, j] <- abs(u[i, j] - u_exact[i, j])
  }
}

max_error <- max(error_abs)
mean_error <- mean(error_abs)
rms_error <- sqrt(mean(error_abs^2))

cat("Сравнение с точным решением")
cat("Максимальная погрешность:", sprintf("%.6e", max_error))
cat("Средняя погрешность:", sprintf("%.6e", mean_error))
cat("Среднеквадратичная:", sprintf("%.6e", rms_error))

par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# Численное решение
image(x, y, t(u), 
      col = heat.colors(50),
      main = "Численное решение u(x,y)",
      xlab = "x", ylab = "y",
      cex.main = 1.2)
contour(x, y, t(u), add = TRUE, col = "black", lwd = 0.5, nlevels = 10)


# Точное решение
image(x, y, t(u_exact), 
      col = heat.colors(50),
      main = "Точное решение",
      xlab = "x", ylab = "y",
      cex.main = 1.2)
contour(x, y, t(u_exact), add = TRUE, col = "black", lwd = 0.5, nlevels = 10)

# Распределение погрешности
image(x, y, t(error_abs), 
      col = terrain.colors(50),
      main = "Абсолютная погрешность",
      xlab = "x", ylab = "y",
      cex.main = 1.2)
contour(x, y, t(error_abs), add = TRUE, col = "black", lwd = 0.5)

# Сходимость итераций
plot(1:length(convergence), convergence,
     type = "l", col = "blue", lwd = 2,
     log = "y",
     xlab = "Номер итерации k",
     ylab = "max |u(k+1) - u(k)| (log)",
     main = "График сходимости метода SOR",
     cex.main = 1.2)
grid(col = "gray", lty = "dotted")
abline(h = eps, col = "red", lty = 2, lwd = 2)
legend("topright", 
       legend = c("Погрешность", paste("Точность ε =", eps)),
       col = c("blue", "red"), 
       lwd = 2, 
       lty = c(1, 2),
       cex = 0.8)

# Значения в характерных точках
x_mid <- which.min(abs(x - 0.5))
y_mid <- which.min(abs(y - 0.5))

cat("Значения решения в характерных точках")
cat(sprintf("  u(0.5, 0.5) численное = %.6f", u[y_mid, x_mid]))
cat(sprintf("  u(0.5, 0.5) точное    = %.6f", u_exact[y_mid, x_mid]))
cat(sprintf("  Погрешность           = %.6e", error_abs[y_mid, x_mid]))

cat("Параметры решения")
cat("Минимум u:", sprintf("%.6f", min(u)))
cat("Максимум u:", sprintf("%.6f", max(u)))