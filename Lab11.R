l <- 1.0          # Длина области [0, l]
a <- 1.0          # Фазовая скорость волны
N <- 50           # Число пространственных шагов
h <- l / N        # Шаг по пространству
tau <- 0.01       # Шаг по времени
T_max <- 2.0      # Конечное время моделирования
M <- 20           # Частота вывода графиков (каждые M шагов)

# Вычисление числа Куранта
c <- a * tau / h

cat("Длина области: l =", l)
cat("Скорость волны: a =", a)
cat("Число узлов: N =", N)
cat("Шаг пространства: h =", h)
cat("Шаг времени: tau =", tau)
cat("Число Куранта: c =", round(c, 4))

# Проверка условия устойчивости
if (abs(c) < 1) {
  cat("Условие устойчивости |c| < 1: выполнено")
} else {
  cat("Условие устойчивости не выполнено")
}
cat("\n")

# Начальное смещение u(x, 0) = μ1(x)
mu1 <- function(x) {
  return(sin(pi * x / l))
}

# Начальная скорость ∂u/∂t(x, 0) = μ2(x)
mu2 <- function(x) {
  return(0.0)
}

# Граничное условие u(0, t) = μ3(t)
mu3 <- function(t) {
  return(0.0)
}

# Граничное условие u(l, t) = μ4(t)
mu4 <- function(t) {
  return(0.0)
}

# Создание пространственной сетки
x <- seq(0, l, length.out = N + 1)

# Создание массивов для трех временных слоев
u_prev <- rep(0, N + 1)  # Слой t - tau (предыдущий)
u_curr <- rep(0, N + 1)  # Слой t (текущий)
u_next <- rep(0, N + 1)  # Слой t + tau (следующий)

# Заполнение начального условия при t = 0
for (j in 1:(N+1)) {
  u_curr[j] <- mu1(x[j])
}

for (j in 2:N) {
  # Вторая производная μ1 методом конечных разностей
  mu1_second_deriv <- (mu1(x[j+1]) - 2*mu1(x[j]) + mu1(x[j-1])) / (h^2)
  
  u_prev[j] <- mu1(x[j]) + tau * mu2(x[j]) + 
    0.5 * tau^2 * a^2 * mu1_second_deriv
}

# Граничные условия для первого шага
u_prev[1] <- mu3(tau)
u_prev[N+1] <- mu4(tau)

# Вычисление максимального числа шагов
k_max <- as.integer(T_max / tau)

# Текущее время
t <- tau
k <- 1

# Настройка окна для графиков
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
plot_counter <- 0
# Основной цикл
while (k <= k_max) {
  for (j in 2:N) {
    u_next[j] <- c^2 * u_curr[j+1] + 
      2 * (1 - c^2) * u_curr[j] + 
      c^2 * u_curr[j-1] - 
      u_prev[j]
  }
  
  # Граничные условия
  u_next[1] <- mu3(t)
  u_next[N+1] <- mu4(t)
  
  if (k %% M == 0 || k == k_max) {
    plot_counter <- plot_counter + 1
    
    if (plot_counter <= 9) {
      plot(x, u_next, type = "l", col = "blue", lwd = 2,
           xlab = "x", ylab = "u(x,t)",
           main = paste("t =", round(t, 3)),
           ylim = c(-1.5, 1.5))
      grid()
    }
    
    cat("Шаг", k, "из", k_max, "  t =", round(t, 4), "\n")
  }
  
  u_prev <- u_curr
  u_curr <- u_next
  
  t <- t + tau
  k <- k + 1
}

cat("Всего выполнено шагов:", k - 1)

t_final <- (k - 1) * tau
u_analytic <- sin(pi * x / l) * cos(pi * a * t_final / l)

# Вычисление погрешности
error <- abs(u_curr - u_analytic)
max_error <- max(error)
mean_error <- mean(error)

cat("Сравнение с аналитическим решением")
cat("Время t =", round(t_final, 4))
cat("Максимальная погрешность:", sprintf("%.6e", max_error))
cat("Средняя погрешность:", sprintf("%.6e", mean_error))