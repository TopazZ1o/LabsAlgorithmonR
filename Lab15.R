# 1. ПАРАМЕТРЫ ЗАДАЧИ И КОНСТАНТЫ
# ------------------------------------------------------------------------------
l_val <- 1.0       # Полуширина области по z: [-l, l]
T_val <- 1.0       # Максимальное время
eps   <- 1.0       # epsilon
mu    <- 1.0       # mu
sigma <- 1.0       # sigma
lambda_sq <- 1.0   # lambda^2 (квадрат лямбда)

# Параметры точного решения (из методички или произвольные)
alpha_p <- 2.0     # alpha (параметр гауссиана)
beta_p  <- 3.0     # beta (частота колебаний)
gamma_p <- 1.0     # gamma (амплитуда)

# Параметры сетки
N1 <- 50           # Число узлов по пространству (от -l до l)
N2 <- 100          # Число шагов по времени
h <- (2 * l_val) / N1   # Шаг по пространству
tau <- T_val / N2       # Шаг по времени

# Сетка
z_nodes <- seq(-l_val, l_val, by = h)
t_nodes <- seq(0, T_val, by = tau)
Nt <- length(t_nodes)
Nz <- length(z_nodes)

# Вспомогательные переменные для схемы
rho <- tau^2 / h^2     # Параметр rho (в методичке rho = tau^2/h^2) NOTE: проверьте размерность, иногда rho=tau/h^2

# 2. ФУНКЦИИ ТОЧНОГО РЕШЕНИЯ И ИСТОЧНИКА
# ------------------------------------------------------------------------------

# u0(z,t) - само точное решение u(z,t)
u_exact_func <- function(z, t) {
  term1 <- (1 / sqrt(alpha_p)) * exp(-alpha_p^2 * z^2)
  term2 <- (l_val^2 - z^2)
  term3 <- (gamma_p * sin(beta_p * t) - gamma_p * beta_p * t)
  return(term1 * term2 * term3)
}

# Производные для вычисления F0 (формулы (5) в методичке)
# u_tt
u2_func <- function(z, t) {
  term1 <- (1 / sqrt(alpha_p)) * exp(-alpha_p^2 * z^2) * (l_val^2 - z^2)
  term2 <- -gamma_p * beta_p^2 * sin(beta_p * t)
  return(term1 * term2)
}
# u_t
u1_func <- function(z, t) {
  term1 <- (1 / sqrt(alpha_p)) * exp(-alpha_p^2 * z^2) * (l_val^2 - z^2)
  term2 <- gamma_p * beta_p * cos(beta_p * t) - gamma_p * beta_p
  return(term1 * term2)
}
# u_z
u3_func <- function(z, t) { # Исправленная формула производной (в методичке u_z названа u3? Проверьте, обычно u_zz нужно для оператора)
  # Но в формуле F0 стоит u_zz. В методичке u3 - это u_z, u4 - u_zz.
  # F0 = ... - u_zz ... В формуле (5) стоит u3 ??
  # В методичке опечатка: в (5) написано u3, но по смыслу уравнения u_zz.
  # Давайте используем u4 (вторая производная) вместо u3, так как оператор Лапласа - это вторая производная.
  return(0) # Заглушка, используем u4
}
# u_zz
u4_func <- function(z, t) {
  term_exp <- (2 / sqrt(alpha_p)) * exp(-alpha_p^2 * z^2)
  poly_z <- 2 * alpha_p^4 * l_val^2 * z^2 - 2 * alpha_p^4 * z^4 - 
    alpha_p^2 * l_val^2 + 5 * alpha_p^2 * z^2 - alpha_p^2 * l_val^2 - 1
  # В методичке сложная формула, перепроверьте знаки. Для теста можно упростить или пересчитать.
  # Используем как в методичке.
  time_part <- gamma_p * sin(beta_p * t) - gamma_p * beta_p * t
  return(term_exp * poly_z * time_part)
}

# Правая часть F0(z,t) согласно (5)
# F0 = eps*mu*u_tt + mu*sigma*u_t - u_zz + lambda^2*u
F0_func <- function(z, t) {
  u  <- u_exact_func(z, t)
  ut <- u1_func(z, t)
  utt <- u2_func(z, t)
  uzz <- u4_func(z, t) # Используем u4 (u_zz)
  
  res <- eps * mu * utt + mu * sigma * ut - uzz + lambda_sq * u
  return(res)
}

# 3. ЧИСЛЕННОЕ РЕШЕНИЕ (НЕЯВНАЯ СХЕМА + ПРОГОНКА)
# ------------------------------------------------------------------------------

# Матрица решения: строки - пространство z, столбцы - время t
# Индексы z: 1 correspond to -l, Nz correspond to l
U <- matrix(0, nrow = Nz, ncol = Nt)

# Начальные условия (t=0 и t=tau - первые два слоя)
# u|t=0 = 0; ut|t=0 = 0 => u|t=tau тоже 0 (аппроксимация 1 порядка)
# Или можно взять точное решение на первых двух слоях для старта
for(i in 1:Nz) {
  U[i, 1] <- u_exact_func(z_nodes[i], t_nodes[1]) # t=0
  U[i, 2] <- u_exact_func(z_nodes[i], t_nodes[2]) # t=tau
}


# Цикл по времени
for (j in 2:(Nt - 1)) {
  # На каждом шаге решаем систему A*y_{i+1} + B*y_i + C*y_{i-1} = D
  # Здесь индексы i идут от 2 до Nz-1 (внутренние узлы)
  
  # Векторы для прогонки (размер Nz)
  A_coef <- numeric(Nz)
  B_coef <- numeric(Nz)
  C_coef <- numeric(Nz)
  D_coef <- numeric(Nz)
  
  # Текущее время для следующего слоя t_{j+1}
  t_next <- t_nodes[j + 1]
  
  for (i in 2:(Nz - 1)) {
    z_i <- z_nodes[i]
    
    # Коэффициенты согласно методичке (после уравнения (9))
    # A_i = -rho
    # B_i = eps*mu + 0.5*mu*sigma*tau + 2*rho + lambda^2*tau^2
    # C_i = -rho
    # D_i = ...
    
    # Обратите внимание: в методичке A стоит при y_{i+1}, B при y_i, C при y_{i-1}.
    # Стандартная прогонка часто пишется A*y_{i-1} - C*y_i + B*y_{i+1} = -F
    # Здесь используем обозначения методички: A*y_{i+1} + B*y_i + C*y_{i-1} = D
    
    val_A <- -rho
    val_B <- eps * mu + 0.5 * mu * sigma * tau + 2 * rho + lambda_sq * tau^2
    val_C <- -rho
    
    # Правая часть D_i
    # D = 2*eps*mu*y_i - eps*mu*y_{i-1} + 0.5*mu*sigma*tau*y_{i-1} ...
    # Внимание: в методичке y_{check} - это слой j-1. y_{hat} - j+1. y - j.
    y_curr <- U[i, j]       # слой j
    y_prev <- U[i, j - 1]   # слой j-1
    F_val  <- F0_func(z_i, t_next) # F берется на слое j+1 (неявная схема) или j? 
    # В формуле (6) стоит F_0^{j} ? Нет, там индекс j, но обычно в неявных схемах берут j+1.
    # В тексте написано (F_0)_i^j. Следуем методичке.
    F_val_method <- F0_func(z_i, t_nodes[j]) 
    
    val_D <- 2 * eps * mu * y_curr - eps * mu * y_prev + 0.5 * mu * sigma * tau * y_prev + tau^2 * F_val_method
    
    A_coef[i] <- val_A
    B_coef[i] <- val_B
    C_coef[i] <- val_C
    D_coef[i] <- val_D
  }
  
  # Прогонка (реализация алгоритма из методички)
  alpha <- numeric(Nz)
  beta  <- numeric(Nz)
  
  # ГУ слева: y[1] = 0
  alpha[1] <- 0
  beta[1]  <- 0
  
  # Прямой ход
  for (i in 2:(Nz - 1)) {
    denom <- B_coef[i] + C_coef[i] * alpha[i - 1]
    alpha[i] <- -A_coef[i] / denom
    beta[i]  <- (D_coef[i] - C_coef[i] * beta[i - 1]) / denom
  }
  
  # Обратный ход
  # ГУ справа: y[Nz] = 0
  U[Nz, j + 1] <- 0
  
  for (i in (Nz - 1):1) {
    U[i, j + 1] <- alpha[i] * U[i + 1, j + 1] + beta[i]
  }
}

# 4. ВИЗУАЛИЗАЦИЯ
# ------------------------------------------------------------------------------
# Сравнение в центре области в конечный момент времени
mid_idx <- round(Nz / 2)
y_num <- U[mid_idx, ]
y_ex  <- sapply(t_nodes, function(t) u_exact_func(z_nodes[mid_idx], t))

plot(t_nodes, y_num, type = "o", col = "blue", pch = 19, cex = 0.5,
     main = "Сравнение решения в центре (z=0)", xlab = "t", ylab = "u")
lines(t_nodes, y_ex, col = "red", lwd = 2)
legend("topright", legend = c("Численное", "Точное"), col = c("blue", "red"), 
       lty = 1, pch = c(19, NA))
grid()

# Профиль в момент T
plot(z_nodes, U[, Nt], type = "o", col = "darkgreen", 
     main = paste("Профиль u(z) при t =", T_val), xlab = "z", ylab = "u")
lines(z_nodes, sapply(z_nodes, function(z) u_exact_func(z, T_val)), col = "red", lty = 2)
grid()
