set.seed(123)

alpha <- rep(-1, 3)    # alpha_i = -1
beta  <- rep(0.125, 3) # beta_i = 0.125
rho   <- rep(0.2, 3)   # rho_i = 0.2
q     <- rep(0.5, 3)   # q_i = 0.5

# Границы интегрирования
a <- -5
b <- 5
volume <- (b - a)^3

# Количество точек Монте-Карло
N <- 1e6

# Генерация N случайных точек (векторно)
x1 <- runif(N, a, b)
x2 <- runif(N, a, b)
x3 <- runif(N, a, b)

# Векторизованная функция p для всех точек
# p(x) = prod_{i=1..3} |x_i - alpha_i|^{beta_i} * exp(-rho_i * |x_i - alpha_i|^{q_i})
# реализуем по координатам и затем перемножаем
abs1 <- abs(x1 - alpha[1])
abs2 <- abs(x2 - alpha[2])
abs3 <- abs(x3 - alpha[3])

term1 <- (abs1^beta[1]) * exp(- rho[1] * (abs1^q[1]))
term2 <- (abs2^beta[2]) * exp(- rho[2] * (abs2^q[2]))
term3 <- (abs3^beta[3]) * exp(- rho[3] * (abs3^q[3]))

values <- term1 * term2 * term3  # длина N

# Оценка интеграла и ошибка
mean_val <- mean(values)
var_val  <- var(values)          # выборочная дисперсия значений p(x)
se_mc    <- sqrt(var_val / N)    # стандартная ошибка для средней p(x)
integral_est <- mean_val * volume
se_integral  <- se_mc * volume

# 95% доверительный интервал (приближенно, по CLT)
z95 <- 1.96
ci_lower <- integral_est - z95 * se_integral
ci_upper <- integral_est + z95 * se_integral

cat("Параметры (векторы):\n")
cat("alpha =", alpha, "\n")
cat("beta  =", beta, "\n")
cat("rho   =", rho, "\n")
cat("q     =", q, "\n\n")

cat("N =", N, "\n")
cat("Оценка интеграла =", integral_est, "\n")
cat("Стандартная ошибка =", se_integral, "\n")
cat(sprintf("95%% CI: [%.6g, %.6g]\n", ci_lower, ci_upper))
cat("Среднее p(x) по точкам =", mean_val, "\n")
cat("Дисперсия p(x) =", var_val, "\n")
