# --- Функция rho(x) для вектора x = (x1, x2, x3) ---
rho3 <- function(x, a, b, p, q) {
  prod(abs(x - a)^b * exp(-p * abs(x - a)^q))
}

# --- Генерация случайных параметров ---
set.seed(42)
a <- runif(3, -5, 5)
b <- runif(3, 0.5, 4)
p <- runif(3, 0.1, 0.5)
q <- runif(3, 0.1, 2)

cat("Случайные параметры:\n")
print(data.frame(a=a, b=b, p=p, q=q))

# --- Область интегрирования: x_i ∈ [-a_max, a_max] ---
low <- rep(-5, 3)
high <- rep(5, 3)
volume <- prod(high - low)

# --- Метод Монте-Карло для тройного интеграла ---
monte_carlo_rho <- function(N, a, b, p, q, low, high) {
  X <- matrix(runif(3 * N, low[1], high[1]), ncol=3)
  vals <- apply(X, 1, function(x) rho3(x, a, b, p, q))
  I_est <- volume * mean(vals)
  return(I_est)
}

# --- Вычисляем интеграл для разных N ---
N_list <- c(1e3, 1e4, 1e5)
cat("\nМетод Монте-Карло для интеграла ∫ρ(x)dx:\n")
for (N in N_list) {
  I_mc <- monte_carlo_rho(N, a, b, p, q, low, high)
  cat(sprintf("  N=%7d   I≈%.6e\n", N, I_mc))
}
