# k(x) - коэффициент теплопроводности
k_func <- function(x) {
  return(1 + 0.5 * x^2) 
}

# q(x) - коэффициент теплообмена
q_func <- function(x) {
  return(1 + x)
}

# f(x) - плотность источников (обычная часть)
f_func <- function(x) {
  return(x + 1)
}

# Если use_delta_source = TRUE, программа добавит источник тепла Q в точку x_source
use_delta_source <- TRUE
x_source <- 0.5  # Точка приложения источника (кси)
Q_source <- 10.0 # Мощность источника

a <- 0.0       # Левая граница
b <- 1.0       # Правая граница
N <- 100       # Число разбиений
h <- (b - a) / N
x_nodes <- seq(a, b, by = h) # Узлы сетки: x_0, x_1, ..., x_N

# Граничные условия 1-го рода: u(0)=u_a, u(1)=u_b
u_a <- 0
u_b <- 0

# кэфы
# Схема имеет вид: a_{i+1} * (y_{i+1}-y_i)/h - a_i * (y_i-y_{i-1})/h - d_i * y_i = -phi_i
# Приводим к виду: A_i*y_{i-1} - C_i*y_i + B_i*y_{i+1} = -F_i

A <- numeric(N + 1)
B <- numeric(N + 1)
C <- numeric(N + 1)
F_vec <- numeric(N + 1)

# Заполняем внутренние узлы (i от 2 до N, т.к. индексы в R с 1)
for (i in 2:N) {
  x_i <- x_nodes[i]
  
  # Вычисление a_i = k(x_{i-0.5})
  # Используем значение в полуцелом узле (середина отрезка)
  k_minus <- k_func(x_i - 0.5 * h)
  k_plus  <- k_func(x_i + 0.5 * h)
  
  # Коэффициенты A и B (деленные на h^2 для удобства уравнения)
  coef_A <- k_minus / (h^2)
  coef_B <- k_plus  / (h^2)
  
  # Вычисление d_i и phi_i методом сумматорных тождеств (через интеграл)
  # В простейшем случае (метод трапеций или средних) это значения в узле:
  d_val   <- q_func(x_i)
  phi_val <- f_func(x_i)
  
  # ОБРАБОТКА СОСРЕДОТОЧЕННОГО ИСТОЧНИКА (как в примере из файла)
  # Если x совпадает с x_source, добавляем Q/h к правой части
  if (use_delta_source && abs(x_i - x_source) < 1e-9) {
    phi_val <- phi_val + Q_source / h
  }
  
  # Коэффициент C (центральный элемент)
  coef_C <- coef_A + coef_B + d_val
  
  # Запись в массивы
  A[i] <- coef_A
  B[i] <- coef_B
  C[i] <- coef_C
  F_vec[i] <- phi_val
}

alpha <- numeric(N + 1)
beta  <- numeric(N + 1)

# Прямой ход (учет левого граничного условия u(0) = u_a)
alpha[2] <- 0
beta[2]  <- u_a

for (i in 2:N) {
  denom <- C[i] - A[i] * alpha[i]
  # Защита от деления на ноль (хотя для данной задачи это редкость)
  if(denom == 0) stop("Деление на ноль в прогонке")
  
  alpha[i + 1] <- B[i] / denom
  beta[i + 1]  <- (F_vec[i] + A[i] * beta[i]) / denom
}

# Обратный ход (находим y)
y <- numeric(N + 1)
y[N + 1] <- u_b # Правое граничное условие

for (i in N:1) {
  y[i] <- alpha[i + 1] * y[i + 1] + beta[i + 1]
}

# Построение графика
plot(x_nodes, y, type = "l", col = "blue", lwd = 2,
     main = paste("Решение задачи (N =", N, ")"),
     xlab = "x", ylab = "u(x)")
grid()

# Если был точечный источник, добавим линию, показывающую его положение
if (use_delta_source) {
  abline(v = x_source, col = "red", lty = 2)
  legend("topright", legend = c("u(x)", "Источник тепла"),
         col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 1))
}

# Вывод таблицы значений (пример)
results <- data.frame(x = x_nodes, u = round(y, 5))
print(head(results))
print(tail(results))
