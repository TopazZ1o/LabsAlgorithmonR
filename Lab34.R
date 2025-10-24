gauss_method <- function(A, b) {
  n <- nrow(A)
  Ab <- cbind(A, b)
  
  # Прямой ход
  for (i in 1:n) {
    # Выбор ведущего элемента
    max_row <- which.max(abs(Ab[i:n, i])) + i - 1
    if (Ab[max_row, i] == 0) stop("Система не имеет единственного решения")
    if (max_row != i) {
      Ab[c(i, max_row), ] <- Ab[c(max_row, i), ]
    }
    
    # Нормализация строки
    Ab[i, ] <- Ab[i, ] / Ab[i, i]
    
    # Обнуление ниже
    for (j in (i+1):n) {
      if (j <= n) {
        Ab[j, ] <- Ab[j, ] - Ab[j, i] * Ab[i, ]
      }
    }
  }
  
  # Обратный ход
  x <- numeric(n)
  for (i in n:1) {
    if (i == n) {
      x[i] <- Ab[i, n+1]
    } else {
      x[i] <- Ab[i, n+1] - sum(Ab[i, (i+1):n] * x[(i+1):n])
    }
  }
  
  return(x)
}

# A <- matrix(c(2,1,-1,
#               -3,-1,2,
#               -2,1,2), nrow=3, byrow=TRUE)
# b <- c(8,-11,-3)

A <- matrix(c(
  1.0, -0.3, -0.2, -0.1, -0.4,
  -0.2,  1.0, -0.1, -0.3, -0.4,
  -0.3, -0.2,  1.0, -0.3, -0.2,
  -0.1, -0.3, -0.4,  1.0, -0.2,
  -0.4, -0.2, -0.1, -0.2,  1.0
), nrow=5, byrow=TRUE)

b <- c(100, 200, -50, 150, 0)

answer<-gauss_method(A, b)
print(answer)