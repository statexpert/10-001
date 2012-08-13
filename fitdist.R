dt <- read.table(file="data.csv", header=TRUE, sep=";")
library("fitdistrplus")
library("boot")
library("parallel")

dists <- c("norm", "logis", "beta") # типы распределений
itr = 10000 # количество итераций

dep <- as.numeric(dt$ln.K.)

## Функция расчитывает AIC для всех поддерживаемых типов распредления
## Возвращает наиболее подходящий (по AIC) тип распределения
define.dist <- function(data) {
  data <- as.numeric(na.omit(data))
  aic <- mclapply(dists, FUN = function(x) fitdist(data = data, distr = x)$aic)
  names(aic) <- dists
  aic <- sort(data.frame(aic))
  return(names(aic[1]))
}

## Функция расчитывает AIC для всех поддерживаемых типов распредления
## Возвращает AIC для всех поддерживаемых типов распредления
aic.fits <- function (data) {
  data <- as.numeric(na.omit(data))
  aic <- mclapply(dists, FUN = function(x) fitdist(data = data, distr = x)$aic)
  names(aic) <- dists
  aic <- sort(data.frame(aic))
  return(aic)
}

## Функция для бутстрапа, расчитывает параметры заданного распределения (mu, sigma, nu, tau)
fit.dist <- function(data, indices, dist) {
  data <- as.numeric(na.omit(data[indices]))
  fit <- fitdist(data = data, distr = dist)
  result <- as.vector(fit$estimate)
  return(result)
}

## Функция бустрапит параметры распределения
## Взвращает вектор: mu, ниж. 95%, верх. 95% (то же для остальных параметров)
boot.dist <- function(x, dist) {
  x <- as.numeric(na.omit(x))
  # Задаём параметры для геренатора (для воспроизводимости результатов)
  set.seed(1234)
  # Расчитываем параметры распределения
  boot.out <- boot(x, statistic=fit.dist, dist = dist, R = itr, parallel = "multicore", ncpus = detectCores())
  # Формируем вывод
  result <- data.frame(boot.out$t0)
  est <- data.frame(boot.out$t0)
  ci <- mclapply(1:length(boot.out$t0), FUN = function(x) boot.ci(boot.out, index = x, type = "bca")$bca[4:5])
  ci.l <- as.numeric(lapply(ci, function(x) x[1]))
  ci.u <- as.numeric(lapply(ci, function(x) x[2]))
  result <- cbind(est, ci.l, ci.u)
  return(result)
}

# Таблица типов распределения
dist.def <- define.dist(dep)

# Таблица AIC для типов распределения
dist.aic <- aic.fits(dep)

# Таблица параметров распределния + доверительные интервалы
dist.params  <- boot.dist(dep, dist = define.dist(dep))

output <- list(dist.def, dist.aic, dist.params)
names(output) <- c("Distributions", "AIC", "Bootstraped.parametrs")

descdist(dep, boot = itr)
print(output)