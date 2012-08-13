dt <- read.table(file="data.csv", header=TRUE, sep=";")
dep <- as.numeric(dt$ln.K.)
library("gamlss")
library("boot")
library("parallel")

itr = 10000 # количество итераций

## Функция расчитывает AIC для всех поддерживаемых типов распредления
## Возвращает наиболее подходящий (по AIC) тип распределения
define.dist <- function(x) {
  x <- as.numeric(na.omit(x))
  fit.out <- fitDist(x, type = "realAll", try.gamlss = TRUE)
  return(fit.out$family[1])
}

## Функция расчитывает AIC для всех поддерживаемых типов распредления
## Возвращает AIC для всех поддерживаемых типов распредления
aic.fits <- function (x) {
  x <- as.numeric(na.omit(x))
  fit.out <- fitDist(x, type = "realAll", try.gamlss = TRUE)
  return(fit.out$fits)
}

## Функция расчитывает параметры заданного типа распределения
details.dist <- function(x, dist) {
  x <- as.numeric(na.omit(x))
  fit.out <- gamlssML(x, family = dist)
  details <- NULL
  params <- fit.out$parameters
  for (i in params) {
    assign(i, eval(parse(text=paste("fit.out$", i, sep=""))))
    details[[i]] <- get(i)
  }
  details$AIC <- AIC(fit.out)
  return(details)
}

## Функция для бутстрапа, расчитывает параметры заданного типа распределения
fit.dist <- function(data, indices, dist) {
  x <- as.numeric(na.omit(data[indices]))
  fit.out <- gamlssML(x, family = dist)
  result <- NULL
  params <- fit.out$parameters
  for (i in params) {
    assign(i, eval(parse(text=paste("fit.out$", i, sep=""))))
    result <- c(result, get(i))
  }
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
  result <- NULL
  params <-  c("mu", "sigma", "nu", "tau")
  for (i in 1:length(boot.out$t0)) {
    assign(params[i], boot.out$t0[i])
    assign(paste(params[i], ".conf", sep = ""), boot.ci(boot.out, index = i, type = "bca")$bca[4:5])
    result <- c(result, get(params[i]), get(paste(params[i], ".conf", sep = "")))
  }
  return(result)
}

# Определяем тип распределения
#dist <- define.dist(x)

# define.dist(x) не всегда даёт лучший варинат, поэтому определяем вручную
dist <- "SHASHo"

histDist(dep, family = dist, density=TRUE)

# Таблица типов распределения
dist.def <- define.dist(dep)

# Таблица AIC для типов распределения
dist.aic <- aic.fits(dep)

# Таблица параметров распределния + доверительные интервалы
dist.params  <- boot.dist(dep, dist = dist)

output <- list(dist.def, dist.aic, dist.params)

dist.params <- data.frame(dist.params)
rownames(dist.params) <- c("mu", "mu.ci.l", "mu.ci.u", "sigma", "sigma.ci.l", 'sigma.ci.u', "nu", "nu.ci.l", "nu.ci.u", "tau", "tau.ci.l", "tau.ci.u")

names(output) <- c("Distributions", "AIC", "Bootstraped.parametrs")

print(output)

fit <- gamlss(ln.K. ~  Cr + Cu + Fe + Fe.Cr + Mo + Nb + Ni + O2+ P + Si + Sn +  Zr, data = dt, family = dist)
fit2 <- stepGAIC(fit, direction="both")