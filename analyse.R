dt <- read.table("data.csv", sep = ";", header = TRUE)

# Проверяем соотвествтие распределения зависимой переменной номральному
hist(dt$ln.K., freq = FALSE)
lines(density(dt$ln.K., bw = 0.5), col = "blue", lwd = 2)
qqnorm(dt$ln.K.)
qqline(dt$ln.K., col = "red", lwd = 2)
shapiro.test(dt$ln.K.)

# Корреляционный аализ зависимой пеерменной с добавками
cor.out <- sapply(dt[2:12], FUN = cor.test, y = dt$ln.K., method = "spearman")[4:3,]
rownames(cor.out) <- c("rho", "p")
# Вывод статистически значимых корреляций
sign.cor <- cor.out[,cor.out[2,] < 0.05]
print(sign.cor)

# Определяем тип распределения
y <- as.numeric(dt$ln.K.)
library("gamlss")
fits <- fitDist(y, type = "realAll", try.gamlss=TRUE)$fits
print(fits)
dist <- "SHASHo"
histDist(y, family = dist, density = TRUE)
fit <- gamlssML(y, family = dist)
#fit <- gamlss(y ~ 1, family = dist)
print(fit)
x <- rSHASHo(n = 1000, mu = fit$mu.coefficients, sigma = exp(fit$sigma.coefficients), nu = fit$nu.coefficients, tau = exp(fit$tau.coefficients))
plot(density(y), col = "red")
lines(density(x), col = "blue")

dist <- "exGAUS"
histDist(y, family = dist, density = TRUE)
fit <- gamlssML(y, family = dist)
fit <- gamlss(y ~ 1, family = dist)
print(fit)
x <- rexGAUS(n = 1000, mu = fit$mu.coefficients, sigma = exp(fit$sigma.coefficients), nu = exp(fit$nu.coefficients))
plot(density(y), col = "red")
lines(density(x), col = "blue")

fit <- lm(ln.K. ~  Cr + Cu + Fe + Fe.Cr + Mo + Nb + Ni + O2+ P + Si + Sn +  Zr, data = dt)
stepfit <- step(fit, direction="both")
summary(stepfit)