library(Hotelling)

dta <- read.delim("oxygen.DAT", header = FALSE, sep = "")
colnames(dta) <- c("X_1", "X_2", "X_3", "X_4", "Gender")
attach(dta)

X_bar_1 <- colMeans(dta[1:25, c(1, 3)])
X_bar_2 <- colMeans(dta[26:50, c(1, 3)])
S_1 <- var(dta[1:25, c(1, 3)])
S_2 <- var(dta[26:50, c(1, 3)])

n_1=n_2=25
p=2

#### testing H0: X_bar_1 - X_bar_2 = 0 
#### (the mean of oxygen volume while resting and strenous exercise between men and women are the same)

## Assuming equal variances
S_po <- ((n_1 - 1) * S_1 + (n_2 - 1) * S_2) / (n_1 + n_2 - 2)

T2 <- (X_bar_1 - X_bar_2) %*% solve(S_po * (1 / n_1 + 1 / n_2)) %*% (X_bar_1 - X_bar_2)
p_value <- 1 - pf((n_1 + n_2 - p - 1) / ((n_1 + n_2 - 2) * p) * T2, p, n_1 + n_2 - p - 1)
qf(0.95, p, n_1 + n_2 - p - 1)
## T2=81.7, p-value=7.17e-11<0.05 -> H0 is rejected.

## Assuming unequal variances, testing H0: X_bar_1 - X_bar_2 = 0 (the mean of oxygen volume while resting and strenous exercise between men and women are the same)
SS_1 <- S_1 %*% solve(S_1 / n_1 + S_2 / n_2) / n_1
SS_2 <- S_2 %*% solve(S_1 / n_1 + S_2 / n_2) / n_2
nu <- (p + p ^ 2) / ((sum(diag(SS_1 ^ 2)) + sum(diag(SS_1)) ^ 2) / n_1 + 
                       (sum(diag(SS_2 ^ 2)) + sum(diag(SS_2)) ^ 2) / n_2)
T2 <- (X_bar_1 - X_bar_2) %*% solve(S_1 / n_1 + S_2 / n_2) %*% (X_bar_1 - X_bar_2)
p_value <- 1 - pf((nu - p + 1) / (nu * p) * T2, p, nu - p + 1)

X_bar_1[1] - X_bar_2[1] + c(-1, 1) * qt(1 - 0.05 / (2 * p), nu) * 
  sqrt(S_1[1, 1] / n_1 + S_2[1, 1] / n_2)
X_bar_1[2] - X_bar_2[2] + c(-1, 1) * qt(1 - 0.05 / (2 * p), nu) * 
  sqrt(S_1[2, 2] / n_1 + S_2[2, 2] / n_2)
## T2=81.7, p-value=9.69e-11<0.05 -> H0 is rejected.

## Checking Normality
qqnorm(dta[1:25, 1])
qqline(dta[1:25, 1])

qqnorm(dta[1:25, 3])
qqline(dta[1:25, 3])

qqnorm(dta[26:50, 1])
qqline(dta[26:50, 1])

qqnorm(dta[26:50, 3])
qqline(dta[26:50, 3])
## Normality assumption seems to be not correct

## Assuming n_1-p and n_2-p are large, normal distribution is not necessary
T2 <- (X_bar_1 - X_bar_2) %*% solve((S_1 * 1 / n_1) + (S_2 * 1 / n_2)) %*% (X_bar_1 - X_bar_2)
p_value <- 1 - pchisq((n_1 + n_2 - p - 1) / ((n_1 + n_2 - 2) * p) * T2, p)
qchisq(0.95, p)
## T2=81.7, p-value=2.06e-09<0.05 -> H0 is rejected.

#### The 3 assumptions end with same result -> H0 is rejected.
## H0 is rejected and the alternative hypothesis, that there is a difference of average oxygen volume while resting and while strenous exercise between men and women, is accepted.