# ---- library ----
library(dplyr)
library(plm)
options(scipen = 9999)

# ---- data preparation ----
#dat = read.csv("https://www.dropbox.com/s/cbvgghuoyd9kpve/IA_corn.csv?dl=1")
dat = readr::read_csv("data/IA_corn.csv") %>%
  select(-X) %>% # drop first unnamed column
  arrange(fips, year) %>%
  as.matrix()


# partition the covariates
y = dat[,"yield", drop=F]
fe = dat[,"fips", drop=F]
parametric = cbind(dat[,"year",drop=F], dat[, grep("tbinneg1",colnames(dat)):grep("tbin40", colnames(dat)), drop=F])
# all of the daily weather variables have the string "jday" in them, likewise "soil" for soil
nonparametric = dat[, grep("jday|soil", colnames(dat)), drop=F]

# generate training and testing indices
tr_idx = which(dat[,"year"] < 2012)
te_idx = which(dat[,"year"] >= 2012)

ytr <- y[tr_idx, , drop = F]
yte <- y[te_idx, , drop = F]
Xtr <- cbind(fe, parametric)[tr_idx,]
Xte <- cbind(fe, parametric)[te_idx,]


Xtr_bar <- as_tibble( Xtr ) %>%
  group_by(fips) %>%
  mutate_at(vars(-group_cols()), mean)

ytr_dmd <- dat[tr_idx,] %>%
  as_tibble() %>%
  #mutate(fips = as.character(fips)) %>%
  group_by(fips) %>%
  summarise(yield_bar = yield - mean(yield)) %>%
  ungroup() %>%
  select(yield_bar) %>%
  as.matrix()

# time-demeaning
Xtr_dmd <- Xtr - as.matrix(Xtr_bar)
Xtr_dmd[,1] <- Xtr[,1] # reset because FE column is 0 because of demeaning

dim(Xtr); dim(Xtr_bar)




# remove columns with no variation (plm() and lm() drops these variables automatically)


Xtr_dmd <- Xtr_dmd[, - which(colnames(Xtr_dmd)=="fips")]
writeLines("remove FE column from demeaned data"); dim(Xtr); dim(Xtr_dmd )
Xtr_bar <- Xtr_bar[, - which(colnames(Xtr_bar)=="fips")]

# ---- fixed-effects regression with plm and lm ----
# regression equation with all X
form <- paste0("yield ~ ", paste(colnames(Xtr), collapse = "+"))

mplm <- plm(form, data = cbind(ytr, Xtr), model="within", index = "fips")

mlm <- lm(ytr ~ . - 1, data = as.data.frame(Xtr_dmd)) # -1 to get rid of intercept


# ---- hand-coded fixed effects regression ----

# calculate coefficients
# note that the X matrix has no columns with 1s
b <-  solve(t(Xtr_dmd) %*% Xtr_dmd) %*% (t(Xtr_dmd) %*% ytr_dmd)

# plm, lm, and hand-coded coefficiets are the same
all.equal(coef(mplm),coef(mlm), b)

# compute average y for each i
ytr_bar_uniq <- aggregate(yield ~ fips, dat[tr_idx, ], function(x) ybar = mean(x))$yield




dim(Xtr_bar); dim(b)
# recover intercepts
# remove distinct() to get matrix with # of rows equal to trainingsmatrix
fe <- ytr_bar_uniq - (as.matrix(distinct(Xtr_bar)) %*% b)

mfe <- fixef(mplm)

# hand-coded and model computed intercepts align
all.equal(as.numeric(fe), as.numeric(mfe))

# balanced panel


fe_expanded <- rep(fe, each = 32)
pred <- Xtr[, -which(colnames(Xtr)=="fips")] %*% b + fe_expanded

# as expected we get the y values with yhat + residuals
plot(ytr, pred + residuals(mplm)); abline(0, 1, col = 'red', lty = 'dashed', lwd = 2)

res <- as.numeric(ytr - pred)
mres <- as.numeric(residuals(mplm))
all.equal(res, mres)









## Define n and k parameters
n <- nrow(fe)
k <- ncol(Xtr_dmd)
t <- count(as_tibble(Xtr), fips) %>% distinct(n) %>% pull() # balanced panel
df <- n*(t-1) - k

# degrees of freedom align
all.equal(df, summary(mplm)$df[2])

## Calculate Variance-Covariance Matrix
VCV = 1/df * as.numeric(t(res)%*%res) * solve(t(Xtr_dmd) %*% Xtr_dmd)
dim(VCV)

## Standard errors of the estimated coefficients
se = sqrt(diag(VCV))
se
t_stat <- as.numeric(abs(b/se))

## Calculate p-value for a t-test of coefficient significance
pvalue = 2*pt(t_stat, df=df,lower.tail= F)

## concatenate into a single data.frame
beta_hat = data.frame(b = as.vector(b),se = se, t_stat = t_stat, pvalue = pvalue)
beta_hat

# standard errors are equal
sem <- summary(mplm)$coefficients[,2]
all.equal(se, sem)

# standard errors p. 271 in wooldrige (econometric analysis of cross section and panel data)


# ---- prediction test data ----
# nte <- as_tibble(Xte) %>% group_by(fips) %>% count()
#
# dd <- as_tibble(fe) %>% bind_cols(nte)
# fe_expanded_te <- data.frame(fe = rep(dd$yield_bar, dd$n),
#                              fips = rep(dd$fips, dd$n))
#
# pred <- Xte[, - c(vars_no_vari, which(colnames(Xte)=="fips"))] %*% b + fe_expanded_te$fe
#
# plot(unlist(pred), as.numeric(yte))
#
# mean((yte - pred)^2)
