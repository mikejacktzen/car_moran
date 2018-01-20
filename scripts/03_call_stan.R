library(rstan)

# rmabd()

load(file=here('/data/proc/stan_dat_moran_tract_la.RData'))
ls()
names(stan_dat_moran)

stan_dat_moran$E %>% summary
stan_dat_moran$Z %>% summary
stan_dat_moran$X %>% summary

stan_dat_moran$E %>% is.na() %>% which
stan_dat_moran$Z %>% is.na() %>% which
stan_dat_moran$X %>% is.na() %>% which


test_moran = stan(file=here("/scripts/bym_moran_op.stan"),
                  data=stan_dat_moran,
                  iter=100,chain=1)

plot(test_moran)
