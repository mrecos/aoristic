library("rstan")

################ retrieve W with some normal noise on W
# works
model_string <- "
data {
  int<lower=0> N;
  vector[N] a;
  vector[N] b;
  int<lower=0> t;
}
parameters {
  vector[N] w;
}
model {
  for (i in 1:N) {
    w[i] ~ normal(t / fabs(b[i] - a[i]), 0.000001);
  }
}
"
model_dat <- list(a = overlap_names$i.begin, 
                  b = overlap_names$i.end, 
                  t = time_step_size, 
                  N = length(overlap_names$i.begin))
fit <- stan(model_code = model_string, data = model_dat, 
            iter = 2000, chains = 2,  warmup=1000)
print(fit)


##### Estimate W with noise on W, a, and b
## Works, but kind of wild
model_string <- "
data {
  int<lower=0> N;
  vector[N] a;
  vector[N] b;
  int<lower=0> t;
}
parameters {
  vector[N] a_n;
  vector[N] b_n;
  vector[N] w;
}
model {
  a_n ~ normal(a, 0.00001);
  b_n ~ normal(b, 0.00001);
  for (i in 1:N) {
    w[i] ~ normal(t / fabs(b_n[i] - a_n[i]), 0.000001);
  }
}
"

model_dat <- list(a = overlap_names$i.begin, 
                  b = overlap_names$i.end, 
                  t = time_step_size, 
                  N = length(overlap_names$i.begin))
fit <- stan(model_code = model_string, data = model_dat, 
            iter = 2000, chains = 2,  warmup=1000)
print(fit)


########## calculate W with no model or params
model_string <- "
data {
int<lower=0> N;
vector[N] a;
vector[N] b;
int<lower=0> t;
}
parameters {

}
model {

}
generated quantities {
vector[N] w;
  for (i in 1:N) {
    w[i] = t / fabs(b[i] - a[i]);
  }
}
"
model_dat <- list(a = overlap_names$i.begin, 
                  b = overlap_names$i.end, 
                  t = time_step_size, 
                  N = length(overlap_names$i.begin))
fit <- stan(model_code = model_string, data = model_dat, 
            iter = 2000, chains = 2,  warmup=1000,
            algorithm="Fixed_param")
print(fit)

########## calculate W with model params WORKS!!!
model_string <- "
data {
int<lower=0> N;
vector[N] a;
vector[N] b;
int<lower=0> t;
}
parameters {
  vector[N] a_n;
  vector[N] b_n;
}
model {
  a_n ~ normal(a, 100);
  b_n ~ normal(b, 100);
}
generated quantities {
  vector[N] w;
  for (i in 1:N) {
    w[i] = t / fabs(b_n[i] - a_n[i]);
  }
}
"
model_dat <- list(a = overlap_names$i.begin, 
                  b = overlap_names$i.end, 
                  t = time_step_size, 
                  N = length(overlap_names$i.begin))
fit <- stan(model_code = model_string, data = model_dat, 
            iter = 2000, chains = 2,  warmup=1000)
print(fit)

########## calculate W with transformed params and then model
# this is a step to the full model... make W than model GP of W
# works
model_string <- "
data {
int<lower=0> N;
vector[N] a;
vector[N] b;
int<lower=0> t;
}
transformed data {
vector[N] w;
  for (i in 1:N) {
    w[i] = t / fabs(b[i] - a[i]);
  }
}
parameters {
  vector[N] w2;
}
model {
  w2 ~ normal(w,0.1);
}
generated quantities {

}
"
model_dat <- list(a = overlap_names$i.begin, 
                  b = overlap_names$i.end, 
                  t = time_step_size, 
                  N = length(overlap_names$i.begin))
fit <- stan(model_code = model_string, data = model_dat, 
            iter = 2000, chains = 2,  warmup=1000)
print(fit)


########## calculate W with transformed params and then model (non trivial)
# this is a step to the full model... make W than model GP of W
# note, W is for each site, not aggregated to span. how to add parameters?
# GP will be across each sites W, not the summed weights for time step.
# WORKING ON THIS
model_string <- "
data {
  //int<lower=0> N;
  int<lower=1> N1;
  int<lower=1> N2;
  vector[N1] a;
  vector[N1] b;
  int<lower=0> t;
  vector[N1] x1;
  vector[N2] x2;
  real sigma_sq;
  real eta_sq;
  real rho_sq;
}
transformed data {
  int<lower=1> N;
  vector[N1] w1;
  vector[N1+N2] x;
  vector[N1+N2] mu;
  cov_matrix[N1+N2] Sigma;
  //vector[N1] y1; // this is now w = calc below
  N = N1 + N2;
  for (i in 1:N1) w1[i] = t / fabs(b[i] - a[i]);
  for (n in 1:N1) x[n]  = x1[n];
  for (n in 1:N2) x[N1 + n] = x2[n];
  for (i in 1:N) mu[i] = 0;
  for (i in 1:N)
    for (j in 1:N)
      //// RBF, aka squared exp kernel, aka Gaussian kernel
      Sigma[i,j] = eta_sq*exp(-rho_sq*pow(x[i] - x[j],2)) +
      if_else(i==j, sigma_sq, 0.0);
}
parameters {
  vector[N2] w2;
}
model {
  vector[N] w;
  for (n in 1:N1) w[n] = w1[n];
  for (n in 1:N2) w[N1 + n] = w2[n];
  w ~ multi_normal(mu,Sigma);
}
generated quantities {

}
"
model_dat <- list(a  = overlap_names$i.begin, 
                  b  = overlap_names$i.end, 
                  t  = time_step_size, 
                  x1 = overlap_names$step, # what should this be?
                  N1 = length(overlap_names$i.begin),
                  x2 = overlap_names$step, # what should this be?
                  N2 = length(overlap_names$i.begin),
                  eta_sq = eta_sq, 
                  rho_sq = rho_sq, 
                  sigma_sq = sigma_sq,
                  iter = iter, chains = chains)
fit <- stan(model_code = model_string, data = model_dat, 
            iter = 2000, chains = 2,  warmup=1000)
print(fit)

sims_gp <- extract(fit, permuted=TRUE)
se_sim_data2 <- adply(sims_gp$w2, 2)
tmp2 <- melt(se_sim_data2)
names(tmp2) <- c("xid", "group", "W")
tmp2 <- mutate(tmp2, x = x1_aor[xid])

fig2b <- ggplot(data = tmp2, aes(x=x, y=W)) +
  geom_line(aes(group=group, color = group), alpha=0.2) +
  theme_bw() +
  geom_point(data=data.frame(x=x1_aor, y=y1_aor), aes(x=x1_aor, y = y1_aor), size = 3) +
  geom_point(data=data.frame(x=x1_aor, y=y1_aor), aes(x=x1_aor, y = y1_aor),
             size = 6, color = "black", shape = 1) +
  labs(caption=paste0("rho_sq = ", rho_sq, ", eta_sq = ", eta_sq, "sigma_sq = ", sigma_sq),
       x = "X", 
       y = "Y") +
  # scale_x_continuous(limits = c(-10,10), breaks = seq(-10,10,by=1), labels = seq(-10,10,by=1)) +
  # scale_y_continuous(limits = c(-3,3), breaks = seq(-3,3,by=1), labels = seq(-3,3,by=1)) +
  theme(
    panel.border = element_rect(colour = "gray90"),
    axis.text.x = element_text(angle = 0, size = 6, family = "Trebuchet MS"),
    axis.text.y = element_text(size = 6, family = "Trebuchet MS"),
    axis.title = element_text(size = 8, family = "Trebuchet MS"),
    plot.caption = element_text(size = 8, hjust=0, 
                                family = "Trebuchet MS"),
    plot.title=element_text(family="TrebuchetMS-Bold"),
    plot.subtitle=element_text(family="TrebuchetMS-Italic"),
    legend.position="none"
  )
plot(fig2b)










