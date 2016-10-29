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


