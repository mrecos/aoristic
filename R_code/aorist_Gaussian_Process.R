############### FUNCTIONS ####################
## Simulate GP based on fixed sigma, rho, and eta
## calls gp-predict.stan
## this sims over same data (Y,X), but uses infered hyperparameters
sim_GP_y <- function(y1_aor, x1_aor, sigma_sq, rho_sq, eta_sq, iter = iter, chains = chains){
  stan_model_loc <- "/Users/mattharris/Documents/R_Local/aoristic/stan_models/"
  sim_fit <- stan(file=paste0(stan_model_loc, "gp-predict_SE.stan"), 
                  data=list(x1=x1_aor, y1=y1_aor, N1=length(x1_aor),
                  x2=x1_aor, N2=length(x1_aor), eta_sq=eta_sq, 
                  rho_sq=rho_sq, sigma_sq=sigma_sq),
                  iter=iter, chains=chains)
  
  sims <- extract(sim_fit, permuted=TRUE)
  return(sims)
}

## takes simulated Y and makes it ready to plot
melt_sim <- function(sim, x_data, param_label){
  tmp <- adply(sim, 2)
  # plot(sim_data[,2], type = "l") same as; plot(sims$y2[1,], type="l")
  # each col of sim_data is each row of sims$y2; 
  # each col sim_data is a complete simulation across all data points x (n = 201)
  # except row 1 is added as a factor of 201 levels for ID
  tmp <- melt(tmp)
  tmp$model <- param_label
  names(tmp) <- c("xid", "group1", "y", "param")
  tmp <- mutate(tmp, x=x_data[xid])
  return(tmp)
}

## like melt_sim, but for two paramaters 
melt_sim_multiparam <- function(sim, x_data, param_label1, param_label2){
  tmp <- adply(sim, 2)
  # plot(sim_data[,2], type = "l") same as; plot(sims$y2[1,], type="l")
  # each col of sim_data is each row of sims$y2; 
  # each col sim_data is a complete simulation across all data points x (n = 201)
  # except row 1 is added as a factor of 201 levels for ID
  tmp <- melt(tmp)
  tmp$param1 <- param_label1
  tmp$param2 <- param_label2
  names(tmp) <- c("xid", "group1", "y", "param1", "param2")
  tmp <- mutate(tmp, x=x_data[xid])
  return(tmp)
}

# simple power function to emulate pow() in stan
r.pow <- function(x,p){x^p}

# http://stackoverflow.com/a/17257422/2259277
number_ticks <- function(n){
  function(limits) pretty(limits, n)
}

############### END FUNCTIONS ####################

############### INCLUDES #########################
## load the required packages
library("rstan")
library("dplyr")
library("forcats")
library("plyr")
# devtools::install_github("hadley/ggplot2")
library("ggplot2")
# devtools::install_github("hrbrmstr/ggalt")
library("ggalt")
# library("reshape2")
library("matrixcalc")
library("cowplot")
library("MASS")
library("lattice")
############### END INCLUDES ######################

############### INITIALIZATION ####################
eta_sq = 1
rho_sq = 0.01
# l = sqrt((1/rho_sq)/2)
# rho_sq_test = 1/(2*l^2)
sigma_sq = 0.05
x <- seq(-10, 10, 0.1)
n <- length(x)
iter = 2000
chains = 3
rho_sq_range = c(0.001, 0.01, 0.1, 1, 10, 100) # for faceted kernel plot
eta_sq_prior = c(0.0001, 0.01, 1)
eta_class <- c("eta 1e-04","eta 0.01","eta 1")
rho_sq_prior = c(0.0001, 0.01, 1)
rho_class <- c("rho 1e-04","rho 0.01","rho 1")
# new data for fit
x1 <- c(-6, -3, -1, 0, 2,4,7,9)
y1 <- c(-2, 0, 1, -1, 0,2.5,0,3)
stan_model_loc <- "/Users/mattharris/Documents/R_Local/aoristic/stan_models/"
plots_loc <- "/Users/mattharris/Documents/R_Local/aoristic/"
############### END INITIALIZATION ####################

############### CODE AND PLOTS ########################

############### Fit STAN model to fixed rho and eta & plot
## Simulate with a few noise-free data points to get posterior.
## still fixed values for eta and rho
## Again pretend the noise is almost zero, but not quite.

x1_aor <- my_aor2$bin.no
y1_aor <- my_aor2$aorist

fit2 <- stan(file=paste0(stan_model_loc, "gp-predict_SE.stan"),
             data=list(x1=x1_aor, y1=y1_aor, N1=length(x1_aor),
             x2=x1_aor, N2=length(x1_aor),eta_sq=eta_sq, 
             rho_sq=rho_sq, sigma_sq=sigma_sq),
             iter=iter, chains=chains)
sims2 <- extract(fit2, permuted=TRUE)

## Rearrange the data for plotting

## NEW dplyr summarise range flow...
post_summarise <- data.frame(sims2$y2) %>%
  summarise_each(funs(`0.50`  = mean,
                      `0.025` = quantile(., p = 0.025),
                      `0.25`  = quantile(., p = 0.25),
                      `0.75`  = quantile(., p = 0.75),
                      `0.975` = quantile(., p = 0.975))) %>%
  gather() %>%
  separate(key, c("key", "stat"), sep = "_") %>%
  mutate(x = rep(1:42,5),
         stat = as.factor(stat))
post_summarise_rect <- spread(post_summarise, stat, value)

ggplot() +
  theme_bw() +
  # geom_point(data=data.frame(x=x1_aor, y=y1_aor), aes(x=x1_aor, y = y1_aor), size = 1) +
  # geom_point(data=data.frame(x=x1_aor, y=y1_aor), aes(x=x1_aor, y = y1_aor),
             # size = 3, color = "black", shape = 1) +
  geom_rect(data = post_summarise_rect, aes(xmin = x, xmax = x + 1,
            ymin = `0.025`, ymax = `0.975`), fill = "goldenrod", color = NA,
            alpha = 0.35) +
  geom_rect(data = post_summarise_rect, aes(xmin = x, xmax = x + 1,
           ymin = `0.25`, ymax = `0.75`), fill = "skyblue", color = NA,
           alpha = 0.75) +
  geom_step(data = post_summarise,
            aes(x = x, y = value, group = stat, 
                color = stat, size = stat)) +
  annotate("text", x = 0.5, y = c(0.65,0.45, 0.23, 0.02, -0.2), 
           label = c("97.5%", "75%", "50%", "25%", "2.75%"), size = 3,
           family="TrebuchetMS-Italic") +
  annotate("text", x = 30, y = -0.41, 
           label = paste0("Parameters: rho_sq = ", rho_sq, ", eta_sq = ", eta_sq, ", sigma_sq = ", sigma_sq),
           family="TrebuchetMS-Italic", size = 3.5) +
  scale_color_manual(values = c("#3B528BFF", "steelblue4", "firebrick3", "steelblue4",
                                "#3B528BFF"), name = "CI") +
  # scale_alpha_manual(values = c(0.5, 0.85, 1, 0.85, 0.5), name = "CI") +
  scale_size_manual(values = c(0.3, 0.3, 1, 0.3, 0.3), name = "CI") +
  scale_x_continuous(breaks = number_ticks(16), limits = c(0,42)) +
  scale_y_continuous(breaks = number_ticks(10)) +
  labs(title = "Aoristic Weight Summed by Time Interval",
       subtitle = "Credible Interval of Gaussian Process with RBF kernel",
       # caption=paste0("Parameters: rho_sq = ", rho_sq, ", eta_sq = ", eta_sq, ", sigma_sq = ", sigma_sq),
       x = "100 yr Time Step", 
       y = "Sum Aoristic Weight") +
  theme(
    axis.text.x = element_text(angle = 0, size = 8, family = "Trebuchet MS"),
    axis.text.y = element_text(size = 8, family = "Trebuchet MS"),
    axis.title = element_text(size = 10, family = "Trebuchet MS"),
    plot.caption = element_text(size = 10, hjust=0, 
                                family = "TrebuchetMS-Bold"),
    plot.title=element_text(size = 12, family="TrebuchetMS-Bold"),
    plot.subtitle=element_text(family="TrebuchetMS-Italic"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none")

# ggsave(filename = paste0(plots_loc, "GP_CI_range_2.png"),
#        width = 7, height = 5)

post_full <- data.frame(sims2$y2) %>%
  data.frame() %>%
  gather() %>%
  mutate(x = rep(1:42,each = nrow(sims2$y2)),
         group = rep(1:nrow(sims2$y2), times = 42))
 
ggplot() +
  geom_step(data = post_full, aes(x = x, group = group, y = value), 
            alpha = 0.05, size = 0.25) +
  geom_step(data = filter(post_summarise, stat == "0.50"),
            aes(x = x, y = value), color = "white", size = 0.5) +
  theme_bw() +
  scale_x_continuous(breaks = number_ticks(16), limits = c(0,42)) +
  scale_y_continuous(breaks = number_ticks(10)) +
  labs(title = "Aoristic Weight Summed by Time Interval",
       subtitle = "5000 Posterior Samples from Gaussian Process with RBF kernel",
       # caption=paste0("Parameters: rho_sq = ", rho_sq, ", eta_sq = ", eta_sq, ", sigma_sq = ", sigma_sq),
       x = "100 yr Time Step", 
       y = "Sum Aoristic Weight") +
  theme(
    axis.text.x = element_text(angle = 0, size = 8, family = "Trebuchet MS"),
    axis.text.y = element_text(size = 8, family = "Trebuchet MS"),
    axis.title = element_text(size = 10, family = "Trebuchet MS"),
    plot.caption = element_text(size = 10, hjust=0, 
                                family = "TrebuchetMS-Bold"),
    plot.title=element_text(size = 12, family="TrebuchetMS-Bold"),
    plot.subtitle=element_text(family="TrebuchetMS-Italic"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none")

# ggsave(filename = paste0(plots_loc, "GP_CI_range_3.png"),
#        width = 7, height = 5)


# plot it - from early summarise plyr methods (deleted)
fig2b <- ggplot(data = tmp2, aes(x=x, y=y)) +
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


############### Fit STAN model to subjective range of rho and eta & plot

## Simulate with a few noise-free data points for posterior over range of eta and rho
## Again pretend the noise is almost zero, but not quite.
param_fits4 <- NULL
for(i in 1:length(eta_sq_prior)){
  for(j in 1:length(rho_sq_prior)){
    ## Parameter value fixed in given example
    fit4 <- stan(file=paste0(stan_model_loc, "gp-predict_SE.stan"),
                 data=list(x1 = x1_aor, y1 = y1_aor, N1=length(x1_aor),
                 x2=x1_aor, N2=length(x1_aor), eta_sq=eta_sq_prior[i], 
                 rho_sq=rho_sq_prior[j], sigma_sq=sigma_sq),
                 iter=iter, chains=chains)
    
    sims4 <- extract(fit4, permuted=TRUE)
    
    ## Rearrange the data and plot it
    data4 <- adply(sims4$y2, 2)
    tmp4 <- melt(data4)
    names(tmp4) <- c("xid", "group1", "y")
    tmp4 <- mutate(tmp4, x=x1_aor[xid])
    tmp4$rho <- rho_sq_prior[j]
    tmp4$eta <- eta_sq_prior[i]
    sigma_sq < sigma_sq
    
    param_fits4 <- rbind(param_fits4, tmp4)
  }
}

param_fits4$rho_label <- paste0("rho ",param_fits4$rho)
param_fits4$eta_label <- paste0("eta ",param_fits4$eta)
param_fits4 <- mutate(param_fits4, 
                      rho_label = fct_relevel(rho_label, rho_class),
                      param_fits4, 
                      eta_label = fct_relevel(eta_label, eta_class))
dp <- data.frame(x=x1_aor, y=y1_aor)
## plot posreriors and data points for range of eta and rho
fig4 <- ggplot() +
  geom_line(data = param_fits4, aes(x=x, y=y, color = group1), alpha=0.1) +
  scale_color_viridis(discrete=TRUE, option = "viridis") +
  geom_point(data = dp, aes(x=x, y = y), size = 1) +
  theme_bw() +
  labs(title="Gaussian Process Fitted to Aoristic Weights",
       subtitle = expression("Squared Exponential Kernel with range of "*rho^2*" and "*eta^2),
       x = "Time Step", 
       y = "Aoristic Weight") +
  facet_grid(rho_label~eta_label) +
  scale_x_continuous(limits = c(0,50), breaks = seq(0,50,by=5), labels = seq(0,50,by=5)) +
  scale_y_continuous(limits = c(0,3.5), breaks = seq(0,3,by=0.5), labels = seq(0,3,by=0.5)) +
  theme(
    # panel.border = element_rect(colour = "gray90"),
    panel.grid = element_blank(),
    strip.background = element_rect(colour = "gray50", fill = "white"),
    strip.text.y = element_text(colour = "black", size = 10, face = "bold", family = "Trebuchet MS"),
    strip.text.x = element_text(colour = "black", size = 10, face = "bold", family = "Trebuchet MS"),
    panel.border = element_rect(colour = "gray90"),
    axis.text.x = element_text(angle = 0, size = 8, family = "Trebuchet MS"),
    axis.text.y = element_text(size = 8, family = "Trebuchet MS"),
    axis.title = element_text(size = 10, family = "Trebuchet MS"),
    plot.caption = element_text(size = 10, hjust=0, 
                                family = "Trebuchet MS"),
    plot.title=element_text(family="TrebuchetMS-Bold"),
    plot.subtitle=element_text(family="TrebuchetMS-Italic"),
    legend.position="none"
  )
plot(fig4)
ggsave(filename = paste0(plots_loc, "GP_SE_simulated_range_rho_eta_facet.png"),
       width = 10, height = 5)

############### Fit STAN model with rho and eta as random variables, estimate, and plot
## Call STAN with our made-up data points, x1_aor and y1_aor
## calls GP_estimate_eta_rho_SE.stan with [x,y] sample data
## sigma_sq is fixed, but rho_sq, eta_sq, and Y are estimated through MCMC
fit5 <- stan(file=paste0(stan_model_loc, "GP_estimate_eta_rho_SE.stan"),
             data=list(x1 = x1_aor, y1 = y1_aor, N1=length(x1_aor), 
             simga_sq = sigma_sq,
             x2=x1_aor, N2=length(x1_aor)),
             iter=iter, chains=chains)

## ESIMATE PARAMETERS
sims5 <- extract(fit5, permuted=TRUE)
## sims5 = 100:201 matrix for each estimated parameter
## 100 rows = 1/2 of interation [iter] (non-warmup samples) for each x value
## 100 plausible values for y at that datapoint x; 
## plot(density(sims5$y2[,1])) # density of estimates for y at x[1]
## 201 cols = each x value (data point)
## plot(sims5$y2[1,], type="l") # response across x for interation 1
## plot median (across all rows) of data point samples (columns)
## xx <- apply(sims5$y2,2,median); plot(xx, type="l")

## testing on known plausible values from model with fixed rho_sq and sigma_sq
## I plotted these as the first of these plots in the blog post
## to do so, just uncomment and assign these known values, then run the rest as usual
# rho_low <- 0.1
# rho_med <- 1.0
# rho_high <- 10
# eta_low <- 0.1
# eta_med <- 1.0
# eta_high <- 10

## derive 2.5%, 50%, and 97.5% quantile from simulated rho_sq and sigma_sq parameters
rho_low <- quantile(sims5$rho_sq, 0.025)
rho_med <- quantile(sims5$rho_sq, 0.50)
rho_high <- quantile(sims5$rho_sq, 0.975)
eta_low <- quantile(sims5$eta_sq, 0.025)
eta_med <- quantile(sims5$eta_sq, 0.50)
eta_high <- quantile(sims5$eta_sq, 0.975)

##### SIMULATE FOR 5%, 50%, and 95% levels for parameter estimate
## SIMULATE NEW RESPONSE BASED ON RANGE OF HYPERPARAMETERS
sim_low <- sim_GP_y(y1_aor, x1_aor, sigma_sq = sigma_sq, rho_low, eta_low, iter=iter, chains=chains)
sim_med <- sim_GP_y(y1_aor, x1_aor, sigma_sq = sigma_sq, rho_med, eta_med, iter=iter, chains=chains)
sim_high <- sim_GP_y(y1_aor, x1_aor, sigma_sq = sigma_sq, rho_high, eta_high, iter=iter, chains=chains)

## PREPARE DATA FOR PLOTTING
tmp_low <- melt_sim(sim_low$y2, x1_aor, "low rho/eta")
tmp_med <- melt_sim(sim_med$y2, x1_aor, "median rho/eta")
tmp_high <- melt_sim(sim_high$y2, x1_aor, "high rho/eta")
#### sample from iterations
tmp_low  <- dplyr::filter(tmp_low, 
            group1 %in% paste0("V",sample(1:((chains*iter)/2),100)))
tmp_med  <- dplyr::filter(tmp_med,  
            group1 %in% paste0("V",sample(1:((chains*iter)/2),100)))
tmp_high <- dplyr::filter(tmp_high, 
            group1 %in% paste0("V",sample(1:((chains*iter)/2),100)))

# bind simulations from various params and order factor for plot
tmp <- rbind(tmp_low, tmp_med, tmp_high)
tmp$param <- factor(tmp$param,levels(factor(tmp$param))[c(2,3,1)])
# create dataframe of training data for plotting
pnt_data <- data.frame(x=rep(x1_aor,3), y=rep(y1_aor,3), 
                       param = rep(c("low rho/eta","median rho/eta","high rho/eta"), each = length(x1_aor)))
pnt_data$param <- factor(pnt_data$param,levels(factor(pnt_data$param))[c(2,3,1)])

##### plot as 5%/5%, 50%/50%, and 95%/95% parameters
Sim_params_plot <- ggplot() +
  geom_line(data = tmp, aes(x=x, y=y, color = group1), alpha=0.1) +
  geom_point(data = pnt_data, aes(x = x, y = y)) +
  theme_bw() +
  facet_grid(param~.) +
  scale_color_viridis(discrete=TRUE) +
  # scale_x_continuous(limits = c(-10,10), breaks = seq(-10,10,by=2), labels = seq(-10,10,by=2)) +
  # scale_y_continuous(limits = c(-3,3), breaks = seq(-3,3,by=1), labels = seq(-3,3,by=1)) +
  labs(title = "Gaussian Process - Hyperparameter Estimation",
       subtitle = expression("Squared Exponential Kernel - Simulating "*rho^2*" and "*eta^2*", Sigma Fixed"),
       x = "X", 
       y = "Y") +
  theme(
    panel.border = element_rect(colour = "gray90"),
    axis.text.x = element_text(angle = 0, size = 6, family = "Trebuchet MS"),
    axis.text.y = element_text(size = 6, family = "Trebuchet MS"),
    axis.title = element_text(size = 8, family = "Trebuchet MS"),
    plot.caption = element_text(size = 8, hjust=0, 
                                family = "Trebuchet MS"),
    plot.title=element_text(family="TrebuchetMS-Bold"),
    plot.subtitle=element_text(family="TrebuchetMS-Italic"),
    legend.position="none",
    strip.background = element_rect(colour = "gray50", fill = "white"),
    strip.text.y = element_text(colour = "black", size = 8, face = "bold", family = "Trebuchet MS")
  )
plot(Sim_params_plot)
# ggsave(filename = "GP_SE_simulated_range_rho_eta_facet.png", width = 7, height = 4)

## Simualte as facted matrix of combinations for 5%, 50%, 95%
rho_sim_HML <- c(rho_low, rho_med, rho_high)
eta_sim_HML <- c(eta_low, eta_med, eta_high)
model_label <- c("low", "median", "high")
sim_hyperparams_plot <- NULL
for(i in 1:length(eta_sim_HML)){
  for(j in 1:length(rho_sim_HML)){
    ## SIMULATE NEW RESPONSE BASED ON RANGE OF HYPERPARAMETERS
    sim_hyperparams <- sim_GP_y(y1_aor, x1_aor, sigma_sq = sigma_sq, 
                                rho_sim_HML[j], 
                                eta_sim_HML[i], 
                                iter=iter, chains=chains)
    ## PREPARE DATA FOR PLOTTING
    label1 <- paste0("rho ", model_label[j])
    label2 <- paste0("eta ", model_label[i])
    tmp_sim_params <- melt_sim_multiparam(sim_hyperparams$y2, x1_aor, label1, label2)
    # bind simulations from various params and order factor for plot
    sim_hyperparams_plot <- rbind(sim_hyperparams_plot, tmp_sim_params)
  }
}

colnames(sim_hyperparams_plot) <- c("xid", "group1", "y", "rho", "eta", "x")
sim_hyperparams_plot$rho <- factor(sim_hyperparams_plot$rho,
                                   levels(factor(sim_hyperparams_plot$rho))[c(3,2,1)])
sim_hyperparams_plot$eta <- factor(sim_hyperparams_plot$eta,
                                   levels(factor(sim_hyperparams_plot$eta))[c(2,3,1)])
# create dataframe of training data for plotting
dp <- data.frame(x=x1_aor, y=y1_aor)
## PLOT
Sim_hyperparams_plot <- ggplot() +
  geom_line(data = sim_hyperparams_plot, aes(x=x, y=y, color = group1), alpha=0.1) +
  geom_point(data = dp, aes(x = x, y = y)) +
  theme_bw() +
  facet_grid(rho~eta) +
  scale_color_viridis(discrete=TRUE) +
  # scale_x_continuous(limits = c(-10,10), breaks = seq(-10,10,by=2), labels = seq(-10,10,by=2)) +
  # scale_y_continuous(limits = c(-3,3), breaks = seq(-3,3,by=1), labels = seq(-3,3,by=1)) +
  labs(title = "Gaussian Process - Hyperparameter Estimation",
       subtitle = expression("Squared Exponential Kernel - Simulating "*rho^2*" and "*eta^2*", Sigma Fixed"),
       x = "X", 
       y = "Y") +
  theme(
    panel.border = element_rect(colour = "gray90"),
    axis.text.x = element_text(angle = 0, size = 6, family = "Trebuchet MS"),
    axis.text.y = element_text(size = 6, family = "Trebuchet MS"),
    axis.title = element_text(size = 8, family = "Trebuchet MS"),
    plot.caption = element_text(size = 8, hjust=0, 
                                family = "Trebuchet MS"),
    plot.title=element_text(family="TrebuchetMS-Bold"),
    plot.subtitle=element_text(family="TrebuchetMS-Italic"),
    legend.position="none",
    strip.background = element_rect(colour = "gray50", fill = "white"),
    strip.text.y = element_text(colour = "black", size = 8, face = "bold", family = "Trebuchet MS")
  )
plot(Sim_hyperparams_plot)
# ggsave(filename = "GP_SE_simulated_facet_rho_eta_facet.png", width = 7, height = 6)


## slab plot of simulated parameter densities
## random sample from our prior (in STAN code) for plotting
## clearly this will differ from actual prior used in STAN code,
## but is same distribution and parameters for the Cauchy()
prior_dist <- rcauchy(nrow(sims5$rho_sq),0,5)
sims_plot_data <- data.frame(value = c(sims5$rho_sq, prior_dist,
                                       sims5$eta_sq, prior_dist), 
                             param = rep(c("rho", "eta"), each = (nrow(sims5$rho_sq))*2),
                             model1 = rep(c("posterior", "prior", "posterior", "prior"),
                                          each = nrow(sims5$rho_sq)))
# plot top plot (rho)
sims_plot_rho <- ggplot(data = dplyr::filter(sims_plot_data, param == "rho"), 
                        aes(value, color = model1, group = model1)) +
  geom_density() +
  xlim(c(0,100)) +
  theme_bw() +
  labs(title = "Hyperparameter Density",
       subtitle = "Prior = halfCauchy(0,5)",
       x = "", 
       y = "Density") +
  scale_color_manual(values = c("dodgerblue2", "firebrick2")) +
  theme(
    panel.border = element_rect(colour = "gray90"),
    axis.text.x = element_text(angle = 0, size = 6, family = "Trebuchet MS"),
    axis.text.y = element_text(size = 6, family = "Trebuchet MS"),
    axis.title = element_text(size = 8, family = "Trebuchet MS"),
    plot.caption = element_text(size = 8, hjust=0, 
                                family = "Trebuchet MS"),
    plot.title=element_text(family="TrebuchetMS-Bold"),
    plot.subtitle=element_text(family="TrebuchetMS-Italic"),
    strip.background = element_rect(colour = "gray50", fill = "white"),
    strip.text.y = element_text(colour = "black", size = 8, face = "bold", family = "Trebuchet MS"),
    legend.position="none"
  )
# plot bottom plot (eta)
sims_plot_eta <- ggplot(data = dplyr::filter(sims_plot_data, param == "eta"), 
                        aes(value, color = model1, group = model1)) +
  geom_density() +
  xlim(c(0,10)) +
  theme_bw() +
  labs(x = "Parameter Value", 
       y = "") +
  scale_color_manual(values = c("dodgerblue2", "firebrick2"), 
                     guide = guide_legend(title = NULL)) +
  theme(
    panel.border = element_rect(colour = "gray90"),
    axis.text.x = element_text(angle = 0, size = 6, family = "Trebuchet MS"),
    axis.text.y = element_text(size = 6, family = "Trebuchet MS"),
    axis.title = element_text(size = 8, family = "Trebuchet MS"),
    plot.caption = element_text(size = 8, hjust=0, 
                                family = "Trebuchet MS"),
    plot.title=element_text(family="TrebuchetMS-Bold"),
    plot.subtitle=element_text(family="TrebuchetMS-Italic"),
    legend.position="bottom",
    strip.background = element_rect(colour = "gray50", fill = "white"),
    strip.text.y = element_text(colour = "black", size = 8, face = "bold", family = "Trebuchet MS")
  )
# use cowplot to combine plots and label
slab_plot <- cowplot::plot_grid(sims_plot_rho, sims_plot_eta, 
                                labels = c("rho^2", "eta^2"),
                                nrow = 2, align = "v", hjust = -0.1, scale = 0.9, vjust = 3)

# cowplot::save_plot("slab_plot_halfcauchy.png", slab_plot,
#           nrow = 2, # and 2 rows
#           base_aspect_ratio = 1,
#           base_width = 6,
#           base_height = 2.5
# )