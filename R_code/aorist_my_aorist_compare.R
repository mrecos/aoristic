###########################################
# Sandbox to compare execution time of functions, in this case archSeries::aorist()
# vs my_aorist2() [my verion with all DT and foverlaps()].  Test can be over a range of sites
# and reports graph of results, summarization of results, intuition of rate per 1,000 sites
# and profvis package to profile functions.
# This will remain pretty much unchangesd.
###########################################

library("ggplot2")
library("profvis")
library("dplyr")
library("tidyr")

time_step_size <- 100 # years
seq_begin <- -5000
seq_end   <- 0

sim_site_count_seq <- c(10,seq(50,50000,50))

execute_time_results <- matrix(nrow = length(sim_site_count_seq), ncol = 3)
colnames(execute_time_results) <- c("site_count", "aorist", "my_aorist")
pb <- txtProgressBar(min = 0, max = length(sim_site_count_seq), style = 3, char="*")
for(i in seq_along(sim_site_count_seq)){
  sim_site_count <- sim_site_count_seq[i]
  sites2 <- data.frame(
    name  = seq(1:sim_site_count),
    begin = sample(seq(-5000,-1000, by = 100), sim_site_count, replace = TRUE)) %>%
    mutate(end = begin + sample(seq(100, 3000, by = 100),
                                length(begin), replace = TRUE),
           end = ifelse(end > 0, 0, end)) %>%
    dplyr::select(begin, end, name) %>%
    mutate(begin = abs(begin), end = abs(end)) %>%
    dplyr::rename(Start = end, End = begin)
  
  aor.start.time <- Sys.time()
    aor <- aorist(sites2, end.date = 5000, bin.width = 100)
  aor.end.time <- Sys.time()
  my.aor.start.time <- Sys.time()
    my_aor <- my_aorist2(sites2, end.date = 5000, bin.width = 100)
  my.aor.end.time <- Sys.time()
  
  aor.execute.time <- aor.end.time - aor.start.time
  my.aor.execute.time <- my.aor.end.time - my.aor.start.time
  execute_time_results[i,1] <- sim_site_count
  execute_time_results[i,2] <- aor.execute.time
  execute_time_results[i,3] <- my.aor.execute.time
  setTxtProgressBar(pb, i)
}
close(pb)

execute_time_results <- data.frame(execute_time_results) %>%
  mutate(diff = my_aorist - aorist)

execute_time_summary <- execute_time_results %>%
  summarise(site_count_min = floor(min(site_count)),
            site_count_max = floor(max(site_count)),
            mean_aorist    = round(mean(aorist),3),
            mean_my_aorist = round(mean(my_aorist),3),
            max_diff       = round(min(diff),3), # reverse b/c negative time
            min_diff       = round(max(diff),4), # reverse b/c negative time
            mean_diff      = round(mean(diff),3)) %>%
  t()
print(execute_time_summary)

plot_dat <- execute_time_results %>%
  dplyr::select(site_count, aorist, my_aorist) %>%
  gather(model, time, -site_count)

ggplot(plot_dat, aes(x = site_count, y = time, group = model, color = model)) +
  geom_line(size = 1) +
  theme_bw()

aor_lm <- lm(aorist ~ site_count, data = execute_time_results)
aor_lm <- (as.numeric(coef(aor_lm))[2]) * 10000
my_aor_lm <- lm(my_aorist ~ site_count, data = execute_time_results)
my_aor_lm <- (as.numeric(coef(my_aor_lm))[2]) * 10000
message(paste0("my_aoristic is ", round(execute_time_results[2,"diff"],3),
  " seconds faster at ", execute_time_results[2,"site_count"], " , and is ",
  round(execute_time_results[201,"diff"],3), " seconds faster at 10,000 sites. For every 1000 sites added: aoristic = ", round(aor_lm,4), " seconds increase; my_aoristic = ",
               round(my_aor_lm,4), " seconds increase in execution time."))

profvis({
  aor_profile <- aorist(sites2, end.date = 5000, bin.width = 100)
  print(aor_profile)
})

profvis({
  my_aor_profile <- my_aorist2(sites2, end.date = 5000, bin.width = 100)
  print(my_aor_profile)
})

