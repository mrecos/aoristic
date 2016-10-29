##########################
# my first attempt at making this work.  This code works, but is not very effiecient
# I'll leave this as is, but the better implimentation is to use the 
# aorist function I wrote with all data.table and data.table::foverlaps()
# It is by far the fastest way to go (without better math or deep parallelization)
########################

# http://stackoverflow.com/a/17257422/2259277
number_ticks <- function(n){
  function(limits) pretty(limits, n)
}

library("ggplot2")
library("dplyr")
library("data.table")
library("tidyr")
library("reshape2")
library("scales")

time_step_size <- 100 # years
seq_begin <- -5000
seq_end   <- 0
sim_site_count <- 250

periods <- data.frame(
  period = c("a","b","c","d"),
  begin  = c(-5000, -4000, -3500, -1500),
  end    = c(-4000, -3500, -1500, -100)) %>%
  mutate(duration = abs(begin - end))

time_steps <- data.frame(
  step  = seq(1:(abs(seq_begin)/time_step_size)),
  begin = seq(seq_begin,(seq_end-time_step_size), by = time_step_size),
  end   = seq((seq_begin+time_step_size),seq_end, by = time_step_size))

# sites <- data.frame(
#   name  = c("Site_1", "Site_2", "Site_3", "Site_4"),
#   begin = c(-3500, -2000, -4800, -3600),
#   end   = c(-1500, -1500, -1200, -1800)) %>%
#   mutate(s_duration = abs(begin - end))

sites <- data.frame(
  name  = seq(1:sim_site_count),
  begin = sample(seq(-5000,-1000, by = 100), sim_site_count, replace = TRUE)) %>%
  mutate(end = begin + sample(seq(100, 3000, by = 100),
                              length(begin), replace = TRUE),
         end = ifelse(end > 0, 0, end))

# join by range using data.table::foverlaps()
time_steps_DT <- time_steps
sites_DT <- sites
setDT(time_steps_DT)
setDT(sites_DT)
setkey(time_steps_DT, begin, end)

# create overlap and report site and time_step
# problem is, this include periods that only overlap by being/end value
# that is: a site dated 1,500 - 1000 will pick up periods on either end 
# of the start/end if they fall on end/start of periods 
# e.g. period 2,500 to 1,500 would be pickedup because the single year of 1,500 overlaps
# but the site is clearly not in the 2,500 to 1,500 period.
# fiter to remove any periods that attach by the single begin/end date overlap
# (could figure out how to do this for index is needed)
overlap_names <- data.table::foverlaps(sites_DT, time_steps_DT, type = "any", which = FALSE) %>%
  filter(i.begin != end & i.end != begin)

# Calculate weights
# summarise weights by step size
# join in time_steps and NA -> 0 for steps with no sites
weights <- overlap_names %>%
  data.frame() %>%
  group_by(name) %>%
  mutate(site_count = n(),
         # W = (1/site_count),
         W = (time_step_size / abs(i.end - i.begin)),
         W = round(W,3)) %>%
  arrange(name, step)

W_site <- weights %>%
  group_by(name) %>%
  summarise(site_W = mean(W))

W_timestep <- weights %>%
  group_by(step) %>%
  summarise(step_W = sum(W),
            median_step_W = median(W),
            site_count = n()) %>%
  right_join(., time_steps) %>%
  replace_na(list(step_W = 0, site_count = 0, median_step_W = 0))

## Sum W to 1 test == PASS!
overlap_names %>%
  data.frame() %>%
  group_by(name) %>%
  mutate(site_count = n(),
         W = (1/site_count),
         W = round(W,3)) %>%
  summarise(w_sum = sum(W)) %>%
  summary()

plot_dat <- W_timestep %>%
  dplyr::select(begin, step_W, median_step_W, site_count) %>%
  rename(`Sum W` = step_W, `Median W` = median_step_W, 
         `Site Count` = site_count) %>%
  melt(id = "begin")

text_df <- data.frame(begin = c(periods$begin + 50, -3000), 
                      value = c(11,11,11,11,13),
                      txt = c("Early", "Middle", "Late", "Classical",""),
                      variable = "Sum W")

ggplot(plot_dat,aes(x = begin, y = value, color = variable)) +
  geom_vline(xintercept = periods$begin, 
             linetype = 2, color = "gray50", size = 0.50) +
  geom_step(size = 1.25) +
  scale_color_manual(values = c("orange", "gray40", "orangered3")) +
  geom_text(data = text_df, color = "black",
            aes(label = txt), hjust = "left",
            family = "TrebuchetMS-Italic") +
  theme_bw() +
  scale_x_continuous(breaks=seq(-5000,0,by = 250)) +
  facet_wrap(~variable, scale = "free_y", 
             nrow = 3, strip.position = "left") +
  labs(title="Distribution of Aoristic Weights for Archaeological Sequence",
       subtitle="Time Interval = 100 years",
       caption="Data: Simulated temporal distribution of 250 sites",
       x = "Years BCE") +
  ylab(NULL) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 12, face = "bold", 
                              family = "Trebuchet MS"),
    axis.text.x = element_text(angle = 90, size = 8, family = "Trebuchet MS"),
    axis.text.y = element_text(size = 8, family = "Trebuchet MS"),
    # axis.ticks.y = element_blank(),
    axis.title = element_text(size = 12, family = "Trebuchet MS"),
    plot.caption = element_text(size = 8, hjust=0, margin=margin(t=5), 
                                family = "Trebuchet MS"),
    plot.title=element_text(family="TrebuchetMS-Bold"),
    plot.subtitle=element_text(family="TrebuchetMS-Italic"),
    legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "gray98"),
    panel.spacing = unit(0.5, "lines")
  )

ggsave("/Users/mattharris/Dropbox/R/gomez_etal_spatio_temporal/sim_aoristic.png", 
       width = 8, height = 5)










