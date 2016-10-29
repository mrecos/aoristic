######################################
# wokring to compare my code to the archSeries::aorist() function()
# below is my original my_aorist() attempt which works, but is slowed by dlpyr
# then is my new my_aorist2() function (will likely change name) which is 100% data.table
# and data.table::foverlaps().  It is twice as fast as archSeries::aorist() and way way 
# faster then my original my_aorist().  I don't think I can squeeze much more out of 
# the new function without serious optimization in math or parallelization.
# This script is probably to remain as is as there is not much more to do with this function.
#######################################

devtools::install_github("davidcorton/archSeries")
library("archSeries")

sites2 <- sites %>%
  dplyr::select(begin, end, name) %>%
  mutate(begin = abs(begin), end = abs(end)) %>%
  # SWITCH BEGIN/END SINCE THIS IS abs() VALUES!
  dplyr::rename(Start = end, End = begin)

data <-  sites2

aor <- aorist(sites2, end.date = 5000, bin.width = 100)
# my_aor <- my_aorist(sites2, end.date = 5000, bin.width = 100)
my_aor2 <- my_aorist2(sites2, end.date = 5000, bin.width = 100)

plot(aor$aorist , my_aor2$aorist)
sqrt(mean((aor$aorist - my_aor2$aorist)^2))

################## My original attempt at aorist function
## works correctly, but profiling showed it did not scale well
## because of the mutate/dplyr stuff.
## I'll moth ball this, and make some additional versions
my_aorist <- function(data, weight = 1, start.date = 0, end.date = 2000, 
                      bin.width = 100, round_int = 4) {
  require(dplyr)
  require(data.table)
  time_steps <- data.frame(
    step  = seq(1:(abs(end.date)/bin.width)), # end not start b/c pos dates
    Start = seq(start.date,(end.date-bin.width), by = bin.width),
    End   = seq((start.date+bin.width),end.date, by = bin.width))
  setDT(time_steps)
  setDT(data)
  setkey(time_steps, Start, End)
  overlap_names <- data.table::foverlaps(data, time_steps, 
                                         type = "any", which = FALSE) %>%
    filter(i.Start != End & i.End != Start)
  aorist <- overlap_names %>%
    data.frame() %>%
    group_by(name) %>%
    mutate(site_count = n(),
           W = (bin.width / (i.End - i.Start)), # SWITCHED b/c of pos dates
           W = round(W,round_int)) %>%
    arrange(name, step) %>%
    group_by(step) %>%
    summarise(step_W = sum(W),
              median_step_W = median(W),
              site_count = n()) %>%
    right_join(., time_steps, by = "step") %>%
    replace_na(list(step_W = 0, site_count = 0, median_step_W = 0)) %>%
    mutate(bin.no = step,
           bin = paste0(Start, "-", End),
           aorist = step_W) %>%
    dplyr::select(bin, bin.no, aorist)
  return(aorist)
}


################## My 2nd attempt at aorist function
## cut down on dplyr/mutate in favor of data.table
## WORKING AND RMSE ~ zero
## TWICE AS FAST as archSeries::aorist() and multitudes faster than my original my_aorist()
my_aorist2 <- function(data, start.date = 0, end.date = 2000, bin.width = 100) {
  require(data.table)
  setDT(data)
  time_steps <- data.table(
    bin.no  = seq(1:(abs(end.date)/bin.width)), # end not start b/c pos dates
    Start = seq(start.date,(end.date-bin.width), by = bin.width),
    End   = seq((start.date+bin.width),end.date, by = bin.width))
  setkey(time_steps, Start, End)
  overlap_names <- data.table::foverlaps(data, time_steps, 
                   type = "any", which = FALSE)
  overlap_names <- overlap_names[i.Start != End & i.End != Start]
  overlap_names[, duration := (i.End - i.Start)] 
  overlap_names[, W := (bin.width / duration)] 
  ov_sum <- overlap_names[, .(aorist = sum(W),
                              median_step_W = median(W),
                              site_count = .N), keyby=.(bin.no)]
  setkey(ov_sum, bin.no)
  setkey(time_steps, bin.no)
  ov_sum2 <- ov_sum[time_steps, nomatch = 0]
  ov_sum2[, bin := paste0(Start, "-", End)]
  return(ov_sum2)
}


############# Annotated archSeries Function below - OG code + comments
function (data, weight = 1, start.date = 0, end.date = 2000, 
          bin.width = 100) 
{
  require(data.table)
  data <- data.table(cbind(data, weight)) # bind 1 weight
  data <- data[End >= start.date & Start <= end.date] # trim to start/end
  data[, `:=`(duration, End - Start)] # calc duration
  data[, `:=`(weight.per.unit, weight/duration)] # calc W per site
  breaks <<- seq(start.date, end.date, bin.width) # make bins
  labels <- numeric(0) # loop to make lables
  for (i in 1:(length(breaks) - 1)) {
    labels[i] <- paste(breaks[i], breaks[i + 1], sep = "-")
  }
  # make DT of time steps
  aorist <- data.table(bin = labels, bin.no = 1:length(labels), 
                       aorist = 0)
  # Loops to great W for time steps
  for (i in 1:length(labels)) {
    bin.1 <- breaks[i] # bin 1
    bin.2 <- breaks[i + 1] # bin 2
    # add lable[i] as column to data
    data[, `:=`(assign("a", labels[i]), 0)] 
    # if Start is in bin, put a weight in label[i] column
    data[Start >= bin.1 & Start < bin.2, `:=`(assign("a", 
        labels[i]), (bin.2 - Start) * weight.per.unit)]
    # if End is in bin, put a weight in label[i] column
    data[End > bin.1 & End <= bin.2, `:=`(assign("a", labels[i]), 
        (End - bin.1) * weight.per.unit)]
    # if Start/End span bin, put a wieght in label[1] column
    data[Start < bin.1 & End > bin.2, `:=`(assign("a", labels[i]), 
        bin.width * weight.per.unit)]
    # if Start/End is contained in bin, put weight in label[i] column
    data[Start >= bin.1 & End <= bin.2, `:=`(assign("a", 
        labels[i]), as.double(weight))]
    #
    aorist$aorist[i] <- sum(data[, get(labels[i])], na.rm = TRUE)
  }
  aorist
}