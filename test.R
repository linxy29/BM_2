# HW8

library(readxl)
library(tidyverse)
library(gee)

data = read_excel("HEALTH.xlsx")
data_no = data %>%
  group_by(ID) %>%
  mutate(base = first(HEALTH),
         HEALTH = ifelse(HEALTH == "Poor", 0, 1)) %>% filter(TIME != 1) %>%
  ungroup() %>%
  mutate(TIME = plyr::mapvalues(TIME, from = c(2, 3, 4),to = c(3, 6, 9)),
         base = factor(base, levels = c("Poor", "Good")))
head(data_no)
