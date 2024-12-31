library(stringr)
library(lubridate)
library(tidyverse)
library(ggplot2)

# consider using SWAPtools


dat_path = "Hupsel.csv"

data = read.csv(dat_path)

data = data%>%
  select(time, station, lat, lon, RH, EV24)%>%
  filter(complete.cases(.))%>%
  mutate(time = ymd(time))

# check if all days are present
data = data%>%
  mutate(diff_days = c(1, diff(time)))
sum(data$diff_days>1, na.rm = TRUE)
# 5 gaps exist of more than 1 day
# find the last complete sequence
data = data%>%
  mutate(in_interval = cumsum(diff_days>1))%>%
  filter(in_interval==max(in_interval))%>%
  select(-in_interval, -diff_days)

data_seasonal = data%>%
  mutate(month = month(time),
         year = year(time))%>%
  mutate(month_cat = ifelse(month%in%c(10,11,12,1,2), 'Winter', 'Summer')%>%
           factor(levels = c('Summer', 'Winter')),
         year = ifelse(month%in%c(1,2), year-1, year))%>%
  group_by(year, month_cat)%>%
  summarise(precip = mean(RH),
            evap   = mean(EV24))
# Eleminating the data that is not available in AMIGO
data_seasonal = data_seasonal%>%
  filter(year>=2004)

data_seasonal%>%
  ggplot()+
  geom_col(aes(year, precip, fill = month_cat),
           position = position_dodge(), alpha = 0.7)+
  geom_line(aes(year, evap, colour = month_cat), size = 1)+
  geom_vline(aes(xintercept = 2008, linetype = "Wet"))+
  geom_vline(aes(xintercept = 2012, linetype = "Average"))+
  geom_vline(aes(xintercept = 2019, linetype = "Dry"))+
  # geom_hline(aes(yintercept = 430, colour = "Mean Summer precip"))+
  # geom_hline(aes(yintercept = 325, colour = "Mean Winter precip"))+
  guides(
    fill = guide_legend(title = "Precipitation"),
    colour = guide_legend(title = "Makkink Evaporation"),
    linetype = guide_legend(title = "")
  )
ggsave("seasonal.png", dpi = 300, width = 6, height = 3)

# Cummulative plot of precipitation - evaporation
data_excess = data%>%
  mutate(precip_excess = RH-EV24,
         month = month(time),
         year = year(time))%>%
  mutate(month_cat = ifelse(month%in%c(10,11,12,1,2), 'Winter', 'Summer')%>%
           factor(levels = c('Summer', 'Winter')))%>%
  # Format the day of year to begin from 1st March
  mutate(day_of_year = yday(time) - yday(dmy(paste0("01-03-", year)))+1)%>%
  mutate(day_of_year = ifelse(day_of_year>0,
                              day_of_year,
                              day_of_year + yday(dmy(paste0("31-12-", year-1)))))%>%
  mutate(year = ifelse(month%in%c(1,2), year-1, year))%>%
  group_by(year)%>%
  summarise(time = time,
            day_of_year = day_of_year, 
            year = year,
            cummulative_excess = cumsum(precip_excess))

mean_data_excess = data_excess%>%
#  filter(year>=2004)%>%
  group_by(day_of_year)%>%
  summarise(median = median(cummulative_excess),
            lower_quantile = quantile(cummulative_excess, 0.25),
            upper_quantile = quantile(cummulative_excess, 0.75))

data_excess%>%
  # filter(year%in%c(2006, 2008, 2010, 2012, 2019, 2018))%>%
  filter(year%in%2011:2012)%>%
  mutate(year = as.factor(year))%>%
  ggplot()+
  geom_line(aes(day_of_year, cummulative_excess, colour = year))+
  geom_ribbon(data = mean_data_excess, aes(day_of_year, ymin = lower_quantile, ymax = upper_quantile), alpha = 0.1)+
  geom_line(data = mean_data_excess, aes(day_of_year, median), colour='red', size = 1.5)+
  geom_vline(aes(xintercept = yday(dmy(paste0("30-09-", year))) - yday(dmy(paste0("01-03-", year)))+1))
ggsave("dryyears.png", dpi = 300, width = 6, height = 3)






