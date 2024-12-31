# create the mete_grid.inp to run the scenario as a comma separated file
# Columns are: 
#   day of year `round(., 2)`
#   Year `int()`
#   RCH - in mm/day
#   ET - in mm/day
#   tempmn grid "NoValue"
#   tempmx grid "NoVakue"
#   temp grid "NoValue"


library(lubridate)
library(data.table)

inp_file = data.table(
  time = seq(ymd("2011-10-01"), ymd("2011-12-03"), by = 1),
  RCH = "Z:\\Sy_scenario\\RCH\\RCH_L1.ASC",
  ET = "Z:\\Sy_scenario\\ET\\ET_L1.ASC",
  tempmn = "NoValue",
  tempmx = "NoValue",
  temp_grid = "NoValue"
)
inp_file[,":="(
  doy = yday(time),
  year = year(time),
  time = NULL
)]

setcolorder(inp_file, c("doy", "year", "RCH", "ET", "tempmn", "tempmx", "temp_grid"))

fwrite(inp_file, "Z:\\Sy_scenario\\METASWAP\\inp\\mete_grid.inp",
       quote = T, col.names = F)





