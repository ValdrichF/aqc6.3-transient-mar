# Run the iMODFLOW models
library(stringr)
library(tidyverse)
library(lubridate)
library(data.table)



setwd("E:/AMIGO_needed/RUNFILES")
imod = "D:/iMOD5/iMOD_V5_5.exe"
############
# Create a .ini file to translate the original.run to original.prj
ini_file_name = "prjtorun.ini"
ini_file = c(
  "function=runfile",
  "runfile_in=AMIGO_v31_NSTAT_25_wegschrijven_per_dag.run",
  "prjfile_out=Original.prj"
)
writeLines(ini_file, ini_file_name)
system(paste(imod, ini_file_name))


##########
# Create a .ini file to translate the original.prj to scenario.run
ini_file_name = "prjtorun.ini"
run_name = "avg_baseline"
start_date = "20111001"
end_date = "20121001"
time_interval = "4" # 3 - daily, 4 - weekly, 5 - 10 days, 6 - 14th and 28th, 7 - monthly

ini_file = c(
  "function=runfile",
  "prjfile_in=Original.prj",
  "runfile_out=run_name.run",
  "output_folder=K:\\MODFLOW_transient_res\\run_name",
  "sim_type=1",
  "iss=1",
  "sdate=start_date",
  "edate=end_date",
  "itt=time_interval",
  "idt=1",
  "saveshd=1,4,5,7",
  "savebnd=1,2,3,4",
  "savedrn=1",
  "saveolf=1",
  "saveriv=1",
  # "window=202275,437050,243450,462275",
  "cellsize=25",
  "ipest=1",
  "SSYSTEM=1"
)
ini_file = str_replace_all(ini_file, "run_name", run_name)
ini_file = str_replace_all(ini_file, "start_date", start_date)
ini_file = str_replace_all(ini_file, "end_date", end_date)
ini_file = str_replace_all(ini_file, "time_interval", time_interval)


writeLines(ini_file, ini_file_name)
system(paste(imod, ini_file_name))

########
# Correct the Run file (Dataset 2,3 and 5 from Original)
run_file = paste0(run_name, ".run")%>%
  readLines()
ds2 = run_file[2]%>%
  str_split("\\s+")%>%
  unlist
ds2_new = c(ds2[1:4], "0          1          0          1        -1")%>%
  paste0(collapse = "          ")

run_file_orj = readLines("AMIGO_v31_NSTAT_25_wegschrijven_per_dag.run", n = 10)
ds3 = run_file_orj[3:5]

new_run = c(run_file[1], 
            ds2_new,
            ds3,
            run_file[-(1:4)]
)

# Use the heads from CHD as as the initial conditions
## Find the previous date available among CHD
# CHDs_all = list.files("E:/AMIGO_needed/CHD/head", full.names = T) # initial
CHDs_all = list.files("E:/AMIGO_needed/SHD", full.names = T) # period
CHDs = CHDs_all%>%
  str_extract("^.*_l")%>%
  unique
# Formatting CHDs path names
CHDs = CHDs%>%
  str_replace_all("/","\\\\\\\\")
CHDs_dates = str_extract(CHDs, "\\d{8}")%>%
  ymd()
SHD_date = ymd(start_date) - CHDs_dates
SHD_date[SHD_date<0] = 9999 # ignore CHDs after the start date
CHD_last = CHDs[which(SHD_date==min(SHD_date, na.rm = T))]

SHD_ind = str_which(new_run, "\\\\SHD\\\\")
new_run[SHD_ind] = str_replace(new_run[SHD_ind], "E:.*_L", CHD_last)

# Save the MetaSWAP outputs
new_run[str_which(new_run, "^1,0\\s.CAP.$")] = "1,1,1 (CAP)"

## if weekly, use the mete_grid_weekly.inp in metaswap
if(F){#time_interval=="4"){
  new_run = new_run%>%
    str_replace("mete_grid", "mete_grid_weekly")
}

# Change the river inputs to use the river stage every month
riv_string = c(
  "  9,(RIV) Rivers",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\RIVIER\\Conductance_ISG_AMIGO.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\LEGGER\\conductance_legger_SEASON.IDF", # Season
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Conductance_TOP10_BREED_BUITEN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Conductance_TOP10_BREED_WRIJ.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Conductance_TOP10_NORMAAL_BUITEN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Conductance_TOP10_NORMAAL_WRIJ.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Conductance_TOP10_SMAL_BUITEN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Conductance_TOP10_SMAL_WRIJ.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TWENTE_KANAAL\\Conductance_TWKAN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\RIVIER\\MAAND\\STAGE_YEARMONTH01.IDF",# YEARMONTH
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\LEGGER\\peil_legger_SEASON.IDF", # SEASON
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Peil_TOP10_BREED_BUITEN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Peil_TOP10_BREED_WRIJ.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Peil_TOP10_NORMAAL_BUITEN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Peil_TOP10_NORMAAL_WRIJ.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Peil_TOP10_SMAL_BUITEN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Peil_TOP10_SMAL_WRIJ.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TWENTE_KANAAL\\Peil_TWKAN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\RIVIER\\Bodemhoogte_ISG_AMIGO.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\LEGGER\\bodemhoogte_legger_SEASON.IDF", # SEASON
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Bodemhoogte_TOP10_BREED_BUITEN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Bodemhoogte_TOP10_BREED_WRIJ.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Bodemhoogte_TOP10_NORMAAL_BUITEN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Bodemhoogte_TOP10_NORMAAL_WRIJ.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Bodemhoogte_TOP10_SMAL_BUITEN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Bodemhoogte_TOP10_SMAL_WRIJ.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TWENTE_KANAAL\\Bodemhoogte_TWKAN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\RIVIER\\Infiltratiefactor_ISG_AMIGO.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\LEGGER\\infiltratiefactor_legger_SEASON.IDF", # SEASON
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Infiltratiefactor_TOP10_BREED_BUITEN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Infiltratiefactor_TOP10_BREED_WRIJ.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Infiltratiefactor_TOP10_NORMAAL_BUITEN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Infiltratiefactor_TOP10_NORMAAL_WRIJ_SEASON.IDF", # SEASON
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Infiltratiefactor_TOP10_SMAL_BUITEN.IDF",
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TOP10\\Infiltratiefactor_TOP10_SMAL_WRIJ_SEASON.IDF", # SEASON
  "     1,1.000,0.0,E:\\AMIGO_needed\\OPPERVLAKTEWATER\\TWENTE_KANAAL\\Infiltratiefactor_TWKAN.IDF"
)

time_index = str_which(new_run, "^\\d{5},")
tim = mapply(
  function(i_start, i_end){
  txt = new_run[i_start:i_end]
  
  riv_index = str_which(txt, "Rivers")
  n_dat = str_extract(txt[riv_index], "-?\\d+")%>%
    as.numeric()
  riv_end_index = ifelse(n_dat==-1, riv_index, riv_index+n_dat*4)
  
  start_date = str_extract(txt[1], "\\d{8}")%>%
    str_sub(1, 6)
  season = ifelse(start_date%>%ym()%>%month()%in%4:9, "zomer", "WINTER")
  
  c(txt[1:(riv_index-1)],
    riv_string%>%
      str_replace("SEASON", season)%>%
      str_replace("YEARMONTH", start_date),
    txt[(riv_end_index+1):length(txt)]
    )
},
time_index,
c(time_index[-1]-1, length(new_run)))%>%
  unlist()
new_run = c(
  new_run[1:(time_index[1]-1)],
  tim
)

# Change the WEL data due to cropping
WEL_string = c(
  " 33,(WEL) Wells",
  "     1,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_1_Bemaling.ipf'",
  "     1,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_1_Landbouw.ipf'",
  "     1,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_1_Onbekend.ipf'",
  "     1,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_1_Overig.ipf'",
  "     2,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_2_Bemaling.ipf'",
  "     2,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_2_Landbouw.ipf'",
  "     2,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_2_Overig.ipf'",
  "     3,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_3_Bemaling.ipf'",
  "     3,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_3_Drinkwaterwinning.ipf'",
  "     3,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_3_Landbouw.ipf'",
  "     3,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_3_Overig.ipf'",
  "     4,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_4_Bemaling.ipf'",
  "     4,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_4_Drinkwaterwinning.ipf'",
  "     4,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_4_Landbouw.ipf'",
  "     4,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_4_Onbekend.ipf'",
  "     4,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_4_Overig.ipf'",
  "     5,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_5_Landbouw.ipf'",
  "     5,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_5_Overig.ipf'",
  "     6,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_6_Bemaling.ipf'",
  "     6,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_6_Drinkwaterwinning.ipf'",
  "     6,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_6_Landbouw.ipf'",
  "     6,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_6_Onbekend.ipf'",
  "     6,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_6_Overig.ipf'",
  "     8,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_8_Bemaling.ipf'",
  "     8,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_8_Drinkwaterwinning.ipf'",
  "     8,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_8_Landbouw.ipf'",
  "     9,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_9_Bemaling.ipf'",
  "     9,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_9_Drinkwaterwinning.ipf'",
  "     9,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_9_Overig.ipf'",
  "    10,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_10_Drinkwaterwinning.ipf'",
  "    11,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_11_Drinkwaterwinning.ipf'",
  "    14,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_14_Overig.ipf'",
  "    14,   1.000000    ,   0.000000    ,'E:\\AMIGO_needed\\ONTTREKKINGEN\\NIET_STATIONAIR_crop\\Locatie_onttrekkingen_AMIGO_14_Landbouw.ipf'"
)

time_index = str_which(new_run, "^\\d{5},")
tim = mapply(
  function(i_start, i_end){
    txt = new_run[i_start:i_end]
    
    wel_index = str_which(txt, "Wells")
    n_dat = str_extract(txt[wel_index], "-?\\d+")%>%
      as.numeric()
    paste(new_run[i_start], n_dat)%>%
      print()
    
    wel_end_index = wel_index+n_dat
    
    start_d = str_extract(new_run[i_start], "\\d{8}")%>%
      ymd()
    
    if((day(start_d)<=7)|(start_d==ymd(start_date))){
      return(
        c(txt[1:(wel_index-1)],
          WEL_string,
          txt[(wel_end_index+1):length(txt)]
        )
      )
    }
    c(txt[1:(wel_index-1)],
      " -1,(WEL) Wells",
      txt[(wel_end_index+1):length(txt)]
    )
  },
  time_index,
  c(time_index[-1]-1, length(new_run)))%>%
  unlist()
new_run = c(
  new_run[1:(time_index[1]-1)],
  tim
)

# Remove repeated files from RUN
previous_size = length(new_run)
time_index = str_which(new_run, "^\\d{5},")%>%
  c(., length(new_run)+1)%>%
  data.table(index = .)
time_index[,":="(
  second_index = c(index[-1], NA),
  third_index = c(index[-(1:2)], NA, NA)
)]
time_index = time_index[complete.cases(time_index)][-1]
tim = mapply(
  function(i_start, i_mid, i_end){
    txt1 = new_run[i_start:(i_mid-1)]
    txt2 = new_run[i_mid:(i_end-1)]
    
    i_packages1 = str_which(txt1, "\\([[:upper:]]{3}\\)")
    i_packages2 = str_which(txt2, "\\([[:upper:]]{3}\\)")
    check_packages2 = str_which(txt2[i_packages2], "^\\s-1", negate = TRUE)
    
    # If all packages reuse previous values, nothing to be changed
    if(length(check_packages2)==0) return(txt2)
    
    is_repeated = rep(F, length(i_packages2))
    
    
    # If the same package is used in txt1, check if they're the same
    for (i_pattern in check_packages2){
      same_package = which(txt1[i_packages1]==txt2[i_packages2[i_pattern]])

      if (length(same_package)==0) next

      package_1 = txt1[i_packages1[same_package]:(c(i_packages1[-1]-1, length(txt1))[same_package])]
      package_2 = txt2[i_packages2[same_package]:(c(i_packages2[-1]-1, length(txt2))[same_package])]
      
      if (identical(package_1, package_2)){
        is_repeated[same_package] = T
      }
    }
    if (!any(is_repeated)) return(txt2)
    # change the number of systems to -1
    txt2[i_packages2[is_repeated]] = txt2[i_packages2[is_repeated]]%>%
      str_replace("^\\s+\\d+", " -1")
    # Remove the data
    i_remove = mapply(seq, 
                      i_packages2[is_repeated]+1,
                      c(i_packages2[-1],length(txt2))[is_repeated]-1)%>%
      unlist
    txt2[i_remove*-1]
  },
  time_index[, index],
  time_index[, second_index],
  time_index[, third_index]
  )%>%
  unlist()

new_run = c(
  new_run[1:(time_index[1, second_index]-1)],
  tim
)
print(previous_size-length(new_run))

writeLines(new_run, paste0(run_name, ".run"))

# Also make the required changes in para_sim.inp
para_sim_dat = readLines("E:/AMIGO_needed/METASWAP/inp/para_sim.inp")
IDBG = start_date%>%
  ymd()%>%
  yday()%>%
  as.character
IYBG = start_date%>%
  ymd()%>%
  year()%>%
  as.character

para_sim_dat = para_sim_dat%>%
  str_replace("(IDBG[:space:]+=[:space:]+)\\d+(.+)$",
              paste0("\\1", IDBG, "\\2"))%>%
  str_replace("(IYBG[:space:]+=[:space:]+)\\d+$",
              paste0("\\1", IYBG))

# Also change the dtgw and dtsw (MODFLOW time-step and fast-iteration timestep)
if(time_interval=="4"){
  para_sim_dat = para_sim_dat%>%
    str_replace("(\\sdtgw\\s.*)\\d\\.", "\\17\\.")
}else{
  para_sim_dat = para_sim_dat%>%
    str_replace("(\\sdtgw\\s.*)\\d\\.", "\\11\\.")
}

para_sim_dat%>%
  writeLines(con = "E:/AMIGO_needed/METASWAP/inp/para_sim.inp")

# Run the model
paste0(dirname(imod), "/bin/iMODFLOW-METASWAP_V5_5.exe ",
# paste0(imod," ",
       run_name, ".run ",
       "> ",
       run_name, ".log ")%>%
  writeLines(paste0(run_name, ".bat "))
