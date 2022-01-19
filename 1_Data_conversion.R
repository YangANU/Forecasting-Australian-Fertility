##############################################
# Conversion of raw data to demography format
#############################################

library(tidyverse)
library(readxl)
library(demography)

## Importing deaths data
Deaths47 = read_excel("./Data/Deaths_47areas_1993_2016 singe year.xlsx", sheet = 1, col_names = T) %>%
  arrange(year, cob, NSD, sex, age)

# Renaming the variables for coding
AgeFct = factor(unique(Deaths47$age), levels = unique(Deaths47$age))
Age_mid_point = c(seq(2, 83, 5), 85)
Deaths47 = mutate(Deaths47, cob = ifelse(cob == "Aus", "AusBorn", "OvsBorn"), 
                  region = paste("R", ifelse(region < 10, paste("0", region, sep = ""), region), sep = ""), 
                  NSD = paste("A", NSD, sep = ""), sex = ifelse(sex == "F", "FEM", "MLE"), 
                  age = factor(age, levels = AgeFct, labels = Age_mid_point))

# Setting global variables
Year = unique(Deaths47$year)
Region = unique(sort(Deaths47$region))
Sex = unique(Deaths47$sex)


# Importing population data
ERP47 = read_excel("./Data/ERP_47areas_1981_2016 single year.xlsx", sheet = 1, col_names = T) %>% 
  filter(year >= Year[1] & year <= Year[length(Year)]) %>% arrange(year, cob, NSD, sex, age) %>%
  mutate(cob = ifelse(cob == "Aus", "AusBorn", "OvsBorn"), 
         region = paste("R", ifelse(region < 10, paste("0", region, sep = ""), region), sep = ""), 
         NSD = paste("A", NSD, sep = ""), sex = ifelse(sex == "F", "FEM", "MLE"), 
         age = factor(age, levels = AgeFct, labels = Age_mid_point))


# Importing NSD_47_code reference table

NSD47 = read_excel("./Data/Geography look-up table.xlsx", sheet = 1, col_names = T) %>%
  arrange(Area, new_region, NSD_47_code, NSD_47_name, new_region_name) %>%
  mutate(new_region = paste("R", ifelse(new_region < 10, paste("0", new_region, sep = ""), new_region), sep = ""), Area = paste("A", ifelse(Area < 10, paste("0", Area, sep = ""), Area), sep = ""))

NSD = unique(Deaths47$NSD)
Area = c("Australia", unique(NSD47$Area))

###############################################
# Collection of subnational areas by Area code
###############################################

# Defining a function to export data and generate demogdata

ConvertDemog_sex <- function(tempDx, tempPop, area)
{
  # Calculate Death rates: Mx = Deaths/Exposure-to-risk
  tempMx = bind_cols(tempDx, tempPop) %>% mutate(Female = Female/Female1, Male = Male/Male1, Total = Total/Total1) %>% select(year, age, Female, Male, Total)
  
  # Export data and generate demogdata 
  datMx = filter(tempMx) %>% select(year, age, Female, Male, Total)
  datPop = filter(tempPop) %>% select(year, age, Female, Male, Total)
  file_Mx = paste("./Data/Australian_demogdata/", area, "_rate.txt", sep = "")
  write.table(datMx, file_Mx, quote = FALSE, row.names = FALSE, col.names = TRUE)
  file_Pop = paste("./Data/Australian_demogdata/", area, "_count.txt", sep = "")
  write.table(datPop, file_Pop, quote = FALSE, row.names = FALSE, col.names = TRUE)
  area_name = filter(NSD47, NSD_47_code == sub('.', '', area)) %>% select(NSD_47_name)
  area_code = filter(NSD47, NSD_47_code == sub('.', '', area)) %>% select(Area)
  if(area == "Australia")
  {
    area_name = "Australia"
    area_code = "Australia"
  }
  assign(as.character(area_code), read.demogdata(file_Mx, file_Pop, type = "mortality", label = area_name, skip = 0, popskip = 0), pos = 1)
}

ConvertDemog_sex_region <- function(tempDx, tempPop, region)
{
  # Calculate Death rates: Mx = Deaths/Exposure-to-risk
  tempMx = bind_cols(tempDx, tempPop) %>% mutate(Female = Female/Female1, Male = Male/Male1, Total = Total/Total1) %>% select(year, age, Female, Male, Total)
  
  # Export data and generate demogdata 
  datMx = filter(tempMx) %>% select(year, age, Female, Male, Total)
  datPop = filter(tempPop) %>% select(year, age, Female, Male, Total)
  file_Mx = paste("./Data/Australian_demogdata/", region, "_rate.txt", sep = "")
  write.table(datMx, file_Mx, quote = FALSE, row.names = FALSE, col.names = TRUE)
  file_Pop = paste("./Data/Australian_demogdata/", region, "_count.txt", sep = "")
  write.table(datPop, file_Pop, quote = FALSE, row.names = FALSE, col.names = TRUE)
  region_name = filter(NSD47, new_region == region) %>% select(new_region_name)
  assign(region, read.demogdata(file_Mx, file_Pop, type = "mortality", label = region_name, skip = 0, popskip = 0), pos = 1)
}


# Obtaining subnational mortality rates for Australia (11 regions containing 47 areas in total)

for (reg in Region)
{
  # Obtain sub-dataset for the region
  tempPop_reg = filter(ERP47, region == reg)
  tempDx_reg = filter(Deaths47, region == reg)
  
  for (area in unique(tempPop_reg$NSD))
  {
    # Calculate Mx and Pop for the sub-regional area(NSD) by birthplace
    tempPop = filter(tempPop_reg, NSD == area)
    tempPop = group_by(tempPop, year, region, sex, age) %>% summarise(ERP = sum(ERP)) %>%
      ungroup() %>% spread(sex, ERP) %>% mutate(Total = FEM + MLE) %>% 
      select(year, age, Female = FEM, Male = MLE, Total)
    
    tempDx = filter(tempDx_reg, NSD == area) 
    tempDx = group_by(tempDx, year, region, sex, age) %>% summarise(deaths = sum(deaths)) %>%
      ungroup() %>% spread(sex, deaths) %>% mutate(Total = FEM + MLE) %>% 
      select(year, age, Female = FEM, Male = MLE, Total)

    # Export data and generate demogdata by birthplace and aggregated across birthplaces
    ConvertDemog_sex(tempDx, tempPop, area)
  
  }
  tempPop_reg = group_by(tempPop_reg, year, region, sex, age) %>% summarise(ERP = sum(ERP)) %>%
    ungroup() %>% spread(sex, ERP) %>% mutate(Total = FEM + MLE) %>% 
    select(year, age, Female = FEM, Male = MLE, Total)
  tempDx_reg = group_by(tempDx_reg, year, region, sex, age) %>% summarise(deaths = sum(deaths)) %>%
    ungroup() %>% spread(sex, deaths) %>% mutate(Total = FEM + MLE) %>% 
    select(year, age, Female = FEM, Male = MLE, Total)
  ConvertDemog_sex_region(tempDx_reg, tempPop_reg, reg)
}

# Aggregated across geography to obtain Australian totals
tempPop_Aus = group_by(ERP47, year, sex, age) %>% summarise(ERP = sum(ERP)) %>%
  ungroup() %>% spread(sex, ERP) %>% mutate(Total = FEM + MLE) %>%
  select(year, age, Female = FEM, Male = MLE, Total)
tempDx_Aus = group_by(Deaths47, year, sex, age) %>% summarise(deaths = sum(deaths)) %>%
  ungroup() %>% spread(sex, deaths) %>% mutate(Total = FEM + MLE) %>%
  select(year, age, Female = FEM, Male = MLE, Total)
ConvertDemog_sex(tempDx_Aus, tempPop_Aus, "Australia")  

# Clean up workspace
rm(tempDx, tempDx_Aus, tempDx_reg, tempPop, tempPop_Aus, tempPop_reg, area, file_Mx, file_Pop, reg, dat, dat_bp)

###################################################################################
# Combine Area level total series (unsmoothed) into a demography time series file
###################################################################################

total_comb = total_comb_pop =  matrix(NA, 18*24, 47)
for(iw in 2:48)
{
  total_comb[,iw-1] = as.numeric(get(Area[iw])$rate$total)
  total_comb_pop[,iw-1] = as.numeric(get(Area[iw])$pop$total)
}
total_comb_v2 = cbind(rep(1993:2016, each=18), rep(seq(0,85,5), 24), total_comb)
total_comb_pop_v2 = cbind(rep(1993:2016, each=18), rep(seq(0,85,5), 24), total_comb_pop)
colnames(total_comb_v2) = colnames(total_comb_pop_v2) = c("Year", "Age", Area[2:48])

temp = 0
for(i in 2:48)
{
  temp = temp + sum(get(Area[i])$rate$total == 0)
}
temp == sum(total_comb == 0)
###################################################################################
# Combine Region level total series (unsmoothed) into a demography time series file
###################################################################################

total_region_comb = total_region_comb_pop = matrix(NA, 18*24, 11)
for(iw in 1:11)
{
  total_region_comb[,iw] = as.numeric(get(Region[iw])$rate$total)
  total_region_comb_pop[,iw] = as.numeric(get(Region[iw])$pop$total)
}
total_region_comb_v2 = cbind(rep(1993:2016, each=18), rep(seq(0,85,5), 24), total_region_comb)
total_region_comb_pop_v2 = cbind(rep(1993:2016, each=18), rep(seq(0,85,5), 24), total_region_comb_pop)
colnames(total_region_comb_v2) = colnames(total_region_comb_pop_v2) = c("Year", "Age", Region)

temp = 0
for(i in 1:11)
{
  temp = temp + sum(get(Region[i])$rate$total == 0)
}
temp == sum(total_region_comb == 0)


##########################################################################################
# smoothed functional curves using penalized regression spline with monotonic constraint
##########################################################################################

## Smoothing Area data
# check the number of missing values
area_missing_value = list()
for(i in 1:48)
{
  area_missing_value[[i]] = lapply(get(Area[i])$rate, function(c)sum(c==0))
}

missing_female = which(sapply(area_missing_value, function(c){c$female > length(Australia$rate$female) * 0.15}))
missing_male = which(sapply(area_missing_value, function(c){c$male > length(Australia$rate$female) * 0.15}))
missing_total = which(sapply(area_missing_value, function(c){c$total > length(Australia$rate$female) * 0.15}))

missing_index = unique(c(missing_female, missing_male, missing_total))

Area_smooth = paste(Area, "smooth", sep = "_")
Area_index = 1:48
smoothing_index = Area_index[! Area_index %in% missing_index]

for(iw in smoothing_index)
{
  assign(Area_smooth[iw], smooth.demogdata(get(Area[iw]), age.grid = 0:85))
}


for(i in missing_index)
{
  demogdata_temp = get(Area[i])
  
  if(i %in% missing_female)
  {
    missing_female_ind = (demogdata_temp$rate$female == 0)
    demogdata_temp$rate$female = replace(demogdata_temp$rate$female, missing_female_ind, get(Area[1])$rate$female[missing_female_ind])
  }
  
  if(i %in% missing_male)
  {
    missing_male_ind = (demogdata_temp$rate$male == 0)
    demogdata_temp$rate$male = replace(demogdata_temp$rate$male, missing_male_ind, get(Area[1])$rate$male[missing_male_ind])
  }
 
  if(i %in% missing_total)
  {
    missing_total_ind = (demogdata_temp$rate$total == 0)
    demogdata_temp$rate$total = replace(demogdata_temp$rate$total, missing_total_ind, get(Area[1])$rate$total[missing_total_ind])
  }
  
  assign(Area_smooth[i], smooth.demogdata(demogdata_temp, age.grid = 0:85))
}

## Smoothing Region data

Region_smooth = paste(Region, "smooth", sep = "_")
for(iw in 1:11)
{
  assign(Region_smooth[iw], smooth.demogdata(get(Region[iw]), age.grid = 0:85))
}


##################################################################################
# Collection of subnational areas by NSD code; splitting data across birth places
##################################################################################

## Not run:

# Define a function to export data and generate demogdata
ConvertDemog = function(tempDx, tempPop, area){
  # Calculate Mx by birthplace: Mx = Dx / Pop
  tempMx = bind_cols(tempDx, tempPop) %>% mutate(Female = Female/Female1, Male = Male/Male1, Total = Total/Total1) %>%
    select(year, cob, age, Female, Male, Total)
  # Export data and generate demogdata for each birthplace
  for (bp in unique(tempPop$cob)){
    datMx = filter(tempMx, cob == bp) %>% select(year, age, Female, Male, Total)
    datPop = filter(tempPop, cob ==bp) %>% select(year, age, Female, Male, Total)
    file_Mx = paste("./Data/Australian_demogdata/Birth_place", area, "_", bp, "_rate.txt", sep = "")
    write.table(datMx, file_Mx, quote = FALSE, row.names = FALSE, col.names = TRUE)
    file_Pop = paste("./Data/Australian_demogdata/Birth_place", area, "_", bp, "_count.txt", sep = "")
    write.table(datPop, file_Pop, quote = FALSE, row.names = FALSE, col.names = TRUE)
    dat = paste(area, bp, sep = "_")
    assign(dat, read.demogdata(file_Mx, file_Pop, type = "mortality", label = paste(area, bp, sep = "_"), skip = 0, popskip = 0), pos = 1)
  }
  # Calculate Mx and Pop aggregated across birthplaces
  tempPop = group_by(tempPop, year, age) %>% 
    summarise(Female = sum(Female), Male = sum(Male), Total = sum(Total)) %>% ungroup()
  tempDx = group_by(tempDx, year, age) %>% 
    summarise(Female = sum(Female), Male = sum(Male), Total = sum(Total)) %>% ungroup()
  tempMx = bind_cols(tempDx, tempPop) %>% mutate(Female = Female/Female1, Male = Male/Male1, Total = Total/Total1) %>%
    select(year, age, Female, Male, Total)
  tempPop = select(tempPop, year, age, Female, Male, Total)
  # Export data and generate demogdata for the area aggerated across birthplaces
  file_Mx = paste("./Data/Australian_demogdata/Birth_place", area, "_rate.txt", sep = "")
  write.table(tempMx, file_Mx, quote = FALSE, row.names = FALSE, col.names = TRUE)
  file_Pop = paste("./Data/Australian_demogdata/Birth_place", area, "_count.txt", sep = "")
  write.table(tempPop, file_Pop, quote = FALSE, row.names = FALSE, col.names = TRUE)
  dat = area
  assign(dat, read.demogdata(file_Mx, file_Pop, type = "mortality", label = dat, skip = 0, popskip = 0), pos = 1)
}


# Obtain sub-dataset by looping across the geography
for (reg in Region){
  # Obtain sub-dataset for the region
  tempPop_reg = filter(ERP47, region == reg)
  tempDx_reg = filter(Deaths47, region == reg)
  for (area in unique(tempPop_reg$NSD)){
    # Calculate Mx and Pop for the sub-regional area(NSD) by birthplace
    tempPop = filter(tempPop_reg, NSD == area) %>% spread(sex, ERP) %>% mutate(Total = FEM + MLE) %>%
      select(year, cob, age, Female = FEM, Male = MLE, Total)
    tempDx = filter(tempDx_reg, NSD == area) %>% spread(sex, deaths) %>% mutate(Total = FEM + MLE) %>%
      select(year, cob, age, Female = FEM, Male = MLE, Total)
    # Export data and generate demogdata by birthplace and aggregated across birthplaces
    ConvertDemog(tempDx, tempPop, area)
  }
  # Aggregated across areas within the region
  tempPop_reg = group_by(tempPop_reg, year, cob, region, sex, age) %>% summarise(ERP = sum(ERP)) %>%
    ungroup() %>% spread(sex, ERP) %>% mutate(Total = FEM + MLE) %>% 
    select(year, cob, age, Female = FEM, Male = MLE, Total)
  tempDx_reg = group_by(tempDx_reg, year, cob, region, sex, age) %>% summarise(deaths = sum(deaths)) %>%
    ungroup() %>% spread(sex, deaths) %>% mutate(Total = FEM + MLE) %>% 
    select(year, cob, age, Female = FEM, Male = MLE, Total)
  ConvertDemog(tempDx_reg, tempPop_reg, reg)
}
# Aggregated across geography to obtain Australian totals
tempPop_Aus = group_by(ERP47, year, cob, sex, age) %>% summarise(ERP = sum(ERP)) %>%
  ungroup() %>% spread(sex, ERP) %>% mutate(Total = FEM + MLE) %>%
  select(year, cob, age, Female = FEM, Male = MLE, Total)
tempDx_Aus = group_by(Deaths47, year, cob, sex, age) %>% summarise(deaths = sum(deaths)) %>%
  ungroup() %>% spread(sex, deaths) %>% mutate(Total = FEM + MLE) %>%
  select(year, cob, age, Female = FEM, Male = MLE, Total)
ConvertDemog(tempDx_Aus, tempPop_Aus, "Australia")  


# Clean up workspace
rm(tempDx, tempDx_Aus, tempDx_reg, tempPop, tempPop_Aus, tempPop_reg, area, bp, file_Mx, file_Pop, reg, dat, dat_bp)

## End(**Not run**)
























