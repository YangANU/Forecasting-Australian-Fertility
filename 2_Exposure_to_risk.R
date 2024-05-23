################################################
# Forecast exposure to risk at the bottom-level 
################################################


###############
# Fixing COB_j
###############

library(ftsa)

# year_horizon: forecast horizon
# top: exposure to risk for a given series at the area level
# bottom: exposure to risk for AUS total series

pop_fore_fix_COB <- function(top, bottom, year_horizon)
{
  pop_fore = array(0, dim = c(year_horizon, 7, year_horizon, 19))

  for(COB in 1:19)
  {
    for(year in 22:31)
    {
      for(j in 1:7)
      {
        ratio = log(top[1:year,j,COB]/bottom[1:year,j,COB] + 1e-6)
        
        pop_fore[,j,year-21,COB] =  exp(forecast(auto.arima(ratio), h = year_horizon)$mean) - 1e-6
      }
    }
  }

  return(pop_fore)
}  

## Area/AUS
pop_ratio_Area_to_AUS = c("pop_ratio_A1_to_AUS", "pop_ratio_A2_to_AUS", "pop_ratio_A3_to_AUS", 
                          "pop_ratio_A4_to_AUS", "pop_ratio_A5_to_AUS", "pop_ratio_A6_to_AUS",
                          "pop_ratio_A7_to_AUS", "pop_ratio_A8_to_AUS", "pop_ratio_A9_to_AUS", 
                          "pop_ratio_A10_to_AUS", "pop_ratio_A11_to_AUS", "pop_ratio_A12_to_AUS", 
                          "pop_ratio_A13_to_AUS", "pop_ratio_A14_to_AUS", "pop_ratio_A15_to_AUS",
                          "pop_ratio_A16_to_AUS", "pop_ratio_A17_to_AUS", "pop_ratio_A18_to_AUS", 
                          "pop_ratio_A19_to_AUS", "pop_ratio_A20_to_AUS", "pop_ratio_A21_to_AUS", 
                          "pop_ratio_A22_to_AUS", "pop_ratio_A23_to_AUS", "pop_ratio_A24_to_AUS",
                          "pop_ratio_A25_to_AUS", "pop_ratio_A26_to_AUS", "pop_ratio_A27_to_AUS", 
                          "pop_ratio_A28_to_AUS", "pop_ratio_A29_to_AUS", "pop_ratio_A30_to_AUS", 
                          "pop_ratio_A31_to_AUS", "pop_ratio_A32_to_AUS", "pop_ratio_A33_to_AUS",
                          "pop_ratio_A34_to_AUS", "pop_ratio_A35_to_AUS", "pop_ratio_A36_to_AUS", 
                          "pop_ratio_A37_to_AUS", "pop_ratio_A38_to_AUS", "pop_ratio_A39_to_AUS", 
                          "pop_ratio_A40_to_AUS", "pop_ratio_A41_to_AUS", "pop_ratio_A42_to_AUS",
                          "pop_ratio_A43_to_AUS", "pop_ratio_A44_to_AUS", "pop_ratio_A45_to_AUS", 
                          "pop_ratio_A46_to_AUS", "pop_ratio_A47_to_AUS")


for(iw in 1:47)
{
  assign(pop_ratio_Area_to_AUS[iw], pop_fore_fix_COB(top = pop_interp_area[[iw]], bottom = pop_interp_national, year_horizon = 10))
  print(iw)
}

for(iw in 1:47)
{
  print(sum(is.na(get(pop_ratio_Area_to_AUS[iw]))))
}

for(iw in 1:47)
{
  print(sum(get(pop_ratio_Area_to_AUS[iw])<0))
}

for(iw in 1:47)
{
  temp = get(pop_ratio_Area_to_AUS[iw])
  temp[get(pop_ratio_Area_to_AUS[iw])<=0] = 0
  assign(pop_ratio_Area_to_AUS[iw], temp)
  rm(temp)
}


## Area/Region
pop_ratio_Area_to_R1 = "pop_ratio_A1_to_R1"
pop_ratio_Area_to_R2 = c("pop_ratio_A2_to_R2", "pop_ratio_A3_to_R2", "pop_ratio_A4_to_R2")
pop_ratio_Area_to_R3 = "pop_ratio_A10_to_R3"
pop_ratio_Area_to_R4 = c("pop_ratio_A11_to_R4", "pop_ratio_A12_to_R5", "pop_ratio_A13_to_R5", "pop_ratio_A16_to_R4", "pop_ratio_A18_to_R4")
pop_ratio_Area_to_R5 = "pop_ratio_A19_to_R5"
pop_ratio_Area_to_R6 = "pop_ratio_A27_to_R6"
pop_ratio_Area_to_R7 = "pop_ratio_A33_to_R7"
pop_ratio_Area_to_R8 = "pop_ratio_A42_to_R8"
pop_ratio_Area_to_R9 = "pop_ratio_A47_to_R9"
pop_ratio_Area_to_R10 = c("pop_ratio_A5_to_R10", "pop_ratio_A6_to_R10", "pop_ratio_A7_to_R10", "pop_ratio_A8_to_R10", "pop_ratio_A14_to_R10", "pop_ratio_A15_to_R10", "pop_ratio_A17_to_R10", "pop_ratio_A20_to_R10", "pop_ratio_A21_to_R10", "pop_ratio_A24_to_R10", "pop_ratio_A28_to_R10", "pop_ratio_A29_to_R10", "pop_ratio_A30_to_R10", "pop_ratio_A34_to_R10", "pop_ratio_A35_to_R10", "pop_ratio_A36_to_R10", "pop_ratio_A37_to_R10", "pop_ratio_A43_to_R10", "pop_ratio_A45_to_R10")
pop_ratio_Area_to_R11 = c("pop_ratio_A9_to_R11", "pop_ratio_A22_to_R11", "pop_ratio_A23_to_R11", "pop_ratio_A25_to_R11", "pop_ratio_A26_to_R11", "pop_ratio_A31_to_R11", "pop_ratio_A32_to_R11", "pop_ratio_A38_to_R11", "pop_ratio_A39_to_R11", "pop_ratio_A40_to_R11", "pop_ratio_A41_to_R11", "pop_ratio_A44_to_R11", "pop_ratio_A46_to_R11")

pop_ratio_Area_to_Region = c(pop_ratio_Area_to_R1, pop_ratio_Area_to_R2, pop_ratio_Area_to_R3,
                             pop_ratio_Area_to_R4, pop_ratio_Area_to_R5, pop_ratio_Area_to_R6,
                             pop_ratio_Area_to_R7, pop_ratio_Area_to_R8, pop_ratio_Area_to_R9,
                             pop_ratio_Area_to_R10, pop_ratio_Area_to_R11)

# R1
assign(pop_ratio_Area_to_R1, pop_fore_fix_COB(top = pop_interp_area[[1]], bottom = pop_interp_region[[1]], year_horizon = 10))
# R2
for(iw in 2:4)
{
  assign(pop_ratio_Area_to_R2[[iw-1]], pop_fore_fix_COB(top = pop_interp_area[[iw]], bottom = pop_interp_region[[2]], year_horizon = 10))
}
# R3
assign(pop_ratio_Area_to_R3, pop_fore_fix_COB(top = pop_interp_area[[10]], bottom = pop_interp_region[[3]], year_horizon = 10))
# R4
for(iw in 1:5)
{
  j = region_list_ind[[4]][iw]
  assign(pop_ratio_Area_to_R4[[iw]], pop_fore_fix_COB(top = pop_interp_area[[j]], bottom = pop_interp_region[[4]], year_horizon = 10))
}
# R5
assign(pop_ratio_Area_to_R5, pop_fore_fix_COB(top = pop_interp_area[[19]], bottom = pop_interp_region[[5]], year_horizon = 10))
# R6
assign(pop_ratio_Area_to_R6, pop_fore_fix_COB(top = pop_interp_area[[27]], bottom = pop_interp_region[[6]], year_horizon = 10))
# R7
assign(pop_ratio_Area_to_R7, pop_fore_fix_COB(top = pop_interp_area[[33]], bottom = pop_interp_region[[7]], year_horizon = 10))
# R8
assign(pop_ratio_Area_to_R8, pop_fore_fix_COB(top = pop_interp_area[[42]], bottom = pop_interp_region[[8]], year_horizon = 10))
# R9
assign(pop_ratio_Area_to_R9, pop_fore_fix_COB(top = pop_interp_area[[47]], bottom = pop_interp_region[[9]], year_horizon = 10))
# R10
for(iw in 1:19)
{
  j = region_list_ind[[10]][iw]
  assign(pop_ratio_Area_to_R10[[iw]], pop_fore_fix_COB(top = pop_interp_area[[j]], bottom = pop_interp_region[[10]], year_horizon = 10))
}
# R11
for(iw in 1:13)
{
  j = region_list_ind[[11]][iw]
  assign(pop_ratio_Area_to_R11[[iw]], pop_fore_fix_COB(top = pop_interp_area[[j]], bottom = pop_interp_region[[11]], year_horizon = 10))
}

for(iw in 1:length(pop_ratio_Area_to_Region))
{
  print(sum(get(pop_ratio_Area_to_Region[iw])<0))
}

for(iw in 1:47)
{
  temp = get(pop_ratio_Area_to_Region[iw])
  temp[get(pop_ratio_Area_to_Region[iw])<=0] = 0
  assign(pop_ratio_Area_to_Region[iw], temp)
  rm(temp)
}



###########################################
# Fixing Area/Region (geographical factor)
###########################################

# year_horizon: forecast horizon
# top: exposure to risk for a given COB series level
# bottom: exposure to risk for COB group total series

pop_fore_fix_geo <- function(top, bottom, year_horizon)
{
  pop_fore = array(0, dim = c(year_horizon, 7, year_horizon, 59))
  
  for(geo in 1:59)
  {
    for(year in 22:31)
    {
      for(j in 1:7)
      {
        ratio = log(top[1:year,j,geo]/bottom[1:year,j,geo] + 1e-6)
        
        pop_fore[,j,year-21,geo] =  exp(forecast(auto.arima(ratio), h = year_horizon)$mean) - 1e-6
      }
    }
  }
  
  
  return(pop_fore)
}


## COB/national
pop_ratio_cob_to_national = c("pop_ratio_cob1_to_national", "pop_ratio_cob2_to_national", 
                              "pop_ratio_cob3_to_national", "pop_ratio_cob4_to_national", 
                              "pop_ratio_cob5_to_national", "pop_ratio_cob6_to_national", 
                              "pop_ratio_cob7_to_national", "pop_ratio_cob8_to_national", 
                              "pop_ratio_cob9_to_national", "pop_ratio_cob10_to_national",
                              "pop_ratio_cob11_to_national", "pop_ratio_cob12_to_national", 
                              "pop_ratio_cob13_to_national", "pop_ratio_cob14_to_national", 
                              "pop_ratio_cob15_to_national", "pop_ratio_cob16_to_national", 
                              "pop_ratio_cob17_to_national", "pop_ratio_cob18_to_national",
                              "pop_ratio_cob19_to_national")



for(iw in 1:19)
{
  assign(pop_ratio_cob_to_national[iw], pop_fore_fix_geo(top = get(paste("pop_interp_cob", iw, sep = "")), bottom = pop_interp_all_cob, year_horizon = 10))
  print(iw)
}


## COB/C_group
pop_ratio_cob_to_C1 = c("pop_ratio_cob1_to_C1", "pop_ratio_cob2_to_C1", "pop_ratio_cob3_to_C1")
pop_ratio_cob_to_C2 = c("pop_ratio_cob4_to_C2", "pop_ratio_cob5_to_C2")
pop_ratio_cob_to_C3 = "pop_ratio_cob6_to_C3"
pop_ratio_cob_to_C4 = "pop_ratio_cob7_to_C4"
pop_ratio_cob_to_C5 = c("pop_ratio_cob8_to_C5", "pop_ratio_cob9_to_C5", "pop_ratio_cob10_to_C5", "pop_ratio_cob11_to_C5", "pop_ratio_cob12_to_C5")
pop_ratio_cob_to_C6 = c("pop_ratio_cob13_to_C6", "pop_ratio_cob14_to_C6")
pop_ratio_cob_to_C7 = c("pop_ratio_cob15_to_C7", "pop_ratio_cob16_to_C7")
pop_ratio_cob_to_C8 = "pop_ratio_cob17_to_C8"
pop_ratio_cob_to_C9 = "pop_ratio_cob18_to_C9"
pop_ratio_cob_to_C10 = "pop_ratio_cob19_to_C10"

# C1
for(iw in 1:3)
{
  j = COB_list[[1]][iw]
  assign(pop_ratio_cob_to_C1[[iw]], pop_fore_fix_geo(top = get(paste("pop_interp_cob", j, sep="")), bottom = pop_interp_C1, year_horizon = 10))
}

# C2
for(iw in 1:2)
{
  j = COB_list[[1]][iw]
  assign(pop_ratio_cob_to_C2[[iw]], pop_fore_fix_geo(top = get(paste("pop_interp_cob", j, sep="")), bottom = pop_interp_C2, year_horizon = 10))
}

# C3
assign(pop_ratio_cob_to_C3, pop_fore_fix_geo(top = pop_interp_cob6, bottom = pop_interp_C3, year_horizon = 10))

# C4
assign(pop_ratio_cob_to_C4, pop_fore_fix_geo(top = pop_interp_cob7, bottom = pop_interp_C4, year_horizon = 10))

# C5
for(iw in 1:5)
{
  j = COB_list[[5]][iw]
  assign(pop_ratio_cob_to_C5[[iw]], pop_fore_fix_geo(top = get(paste("pop_interp_cob", j, sep="")), bottom = pop_interp_C5, year_horizon = 10))
}

# C6
for(iw in 1:2)
{
  j = COB_list[[6]][iw]
  assign(pop_ratio_cob_to_C6[[iw]], pop_fore_fix_geo(top = get(paste("pop_interp_cob", j, sep="")), bottom = pop_interp_C6, year_horizon = 10))
}

# C7
for(iw in 1:2)
{
  j = COB_list[[7]][iw]
  assign(pop_ratio_cob_to_C7[[iw]], pop_fore_fix_geo(top = get(paste("pop_interp_cob", j, sep="")), bottom = pop_interp_C7, year_horizon = 10))
}

# C8
assign(pop_ratio_cob_to_C8, pop_fore_fix_geo(top = pop_interp_cob17, bottom = pop_interp_C8, year_horizon = 10))

# C9
assign(pop_ratio_cob_to_C9, pop_fore_fix_geo(top = pop_interp_cob18, bottom = pop_interp_C9, year_horizon = 10))

# C10
assign(pop_ratio_cob_to_C10, pop_fore_fix_geo(top = pop_interp_cob19, bottom = pop_interp_C10, year_horizon = 10))
