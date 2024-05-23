####################################
# Australian sub-national fertility 
####################################

library(tidyverse)
library(demography)
library(Hmisc)
library(forecast)

fertility_all = read.csv(file = "ASFR_11region_1981_2016.csv", header = TRUE)
population_all = read.csv(file = "ERP_1981_2016.csv", header = TRUE)
colnames(population_all) = c("year", "cob", "NSD", "sex", "age", "ERP")

region_nsd = read.csv(file = "Geography.csv", header = TRUE)
colnames(region_nsd) = c("NSD", "NSD name", "Region", "Region name")
region_nsd = arrange(region_nsd, Region)

births_region_obs = read.csv(file = "births_11_region.csv")
births_area_obs = read.csv(file = "births_47_area.csv")

region_list = list()
for(i in 1:length(unique(region_nsd$Region)))
{
  region_ind = unique(region_nsd$Region)[i]
  
  region_list[[i]] = as.numeric(unlist(filter(region_nsd, Region %in% region_ind) %>% select(NSD)))
}

region_list_ind = list()
for(i in 1:length(region_list))
{
  region_list_ind[[i]] = match(region_list[[i]], sort(unique(region_nsd$NSD)))
}

for(i in 1:11)
{
  print(paste("region ", i))
  print(match(region_list[[i]], nsd_all))
}


COB = c("Australia", "New Zealand", "Other Oceania", "United Kindom", "Other Europe", "South-East Europe", "Middle East", "Vietnam", "Philippines", "Malaysia", "Indonesia", "Other South-East Asia", "China", "Other North-East Asia", "India", "Other Southern Asia", "North America", "South America", "Sub-Saharan Africa")


###########################
# Extract number of births
###########################

region_birth = list()
for(i in 1:11)
{
  region_birth[[i]] = array(0, dim = c(7,7,5,19)) # dimension = year, age, year_select, cob
  for(j in 1:7)
  {
    year_select = (year_all[j]:(year_all[j]+4))
    
    for(k in 1:19)
    {
      # sum up briths in 5 consecutive years
      births_select = filter(births_region_obs, year %in% year_select) %>% filter(region == i) %>% filter(cob == k) %>% select(age, births)
      region_birth[[i]][j,,,k] = matrix(births_select$births, nrow = 7, ncol = 5, byrow = FALSE)
      rownames(region_birth[[i]][j,,,k]) = unique(births_select$age)
      colnames(region_birth[[i]][j,,,k]) = year_select
    }
  }
}

area_birth = list()
for(i in 1:47)
{
  area_birth[[i]] = array(0, dim = c(7,7,5,19)) # dimension = year, age, year, cob
  for(j in 1:7)
  {
    year_select = (year_all[j]:(year_all[j]+4))
    
    for(k in 1:19)
    {
      # sum up briths in 5 consecutive years
      births_select = filter(births_area_obs, year %in% year_select) %>% filter(NSD == nsd_all[i]) %>% arrange(year) %>% filter(cob == k) %>% select(age, Birth)
      area_birth[[i]][j,,,k] = matrix(births_select$Birth, nrow = 7, ncol = 5, byrow = FALSE)
      rownames(area_birth[[i]][j,,,k]) = unique(births_select$age)
      colnames(area_birth[[i]][j,,,k]) = year_select
    }
  }
}

##############
# Extract ERP
##############

year_erp = unique(population_all$year)
age_all = unique(fertility_all$age)

ERP_area = list()
for(i in 1:47)
{
  ERP_area[[i]] = array(0, dim = c(8,7,19)) # dimension = year, age, cob
  for(j in 1:8)
  {
    year_select = year_erp[j]
    for(k in 1:19)
    {
      pop_select = filter(population_all, year == year_select) %>% filter(sex == "F") %>% filter(NSD == nsd_all[i]) %>% filter(cob == k) %>% filter(age %in% age_all) %>% select(age, ERP)
      ERP_area[[i]][j,,k] = pop_select$ERP
    }
  }
}

region_list = region_list_ind = list()
for(i in 1:length(unique(region_nsd$Region)))
{
  region_ind = unique(region_nsd$Region)[i]
  
  region_list[[i]] = as.numeric(unlist(filter(region_nsd, Region %in% region_ind) %>% select(NSD)))
  region_list_ind[[i]] = which(nsd_all %in% region_list[[i]])
}

ERP_region = list()
for(i in 1:11)
{
  nsd_id = match(region_list[[i]], nsd_all)
  ERP_region[[i]] = array(0, dim = c(8,7,19))
  for(j in 1:length(nsd_id))
  {
    ERP_region[[i]] = ERP_region[[i]] + ERP_area[[nsd_id[j]]]
  }
}

ERP_national = apply(Reduce("+", ERP_region), c(1,2,3), sum)


pop_region = list()
for(i in 1:11)
{
  pop_region[[i]] = array(0, dim = c(7,7,19)) # dimension = (year, age, cob)
  for(j in 1:7)
  {
    pop_region[[i]][j,,] = (ERP_region[[i]][j,,] + ERP_region[[i]][j+1,,])/2
  }
}


for(i in 1:11)
{
  print(sum(pop_region[[i]] == 0))
  
}

pop_area = list()
for(i in 1:47)
{
  pop_area[[i]] = array(0, dim = c(7,7,19)) # dimension = (year, age, cob)
  for(j in 1:7)
  {
    pop_area[[i]][j,,] = (ERP_area[[i]][j,,] + ERP_area[[i]][j+1,,])/2
  }
}

for(i in 1:47)
{
  print(sum(pop_area[[i]] == 0))
}


#####################################################################
# fill in missing birth observations at the most disaggregated level
####################################################################

birth_total_region = birth_total_area = list()
for(i in 1:11)
{
  birth_total_region[[i]] = apply(region_birth[[i]], c(1,2,4), sum)/5 # dimension = (year, age, cob)
}

for(i in 1:47)
{
  birth_total_area[[i]] = apply(area_birth[[i]], c(1,2,4), sum)/5 # dimension = (year, age, cob)
}

## For a particular age, fill in missing values via na.interp() if the number of missing values is smaller than 4.
## Otherwise, extrapolate in Age direction to fill in missings via approxExtrap.

fill_in_missing <- function(dat_mat)
{
  nT = nrow(dat_mat)
  nA = ncol(dat_mat)
  A_ind = 1:nA
  
  data_mat_na = dat_mat
  data_mat_na[dat_mat == 0] = NA # dimension = (year, age)
  
  update_ind = is.na(data_mat_na)
  
  if(all(is.na(data_mat_na)))
  {
    data_mat_update = dat_mat+1e-6
  } else {
    # Step 1: Intrapolate in Time direction
    data_mat_update = data_mat_na
    for(j in 1:nT)
    {
      if((sum(is.na(data_mat_na[,j])) > 0) && (sum(is.na(data_mat_na[,j])) <= 3))
      {
        fill_value = as.numeric(na.interp(data_mat_na[,j]))
        if(sum(fill_value<=0) > 0)
        {
          negative_ind = which(fill_value<=0)
          fill_value[negative_ind] = NA
        }
        
        data_mat_update[,j] = fill_value
      }
    }
    
    # Step 2: Extrapolate in Age direction
    if(sum(is.na(data_mat_update))>0)
    {
      for(i in 1:nA)
      {
        if(sum(is.na(data_mat_update[i,])) > 0)
        {
          na_ind = A_ind[is.na(data_mat_update[i,])]
          if(sum(is.na(data_mat_update[i,])) == 7)
          {
            # no observation for the year
            mean_all_year = apply(data_mat_update, 1, mean, na.rm = TRUE)
            year_with_observation = which(!is.nan(mean_all_year))
            data_mat_update[i,] = mean_all_year[year_with_observation[which.min(abs(year_with_observation - i))]]
          } else {
            if(sum(is.na(data_mat_update[i,])) >= 5)
            {
              # one or two observations for the year
              data_mat_update[i,na_ind] = mean(data_mat_update[i,], na.rm = TRUE)
            } else {
              fill_value_first_half = data_mat_update[i,1:4]
              first_half_na_ind = is.na(fill_value_first_half)
              
              #############################
              ## first half of age vectors
              #############################
              
              if(sum(is.na(data_mat_update[i,1:4])) == 4)
              {
                # Missing the first 4 observations
                fill_value_first_half[first_half_na_ind] = rep(data_mat_update[i,5], sum(first_half_na_ind))
              }
              
              if(sum(is.na(data_mat_update[i,1:4])) == 3)
              {
                # Missing the first 3 observations
                fill_value_first_half[first_half_na_ind] = rep(mean(data_mat_update[i,1:4], na.rm = TRUE), sum(first_half_na_ind))
              }
              
              if(sum(is.na(data_mat_update[i,1:4])) == 2)
              {
                fill_value_first_half[first_half_na_ind] = mean(data_mat_update[i,1:4], na.rm = TRUE)
              }
              
              if(sum(is.na(data_mat_update[i,1:4])) == 1)
              {
                if(is.na(data_mat_update[i,1]))
                {
                  if(data_mat_update[i,2] > data_mat_update[i,3])
                  {
                    fill_value_first_half[first_half_na_ind] = c(rep(data_mat_update[i,2], 2), data_mat_update[i,3:4])[first_half_na_ind]
                  } else {
                    fill_age_1 = 15
                    
                    fill_value_first_half[first_half_na_ind] = rev(Hmisc::approxExtrap(x = rev(age_all[2:4]), y = rev(data_mat_update[i,2:4]), xout = rev(c(fill_age_1, age_all[2:4])))$y)[first_half_na_ind]
                    while(fill_value_first_half[1] < 0.01)
                    {
                      fill_age_1 = fill_age_1 + 0.5
                      fill_value_first_half[first_half_na_ind] = rev(Hmisc::approxExtrap(x = rev(age_all[2:4]), y = rev(data_mat_update[i,2:4]), xout = rev(c(fill_age_1, age_all[2:4])))$y)[first_half_na_ind]
                      
                      if(fill_age_1 > 19.5)
                      {
                        break
                      }
                    }
                  }
                }
                
                if(!is.na(data_mat_update[i,1]))
                {
                  fill_value_first_half[first_half_na_ind] = rev(Hmisc::approxExtrap(x = rev(age_all[1:4]), y = rev(data_mat_update[i,1:4]), xout = rev(age_all[1:4]))$y)[first_half_na_ind]
                }
              }
              
              ##############################
              ## second half of age vectors
              ##############################
              
              fill_value_second_half = data_mat_update[i,5:7]
              second_half_na_ind = is.na(fill_value_second_half)
              
              if(sum(is.na(data_mat_update[i,5:7])) == 3)
              {
                # Missing the last 3 observations
                fill_value_second_half[second_half_na_ind] = rep(rev(fill_value_first_half)[4], 3)
              }
              
              if(sum(is.na(data_mat_update[i,5:7])) == 2)
              {
                if(is.na(data_mat_update[i,6]) & is.na(data_mat_update[i,7]))
                {
                  # Missing last 2 observations
                  if(is.na(data_mat_update[i,4]))
                  {
                    fill_value_second_half[second_half_na_ind] = rep(mean(data_mat_update[i,], na.rm = TRUE), 2)
                  } else {
                    if(data_mat_update[i,5] < 0.05)
                    {
                      fill_value_second_half[second_half_na_ind] = rep(data_mat_update[i,5], 2)
                    } else {
                      fill_age_2 = c(40,43)
                      fill_value_second_half[second_half_na_ind] = (Hmisc::approxExtrap(x = age_all[4:5], y = data_mat_update[i,4:5], xout = c(age_all[5], fill_age_2))$y)[second_half_na_ind]
                      while(fill_value_second_half[2] < 0.01 || fill_value_second_half[3] < 0.01)
                      {
                        fill_age_2[1] = fill_age_2[1] - 0.2
                        fill_age_2[2] = fill_age_2[2] - 0.2
                        fill_value_second_half[second_half_na_ind] = (Hmisc::approxExtrap(x = age_all[4:5], y = data_mat_update[i,4:5], xout = c(age_all[5], fill_age_2))$y)[second_half_na_ind]
                        if(fill_age_2[2] < 41)
                        {
                          if(fill_value_second_half[2] > 0)
                          {
                            fill_value_second_half[second_half_na_ind] = rep(fill_value_second_half[2], 2)
                          } else {
                            fill_value_second_half[second_half_na_ind] = rep(fill_value_second_half[1], 2)
                          }
                          break
                        }
                      }
                    }
                  }
                }
                
                if(is.na(data_mat_update[i,5]) & is.na(data_mat_update[i,7]))
                {
                  # Missing the 5th and 7th observations
                  fill_value_second_half_1 = Hmisc::approxExtrap(x = age_all[4:6], y = c(fill_value_first_half[4], data_mat_update[i,5:6]), xout = age_all[4:6])$y
                  if(fill_value_second_half_1[2] < fill_value_second_half_1[3])
                  {
                    fill_value_second_half[second_half_na_ind] = c(fill_value_second_half_1[-1], fill_value_second_half_1[2])[second_half_na_ind]
                    # prevent strictly increasing fertility at 35, 40, 45
                  } else {
                    fill_age_2 = 45
                    fill_value_second_half[second_half_na_ind] = (Hmisc::approxExtrap(x = age_all[4:6], y = fill_value_second_half_1, xout = c(age_all[4:6], fill_age_2))$y)[c(FALSE, second_half_na_ind)]
                    
                    while(fill_value_second_half[3] < 0.01)
                    {
                      fill_age_2 = fill_age_2 - 0.5
                      fill_value_second_half[second_half_na_ind] = (Hmisc::approxExtrap(x = age_all[4:6], y = fill_value_second_half_1, xout = c(age_all[4:6], fill_age_2))$y)[c(FALSE, second_half_na_ind)]
                    }
                  }
                }
                
                if(is.na(data_mat_update[i,5]) & is.na(data_mat_update[i,6]))
                {
                  # Missing the 5th and 6th observations
                  fill_value_second_half[second_half_na_ind] = (Hmisc::approxExtrap(x = age_all[4:7], y = c(fill_value_first_half[4], data_mat_update[i,5:7]), xout = age_all[4:7])$y)[c(FALSE, second_half_na_ind)]
                }
              }
              
              if(sum(is.na(data_mat_update[i,5:7])) == 1)
              {
                if(is.na(data_mat_update[i,7]))
                {
                  fill_age_2 = 45
                  fill_value_second_half[second_half_na_ind] = (Hmisc::approxExtrap(x = age_all[4:6], y = data_mat_update[i,4:6], xout = c(age_all[4:6], fill_age_2))$y)[c(FALSE, second_half_na_ind)]
                  
                  while(fill_value_second_half[3] < 0.01)
                  {
                    fill_age_2 = fill_age_2 - 0.5
                    fill_value_second_half[second_half_na_ind] = (Hmisc::approxExtrap(x = age_all[4:6], y = data_mat_update[i,4:6], xout = c(age_all[4:6], fill_age_2))$y)[c(FALSE, second_half_na_ind)]
                    if(fill_age_2 == 40.5)
                    {
                      fill_value_second_half[second_half_na_ind] = data_mat_update[i,6]
                      break
                    }
                  }
                } else {
                  fill_value_second_half[second_half_na_ind] = (Hmisc::approxExtrap(x = age_all[4:7], y = data_mat_update[i,4:7], xout = age_all[4:7])$y)[c(FALSE, second_half_na_ind)]
                }
              }
              
              
              # combine both halves
              fill_value_both = c(fill_value_first_half, fill_value_second_half)
              
              data_mat_update[i,na_ind] = fill_value_both[na_ind]
              
              rm(fill_value_first_half, fill_value_second_half, first_half_na_ind, second_half_na_ind)
            }
          }
        }
      }
    }
  }
  
  return(list(data_mat_update = data_mat_update, update_ind = update_ind))
}


birth_update_area = birth_update_ind_area = list()
for(i in 1:47)
{
  birth_update_area[[i]] = birth_update_ind_area[[i]] = array(0, dim = c(7,7,19))
  for(j in 1:19)
  {
    birth_update_temp = fill_in_missing(birth_total_area[[i]][,,j])
    
    birth_update_area[[i]][,,j] = birth_update_temp$data_mat_update
    birth_update_ind_area[[i]][,,j] = birth_update_temp$update_ind
    print(paste("Area ", i, " COB ", j, " complete", sep = "")) 
  }
}

sum(is.na(birth_update_area)) # return 0

#########################################
# Update region birth and national birth
#########################################

birth_update_region = birth_update_ind_region = list()
for(i in 1:11)
{
  # Birth
  birth_update_region[[i]] = array(0, dim = c(7,7,19))
  for(k in 1:length(region_list_ind[[i]]))
  {
    ind_select = region_list_ind[[i]][k]
    birth_update_select = birth_update_area[[ind_select]]
    birth_update_region[[i]] = birth_update_region[[i]] + birth_update_select
  }
  
  # Update index
  update_ind_region = array(0, dim = c(7,7,19))
  for(k in 1:length(region_list_ind[[i]]))
  {
    ind_select = region_list_ind[[i]][k]
    birth_update_ind_select = birth_update_ind_area[[ind_select]]
    update_ind_region = update_ind_region + birth_update_ind_select
  }
  birth_update_ind_region[[i]] = (update_ind_region != 0)
  rm(update_ind_region)
}

birth_update_national = Reduce('+', birth_update_area)

##################################################
# Update Population at the area and region levels
##################################################

pop_update_area = pop_update_region = pop_update_ind_area = pop_update_ind_region = list()

for(i in 1:47)
{
  pop_update_area[[i]] = pop_update_ind_area[[i]] = array(0, dim = c(7,7,19))
  for(j in 1:19)
  {
    pop_update_temp = fill_in_missing(pop_area[[i]][,,j])
  
    # check if the filled population value is less than the smallest observation
    for(k in 1:7)
    {
      if(sum(pop_update_temp$update_ind[k,]) > 0 & sum(pop_update_temp$update_ind[k,]) < 6)
      {
        small_fill_ind = pop_update_temp$data_mat_update[k,] < min(pop_area[[i]][k,,j][!pop_update_temp$update_ind[k,]])
        if(sum(small_fill_ind) > 0)
        {
          pop_update_temp$data_mat_update[k,small_fill_ind] = min(pop_area[[i]][k,,j][!pop_update_temp$update_ind[k,]])
        }
      }
    }

    # check if population is too small compared to birth number
    large_birth_index = (birth_update_area[[i]][,,j]/pop_update_temp$data_mat_update > 1)
      
    if(sum(large_birth_index) > 0)
    {
      if(sum(birth_total_area[[i]][,,j][large_birth_index]) == 0)
      {
        row_ind = ifelse(which(large_birth_index==TRUE)%%7 > 0, which(large_birth_index==TRUE)%%7, 1)
        col_ind = ifelse(which(large_birth_index==TRUE)%/%7 + 1 < 8, which(large_birth_index==TRUE)%/%7 + 1, 7)
        
        for(ik in 1:length(col_ind))
        {
          ri = row_ind[ik]
          ci = col_ind[ik]
          
          birth_update_area[[i]][ri,ci,j] = mean(birth_update_area[[i]][-ri,ci,j], na.rm = TRUE)
        }
      } else {
        if(sum(pop_update_temp$data_mat_update[large_birth_index]) == 0)
        {
          pop_update_temp$data_mat_update[large_birth_index] = birth_total_area[[i]][,,j][large_birth_index]
        }
      }

      warning(paste("Area ", i, " COB ", j, " Pop too small. Filled birth number changed.", sep = ""))
    }
    
    pop_update_area[[i]][,,j] = pop_update_temp$data_mat_update
    pop_update_ind_area[[i]][,,j] = pop_update_temp$update_ind
    print(paste("Area ", i, " COB ", j, " complete", sep = "")) 
  }
}


pop_update_region = pop_update_ind_region = list()
for(i in 1:11)
{
  # population
  pop_update_region[[i]] = array(0, dim = c(7,7,19))
  for(k in 1:length(region_list_ind[[i]]))
  {
    ind_select = region_list_ind[[i]][k]
    pop_update_select = pop_update_area[[ind_select]]
    pop_update_region[[i]] = pop_update_region[[i]] + pop_update_select
  }
  
  # Update index
  update_ind_region = array(0, dim = c(7,7,19))
  for(k in 1:length(region_list_ind[[i]]))
  {
    ind_select = region_list_ind[[i]][k]
    pop_update_ind_select = pop_update_ind_area[[ind_select]]
    update_ind_region = update_ind_region + pop_update_ind_select
  }
  pop_update_ind_region[[i]] = (update_ind_region != 0) # 1 means updated value
  rm(update_ind_region)
}

pop_update_national = Reduce('+', pop_update_area)

###################################################################################
# Interpolate in the time direction via na.spline (i.e., polynomial interpolation)
###################################################################################

time_interp <- function(dat_list)
{
  if(!is.list(dat_list))
  {
    n_age = dim(dat_list)[2]
    n_cob = dim(dat_list)[3]
    
    dat_interp = array(0, dim = c(31, n_age, n_cob))
    for(k in 1:n_cob)
    {
      for(j in 1:n_age)
      {
        y_missing = rep(NA, 31)
        y_missing[c(1,6,11,16,21,26,31)] = dat_list[,j,k]
        
        dat_interp[,j,k] = exp(na.spline(log(y_missing)))
      }
    }
  } else {
    L = length(dat_list)
    n_year = dim(dat_list[[1]])[1]
    if(n_year == 8)
    {
      n_age = dim(dat_list[[1]])[2]
      n_cob = dim(dat_list[[1]])[3]
      
      dat_interp = list()
      for(l in 1:L)
      {
        dat_filled = array(0, dim = c(31, n_age, n_cob))
        for(k in 1:n_cob)
        {
          for(j in 1:n_age)
          {
            y_missing = rep(NA, 31)
            y_missing[c(1,6,11,16,21,26,31)] = dat_list[[l]][1:7,j,k]
            
            dat_filled[,j,k] = exp(na.spline(log(y_missing)))
          }
        }
        dat_interp[[l]] = dat_filled
      }
    } else {
      n_age = dim(dat_list[[1]])[2]
      n_cob = dim(dat_list[[1]])[3]
      
      dat_interp = list()
      for(l in 1:L)
      {
        dat_filled = array(0, dim = c(31, n_age, n_cob))
        for(k in 1:n_cob)
        {
          for(j in 1:n_age)
          {
            y_missing = rep(NA, 31)
            y_missing[c(1,6,11,16,21,26,31)] = dat_list[[l]][,j,k]
            
            dat_filled[,j,k] = exp(na.spline(log(y_missing)))
          }
        }
        dat_interp[[l]] = dat_filled
      }
    }
  }
  
  return(dat_interp)
}

# Interpolate birth data
birth_interp_national = time_interp(birth_update_national)
birth_interp_region = time_interp(birth_update_region)
birth_interp_area = time_interp(birth_update_area)

## check
sum(birth_interp_national <= 0)

for(i in 1:11)
{
  print(sum(birth_interp_region[[i]] <= 0))
}

for(i in 1:47)
{
  print(sum(birth_interp_area[[i]] <= 0))
}

# Interpolate population data
pop_interp_national = time_interp(pop_update_national)
pop_interp_region = time_interp(pop_update_region)
pop_interp_area = time_interp(pop_update_area)

## check
sum(pop_interp_national <= 0)

for(i in 1:11)
{
  print(sum(pop_interp_region[[i]] <= 0))
}

for(i in 1:47)
{
  print(sum(pop_interp_area[[i]] <= 0))
}

#####################################################
# Compute age-specific fertility for region and area
#####################################################

# '%notin%' <- Negate('%in%')
# 
# # national level
# fertility_interp_national = birth_interp_national/pop_interp_national
# 
# # region level
# fertility_interp_region_fix_cob = list()
# for(i in 1:11)
# {
#   fertility_interp_region_fix_cob[[i]] = birth_interp_region[[i]]/pop_interp_region[[i]] # dimension = (year, age, cob)
# }
# 
# # area level
# fertility_interp_area_fix_cob = list()
# for(i in 1:47)
# {
#   fertility_interp_area_fix_cob[[i]] = birth_interp_area[[i]]/pop_interp_area[[i]] # dimension = (year, age, cob)
# }

#########################################
# Compute age-specific fertility for COB 
#########################################

COB_list = list()
COB_list[[1]] = 1:3
COB_list[[2]] = 4:5
COB_list[[3]] = 6
COB_list[[4]] = 7
COB_list[[5]] = 8:12
COB_list[[6]] = 13:14
COB_list[[7]] = 15:16
COB_list[[8]] = 17
COB_list[[9]] = 18
COB_list[[10]] = 19

# COB level (base level)
for(i in 1:19)
{
  # Population
  pop_array = array(0, dim = c(31, 7, 59))
  for(j in 1:47)
  {
    pop_array[,,j] = pop_interp_area[[j]][,,i]
  }
  
  for(k in 48:58)
  {
    pop_array[,,k] = pop_interp_region[[k-47]][,,i]
  }
  
  pop_array[,,59] = pop_interp_national[,,i]
  
  assign(paste("pop_interp_cob", i, sep=""),  pop_array)
  
  # Births
  birth_array = array(0, dim = c(31, 7, 59))
  for(j in 1:47)
  {
    birth_array[,,j] = birth_interp_area[[j]][,,i]
  }
  
  for(k in 48:58)
  {
    birth_array[,,k] = birth_interp_region[[k-47]][,,i]
  }
  
  birth_array[,,59] = birth_interp_national[,,i]
  
  assign(paste("birth_interp_cob", i, sep=""), birth_array)
  
  # Fertility
  fertility_array = array(0, dim = c(31, 7, 59))
  for(j in 1:47)
  {
    fertility_array[,,j] = fertility_interp_area[[j]][,,i]
  }
  
  for(k in 48:58)
  {
    fertility_array[,,k] = fertility_interp_region[[k-47]][,,i]
  }
  
  fertility_array[,,59] = fertility_interp_national[,,i]
  
  assign(paste("fertility_interp_cob", i, sep=""),  fertility_array)
  
  # clean memory
  rm(pop_array, birth_array, fertility_array)
}

for(i in 1:19)
{
  print(sum(get(paste("fertility_interp_cob",i,sep=""))==0))
}

# C_group level (middle level)

for(i in 1:10)
{
  cob_ind = COB_list[[i]]
  # Population
  pop_array = array(0, dim = c(31, 7, 59))
  for(j in 1:59)
  {
    for(iw in 1:length(cob_ind))
    {
      pop_array[,,j] = pop_array[,,j] + get(paste("pop_interp_cob", cob_ind[iw], sep = ""))[,,j]
    }
  }
  
  assign(paste("pop_interp_C", i, sep = ""),  pop_array)
  
  # Births
  birth_array = array(0, dim = c(31, 7, 59))
  for(j in 1:59)
  {
    for(iw in 1:length(cob_ind))
    {
      birth_array[,,j] = birth_array[,,j] + get(paste("birth_interp_cob", cob_ind[iw], sep = ""))[,,j]
    }
  }
  
  assign(paste("birth_interp_C", i, sep = ""),  birth_array)
  
  # Fertility
  fertility_array = array(0, dim = c(31, 7, 59))
  for(j in 1:59)
  {
    for(iw in 1:length(cob_ind))
    {
      fertility_array[,,j] = fertility_array[,,j] + get(paste("fertility_interp_cob", cob_ind[iw], sep = ""))[,,j]
    }
  }
  
  assign(paste("fertility_interp_C", i, sep = ""),  fertility_array)
  
  # clean memory
  rm(pop_array, birth_array, fertility_array)
}

for(i in 1:10)
{
  print(sum(get(paste("fertility_interp_C", i, sep = "")) <= 0))
}

# total level (top level)
pop_interp_all_cob = array(0, dim = c(31, 7, 59))
for(j in 1:59)
{
  for(iw in 1:19)
  {
    pop_interp_all_cob[,,j] = pop_interp_all_cob[,,j] + get(paste("pop_interp_cob", iw, sep = ""))[,,j]
  }
}

birth_interp_all_cob = array(0, dim = c(31, 7, 59))
for(j in 1:59)
{
  for(iw in 1:19)
  {
    birth_interp_all_cob[,,j] = birth_interp_all_cob[,,j] + get(paste("birth_interp_cob", iw, sep = ""))[,,j]
  }
}

fertility_interp_all_cob = array(0, dim = c(31, 7, 59))
for(j in 1:59)
{
  # for(iw in 1:19)
  # {
  #   fertility_interp_all_cob[,,j] = fertility_interp_all_cob[,,j] + get(paste("fertility_interp_cob", iw, sep = ""))[,,j]
  # }
  
  fertility_interp_all_cob = birth_interp_all_cob/pop_interp_all_cob
}

sum(fertility_interp_all_cob <= 0)

###################################################
# Create FTS objects after nonparametric smoothing
###################################################

library(demography)

##########
# Fix COB
##########

year_interp = 1981:2011

# national level

fts_national_fix_cob = fts_national_fix_cob_smooth = vector("character", length = 19)
for(i in 1:19)
{
  fts_national_fix_cob[i] = paste("fts_national_cob", i, sep = "")
  fts_national_fix_cob_smooth[i] = paste("fts_national_smooth_cob", i, sep = "")
}

for(i in 1:19)
{
  demo = demogdata(data = t(fertility_interp_national[,,i]), pop = t(pop_interp_national[,,i]), type = "fertility", ages = age_all, years = year_interp, label = paste("Australia COB", i, sep = " "), name = "female")
  demo_smooth = smooth.demogdata(demo)
  assign(fts_national_fix_cob[i], demo)
  assign(fts_national_fix_cob_smooth[i], demo_smooth)
  rm(demo,demo_smooth)
}

# region level

fts_region_fix_cob = fts_region_fix_cob_smooth = list()
for(j in 1:11)
{
  fts_region_fix_cob[[j]] = vector("character", length = 19)
  fts_region_fix_cob_smooth[[j]] = vector("character", length = 19)
  
  for(i in 1:19)
  {
    fts_region_fix_cob[[j]][i] = paste("fts_region_", j, "_cob", i, sep = "")
    fts_region_fix_cob_smooth[[j]][i] = paste("fts_region_smooth_", j, "_cob", i, sep = "")
    
    demo = demogdata(data = t(fertility_interp_region[[j]][,,i]), pop = t(pop_interp_region[[j]][,,i]), type = "fertility", ages = age_all, years = year_interp, label = paste("Region ", j, " COB ", i, sep = " "), name = "female")
    demo_smooth = smooth.demogdata(demo)
    
    assign(fts_region_fix_cob[[j]][i], demo)
    assign(fts_region_fix_cob_smooth[[j]][i], demo_smooth)
    
    rm(demo, demo_smooth)
  }
}

# area level

fts_area_fix_cob = fts_area_fix_cob_smooth = list()
for(j in 1:47)
{
  fts_area_fix_cob[[j]] = vector("character", length = 19)
  fts_area_fix_cob_smooth[[j]] = vector("character", length = 19)
  
  for(i in 1:19)
  {
    fts_area_fix_cob[[j]][i] = paste("fts_area_", j, "_cob", i, sep = "")
    fts_area_fix_cob_smooth[[j]][i] = paste("fts_area_smooth_", j, "_cob", i, sep = "")
    
    demo = demogdata(data = t(fertility_interp_area[[j]][,,i]), pop = t(pop_interp_area[[j]][,,i]), type = "fertility", ages = age_all, years = year_interp, label = paste("Area ", j, " COB ", i, sep = " "), name = "female")
    demo_smooth = smooth.demogdata(demo)
    
    assign(fts_area_fix_cob[[j]][i], demo)
    assign(fts_area_fix_cob_smooth[[j]][i], demo_smooth)
    
    rm(demo, demo_smooth)
  }
}


for(i in 1:47)
{
  for(j in 1:19)
  {
    print(paste("area",i,"cob",j))
    print(sum(get(fts_area_fix_cob[[i]][j])$rate$female>1))
  }
}

fts_area_23_cob14 = fts_area_23_cob14 = NULL


##########
# Fix geo
##########

# national level

fts_national_fix_geo = fts_national_fix_geo_smooth = vector("character", length = 59) # Area, Region, National
for(i in 1:59)
{
  fts_national_fix_geo[i] = paste("fts_cob_all_geo", i, sep = "")
  fts_national_fix_geo_smooth[i] = paste("fts_cob_all_smooth_geo", i, sep = "")
}

for(i in 1:59)
{
  demo = demogdata(data = t(fertility_interp_all_cob[,,i]), pop = t(pop_interp_all_cob[,,i]), type = "fertility", ages = age_all, years = year_interp, label = paste("Australia GEO", i, sep = " "), name = "female")
  demo_smooth = smooth.demogdata(demo)
  
  assign(fts_national_fix_geo[i], demo)
  assign(fts_national_fix_geo_smooth[i], demo_smooth)
  
  rm(demo,demo_smooth)
}


# region level

fts_region_fix_geo = fts_region_fix_geo_smooth = list()
for(j in 6:10)
{
  fts_region_fix_geo[[j]] = vector("character", length = 59)
  fts_region_fix_geo_smooth[[j]] = vector("character", length = 59)
  
  for(i in 1:59)
  {
    if (j == 6 & i == 23)
    {
      fts_region_fix_geo[[j]][i] = paste("fts_C_", j, "_geo", i, sep = "")
      fts_region_fix_geo_smooth[[j]][i] = paste("fts_smooth_C", j, "_geo", i, sep = "")
      
      assign(fts_region_fix_geo[[j]][i], NA)
      assign(fts_region_fix_geo_smooth[[j]][i], NA)
    } else {
      fts_region_fix_geo[[j]][i] = paste("fts_C_", j, "_geo", i, sep = "")
      fts_region_fix_geo_smooth[[j]][i] = paste("fts_smooth_C", j, "_geo", i, sep = "")
      
      demo = demogdata(data = t(get(paste("fertility_interp_C", j, sep = ""))[,,i]), pop = t(get(paste("pop_interp_C", j, sep = ""))[,,i]), type = "fertility", ages = age_all, years = year_interp, label = paste("C ", j, " GEO ", i, sep = " "), name = "female")
      demo_smooth = smooth.demogdata(demo)
      assign(fts_region_fix_geo[[j]][i], demo)
      assign(fts_region_fix_geo_smooth[[j]][i], demo_smooth)
      
      rm(demo, demo_smooth)
    }
  }
}


# area level

fts_area_fix_geo = fts_area_fix_geo_smooth = list()
for(j in 1:19)
{
  fts_area_fix_geo[[j]] = vector("character", length = 59)
  fts_area_fix_geo_smooth[[j]] = vector("character", length = 59)
  
  for(i in 1:59)
  {
    fts_area_fix_geo[[j]][i] = paste("fts_cob_", j, "_geo", i, sep = "")
    fts_area_fix_geo_smooth[[j]][i] = paste("fts_smooth_cob", j, "_geo", i, sep = "")
    
    demo = demogdata(data = t(get(paste("fertility_interp_cob", j, sep = ""))[,,i]), pop = t(get(paste("pop_interp_cob", j, sep = ""))[,,i]), type = "fertility", ages = age_all, years = year_interp, label = paste("COB ", j, " GEO ", i, sep = " "), name = "female")
    demo_smooth = smooth.demogdata(demo)
    
    assign(fts_area_fix_geo[[j]][i], demo)
    assign(fts_area_fix_geo_smooth[[j]][i], demo_smooth)
    
    rm(demo, demo_smooth)
  }
}


## Make rainbow plots

# National raw & smoothed

national_raw_demog = demogdata(data = t(apply(Reduce("+", birth_total_area), c(1,2), sum)/apply(Reduce("+", pop_area), c(1,2), sum)), pop = t(apply(Reduce("+", pop_area), c(1,2), sum)), type = "fertility", ages = age_all, years = year_all, label = "Australia", name = "raw")

savepdf("Australia_raw", width = 24, height = 24, toplines = 0.9)
par(mar = c(5,7,4,2) + 0.1)
plot(national_raw_demog, ylab = "", xlab = "", cex.main = 2,  xaxt = "n", yaxt = "n", ylim = c(0, 0.16))
axis(side = 1, at = age_all, cex.axis = 2, mgp = c(3, 1, 0))
axis(side = 2, at = seq(0, 0.16, 0.04), cex.axis = 2, las = 1)
title(ylab = "Age-specific fertility rate", line = 5, cex.lab = 2)
title(xlab = "Age", line = 3, cex.lab = 2)
legend("topright", lty = c(1,1), col = c(2,"purple"), c("1981", "2011"), cex = 2)
dev.off()


national_interp_demog = demogdata(data = t(apply(fertility_interp_national*pop_interp_national, c(1,2), sum)/apply(pop_interp_national, c(1,2), sum)), pop = t(apply(pop_interp_national, c(1,2), sum)), type = "fertility", ages = age_all, years = year_interp, label = "Australia", name = "smoothed")
national_interp_smooth_demog = smooth.demogdata(national_interp_demog)

savepdf("Australia_interp", width = 24, height = 24, toplines = 0.9)
par(mar = c(5,7,4,2) + 0.1)
plot(national_interp_smooth_demog, ylab = "", xlab = "", yaxt = "n", cex.main = 2,  xaxt = "n", ylim = c(0, 0.16))
axis(side = 1, at = age_all, cex.axis = 2, mgp = c(3, 1, 0))
axis(side = 2, at = seq(0, 0.16, 0.04), cex.axis = 2, las = 1)
title(xlab = "Age", line = 3, cex.lab = 2)
dev.off()

# Area 1 total

area_1_interp_demog = demogdata(data = t(apply(fertility_interp_area[[1]]*pop_interp_area[[1]], c(1,2), sum)/apply(pop_interp_area[[1]], c(1,2), sum)), pop = t(apply(pop_interp_area[[1]], c(1,2), sum)), type = "fertility", ages = age_all, years = year_interp, label = "Sydney", name = "smoothed")
area_1_interp_smooth_demog = smooth.demogdata(area_1_interp_demog)

savepdf("A1_interp", width = 24, height = 24, toplines = 0.9)
par(mar = c(5,7,4,2) + 0.1)
plot(area_1_interp_smooth_demog, cex.main = 2, ylab = "", xlab = "", xaxt = "n", yaxt = "n", ylim = c(0, 0.16))
axis(side = 1, at = age_all, cex.axis = 2, mgp = c(3, 1, 0))
axis(side = 2, at = seq(0, 0.16, 0.04), cex.axis = 2, las = 1)
title(ylab = "Age-specific fertility rate", line = 5, cex.lab = 2)
title(xlab = "Age", line = 3, cex.lab = 2)
dev.off()

# Area 11 total

area_11_interp_demog = demogdata(data = t(apply(fertility_interp_area[[11]]*pop_interp_area[[11]], c(1,2), sum)/apply(pop_interp_area[[11]], c(1,2), sum)), pop = t(apply(pop_interp_area[[11]], c(1,2), sum)), type = "fertility", ages = age_all, years = year_interp, label = "Remote Australia", name = "smoothed")
area_11_interp_smooth_demog = smooth.demogdata(area_11_interp_demog)

savepdf("A11_interp", width = 24, height = 24, toplines = 0.9)
par(mar = c(5,7,4,2) + 0.1)
plot(area_11_interp_smooth_demog, cex.main = 2, ylab = "", xlab = "", xaxt = "n", yaxt = "n", ylim = c(0, 0.165))
axis(side = 1, at = age_all, cex.axis = 2, mgp = c(3, 1, 0))
axis(side = 2, at = seq(0, 0.165, 0.04), cex.axis = 2, las = 1)
title(ylab = "Age-specific fertility rate", line = 5, cex.lab = 2)
title(xlab = "Age", line = 3, cex.lab = 2)
dev.off()


# Zone 1 total
cob_1_interp_demog = demogdata(data = t(apply(apply(fertility_interp_cob1*pop_interp_cob1, c(1,2), sum)/apply(pop_interp_cob1, c(1,2), sum), c(1,2), sum)), pop = t(apply(pop_interp_cob1, c(1,2), sum)), type = "fertility", ages = age_all, years = year_interp, label = "Birthplace Australia", name = "smoothed")
cob_1_interp_smooth_demog = smooth.demogdata(cob_1_interp_demog)

savepdf("COB_1_interp", width = 24, height = 24, toplines = 0.9)
par(mar = c(5,7,4,2) + 0.1)
plot(cob_1_interp_smooth_demog, cex.main = 2,  ylab = "", xlab = "", xaxt = "n", yaxt = "n", ylim = c(0, 0.16))
axis(side = 1, at = age_all, cex.axis = 2, mgp = c(3, 1, 0))
axis(side = 2, at = seq(0, 0.16, 0.04), cex.axis = 2, las = 1)
title(ylab = "Age-specific fertility rate", line = 5, cex.lab = 2)
title(xlab = "Age", line = 3, cex.lab = 2)
dev.off()


# Zone 4 total

cob_4_interp_demog = demogdata(data = t(apply(apply(fertility_interp_cob4*pop_interp_cob4, c(1,2), sum)/apply(pop_interp_cob4, c(1,2), sum), c(1,2), sum)), pop = t(apply(pop_interp_cob4, c(1,2), sum)), type = "fertility", ages = age_all, years = year_interp, label = "Birthplace UK", name = "smoothed")
cob_4_interp_smooth_demog = smooth.demogdata(cob_4_interp_demog)

savepdf("COB_4_interp", width = 24, height = 24, toplines = 0.9)
par(mar = c(5,7,4,2) + 0.1)
plot(cob_4_interp_smooth_demog, cex.main = 2,  ylab = "", xlab = "", xaxt = "n", yaxt = "n", ylim = c(0, 0.16))
axis(side = 1, at = age_all, cex.axis = 2, mgp = c(3, 1, 0))
axis(side = 2, at = seq(0, 0.16, 0.04), cex.axis = 2, las = 1)
title(ylab = "Age-specific fertility rate", line = 5, cex.lab = 2)
title(xlab = "Age", line = 3, cex.lab = 2)
dev.off()


##########
# Summary
##########


total_birth_cob = rep(0, 19)

for(i in 1:19)
{
  temp = 0
  for(j in 1:11)
  {
    temp = temp + sum(birth_total_region[[j]][,,i])
  }
  
  total_birth_cob[i] = temp
  
  rm(temp)
}

round(total_birth_cob)


total_birth_region = rep(0, 11)
for(j in 1:11)
{
  total_birth_region[j] = sum(birth_total_region[[j]])
}

round(total_birth_region)

pop_2011 = 0
for(i in 1:11)
{
  pop_2011 = pop_2011 + pop_interp_region[[i]][,7,]
}

(sum(pop_interp_region[[1]][,7,]) + sum(pop_interp_region[[2]][,7,]) + sum(pop_interp_region[[3]][,7,]) + sum(pop_interp_region[[5]][,7,]) + sum(pop_interp_region[[6]][,7,]) + sum(pop_interp_region[[7]][,7,]) + sum(pop_interp_region[[9]][,7,]))/sum(pop_2011)


sum(pop_interp_region[[11]][,7,])/sum(pop_2011)


###################
# compute variance
###################

age_ind = match(age_all, age_interp)

fertility_all_cob = list()
for(i in 1:47)
{
  fertility_all_cob[[i]] = apply(birth_interp_area[[i]][age_ind,,], c(1,2), sum)/apply(pop_interp_area[[i]][age_ind,,], c(1,2), sum)
}

all_cob_mad = rep(0, 47)
for(i in 1:47)
{
  all_cob_mad[i] = mean(apply(fertility_all_cob[[i]], 2, mad))
}

round(all_cob_mad, 4)

all_cob_mad[which.max(all_cob_mad)]/all_cob_mad[which.min(all_cob_mad)]






