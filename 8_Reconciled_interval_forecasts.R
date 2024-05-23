#########################################
# Base and reconciled interval forecasts
#########################################

interval_score = function(PI_val, data_series, fh, alpha = 0.8)
{
  require(demography)
  
  test_val = extract.years(data_series, ((2001+fh):2011))$rate$female
  boot_sample = PI_val[[fh]]
  
  # replace negative values
  boot_index_negative = which(boot_sample < 0)
  if(length(boot_index_negative) > 0)
  {
    boot_sample_v1 = replace(boot_sample, boot_index_negative, 0)
  } else {
    boot_sample_v1 = boot_sample
  }
  
  boot_index_positive = which(boot_sample_v1 > 1)
  if(length(boot_index_positive) > 0)
  {
    boot_sample_v2 = replace(boot_sample_v1, boot_index_positive, 1)
  } else {
    boot_sample_v2 = boot_sample_v1
  }
  
  # compute lower and upper bounds
  if(fh == 10)
  {
    dummy_boot = apply(boot_sample_v2, 1, quantile, c((1 - alpha)/2, (1 + alpha)/2), na.rm = TRUE)
    PI_lb_val = dummy_boot[1,]
    PI_ub_val  = dummy_boot[2,]
  } else {
    dummy_boot = apply(boot_sample_v2, c(1,3), quantile, c((1 - alpha)/2, (1 + alpha)/2), na.rm = TRUE)
    PI_lb_val = dummy_boot[1,,]
    PI_ub_val  = dummy_boot[2,,]
  }
  
  lb_ind = ifelse(test_val < PI_lb_val, 1, 0)
  ub_ind = ifelse(test_val > PI_ub_val, 1, 0)
  
  score = mean((PI_ub_val - PI_lb_val) + 2/(1 - alpha) * (PI_lb_val - test_val) * lb_ind + 2/(1 - alpha) * (test_val - PI_ub_val) * ub_ind)
  
  return(list(score = score))
}

interval_score_recon = function(PI_val, data_series, fh, alpha = 0.8)
{
  require(demography)
  
  test_val = extract.years(data_series, ((2001+fh):2011))$rate$female
  boot_sample = PI_val
  
  # replace negative values
  boot_index_negative = which(boot_sample < 0)
  if(length(boot_index_negative) > 0)
  {
    boot_sample_v1 = replace(boot_sample, boot_index_negative, 0)
  } else {
    boot_sample_v1 = boot_sample
  }
  
  boot_index_positive = which(boot_sample_v1 > 1)
  if(length(boot_index_positive) > 0)
  {
    boot_sample_v2 = replace(boot_sample_v1, boot_index_positive, 1)
  } else {
    boot_sample_v2 = boot_sample_v1
  }
  
  
  # compute lower and upper bounds
  dummy_boot = apply(boot_sample_v2, 1, quantile, c((1 - alpha)/2, (1 + alpha)/2), na.rm = TRUE)
  PI_lb_val = dummy_boot[1,]
  PI_ub_val  = dummy_boot[2,]
  
  lb_ind = ifelse(test_val < PI_lb_val, 1, 0)
  ub_ind = ifelse(test_val > PI_ub_val, 1, 0)
  
  score = mean((PI_ub_val - PI_lb_val) + 2/(1 - alpha) * (PI_lb_val - test_val) * lb_ind + 2/(1 - alpha) * (test_val - PI_ub_val) * ub_ind)
  
  return(list(score = score))
}


##########
# Fix COB
##########

# Area level
area_interval_score_all = vector("character", length = 47)
for(i in 1:47)
{
  area_interval_score_all[i] = paste("A", i, "_interval_score_fix_COB", sep = "")
}

for(ij in 1:47)
{
  interval_score_temp = matrix(0, nrow = 10, ncol = 19)
  for(k in 1:19)
  {
    if(ij  == 23 && k == 14)
    {
      interval_score_temp[,k] = rep(0, 10)
    } else {
      for(h in 1:10)
      {
        interval_score_temp[h,k] = interval_score(PI_val = get(area_PI_fh_all[ij])[[k]], data_series = get(fts_area_fix_cob[[ij]][k]), fh = h, alpha = 0.8)$score
      }
    }
  }
  assign(area_interval_score_all[ij], interval_score_temp)
  
  rm(interval_score_temp)
}

# Region level
region_interval_score_all = vector("character", length = 11)
for(i in 1:11)
{
  region_interval_score_all[i] = paste("R", i, "_interval_score_fix_COB", sep = "")
}

for(ij in 1:11)
{
  interval_score_temp = matrix(0, nrow = 10, ncol = 19)
  for(k in 1:19)
  {
    for(h in 1:10)
    {
      interval_score_temp[h,k] = interval_score(PI_val = get(region_PI_fh_all[ij])[[k]], data_series = get(fts_region_fix_cob[[ij]][k]), fh = h, alpha = 0.8)$score
    }
  }
  assign(region_interval_score_all[ij], interval_score_temp)
  
  rm(interval_score_temp)
}

# summary of interval scores

IS_area_fix_cob = array(0, dim = c(10,47,19))
for(ij in 1:47)
{
  for(cob in 1:19)
  {
    IS_area_fix_cob[,ij,cob] = get(area_interval_score_all[ij])[,cob]
  }
}

IS_region_fix_cob = array(0, dim = c(10, 11, 19))
for(ij in 1:11)
{
  for(cob in 1:19)
  {
    IS_region_fix_cob[,ij,cob] = get(region_interval_score_all[ij])[,cob]
  }
}

IS_fix_COB_all = array(0, dim = c(10,3,19))
for(cob in 1:19)
{
  IS_fix_COB_all[,1,cob] = ind_national_interval_score_fix_COB[,cob]
  IS_fix_COB_all[,2,cob] = rowMeans(IS_region_fix_cob[,,cob])
  IS_fix_COB_all[,3,cob] = rowMeans(IS_area_fix_cob[,,cob])
}

base_IS_fix_COB_averaged = apply(IS_fix_COB_all, c(1,2), mean)

##########
# Fix GEO
##########

# COB level

cob_interval_score_all = vector("character", length = 19)
for(i in 1:19)
{
  cob_interval_score_all[i] = paste("COB_", i, "_interval_score_fix_GEO", sep = "")
}

for(ij in 1:19)
{
  interval_score_temp = matrix(0, nrow = 10, ncol = 59)
  for(k in 1:59)
  {
    if(ij  == 14 && k == 23)
    {
      interval_score_temp[,k] = rep(0, 10)
    } else {
      for(h in 1:10)
      {
        interval_score_temp[h,k] = interval_score(PI_val = get(cob_PI_fh_all[ij])[[k]], data_series = get(fts_area_fix_geo[[ij]][k]), fh = h, alpha = 0.8)$score
      }
    }
  }
  assign(cob_interval_score_all[ij], interval_score_temp)
  
  rm(interval_score_temp)
}

# C group level

c_group_interval_score_all = vector("character", length = 10)
for(i in 1:10)
{
  c_group_interval_score_all[i] = paste("C_group_", i, "_interval_score_fix_GEO", sep = "")
}

for(ij in 1:10)
{
  interval_score_temp = matrix(0, nrow = 10, ncol = 59)
  for(k in 1:59)
  {
    if(ij  == 6 && k == 23)
    {
      interval_score_temp[,k] = rep(0, 10)
    } else {
      for(h in 1:10)
      {
        interval_score_temp[h,k] = interval_score(PI_val = get(c_group_PI_fh_all[ij])[[k]], data_series = get(fts_region_fix_geo[[ij]][k]), fh = h, alpha = 0.8)$score
      }
    }
  }
  assign(c_group_interval_score_all[ij], interval_score_temp)
  
  rm(interval_score_temp)
}

# summary of interval scores

IS_cob_fix_geo = array(0, dim = c(10,19,59))
for(ij in 1:19)
{
  for(geo in 1:59)
  {
    IS_cob_fix_geo[,ij,geo] = get(cob_interval_score_all[ij])[,geo]
  }
}

IS_c_group_fix_geo = array(0, dim = c(10, 10, 59))
for(ij in 1:10)
{
  for(geo in 1:59)
  {
    IS_c_group_fix_geo[,ij,geo] = get(c_group_interval_score_all[ij])[,geo]
  }
}

IS_fix_GEO_all = array(0, dim = c(10,3,59))
for(geo in 1:59)
{
  IS_fix_GEO_all[,1,geo] = ind_global_interval_score_fix_GEO[,geo]
  IS_fix_GEO_all[,2,geo] = rowMeans(IS_cob_fix_geo[,,geo])
  IS_fix_GEO_all[,3,geo] = rowMeans(IS_c_group_fix_geo[,,geo])
}

base_IS_fix_GEO_averaged = apply(IS_fix_GEO_all, c(1,2), mean)

#########################################################################################
# All-level bootstrapped base forecasts (B = 1000) using multivariate forecasting method
#########################################################################################

##########
# Fix COB
##########

A23_PI_fh_fix_COB[[14]] = list()
for(h in 1:9)
{
  A23_PI_fh_fix_COB[[14]][[h]] = array(rep(1e-6, 7000*(11-h)), dim = c(7, 1000, (11-h)))
}
A23_PI_fh_fix_COB[[14]][[10]] = matrix(rep(1e-6, 7000), nrow = 7, ncol = 1000)


PI_all_level_fix_cob <- function(fh, cob)
{
  all_level_fh_PI_fix_cob = array(0, dim = c(7, (11-fh), 59, 1000))
  
  if(fh == 10)
  {
    for(ij in 1:1000)
    {
      all_level_fh_PI_fix_cob[,,,ij] = cbind(ind_national_PI_all_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                             
                                             R1_PI_fh_fix_COB[[cob]][[fh]][,ij], R2_PI_fh_fix_COB[[cob]][[fh]][,ij], R3_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             R4_PI_fh_fix_COB[[cob]][[fh]][,ij], R5_PI_fh_fix_COB[[cob]][[fh]][,ij], R6_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             R7_PI_fh_fix_COB[[cob]][[fh]][,ij], R8_PI_fh_fix_COB[[cob]][[fh]][,ij], R9_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             R10_PI_fh_fix_COB[[cob]][[fh]][,ij], R11_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             
                                             A1_PI_fh_fix_COB[[cob]][[fh]][,ij], A2_PI_fh_fix_COB[[cob]][[fh]][,ij], A3_PI_fh_fix_COB[[cob]][[fh]][,ij], A4_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             A5_PI_fh_fix_COB[[cob]][[fh]][,ij], A6_PI_fh_fix_COB[[cob]][[fh]][,ij], A7_PI_fh_fix_COB[[cob]][[fh]][,ij], A8_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             A9_PI_fh_fix_COB[[cob]][[fh]][,ij], A10_PI_fh_fix_COB[[cob]][[fh]][,ij], A11_PI_fh_fix_COB[[cob]][[fh]][,ij], A12_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             A13_PI_fh_fix_COB[[cob]][[fh]][,ij], A14_PI_fh_fix_COB[[cob]][[fh]][,ij], A15_PI_fh_fix_COB[[cob]][[fh]][,ij], A16_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             A17_PI_fh_fix_COB[[cob]][[fh]][,ij], A18_PI_fh_fix_COB[[cob]][[fh]][,ij], A19_PI_fh_fix_COB[[cob]][[fh]][,ij], A20_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             A21_PI_fh_fix_COB[[cob]][[fh]][,ij], A22_PI_fh_fix_COB[[cob]][[fh]][,ij], A23_PI_fh_fix_COB[[cob]][[fh]][,ij], A24_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             A25_PI_fh_fix_COB[[cob]][[fh]][,ij], A26_PI_fh_fix_COB[[cob]][[fh]][,ij], A27_PI_fh_fix_COB[[cob]][[fh]][,ij], A28_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             A29_PI_fh_fix_COB[[cob]][[fh]][,ij], A30_PI_fh_fix_COB[[cob]][[fh]][,ij], A31_PI_fh_fix_COB[[cob]][[fh]][,ij], A32_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             A33_PI_fh_fix_COB[[cob]][[fh]][,ij], A34_PI_fh_fix_COB[[cob]][[fh]][,ij], A35_PI_fh_fix_COB[[cob]][[fh]][,ij], A36_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             A37_PI_fh_fix_COB[[cob]][[fh]][,ij], A38_PI_fh_fix_COB[[cob]][[fh]][,ij], A39_PI_fh_fix_COB[[cob]][[fh]][,ij], A40_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             A41_PI_fh_fix_COB[[cob]][[fh]][,ij], A42_PI_fh_fix_COB[[cob]][[fh]][,ij], A43_PI_fh_fix_COB[[cob]][[fh]][,ij], A44_PI_fh_fix_COB[[cob]][[fh]][,ij],
                                             A45_PI_fh_fix_COB[[cob]][[fh]][,ij], A46_PI_fh_fix_COB[[cob]][[fh]][,ij], A47_PI_fh_fix_COB[[cob]][[fh]][,ij])
    }
  } else {
    for(ij in 1:1000)
    {
      all_level_fh_PI_fix_cob[,,,ij] = cbind(ind_national_PI_all_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                             
                                             R1_PI_fh_fix_COB[[cob]][[fh]][,ij,], R2_PI_fh_fix_COB[[cob]][[fh]][,ij,], R3_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             R4_PI_fh_fix_COB[[cob]][[fh]][,ij,], R5_PI_fh_fix_COB[[cob]][[fh]][,ij,], R6_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             R7_PI_fh_fix_COB[[cob]][[fh]][,ij,], R8_PI_fh_fix_COB[[cob]][[fh]][,ij,], R9_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             R10_PI_fh_fix_COB[[cob]][[fh]][,ij,], R11_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             
                                             A1_PI_fh_fix_COB[[cob]][[fh]][,ij,], A2_PI_fh_fix_COB[[cob]][[fh]][,ij,], A3_PI_fh_fix_COB[[cob]][[fh]][,ij,], A4_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             A5_PI_fh_fix_COB[[cob]][[fh]][,ij,], A6_PI_fh_fix_COB[[cob]][[fh]][,ij,], A7_PI_fh_fix_COB[[cob]][[fh]][,ij,], A8_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             A9_PI_fh_fix_COB[[cob]][[fh]][,ij,], A10_PI_fh_fix_COB[[cob]][[fh]][,ij,], A11_PI_fh_fix_COB[[cob]][[fh]][,ij,], A12_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             A13_PI_fh_fix_COB[[cob]][[fh]][,ij,], A14_PI_fh_fix_COB[[cob]][[fh]][,ij,], A15_PI_fh_fix_COB[[cob]][[fh]][,ij,], A16_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             A17_PI_fh_fix_COB[[cob]][[fh]][,ij,], A18_PI_fh_fix_COB[[cob]][[fh]][,ij,], A19_PI_fh_fix_COB[[cob]][[fh]][,ij,], A20_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             A21_PI_fh_fix_COB[[cob]][[fh]][,ij,], A22_PI_fh_fix_COB[[cob]][[fh]][,ij,], A23_PI_fh_fix_COB[[cob]][[fh]][,ij,], A24_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             A25_PI_fh_fix_COB[[cob]][[fh]][,ij,], A26_PI_fh_fix_COB[[cob]][[fh]][,ij,], A27_PI_fh_fix_COB[[cob]][[fh]][,ij,], A28_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             A29_PI_fh_fix_COB[[cob]][[fh]][,ij,], A30_PI_fh_fix_COB[[cob]][[fh]][,ij,], A31_PI_fh_fix_COB[[cob]][[fh]][,ij,], A32_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             A33_PI_fh_fix_COB[[cob]][[fh]][,ij,], A34_PI_fh_fix_COB[[cob]][[fh]][,ij,], A35_PI_fh_fix_COB[[cob]][[fh]][,ij,], A36_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             A37_PI_fh_fix_COB[[cob]][[fh]][,ij,], A38_PI_fh_fix_COB[[cob]][[fh]][,ij,], A39_PI_fh_fix_COB[[cob]][[fh]][,ij,], A40_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             A41_PI_fh_fix_COB[[cob]][[fh]][,ij,], A42_PI_fh_fix_COB[[cob]][[fh]][,ij,], A43_PI_fh_fix_COB[[cob]][[fh]][,ij,], A44_PI_fh_fix_COB[[cob]][[fh]][,ij,],
                                             A45_PI_fh_fix_COB[[cob]][[fh]][,ij,], A46_PI_fh_fix_COB[[cob]][[fh]][,ij,], A47_PI_fh_fix_COB[[cob]][[fh]][,ij,])
    }
  }
  
  return(all_level_fh_PI_fix_cob)
}

test = PI_all_level_fix_cob(fh = 1, cob = 14)

# one- to 10-step-ahead bootstrapped base forecasts (B = 1000)

PI_base_fix_cob = vector("character", length = 19)
for(i in 1:19)
{
  PI_base_fix_cob[i] = paste("PI_base_", i, "_fix_COB", sep = "")
}

for(cob in 1:19)
{
  temp_list = list()
  for(h in 1:10)
  {
    temp_list[[h]] = PI_all_level_fix_cob(fh = h, cob = cob)
  }
  
  assign(PI_base_fix_cob[cob], temp_list)
  rm(temp_list)
}


# All-level bootstrapped grouped forecasts

recon_PI_fix_cob <- function(kj, age, hier_method = c("BU", "comb_OLS", "mint"), cob_ind)
{
  hier_method = match.arg(hier_method)
  
  age_ind = match(age, age_all)
  
  hier_fore = array(NA, dim = c(59, (11-kj), 1000)) # in total 59 series
  
  # Summing matrix for kj horizon at a given age
  summ_mat = Smat_fun_fix_cob(kj = kj, age = age, no_area = 47, cob = cob_ind)
  
  # ik: number of years in the forecasting period, it depends on kj
  for(ik in 1:(11-kj))
  {
    hier = summ_mat[,,ik]
    
    if(hier_method == "BU")
    {
      for(ij in 1:1000)
      {
        hier_fore[,ik,ij] = hier %*% (get(paste("PI_base_", cob_ind, "_fix_COB", sep = "")))[[kj]][age_ind,ik,13:59,ij]
      }
    }
    
    if(hier_method == "comb_OLS")
    {
      for(ij in 1:1000)
      {
        hier_fore[,ik,ij] = hier %*% ginv(t(hier) %*% hier) %*% t(hier) %*% (get(paste("PI_base_", cob_ind, "_fix_COB", sep = "")))[[kj]][age_ind,ik,,ij]
      }
    }
    
    if(hier_method == "mint")
    {
      ginv = MASS:::ginv
      wh = wh_fun_fix_cob(kj = kj, age = age, cob_ind = cob_ind)
      for(ij in 1:1000)
      {
        hier_fore[,ik,ij] = hier %*% ginv(t(hier) %*% ginv(wh) %*% hier) %*% t(hier) %*% ginv(wh) %*% (get(paste("PI_base_", cob_ind, "_fix_COB", sep = "")))[[kj]][age_ind,ik,,ij]
      }
    }
  }
  return(hier_fore)
}

test = recon_PI_fix_cob(kj = 10, age = age_all[4], hier_method = "mint", cob_ind = 14)

# Calculation of reconciled interval predictions (hier_method = "BU")
BU_recon_PI_fix_cob = list()

for(cob in 1:19)
{
  BU_recon_PI_fix_cob[[cob]] = list()
  for(h in 1:10)
  {
    BU_recon_PI_fix_cob[[cob]][[h]] = array(NA, dim = c(7, 59, (11-h), 1000))
  
    for(ij in 1:7)
    {
      BU_recon_PI_fix_cob[[cob]][[h]][ij,,,] = recon_PI_fix_cob(kj = h, age = age_all[ij], hier_method = "BU", cob_ind = cob)
    }
  }
}

# Calculation of reconciled interval predictions (hier_method = "comb_OLS")
OLS_recon_PI_fix_cob = list()

for(cob in 1:19)
{
  OLS_recon_PI_fix_cob[[cob]] = list()
  for(h in 1:10)
  {
    OLS_recon_PI_fix_cob[[cob]][[h]] = array(NA, dim = c(7, 59, (11-h), 1000))
    
    for(ij in 1:7)
    {
      OLS_recon_PI_fix_cob[[cob]][[h]][ij,,,] = recon_PI_fix_cob(kj = h, age = age_all[ij], hier_method = "comb_OLS", cob_ind = cob)
    }
  }
}

# Calculation of reconciled interval predictions (hier_method = "mint")
MINT_recon_PI_fix_cob = list()

for(cob in 1:19)
{
  MINT_recon_PI_fix_cob[[cob]] = list()
  for(h in 1:10)
  {
    MINT_recon_PI_fix_cob[[cob]][[h]] = array(NA, dim = c(7, 59, (11-h), 1000))
    
    for(ij in 1:7)
    {
      MINT_recon_PI_fix_cob[[cob]][[h]][ij,,,] = recon_PI_fix_cob(kj = h, age = age_all[ij], hier_method = "mint", cob_ind = cob)
    }
  }
}

# function for computing reconciled interval scores
interval_score_recon_fix_cob <- function(PI_val_list, cob_ind)
{
  interval_score_all_level_temp = matrix(0, nrow = 10, ncol = 59)
  
  # national level
  for(h in 1:10)
  {
    interval_score_all_level_temp[h,1] = interval_score_recon(PI_val = PI_val_list[[cob_ind]][[h]][,1,,], data_series = get(fts_national_fix_cob[cob_ind]), fh = h, alpha = 0.8)$score
  }
  
  # region level
  for(k in 1:11)
  {
    for(h in 1:10)
    {
      interval_score_all_level_temp[h,k+1] = interval_score_recon(PI_val = PI_val_list[[cob_ind]][[h]][,k+1,,], data_series = get(fts_region_fix_cob[[k]][cob_ind]), fh = h, alpha = 0.8)$score
    }
  }
  
  # area level
  for(j in 1:47)
  {
    if(j == 23 && cob_ind == 14)
    {
      interval_score_all_level_temp[,(j+12)] = rep(0, 10)
    } else {
      for(h in 1:10)
      {
        interval_score_all_level_temp[h,(j+12)] = interval_score_recon(PI_val = PI_val_list[[cob_ind]][[h]][,j+12,,], data_series = get(fts_area_fix_cob[[j]][cob_ind]), fh = h, alpha = 0.8)$score
      }
    }
  }
  return(list(interval_score_all = interval_score_all_level_temp))
}

BU_recon_interval_score_fix_cob = OLS_recon_interval_score_fix_cob = MINT_recon_interval_score_fix_cob = list()
for(cob in 1:19)
{
  # BU
  BU_recon_interval_score_fix_cob[[cob]] = interval_score_recon_fix_cob(PI_val_list = BU_recon_PI_fix_cob, cob_ind = cob)$interval_score_all
    
  # OLS
  OLS_recon_interval_score_fix_cob[[cob]] = interval_score_recon_fix_cob(PI_val_list = OLS_recon_PI_fix_cob, cob_ind = cob)$interval_score_all
  
  # MINT
  MINT_recon_interval_score_fix_cob[[cob]] = interval_score_recon_fix_cob(PI_val_list = MINT_recon_PI_fix_cob, cob_ind = cob)$interval_score_all
}


# save(BU_recon_PI_fix_cob, file = "BU_recon_PI_fix_cob.RData")
# save(OLS_recon_PI_fix_cob, file = "OLS_recon_PI_fix_cob.RData")
# save(MINT_recon_PI_fix_cob, file = "MINT_recon_PI_fix_cob.RData")

rm(BU_recon_PI_fix_cob, OLS_recon_PI_fix_cob, MINT_recon_PI_fix_cob)

load(file = "BU_recon_PI_fix_cob.RData")
load(file = "OLS_recon_PI_fix_cob.RData")
load(file = "MINT_recon_PI_fix_cob.RData")

# summary of interval scores

BU_IS_fix_cob = OLS_IS_fix_cob = MINT_IS_fix_cob = array(0, dim = c(10, 3, 19))
for(cob in 1:19)
{
  # BU
  BU_IS_fix_cob[,1,cob] = BU_recon_interval_score_fix_cob[[cob]][,1]
  BU_IS_fix_cob[,2,cob] = rowMeans(BU_recon_interval_score_fix_cob[[cob]][,2:12])
  BU_IS_fix_cob[,3,cob] = rowMeans(BU_recon_interval_score_fix_cob[[cob]][,13:59])
  
  # OP
  OLS_IS_fix_cob[,1,cob] = OLS_recon_interval_score_fix_cob[[cob]][,1]
  OLS_IS_fix_cob[,2,cob] = rowMeans(OLS_recon_interval_score_fix_cob[[cob]][,2:12])
  OLS_IS_fix_cob[,3,cob] = rowMeans(OLS_recon_interval_score_fix_cob[[cob]][,13:59])

  # MINT
  MINT_IS_fix_cob[,1,cob] = MINT_recon_interval_score_fix_cob[[cob]][,1]
  MINT_IS_fix_cob[,2,cob] = rowMeans(MINT_recon_interval_score_fix_cob[[cob]][,2:12])
  MINT_IS_fix_cob[,3,cob] = rowMeans(MINT_recon_interval_score_fix_cob[[cob]][,13:59])
}

BU_IS_fix_cob_averaged = apply(BU_IS_fix_cob, c(1,2), mean)
OLS_IS_fix_cob_averaged = apply(OLS_IS_fix_cob, c(1,2), mean)
MINT_IS_fix_cob_averaged = apply(MINT_IS_fix_cob, c(1,2), mean)


##########
# Fix GEO
##########

COB_14_PI_fh_fix_GEO[[23]] = C_group_6_PI_fh_fix_GEO[[23]] = list()
for(h in 1:9)
{
  COB_14_PI_fh_fix_GEO[[23]][[h]] = C_group_6_PI_fh_fix_GEO[[23]][[h]] = array(rep(1e-6, 7000*(11-h)), dim = c(7, 1000, (11-h)))
}
COB_14_PI_fh_fix_GEO[[23]][[10]] = C_group_6_PI_fh_fix_GEO[[23]][[10]] = matrix(rep(1e-6, 7000), nrow = 7, ncol = 1000)

PI_all_level_fix_geo <- function(fh, geo)
{
  all_level_fh_PI_fix_geo = array(0, dim = c(7, (11-fh), 30, 1000))
  
  if(fh == 10)
  {
    for(ij in 1:1000)
    {
      all_level_fh_PI_fix_geo[,,,ij] = cbind(ind_global_PI_all_fix_GEO[[geo]][[fh]]$PI_boot[,ij,],
                                             
                                             COB_1_PI_fh_fix_GEO[[geo]][[fh]][,ij], COB_2_PI_fh_fix_GEO[[geo]][[fh]][,ij], COB_3_PI_fh_fix_GEO[[geo]][[fh]][,ij], 
                                             COB_4_PI_fh_fix_GEO[[geo]][[fh]][,ij], COB_5_PI_fh_fix_GEO[[geo]][[fh]][,ij], COB_6_PI_fh_fix_GEO[[geo]][[fh]][,ij],
                                             COB_7_PI_fh_fix_GEO[[geo]][[fh]][,ij], COB_8_PI_fh_fix_GEO[[geo]][[fh]][,ij], COB_9_PI_fh_fix_GEO[[geo]][[fh]][,ij],
                                             COB_10_PI_fh_fix_GEO[[geo]][[fh]][,ij], COB_11_PI_fh_fix_GEO[[geo]][[fh]][,ij], COB_12_PI_fh_fix_GEO[[geo]][[fh]][,ij],
                                             COB_13_PI_fh_fix_GEO[[geo]][[fh]][,ij], COB_14_PI_fh_fix_GEO[[geo]][[fh]][,ij], COB_15_PI_fh_fix_GEO[[geo]][[fh]][,ij],
                                             COB_16_PI_fh_fix_GEO[[geo]][[fh]][,ij], COB_17_PI_fh_fix_GEO[[geo]][[fh]][,ij], COB_18_PI_fh_fix_GEO[[geo]][[fh]][,ij],
                                             COB_19_PI_fh_fix_GEO[[geo]][[fh]][,ij],
                                             
                                             C_group_1_PI_fh_fix_GEO[[geo]][[fh]][,ij], C_group_2_PI_fh_fix_GEO[[geo]][[fh]][,ij], C_group_3_PI_fh_fix_GEO[[geo]][[fh]][,ij],
                                             C_group_4_PI_fh_fix_GEO[[geo]][[fh]][,ij], C_group_5_PI_fh_fix_GEO[[geo]][[fh]][,ij], C_group_6_PI_fh_fix_GEO[[geo]][[fh]][,ij],
                                             C_group_7_PI_fh_fix_GEO[[geo]][[fh]][,ij], C_group_8_PI_fh_fix_GEO[[geo]][[fh]][,ij], C_group_9_PI_fh_fix_GEO[[geo]][[fh]][,ij],
                                             C_group_10_PI_fh_fix_GEO[[geo]][[fh]][,ij])
    }
  } else {
    for(ij in 1:1000)
    {
      all_level_fh_PI_fix_geo[,,,ij] = cbind(ind_global_PI_all_fix_GEO[[geo]][[fh]]$PI_boot[,ij,],
                                             
                                             COB_1_PI_fh_fix_GEO[[geo]][[fh]][,ij,], COB_2_PI_fh_fix_GEO[[geo]][[fh]][,ij,], COB_3_PI_fh_fix_GEO[[geo]][[fh]][,ij,], 
                                             COB_4_PI_fh_fix_GEO[[geo]][[fh]][,ij,], COB_5_PI_fh_fix_GEO[[geo]][[fh]][,ij,], COB_6_PI_fh_fix_GEO[[geo]][[fh]][,ij,],
                                             COB_7_PI_fh_fix_GEO[[geo]][[fh]][,ij,], COB_8_PI_fh_fix_GEO[[geo]][[fh]][,ij,], COB_9_PI_fh_fix_GEO[[geo]][[fh]][,ij,],
                                             COB_10_PI_fh_fix_GEO[[geo]][[fh]][,ij,], COB_11_PI_fh_fix_GEO[[geo]][[fh]][,ij,], COB_12_PI_fh_fix_GEO[[geo]][[fh]][,ij,],
                                             COB_13_PI_fh_fix_GEO[[geo]][[fh]][,ij,], COB_14_PI_fh_fix_GEO[[geo]][[fh]][,ij,], COB_15_PI_fh_fix_GEO[[geo]][[fh]][,ij,],
                                             COB_16_PI_fh_fix_GEO[[geo]][[fh]][,ij,], COB_17_PI_fh_fix_GEO[[geo]][[fh]][,ij,], COB_18_PI_fh_fix_GEO[[geo]][[fh]][,ij,],
                                             COB_19_PI_fh_fix_GEO[[geo]][[fh]][,ij,],
                                             
                                             C_group_1_PI_fh_fix_GEO[[geo]][[fh]][,ij,], C_group_2_PI_fh_fix_GEO[[geo]][[fh]][,ij,], C_group_3_PI_fh_fix_GEO[[geo]][[fh]][,ij,],
                                             C_group_4_PI_fh_fix_GEO[[geo]][[fh]][,ij,], C_group_5_PI_fh_fix_GEO[[geo]][[fh]][,ij,], C_group_6_PI_fh_fix_GEO[[geo]][[fh]][,ij,],
                                             C_group_7_PI_fh_fix_GEO[[geo]][[fh]][,ij,], C_group_8_PI_fh_fix_GEO[[geo]][[fh]][,ij,], C_group_9_PI_fh_fix_GEO[[geo]][[fh]][,ij,],
                                             C_group_10_PI_fh_fix_GEO[[geo]][[fh]][,ij,])
    }
  }
  
  return(all_level_fh_PI_fix_geo)
}


test = PI_all_level_fix_geo(fh = 1, geo = 23)

# one- to 10-step-ahead bootstrapped base forecasts (B = 1000)

PI_base_fix_geo = vector("character", length = 59)
for(i in 1:59)
{
  PI_base_fix_geo[i] = paste("PI_base_", i, "_fix_GEO", sep = "")
}

for(geo in 1:59)
{
  temp_list = list()
  for(h in 1:10)
  {
    temp_list[[h]] = PI_all_level_fix_geo(fh = h, geo = geo)
  }
  
  assign(PI_base_fix_geo[geo], temp_list)
  rm(temp_list)
}

# All-level bootstrapped grouped forecasts

recon_PI_fix_geo <- function(kj, age, hier_method = c("BU", "comb_OLS", "mint"), geo_ind)
{
  hier_method = match.arg(hier_method)
  
  age_ind = match(age, age_all)
  
  hier_fore = array(NA, dim = c(30, (11-kj), 1000)) # in total 30 series
  
  # Summing matrix for kj horizon at a given age
  summ_mat = Smat_fun_fix_geo(kj = kj, age = age, no_cob = 19, geo = geo_ind)
  
  # ik: number of years in the forecasting period, it depends on kj
  for(ik in 1:(11-kj))
  {
    hier = summ_mat[,,ik]
    
    if(hier_method == "BU")
    {
      for(ij in 1:1000)
      {
        hier_fore[,ik,ij] = hier %*% (get(paste("PI_base_", geo_ind, "_fix_GEO", sep = "")))[[kj]][age_ind,ik,12:30,ij]
      }
    }
    
    if(hier_method == "comb_OLS")
    {
      for(ij in 1:1000)
      {
        hier_fore[,ik,ij] = hier %*% ginv(t(hier) %*% hier) %*% t(hier) %*% (get(paste("PI_base_", geo_ind, "_fix_GEO", sep = "")))[[kj]][age_ind,ik,,ij]
      }
    }
    
    if(hier_method == "mint")
    {
      ginv = MASS:::ginv
      wh = wh_fun_fix_geo(kj = kj, age = age, geo_ind = geo_ind)
      for(ij in 1:1000)
      {
        hier_fore[,ik,ij] = hier %*% ginv(t(hier) %*% ginv(wh) %*% hier) %*% t(hier) %*% ginv(wh) %*% (get(paste("PI_base_", geo_ind, "_fix_GEO", sep = "")))[[kj]][age_ind,ik,,ij]
      }
    }
  }
  return(hier_fore)
}

test = recon_PI_fix_geo(kj = 10, age = age_all[4], hier_method = "mint", geo_ind = 23)


# Calculation of reconciled interval predictions (hier_method = "BU")
BU_recon_PI_fix_geo = list()

for(geo in 1:59)
{
  BU_recon_PI_fix_geo[[geo]] = list()
  for(h in 1:10)
  {
    BU_recon_PI_fix_geo[[geo]][[h]] = array(NA, dim = c(7, 30, (11-h), 1000))
    
    for(ij in 1:7)
    {
      BU_recon_PI_fix_geo[[geo]][[h]][ij,,,] = recon_PI_fix_geo(kj = h, age = age_all[ij], hier_method = "BU", geo_ind = geo)
    }
  }
}


# Calculation of reconciled interval predictions (hier_method = "comb_OLS")
OLS_recon_PI_fix_geo = list()

for(geo in 1:59)
{
  OLS_recon_PI_fix_geo[[geo]] = list()
  for(h in 1:10)
  {
    OLS_recon_PI_fix_geo[[geo]][[h]] = array(NA, dim = c(7, 30, (11-h), 1000))
    
    for(ij in 1:7)
    {
      OLS_recon_PI_fix_geo[[geo]][[h]][ij,,,] = recon_PI_fix_geo(kj = h, age = age_all[ij], hier_method = "comb_OLS", geo_ind = geo)
    }
  }
}

# Calculation of reconciled interval predictions (hier_method = "mint")
MINT_recon_PI_fix_geo = list()

for(geo in 1:59)
{
  MINT_recon_PI_fix_geo[[geo]] = list()
  for(h in 1:10)
  {
    MINT_recon_PI_fix_geo[[geo]][[h]] = array(NA, dim = c(7, 30, (11-h), 1000))
    
    for(ij in 1:7)
    {
      MINT_recon_PI_fix_geo[[geo]][[h]][ij,,,] = recon_PI_fix_geo(kj = h, age = age_all[ij], hier_method = "mint", geo_ind = geo)
    }
  }
}


# function for computing reconciled interval scores
interval_score_recon_fix_geo <- function(PI_val_list, geo_ind)
{
  interval_score_all_level_temp = matrix(0, nrow = 10, ncol = 30)
  
  # Global level
  for(h in 1:10)
  {
    interval_score_all_level_temp[h,1] = interval_score_recon(PI_val = PI_val_list[[geo_ind]][[h]][,1,,], data_series = get(fts_national_fix_geo[geo_ind]), fh = h, alpha = 0.8)$score
  }
  
  
  # C_group level
  for(k in 1:10)
  {
    if(k == 6 && geo == 23)
    {
      interval_score_all_level_temp[,k+1] = rep(0,10)
    } else {
      for(h in 1:10)
      {
        interval_score_all_level_temp[h,k+1] = interval_score_recon(PI_val = PI_val_list[[geo_ind]][[h]][,k+1,,], data_series = get(fts_region_fix_geo[[k]][geo_ind]), fh = h, alpha = 0.8)$score
      }
    }
  }
  
  # COB level
  for(j in 1:19)
  {
    if(j == 14 && geo_ind == 23)
    {
      interval_score_all_level_temp[,(j+11)] = rep(0, 10)
    } else {
      for(h in 1:10)
      {
        interval_score_all_level_temp[h,(j+11)] = interval_score_recon(PI_val = PI_val_list[[geo_ind]][[h]][,j+11,,], data_series = get(fts_area_fix_geo[[j]][geo_ind]), fh = h, alpha = 0.8)$score
      }
    }
  }
  return(list(interval_score_all = interval_score_all_level_temp))
}

BU_recon_interval_score_fix_geo = OLS_recon_interval_score_fix_geo = MINT_recon_interval_score_fix_geo = list()
for(geo in 1:59)
{
  # BU
  BU_recon_interval_score_fix_geo[[geo]] = interval_score_recon_fix_geo(PI_val_list = BU_recon_PI_fix_geo, geo_ind = geo)$interval_score_all
  
  # OLS
  OLS_recon_interval_score_fix_geo[[geo]] = interval_score_recon_fix_geo(PI_val_list = OLS_recon_PI_fix_geo, geo_ind = geo)$interval_score_all
  
  # MINT
  MINT_recon_interval_score_fix_geo[[geo]] = interval_score_recon_fix_geo(PI_val_list = MINT_recon_PI_fix_geo, geo_ind = geo)$interval_score_all
}


# save(BU_recon_PI_fix_geo, file = "BU_recon_PI_fix_geo.RData")
# save(OLS_recon_PI_fix_geo, file = "OLS_recon_PI_fix_geo.RData")
# save(MINT_recon_PI_fix_geo, file = "MINT_recon_PI_fix_geo.RData")

rm(BU_recon_PI_fix_geo, OLS_recon_PI_fix_geo, MINT_recon_PI_fix_geo)

load(file = "BU_recon_PI_fix_geo.RData")
load(file = "OLS_recon_PI_fix_geo.RData")
load(file = "MINT_recon_PI_fix_geo.RData")



# summary of interval scores
BU_IS_fix_geo = OLS_IS_fix_geo = MINT_IS_fix_geo = array(0, dim = c(10, 3, 59))
for(geo in 1:59)
{
  # BU
  BU_IS_fix_geo[,1,geo] = BU_recon_interval_score_fix_geo[[geo]][,1]
  BU_IS_fix_geo[,2,geo] = rowMeans(BU_recon_interval_score_fix_geo[[geo]][,2:11])
  BU_IS_fix_geo[,3,geo] = rowMeans(BU_recon_interval_score_fix_geo[[geo]][,12:30])
  
  # OP
  OLS_IS_fix_geo[,1,geo] = OLS_recon_interval_score_fix_geo[[geo]][,1]
  OLS_IS_fix_geo[,2,geo] = rowMeans(OLS_recon_interval_score_fix_geo[[geo]][,2:11])
  OLS_IS_fix_geo[,3,geo] = rowMeans(OLS_recon_interval_score_fix_geo[[geo]][,12:30])

  # MINT
  MINT_IS_fix_geo[,1,geo] = MINT_recon_interval_score_fix_geo[[geo]][,1]
  MINT_IS_fix_geo[,2,geo] = rowMeans(MINT_recon_interval_score_fix_geo[[geo]][,2:11])
  MINT_IS_fix_geo[,3,geo] = rowMeans(MINT_recon_interval_score_fix_geo[[geo]][,12:30])
}

BU_IS_fix_geo_averaged = apply(BU_IS_fix_geo, c(1,2), mean)
OLS_IS_fix_geo_averaged = apply(OLS_IS_fix_geo, c(1,2), mean)
MINT_IS_fix_geo_averaged = apply(MINT_IS_fix_geo, c(1,2), mean)
