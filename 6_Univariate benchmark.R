######################################################################
# Univariate point forecast and interval forecasts (benchmark method)
#####################################################################

###############################################
# Fix COB; disaggregate by geographical factor
###############################################

ind_recon_hier_fix_cob <- function(national_forc, region_forc, area_forc, age, kj, hier_method = c("BU", "comb_OLS", "comb_GLS", "mint"), cob_ind)
{
  hier_method = match.arg(hier_method)
  
  age_ind = match(age, age_interp)
  
  hier_rate = matrix(NA,59,(11-kj)) # in total 59 series
  
  # Level_1 national
  
  hier_rate[1,] = get(national_forc[1])[[cob_ind]][kj,age_ind,1:(11-kj)] #  arguments: horizon, age, 11- horzon; length of list == cob
  
  
  # Level 2 disaggregate by region
  for(iw in 1:11)
  {
    hier_rate[iw+1,]  = get(region_forc[iw])[[cob_ind]][kj,age_ind,1:(11-kj)]  
  }
  
  # Level 3 disaggregate by area
  
  for(iw in 1:47)
  {
    if(iw == 23)
    {
      hier_rate[iw+12,] = rep(0, 11-kj)
    } else {
      hier_rate[iw+12,] = get(area_forc[iw])[[cob_ind]][kj,age_ind,1:(11-kj)]
    }
  }  
  
  # forecast reconciliation via bottom-up, optimal combination or MinT methods
  hier_fore = matrix(NA, 59, (11-kj))
  summing_mat = Smat_fun_fix_cob(kj = kj, age = age, no_area = 47, cob = cob_ind)
  for(ik in 1:(11-kj))
  {
    hier = summing_mat[,,ik]
    if(hier_method == "BU")
    {
      hier_fore[,ik] = (hier %*% hier_rate[13:59,ik])
    }
    if(hier_method == "comb_OLS")
    {
      hier_fore[,ik] = hier %*% ginv(t(hier) %*% hier) %*% t(hier) %*% hier_rate[,ik]
    }
    if(hier_method == "comb_GLS")
    {
      hier_fore[,ik] = hier %*% ginv(t(hier) %*% ginv(diag(rowSums(hier^2))) %*% hier) %*% t(hier) %*% ginv(diag(rowSums(hier^2))) %*% hier_rate[,ik]
    }
    if(hier_method == "mint")
    {
      ginv = MASS:::ginv
      wh = ind_wh_fun_fix_cob(kj = kj, age = age, cob_ind = cob_ind)
      hier_fore[,ik] = hier %*% ginv(t(hier) %*% ginv(wh) %*% hier) %*% t(hier) %*% ginv(wh) %*% hier_rate[,ik]
    }
  }
  return(hier_fore)
}

ind_BU_optim_err_fix_cob <- function(ik, hier_method, national_forc, region_forc, area_forc, cob_ind)
{
  age_all = c(15, 20, 25, 30, 35, 40, 45)
  
  me = ftsa:::me; mae = ftsa:::mae; rmse = ftsa:::rmse; mase = ftsa:::mase
  BU_optim_hier_comb = array(NA, dim = c(7, 59,(11-ik)))
  
  for(i in 1:7)
  {
    BU_optim_hier_comb[i,,] = ind_recon_hier_fix_cob(national_forc, region_forc, area_forc,      
                                                 age = age_all[i], kj = ik, hier_method = hier_method,
                                                 cob_ind = cob_ind)
  }
  
  #####################################
  # Errors, including ME, MAE and RMSE
  #####################################
  
  # Level 1 (Total)
  mae_national = mae(BU_optim_hier_comb[,1,], extract.years(get(fts_national_fix_cob[cob_ind]), years = (2001+ik):2011)$rate$female)
  mase_national = mase(BU_optim_hier_comb[,1,], extract.years(get(fts_national_fix_cob[cob_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_national_fix_cob[cob_ind]), years = 1986:(2001+ik-1))$rate$female)
  rmse_national = rmse(BU_optim_hier_comb[,1,], extract.years(get(fts_national_fix_cob[cob_ind]), years = (2001+ik):2011)$rate$female)
  
  # Level 2 (Region)
  mae_region = mase_region= rmse_region = rep(0, 11)
  for(iw in 1:11)
  {
    mae_region[iw]  = mae(BU_optim_hier_comb[,(iw+1),],  extract.years(get(fts_region_fix_cob[[iw]][cob_ind]), years = (2001+ik):2011)$rate$female)
    mase_region[iw]  = mase(BU_optim_hier_comb[,(iw+1),],  extract.years(get(fts_region_fix_cob[[iw]][cob_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_region_fix_cob[[iw]][cob_ind]), years = 1986:(2001+ik-1))$rate$female)
    rmse_region[iw] = rmse(BU_optim_hier_comb[,(iw+1),], extract.years(get(fts_region_fix_cob[[iw]][cob_ind]), years = (2001+ik):2011)$rate$female)
  }
  
  # Level 3 (Area)
  mae_area = mase_area = mase_area = rmse_area = rep(0, 47)
  
  for(iwk in 1:47)
  {
    if(iwk != 23)
    {
      mae_area[iwk]  = mae(BU_optim_hier_comb[,(iwk+12),],   extract.years(get(fts_area_fix_cob[[iwk]][cob_ind]), years = (2001+ik):2011)$rate$female)
      mase_area[iwk]  = mase(BU_optim_hier_comb[,(iwk+12),],   extract.years(get(fts_area_fix_cob[[iwk]][cob_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_area_fix_cob[[iwk]][cob_ind]), years = 1986:(2001+ik-1))$rate$female)
      rmse_area[iwk] = rmse(BU_optim_hier_comb[,(iwk+12),],  extract.years(get(fts_area_fix_cob[[iwk]][cob_ind]), years = (2001+ik):2011)$rate$female)
    }
  }
  
  
  return(list(mae_national = mae_national, 
              mae_region = mae_region, 
              mae_area  = mae_area, 
              
              mase_national = mase_national, 
              mase_region = mase_region, 
              mase_area  = mase_area, 
              
              rmse_national = rmse_national, 
              rmse_region = rmse_region, 
              rmse_area  = rmse_area))
}



# Area level 

ind_area_train_residual_fix_COB = vector("character", length = 47)
for(i in 1:47)
{
  ind_area_train_residual_fix_COB[i] = paste("ind_A", i, "_train_residual_fix_COB", sep = "")
}

ind_area_forecast_fix_COB = vector("character", length = 47)
for(i in 1:47)
{
  ind_area_forecast_fix_COB[i] = paste("ind_A", i, "_forecast_fix_COB", sep = "")
}

ind_area_mae_fix_COB = vector("character", length = 47)
for(i in 1:47)
{
  ind_area_mae_fix_COB[i] = paste("ind_A", i, "_mae_fix_COB", sep = "")
}

ind_area_mase_fix_COB = vector("character", length = 47)
for(i in 1:47)
{
  ind_area_mase_fix_COB[i] = paste("ind_A", i, "_mase_fix_COB", sep = "")
}


ind_area_rmse_fix_COB = vector("character", length = 47)
for(i in 1:47)
{
  ind_area_rmse_fix_COB[i] = paste("ind_A", i, "_rmse_fix_COB", sep = "")
}

# Region level

ind_region_train_residual_fix_COB = vector("character", length = 11)
for(i in 1:11)
{
  ind_region_train_residual_fix_COB[i] = paste("ind_R", i, "_train_residual_fix_COB", sep = "")
}

ind_region_forecast_fix_COB = vector("character", length = 11)
for(i in 1:11)
{
  ind_region_forecast_fix_COB[i] = paste("ind_R", i, "_forecast_fix_COB", sep = "")
}

ind_region_mae_fix_COB = vector("character", length = 11)
for(i in 1:11)
{
  ind_region_mae_fix_COB[i] = paste("ind_R", i, "_mae_fix_COB", sep = "")
}

ind_region_mase_fix_COB = vector("character", length = 11)
for(i in 1:11)
{
  ind_region_mase_fix_COB[i] = paste("ind_R", i, "_mase_fix_COB", sep = "")
}

ind_region_rmse_fix_COB = vector("character", length = 11)
for(i in 1:11)
{
  ind_region_rmse_fix_COB[i] = paste("ind_R", i, "_rmse_fix_COB", sep = "")
}

# National level

ind_national_train_residual_fix_COB = paste("ind_National", "_train_residual_fix_COB", sep = "")
ind_national_forecast_fix_COB = paste("ind_National", "_forcast_fix_COB", sep = "")
ind_national_mae_fix_COB = paste("ind_National", "_mae_fix_COB", sep = "")
ind_national_mase_fix_COB = paste("ind_National", "_mase_fix_COB", sep = "")
ind_national_rmse_fix_COB = paste("ind_National", "_rmse_fix_COB", sep = "")

# Area independent forecasts

ind_area_fix_COB <- function(area_index = 1,  fmethod = c("classical", "M"), year_horizon, cob_ind)
{
  require(ftsa)
  require(demography)
  fmethod = match.arg(fmethod)
  
  n_year = 31 - (year_horizon + 1)
  
  if(area_index == 23 && cob_ind == 14)
  {
    res_area = train_residual = area_mae = area_mase = area_rmse = 0  # no observation at all for Area 23 COB 14
  } else {
    area_comb = get(fts_area_fix_cob_smooth[[area_index]][cob_ind])$rate$female[3:33,]
    colnames(area_comb) = year_interp
    rownames(area_comb) = age_interp
    
    res_area = array(NA, dim = c(year_horizon, 31, year_horizon))
    train_residual = list()
    for(ij in 1:year_horizon)
    {
      ind_dat = fts(x = age_interp, y = area_comb[,1:(n_year+ij)], xname = "age", yname = "fertility")
      ftsm_order = head(which(round(cumsum(ftsm(ind_dat, method = fmethod, lambda = 3, order = ncol(ind_dat$y))$varprop),3)>=0.95),1)
      ftsm_area = ftsm(ind_dat, order = ftsm_order, method = fmethod, lambda = 3)
      
      fun_forc  = forecast(ftsm_area, h = year_horizon)
      
      res_area[,,ij] = t(fun_forc$mean$y)
      
      train_residual[[ij]] = fun_forc$error$y
    }
    
    # Errors
    
    area_mae = area_mase = area_rmse = matrix(NA, year_horizon, 1)
    age_ind = match(age_all, age_interp) 
    
    for(ik in 1:year_horizon)
    {
      area_mae[ik,1] = ftsa:::mae(res_area[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_area_fix_cob[[area_index]][cob_ind]), years = (2001+ik):2011)$rate$female)
      area_mase[ik,1] = ftsa:::mase(res_area[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_area_fix_cob[[area_index]][cob_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_area_fix_cob[[area_index]][cob_ind]), years = 1986:(2001+ik-1))$rate$female)
      area_rmse[ik,1] = ftsa:::rmse(res_area[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_area_fix_cob[[area_index]][cob_ind]), years = (2001+ik):2011)$rate$female)
      
    }
  }
  
  return(list(res_area = res_area, train_residual = train_residual, area_mae = area_mae, area_rmse = area_rmse, area_mase = area_mase))    
}

test = ind_area_fix_COB(area_index = 5, fmethod = "classical", year_horizon = 10, cob_ind = 1)

for(area_ind in 1:47)
{
  temp_res_forc_list = temp_train_residual_list = temp_mae_list = temp_mase_list = temp_rmse_list = list()
  for(cob in 1:19)
  {
    temp = ind_area_fix_COB(area_index = area_ind, fmethod = "classical", year_horizon = 10, cob_ind = cob)
    temp_res_forc_list[[cob]] = temp$res_area
    temp_train_residual_list[[cob]] = temp$train_residual
    temp_mae_list[[cob]] = temp$area_mae
    temp_mase_list[[cob]] = temp$area_mase
    temp_rmse_list[[cob]] = temp$area_rmse
  }
  
  assign(ind_area_forecast_fix_COB[area_ind], temp_res_forc_list)
  assign(ind_area_train_residual_fix_COB[area_ind], temp_train_residual_list)
  assign(ind_area_mae_fix_COB[area_ind], temp_mae_list)
  assign(ind_area_mase_fix_COB[area_ind], temp_mase_list)
  assign(ind_area_rmse_fix_COB[area_ind], temp_rmse_list)
  
  rm(temp, temp_res_forc_list, temp_train_residual_list, temp_mae_list, temp_mase_list, temp_rmse_list)
  
  print(area_ind)
}


ind_fix_cob_all_area_mae = array(0, dim = c(10, 47, 19))
for(ik in 1:47)
{
  temp = get(ind_area_mae_fix_COB[ik])
  for(ij in 1:19)
  {
    ind_fix_cob_all_area_mae[,ik,ij] = temp[[ij]]
  }
}


ind_fix_cob_all_area_mase = array(0, dim = c(10, 47, 19))
for(ik in 1:47)
{
  temp = get(ind_area_mase_fix_COB[ik])
  for(ij in 1:19)
  {
    ind_fix_cob_all_area_mase[,ik,ij] = temp[[ij]]
  }
}


ind_fix_cob_all_area_rmse = array(0, dim = c(10, 47, 19))
for(ik in 1:47)
{
  temp = get(ind_area_rmse_fix_COB[ik])
  for(ij in 1:19)
  {
    ind_fix_cob_all_area_rmse[,ik,ij] = temp[[ij]]
  }
}


# Region independent forecast

ind_region_fix_COB <- function(region_index = 1,  fmethod = c("classical", "M"), year_horizon, cob_ind)
{
  require(ftsa)
  require(demography)
  fmethod = match.arg(fmethod)
  
  n_year = 31 - (year_horizon + 1)
  
  region_comb = get(fts_region_fix_cob_smooth[[region_index]][cob_ind])$rate$female[3:33,]
  colnames(region_comb) = year_interp
  rownames(region_comb) = age_interp
  
  res_region = array(NA, dim = c(year_horizon, 31, year_horizon))
  train_residual = list()
  for(ij in 1:year_horizon)
  {
    ind_dat = fts(x = age_interp, y = region_comb[,1:(n_year+ij)], xname = "age", yname = "fertility")
    ftsm_order = head(which(round(cumsum(ftsm(ind_dat, method = fmethod, lambda = 3, order = ncol(ind_dat$y))$varprop),3)>=0.95),1)
    ftsm_region = ftsm(ind_dat, order = ftsm_order, method = fmethod, lambda = 3)
    
    fun_forc  = forecast(ftsm_region, h = year_horizon)
    
    res_region[,,ij] = t(fun_forc$mean$y)
    
    train_residual[[ij]] = fun_forc$error$y
  }
  
  # Errors
  
  region_mae = region_mase = region_rmse = matrix(NA, year_horizon, 1)
  age_ind = match(age_all, age_interp) 
  
  for(ik in 1:year_horizon)
  {
    region_mae[ik,1] = ftsa:::mae(res_region[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_region_fix_cob[[region_index]][cob_ind]), years = (2001+ik):2011)$rate$female)
    region_mase[ik,1] = ftsa:::mase(res_region[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_region_fix_cob[[region_index]][cob_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_region_fix_cob[[region_index]][cob_ind]), years = 1986:(2001+ik-1))$rate$female)
    region_rmse[ik,1] = ftsa:::rmse(res_region[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_region_fix_cob[[region_index]][cob_ind]), years = (2001+ik):2011)$rate$female)
  }
  
  return(list(res_region = res_region, train_residual = train_residual, region_mae = region_mae, region_mase = region_mase, region_rmse = region_rmse))    
}

test = ind_region_fix_COB(region_index = 5, fmethod = "classical", year_horizon = 10, cob_ind = 1)

for(region_ind in 1:11)
{
  temp_res_forc_list = temp_train_residual_list = temp_mae_list = temp_mase_list = temp_rmse_list = list()
  for(cob in 1:19)
  {
    temp = ind_region_fix_COB(region_index = region_ind, fmethod = "classical", year_horizon = 10, cob_ind = cob)
    temp_res_forc_list[[cob]] = temp$res_region
    temp_train_residual_list[[cob]] = temp$train_residual
    temp_mae_list[[cob]] = temp$region_mae
    temp_mase_list[[cob]] = temp$region_mase
    temp_rmse_list[[cob]] = temp$region_rmse
  }
  
  assign(ind_region_forecast_fix_COB[region_ind], temp_res_forc_list)
  assign(ind_region_train_residual_fix_COB[region_ind], temp_train_residual_list)
  assign(ind_region_mae_fix_COB[region_ind], temp_mae_list)
  assign(ind_region_mase_fix_COB[region_ind], temp_mase_list)
  assign(ind_region_rmse_fix_COB[region_ind], temp_rmse_list)
  
  rm(temp, temp_res_forc_list, temp_train_residual_list, temp_mae_list, temp_mase_list, temp_rmse_list)
}

ind_fix_cob_all_region_mae = array(0, dim = c(10, 11, 19))
for(ik in 1:11)
{
  temp = get(ind_region_mae_fix_COB[ik])
  for(ij in 1:19)
  {
    ind_fix_cob_all_region_mae[,ik,ij] = temp[[ij]]
  }
}

ind_fix_cob_all_region_mase = array(0, dim = c(10, 11, 19))
for(ik in 1:11)
{
  temp = get(ind_region_mase_fix_COB[ik])
  for(ij in 1:19)
  {
    ind_fix_cob_all_region_mase[,ik,ij] = temp[[ij]]
  }
}

ind_fix_cob_all_region_rmse = array(0, dim = c(10, 11, 19))
for(ik in 1:11)
{
  temp = get(ind_region_rmse_fix_COB[ik])
  for(ij in 1:19)
  {
    ind_fix_cob_all_region_rmse[,ik,ij] = temp[[ij]]
  }
}

# National independent forecast

ind_national_fix_COB <- function(fmethod = c("classical", "M"), year_horizon, cob_ind)
{
  require(ftsa)
  require(demography)
  fmethod = match.arg(fmethod)
  
  n_year = 31 - (year_horizon + 1)
  
  national_comb = get(fts_national_fix_cob_smooth[cob_ind])$rate$female[3:33,]
  colnames(national_comb) = year_interp
  rownames(national_comb) = age_interp
  
  res_national = array(NA, dim = c(year_horizon, 31, year_horizon))
  train_residual = list()
  for(ij in 1:year_horizon)
  {
    ind_dat = fts(x = age_interp, y = national_comb[,1:(n_year+ij)], xname = "age", yname = "fertility")
    ftsm_order = head(which(round(cumsum(ftsm(ind_dat, method = fmethod, lambda = 3, order = ncol(ind_dat$y))$varprop),3)>=0.95),1)
    ftsm_national = ftsm(ind_dat, order = ftsm_order, method = fmethod, lambda = 3)
    
    fun_forc  = forecast(ftsm_national, h = year_horizon)
    
    res_national[,,ij] = t(fun_forc$mean$y)
    
    train_residual[[ij]] = fun_forc$error$y
  }
  
  # Errors
  
  national_mae = national_mase = national_rmse = matrix(NA, year_horizon, 1)
  age_ind = match(age_all, age_interp) 
  
  for(ik in 1:year_horizon)
  {
    national_mae[ik,1] = ftsa:::mae(res_national[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_national_fix_cob[cob_ind]), years = (2001+ik):2011)$rate$female)
    national_mase[ik,1] = ftsa:::mase(res_national[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_national_fix_cob[cob_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_national_fix_cob[cob_ind]), years = 1986:(2001+ik-1))$rate$female)
    national_rmse[ik,1] = ftsa:::rmse(res_national[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_national_fix_cob[cob_ind]), years = (2001+ik):2011)$rate$female)
    
  }
  
  return(list(res_national = res_national, train_residual = train_residual, national_mae = national_mae, national_mase = national_mase, national_rmse = national_rmse))    
}


ind_national_train_residual_list = ind_national_forcast_list = ind_national_mae_list = ind_national_mase_list= ind_national_rmse_list = list()
for(cob in 1:19)
{
  temp = ind_national_fix_COB(fmethod = "classical", year_horizon = 10, cob_ind = cob)
  ind_national_forcast_list[[cob]] = temp$res_national
  ind_national_train_residual_list[[cob]] = temp$train_residual
  ind_national_mae_list[[cob]] = temp$national_mae
  ind_national_mase_list[[cob]] = temp$national_mase
  ind_national_rmse_list[[cob]] = temp$national_rmse
}
rm(temp)

assign(ind_national_forecast_fix_COB[1], ind_national_forcast_list)
assign(ind_national_train_residual_fix_COB[1], ind_national_train_residual_list)
assign(ind_national_mae_fix_COB[1], ind_national_mae_list)
assign(ind_national_mase_fix_COB[1], ind_national_mase_list)
assign(ind_national_rmse_fix_COB[1], ind_national_rmse_list)

ind_fix_cob_all_national_mae = matrix(0, nrow = 10, ncol = 19)
for(ij in 1:19)
{
  ind_fix_cob_all_national_mae[,ij] = get(ind_national_mae_fix_COB[1])[[ij]]
}

ind_fix_cob_all_national_mase = matrix(0, nrow = 10, ncol = 19)
for(ij in 1:19)
{
  ind_fix_cob_all_national_mase[,ij] = get(ind_national_mase_fix_COB[1])[[ij]]
}

ind_fix_cob_all_national_rmse = matrix(0, nrow = 10, ncol = 19)
for(ij in 1:19)
{
  ind_fix_cob_all_national_rmse[,ij] = get(ind_national_rmse_fix_COB[1])[[ij]]
}

# combine all point forecasts into an array

ind_fix_cob_all_level_mae = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  ind_fix_cob_all_level_mae[,,i] = cbind(ind_fix_cob_all_national_mae[,i], rowMeans(ind_fix_cob_all_region_mae[,,i]), rowMeans(ind_fix_cob_all_area_mae[,,i]))
}
ind_fix_cob_averaged_mae = apply(ind_fix_cob_all_level_mae, c(1,2), mean)


###
ind_fix_cob_all_area_mase_copy = ind_fix_cob_all_area_mase
for(i in 1:19)
{
  ind_fix_cob_all_area_mase[,,i][,which(colMeans(ind_fix_cob_all_area_mase[,,i]) > 300)] = 0
}
###


ind_fix_cob_all_level_mase = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  ind_fix_cob_all_level_mase[,,i] = cbind(ind_fix_cob_all_national_mase[,i], rowMeans(ind_fix_cob_all_region_mase[,,i]), rowMeans(ind_fix_cob_all_area_mase[,,i]))
}
ind_fix_cob_averaged_mase = apply(ind_fix_cob_all_level_mase, c(1,2), mean)


ind_fix_cob_all_level_rmse = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  ind_fix_cob_all_level_rmse[,,i] = cbind(ind_fix_cob_all_national_rmse[,i], rowMeans(ind_fix_cob_all_region_rmse[,,i]), rowMeans(ind_fix_cob_all_area_rmse[,,i]))
}
ind_fix_cob_averaged_rmse = apply(ind_fix_cob_all_level_rmse, c(1,2), mean)


###################################################
# Fix GEO; disaggregate by Country/Region of birth
###################################################


# COB level 

ind_cob_train_residual_fix_GEO = vector("character", length = 19)
for(i in 1:19)
{
  ind_cob_train_residual_fix_GEO[i] = paste("ind_cob", i, "_train_residual_fix_GEO", sep = "")
}

ind_cob_forecast_fix_GEO = vector("character", length = 19)
for(i in 1:19)
{
  ind_cob_forecast_fix_GEO[i] = paste("ind_cob", i, "_forecast_fix_GEO", sep = "")
}

ind_cob_mae_fix_GEO = vector("character", length = 19)
for(i in 1:19)
{
  ind_cob_mae_fix_GEO[i] = paste("ind_cob", i, "_mae_fix_GEO", sep = "")
}

ind_cob_mase_fix_GEO = vector("character", length = 19)
for(i in 1:19)
{
  ind_cob_mase_fix_GEO[i] = paste("ind_cob", i, "_mase_fix_GEO", sep = "")
}

ind_cob_rmse_fix_GEO = vector("character", length = 19)
for(i in 1:19)
{
  ind_cob_rmse_fix_GEO[i] = paste("ind_cob", i, "_rmse_fix_GEO", sep = "")
}

# C_group level

ind_c_group_train_residual_fix_GEO = vector("character", length = 10)
for(i in 1:10)
{
  ind_c_group_train_residual_fix_GEO[i] = paste("ind_C", i, "_train_residual_fix_GEO", sep = "")
}

ind_c_group_forecast_fix_GEO = vector("character", length = 10)
for(i in 1:10)
{
  ind_c_group_forecast_fix_GEO[i] = paste("ind_C", i, "_forecast_fix_GEO", sep = "")
}

ind_c_group_mae_fix_GEO = vector("character", length = 10)
for(i in 1:10)
{
  ind_c_group_mae_fix_GEO[i] = paste("ind_C", i, "_mae_fix_GEO", sep = "")
}

ind_c_group_mase_fix_GEO = vector("character", length = 10)
for(i in 1:10)
{
  ind_c_group_mase_fix_GEO[i] = paste("ind_C", i, "_mase_fix_GEO", sep = "")
}

ind_c_group_rmse_fix_GEO = vector("character", length = 10)
for(i in 1:10)
{
  ind_c_group_rmse_fix_GEO[i] = paste("ind_C", i, "_rmse_fix_GEO", sep = "")
}

# global

ind_global_train_residual_fix_GEO = paste("ind_Global", "_train_residual_fix_GEO", sep = "")
ind_global_forecast_fix_GEO = paste("ind_Global", "_forecast_fix_GEO", sep = "")
ind_global_mae_fix_GEO = paste("ind_Global", "_mae_fix_GEO", sep = "")
ind_global_mase_fix_GEO = paste("ind_Global", "_mase_fix_GEO", sep = "")
ind_global_rmse_fix_GEO = paste("ind_Global", "_rmse_fix_GEO", sep = "")


# COB independent forecasts

ind_cob_fix_GEO <- function(cob_index = 1,  fmethod = c("classical", "M"), year_horizon, geo_ind)
{
  require(ftsa)
  require(demography)
  fmethod = match.arg(fmethod)
  
  n_year = 31 - (year_horizon + 1)
  
  if(cob_index == 14 && geo_ind == 23)
  {
    res_cob = train_residual = cob_mae = cob_mase = cob_rmse = 0  # no observation at all for cob 14 area 23
  } else {
    cob_comb = get(fts_area_fix_geo_smooth[[cob_index]][geo_ind])$rate$female[3:33,]
    colnames(cob_comb) = year_interp
    rownames(cob_comb) = age_interp
    
    res_cob = array(NA, dim = c(year_horizon, 31, year_horizon))
    train_residual = list()
    for(ij in 1:year_horizon)
    {
      ind_dat = fts(x = age_interp, y = cob_comb[,1:(n_year+ij)], xname = "age", yname = "fertility")
      ftsm_order = head(which(round(cumsum(ftsm(ind_dat, method = fmethod, lambda = 3, order = ncol(ind_dat$y))$varprop),3)>=0.95),1)
      ftsm_cob = ftsm(ind_dat, order = ftsm_order, method = fmethod, lambda = 3)
      
      fun_forc  = forecast(ftsm_cob, h = year_horizon)
      
      res_cob[,,ij] = t(fun_forc$mean$y)
      
      train_residual[[ij]] = fun_forc$error$y
    }
    
    # Errors
    
    cob_mae = cob_mase = cob_rmse = matrix(NA, year_horizon, 1)
    age_ind = match(age_all, age_interp) 
    
    for(ik in 1:year_horizon)
    {
      cob_mae[ik,1] = ftsa:::mae(res_cob[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_area_fix_geo[[cob_index]][geo_ind]), years = (2001+ik):2011)$rate$female)
      cob_mase[ik,1] = ftsa:::mase(res_cob[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_area_fix_geo[[cob_index]][geo_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_area_fix_geo[[cob_index]][geo_ind]), years = 1986:(2001+ik-1))$rate$female)
      cob_rmse[ik,1] = ftsa:::rmse(res_cob[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_area_fix_geo[[cob_index]][geo_ind]), years = (2001+ik):2011)$rate$female)
    }
  }
  
  return(list(res_cob = res_cob, train_residual = train_residual, cob_mae = cob_mae, cob_mase = cob_mase, cob_rmse = cob_rmse))    
}

test = ind_cob_fix_GEO(cob_index = 1, fmethod = "classical", year_horizon = 10, geo_ind = 59)

for(cob_ind in 1:19)
{
  temp_res_forc_list = temp_train_residual_list = temp_mae_list = temp_mase_list = temp_rmse_list = list()
  for(geo in 1:59)
  {
    temp = ind_cob_fix_GEO(cob_index = cob_ind, fmethod = "classical", year_horizon = 10, geo_ind = geo)
    temp_res_forc_list[[geo]] = temp$res_cob
    temp_train_residual_list[[geo]] = temp$train_residual
    temp_mae_list[[geo]] = temp$cob_mae
    temp_mase_list[[geo]] = temp$cob_mase
    temp_rmse_list[[geo]] = temp$cob_rmse
  }
  
  assign(ind_cob_forecast_fix_GEO[cob_ind], temp_res_forc_list)
  assign(ind_cob_train_residual_fix_GEO[cob_ind], temp_train_residual_list)
  assign(ind_cob_mae_fix_GEO[cob_ind], temp_mae_list)
  assign(ind_cob_mase_fix_GEO[cob_ind], temp_mase_list)
  assign(ind_cob_rmse_fix_GEO[cob_ind], temp_rmse_list)
  
  rm(temp, temp_res_forc_list, temp_train_residual_list, temp_mae_list, temp_mase_list, temp_rmse_list)
  
  print(cob_ind)
}

# summary of results

ind_fix_geo_all_cob_mae = array(0, dim = c(10, 19, 59))
for(ik in 1:19)
{
  temp = get(ind_cob_mae_fix_GEO[ik])
  for(ij in 1:59)
  {
    if(ik == 14 && ij == 23)
    {
      ind_fix_geo_all_cob_mae[,ik,ij] = rep(0, 10)
    } else {
      ind_fix_geo_all_cob_mae[,ik,ij] = temp[[ij]]
    }
  }
}

ind_fix_geo_all_cob_mase = array(0, dim = c(10, 19, 59))
for(ik in 1:19)
{
  temp = get(ind_cob_mase_fix_GEO[ik])
  for(ij in 1:59)
  {
    if(ik == 14 && ij == 23)
    {
      ind_fix_geo_all_cob_mase[,ik,ij] = rep(0, 10)
    } else {
      ind_fix_geo_all_cob_mase[,ik,ij] = temp[[ij]]
    }
  }
}

ind_fix_geo_all_cob_rmse = array(0, dim = c(10, 19, 59))
for(ik in 1:19)
{
  temp = get(ind_cob_rmse_fix_GEO[ik])
  for(ij in 1:59)
  {
    if(ik == 14 && ij == 23)
    {
      ind_fix_geo_all_cob_rmse[,ik,ij] = rep(0, 10)
    } else {
      ind_fix_geo_all_cob_rmse[,ik,ij] = temp[[ij]]
    }
  }
}


# C_group independent forecasts

ind_c_group_fix_GEO <- function(c_group_index = 1,  fmethod = c("classical", "M"), year_horizon, geo_ind)
{
  require(ftsa)
  require(demography)
  fmethod = match.arg(fmethod)
  
  n_year = 31 - (year_horizon + 1)
  
  if(c_group_index == 6 && geo_ind == 23)
  {
    res_c_group = train_residual = c_group_mae = c_group_mase = c_group_rmse = 0  # no observation at all for cob 14 area 23
  } else {
    c_group_comb = get(fts_region_fix_geo_smooth[[c_group_index]][geo_ind])$rate$female[3:33,]
    colnames(c_group_comb) = year_interp
    rownames(c_group_comb) = age_interp
    
    res_c_group = array(NA, dim = c(year_horizon, 31, year_horizon))
    train_residual = list()
    for(ij in 1:year_horizon)
    {
      ind_dat = fts(x = age_interp, y = c_group_comb[,1:(n_year+ij)], xname = "age", yname = "fertility")
      ftsm_order = head(which(round(cumsum(ftsm(ind_dat, method = fmethod, lambda = 3, order = ncol(ind_dat$y))$varprop),3)>=0.95),1)
      ftsm_c_group = ftsm(ind_dat, order = ftsm_order, method = fmethod, lambda = 3)
      
      fun_forc  = forecast(ftsm_c_group, h = year_horizon)
      
      res_c_group[,,ij] = t(fun_forc$mean$y)
      
      train_residual[[ij]] = fun_forc$error$y
    }
    
    # Errors
    
    c_group_mae = c_group_mase = c_group_rmse = matrix(NA, year_horizon, 1)
    age_ind = match(age_all, age_interp) 
    
    for(ik in 1:year_horizon)
    {
      c_group_mae[ik,1] = ftsa:::mae(res_c_group[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_region_fix_geo[[c_group_index]][geo_ind]), years = (2001+ik):2011)$rate$female)
      c_group_mase[ik,1] = ftsa:::mase(res_c_group[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_region_fix_geo[[c_group_index]][geo_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_region_fix_geo[[c_group_index]][geo_ind]), years = 1986:(2001+ik-1))$rate$female)
      c_group_rmse[ik,1] = ftsa:::rmse(res_c_group[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_region_fix_geo[[c_group_index]][geo_ind]), years = (2001+ik):2011)$rate$female)
    }
  }

  return(list(res_c_group = res_c_group, train_residual = train_residual, c_group_mae = c_group_mae, c_group_mase = c_group_mase, c_group_rmse = c_group_rmse))    
}

test = ind_c_group_fix_GEO(c_group_index = 6, fmethod = "classical", year_horizon = 10, geo_ind = 24)

for(c_group_ind in 1:10)
{
  temp_res_forc_list = temp_train_residual_list = temp_mae_list = temp_mase_list = temp_rmse_list = list()
  for(geo in 1:59)
  {
    temp = ind_c_group_fix_GEO(c_group_index = c_group_ind, fmethod = "classical", year_horizon = 10, geo_ind = geo)
    temp_res_forc_list[[geo]] = temp$res_c_group
    temp_train_residual_list[[geo]] = temp$train_residual
    temp_mae_list[[geo]] = temp$c_group_mae
    temp_mase_list[[geo]] = temp$c_group_mase
    temp_rmse_list[[geo]] = temp$c_group_rmse
  }
  
  assign(ind_c_group_forecast_fix_GEO[c_group_ind], temp_res_forc_list)
  assign(ind_c_group_train_residual_fix_GEO[c_group_ind], temp_train_residual_list)
  assign(ind_c_group_mae_fix_GEO[c_group_ind], temp_mae_list)
  assign(ind_c_group_mase_fix_GEO[c_group_ind], temp_mase_list)
  assign(ind_c_group_rmse_fix_GEO[c_group_ind], temp_rmse_list)
  
  rm(temp, temp_res_forc_list, temp_train_residual_list, temp_mae_list, temp_mase_list, temp_rmse_list)
  
  print(c_group_ind)
}

# summary of results

ind_fix_geo_all_c_group_mae = array(0, dim = c(10, 10, 59))
for(ik in 1:10)
{
  temp = get(ind_c_group_mae_fix_GEO[ik])
  for(ij in 1:59)
  {
    if(ik == 6 && ij == 23)
    {
      ind_fix_geo_all_c_group_mae[,ik,ij] = rep(0, 10)
    } else {
      ind_fix_geo_all_c_group_mae[,ik,ij] = temp[[ij]]
    }
  }
}

ind_fix_geo_all_c_group_mase = array(0, dim = c(10, 10, 59))
for(ik in 1:10)
{
  temp = get(ind_c_group_mase_fix_GEO[ik])
  for(ij in 1:59)
  {
    if(ik == 6 && ij == 23)
    {
      ind_fix_geo_all_c_group_mase[,ik,ij] = rep(0, 10)
    } else {
      ind_fix_geo_all_c_group_mase[,ik,ij] = temp[[ij]]
    }
  }
}


ind_fix_geo_all_c_group_rmse = array(0, dim = c(10, 10, 59))
for(ik in 1:10)
{
  temp = get(ind_c_group_rmse_fix_GEO[ik])
  for(ij in 1:59)
  {
    if(ik == 6 && ij == 23)
    {
      ind_fix_geo_all_c_group_rmse[,ik,ij] = rep(0, 10)
    } else {
      ind_fix_geo_all_c_group_rmse[,ik,ij] = temp[[ij]]
    }
  }
}

# All birthplace independent forecasts

ind_c_global_fix_GEO <- function(fmethod = c("classical", "M"), year_horizon, geo_ind)
{
  require(ftsa)
  require(demography)
  fmethod = match.arg(fmethod)
  
  n_year = 31 - (year_horizon + 1)
  
  c_global_comb = get(fts_national_fix_geo_smooth[geo_ind])$rate$female[3:33,]
  colnames(c_global_comb) = year_interp
  rownames(c_global_comb) = age_interp
  
  res_c_global = array(NA, dim = c(year_horizon, 31, year_horizon))
  train_residual = list()
  for(ij in 1:year_horizon)
  {
    ind_dat = fts(x = age_interp, y = c_global_comb[,1:(n_year+ij)], xname = "age", yname = "fertility")
    ftsm_order = head(which(round(cumsum(ftsm(ind_dat, method = fmethod, lambda = 3, order = ncol(ind_dat$y))$varprop),3)>=0.95),1)
    ftsm_c_global = ftsm(ind_dat, order = ftsm_order, method = fmethod, lambda = 3)
    
    fun_forc  = forecast(ftsm_c_global, h = year_horizon)
    
    res_c_global[,,ij] = t(fun_forc$mean$y)
    
    train_residual[[ij]] = fun_forc$error$y
  }
  
  # Errors
  
  c_global_mae = c_global_mase = c_global_rmse = matrix(NA, year_horizon, 1)
  age_ind = match(age_all, age_interp) 
  
  for(ik in 1:year_horizon)
  {
    c_global_mae[ik,1] = ftsa:::mae(res_c_global[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_national_fix_geo[geo_ind]), years = (2001+ik):2011)$rate$female)
    c_global_mase[ik,1] = ftsa:::mase(res_c_global[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_national_fix_geo[geo_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_national_fix_geo[geo_ind]), years = 1986:(2001+ik-1))$rate$female)
    c_global_rmse[ik,1] = ftsa:::rmse(res_c_global[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_national_fix_geo[geo_ind]), years = (2001+ik):2011)$rate$female)
  }
  
  return(list(res_c_global = res_c_global, train_residual = train_residual, c_global_mae = c_global_mae, c_global_mase = c_global_mase, c_global_rmse = c_global_rmse))    
}

ind_global_train_residual_list = ind_global_forcast_list = ind_global_mae_list = ind_global_mase_list = ind_global_rmse_list = list()
for(geo in 1:59)
{
  temp = ind_c_global_fix_GEO(fmethod = "classical", year_horizon = 10, geo_ind = geo)
  ind_global_forcast_list[[geo]] = temp$res_c_global
  ind_global_train_residual_list[[geo]] = temp$train_residual
  ind_global_mae_list[[geo]] = temp$c_global_mae
  ind_global_mase_list[[geo]] = temp$c_global_mase
  ind_global_rmse_list[[geo]] = temp$c_global_rmse
  rm(temp)
}

assign(ind_global_forecast_fix_GEO[1], ind_global_forcast_list)
assign(ind_global_train_residual_fix_GEO[1], ind_global_train_residual_list)
assign(ind_global_mae_fix_GEO[1], ind_global_mae_list)
assign(ind_global_mase_fix_GEO[1], ind_global_mase_list)
assign(ind_global_rmse_fix_GEO[1], ind_global_rmse_list)

# summary of results

ind_fix_geo_all_global_mae = matrix(0, nrow = 10, ncol = 59)
for(ij in 1:59)
{
  ind_fix_geo_all_global_mae[,ij] = get(ind_global_mae_fix_GEO[1])[[ij]]
}

ind_fix_geo_all_global_mase = matrix(0, nrow = 10, ncol = 59)
for(ij in 1:59)
{
  ind_fix_geo_all_global_mase[,ij] = get(ind_global_mase_fix_GEO[1])[[ij]]
}

ind_fix_geo_all_global_rmse = matrix(0, nrow = 10, ncol = 59)
for(ij in 1:59)
{
  ind_fix_geo_all_global_rmse[,ij] = get(ind_global_rmse_fix_GEO[1])[[ij]]
}

# combine all point forecasts into an array

ind_fix_geo_all_level_mae = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  ind_fix_geo_all_level_mae[,,i] = cbind(ind_fix_geo_all_global_mae[,i], rowMeans(ind_fix_geo_all_c_group_mae[,,i]), rowMeans(ind_fix_geo_all_cob_mae[,,i]))
}
ind_fix_geo_mae_averaged = apply(ind_fix_geo_all_level_mae, c(1,2), mean)

###
ind_fix_geo_all_c_group_mase_copy = ind_fix_geo_all_c_group_mase
for(i in 1:59)
{
  ind_fix_geo_all_c_group_mase[,,i][,which(colMeans(ind_fix_geo_all_c_group_mase[,,i]) > 100)] = 0
}

ind_fix_geo_all_cob_mase_copy = ind_fix_geo_all_cob_mase
for(i in 1:59)
{
  ind_fix_geo_all_cob_mase[,,i][,which(colMeans(ind_fix_geo_all_cob_mase[,,i]) > 300)] = 0
}
###

ind_fix_geo_all_level_mase = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  ind_fix_geo_all_level_mase[,,i] = cbind(ind_fix_geo_all_global_mase[,i], rowMeans(ind_fix_geo_all_c_group_mase[,,i]), rowMeans(ind_fix_geo_all_cob_mase[,,i]))
}
ind_fix_geo_mase_averaged = apply(ind_fix_geo_all_level_mase, c(1,2), mean)


ind_fix_geo_all_level_rmse = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  ind_fix_geo_all_level_rmse[,,i] = cbind(ind_fix_geo_all_global_rmse[,i], rowMeans(ind_fix_geo_all_c_group_rmse[,,i]), rowMeans(ind_fix_geo_all_cob_rmse[,,i]))
}
ind_fix_geo_rmse_averaged = apply(ind_fix_geo_all_level_rmse, c(1,2), mean)

#########################################
# Reconciled independent point forecasts
#########################################

##########
# Fix COB
##########

########################################
# Point forecast errors; fmethod = "BU" 
########################################

ind_national_mae_BU = ind_national_mase_BU = ind_national_rmse_BU = array(0, dim = c(10, 1, 19))
ind_region_mae_BU = ind_region_mase_BU = ind_region_rmse_BU = array(0, dim = c(10, 11, 19))
ind_area_mae_BU  = ind_area_mase_BU = ind_area_rmse_BU = array(0, dim = c(10, 47, 19))

for(cob in 1:19)
{
  for(ikw in 1:10)
  {
    dum = BU_optim_err_fix_cob(ik = ikw, hier_method = "BU", 
                               national_forc = ind_national_forecast_fix_COB, 
                               region_forc = ind_region_forecast_fix_COB, 
                               area_forc = ind_area_forecast_fix_COB, 
                               cob_ind = cob)
    
    # MAE
    ind_national_mae_BU[ikw,,cob] = dum$mae_national
    ind_region_mae_BU[ikw,,cob] = dum$mae_region
    ind_area_mae_BU[ikw,,cob] = dum$mae_area
    
    # MAE
    ind_national_mase_BU[ikw,,cob] = dum$mase_national
    ind_region_mase_BU[ikw,,cob] = dum$mase_region
    ind_area_mase_BU[ikw,,cob] = dum$mase_area
    
    # RMSE
    ind_national_rmse_BU[ikw,,cob] = dum$rmse_national
    ind_region_rmse_BU[ikw,,cob] = dum$rmse_region
    ind_area_rmse_BU[ikw,,cob] = dum$rmse_area
    
  }
}

###
ind_area_mase_BU_copy = ind_area_mase_BU
for(i in 1:19)
{
  ind_area_mase_BU[,,i][,which(colMeans(ind_area_mase_BU[,,i]) > 300)] = 0
}
###



# Summary of results

ind_fix_cob_BU_mae = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  ind_fix_cob_BU_mae[,,i] = cbind(ind_national_mae_BU[,,i], rowMeans(ind_region_mae_BU[,,i]), rowMeans(ind_area_mae_BU[,,i]))
}
ind_fix_cob_BU_mae_averaged = apply(ind_fix_cob_BU_mae, c(1,2), mean)

ind_fix_cob_BU_mase = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  ind_fix_cob_BU_mase[,,i] = cbind(ind_national_mase_BU[,,i], rowMeans(ind_region_mase_BU[,,i]), rowMeans(ind_area_mase_BU[,,i]))
}
ind_fix_cob_BU_mase_averaged = apply(ind_fix_cob_BU_mase, c(1,2), mean)


ind_fix_cob_BU_rmse = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  ind_fix_cob_BU_rmse[,,i] = cbind(ind_national_rmse_BU[,,i], rowMeans(ind_region_rmse_BU[,,i]), rowMeans(ind_area_rmse_BU[,,i]))
}
ind_fix_cob_BU_rmse_averaged = apply(ind_fix_cob_BU_rmse, c(1,2), mean)

##############################################
# Point forecast errors; fmethod = "comb_OLS" 
##############################################

ind_national_mae_OP = ind_national_mase_OP = ind_national_rmse_OP = array(0, dim = c(10, 1, 19))
ind_region_mae_OP = ind_region_mase_OP = ind_region_rmse_OP = array(0, dim = c(10, 11, 19))
ind_area_mae_OP = ind_area_mase_OP = ind_area_rmse_OP = array(0, dim = c(10, 47, 19))

for(cob in 1:19)
{
  for(ikw in 1:10)
  {
    dum = BU_optim_err_fix_cob(ik = ikw, hier_method = "comb_OLS", 
                               national_forc = ind_national_forecast_fix_COB, 
                               region_forc = ind_region_forecast_fix_COB, 
                               area_forc = ind_area_forecast_fix_COB, 
                               cob_ind = cob)
    
    # MAE
    ind_national_mae_OP[ikw,,cob] = dum$mae_national
    ind_region_mae_OP[ikw,,cob] = dum$mae_region
    ind_area_mae_OP[ikw,,cob] = dum$mae_area
    
    # MASE
    ind_national_mase_OP[ikw,,cob] = dum$mase_national
    ind_region_mase_OP[ikw,,cob] = dum$mase_region
    ind_area_mase_OP[ikw,,cob] = dum$mase_area
    
    # RMSE
    ind_national_rmse_OP[ikw,,cob] = dum$rmse_national
    ind_region_rmse_OP[ikw,,cob] = dum$rmse_region
    ind_area_rmse_OP[ikw,,cob] = dum$rmse_area
    
  }
}

###
ind_area_mase_OP_copy = ind_area_mase_OP
for(i in 1:19)
{
  ind_area_mase_OP[,,i][,which(colMeans(ind_area_mase_OP[,,i]) > 300)] = 0
}
###


# Summary of results

ind_fix_cob_OP_mae = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  ind_fix_cob_OP_mae[,,i] = cbind(ind_national_mae_OP[,,i], rowMeans(ind_region_mae_OP[,,i]), rowMeans(ind_area_mae_OP[,,i]))
}
ind_fix_cob_OP_mae_averaged = apply(ind_fix_cob_OP_mae, c(1,2), mean)


ind_fix_cob_OP_mase = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  ind_fix_cob_OP_mase[,,i] = cbind(ind_national_mase_OP[,,i], rowMeans(ind_region_mase_OP[,,i]), rowMeans(ind_area_mase_OP[,,i]))
}
ind_fix_cob_OP_mase_averaged = apply(ind_fix_cob_OP_mase, c(1,2), mean)


ind_fix_cob_OP_rmse = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  ind_fix_cob_OP_rmse[,,i] = cbind(ind_national_rmse_OP[,,i], rowMeans(ind_region_rmse_OP[,,i]), rowMeans(ind_area_rmse_OP[,,i]))
}
ind_fix_cob_OP_rmse_averaged = apply(ind_fix_cob_OP_rmse, c(1,2), mean)

##########################################
# Point forecast errors; fmethod = "mint" 
##########################################

ind_national_mae_MINT = ind_national_mase_MINT = ind_national_rmse_MINT = array(0, dim = c(10, 1, 19))
ind_region_mae_MINT = ind_region_mase_MINT = ind_region_rmse_MINT = array(0, dim = c(10, 11, 19))
ind_area_mae_MINT = ind_area_mase_MINT = ind_area_rmse_MINT = array(0, dim = c(10, 47, 19))

for(cob in 1:19)
{
  for(ikw in 1:10)
  {
    dum = ind_BU_optim_err_fix_cob(ik = ikw, hier_method = "mint", 
                                   national_forc = ind_national_forecast_fix_COB, 
                                   region_forc = ind_region_forecast_fix_COB, 
                                   area_forc = ind_area_forecast_fix_COB, 
                                   cob_ind = cob)
    
    # MAE
    ind_national_mae_MINT[ikw,,cob] = dum$mae_national
    ind_region_mae_MINT[ikw,,cob] = dum$mae_region
    ind_area_mae_MINT[ikw,,cob] = dum$mae_area
    
    # MASE
    ind_national_mase_MINT[ikw,,cob] = dum$mase_national
    ind_region_mase_MINT[ikw,,cob] = dum$mase_region
    ind_area_mase_MINT[ikw,,cob] = dum$mase_area
    
    # RMSE
    ind_national_rmse_MINT[ikw,,cob] = dum$rmse_national
    ind_region_rmse_MINT[ikw,,cob] = dum$rmse_region
    ind_area_rmse_MINT[ikw,,cob] = dum$rmse_area
    
  }
}

###
ind_area_mase_MINT_copy = ind_area_mase_MINT
for(i in 1:19)
{
  ind_area_mase_MINT[,,i][,which(colMeans(ind_area_mase_MINT[,,i]) > 300)] = 0
}
###

# Summary of results

ind_fix_cob_MINT_mae = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  ind_fix_cob_MINT_mae[,,i] = cbind(ind_national_mae_MINT[,,i], rowMeans(ind_region_mae_MINT[,,i]), rowMeans(ind_area_mae_MINT[,,i]))
}
ind_fix_cob_MINT_mae_averaged = apply(ind_fix_cob_MINT_mae, c(1,2), mean)


ind_fix_cob_MINT_mase = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  ind_fix_cob_MINT_mase[,,i] = cbind(ind_national_mase_MINT[,,i], rowMeans(ind_region_mase_MINT[,,i]), rowMeans(ind_area_mase_MINT[,,i]))
}
ind_fix_cob_MINT_mase_averaged = apply(ind_fix_cob_MINT_mase, c(1,2), mean)


ind_fix_cob_MINT_rmse = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  ind_fix_cob_MINT_rmse[,,i] = cbind(ind_national_rmse_MINT[,,i], rowMeans(ind_region_rmse_MINT[,,i]), rowMeans(ind_area_rmse_MINT[,,i]))
}
ind_fix_cob_MINT_rmse_averaged = apply(ind_fix_cob_MINT_rmse, c(1,2), mean)


##########
# Fix GEO
##########

ind_recon_hier_fix_geo <- function(global_forc, c_group_forc, cob_forc, age, kj, hier_method = c("BU", "comb_OLS", "comb_GLS", "mint"), geo_ind)
{
  hier_method = match.arg(hier_method)
  
  age_ind = match(age, age_interp)
  
  hier_rate = matrix(NA,30,(11-kj)) # in total 59 series
  
  # Level_1 global
  
  hier_rate[1,] = get(global_forc[1])[[geo_ind]][kj,age_ind,1:(11-kj)] #  arguments: horizon, age, 11- horzon; length of list == geo
  
  
  # Level 2 disaggregate by C group
  for(iw in 1:10)
  {
    if(iw == 6 && geo_ind == 23)
    {
      hier_rate[iw+1,] = rep(0, 11-kj)
    } else {
      hier_rate[iw+1,] = get(c_group_forc[iw])[[geo_ind]][kj,age_ind,1:(11-kj)]  
    }
  }
  
  # Level 3 disaggregate by cob
  
  for(iw in 1:19)
  {
    if(iw == 14 && geo_ind == 23)
    {
      hier_rate[iw+11,] = rep(0, 11-kj)
    } else {
      hier_rate[iw+11,] = get(cob_forc[iw])[[geo_ind]][kj,age_ind,1:(11-kj)]
    }
  }  
  
  # forecast reconciliation via bottom-up, optimal combination or MinT methods
  hier_fore = matrix(NA, 30, (11-kj))
  summing_mat = Smat_fun_fix_geo(kj = kj, age = age, no_cob = 19, geo = geo_ind)
  for(ik in 1:(11-kj))
  {
    hier = summing_mat[,,ik]
    if(hier_method == "BU")
    {
      hier_fore[,ik] = (hier %*% hier_rate[12:30,ik])
    }
    if(hier_method == "comb_OLS")
    {
      hier_fore[,ik] = hier %*% ginv(t(hier) %*% hier) %*% t(hier) %*% hier_rate[,ik]
    }
    if(hier_method == "comb_GLS")
    {
      hier_fore[,ik] = hier %*% ginv(t(hier) %*% ginv(diag(rowSums(hier^2))) %*% hier) %*% t(hier) %*% ginv(diag(rowSums(hier^2))) %*% hier_rate[,ik]
    }
    if(hier_method == "mint")
    {
      ginv = MASS:::ginv
      wh = ind_wh_fun_fix_geo(kj = kj, age = age, geo_ind = geo_ind)
      hier_fore[,ik] = hier %*% ginv(t(hier) %*% ginv(wh) %*% hier) %*% t(hier) %*% ginv(wh) %*% hier_rate[,ik]
    }
  }
  return(hier_fore)
}


# point forecast errors function (shared by both univariate and multivariate functional time series forecasting methods)

ind_BU_optim_err_fix_geo <- function(ik, hier_method, global_forc, c_group_forc, cob_forc, geo_ind)
{
  age_all = c(15, 20, 25, 30, 35, 40, 45)
  
  me = ftsa:::me; mae = ftsa:::mae; rmse = ftsa:::rmse; mase = ftsa:::mase
  BU_optim_hier_comb = array(NA, dim = c(7, 30,(11-ik)))
  
  for(i in 1:7)
  {
    BU_optim_hier_comb[i,,] = ind_recon_hier_fix_geo(global_forc, c_group_forc, cob_forc,      
                                                     age = age_all[i], kj = ik, hier_method = hier_method,
                                                     geo_ind = geo_ind)
  }
  
  #####################################
  # Errors, including ME, MAE and RMSE
  #####################################
  
  # Level 1 (Global)
  mae_global = mae(BU_optim_hier_comb[,1,], extract.years(get(fts_national_fix_geo[geo_ind]), years = (2001+ik):2011)$rate$female)
  mase_global = mase(BU_optim_hier_comb[,1,], extract.years(get(fts_national_fix_geo[geo_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_national_fix_geo[geo_ind]), years = 1986:(2001+ik-1))$rate$female)
  rmse_global = rmse(BU_optim_hier_comb[,1,], extract.years(get(fts_national_fix_geo[geo_ind]), years = (2001+ik):2011)$rate$female)
  
  # Level 2 (C_group)
  mae_c_group = mase_c_group = rmse_c_group = rep(0, 10)
  for(iw in 1:10)
  {
    if(iw ==6 && geo_ind == 23)
    {
      mae_c_group[iw] = mase_c_group[iw] = rmse_c_group[iw] = 0
    } else {
      mae_c_group[iw]  = mae(BU_optim_hier_comb[,(iw+1),],  extract.years(get(fts_region_fix_geo[[iw]][geo_ind]), years = (2001+ik):2011)$rate$female)
      mase_c_group[iw]  = mase(BU_optim_hier_comb[,(iw+1),],  extract.years(get(fts_region_fix_geo[[iw]][geo_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_region_fix_geo[[iw]][geo_ind]), years = 1986:(2001+ik-1))$rate$female)
      rmse_c_group[iw] = rmse(BU_optim_hier_comb[,(iw+1),], extract.years(get(fts_region_fix_geo[[iw]][geo_ind]), years = (2001+ik):2011)$rate$female)
    }
  }
  
  # Level 3 (COB)
  mae_cob  = mase_cob = rmse_cob = rep(0, 10)
  
  for(iwk in 1:19)
  {
    if(iwk == 14 && geo_ind == 23)
    {
      mae_cob[iwk] = mase_cob[iwk] = rmse_cob[iwk] = 0
    } else {
      mae_cob[iwk]  = mae(BU_optim_hier_comb[,(iwk+11),],   extract.years(get(fts_area_fix_geo[[iwk]][geo_ind]), years = (2001+ik):2011)$rate$female)
      mase_cob[iwk]  = mase(BU_optim_hier_comb[,(iwk+11),],   extract.years(get(fts_area_fix_geo[[iwk]][geo_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_area_fix_geo[[iwk]][geo_ind]), years = 1986:(2001+ik-1))$rate$female)
      rmse_cob[iwk] = rmse(BU_optim_hier_comb[,(iwk+11),],  extract.years(get(fts_area_fix_geo[[iwk]][geo_ind]), years = (2001+ik):2011)$rate$female)
    }
  }
  
  
  return(list(mae_global = mae_global, 
              mae_c_group = mae_c_group, 
              mae_cob  = mae_cob, 
              
              mase_global = mase_global, 
              mase_c_group = mase_c_group, 
              mase_cob  = mase_cob, 
              
              rmse_global = rmse_global, 
              rmse_c_group = rmse_c_group, 
              rmse_cob  = rmse_cob))
}


########################################
# Point forecast errors; fmethod = "BU" 
########################################

ind_global_mae_BU  = ind_global_mase_BU = ind_global_rmse_BU = array(0, dim = c(10, 1, 59))
ind_c_group_mae_BU = ind_c_group_mase_BU = ind_c_group_rmse_BU = array(0, dim = c(10, 10, 59))
ind_cob_mae_BU = ind_cob_mase_BU = ind_cob_rmse_BU = array(0, dim = c(10, 19, 59))

for(geo in 1:59)
{
  for(ikw in 1:10)
  {
    dum = BU_optim_err_fix_geo(ik = ikw, hier_method = "BU", 
                               global_forc = ind_global_forecast_fix_GEO, 
                               c_group_forc = ind_c_group_forecast_fix_GEO, 
                               cob_forc = ind_cob_forecast_fix_GEO, 
                               geo_ind = geo)
    
    # MAE
    ind_global_mae_BU[ikw,,geo] = dum$mae_global
    ind_c_group_mae_BU[ikw,,geo] = dum$mae_c_group
    ind_cob_mae_BU[ikw,,geo] = dum$mae_cob
    
    # MASE
    ind_global_mase_BU[ikw,,geo] = dum$mase_global
    ind_c_group_mase_BU[ikw,,geo] = dum$mase_c_group
    ind_cob_mase_BU[ikw,,geo] = dum$mase_cob
    
    # RMSE
    ind_global_rmse_BU[ikw,,geo] = dum$rmse_global
    ind_c_group_rmse_BU[ikw,,geo] = dum$rmse_c_group
    ind_cob_rmse_BU[ikw,,geo] = dum$rmse_cob
    
  }
}

# Summary of results

ind_fix_geo_BU_mae = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  ind_fix_geo_BU_mae[,,i] = cbind(ind_global_mae_BU[,,i], rowMeans(ind_c_group_mae_BU[,,i]), rowMeans(ind_cob_mae_BU[,,i]))
}
ind_fix_geo_BU_mae_averaged = apply(ind_fix_geo_BU_mae, c(1,2), mean)

###
ind_c_group_mase_BU_copy = ind_c_group_mase_BU
for(i in 1:59)
{
  ind_c_group_mase_BU[,,i][,which(colMeans(ind_c_group_mase_BU[,,i]) > quantile(ind_c_group_mase_BU, 0.99) )] = 0
}

ind_cob_mase_BU_copy = ind_cob_mase_BU
for(i in 1:59)
{
  ind_cob_mase_BU[,,i][,which(colMeans(ind_cob_mase_BU[,,i]) > 300)] = 0
}
###

ind_fix_geo_BU_mase = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  ind_fix_geo_BU_mase[,,i] = cbind(ind_global_mase_BU[,,i], rowMeans(ind_c_group_mase_BU[,,i]), rowMeans(ind_cob_mase_BU[,,i]))
}
ind_fix_geo_BU_mase_averaged = apply(ind_fix_geo_BU_mase, c(1,2), mean)


ind_fix_geo_BU_rmse = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  ind_fix_geo_BU_rmse[,,i] = cbind(ind_global_rmse_BU[,,i], rowMeans(ind_c_group_rmse_BU[,,i]), rowMeans(ind_cob_rmse_BU[,,i]))
}
ind_fix_geo_BU_rmse_averaged = apply(ind_fix_geo_BU_rmse, c(1,2), mean)

##############################################
# Point forecast errors; fmethod = "comb_OLS" 
##############################################

ind_global_mae_OP = ind_global_mase_OP = ind_global_rmse_OP = array(0, dim = c(10, 1, 59))
ind_c_group_mae_OP = ind_c_group_mase_OP = ind_c_group_rmse_OP = array(0, dim = c(10, 10, 59))
ind_cob_mae_OP = ind_cob_mase_OP = ind_cob_rmse_OP = array(0, dim = c(10, 19, 59))

for(geo in 1:59)
{
  for(ikw in 1:10)
  {
    dum = BU_optim_err_fix_geo(ik = ikw, hier_method = "comb_OLS", 
                               global_forc = ind_global_forecast_fix_GEO, 
                               c_group_forc = ind_c_group_forecast_fix_GEO, 
                               cob_forc = ind_cob_forecast_fix_GEO, 
                               geo_ind = geo)
    
    # MAE
    ind_global_mae_OP[ikw,,geo] = dum$mae_global
    ind_c_group_mae_OP[ikw,,geo] = dum$mae_c_group
    ind_cob_mae_OP[ikw,,geo] = dum$mae_cob
    
    # MASE
    ind_global_mase_OP[ikw,,geo] = dum$mase_global
    ind_c_group_mase_OP[ikw,,geo] = dum$mase_c_group
    ind_cob_mase_OP[ikw,,geo] = dum$mase_cob
    
    # RMSE
    ind_global_rmse_OP[ikw,,geo] = dum$rmse_global
    ind_c_group_rmse_OP[ikw,,geo] = dum$rmse_c_group
    ind_cob_rmse_OP[ikw,,geo] = dum$rmse_cob
    
  }
}

ind_fix_geo_OP_mae = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  ind_fix_geo_OP_mae[,,i] = cbind(ind_global_mae_OP[,,i], rowMeans(ind_c_group_mae_OP[,,i]), rowMeans(ind_cob_mae_OP[,,i]))
}
ind_fix_geo_OP_mae_averaged = apply(ind_fix_geo_OP_mae, c(1,2), mean)

###
ind_c_group_mase_OP_copy = ind_c_group_mase_OP
for(i in 1:59)
{
  ind_c_group_mase_OP[,,i][,which(colMeans(ind_c_group_mase_OP[,,i]) > 400)] = 0
}

ind_cob_mase_OP_copy = ind_cob_mase_OP
for(i in 1:59)
{
  ind_cob_mase_OP[,,i][,which(colMeans(ind_cob_mase_OP[,,i]) > 200)] = 0
}
###

ind_fix_geo_OP_mase = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  ind_fix_geo_OP_mase[,,i] = cbind(ind_global_mase_OP[,,i], rowMeans(ind_c_group_mase_OP[,,i]), rowMeans(ind_cob_mase_OP[,,i]))
}
ind_fix_geo_OP_mase_averaged = apply(ind_fix_geo_OP_mase, c(1,2), mean)


ind_fix_geo_OP_rmse = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  ind_fix_geo_OP_rmse[,,i] = cbind(ind_global_rmse_OP[,,i], rowMeans(ind_c_group_rmse_OP[,,i]), rowMeans(ind_cob_rmse_OP[,,i]))
}
ind_fix_geo_OP_rmse_averaged = apply(ind_fix_geo_OP_rmse, c(1,2), mean)

##########################################
# Point forecast errors; fmethod = "mint" 
##########################################

ind_global_mae_MINT = ind_global_mase_MINT = ind_global_rmse_MINT = array(0, dim = c(10, 1, 59))
ind_c_group_mae_MINT = ind_c_group_mase_MINT = ind_c_group_rmse_MINT = array(0, dim = c(10, 10, 59))
ind_cob_mae_MINT = ind_cob_mase_MINT = ind_cob_rmse_MINT = array(0, dim = c(10, 19, 59))

for(geo in 1:59)
{
  for(ikw in 1:10)
  {
    dum = ind_BU_optim_err_fix_geo(ik = ikw, hier_method = "mint", 
                                   global_forc = ind_global_forecast_fix_GEO, 
                                   c_group_forc = ind_c_group_forecast_fix_GEO, 
                                   cob_forc = ind_cob_forecast_fix_GEO,
                                   geo_ind = geo)
    
    # MAE
    ind_global_mae_MINT[ikw,,geo] = dum$mae_global
    ind_c_group_mae_MINT[ikw,,geo] = dum$mae_c_group
    ind_cob_mae_MINT[ikw,,geo] = dum$mae_cob
    
    # MASE
    ind_global_mase_MINT[ikw,,geo] = dum$mase_global
    ind_c_group_mase_MINT[ikw,,geo] = dum$mase_c_group
    ind_cob_mase_MINT[ikw,,geo] = dum$mase_cob
    
    # RMSE
    ind_global_rmse_MINT[ikw,,geo] = dum$rmse_global
    ind_c_group_rmse_MINT[ikw,,geo] = dum$rmse_c_group
    ind_cob_rmse_MINT[ikw,,geo] = dum$rmse_cob
    
  }
}

ind_fix_geo_MINT_mae = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  ind_fix_geo_MINT_mae[,,i] = cbind(ind_global_mae_MINT[,,i], rowMeans(ind_c_group_mae_MINT[,,i]), rowMeans(ind_cob_mae_MINT[,,i]))
}
ind_fix_geo_MINT_mae_averaged = apply(ind_fix_geo_MINT_mae, c(1,2), mean)

###
ind_c_group_mase_MINT_copy = ind_c_group_mase_MINT
for(i in 1:59)
{
  ind_c_group_mase_MINT[,,i][,which(colMeans(ind_c_group_mase_MINT[,,i]) > 300)] = 0
}

ind_cob_mase_MINT_copy = ind_cob_mase_MINT
for(i in 1:59)
{
  ind_cob_mase_MINT[,,i][,which(colMeans(ind_cob_mase_MINT[,,i]) > 300)] = 0
}
###

ind_fix_geo_MINT_mase = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  ind_fix_geo_MINT_mase[,,i] = cbind(ind_global_mase_MINT[,,i], rowMeans(ind_c_group_mase_MINT[,,i]), rowMeans(ind_cob_mase_MINT[,,i]))
}
ind_fix_geo_MINT_mase_averaged = apply(ind_fix_geo_MINT_mase, c(1,2), mean)


ind_fix_geo_MINT_rmse = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  ind_fix_geo_MINT_rmse[,,i] = cbind(ind_global_rmse_MINT[,,i], rowMeans(ind_c_group_rmse_MINT[,,i]), rowMeans(ind_cob_rmse_MINT[,,i]))
}
ind_fix_geo_MINT_rmse_averaged = apply(ind_fix_geo_MINT_rmse, c(1,2), mean)


##########################################################
# constructing independent pointwise prediction intervals
##########################################################

# Area
ind_PI <- function(dat, dat_unsmooth, fh = 1, nboot = 1000, alpha = 0.8)
{
  # select ncomp_comb based on 95% of total variation; use all available data
  comb_object = t(scale(t(log(dat$rate[[1]][3:33,])), center = TRUE, scale = FALSE))
  comb_mean = rowMeans(log(dat$rate[[1]][3:33,]))
  rownames(comb_object) = dat$age[3:33]
  
  C_0 = long_run_covariance_estimation(comb_object, H = 3, C0 = 3)
  eigen_decomp = eigen(C_0)
  ncomp_comb = min(head(which(cumsum(eigen_decomp$values)/sum(eigen_decomp$values) >= 0.95),1), 2)
  
  dynamic_basis = as.matrix(eigen_decomp$vectors[,1:ncomp_comb])
  dynamic_scores = t(dynamic_basis) %*% comb_object
  
  # calculate in-sample forecast curves
  fore_curve_total = matrix(NA, nrow = 31, ncol = length(dat$year) - ncomp_comb - fh + 1)
  
  for(ij in 1:(length(dat$year) - ncomp_comb - fh + 1))
  {
    # compute mean and standard deviation functions
    dat_one_step = as.data.frame(log(dat$rate[[1]][3:33, 1:(ncomp_comb + ij - 1)]))
    
    scores_fit = scores_fore = list()
    fore_ftsm_dyn = matrix(NA, nrow = nrow(dat_one_step), ncol = fh)
    
    for(ik in 1:ncomp_comb)
    {
      scores_fit[[ik]] = auto.arima(dynamic_scores[ik,])
      scores_fore[[ik]] = forecast(scores_fit[[ik]], h = fh)$mean
    }
    
    for(ih in 1:fh)
    {
      fore_ftsm_dyn[,ih] = dynamic_basis %*% unlist(lapply(scores_fore,`[[`,ih))
    }
    
    fore_res = exp(fore_ftsm_dyn + comb_mean)
    
    # fill in gender specific in-sample forecast curves
    fore_curve_total[,ij] = t(fore_res[,fh])
  }
  
  # holdout data samples
  
  true_dat = extract.years(dat_unsmooth, dat_unsmooth$year[(ncomp_comb + fh):length(dat_unsmooth$year)])
  true_dat_smooth = extract.years(dat, dat$year[(ncomp_comb + fh):length(dat$year)])
  
  holdout_val = holdout_smooth = array(NA, dim = c(length(dat_unsmooth$age), length(dat_unsmooth$year) - ncomp_comb - fh + 1, n_pop))
  age_ind = match(dat_unsmooth$age, dat$age[3:33])
  
  holdout_val = (true_dat$rate[[1]])
  holdout_smooth = (true_dat_smooth$rate[[1]][age_ind,])
  holdout_val_index = which(!is.finite(holdout_val))
  
  if(length(holdout_val_index) > 0)
  {
    holdout_val = replace(holdout_val, holdout_val_index, holdout_smooth[holdout_val_index])
  }
  
  # Warnings if infinite values detected
  if(length(which(!is.finite(holdout_val))) > 0)
  {
    warning("some holdout data have Inf")
  }
  
  # region level total series errors
  err_total = holdout_val - fore_curve_total[age_ind,,drop = FALSE]
  
  
  # bootstrap error function
  err_boot_total = matrix(NA, nrow = nrow(err_total), ncol = nboot)
  for(ij in 1:nboot)
  {
    err_boot_total[,ij] = err_total[, sample(1:ncol(err_total), 1, replace = TRUE)]
  }
  # constructing PI
  
  fore_ftsm = forecast(ftsm(fts(1:nrow(comb_object), comb_object), order = ncomp_comb), h = fh)
  fore_res  = exp(fore_ftsm$mean$y + comb_mean)
  fore_mfts_total = cbind(fore_res[,fh])
  
  boot_PI_total = matrix(NA, nrow = length(dat_unsmooth$age), ncol = nboot)
  boot_PI_total = err_boot_total + matrix(rep(fore_mfts_total[age_ind,], nboot), nrow = length(dat_unsmooth$age), ncol = nboot)
  
  
  return(list(boot_sample = boot_PI_total))
}


##########
# Fix COB
##########

ind_area_PI_fix_COB <- function(data_series = get(fts_area_fix_cob_smooth[[area_index]][cob_ind]),  data_series_unsmooth = get(fts_area_fix_cob[[area_index]][cob_ind]), fh, nboot = 1000)
{
  require(demography)
  
  PI_boot_area = array(0, dim = c(7, nboot, (11-fh)))
  
  for(ij in 1:(11-fh))
  {
    dum_pointwise = ind_PI(dat = extract.years(data_series, 1981:(2000+ij)), dat_unsmooth = extract.years(data_series_unsmooth, 1981:(2000+ij)), fh = fh, nboot = nboot)
    PI_boot_area[,,ij]  = dum_pointwise$boot_sample
  }
  
  return(list(PI_boot = PI_boot_area))
}

test = ind_area_PI_fix_COB(data_series = get(fts_area_fix_cob_smooth[[1]][1]),  data_series_unsmooth = get(fts_area_fix_cob[[1]][1]), fh = 1, nboot = 20)

library(doParallel)

ind_area_PI_all_fix_COB = vector("character", length = 47)
for(i in 1:47)
{
  ind_area_PI_all_fix_COB[[i]] = paste("ind_A", i, "_PI_fix_COB", sep = "")
  temp_list = list()

  cl <- makeCluster(20) 
  registerDoParallel(cl)
  
  for(cob in 1:19)
  {
    if(i == 23 && cob == 14)
    {
      temp_list[[cob]] = NA
    } else {
      clusterExport(cl,list(fts_area_fix_cob_smooth[[i]][cob], fts_area_fix_cob[[i]][cob]))
      
      temp_list[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% ind_area_PI_fix_COB(data_series = get(fts_area_fix_cob_smooth[[i]][cob]), data_series_unsmooth = get(fts_area_fix_cob[[i]][cob]) , fh = ik, nboot = 1000)
    }
  }
  
  stopCluster(cl)
  rm(cl)
  
  assign(ind_area_PI_all_fix_COB[[i]], temp_list)
  rm(temp_list)
  
  print(i)
}

# Region

ind_region_PI_fix_COB <- function(data_series = get(fts_region_fix_cob_smooth[[region_index]][cob_ind]),  data_series_unsmooth = get(fts_region_fix_cob[[region_index]][cob_ind]), fh, nboot = 1000)
{
  require(demography)
  
  PI_boot_region = array(0, dim = c(7, nboot, (11-fh)))
  
  for(ij in 1:(11-fh))
  {
    dum_pointwise = ind_PI(dat = extract.years(data_series, 1981:(2000+ij)), dat_unsmooth = extract.years(data_series_unsmooth, 1981:(2000+ij)), fh = fh, nboot = nboot)
    PI_boot_region[,,ij]  = dum_pointwise$boot_sample
  }
  
  return(list(PI_boot = PI_boot_region))
}

test = ind_region_PI_fix_COB(data_series = get(fts_region_fix_cob_smooth[[11]][14]),  data_series_unsmooth = get(fts_region_fix_cob[[11]][14]), fh = 1, nboot = 20)


ind_region_PI_all_fix_COB = vector("character", length = 11)
for(i in 1:11)
{
  ind_region_PI_all_fix_COB[[i]] = paste("ind_R", i, "_PI_fix_COB", sep = "")
  temp_list = list()
  
  cl <- makeCluster(20) 
  registerDoParallel(cl)
  
  for(cob in 1:19)
  {
    clusterExport(cl,list(fts_region_fix_cob_smooth[[i]][cob], fts_region_fix_cob[[i]][cob]))
    
    temp_list[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% ind_region_PI_fix_COB(data_series = get(fts_region_fix_cob_smooth[[i]][cob]), data_series_unsmooth = get(fts_region_fix_cob[[i]][cob]) , fh = ik, nboot = 1000)
  }
  
  stopCluster(cl)
  rm(cl)
  
  assign(ind_region_PI_all_fix_COB[[i]], temp_list)
  rm(temp_list)
  
  print(i)
}

# National

ind_national_PI_fix_COB <- function(data_series = get(fts_national_fix_cob_smooth[cob_ind]),  data_series_unsmooth = get(fts_national_fix_cob[cob_ind]), fh, nboot = 1000)
{
  require(demography)
  require(ftsa)
  
  PI_boot_national = array(0, dim = c(7, nboot, (11-fh)))
  
  for(ij in 1:(11-fh))
  {
    dum_pointwise = ind_PI(dat = extract.years(data_series, 1981:(2000+ij)), dat_unsmooth = extract.years(data_series_unsmooth, 1981:(2000+ij)), fh = fh, nboot = nboot)
    PI_boot_national[,,ij]  = dum_pointwise$boot_sample
  }
  
  return(list(PI_boot = PI_boot_national))
}

test = ind_national_PI_fix_COB(data_series = get(fts_national_fix_cob_smooth[14]),  data_series_unsmooth = get(fts_national_fix_cob[14]), fh = 1, nboot = 20)

cl <- makeCluster(20) 
registerDoParallel(cl)

ind_national_PI_all_fix_COB = list()

for(cob in 1:19)
{
  clusterExport(cl,list(fts_national_fix_cob_smooth[cob], fts_national_fix_cob[cob]))
  
  ind_national_PI_all_fix_COB[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% ind_national_PI_fix_COB(data_series = get(fts_national_fix_cob_smooth[cob]), data_series_unsmooth = get(fts_national_fix_cob[cob]) , fh = ik, nboot = 1000)
}

stopCluster(cl)
rm(cl)


##########
# Fix GEO
##########

# COB
ind_cob_PI_fix_GEO <- function(data_series = get(fts_area_fix_geo_smooth[[cob_index]][geo_ind]),  data_series_unsmooth = get(fts_area_fix_geo[[cob_index]][geo_ind]), fh, nboot = 1000)
{
  require(demography)
  
  PI_boot_cob = array(0, dim = c(7, nboot, (11-fh)))
  
  for(ij in 1:(11-fh))
  {
    dum_pointwise = ind_PI(dat = extract.years(data_series, 1981:(2000+ij)), dat_unsmooth = extract.years(data_series_unsmooth, 1981:(2000+ij)), fh = fh, nboot = nboot)
    PI_boot_cob[,,ij]  = dum_pointwise$boot_sample
  }
  
  return(list(PI_boot = PI_boot_cob))
}


test = ind_cob_PI_fix_GEO(data_series = get(fts_area_fix_geo_smooth[[1]][1]),  data_series_unsmooth = get(fts_area_fix_geo[[1]][1]), fh = 1, nboot = 20)

library(doParallel)

ind_cob_PI_all_fix_GEO = vector("character", length = 19)
for(i in 1:19)
{
  ind_cob_PI_all_fix_GEO[[i]] = paste("ind_COB_", i, "_PI_fix_GEO", sep = "")
  temp_list = list()
  
  cl <- makeCluster(10) 
  registerDoParallel(cl)
  
  for(geo in 1:59)
  {
    if(i == 14 && geo == 23)
    {
      temp_list[[geo]] = NA
    } else {
      clusterExport(cl,list(fts_area_fix_geo_smooth[[i]][geo], fts_area_fix_geo[[i]][geo]))
      
      temp_list[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% ind_cob_PI_fix_GEO(data_series = get(fts_area_fix_geo_smooth[[i]][geo]), data_series_unsmooth = get(fts_area_fix_geo[[i]][geo]) , fh = ik, nboot = 1000)
    }
  }
  
  stopCluster(cl)
  rm(cl)
  
  assign(ind_cob_PI_all_fix_GEO[[i]], temp_list)
  rm(temp_list)
  
  print(i)
}

# C group

ind_C_group_PI_fix_GEO <- function(data_series = get(fts_region_fix_geo_smooth[[c_index]][geo_ind]),  data_series_unsmooth = get(fts_region_fix_geo[[c_index]][geo_ind]), fh, nboot = 1000)
{
  require(demography)
  
  PI_boot_c = array(0, dim = c(7, nboot, (11-fh)))
  
  for(ij in 1:(11-fh))
  {
    dum_pointwise = ind_PI(dat = extract.years(data_series, 1981:(2000+ij)), dat_unsmooth = extract.years(data_series_unsmooth, 1981:(2000+ij)), fh = fh, nboot = nboot)
    PI_boot_c[,,ij]  = dum_pointwise$boot_sample
  }
  
  return(list(PI_boot = PI_boot_c))
}

test = ind_C_group_PI_fix_GEO(data_series = get(fts_area_fix_geo_smooth[[1]][1]),  data_series_unsmooth = get(fts_area_fix_geo[[1]][1]), fh = 1, nboot = 20)

library(doParallel)

ind_C_group_PI_all_fix_GEO = vector("character", length = 10)
for(i in 1:10)
{
  ind_C_group_PI_all_fix_GEO[[i]] = paste("ind_C_group_", i, "_PI_fix_GEO", sep = "")
  temp_list = list()
  
  cl <- makeCluster(10) 
  registerDoParallel(cl)
  
  for(geo in 1:59)
  {
    if(i == 6 && geo == 23)
    {
      temp_list[[geo]] = NA
    } else {
      clusterExport(cl,list(fts_region_fix_geo_smooth[[i]][geo], fts_region_fix_geo[[i]][geo]))
      
      temp_list[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% ind_C_group_PI_fix_GEO(data_series = get(fts_region_fix_geo_smooth[[i]][geo]), data_series_unsmooth = get(fts_region_fix_geo[[i]][geo]) , fh = ik, nboot = 1000)
    }
  }
  
  stopCluster(cl)
  rm(cl)
  
  assign(ind_C_group_PI_all_fix_GEO[[i]], temp_list)
  rm(temp_list)
  
  print(i)
}

# Global

ind_global_PI_fix_GEO <- function(data_series = get(fts_national_fix_geo_smooth[geo_ind]),  data_series_unsmooth = get(fts_national_fix_geo[geo_ind]), fh, nboot = 1000)
{
  require(demography)
  
  PI_boot_c = array(0, dim = c(7, nboot, (11-fh)))
  
  for(ij in 1:(11-fh))
  {
    dum_pointwise = ind_PI(dat = extract.years(data_series, 1981:(2000+ij)), dat_unsmooth = extract.years(data_series_unsmooth, 1981:(2000+ij)), fh = fh, nboot = nboot)
    PI_boot_c[,,ij]  = dum_pointwise$boot_sample
  }
  
  return(list(PI_boot = PI_boot_c))
}

cl <- makeCluster(10) 
registerDoParallel(cl)

ind_global_PI_all_fix_GEO = list()

for(geo in 1:59)
{
  clusterExport(cl,list(fts_national_fix_geo_smooth[geo], fts_national_fix_geo[geo]))
  
  ind_global_PI_all_fix_GEO[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% ind_global_PI_fix_GEO(data_series = get(fts_national_fix_geo_smooth[geo]), data_series_unsmooth = get(fts_national_fix_geo[geo]) , fh = ik, nboot = 1000)
}


stopCluster(cl)
rm(cl)


####################################
# Interval forecasts reconciliation
####################################

ind_interval_score = function(PI_val, data_series, fh, alpha = 0.8)
{
  require(demography)
  
  test_val = extract.years(data_series, ((2001+fh):2011))$rate$female
  boot_sample = PI_val[[fh]]$PI_boot
  
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



##########
# Fix COB
##########

# Area level
ind_area_interval_score_all = vector("character", length = 47)
for(i in 1:47)
{
  ind_area_interval_score_all[i] = paste("ind_A_", i, "_interval_score_fix_COB", sep = "")
}

for(ij in 1:47)
{
  ind_interval_score_temp = matrix(0, nrow = 10, ncol = 19)
  for(k in 1:19)
  {
    if(ij  == 23 && k == 14)
    {
      ind_interval_score_temp[,k] = rep(0, 10)
    } else {
      for(h in 1:10)
      {
        ind_interval_score_temp[h,k] = ind_interval_score(PI_val = get(ind_area_PI_all_fix_COB[ij])[[k]], data_series = get(fts_area_fix_cob[[ij]][k]), fh = h, alpha = 0.8)$score
      }
    }
  }
  assign(ind_area_interval_score_all[ij], ind_interval_score_temp)
  
  rm(ind_interval_score_temp)
}

# Region level

ind_region_interval_score_all = vector("character", length = 11)
for(i in 1:11)
{
  ind_region_interval_score_all[i] = paste("ind_R_", i, "_interval_score_fix_COB", sep = "")
}

for(ij in 1:11)
{
  ind_interval_score_temp = matrix(0, nrow = 10, ncol = 19)
  for(k in 1:19)
  {
    for(h in 1:10)
    {
      ind_interval_score_temp[h,k] = ind_interval_score(PI_val = get(ind_region_PI_all_fix_COB[ij])[[k]], data_series = get(fts_region_fix_cob[[ij]][k]), fh = h, alpha = 0.8)$score
    }
  }
  assign(ind_region_interval_score_all[ij], ind_interval_score_temp)
  
  rm(ind_interval_score_temp)
}

# National level

ind_national_interval_score_fix_COB = matrix(0, nrow = 10, ncol = 19)
for(k in 1:19)
{
  for(h in 1:10)
  {
    ind_national_interval_score_fix_COB[h,k] = ind_interval_score(PI_val = ind_national_PI_all_fix_COB[[k]], data_series = get(fts_national_fix_cob[k]), fh = h, alpha = 0.8)$score
  }
}

##########
# Fix GEO
##########

# COB level

ind_cob_interval_score_all = vector("character", length = 19)
for(i in 1:19)
{
  ind_cob_interval_score_all[i] = paste("ind_COB_", i, "_interval_score_fix_GEO", sep = "")
}

for(ij in 1:19)
{
  ind_interval_score_temp = matrix(0, nrow = 10, ncol = 59)
  for(k in 1:59)
  {
    if(ij  == 14 && k == 23)
    {
      ind_interval_score_temp[,k] = rep(0, 10)
    } else {
      for(h in 1:10)
      {
        ind_interval_score_temp[h,k] = ind_interval_score(PI_val = get(ind_cob_PI_all_fix_GEO[ij])[[k]], data_series = get(fts_area_fix_geo[[ij]][k]), fh = h, alpha = 0.8)$score
      }
    }
  }
  assign(ind_cob_interval_score_all[ij], ind_interval_score_temp)
  
  rm(ind_interval_score_temp)
}

# C group level

ind_c_group_interval_score_all = vector("character", length = 10)
for(i in 1:10)
{
  ind_c_group_interval_score_all[i] = paste("ind_C_group_", i, "_interval_score_fix_GEO", sep = "")
}

for(ij in 1:10)
{
  ind_interval_score_temp = matrix(0, nrow = 10, ncol = 59)
  for(k in 1:59)
  {
    if(ij  == 6 && k == 23)
    {
      ind_interval_score_temp[,k] = rep(0, 10)
    } else {
      for(h in 1:10)
      {
        ind_interval_score_temp[h,k] = ind_interval_score(PI_val = get(ind_C_group_PI_all_fix_GEO[ij])[[k]], data_series = get(fts_region_fix_geo[[ij]][k]), fh = h, alpha = 0.8)$score
      }
    }
  }
  assign(ind_c_group_interval_score_all[ij], ind_interval_score_temp)
  
  rm(ind_interval_score_temp)
}

# Global level

ind_global_interval_score_fix_GEO = matrix(0, nrow = 10, ncol = 59)
for(k in 1:59)
{
  for(h in 1:10)
  {
    ind_global_interval_score_fix_GEO[h,k] = ind_interval_score(PI_val = ind_global_PI_all_fix_GEO[[k]], data_series = get(fts_national_fix_geo[k]), fh = h, alpha = 0.8)$score
  }
}


#########################################################################################
# All-level bootstrapped base forecasts (B = 1000) using multivariate forecasting method
#########################################################################################

##########
# Fix COB
##########

ind_A23_PI_fix_COB[[14]] = list()
for(h in 1:9)
{
  ind_A23_PI_fix_COB[[14]][[h]] = list(PI_boot = array(rep(1e-6, 7000*(11-h)), dim = c(7, 1000, (11-h))))
}
ind_A23_PI_fix_COB[[14]][[10]] = list(PI_boot = array(rep(1e-6, 7000*(1)), dim = c(7, 1000, 1)))

ind_PI_all_level_fix_cob <- function(fh, cob)
{
  ind_all_level_fh_PI_fix_cob = array(0, dim = c(7, (11-fh), 59, 1000))

  for(ij in 1:1000)
  {
    ind_all_level_fh_PI_fix_cob[,,,ij] = cbind(ind_national_PI_all_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               
                                               ind_R1_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_R2_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_R3_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_R4_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_R5_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_R6_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_R7_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_R8_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_R9_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_R10_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_R11_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               
                                               ind_A1_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A2_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A3_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A4_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_A5_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A6_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A7_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A8_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_A9_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A10_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A11_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A12_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_A13_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A14_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A15_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A16_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_A17_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A18_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A19_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A20_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_A21_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A22_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A23_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A24_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_A25_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A26_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A27_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A28_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_A29_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A30_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A31_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A32_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_A33_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A34_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A35_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A36_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_A37_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A38_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A39_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A40_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_A41_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A42_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A43_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A44_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,],
                                               ind_A45_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A46_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,], ind_A47_PI_fix_COB[[cob]][[fh]]$PI_boot[,ij,])
  }
  
  return(ind_all_level_fh_PI_fix_cob)
}

test = ind_PI_all_level_fix_cob(fh = 10, cob = 14)


# one- to 10-step-ahead bootstrapped base forecasts (B = 1000)

ind_PI_base_fix_cob = vector("character", length = 19)
for(i in 1:19)
{
  ind_PI_base_fix_cob[i] = paste("ind_PI_base_", i, "_fix_COB", sep = "")
}

for(cob in 1:19)
{
  temp_list = list()
  for(h in 1:10)
  {
    temp_list[[h]] = ind_PI_all_level_fix_cob(fh = h, cob = cob)
  }
  
  assign(ind_PI_base_fix_cob[cob], temp_list)
  rm(temp_list)
}

# All-level bootstrapped grouped forecasts

ind_recon_PI_fix_cob <- function(kj, age, hier_method = c("BU", "comb_OLS", "mint"), cob_ind)
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
        hier_fore[,ik,ij] = hier %*% (get(paste("ind_PI_base_", cob_ind, "_fix_COB", sep = "")))[[kj]][age_ind,ik,13:59,ij]
      }
    }
    
    if(hier_method == "comb_OLS")
    {
      for(ij in 1:1000)
      {
        hier_fore[,ik,ij] = hier %*% ginv(t(hier) %*% hier) %*% t(hier) %*% (get(paste("ind_PI_base_", cob_ind, "_fix_COB", sep = "")))[[kj]][age_ind,ik,,ij]
      }
    }
    
    if(hier_method == "mint")
    {
      ginv = MASS:::ginv
      wh = ind_wh_fun_fix_cob(kj = kj, age = age, cob_ind = cob_ind)
      for(ij in 1:1000)
      {
        hier_fore[,ik,ij] = hier %*% ginv(t(hier) %*% ginv(wh) %*% hier) %*% t(hier) %*% ginv(wh) %*% (get(paste("ind_PI_base_", cob_ind, "_fix_COB", sep = "")))[[kj]][age_ind,ik,,ij]
      }
    }
  }
  return(hier_fore)
}

test = ind_recon_PI_fix_cob(kj = 10, age = age_all[4], hier_method = "comb_OLS", cob_ind = 14)

# Calculation of reconciled interval predictions (hier_method = "BU")
ind_BU_recon_PI_fix_cob = list()

for(cob in 1:19)
{
  ind_BU_recon_PI_fix_cob[[cob]] = list()
  for(h in 1:10)
  {
    ind_BU_recon_PI_fix_cob[[cob]][[h]] = array(NA, dim = c(7, 59, (11-h), 1000))
    
    for(ij in 1:7)
    {
      ind_BU_recon_PI_fix_cob[[cob]][[h]][ij,,,] = ind_recon_PI_fix_cob(kj = h, age = age_all[ij], hier_method = "BU", cob_ind = cob)
    }
  }
}


# Calculation of reconciled interval predictions (hier_method = "comb_OLS")
ind_OLS_recon_PI_fix_cob = list()

for(cob in 1:19)
{
  ind_OLS_recon_PI_fix_cob[[cob]] = list()
  for(h in 1:10)
  {
    ind_OLS_recon_PI_fix_cob[[cob]][[h]] = array(NA, dim = c(7, 59, (11-h), 1000))
    
    for(ij in 1:7)
    {
      ind_OLS_recon_PI_fix_cob[[cob]][[h]][ij,,,] = ind_recon_PI_fix_cob(kj = h, age = age_all[ij], hier_method = "comb_OLS", cob_ind = cob)
    }
  }
}

# Calculation of reconciled interval predictions (hier_method = "mint")
ind_MINT_recon_PI_fix_cob = list()

for(cob in 1:19)
{
  ind_MINT_recon_PI_fix_cob[[cob]] = list()
  for(h in 1:10)
  {
    ind_MINT_recon_PI_fix_cob[[cob]][[h]] = array(NA, dim = c(7, 59, (11-h), 1000))
    
    for(ij in 1:7)
    {
      ind_MINT_recon_PI_fix_cob[[cob]][[h]][ij,,,] = ind_recon_PI_fix_cob(kj = h, age = age_all[ij], hier_method = "mint", cob_ind = cob)
    }
  }
}

ind_BU_recon_interval_score_fix_cob = ind_OLS_recon_interval_score_fix_cob = ind_MINT_recon_interval_score_fix_cob = list()
for(cob in 1:19)
{
  # BU
  ind_BU_recon_interval_score_fix_cob[[cob]] = interval_score_recon_fix_cob(PI_val_list = ind_BU_recon_PI_fix_cob, cob_ind = cob)$interval_score_all
  
  # OLS
  ind_OLS_recon_interval_score_fix_cob[[cob]] = interval_score_recon_fix_cob(PI_val_list = ind_OLS_recon_PI_fix_cob, cob_ind = cob)$interval_score_all
  
  # MINT
  ind_MINT_recon_interval_score_fix_cob[[cob]] = interval_score_recon_fix_cob(PI_val_list = ind_MINT_recon_PI_fix_cob, cob_ind = cob)$interval_score_all
}

# save(ind_BU_recon_PI_fix_cob, file = "ind_BU_recon_PI_fix_cob.RData")
# save(ind_OLS_recon_PI_fix_cob, file = "ind_OLS_recon_PI_fix_cob.RData")
# save(ind_MINT_recon_PI_fix_cob, file = "ind_MINT_recon_PI_fix_cob.RData")

rm(ind_BU_recon_PI_fix_cob, ind_OLS_recon_PI_fix_cob, ind_MINT_recon_PI_fix_cob)

load(file = "ind_BU_recon_PI_fix_cob.RData")
load(file = "ind_OLS_recon_PI_fix_cob.RData")
load(file = "ind_MINT_recon_PI_fix_cob.RData")


##########
# Fix GEO
##########

ind_COB_14_PI_fix_GEO[[23]] = ind_C_group_6_PI_fix_GEO[[23]] = list()
for(h in 1:10)
{
  ind_COB_14_PI_fix_GEO[[23]][[h]] = ind_C_group_6_PI_fix_GEO[[23]][[h]] = list(PI_boot = array(rep(1e-6, 7000*(11-h)), dim = c(7, 1000, (11-h))))
}

ind_PI_all_level_fix_geo <- function(fh, geo)
{
  ind_all_level_fh_PI_fix_geo = array(0, dim = c(7, (11-fh), 30, 1000))
  
  for(ij in 1:1000)
  {
    ind_all_level_fh_PI_fix_geo[,,,ij] = cbind(ind_global_PI_all_fix_GEO[[geo]][[fh]]$PI_boot[,ij,],
                                               
                                               ind_COB_1_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_COB_2_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_COB_3_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], 
                                               ind_COB_4_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_COB_5_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_COB_6_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,],
                                               ind_COB_7_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_COB_8_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_COB_9_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,],
                                               ind_COB_10_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_COB_11_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_COB_12_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,],
                                               ind_COB_13_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_COB_14_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_COB_15_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,],
                                               ind_COB_16_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_COB_17_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_COB_18_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,],
                                               ind_COB_19_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,],
                                               
                                               ind_C_group_1_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_C_group_2_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_C_group_3_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,],
                                               ind_C_group_4_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_C_group_5_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_C_group_6_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,],
                                               ind_C_group_7_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_C_group_8_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,], ind_C_group_9_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,],
                                               ind_C_group_10_PI_fix_GEO[[geo]][[fh]]$PI_boot[,ij,])
  }
  
  return(ind_all_level_fh_PI_fix_geo)
}

test = ind_PI_all_level_fix_geo(fh = 1, geo = 23)

# one- to 10-step-ahead bootstrapped base forecasts (B = 1000)

ind_PI_base_fix_geo = vector("character", length = 59)
for(i in 1:59)
{
  ind_PI_base_fix_geo[i] = paste("ind_PI_base_", i, "_fix_GEO", sep = "")
}

for(geo in 1:59)
{
  temp_list = list()
  for(h in 1:10)
  {
    temp_list[[h]] = ind_PI_all_level_fix_geo(fh = h, geo = geo)
  }
  
  assign(ind_PI_base_fix_geo[geo], temp_list)
  rm(temp_list)
}

# All-level bootstrapped grouped forecasts

ind_recon_PI_fix_geo <- function(kj, age, hier_method = c("BU", "comb_OLS", "mint"), geo_ind)
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
        hier_fore[,ik,ij] = hier %*% (get(paste("ind_PI_base_", geo_ind, "_fix_GEO", sep = "")))[[kj]][age_ind,ik,12:30,ij]
      }
    }
    
    if(hier_method == "comb_OLS")
    {
      for(ij in 1:1000)
      {
        hier_fore[,ik,ij] = hier %*% ginv(t(hier) %*% hier) %*% t(hier) %*% (get(paste("ind_PI_base_", geo_ind, "_fix_GEO", sep = "")))[[kj]][age_ind,ik,,ij]
      }
    }
    
    if(hier_method == "mint")
    {
      ginv = MASS:::ginv
      wh = ind_wh_fun_fix_geo(kj = kj, age = age, geo_ind = geo_ind)
      for(ij in 1:1000)
      {
        hier_fore[,ik,ij] = hier %*% ginv(t(hier) %*% ginv(wh) %*% hier) %*% t(hier) %*% ginv(wh) %*% (get(paste("ind_PI_base_", geo_ind, "_fix_GEO", sep = "")))[[kj]][age_ind,ik,,ij]
      }
    }
  }
  return(hier_fore)
}

test = ind_recon_PI_fix_geo(kj = 10, age = age_all[4], hier_method = "mint", geo_ind = 23)

# Calculation of reconciled interval predictions (hier_method = "BU")
ind_BU_recon_PI_fix_geo = list()

for(geo in 1:59)
{
  ind_BU_recon_PI_fix_geo[[geo]] = list()
  for(h in 1:10)
  {
    ind_BU_recon_PI_fix_geo[[geo]][[h]] = array(NA, dim = c(7, 30, (11-h), 1000))
    
    for(ij in 1:7)
    {
      ind_BU_recon_PI_fix_geo[[geo]][[h]][ij,,,] = ind_recon_PI_fix_geo(kj = h, age = age_all[ij], hier_method = "BU", geo_ind = geo)
    }
  }
}


# Calculation of reconciled interval predictions (hier_method = "comb_OLS")
ind_OLS_recon_PI_fix_geo = list()

for(geo in 1:59)
{
  ind_OLS_recon_PI_fix_geo[[geo]] = list()
  for(h in 1:10)
  {
    ind_OLS_recon_PI_fix_geo[[geo]][[h]] = array(NA, dim = c(7, 30, (11-h), 1000))
    
    for(ij in 1:7)
    {
      ind_OLS_recon_PI_fix_geo[[geo]][[h]][ij,,,] = ind_recon_PI_fix_geo(kj = h, age = age_all[ij], hier_method = "comb_OLS", geo_ind = geo)
    }
  }
}

# Calculation of reconciled interval predictions (hier_method = "mint")
ind_MINT_recon_PI_fix_geo = list()

for(geo in 1:59)
{
  ind_MINT_recon_PI_fix_geo[[geo]] = list()
  for(h in 1:10)
  {
    ind_MINT_recon_PI_fix_geo[[geo]][[h]] = array(NA, dim = c(7, 30, (11-h), 1000))
    
    for(ij in 1:7)
    {
      ind_MINT_recon_PI_fix_geo[[geo]][[h]][ij,,,] = ind_recon_PI_fix_geo(kj = h, age = age_all[ij], hier_method = "mint", geo_ind = geo)
    }
  }
}


ind_BU_recon_interval_score_fix_geo = ind_OLS_recon_interval_score_fix_geo = ind_MINT_recon_interval_score_fix_geo = list()
for(geo in 1:59)
{
  # BU
  ind_BU_recon_interval_score_fix_geo[[geo]] = interval_score_recon_fix_geo(PI_val_list = ind_BU_recon_PI_fix_geo, geo_ind = geo)$interval_score_all
  
  # OLS
  ind_OLS_recon_interval_score_fix_geo[[geo]] = interval_score_recon_fix_geo(PI_val_list = ind_OLS_recon_PI_fix_geo, geo_ind = geo)$interval_score_all
  
  # MINT
  ind_MINT_recon_interval_score_fix_geo[[geo]] = interval_score_recon_fix_geo(PI_val_list = ind_MINT_recon_PI_fix_geo, geo_ind = geo)$interval_score_all
}

# save(ind_BU_recon_PI_fix_geo, file = "ind_BU_recon_PI_fix_geo.RData")
# save(ind_OLS_recon_PI_fix_geo, file = "ind_OLS_recon_PI_fix_geo.RData")
# save(ind_MINT_recon_PI_fix_geo, file = "ind_MINT_recon_PI_fix_geo.RData")

rm(ind_BU_recon_PI_fix_geo, ind_OLS_recon_PI_fix_geo, ind_MINT_recon_PI_fix_geo)

load(file = "ind_BU_recon_PI_fix_geo.RData")
load(file = "ind_OLS_recon_PI_fix_geo.RData")
load(file = "ind_MINT_recon_PI_fix_geo.RData")


#############################
# Summary of interval scores
#############################

# Fix cob, independent forecasts

ind_IS_area_fix_cob = array(0, dim = c(10,47,19))
for(ij in 1:47)
{
  for(cob in 1:19)
  {
    ind_IS_area_fix_cob[,ij,cob] = get(ind_area_interval_score_all[ij])[,cob]
  }
}

ind_IS_region_fix_cob = array(0, dim = c(10, 11, 19))
for(ij in 1:11)
{
  for(cob in 1:19)
  {
    ind_IS_region_fix_cob[,ij,cob] = get(ind_region_interval_score_all[ij])[,cob]
  }
}

ind_IS_national_fix_cob = array(0, dim = c(10, 1, 19))
for(cob in 1:19)
{
  ind_IS_national_fix_cob[,1,cob] = ind_national_interval_score_fix_COB[,cob]
}

ind_IS_fix_COB_all = array(0, dim = c(10,3,19))
for(cob in 1:19)
{
  ind_IS_fix_COB_all[,1,cob] = ind_IS_national_fix_cob[,,cob]
  ind_IS_fix_COB_all[,2,cob] = rowMeans(ind_IS_region_fix_cob[,,cob])
  ind_IS_fix_COB_all[,3,cob] = rowMeans(ind_IS_area_fix_cob[,,cob])
}

base_ind_IS_fix_COB_averaged = apply(ind_IS_fix_COB_all, c(1,2), mean)

# Fix geo, independent forecasts

ind_IS_cob_fix_geo = array(0, dim = c(10,19,59))
for(ij in 1:19)
{
  for(geo in 1:59)
  {
    ind_IS_cob_fix_geo[,ij,geo] = get(ind_cob_interval_score_all[ij])[,geo]
  }
}

ind_IS_c_group_fix_geo = array(0, dim = c(10, 10, 59))
for(ij in 1:10)
{
  for(geo in 1:59)
  {
    ind_IS_c_group_fix_geo[,ij,geo] = get(ind_c_group_interval_score_all[ij])[,geo]
  }
}

ind_IS_fix_GEO_all = array(0, dim = c(10,3,59))
for(geo in 1:59)
{
  ind_IS_fix_GEO_all[,1,geo] = ind_global_interval_score_fix_GEO[,geo]
  ind_IS_fix_GEO_all[,2,geo] = rowMeans(ind_IS_cob_fix_geo[,,geo])
  ind_IS_fix_GEO_all[,3,geo] = rowMeans(ind_IS_c_group_fix_geo[,,geo])
}

base_ind_IS_fix_GEO_averaged = apply(ind_IS_fix_GEO_all, c(1,2), mean)

# fix cob; reconciled

ind_BU_IS_fix_cob = ind_OLS_IS_fix_cob = ind_MINT_IS_fix_cob = array(0, dim = c(10, 3, 19))
for(cob in 1:19)
{
  # BU
  ind_BU_IS_fix_cob[,1,cob] = ind_BU_recon_interval_score_fix_cob[[cob]][,1]
  ind_BU_IS_fix_cob[,2,cob] = rowMeans(ind_BU_recon_interval_score_fix_cob[[cob]][,2:12])
  ind_BU_IS_fix_cob[,3,cob] = rowMeans(ind_BU_recon_interval_score_fix_cob[[cob]][,13:59])
  
  # OP
  ind_OLS_IS_fix_cob[,1,cob] = ind_OLS_recon_interval_score_fix_cob[[cob]][,1]
  ind_OLS_IS_fix_cob[,2,cob] = rowMeans(ind_OLS_recon_interval_score_fix_cob[[cob]][,2:12])
  ind_OLS_IS_fix_cob[,3,cob] = rowMeans(ind_OLS_recon_interval_score_fix_cob[[cob]][,13:59])
  
  # MINT
  ind_MINT_IS_fix_cob[,1,cob] = ind_MINT_recon_interval_score_fix_cob[[cob]][,1]
  ind_MINT_IS_fix_cob[,2,cob] = rowMeans(ind_MINT_recon_interval_score_fix_cob[[cob]][,2:12])
  ind_MINT_IS_fix_cob[,3,cob] = rowMeans(ind_MINT_recon_interval_score_fix_cob[[cob]][,13:59])
}

ind_BU_IS_fix_cob_averaged = apply(ind_BU_IS_fix_cob, c(1,2), mean)
ind_OLS_IS_fix_cob_averaged = apply(ind_OLS_IS_fix_cob, c(1,2), mean)
ind_MINT_IS_fix_cob_averaged = apply(ind_MINT_IS_fix_cob, c(1,2), mean)


# fix geo; reconciled
ind_BU_IS_fix_geo = ind_OLS_IS_fix_geo = ind_MINT_IS_fix_geo = array(0, dim = c(10, 3, 59))
for(geo in 1:59)
{
  # BU
  ind_BU_IS_fix_geo[,1,geo] = ind_BU_recon_interval_score_fix_geo[[geo]][,1]
  ind_BU_IS_fix_geo[,2,geo] = rowMeans(ind_BU_recon_interval_score_fix_geo[[geo]][,2:11])
  ind_BU_IS_fix_geo[,3,geo] = rowMeans(ind_BU_recon_interval_score_fix_geo[[geo]][,12:30])
  
  # OP
  ind_OLS_IS_fix_geo[,1,geo] = ind_OLS_recon_interval_score_fix_geo[[geo]][,1]
  ind_OLS_IS_fix_geo[,2,geo] = rowMeans(ind_OLS_recon_interval_score_fix_geo[[geo]][,2:11])
  ind_OLS_IS_fix_geo[,3,geo] = rowMeans(ind_OLS_recon_interval_score_fix_geo[[geo]][,12:30])
  
  # MINT
  ind_MINT_IS_fix_geo[,1,geo] = ind_MINT_recon_interval_score_fix_geo[[geo]][,1]
  ind_MINT_IS_fix_geo[,2,geo] = rowMeans(ind_MINT_recon_interval_score_fix_geo[[geo]][,2:11])
  ind_MINT_IS_fix_geo[,3,geo] = rowMeans(ind_MINT_recon_interval_score_fix_geo[[geo]][,12:30])
}

ind_BU_IS_fix_geo_averaged = apply(ind_BU_IS_fix_geo, c(1,2), mean)
ind_OLS_IS_fix_geo_averaged = apply(ind_OLS_IS_fix_geo, c(1,2), mean)
ind_MINT_IS_fix_geo_averaged = apply(ind_MINT_IS_fix_geo, c(1,2), mean)
