##########################################################
# Reconciliation of point forecasts using various methods
##########################################################

##########
# Fix COB
##########

recon_hier_fix_cob <- function(national_forc, region_forc, area_forc, age, kj, hier_method = c("BU", "comb_OLS", "comb_GLS", "mint"), cob_ind)
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
      wh = wh_fun_fix_cob(kj = kj, age = age, cob_ind = cob_ind)
      hier_fore[,ik] = hier %*% ginv(t(hier) %*% ginv(wh) %*% hier) %*% t(hier) %*% ginv(wh) %*% hier_rate[,ik]
    }
  }
  return(hier_fore)
}


# point forecast errors function (shared by both univariate and multivariate functional time series forecasting methods)

BU_optim_err_fix_cob <- function(ik, hier_method, national_forc, region_forc, area_forc, cob_ind)
{
  age_all = c(15, 20, 25, 30, 35, 40, 45)
  
  me = ftsa:::me; mae = ftsa:::mae; rmse = ftsa:::rmse; mase = ftsa:::mase
  BU_optim_hier_comb = array(NA, dim = c(7, 59,(11-ik)))
  
  for(i in 1:7)
  {
    BU_optim_hier_comb[i,,] = recon_hier_fix_cob(national_forc, region_forc, area_forc,      
                                                      age = age_all[i], kj = ik, hier_method = hier_method,
                                                      cob_ind = cob_ind)
  }
  
  #####################################
  # Errors, including ME, MAE and RMSE
  #####################################
  
  # Level 1 (Total)
  mae_national = mae(BU_optim_hier_comb[,1,], extract.years(get(fts_national_fix_cob[cob_ind]), years = (2001+ik):2011)$rate$female)
  mase_national = mase(forecast = BU_optim_hier_comb[,1,], outsampletrue = extract.years(get(fts_national_fix_cob[cob_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_national_fix_cob[cob_ind]), years = 1986:(2001+ik-1))$rate$female)
  rmse_national = rmse(BU_optim_hier_comb[,1,], extract.years(get(fts_national_fix_cob[cob_ind]), years = (2001+ik):2011)$rate$female)
  
  # Level 2 (Region)
 mae_region = mase_region = rmse_region = rep(0, 11)
  for(iw in 1:11)
  {
    mae_region[iw]  = mae(BU_optim_hier_comb[,(iw+1),],  extract.years(get(fts_region_fix_cob[[iw]][cob_ind]), years = (2001+ik):2011)$rate$female)
    mase_region[iw]  = mase(forecast = BU_optim_hier_comb[,(iw+1),], outsampletrue = extract.years(get(fts_region_fix_cob[[iw]][cob_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_region_fix_cob[[iw]][cob_ind]), years = 1986:(2001+ik-1))$rate$female)
    rmse_region[iw] = rmse(BU_optim_hier_comb[,(iw+1),], extract.years(get(fts_region_fix_cob[[iw]][cob_ind]), years = (2001+ik):2011)$rate$female)
  }
  
  # Level 3 (Area)
  mae_area = mase_area= rmse_area = rep(0, 47)
  
  for(iwk in 1:47)
  {
    if(iwk != 23)
    {
      mae_area[iwk]  = mae(BU_optim_hier_comb[,(iwk+12),],   extract.years(get(fts_area_fix_cob[[iwk]][cob_ind]), years = (2001+ik):2011)$rate$female)
      mase_area[iwk]  = mase(forecast = BU_optim_hier_comb[,(iwk+12),], outsampletrue = extract.years(get(fts_area_fix_cob[[iwk]][cob_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_area_fix_cob[[iwk]][cob_ind]), years = 1986:(2001+ik-1))$rate$female)
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


########################################
# Point forecast errors; fmethod = "BU" 
########################################

national_mae_BU = national_mase_BU = national_rmse_BU = array(0, dim = c(10, 1, 19))
region_mae_BU = region_mase_BU = region_rmse_BU = array(0, dim = c(10, 11, 19))
area_mae_BU = area_mase_BU = area_rmse_BU = array(0, dim = c(10, 47, 19))

for(cob in 1:19)
{
  for(ikw in 1:10)
  {
    dum = BU_optim_err_fix_cob(ik = ikw, hier_method = "BU", 
                               national_forc = ind_national_forecast_fix_COB, 
                               region_forc = region_forecast_fix_COB, 
                               area_forc = area_forecast_fix_COB, 
                               cob_ind = cob)
    
    # MAE
    national_mae_BU[ikw,,cob] = dum$mae_national
    region_mae_BU[ikw,,cob] = dum$mae_region
    area_mae_BU[ikw,,cob] = dum$mae_area
    
    # MASE
    national_mase_BU[ikw,,cob] = dum$mase_national
    region_mase_BU[ikw,,cob] = dum$mase_region
    area_mase_BU[ikw,,cob] = dum$mase_area
    
    # RMSE
    national_rmse_BU[ikw,,cob] = dum$rmse_national
    region_rmse_BU[ikw,,cob] = dum$rmse_region
    area_rmse_BU[ikw,,cob] = dum$rmse_area
    
  }
}

# Summary of results

fix_cob_BU_mae = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  fix_cob_BU_mae[,,i] = cbind(national_mae_BU[,,i], rowMeans(region_mae_BU[,,i]), rowMeans(area_mae_BU[,,i]))
}
fix_cob_BU_mae_averaged = apply(fix_cob_BU_mae, c(1,2), mean)

###
area_mase_BU_copy = area_mase_BU
for(i in 1:19)
{
  area_mase_BU[,,i][,which(colMeans(area_mase_BU[,,i]) > 300)] = 0
}
###

fix_cob_BU_mase = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  fix_cob_BU_mase[,,i] = cbind(national_mase_BU[,,i], rowMeans(region_mase_BU[,,i]), rowMeans(area_mase_BU[,,i]))
}
fix_cob_BU_mase_averaged = apply(fix_cob_BU_mase, c(1,2), mean)


fix_cob_BU_rmse = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  fix_cob_BU_rmse[,,i] = cbind(national_rmse_BU[,,i], rowMeans(region_rmse_BU[,,i]), rowMeans(area_rmse_BU[,,i]))
}
fix_cob_BU_rmse_averaged = apply(fix_cob_BU_rmse, c(1,2), mean)

##############################################
# Point forecast errors; fmethod = "comb_OLS" 
##############################################

national_mae_OP = national_mase_OP = national_rmse_OP = array(0, dim = c(10, 1, 19))
region_mae_OP = region_mase_OP = region_rmse_OP = array(0, dim = c(10, 11, 19))
area_mae_OP = area_mase_OP = area_rmse_OP = array(0, dim = c(10, 47, 19))

for(cob in 1:19)
{
  for(ikw in 1:10)
  {
    dum = BU_optim_err_fix_cob(ik = ikw, hier_method = "comb_OLS", 
                               national_forc = ind_national_forecast_fix_COB, 
                               region_forc = region_forecast_fix_COB, 
                               area_forc = area_forecast_fix_COB, 
                               cob_ind = cob)
    
    # MAE
    national_mae_OP[ikw,,cob] = dum$mae_national
    region_mae_OP[ikw,,cob] = dum$mae_region
    area_mae_OP[ikw,,cob] = dum$mae_area
    
    # MASE
    national_mase_OP[ikw,,cob] = dum$mase_national
    region_mase_OP[ikw,,cob] = dum$mase_region
    area_mase_OP[ikw,,cob] = dum$mase_area
    
    # RMSE
    national_rmse_OP[ikw,,cob] = dum$rmse_national
    region_rmse_OP[ikw,,cob] = dum$rmse_region
    area_rmse_OP[ikw,,cob] = dum$rmse_area
    
  }
}

# Summary of results

fix_cob_OP_mae = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  fix_cob_OP_mae[,,i] = cbind(national_mae_OP[,,i], rowMeans(region_mae_OP[,,i]), rowMeans(area_mae_OP[,,i]))
}
fix_cob_OP_mae_averaged = apply(fix_cob_OP_mae, c(1,2), mean)

###
area_mase_OP_copy = area_mase_OP
for(i in 1:19)
{
  area_mase_OP[,,i][,which(colMeans(area_mase_OP[,,i]) > 300)] = 0
}
###

fix_cob_OP_mase = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  fix_cob_OP_mase[,,i] = cbind(national_mase_OP[,,i], rowMeans(region_mase_OP[,,i]), rowMeans(area_mase_OP[,,i]))
}
fix_cob_OP_mase_averaged = apply(fix_cob_OP_mase, c(1,2), mean)


fix_cob_OP_rmse = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  fix_cob_OP_rmse[,,i] = cbind(national_rmse_OP[,,i], rowMeans(region_rmse_OP[,,i]), rowMeans(area_rmse_OP[,,i]))
}
fix_cob_OP_rmse_averaged = apply(fix_cob_OP_rmse, c(1,2), mean)


##########################################
# Point forecast errors; fmethod = "mint" 
##########################################

national_mae_MINT = national_mase_MINT = national_rmse_MINT = array(0, dim = c(10, 1, 19))
region_mae_MINT = region_mase_MINT = region_rmse_MINT = array(0, dim = c(10, 11, 19))
area_mae_MINT = area_mase_MINT = area_rmse_MINT = array(0, dim = c(10, 47, 19))

for(cob in 1:19)
{
  for(ikw in 1:10)
  {
    dum = BU_optim_err_fix_cob(ik = ikw, hier_method = "mint", 
                               national_forc = ind_national_forecast_fix_COB, 
                               region_forc = region_forecast_fix_COB, 
                               area_forc = area_forecast_fix_COB, 
                               cob_ind = cob)
    
    # MAE
    national_mae_MINT[ikw,,cob] = dum$mae_national
    region_mae_MINT[ikw,,cob] = dum$mae_region
    area_mae_MINT[ikw,,cob] = dum$mae_area
    
    # MAE
    national_mase_MINT[ikw,,cob] = dum$mase_national
    region_mase_MINT[ikw,,cob] = dum$mase_region
    area_mase_MINT[ikw,,cob] = dum$mase_area
    
    # RMSE
    national_rmse_MINT[ikw,,cob] = dum$rmse_national
    region_rmse_MINT[ikw,,cob] = dum$rmse_region
    area_rmse_MINT[ikw,,cob] = dum$rmse_area
    
  }
}

fix_cob_MINT_mae = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  fix_cob_MINT_mae[,,i] = cbind(national_mae_MINT[,,i], rowMeans(region_mae_MINT[,,i]), rowMeans(area_mae_MINT[,,i]))
}
fix_cob_MINT_mae_averaged = apply(fix_cob_MINT_mae, c(1,2), mean)

###
area_mase_MINT_copy = area_mase_MINT
for(i in 1:19)
{
  area_mase_MINT[,,i][,which(colMeans(area_mase_MINT[,,i]) > 300)] = 0
}
###

fix_cob_MINT_mase = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  fix_cob_MINT_mase[,,i] = cbind(national_mase_MINT[,,i], rowMeans(region_mase_MINT[,,i]), rowMeans(area_mase_MINT[,,i]))
}
fix_cob_MINT_mase_averaged = apply(fix_cob_MINT_mase, c(1,2), mean)


fix_cob_MINT_rmse = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  fix_cob_MINT_rmse[,,i] = cbind(national_rmse_MINT[,,i], rowMeans(region_rmse_MINT[,,i]), rowMeans(area_rmse_MINT[,,i]))
}
fix_cob_MINT_rmse_averaged = apply(fix_cob_MINT_rmse, c(1,2), mean)

##########
# Fix GEO
##########

recon_hier_fix_geo <- function(global_forc, c_group_forc, cob_forc, age, kj, hier_method = c("BU", "comb_OLS", "comb_GLS", "mint"), geo_ind)
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
      wh = wh_fun_fix_geo(kj = kj, age = age, geo_ind = geo_ind)
      hier_fore[,ik] = hier %*% ginv(t(hier) %*% ginv(wh) %*% hier) %*% t(hier) %*% ginv(wh) %*% hier_rate[,ik]
    }
  }
  return(hier_fore)
}


# point forecast errors function (shared by both univariate and multivariate functional time series forecasting methods)

BU_optim_err_fix_geo <- function(ik, hier_method, global_forc, c_group_forc, cob_forc, geo_ind)
{
  age_all = c(15, 20, 25, 30, 35, 40, 45)
  
  me = ftsa:::me; mae = ftsa:::mae; rmse = ftsa:::rmse; mase = ftsa:::mase
  BU_optim_hier_comb = array(NA, dim = c(7, 30,(11-ik)))
  
  for(i in 1:7)
  {
    BU_optim_hier_comb[i,,] = recon_hier_fix_geo(global_forc, c_group_forc, cob_forc,      
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
  mae_cob = mase_cob = rmse_cob = rep(0, 10)
  
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

global_mae_BU = global_mase_BU = global_rmse_BU = array(0, dim = c(10, 1, 59))
c_group_mae_BU = c_group_mase_BU = c_group_rmse_BU = array(0, dim = c(10, 10, 59))
cob_mae_BU = cob_mase_BU = cob_rmse_BU = array(0, dim = c(10, 19, 59))

for(geo in 1:59)
{
  for(ikw in 1:10)
  {
    dum = BU_optim_err_fix_geo(ik = ikw, hier_method = "BU", 
                               global_forc = global_forecast_fix_GEO, 
                               c_group_forc = C_group_forecast_fix_GEO, 
                               cob_forc = cob_forecast_fix_GEO, 
                               geo_ind = geo)
    
    # MAE
    global_mae_BU[ikw,,geo] = dum$mae_global
    c_group_mae_BU[ikw,,geo] = dum$mae_c_group
    cob_mae_BU[ikw,,geo] = dum$mae_cob
    
    # MAE
    global_mase_BU[ikw,,geo] = dum$mase_global
    c_group_mase_BU[ikw,,geo] = dum$mase_c_group
    cob_mase_BU[ikw,,geo] = dum$mase_cob
    
    # RMSE
    global_rmse_BU[ikw,,geo] = dum$rmse_global
    c_group_rmse_BU[ikw,,geo] = dum$rmse_c_group
    cob_rmse_BU[ikw,,geo] = dum$rmse_cob
    
  }
}

# Summary of results

fix_geo_BU_mae = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  fix_geo_BU_mae[,,i] = cbind(global_mae_BU[,,i], rowMeans(c_group_mae_BU[,,i]), rowMeans(cob_mae_BU[,,i]))
}
fix_geo_BU_mae_averaged = apply(fix_geo_BU_mae, c(1,2), mean)

###
c_group_mase_BU_copy = c_group_mase_BU
for(i in 1:59)
{
  c_group_mase_BU[,,i][,which(colMeans(c_group_mase_BU[,,i]) > 35)] = 0
}

cob_mase_BU_copy = cob_mase_BU
for(i in 1:59)
{
  cob_mase_BU[,,i][,which(colMeans(cob_mase_BU[,,i]) > 300)] = 0
}
###

fix_geo_BU_mase = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  fix_geo_BU_mase[,,i] = cbind(global_mase_BU[,,i], rowMeans(c_group_mase_BU[,,i]), rowMeans(cob_mase_BU[,,i]))
}
fix_geo_BU_mase_averaged = apply(fix_geo_BU_mase, c(1,2), mean, na.rm = TRUE)


fix_geo_BU_rmse = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  fix_geo_BU_rmse[,,i] = cbind(global_rmse_BU[,,i], rowMeans(c_group_rmse_BU[,,i]), rowMeans(cob_rmse_BU[,,i]))
}
fix_geo_BU_rmse_averaged = apply(fix_geo_BU_rmse, c(1,2), mean)

##############################################
# Point forecast errors; fmethod = "comb_OLS" 
##############################################

global_mae_OP = global_mase_OP = global_rmse_OP = array(0, dim = c(10, 1, 59))
c_group_mae_OP = c_group_mase_OP = c_group_rmse_OP = array(0, dim = c(10, 10, 59))
cob_mae_OP = cob_mase_OP = cob_rmse_OP = array(0, dim = c(10, 19, 59))

for(geo in 1:59)
{
  for(ikw in 1:10)
  {
    dum = BU_optim_err_fix_geo(ik = ikw, hier_method = "comb_OLS", 
                               global_forc = global_forecast_fix_GEO, 
                               c_group_forc = C_group_forecast_fix_GEO, 
                               cob_forc = cob_forecast_fix_GEO, 
                               geo_ind = geo)
    
    # MAE
    global_mae_OP[ikw,,geo] = dum$mae_global
    c_group_mae_OP[ikw,,geo] = dum$mae_c_group
    cob_mae_OP[ikw,,geo] = dum$mae_cob
    
    # MAE
    global_mase_OP[ikw,,geo] = dum$mase_global
    c_group_mase_OP[ikw,,geo] = dum$mase_c_group
    cob_mase_OP[ikw,,geo] = dum$mase_cob
    
    # RMSE
    global_rmse_OP[ikw,,geo] = dum$rmse_global
    c_group_rmse_OP[ikw,,geo] = dum$rmse_c_group
    cob_rmse_OP[ikw,,geo] = dum$rmse_cob
    
  }
}

fix_geo_OP_mae = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  fix_geo_OP_mae[,,i] = cbind(global_mae_OP[,,i], rowMeans(c_group_mae_OP[,,i]), rowMeans(cob_mae_OP[,,i]))
}
fix_geo_OP_mae_averaged = apply(fix_geo_OP_mae, c(1,2), mean)

###
c_group_mase_OP_copy = c_group_mase_OP
for(i in 1:59)
{
  c_group_mase_OP[,,i][,which(colMeans(c_group_mase_OP[,,i]) > 300)] = 0
}

cob_mase_OP_copy = cob_mase_OP
for(i in 1:59)
{
  cob_mase_OP[,,i][,which(colMeans(cob_mase_OP[,,i]) > 300)] = 0
}
###


fix_geo_OP_mase = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  fix_geo_OP_mase[,,i] = cbind(global_mase_OP[,,i], rowMeans(c_group_mase_OP[,,i]), rowMeans(cob_mase_OP[,,i]))
}
fix_geo_OP_mase_averaged = apply(fix_geo_OP_mase, c(1,2), mean)


fix_geo_OP_rmse = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  fix_geo_OP_rmse[,,i] = cbind(global_rmse_OP[,,i], rowMeans(c_group_rmse_OP[,,i]), rowMeans(cob_rmse_OP[,,i]))
}
fix_geo_OP_rmse_averaged = apply(fix_geo_OP_rmse, c(1,2), mean)


##########################################
# Point forecast errors; fmethod = "mint" 
##########################################

global_mae_MINT = global_mase_MINT = global_rmse_MINT = array(0, dim = c(10, 1, 59))
c_group_mae_MINT = c_group_mase_MINT = c_group_rmse_MINT = array(0, dim = c(10, 10, 59))
cob_mae_MINT = cob_mase_MINT = cob_rmse_MINT = array(0, dim = c(10, 19, 59))

for(geo in 1:59)
{
  for(ikw in 1:10)
  {
    dum = BU_optim_err_fix_geo(ik = ikw, hier_method = "mint", 
                               global_forc = global_forecast_fix_GEO, 
                               c_group_forc = C_group_forecast_fix_GEO, 
                               cob_forc = cob_forecast_fix_GEO,
                               geo_ind = geo)
    
    # MAE
    global_mae_MINT[ikw,,geo] = dum$mae_global
    c_group_mae_MINT[ikw,,geo] = dum$mae_c_group
    cob_mae_MINT[ikw,,geo] = dum$mae_cob
    
    # MAsE
    global_mase_MINT[ikw,,geo] = dum$mase_global
    c_group_mase_MINT[ikw,,geo] = dum$mase_c_group
    cob_mase_MINT[ikw,,geo] = dum$mase_cob
    
    # RMSE
    global_rmse_MINT[ikw,,geo] = dum$rmse_global
    c_group_rmse_MINT[ikw,,geo] = dum$rmse_c_group
    cob_rmse_MINT[ikw,,geo] = dum$rmse_cob
    
  }
}

fix_geo_MINT_mae = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  fix_geo_MINT_mae[,,i] = cbind(global_mae_MINT[,,i], rowMeans(c_group_mae_MINT[,,i]), rowMeans(cob_mae_MINT[,,i]))
}
fix_geo_MINT_mae_averaged = apply(fix_geo_MINT_mae, c(1,2), mean)

###
c_group_mase_MINT_copy = c_group_mase_MINT
for(i in 1:59)
{
  c_group_mase_MINT[,,i][,which(colMeans(c_group_mase_MINT[,,i]) > 300)] = 0
}

cob_mase_MINT_copy = cob_mase_MINT
for(i in 1:59)
{
  cob_mase_MINT[,,i][,which(colMeans(cob_mase_MINT[,,i]) > 300)] = 0
}
###

fix_geo_MINT_mase = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  fix_geo_MINT_mase[,,i] = cbind(global_mase_MINT[,,i], rowMeans(c_group_mase_MINT[,,i]), rowMeans(cob_mase_MINT[,,i]))
}
fix_geo_MINT_mase_averaged = apply(fix_geo_MINT_mase, c(1,2), mean)


fix_geo_MINT_rmse = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  fix_geo_MINT_rmse[,,i] = cbind(global_rmse_MINT[,,i], rowMeans(c_group_rmse_MINT[,,i]), rowMeans(cob_rmse_MINT[,,i]))
}
fix_geo_MINT_rmse_averaged = apply(fix_geo_MINT_rmse, c(1,2), mean)



















