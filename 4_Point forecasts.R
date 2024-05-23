######################################################
# Multivariate functional time series point forecasts 
######################################################

###############################################
# Fix COB; disaggregate by geographical factor
###############################################

# Area level 

area_train_residual_fix_COB = vector("character", length = 47)
for(i in 1:47)
{
  area_train_residual_fix_COB[i] = paste("A", i, "_train_residual_fix_COB", sep = "")
}

area_forecast_fix_COB = vector("character", length = 47)
for(i in 1:47)
{
  area_forecast_fix_COB[i] = paste("A", i, "_forecast_fix_COB", sep = "")
}

area_mae_fix_COB = vector("character", length = 47)
for(i in 1:47)
{
  area_mae_fix_COB[i] = paste("A", i, "_mae_fix_COB", sep = "")
}

area_mase_fix_COB = vector("character", length = 47)
for(i in 1:47)
{
  area_mase_fix_COB[i] = paste("A", i, "_mase_fix_COB", sep = "")
}

area_rmse_fix_COB = vector("character", length = 47)
for(i in 1:47)
{
  area_rmse_fix_COB[i] = paste("A", i, "_rmse_fix_COB", sep = "")
}

# Region level

region_train_residual_fix_COB = vector("character", length = 11)
for(i in 1:11)
{
  region_train_residual_fix_COB[i] = paste("R", i, "_train_residual_fix_COB", sep = "")
}

region_forecast_fix_COB = vector("character", length = 11)
for(i in 1:11)
{
  region_forecast_fix_COB[i] = paste("R", i, "_forecast_fix_COB", sep = "")
}

region_mae_fix_COB = vector("character", length = 11)
for(i in 1:11)
{
  region_mae_fix_COB[i] = paste("R", i, "_mae_fix_COB", sep = "")
}

region_mase_fix_COB = vector("character", length = 11)
for(i in 1:11)
{
  region_mase_fix_COB[i] = paste("R", i, "_mase_fix_COB", sep = "")
}

region_rmse_fix_COB = vector("character", length = 11)
for(i in 1:11)
{
  region_rmse_fix_COB[i] = paste("R", i, "_rmse_fix_COB", sep = "")
}

# National level

national_train_residual_fix_COB = paste("National", "_train_residual_fix_COB", sep = "")
national_forecast_fix_COB = paste("National", "_forcast_fix_COB", sep = "")
national_mae_fix_COB = paste("National", "_mae_fix_COB", sep = "")
national_mase_fix_COB = paste("National", "_mase_fix_COB", sep = "")
national_rmse_fix_COB = paste("National", "_rmse_fix_COB", sep = "")

# function for multivariate forecasting

mfts <- function(dat, pcamethod = c("static", "dynamic"), year_horizon)
{
  require(ftsa)
  n = length(names(dat$rate))
  # if the data contains female, male and total, then we consider only female and male
  # if(n == 3)
  # {
  #   n_pop = 2
  # }else
  # {
  #   n_pop = n
  # }
  
  n_pop = n
  
  rowmeans_object = sd_object = decenter_object = list()
  for(ik in 1:n_pop)
  {
    # compute mean and standard deviation functions
    rowmeans_object[[ik]] = rowMeans(dat$rate[[ik]], na.rm=TRUE)
    sd_object[[ik]] = apply(dat$rate[[ik]], 1, sd, na.rm=TRUE)
    
    # de-center functional data
    decenter_object[[ik]] = t(scale(t(dat$rate[[ik]]), center = TRUE, scale = TRUE))
  }
  comb_object = do.call(rbind, decenter_object)
  
  pcamethod = match.arg(pcamethod)
  
  if (pcamethod == "static")
  {
    ncomp_comb = head(which(cumsum(ftsm(fts(1:nrow(comb_object), comb_object), order = 20)$varprop) >= 0.95), 1)
    fit_ftsm = ftsm(fts(1:nrow(comb_object), comb_object), order = ncomp_comb)
    train_residual = (comb_object * do.call(c,sd_object) + do.call(c, rowmeans_object)) - (fit_ftsm$fitted$y * do.call(c,sd_object) + do.call(c, rowmeans_object))
    fore_ftsm = forecast(fit_ftsm, h = year_horizon)
    fore_res = fore_ftsm$mean$y * do.call(c,sd_object) + do.call(c, rowmeans_object)
  }
  
  
  if (pcamethod == "dynamic")
  {
    data_dum = comb_object
    
    C_0 = long_run_covariance_estimation(data_dum, H = 3, C0 = 3)
    eigen_decomp = eigen(C_0)
    dynamic_order = head(which(cumsum(eigen_decomp$values)/sum(eigen_decomp$values) >= 0.95),1)
    dynamic_basis = as.matrix(eigen_decomp$vectors[,1:dynamic_order])
    dynamic_scores = t(dynamic_basis) %*% data_dum
    
    train_residual = (comb_object * do.call(c,sd_object) + do.call(c, rowmeans_object)) - (dynamic_basis %*%dynamic_scores * do.call(c,sd_object) + do.call(c, rowmeans_object))
    
    # making forecasts
    scores_fit = scores_fore = list()
    fore_ftsm_dyn = matrix(NA, nrow = nrow(data_dum), ncol = year_horizon)
    
    for(ik in 1:dynamic_order)
    {
      scores_fit[[ik]] = auto.arima(dynamic_scores[ik,])
      scores_fore[[ik]] = forecast(scores_fit[[ik]], h = year_horizon)$mean
    }
    
    for(ih in 1:year_horizon)
    {
      fore_ftsm_dyn[,ih] = dynamic_basis%*% unlist(lapply(scores_fore,`[[`,ih))
    }
    
    fore_res = (fore_ftsm_dyn * do.call(c,sd_object) + do.call(c, rowmeans_object))
  }
  
  return(list(fore_res = fore_res, train_residual = train_residual))
}

# Collectively forecast all 47 areas

mfts_area_fix_COB <- function(region_index = 1,  pcamethod = c("static", "dynamic"), year_horizon, cob_ind)
{
  require(ftsa)
  require(demography)
  n_year = 2011 - (year_horizon + 1)
  
  area_ind = region_list_ind[[region_index]]
  
  if(region_index == 11 && cob_ind == 14)
  {
    area_ind = area_ind[-3] # no observation at all for Area 23 COB 14
  }
  
  # region_comb = region_comb_pop = matrix(NA, 31*31, length(area_ind))
  # for(iw in 1:length(area_ind))
  # {
  #   region_comb[,iw] = as.numeric(get(fts_area_fix_cob_smooth[[area_ind[iw]]][cob_ind])$rate$female[3:33,])
  #   region_comb_pop[,iw] = as.numeric(get(fts_area_fix_cob_smooth[[area_ind[iw]]][cob_ind])$pop$female[3:33,])
  # }
  # region_comb_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), region_comb)
  # region_comb_pop_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), region_comb_pop)
  # 
  # colnames(region_comb_v2) = colnames(region_comb_pop_v2) = c("Year", "Age", paste("Area_", area_ind, sep = ""))
  # 
  # write.table(region_comb_v2, file = paste("R_", region_index, "_COB_", cob_ind, "_region_comb.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)  
  # write.table(region_comb_pop_v2, file =  paste("R_", region_index, "_COB_", cob_ind, "_region_comb_pop.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  region_comb_demogdata = read.demogdata(paste("R_", region_index, "_COB_", cob_ind, "_region_comb.txt", sep = ""), paste("R_", region_index, "_COB_", cob_ind, "_region_comb_pop.txt", sep = ""), type = "fertility", label = paste("R_", region_index, sep = ""), skip = 0)
  
  res_region = array(NA, dim = c(year_horizon, 31, year_horizon, length(area_ind)))
  train_residual = list()
  for(ij in 1:year_horizon)
  {
    ind_dat = extract.years(region_comb_demogdata, years = 1981:(n_year+ij))
    # pcamethod = match.arg(pcamethod)
    fun_forc = mfts(ind_dat, pcamethod = pcamethod, year_horizon = year_horizon)
    
    train_residual[[ij]] = array(NA, dim = c(31, (20+ij), length(area_ind)))
    for(iwk in 1:length(area_ind))
    {
      res_region[,,ij,iwk] = t(fun_forc$fore_res[(31*(iwk-1)+1):(31*iwk),])
      
      train_residual[[ij]][,,iwk] = fun_forc$train_residual[(31*(iwk-1)+1):(31*iwk),]
    }
  }
  
  # Errors
  
  region_mae = region_mase = region_rmse = matrix(NA, year_horizon, length(area_ind))
  age_ind = match(age_all, age_interp) 

  for(iw in 1:length(area_ind))
  {
    for(ik in 1:year_horizon)
    {
      region_mae[ik,iw] = ftsa:::mae(res_region[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_area_fix_cob[[area_ind[iw]]][cob_ind]), years = (2001+ik):2011)$rate$female)
      region_mase[ik,iw] = ftsa:::mase(res_region[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_area_fix_cob[[area_ind[iw]]][cob_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_area_fix_cob[[area_ind[iw]]][cob_ind]), years = 1981:(2001+ik-1))$rate$female)
      region_rmse[ik,iw] = ftsa:::rmse(res_region[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_area_fix_cob[[area_ind[iw]]][cob_ind]), years = (2001+ik):2011)$rate$female)
      
    }
  }
  return(list(res_region = res_region, train_residual = train_residual, region_mae = region_mae, region_mase = region_mase,
              region_rmse = region_rmse))    
  
}

test = mfts_area_fix_COB(region_index = 2,  pcamethod = c("dynamic"), year_horizon = 10, cob_ind = 1)

region_all_forc = list()
for(r_ind in 1:11)
{
  region_all_forc[[r_ind]] = list()
  for(cob in 1:19)
  {
    region_all_forc[[r_ind]][[cob]] = mfts_area_fix_COB(region_index = r_ind,  pcamethod = c("dynamic"), year_horizon = 10, cob_ind = cob)
  }
  print(r_ind)
}


for(r_ind in 1:11)
{
  area_ind_all = region_list_ind[[r_ind]]
  if(r_ind == 11)
  {
    area_ind_all = region_list_ind[[r_ind]][-3]
  }
  
  
  for(ik in 1:length(area_ind_all))
  {
    temp_res_forc_list = temp_train_residual_list = temp_mae_list = temp_mase_list = temp_rmse_list = list()
    for(cob in 1:19)
    {
      temp_res_forc_list[[cob]] = region_all_forc[[r_ind]][[cob]]$res_region[,,,ik]
      temp_train_residual_list[[cob]] = lapply(region_all_forc[[r_ind]][[cob]]$train_residual, "[", ,,ik)
      temp_mae_list[[cob]] = region_all_forc[[r_ind]][[cob]]$region_mae[,ik]
      temp_mase_list[[cob]] = region_all_forc[[r_ind]][[cob]]$region_mase[,ik]
      temp_rmse_list[[cob]] = region_all_forc[[r_ind]][[cob]]$region_rmse[,ik]
    }
    
    assign(area_forecast_fix_COB[area_ind_all[ik]], temp_res_forc_list)
    assign(area_train_residual_fix_COB[area_ind_all[ik]], temp_train_residual_list)
    assign(area_mae_fix_COB[area_ind_all[ik]], temp_mae_list)
    assign(area_mase_fix_COB[area_ind_all[ik]], temp_mase_list)
    assign(area_rmse_fix_COB[area_ind_all[ik]], temp_rmse_list)
  
    rm(temp_res_forc_list, temp_train_residual_list, temp_mae_list, temp_mase_list, temp_rmse_list)
  }
}

fix_cob_all_area_mae = array(0, dim = c(10, 47, 19))
for(ik in 1:47)
{
  temp = get(area_mae_fix_COB[ik])
  for(ij in 1:19)
  {
    if(ik == 23 && ij == 14)
    {
      fix_cob_all_area_mae[,ik,ij] = rep(0, 10)
    } else {
      fix_cob_all_area_mae[,ik,ij] = temp[[ij]]
    }
  }
}


fix_cob_all_area_mase = array(0, dim = c(10, 47, 19))
for(ik in 1:47)
{
  for(ij in 1:19)
  {
    if(ik == 23)
    {
      fix_cob_all_area_mase[,ik,ij] = rep(0, 10)
    } else {
      temp = get(area_mase_fix_COB[ik])
      fix_cob_all_area_mase[,ik,ij] = temp[[ij]]
    }
  }
}



fix_cob_all_area_rmse = array(0, dim = c(10, 47, 19))
for(ik in 1:47)
{
  temp = get(area_rmse_fix_COB[ik])
  for(ij in 1:19)
  {
    if(ik == 23 && ij == 14)
    {
      fix_cob_all_area_rmse[,ik,ij] = rep(0, 10)
    } else {
      fix_cob_all_area_rmse[,ik,ij] = temp[[ij]]
    }
  }
}


# Collectively forecast all 11 regions 

mfts_region_fix_COB <- function(pcamethod = c("static", "dynamic"), year_horizon, cob_ind)
{
  require(ftsa)
  require(demography)
  n_year = 2011 - (year_horizon + 1)
  
  # national_comb = national_comb_pop = matrix(NA, 31*31, 11)
  # for(iw in 1:11)
  # {
  #   national_comb[,iw] = as.numeric(get(fts_region_fix_cob_smooth[[iw]][cob_ind])$rate$female[3:33,])
  #   national_comb_pop[,iw] = as.numeric(get(fts_region_fix_cob_smooth[[iw]][cob_ind])$pop$female[3:33,])
  # }
  # national_comb_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), national_comb)
  # national_comb_pop_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), national_comb_pop)
  # 
  # colnames(national_comb_v2) = colnames(national_comb_pop_v2) = c("Year", "Age", paste("Region_", 1:11, sep = ""))
  # 
  # write.table(national_comb_v2, file = paste("national_COB_", cob_ind, "_comb.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)  
  # write.table(national_comb_pop_v2, file =  paste("national_COB_", cob_ind, "_comb_pop.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  national_comb_demogdata = read.demogdata(paste("national_COB_", cob_ind, "_comb.txt", sep = ""), paste("national_COB_", cob_ind, "_comb_pop.txt", sep = ""), type = "fertility", label = "national", skip = 0)
  
  res_national = array(NA, dim = c(year_horizon, 31, year_horizon, 11))
  train_residual = list()
  for(ij in 1:year_horizon)
  {
    ind_dat = extract.years(national_comb_demogdata, years = 1981:(n_year+ij))
    # pcamethod = match.arg(pcamethod)
    fun_forc = mfts(ind_dat, pcamethod = pcamethod, year_horizon = year_horizon)
    
    train_residual[[ij]] = array(NA, dim = c(31, (20+ij), 11))
    for(iwk in 1:11)
    {
      res_national[,,ij,iwk] = t(fun_forc$fore_res[(31*(iwk-1)+1):(31*iwk),])
      
      train_residual[[ij]][,,iwk] = fun_forc$train_residual[(31*(iwk-1)+1):(31*iwk),]
    }
  }
  
  # Errors
  
  national_mae = national_mase = national_rmse = matrix(NA, year_horizon, 11)
  age_ind = match(age_all, age_interp) 
  
  for(iw in 1:11)
  {
    for(ik in 1:year_horizon)
    {
      national_mae[ik,iw] = ftsa:::mae(res_national[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_region_fix_cob[[iw]][cob_ind]), years = (2001+ik):2011)$rate$female)
      national_mase[ik,iw] = ftsa:::mase(res_national[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_region_fix_cob[[iw]][cob_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_region_fix_cob[[iw]][cob_ind]), years = 1981:(2001+ik-1))$rate$female)
      national_rmse[ik,iw] = ftsa:::rmse(res_national[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_region_fix_cob[[iw]][cob_ind]), years = (2001+ik):2011)$rate$female)
    }
  }
  return(list(res_national = res_national, train_residual = train_residual, national_mae = national_mae, national_mase = national_mase, national_rmse = national_rmse))  
}

region_all_forc_fix_cob = list()
for(cob in 1:19)
{
  region_all_forc_fix_cob[[cob]] = mfts_region_fix_COB(pcamethod = c("dynamic"), year_horizon = 10, cob_ind = cob)
  print(cob)
}


for(ik in 1:11)
{
  temp_res_forc_list = temp_train_residual_list = temp_mae_list = temp_mase_list = temp_rmse_list = list()
  for(cob in 1:19)
  {
    temp_res_forc_list[[cob]] = region_all_forc_fix_cob[[cob]]$res_national[,,,ik]
    temp_train_residual_list[[cob]] = lapply(region_all_forc_fix_cob[[cob]]$train_residual, "[", ,,ik)
    temp_mae_list[[cob]] = region_all_forc_fix_cob[[cob]]$national_mae[,ik]
    temp_mase_list[[cob]] = region_all_forc_fix_cob[[cob]]$national_mase[,ik]
    temp_rmse_list[[cob]] = region_all_forc_fix_cob[[cob]]$national_rmse[,ik]
  }
  
  assign(region_forecast_fix_COB[ik], temp_res_forc_list)
  assign(region_train_residual_fix_COB[ik], temp_train_residual_list)
  assign(region_mae_fix_COB[ik], temp_mae_list)
  assign(region_mase_fix_COB[ik], temp_mase_list)
  assign(region_rmse_fix_COB[ik], temp_rmse_list)
  
  rm(temp_res_forc_list, temp_train_residual_list, temp_mae_list, temp_mase_list, temp_rmse_list)
}

fix_cob_all_region_mae = array(0, dim = c(10, 11, 19))
for(ik in 1:11)
{
  temp = get(region_mae_fix_COB[ik])
  for(ij in 1:19)
  {
    fix_cob_all_region_mae[,ik,ij] = temp[[ij]]
  }
}

fix_cob_all_region_mase = array(0, dim = c(10, 11, 19))
for(ik in 1:11)
{
  temp = get(region_mase_fix_COB[ik])
  for(ij in 1:19)
  {
    fix_cob_all_region_mase[,ik,ij] = temp[[ij]]
  }
}

fix_cob_all_region_rmse = array(0, dim = c(10, 11, 19))
for(ik in 1:11)
{
  temp = get(region_rmse_fix_COB[ik])
  for(ij in 1:19)
  {
    fix_cob_all_region_rmse[,ik,ij] = temp[[ij]]
  }
}

# National series by itself

# mfts_national_fix_COB <- function(pcamethod = c("static", "dynamic"), year_horizon, cob_ind)
# {
#   require(ftsa)
#   require(demography)
#   n_year = 2011 - (year_horizon + 1)
#   
#   national_comb = as.numeric(get(fts_national_fix_cob_smooth[cob_ind])$rate$female[3:33,])
#   national_comb_pop = as.numeric(get(fts_national_fix_cob_smooth[cob_ind])$pop$female[3:33,])
#   
#   national_comb_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), national_comb)
#   national_comb_pop_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), national_comb_pop)
#   
#   colnames(national_comb_v2) = colnames(national_comb_pop_v2) = c("Year", "Age", "National")
#   
#   write.table(national_comb_v2, file = paste("AUS_COB_", cob_ind, "_comb.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)  
#   write.table(national_comb_pop_v2, file =  paste("AUS_COB_", cob_ind, "_comb_pop.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
#   
#   national_comb_demogdata = read.demogdata(paste("AUS_COB_", cob_ind, "_comb.txt", sep = ""), paste("AUS_COB_", cob_ind, "_comb_pop.txt", sep = ""), type = "fertility", label = "national", skip = 0)
#   
#   res_national = array(NA, dim = c(year_horizon, 31, year_horizon))
#   train_residual = list()
#   for(ij in 1:year_horizon)
#   {
#     ind_dat = extract.years(national_comb_demogdata, years = 1981:(n_year+ij))
#     # pcamethod = match.arg(pcamethod)
#     fun_forc = mfts(ind_dat, pcamethod = pcamethod, year_horizon = year_horizon)
#     
#     res_national[,,ij] = t(fun_forc$fore_res)
#     
#     train_residual[[ij]] = fun_forc$train_residual
#   }
#   
#   # Errors
#   
#   national_mae = national_rmse = matrix(NA, year_horizon, 1)
#   age_ind = match(age_all, age_interp) 
#   
#   for(ik in 1:year_horizon)
#   {
#     national_mae[ik,1] = ftsa:::mae(res_national[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_national_fix_cob[cob_ind]), years = (2001+ik):2011)$rate$female)
#     national_rmse[ik,1] = ftsa:::rmse(res_national[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_national_fix_cob[cob_ind]), years = (2001+ik):2011)$rate$female)
#   }
#   
#   return(list(res_national = res_national, train_residual = train_residual, 
#               national_mae = national_mae, national_rmse = national_rmse))    
# }
# 
# 
# national_forc_fix_cob = list()
# for(cob in 1:19)
# {
#   national_forc_fix_cob[[cob]] = mfts_national_fix_COB(pcamethod = c("dynamic"), year_horizon = 10, cob_ind = cob)
# }
# 
# national_forecast_list = national_train_residual_list = national_mae_list = national_rmse_list = list()
# 
# for(cob in 1:19)
# {
#   national_forecast_list[[cob]] = national_forc_fix_cob[[cob]]$res_national
#   national_train_residual_list[[cob]] = national_forc_fix_cob[[cob]]$train_residual
#   national_mae_list[[cob]] = national_forc_fix_cob[[cob]]$national_mae
#   national_rmse_list[[cob]] = national_forc_fix_cob[[cob]]$national_rmse
# }
# 
# assign(national_forecast_fix_COB[1], national_forecast_list)
# assign(national_train_residual_fix_COB[1], national_train_residual_list)
# assign(national_mae_fix_COB[1], national_mae_list)
# assign(national_rmse_fix_COB[1], national_rmse_list)

fix_cob_all_national_mae = matrix(0, nrow = 10, ncol = 19)
for(ij in 1:19)
{
  fix_cob_all_national_mae[,ij] = get(ind_national_mae_fix_COB[1])[[ij]]
}

fix_cob_all_national_mase = matrix(0, nrow = 10, ncol = 19)
for(ij in 1:19)
{
  fix_cob_all_national_mase[,ij] = get(ind_national_mase_fix_COB[1])[[ij]]
}


fix_cob_all_national_rmse = matrix(0, nrow = 10, ncol = 19)
for(ij in 1:19)
{
  fix_cob_all_national_rmse[,ij] = get(ind_national_rmse_fix_COB[1])[[ij]]
}

# combine all point forecasts into an array

fix_cob_all_level_mae = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  fix_cob_all_level_mae[,,i] = cbind(fix_cob_all_national_mae[,i], rowMeans(fix_cob_all_region_mae[,,i]), rowMeans(fix_cob_all_area_mae[,,i]))
}
fix_cob_averaged_mae = apply(fix_cob_all_level_mae, c(1,2), mean)

###
fix_cob_all_area_mase_copy = fix_cob_all_area_mase
for(i in 1:19)
{
  fix_cob_all_area_mase[,,i][,which(colMeans(fix_cob_all_area_mase[,,i]) > 300)] = 0
}
###

fix_cob_all_level_mase = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  fix_cob_all_level_mase[,,i] = cbind(fix_cob_all_national_mase[,i], rowMeans(fix_cob_all_region_mase[,,i]), rowMeans(fix_cob_all_area_mase[,,i]))
}
fix_cob_averaged_mase = apply(fix_cob_all_level_mase, c(1,2), mean)

fix_cob_all_level_rmse = array(0, dim = c(10,3,19))
for(i in 1:19)
{
  fix_cob_all_level_rmse[,,i] = cbind(fix_cob_all_national_rmse[,i], rowMeans(fix_cob_all_region_rmse[,,i]), rowMeans(fix_cob_all_area_rmse[,,i]))
}
fix_cob_averaged_rmse = apply(fix_cob_all_level_rmse, c(1,2), mean)


###################################################
# Fix GEO; disaggregate by Country/Region of birth
###################################################

# COB level 

cob_train_residual_fix_GEO = vector("character", length = 19)
for(i in 1:19)
{
  cob_train_residual_fix_GEO[i] = paste("COB", i, "_train_residual_fix_GEO", sep = "")
}

cob_forecast_fix_GEO = vector("character", length = 19)
for(i in 1:19)
{
  cob_forecast_fix_GEO[i] = paste("COB", i, "_forecast_fix_GEO", sep = "")
}

cob_mae_fix_GEO = vector("character", length = 19)
for(i in 1:19)
{
  cob_mae_fix_GEO[i] = paste("COB", i, "_mae_fix_GEO", sep = "")
}

cob_mase_fix_GEO = vector("character", length = 19)
for(i in 1:19)
{
  cob_mase_fix_GEO[i] = paste("COB", i, "_mase_fix_GEO", sep = "")
}

cob_rmse_fix_GEO = vector("character", length = 19)
for(i in 1:19)
{
  cob_rmse_fix_GEO[i] = paste("COB", i, "_rmse_fix_GEO", sep = "")
}

# C_group level

C_group_train_residual_fix_GEO = vector("character", length = 10)
for(i in 1:10)
{
  C_group_train_residual_fix_GEO[i] = paste("C", i, "_train_residual_fix_GEO", sep = "")
}

C_group_forecast_fix_GEO = vector("character", length = 10)
for(i in 1:10)
{
  C_group_forecast_fix_GEO[i] = paste("C", i, "_forecast_fix_GEO", sep = "")
}

C_group_mae_fix_GEO = vector("character", length = 10)
for(i in 1:10)
{
  C_group_mae_fix_GEO[i] = paste("C", i, "_mae_fix_GEO", sep = "")
}

C_group_mase_fix_GEO = vector("character", length = 10)
for(i in 1:10)
{
  C_group_mase_fix_GEO[i] = paste("C", i, "_mase_fix_GEO", sep = "")
}

C_group_rmse_fix_GEO = vector("character", length = 10)
for(i in 1:10)
{
  C_group_rmse_fix_GEO[i] = paste("C", i, "_rmse_fix_GEO", sep = "")
}

# global level

global_train_residual_fix_GEO = paste("Global", "_train_residual_fix_GEO", sep = "")
global_forecast_fix_GEO = paste("Global", "_forecast_fix_GEO", sep = "")
global_mae_fix_GEO = paste("Global", "_mae_fix_GEO", sep = "")
global_mase_fix_GEO = paste("Global", "_mase_fix_GEO", sep = "")
global_rmse_fix_GEO = paste("Global", "_rmse_fix_GEO", sep = "")


# Collectively forecast all 19 COB places

mfts_area_fix_GEO <- function(C_group = 1,  pcamethod = c("static", "dynamic"), year_horizon, geo_ind = 1)
{
  require(ftsa)
  require(demography)
  
  n_year = 2011 - (year_horizon + 1)
  
  cob_ind = COB_list[[C_group]]
  
  if(C_group == 6 && geo_ind == 23)
  {
    cob_ind = cob_ind[-2] # no observation at all for Area 23 COB 14
  }
  
  # C_comb = C_comb_pop = matrix(NA, 31*31, length(cob_ind))
  # for(iw in 1:length(cob_ind))
  # {
  #   C_comb[,iw] = as.numeric(get(fts_area_fix_geo_smooth[[cob_ind[iw]]][geo_ind])$rate$female[3:33,])
  #   C_comb_pop[,iw] = as.numeric(get(fts_area_fix_geo_smooth[[cob_ind[iw]]][geo_ind])$pop$female[3:33,])
  # }
  # 
  # C_comb_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), C_comb)
  # C_comb_pop_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), C_comb_pop)
  # 
  # colnames(C_comb_v2) = colnames(C_comb_pop_v2) = c("Year", "Age", paste("COB_", cob_ind, sep = ""))
  # 
  # write.table(C_comb_v2, file = paste("C_", C_group, "_GEO_", geo_ind, "_comb.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)  
  # write.table(C_comb_pop_v2, file =  paste("C_", C_group, "_GEO_", geo_ind, "_comb_pop.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
  # 
  C_comb_demogdata = read.demogdata(paste("C_", C_group, "_GEO_", geo_ind, "_comb.txt", sep = ""), paste("C_", C_group, "_GEO_", geo_ind, "_comb_pop.txt", sep = ""), type = "fertility", label = paste("C_", C_group, sep = ""), skip = 0)
  
  res_C = array(NA, dim = c(year_horizon, 31, year_horizon, length(cob_ind)))
  train_residual = list()
  for(ij in 1:year_horizon)
  {
    ind_dat = extract.years(C_comb_demogdata, years = 1981:(n_year+ij))
    # pcamethod = match.arg(pcamethod)
    fun_forc = mfts(ind_dat, pcamethod = pcamethod, year_horizon = year_horizon)
    
    train_residual[[ij]] = array(NA, dim = c(31, (20+ij), length(cob_ind)))
    for(iwk in 1:length(cob_ind))
    {
      res_C[,,ij,iwk] = t(fun_forc$fore_res[(31*(iwk-1)+1):(31*iwk),])
      
      train_residual[[ij]][,,iwk] = fun_forc$train_residual[(31*(iwk-1)+1):(31*iwk),]
    }
  }
  
  # Errors
  
  C_mae = C_mase = C_rmse = matrix(NA, year_horizon, length(cob_ind))
  age_ind = match(age_all, age_interp) 
  
  for(iw in 1:length(cob_ind))
  {
    for(ik in 1:year_horizon)
    {
      C_mae[ik,iw] = ftsa:::mae(res_C[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_area_fix_geo[[cob_ind[iw]]][geo_ind]), years = (2001+ik):2011)$rate$female)
      C_mase[ik,iw] = ftsa:::mase(res_C[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_area_fix_geo[[cob_ind[iw]]][geo_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_area_fix_geo[[cob_ind[iw]]][geo_ind]), years = 1986:(2001+ik-1))$rate$female)
      C_rmse[ik,iw] = ftsa:::rmse(res_C[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_area_fix_geo[[cob_ind[iw]]][geo_ind]), years = (2001+ik):2011)$rate$female)
      
    }
  }
  return(list(res_C = res_C, train_residual = train_residual, C_mae = C_mae, C_mase = C_mase, C_rmse = C_rmse))    
  
}

test = mfts_area_fix_GEO(C_group = 2,  pcamethod = c("dynamic"), year_horizon = 10, geo_ind = 59)

C_group_all_forc = list()
for(C_ind in 1:10)
{
  C_group_all_forc[[C_ind]] = list()
  for(geo in 1:59)
  {
    C_group_all_forc[[C_ind]][[geo]] = mfts_area_fix_GEO(C_group = C_ind,  pcamethod = c("dynamic"), year_horizon = 10, geo_ind = geo)
  }
  print(C_ind)
}


for(C_ind in 1:10)
{
  cob_ind_all = COB_list[[C_ind]]
  for(ik in 1:length(cob_ind_all))
  {
    temp_res_forc_list = temp_train_residual_list = temp_mae_list = temp_mase_list = temp_rmse_list = list()
    if(cob_ind_all[ik] == 14)
    {
      no_23 = c(1:22,24:59)
      for(geo in no_23)
      {
        temp_res_forc_list[[geo]] = C_group_all_forc[[C_ind]][[geo]]$res_C[,,,ik]
        temp_train_residual_list[[geo]] = lapply(C_group_all_forc[[C_ind]][[geo]]$train_residual, "[", ,,ik)
        temp_mae_list[[geo]] = C_group_all_forc[[C_ind]][[geo]]$C_mae[,ik]
        temp_mase_list[[geo]] = C_group_all_forc[[C_ind]][[geo]]$C_mase[,ik]
        temp_rmse_list[[geo]] = C_group_all_forc[[C_ind]][[geo]]$C_rmse[,ik]
      }
      
      temp_res_forc_list[[23]] = temp_train_residual_list[[23]] = temp_mae_list[[23]] = temp_mase_list[[23]] = temp_rmse_list[[23]] = NA
    } else {
      for(geo in 1:59)
      {
        temp_res_forc_list[[geo]] = C_group_all_forc[[C_ind]][[geo]]$res_C[,,,ik]
        temp_train_residual_list[[geo]] = lapply(C_group_all_forc[[C_ind]][[geo]]$train_residual, "[", ,,ik)
        temp_mae_list[[geo]] = C_group_all_forc[[C_ind]][[geo]]$C_mae[,ik]
        temp_mase_list[[geo]] = C_group_all_forc[[C_ind]][[geo]]$C_mase[,ik]
        temp_rmse_list[[geo]] = C_group_all_forc[[C_ind]][[geo]]$C_rmse[,ik]
      }
    }
    
    assign(cob_forecast_fix_GEO[cob_ind_all[ik]], temp_res_forc_list)
    assign(cob_train_residual_fix_GEO[cob_ind_all[ik]], temp_train_residual_list)
    assign(cob_mae_fix_GEO[cob_ind_all[ik]], temp_mae_list)
    assign(cob_mase_fix_GEO[cob_ind_all[ik]], temp_mase_list)
    assign(cob_rmse_fix_GEO[cob_ind_all[ik]], temp_rmse_list)
    
    rm(temp_res_forc_list, temp_train_residual_list, temp_mae_list, temp_mase_list, temp_rmse_list)
  }
}

# summary of results

fix_geo_all_cob_mae = array(0, dim = c(10, 19, 59))
for(ik in 1:19)
{
  temp = get(cob_mae_fix_GEO[ik])
  for(ij in 1:59)
  {
    if(ik == 14 && ij == 23)
    {
      fix_geo_all_cob_mae[,ik,ij] = rep(0, 10)
    } else {
      fix_geo_all_cob_mae[,ik,ij] = temp[[ij]]
    }
  }
}

fix_geo_all_cob_mase = array(0, dim = c(10, 19, 59))
for(ik in 1:19)
{
  temp = get(cob_mase_fix_GEO[ik])
  for(ij in 1:59)
  {
    if(ik == 14 && ij == 23)
    {
      fix_geo_all_cob_mase[,ik,ij] = rep(0, 10)
    } else {
      fix_geo_all_cob_mase[,ik,ij] = temp[[ij]]
    }
  }
}

fix_geo_all_cob_rmse = array(0, dim = c(10, 19, 59))
for(ik in 1:19)
{
  temp = get(cob_rmse_fix_GEO[ik])
  for(ij in 1:59)
  {
    if(ik == 14 && ij == 23)
    {
      fix_geo_all_cob_rmse[,ik,ij] = rep(0, 10)
    } else {
      fix_geo_all_cob_rmse[,ik,ij] = temp[[ij]]
    }
  }
}


# Collectively forecast 10 C_groups

mfts_C_fix_GEO <- function(pcamethod = c("static", "dynamic"), year_horizon, geo_ind = 1)
{
  require(ftsa)
  require(demography)
  
  n_year = 2011 - (year_horizon + 1)
  
  if(geo_ind == 23)
  {
    # national_comb = national_comb_pop = matrix(NA, 31*31, 9)
    # c_sub = c(1:5,7:10)
    # for(iw in 1:9)
    # {
    #   national_comb[,iw] = as.numeric(get(fts_region_fix_geo_smooth[[c_sub[iw]]][geo_ind])$rate$female[3:33,])
    #   national_comb_pop[,iw] = as.numeric(get(fts_region_fix_geo_smooth[[c_sub[iw]]][geo_ind])$pop$female[3:33,])
    # }
    # national_comb_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), national_comb)
    # national_comb_pop_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), national_comb_pop)
    # 
    # colnames(national_comb_v2) = colnames(national_comb_pop_v2) = c("Year", "Age", paste("C_", c_sub, sep = ""))
    # 
    # write.table(national_comb_v2, file = paste("C_total_GEO_", geo_ind, "_comb.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)  
    # write.table(national_comb_pop_v2, file =  paste("C_total_GEO_", geo_ind, "_comb_pop.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
    
    national_comb_demogdata = read.demogdata(paste("C_total_GEO_", geo_ind, "_comb.txt", sep = ""), paste("C_total_GEO_", geo_ind, "_comb_pop.txt", sep = ""), type = "fertility", label = "C_total", skip = 0)
    
    res_national = array(NA, dim = c(year_horizon, 31, year_horizon, 9))
    train_residual = list()
    for(ij in 1:year_horizon)
    {
      ind_dat = extract.years(national_comb_demogdata, years = 1981:(n_year+ij))
      # pcamethod = match.arg(pcamethod)
      fun_forc = mfts(ind_dat, pcamethod = pcamethod, year_horizon = year_horizon)
      
      train_residual[[ij]] = array(NA, dim = c(31, (20+ij), 9))
      for(iwk in 1:9)
      {
        res_national[,,ij,iwk] = t(fun_forc$fore_res[(31*(iwk-1)+1):(31*iwk),])
        
        train_residual[[ij]][,,iwk] = fun_forc$train_residual[(31*(iwk-1)+1):(31*iwk),]
      }
    }
    
    # Errors
    
    national_mae = national_mase = national_rmse = matrix(NA, year_horizon, 9)
    age_ind = match(age_all, age_interp) 
    
    for(iw in 1:9)
    {
      for(ik in 1:year_horizon)
      {
        national_mae[ik,iw] = ftsa:::mae(res_national[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_region_fix_geo[[c_sub[iw]]][geo_ind]), years = (2001+ik):2011)$rate$female)
        national_mase[ik,iw] = ftsa:::mase(res_national[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_region_fix_geo[[c_sub[iw]]][geo_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_region_fix_geo[[c_sub[iw]]][geo_ind]), years = 1986:(2001+ik-1))$rate$female)
        national_rmse[ik,iw] = ftsa:::rmse(res_national[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_region_fix_geo[[c_sub[iw]]][geo_ind]), years = (2001+ik):2011)$rate$female)
      }
    }
    
  } else {
    # national_comb = national_comb_pop = matrix(NA, 31*31, 10)
    # for(iw in 1:10)
    # {
    #   national_comb[,iw] = as.numeric(get(fts_region_fix_geo_smooth[[iw]][geo_ind])$rate$female[3:33,])
    #   national_comb_pop[,iw] = as.numeric(get(fts_region_fix_geo_smooth[[iw]][geo_ind])$pop$female[3:33,])
    # }
    # national_comb_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), national_comb)
    # national_comb_pop_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), national_comb_pop)
    # 
    # colnames(national_comb_v2) = colnames(national_comb_pop_v2) = c("Year", "Age", paste("C_", 1:10, sep = ""))
    # 
    # write.table(national_comb_v2, file = paste("C_total_GEO_", geo_ind, "_comb.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)  
    # write.table(national_comb_pop_v2, file =  paste("C_total_GEO_", geo_ind, "_comb_pop.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
    
    national_comb_demogdata = read.demogdata(paste("C_total_GEO_", geo_ind, "_comb.txt", sep = ""), paste("C_total_GEO_", geo_ind, "_comb_pop.txt", sep = ""), type = "fertility", label = "C_total", skip = 0)
    
    res_national = array(NA, dim = c(year_horizon, 31, year_horizon, 10))
    train_residual = list()
    for(ij in 1:year_horizon)
    {
      ind_dat = extract.years(national_comb_demogdata, years = 1981:(n_year+ij))
      # pcamethod = match.arg(pcamethod)
      fun_forc = mfts(ind_dat, pcamethod = pcamethod, year_horizon = year_horizon)
      
      train_residual[[ij]] = array(NA, dim = c(31, (20+ij), 10))
      for(iwk in 1:10)
      {
        res_national[,,ij,iwk] = t(fun_forc$fore_res[(31*(iwk-1)+1):(31*iwk),])
        
        train_residual[[ij]][,,iwk] = fun_forc$train_residual[(31*(iwk-1)+1):(31*iwk),]
      }
    }
    
    # Errors
    
    national_mae = national_mase = national_rmse = matrix(NA, year_horizon, 10)
    age_ind = match(age_all, age_interp) 
    
    for(iw in 1:10)
    {
      for(ik in 1:year_horizon)
      {
        national_mae[ik,iw] = ftsa:::mae(res_national[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_region_fix_geo[[iw]][geo_ind]), years = (2001+ik):2011)$rate$female)
        national_mase[ik,iw] = ftsa:::mase(res_national[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_region_fix_geo[[iw]][geo_ind]), years = (2001+ik):2011)$rate$female, insampletrue = extract.years(get(fts_region_fix_geo[[iw]][geo_ind]), years = 1986:(2001+ik-1))$rate$female)
        national_rmse[ik,iw] = ftsa:::rmse(res_national[ik,age_ind,1:((year_horizon + 1)-ik),iw], extract.years(get(fts_region_fix_geo[[iw]][geo_ind]), years = (2001+ik):2011)$rate$female)
      }
    }
  }

  return(list(res_C_total = res_national, train_residual = train_residual, C_total_mae = national_mae,  C_total_mase = national_mase, C_total_rmse = national_rmse))    
}


All_COB_forc = list()
for(geo in 1:59)
{
  All_COB_forc[[geo]] = mfts_C_fix_GEO(pcamethod = c("dynamic"), year_horizon = 10, geo_ind = geo)
  
  print(geo)
}


for(ik in 1:10)
{
  temp_res_forc_list = temp_train_residual_list = temp_mae_list = temp_mase_list = temp_rmse_list = list()
  for(geo in 1:59)
  {
    if(geo == 23)
    {
      if(ik %in% 1:6)
      {
        temp_res_forc_list[[geo]] = All_COB_forc[[geo]]$res_C_total[,,,ik]
        temp_train_residual_list[[geo]] = lapply(All_COB_forc[[geo]]$train_residual, "[", ,,ik)
        temp_mae_list[[geo]] = All_COB_forc[[geo]]$C_total_mae[,ik]
        temp_mase_list[[geo]] = All_COB_forc[[geo]]$C_total_mase[,ik]
        temp_rmse_list[[geo]] = All_COB_forc[[geo]]$C_total_rmse[,ik]
      }
      
      if(ik %in% 7:10)
      {
        temp_res_forc_list[[geo]] = All_COB_forc[[geo]]$res_C_total[,,,ik-1]
        temp_train_residual_list[[geo]] = lapply(All_COB_forc[[geo]]$train_residual, "[", ,,ik-1)
        temp_mae_list[[geo]] = All_COB_forc[[geo]]$C_total_mae[,ik-1]
        temp_mase_list[[geo]] = All_COB_forc[[geo]]$C_total_mase[,ik-1]
        temp_rmse_list[[geo]] = All_COB_forc[[geo]]$C_total_rmse[,ik-1]
      }
      
      if(ik == 6)
      {
        temp_res_forc_list[[geo]] = temp_train_residual_list[[geo]] = temp_mae_list[[geo]] = temp_mase_list[[geo]] = temp_rmse_list[[geo]] = NA
      }
      
    } else {
      temp_res_forc_list[[geo]] = All_COB_forc[[geo]]$res_C_total[,,,ik]
      temp_train_residual_list[[geo]] = lapply(All_COB_forc[[geo]]$train_residual, "[", ,,ik)
      temp_mae_list[[geo]] = All_COB_forc[[geo]]$C_total_mae[,ik]
      temp_mase_list[[geo]] = All_COB_forc[[geo]]$C_total_mase[,ik]
      temp_rmse_list[[geo]] = All_COB_forc[[geo]]$C_total_rmse[,ik]
    }
  }
  
  assign(C_group_forecast_fix_GEO[ik], temp_res_forc_list)
  assign(C_group_train_residual_fix_GEO[ik], temp_train_residual_list)
  assign(C_group_mae_fix_GEO[ik], temp_mae_list)
  assign(C_group_mase_fix_GEO[ik], temp_mase_list)
  assign(C_group_rmse_fix_GEO[ik], temp_rmse_list)
  
  rm(temp_res_forc_list, temp_train_residual_list, temp_mae_list, temp_mase_list, temp_rmse_list)
}

# summary of results

fix_geo_all_c_group_mae = array(0, dim = c(10, 10, 59))
for(ik in 1:10)
{
  temp = get(C_group_mae_fix_GEO[ik])
  for(ij in 1:59)
  {
    if(ik == 6 && ij == 23)
    {
      fix_geo_all_c_group_mae[,ik,ij] = rep(0, 10)
    } else {
      fix_geo_all_c_group_mae[,ik,ij] = temp[[ij]]
    }
  }
}

fix_geo_all_c_group_mase = array(0, dim = c(10, 10, 59))
for(ik in 1:10)
{
  temp = get(C_group_mase_fix_GEO[ik])
  for(ij in 1:59)
  {
    if(ik == 6 && ij == 23)
    {
      fix_geo_all_c_group_mase[,ik,ij] = rep(0, 10)
    } else {
      fix_geo_all_c_group_mase[,ik,ij] = temp[[ij]]
    }
  }
}

fix_geo_all_c_group_rmse = array(0, dim = c(10, 10, 59))
for(ik in 1:10)
{
  temp = get(C_group_rmse_fix_GEO[ik])
  for(ij in 1:59)
  {
    if(ik == 6 && ij == 23)
    {
      fix_geo_all_c_group_rmse[,ik,ij] = rep(0, 10)
    } else {
      fix_geo_all_c_group_rmse[,ik,ij] = temp[[ij]]
    }
  }
}

# global

# mfts_AUS_fix_GEO <- function(pcamethod = c("static", "dynamic"), year_horizon, geo_ind = 1)
# {
#   require(ftsa)
#   require(demography)
#   
#   n_year = 2011 - (year_horizon + 1)
#   
#   national_comb = as.numeric(get(fts_national_fix_geo_smooth[geo_ind])$rate$female[3:33,])
#   national_comb_pop = as.numeric(get(fts_national_fix_geo_smooth[geo_ind])$pop$female[3:33,])
#   
#   national_comb_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), national_comb)
#   national_comb_pop_v2 = cbind(rep(year_interp, each=31), rep(age_interp, 31), national_comb_pop)
#   
#   colnames(national_comb_v2) = colnames(national_comb_pop_v2) = c("Year", "Age", "National")
#   
#   write.table(national_comb_v2, file = paste("AUS_GEO_", geo_ind, "_comb.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)  
#   write.table(national_comb_pop_v2, file =  paste("AUS_GEO_", geo_ind, "_comb_pop.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
#   
#   national_comb_demogdata = read.demogdata(paste("AUS_GEO_", geo_ind, "_comb.txt", sep = ""), paste("AUS_GEO_", geo_ind, "_comb_pop.txt", sep = ""), type = "fertility", label = "AUS", skip = 0)
#   
#   res_national = array(NA, dim = c(year_horizon, 31, year_horizon))
#   train_residual = list()
#   for(ij in 1:year_horizon)
#   {
#     ind_dat = extract.years(national_comb_demogdata, years = 1981:(n_year+ij))
#     # pcamethod = match.arg(pcamethod)
#     fun_forc = mfts(ind_dat, pcamethod = pcamethod, year_horizon = year_horizon)
#     
#     res_national[,,ij] = t(fun_forc$fore_res)
#     
#     train_residual[[ij]] = fun_forc$train_residual
#   }
#   
#   # Errors
#   
#   national_mae = national_rmse = matrix(NA, year_horizon, 1)
#   age_ind = match(age_all, age_interp) 
#   
#   for(ik in 1:year_horizon)
#   {
#     national_mae[ik,1] = ftsa:::mae(res_national[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_national_fix_geo[geo_ind]), years = (2001+ik):2011)$rate$female)
#     national_rmse[ik,1] = ftsa:::rmse(res_national[ik,age_ind,1:((year_horizon + 1)-ik)], extract.years(get(fts_national_fix_geo[geo_ind]), years = (2001+ik):2011)$rate$female)
#   }
#   
#   return(list(res_C_total = res_national, train_residual = train_residual, C_total_mae = national_mae, C_total_rmse = national_rmse))    
# }
# 
# national_forc_fix_geo = list()
# for(geo in 1:59)
# {
#   national_forc_fix_geo[[geo]] = mfts_AUS_fix_GEO(pcamethod = c("dynamic"), year_horizon = 10, geo_ind = geo)
# }
# 
# national_forecast_list = national_train_residual_list = national_mae_list = national_rmse_list = list()
# for(geo in 1:59)
# {
#   national_forecast_list[[geo]] = national_forc_fix_geo[[geo]]$res_C_total
#   national_train_residual_list[[geo]] = national_forc_fix_geo[[geo]]$train_residual
#   national_mae_list[[geo]] = national_forc_fix_geo[[geo]]$C_total_mae
#   national_rmse_list[[geo]] = national_forc_fix_geo[[geo]]$C_total_rmse
# }
# 
# assign(global_forecast_fix_GEO[1], national_forecast_list)
# assign(global_train_residual_fix_GEO[1], national_train_residual_list)
# assign(global_mae_fix_GEO[1], national_mae_list)
# assign(global_rmse_fix_GEO[1], national_rmse_list)

fix_geo_all_global_mae = ind_fix_geo_all_global_mae
fix_geo_all_global_mase = ind_fix_geo_all_global_mase
fix_geo_all_global_rmse = ind_fix_geo_all_global_rmse


# combine all point forecasts into an array

fix_geo_all_level_mae = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  fix_geo_all_level_mae[,,i] = cbind(fix_geo_all_global_mae[,i], rowMeans(fix_geo_all_c_group_mae[,,i]), rowMeans(fix_geo_all_cob_mae[,,i]))
}
fix_geo_mae_averaged = apply(fix_geo_all_level_mae, c(1,2), mean)

###
fix_geo_all_c_group_mase_copy = fix_geo_all_c_group_mase
for(i in 1:59)
{
  fix_geo_all_c_group_mase[,,i][,which(colMeans(fix_geo_all_c_group_mase[,,i]) > 300)] = 0
}

fix_geo_all_cob_mase_copy = fix_geo_all_cob_mase
for(i in 1:59)
{
  fix_geo_all_cob_mase[,,i][,which(colMeans(fix_geo_all_cob_mase[,,i]) > 300)] = 0
}
###

fix_geo_all_level_mase = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  fix_geo_all_level_mase[,,i] = cbind(fix_geo_all_global_mase[,i], rowMeans(fix_geo_all_c_group_mase[,,i]), rowMeans(fix_geo_all_cob_mase[,,i]))
}
fix_geo_mase_averaged = apply(fix_geo_all_level_mase, c(1,2), mean)


fix_geo_all_level_rmse = array(0, dim = c(10,3,59))
for(i in 1:59)
{
  fix_geo_all_level_rmse[,,i] = cbind(fix_geo_all_global_rmse[,i], rowMeans(fix_geo_all_c_group_rmse[,,i]), rowMeans(fix_geo_all_cob_rmse[,,i]))
}
fix_geo_rmse_averaged = apply(fix_geo_all_level_rmse, c(1,2), mean)



