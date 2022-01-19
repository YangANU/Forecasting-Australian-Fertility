#####################################################
# Multivariate functional time series point forecast 
#####################################################

# Area level mfts training residuals 
mfts_Area_train_residual_female = paste("mfts", Area, "train_residual_female", sep = "_")
mfts_Area_train_residual_male   = paste("mfts", Area, "train_residual_male", sep = "_")
mfts_Area_train_residual_total  = paste("mfts", Area, "train_residual_total", sep = "_")

# Area level mfts forecasts
mfts_Area_forc_female = paste("mfts", Area, "forc_female", sep = "_")
mfts_Area_forc_male = paste("mfts", Area, "forc_male", sep = "_")
mfts_Area_forc_total = paste("mfts", Area, "forc_total", sep = "_")

# Mean error of Area level mfts forecasts
mfts_Area_me_female = paste("mfts", Area, "me_female", sep = "_")
mfts_Area_me_male = paste("mfts", Area, "me_male", sep = "_")
mfts_Area_me_total = paste("mfts", Area, "me_total", sep = "_")

# Mean absolute error of Area level mfts forecasts
mfts_Area_mae_female = paste("mfts", Area, "mae_female", sep = "_")
mfts_Area_mae_male = paste("mfts", Area, "mae_male", sep = "_")
mfts_Area_mae_total = paste("mfts", Area, "mae_total", sep = "_")

# Root mean square error of Area level mfts forecasts
mfts_Area_rmse_female = paste("mfts", Area, "rmse_female", sep = "_")
mfts_Area_rmse_male = paste("mfts", Area, "rmse_male", sep = "_")
mfts_Area_rmse_total = paste("mfts", Area, "rmse_total", sep = "_")

# Region level mfts training residuals 
mfts_Region_train_residual_female = paste("mfts", Region, "train_residual_female", sep = "_")
mfts_Region_train_residual_male   = paste("mfts", Region, "train_residual_male", sep = "_")
mfts_Region_train_residual_total  = paste("mfts", Region, "train_residual_total", sep = "_")

# Region level mfts forecasts
mfts_Region_forc_female = paste("mfts", Region, "forc_female", sep = "_")
mfts_Region_forc_male = paste("mfts", Region, "forc_male", sep = "_")
mfts_Region_forc_total = paste("mfts", Region, "forc_total", sep = "_")

# Mean error of Region level mfts forecasts
mfts_Region_me_female = paste("mfts", Region, "me_female", sep = "_")
mfts_Region_me_male = paste("mfts", Region, "me_male", sep = "_")
mfts_Region_me_total = paste("mfts", Region, "me_total", sep = "_")

# Mean absolute error of Region level mfts forecasts
mfts_Region_mae_female = paste("mfts", Region, "mae_female", sep = "_")
mfts_Region_mae_male = paste("mfts", Region, "mae_male", sep = "_")
mfts_Region_mae_total = paste("mfts", Region, "mae_total", sep = "_")

# Root mean square error of Region level mfts forecasts
mfts_Region_rmse_female = paste("mfts", Region, "rmse_female", sep = "_")
mfts_Region_rmse_male = paste("mfts", Region, "rmse_male", sep = "_")
mfts_Region_rmse_total = paste("mfts", Region, "rmse_total", sep = "_")


# multivarite functional time series model function

mfts <- function(dat, pcamethod = c("static", "dynamic"), year_horizon)
{
  n = length(names(dat$rate))
  # if the data contains female, male and total, then we consider only female and male
  if(n == 3)
  {
    n_pop = 2
  }else
  {
    n_pop = n
  }
  
  rowmeans_object = sd_object = decenter_object = list()
  for(ik in 1:n_pop)
  {
    # compute mean and standard deviation functions
    rowmeans_object[[ik]] = rowMeans(log(dat$rate[[ik]]), na.rm=TRUE)
    sd_object[[ik]] = apply(log(dat$rate[[ik]]), 1, sd, na.rm=TRUE)
    
    # de-center functional data
    decenter_object[[ik]] = t(scale(t(log(dat$rate[[ik]])), center = TRUE, scale = TRUE))
  }
  comb_object = do.call(rbind, decenter_object)
  
  pcamethod = match.arg(pcamethod)
  
  if (pcamethod == "static")
  {
    ncomp_comb = head(which(cumsum(ftsm(fts(1:nrow(comb_object), comb_object), order = 19)$varprop) >= 0.95), 1)
    fit_ftsm = ftsm(fts(1:nrow(comb_object), comb_object), order = ncomp_comb)
    train_residual = exp(comb_object * do.call(c,sd_object) + do.call(c, rowmeans_object)) - exp(fit_ftsm$fitted$y * do.call(c,sd_object) + do.call(c, rowmeans_object))
    fore_ftsm = forecast(fit_ftsm, h = year_horizon)
    fore_res = exp(fore_ftsm$mean$y * do.call(c,sd_object) + do.call(c, rowmeans_object))
  }
  
  
  if (pcamethod == "dynamic")
  {
    data_dum = comb_object
    
    C_0 = long_run_covariance_estimation(data_dum, H = 3, C0 = 3)
    eigen_decomp = eigen(C_0$BT_FT_fix_C0)
    dynamic_order = head(which(cumsum(eigen_decomp$values)/sum(eigen_decomp$values) >= 0.95),1)
    dynamic_basis = as.matrix(eigen_decomp$vectors[,1:dynamic_order])
    dynamic_scores = t(dynamic_basis) %*% data_dum
    
    train_residual = exp(comb_object * do.call(c,sd_object) + do.call(c, rowmeans_object)) - exp(dynamic_basis %*%dynamic_scores * do.call(c,sd_object) + do.call(c, rowmeans_object))
    
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
    
    fore_res = exp(fore_ftsm_dyn * do.call(c,sd_object) + do.call(c, rowmeans_object))
  }
  
  return(list(fore_res = fore_res, train_residual = train_residual))
}

# function for Area forecast

me   = ftsa:::me
mae  = ftsa:::mae
rmse = ftsa:::rmse

mfts_back_test_Area <- function(iw,  pcamethod = c("static", "dynamic"), fmethod = "classical", year_horizon)
{
  res_male = res_female = array(NA, dim = c(year_horizon, 86, year_horizon))
  train_residual_female = train_residual_male = list()
  
  pcamethod = match.arg(pcamethod)
  
  for(ij in 1:year_horizon)
  {
    ind_dat = extract.years(get(Area_smooth[iw]), years = 1993:(2010+ij))
    fun_forc = mfts(ind_dat, pcamethod = pcamethod, year_horizon = year_horizon)
    
    train_residual_female[[ij]] = fun_forc$train_residual[1:86,]
    train_residual_male[[ij]] =  fun_forc$train_residual[87:172,]
    
    res_female[,,ij] = t(fun_forc$fore_res[1:86,])
    res_male[,,ij]   = t(fun_forc$fore_res[87:172,])
    
  }
  # MAE & RMSE
  female_me = male_me = female_mae = male_mae = female_rmse = male_rmse = vector("numeric",year_horizon)
  for(ik in 1:year_horizon)
  {
    age_forc = 0:85
    age_obs  = get(Area[iw])$age
    age_ind  = age_forc %in% age_obs
    
    data_forc_female    = res_female[ik,age_ind,1:(6-ik)]
    data_testing_female = extract.years(get(Area[iw]), years = (2011+ik):2016)$rate$female
    testing_ind_female  = data_testing_female != 0
    
    data_forc_male    = res_male[ik,age_ind,1:(6-ik)]
    data_testing_male = extract.years(get(Area[iw]), years = (2011+ik):2016)$rate$male
    testing_ind_male  = data_testing_male != 0
    
    # me 
    female_me[ik] = me(data_forc_female[testing_ind_female], data_testing_female[testing_ind_female])
    male_me[ik]   = me(data_forc_male[testing_ind_male],   data_testing_male[testing_ind_male])
    
    # mae
    female_mae[ik] = mae(data_forc_female[testing_ind_female], data_testing_female[testing_ind_female])
    male_mae[ik]   = mae(data_forc_male[testing_ind_male],   data_testing_male[testing_ind_male])
    
    # rmse
    female_rmse[ik] = rmse(data_forc_female[testing_ind_female], data_testing_female[testing_ind_female])
    male_rmse[ik]   = rmse(data_forc_male[testing_ind_male],   data_testing_male[testing_ind_male])
    
  }
  return(list(res_female = res_female, res_male = res_male,
              train_residual_female = train_residual_female, 
              train_residual_male = train_residual_male,
              female_me = female_me, male_me = male_me, 
              female_mae = female_mae, male_mae = male_mae,
              female_rmse = female_rmse, male_rmse = male_rmse))
}

# function for Region forecast

mfts_back_test_Region <- function(iw, pcamethod = c("static", "dynamic"), year_horizon, fmethod = "classical")
{
  res_male = res_female = array(NA, dim = c(year_horizon, 86, year_horizon))
  
  train_residual_female = train_residual_male = list()
  
  pcamethod = match.arg(pcamethod)
  
  for(ij in 1:year_horizon)
  {
    ind_dat = extract.years(get(Region_smooth[iw]), years = 1993:(2010+ij))
    fun_forc = mfts(ind_dat, pcamethod = pcamethod, year_horizon = year_horizon)
    
    train_residual_female[[ij]] = fun_forc$train_residual[1:86,]
    train_residual_male[[ij]] =  fun_forc$train_residual[87:172,]
    
    res_female[,,ij] = t(fun_forc$fore_res[1:86,])
    res_male[,,ij]   = t(fun_forc$fore_res[87:172,])
  }
  # MAE
  female_me = male_me = female_mae = male_mae = female_rmse = male_rmse = vector("numeric", year_horizon)
  for(ik in 1:year_horizon)
  {
    age_forc = 0:85
    age_obs  = get(Region[iw])$age
    age_ind  = age_forc %in% age_obs
    
    data_forc_female    = res_female[ik,age_ind,1:(6-ik)]
    data_testing_female = extract.years(get(Region[iw]), years = (2011+ik):2016)$rate$female
    testing_ind_female  = data_testing_female != 0
    
    data_forc_male    = res_male[ik,age_ind,1:(6-ik)]
    data_testing_male = extract.years(get(Region[iw]), years = (2011+ik):2016)$rate$male
    testing_ind_male  = data_testing_male != 0
    
    female_me[ik] = me(data_forc_female[testing_ind_female], data_testing_female[testing_ind_female])
    male_me[ik]   = me(data_forc_male[testing_ind_male],     data_testing_male[testing_ind_male])
    
    female_mae[ik] = mae(data_forc_female[testing_ind_female], data_testing_female[testing_ind_female])
    male_mae[ik]   = mae(data_forc_male[testing_ind_male],     data_testing_male[testing_ind_male])
    
    female_rmse[ik] = rmse(data_forc_female[testing_ind_female], data_testing_female[testing_ind_female])
    male_rmse[ik]   = rmse(data_forc_male[testing_ind_male],     data_testing_male[testing_ind_male])
  }
  return(list(res_female = res_female, res_male = res_male,
              train_residual_female = train_residual_female, 
              train_residual_male = train_residual_male,
              female_me = female_me, male_me = male_me,
              female_mae = female_mae, male_mae = male_mae,
              female_rmse = female_rmse, male_rmse = male_rmse))
}


################################################################
# Multivariate forecast of age-specific mortality rates by Area
################################################################

for(ik in 1:48)
{
  dum = mfts_back_test_Area(iw = ik, pcamethod = c("static"), fmethod = "classical", year_horizon = 5)
  
  ## mortality rate
  assign(mfts_Area_forc_female[ik], dum$res_female)
  assign(mfts_Area_forc_male[ik],   dum$res_male)
  
  assign(mfts_Area_train_residual_female[ik], dum$train_residual_female)
  assign(mfts_Area_train_residual_male[ik], dum$train_residual_male)
  
  assign(mfts_Area_me_female[ik], dum$female_me)
  assign(mfts_Area_me_male[ik],   dum$male_me)
  
  assign(mfts_Area_mae_female[ik], dum$female_mae)
  assign(mfts_Area_mae_male[ik],   dum$male_mae)
  
  assign(mfts_Area_rmse_female[ik], dum$female_rmse)
  assign(mfts_Area_rmse_male[ik],   dum$male_rmse)
  
  rm(dum)
}

###################################################################
# Multivariate forecast of age-specific mortality rates by Regions
###################################################################

for(ik in 1:11)
{
  dum = mfts_back_test_Region(iw = ik, pcamethod = c("static"), year_horizon = 5)
  
  # mortality rate
  assign(mfts_Region_forc_female[ik], dum$res_female)
  assign(mfts_Region_forc_male[ik],   dum$res_male)
  
  assign(mfts_Region_train_residual_female[ik], dum$train_residual_female)
  assign(mfts_Region_train_residual_male[ik], dum$train_residual_male)
  
  assign(mfts_Region_me_female[ik], dum$female_me)
  assign(mfts_Region_me_male[ik],   dum$male_me)
  
  assign(mfts_Region_mae_female[ik], dum$female_mae)
  assign(mfts_Region_mae_male[ik],   dum$male_mae)
  
  assign(mfts_Region_rmse_female[ik], dum$female_rmse)
  assign(mfts_Region_rmse_male[ik],   dum$male_rmse)

  rm(dum)
}

###################################################################
# Multivariate forecast of Area Total age-specific mortality rates
###################################################################

mfts_back_test_Area_total <- function(fmethod = "classical", pcamethod = c("static", "dynamic"), year_horizon)
{
  total_comb = total_comb_pop =  matrix(NA, 86*24, 47)
  for(iw in 2:48)
  {
    total_comb[,iw-1] = as.numeric(get(Area_smooth[iw])$rate$total)
    total_comb_pop[,iw-1] = as.numeric(get(Area_smooth[iw])$pop$total)
  }

  total_comb_v2 = cbind(rep(1993:2016, each=86), rep(0:85, 24), total_comb)
  total_comb_pop_v2 = cbind(rep(1993:2016, each=86), rep(0:85, 24), total_comb_pop)
  colnames(total_comb_v2) = colnames(total_comb_pop_v2) = c("Year", "Age", Area[2:48])
  
  write.table(total_comb_v2, file = "total_comb_Area.txt", quote = FALSE, row.names = TRUE, col.names=TRUE)
  write.table(total_comb_pop_v2, file = "total_comb_pop_Area.txt", quote = FALSE, row.names = TRUE, col.names=TRUE)
  
  total_comb_demogdata = read.demogdata("total_comb_Area.txt", "total_comb_pop_Area.txt", 
                                        type = "mortality", label="total_comb_Area", skip = 0)
  
  res_Area = array(NA, dim = c(year_horizon, 86, year_horizon, 47))
  train_residual_total = list()
  for(ij in 1:year_horizon)
  {
    ind_dat = extract.years(total_comb_demogdata, years = 1993:(2010+ij))
    fun_forc = mfts(ind_dat, pcamethod = pcamethod, year_horizon = year_horizon)
    
    train_residual_total[[ij]] = array(NA, dim = c(86, (18+ij), 47))
    for(iwk in 1:47)
    {
      res_Area[,,ij,iwk] = t(fun_forc$fore_res[(86*(iwk-1)+1):(86*iwk),])
      train_residual_total[[ij]][,,iwk] = fun_forc$train_residual[(86*(iwk-1)+1):(86*iwk),]
    }
    
  }
  
  # MAE & RMSE
  total_Area_me = total_Area_mae = total_Area_rmse = matrix(NA,year_horizon,47)
  for(iw in 1:47)
  {
    for(ik in 1:year_horizon)
    {
      age_forc = 0:85
      age_obs  = get(Area[iw])$age
      age_ind  = age_forc %in% age_obs
      
      data_forc_total    = res_Area[ik,age_ind,1:(6-ik),iw]
      data_testing_total = extract.years(get(Area[iw+1]), years = (2011+ik):2016)$rate$total
      testing_ind_total  = data_testing_total != 0
      
      # me
      total_Area_me[ik,iw] = me(data_forc_total[testing_ind_total],     data_testing_total[testing_ind_total])
      
      # mae
      total_Area_mae[ik,iw] = mae(data_forc_total[testing_ind_total],   data_testing_total[testing_ind_total])
      
      # rmse
      total_Area_rmse[ik,iw] = rmse(data_forc_total[testing_ind_total], data_testing_total[testing_ind_total])
    }
  }
  return(list(res_Area = res_Area, 
              total_Area_me = total_Area_me, 
              total_Area_mae = total_Area_mae, 
              total_Area_rmse = total_Area_rmse, 
              train_residual_total = train_residual_total))
}

# forecasting mortality rates of Area total

dum_Area = mfts_back_test_Area_total(pcamethod = c("static"), year_horizon = 5)
mfts_Area_me_total = dum_Area$total_Area_me
mfts_Area_mae_total = dum_Area$total_Area_mae
mfts_Area_rmse_total = dum_Area$total_Area_rmse
res_Area = dum_Area$res_Area
mfts_Area_residual = dum_Area$train_residual_total
colnames(mfts_Area_me_total) = colnames(mfts_Area_mae_total) = colnames(mfts_Area_rmse_total) = Area[2:48]

for(ik in 1:47)
{
  assign(mfts_Area_forc_total[ik+1], res_Area[,,,ik])
}

for(ik in 1:47)
{
  tmp = list()
  for(ih in 1:5)
  {
    tmp[[ih]] = mfts_Area_residual[[ih]][,,ik]
  }
  assign(mfts_Area_train_residual_total[ik+1], tmp)
}

mfts_Australia_forc_total  = ind_Australia_forc_total
mfts_Australia_forc_female = ind_Australia_forc_female
mfts_Australia_forc_male   = ind_Australia_forc_male

mfts_Australia_train_residual_total = ind_Australia_train_residual_total
mfts_Australia_train_residual_female = ind_Australia_train_residual_female
mfts_Australia_train_residual_male = ind_Australia_train_residual_male

######################################################################
# Multivariate forecast of Region total  age-specific mortality rates
######################################################################

mfts_back_test_Region_total = function(fmethod = "classical", pcamethod = c("static", "dynamic"), year_horizon)
{
  total_Region_comb = total_Region_comb_pop = matrix(NA, 86*24, 11)
  for(iw in 1:11)
  {
    total_Region_comb[,iw] = as.numeric(get(Region_smooth[iw])$rate$total)
    total_Region_comb_pop[,iw] = as.numeric(get(Region_smooth[iw])$pop$total)
  }
  total_Region_comb_v2 = cbind(rep(1993:2016, each=86), rep(0:85, 24), total_Region_comb)
  total_Region_comb_pop_v2 = cbind(rep(1993:2016, each=86), rep(0:85, 24), total_Region_comb_pop)
  colnames(total_Region_comb_v2) = colnames(total_Region_comb_pop_v2) = c("Year", "Age", Region)
  
  write.table(total_Region_comb_v2, file = "total_Region_comb.txt", quote = FALSE, row.names = TRUE,
              col.names = TRUE)  
  write.table(total_Region_comb_pop_v2, file = "total_Region_comb_pop.txt", quote = FALSE, row.names = TRUE,
              col.names = TRUE)
  
  total_Region_comb_demogdata = read.demogdata("total_Region_comb.txt", "total_Region_comb_pop.txt",
                                               type = "mortality", label = "total_Region_comb", skip = 0)
  
  res_Region = array(NA, dim = c(year_horizon, 86, year_horizon, 11))
  
  train_residual_total = list()
  for(ij in 1:year_horizon)
  {
    ind_dat = extract.years(total_Region_comb_demogdata, years = 1993:(2010+ij))
    fun_forc = mfts(ind_dat, pcamethod = pcamethod, year_horizon = year_horizon)
    
    train_residual_total[[ij]] = array(NA, dim = c(86, (18+ij), 11))
    for(iwk in 1:11)
    {
      res_Region[,,ij,iwk] = t(fun_forc$fore_res[(86*(iwk-1)+1):(86*iwk),])
      
      train_residual_total[[ij]][,,iwk] = fun_forc$train_residual[(86*(iwk-1)+1):(86*iwk),]
    }
  }
  # Errors
  total_Region_me = total_Region_mae = total_Region_rmse = matrix(NA, year_horizon, 11)
  for(iw in 1:11)
  {
    for(ik in 1:year_horizon)
    {
      age_forc = 0:85
      age_obs  = get(Region[iw])$age
      age_ind  = age_forc %in% age_obs
      
      data_forc_total    = res_Region[ik,age_ind,1:(6-ik),iw]
      data_testing_total = extract.years(get(Region[iw]), years = (2011+ik):2016)$rate$total
      testing_ind_total  = data_testing_total != 0
      
      total_Region_me[ik,iw]   = me(data_forc_total[testing_ind_total],   data_testing_total[testing_ind_total])
      total_Region_mae[ik,iw]  = mae(data_forc_total[testing_ind_total],  data_testing_total[testing_ind_total])
      total_Region_rmse[ik,iw] = rmse(data_forc_total[testing_ind_total], data_testing_total[testing_ind_total])
    }
  }
  return(list(res_Region = res_Region, 
              train_residual_total = train_residual_total,
              total_Region_me = total_Region_me, 
              total_Region_mae = total_Region_mae,
              total_Region_rmse = total_Region_rmse))    
}

# forecasting mortality rates of Region total

dum_Region = mfts_back_test_Region_total(pcamethod = "static", year_horizon = 5)
mfts_Region_me_total = dum_Region$total_Region_me
mfts_Region_mae_total = dum_Region$total_Region_mae
mfts_Region_rmse_total = dum_Region$total_Region_rmse
mfts_Region_residual = dum_Region$train_residual_total
colnames(mfts_Region_me_total) = colnames(mfts_Region_rmse_total) = Region

for(ik in 1:11)
{
  assign(mfts_Region_forc_total[ik], dum_Region$res_Region[,,,ik])
  
}

for(ik in 1:11)
{
  tmp = list()
  for(ih in 1:5)
  {
    tmp[[ih]] = mfts_Region_residual[[ih]][,,ik]
  }
  assign(mfts_Region_train_residual_total[ik], tmp)
}


############################################
# Calculation of mfts point forecast errors
############################################

#  Area level

mfts_me_female_mean_overall = mfts_me_male_mean_overall =
  mfts_mae_female_mean_overall = mfts_mae_male_mean_overall = 
  mfts_rmse_female_mean_overall = mfts_rmse_male_mean_overall = NULL
for(ik in 1:48)
{
  # me
  mfts_me_female_mean_overall  = cbind(mfts_me_female_mean_overall, get(mfts_Area_me_female[ik]))
  mfts_me_male_mean_overall    = cbind(mfts_me_male_mean_overall,   get(mfts_Area_me_male[ik]))
  
  # me
  mfts_mae_female_mean_overall  = cbind(mfts_mae_female_mean_overall, get(mfts_Area_mae_female[ik]))
  mfts_mae_male_mean_overall    = cbind(mfts_mae_male_mean_overall,   get(mfts_Area_mae_male[ik]))
  
  # rmse
  mfts_rmse_female_mean_overall  = cbind(mfts_rmse_female_mean_overall, get(mfts_Area_rmse_female[ik]))
  mfts_rmse_male_mean_overall    = cbind(mfts_rmse_male_mean_overall,   get(mfts_Area_rmse_male[ik]))
}

# take an average

rowmeans_mfts_me_female_mean_overall = apply(mfts_me_female_mean_overall[,2:48], 1, mean)
rowmeans_mfts_me_male_mean_overall   = apply(mfts_me_male_mean_overall[,2:48], 1, mean)
rowmeans_mfts_me_total_mean_overall  = apply(mfts_Area_me_total, 1, mean)

rowmeans_mfts_mae_female_mean_overall = apply(mfts_mae_female_mean_overall[,2:48], 1, mean)
rowmeans_mfts_mae_male_mean_overall   = apply(mfts_mae_male_mean_overall[,2:48], 1, mean)
rowmeans_mfts_mae_total_mean_overall  = apply(mfts_Area_mae_total, 1, mean)

rowmeans_mfts_rmse_female_mean_overall = apply(mfts_rmse_female_mean_overall[,2:48], 1, mean)
rowmeans_mfts_rmse_male_mean_overall   = apply(mfts_rmse_male_mean_overall[,2:48], 1, mean)
rowmeans_mfts_rmse_total_mean_overall  = apply(mfts_Area_rmse_total, 1, mean)

# Region level

mfts_Region_me_female_mean_overall  = mfts_Region_me_male_mean_overall  =  
  mfts_Region_mae_female_mean_overall  = mfts_Region_mae_male_mean_overall  = 
  mfts_Region_rmse_female_mean_overall = mfts_Region_rmse_male_mean_overall = NULL
for(ik in 1:11)
{
  # me
  mfts_Region_me_female_mean_overall = cbind(mfts_Region_me_female_mean_overall, get(mfts_Region_me_female[ik]))
  mfts_Region_me_male_mean_overall   = cbind(mfts_Region_me_male_mean_overall,   get(mfts_Region_me_male[ik]))
  
  # mae
  mfts_Region_mae_female_mean_overall = cbind(mfts_Region_mae_female_mean_overall, get(mfts_Region_mae_female[ik]))
  mfts_Region_mae_male_mean_overall   = cbind(mfts_Region_mae_male_mean_overall,   get(mfts_Region_mae_male[ik]))
  
  # rmse
  mfts_Region_rmse_female_mean_overall = cbind(mfts_Region_rmse_female_mean_overall, get(mfts_Region_rmse_female[ik]))
  mfts_Region_rmse_male_mean_overall   = cbind(mfts_Region_rmse_male_mean_overall,   get(mfts_Region_rmse_male[ik]))
}

# take an average

rowmeans_mfts_Region_me_female_mean_overall = apply(mfts_Region_me_female_mean_overall[,1:11], 1, mean)
rowmeans_mfts_Region_me_male_mean_overall   = apply(mfts_Region_me_male_mean_overall[,1:11], 1, mean)
rowmeans_mfts_Region_me_total_mean_overall  = apply(mfts_Region_me_total, 1, mean)

rowmeans_mfts_Region_mae_female_mean_overall = apply(mfts_Region_mae_female_mean_overall[,1:11], 1, mean)
rowmeans_mfts_Region_mae_male_mean_overall   = apply(mfts_Region_mae_male_mean_overall[,1:11], 1, mean)
rowmeans_mfts_Region_mae_total_mean_overall  = apply(mfts_Region_mae_total, 1, mean)

rowmeans_mfts_Region_rmse_female_mean_overall = apply(mfts_Region_rmse_female_mean_overall[,1:11], 1, mean)
rowmeans_mfts_Region_rmse_male_mean_overall   = apply(mfts_Region_rmse_male_mean_overall[,1:11], 1, mean)
rowmeans_mfts_Region_rmse_total_mean_overall  = apply(mfts_Region_rmse_total, 1, mean)

#####################
# Summary of results
#####################

# Mean Errors

mfts_me_total_mean_overall =  ind_Australia_me_total 

mfts_all_level_err_me = cbind(mfts_me_total_mean_overall, 
                              apply(cbind(mfts_me_female_mean_overall[,1], mfts_me_male_mean_overall[,1]), 1,mean),
                              apply(mfts_Region_me_total, 1, mean), 
                              apply(cbind(mfts_Region_me_female_mean_overall, mfts_Region_me_male_mean_overall),1,mean),
                              rowmeans_mfts_me_total_mean_overall, 
                              apply(cbind(rowmeans_mfts_me_female_mean_overall, rowmeans_mfts_me_male_mean_overall),1,mean))

mfts_all_level_err_me_all = rbind(mfts_all_level_err_me, colMeans(mfts_all_level_err_me), apply(mfts_all_level_err_me, 2, median))

# Mean Absolute Errors

mfts_mae_total_mean_overall =  ind_Australia_mae_total 

mfts_all_level_err_mae = cbind(mfts_mae_total_mean_overall, 
                               apply(cbind(mfts_mae_female_mean_overall[,1], mfts_mae_male_mean_overall[,1]), 1,mean),
                               apply(mfts_Region_mae_total, 1, mean), 
                               apply(cbind(mfts_Region_mae_female_mean_overall, mfts_Region_mae_male_mean_overall),1,mean),
                               rowmeans_mfts_mae_total_mean_overall, 
                               apply(cbind(rowmeans_mfts_mae_female_mean_overall, rowmeans_mfts_mae_male_mean_overall),1,mean))

mfts_all_level_err_mae_all = rbind(mfts_all_level_err_mae, colMeans(mfts_all_level_err_mae), apply(mfts_all_level_err_mae, 2, median))

# Root Mean Square Errors

mfts_rmse_total_mean_overall = ind_Australia_rmse_total 

mfts_all_level_err_rmse = cbind(mfts_rmse_total_mean_overall, 
                                apply(cbind(mfts_rmse_female_mean_overall[,1], mfts_rmse_male_mean_overall[,1]), 1,mean),
                                apply(mfts_Region_rmse_total, 1, mean), 
                                apply(cbind(mfts_Region_rmse_female_mean_overall, mfts_Region_rmse_male_mean_overall),1,mean),
                                rowmeans_mfts_rmse_total_mean_overall, 
                                apply(cbind(rowmeans_mfts_rmse_female_mean_overall, rowmeans_mfts_rmse_male_mean_overall),1,mean))

mfts_all_level_err_rmse_all = rbind(mfts_all_level_err_rmse, colMeans(mfts_all_level_err_rmse), apply(mfts_all_level_err_rmse, 2, median))

colnames(mfts_all_level_err_me) = colnames(mfts_all_level_err_mae) = colnames(mfts_all_level_err_rmse) = c("Total", "Sex", "Region", "Region + Sex", "Area", "Area + Sex")
colnames(mfts_all_level_err_me_all) = colnames(mfts_all_level_err_mae_all) = colnames(mfts_all_level_err_rmse_all) = c("Total", "Sex", "Region", "Region + Sex", "Area", "Area + Sex")
rownames(mfts_all_level_err_me_all) = rownames(mfts_all_level_err_mae_all) = rownames(mfts_all_level_err_rmse_all) = c(1:5, "Mean", "Median")













