#################################################################################
# Univariate functional time series point forecast of subnatioal mortality rates
#################################################################################


library(demography)
library(ftsa)

##################
# independent FDM
##################

# Area level training residuals

ind_Area_train_residual_female = paste("ind", Area, "train_residual_female", sep = "_")
ind_Area_train_residual_male = paste("ind", Area, "train_residual_male", sep = "_")
ind_Area_train_residual_total = paste("ind", Area, "train_residual_total", sep = "_")

# Area level forecasts

ind_Area_forc_female = paste("ind", Area, "forc_female", sep = "_")
ind_Area_forc_male = paste("ind", Area, "forc_male", sep = "_")
ind_Area_forc_total = paste("ind", Area, "forc_total", sep = "_")

# Mean errors at Area level

ind_Area_me_female = paste("ind", Area, "me_female", sep = "_")
ind_Area_me_male = paste("ind", Area, "me_male", sep = "_")
ind_Area_me_total = paste("ind", Area, "me_total", sep = "_")

# Mean absolute errors at Area level

ind_Area_mae_female = paste("ind", Area, "mae_female", sep = "_")
ind_Area_mae_male = paste("ind", Area, "mae_male", sep = "_")
ind_Area_mae_total = paste("ind", Area, "mae_total", sep = "_")

# Root mean square errors at Area level

ind_Area_rmse_female = paste("ind", Area, "rmse_female", sep = "_")
ind_Area_rmse_male = paste("ind", Area, "rmse_male", sep = "_")
ind_Area_rmse_total = paste("ind", Area, "rmse_total", sep = "_")

# Region level training residuals

ind_Region_train_residual_female = paste("ind", Region, "train_residual_female", sep = "_")
ind_Region_train_residual_male = paste("ind", Region, "train_residual_male", sep = "_")
ind_Region_train_residual_total = paste("ind", Region, "train_residual_total", sep = "_")

# Region level forecasts

ind_Region_forc_female = paste("ind", Region, "forc_female", sep = "_")
ind_Region_forc_male = paste("ind", Region, "forc_male", sep = "_")
ind_Region_forc_total = paste("ind", Region, "forc_total", sep = "_")

# Mean errors at Region level

ind_Region_me_female = paste("ind", Region, "me_female", sep = "_")
ind_Region_me_male = paste("ind", Region, "me_male", sep = "_")
ind_Region_me_total = paste("ind", Region, "me_total", sep = "_")

# Mean absolute errors at Region level

ind_Region_mae_female = paste("ind", Region, "mae_female", sep = "_")
ind_Region_mae_male = paste("ind", Region, "mae_male", sep = "_")
ind_Region_mae_total = paste("ind", Region, "mae_total", sep = "_")

# Root mean square errors at Region level

ind_Region_rmse_female = paste("ind", Region, "rmse_female", sep = "_")
ind_Region_rmse_male = paste("ind", Region, "rmse_male", sep = "_")
ind_Region_rmse_total = paste("ind", Region, "rmse_total", sep = "_")


# function for Region forecast

me   = ftsa:::me
mae  = ftsa:::mae
rmse = ftsa:::rmse

ind_back_test_Region <- function(iw, year_horizon, fmethod = c("classical", "M"), alpha = 0.2)
{
  fmethod = match.arg(fmethod)
  train_residual_female = train_residual_male = train_residual_total = list()
  
  res_male = res_male_lb = res_male_ub = 
    res_female = res_female_lb = res_female_ub = 
    res_total = res_total_lb = res_total_ub = array(NA, dim = c(year_horizon, 86, year_horizon))
  
  for(j in 1:year_horizon)
  {
    ind_dat = extract.years(get(Region_smooth[iw]), years = 1993:(2010+j))
    
    fdm_female_order = head(which(round(cumsum(fdm(ind_dat, series = "female", method = fmethod, lambda = 2.33, order = length(ind_dat$year))$varprop),3)>=0.95),1)
    fdm_female = fdm(ind_dat, series = "female", order = fdm_female_order, method = fmethod, lambda = 2.33)
    train_residual_female[[j]] = ind_dat$rate$female - exp(fdm_female$fitted$y)
    fun_forc_female  = forecast(fdm_female, h = year_horizon)        
    
    fdm_male_order = head(which(round(cumsum(fdm(ind_dat, series = "male", method = fmethod, lambda = 2.33, order = length(ind_dat$year))$varprop),3)>=0.95),1)
    fdm_male = fdm(ind_dat, series = "male", order = fdm_male_order, method = fmethod, lambda = 2.33)
    train_residual_male[[j]] = ind_dat$rate$male - exp(fdm_male$fitted$y)
    fun_forc_male  = forecast(fdm_male, h = year_horizon)
    
    fdm_total_order = head(which(round(cumsum(fdm(ind_dat, series = "total", method = fmethod, lambda = 2.33, order = length(ind_dat$year))$varprop),3)>=0.95),1)
    fdm_total = fdm(ind_dat, series = "total", order = fdm_total_order, method = fmethod, lambda = 2.33)
    train_residual_total[[j]] = ind_dat$rate$total - exp(fdm_total$fitted$y)
    fun_forc_total  = forecast(fdm_total, h = year_horizon)
    
    res_female[,,j]    = t(fun_forc_female$rate$female)
    res_female_lb[,,j] = t(fun_forc_female$rate$lower)
    res_female_ub[,,j] = t(fun_forc_female$rate$upper)
    
    res_male[,,j]    = t(fun_forc_male$rate$male)
    res_male_lb[,,j] = t(fun_forc_male$rate$lower)
    res_male_ub[,,j] = t(fun_forc_male$rate$upper)
    
    res_total[,,j]    = t(fun_forc_total$rate$total)
    res_total_lb[,,j] = t(fun_forc_total$rate$lower)
    res_total_ub[,,j] = t(fun_forc_total$rate$upper)
  }
  
  
  # Errors
  female_me = male_me = total_me = female_mae = male_mae = total_mae = female_rmse = male_rmse = total_rmse = vector("numeric",year_horizon)
  for(k in 1:year_horizon)
  {
    age_forc = 0:85
    age_obs  = get(Region[iw])$age
    age_ind  = age_forc %in% age_obs
    
    data_forc_female    = res_female[k,age_ind,1:(6-k)]
    data_testing_female = extract.years(get(Region[iw]), years = (2011+k):2016)$rate$female
    testing_ind_female  = data_testing_female != 0
    
    data_forc_male    = res_male[k,age_ind,1:(6-k)]
    data_testing_male = extract.years(get(Region[iw]), years = (2011+k):2016)$rate$male
    testing_ind_male  = data_testing_male != 0
    
    data_forc_total    = res_total[k,age_ind,1:(6-k)]
    data_testing_total = extract.years(get(Region[iw]), years = (2011+k):2016)$rate$total
    testing_ind_total  = data_testing_total != 0
    
    female_me[k]  = me(data_forc_female[testing_ind_female], data_testing_female[testing_ind_female])
    male_me[k]    = me(data_forc_male[testing_ind_male],     data_testing_male[testing_ind_male])
    total_me[k]   = me(data_forc_total[testing_ind_total],  data_testing_total[testing_ind_total])
    
    female_mae[k]  = mae(data_forc_female[testing_ind_female],  data_testing_female[testing_ind_female])
    male_mae[k]    = mae(data_forc_male[testing_ind_male],      data_testing_male[testing_ind_male])
    total_mae[k]   = mae(data_forc_total[testing_ind_total],   data_testing_total[testing_ind_total])
    
    female_rmse[k]  = rmse(data_forc_female[testing_ind_female], data_testing_female[testing_ind_female])
    male_rmse[k]    = rmse(data_forc_male[testing_ind_male],     data_testing_male[testing_ind_male])
    total_rmse[k]   = rmse(data_forc_total[testing_ind_total],  data_testing_total[testing_ind_total])
  }	
  return(list(res_female = res_female, res_male = res_male, res_total = res_total,
              train_residual_female = train_residual_female, train_residual_male = train_residual_male, train_residual_total = train_residual_total,
              female_me = female_me, male_me = male_me, total_me = total_me,
              female_mae = female_mae, male_mae = male_mae, total_mae = total_mae,
              female_rmse = female_rmse, male_rmse = male_rmse, total_rmse = total_rmse))
}

# function for Area forecast

ind_back_test_Area <- function(iw, year_horizon, fmethod = c("classical", "M"), alpha = 0.2)
{
  fmethod = match.arg(fmethod)
  train_residual_female = train_residual_male = train_residual_total = list()
  
  res_male  = res_male_lb  = res_male_ub  = 
    res_female = res_female_lb = res_female_ub = 
    res_total = res_total_lb = res_total_ub = array(NA, dim = c(year_horizon, 86, year_horizon))
  
  for(j in 1:year_horizon)
  {
    ind_dat = extract.years(get(Area_smooth[iw]), years = 1993:(2010+j))
    
    fdm_female_order = head(which(round(cumsum(fdm(ind_dat, series = "female", method = fmethod, lambda = 2.33, order = length(ind_dat$year))$varprop),3)>=0.95),1)
    fdm_female = fdm(ind_dat, series = "female", order = fdm_female_order, method = fmethod, lambda = 2.33)
    train_residual_female[[j]] = ind_dat$rate$female - exp(fdm_female$fitted$y)
    fun_forc_female  = forecast(fdm_female, h = year_horizon)   
    
    fdm_male_order = head(which(round(cumsum(fdm(ind_dat, series = "male", method = fmethod, lambda = 2.33, order = length(ind_dat$year))$varprop),3)>=0.95),1)
    fdm_male = fdm(ind_dat, series = "male", order = fdm_male_order, method = fmethod, lambda = 2.33)
    train_residual_male[[j]] = ind_dat$rate$male - exp(fdm_male$fitted$y)
    fun_forc_male  = forecast(fdm_male, h = year_horizon)
    
    fdm_total_order = head(which(round(cumsum(fdm(ind_dat, series = "total", method = fmethod, lambda = 2.33, order = length(ind_dat$year))$varprop),3)>=0.95),1)
    fdm_total = fdm(ind_dat, series = "total", order = fdm_total_order, method = fmethod, lambda = 2.33)
    train_residual_total[[j]] = ind_dat$rate$total - exp(fdm_total$fitted$y)
    fun_forc_total  = forecast(fdm_total, h = year_horizon)
    
    res_female[,,j] = t(fun_forc_female$rate$female)
    res_female_lb[,,j] = t(fun_forc_female$rate$lower)
    res_female_ub[,,j] = t(fun_forc_female$rate$upper)
    
    res_male[,,j] = t(fun_forc_male$rate$male)
    res_male_lb[,,j] = t(fun_forc_male$rate$lower)
    res_male_ub[,,j] = t(fun_forc_male$rate$upper)
    
    res_total[,,j] = t(fun_forc_total$rate$total)
    res_total_lb[,,j] = t(fun_forc_total$rate$lower)
    res_total_ub[,,j] = t(fun_forc_total$rate$upper)
  }
  # Errors
  female_me = male_me = total_me = female_mae = male_mae = total_mae = female_rmse = male_rmse = total_rmse = vector("numeric",year_horizon)
  for(k in 1:year_horizon)
  {    
    age_forc = 0:85
    age_obs  = get(Area[iw])$age
    age_ind  = age_forc %in% age_obs
    
    data_forc_female    = res_female[k,age_ind,1:(6-k)]
    data_testing_female = extract.years(get(Area[iw]), years = (2011+k):2016)$rate$female
    testing_ind_female  = data_testing_female != 0
    
    data_forc_male    = res_male[k,age_ind,1:(6-k)]
    data_testing_male = extract.years(get(Area[iw]), years = (2011+k):2016)$rate$male
    testing_ind_male  = data_testing_male != 0
    
    data_forc_total    = res_total[k,age_ind,1:(6-k)]
    data_testing_total = extract.years(get(Area[iw]), years = (2011+k):2016)$rate$total
    testing_ind_total  = data_testing_total != 0
    
    female_me[k]  = me(data_forc_female[testing_ind_female],  data_testing_female[testing_ind_female])
    male_me[k]    = me(data_forc_male[testing_ind_male],      data_testing_male[testing_ind_male])
    total_me[k]   = me(data_forc_total[testing_ind_total],    data_testing_total[testing_ind_total])
    
    female_mae[k]  = mae(data_forc_female[testing_ind_female],  data_testing_female[testing_ind_female])
    male_mae[k]    = mae(data_forc_male[testing_ind_male],      data_testing_male[testing_ind_male])
    total_mae[k]   = mae(data_forc_total[testing_ind_total],    data_testing_total[testing_ind_total])
    
    female_rmse[k]  = rmse(data_forc_female[testing_ind_female],  data_testing_female[testing_ind_female])
    male_rmse[k]    = rmse(data_forc_male[testing_ind_male],      data_testing_male[testing_ind_male])
    total_rmse[k]   = rmse(data_forc_total[testing_ind_total],    data_testing_total[testing_ind_total])
  }	
  
  # Interval score	  
  female_score = male_score = total_score = vector("numeric",year_horizon)
  for(k in 1:year_horizon)
  {
    age_forc = 0:85
    age_obs  = get(Area[iw])$age
    age_ind  = age_forc %in% age_obs

    # female
    
    test_female = extract.years(get(Area[iw]), years = (2011+k):2016)$rate$female	  
    test_female_notmissing = test_female != 0
    
    lb_ind = ifelse(test_female < res_female_lb[k,age_ind,1:(6-k)], 1, 0)
    ub_ind = ifelse(test_female > res_female_ub[k,age_ind,1:(6-k)], 1, 0)
    female_score[k] = mean((res_female_ub[k,age_ind,1:(6-k)] - res_female_lb[k,age_ind,1:(6-k)])[test_female_notmissing] 
                           + (2/alpha * (res_female_lb[k,age_ind,1:(6-k)] - test_female) * lb_ind)[test_female_notmissing] 
                           + (2/alpha * (test_female - res_female_ub[k,age_ind,1:(6-k)]) * ub_ind)[test_female_notmissing])
    
    # male
    
    test_male = extract.years(get(Area[iw]), years = (2011+k):2016)$rate$male	  
    test_male_notmissing = test_male != 0
    
    lb_ind = ifelse(test_male < res_male_lb[k,age_ind,1:(6-k)], 1, 0)
    ub_ind = ifelse(test_male > res_male_ub[k,age_ind,1:(6-k)], 1, 0)
    male_score[k] = mean((res_male_ub[k,age_ind,1:(6-k)] - res_male_lb[k,age_ind,1:(6-k)])[test_male_notmissing] 
                         + (2/alpha * (res_male_lb[k,age_ind,1:(6-k)] - test_male) * lb_ind)[test_male_notmissing] 
                         + (2/alpha * (test_male - res_male_ub[k,age_ind,1:(6-k)]) * ub_ind)[test_male_notmissing])
    
    # total
    
    test_total = extract.years(get(Area[iw]), years = (2011+k):2016)$rate$total	 
    test_total_notmissing = test_total != 0
    
    lb_ind = ifelse(test_total < res_total_lb[k,age_ind,1:(6-k)], 1, 0)
    ub_ind = ifelse(test_total > res_total_ub[k,age_ind,1:(6-k)], 1, 0)
    total_score[k] = mean((res_total_ub[k,age_ind,1:(6-k)] - res_total_lb[k,age_ind,1:(6-k)])[test_total_notmissing]
                          + (2/alpha * (res_total_lb[k,age_ind,1:(6-k)] - test_total) * lb_ind)[test_total_notmissing] 
                          + (2/alpha * (test_total - res_total_ub[k,age_ind,1:(6-k)]) * ub_ind)[test_total_notmissing])
  }
  return(list(res_female = res_female, res_male = res_male, res_total = res_total,
              train_residual_female = train_residual_female, train_residual_male = train_residual_male, train_residual_total = train_residual_total,
              female_me = female_me, male_me = male_me, total_me = total_me,
              female_mae = female_mae, male_mae = male_mae, total_mae = total_mae,
              female_rmse = female_rmse, male_rmse = male_rmse, total_rmse = total_rmse,
              female_score = female_score, male_score = male_score, total_score = total_score))
}

##############################################################
# Univariate forecast of age-specific mortality rates by Area
##############################################################

for(ik in 1:48)
{
  dum = ind_back_test_Area(iw = ik, year_horizon = 5, fmethod = "classical")
  assign(ind_Area_forc_female[ik], dum$res_female)
  assign(ind_Area_forc_male[ik],   dum$res_male)
  assign(ind_Area_forc_total[ik],  dum$res_total)
  
  assign(ind_Area_train_residual_female[ik], dum$train_residual_female)
  assign(ind_Area_train_residual_male[ik], dum$train_residual_male)
  assign(ind_Area_train_residual_total[ik], dum$train_residual_total)
  
  assign(ind_Area_me_female[ik], dum$female_me)
  assign(ind_Area_me_male[ik],   dum$male_me)
  assign(ind_Area_me_total[ik],  dum$total_me) 
  
  assign(ind_Area_mae_female[ik], dum$female_mae)
  assign(ind_Area_mae_male[ik],   dum$male_mae)
  assign(ind_Area_mae_total[ik],  dum$total_mae)   
  
  assign(ind_Area_rmse_female[ik], dum$female_rmse)
  assign(ind_Area_rmse_male[ik],   dum$male_rmse)
  assign(ind_Area_rmse_total[ik],  dum$total_rmse)   

  rm(dum)
}


###########################################################
# Univariate forecast of age-specific mortality by regions
###########################################################

for(ik in 1:11)
{
  dum = ind_back_test_Region(iw = ik, year_horizon = 5)
  assign(ind_Region_forc_female[ik], dum$res_female)
  assign(ind_Region_forc_male[ik], dum$res_male)
  assign(ind_Region_forc_total[ik], dum$res_total)
  
  assign(ind_Region_train_residual_female[ik], dum$train_residual_female)
  assign(ind_Region_train_residual_male[ik], dum$train_residual_male)
  assign(ind_Region_train_residual_total[ik], dum$train_residual_total)
  
  assign(ind_Region_me_female[ik], dum$female_me)
  assign(ind_Region_me_male[ik],   dum$male_me)
  assign(ind_Region_me_total[ik],  dum$total_me)
  
  assign(ind_Region_mae_female[ik], dum$female_mae)
  assign(ind_Region_mae_male[ik],   dum$male_mae)
  assign(ind_Region_mae_total[ik],  dum$total_mae)
  
  assign(ind_Region_rmse_female[ik], dum$female_rmse)
  assign(ind_Region_rmse_male[ik],   dum$male_rmse)
  assign(ind_Region_rmse_total[ik],  dum$total_rmse)

  rm(dum)
}

#####################
# Summary of results
#####################

# Region + Sex

ind_Region_me_female_mean_overall = ind_Region_me_male_mean_overall = ind_Region_me_total_mean_overall =
  ind_Region_mae_female_mean_overall = ind_Region_mae_male_mean_overall = ind_Region_mae_total_mean_overall = 
  ind_Region_rmse_female_mean_overall = ind_Region_rmse_male_mean_overall = ind_Region_rmse_total_mean_overall = NULL

for(ik in 1:11)  
{
  ind_Region_me_female_mean_overall = cbind(ind_Region_me_female_mean_overall, get(ind_Region_me_female[ik]))
  ind_Region_me_male_mean_overall   = cbind(ind_Region_me_male_mean_overall,   get(ind_Region_me_male[ik]))
  ind_Region_me_total_mean_overall  = cbind(ind_Region_me_total_mean_overall,   get(ind_Region_me_total[ik]))
  
  ind_Region_mae_female_mean_overall = cbind(ind_Region_mae_female_mean_overall, get(ind_Region_mae_female[ik]))
  ind_Region_mae_male_mean_overall   = cbind(ind_Region_mae_male_mean_overall,   get(ind_Region_mae_male[ik]))
  ind_Region_mae_total_mean_overall  = cbind(ind_Region_mae_total_mean_overall,   get(ind_Region_mae_total[ik]))
  
  ind_Region_rmse_female_mean_overall = cbind(ind_Region_rmse_female_mean_overall, get(ind_Region_rmse_female[ik]))
  ind_Region_rmse_male_mean_overall   = cbind(ind_Region_rmse_male_mean_overall,   get(ind_Region_rmse_male[ik]))
  ind_Region_rmse_total_mean_overall  = cbind(ind_Region_rmse_total_mean_overall,   get(ind_Region_rmse_total[ik]))
}

# Area + Sex

ind_me_female_mean_overall = ind_me_male_mean_overall = ind_me_total_mean_overall =
  ind_mae_female_mean_overall = ind_mae_male_mean_overall = ind_mae_total_mean_overall = 
  ind_rmse_female_mean_overall = ind_rmse_male_mean_overall = ind_rmse_total_mean_overall = NULL

for(ik in 1:48)
{
  ind_me_female_mean_overall  = cbind(ind_me_female_mean_overall, get(ind_Area_me_female[ik]))
  ind_me_male_mean_overall    = cbind(ind_me_male_mean_overall,   get(ind_Area_me_male[ik]))
  ind_me_total_mean_overall   = cbind(ind_me_total_mean_overall,   get(ind_Area_me_total[ik]))
  
  ind_mae_female_mean_overall  = cbind(ind_mae_female_mean_overall, get(ind_Area_mae_female[ik]))
  ind_mae_male_mean_overall    = cbind(ind_mae_male_mean_overall,   get(ind_Area_mae_male[ik]))
  ind_mae_total_mean_overall   = cbind(ind_mae_total_mean_overall,   get(ind_Area_mae_total[ik]))
  
  ind_rmse_female_mean_overall  = cbind(ind_rmse_female_mean_overall, get(ind_Area_rmse_female[ik]))
  ind_rmse_male_mean_overall    = cbind(ind_rmse_male_mean_overall,   get(ind_Area_rmse_male[ik]))
  ind_rmse_total_mean_overall   = cbind(ind_rmse_total_mean_overall,   get(ind_Area_rmse_total[ik]))
}

##########################
# Averaging over 47 Areas
##################$#######

rowmeans_ind_me_female_mean_overall = apply(ind_me_female_mean_overall[,2:48], 1, mean)
rowmeans_ind_me_male_mean_overall   = apply(ind_me_male_mean_overall[,2:48], 1, mean)
rowmeans_ind_me_total_mean_overall  = apply(ind_me_total_mean_overall[,2:48], 1, mean)

rowmeans_ind_mae_female_mean_overall = apply(ind_mae_female_mean_overall[,2:48], 1, mean)
rowmeans_ind_mae_male_mean_overall   = apply(ind_mae_male_mean_overall[,2:48], 1, mean)
rowmeans_ind_mae_total_mean_overall  = apply(ind_mae_total_mean_overall[,2:48], 1, mean)

rowmeans_ind_rmse_female_mean_overall = apply(ind_rmse_female_mean_overall[,2:48], 1, mean)
rowmeans_ind_rmse_male_mean_overall   = apply(ind_rmse_male_mean_overall[,2:48], 1, mean)
rowmeans_ind_rmse_total_mean_overall  = apply(ind_rmse_total_mean_overall[,2:48], 1, mean)

#####################
# Summary of results
#####################

ind_all_level_err_me = cbind(ind_me_total_mean_overall[,1], apply(cbind(ind_me_female_mean_overall[,1], ind_me_male_mean_overall[,1]), 1,mean),
                             apply(ind_Region_me_total_mean_overall, 1, mean), apply(cbind(ind_Region_me_female_mean_overall, ind_Region_me_male_mean_overall),1,mean),
                             rowmeans_ind_me_total_mean_overall, apply(cbind(rowmeans_ind_me_female_mean_overall, rowmeans_ind_me_male_mean_overall),1,mean))
ind_all_level_err_me_all = rbind(ind_all_level_err_me, colMeans(ind_all_level_err_me), apply(ind_all_level_err_me, 2, median))


ind_all_level_err_mae = cbind(ind_mae_total_mean_overall[,1], apply(cbind(ind_mae_female_mean_overall[,1], ind_mae_male_mean_overall[,1]), 1,mean),
                              apply(ind_Region_mae_total_mean_overall, 1, mean), apply(cbind(ind_Region_mae_female_mean_overall, ind_Region_mae_male_mean_overall),1,mean),
                              rowmeans_ind_mae_total_mean_overall, apply(cbind(rowmeans_ind_mae_female_mean_overall, rowmeans_ind_mae_male_mean_overall),1,mean))
ind_all_level_err_mae_all = rbind(ind_all_level_err_mae, colMeans(ind_all_level_err_mae), apply(ind_all_level_err_mae, 2, median))


ind_all_level_err_rmse = cbind(ind_rmse_total_mean_overall[,1], apply(cbind(ind_rmse_female_mean_overall[,1], ind_rmse_male_mean_overall[,1]), 1,mean),
                               apply(ind_Region_rmse_total_mean_overall, 1, mean), apply(cbind(ind_Region_rmse_female_mean_overall, ind_Region_rmse_male_mean_overall),1,mean),
                               rowmeans_ind_rmse_total_mean_overall, apply(cbind(rowmeans_ind_rmse_female_mean_overall, rowmeans_ind_rmse_male_mean_overall),1,mean))
ind_all_level_err_rmse_all = rbind(ind_all_level_err_rmse, colMeans(ind_all_level_err_rmse), apply(ind_all_level_err_rmse, 2, median))

colnames(ind_all_level_err_rmse) = colnames(ind_all_level_err_mae) = colnames(ind_all_level_err_me)  =  c("Total", "Sex", "Region", "Region + Sex", "Area", "Area + Sex")

colnames(ind_all_level_err_rmse_all) = colnames(ind_all_level_err_mae_all) = colnames(ind_all_level_err_me_all)  =  c("Total", "Sex", "Region", "Region + Sex", "Area", "Area + Sex")
rownames(ind_all_level_err_me_all) = rownames(ind_all_level_err_mae_all) = rownames(ind_all_level_err_rmse_all)=c(1:5, "Mean", "Median")








