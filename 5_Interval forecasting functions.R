###########################################################
# constructing multivariate pointwise prediction intervals
###########################################################

# data_series: specific data series
# fh: forecast horizon
# nboot: number of bootstrap replication
# alpha: nominal coverage probability

##########
# Fix COB
##########

PI_fix <- function(dat, dat_unsmooth, pcamethod = c("static", "dynamic"), fh = 1, nboot = 1000, alpha = 0.8)
{
  # select ncomp_comb based on 95% of total variation; use all available data
  n_pop = length(names(dat$rate))
  
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
  
  if (pcamethod == "static")
  {
    ncomp_comb = min(head(which(cumsum(ftsm(fts(1:nrow(comb_object), comb_object), order = 20)$varprop) >= 0.95), 1), 9)
    
    if(ncomp_comb> 9 || ncomp_comb == 0)
    {
      warning("check ncomp_comb")
    }
    
    # calculate in-sample forecast curves
    fore_curve_region_total = array(NA, dim = c(length(dat$age), length(dat$year) - ncomp_comb - fh + 1, n_pop))
    
    for(ij in 1:(length(dat$year) - ncomp_comb - fh + 1))
    {
      
      rowmeans_object_one_step = sd_object_one_step = decenter_object_one_step = list()
      for(ik in 1:n_pop)
      {
        # compute mean and standard deviation functions
        dat_one_step = as.data.frame(log(dat$rate[[ik]][, 1:(ncomp_comb + ij - 1)]))
        rowmeans_object_one_step[[ik]] = rowMeans(dat_one_step, na.rm=TRUE)
        sd_object_one_step[[ik]] = apply(dat_one_step, 1, sd, na.rm=TRUE)
        
        # de-center functional data
        decenter_object_one_step[[ik]] = t(scale(t(dat_one_step), center = TRUE, scale = TRUE))
      }
      
      # h-step-ahead at each time
      comb_object_one_step = do.call(rbind, decenter_object_one_step)
      
      fore_ftsm = forecast(ftsm(fts(1:nrow(comb_object_one_step), comb_object_one_step), order = ncomp_comb), h = fh)
      fore_res  = exp(fore_ftsm$mean$y * do.call(c,sd_object_one_step) + do.call(c, rowmeans_object_one_step))
      
      # fill in in-sample forecast curves
      for(iwk in 1:n_pop)
      {
        fore_curve_region_total[,,iwk] = t(fore_res[(length(dat$age)*(iwk-1)+1):(length(dat$age)*iwk),fh])
      }
    }
  }
  
  if (pcamethod == "dynamic")
  {
    data_dum = comb_object
    
    C_0 = long_run_covariance_estimation(data_dum, H = 3, C0 = 3)
    eigen_decomp = eigen(C_0)
    ncomp_comb = min(head(which(cumsum(eigen_decomp$values)/sum(eigen_decomp$values) >= 0.95),1), 2)
    
    dynamic_basis = as.matrix(eigen_decomp$vectors[,1:ncomp_comb])
    dynamic_scores = t(dynamic_basis) %*% data_dum
    
    # calculate in-sample forecast curves
    fore_curve_region_total = array(NA, dim = c(length(dat$age), length(dat$year) - ncomp_comb - fh + 1, n_pop))
    for(ij in 1:(length(dat$year) - ncomp_comb - fh + 1))
    {
      rowmeans_object_one_step = sd_object_one_step = decenter_object_one_step = list()
      for(ik in 1:n_pop)
      {
        # compute mean and standard deviation functions
        dat_one_step = as.data.frame(log(dat$rate[[ik]][, 1:(ncomp_comb + ij - 1)]))
        rowmeans_object_one_step[[ik]] = rowMeans(dat_one_step, na.rm=TRUE)
        sd_object_one_step[[ik]] = apply(dat_one_step, 1, sd, na.rm=TRUE)
        
        # de-center functional data
        decenter_object_one_step[[ik]] = t(scale(t(dat_one_step), center = TRUE, scale = TRUE))
      }
      
      # h-step-ahead at each time
      comb_object_one_step = do.call(rbind, decenter_object_one_step)
      
      scores_fit = scores_fore = list()
      fore_ftsm_dyn = matrix(NA, nrow = nrow(comb_object_one_step), ncol = fh)
      
      for(ik in 1:ncomp_comb)
      {
        scores_fit[[ik]] = auto.arima(dynamic_scores[ik,])
        scores_fore[[ik]] = forecast(scores_fit[[ik]], h = fh)$mean
      }
      
      for(ih in 1:fh)
      {
        fore_ftsm_dyn[,ih] = dynamic_basis %*% unlist(lapply(scores_fore,`[[`,ih))
      }
      
      fore_res = exp(fore_ftsm_dyn * do.call(c,sd_object_one_step) + do.call(c, rowmeans_object_one_step))
      
      # fill in gender specific in-sample forecast curves
      for(iwk in 1:n_pop)
      {
        fore_curve_region_total[,,iwk] = t(fore_res[(length(dat$age)*(iwk-1)+1):(length(dat$age)*iwk),fh])
      }
    }
  }
  
  # holdout data samples
  
  true_dat = extract.years(dat_unsmooth, dat_unsmooth$year[(ncomp_comb + fh):length(dat_unsmooth$year)])
  true_dat_smooth = extract.years(dat, dat$year[(ncomp_comb + fh):length(dat$year)])
  
  holdout_val = holdout_smooth = array(NA, dim = c(length(dat_unsmooth$age), length(dat_unsmooth$year) - ncomp_comb - fh + 1, n_pop))
  age_ind = match(dat_unsmooth$age, dat$age)
  
  for(iwk in 1:n_pop)
  {
    holdout_val[,,iwk] = (true_dat$rate[[iwk]])
    holdout_smooth[,,iwk] = (true_dat_smooth$rate[[iwk]][age_ind,])
  }
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
  err_region_total = holdout_val - fore_curve_region_total[age_ind,,,drop = FALSE]
  
  
  # bootstrap error function
  err_boot_region_total = array(NA, dim = c(nrow(err_region_total), nboot, n_pop))
  for(iwk in 1:n_pop)
  {
    for(ij in 1:nboot)
    {
      err_boot_region_total[,ij,iwk] = err_region_total[, sample(1:ncol(err_region_total), 1, replace = TRUE), iwk]
    }
  }
  
  # constructing PI
  
  fore_ftsm = forecast(ftsm(fts(1:nrow(comb_object), comb_object), order = ncomp_comb), h = fh)
  fore_res  = exp(fore_ftsm$mean$y * do.call(c,sd_object) + do.call(c, rowmeans_object))
  fore_mfts_region_total = array(NA, dim = c(length(dat$age), 1, n_pop))
  for(iwk in 1:n_pop)
  {
    fore_mfts_region_total[,,iwk] = t(fore_res[(length(dat$age)*(iwk-1)+1):(length(dat$age)*iwk),fh])
  }
  
  boot_PI_region_total = array(NA, dim = c(length(dat_unsmooth$age), nboot, n_pop))
  for(iwk in 1:n_pop)
  {
    boot_PI_region_total[,,iwk] = err_boot_region_total[,,iwk] + matrix(rep(fore_mfts_region_total[age_ind,,iwk], nboot), nrow = length(dat_unsmooth$age), ncol = nboot)
  }
  
  
  return(list( boot_sample = boot_PI_region_total))
}


# for(region_index in 1:11)
# {
#   for(cob_ind in 1:19)
#   {
#     area_ind = region_list_ind[[region_index]]
#     
#     if(region_index == 11 && cob_ind == 14)
#     {
#       area_ind = area_ind[-3] # no observation at all for Area 23 COB 14
#     }
#     
#     region_comb = region_comb_pop = matrix(NA, 7*31, length(area_ind))
#     for(iw in 1:length(area_ind))
#     {
#       region_comb[,iw] = as.numeric(get(fts_area_fix_cob[[area_ind[iw]]][cob_ind])$rate$female)
#       region_comb_pop[,iw] = as.numeric(get(fts_area_fix_cob[[area_ind[iw]]][cob_ind])$pop$female)
#     }
#     region_comb_v2 = cbind(rep(year_interp, each=7), rep(age_all, 31), region_comb)
#     region_comb_pop_v2 = cbind(rep(year_interp, each=7), rep(age_all, 31), region_comb_pop)
#     
#     colnames(region_comb_v2) = colnames(region_comb_pop_v2) = c("Year", "Age", paste("Area_", area_ind, sep = ""))
#     
#     write.table(region_comb_v2, file = paste("R_", region_index, "_COB_", cob_ind, "_region_comb_unsmooth.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)  
#     write.table(region_comb_pop_v2, file =  paste("R_", region_index, "_COB_", cob_ind, "_region_comb_pop_unsmooth.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
#   }
# }


area_PI_fh_mfts <- function(region_index, fh, pcamethod = c("static", "dynamic"), nboot = 1000, cob_ind)
{
  if(region_index == 11 && cob_ind == 14)
  {
    n_area = length(region_list[[region_index]])-1 # no observation at all for Area 23 COB 14
  } else {
    n_area = length(region_list[[region_index]])
  }

  PI_boot_mfts_area = array(0, dim = c(7, nboot, (11-fh), n_area))

  data_series = read.demogdata(paste("R_", region_index, "_COB_", cob_ind, "_region_comb.txt", sep = ""), paste("R_", region_index, "_COB_", cob_ind, "_region_comb_pop.txt", sep = ""), type = "fertility", label = paste("R_", region_index, sep = ""), skip = 0)
  
  data_series_unsmooth = read.demogdata(paste("R_", region_index, "_COB_", cob_ind, "_region_comb_unsmooth.txt", sep = ""), paste("R_", region_index, "_COB_", cob_ind, "_region_comb_pop_unsmooth.txt", sep = ""), type = "fertility", label = paste("R_", region_index, "_unsmooth", sep = ""), skip = 0)
  
  for(ij in 1:(11-fh))
  {
    dum_pointwise = PI_fix(dat = extract.years(data_series, 1981:(2000+ij)), dat_unsmooth = extract.years(data_series_unsmooth, 1981:(2000+ij)), pcamethod = pcamethod, fh = fh, nboot = nboot)
    PI_boot_mfts_area[,,ij,]  = dum_pointwise$boot_sample
  }
  return(list(PI_boot_mfts = PI_boot_mfts_area))
}

test = area_PI_fh_mfts(region_index = 2, fh = 1, pcamethod = "dynamic", nboot = 20, cob_ind = 19)

# Area level

library(doParallel)
# R1
cl <- makeCluster(10) 
registerDoParallel(cl)

R1_PI_all_fix_cob = list()
for(cob in 1:19)
{
  R1_PI_all_fix_cob[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% area_PI_fh_mfts(region_index = 1,  pcamethod = "dynamic", fh = ik, cob_ind = cob, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# R2
cl <- makeCluster(10) 
registerDoParallel(cl)

R2_PI_all_fix_cob = list()
for(cob in 1:19)
{
  R2_PI_all_fix_cob[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% area_PI_fh_mfts(region_index = 2,  pcamethod = "dynamic", fh = ik, cob_ind = cob, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# R3
cl <- makeCluster(10) 
registerDoParallel(cl)

R3_PI_all_fix_cob = list()
for(cob in 1:19)
{
  R3_PI_all_fix_cob[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% area_PI_fh_mfts(region_index = 3,  pcamethod = "dynamic", fh = ik, cob_ind = cob, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# R4
cl <- makeCluster(10) 
registerDoParallel(cl)

R4_PI_all_fix_cob = list()
for(cob in 1:19)
{
  R4_PI_all_fix_cob[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% area_PI_fh_mfts(region_index = 4,  pcamethod = "dynamic", fh = ik, cob_ind = cob, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# R5
cl <- makeCluster(10) 
registerDoParallel(cl)

R5_PI_all_fix_cob = list()
for(cob in 1:19)
{
  R5_PI_all_fix_cob[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% area_PI_fh_mfts(region_index = 5,  pcamethod = "dynamic", fh = ik, cob_ind = cob, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# R6
cl <- makeCluster(10) 
registerDoParallel(cl)

R6_PI_all_fix_cob = list()
for(cob in 1:19)
{
  R6_PI_all_fix_cob[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% area_PI_fh_mfts(region_index = 6,  pcamethod = "dynamic", fh = ik, cob_ind = cob, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# R7
cl <- makeCluster(10) 
registerDoParallel(cl)

R7_PI_all_fix_cob = list()
for(cob in 1:19)
{
  R7_PI_all_fix_cob[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% area_PI_fh_mfts(region_index = 7,  pcamethod = "dynamic", fh = ik, cob_ind = cob, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# R8
cl <- makeCluster(10) 
registerDoParallel(cl)

R8_PI_all_fix_cob = list()
for(cob in 1:19)
{
  R8_PI_all_fix_cob[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% area_PI_fh_mfts(region_index = 8,  pcamethod = "dynamic", fh = ik, cob_ind = cob, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# R9
cl <- makeCluster(10) 
registerDoParallel(cl)

R9_PI_all_fix_cob = list()
for(cob in 1:19)
{
  R9_PI_all_fix_cob[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% area_PI_fh_mfts(region_index = 9,  pcamethod = "dynamic", fh = ik, cob_ind = cob, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# R10
cl <- makeCluster(10) 
registerDoParallel(cl)

R10_PI_all_fix_cob = list()
for(cob in 1:19)
{
  R10_PI_all_fix_cob[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% area_PI_fh_mfts(region_index = 10,  pcamethod = "dynamic", fh = ik, cob_ind = cob, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# R11
cl <- makeCluster(10) 
registerDoParallel(cl)

R11_PI_all_fix_cob = list()
for(cob in 15:19)
{
  R11_PI_all_fix_cob[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% area_PI_fh_mfts(region_index = 11,  pcamethod = "dynamic", fh = ik, cob_ind = cob, nboot = 1000)
}

stopCluster(cl)
rm(cl)

area_PI_fh_all = vector("character", length = 47)
for(i in 1:47)
{
  area_PI_fh_all[i] = paste("A", i, "_PI_fh_fix_COB", sep = "")
}

for(j in 1:11)
{
  area_ind = region_list_ind[[j]]
  nk = length(area_ind)
  for(k in 1:nk)
  {
    area_select = area_ind[k]
    temp_list = list()
    for(cob in 1:19)
    {
      if(area_select == 23 && cob == 14)
      {
        temp_list[[cob]] = NA
      } else {
        temp_list[[cob]] = list()
        
        if(j == 11 && k > 3)
        {
          for(h in 1:10)
          {
            temp_list[[cob]][[h]] = get(paste("R", j, "_PI_all_fix_cob", sep = ""))[[cob]][[h]]$PI_boot_mfts[,,,k-1]
          }
        } else {
          for(h in 1:10)
          {
            temp_list[[cob]][[h]] = get(paste("R", j, "_PI_all_fix_cob", sep = ""))[[cob]][[h]]$PI_boot_mfts[,,,k]
          }
        }
      }
    }
    assign(area_PI_fh_all[area_select], temp_list)
  }
}

# Region level

for(cob_ind in 1:19)
{
  region_comb = region_comb_pop = matrix(NA, 7*31, 11)
  for(iw in 1:11)
  {
    region_comb[,iw] = as.numeric(get(fts_region_fix_cob[[iw]][cob_ind])$rate$female)
    region_comb_pop[,iw] = as.numeric(get(fts_region_fix_cob[[iw]][cob_ind])$pop$female)
  }
  region_comb_v2 = cbind(rep(year_interp, each=7), rep(age_all, 31), region_comb)
  region_comb_pop_v2 = cbind(rep(year_interp, each=7), rep(age_all, 31), region_comb_pop)
  
  colnames(region_comb_v2) = colnames(region_comb_pop_v2) = c("Year", "Age", paste("Region_", 1:11, sep = ""))
  
  write.table(region_comb_v2, file = paste("national", "_COB_", cob_ind, "_comb_unsmooth.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(region_comb_pop_v2, file =  paste("national", "_COB_", cob_ind, "_comb_pop_unsmooth.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
}

region_PI_fh_mfts <- function(fh, pcamethod = c("static", "dynamic"), nboot = 1000, cob_ind)
{

  PI_boot_mfts_region = array(0, dim = c(7, nboot, (11-fh), 11))
  
  data_series = read.demogdata(paste("national", "_COB_", cob_ind, "_comb.txt", sep = ""), paste("national", "_COB_", cob_ind, "_comb_pop.txt", sep = ""), type = "fertility", label = paste("national", sep = ""), skip = 0)
  
  data_series_unsmooth = read.demogdata(paste("national", "_COB_", cob_ind, "_comb_unsmooth.txt", sep = ""), paste("national", "_COB_", cob_ind, "_comb_pop_unsmooth.txt", sep = ""), type = "fertility", label = paste("national", "_unsmooth", sep = ""), skip = 0)
  
  for(ij in 1:(11-fh))
  {
    dum_pointwise = PI_fix(dat = extract.years(data_series, 1981:(2000+ij)), dat_unsmooth = extract.years(data_series_unsmooth, 1981:(2000+ij)), pcamethod = pcamethod, fh = fh, nboot = nboot)
    PI_boot_mfts_region[,,ij,]  = dum_pointwise$boot_sample
  }
  return(list(PI_boot_mfts = PI_boot_mfts_region))
}


cl <- makeCluster(10)
registerDoParallel(cl)

national_PI_all_fix_cob = list()
for(cob in 1:19)
{
  national_PI_all_fix_cob[[cob]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% region_PI_fh_mfts(pcamethod = "dynamic", fh = ik, cob_ind = cob, nboot = 1000)
}

stopCluster(cl)
rm(cl)


region_PI_fh_all = vector("character", length = 11)
for(i in 1:11)
{
  region_PI_fh_all[i] = paste("R", i, "_PI_fh_fix_COB", sep = "")
}

for(j in 1:11)
{
  temp_list = list()
  for(cob in 1:19)
  {
    temp_list[[cob]] = list()
    for(h in 1:10)
    {
      temp_list[[cob]][[h]] = national_PI_all_fix_cob[[cob]][[h]]$PI_boot_mfts[,,,j]
    }
  }
  assign(region_PI_fh_all[j], temp_list)
  rm(temp_list)
}


##########
# Fix GEO
##########

# COB level

for(C_index in 1:10)
{
  for(geo_ind in 1:59)
  {
    cob_ind = COB_list[[C_index]]
    
    if(C_index == 6 && geo_ind == 23)
    {
      cob_ind = cob_ind[1] # no observation at all for Area 23 COB 14
    }

    region_comb = region_comb_pop = matrix(NA, 7*31, length(cob_ind))
    for(iw in 1:length(cob_ind))
    {
      region_comb[,iw] = as.numeric(get(fts_area_fix_geo[[cob_ind[iw]]][geo_ind])$rate$female)
      region_comb_pop[,iw] = as.numeric(get(fts_area_fix_geo[[cob_ind[iw]]][geo_ind])$pop$female)
    }
    region_comb_v2 = cbind(rep(year_interp, each=7), rep(age_all, 31), region_comb)
    region_comb_pop_v2 = cbind(rep(year_interp, each=7), rep(age_all, 31), region_comb_pop)

    colnames(region_comb_v2) = colnames(region_comb_pop_v2) = c("Year", "Age", paste("COB_", cob_ind, sep = ""))

    write.table(region_comb_v2, file = paste("C_", C_index, "_GEO_", geo_ind, "_comb_unsmooth.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
    write.table(region_comb_pop_v2, file =  paste("C_", C_index, "_GEO_", geo_ind, "_comb_pop_unsmooth.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
  }
}


COB_PI_fh_mfts <- function(C_index, fh, pcamethod = c("static", "dynamic"), nboot = 1000, geo_ind)
{
  if(C_index == 6 && geo_ind == 23)
  {
    cob_ind = COB_list[[C_index]][1] # no observation at all for C_index 6 COB 14
  } else {
    cob_ind = COB_list[[C_index]]
  }
  
  n_area = length(cob_ind)

  PI_boot_mfts_area = array(0, dim = c(7, nboot, (11-fh), n_area))
  
  data_series = read.demogdata(paste("C_", C_index, "_GEO_", geo_ind, "_comb.txt", sep = ""), paste("C_", C_index, "_GEO_", geo_ind, "_comb_pop.txt", sep = ""), type = "fertility", label = paste("C_", C_index, sep = ""), skip = 0)
  
  data_series_unsmooth = read.demogdata(paste("C_", C_index, "_GEO_", geo_ind, "_comb_unsmooth.txt", sep = ""), paste("C_", C_index, "_GEO_", geo_ind, "_comb_pop_unsmooth.txt", sep = ""), type = "fertility", label = paste("C_", C_index, "_unsmooth", sep = ""), skip = 0)
  
  for(ij in 1:(11-fh))
  {
    dum_pointwise = PI_fix(dat = extract.years(data_series, 1981:(2000+ij)), dat_unsmooth = extract.years(data_series_unsmooth, 1981:(2000+ij)), pcamethod = pcamethod, fh = fh, nboot = nboot)
    PI_boot_mfts_area[,,ij,]  = dum_pointwise$boot_sample
  }
  return(list(PI_boot_mfts = PI_boot_mfts_area))
}

test = COB_PI_fh_mfts(C_index = 5, fh = 1, pcamethod = "dynamic", nboot = 20, geo_ind = 1)

library(doParallel)

# C1
cl <- makeCluster(10) 
registerDoParallel(cl)

C1_PI_all_fix_geo = list()
for(geo in 1:59)
{
  C1_PI_all_fix_geo[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% COB_PI_fh_mfts(C_index = 1,  pcamethod = "dynamic", fh = ik, geo_ind = geo, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# C2 
cl <- makeCluster(10) 
registerDoParallel(cl)

C2_PI_all_fix_geo = list()
for(geo in 1:59)
{
  C2_PI_all_fix_geo[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% COB_PI_fh_mfts(C_index = 2,  pcamethod = "dynamic", fh = ik, geo_ind = geo, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# C3
cl <- makeCluster(10) 
registerDoParallel(cl)

C3_PI_all_fix_geo = list()
for(geo in 1:59)
{
  C3_PI_all_fix_geo[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% COB_PI_fh_mfts(C_index = 3,  pcamethod = "dynamic", fh = ik, geo_ind = geo, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# C4
cl <- makeCluster(10) 
registerDoParallel(cl)

C4_PI_all_fix_geo = list()
for(geo in 1:59)
{
  C4_PI_all_fix_geo[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% COB_PI_fh_mfts(C_index = 4,  pcamethod = "dynamic", fh = ik, geo_ind = geo, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# C5
cl <- makeCluster(10) 
registerDoParallel(cl)

C5_PI_all_fix_geo = list()
for(geo in 1:59)
{
  C5_PI_all_fix_geo[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% COB_PI_fh_mfts(C_index = 5,  pcamethod = "dynamic", fh = ik, geo_ind = geo, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# C6
cl <- makeCluster(10) 
registerDoParallel(cl)

C6_PI_all_fix_geo = list()
for(geo in 1:59)
{
  C6_PI_all_fix_geo[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% COB_PI_fh_mfts(C_index = 6,  pcamethod = "dynamic", fh = ik, geo_ind = geo, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# C7
cl <- makeCluster(10) 
registerDoParallel(cl)

C7_PI_all_fix_geo = list()
for(geo in 1:59)
{
  C7_PI_all_fix_geo[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% COB_PI_fh_mfts(C_index = 7,  pcamethod = "dynamic", fh = ik, geo_ind = geo, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# C8
cl <- makeCluster(10) 
registerDoParallel(cl)

C8_PI_all_fix_geo = list()
for(geo in 1:59)
{
  C8_PI_all_fix_geo[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% COB_PI_fh_mfts(C_index = 8,  pcamethod = "dynamic", fh = ik, geo_ind = geo, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# C9
cl <- makeCluster(10) 
registerDoParallel(cl)

C9_PI_all_fix_geo = list()
for(geo in 1:59)
{
  C9_PI_all_fix_geo[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% COB_PI_fh_mfts(C_index = 9,  pcamethod = "dynamic", fh = ik, geo_ind = geo, nboot = 1000)
}

stopCluster(cl)
rm(cl)

# C10
cl <- makeCluster(10) 
registerDoParallel(cl)

C10_PI_all_fix_geo = list()
for(geo in 1:59)
{
  C10_PI_all_fix_geo[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% COB_PI_fh_mfts(C_index = 10,  pcamethod = "dynamic", fh = ik, geo_ind = geo, nboot = 1000)
}

stopCluster(cl)
rm(cl)

cob_PI_fh_all = vector("character", length = 19)
for(i in 1:19)
{
  cob_PI_fh_all[i] = paste("COB_", i, "_PI_fh_fix_GEO", sep = "")
}

for(j in 1:10)
{
  cob_ind = COB_list[[j]]
  nk = length(cob_ind)
  for(k in 1:nk)
  {
    cob_select = cob_ind[k]
    temp_list = list()
    for(geo in 1:59)
    {
      if(cob_select == 14 && geo == 23)
      {
        temp_list[[geo]] = NA
      } else {
        temp_list[[geo]] = list()
        
        if(j == 6 && k > 1)
        {
          for(h in 1:10)
          {
            temp_list[[geo]][[h]] = get(paste("C", j, "_PI_all_fix_geo", sep = ""))[[geo]][[h]]$PI_boot_mfts[,,,k-1]
          }
        } else {
          for(h in 1:10)
          {
            temp_list[[geo]][[h]] = get(paste("C", j, "_PI_all_fix_geo", sep = ""))[[geo]][[h]]$PI_boot_mfts[,,,k]
          }
        }
      }
    }
    assign(cob_PI_fh_all[cob_select], temp_list)
  }
}


# C_group level

for(geo_ind in 1:59)
{

  if(geo_ind == 23)
  {
    cob_ind = c(1:5, 7:10) # no observation at all for Area 23 COB 14 C_group 6
  } else {
    cob_ind = 1:10
  }
  
  region_comb = region_comb_pop = matrix(NA, 7*31, length(cob_ind))
  for(iw in 1:length(cob_ind))
  {
    region_comb[,iw] = as.numeric(get(fts_region_fix_geo[[cob_ind[iw]]][geo_ind])$rate$female)
    region_comb_pop[,iw] = as.numeric(get(fts_region_fix_geo[[cob_ind[iw]]][geo_ind])$pop$female)
  }
  region_comb_v2 = cbind(rep(year_interp, each=7), rep(age_all, 31), region_comb)
  region_comb_pop_v2 = cbind(rep(year_interp, each=7), rep(age_all, 31), region_comb_pop)
  
  colnames(region_comb_v2) = colnames(region_comb_pop_v2) = c("Year", "Age", paste("C_", cob_ind, sep = ""))
  
  write.table(region_comb_v2, file = paste("national", "_GEO_", geo_ind, "_comb_unsmooth.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(region_comb_pop_v2, file =  paste("national", "_GEO_", geo_ind, "_comb_pop_unsmooth.txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE)
}

C_PI_fh_mfts <- function(fh, pcamethod = c("static", "dynamic"), nboot = 1000, geo_ind)
{
  if(geo_ind == 23)
  {
    n_c = 9 # no observation at all for C_index 6 COB 14
  } else {
    n_c = 10
  }
  
  PI_boot_mfts_c = array(0, dim = c(7, nboot, (11-fh), n_c))
  
  data_series = read.demogdata(paste("C_total_GEO_", geo_ind, "_comb.txt", sep = ""), paste("C_total_GEO_", geo_ind, "_comb_pop.txt", sep = ""), type = "fertility", label = paste("national", sep = ""), skip = 0)
  
  data_series_unsmooth = read.demogdata(paste("national", "_GEO_", geo_ind, "_comb_unsmooth.txt", sep = ""), paste("national", "_GEO_", geo_ind, "_comb_pop_unsmooth.txt", sep = ""), type = "fertility", label = paste("national", "_unsmooth", sep = ""), skip = 0)
  
  for(ij in 1:(11-fh))
  {
    dum_pointwise = PI_fix(dat = extract.years(data_series, 1981:(2000+ij)), dat_unsmooth = extract.years(data_series_unsmooth, 1981:(2000+ij)), pcamethod = pcamethod, fh = fh, nboot = nboot)
    PI_boot_mfts_c[,,ij,]  = dum_pointwise$boot_sample
  }
  return(list(PI_boot_mfts = PI_boot_mfts_c))
}

test = C_PI_fh_mfts(pcamethod = "dynamic", fh = ik, geo_ind = geo, nboot = 10)

cl <- makeCluster(10) 
registerDoParallel(cl)

ALL_C_PI_all_fix_geo = list()
for(geo in 1:59)
{
  ALL_C_PI_all_fix_geo[[geo]] = foreach(ik = 1:10, .packages = c("demography", "ftsa")) %dopar% C_PI_fh_mfts(pcamethod = "dynamic", fh = ik, geo_ind = geo, nboot = 1000)
}

stopCluster(cl)
rm(cl)

c_group_PI_fh_all = vector("character", length = 10)
for(i in 1:10)
{
  c_group_PI_fh_all[i] = paste("C_group_", i, "_PI_fh_fix_GEO", sep = "")
}

for(j in 1:10)
{
  temp_list = list()
  for(geo in 1:59)
  {
    temp_list[[geo]] = list()
    if(geo == 23 && j == 6)
    {
      for(h in 1:10)
      {
        temp_list[[geo]][[h]] = NA
      }
    }
    
    if(geo == 23 && j > 6)
    {
      for(h in 1:10)
      {
        temp_list[[geo]][[h]] = ALL_C_PI_all_fix_geo[[geo]][[h]]$PI_boot_mfts[,,,j-1]
      }
    } 
    
    if(geo == 23 && j < 6)
    {
      for(h in 1:10)
      {
        temp_list[[geo]][[h]] = ALL_C_PI_all_fix_geo[[geo]][[h]]$PI_boot_mfts[,,,j]
      }
    }
    
    if(geo != 23)
    {
      for(h in 1:10)
      {
        temp_list[[geo]][[h]] = ALL_C_PI_all_fix_geo[[geo]][[h]]$PI_boot_mfts[,,,j]
      }
    }
    
  }
  assign(c_group_PI_fh_all[j], temp_list)
  rm(temp_list)
}
