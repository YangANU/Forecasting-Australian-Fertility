################
# Rainbow plots
################

savepdf("Aus_female", width = 12, height = 10, toplines = 0.9)
plot(extract.years(Australia, 1993:2016), series = "female", lwd = 1, main="Australia: female death rates (1993-2016)")
dev.off()

savepdf("Aus_male", width = 12, height = 10, toplines = 0.9)
plot(extract.years(Australia, 1993:2016), series = "male", lwd = 1, main="Australia: male death rates (1993-2016)")
dev.off()

savepdf("Aus_female_smooth", width = 12, height = 10, toplines = 0.9)
plot(extract.years(Australia_smooth, 1993:2016), series = "female", lwd = 1, main="")
dev.off()

savepdf("Aus_male_smooth", width = 12, height = 10, toplines = 0.9)
plot(extract.years(Australia_smooth, 1993:2016), series = "male", lwd = 1, main = "")
dev.off()

#############################
# ggplot for point forecasts
#############################

library(forcats)
library(ggplot2)
library(tidyverse)

point_forecasts_mae_series <- map2_dfr(1:6, c("Australia", "Sex", "Region", "Region + Sex", "Area", "Area + Sex"),
                                        ~tibble(
                                          variety = forcats::as_factor(rep(c("Base", "BU", "OP", "MinT", "Comb_av", "Base", "BU", "OP", "MinT", "Comb_av"), each = 5)),
                                          Method = forcats::as_factor(rep(c("FPCA", "MFPCA"), each = 25)),
                                          series = c(ind_all_level_err_mae[1:5, .x]*100,
                                                     BU_all_level_err_mae_ind[1:5, .x]*100,
                                                     optim_all_level_err_mae_ind[1:5, .x]*100,
                                                     optim_mint_all_level_err_mae_ind[1:5, .x]*100,
                                                     Comb_Av_all_level_err_mae_ind[1:5, .x]*100,
                                                     mfts_all_level_err_mae[1:5, .x]*100,
                                                     BU_all_level_err_mae_mfts[1:5, .x]*100,
                                                     optim_all_level_err_mae_mfts[1:5, .x]*100,
                                                     optim_mint_all_level_err_mae_mfts[1:5, .x]*100,
                                                     Comb_Av_all_level_err_mae_mfts[1:5, .x]*100),
                                          type = factor(paste(.y, "series"), 
                                                        levels = c("Australia series", "Region series", "Area series", "Sex series", "Region + Sex series", "Area + Sex series"))))

ggplot(point_forecasts_mae_series, aes(x = variety, y = series, fill = Method)) +
  geom_boxplot() + ylab("MAFE (x 100)") + xlab("") + facet_wrap(~type) + scale_y_log10()
ggsave("point_forecasts_mae.pdf", height = 8.27, width = 11.69, units = "in")


point_forecasts_rmse_series <- map2_dfr(1:6, c("Australia", "Sex", "Region", "Region + Sex", "Area", "Area + Sex"),
                                ~tibble(
                                  variety = forcats::as_factor(rep(c("Base", "BU", "OP", "MinT", "Comb_av", "Base", "BU", "OP", "MinT", "Comb_av"), each = 5)),
                                  Method = forcats::as_factor(rep(c("FPCA", "MFPCA"), each = 25)),
                                  series = c(ind_all_level_err_rmse[1:5, .x]*100,
                                             BU_all_level_err_rmse_ind[1:5, .x]*100,
                                             optim_all_level_err_rmse_ind[1:5, .x]*100,
                                             optim_mint_all_level_err_rmse_ind[1:5, .x]*100,
                                             Comb_Av_all_level_err_rmse_ind[1:5, .x]*100,
                                             mfts_all_level_err_rmse[1:5, .x]*100,
                                             BU_all_level_err_rmse_mfts[1:5, .x]*100,
                                             optim_all_level_err_rmse_mfts[1:5, .x]*100,
                                             optim_mint_all_level_err_rmse_mfts[1:5, .x]*100,
                                             Comb_Av_all_level_err_rmse_mfts[1:5, .x]*100),
                                  type = factor(paste(.y, "series"), 
                                  levels = c("Australia series", "Region series", "Area series", "Sex series", "Region + Sex series", "Area + Sex series"))))

ggplot(point_forecasts_rmse_series, aes(x = variety, y = series, fill = Method)) +
  geom_boxplot() + ylab("RMSFE (x 100)") + xlab("") + facet_wrap(~type) + scale_y_log10()
ggsave("point_forecasts_rmse.pdf", height = 8.27, width = 11.69, units = "in")


#################################################################

# ratio of each Area to the whole country
Area_total = Area_female = Area_male = array(NA, dim = c(18, 24, 47))
for(iw in 2:48)
{
  gettotal <- get(Area[iw])$rate$total
  gettotal[gettotal==0] <- NA
  getmale <- get(Area[iw])$rate$male
  getmale[getmale==0] <- NA
  getfemale <- get(Area[iw])$rate$female
  getfemale[getfemale==0] <- NA
  Area_total[,,iw-1]  = log(gettotal/Australia$rate$total)
  Area_female[,,iw-1] = log(getfemale/Australia$rate$female) 
  Area_male[,,iw-1]   = log(getmale/Australia$rate$male)
}  

# raw (averaged over age and Area)
total_rate_raw = apply(Area_total, c(1,3), mean, na.rm = TRUE)
female_rate_raw = apply(Area_female, c(1,3), mean, na.rm = TRUE)
male_rate_raw = apply(Area_male, c(1,3), mean, na.rm = TRUE)
colnames(male_rate_raw) = colnames(female_rate_raw) = colnames(total_rate_raw) = Area[2:48]

total_rate_raw_arrange = total_rate_raw
female_rate_raw_arrange = female_rate_raw
male_rate_raw_arrange = male_rate_raw


# raw (averaged over year and Area)
total_rate_raw_year  = apply(Area_total, c(2,3), mean, na.rm = TRUE)
female_rate_raw_year = apply(Area_female, c(2,3), mean, na.rm = TRUE)
male_rate_raw_year   = apply(Area_male, c(2,3), mean, na.rm = TRUE)
colnames(male_rate_raw_year) = colnames(female_rate_raw_year) = colnames(total_rate_raw_year) = Area[2:48]

total_rate_raw_arrange_year = total_rate_raw_year
female_rate_raw_arrange_year = female_rate_raw_year
male_rate_raw_arrange_year = male_rate_raw_year


# heatmap of log ratios; averaged over age
require(RColorBrewer)


# two graphs together
age_obs

savepdf("image_plot",width=20,height=14,toplines=.8,pointsize=12)
require(RColorBrewer)
par(mfrow=c(2,3))
image(age_obs, 1:47, total_rate_raw_arrange, col = brewer.pal(7, "RdBu"), xlab = "Age", ylab = "Area", main = "Total",zlim=c(-1,1))
box()
image(age_obs, 1:47, female_rate_raw_arrange, col = brewer.pal(7, "RdBu"), xlab = "Age", ylab = "", main = "Female",zlim=c(-1,1))
box()
image(age_obs, 1:47, male_rate_raw_arrange, col = brewer.pal(7, "RdBu"), xlab = "Age", ylab = "", main = "Male",zlim=c(-1,1))
box()

image(1993:2016, 1:47, total_rate_raw_arrange_year, col = brewer.pal(7, "RdBu"), xlab = "Year", ylab = "Area",zlim=c(-1,1))
box()
image(1993:2016, 1:47, female_rate_raw_arrange_year, col = brewer.pal(7, "RdBu"), xlab = "Year", ylab = "",zlim=c(-1,1))
box()
image(1993:2016, 1:47, male_rate_raw_arrange_year, col = brewer.pal(7, "RdBu"), xlab = "Year", ylab = "",zlim=c(-1,1))
box()
dev.off()









