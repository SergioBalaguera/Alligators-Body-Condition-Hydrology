###############################################################################X
###############   Alligator body condition and hydrology    ###################X
###############################################################################X
###############   Sergio A. Balaguera-Reina Oct 1st, 2024    ##################X
###############################################################################X

rm(list = ls())
dev.off()

## Libraries----
library(tidyverse)
library(scales)
library(MASS)
library(emmeans)
library(ggpmisc)
library (ncdf4)
library(reshape)
library(lubridate)
library(mgcv)
library(Boruta)
library(gratia) 
library(patchwork)

## 1. Body condition----
## 1.1 Call data----
Gators <- read.csv("Gators.csv")
Gators <- Gators %>%
  rename(Capture.Date = date) %>%
  mutate(Capture.Date = as.Date(Capture.Date, "%m/%d/%Y"))

## 1.2 Summarize data (create and export table 1)----
Gators %>%
  group_by(Sex) %>%
  summarise(
    n()
  )

Table1 <- Gators %>%
  group_by(Route) %>%
  summarise(
    `Number of alligators` = n(),
    Spring = length(which(Season %in% "Spring")),
    Fall = length(which(Season %in% "Fall")),
    Male = length(which(Sex %in% "Male")),
    Female = length(which(Sex %in% "Female")),
    Unknown = length(which(Sex %in% "U")),
    Subadult = length(which(SizeClass %in% "Subadult")),
    Adult = length(which(SizeClass %in% "Adult")),
    WY = paste(min(WY), max(WY), sep = "-")
  )

write.csv(Table1, "./Tables/Table1.csv")

sum(Table1$Subadult)
sum(Table1$Adult)

sum(Table1$Fall)
sum(Table1$Spring)

Appendix1 <- Gators %>%
  group_by(WY, Season, Route) %>%
  summarise(
    `Number of alligators` = n(),
    Male = length(which(Sex %in% "Male")),
    Female = length(which(Sex %in% "Female")),
    Unknown = length(which(Sex %in% "U")),
    Subadult = length(which(SizeClass %in% "Subadult")),
    Adult = length(which(SizeClass %in% "Adult")),
  )

Appendix1 %>%
  group_by(Route) %>%
  summarise(
    n = sum(`Number of alligators`)
  ) %>%
  arrange(n)

a <- Appendix1 %>%
  group_by(WY) %>%
  summarise(
    n = sum(`Number of alligators`)
  ) %>%
  arrange(n)

Appendix1 %>%
  group_by(Route, Season) %>%
  summarise(
    n = sum(`Number of alligators`)
  ) %>%
  subset(Season %in% "Fall") %>%
  arrange(n)

Appendix1 %>%
  group_by(Route, Season) %>%
  summarise(
    n = sum(`Number of alligators`)
  ) %>%
  subset(Season %in% "Spring") %>%
  arrange(n)

range(Gators$Total.Length)
mean(Gators$Total.Length)
range(Gators$SV.length)
mean(Gators$SV.length)
range(Gators$Weight)
mean(Gators$Weight)

## 1.3 Define allometric coefficients----
#Test for outliers
W_SVL <- lm(lnWeight ~ lnSVL, data = Gators)
summary(W_SVL)
cor.test(x = W_SVL$residuals, y = W_SVL$fitted.values, # Test for heteroscedasticity 
         method = "spearman") 
stud_resid <- studres(W_SVL) # identify outliers 
plot(stud_resid ~ Gators$lnSVL, ylab = 'Studentized Residuals', 
     xlab = 'Displacement')
abline(0,0)
Gators <- cbind(Gators, stud_resid)
Outliers <- subset(Gators, stud_resid < -3.0 | stud_resid > 3.0)
Outliers %>%
  group_by(Route) %>%
  summarise(
    n()
  )
Gators_NO <- subset(Gators, !stud_resid < -3.0 & !stud_resid > 3.0) # 57 outliers (larger than 3 and lower than -3)
W_SVL <- lm(lnWeight ~ lnSVL, data = Gators_NO)
summary(W_SVL)
cor.test(x = W_SVL$residuals, y = W_SVL$fitted.values, 
         method = "spearman")
coefficients(W_SVL)
confint(W_SVL, level=0.95)

#Assess effects data without outliers
Appendix2 <- data.frame(contrast = NULL,
                        estimate = NULL,
                        SE = NULL,
                        df = NULL,
                        t.ratio = NULL,
                        p.value = NULL)
#By Sex
Gators_NOS <- subset(Gators_NO, !Sex %in% "U")
a <- lm(lnWeight ~ lnSVL*Sex, data = Gators_NOS)
summary(a)
anova(a)
a <- lstrends(a, "Sex", var = "lnSVL")
a <- data.frame(pairs(a)) # Moderate evidence by sex
Appendix2 <- rbind(Appendix2, a)

#By Size Class
a <- lm(lnWeight ~ lnSVL*SizeClass, data = Gators_NO)
summary(a)
anova(a)
a <- lstrends(a, "SizeClass", var = "lnSVL")
a <- data.frame(pairs(a)) # Strong evidence by Size Class
Appendix2 <- rbind(Appendix2, a)

#By season
a <- lm(lnWeight ~ lnSVL*Season, data = Gators_NO)
summary(a)
anova(a)
a <- lstrends(a, "Season", var = "lnSVL")
a <- data.frame(pairs(a)) # No evidence of an effect by season
Appendix2 <- rbind(Appendix2, a)

#By Route
a <- lm(lnWeight ~ lnSVL*Route, data = Gators_NO)
summary(a)
anova(a)
confint(lstrends(a, "Route", var = "lnSVL"), adjust = "Bonferroni")
a <- test(pairs(lstrends(a, "Route", var = "lnSVL"), adjust = "Bonferroni"))
Appendix2 <- rbind(Appendix2, a) # Strong evidence of an effect by route

#By WY
a <- lm(lnWeight ~ lnSVL*as.factor(WY), data = Gators_NO)
summary(a)
anova(a)
confint(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")
a <- test(pairs(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")) 
Appendix2 <- rbind(Appendix2, a) # Strong evidence of an effect by route

Appendix2 <- Appendix2 %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3)))

#Remove most different routes
Gators_1 <- subset(Gators_NO, !Route %in% c("Lox_Marsh", "BICY_Marsh", "ENP-SS"))

#Assess effects data without most different routes (Central)
Appendix2.1 <- data.frame(contrast = NULL,
                          estimate = NULL,
                          SE = NULL,
                          df = NULL,
                          t.ratio = NULL,
                          p.value = NULL)
#By Sex
Gators_1S <- subset(Gators_1, !Sex %in% "U")
a <- lm(lnWeight ~ lnSVL*Sex, data = Gators_1S)
summary(a)
anova(a)
a <- lstrends(a, "Sex", var = "lnSVL")
a <- data.frame(pairs(a)) # No difference by sex
Appendix2.1 <- rbind(Appendix2.1, a)

#By Size Class
a <- lm(lnWeight ~ lnSVL*SizeClass, data = Gators_1)
summary(a)
anova(a)
a <- lstrends(a, "SizeClass", var = "lnSVL")
a <- data.frame(pairs(a)) # moderate evidence by Size Class

#By season
a <- lm(lnWeight ~ lnSVL*Season, data = Gators_1)
summary(a)
anova(a)
a <- lstrends(a, "Season", var = "lnSVL")
a <- data.frame(pairs(a)) # No differences by season
Appendix2.1 <- rbind(Appendix2.1, a)

#By Route
a <- lm(lnWeight ~ lnSVL*Route, data = Gators_1)
summary(a)
anova(a)
confint(lstrends(a, "Route", var = "lnSVL"), adjust = "Bonferroni")
a <- test(pairs(lstrends(a, "Route", var = "lnSVL"), adjust = "Bonferroni")) #Weak differences 
Appendix2.1 <- rbind(Appendix2.1, a)

#By WY
a <- lm(lnWeight ~ lnSVL*as.factor(WY), data = Gators_1)
summary(a)
anova(a)
b <- confint(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")
b$Group <- "Central"
Slope_Groups <- b
a <- test(pairs(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")) #Strong differences
Appendix2.1 <- rbind(Appendix2.1, a)

#No covariates
a <- lm(lnWeight ~ lnSVL, data = Gators_1)
summary(a)
coefficients(a)
confint(a)

Appendix2.1 <- Appendix2.1 %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3)))

## 1.4 Group sampling areas based on evidence----
#Group data based on slopes
Gators <- Gators_NO %>%
  mutate(Groups = case_when(Route %in% c("ENP-FC", "WCA3A-N41_Marsh", "WCA3B_Marsh", 
                                         "WCA3A-HD_Marsh", "WCA3A-TW_Marsh", "WCA2A_Marsh") ~ "Central",
                            Route %in% c("Lox_Marsh") ~ "Northeast",
                            Route %in% c("ENP-SS") ~ "South",
                            Route %in% c("BICY_Marsh") ~ "West"))

Groups <- split(Gators, ~ Gators$Groups)

#Assess effects data grouped
#South
Appendix2.2 <- data.frame(contrast = NULL,
                          estimate = NULL,
                          SE = NULL,
                          df = NULL,
                          t.ratio = NULL,
                          p.value = NULL)

#By Sex
Gators_S <- subset(Groups$South, !Sex %in% "U")
a <- lm(lnWeight ~ lnSVL*Sex, data = Gators_S)
summary(a)
anova(a)
a <- lstrends(a, "Sex", var = "lnSVL")
a <- data.frame(pairs(a)) # Weak difference by sex
Appendix2.2 <- rbind(Appendix2.2, a)

#By Size Class
a <- lm(lnWeight ~ lnSVL*SizeClass, data = Groups$South)
summary(a)
anova(a)
a <- lstrends(a, "SizeClass", var = "lnSVL")
a <- data.frame(pairs(a)) # No evidence by Size Class

#By season
a <- lm(lnWeight ~ lnSVL*Season, data = Groups$South)
summary(a)
anova(a)
a <- lstrends(a, "Season", var = "lnSVL")
a <- data.frame(pairs(a)) # No differences by season
Appendix2.2 <- rbind(Appendix2.2, a)

#By WY
a <- lm(lnWeight ~ lnSVL*as.factor(WY), data =Groups$South)
summary(a)
anova(a)
b <- confint(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")
b$Group <- "South"
Slope_Groups <- rbind(Slope_Groups, b)
a <- test(pairs(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")) 
Appendix2.2 <- rbind(Appendix2.2, a)

Appendix2.2 <- Appendix2.2 %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3)))

#No covariates
a <- lm(lnWeight ~ lnSVL, data = Groups$South)
summary(a)
coefficients(a)
confint(a)

#Northeast
Appendix2.3 <- data.frame(contrast = NULL,
                          estimate = NULL,
                          SE = NULL,
                          df = NULL,
                          t.ratio = NULL,
                          p.value = NULL)

#By Sex
Gators_S <- subset(Groups$Northeast, !Sex %in% "U")
a <- lm(lnWeight ~ lnSVL*Sex, data = Gators_S)
summary(a)
anova(a)
a <- lstrends(a, "Sex", var = "lnSVL")
a <- data.frame(pairs(a)) # No difference by sex
Appendix2.3 <- rbind(Appendix2.3, a)

#By Size Class
a <- lm(lnWeight ~ lnSVL*SizeClass, data = Groups$Northeast)
summary(a)
anova(a)
a <- lstrends(a, "SizeClass", var = "lnSVL")
a <- data.frame(pairs(a)) # No evidence by Size Class

#By season
a <- lm(lnWeight ~ lnSVL*Season, data = Groups$Northeast)
summary(a)
anova(a)
a <- lstrends(a, "Season", var = "lnSVL")
a <- data.frame(pairs(a)) # No differences by season
Appendix2.3 <- rbind(Appendix2.3, a)

#By WY
a <- lm(lnWeight ~ lnSVL*as.factor(WY), data =Groups$Northeast)
summary(a)
anova(a)
b <- confint(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")
b$Group <- "Northeast"
Slope_Groups <- rbind(Slope_Groups, b)
a <- test(pairs(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")) 
Appendix2.3 <- rbind(Appendix2.3, a)

Appendix2.3 <- Appendix2.3 %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3)))

#No covariates
a <- lm(lnWeight ~ lnSVL, data = Groups$Northeast)
summary(a)
coefficients(a)
confint(a)

#West
Appendix2.4 <- data.frame(contrast = NULL,
                          estimate = NULL,
                          SE = NULL,
                          df = NULL,
                          t.ratio = NULL,
                          p.value = NULL)

#By Sex
Gators_S <- subset(Groups$West, !Sex %in% "U")
a <- lm(lnWeight ~ lnSVL*Sex, data = Gators_S)
summary(a)
anova(a)
a <- lstrends(a, "Sex", var = "lnSVL")
a <- data.frame(pairs(a)) # No difference by sex
Appendix2.4 <- rbind(Appendix2.4, a)

#By Size Class
a <- lm(lnWeight ~ lnSVL*SizeClass, data = Groups$West)
summary(a)
anova(a)
a <- lstrends(a, "SizeClass", var = "lnSVL")
a <- data.frame(pairs(a)) # No evidence by Size Class

#By season
a <- lm(lnWeight ~ lnSVL*Season, data = Groups$West)
summary(a)
anova(a)
a <- lstrends(a, "Season", var = "lnSVL")
a <- data.frame(pairs(a)) # Weak differences by season
Appendix2.4 <- rbind(Appendix2.4, a)

#By WY
a <- lm(lnWeight ~ lnSVL*as.factor(WY), data =Groups$West)
summary(a)
anova(a)
b <- confint(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")
b$Group <- "West"
Slope_Groups <- rbind(Slope_Groups, b)
a <- test(pairs(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")) 
Appendix2.4 <- rbind(Appendix2.4, a)

#No covariates
a <- lm(lnWeight ~ lnSVL, data = Groups$West)
summary(a)
coefficients(a)
confint(a)

Appendix2.4 <- Appendix2.4 %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3)))

## 1.5 Export the supplementary material table and produce Figure 1----
SupMat <- list(Appendix1, Appendix2, Appendix2.1, Appendix2.2, Appendix2.3, Appendix2.4)
names(SupMat) <- c("Summary Gators", "All data no outliers", "Central", "South", "Northeast", "West")
writexl::write_xlsx(SupMat, "./Supplementary material/SupMat.xlsx")

#Figure 2
ggplot(Gators, aes(x = lnSVL, y = lnWeight, color = Groups)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_point(alpha = 0.05) +
  stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = T, coef.digits = 3, f.digits = 3, p.digits = 3) +
  annotate(geom = "text", x = 4.6, y = 5.05, label = "CI = 3.14 - 3.19", col = "#F8766D") +
  annotate(geom = "text", x = 4.6, y = 4.84, label = "CI = 2.86 - 2.95", col = "#7CAE00") +
  annotate(geom = "text", x = 4.6, y = 4.63, label = "CI = 3.00 - 3.12", col = "#00BFC4") +
  annotate(geom = "text", x = 4.6, y = 4.42, label = "CI = 3.19 - 3.33", col = "#C77CFF") +
  labs(y = "ln Weight", x = "ln Snout-Vent Length") +
  theme_bw() +
  theme(legend.position = c(0.9, 0.1))

ggsave("./Figures/Figure2.jpeg", plot = last_plot(), height = 3300, width = 2550, dpi = 370, units = "px")

## 2. Body condition ~ hydrology analysis ####
#### 2.1 Get the hydro variables ####
#Set up gator Data
cgators <- data.frame(Consecutive = as.integer(paste(Gators$Consecutive)),
                      Capture.Date = paste(Gators$Capture.Date),
                      Route = paste(Gators$Route),
                      UTM_Easting = as.numeric(paste(Gators$UTM_Easting)),
                      UTM_Northing = as.numeric(paste(Gators$UTM_Northing)))
cgators <- cgators[complete.cases(cgators$UTM_Easting),]

cgators$Capture.Date <- as.Date(cgators$Capture.Date)
cgators <- cgators[order(cgators$Consecutive, decreasing = F),]
cdate <- sort(unique(cgators$Capture.Date))

# Set up DEM file
if (!file.exists("./input/eden_dem_cm_oc11.nc")) {
  err <- try (download.file("https://sofia.usgs.gov/eden/data/dem/eden_dem_cm_oc11.zip", 
                            "./input/eden_dem_cm_oc11.zip"))
  unzip("./input/eden_dem_cm_oc11.zip", "eden_dem_cm_oc11.nc", exdir = "./input")
}

dem.nc <- nc_open("./input/eden_dem_cm_oc11.nc")
x <- ncvar_get(dem.nc, "x")
y <- ncvar_get(dem.nc, "y")
dem <- ncvar_get(dem.nc, "dem")
nc_close(dem.nc)

# Set up WL files
st_date <- seq(cdate[1], length = 2, by = "-3 years")[2] + 1
en_date <- rev(cdate)[1]
st_yr <- as.POSIXlt(st_date)$year + 1900
en_yr <- as.POSIXlt(en_date)$year + 1900
st_qtr <- substr(quarters(st_date), 2, 2)
en_qtr <- substr(quarters(en_date), 2, 2)

surf_files <- tdate <- NULL
options(timeout = max(300, getOption("timeout")))
for (i in st_yr:en_yr) {
  first_q <- ifelse (i == st_yr, st_qtr, 1)
  last_q <- ifelse (i == en_yr, en_qtr, 4)
  for (j in first_q:last_q) {
    f <- paste0("./surfaces/", i, "_q", j, ".nc")
    if (!file.exists(f))
      download.file(paste0("https://sflthredds.er.usgs.gov/thredds/fileServer/eden/", f), mode = "wb", f)
    surf_files <- c(surf_files, f)
    t.nc <- nc_open(f)
    d <- ncvar_get(t.nc, "time")
    tdate <- as.Date(c(tdate, as.Date(t.nc$dim$time$units, format = "days since %Y-%m-%dT%H:%M:%S") + 
                         d, recursive = T), origin = "1970/1/1")
    nc_close(t.nc)
  }
}

#Calculate water depth
depth_df <- data.frame(date = tdate)
system.time(
  for (i in 1:dim(cgators)[1]) {
    if (i %% 100 == 0) print(i)
    # Find EDEN pixel for each capture
    s <- cgators[i, c("UTM_Easting", "UTM_Northing")]
    s$grid_x <- which.min(abs(x - s$UTM_Easting))
    s$grid_y <- which.min(abs(y - s$UTM_Northing))
    stage <- NULL
    for (k in 1:length(surf_files)) {
      t.nc <- nc_open(surf_files[k])
      t <- ncvar_get(t.nc, "stage", c(s$grid_x, s$grid_y, 1), c(1, 1, -1))
      stage <- c(stage, t)
      nc_close(t.nc)
    }
    depth_df <- cbind(depth_df, stage - dem[s$grid_x, s$grid_y])
    names(depth_df)[ncol(depth_df)] <- cgators$Consecutive[i]
  }
)

#Calculate day since dry
wet <- ifelse(depth_df[, 2:ncol(depth_df)] >= 15, 1, 0) # wet == 1
dsd <- wet
for (i in 2:nrow(wet))
  dsd[i, ] <- (dsd[i - 1, ] + wet[i, ]) * wet[i, ] # resets if 0
wet <- data.frame(date = tdate, wet)
dsd <- data.frame(date = tdate, dsd)
for (i in 2:ncol(dsd)) {
  d <- which(dsd[, i] == 0)[1] - 1
  if (is.na(d)) d <- dim(dsd)[1] 
  if (d != 0) dsd[1:d, i] <- Inf
}

#Calculate median, max, min, and ranges
min_30 <- 
  med_30 <- 
  max_30 <- 
  range_30 <- 
  min_90 <- 
  med_90 <- 
  max_90 <- 
  range_90 <- 
  min_180 <- 
  med_180 <- 
  max_180 <- 
  range_180 <-
  min_365 <- 
  med_365 <- 
  max_365 <- 
  range_365 <- depth_df

min_30 <- as.matrix(min_30[min_30$date %in% cdate, 2:ncol(min_30)])
med_30 <- as.matrix(med_30[med_30$date %in% cdate, 2:ncol(med_30)])
max_30 <- as.matrix(max_30[max_30$date %in% cdate, 2:ncol(max_30)])
range_30 <- as.matrix(range_30[range_30$date %in% cdate, 2:ncol(range_30)])
min_90 <- as.matrix(min_90[min_90$date %in% cdate, 2:ncol(min_90)])
med_90 <- as.matrix(med_90[med_90$date %in% cdate, 2:ncol(med_90)])
max_90 <- as.matrix(max_90[max_90$date %in% cdate, 2:ncol(max_90)])
range_90 <- as.matrix(range_90[range_90$date %in% cdate, 2:ncol(range_90)])
min_180 <- as.matrix(min_180[min_180$date %in% cdate, 2:ncol(min_180)])
med_180 <- as.matrix(med_180[med_180$date %in% cdate, 2:ncol(med_180)])
max_180 <- as.matrix(max_180[max_180$date %in% cdate, 2:ncol(max_180)])
range_180 <- as.matrix(range_180[range_180$date %in% cdate, 2:ncol(range_180)])
min_365 <- as.matrix(min_365[min_365$date %in% cdate, 2:ncol(min_365)])
med_365 <- as.matrix(med_365[med_365$date %in% cdate, 2:ncol(med_365)])
max_365 <- as.matrix(max_365[max_365$date %in% cdate, 2:ncol(max_365)])
range_365 <- as.matrix(range_365[range_365$date %in% cdate, 2:ncol(range_365)])
for (i in 1:length(cdate)) {
  print(i)
  j <- which(depth_df$date == cdate[i])
  min_30[i, ] <- apply(depth_df[(j - 29):j, 2:ncol(depth_df)], 2, min)
  med_30[i, ] <- apply(depth_df[(j - 29):j, 2:ncol(depth_df)], 2, median)
  max_30[i, ] <- apply(depth_df[(j - 29):j, 2:ncol(depth_df)], 2, max)
  min_90[i, ] <- apply(depth_df[(j - 89):j, 2:ncol(depth_df)], 2, min)
  med_90[i, ] <- apply(depth_df[(j - 89):j, 2:ncol(depth_df)], 2, median)
  max_90[i, ] <- apply(depth_df[(j - 89):j, 2:ncol(depth_df)], 2, max)
  min_180[i, ] <- apply(depth_df[(j - 179):j, 2:ncol(depth_df)], 2, min)
  med_180[i, ] <- apply(depth_df[(j - 179):j, 2:ncol(depth_df)], 2, median)
  max_180[i, ] <- apply(depth_df[(j - 179):j, 2:ncol(depth_df)], 2, max)
  min_365[i, ] <- apply(depth_df[(j - 364):j, 2:ncol(depth_df)], 2, min)
  med_365[i, ] <- apply(depth_df[(j - 364):j, 2:ncol(depth_df)], 2, median)
  max_365[i, ] <- apply(depth_df[(j - 364):j, 2:ncol(depth_df)], 2, max)
}

min_30 <- data.frame(cdate, min_30)
med_30 <- data.frame(cdate, med_30)
max_30 <- data.frame(cdate, max_30)
range_30 <- max_30 - min_30
range_30 <- data.frame(cdate, range_30)

min_90 <- data.frame(cdate, min_90)
med_90 <- data.frame(cdate, med_90)
max_90 <- data.frame(cdate, max_90)
range_90 <- max_90 - min_90
range_90 <- data.frame(cdate, range_90)

min_180 <- data.frame(cdate, min_180)
med_180 <- data.frame(cdate, med_180)
max_180 <- data.frame(cdate, max_180)
range_180 <- max_180 - min_180
range_180 <- data.frame(cdate, range_180)

min_365 <- data.frame(cdate, min_365)
med_365 <- data.frame(cdate, med_365)
max_365 <- data.frame(cdate, max_365)
range_365 <- max_365 - min_365
range_365 <- data.frame(cdate, range_365)

#Calculate average hydroperiod
ahp <- data.frame(year = 2000:2024)
for (i in 2000:2024) {
  st <- which(tdate == paste0(i - 1, "-05-01"))
  en <- which(tdate == paste0(i, "-04-30"))
  ahp[i - 1999, 2:ncol(wet)] <- colSums(wet[st:en, 2:ncol(wet)])
}
names(ahp) <- names(wet)

# Restrict to dates
depth_df <- depth_df[depth_df$date %in% cdate, ]
dsd <- dsd[dsd$date %in% cdate, ]

#Convert tables for analysis
depth_df2 <- data.frame(capture.date = NULL, consecutive = NULL, depth = NULL)
for (i in 1:dim(cgators)[1]) {
  d <- which(names(depth_df) == paste0(cgators$Consecutive[i]))
  c <- which(depth_df$date == cgators$Capture.Date[i])
  if (length(d))
    depth_df2 <- rbind(depth_df2, data.frame(capture.date = cgators$Capture.Date[i], 
                                             consecutive = cgators$Consecutive[i], depth = depth_df[c, d]))
}
names(depth_df2) <- c("date", "consecutive", "depth")
head(depth_df2)

dsd2 <- data.frame(capture.date = NULL, consecutive = NULL, dsd = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(dsd$date == cgators$Capture.Date[i])
  d <- which(names(dsd) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    dsd2 <- rbind(dsd2, data.frame(capture.date = cgators$Capture.Date[i], consecutive = cgators$Consecutive[i], dsd = dsd[c, d]))
}
names(dsd2) <- c("date", "consecutive", "dsd")
head(dsd2)

min_30_2 <- data.frame(capture.date = NULL, consecutive = NULL, min_30 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(min_30$cdate == cgators$Capture.Date[i])
  d <- which(names(min_30) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    min_30_2 <- rbind(min_30_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                           consecutive = cgators$Consecutive[i], 
                                           min_30 = min_30[c, d]))
}
names(min_30_2) <- c("date", "consecutive", "min_30")
head(min_30_2)

med_30_2 <- data.frame(capture.date = NULL, consecutive = NULL, med_30 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(med_30$cdate == cgators$Capture.Date[i])
  d <- which(names(med_30) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    med_30_2 <- rbind(med_30_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                           consecutive = cgators$Consecutive[i], 
                                           med_30 = med_30[c, d]))
}
names(med_30_2) <- c("date", "consecutive", "med_30")
head(med_30_2)

max_30_2 <- data.frame(capture.date = NULL, consecutive = NULL, max_30 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(max_30$cdate == cgators$Capture.Date[i])
  d <- which(names(max_30) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    max_30_2 <- rbind(max_30_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                           consecutive = cgators$Consecutive[i], 
                                           max_30 = max_30[c, d]))
}
names(max_30_2) <- c("date", "consecutive", "max_30")
head(max_30_2)

range_30_2 <- data.frame(capture.date = NULL, consecutive = NULL, range_30 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(range_30$cdate == cgators$Capture.Date[i])
  d <- which(names(range_30) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    range_30_2 <- rbind(range_30_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                               consecutive = cgators$Consecutive[i], 
                                               range_30 = range_30[c, d]))
}
names(range_30_2) <- c("date", "consecutive", "range_30")
head(range_30_2)

min_90_2 <- data.frame(capture.date = NULL, consecutive = NULL, min_90 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(min_90$cdate == cgators$Capture.Date[i])
  d <- which(names(min_90) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    min_90_2 <- rbind(min_90_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                           consecutive = cgators$Consecutive[i], 
                                           min_90 = min_90[c, d]))
}
names(min_90_2) <- c("date", "consecutive", "min_90")
head(min_90_2)

med_90_2 <- data.frame(capture.date = NULL, consecutive = NULL, med_90 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(med_90$cdate == cgators$Capture.Date[i])
  d <- which(names(med_90) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    med_90_2 <- rbind(med_90_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                           consecutive = cgators$Consecutive[i], 
                                           med_90 = med_90[c, d]))
}
names(med_90_2) <- c("date", "consecutive", "med_90")
head(med_90_2)

max_90_2 <- data.frame(capture.date = NULL, consecutive = NULL, max_90 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(max_90$cdate == cgators$Capture.Date[i])
  d <- which(names(max_90) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    max_90_2 <- rbind(max_90_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                           consecutive = cgators$Consecutive[i], 
                                           max_90 = max_90[c, d]))
}
names(max_90_2) <- c("date", "consecutive", "max_90")
head(max_90_2)

range_90_2 <- data.frame(capture.date = NULL, consecutive = NULL, range_90 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(range_90$cdate == cgators$Capture.Date[i])
  d <- which(names(range_90) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    range_90_2 <- rbind(range_90_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                               consecutive = cgators$Consecutive[i], 
                                               range_90 = range_90[c, d]))
}
names(range_90_2) <- c("date", "consecutive", "range_90")
head(range_90_2)

min_180_2 <- data.frame(capture.date = NULL, consecutive = NULL, min_180 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(min_180$cdate == cgators$Capture.Date[i])
  d <- which(names(min_180) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    min_180_2 <- rbind(min_180_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                             consecutive = cgators$Consecutive[i], 
                                             min_180 = min_180[c, d]))
}
names(min_180_2) <- c("date", "consecutive", "min_180")
head(min_180_2)

med_180_2 <- data.frame(capture.date = NULL, consecutive = NULL, med_180 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(med_180$cdate == cgators$Capture.Date[i])
  d <- which(names(med_180) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    med_180_2 <- rbind(med_180_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                             consecutive = cgators$Consecutive[i], 
                                             med_180 = med_180[c, d]))
}
names(med_180_2) <- c("date", "consecutive", "med_180")
head(med_180_2)

max_180_2 <- data.frame(capture.date = NULL, consecutive = NULL, max_180 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(max_180$cdate == cgators$Capture.Date[i])
  d <- which(names(max_180) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    max_180_2 <- rbind(max_180_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                             consecutive = cgators$Consecutive[i], 
                                             max_180 = max_180[c, d]))
}
names(max_180_2) <- c("date", "consecutive", "max_180")
head(max_180_2)

range_180_2 <- data.frame(capture.date = NULL, consecutive = NULL, range_180 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(range_180$cdate == cgators$Capture.Date[i])
  d <- which(names(range_180) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    range_180_2 <- rbind(range_180_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                                 consecutive = cgators$Consecutive[i], 
                                                 range_180 = range_180[c, d]))
}
names(range_180_2) <- c("date", "consecutive", "range_180")
head(range_180_2)

min_365_2 <- data.frame(capture.date = NULL, consecutive = NULL, min_365 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(min_365$cdate == cgators$Capture.Date[i])
  d <- which(names(min_365) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    min_365_2 <- rbind(min_365_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                             consecutive = cgators$Consecutive[i], 
                                             min_365 = min_365[c, d]))
}
names(min_365_2) <- c("date", "consecutive", "min_365")
head(min_365_2)

med_365_2 <- data.frame(capture.date = NULL, consecutive = NULL, med_365 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(med_365$cdate == cgators$Capture.Date[i])
  d <- which(names(med_365) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    med_365_2 <- rbind(med_365_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                             consecutive = cgators$Consecutive[i], 
                                             med_365 = med_365[c, d]))
}
names(med_365_2) <- c("date", "consecutive", "med_365")
head(med_365_2)

max_365_2 <- data.frame(capture.date = NULL, consecutive = NULL, max_365 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(max_365$cdate == cgators$Capture.Date[i])
  d <- which(names(max_365) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    max_365_2 <- rbind(max_365_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                             consecutive = cgators$Consecutive[i], 
                                             max_365 = max_365[c, d]))
}
names(max_365_2) <- c("date", "consecutive", "max_365")
head(max_365_2)

range_365_2 <- data.frame(capture.date = NULL, consecutive = NULL, range_365 = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(range_365$cdate == cgators$Capture.Date[i])
  d <- which(names(range_365) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    range_365_2 <- rbind(range_365_2, data.frame(capture.date = cgators$Capture.Date[i], 
                                                 consecutive = cgators$Consecutive[i], 
                                                 range_365 = range_365[c, d]))
}
names(range_365_2) <- c("date", "consecutive", "range_365")
head(range_365_2)

cgators <- cgators %>%
  mutate(year = as.integer(format(Capture.Date, "%Y"))) %>%
  mutate(month = as.integer(format(Capture.Date, "%m")))

cgators$wy <- with(cgators, 
                   ifelse(month >= 6 & month <= 12, year + 1, year))

ahp_2 <- data.frame(wy = NULL, consecutive = NULL, ahp = NULL)
for (i in 1:dim(cgators)[1]) {
  c <- which(ahp$date == cgators$wy[i])
  d <- which(names(ahp) == paste0("X", cgators$Consecutive[i]))
  if (length(d))
    ahp_2 <- rbind(ahp_2, data.frame(wy = cgators$wy[i], 
                                     consecutive = cgators$Consecutive[i], 
                                     ahp = ahp[c, d]))
}
names(ahp_2) <- c("wy", "consecutive", "ahp")
head(ahp_2)

#Get on spot water depth
onSpot <- data.frame(date = as.Date(Gators$Capture.Date),
                     consecutive = as.integer(Gators$Consecutive),
                     OS_water_depth = as.numeric(Gators$Water.Depth),
                     Route = Gators$Route)
summary(onSpot)
onSpot$OS_water_depth[onSpot$OS_water_depth == 999.0] = NA

#Consolidate variables in one table
covariates <- plyr::join_all(list(onSpot, depth_df2, min_30_2, med_30_2, max_30_2,
                                  range_30_2, min_90_2, med_90_2, max_90_2, range_90_2, 
                                  min_180_2, med_180_2, max_180_2, range_180_2, min_365_2, 
                                  med_365_2, max_365_2, range_365_2), 
                             by = c("date", "consecutive"),
                             type = "left")
covariates <- covariates %>%
  mutate(year = as.integer(format(date, "%Y"))) %>%
  mutate(month = as.integer(format(date, "%m")))

covariates$wy <- with(covariates,
                      ifelse(month >=  6 & month <= 12, year + 1, year))

covariates <- left_join(covariates, ahp_2, 
                        by = c("wy", "consecutive"))

covariates <- covariates %>%
  dplyr::select(-wy, -year, -month) %>%
  relocate(Route, .before = OS_water_depth)

#### 2.2 merge body condition with covariates and run models ####
covariates <- covariates %>%
  mutate(date = as.Date(date))

Gators <- dplyr::rename(Gators,
                        "date" = "Capture.Date",
                        "consecutive" = "Consecutive")

Gators_Covariates <- left_join(Gators, covariates, 
                               by = c("date", "consecutive", "Route"))

Gators_Covariates <- Gators_Covariates %>%
  dplyr::select(consecutive, date, Route, UTM_Easting, UTM_Northing, SV.length,
                Weight, Sex, lnWeight, lnSVL, 
                SizeClass, month, Season, year, WY, stud_resid, Groups, OS_water_depth,
                depth, min_30, med_30, max_30, range_30, min_90, med_90, max_90,
                range_90, min_180, med_180, max_180, range_180, min_365, med_365,
                max_365, range_365, ahp) %>%
  mutate(Relativek = case_when(Groups %in% "Central" ~ (Weight/SV.length^3.17)*10^5,
                               Groups %in% "Northeast" ~ (Weight/SV.length^2.83)*10^5,
                               Groups %in% "South" ~ (Weight/SV.length^3.03)*10^5,
                               Groups %in% "West" ~ (Weight/SV.length^3.27)*10^5)) %>%
  mutate(Groups = as.factor(Groups)) %>%
  mutate(WY = as.factor(WY)) %>%
  mutate(Season = as.factor(Season)) %>%
  mutate(Season_num = case_when(Season %in% "Fall" ~ 1,
                                Season %in% "Spring" ~ 2)) %>%
  mutate(Sex = as.factor(Sex)) %>%
  mutate(month = as.factor(month)) %>%
  mutate(SizeClass = as.factor(SizeClass)) %>%
  mutate(Route = as.factor(Route)) %>%
  mutate(Julian = yday(date))

summary(Gators_Covariates) 

#Body condition Central
Gators_Covariates %>%
  subset(Groups %in% "Central") %>%
  group_by(Route) %>%
  summarise(
    mean(Relativek),
    sd(Relativek)
  )

## 3. GAM----
## 3.1 Mean relative K per route, season, year----
Gators_mean <- Gators_Covariates %>%
  group_by(Route, Season, WY) %>%
  summarise(across(c(Relativek, depth, min_30, med_30, max_30, range_30,
                     min_90, med_90, max_90, range_90, min_180, med_180, 
                     max_180, range_180, min_365, med_365, max_365, range_365, 
                     ahp), mean), .groups = 'drop') %>%
  mutate(Groups = case_when(Route %in% c("ENP-FC", "WCA3A-N41_Marsh", "WCA3B_Marsh", 
                                         "WCA3A-HD_Marsh", "WCA3A-TW_Marsh", "WCA2A_Marsh") ~ "Central",
                            Route %in% c("Lox_Marsh") ~ "Northeast",
                            Route %in% c("ENP-SS") ~ "South",
                            Route %in% c("BICY_Marsh") ~ "West")) %>%
  as.data.frame()

## 3.1.1 Variable selection----
set.seed(123)
Boruta_Analysis_mean_Central <- Gators_mean %>%
  subset(Groups %in% "Central") %>%
  dplyr::select(4:22) %>%
  drop_na() %>%
  Boruta(Relativek ~ ., data = ., doTrace = 2, maxRuns = 1000)
plot(Boruta_Analysis_mean_Central)
attStats(Boruta_Analysis_mean_Central) %>%
  as.data.frame() %>%
  rownames_to_column('var') %>%
  dplyr::select(var, decision, meanImp) %>%
  arrange(desc(meanImp))

set.seed(123)
Boruta_Analysis_mean_Northeast <- Gators_mean %>%
  subset(Groups %in% "Northeast") %>%
  dplyr::select(4:22) %>%
  drop_na() %>%
  Boruta(Relativek ~ ., data = ., doTrace = 2, maxRuns = 1000)
plot(Boruta_Analysis_mean_Northeast)
attStats(Boruta_Analysis_mean_Northeast) %>%
  as.data.frame() %>%
  rownames_to_column('var') %>%
  dplyr::select(var, decision, meanImp) %>%
  arrange(desc(meanImp)) 

set.seed(123)
Boruta_Analysis_mean_South <- Gators_mean %>%
  subset(Groups %in% "South") %>%
  dplyr::select(4:22) %>%
  drop_na() %>%
  Boruta(Relativek ~ ., data = ., doTrace = 2, maxRuns = 1000)
plot(Boruta_Analysis_mean_South)
attStats(Boruta_Analysis_mean_South) %>%
  as.data.frame() %>%
  rownames_to_column('var') %>%
  dplyr::select(var, decision, meanImp) %>%
  arrange(desc(meanImp))

set.seed(123)
Boruta_Analysis_mean_West <- Gators_mean %>%
  subset(Groups %in% "West") %>%
  dplyr::select(4:22) %>%
  drop_na() %>%
  Boruta(Relativek ~ ., data = ., doTrace = 2, maxRuns = 1000)
plot(Boruta_Analysis_mean_West)
attStats(Boruta_Analysis_mean_West) %>%
  as.data.frame() %>%
  rownames_to_column('var') %>%
  dplyr::select(var, decision, meanImp) %>%
  arrange(desc(meanImp))

## 3.1.2 Models----
## 3.1.2.1 Central----
#Gaussian
system.time(
  m1 <- bam(Relativek ~ #model chosen
              s(max_365, Route, bs = 'fs') + #Model using only boruta confirmed selection.
              s(range_365, Route, bs = 'fs') +
              s(med_365, Route, bs = 'fs') +
              s(ahp, Route, bs = 'fs') +
              s(max_180, Route, bs = 'fs') +
              s(med_180, Route, bs = 'fs') +
              s(range_30, Route, bs = 'fs') +
              s(min_365, Route, bs = 'fs') +
              s(max_30, Route, bs = 'fs') +
              s(max_90, Route, bs = 'fs') +
              s(min_30, Route, bs = 'fs') +
              s(med_30, Route, bs = 'fs') +
              s(depth, Route, bs = 'fs') +
              s(min_90, Route, bs = 'fs') +
              s(med_90, Route, bs = 'fs') +
              s(range_180, Route, bs = 'fs') +
              s(range_90, Route, bs = 'fs') + 
              s(min_180, Route, bs= 'fs') + 
              ti(max_365, range_365, Route, bs = 'fs') +
              ti(max_365, med_365, Route, bs = 'fs') +
              ti(max_365, ahp, Route, bs = 'fs') + 
              ti(max_365, max_180, Route, bs = 'fs') + 
              ti(max_365, med_180, Route, bs = 'fs') + 
              ti(max_365, range_30, Route, bs = 'fs') + 
              ti(max_365, min_365, Route, bs = 'fs') + 
              ti(max_365, max_30, Route, bs = 'fs') + 
              ti(max_365, max_90, Route, bs = 'fs') + 
              ti(max_365, min_30, Route, bs = 'fs') + 
              ti(max_365, med_30, Route, bs = 'fs') + 
              ti(max_365, depth, Route, bs = 'fs') + 
              ti(max_365, min_90, Route, bs = 'fs') + 
              ti(max_365, med_90, Route, bs = 'fs') + 
              ti(max_365, range_180, Route, bs = 'fs') +
              ti(max_365, range_90, Route, bs = 'fs') + 
              ti(max_365, min_180, Route, bs = 'fs') + 
              s(WY, bs = 're'),
            data = subset(Gators_mean, Groups %in% "Central"), method = 'fREML',
            nthreads = c(4,1), select = T, discrete = T)
)

m1_sum <- summary(m1) # Deviance explained 71.8%
k.check(m1) %>% # Ensure sufficient complexity was allocated to splines and check significant predictors  
  as.data.frame() %>%
  arrange(`p-value`) %>%
  filter(`p-value` < 0.1) %>% 
  rownames_to_column('smooth') %>% 
  dplyr::select(c(1:5)) %>% 
  mutate(edf = round(edf, 1))
appraise(m1, method = 'simulate') # It's not perfect but looks good

#Gamma - inverse 
system.time(
  m2 <- bam(Relativek ~
              s(max_365, Route, bs = 'fs') +
              s(range_365, Route, bs = 'fs') +
              s(med_365, Route, bs = 'fs') +
              s(ahp, Route, bs = 'fs') +
              s(max_180, Route, bs = 'fs') +
              s(med_180, Route, bs = 'fs') +
              s(range_30, Route, bs = 'fs') +
              s(min_365, Route, bs = 'fs') +
              s(max_30, Route, bs = 'fs') +
              s(max_90, Route, bs = 'fs') +
              s(min_30, Route, bs = 'fs') +
              s(med_30, Route, bs = 'fs') +
              s(depth, Route, bs = 'fs') +
              s(min_90, Route, bs = 'fs') +
              s(med_90, Route, bs = 'fs') +
              s(range_180, Route, bs = 'fs') +
              s(range_90, Route, bs = 'fs') + 
              s(min_180, Route, bs= 'fs') + 
              ti(max_365, range_365, Route, bs = 'fs') +
              ti(max_365, med_365, Route, bs = 'fs') +
              ti(max_365, ahp, Route, bs = 'fs') + 
              ti(max_365, max_180, Route, bs = 'fs') + 
              ti(max_365, med_180, Route, bs = 'fs') + 
              ti(max_365, range_30, Route, bs = 'fs') + 
              ti(max_365, min_365, Route, bs = 'fs') + 
              ti(max_365, max_30, Route, bs = 'fs') + 
              ti(max_365, max_90, Route, bs = 'fs') + 
              ti(max_365, min_30, Route, bs = 'fs') + 
              ti(max_365, med_30, Route, bs = 'fs') + 
              ti(max_365, depth, Route, bs = 'fs') + 
              ti(max_365, min_90, Route, bs = 'fs') + 
              ti(max_365, med_90, Route, bs = 'fs') + 
              ti(max_365, range_180, Route, bs = 'fs') +
              ti(max_365, range_90, Route, bs = 'fs') + 
              ti(max_365, min_180, Route, bs = 'fs') + 
              s(WY, bs = 're'),
            data = subset(Gators_mean, Groups %in% "Central"), method = 'fREML',
            nthreads = c(4,1), family = Gamma, select = T, discrete = T)
)

m2_sum <- summary(m2) # Deviance explained 76.5%
k.check(m2) %>% # Ensure sufficient complexity was allocated to splines and check significant predictors  
  as.data.frame() %>%
  arrange(`p-value`) %>%
  filter(`p-value` < 0.1) %>% 
  rownames_to_column('smooth') %>% 
  dplyr::select(c(1:3)) %>% 
  mutate(edf = round(edf, 1))
appraise(m2, method = 'simulate') # No good looking

#Gamma log
system.time(
  m3 <- bam(Relativek ~
              s(max_365, Route, bs = 'fs') +
              s(range_365, Route, bs = 'fs') +
              s(med_365, Route, bs = 'fs') +
              s(ahp, Route, bs = 'fs') +
              s(max_180, Route, bs = 'fs') +
              s(med_180, Route, bs = 'fs') +
              s(range_30, Route, bs = 'fs') +
              s(min_365, Route, bs = 'fs') +
              s(max_30, Route, bs = 'fs') +
              s(max_90, Route, bs = 'fs') +
              s(min_30, Route, bs = 'fs') +
              s(med_30, Route, bs = 'fs') +
              s(depth, Route, bs = 'fs') +
              s(min_90, Route, bs = 'fs') +
              s(med_90, Route, bs = 'fs') +
              s(range_180, Route, bs = 'fs') +
              s(range_90, Route, bs = 'fs') + 
              s(min_180, Route, bs= 'fs') + 
              ti(max_365, range_365, Route, bs = 'fs') +
              ti(max_365, med_365, Route, bs = 'fs') +
              ti(max_365, ahp, Route, bs = 'fs') + 
              ti(max_365, max_180, Route, bs = 'fs') + 
              ti(max_365, med_180, Route, bs = 'fs') + 
              ti(max_365, range_30, Route, bs = 'fs') + 
              ti(max_365, min_365, Route, bs = 'fs') + 
              ti(max_365, max_30, Route, bs = 'fs') + 
              ti(max_365, max_90, Route, bs = 'fs') + 
              ti(max_365, min_30, Route, bs = 'fs') + 
              ti(max_365, med_30, Route, bs = 'fs') + 
              ti(max_365, depth, Route, bs = 'fs') + 
              ti(max_365, min_90, Route, bs = 'fs') + 
              ti(max_365, med_90, Route, bs = 'fs') + 
              ti(max_365, range_180, Route, bs = 'fs') +
              ti(max_365, range_90, Route, bs = 'fs') + 
              ti(max_365, min_180, Route, bs = 'fs') + 
              s(WY, bs = 're'),
            data = subset(Gators_mean, Groups %in% "Central"), method = 'fREML',
            nthreads = c(4,1), family = Gamma(link = log), select = T, discrete = T)
)

m3_sum <- summary(m3) # Deviance explained 64.8%
k.check(m3) %>% # Ensure sufficient complexity was allocated to splines and check significant predictors  
  as.data.frame() %>%
  arrange(`p-value`) %>%
  filter(`p-value` < 0.1) %>% 
  rownames_to_column('smooth') %>% 
  dplyr::select(c(1:3)) %>% 
  mutate(edf = round(edf, 1))
appraise(m3, method = 'simulate') #No good looking

#Test models
AIC(m1, m2, m3) # m1 has the lowest AIC
anova(m1, m2, m3, test = 'F') # no difference, choose either

draw(m1, select = c("s(range_365,Route)",
                    "s(ahp,Route)",
                    "s(min_365,Route)",
                    "s(max_30,Route)",
                    "s(max_90,Route)",
                    "s(med_30,Route)",
                    "s(depth,Route)",
                    "s(min_90,Route)",
                    "s(med_90,Route)"), scales = 'free') +
  plot_layout(guides = "collect") & theme_bw()

ggsave("./Figures/GAM_Univariate_Central_Mean.jpeg", plot = last_plot(), height = 1500, width = 2500, 
       units = "px", dpi = 250)


draw(m1, ncol = 1, select = c("ti(max_365,med_180,Route)",
                              "ti(max_365,ahp,Route)"), scales = 'free') +
  theme_bw()

ggsave("./Figures/GAM_Bivariate_Central_Mean.jpeg", plot = last_plot(), height = 2000, width = 1500, 
       units = "px", dpi = 250)

## 3.1.2.2 Northeast----
#None model works kept it here for reference of what I did.
#Gaussian
system.time(
  m4 <- gam(Relativek ~ #Model using only boruta confirmed selection.
              s(range_365, Season, bs = 'fs') +
              s(min_365, Season, bs = 'fs') +
              s(range_30, Season, bs = 'fs') +
              s(max_30, Season, bs = 'fs') +
              s(range_180, Season, bs = 'fs') +
              ti(range_365, min_365, Season, bs = 'fs') +
              ti(range_365, range_30, Season, bs = 'fs') +
              ti(range_365, max_30, Season, bs = 'fs') + 
              ti(range_365, range_180, Season, bs = 'fs') + 
              s(WY, bs = 're'),
            data = subset(Gators_mean, Groups %in% "Northeast"), method = 'REML',
            nthreads = c(4,1), select = F)
)

m4_sum <- summary(m4) # Deviance explained 72.6% However, none predictors came significant, not useful
k.check(m4) %>% 
  as.data.frame() %>%
  arrange(`p-value`) %>%
  filter(`p-value` < 0.1) %>% 
  rownames_to_column('smooth') %>% 
  dplyr::select(c(1:5)) %>% 
  mutate(edf = round(edf, 1))

system.time(
  m5 <- bam(Relativek ~ #Model including Boruta relevance order (including rejected) plus season as factor smooth
              s(range_365, Season, bs = 'fs') +
              s(min_365, Season, bs = 'fs') +
              s(range_30, Season, bs = 'fs') +
              s(max_30, Season, bs = 'fs') +
              s(range_180, Season, bs = 'fs') +
              s(ahp, Season, bs = 'fs') +
              s(med_90, Season, bs = 'fs') +
              s(min_30, Season, bs = 'fs') +
              s(max_90, Season, bs = 'fs') +
              s(depth, Season, bs = 'fs') +
              s(med_30, Season, bs = 'fs') +
              s(max_365, Season, bs = 'fs') +
              s(min_90, Season, bs = 'fs') +
              s(max_180, Season, bs = 'fs') +
              s(med_365, Season, bs = 'fs') +
              s(range_90, Season, bs = 'fs') +
              s(med_180, Season, bs = 'fs') +
              s(min_180, Season, bs = 'fs') +
              ti(range_365, min_365, Season, bs = 'fs') +
              ti(range_365, range_30, Season, bs = 'fs') +
              ti(range_365, max_30, Season, bs = 'fs') + 
              ti(range_365, range_180, Season, bs = 'fs') + 
              ti(range_365, ahp, Season, bs = 'fs') +
              ti(range_365, med_90, Season, bs = 'fs') +
              ti(range_365, min_30, Season, bs = 'fs') + 
              ti(range_365, max_90, Season, bs = 'fs') + 
              ti(range_365, depth, Season, bs = 'fs') + 
              ti(range_365, med_30, Season, bs = 'fs') + 
              ti(range_365, max_365, Season, bs = 'fs') + 
              ti(range_365, min_90, Season, bs = 'fs') + 
              ti(range_365, max_180, Season, bs = 'fs') + 
              ti(range_365, med_365, Season, bs = 'fs') + 
              ti(range_365, range_90, Season, bs = 'fs') + 
              ti(range_365, med_180, Season, bs = 'fs') + 
              ti(range_365, min_180, Season, bs = 'fs') + 
              s(WY, bs = 're'),
            data = subset(Gators_mean, Groups %in% "Northeast"), method = 'fREML',
            nthreads = c(4,1), select = T, discrete = T)
)

m5_sum <- summary(m5) # Deviance explained 99.8%
k.check(m5) %>% 
  as.data.frame() %>%
  arrange(`p-value`) %>%
  filter(`p-value` < 0.1) %>% 
  rownames_to_column('smooth') %>% 
  dplyr::select(c(1:3)) %>% 
  mutate(edf = round(edf, 1))
appraise(m5, method = 'simulate') #Complete unfit in QQplot.Do not use this model

##Gamma log
system.time(
  m6 <- bam(Relativek ~ #Model including Boruta relevance order (including rejected) plus season as factor smooth
              s(range_365, Season, bs = 'fs') +
              s(min_365, Season, bs = 'fs') +
              s(range_30, Season, bs = 'fs') +
              s(max_30, Season, bs = 'fs') +
              s(range_180, Season, bs = 'fs') +
              s(ahp, Season, bs = 'fs') +
              s(med_90, Season, bs = 'fs') +
              s(min_30, Season, bs = 'fs') +
              s(max_90, Season, bs = 'fs') +
              s(depth, Season, bs = 'fs') +
              s(med_30, Season, bs = 'fs') +
              s(max_365, Season, bs = 'fs') +
              s(min_90, Season, bs = 'fs') +
              s(max_180, Season, bs = 'fs') +
              s(med_365, Season, bs = 'fs') +
              s(range_90, Season, bs = 'fs') +
              s(med_180, Season, bs = 'fs') +
              s(min_180, Season, bs = 'fs') +
              ti(range_365, min_365, Season, bs = 'fs') +
              ti(range_365, range_30, Season, bs = 'fs') +
              ti(range_365, max_30, Season, bs = 'fs') + 
              ti(range_365, range_180, Season, bs = 'fs') + 
              ti(range_365, ahp, Season, bs = 'fs') +
              ti(range_365, med_90, Season, bs = 'fs') +
              ti(range_365, min_30, Season, bs = 'fs') + 
              ti(range_365, max_90, Season, bs = 'fs') + 
              ti(range_365, depth, Season, bs = 'fs') + 
              ti(range_365, med_30, Season, bs = 'fs') + 
              ti(range_365, max_365, Season, bs = 'fs') + 
              ti(range_365, min_90, Season, bs = 'fs') + 
              ti(range_365, max_180, Season, bs = 'fs') + 
              ti(range_365, med_365, Season, bs = 'fs') + 
              ti(range_365, range_90, Season, bs = 'fs') + 
              ti(range_365, med_180, Season, bs = 'fs') + 
              ti(range_365, min_180, Season, bs = 'fs') + 
              s(WY, bs = 're'),
            data = subset(Gators_mean, Groups %in% "Northeast"), method = 'fREML',
            nthreads = c(4,1), family = Gamma(link = log), select = T, discrete = T)
)

m6_sum <- summary(m6) # Warning, do not run, did not work, gives null, leave here to proof it was ran but it 

## 3.1.2.3 South----
system.time(
  m7 <- gam(Relativek ~ #Model using only boruta confirmed selection.
              s(range_365, Season, bs = 'fs') +
              s(WY, bs = 're'),
            data = subset(Gators_mean, Groups %in% "South"), method = 'REML',
            nthreads = c(4,1), select = F)
)

m7_sum <- summary(m7) # Deviance explained 86.2%
appraise(m7, method = 'simulate') # It is not perfect but doesn't look too off

draw(m7, select = "s(range_365,Season)", scale = 'free') & theme_bw()

ggsave("./Figures/GAM_South_Mean.jpeg", plot = last_plot(), height = 1500, width = 2000,
       dpi = 500, units = 'px')

## 3.1.2.4 West----
#Gaussian
#None model works kept it here for reference of what I did.
system.time(
  m8 <- gam(Relativek ~
              s(range_180) +
              s(range_365) +
              s(min_180) +
              s(min_30) +
              ti(range_180, range_365) +
              ti(range_180, min_180) +
              ti(range_180, min_30) + 
              s(Season, bs = 're') +
              s(WY, bs = 're'),
            data = subset(Gators_mean, Groups %in% "West"), method = 'REML',
            nthreads = c(4,1), select = F)
)

m8_sum <- summary(m8) # Deviance explained 63.5% None of predictors significant
appraise(m8, method = 'simulate') # No perfect but works

system.time(
  m9 <- gam(Relativek ~ 
              s(range_180, Season, bs = 'fs') +
              s(range_365, Season, bs = 'fs') +
              s(min_180, Season, bs = 'fs') +
              s(min_30, Season, bs = 'fs') +
              ti(range_180, range_365, Season, bs = 'fs') +
              ti(range_180, min_180, Season, bs = 'fs') +
              ti(range_180, min_30, Season, bs = 'fs') + 
              s(WY, bs = 're'),
            data = subset(Gators_mean, Groups %in% "West"), method = 'REML',
            nthreads = c(4,1), select = F)
)
m9_sum <- summary(m9) #Deviance explained 56.9%
appraise(m9, method = 'simulate')

AIC(m8, m9) # Model 9 de best

draw(m9, select = ("s(range_365,Season)"), scale = 'free') & theme_bw()
ggsave("./Figures/GAM_West_Mean.jpeg", plot = last_plot(), height = 1500, width = 2000,
       dpi = 500, units = 'px')

## 3.2 Individual by group----
## 3.2.1 Variable selection----
set.seed(123)
Boruta_Analysis_Central <- Gators_Covariates %>%
  subset(Groups %in% "Central") %>%
  dplyr::select(22, 19:37) %>%
  drop_na() %>%
  Boruta(Relativek ~ ., data = ., doTrace = 2, maxRuns = 1000)
plot(Boruta_Analysis_Central)
attStats(Boruta_Analysis_Central) %>%
  as.data.frame() %>%
  rownames_to_column('var') %>%
  dplyr::select(var, decision, meanImp) %>%
  arrange(desc(meanImp))

set.seed(123)
Boruta_Analysis_Northeast <- Gators_Covariates %>%
  subset(Groups %in% "Northeast") %>%
  dplyr::select(19:37) %>%
  drop_na() %>%
  Boruta(Relativek ~ ., data = ., doTrace = 2, maxRuns = 1000)
plot(Boruta_Analysis_Northeast)
attStats(Boruta_Analysis_Northeast) %>%
  as.data.frame() %>%
  rownames_to_column('var') %>%
  dplyr::select(var, decision, meanImp) %>%
  arrange(desc(meanImp)) 

set.seed(123)
Boruta_Analysis_South <- Gators_Covariates %>%
  subset(Groups %in% "South") %>%
  dplyr::select(19:37) %>%
  drop_na() %>%
  Boruta(Relativek ~ ., data = ., doTrace = 2, maxRuns = 1000)
plot(Boruta_Analysis_South)
attStats(Boruta_Analysis_South) %>%
  as.data.frame() %>%
  rownames_to_column('var') %>%
  dplyr::select(var, decision, meanImp) %>%
  arrange(desc(meanImp))

set.seed(123)
Boruta_Analysis_West <- Gators_Covariates %>%
  subset(Groups %in% "West") %>%
  dplyr::select(19:37) %>%
  drop_na() %>%
  Boruta(Relativek ~ ., data = ., doTrace = 2, maxRuns = 1000)
plot(Boruta_Analysis_West)
attStats(Boruta_Analysis_West) %>%
  as.data.frame() %>%
  rownames_to_column('var') %>%
  dplyr::select(var, decision, meanImp) %>%
  arrange(desc(meanImp))

## 3.2.2 Models----
## 3.2.2.1 Central----
system.time(
  m10 <- bam(Relativek ~ #Model using only boruta confirmed selection. No rejected variables
               Route +
               s(range_365, by = Route, bs = 'ts') +
               s(range_180, by = Route, bs = 'ts') +
               s(max_365, by = Route, bs = 'ts') +
               s(med_365, by = Route, bs = 'ts') +
               s(min_365, by = Route, bs = 'ts') +
               s(ahp, by = Route, bs = 'ts') +
               s(med_180, by = Route, bs = 'ts') +
               s(range_30, by = Route, bs = 'ts') +
               s(range_90, by = Route, bs = 'ts') + 
               s(max_30, by = Route, bs = 'ts') +
               s(max_90, by = Route, bs = 'ts') +
               s(med_90, by = Route, bs = 'ts') +
               s(depth, by = Route, bs = 'ts') +
               s(med_30, by = Route, bs = 'ts') +
               s(max_180, by = Route, bs = 'ts') +
               s(min_180, by = Route, bs = 'ts') + 
               s(min_30, by = Route, bs = 'ts') +
               s(min_90, by = Route, bs = 'ts') +
               ti(range_365, range_180, by = Route, bs = c('ts', 'ts')) +
               ti(range_365, max_365, by = Route, bs = c('ts', 'ts')) +
               ti(range_365, med_365, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, min_365, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, ahp, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, med_180, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, range_30, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, range_90, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, max_30, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, max_90, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, med_90, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, depth, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, med_30, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, max_180, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, min_180, by = Route, bs = c('ts', 'ts')) +
               ti(range_365, min_30, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, min_90, by = Route, bs = c('ts', 'ts')) + 
               s(WY, bs = 're') +
               s(Sex, bs = 're') +
               s(month, bs = 're') +
               s(SizeClass, bs = 're'), 
             data = subset(Gators_Covariates, Groups %in% "Central"), method = 'fREML',
             nthreads = c(4,1), select = T, discrete = T)
)

m10_sum <- summary(m10) # Deviance explained 20.3%
k.check(m10) %>%  
  as.data.frame() %>%
  arrange(`p-value`) %>%
  filter(`p-value` < 0.1) %>% 
  rownames_to_column('smooth') %>% 
  dplyr::select(c(1:5)) %>% 
  mutate(edf = round(edf, 1))
appraise(m10, method = 'simulate') # It's not perfect but looks good

#Gamma log
system.time(
  m11 <- bam(Relativek ~ #Model using only boruta confirmed selection. No rejected variables
               Route +
               s(range_365, by = Route, bs = 'ts') +
               s(range_180, by = Route, bs = 'ts') +
               s(max_365, by = Route, bs = 'ts') +
               s(med_365, by = Route, bs = 'ts') +
               s(min_365, by = Route, bs = 'ts') +
               s(ahp, by = Route, bs = 'ts') +
               s(med_180, by = Route, bs = 'ts') +
               s(range_30, by = Route, bs = 'ts') +
               s(range_90, by = Route, bs = 'ts') + 
               s(max_30, by = Route, bs = 'ts') +
               s(max_90, by = Route, bs = 'ts') +
               s(med_90, by = Route, bs = 'ts') +
               s(depth, by = Route, bs = 'ts') +
               s(med_30, by = Route, bs = 'ts') +
               s(max_180, by = Route, bs = 'ts') +
               s(min_180, by = Route, bs = 'ts') + 
               s(min_30, by = Route, bs = 'ts') +
               s(min_90, by = Route, bs = 'ts') +
               ti(range_365, range_180, by = Route, bs = c('ts', 'ts')) +
               ti(range_365, max_365, by = Route, bs = c('ts', 'ts')) +
               ti(range_365, med_365, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, min_365, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, ahp, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, med_180, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, range_30, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, range_90, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, max_30, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, max_90, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, med_90, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, depth, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, med_30, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, max_180, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, min_180, by = Route, bs = c('ts', 'ts')) +
               ti(range_365, min_30, by = Route, bs = c('ts', 'ts')) + 
               ti(range_365, min_90, by = Route, bs = c('ts', 'ts')) + 
               s(WY, bs = 're') +
               s(Sex, bs = 're') +
               s(month, bs = 're') +
               s(SizeClass, bs = 're'), 
             data = subset(Gators_Covariates, Groups %in% "Central"), method = 'fREML',
             nthreads = c(4,1), family = Gamma(link = log), select = T, discrete = T)
)

m11_sum <- summary(m11) # Deviance explained 20.2%

k.check(m11) %>%  
  as.data.frame() %>%
  arrange(`p-value`) %>%
  filter(`p-value` < 0.1) %>% 
  rownames_to_column('smooth') %>% 
  dplyr::select(c(1:5)) %>% 
  mutate(edf = round(edf, 1))
appraise(m11, method = 'simulate') # It's not perfect but looks good

#Test models
AIC(m10, m11) # m11 has the lowest AIC. Model selected

#univariate
draw(m11, select = c("s(max_365):RouteWCA3A-TW_Marsh",
                     "s(med_365):RouteWCA3A-N41_Marsh",
                     "s(med_90):RouteWCA3A-N41_Marsh",
                     "s(med_365):RouteWCA3B_Marsh",
                     "s(min_180):RouteWCA3B_Marsh",
                     "s(min_365):RouteENP-FC",
                     "s(med_180):RouteENP-FC",
                     "s(max_180):RouteENP-FC",
                     "s(range_180):RouteWCA3A-HD_Marsh",
                     "s(min_30):RouteWCA3A-HD_Marsh",
                     "s(min_30):RouteWCA2A_Marsh",
                     "s(min_90):RouteWCA2A_Marsh"), scales = 'free') +
  plot_layout(ncol = 3, nrow = 4, guides = "collect") & theme_bw()

ggsave("./Figures/GAM_Central_Ind.jpeg", plot = last_plot(), height = 2500, 
       width = 2000, units = "px", dpi = 250)

#bivariate
draw(m11, select = c("ti(range_365,max_365):RouteWCA3A-N41_Marsh",
                     "ti(range_365,med_365):RouteWCA3A-N41_Marsh",
                     "ti(range_365,med_365):RouteWCA3B_Marsh",
                     "ti(range_365,med_180):RouteWCA3A-N41_Marsh",
                     "ti(range_365,med_180):RouteWCA3A-HD_Marsh",
                     "ti(range_365,range_30):RouteWCA2A_Marsh",
                     "ti(range_365,range_90):RouteWCA3A-HD_Marsh",
                     "ti(range_365,range_90):RouteWCA3B_Marsh",
                     "ti(range_365,med_90):RouteWCA3B_Marsh",
                     "ti(range_365,max_30):RouteWCA3B_Marsh",
                     "ti(range_365,max_180):RouteWCA3A-TW_Marsh",
                     "ti(range_365,max_180):RouteWCA3B_Marsh", 
                     "ti(range_365,min_180):RouteENP-FC",
                     "ti(range_365,min_180):RouteWCA3A-TW_Marsh",
                     "ti(range_365,min_30):RouteENP-FC"), scales = 'free') +
  theme_bw()

#filtered
draw(m11, select = c("ti(range_365,med_180):RouteWCA3A-N41_Marsh",
                     "ti(range_365,max_30):RouteWCA3B_Marsh",
                     "ti(range_365,med_90):RouteWCA3B_Marsh",
                     "ti(range_365,min_180):RouteWCA3A-TW_Marsh"), scales = 'free') +
  theme_bw()


ggsave("./Figures/GAM_Central_Ind_bivariate.jpeg", plot = last_plot(), height = 1500, 
       width = 2000, units = "px", dpi = 250)

## 3.3.2.2 Northeast----
system.time(
  m12 <-bam(Relativek ~ 
              Season +
              s(range_365, by = Season, bs = 'ts') +
              s(min_365, by = Season, bs = 'ts') +
              s(range_180, by = Season, bs = 'ts') + 
              s(max_180, by = Season, bs = 'ts') +
              s(range_90, by = Season, bs = 'ts') +
              s(max_90, by = Season, bs = 'ts') +
              s(med_180, by = Season, bs = 'ts') +
              s(med_90, by = Season, bs = 'ts') +
              s(range_30, by = Season, bs = 'ts') +
              s(min_180, by = Season, bs = 'ts') +
              s(max_30, by = Season, bs = 'ts') +
              s(ahp, by = Season, bs = 'ts') +
              s(depth, by = Season, bs = 'ts') +
              s(min_30, by = Season, bs = 'ts') +
              s(med_30, by = Season, bs = 'ts') +
              s(max_365, by = Season, bs = 'ts') +
              s(min_90, by = Season, bs = 'ts') +
              s(med_365, by = Season, bs = 'ts') +
              ti(range_365, min_365, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, range_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, range_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, range_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, ahp, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, depth, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_365, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_365, by = Season, bs = c('ts', 'ts')) +
              s(WY, bs = 're') +
              s(Sex, bs = 're') +
              s(month, bs = 're') +
              s(SizeClass, bs = 're'),
            data = subset(Gators_Covariates, Groups %in% "Northeast"), method = 'fREML', 
            nthreads = c(4, 1), select = T, discrete = T)
)
m12_sum <- summary(m12) #Check deviance explained
appraise(m12, method = 'simulate') # Ensure distribution model fit data adequately  

system.time(
  m13 <-bam(Relativek ~ 
              Season +
              s(range_365, by = Season, bs = 'ts') +
              s(min_365, by = Season, bs = 'ts') +
              s(range_180, by = Season, bs = 'ts') + 
              s(max_180, by = Season, bs = 'ts') +
              s(range_90, by = Season, bs = 'ts') +
              s(max_90, by = Season, bs = 'ts') +
              s(med_180, by = Season, bs = 'ts') +
              s(med_90, by = Season, bs = 'ts') +
              s(range_30, by = Season, bs = 'ts') +
              s(min_180, by = Season, bs = 'ts') +
              s(max_30, by = Season, bs = 'ts') +
              s(ahp, by = Season, bs = 'ts') +
              s(depth, by = Season, bs = 'ts') +
              s(min_30, by = Season, bs = 'ts') +
              s(med_30, by = Season, bs = 'ts') +
              s(max_365, by = Season, bs = 'ts') +
              s(min_90, by = Season, bs = 'ts') +
              s(med_365, by = Season, bs = 'ts') +
              ti(range_365, min_365, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, range_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, range_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, range_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, ahp, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, depth, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_365, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_365, by = Season, bs = c('ts', 'ts')) +
              s(WY, bs = 're') +
              s(Sex, bs = 're') +
              s(month, bs = 're') +
              s(SizeClass, bs = 're'),
            data = subset(Gators_Covariates, Groups %in% "Northeast"), method = 'fREML', 
            nthreads = c(4, 1), family = Gamma(link = log), select = T, discrete = T)
)

m13_sum <- summary(m13)
appraise(m13, method = 'simulate')

AIC(m12, m13) # m13 best model. model chosen.

draw(m13, select = c("s(min_365):SeasonFall",
                     "s(max_365):SeasonSpring")) + 
  theme_bw() 

ggsave("./Figures/GAM_Northeast_Ind_1.jpeg", plot = last_plot(), height = 1500, 
       width = 2500, units = "px", dpi = 450)

draw(m13, nrow = 1, select = c("ti(range_365,range_180):SeasonFall",
                               "ti(range_365,max_180):SeasonFall",
                               "ti(range_365,min_30):SeasonSpring")) + 
  theme_bw() 

ggsave("./Figures/GAM_Northeast_Ind_2.jpeg", plot = last_plot(), height = 1500, 
       width = 4000, units = "px", dpi = 450)

## 3.3.2.3 South----
system.time(
  m14 <-bam(Relativek ~ 
              Season +
              s(range_365, by = Season, bs = 'ts') +
              s(max_365, by = Season, bs = 'ts') +
              s(min_30, by = Season, bs = 'ts') + 
              s(max_180, by = Season, bs = 'ts') +
              s(max_90, by = Season, bs = 'ts') +
              s(med_30, by = Season, bs = 'ts') +
              s(min_365, by = Season, bs = 'ts') +
              s(ahp, by = Season, bs = 'ts') +
              s(med_90, by = Season, bs = 'ts') +
              s(min_90, by = Season, bs = 'ts') +
              s(max_30, by = Season, bs = 'ts') +
              s(depth, by = Season, bs = 'ts') +
              s(range_180, by = Season, bs = 'ts') +
              s(min_180, by = Season, bs = 'ts') +
              s(med_180, by = Season, bs = 'ts') +
              s(range_30, by = Season, bs = 'ts') +
              s(range_90, by = Season, bs = 'ts') +
              s(med_365, by = Season, bs = 'ts') +
              ti(range_365, max_365, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_365, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, ahp, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, depth, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, range_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, range_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, range_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_365, by = Season, bs = c('ts', 'ts')) +
              s(WY, bs = 're') +
              s(Sex, bs = 're') +
              s(month, bs = 're') +
              s(SizeClass, bs = 're'),
            data = subset(Gators_Covariates, Groups %in% "South"), method = 'fREML', 
            nthreads = c(4, 1), select = T, discrete = T)
)
m14_sum <- summary(m14) #Check deviance explained (34.4%)
appraise(m14, method = 'simulate') # Good  

system.time(
  m15 <-bam(Relativek ~ 
              Season +
              s(range_365, by = Season, bs = 'ts') +
              s(max_365, by = Season, bs = 'ts') +
              s(min_30, by = Season, bs = 'ts') + 
              s(max_180, by = Season, bs = 'ts') +
              s(max_90, by = Season, bs = 'ts') +
              s(med_30, by = Season, bs = 'ts') +
              s(min_365, by = Season, bs = 'ts') +
              s(ahp, by = Season, bs = 'ts') +
              s(med_90, by = Season, bs = 'ts') +
              s(min_90, by = Season, bs = 'ts') +
              s(max_30, by = Season, bs = 'ts') +
              s(depth, by = Season, bs = 'ts') +
              s(range_180, by = Season, bs = 'ts') +
              s(min_180, by = Season, bs = 'ts') +
              s(med_180, by = Season, bs = 'ts') +
              s(range_30, by = Season, bs = 'ts') +
              s(range_90, by = Season, bs = 'ts') +
              s(med_365, by = Season, bs = 'ts') +
              ti(range_365, max_365, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_365, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, ahp, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, max_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, depth, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, range_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, min_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_180, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, range_30, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, range_90, by = Season, bs = c('ts', 'ts')) +
              ti(range_365, med_365, by = Season, bs = c('ts', 'ts')) +
              s(WY, bs = 're') +
              s(Sex, bs = 're') +
              s(month, bs = 're') +
              s(SizeClass, bs = 're'),
            data = subset(Gators_Covariates, Groups %in% "South"), method = 'fREML', 
            nthreads = c(4, 1), family = Gamma(link = log), select = T, discrete = T)
)
m15_sum <- summary(m15) #Check deviance explained (31.7%)
appraise(m15, method = 'simulate') # Ensure distribution model fit data adequately

AIC(m14, m15) # m14  best model. model chosen.

draw(m14, select = c("s(range_30):SeasonFall")) + 
  theme_bw() 

ggsave("./Figures/GAM_South_Ind_1.jpeg", plot = last_plot(), height = 1500, 
       width = 1250, units = "px", dpi = 450)

draw(m14, nrow = 1, select = c("ti(range_365,max_180):SeasonSpring",
                               "ti(range_365,med_90):SeasonSpring",
                               "ti(range_365,depth):SeasonFall",
                               "ti(range_365,range_30):SeasonSpring")) + 
  theme_bw() 

ggsave("./Figures/GAM_South_Ind_2.jpeg", plot = last_plot(), height = 1500, 
       width = 4000, units = "px", dpi = 450)

## 3.3.2.4 West----
system.time(
  m16 <-bam(Relativek ~ 
              Season +
              s(med_180, by = Season, bs = 'ts') +
              s(min_180, by = Season, bs = 'ts') +
              s(range_180, by = Season, bs = 'ts') + 
              s(med_365, by = Season, bs = 'ts') +
              s(range_365, by = Season, bs = 'ts') +
              s(ahp, by = Season, bs = 'ts') +
              s(min_365, by = Season, bs = 'ts') +
              s(min_90, by = Season, bs = 'ts') +
              s(med_90, by = Season, bs = 'ts') +
              s(max_180, by = Season, bs = 'ts') +
              s(max_365, by = Season, bs = 'ts') +
              s(max_30, by = Season, bs = 'ts') +
              s(min_30, by = Season, bs = 'ts') +
              s(max_90, by = Season, bs = 'ts') +
              s(depth, by = Season, bs = 'ts') +
              s(med_30, by = Season, bs = 'ts') +
              ti(med_180, min_180, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, range_180, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, med_365, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, range_365, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, ahp, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, min_365, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, min_90, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, med_90, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, max_180, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, max_365, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, max_30, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, min_30, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, max_90, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, depth, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, med_30, by = Season, bs = c('ts', 'ts')) +
              s(WY, bs = 're') +
              s(Sex, bs = 're') +
              s(month, bs = 're') +
              s(SizeClass, bs = 're'),
            data = subset(Gators_Covariates, Groups %in% "West"), method = 'fREML', 
            nthreads = c(4, 1), select = T, discrete = T)
)
m16_sum <- summary(m16) #Check deviance explained (25.7%)
appraise(m16, method = 'simulate') # A bit of in the upper extreme  

system.time(
  m17 <-bam(Relativek ~ 
              Season +
              s(med_180, by = Season, bs = 'ts') +
              s(min_180, by = Season, bs = 'ts') +
              s(range_180, by = Season, bs = 'ts') + 
              s(med_365, by = Season, bs = 'ts') +
              s(range_365, by = Season, bs = 'ts') +
              s(ahp, by = Season, bs = 'ts') +
              s(min_365, by = Season, bs = 'ts') +
              s(min_90, by = Season, bs = 'ts') +
              s(med_90, by = Season, bs = 'ts') +
              s(max_180, by = Season, bs = 'ts') +
              s(max_365, by = Season, bs = 'ts') +
              s(max_30, by = Season, bs = 'ts') +
              s(min_30, by = Season, bs = 'ts') +
              s(max_90, by = Season, bs = 'ts') +
              s(depth, by = Season, bs = 'ts') +
              s(med_30, by = Season, bs = 'ts') +
              ti(med_180, min_180, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, range_180, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, med_365, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, range_365, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, ahp, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, min_365, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, min_90, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, med_90, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, max_180, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, max_365, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, max_30, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, min_30, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, max_90, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, depth, by = Season, bs = c('ts', 'ts')) +
              ti(med_180, med_30, by = Season, bs = c('ts', 'ts')) +
              s(WY, bs = 're') +
              s(Sex, bs = 're') +
              s(month, bs = 're') +
              s(SizeClass, bs = 're'),
            data = subset(Gators_Covariates, Groups %in% "West"), method = 'fREML', 
            nthreads = c(4, 1), family = Gamma(link = log), select = T, discrete = T)
)
m17_sum <- summary(m17) #Check deviance explained (25.3%)
appraise(m17, method = 'simulate') # Good  

AIC(m16, m17) # m17  best model. model chosen.

draw(m17, nrow = 1, select = c("s(min_180):SeasonSpring",
                               "s(min_30):SeasonSpring",
                               "s(med_30):SeasonSpring")) + 
  theme_bw() 

ggsave("./Figures/GAM_West_Ind_1.jpeg", plot = last_plot(), height = 1500, 
       width = 3750, units = "px", dpi = 450)

draw(m17, select = c("ti(med_180,min_30):SeasonFall",
                     "ti(med_180,max_90):SeasonFall")) + 
  theme_bw() 

ggsave("./Figures/GAM_West_Ind_2.jpeg", plot = last_plot(), height = 1500, 
       width = 3750, units = "px", dpi = 450)

