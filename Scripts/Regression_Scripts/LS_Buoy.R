library(tidyverse)
library(lubridate)
library(corrplot)
library(SARP.compo)

#import buoy data
buoy <- read.csv("Data/Buoy/daily_buoy_tempcorr_clean.csv") %>%  
  mutate(sampledate = ymd(sampledate)) %>% 
  select(-X)

# read LS data

ls8 <- read.csv("Data/Landsat/LS8_Buoy100m.csv") %>% 
  filter(pixelCount == 32) %>% 
  mutate(DATE_ACQUIRED = ymd(DATE_ACQUIRED)) %>% 
  dplyr::rename('Aerosol' = SR_B1,
         'Blue' = SR_B2,
         'Green' = SR_B3,
         'Red' = SR_B4,
         'NIR' = SR_B5,
         'SWIR1' = SR_B6,
         'SWIR2' = SR_B7,
         'Thermal' = ST_B10) %>% 
  filter(Aerosol>0,Blue>0,Green>0,Red>0,NIR>0,SWIR1>0,SWIR2>0) 

band_ts <- ls8 %>% 
  pivot_longer(cols=2:9,names_to='band',values_to='value')

ggplot(band_ts,aes(DATE_ACQUIRED,value))+
  geom_point()+
  geom_line()+
  facet_wrap(~band,scales='free')

#merge datasets
combined <- merge(x = buoy, y = ls8,
                  by.x = 'sampledate', by.y = 'DATE_ACQUIRED',
                  all.x = F, all.y = T)

#seasonal distribution of available data
combined %>%
  filter(!is.na(avg_chlor_rfu_corr)) %>%
  ggplot()+
  geom_histogram(aes(x=month(sampledate)))

#create all unique band ratios
band_combos <- combined %>% 
  select(sampledate,8:15)

names <- c('Aerosol','Blue','Green','Red','NIR','SWIR1','SWIR2','Thermal')

all_ratios <- calc.rapports(band_combos, noms=names, log = FALSE, isoler = FALSE )

ar_merge <- merge(x = buoy, y = all_ratios,
                  by.x = 'sampledate', by.y = 'sampledate',
                  all.x = F, all.y = T)

#filter out days with ice cover
ice <- read.csv('Data/IceCover/cleaned_ice_cover.csv') %>% 
  select(-X) %>% 
  mutate(Date = ymd(Date))

ice_filt <- unique(merge(x = ar_merge, y = ice,
                         by.x = 'sampledate', by.y = 'Date',
                         all.x = T, all.y = F)) %>% 
  filter(Ice == 'No')
#removes 14 rows

allratioreg <- ice_filt %>% 
  pivot_longer(cols = 16:43, names_to = 'ratio', values_to = 'value')

corrplot(cor(ar_merge[,-1],use='pairwise'),type='lower')

cors <- cor(ar_merge[,-1],use='pairwise')

#Chl--------------------------------
#make corrplot for chl-a
chlcors = cor(ar_merge[,c(2,18,20:55)], use = 'pairwise')
corrplot(chlcors, type = 'lower', tl.cex = 0.7)

as_tibble(chlcors, rownames = NA) %>% 
  rownames_to_column() %>% 
  select(1:3) %>% 
  mutate(abscor = abs(avg_chlor_rfu),
         abs_scale_cor = abs(scaled_chlor_rfu)) %>% 
  arrange(-abscor) %>% 
  head(10)
as_tibble(chlcors, rownames = NA) %>% 
  rownames_to_column() %>% 
  select(1:3) %>% 
  mutate(abscor = abs(avg_chlor_rfu),
         abs_scale_cor = abs(scaled_chlor_rfu)) %>% 
  arrange(-abs_scale_cor) %>% 
  head(10)

#unscaled
ggplot(allratioreg%>% filter(!is.na(avg_chlor_rfu)),aes(avg_chlor_rfu,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#scaled
ggplot(allratioreg%>% filter(!is.na(scaled_chlor_rfu)),aes(scaled_chlor_rfu,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#run linear regressions
chl_reg <- ice_filt %>% 
  select(avg_chlor_rfu,20:55) %>% 
  filter(!is.na(avg_chlor_rfu))

bands <- colnames(chl_reg)
chl_lm <- data.frame(band = rep(NA,37),
                     r2 = rep(NA,37),
                     p = rep(NA,37),
                     slope = rep(NA,37),
                     int = rep(NA,37))

for(k in 2:length(bands)){
  df = chl_reg %>% 
    select(avg_chlor_rfu,bands[k])
  
  chl_lm$band[k-1] = bands[k]
  chl_lm$r2[k-1] = summary(lm(avg_chlor_rfu~.,data=df))$adj.r.squared
  chl_lm$p[k-1] = summary(lm(avg_chlor_rfu~.,data=df))$coefficients[2,4]
  chl_lm$slope[k-1] = summary(lm(avg_chlor_rfu~.,data=df))$coefficients[2]
  chl_lm$int[k-1] = summary(lm(avg_chlor_rfu~.,data=df))$coefficients[1]
}

sig_chl_lm <- chl_lm %>% 
  filter(p <= 0.05)

#chl does not seem to have a meaningful relationship with any band ratios

#Phyco------------------------------

#unscaled
ggplot(allratioreg%>% filter(!is.na(avg_phyco_rfu)),aes(avg_phyco_rfu,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#scaled
ggplot(allratioreg%>% filter(!is.na(scaled_phyco_rfu)),aes(scaled_phyco_rfu,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#run linear regressions
phyco_reg <- ice_filt %>% 
  select(avg_phyco_rfu,20:55) %>% 
  filter(!is.na(avg_phyco_rfu))

phyco_lm <- data.frame(band = rep(NA,37),
                       r2 = rep(NA,37),
                       p = rep(NA,37),
                       slope = rep(NA,37),
                       int = rep(NA,37))

for(k in 2:length(bands)){
  df = phyco_reg %>% 
    select(avg_phyco_rfu,bands[k])
  
  phyco_lm$band[k-1] = bands[k]
  phyco_lm$r2[k-1] = summary(lm(avg_phyco_rfu~.,data=df))$adj.r.squared
  phyco_lm$p[k-1] = summary(lm(avg_phyco_rfu~.,data=df))$coefficients[2,4]
  phyco_lm$slope[k-1] = summary(lm(avg_phyco_rfu~.,data=df))$coefficients[2]
  phyco_lm$int[k-1] = summary(lm(avg_phyco_rfu~.,data=df))$coefficients[1]
}

sig_phyco_lm <- phyco_lm %>% 
  filter(p <= 0.05)

#not very good

#fdom------------------------------

ggplot(allratioreg%>% filter(!is.na(avg_fdom)),aes(avg_fdom,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#run linear regressions
fdom_reg <- ice_filt %>% 
  select(avg_fdom,20:55) %>% 
  filter(!is.na(avg_fdom))

fdom_lm <- data.frame(band = rep(NA,37),
                      r2 = rep(NA,37),
                      p = rep(NA,37),
                      slope = rep(NA,37),
                      int = rep(NA,37))

for(k in 2:length(bands)){
  df = fdom_reg %>% 
    select(avg_fdom,bands[k])
  
  fdom_lm$band[k-1] = bands[k]
  fdom_lm$r2[k-1] = summary(lm(avg_fdom~.,data=df))$adj.r.squared
  fdom_lm$p[k-1] = summary(lm(avg_fdom~.,data=df))$coefficients[2,4]
  fdom_lm$slope[k-1] = summary(lm(avg_fdom~.,data=df))$coefficients[2]
  fdom_lm$int[k-1] = summary(lm(avg_fdom~.,data=df))$coefficients[1]
}

sig_fdom_lm <- fdom_lm %>% 
  filter(p <= 0.05)

#also not good

#turbidity------------------------------
allratioreg%>% filter(!is.na(avg_turbidity)) %>% 
  filter(avg_turbidity < 7) %>% #remove 1 outlier?
ggplot(aes(avg_turbidity,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

ggplot(allratioreg%>% filter(!is.na(avg_turbidity)),aes(avg_turbidity,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#run linear regressions
turb_reg <- ice_filt %>% 
  select(avg_turbidity,20:55) %>% 
  filter(!is.na(avg_turbidity), avg_turbidity < 7)

turb_lm <- data.frame(band = rep(NA,37),
                      r2 = rep(NA,37),
                      p = rep(NA,37),
                      slope = rep(NA,37),
                      int = rep(NA,37))

for(k in 2:length(bands)){
  df = turb_reg %>% 
    select(avg_turbidity,bands[k])
  
  turb_lm$band[k-1] = bands[k]
  turb_lm$r2[k-1] = summary(lm(avg_turbidity~.,data=df))$adj.r.squared
  turb_lm$p[k-1] = summary(lm(avg_turbidity~.,data=df))$coefficients[2,4]
  turb_lm$slope[k-1] = summary(lm(avg_turbidity~.,data=df))$coefficients[2]
  turb_lm$int[k-1] = summary(lm(avg_turbidity~.,data=df))$coefficients[1]
}

sig_turb_lm <- turb_lm %>% 
  filter(p <= 0.05)

#lots of good regressions





