library(tidyverse)

#designate seasons-------------
empty <- data.frame(date = seq(as.Date("01-01-2000", format = "%m-%d-%Y"), as.Date("12-31-2022", format = "%m-%d-%Y"), by = "day")) %>%
  mutate(doy = yday(date),
         year = year(date))
seasons <- read.csv('MetabolismModel/seasons_from_Robins_paper_long_withplaceholders.csv') %>%
  select(-placeholder) %>%
  mutate(date = as.Date(date,format = "%m/%d/%y"))
DATA <- left_join(empty, seasons) %>%
  fill(season, .direction="down")

DATA$season <- factor(DATA$season, levels = c("Ice-on", "Spring", "Clearwater", "Early Summer", "Late Summer", "Fall"))

# new <- merge(x=pca_new,y=DATA,
#              by.x='date',by.y='date') %>% 
#   mutate(img = case_when(season == 'Spring' ~ 1,
#                          season == 'Clearwater' ~ 1,
#                          season == 'Early Summer' ~ 1,
#                          season == 'Late Summer' & date == '2021-09-05' ~ 1,
#                          season == 'Late Summer' & date == '2022-07-27' ~ 2,
#                          season == 'Fall' & date == '2021-11-02' ~ 1,
#                          season == 'Fall' & date == '2021-12-13' ~ 2))
# 
# scale_color_manual(values=c("#73456D",  "#B2CCF1","#EE914A", "#9AD67A","#138E90"),name='Season') 
  


#buoy and metabolism regression (NPP vs phyco)-----------------------------------

buoy <- read.csv("Data/Buoy/daily_buoy_tempcorr_clean.csv") %>%  
  mutate(sampledate = ymd(sampledate)) %>% 
  select(-X)

metab <- read.csv('Data/Metabolism_Model_Output/SimResultsMatrix_MetabData_run05oct23.csv') %>% 
  mutate(SimDate = ymd(SimDate)) 


buoy_metab <- merge(x=buoy,y=metab,by.x='sampledate',by.y='SimDate')

buoy_metab_plot <- merge(x=buoy_metab,y=DATA,
                         by.x='sampledate',by.y='date') %>% 
  filter(season != 'Ice-on') %>% 
  filter(EpiNPP_mgC_L < 0.4)

p1 <- ggplot(buoy_metab_plot,aes(avg_phyco_rfu_corr,EpiNPP_mgC_L,color=season))+
  geom_point()+
  scale_color_manual(values=c("#73456D",  "#B2CCF1","#EE914A", "#9AD67A","#138E90"),name='Season') +
  scale_x_log10()+
  theme_bw()+
  geom_smooth(color='grey',method='lm')+
  labs(y='NPP (mgC/L/d)',
       x='Phycocyanin (RFU)')+
  theme(legend.position='none')
p1

#landsat and metabolism regression-------------------------------------

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
  filter(Aerosol>0,Blue>0,Green>0,Red>0,NIR>0,SWIR1>0,SWIR2>0) %>% 
  select(DATE_ACQUIRED,Thermal)

ls_metab <- merge(x=metab,y=ls8,
                  by.x='SimDate',by.y='DATE_ACQUIRED')

ls_metab_plot <- merge(x=ls_metab,y=DATA,
                       by.x='SimDate',by.y='date') %>% 
  select(Thermal,EpiDO_mgO2_L,EpiDO_sat,EpiNPP_mgC_L,EpiR_mgC_L,season) %>% 
  filter(season != 'Ice-on') %>% 
  pivot_longer(2:5,names_to='var',values_to='val') %>% 
  mutate(var = case_when(var == 'EpiDO_mgO2_L' ~ 'Dissolved Oxygen (mg/L)',
                         var == 'EpiDO_sat' ~ 'Dissolved Oxygen (%)',
                         var == 'EpiNPP_mgC_L' ~ 'NPP (mgC/L/d)',
                         var == 'EpiR_mgC_L' ~ 'Respiration (mgC/L/d)'))
  
p2 <- ggplot(ls_metab_plot,aes(val,Thermal,color=season))+
  geom_point()+
  scale_color_manual(values=c("#73456D",  "#B2CCF1","#EE914A", "#9AD67A","#138E90"),name='Season') +
  facet_wrap(~var,scales='free')+
  geom_smooth(color='grey',method='lm')+
  theme_bw()+
  labs(y='Landsat 8 Surface Temperature (K)',
       x='Value')+
  theme(legend.position = 'none')
p2

#sentinel and buoy regression------------------------  

s2 <- read.csv("Data/Sentinel2/S2_Buoy100mMean.csv") %>% 
  mutate(Date = substr(DATATAKE_IDENTIFIER, 6, 13),
         Year = substr(Date,0,4),
         Month = substr(Date,5,6),
         Day = substr(Date,7,8),
         Date = ymd(paste(Year,Month,Day,sep='-'))) %>% 
  select(-DATATAKE_IDENTIFIER,-Year,-Month,-Day) %>% 
  relocate(Date) %>% 
  rename('Aerosol' = B1,
         'Blue' = B2,
         'Green' = B3,
         'Red' = B4,
         'RedEdge1' = B5,
         'RedEdge2' = B6,
         'RedEdge3' = B7,
         'NIR' = B8,
         'RedEdge4' = B8A,
         'SWIR1' = B11,
         'SWIR2' = B12) %>% 
  filter(Aerosol>0,Blue>0,Green>0,Red>0,NIR>0,SWIR1>0,SWIR2>0) %>% 
  filter(pixelCount >300) %>% 
  filter(MGRS_TILE == '15TYH') %>% 
  mutate(ratio = Red/RedEdge1) %>% 
  select(Date,ratio)


s2_buoy <- merge(x=s2,y=buoy,by.x='Date',by.y='sampledate') %>% 
  filter(!is.na(avg_phyco_rfu_corr))


s2_buoy_plot <- merge(x=s2_buoy,y=DATA,by.x='Date',by.y='date') %>% 
  filter(season != 'Ice-on')

p3 <- ggplot(s2_buoy_plot,aes(avg_phyco_rfu_corr,ratio,color=season))+
  geom_point()+
  geom_smooth(color='grey',method='lm')+
  scale_color_manual(values=c("#73456D",  "#B2CCF1","#EE914A", "#9AD67A","#138E90"),name='Season') +
  labs(y='Sentinel 2 Red/Red Edge 1',
       x='Phycocyanin (RFU)')+
  theme_bw()+
  theme(legend.position='none')
p3

#arrange plots----------------------
library(ggpubr)

arranged <- ggarrange(p1,p2,p3,ncol=1,nrow=3,labels=c('a','b','c'),common.legend=TRUE,legend='right')
arranged

ggsave(filename='Figures/fig3.jpg',
       plot=arranged,
       width=5.5,height=10,
       units='in',dpi=320)

