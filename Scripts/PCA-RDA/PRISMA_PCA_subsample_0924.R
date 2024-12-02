# PRISMA PCA comparisons -----------------------------------------------------

# Load libraries ----------------------------------------------
library(sp)
#library(raster)
library(tidyverse)
library(plyr)
library(corrplot)
library(caret)
library(rpart)
library(rattle)
library(ade4)
library(plyr)

prisma <- read.csv("Data/PRISMA/L2W_vnorm/vnorm_prisma_lake_500clip_noBin.csv")
prisma_wide <- pivot_wider(prisma, 
                           id_cols = c(index,date,x,y), 
                           names_from = wavelength,
                           values_from = norm)
# wav <- as.numeric(colnames(prisma_wide[,-c(1:4)]))
# matplot(x = wav,
#         y = t(prisma_wide[,-c(1:4)]),
#         type = "l",
#         xlab = "Wavelength (nm)",
#         ylab = "Reflectance")

library(prospectr)
# https://rdrr.io/cran/prospectr/man/savitzkyGolay.html
# Savitzy-Golay Smoothing Filter with SMASH parameters for DESIS:
# m	 = an integer indicating the differentiation order       = 2
# p	 = an integer indicating the polynomial order           = 3
# w	 = an integer indicating the window size (must be odd)  = 7
index <- prisma_wide[,1:4]
prisma_smooth <- savitzkyGolay(X = prisma_wide[,-c(1:4)], m=2 , p=3 , w=7 )
# wav <- as.numeric(colnames(prisma_smooth))
# matplot(x = wav,
#         y = t(prisma_smooth),
#         type = "l",
#         xlab = "Wavelength (nm)",
#         ylab = "Reflectance")

# Load list of scenes

flist = dir("Data/PRISMA/L2W_vnorm/lake")
#Initialize storage for variances and components
vars = tibble()
lds = tibble()

# complying done --------------------------------------
## ONLY NEED TO RUN ONCE AND THEN CAN LOAD COMPILED CSV IN LINE BELOW (COMMENT THE REST OUT)
# compile PRISMA scenes to one common file
# PRISMA lake clipped data
# 
# # #Load first file and add to PRISMA_lake that will eventually contain all PRISMA lake data
prisma <- read.csv(file.path("Data/PRISMA/L2W_vnorm/lake",flist[1]))
 
date = flist[1]
PRISMA_lake <- cbind(date,prisma)

 for(i in 2:length(flist)){
   prisma <- read.csv(file.path("Data/PRISMA/L2W_vnorm/lake",flist[i]))
   date = flist[i]
  temp <- cbind(date,prisma)
  PRISMA_lake <- plyr::rbind.fill(PRISMA_lake,temp)
 }
#
 PRISMA_lake <- PRISMA_lake %>% mutate(date = substr(date,8,17))
#
# # #export summary file
# write.csv(PRISMA_lake,'Data/PRISMA/Vnorm/PRISMA_lake_compiled.csv')

# Load compiled csv ---------------------------------------------------
# MAKE SURE YOU PULL PRISMA_lake_compiled off Google Drive: 
# https://drive.google.com/drive/u/1/folders/1z6HNn4lAW20FnxJWGRZyuyDgpvSD1sYu
# And save in Data/PRISMA
PRISMA_lake <- read.csv('Data/PRISMA/L2W_Vnorm/vnorm_prisma_lake_500clip_noBin.csv')


##### Subsample
set.seed(1234) # to make runs consistent
# subsample 10% PRISMA data because 
# Critical Scale of Variability ~100m
# PRISMA pixels are 30m, so 3x3 grid of PRISMA pixels is ~100x100m 
# Which means we subsample 1/9 = 11.1% so 10% is fine
prsma_sub <- prisma_wide %>%
  rownames_to_column("X") %>% 
  # dplyr::select(-"765") %>% #drop 765 as over half are missing
  group_by(date) %>% 
  dplyr::select(index, c("505":"796")) %>%
  drop_na() %>%
slice_sample(prop=0.10)
# run pca
prsm_s2 <- prsma_sub %>% ungroup(date) %>% select(-date, -index)

pca_prsm <- princomp(prsm_s2,cor = TRUE)
  
summary(pca_prsm)
biplot(pca_prsm)

screeplot(pca_prsm)#shows variance represented by each component

#save PCA output
pca_prsm_scores = pca_prsm$scores
rownames(pca_prsm_scores) = prsma_sub$index
pca_prsm_ldgs = pca_prsm$loadings
saveRDS(pca_prsm_ldgs, "Outputs/PCA_outputs/PRISMA/pca_prsm_ldgs.rds")
#saveRDS(pca_prsm_scores, "Outputs/PCA_outputs/PRISMA/pca_prsm_scores.rds")
#write_csv(PRISMA_lake[,c(2,4,3,1,5:8)],"Outputs/PCA_outputs/PRISMA/prisma_geo.csv" )

write_csv(rownames_to_column(PRISMA_lake, "X")[,c(1,2,158,159,162)],"Outputs/PCA_outputs/PRISMA/prisma_geo.csv")

#improve biplot
library(ggfortify)
peu <- autoplot(pca_prsm,loadings=TRUE,loadings.label=T,
                loadings.label.repel=T, loadings.colour="black",loadings.label.colour="black")+
  theme_bw()
peu
# or to include date
#install.packages('ggord')
options(repos = c(
  fawda123 = 'https://fawda123.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

library(factoextra)
library(ggord)
# Perform PCA
res.pca <- prcomp(your_data)
# Extract loadings
loadings <- res.pca$rotation

loadings <- data.frame(pca_prsm_ldgs) %>% 
  select(1:2)
# Calculate absolute values of loadings for each PC
abs_loadings <- abs(pca_prsm_ldgs)
# Set a threshold (e.g., top 20% most influential variables)
threshold <- quantile(abs_loadings, 0.8)
# Create the biplot with filtering
ggord(pca_prsm, scale = ifelse(abs_loadings > threshold, 1, 0))

# Install ggord
install.packages('ggord')
require(ggord)

psm_gg <- ggord(pca_prsm, prsma_sub$date,
              cols=c("firebrick4","orange","yellowgreen","cornflowerblue",
                     "purple3","pink3","deeppink2"),
              vec_ext=4.5,vectyp="solid",xlims=c(-4,3),ylims=c(-3,4),
              repel=TRUE,grp_title="Bout",size=3,max.overlaps=50)
psm_gg

pca_new <- data.frame(pca_prsm$scores)
pca_new$date <- prsma_sub$date

pca_new <- pca_new %>% 
  filter(Comp.1 >= -4 & Comp.1 <= 4 & Comp.2 >= -4 & Comp.2 <= 4) %>% 
  group_by(date) %>% 
  mutate(mean1 = mean(Comp.1),
         mean2 = mean(Comp.2))

date_cols <- unique(pca_new$date)



ggplot(pca_new,aes(Comp.1,Comp.2,color=date))+
  geom_point(alpha=0.01)+
  geom_point(data=pca_new %>% 
               group_by(date) %>% 
               summarise_at(c('Comp.1','Comp.2'), mean),
             size=2)+
  scale_y_continuous(limits=c(-4,4))+
  scale_x_continuous(limits=c(-4,4))+
  scale_color_manual(values=c("firebrick4","orange","yellowgreen","cornflowerblue","purple3","pink3","deeppink2"))+
  stat_ellipse(aes())+
  scale_fill_manual(values=c("firebrick4","orange","yellowgreen","cornflowerblue","purple3","pink3","deeppink2"))+
  theme_bw()


empty <- data.frame(date = seq(as.Date("01-01-2000", format = "%m-%d-%Y"), as.Date("12-31-2022", format = "%m-%d-%Y"), by = "day")) %>%
  mutate(doy = yday(date),
         year = year(date))
seasons <- read.csv('MetabolismModel/seasons_from_Robins_paper_long_withplaceholders.csv') %>%
  select(-placeholder) %>%
  mutate(date = as.Date(date,format = "%m/%d/%y"))
DATA <- left_join(empty, seasons) %>%
  fill(season, .direction="down")

DATA$season <- factor(DATA$season, levels = c("Ice-on", "Spring", "Clearwater", "Early Summer", "Late Summer", "Fall"))

new <- merge(x=pca_new,y=DATA,
             by.x='date',by.y='date') %>% 
  mutate(img = case_when(season == 'Spring' ~ 1,
                         season == 'Clearwater' ~ 1,
                         season == 'Early Summer' ~ 1,
                         season == 'Late Summer' & date == '2021-09-05' ~ 1,
                         season == 'Late Summer' & date == '2022-07-27' ~ 2,
                         season == 'Fall' & date == '2021-11-02' ~ 1,
                         season == 'Fall' & date == '2021-12-13' ~ 2))

# ggplot(new,aes(wavelength,normal,color=season,fill=season,group=date))+
#   geom_line(aes(linetype=factor(img)),size=1)+
#   scale_color_manual(values=c("#73456D",  "#B2CCF1","#EE914A", "#9AD67A","#138E90")) +
#   scale_fill_manual(values=c( "#73456D",  "#B2CCF1","#EE914A", "#9AD67A","#138E90")) +
#   facet_wrap(~season)+
#   # scale_x_log10()+
#   labs(y='Normalized Rrs',
#        x='Wavelength (nm)')+
#   theme(legend.position = 'none')

ggplot(new,aes(Comp.1,Comp.2,color=season))+
  geom_point(alpha=0.01)+
  geom_point(data=new %>%
               group_by(date,season) %>%
               summarise_at(c('Comp.1','Comp.2'), mean),
             size=2)+
  scale_y_continuous(limits=c(-4,4))+
  scale_x_continuous(limits=c(-4,4))+
  scale_color_manual(values=c("#73456D",  "#B2CCF1","#EE914A", "#9AD67A","#138E90"),name='Season') +
  stat_ellipse(aes(lty=factor(img)),size=1)+
  scale_fill_manual(values=c( "#73456D",  "#B2CCF1","#EE914A", "#9AD67A","#138E90")) +
  theme_bw()+
  # theme(legend.position='none')+
  labs(y='Component 2',
       x='Component 1')+
  guides(linetype=guide_legend(title="Image"))

  

# get variances represented by each PC ----------------------------
  variances <- as_tibble(pca_prsm$sdev^2, rownames = "component") %>% #convert SD to variance
    mutate(proportion = value/sum(value) #proportion of total variance represented
           ) #add image information in a separate column ? date = prsma_sub$date
# decide how many components to retain
  vari10 <- variances %>% 
    filter(component %in% c("Comp.1","Comp.2","Comp.3","Comp.4","Comp.5","Comp.6",
                                                 "Comp.7","Comp.8","Comp.9","Comp.10"))
  
  #variances = variances %>% 
  # filter(proportion >= 0.01) #just keep those representing at least 1% of variance
  #slice_head(n = (nrow(variances[variances$proportion >= 0.01,])+1)) #or keep those plus one
  
  # get loadings for top 10 components
  loadings = pca_prsm$loadings[,1:10] %>% #only keep the ones that have been retained
    as_tibble(rownames = "band") %>% 
    janitor::clean_names() %>%
    pivot_longer(cols = starts_with("comp")) %>% #this will cause an error if only retaining one PC - just comment out this line and fix it manually for that case
    mutate(band = as.numeric(str_remove(band, "Rrs_")), #format wavelengths for plotting
           img_date = prsma_sub$date) #tag with image date
  
# plot loadings  
loadings %>% filter(name %in%c("comp_1","comp_2","comp_3","comp_4","comp_5")) %>%
  ggplot(aes(x = band, y = value))+
    theme_bw()+
    geom_line(aes(group = name, color = name))+
    scale_color_discrete(name = "Component")+
    labs(x = "Wavelength", y = "Loading")


loadings %>% filter(name %in%c("comp_1","comp_2","comp_3","comp_4", "comp_5")) %>%
  ggplot(aes(x = band, y = value))+
  theme_bw()+
  geom_line(aes(group = name, color = name))+
  scale_color_discrete(name = "Component")+
  labs(x = "Wavelength", y = "Loading")+
  facet_wrap(~name)
# save information from this image -NOT WORKING
   vars = bind_rows(vars, variances)
   lds = bind_rows(lds, loadings)


# save information to files
write_csv(loadings, "Data/PRISMA/prisma_PCAtop10_loadings.csv")
write_csv(vari10, "Data/PRISMA/prisma_PCAtop10_variances.csv")

loadings <- read.csv("Data/PRISMA/prisma_PCAtop10_loadings_092724.csv")

#plot all PC loadings by date ------------------------------------------------------

lds = lds %>% 
  filter(img_date!="orrected.z")
lds %>% 
  mutate(grp = paste(name, img_date),
         #following line is to group the panels
         cat = case_when(name %in% c("comp_1", "comp_2") ~ "PC 1, 2",
                         name %in% c("comp_4", "comp_5") ~ "PC 4, 5",
                         name %in% c("comp_6", "comp_7") ~ "PC 6, 7",
                         T ~ "PC 3")) %>% 
  ggplot(aes(x = band, y = value))+
  geom_line(aes(group = grp, color= img_date, lty = name))+
  scale_linetype_manual(values = c(1,2,1,2,1,2,1), guide = "none")+
#  facet_wrap(~ cat)+
  facet_wrap(~name)+
  scale_color_discrete(name = "Image Date")+
  labs(x = "Wavelength", y = "Loading", title = "PRISMA: all PCs over 1% variance")+
  theme_bw()
  
ggsave("Outputs/PCA_outputs/PRISMA/loadings_by_component.png")

lds %>% 
  mutate(grp = paste(name, img_date)) %>% 
  ggplot(aes(x = band, y = value))+
  geom_line(aes(group = grp, color= name))+
 # scale_linetype_manual(values = c(1,2,1,2,1,2,1), guide = "none")+
  #  facet_wrap(~ cat)+
  facet_wrap(~img_date)+
  scale_color_discrete(name = "PC")+
  labs(x = "Wavelength", y = "Loading", title = "PRISMA: all PCs over 1% variance")+
  theme_bw()

ggsave("Outputs/PCA_outputs/PRISMA/loadings_by_date.png")

#denoise
lds %>% 
  filter(name != "comp_5",
         name != "comp_4" | img_date %in% c("2022_06_16", "2022_07_27"),
         name != "comp_6" | img_date != "2022_06_16") %>% 
  mutate(grp = paste(name, img_date)) %>% 
  ggplot(aes(x = band, y = value))+
  geom_line(aes(group = grp, color= name))+
  # scale_linetype_manual(values = c(1,2,1,2,1,2,1), guide = "none")+
  #  facet_wrap(~ cat)+
  facet_wrap(~img_date)+
  scale_color_discrete(name = "PC")+
  labs(x = "Wavelength", y = "Loading", title = "PRISMA: all PCs over 1% variance")+
  theme_bw()

ggsave("Outputs/PCA_outputs/PRISMA/loadings_bydate_denoised.png")

## scratch -- transparancy by weight
vars = read_csv("Outputs/PCA_outputs/PRISMA/prisma_PCA_variances.csv") %>% 
  mutate(component = component %>% str_replace("Comp.", "comp_")) %>% 
  dplyr::select(-value)

data = lds %>% 
  left_join(vars, by = c("name" = "component", "img_date"))

data %>% 
  mutate(grp = paste(name, img_date),
         #following line is to group the panels by shape
         cat = case_when(name %in% c("comp_1", "comp_2") ~ "PC 1, 2",
                         name %in% c("comp_4", "comp_5") ~ "PC 4, 5",
                         name %in% c("comp_6", "comp_7") ~ "PC 6, 7",
                         T ~ "PC 3"))%>% 
  ggplot(aes(x = band, y = value))+
  geom_line(aes(group = grp, color= img_date, alpha = sqrt(sqrt(proportion))))+
  facet_wrap(~ name)+
  scale_color_discrete(name = "Image Date")+
  scale_alpha(name = "Prop. Var \nExplained", range = c(0, 1))+
  labs(x = "Wavelength", y = "Loading", title = "PRISMA: all PCs at/over 1% variance")+
  theme_bw()
