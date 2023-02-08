#Code for 
# Municipal compost public health, waste management, and urban agriculture: 
#A decadal study of fugitive Pb in City of Boston, Massachusetts, USA 

#####Figure 2: Pb concentrations in compost from City (City of Boston) and non-city sources.#####
knitr::opts_chunk$set(echo = TRUE)
library(janitor)
library(readr)
library(tidyverse)

#load both datasets (raw)

workingdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workingdir)
getwd()

# compost_dat <- read_excel("Compost_2022_djb_2.xlsx")

non_city_raw <- read_csv("non_city_raw.csv")%>%clean_names()

city_vendor_raw <- read_csv("city_vendor_raw.csv")%>%clean_names()

#select Pb as variable of interest
non_city_label <- rep(c("Non-city", "Non-city"), times=48)

raw_Pb_non_city <- non_city_raw %>% 
  dplyr::select(sample_name, year_collected, pb_ug_g)%>% 
  na.omit() %>% 
  mutate(Source=non_city_label)

city_label <- rep(c("City", "City"), times=57)

raw_Pb_vendor <- city_vendor_raw  %>% 
  dplyr::select(sample_name, year_collected, pb_ug_g) %>% 
  na.omit()%>% 
  mutate(Source=city_label)

average_city_label <- rep(c("City"), times=9)
average_noncity_label <- rep(c("Non-city"), times=9)


#Calculate global average Pb values for City and non-city samples
global_average_city <- raw_Pb_vendor %>% 
  group_by(year_collected) %>% 
  summarise(ave=mean(pb_ug_g))%>% 
  mutate(Source=average_city_label) 
global_average_non_city <- raw_Pb_non_city %>% 
  group_by(year_collected) %>% 
  summarise(ave=mean(pb_ug_g)) %>% 
  mutate(Source=average_noncity_label)

global_average_compare <- rbind(global_average_city, global_average_non_city)

t.test(ave~Source ,data=global_average_compare)


fig1_raw <- rbind(raw_Pb_non_city, raw_Pb_vendor)

#run summary stats
fig1_summ <- fig1_raw %>% 
  group_by(year_collected, Source) %>% 
  summarise(Pb_mean=mean(pb_ug_g), Pb_sd=sd(pb_ug_g))

global_mean <- fig1_summ %>% 
  group_by(Source) %>% 
  summarise(Pb_global=mean(Pb_mean))

##FIGURE2: color coded bar graph (city vs non city )  

##Figures 2B and 2C
fig2 <- ggplot(data=fig1_raw, aes(x=year_collected, y=pb_ug_g, group=year_collected))
fig2 + geom_boxplot(aes(fill=Source))+ 
  facet_grid(. ~  Source)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size = 16), axis.text=element_text(size=12))+ 
  labs(x = "Year Collected", y= "Pb (ppm)")



######Figure 3: PCA Whole matrix geochemical characterization##### 
##Load Required Packages
library(readxl)
library(ggplot2)
library(readr)
library(tidyverse)
library(stringr)
library(corrplot)
library(Hmisc)
library(factoextra)
library(ggfortify)
library(colorspace)
library(RColorBrewer)


#Load data
workingdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workingdir)
getwd()

compost_dat <- read_excel("Compost_2022_djb_2.xlsx")


#Create dataframe with only continuous variables
compost_dat_cont <- select(compost_dat, c(Mg, Al, Si, 
                                     P, S, Cl,
                                     K, Ca, Ti, 
                                     V, Cr, Mn,
                                     Fe, Co, Ni,
                                     Cu, Zn, As, 
                                     Br, Rb, Sr,
                                     Pb))
  

#Convert all values to numeric
compost_dat_ContNum <- as.data.frame(lapply(compost_dat_cont, as.numeric))


##Calculate the Correlation Matrix (including p-values)

compost_dat_cor <- rcorr(as.matrix(compost_dat_ContNum), type="pearson")
compost_dat_cor$r # look at correlation coefficients (r)
compost_dat_cor$P # look at p-values for each pairwise comparison

#Make the correlogram!
corrplot(compost_dat_cor$r, type = "upper", 
         tl.col = "black", 
         tl.srt = 45, # tilt top row to 45-degree angle
         addCoef.col = "black",
         number.cex=0.75) 

#PCA prep
#Dataframe with both continuous and categorical variables
compost_dat

# 1. Specify categorical variables to be their own objects and remove
#    them from the data frame

#Drop NAs from the dataframe
compost_pca <- compost_dat %>%
  drop_na() 

Year_vector <- as.vector(compost_pca$Year)
Source_vector <- as.vector(compost_pca$'Waste Stream Source')

# 2. Create a dataframe with only continuous variable columns for the PCA

compost_pca <- select(compost_pca, c(Mg, Al, Si, 
                                     P, S, Cl,
                                     K, Ca, Ti, 
                                     V, Cr, Mn,
                                     Fe, Co, Ni,
                                     Cu, Zn, As, 
                                     Br, Rb, Sr,
                                     Pb))

#See if you need to log transform the data

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Mg))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Al))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Si))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = P))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = S))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Cl))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = K))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Ca))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Ti))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = V))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Cr))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Mn))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Fe))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Co))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Ni))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Cu))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Zn))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = As))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Br))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Rb))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Sr))

ggplot(data = compost_pca) + 
  geom_freqpoly(aes(x = Pb))

#OK without transforming 
# Estimate the PCA Model

compost_pca_est <- prcomp(compost_pca, center= TRUE, scale.=TRUE)

# Look at summary output of the principal components

summary(compost_pca_est)

# 6. Summary plot showing percent variation in each principal component

fviz_eig(compost_pca_est)

# PCA plot mapping variables in relation to each other on PC1 & PC2
fviz_pca_var(compost_pca_est,
             geom.ind = "point", # show points only (but not "text")
             mean.point = FALSE, # Remove point that represents the mean of each group
             addEllipses = FALSE, # add ellipses
             col.var = "black") + # make variables & arrows black (default is blue)
  theme_bw()

# You can plot PCA also grouped based on one of your categorical variables
fviz_pca_biplot(compost_pca_est,
                geom.ind = "point", # show points only (but not "text")
                labelsize = 8,
                pointsize = 3,
                col.ind = Source_vector, # color by categorical variable
                mean.point = FALSE, # Remove point that represents the mean of each group
                addEllipses = TRUE, # add ellipses
                col.var = "black", # make variables & arrows black (default is blue)
                legend.title = "Waste Stream Source")  +
  theme_bw(base_size = 20)


fviz_pca_biplot(compost_pca_est,
                geom.ind = "point",
                col.ind = Year_vector,
                col.var = "black",
                mean.point = FALSE, # Remove point that represents the mean of each group
                addEllipses = TRUE, # add ellipses
                gradient.cols = "Spectral",
                legend.title = "Year")  +
  theme_bw()

fviz_pca_biplot(compost_pca_est,
                geom.ind = "point",
                col.ind = Year_vector,
                col.var = "black",
                mean.point = FALSE, # Remove point that represents the mean of each group
                addEllipses = FALSE, # add ellipses
                gradient.cols = "Spectral",
                legend.title = "Year")  +
  theme_bw()



#####Figure 4: Geochemical fingerprinting by comparing: #####
#A) City of Boston vendor sources versus non-city sources, 
#B) Change in Pb/Ti and Si/Ca ratios over time for city compost and 
#C) Waste streams of non-urban source compost.

###Figure 4a
##select Si and Ca as variables of interest from dataframes used to create Figure 2
si_ca_non <- non_city_raw %>% 
  dplyr::select(year_collected, si_percent, ca_percent, pb_ug_g, ti_ug_g, loi)%>% 
  mutate(si_ca_ratio=si_percent/ca_percent, pb_ti_ratio=pb_ug_g/ti_ug_g)%>%  
  mutate(Source=non_city_label) %>% 
  na.omit

city_label_si<- rep(c("City"), times=133)

si_ca_city <- city_vendor_raw %>% 
  dplyr::select(year_collected, si_percent, ca_percent, pb_ug_g, ti_ug_g, loi)%>%
  mutate(si_ca_ratio=si_percent/ca_percent, pb_ti_ratio=pb_ug_g/ti_ug_g)%>% 
  mutate(Source=city_label_si) %>% 
  na.omit()

#combine dataframes 
si_ca <- rbind(si_ca_non, si_ca_city) %>% 
  mutate(si_ca_ratio=si_percent/ca_percent, pb_ti_ratio=pb_ug_g/ti_ug_g)

#Si/Ca Plot
gradient <- c("#99999", "#E69F00")

fig4a <- ggplot(data=si_ca, aes(x=si_ca_ratio, y=pb_ti_ratio, color=Source, shape=Source))
fig4a + geom_jitter()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  labs(x = "Si/Ca", y= "Pb/Ti") + 
  scale_fill_manual(values=gradient)

###Figure 4b
#expand color range (two color gradient)
#City ratio plots by year
library(RColorBrewer)
library(scales)

colnames(si_ca_city)[1] <- "Year Collected"

ggplot()+ geom_point(data=si_ca_city, aes(x=si_ca_ratio, y=pb_ti_ratio, colour=`Year Collected`))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size = 18), axis.text=element_text(size=14))+ labs(x = "Si/Ca ", y= "Pb/Ti", fill="Year Collected")+
  scale_colour_gradientn(colors=c("#B2182B", "#EF8A62", "#FDDBC7", "#F7F7F7" ,"#D1E5F0", "#67A9CF","#2166AC"), values = rescale(si_ca_city$year_collected, c(2006, 2008, 2009, 2010, 2011, 2021, 2022)), guide = "colorbar", limits=c(2006, 2022))

###Figure 4c
si_ca_non_BINN <- non_city_raw %>% 
  dplyr::select(year_collected, waste_stream_source, si_percent, ca_percent, pb_ug_g, ti_ug_g, loi)%>% 
  mutate(si_ca_ratio=si_percent/ca_percent, pb_ti_ratio=pb_ug_g/ti_ug_g)%>%  
  mutate(Source=non_city_label) %>% 
  na.omit

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

legend_title <- "Waste Stream Source"
fig4c <- ggplot(data=si_ca_non_BINN, aes(x=si_ca_ratio, y=pb_ti_ratio, shape=waste_stream_source, color=waste_stream_source))+ geom_jitter()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size = 18), axis.text=element_text(size=14))+ labs(x = "Si/Ca", y= "Pb/Ti")+ scale_y_continuous(limits = c(0, 0.3)) + scale_x_continuous(limits=c(0,25)) + guides(shape = "none") + scale_shape_manual(values=c(2, 8, 1)) + scale_color_discrete(legend_title)

fig4c



#####Figure 5: Pb Concentrations in Size Fractionated Compost #####
#(A) from 2006 to 2022 
#(B) in City and Suburban settings in 2022

compost_size_dat <- read_excel("Compost_2022_size_for_annual_bar_graph.xlsx")

compost_elements_vs_size <- compost_size_dat %>%
  select('Year collected', 'Size (microns)', Pb) %>%
  rename(size = 'Size (microns)',
         year = 'Year collected')

compost_elements_vs_size$size <- factor(compost_elements_vs_size$size, levels = 
                                     c(">2000", "2000-125", "125-63", 
                                       "63-37", "<37"))

###Figure 5a
ggplot(data = compost_elements_vs_size, x = size, y = Pb, group = year)+
  geom_bar(aes(x = size, y = Pb), stat = "identity", width = 0.5, fill = "#080707")+
  facet_wrap(~year)+ 
  theme_bw(base_size = 16,) +
  labs(y = "Pb Concentration (ug/g)", x = "Grain Size") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 700))


#Suburban vs city data from 2021

compost_size_dat_22 <- read_csv("C:/Users/Claire/Desktop/Wellesley/Compost MS 2023/2022 Compost size.csv")


compost_size_dat_22$Size <- factor(compost_size_dat_22$Size, levels = 
                                          c("BULK", ">2000", "2000-250", "250-125", "125-63", 
                                            "63-37", "<37"))

###Figure 5b
ggplot(data = compost_size_dat_22, x = Size, y = Pb, group = Source)+
  geom_bar(aes(x = Size, y = Pb), stat = "identity", width = 0.5, fill = "#080707")+
  scale_fill_manual(values = c("#808080", "#080707", "#080707", "#080707", "#080707", "#080707", "#080707"))+
  facet_grid(~Source)+ 
  theme_bw(base_size = 16,) +
  labs(y = "Pb Concentration (ug/g)", x = "Grain Size") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 400))









######Grain Size PCA##########

#Load data
workingdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workingdir)
getwd()

compost_size_dat <- read_xlsx("Compost size for Pca.xlsx")



#Create dataframe with only continuous variables
compost_size_dat_cont <- select(compost_size_dat, c(Mg, Al, Si, 
                                                    P, S, Cl,
                                                    K, Ca, Ti, 
                                                    Cr, Mn,
                                                    Fe, Ni,
                                                    Cu, Zn, 
                                                    Br, Rb, Sr,
                                                    Pb))

#to remove use -c(...)
#Convert all values to numeric
compost_size_dat_ContNum <- as.data.frame(lapply(compost_size_dat_cont, as.numeric))


##Step 1: Calculate the Correlation Matrix (including p-values)

compost_size_dat_cor <- rcorr(as.matrix(compost_size_dat_ContNum), type="pearson")
compost_size_dat_cor$r # look at correlation coefficients (r)
compost_size_dat_cor$P # look at p-values for each pairwise comparison

#Make the correlogram!
corrplot(compost_size_dat_cor$r, type = "upper", 
         tl.col = "black", 
         tl.srt = 45, # tilt top row to 45-degree angle
         addCoef.col = "black",
         number.cex=0.75) 

#PCA time bae-bee
#Dataframe with both contin and categorical variables
compost_size_dat

# 1. Specify categorical variables to be their own objects and remove
#    them from the data frame

#Drop NAs from the dataframe
compost_size_pca <- compost_size_dat %>%
  drop_na() 

Size_vector <- as.vector(compost_size_pca$"Size (microns)")


# 2. Create a dataframe with only continuous variable columns for the PCA

compost_size_pca <- select(compost_size_pca, c(Mg, Al, Si, 
                                               P, S, Cl,
                                               K, Ca, Ti, 
                                               Cr, Mn,
                                               Fe, Ni,
                                               Cu, Zn, 
                                               Br, Rb, Sr,
                                               Pb))

#See if you need to log transform the data lmao

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Mg))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Al))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Si))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = P))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = S))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Cl))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = K))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Ca))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Ti))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = V))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Cr))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Mn))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Fe))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Co))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Ni))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Cu))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Zn))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = As))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Br))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Rb))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Sr))

ggplot(data = compost_size_pca) + 
  geom_freqpoly(aes(x = Pb))

#OK without transforming - here is code if we would like to from Tar Creek code
#Log transform the data for the elements that are skewed
# TarCreek_pcaLog <- TarCreek_pca %>%
#   mutate(Ag_log = log(Ag), Ga_log = log(Ga),
#          Pb_log = log(Pb), Si_log = log(Si),
#          As_log = log(As), Cd_log = log(Cd),
#          Ba_log = log(Ba), Pd_log = log(Pd),
#          Sn_log = log(Sn), Ge_log = log(Ge),
#          Zn_log = log(Zn), Fe_log = log(Fe),
#          Ca_log = log(Ca), S_log = log(S),
#          Se_log = log(Se),
#   )

# Estimate the PCA Model

compost_size_pca_est <- prcomp(compost_size_pca, center= TRUE, scale.=TRUE)

# Look at summary output of the principal components

summary(compost_size_pca_est)

# 6. Summary plot showing percent variation in each principal component

fviz_eig(compost_size_pca_est)

# PCA plot mapping variables in relation to each other on PC1 & PC2
fviz_pca_var(compost_size_pca_est,
             geom.ind = "point", # show points only (but not "text")
             mean.point = FALSE, # Remove point that represents the mean of each group
             addEllipses = FALSE, # add ellipses
             col.var = "black") + # make variables & arrows black (default is blue)
  theme_bw()

# You can plot PCA also grouped based on one of your categorical variables
fviz_pca_biplot(compost_size_pca_est,
                geom.ind = "point", # show points only (but not "text")
                labelsize = 8,
                pointsize = 3,
                col.ind = Size_vector, # color by categorical variable
                mean.point = FALSE, # Remove point that represents the mean of each group
                addEllipses = TRUE, # add ellipses
                col.var = "black", # make variables & arrows black (default is blue)
                legend.title = "Size")  +
  theme_bw(base_size = 20)
