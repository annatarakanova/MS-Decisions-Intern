rm(list=ls())
library(ospsuite)
library(ggplot2)
library(dplyr)
library(cowplot)
############# Read PK data were digitized from Sandhu et al. 2013 paper: 22 QD dosing events
pkdat <- read.csv("Data/Olaparib_observed_data2.csv")

ggplot(pkdat, aes(x=TIME, y= DV, col=as.factor(ID))) +
  geom_line()+
  theme_bw()  +  theme(legend.position='right', legend.direction='vertical')+
  scale_x_continuous("Time, h", limits = c(0,100)) + scale_y_log10("Drug conc., uM")



#########  Step 1:To translate data into Monolix format


pkdat_mlx <- pkdat 
pkdat_mlx_adm <-  pkdat_mlx %>% filter(!duplicated(ID))
pkdat_mlx_adm$TIME <- 0
pkdat_mlx_adm$EVID <- 1


dosing_times_h <- seq(0,2*24,by=12)
dosing_times_h 

pkdat_200bid <-pkdat_mlx_adm %>% filter(ID == "200_Gao_2023")
pkdat_200bid_adm <- data.frame()
for(i in 2:length(dosing_times_h )){
  pkdat_200bid$TIME <- dosing_times_h[i] 
  pkdat_200bid_adm <- rbind(pkdat_200bid_adm, pkdat_200bid)
  
}
pkdat_400bid <- pkdat_mlx_adm %>% filter(ID == "400_Gao_2023")
pkdat_400bid_adm <- data.frame()
for(i in 2:length(dosing_times_h )){
  pkdat_400bid$TIME <- dosing_times_h[i] 
  pkdat_400bid_adm <- rbind(pkdat_400bid_adm, pkdat_400bid)
}
pkdat_mlx_adm_sum <- rbind(pkdat_mlx_adm, pkdat_200bid_adm, pkdat_400bid_adm)

# calculate expected amount of the drug (umol) in urine at the end
pkdat_mlx_dvid2 <-  pkdat_mlx %>% filter(!duplicated(ID))
pkdat_mlx_dvid2$EVID <- 0
pkdat_mlx_dvid2$DVID <- 2

MW_olaparib <- 434
pkdat_mlx_dvid2 <- pkdat_mlx_dvid2 %>% 
    mutate(DV=if_else(ID %in% c("200_Gao_2023", "400_Gao_2023"), 
                      length(dosing_times_h)*AMT*0.15*1000/MW_olaparib, AMT*0.15*1000/MW_olaparib )) %>%
    mutate(TIME= 100, DVU="umol")
                    


###
pkdat_mlx_all <- rbind(pkdat_mlx_adm_sum, pkdat_mlx, pkdat_mlx_dvid2)
write.table(pkdat_mlx_all,"Data/Olaparib_PK_data_mlx2.csv", row.names = F, sep = ",") 

