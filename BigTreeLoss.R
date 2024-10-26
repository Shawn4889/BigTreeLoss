#Chapater SAR biomass mapping
#Xiaoxuan Li 
library(lidR)
library(tools) 
library(raster)
library(ggplot2)
library(ggpointdensity)
library(viridis)
library(ggcorrplot)
library(ggpmisc)
library(dplyr)
library(tidyverse) 
library(rsq)
library(readxl)
library(grid)
library(mltools)
library(Metrics)
library(devtools)
library(ggsignif)
library(igraph)
library(ff)
library(corrplot)
library(psych)
library(xlsx)
library(car)
library(MASS)
library(ggrepel)
library(rlist)
library(data.table)
library(lemon)
library(lmtest)
library(reshape2)
library(ggpubr)
library(tidyr)
library(stats)
library(igraph)
library(caret)
library(vegan)
library(rasterdiv)
library(snow)

#The following codes are for LiDAR batch processing and big tree loss detection in South Africa
# ------------------------------------------------------------------------------------------------ #
ras_2018_dir <- "E:\\ChangMap\\CHM\\DB_chm\\Noshurb_pitfree\\Site\\CHM_2008_Justicia.tif"
ras <- raster(ras_2018_dir)
values(ras)[values(ras) < 5] = NA
#ras_smooth <- focal(ras, w=matrix(1,3,3), fun=mean,na.rm=TRUE)

#f_tree <- function(x) {ifelse(x < 10, x + 3, x + 3)}
#f_tree <- function(x) {x+5}

#tt <- find_trees(ras, lmf(f_tree, 5, shape = "circular"))
#tt <- find_trees(ras, lmf(f_tree, 9, shape = "circular"))
tt <- find_trees(ras, lmf(ws = 15, hmin = 7, shape = "circular"))

#plot(tt,add=T)
#save detected trees 
shapefile(tt, filename='E:\\Ghermay\\04212024\\tt_2008_Justicia.shp',overwrite=TRUE)
crowns = dalponte2016(ras, tt, th_tree = 5, th_seed = 0.7, th_cr = 0.4, max_cr = 15)()
writeRaster(crowns, "E:\\Ghermay\\04212024\\Crowns_2008_Justicia.tif",datatype="INT2U",overwrite=TRUE)



# ------------------------------------------------------------------------------------------------ #
#batch preprocessing lidar point cloud to chm
dir <- "E:\\NLT\\DC_Tree\\06232024\\LiDAR\\ALS\\2020"         
output = "E:\\NLT\\DC_Tree\\06232024\\LiDAR\\CHM\\2020"
las_files = list.files(dir,pattern = glob2rx('*.las'),full.names = T)
for (i in las_files){
  out = file.path(output,paste0(basename(file_path_sans_ext(i)), ".tif"))
  if (!file.exists(out)) 
  {
    print(i)
    print("read ALS")
    las <- readLAS(i,select = "i,c,r")
    print("normalize ALS")
    las = normalize_height(las, tin(),na.rm = TRUE)
    print("filter out ALS point cloud > 100m")
    las <- filter_poi(las, Z >= 0, Z <= 100)
    print("1m CHM generation")
    chm <- grid_canopy(las, 1, p2r())
    print("Output chm")
    writeRaster(chm,out,options=c('TFW=YES'),overwrite=TRUE)
  }
}


# ------------------------------------------------------------------------------------------------ #
#04212024 - batch process
#tree local maximum polygons and tree crown raster
ras_dir <- "E:\\NLT\\DC_Tree\\06232024\\LiDAR\\CHM\\2020"
out_dir <- "E:\\NLT\\DC_Tree\\06232024\\LiDAR\\CHM"       
ras_files <- list.files(ras_dir, pattern = "\\.tif$", full.names = TRUE)
for (i in ras_files){
  print(i)
  ras <- raster(i)
  values(ras)[values(ras) < 15] = NA
  tt <- find_trees(ras, lmf(ws = 15, hmin = 7, shape = "circular"))
  shapefile(tt, 
            filename=file.path(out_dir,paste0(basename(file_path_sans_ext(i)),"_tt.shp")),
            overwrite=TRUE)
  crowns = dalponte2016(ras, tt, th_tree = 5, th_seed = 0.7, th_cr = 0.4, max_cr = 15)()
  writeRaster(crowns, 
              file.path(out_dir,paste0(basename(file_path_sans_ext(i)),"_crown.tif")),
              datatype="INT2U",
              overwrite=TRUE)}


#python focal stats to do tree crown pit-filling 
#zonal stats using 2008/2010 local maxima and 2018 tree crown 


# ------------------------------------------------------------------------------------------------ #
#batch process zonal stats results and do difference
sur_list_1 <- c("Agincourt_2008","Ireagh_1_2008","Ireagh_2_2008","Andover_2010",
                "Welverdiendt_2010","Justicia_2008","Justicia_2010")
sur_list_2 <- c("Agincourt_2018","Ireagh_1_2018","Ireagh_2_2018","Andover_2018",
                "Welverdiendt_2018","Justicia_2018","Justicia_2018")
sur_list_3 <- c("Agincourt_2018_2008","Ireagh_1_2018_2008","Ireagh_2_2018_2008","Andover_2018_2010",
                "Welverdiendt_2018_2010","Justicia_2018_2008","Justicia_2018_2010")

dir <- "E:\\NLT\\Ghermay\\04212024\\TT_zonal_table_mean"

for (i in 1:length(sur_list_1)){
  print(sur_list_1[i])
  print(sur_list_2[i])
  dir_1 <- file.path(dir,paste0(sur_list_1[i],"_Zonal.csv"))
  dir_2 <- file.path(dir,paste0(sur_list_2[i],"_Zonal.csv"))
  Data_1 = read.csv(dir_1)
  Data_2 = read.csv(dir_2)
  Data_1 <- Data_1[Data_1$MEAN > 7,]
  Data_index <- merge(Data_1, Data_2, by.x = "FID", by.y = "FID")
  Data_index$Diff_percent <- (Data_index$MEAN.y-Data_index$MEAN.x)/Data_index$MEAN.x*100
  write.csv(Data_index, file.path(dir,paste0(sur_list_3[i],".csv")), row.names=FALSE)
}




# ------------------------------------------------------------------------------------------------ #
#04212024 - batch process
dir <- "E:\\NLT\\Ghermay\\05042024"
sur_list_3 <- c("Agincourt_2018_2008", "Ireagh_Sabi_2018_2008", "Ireagh_Comm_2018_2008",
                "Andover_2018_2010", "Justicia_Sabi_2018_2008", "Justicia_Comm_2018_2008",
                "Welverdiendt_2018_2010","Justicia_Sabi_2018_2010", "Justicia_Comm_2018_2010")
plot_list <- list()
for (i in 1:length(sur_list_3)){
  dir_csv <- file.path(dir,paste0(sur_list_3[i],".csv"))
  Data_index = read.csv(dir_csv)
  Data_index <- na.omit(Data_index)
  p <- ggplot(Data_index, aes(x=Diff_percent)) +
    geom_histogram(binwidth = 5, color="black", fill="#9AC8CD") +
    theme_minimal()+
    labs(title = paste0("Crown feature, p98, ", basename(file_path_sans_ext(dir_csv))), y="Frequency", 
         x="Big tree loss/gain (%) between 2018 and 2008/2010 ")+
    theme(legend.position = "none")+
    theme(legend.title = element_blank())+
    coord_cartesian(xlim =  c(-100, 100))+
    theme(text=element_text(size=25))+
    theme(plot.title = element_text(hjust = 0.5))
  plot_list[[i]] <- p
  print(sur_list_3[i])
  #print(summary(Data_index$Diff_percent))
  print(sum(Data_index$Diff_percent <= 0)/nrow(Data_index)*100)
  print(sum(Data_index$Diff_percent <= -60)/nrow(Data_index)*100)
  print(sum(Data_index$Diff_percent <= -75)/nrow(Data_index)*100)
}


ggarrange(plotlist = plot_list, ncol = 3, nrow = 3)
out = "E:\\NLT\\Ghermay\\05042024\\Crown_zonal_table_p98.jpg"
ggsave(out,height=15, width=30, dpi=600)






#The following codes are for DC tree detection
# ------------------------------------------------------------------------------------------------ #
#04212024 - batch process - DC
#tree local maximum polygons and tree crown raster
ras_dir <- "E:\\NLT\\DC_Tree\\06232024\\LiDAR\\CHM\\DC_CHM_2015.tif"
out_dir <- "E:\\NLT\\DC_Tree\\06232024\\LiDAR\\CHM"       

ras <- raster(ras_dir)

values(ras)[values(ras) < 15] = NA
tt <- find_trees(ras, lmf(ws = 20, hmin = 15, shape = "circular"))

crowns = dalponte2016(ras, tt, th_tree = 12, th_seed = 0.7, th_cr = 0.4, max_cr = 25)()

shapefile(tt, 
          filename=file.path(out_dir,"Local_maxima.shp"),
          overwrite=TRUE)

writeRaster(crowns, 
            file.path(out_dir,"Tree_crown.tif"),
            datatype="INT2U",
            overwrite=TRUE)

# ------------------------------------------------------------------------------------------------ #
#batch process zonal stats results and do difference -DC
dir <- "E:\\NLT\\DC_Tree\\07302024\\Result"
dir_1 <- file.path(dir,"DC_CHM_2015_zonal.csv")
dir_2 <- file.path(dir,"DC_CHM_2020_zonal.csv")
Data_1 = read.csv(dir_1)
Data_2 = read.csv(dir_2)
#Data_1 <- Data_1[Data_1$MEAN > 7,]
Data_index <- merge(Data_1, Data_2, by.x = "Id", by.y = "Id")
Data_index$Diff_percent <- (Data_index$PCT98.y-Data_index$PCT98.x)/Data_index$PCT98.x*100
write.csv(Data_index, file.path(dir,"Merge_2015_2020.csv"), row.names=FALSE)


# ------------------------------------------------------------------------------------------------ #
#04212024 - batch process
dir <- "E:\\NLT\\DC_Tree\\07302024\\Result\\Merge_2015_2020.csv"
Data_index = read.csv(dir)
Data_index <- na.omit(Data_index)
ggplot(Data_index, aes(x=Diff_percent)) +
  geom_histogram(binwidth = 5, color="black", fill="#9AC8CD") +
  theme_minimal()+
  labs(title = paste0("Big tree change in Washington D.C. ward 7"), y="Frequency", 
       x="Big tree loss/gain (%) between 2015 and 2020")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  coord_cartesian(xlim =  c(-100, 100))+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))





# ------------------------------------------------------------------------------------------------ #
#07042024 - batch process - DC
#tree local maximum polygons and tree crown raster
ras_dir <- "E:\\NLT\\DC_Tree\\06232024\\LiDAR\\CHM\\DC_CHM_2015.tif"
out_dir <- "E:\\NLT\\DC_Tree\\07032024\\Crown"       

ras <- raster(ras_dir)

values(ras)[values(ras) < 5] = NA
tt <- find_trees(ras, lmf(ws = 12, hmin = 5, shape = "circular"))
crowns = dalponte2016(ras, tt, th_tree = 5, th_seed = 0.7, th_cr = 0.4, max_cr = 12)()

shapefile(tt, 
          filename=file.path(out_dir,"Local_maxima.shp"),
          overwrite=TRUE)

writeRaster(crowns, 
            file.path(out_dir,"Tree_crown.tif"),
            datatype="INT2U",
            overwrite=TRUE)



# ------------------------------------------------------------------------------------------------ #
#07042024 - batch process
dir <- "E:\\NLT\\DC_Tree\\07302024\\Result\\stats_2015_2020_5m.csv"
Data = read.csv(dir)

reg1 <- lm(PCT98_2015 ~ Ref_HEIGHT_m, data = Data)
r2_1 <- round(summary(reg1)$adj.r.squared,2)
MD_1 <- mean(Data$PCT98_2015 - Data$Ref_HEIGHT_m)

reg2 <- lm(PCT98_2020 ~ Ref_HEIGHT_m, data = Data)
r2_2 <- round(summary(reg2)$adj.r.squared,2)
MD_2 <- mean(Data$PCT98_2020 - Data$Ref_HEIGHT_m)


ggplot(Data)+ 
  theme_bw()+
  coord_cartesian(xlim = c(5, 30),ylim =  c(5, 30))+
  geom_point(aes(x = Ref_HEIGHT_m, y = PCT98_2015), color = 'blue') + 
  geom_point(aes(x = Ref_HEIGHT_m, y = PCT98_2020), color = 'red') +  
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+ 
  annotate("text", x=17, y=30, 
           label = paste("R2 of 2015: ",round(r2_1,3)), color = "blue", size = 8) +
  annotate("text", x=17, y=29, 
           label = paste("R2 of 2020: ",round(r2_2,3)), color = "red", size = 8)+
  annotate("text", x=17, y=28, 
           label = "n = 8392", color = "black", size = 8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(size=35)) +
  labs(x="Reference height (m)", 
       y="LiDAR-based height (m)")+
  theme(plot.title = element_text(hjust = 0.5))

out = "E:\\NLT\\DC_Tree\\07302024\\LiDAR_2015_2020_5m.jpg"
ggsave(out,height=15, width=15, dpi=600)



#histogram of height change 
dir <- "E:\\NLT\\DC_Tree\\07302024\\Result\\stats_2015_2020_5m.csv"
Data = read.csv(dir)

Data$Diff_2023_2015 <- (Data$Ref_HEIGHT_m - Data$PCT98_2015)/Data$PCT98_2015*100
Data$Diff_2023_2020 <- (Data$Ref_HEIGHT_m - Data$PCT98_2020)/Data$PCT98_2020*100

Data_sel <- Data[, c("Id", "Diff_2023_2015", "Diff_2023_2020")]
Data_melt <- melt(Data_sel,id = "Id", variable.name = "Year", value.name = "Value")

ggplot(Data_melt, aes(x = Value, fill = Year)) +
  geom_boxplot() +
  theme_minimal()+
  labs(title = paste0("Big tree change in Washington, D.C. ward 7, n = 8392"), y="Frequency", 
       x="Big tree loss/gain (%) between 2015/2020 and 2023")+
  theme(legend.position = c(0.2,0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(xlim =  c(-100, 100))+
  scale_y_continuous("Frequency (%)", labels=scales::percent) +
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))





# ------------------------------------------------------------------------------------------------ #
#07182024 - 3m DC known trees vs. CHM 
dir <- "E:\\NLT\\DC_Tree\\07182024\\DC_trees_2015_2020_3m.csv"
Data = read.csv(dir)
Data <- Data[Data$HEIGHT_2015_p98 >= 2,]
Data <- Data[Data$HEIGHT_2020_p98 >= 2,]
Data <- Data[Data$HEIGHT_2023_ref>= 2,]


reg1 <- lm(HEIGHT_2015_p98 ~ HEIGHT_2023_ref, data = Data)
r2_1 <- round(summary(reg1)$adj.r.squared,2)
MD_1 <- mean(Data$HEIGHT_2015_p98 - Data$HEIGHT_2023_ref)

reg2 <- lm(HEIGHT_2020_p98 ~ HEIGHT_2023_ref, data = Data)
r2_2 <- round(summary(reg2)$adj.r.squared,2)
MD_2 <- mean(Data$HEIGHT_2020_p98 - Data$HEIGHT_2023_ref)

x_axis <- c(0,5,10,15,20,25,30)


p1 <- ggplot(Data)+ 
  theme_bw()+
  coord_cartesian(xlim = c(2, 30),ylim =  c(2, 30))+
  geom_point(aes(x = HEIGHT_2023_ref, y = HEIGHT_2015_p98), color = 'red') + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+ 
  annotate("text", x=5, y=30, 
           label = paste("R2 of 2015 LiDAR: ",round(r2_1,3)), color = "black", size = 12, hjust = 0) +
  annotate("text", x=5, y=28, 
           label = paste("Bias of 2015 LiDAR: ",round(MD_1,3)), color = "black", size = 12, hjust = 0)+
  annotate("text", x=5, y=26, 
           label = "n = 31934", color = "black", size = 12, hjust = 0)+
  scale_y_continuous(minor_breaks = x_axis,breaks =  x_axis)+
  scale_x_continuous(minor_breaks = round(x_axis,digits = 1),
                     breaks = round(x_axis,digits = 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(size=40)) +
  labs(x="2023 DC tree reference height (m)", 
       y="2015 LiDAR-based height (m)")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(Data)+ 
  theme_bw()+
  coord_cartesian(xlim = c(2, 30),ylim =  c(2, 30))+
  geom_point(aes(x = HEIGHT_2023_ref, y = HEIGHT_2020_p98), color = 'blue') +  
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+ 
  annotate("text", x=5, y=30, 
           label = paste("R2 of 2020 LiDAR: ",round(r2_2,3)), color = "black", size = 12, hjust = 0) +
  annotate("text", x=5, y=28, 
           label = paste("Bias of 2020 LiDAR: ",round(MD_2,3)), color = "black", size = 12, hjust = 0)+
  annotate("text", x=5, y=26, 
           label = "n = 31934", color = "black", size = 12, hjust = 0)+
  scale_y_continuous(minor_breaks = x_axis,breaks =  x_axis)+
  scale_x_continuous(minor_breaks = round(x_axis,digits = 1),
                     breaks = round(x_axis,digits = 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(size=40)) +
  labs(x="2023 DC tree reference height (m)", 
       y="2020 LiDAR-based height (m)")+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(p1,p2)
out = "E:\\NLT\\DC_Tree\\08082024\\DC_trees_2015_2020_3m.jpg"
ggsave(out,height=10, width=20, dpi=600)




#histogram of height change 
dir <- "E:\\NLT\\DC_Tree\\07182024\\DC_trees_2015_2020_3m.csv"
Data = read.csv(dir)
Data <- Data[Data$HEIGHT_2015_p98 >= 2,]
Data <- Data[Data$HEIGHT_2020_p98 >= 2,]
Data <- Data[Data$HEIGHT_2023_ref>= 2,]


Data_sel <- Data[, c("ID", "Diff_2023_2015", "Diff_2023_2020")]
Data_melt <- melt(Data_sel,id = "ID", variable.name = "Year", value.name = "Value")

ggplot(Data_melt, aes(x = Value, fill = Year)) +
  geom_density(alpha = 0.5) +
  theme_minimal()+
  labs(title = paste0("Big tree change in Washington, D.C. ward 7, n = 31932"), y="Frequency", 
       x="Big tree loss/gain (%) between 2023 and 2015/2020")+
  theme(legend.position = c(0.2,0.7))+
  theme(legend.title = element_blank())+
  coord_cartesian(xlim =  c(-100, 100))+
  scale_y_continuous("Frequency (%)", labels=scales::percent) +
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))

out = "E:\\NLT\\DC_Tree\\08082024\\DC_trees_2015_2020_3m_histogram.jpg"
ggsave(out,height=10, width=20, dpi=600)


# ------------------------------------------------------------------------------------------------ #
#08082024 - batch process - DC
#tree local maximum polygons and tree crown raster
ras_dir <- "E:\\NLT\\DC_Tree\\06232024\\LiDAR\\CHM\\DC_CHM_2015.tif"
out_dir <- "E:\\NLT\\DC_Tree\\08082024\\Crown"       

ras <- raster(ras_dir)

values(ras)[values(ras) < 2] = NA
tt <- find_trees(ras, lmf(ws = 8, hmin = 2, shape = "circular"))
crowns = dalponte2016(ras, tt, th_tree = 2, th_seed = 0.7, th_cr = 0.5, max_cr = 4)()

shapefile(tt, 
          filename=file.path(out_dir,"Local_maxima_2015_2m.shp"),
          overwrite=TRUE)

writeRaster(crowns, 
            file.path(out_dir,"Tree_crown_2015_2m.tif"),
            datatype="INT4U",
            overwrite=TRUE)



# ------------------------------------------------------------------------------------------------ #
#08082024 - batch process
dir <- "E:\\NLT\\DC_Tree\\08082024\\Result\\DC_CHM_2015_zonal_3m.csv"
Data1 = read.csv(dir)
Data1 <- Data1[Data1$Ref_height >= 2,]
Data1 <- Data1[Data1$CHM_2015 >= 2,]
reg1 <- lm(CHM_2015 ~ Ref_height, data = Data1)
r2_1 <- round(summary(reg1)$adj.r.squared,2)
MD_1 <- mean(Data1$CHM_2015 - Data1$Ref_height)

dir <- "E:\\NLT\\DC_Tree\\08082024\\Result\\DC_CHM_2020_zonal_3m.csv"
Data2 = read.csv(dir)
Data2 <- Data2[Data2$Ref_height >= 2,]
Data2 <- Data2[Data2$CHM_2020 >= 2,]
reg2 <- lm(CHM_2020 ~ Ref_height, data = Data2)
r2_2 <- round(summary(reg2)$adj.r.squared,2)
MD_2 <- mean(Data2$CHM_2020 - Data2$Ref_height)

x_axis <- c(0,5,10,15,20,25,30)

p1<- ggplot(Data1)+ 
  theme_bw()+
  coord_cartesian(xlim = c(2, 30),ylim =  c(2, 30))+
  geom_point(aes(x = Ref_height, y = CHM_2015), color = 'red') + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+ 
  annotate("text", x=6, y=30, 
           label = paste("R2 of 2015 LiDAR: ",round(r2_1,3)), color = "black", size = 12, hjust = 0) +
  annotate("text", x=6, y=28, 
           label = paste("Bias of 2015 LiDAR: ",round(MD_1,3)), color = "black", size = 12, hjust = 0) +
  annotate("text", x=6, y=26, 
           label = "n = 10224 (32%)", color = "black", size = 12, hjust = 0)+
  scale_y_continuous(minor_breaks = x_axis,breaks =  x_axis)+
  scale_x_continuous(minor_breaks = round(x_axis,digits = 1),
                     breaks = round(x_axis,digits = 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(size=40)) +
  labs(x="2023 DC tree reference height (m)", 
       y="2015 LiDAR-based height (m)")+
  theme(plot.title = element_text(hjust = 0.5))

p2<- ggplot(Data2)+ 
  theme_bw()+
  coord_cartesian(xlim = c(2, 30),ylim =  c(2, 30))+
  geom_point(aes(x = Ref_height, y = CHM_2020), color = 'blue') +  
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+ 
  annotate("text", x=6, y=30, 
           label = paste("R2 of 2020 LiDAR: ",round(r2_2,3)), color = "black", size = 12, hjust = 0)+
  annotate("text", x=6, y=28, 
           label = paste("Bias of 2020 LiDAR: ",round(MD_2,3)), color = "black", size = 12, hjust = 0) +
  annotate("text", x=6, y=26, 
           label = "n = 15200 (47.6%)", color = "black", size = 12, hjust = 0)+
  scale_y_continuous(minor_breaks = x_axis,breaks =  x_axis)+
  scale_x_continuous(minor_breaks = round(x_axis,digits = 1),
                     breaks = round(x_axis,digits = 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(size=40)) +
  labs(x="2023 DC tree reference height (m)", 
       y="2020 LiDAR-based height (m)")+
  theme(plot.title = element_text(hjust = 0.5))


ggarrange(p1,p2)
out = "E:\\NLT\\DC_Tree\\08082024\\LiDAR_2015_2020_2m.jpg"
ggsave(out,height=10, width=20, dpi=600)




#histogram of height change 
dir <- "E:\\NLT\\DC_Tree\\08082024\\Result\\DC_CHM_2015_zonal_3m.csv"
Data1 = read.csv(dir)
Data1 <- Data1[Data1$Ref_height >= 2,]
Data1 <- Data1[Data1$CHM_2015 >= 2,]
Data1$Diff_2023_2015 <- (Data1$Ref_height - Data1$CHM_2015)/Data1$CHM_2015*100

dir <- "E:\\NLT\\DC_Tree\\08082024\\Result\\DC_CHM_2020_zonal_3m.csv"
Data2 = read.csv(dir)
Data2 <- Data2[Data2$Ref_height >= 2,]
Data2 <- Data2[Data2$CHM_2020 >= 2,]
Data2$Diff_2023_2020 <- (Data2$Ref_height - Data2$CHM_2020)/Data2$CHM_2020*100

breaks <- seq(-80,80,10)
Data1 <- Data1 %>% group_by(gr=cut(Diff_2023_2015, breaks= breaks)) %>% 
  arrange(as.numeric(gr))
Data1 <- na.omit(Data1)
num1 <- table(Data1$gr)/nrow(Data1)*100 
names_list1 <- names(num1)
List1 <- cbind(names_list1,round(as.numeric(num1),3))

Data2 <- Data2 %>% group_by(gr=cut(Diff_2023_2020, breaks= breaks)) %>% 
  arrange(as.numeric(gr))
Data2 <- na.omit(Data2)
num2 <- table(Data2$gr)/nrow(Data2)*100
List2 <- cbind(List1,round(as.numeric(num2),3))
colnames(List2) <- c("range", "Diff_2023_2015", "Diff_2023_2020")


Data_melt <- melt(data.frame(List2),id = "range")
Data_melt$value <- as.numeric(Data_melt$value)

levels_order <- c("(-100,-90]", "(-90,-80]", "(-80,-70]", "(-70,-60]", "(-60,-50]",
                  "(-50,-40]", "(-40,-30]", "(-30,-20]", "(-20,-10]", "(-10,0]",
                  "(0,10]", "(10,20]", "(20,30]", "(30,40]", "(40,50]", "(50,60]",
                  "(60,70]", "(70,80]", "(80,90]", "(90,100]")
Data_melt$range <- factor(Data_melt$range, levels = levels_order)

colors <- c("Diff_2023_2015" = "#BB5566", "Diff_2023_2020" = "#004488")
ggplot(Data_melt, aes(x = range, y = value, fill = variable)) +
  theme_bw()+
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = colors, 
                    labels = c("Diff_2023_2015","Diff_2023_2020")) +
  labs(
    x = "Difference height between 2023 and 2015/2020 (Mg/ha)",
    y = "Frequency (%)"
  ) +
  theme(legend.position=c(.1,.8))+
  theme(text=element_text(size=40, color = "black"))+
  theme(text=element_text(size=40, color = "black"),
        axis.text=element_text(size=40, color = "black"))+
  theme(legend.title = element_blank(),
        legend.justification = "left",
        legend.text.align = 0)

out = "E:\\NLT\\DC_Tree\\08082024\\Diff_2m.jpg"
ggsave(out,height=10, width=40, dpi=600)
