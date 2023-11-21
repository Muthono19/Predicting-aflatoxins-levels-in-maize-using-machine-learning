#### Script to train a random forest model with cross-validation
### to predict Aflatoxin prevalence in maize in South and Eastern Africa
### Two approaches used "d.grid" = aggregated data 10km grids; d.farm= original XY
### .temp = withtempral matching data
### .gap = CHELSA tempral data gap filled with TERRACLIMATE for 2021/2022
#### To reduce the model uncertainity, we fit the model using log(aflatoxin+1)
### Created by Gachoki S
### Date 24th october 2023

######^^^^^^^^^^^^**********************************************################
###Load all the required libraries
library(tidyr)
library(ggplot2)
library(party)
library(caret)
library(ranger)
library(randomForest)
library(VSURF)
library(rgdal)
library(raster)
library(CAST)
library(reshape2)
library(viridis)
library(car)
library(corrplot)
library(gbm)
library(gridExtra)
library(terra)
library(pdp)
library(tidyverse)
library(tidyquant)
library(ggdist)
library(ggthemes)

#### set working directory #####
setwd("~/South Eastern Africa Aflatoxins")

#### Load the two dataset; at farm level and lowest admin level
#### Original XY
d.farm <- read.csv("dBaseSEA_afla_FinalR.csv")
d.farm$aflalog <- log(d.farm$Aflatoxin + 1)
names(d.farm)
#### aggregated at 10km grids
d.grid <- read.csv("dBaseSEA_afla_Final_Grid10km.csv")
d.grid$aflalog <- log(d.grid$Aflatoxin + 1)
names(d.grid)

### Make the relevant columns categorical
d.farm$cluster <- as.factor(d.farm$cluster)
d.farm$Time <- as.factor(d.farm$Time)
d.farm$Treat <- as.factor(d.farm$Treat)
d.grid$cluster <- as.factor(d.grid$cluster)
d.grid$Time <- as.factor(d.grid$Time)
d.grid$Treat <- as.factor(d.grid$Treat)

#### Exploratory data analysis: Rain clouds plots of the farm and admin dataset####
rn.farm <- ggplot(d.farm, aes(x = factor(cluster), y = aflalog, fill = factor(cluster))) +
  stat_halfeye(adjust = 0.5,justification = -0.2, .width = 0,point_colour = NA) +
  geom_boxplot(width = 0.2,outlier.color = NA,alpha = NA) +
  stat_dots(side = "left", justification = 1.1,binwidth = 0.03,dotsize=1, overlaps="nudge") +
  scale_fill_tq() +theme_tq() +
  labs(title = "Farm level",x = "",y = "log(aflatoxins+1)",fill = "Clusters") +
  theme(legend.position="none", text = element_text(size = 16, color="black"), 
        axis.title = element_text(size = 16, color = "black"),plot.title = element_text(size = 20, hjust = 0.5, face="bold"))

rn.admin <- ggplot(d.grid, aes(x = factor(cluster), y = aflalog, fill = factor(cluster))) +
  stat_halfeye(adjust = 0.5,justification = -0.2, .width = 0,point_colour = NA) +
  geom_boxplot(width = 0.2,outlier.color = NA,alpha = NA) +
  stat_dots(side = "left", justification = 1.1,binwidth = 0.03,dotsize=1, overlaps="nudge") +
  scale_fill_tq() +theme_tq() +
  labs(title = "10km grid level",x = "Long-term precipitation clusters",y = "log(aflatoxins+1)",fill = "Clusters") +
  theme(legend.position="none", text = element_text(size = 16, color="black"), 
        axis.title = element_text(size = 16, color = "black"),plot.title = element_text(size = 20, hjust = 0.5, face="bold"))

png(file = "ClusterRaincloudplots.png", width = 7400, height = 7400, units = "px", res = 800, type = "cairo")
grid.arrange(rn.farm,rn.admin,nrow=2)
dev.off()

### Plot correlation matrix for the continous predictor variables #####
cor.farm <- cor(d.farm[,c(12:33)])
cor.grid <- cor(d.grid[,c(7:28)])

png(file = "CorrelationMatrix.png", width = 14000, height =8500, units = "px", res = 650, type = "cairo")
par(mfrow=c(1,2))
corrplot(cor.farm, type="lower", tl.cex=2.5,cl.cex = 2.5,tl.col="black")
mtext("Farm level", side = 3, line = -1, font = 2, cex = 2.5)
corrplot(cor.grid, type="lower", tl.cex=2.5,cl.cex = 2.5,tl.col="black")
mtext("10-km grid level", side = 3, line = -1, font = 2, cex = 2.5)
dev.off()

###### *******************Model training initialization******************* ####

### Subset the various set of datasets: the Y and Predictors
#### %%%%% FARM LEVEL %%%%%********
###temporally matching satellite data
d.farm.temp <- d.farm[,c(34,11:19,24:33)]
names(d.farm.temp)
pred.farm.temp <- d.farm[,c(11:19,24:33)]
names(pred.farm.temp)
###gap-filled RS data
d.farm.gap <- d.farm[,c(34,11,14,16:18,20:33)]
names(d.farm.gap)
pred.farm.gap <- d.farm[,c(11,14,16:18,20:33)]
names(pred.farm.gap)

####### %%%%%%%%% GRID LEVEL #********
###temporally matching satellite adat
d.grid.temp <- d.grid[,c(29,6:14,19:28)]
names(d.grid.temp)
pred.grid.temp <- d.grid[,c(6:14,19:28)]
names(pred.grid.temp)
###gap-filled RS data
d.grid.gap <- d.grid[,c(29,6,9,11:13,15:28)]
names(d.grid.gap)
pred.grid.gap <- d.grid[,c(6,9,11:13,15:28)]
names(pred.grid.gap)

#### Feature elimination
##### calculate the best mtry
bestmtry.farm.temp <- tuneRF(pred.farm.temp,d.farm.temp$aflalog,stepFactor=1.5, improve=1e-5,ntree=300)
bestmtry.farm.gap <- tuneRF(pred.farm.gap,d.farm.gap$aflalog,stepFactor=1.5, improve=1e-5,ntree=300)
bestmtry.grid.temp <- tuneRF(pred.grid.temp,d.grid.temp$aflalog,stepFactor=1.5, improve=1e-5,ntree=300)
bestmtry.grid.gap <- tuneRF(pred.grid.gap,d.grid.gap$aflalog,stepFactor=1.5, improve=1e-5,ntree=300)
##### VSURF feature elimination
varvsurf.farm.temp <- VSURF(pred.farm.temp,d.farm.temp$aflalog,data=d.farm.temp,mtry=2,ntree=300,ntree.thres = 30,nfor.thres = 10,
                  ntree.interp = 10,nfor.interp = 10)
varvsurf.farm.temp$varselect.thres
varvsurf.farm.gap <- VSURF(pred.farm.gap,d.farm.gap$aflalog,data=d.farm.gap,mtry=3,ntree=300,ntree.thres = 30,nfor.thres = 10,
                            ntree.interp = 10,nfor.interp = 10)
varvsurf.farm.gap$varselect.thres

varvsurf.grid.temp <- VSURF(pred.grid.temp,d.grid.temp$aflalog,data=d.grid.temp,mtry=6,ntree=300,ntree.thres = 30,nfor.thres = 10,
                            ntree.interp = 10,nfor.interp = 10)
varvsurf.grid.temp$varselect.thres
varvsurf.grid.gap <- VSURF(pred.grid.gap,d.grid.gap$aflalog,data=d.grid.gap,mtry=4,ntree=300,ntree.thres = 30,nfor.thres = 10,
                           ntree.interp = 10,nfor.interp = 10)
varvsurf.grid.gap$varselect.thres

##thres.var <- data.frame(pred.grid.gap[,c( 18,  6,  2,  5, 19,  7, 17,  8,  9, 16, 15,  3, 10,  4,  1, 13, 14, 12, 11)])
##names(thres.var)

##### Define the formulas
fm.farm.temp <- aflalog ~ Longitudes+soiltemp+Latitudes+mintempTERRA+dem+pdsi+precTERRA+maxtempTERRA+vpd+orgcarb+
  evi+soilbulk+ndvi+sand+clay+CEC+PH+silt
fm.farm.gap <- aflalog ~ soiltemp+Longitudes+Latitudes+mintempCHELSA+maxtempCHELSA+pdsi+precCHELSA+dem+humidity+orgcarb+
  ndvi+evi+CEC+sand+soilbulk+PH+clay+silt

fm.grid.temp <- aflalog ~ soiltemp+Longitudes+dem+pdsi+Latitudes+precTERRA+mintempTERRA+vpd+maxtempTERRA+soilbulk+orgcarb+
  ndvi+silt+sand+evi+clay+PH+CEC
fm.grid.gap <- aflalog ~ Latitudes+maxtempCHELSA+pdsi+soiltemp+Longitudes+mintempCHELSA+dem+precCHELSA+humidity+orgcarb+
  soilbulk+ndvi+sand+evi+clay+CEC+silt+PH 

### create space-time folds for cross validation
set.seed(10)
fold.farm <- CreateSpacetimeFolds(d.farm,spacevar="cluster",timevar="Time", k=5)
fold.grid <- CreateSpacetimeFolds(d.grid,spacevar="cluster",timevar="Time", k=5)

## Tune grid and control parameters
ranger_grid.farm <- expand.grid(.mtry = 2,.splitrule = "variance",.min.node.size = 30)
ranger_grid.grid <- expand.grid(.mtry = 5,.splitrule = "variance",.min.node.size = 30)

trcntrl.farm <- trainControl(method="cv",savePredictions = "all",index=fold.farm$index,indexOut=fold.farm$indexOut)
trcntrl.grid <- trainControl(method="cv",savePredictions = "all",index=fold.grid$index,indexOut=fold.grid$indexOut)

### fit the model using caret Ranger method
#######%%%%%%%^^^^^^^^^^FARM LEVEL &&&^^^^^^^^^
model.farm.temp <- caret::train(fm.farm.temp,d.farm.temp,method="ranger",metric="RMSE",importance="permutation", 
                                tuneGrid = ranger_grid.farm,trControl = trcntrl.farm,num.trees=300)
model.farm.temp$finalModel
exp(sqrt(model.farm.temp$finalModel$prediction.error)) - 1

model.farm.gap <- caret::train(fm.farm.gap,d.farm.gap,method="ranger",metric="RMSE",importance="permutation", 
                               tuneGrid = ranger_grid.farm,trControl = trcntrl.farm,num.trees=300)
model.farm.gap$finalModel
exp(sqrt(model.farm.temp$finalModel$prediction.error)) - 1

#######%%%%%%%^^^^^^^^^^10KM GRID LEVEL &&&^^^^^^^^^
model.grid.temp <- caret::train(fm.grid.temp,d.grid.temp,method="ranger",metric="RMSE",importance="permutation", 
                                tuneGrid = ranger_grid.grid,trControl = trcntrl.grid,num.trees=300)
model.grid.temp$finalModel
exp(sqrt(model.grid.temp$finalModel$prediction.error)) - 1
mean(d.grid$Aflatoxin)

model.grid.gap <- caret::train(fm.grid.gap,d.grid.gap,method="ranger",metric="RMSE",importance="permutation", 
                               tuneGrid = ranger_grid.grid,trControl = trcntrl.grid,num.trees=300)
model.grid.gap$finalModel
exp(sqrt(model.grid.temp$finalModel$prediction.error)) - 1

### Plot predicted vs observed
### add the various predictions to the main database
d.farm$temppred <- model.farm.temp$finalModel$predictions
d.farm$gappred <- model.farm.gap$finalModel$predictions

d.grid$temppred <- model.grid.temp$finalModel$predictions
d.grid$gappred <- model.grid.gap$finalModel$predictions

obs.farm.temp <- ggplot(d.farm, aes(x = temppred, y = aflalog, color = Treat)) +geom_point(size = 2, pch=8) +
  labs(x = "Predicted log(aflatoxins+1)", y = "Observed log(aflatoxins+1)") +
  ggtitle("D1-farm") +scale_color_manual(values = c("1" = "magenta", "2" = "brown", "3" = "green"),
    labels = c("1" = "Unknown", "2" = "Treated", "3" = "Control")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth=2)+
  annotate("text", x = min(d.farm$temppred), y = max(d.farm$aflalog), 
           label = paste("r =", 0.6), hjust = 0, vjust = 1, size=8, col="red") +
  labs(color = "Treatment") +theme_base()+ theme(legend.position ="none",plot.title = element_text(hjust = 0.5, vjust = -7))

obs.farm.gap <- ggplot(d.farm, aes(x = gappred, y = aflalog, color = Treat)) +geom_point(size = 2, pch=8) +
  labs(x = "Predicted log(aflatoxins+1)", y = "Observed log(aflatoxins+1)") +
  ggtitle("D2-farm") +scale_color_manual(values = c("1" = "magenta", "2" = "brown", "3" = "green"),
                                  labels = c("1" = "Unknown", "2" = "Treated", "3" = "Control")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth=2)+
  annotate("text", x = min(d.farm$temppred), y = max(d.farm$aflalog), 
           label = paste("r =", 0.6), hjust = 0, vjust = 1, size=8, col="red") +
  labs(color = "Treatment") +theme_base()+ theme(legend.position ="none",plot.title = element_text(hjust = 0.5, vjust = -7))

obs.grid.temp <- ggplot(d.grid, aes(x = temppred, y = aflalog, color = Treat)) +geom_point(size = 2, pch=8) +
  labs(x = "Predicted log(aflatoxins+1)", y = "Observed log(aflatoxins+1)") +
  ggtitle("D1-10km") +scale_color_manual(values = c("1" = "magenta", "2" = "brown", "3" = "green"),
                                  labels = c("1" = "Unknown", "2" = "Treated", "3" = "Control")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth=2)+
  annotate("text", x = min(d.farm$temppred), y = max(d.farm$aflalog), 
           label = paste("r =", 0.54), hjust = 0, vjust = 1, size=8, col="red") +
  labs(color = "Treatment") +theme_base()+ theme(legend.position ="none",plot.title = element_text(hjust = 0.5, vjust = -7))

obs.grid.gap <- ggplot(d.grid, aes(x = gappred, y = aflalog, color = Treat)) +geom_point(size = 2, pch=8) +
  labs(x = "Predicted log(aflatoxins+1)", y = "Observed log(aflatoxins+1)") +
  ggtitle("D2-10km") +scale_color_manual(values = c("1" = "magenta", "2" = "brown", "3" = "green"),
                                  labels = c("1" = "Unknown", "2" = "Treated", "3" = "Control")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth=2)+
  annotate("text", x = min(d.farm$temppred), y = max(d.farm$aflalog), 
           label = paste("r =", 0.53), hjust = 0, vjust = 1, size=8, col="red") +
  labs(color = "Treatment") +theme_base()+ theme(legend.position ="none",plot.title = element_text(hjust = 0.5, vjust = -7))

png(file = "predvsobs_all.png", width = 7500, height = 6500, units = "px", res = 600, type = "cairo")
grid.arrange(obs.farm.temp,obs.farm.gap,obs.grid.temp,obs.grid.gap, ncol=2, nrow=2)
dev.off()

#### Variable importance and PDPs
### with temporal matching variables
imp.farm.temp <- varImp(model.farm.temp)
imp.farm.temp2 <- as.data.frame(imp.farm.temp$importance)
imp.farm.temp2$varnames <- rownames(imp.farm.temp2)
impplt.farm.temp <- ggplot(imp.farm.temp2, aes(x=reorder(varnames, Overall), y=Overall)) +  geom_point(color="blue",size=4)+
  ggtitle("D1-farm")+ xlab("") + ylab("")+ coord_flip()+theme_tq() + theme(plot.title = element_text(size=16, hjust = 0.5, color="black"),
                                                                    text = element_text(size = 16, face="bold",color = "black"))
imp.farm.gap <- varImp(model.farm.gap)
imp.farm.gap2 <- as.data.frame(imp.farm.gap$importance)
imp.farm.gap2$varnames <- rownames(imp.farm.gap2)
impplt.farm.gap <-ggplot(imp.farm.gap2, aes(x=reorder(varnames, Overall), y=Overall)) +  geom_point(color="blue",size=4)+
  ggtitle("D2-farm")+ xlab("") + ylab("")+ coord_flip()+theme_tq() + theme(plot.title = element_text(size=16, hjust = 0.5, color="black"),
                                                                      text = element_text(size = 16, face="bold",color = "black"))

imp.grid.temp <- varImp(model.grid.temp)
imp.grid.temp2 <- as.data.frame(imp.grid.temp$importance)
imp.grid.temp2$varnames <- rownames(imp.grid.temp2)
impplt.grid.temp <- ggplot(imp.grid.temp2, aes(x=reorder(varnames, Overall), y=Overall)) +  geom_point(color="blue",size=4)+
  ggtitle("D1-10km")+ xlab("") + ylab("")+ coord_flip()+theme_tq() + theme(plot.title = element_text(size=16, hjust = 0.5, color="black"),
                                                                           text = element_text(size = 16, face="bold",color = "black"))
imp.grid.gap <- varImp(model.grid.gap)
imp.grid.gap2 <- as.data.frame(imp.grid.gap$importance)
imp.grid.gap2$varnames <- rownames(imp.grid.gap2)
impplt.grid.gap <-ggplot(imp.grid.gap2, aes(x=reorder(varnames, Overall), y=Overall)) +  geom_point(color="blue",size=4)+
  ggtitle("D2-10km")+ xlab("") + ylab("")+ coord_flip()+theme_tq() + theme(plot.title = element_text(size=16, hjust = 0.5, color="black"),
                                                                      text = element_text(size = 16, face="bold",color = "black"))

png(file = "VarImp_all.png", width = 14000, height = 7500, units = "px", res = 750, type = "cairo")
grid.arrange(impplt.farm.temp,impplt.farm.gap,impplt.grid.temp,impplt.grid.temp,ncol=4)
dev.off()

#####Partial plots
library(pdp)
all.farm.temp = topPredictors(model.farm.temp,n=6)
pd.farm.temp <- NULL
for (i in all.farm.temp) {
  tmp <- partial(model.farm.temp, pred.var = i, data = d.farm.temp)
  names(tmp) <- c("x", "y")
  pd.farm.temp <- rbind(pd.farm.temp, cbind(tmp, predictor = i))
}
pd.farm.temp$predictor <- factor(pd.farm.temp$predictor, levels = unique(pd.farm.temp$predictor))
pp.farm.temp <- ggplot(pd.farm.temp, aes(x, y)) + geom_line(linewidth=1) + theme_classic() +theme(text = element_text(size = 25, face="bold",
          color = "black"),axis.text.y = element_text(size=18, face="bold",color = "black"),
          plot.title = element_text(hjust = 0.5,color="red"),axis.text = element_text(size = 18,color = "black"))+
  ggtitle("D1 farm")+ ylab("log(aflatoxins+1)") +xlab("") + facet_wrap(~ predictor, scales = "free")

### with gap filled variables
all.farm.gap = topPredictors(model.farm.gap,n=6)
pd.farm.gap <- NULL
for (i in all.farm.gap) {
  tmp <- partial(model.farm.gap, pred.var = i)
  names(tmp) <- c("x", "y")
  pd.farm.gap <- rbind(pd.farm.gap, cbind(tmp, predictor = i))
}
pd.farm.gap$predictor <- factor(pd.farm.gap$predictor, levels = unique(pd.farm.gap$predictor))
pp.farm.gap <- ggplot(pd.farm.gap, aes(x, y)) + geom_line(linewidth=1) + theme_classic() +
  theme(text = element_text(size = 25, face="bold",color = "black"),axis.text.y = element_text(size=18, face="bold",color = "black"),
  plot.title = element_text(hjust = 0.5,color="red"),axis.text = element_text(size = 18,color = "black"))+
  ggtitle("D2 farm")+ ylab("log(aflatoxins+1)") +xlab("") + facet_wrap(~ predictor, scales = "free")

all.grid.temp = topPredictors(model.grid.temp,n=6)
pd.grid.temp <- NULL
for (i in all.grid.temp) {
  tmp <- partial(model.grid.temp, pred.var = i, data = d.grid.temp)
  names(tmp) <- c("x", "y")
  pd.grid.temp <- rbind(pd.grid.temp, cbind(tmp, predictor = i))
}
pd.grid.temp$predictor <- factor(pd.grid.temp$predictor, levels = unique(pd.grid.temp$predictor))
pp.grid.temp <- ggplot(pd.grid.temp, aes(x, y)) + geom_line(linewidth=1) +
  theme_classic() +theme(text = element_text(size = 25, face="bold",
  color = "black"),axis.text.y = element_text(size=18, face="bold",color = "black"),
plot.title = element_text(hjust = 0.5,color="red"),axis.text = element_text(size = 18,color = "black"))+
  ggtitle("D1 10km")+ ylab("log(aflatoxins+1)") +xlab("") + facet_wrap(~ predictor, scales = "free")

### with gap filled variables
all.grid.gap = topPredictors(model.grid.gap,n=6)
pd.grid.gap <- NULL
for (i in all.grid.gap) {
  tmp <- partial(model.grid.gap, pred.var = i)
  names(tmp) <- c("x", "y")
  pd.grid.gap <- rbind(pd.grid.gap, cbind(tmp, predictor = i))
}
pd.grid.gap$predictor <- factor(pd.grid.gap$predictor, levels = unique(pd.grid.gap$predictor))
pp.grid.gap <- ggplot(pd.grid.gap, aes(x, y)) + geom_line(linewidth=1) + theme_classic() +
  theme(text = element_text(size = 25, face="bold",
  color = "black"),axis.text.y = element_text(size=18, face="bold",color = "black"),
  plot.title = element_text(hjust = 0.5,color="red"),axis.text = element_text(size = 18,color = "black"))+
  ggtitle("D2 10km")+ ylab("log(aflatoxins+1)") +xlab("") + facet_wrap(~ predictor, scales = "free")

png(file = "partialplots_all.png", width = 14000, height =10000, units = "px", res = 600, type = "cairo")
grid.arrange(pp.farm.temp,pp.farm.gap,pp.grid.temp,pp.grid.gap,nrow=2, ncol=2)
dev.off()

#### Spatial predictions
### with Wet year 2020
farm.temp.ras20 <- stack("Longitudes.tif","stemp2020.tif","Latitudes.tif","mintemp2020.tif","dem.tif","pdsi2020.tif",
"prec2020.tif","maxtemp2020.tif","vpd2020.tif","orgcarb.tif","evi2020.tif","soilbulk.tif","ndvi2020.tif",
"sand.tif","clay.tif","CEC.tif","PH.tif","silt.tif")
names(farm.temp.ras20)<- c("Longitudes","soiltemp","Latitudes","mintempTERRA","dem","pdsi","precTERRA","maxtempTERRA","vpd",
"orgcarb","evi","soilbulk","ndvi","sand","clay","CEC","PH","silt")

farm.temp.rast20 <- rast(farm.temp.ras20)
gc()
pred.farm.temp20 <- predict(farm.temp.rast20,model.farm.temp,na.rm=T)
#pred.farm.temp2 <- exp(pred.farm.temp) - 1
plot(pred.farm.temp20)
#aoa.farm.temp <- aoa(farm.temp.rast, model.farm.temp, useWeight = TRUE)
#plot(aoa.farm.temp$AOA)#, col=c("grey","transparent"), add=T)
writeRaster(pred.farm.temp20,"prediction_Farm_Temp_Wet.tif")

### with dry year 2021
farm.temp.ras21 <- stack("Longitudes.tif","stemp2021.tif","Latitudes.tif","mintemp2021.tif","dem.tif","pdsi2021.tif",
                         "prec2021.tif","maxtemp2021.tif","vpd2021.tif","orgcarb.tif","evi2021.tif","soilbulk.tif","ndvi2020.tif",
                         "sand.tif","clay.tif","CEC.tif","PH.tif","silt.tif")
names(farm.temp.ras21)<- c("Longitudes","soiltemp","Latitudes","mintempTERRA","dem","pdsi","precTERRA","maxtempTERRA","vpd",
                           "orgcarb","evi","soilbulk","ndvi","sand","clay","CEC","PH","silt")

farm.temp.rast21 <- rast(farm.temp.ras21)
gc()
pred.farm.temp21 <- predict(farm.temp.rast21,model.farm.temp,na.rm=T)
#pred.farm.temp2 <- exp(pred.farm.temp) - 1
plot(pred.farm.temp21)
#aoa.farm.temp <- aoa(farm.temp.rast, model.farm.temp, useWeight = TRUE)
#plot(aoa.farm.temp$AOA)#, col=c("grey","transparent"), add=T)
writeRaster(pred.farm.temp21,"prediction_Farm_Temp_Dry.tif",overwrite=TRUE)

###wet year future projections
farm.temp.ras20F <- stack("Longitudes.tif","stemp2020.tif","Latitudes.tif","mintemp70.tif","dem.tif","pdsi2020.tif",
                         "prec70.tif","maxtemp70.tif","vpd2020.tif","orgcarb.tif","evi2020.tif","soilbulk.tif","ndvi2020.tif",
                         "sand.tif","clay.tif","CEC.tif","PH.tif","silt.tif")
names(farm.temp.ras20F)<- c("Longitudes","soiltemp","Latitudes","mintempTERRA","dem","pdsi","precTERRA","maxtempTERRA","vpd",
                           "orgcarb","evi","soilbulk","ndvi","sand","clay","CEC","PH","silt")

farm.temp.rast20F <- rast(farm.temp.ras20F)
gc()
pred.farm.temp20F <- predict(farm.temp.rast20F,model.farm.temp,na.rm=T)
#pred.farm.temp2 <- exp(pred.farm.temp) - 1
plot(pred.farm.temp20F)
#aoa.farm.temp <- aoa(farm.temp.rast, model.farm.temp, useWeight = TRUE)
#plot(aoa.farm.temp$AOA)#, col=c("grey","transparent"), add=T)
writeRaster(pred.farm.temp20F,"prediction_Farm_Temp_WetF.tif")

### with dry year 2021 future projections
farm.temp.ras21F <- stack("Longitudes.tif","stemp2021.tif","Latitudes.tif","mintemp70.tif","dem.tif","pdsi2021.tif",
                         "prec70.tif","maxtemp70.tif","vpd2021.tif","orgcarb.tif","evi2021.tif","soilbulk.tif","ndvi2020.tif",
                         "sand.tif","clay.tif","CEC.tif","PH.tif","silt.tif")
names(farm.temp.ras21F)<- c("Longitudes","soiltemp","Latitudes","mintempTERRA","dem","pdsi","precTERRA","maxtempTERRA","vpd",
                           "orgcarb","evi","soilbulk","ndvi","sand","clay","CEC","PH","silt")

farm.temp.rast21F <- rast(farm.temp.ras21F)
gc()
pred.farm.temp21F <- predict(farm.temp.rast21F,model.farm.temp,na.rm=T)
#pred.farm.temp2 <- exp(pred.farm.temp) - 1
plot(pred.farm.temp21F)
#aoa.farm.temp <- aoa(farm.temp.rast, model.farm.temp, useWeight = TRUE)
#plot(aoa.farm.temp$AOA)#, col=c("grey","transparent"), add=T)
writeRaster(pred.farm.temp21F,"prediction_Farm_Temp_DryF.tif",overwrite=TRUE)

#### lEVELPLOTS FOR TH PREDICTIONS
### Load the national boundaries
aoi <- readOGR("AOI_BoundariesPRJ.shp")

## Load the predicted rasters
afla.pred <- stack("prediction_Farm_Temp_Dry.tif","prediction_Farm_Temp_DryF.tif",
                   "prediction_Farm_Temp_Wet.tif","prediction_Farm_Temp_WetF.tif")
afla.pred2 <- exp(afla.pred)-1

plot(afla.pred2)
names(afla.pred) <- c("","","","")

allnames <- c("Driest 2021","Driest 2070","Wettest 2020","Wettest 2070")
library(RColorBrewer)

myInverseColors <- rev(brewer.pal(11, "RdYlGn"))

myTheme=rasterTheme(region=myInverseColors)

png(file = "predictionsAll.png", width = 4500, height = 4500, units = "px", res = 600, type = "cairo")
levelplot(afla.pred2,par.settings=myTheme, names.attr=allnames, layout=c(2,2), margin=FALSE,ylab="",
          panel=function(...) {
            panel.levelplot(...)
            sp.polygons(aoi, fill = NA, col = "black", lwd=0.3)
          })
dev.off()






### with gap filled variables
farm.gap.ras <- stack("maxtempCHELSA.tif","long.tif","pdsi.tif","lat.tif","humidity.tif","mintempCHELSA.tif","soiltemp.tif",
"dem.tif","orgcarb.tif","evi.tif","soilbulk.tif","precCHELSA.tif","clay.tif","sand.tif","ndvi.tif","PH.tif","silt.tif","CEC.tif")
farm.gap.rast <- rast(farm.gap.ras)
gc()
pred.farm.gap <- predict(farm.gap.rast,model.farm.gap,na.rm=T)
pred.farm.gap.crp <- mask(pred.farm.gap,cropmask)
plot(pred.farm.gap.crp)
writeRaster(pred.farm.gap.crp,"prediction_Farm_Gap.tif")

### with admin variables
admin.temp.ras <- stack("dem_a.tif","evi_a.tif","ndvi_a.tif","mintemp_a.tif","soiltemp_a.tif","sand_a.tif","pdsi_a.tif","lat.tif",
"maxtemp_a.tif","PH_a.tif","soilbulk_a.tif","long.tif","silt_a.tif","prec_a.tif","orgcarb_a.tif","vpd_a.tif","clay_a.tif","CEC_a.tif")
admin.temp.rast <- rast(admin.temp.ras)
gc()
pred.admin.temp <- predict(admin.temp.rast,model.admin.temp,na.rm=T)
pred.admin.temp.crp <- mask(pred.admin.temp,cropmask)
plot(pred.admin.temp.crp)
writeRaster(pred.admin.temp.crp,"prediction_admin_Temp.tif")

#### Extrapolating the temporal matching model to future climatic projections
fut.temp.ras <- stack("mintemp70.tif","long.tif","pdsi.tif","lat.tif","dem.tif","soiltemp.tif","maxtemp70.tif",
                       "vpd.tif","sand.tif","soilbulk.tif","orgcarb.tif","prec70.tif","evi.tif","clay.tif","ndvi.tif","PH.tif","silt.tif","CEC.tif")
names(fut.temp.ras) <- c("mintempTERRA","long","pdsi","lat","dem","soiltemp","maxtempTERRA","vpd","sand","soilbulk","orgcarb","precTERRA","evi","clay","ndvi","PH","silt","CEC")
fut.temp.rast <- rast(fut.temp.ras)
gc()
pred.fut.temp <- predict(fut.temp.rast,model.farm.temp,na.rm=T)
pred.fut.temp.crp <- mask(pred.fut.temp,cropmask)
plot(pred.fut.temp.crp)
writeRaster(pred.fut.temp.crp,"prediction_Future_Temp.tif")


##### *********************************The END ************************************ #####


