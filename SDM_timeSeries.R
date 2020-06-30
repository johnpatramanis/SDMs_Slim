library(data.table)
library(raster)
library(randomForest)
library(lattice)
library(RColorBrewer)
library(PresenceAbsence)

# this tutorial was lifted from:
# tutorial in: https://damariszurell.github.io/SDM-Intro/
# the tutorial dowloads large file, so better load the outputs:
#"bio_analog_df.txt"
#"bio_novel_df_df.txt"

# at the end of the tutorial I added how to approximate an interpolation between a current and future SDMs for a time series of 50 years
# starting from these SDM outputs

avi_dat <- read.table('Data_SwissBreedingBirds.csv', header=T, sep=',')
summary(avi_dat)
names(avi_dat)
avi_cols <- c('Turdus_torquatus', 'bio_5', 'bio_2', 'bio_14', 'std', 'rad', 'blockCV_tile')
avi_df <- data.frame(avi_dat)[avi_cols]
head(avi_df)
# Please note that you have to set download=T if you haven't downloaded the data before:
bio_curr <- getData('worldclim', var='bio', res=0.5, lon=5.5, lat=45.5,path = "./",download = F)[[c(2,5,14)]]
# Please note that you have to set download=T if you haven't downloaded the data before:
bio_fut <- getData('CMIP5', var='bio', res=0.5, lon=5.5, lat=45.5, rcp=45, model='NO', year=50, path = "./", download=F)[[c(2,5,14)]]

# A spatial mask of Switzerland in Swiss coordinates
bg <- raster('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')

# The spatial extent of Switzerland in Lon/Lat coordinates is roughly:
ch_ext <- c(5, 11, 45, 48)

# Crop the climate layers to the extent of Switzerland
bio_curr <- crop(bio_curr, ch_ext)

# Re-project to Swiss coordinate system and clip to Swiss political bounday
bio_curr <- projectRaster(bio_curr, bg)
bio_curr <- resample(bio_curr, bg)
bio_curr <- mask(bio_curr, bg)
names(bio_curr) <- c('bio_2', 'bio_5', 'bio_14')

# For storage reasons the temperature values in worldclim are multiplied by 10. For easier interpretability, we change it back to °C.
bio_curr[[1]] <- bio_curr[[1]]/10
bio_curr[[2]] <- bio_curr[[2]]/10

# Repeat above steps for future climate layers
bio_fut <- crop(bio_fut, ch_ext)
bio_fut <- projectRaster(bio_fut, bg)
bio_fut <- resample(bio_fut, bg)
bio_fut <- mask(bio_fut, bg)
names(bio_fut) <- c('bio_2', 'bio_5', 'bio_14')
bio_fut[[1]] <- bio_fut[[1]]/10
bio_fut[[2]] <- bio_fut[[2]]/10

plot(bio_curr)

# I recommend to always store the scaling coefficients:
scale_attrib <- attributes(scale(avi_df[,2:6]))[3:4]

# Make new dataset with scaled predictors
avi_dfst <- avi_df
avi_dfst[,2:6] <- scale(avi_df[,2:6])
values(bio_curr) <- scale(values(bio_curr),center = scale_attrib$`scaled:center`[names(bio_curr)], scale = scale_attrib$`scaled:scale`[names(bio_curr)])
values(bio_fut) <- scale(values(bio_fut),center = scale_attrib$`scaled:center`[names(bio_fut)], scale = scale_attrib$`scaled:scale`[names(bio_fut)])

# Fit GLM
m_glm <- glm( Turdus_torquatus ~ bio_2 + I(bio_2^2) + bio_5 + I(bio_5^2) + bio_14 + I(bio_14^2), family='binomial', data=avi_dfst)
summary(m_glm)
pred <- c('bio_2', 'bio_5', 'bio_14')
par(mfrow=c(1,3)) 
# Loop over all predictors:
for (i in 1:3) {
  # Make new dummy data set with all predictors kept constant at their mean
  xz <- data.frame(sapply(colMeans(avi_dfst[,pred]),rep,each=50))
  
  # Let predictor i be a equal-spaced sequence from lowest to highest observed value
  xz[,pred[i]] <- seq(min(avi_dfst[,pred[i]]),max(avi_dfst[,pred[i]]),length=50)
  
  # Make predictions to dummy data set
  xz$z <- predict(m_glm, newdata=xz, type='response')
  
  # Plot response along predictor i
  plot(xz[,i],xz$z,type='l', xlab=pred[i], ylab='Occurrence probability')
}
# We prepare the response surface by making a dummy data set where two predictor variables range from their minimum to maximum value, and the remaining predictor is kept constant at its mean:
xyz <- data.frame(expand.grid(seq(min(avi_dfst[,pred[1]]),max(avi_dfst[,pred[1]]),length=50), seq(min(avi_dfst[,pred[2]]),max(avi_dfst[,pred[2]]),length=50)), mean(avi_dfst[,pred[3]]))
names(xyz) <- pred
# Make predictions
xyz$z <- predict(m_glm, xyz, type='response')
summary(xyz)
# Make a colour scale
cls <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(100)
# plot 3D-surface
wireframe(z ~ bio_2 + bio_5, data = xyz, zlab = list("Occurrence prob.", rot=90),
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), zlim = c(0, 1), 
          main='GLM', xlab='bio_2', ylab='bio_5', screen=list(z = 120, x = -70, y = 3))

#Codes for plotting inflated response curves are available from github and we need to source these before plotting. You can best imagine inflated response curves as multiple slices through above 3D-surface.
# Get source code for inflated response curves:
script<-readLines("https://raw.githubusercontent.com/damariszurell/Rcodes_MapNovelEnvironments_SDMs/master/appendixS1_functions.r")
eval(parse(text  = script)) 
# Plot inflated response curves:
par(mfrow=c(1,3)) 
inflated.response(m_glm, predictors = avi_dfst[,pred], method = "stat6", lwd = 3, main='GLM') 

# Fit RF
(m_rf <- randomForest( x=avi_dfst[,2:4], y=avi_dfst[,1], ntree=1000, nodesize=10, importance =T))
# Variable importance:
importance(m_rf,type=1)
varImpPlot(m_rf)
# Look at single trees:
head(getTree(m_rf,1,T))
# Now, we plot response curves in the same way as we did for GLMs above:
par(mfrow=c(1,3)) 
for (i in 1:3) {
  xz <- data.frame(sapply(colMeans(avi_dfst[,pred]),rep,each=50))
  xz[,pred[i]] <- seq(min(avi_dfst[,pred[i]]),max(avi_dfst[,pred[i]]),length=50)
  xz$z <- predict(m_rf, newdata=xz)
  plot(xz[,i],xz$z,type='l', xlab=pred[i], 
       ylab='Occurrence probability',ylim=c(0,1), main='Random forest')
}
# Plot the response surface:
xyz$z <- predict(m_rf, xyz)
wireframe(z ~ bio_2 + bio_5, data = xyz, zlab = list("Occurrence prob.", rot=90),
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), zlim = c(0, 1), 
          main='RF', xlab='bio_2', ylab='bio_5', screen=list(z = 120, x = -70, y = 3))
# Plot inflated response curves:
par(mfrow=c(1,3)) 
inflated.response(m_rf, predictors = avi_dfst[,pred], method = "stat6", lwd = 3, main='RF') 

#2.4.1 Ring Ouzel: Model assessment
#We assess model predictive performance using (spatial block) cross-validation. Also, we select thresholds for binarising the predicted occurrence probabilities based on the cross-validated predictions.
#This time, we will need a number of helper functions for making predictions, for making cross-validated predictions, and for assessing model performance.

# Function to make predictions:
make.preds <- function(model, newdata) {
  switch(class(model)[1],
         glm = predict(model, newdata, type='response'),
         randomForest = predict(model, newdata))
}

# Function to make cross-validated predictions of the folds ("cv.tiles") are known. Remember that our data set already contains these folds in the column "blockCV_tile":
crossval.preds <- function(model, dat, pred, cv.tiles) {
  
  # How many folds are the data split in?
  kfold <- length(unique(cv.tiles))
  
  # Prepare data frame to store the cross-validated predictions
  cross.val.preds <- data.frame(row = row.names(dat), 
                                cross.val.preds = numeric(length = nrow(dat))) 
  
  # Loop through all folds
  for(i in seq_len(kfold)){
    # All data rows except those belonging to fold i are used as training data:
    cv.train <- dat[cv.tiles!=i,]
    # All data rows belonging to fold i are used as test data:
    cv.test <- dat[cv.tiles ==i,]
    
    # We update the model for the new training data. We use the switch statement to check the model class. 
    modtmp <- switch(class(model)[1],
                     glm = update(model, data=cv.train),
                     randomForest = update(model, data=cv.train))
    
    # We make predictions for the test data (fold i):
    cross.val.preds[which(cv.tiles ==i),2] <- make.preds(modtmp, cv.test[, pred])
  }
  
  return(cross.val.preds[,2])
}

# Function to calculate model performance based on several threshold-depedent and threshold-independent measures:
calc.eval <- function(obs, predictions, thresh.method='MaxSens+Spec'){
  require(PresenceAbsence)
  
  # Helper functions:
  # True Skill Statistic:
  TSS = function(cmx){
    PresenceAbsence::sensitivity(cmx, st.dev=F) + 
      PresenceAbsence::specificity(cmx, st.dev=F) - 1
  }
  
  # Prepare data frame for optimising the threshold for binarising the predictions:
  thresh.dat <- data.frame(ID=length(obs), 
                           obs = obs,
                           pred = predictions)
  
  # Find the optimal threshold:
  thresh <- optimal.thresholds(DATA= thresh.dat)
  
  # Make a confusion matrix using the optimal threshold:
  cmx.maxSSS <- cmx(DATA= thresh.dat, threshold=thresh[thresh$Method==thresh.method,2])
  
  # Calculate and return AUC, TSS, Sensitivity, and Specificity
  data.frame(AUC = PresenceAbsence::auc(thresh.dat, st.dev=F),
             TSS = TSS(cmx.maxSSS), 
             Sens = PresenceAbsence::sensitivity(cmx.maxSSS, st.dev=F),
             Spec = PresenceAbsence::specificity(cmx.maxSSS, st.dev=F),
             thresh = thresh[thresh$Method==thresh.method,2])
}

# Make cross-validated predictions for GLM:
crosspred_glm <- crossval.preds(m_glm, dat= avi_dfst[!is.na(avi_dfst$blockCV_tile),], pred=pred, cv.tiles= avi_dfst[!is.na(avi_dfst$blockCV_tile),'blockCV_tile'])
# Make cross-validated predictions for RF:
crosspred_rf <- crossval.preds(m_rf, dat= avi_dfst[!is.na(avi_dfst$blockCV_tile),], pred=pred, cv.tiles= avi_dfst[!is.na(avi_dfst$blockCV_tile),'blockCV_tile'])
# Look at correlation between GLM and RF predictions:
dev.off()
plot(crosspred_glm, crosspred_rf, pch=19, col='grey35')
#Next, we assess cross-validated model performance. We inspect different measures: AUC, the area under the receiver operating characteristic (ROC) curve (Hosmer and Lemeshow 2013); TSS, the true skill statistic (Allouche, Tsoar, and Kadmon 2006); sensitivity, the true positive rate; and specificity, the true negative rate. Simultaneously, we estimate the optimal threshold for making binary predictions. For this, we use a threshold that maximises TSS or maximised the sum of sensitivity and specificity (Liu et al. 2005).
(eval_glm <- calc.eval(obs = avi_dfst[!is.na(avi_dfst$blockCV_tile),1], predictions = crosspred_glm))
(eval_rf <- calc.eval(obs = avi_dfst[!is.na(avi_dfst$blockCV_tile),1], predictions = crosspred_rf))
#We can also combine predictions from the two SDM algorithms and make an ensemble prediction, for example by taking the median.
# Derive median predictions:
crosspred_ens <- apply(data.frame(crosspred_glm, crosspred_rf),1,median)
# Evaluate ensemble predictions
(eval_ens <- calc.eval(obs = avi_dfst[!is.na(avi_dfst$blockCV_tile),1], predictions = crosspred_ens))

#2.5 Predictions
#Now that we carefully fitted the SDMs, inspected model and extrapolation behaviour, and assessed predictive performance, it is finally time to make predictions in space and time. Importance points to consider here are quantification of uncertainty due to input data, algorithms, model complexity and boundary conditions (e.g. climate scenarios)(Araújo et al. 2019; Thuiller et al. 2019). When transferring the model to a different geographic area or time period, it is also recommended to quantify uncertainty due to novel environments (Zurell, Elith, and Schroeder 2012).
# Make predictions to current climate:
bio_curr_df <- data.frame(rasterToPoints(bio_curr))
bio_curr_df$pred_glm <- make.preds(m_glm, bio_curr_df)
bio_curr_df$pred_rf <- make.preds(m_rf, bio_curr_df)
bio_curr_df$pred_ens <- apply(bio_curr_df[,-c(1:5)],1,median)
# Make binary predictions:
bio_curr_df$bin_glm <- ifelse(bio_curr_df$pred_glm > eval_glm$thresh, 1, 0)
bio_curr_df$bin_rf <- ifelse(bio_curr_df$pred_rf > eval_rf$thresh, 1, 0)
bio_curr_df$bin_ens <- ifelse(bio_curr_df$pred_ens > eval_ens$thresh, 1, 0)
# Make raster stack of predictions:
r_pred_curr <- rasterFromXYZ(bio_curr_df[,-c(3:5)])
plot(r_pred_curr)

# Assess novel environments in future climate layer:
bio_fut_df <- data.frame(rasterToPoints(bio_fut))
# Values of 1 in the eo.mask will indicate novel environmental conditions
bio_fut_df$eo.mask <- eo.mask(avi_dfst[,pred], bio_fut_df[,pred])
plot(rasterFromXYZ(bio_fut_df[,-c(3:5)]), main='Environmental novelty')
# Make predictions to future climate:
bio_fut_df$pred_glm <- make.preds(m_glm, bio_fut_df)
bio_fut_df$pred_rf <- make.preds(m_rf, bio_fut_df)
bio_fut_df$pred_ens <- apply(bio_fut_df[,-c(1:5)],1,median)
# Make binary predictions:
bio_fut_df$bin_glm <- ifelse(bio_fut_df$pred_glm > eval_glm$thresh, 1, 0)
bio_fut_df$bin_rf <- ifelse(bio_fut_df$pred_rf > eval_rf$thresh, 1, 0)
bio_fut_df$bin_ens <- ifelse(bio_fut_df$pred_ens > eval_ens$thresh, 1, 0)
# Make raster stack of predictions:
r_pred_fut <- rasterFromXYZ(bio_fut_df[,-c(3:5)])
plot(r_pred_fut[[-1]])

# Predictions to analogous climates:
bio_analog_df <- bio_fut_df[,c('x','y','pred_glm','pred_rf')]
bio_analog_df[bio_fut_df$eo.mask>0,c('pred_glm','pred_rf')] <- NA
plot(rasterFromXYZ(bio_analog_df))
# Predictions to novel climates:
bio_novel_df <- bio_fut_df[,c('x','y','pred_glm','pred_rf')]
bio_novel_df[bio_fut_df$eo.mask==0,c('pred_glm','pred_rf')] <- NA
plot(rasterFromXYZ(bio_novel_df))
write.table(bio_analog_df[,c(1:2,4)],"bio_analog_df.txt",row.names = F,sep="\t",quote = F)
write.table(bio_novel_df[,c(1:2,4)],"bio_novel_df_df.txt",row.names = F,sep="\t",quote = F)

######### 
# Load outputs and project to future
#########

SDM_current<-fread("bio_analog_df.txt")
SDM_current=rasterFromXYZ(SDM_current)
plot(SDM_current)
SDM_future<-fread("bio_novel_df_df.txt")
SDM_future=rasterFromXYZ(SDM_future)
plot(SDM_future)

# Interpolate values between current and future
# the future layers are for year 2050, we can get data further in time up to year 2100 (I think)
# current data is approximately for condition around 2000
# 50 year series
# raster of extent = predictions
r <- SDM_current
# all values NA
values(r) <- NA
# 50 rasters (here we can generate N rasters)
x <- sapply(1:50, function(...) r)
# assign 1st raster as current and last raster as future
values(x[[1]]) <- values(SDM_current)
# convert NA values to 0
values(x[[1]])[is.na(values(x[[1]]))]=0
# same for future projection which is stored in the LAST layer
values(x[[50]]) <- values(SDM_future)
values(x[[50]])[is.na(values(x[[50]]))]=0
#plot(x[[1]])
#plot(x[[50]])
# stack them
s <- stack(x)
#plot(s)
# interpolate NA's
# this is brute force and I'm not sure how this function operates, we might need a more sophisticated approach
z <- approxNA(s)
#plot(z)
dim(z)
zValues<-data.table(values(z))
write.table(zValues,"timeSeries_50years_nrow221_ncol349.txt",row.names = F,quote = F,sep = "\t")


