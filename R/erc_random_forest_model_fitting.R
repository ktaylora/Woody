cat("## Eastern Redcedar Landcover Builder ###\n");

suppressMessages(require(rgdal))
suppressMessages(require(randomForest))
suppressMessages(require(raster))
suppressMessages(require(pROC))

setwd("/home/ktaylora/Incoming/offyear_study/")

#
# RUNTIME ARGUMENTS
#

ARGV = commandArgs(trailingOnly=T)

SCALE = 10;
TIME_STEP = 5 # number of years between the NAIP imagery segments used for our trend models
ACCEPTABLE_TYPES = c("grass","forest","water","row_crop","small_grains","erc_high","erc_low")
FAIL_THRESHOLD = 0.4 # proportion of training records containing bunk values
N_BS_REPLICATES = 999
DO_LOCAL_EXTRACTION = FALSE

if( length(ARGV) > 0 ){
  if( file.exists(ARGV[1]) ) TRAINING_DATA = ARGV[1]
} else {
  TRAINING_DATA = paste("~/Incoming/offyear_study/vector/erc_training_data/erc_balanced_training_pts_2015_2017_",SCALE,"m.geojson",sep="")
}

if(!file.exists(TRAINING_DATA)){
  cat("DEBUG: Fatal -- Training data file path not found\n")
}

# add latitude and longitude to our stack
# EXPLANATORY_COVS <- c(EE_RASTER_STACK_NAMES,"lon","lat")

#
# Local Functions
#
est_runtime <- function(expression){
    elapsed = 3;
    return(as.vector(summary(system.time(
        expression
    )))[elapsed])
}

bs_rf_fit_df <- function(df=NULL, n=N_BS_REPLICATES, m_rf=NULL, target_rf_classes=NULL){
    cl <- parallel::makeCluster(parallel::detectCores()*0.75)
    parallel::clusterExport(cl, varlist=c('df','n','m_rf','target_rf_classes'), envir=environment())
    df <- do.call(rbind, parallel::parLapply(
        cl=cl,
        X=1:n,
        fun=function(n){
            library('randomForest')
            p <- predict(
                m_rf, 
                df[sample(1:nrow(df), replace=T),], 
                type="response"
            )
            p <- as.vector(table(as.vector(p)))[target_rf_classes] 
            return( round(p / nrow(df),4) )
    }));
    parallel::stopCluster(cl); rm(cl); 
    return(df)
}

#
# MAIN
#

if( DO_LOCAL_EXTRACTION ){
  EE_RASTER_STACK_NAMES <- c("BLAH!")
  cat("DEBUG: Reading explanatory variables as raster stack\n")
  explanatories <- raster::stack(paste("/home/ktaylora/erc_study_naip_ndvi_2015_2017_",SCALE,"m.tif",sep=""))
  explanatories <- raster::stack(
    explanatories,
    raster::init(explanatories, 'x'),
    raster::init(explanatories, 'y')
  )
  names(explanatories) <- c(EE_RASTER_STACK_NAMES, "lon", "lat")
  cat("DEBUG: Performing raster extraction and QA checks for non-sense values\n")
  training_data <- cbind(
    type=training_data$type, 
    raster::extract(x=explanatories, y=training_data, df=T)
  )
} else {
  cat("DEBUG: Reading model training data from disc\n")

  training_data <- suppressMessages(suppressWarnings(readOGR(
    TRAINING_DATA, 
    stringsAsFactors=F, 
    verbose=F,
    require_geomType='wkbPoint'
  )))

  cat("DEBUG: Number of training records:\n")
  print(table(training_data$type))
  
  cat("DEBUG: Cleaning column names and affixing latitude & longitude (if missing)\n")
  colnames(training_data@data) <- tolower(colnames(training_data@data))

  if(sum(grepl(colnames(training_data@data), pattern="lon|lat")) < 2){
    training_data@data <- cbind(
      training_data@data, 
      lon=training_data@coords[,1], 
      lat=training_data@coords[,2]
    )
  }

  cat("DEBUG: Dropping columns that don't look useful\n")

  training_data@data <- training_data@data[,!grepl(tolower(colnames(training_data@data)), pattern="nlcd|nass|fid|id|row_crop|small")]
  training_data <- as.data.frame(training_data)

  # figure out our greeness index
  if( sum(grepl(colnames(training_data), pattern="ndvi")) > 0 ){
    cat("DEBUG: Looks like we are using NDVI as our greenness index\n")
    GREENESS_INDEX <- "ndvi"
  } else if ( sum(grepl(colnames(training_data), pattern="g_")) > 0 ) {
    cat("DEBUG: Looks like we are using an alternative greeness index specified by the user (e.g., as 'green_') \n")
    GREENESS_INDEX <- "gf"
  } else {
    cat("DEBUG: Couldn't find a suitable greeness index in input data -- this shouldn't happen\n")
    stop()
  }

  # Drop any lurking non-sense fields in our training data
  EXPLANATORY_COVS <- colnames(training_data)
    EXPLANATORY_COVS <- EXPLANATORY_COVS[grepl(EXPLANATORY_COVS, pattern="type|ndvi|g_|^r$|^g$|^b$|^n$|elev|aspect|slope|lon|lat")]
  cat("DEBUG: Pool of explanatory/response variables available:", paste(EXPLANATORY_COVS, collapse=","), "\n")
  training_data <- training_data[,EXPLANATORY_COVS]
}

cat("DEBUG: Checking for NA inflation in our training data\n")
NA_EXTRACTION_ERR <- sum(is.na(training_data[,GREENESS_INDEX]))/length(training_data[,GREENESS_INDEX])
if(NA_EXTRACTION_ERR > 0){
  warning(
    round(NA_EXTRACTION_ERR,3),
    "% of records contained NA values; Dropping. Check to make sure your",
    "training data don't fall outside of the boundary of your project extent."
  )
  training_data <- training_data[!is.na(training_data[,GREENESS_INDEX]),]
}

# Maybe we have an acceptable number of training points that are technically 
# outside of the project region -- let's drop them from consideration here
WARN_RATIO <- round( sum( training_data$b == 0 & training_data[,GREENESS_INDEX] == 0 ) / length( training_data[,GREENESS_INDEX] ) , 2)
if( WARN_RATIO > 0 ){
  cat("DEBUG: Warning:", WARN_RATIO*100, "% of our training data contained useless values that we are dropping.",
      "Please review the source training data for errors.\n")
  training_data <- training_data[ which( !( training_data$b == 0 & training_data[,GREENESS_INDEX] == 0 ) ) , ]
}

# drop any lurking "type" values that aren't clearly coded
training_data <- training_data[grepl(training_data$type, pattern=paste(ACCEPTABLE_TYPES, collapse="|")),]
training_data <- training_data[,EXPLANATORY_COVS]

# split our dataset into training / testing components for ROC and evaluation purposes -- we will
# merge these back in to a full dataset for final model fitting later

testing <- sample(1:nrow(training_data), size=0.2*nrow(training_data))
training <- which(!(1:nrow(training_data) %in% testing))

cat("DEBUG: Fitting random forests\n")

# fit a throw-away random forests model 

test <- est_runtime(
  m_rf <- randomForest::randomForest(
    as.factor(type)~., 
    data=training_data[training,], 
    do.trace=F, 
    norm.votes=F, 
    importance=F, 
    ntree=3000
  )
)

cat("DEBUG: Finished RF fitting. Took (minutes): ", round(test/60), "\n")
cat("DEBUG: Checking for OOB error convergence. See: 'figures/erc_rf_oob_error_convergence.png'\n")

oob_err_rate_vs_trees <- round(
  diff(diff(m_rf$err.rate[,1])), 
  5 # just-in-case, only accept 5 decimal digits of precision
)

RF_BURNIN_THRESHOLD <- sd(oob_err_rate_vs_trees) 

# take an arbitrary value in the right tail of the 
# err distribution as our cut-off and re-fit our RF

N_TREES <- as.vector(
  round(
    quantile(
      which(abs(oob_err_rate_vs_trees) > RF_BURNIN_THRESHOLD), 
      probs=0.81)
  )
)

png("figures/erc_rf_oob_error_convergence.png", width=800, height=600)
  plot(
    oob_err_rate_vs_trees, 
    type="l", 
    col="DarkGrey",
    xlab="Number of Trees",
    ylab="OOB Error Rate",
    cex.lab=1.75,
    cex.axis=1.75,
    cex=1.75,
    ylim=c(mean(oob_err_rate_vs_trees)-(4*RF_BURNIN_THRESHOLD), mean(oob_err_rate_vs_trees)+(4*RF_BURNIN_THRESHOLD))
  )

  grid();
  
  abline(
    h=mean(oob_err_rate_vs_trees), 
    col="Black", 
    lty=2
  )
  
  abline(
    h=mean(oob_err_rate_vs_trees) + RF_BURNIN_THRESHOLD, 
    col="DarkRed", 
    lty=2
  )
  
  abline(
    h=mean(oob_err_rate_vs_trees) - RF_BURNIN_THRESHOLD, 
    col="DarkRed", 
    lty=2
  )
  
  abline(
    v=N_TREES, 
    col="DarkGreen", 
  )

  suppressMessages(dev.off())

cat("DEBUG: Refitting RF using", N_TREES, "trees to control for over-fit\n")

m_rf <- randomForest::randomForest(
    as.factor(type)~., 
    data=training_data[training,], 
    do.trace=F, 
    norm.votes=T, 
    importance=T, 
    ntree=N_TREES
)

test <- as.matrix(predict(m_rf, newdata=training_data[testing,], type="prob", norm.votes=T))
   n <- colnames(test)

test <- do.call(rbind, lapply(
  1:nrow(test),
  FUN=function(row){
    row <- matrix(round(as.numeric(test[row,]), 4), nrow=1)
    colnames(row) <- n 
    row <- cbind(row, sum(row[,grepl(colnames(row), pattern="erc_")]))
    row <- cbind(row, sum(row[,grepl(colnames(row), pattern="forest|grass|row_crop|small_grains|water")]))
    colnames(row) <- c(n, "erc", "other")
    return(row)
  }
))

test <- as.data.frame(cbind(
  test,
  predicted=as.vector(sapply(X=1:nrow(test), FUN=function(i) names(test[i,])[ which.max(test[i,])])),
  observed=training_data[testing,'type']
))

cat("DEBUG: Confusion Matrix for RF Model:\n");
print(round(m_rf$confusion,2));
cat("\n");

# test$predicted <- ifelse(grepl(test$predicted, pattern="erc"), yes="erc", no="other")

cat("DEBUG: Building ROC curves. See: 'figures/erc_roc_curve.png'\n")

# test$observed <- ifelse(grepl(test$observed, pattern="erc"), yes="erc", no="other")

# build a table of our per-class AUC values
roc_auc_values <- round(data.frame(
  erc=suppressMessages(as.numeric(auc(roc(ifelse(grepl(test$observed, pattern="erc"), yes="erc", no="other"), as.numeric(test$erc))))),
  grass=suppressMessages(as.numeric(auc(roc(ifelse(grepl(test$observed, pattern="grass"), yes="grass", no="other"), as.numeric(test$grass))))),
  small_grains=suppressMessages(as.numeric(auc(roc(ifelse(grepl(test$observed, pattern="small_grains"), yes="small_grains", no="other"), as.numeric(test$small_grains))))),
  row_crop=suppressMessages(as.numeric(auc(roc(ifelse(grepl(test$observed, pattern="row_crop"), yes="row_crop", no="other"), as.numeric(test$row_crop)))))
), 2)

png("figures/erc_roc_curve.png", width=800, height=600)
  
  plot(
    roc(ifelse(grepl(test$observed, pattern="erc"), yes="erc", no="other"), as.numeric(test$erc)),
    col="#86b300", # olive
    lwd=1.25,
    cex.lab=1.75,
    cex.axis=1.75,
    cex=1.75
  );

  grid(); grid();

  lines(
    suppressMessages(roc(ifelse(grepl(test$observed, pattern="grass"), yes="grass", no="other"), as.numeric(test$grass))),
    lwd=1.25,
    col="#0040ff" # blue
  );

  lines(
    suppressMessages(roc(ifelse(grepl(test$observed, pattern="row_crop"), yes="row_crop", no="other"), as.numeric(test$row_crop))),
    lwd=1.25,
    col="#c2f0c2" # light green
  );

  lines(
    suppressMessages(roc(ifelse(grepl(test$observed, pattern="small_grains"), yes="small_grains", no="other"), as.numeric(test$small_grains))),
    col="#ffdb4d" # mustard
  );

  legend(
    "bottomright", 
    pch = c(19, 19, 19, 19), 
    col = c("#86b300","#0040ff", "#c2f0c2", "#ffdb4d"), 
    legend = c(
      paste("eastern redcedar (AUC=",roc_auc_values$erc,")", sep=""), 
      paste("grass (AUC=",roc_auc_values$grass,")", sep=""), 
      paste("row crop (AUC=",roc_auc_values$row_crop,")", sep=""), 
      paste("small grains (AUC=",roc_auc_values$small_grains,")", sep="")
    ),
    cex=1.75
  )

  suppressMessages(dev.off())

cat("DEBUG: Refitting using full training dataset\n")
m_rf <- randomForest::randomForest(
    as.factor(type)~., 
    data=training_data, 
    do.trace=F, 
    norm.votes=T, 
    importance=T, 
    ntree=N_TREES
)
cat("DEBUG: Flushing confusion matrix to disk\n")

write.csv(
  as.data.frame((m_rf$confusion)), 
  "products/rf_erc_confusion_matrix.csv",
  row.names=F
)

# flush our session to disk for future review

cat("DEBUG: Saving workspace to disc for review\n")
r_data_file <- tolower(paste(
  "workflows/",
  "rf_erc_modeling_workflow_",
  gsub(format(Sys.time(), "%b %d %Y"), pattern = " ", replacement = "_"),
  ".rdata",
  sep = ""
))

save(
  compress = T,
  list = ls()[!grepl(ls(), pattern="explanatories")],
  file = r_data_file
)
