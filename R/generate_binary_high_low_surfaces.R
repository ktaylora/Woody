suppressMessages(require(parallel))

ARGV <- commandArgs(trailingOnly=T)

if(length(ARGV)>0) {
  MERGED_PREDICTION_RASTERS_DIR = ARGV[1]
} else {
  MERGED_PREDICTION_RASTERS_DIR = '.'
}

cat("DEBUG: Looking for rf_prediction raster products in : ", MERGED_PREDICTION_RASTERS_DIR, "\n")

prediction_rasters <- list.files(
  MERGED_PREDICTION_RASTERS_DIR, 
  pattern="rf_prediction_.*.tif", 
  full.names=T
)

prediction_rasters <- prediction_rasters[!grepl(prediction_rasters, pattern="high|low")]
prediction_rasters <- prediction_rasters[grepl(prediction_rasters, pattern="[.]tif$")]

if( length(prediction_rasters) == 0 ){
  cat("DEBUG: No rf_predction raster products found -- quitting\n")
  quit("no", 1)
}

cl <- parallel::makeCluster(parallel::detectCores()*0.8, outfile='')

result <- parallel::parLapply(
  cl=cl,
  X=prediction_rasters,
  fun=function(r){
    suppressMessages(require(raster))
    cat("DEBUG: Processing image: ", r, "\n")
    if(grepl(r, pattern="mesq")){
      high_class <- 3
      low_class <- 4
      forest_class <- 1
    } else {
      high_class <- 1
      low_class <- 2
      forest_class <- 3 
    }

    dst_high_file <- gsub(r, pattern="rf_prediction", replacement="high_rf_prediction")
    dst_low_file <- gsub(r, pattern="rf_prediction", replacement="low_rf_prediction")
    dst_forest_file <- gsub(r, pattern="rf_prediction", replacement="mixed_forest_rf_prediction")
    
    if(file.exists(dst_high_file)){
      cat("DEBUG: Skipping found file:", dst_high_file,"\n");
    } else {
      raster::writeRaster(raster(r) == high_class, dst_high_file, progress='text', datatype='INT1U', overwrite=T)
    } 
    
    if(file.exists(dst_low_file)){
      cat("DEBUG: Skipping found file:", dst_low_file,"\n");
    } else {
      raster::writeRaster(raster(r) == low_class, dst_low_file, progress='text', datatype='INT1U', overwrite=T)
    }
    
    if(file.exists(dst_forest_file)){
      cat("DEBUG: Skipping found file:", dst_forest_file,"\n");
    } else {
      raster::writeRaster(raster(r) == forest_class, dst_forest_file, progress='text', datatype='INT1U', overwrite=T)
    }
  }
)

cat("DEBUG: Cleaning-up cluster interface\n")
parallel::stopCluster(cl); rm(cl); gc()

cat("DEBUG: Done")
quit("no", 0)
