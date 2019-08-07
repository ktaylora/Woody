suppressMessages(require(parallel))

ARGV <- commandArgs(trailingOnly=T)

if(length(ARGV)>0) {
  MERGED_PREDICTION_RASTERS_DIR = ARGV[1]
} else {
  MERGED_PREDICTION_RASTERS_DIR = '.'
}

cat("DEBUG: Looking for rf_prediction raster products in : ", MERGED_PREDICTION_RASTERS_DIR, "\n")

prediction_rasters <- list.files(
  MERGED_PREDICTION_RASTERS, 
  pattern="rf_prediction_.*.tif", 
  full.names=T
)

if( length(predction_rasters) == 0 ){
  cat("DEBUG: No rf_predction raster products found -- quitting\n")
  quit("no", 1)
}

cl <- parallel::makeCluster(parallel::detectCores()*0.8, outfile='')

result <- parallel::parLapply(
  cl=cl,
  fun=function(r){
    require(raster)
    cat("DEBUG: Processing image: ", r, "\n")
    dst_high_file <- gsub(r, pattern="rf_prediction", replacement="high_rf_prediction")
     dst_low_file <- gsub(r, pattern="rf_prediction", replacement="low_rf_prediction")
    raster::writeRaster(raster(r) == 1, dst_high_file, progress='text', datatype='INT1U', overwrite=T)
    raster::writeRaster(raster(r) == 2, dst_low_file, progress='text', datatype='INT1U', overwrite=T)
  }
)

cat("DEBUG: Cleaning-up cluster interface\n")
parallel::stopCluster(cl); rm(cl); gc()

cat("DEBUG: one")
quit("no", 0)
