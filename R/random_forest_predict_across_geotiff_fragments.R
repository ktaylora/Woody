cat("### RF Predictor Script\n")

suppressMessages(require(parallel))

argv <- commandArgs(trailingOnly=T)

# Runtime Options
GEOTIFF_EE_IMAGE_SEGMENTS_DIR = NA
RF_MODEL_OBJECT_PATH = NA

if( is.na(GEOTIFF_EE_IMAGE_SEGMENTS_DIR) ){
  if( !is.na(argv[1]) > 0 ){
    GEOTIFF_EE_IMAGE_SEGMENTS_DIR = argv[1]
  } else {
    GEOTIFF_EE_IMAGE_SEGMENTS_DIR = '.'
  }
}

if( is.na(RF_MODEL_OBJECT_PATH) ){
  if( !is.na(argv[2]) ){
    RF_MODEL_OBJECT_PATH = argv[2]
  } else {
    RF_MODEL_OBJECT_PATH = '.'
  }
}

cat("DEBUG: Loading RF model object from disc\n")

load(
  list.files(
    RF_MODEL_OBJECT_PATH,
    pattern="rf.*.modeling_workflow",
    full.names=T)[1]
)

cat("DEBUG: Reading GeoTIFF raster segments from disc\n")

geotiff_segments <- list.files(
    GEOTIFF_EE_IMAGE_SEGMENTS_DIR, 
    pattern="tif$", 
    full.names=T
) 

# don't try and use any lurking _prediction rasters if they are in the CWD
geotiff_segments <- geotiff_segments[!grepl(geotiff_segments, pattern="_pred")]

# Figure out the band names using gdalinfo from our original EE GeoTIFF files
# Otherwise, just guess the band names using sane defaults and see if any errors 
# are raised

band_names <- list.files(paste(GEOTIFF_EE_IMAGE_SEGMENTS_DIR, "..", sep="/"), pattern="tif$", recursive=F, full.names=T)[1]

if( length(band_names) > 0 ) {
  band_names <- system(paste('gdalinfo', band_names),  intern = TRUE)
  band_names <- sapply(strsplit(band_names[grepl(band_names, pattern="Desc")], split=" = "), function(i) i[2])
} 

if (  length(band_names) > 0 ){
  cat("DEBUG: Band names available:",paste(band_names, collapse=", "),"\n")
} else {
  cat("DEBUG: Didn't have any source rasters available to determine band names -- quitting\n")
  quit("no", 1)
}

# Determine if we need to re-fit our random forests to accomodate shifts in bands for this
# NAIP interval
if ( sum(grepl(band_names, pattern="G_|GF_")) > 0 ) {
  cat("DEBUG: Detected alternate green fraction used in input bands -- refitting RF model\n")
  training_data_alt_gf <- strsplit(TRAINING_DATA, split="/")[[1]]
  training_data_alt_gf <- paste(training_data_alt_gf[1:(length(training_data_alt_gf)-1)], collapse="/")
  training_data_alt_gf <- list.files(training_data_alt_gf, pattern=".*alt_green.*.[.]geojson", full.names=T)
  training_data_alt_gf <- suppressWarnings(suppressMessages(rgdal::readOGR(training_data_alt_gf, require_geomType='wkbPoint', verbose=F)))
  
  coords <- sp::spTransform(training_data_alt_gf, sp::CRS(raster::projection("+init=epsg:4326")))@coords
    colnames(coords) <- c("lon","lat")
  training_data_alt_gf@data <- cbind(training_data_alt_gf@data, coords)
  
  colnames(training_data_alt_gf@data) <- tolower(colnames(training_data_alt_gf@data))
  colnames(training_data_alt_gf@data) <- gsub(colnames(training_data_alt_gf@data), pattern="g_", replacement="gf_")
  
  m_rf <- randomForest::randomForest(
      as.factor(type)~r + b + g + gf + gf_3 + gf_18 + gf_36 + gf_sd_3 + gf_sd_9 + gf_sd_18 + gf_sd_36 + lat + lon, 
      data=training_data_alt_gf@data,
      do.trace=F, 
      norm.votes=T, 
      importance=T, 
      ntree=N_TREES
  )
}

cat("DEBUG: Starting workers...\n")

cl <- parallel::makeCluster(
  floor(parallel::detectCores()*0.6), 
  outfile=''
)

parallel::clusterExport(
  cl, 
  varlist=c("m_rf", "geotiff_segments", "band_names", "GEOTIFF_EE_IMAGE_SEGMENTS_DIR")
)

result <- try(parallel::parLapply(
    cl=cl,
    X=geotiff_segments,
    fun=function(segment){
    
      if( file.exists(gsub(segment, pattern="[.]tif$", replacement="_prediction.tif")) ){
        cat(
          "DEBUG: Skipping segment (prediction already found in CWD) --",
          gsub(segment, pattern="[.]tif$", replacement="_prediction.tif"),
          "\n"
        )
        return(FALSE)
      }

      suppressMessages(require(raster))
      suppressMessages(require(randomForest))

      cat("DEBUG: Processing segment:", segment, "\n")
      w <- which(geotiff_segments %in% segment)[1]

      prediction_surface <- suppressMessages(try(raster::stack(segment), silent=TRUE))

      if(class(prediction_surface) == "try-error"){
        cat(
          "DEBUG: Error reading segment from disc (skipping):",
          segment,
          "\n"
        )
        return(FALSE)
      }

      if( raster::nlayers(prediction_surface) != length(band_names) ){
        cat(
          "DEBUG: Error mis-match between number of layers and band names (quitting):",
          segment,
          "\n"
        )
        return(FALSE)
      }

      cat("DEBUG: Calculating lat/lon for cell centroids for segment",w,"\n")
      coords <- sp::spTransform(
          raster::rasterToPoints(prediction_surface, spatial=T),
          sp::CRS(raster::projection("+init=epsg:4326"))
      )@coords
      
      cat("DEBUG: Building lat/lon raster surface for segment",w,"\n")
      
      template <- prediction_surface[[1]]
        values(template) <- coords[,1]

      prediction_surface <- addLayer(
          prediction_surface,
          template
      )

      template <- prediction_surface[[1]]
        values(template) <- coords[,2]

      prediction_surface <- addLayer(
          prediction_surface,
          template
      )

      rm(template)

      cat("DEBUG: Setting band names for segment",w,"\n")

      band_names <- c(band_names, "lon", "lat")
      names(prediction_surface) <- tolower(band_names)

      cat("DEBUG: Running random forests against segment",w,"\n")
      raster::predict(
          prediction_surface, 
          m_rf, 
          filename=gsub(segment, pattern="[.]tif$", replacement="_prediction.tif"),
          dataType='INT1U',
          format='GTiff'
      )
      cat("DEBUG: Finished random forests run for segment",w,"\n")
      w <- list.files(GEOTIFF_EE_IMAGE_SEGMENTS_DIR, pattern="tif$")
      cat(
        "DEBUG: Percent complete -- ", 
        round( length(w[grepl(w, pattern="_prediction")]) / length(w[!grepl(w, pattern="_prediction")]) , 2),
        "% \n"
      )
      return(TRUE)
    }
))

cat("DEBUG: Cleaning up our 'R' cluster interface\n")
parallel::stopCluster(cl); rm(cl); gc()

cat("DEBUG: DONE\n")

if( class(result) == "try-error" ){
  cat("DEBUG: (But we encountered some problems)\n")
  quit(save = "no", status = TRUE)
} else {
  quit(save = "no", status = FALSE)
}



