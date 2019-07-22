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

band_names  <- c(
  'R','G','B','N','NDVI','NDVI_3','NDVI_9',
  'NDVI_18','NDVI_36','NDVI_SD_3','NDVI_SD_9',
  'NDVI_SD_18','NDVI_SD_36','ELEV','ELEV_SD_11',
  'ASPECT','SLOPE'
)  

cat("DEBUG: Using band names:",paste(band_names, collapse=", "),"\n")

cat("DEBUG: Starting workers...\n")

cl <- parallel::makeCluster(
  floor(parallel::detectCores()*0.8), 
  outfile=''
)

parallel::clusterExport(
  cl, 
  varlist=c("m_rf", "geotiff_segments", "band_names")
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
      w <- list.files('.', pattern="tif$")
      cat(
        "DEBUG: Percent complete -- ", 
        round( length(w[grepl(w, pattern="_prediction")]) / length(w[!grepl(w, pattern="_prediction")]) , 2),
        "% \n"
      )
      return(TRUE)
    }
))

cat("DEBUG: DONE\n")

if( class(result) == "try-error" ){
  cat("DEBUG: (But we encountered some problems)\n")
  quit(save = "no", status = FALSE)
} else {
  quit(save = "no", status = TRUE)
}



