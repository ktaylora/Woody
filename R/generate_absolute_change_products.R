
suppressMessages(require(rgdal))
suppressMessages(require(raster))
suppressMessages(require(rgeos))

#
# Runtime parameters
#

CELL_RESOLUTION     = 10
TIMESTEP_FOR_SERIES = 3
DECIMAL_PRECISION   = 4
LOGGING = T

#
# Local functions
#
logger <- function(x) if(LOGGING) cat(paste("DEBUG:", x, "\n", sep=" "))

est_timeseries_rate_of_change <- function(current=NA, past=NA){
    rate_of_change <- current

    logger("Estimating area of each USNG unit")
    area <- rgeos::gArea(spTransform(past, CRS(projection("+init=epsg:2163"))), byid=T)

    logger("Estimating rate-of-change for first time-series")
    result <- ( (CELL_RESOLUTION^2) * (current@data - past@data) ) / area
        result <- (result / TIMESTEP_FOR_SERIES ) # annualize
        result <- round(result, DECIMAL_PRECISION)

    rate_of_change@data <- result
    return(rate_of_change)
}

#
# Main
#

mesq_products <- list(
    bcr_18=list( c("/home/ktaylora/Incoming/offyear_study/products/","mesq_bcr_18_2015_2017"), c("/home/ktaylora/Incoming/offyear_study/products/","mesq_bcr_18_2012_2014") ),
    bcr_19=list( c("/home/ktaylora/Incoming/offyear_study/products/","mesq_bcr_19_2015_2017"), c("/home/ktaylora/Incoming/offyear_study/products/","mesq_bcr_19_2012_2014") )    
)

logger("Reading GeoJSON files from disc")

result <- lapply(
    X=mesq_products,
    FUN=function(bcr){
        current <- readOGR(bcr[[1]][1], bcr[[1]][2], stringsAsFactors=F, verbose=F)
          current$both <- rowSums(current@data)
        past <- readOGR(bcr[[2]][1], bcr[[2]][2], stringsAsFactors=F, verbose=F)
          past$both <- rowSums(past@data)
        return(est_timeseries_rate_of_change(current, past))
    }
)

names(result) <- c("bcr_18","bcr_19")

logger("Writing Geometries to disc as shapefiles")
writeOGR(result$bcr_18, "/home/ktaylora/Incoming/offyear_study/products/","mesq_bcr_18_rate_of_change_2012_2017", driver="ESRI Shapefile", overwrite=T)
writeOGR(result$bcr_19, "/home/ktaylora/Incoming/offyear_study/products/","mesq_bcr_19_rate_of_change_2012_2017", driver="ESRI Shapefile", overwrite=T)

