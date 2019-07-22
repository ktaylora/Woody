#!/bin/bash
echo ""
echo "## Random Forests Model Predictor with Raster Tiles Comprehension"
echo ""

if [ $1 ]; then
  DESTINATION_FILENAME=$1
else
  DESTINATION_FILENAME="rf_prediction.tif"
fi

if [ ! -d splits ]; then
  echo  "DEBUG: Generating tiles from Earth Engine source rasters (this could take awhile)"
  mkdir splits;
  bash ./gdal_splitter;
else
  echo "DEBUG: Existing ./splits/ directory found -- keeping contents as-is and running RF"
fi

echo "DEBUG: Staging scripts in ./splits directory"

mv *.rdata ./splits/
mv random_forest_predict*.R ./splits/

cd ./splits;

echo "DEBUG: Launching Rscript"
echo ""

Rscript random_forest_predict_across_geotiff_fragments.R $PWD;

if [ $? ]; then
  echo "DEBUG: Merging prediction fragments as : ../$DESTINATION_FILENAME"
  echo ""
  gdal_merge.py -o "../$DESTINATION_FILENAME" *_prediction.tif
else:
  echo "DEBUG: 'R' failed. Skipping gdal merge operation"
fi

echo "DEBUG: Cleaning-up"

mv random_forest*.R ../;
mv *.rdata ../;

cd ..;

echo "DEBUG: Done"
