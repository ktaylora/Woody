#!/usr/bin/env python3

"""
Random Forest Predict Across GeoTIFF Fragments

This script will accept a directory containing a series of GeoTIFF files
produced by Earth Engine (i.e., from a raster export to Google Drive). 
It will split each of the GeoTIFF files into segments that are digestable
by R and then call 'R' using an RData file to predict woody species hotspots
for each segment (saved as a raster file). The script will then merge the 
results together into a single (small) GeoTIFF file with byte precision.
"""

import sys
import os
import subprocess

import glob

from beatbox.vector import Vector
from beatbox.raster import Raster

#
# Be verbose by default
#

import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

#
#
# Runtime arguments
#

import argparse

parser = argparse.ArgumentParser(description='Use a Random Forests classifier to predict across Earth Engine GeoTIFF fragments')
parser.add_argument('--dir', metavar="d", dest='geotiff_dir_path', default='.',help='specify the full path to a directory containing EE GeoTIFF files')
parser.add_argument('--target_dir', metavar="t", dest='splits_target_dir', default='splits',help='specify the local path to a target directory to house our split EE GeoTIFF files')
parser.add_argument('--outfile', metavar="o", dest='outfile', default='rf_prediction.tif', help="specify the filename for our merged random forests prediction raster")

args = parser.parse_args()

#
# Local functions
#

def dig_path(**kwargs):
  PATH = kwargs.get('path', '.')
  EXTENSION_PATTERN = kwargs.get('extension_pattern', '*.tif')
  return glob.glob(os.path.join(PATH, EXTENSION_PATTERN), recursive=False)

def gdal_retile(**kwargs):
  """
  Wrapper for gdal_retile.py
  """
  GEOTIFF_FILES = kwargs.get('geotiff_files', None)
  TARGET_DIR = kwargs.get('target_dir', args.splits_target_dir)
  PIXEL_SIZE = kwargs.get('pixel_size', 2500)

  if GEOTIFF_FILES is None:
    raise ValueError("geotiff_files= argument cannot be empty")
  elif len(GEOTIFF_FILES) == 0:
    raise ValueError("geotiff_files= argument cannot be empty")

  if os.path.exists(TARGET_DIR):
    logger.debug("Existing gdal_retile target directory found : " + \
      TARGET_DIR)
    logger.debug("Assuming you want to use existing tiles and " + \
      "skipping GeoTIFF file splitting operations")
    return
  else:
    os.mkdir(TARGET_DIR)

  progress = 0

  for segment in GEOTIFF_FILES:
    subprocess.run(["gdal_retile.py", "-ps", str(PIXEL_SIZE), str(PIXEL_SIZE), "-targetDir", TARGET_DIR, segment])
    progress += 1
    logger.debug("Overall % complete: "+str(round(progress/len(GEOTIFF_FILES), 2)))

def rf_predict(**kwargs):
  """
  Wrapper for Rscript
  """
  R_SCRIPT_SOURCE_PATH = kwargs.get('r_script_path',
    os.path.join('.','random_forest_predict_across_geotiff_fragments.R'))
  GEOTIFF_SEGMENTS_DIR = kwargs.get('geotiff_segments_path',
    os.path.join('.',args.splits_target_dir))
  logger.debug("Calling: Rscript " + R_SCRIPT_SOURCE_PATH + ' ' + GEOTIFF_SEGMENTS_DIR)
  subprocess.run(['Rscript', R_SCRIPT_SOURCE_PATH,
    GEOTIFF_SEGMENTS_DIR])

def gdal_merge(**kwargs):
  """
  Wrapper for gdal_merge
  """
  GEOTIFF_SEGMENTS_DIR = kwargs.get('geotiff_segments_path', args.splits_target_dir)
  
  logger.debug("Calling gdal_merge.py -ot Byte -o " + args.outfile + " using our prediction segments in : "+ GEOTIFF_SEGMENTS_DIR)
  
  merge_runner = ["gdal_merge.py",'-ot','Byte', '-o', args.outfile]
  merge_runner = merge_runner + glob.glob(os.path.join(GEOTIFF_SEGMENTS_DIR,'*_prediction.tif'))

  subprocess.run(merge_runner)

def gdal_aggregate_sum_by_usng_unit(**kwargs):
  """
  Wrapper for gdal/beatbox that 
  """
  VECTOR_GEOMETRIES = kwargs.get('usng_units', None)
  GEOTIFF_PREDICTION_RASTER = kwargs.get('geotiff_segment', None)
  TARGET_VALUES = kwargs.get('target_values', [1,2])
  
  vector_units = Vector(VECTOR_GEOMETRIES).to_geodataframe()
  prediction_raster = Raster(GEOTIFF_PREDICTION_RASTER).to_georaster()

#
# Main
#

if __name__ == "__main__":

  print("\n### Random Forests Model Predictor with Raster Tiles Comprehension\n")

  logger.debug("Parsing input directory for files: " + os.path.join(args.geotiff_dir_path,'*.tif'))
  
  geotiff_files = dig_path(path=args.geotiff_dir_path)
  logger.debug("Found "+str(len(geotiff_files))+" raster files")
  logger.debug("Wrapping gdal_retile.py using Earth Engine segments")

  try:
    gdal_retile(geotiff_files=geotiff_files,
      target_dir=os.path.join(args.geotiff_dir_path, args.splits_target_dir))
  except ValueError as e:
    logger.debug("Cannot pass gdal_retile.py an empty list of GeoTIFFs -- this shouldn't happen; quitting")
    sys.exit(1)

  if len(dig_path(path=os.path.join(args.geotiff_dir_path, args.splits_target_dir))) == 0:
    logger.debug("gdal_retile target 'splits' directory didn't appear to contain any GeoTIFF files -- this shouldn't happen; quitting")
    sys.exit(1)
  
  rf_predict(geotiff_segments_path=os.path.join(args.geotiff_dir_path, args.splits_target_dir))

  if len(dig_path(path=os.path.join(args.geotiff_dir_path, args.splits_target_dir), extension_pattern='*_prediction.tif')) < 1 :
    logger.debug("Couldn't find any *_prediction.tif products in the 'splits' directory -- this shouldn't happen; quitting")
    sys.exit(1)
  else:
    gdal_merge(geotiff_segments_path=os.path.join(args.geotiff_dir_path, args.splits_target_dir))

  logger.debug("DONE")
  sys.exit(0)
