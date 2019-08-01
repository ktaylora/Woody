#!/bin/bash
ls -1 *.tif | awk '{ print "gdal_retile.py -ps 2500 2500 -targetDir splits/ " $1 }' | /bin/bash
