#! /bin/sh

# Arguments.
xyzfile=$1
grdfile=$2
incr=$3

# Set parameters.
region=-R170.0/180.0/-45.0/-36.0
incrarg=-I$incr

# Grid data.
xyz2grd $xyzfile -H1 -G$grdfile $incrarg $region -V
