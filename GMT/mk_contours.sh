#! /bin/sh

# Arguments.
grdfile=$1
outfile=$2

# Set parameters.
region=-R170.0/180.0/-45.0/-36.0
proj=-JM175.0/-40.5/6i
annot=-Ba2f1g1
clipfile=bounds_wgs84.txt
dry=grey90
clabels=-Gl176/-36/179/-40,173/-39/177/-42,171/-42/173/-44.5
cint=-Ccontoursb.txt
cannot=-A40+f28
clines=-Wcthick,black,solid
clevels=-L-300/0
blines=-W0.75p,red,solid

# Set paper size.
gmtset PAPER_MEDIA=a4+
gmtset GRID_PEN_PRIMARY=0.25p/grey
gmtset GRID_PEN_SECONDARY=0.25p/grey

# Draw map and coastline.
pscoast $annot $proj $region -Df -Ia -Na -K -P -V -G$dry -W > $outfile

# Generate contours.
grdcontour $grdfile $proj $region $cint $annot $cannot $clevels $clines $clabels -K -O -P -V >> $outfile
grdcontour $grdfile $cint $clevels -Dhiku_cont.txt

awk '{if ($3 < 0) print $0}' hiku_cont.txt > tmp
mv tmp hiku_cont.txt

cat hiku_cont.txt hik_trench_surf_wgs84_dec.txt > hiku_cont_wtrench.txt

# Draw boundary of valid points.
psxy $clipfile $proj $region $annot -A -O -P -V $blines >> $outfile
