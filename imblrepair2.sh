#!/bin/csh -f

#WARNING: DO NOT USE if zeros are likely to be real data!
#This function assumes zeros are targets for replacement, &
# that bad pixels are either isolated or in single-file lines.
# Not be used for inpainting large groups of contiguous blanks.

set ifile = $1
set ofile = $ifile #default is to replace file input file
set sweep = 'n' #if y and $ofile exists, remove $ifile
set axis = 'x0' #choices are x0, y0, and xy or yx

foreach par ( $argv[*] )
    set check=`echo $par | awk -F= '{print NF}'`
    if ( "$check" >= 2 ) set $par
end
echo $argv[*]

goto start
start:
set tempr=$ifile:r
echo "imblrepair.sh: Beginning repairs to $ifile..."

rm -r $tempr.??x
rm -r $tempr.??y
rm -r $tempr.???.mvd
rm -r $tempr.avg*
rm -r $tempr.?avg
rm -r $tempr.*zed
set norz = `immask in="$ifile"`
if ($norz[5] < $norz[8]) then
    imblr in="$ifile" out="$tempr.zed"
else
    cp -r $ifile $tempr.zed
endif

if ($axis == "x0" || $axis == "xy" || $axis == "yx") then
    cp -r $ifile $tempr.mnx
    cp -r $ifile $tempr.plx
    set xrefpix=`gethd in=$tempr.zed/crpix1`
    puthd in="$tempr.mnx/crpix1" value=`calc "$xrefpix-1"` type=double
    puthd in="$tempr.plx/crpix1" value=`calc "$xrefpix+1"` type=double
    regrid in="$tempr.mnx" tin="$tempr.zed" out="$tempr.mnx.mvd"
    regrid in="$tempr.plx" tin="$tempr.zed" out="$tempr.plx.mvd"
    maths exp="(<$tempr.plx.mvd>+<$tempr.mnx.mvd>)/2" mask="<$tempr.zed>.eq.0" out="$tempr.avgx"
endif

if ($axis == "y0" || $axis == "xy" || $axis == "yx") then
    cp -r $ifile $tempr.mny
    cp -r $ifile $tempr.ply
    set yrefpix=`gethd in=$tempr.zed/crpix2`
    puthd in="$tempr.mny/crpix2" value=`calc "$yrefpix-1"` type=double
    puthd in="$tempr.ply/crpix2" value=`calc "$yrefpix+1"` type=double
    regrid in="$tempr.mny" tin="$tempr.zed" out="$tempr.mny.mvd"
    regrid in="$tempr.ply" tin="$tempr.zed" out="$tempr.ply.mvd"
    maths exp="(<$tempr.ply.mvd>+<$tempr.mny.mvd>)/2" mask="<$tempr.zed>.eq.0" out="$tempr.avgy"
endif

if ($axis == "x0") then
    imblr in="$tempr.avgx" out="$tempr.avg.zed"

else if ($axis == "y0") then
    imblr in="$tempr.avgy" out="$tempr.avg.zed"

else if ($axis == "xy" || $axis == "xy") then
    imcomb in="$tempr.avgx $tempr.avgy" out="$tempr.navg"
    imblr in="$tempr.navg" out="$tempr.avg.zed"
endif

if ($ofile == $ifile) echo "overwriting $ifile..."
rm -r $ofile
#the following only works b/c the 2 images are complements
#pixels to be replaced = 0 in $tempr.zed, & all pixels
# EXCEPT replacement values = 0 in $tempr.avg.zed. x+0=x
maths exp="<$tempr.zed>+<$tempr.avg.zed>" out="$ofile"

#echo "Begin clean-up..."
rm -r $tempr.??x
rm -r $tempr.??y
rm -r $tempr.???.mvd
rm -r $tempr.avg*
rm -r $tempr.?avg
rm -r $tempr.*zed

if ($sweep == "y" && $ifile != $ofile) then
    echo "Sweep = y: deleting $ifile"
    rm -r $ifile
endif

set mssg2 = `immask in="$ofile"`
if ($mssg2[5] < $mssg2[8]) then
    echo "ERROR: Blank pixel repair failed"
    exit 1
else
    echo "Blank pixel repair completed."
endif
end:
