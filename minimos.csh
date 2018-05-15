#!/bin/csh -f
# Combine multiple map/cube files.  Adapted from SOD's original Merge.csh for 
# just 3 cubes.  Assumes axes are in lbv order.
# 2009Jun17: Original -Stefan O'Dougherty
# 2009Aug13: Fixed template descriptor calculations, other minor tweaks -PJB
# 2011Oct12: Corrected newLonPixels calculation for cosine(Dec or Lat)
#	and the fact that default angle units in regrid are now degrees -PJB
# 2014May02: Rewrote for arbitrary list of files -PJB
# 2014Sep14: Added merging on z axis, allowed for any miriad-style suffix 
#	in file names -PJB
# 2015Jul14: Resolved longstanding imcomb error re common filename parts,
#	and also z-axis/imstat issues traced to mean pixel error -PJB

set datadir = .			# Where the files are located; '.' =  cwd
# Will have to set datadir manually either at command line or here 
#if (! -e maplist) then		# If no list given, merge all fits files
#    ls *.fits >! maplist
#endif
set files = `cat ./maplist` #maplist is a file, not a dir
# Note: cat command followed by file name prints contents of that file
# so this assigns $files as the variable name of the list in maplist 
set filect = 0                  #count of usable files
# Warnings: $out1 MUST be part of the beginning of each $filer below,
# including any assumed punctuation, but NOT TRAILING DOT!  $out2 MUST
# INCLUDE LEADING DOT if required.
set out1 = dr3			#\ common parts of output filenames, but NOT
set out2 = .12co.mom0		#/ incl fits or mir suffixes!
set delmir = y			#\
set delregr = y			# > Cleanup options
set delcomb = y			#/

#Notes 1:
# set assigns variables local to this shell
# setenv assigns variables for this shell and all subshells

###OVERRIDE ABOVE PARAMETERS IF GIVEN ON COMMAND LINE
foreach par ( $argv[*] )
    set check=`echo $par | awk -F= '{print NF}'`
    if ( "$check" >= 2 ) set $par
end

#Notes 2:
# Awk is like a stripped-down numpy - it handles floats & their calculations
# see http://www.grymoire.com/Unix/Awk.html for syntax & built-ins  

goto start
start:
#check if both start & end actually need colons after them

#Notes 3:
# colons are string-editing operators.
# *see https://en.wikibooks.org/wiki/C_Shell_Scripting/Modifiers
# Comparison operators are similar to python.
# Logical operators include || (OR) and && (AND).
# & is bitwise AND.
# CSH also has inquiry operators exclusively for files & their names.
# *see http://docstore.mik.ua/orelly/unix/upt/ch47_04.htm, section 47.4.3.5

echo ""
echo "===================================="
echo "Start with fits to miriad conversion"

#Notes 4:
# here, fits is a miriad command to convert b/w fits & miriad image formats
# *see http://www.atnf.csiro.au/computing/software/miriad/doc/fits.html
# Unix treats miriad files as directories, each containing:
#	-1 history file (plain text)
#	-1 header file (binary)
#	-1 mask file (binary)
#	-1 image/data file (binary)
# but miriad treats whole directory as 1 file.
# Commands in `` force in-line execution
# *see http://www.atnf.csiro.au/computing/software/miriad/taskindex.html

foreach file ($files) #recall $files = `cat ./maplist` 
    set filee = $file:e #extension
    set filer = $file:r #file name root
    if ($filee == 'fits' || -d $file) then
    #if $file extension is fits OR $file is a directory  
	if (-e $filer.fits && ! -d $file) then
	#if $file of given name exists AND is not a directory
	    rm -rf $filer.mir #remove existing version of $file
	    fits in=$filer.fits out=$filer.mir op=xyin #miriad fits conversion
	    set starttype = fits 
	else
	    if (-d $file) then #why not test if "-e $filer.mir" instead?
		echo "$file is an existing directory, assuming miriad-style"
		set starttype = mir
		set delmir=n	# don't delete it!!!
	    endif
	endif
	set filect = `calc -i "($filect)+(1)"` #increments $filect
    else
	echo "$file must be *.fits or miriad-style, skipping"
    endif
end
echo "Found $filect usable map files" #prints count of usable files


echo ""
echo "=========================="
echo "Now set up the mosaic area"
#will need this to combine areas of mopra files in Regions/R*
foreach file ($files)
    set filer = $file:r #assign template its root name
    set filee = $file:e #assign template extension
    if ($starttype == 'fits') set filee = mir #set output file type to miriad
    set first = `echo "$file $files" | awk '{if ($1 == $2) print 1; else print 0}'` 
    # $1 refs 1st parameter ($file), $2 is 2nd param ($files), etc
    # 
    set filedim = `gethd in=$filer.$filee/naxis` #get image dims from header
    #echo ""
    #echo $filer

    #Accumulate X axis descriptors (I know all these terms)
    set ctyp1 = `gethd in=$filer.$filee/ctype1`
    set npix1 = `gethd in=$filer.$filee/naxis1`
    set delt1 = `gethd in=$filer.$filee/cdelt1`
    set rpix1 = `gethd in=$filer.$filee/crpix1`
    set rval1 = `gethd in=$filer.$filee/crval1`
    #echo $ctyp1 $npix1 $delt1 $rpix1 $rval1
    set lmaxf = `calc "($rval1)+($delt1)*(($rpix1)-($npix1))"`
    set lminf = `calc "($rval1)+($delt1)*(($rpix1)-(1))"`
    #echo $lmaxf $lminf #$files
    if ($first) then
	set lmax = $lmaxf
	set lmin = $lminf
    else
	set lmax = `echo "$lmax $lmaxf" | awk '{if ($1 > $2) print $1; else print $2}'`
	set lmin = `echo "$lmin $lminf" | awk '{if ($1 < $2) print $1; else print $2}'`
    endif
    #echo $lmax $lmin

    #Accumulate Y axis descriptors
    set ctyp2 = `gethd in=$filer.$filee/ctype2`
    set npix2 = `gethd in=$filer.$filee/naxis2`
    set delt2 = `gethd in=$filer.$filee/cdelt2`
    set rpix2 = `gethd in=$filer.$filee/crpix2`
    set rval2 = `gethd in=$filer.$filee/crval2`
    #echo $ctyp2 $npix2 $delt2 $rpix2 $rval2
    set bmaxf = `calc "($rval2)-($delt2)*(($rpix2)-($npix2))"`
    set bminf = `calc "($rval2)-($delt2)*(($rpix2)-(1))"`
    #echo $bmaxf $bminf
    if ($first) then
	set bmax = $bmaxf
	set bmin = $bminf
    else
	set bmax=`echo "$bmax $bmaxf" | awk '{if ($1 > $2) print $1; else print $2}'`
	set bmin=`echo "$bmin $bminf" | awk '{if ($1 < $2) print $1; else print $2}'`
    endif
    #echo $bmax $bmin
end

#Use last map to define coord types
set cbase1 = `echo "$ctyp1" | cut -d - -f 1 -`
set cbase2 = `echo "$ctyp2" | cut -d - -f 1 -`
if ($cbase1 == 'RA' || $cbase1 == 'DEC' || $cbase1 == 'GLON' || $cbase1 == 'GLAT') then
    set xunit = radians
endif
if ($cbase1 == 'VELO') then
    set xunit = km/sec
endif
if ($cbase2 == 'RA' || $cbase2 == 'DEC' || $cbase2 == 'GLON' || $cbase2 == 'GLAT') then
    set yunit = radians
endif
if ($cbase2 == 'VELO') then
    set yunit = km/sec
endif
# Insert any other unit definitions here
echo "Maximum x coord found = $lmax $xunit"
echo "Minimum x coord found = $lmin $xunit"
echo "Maximum y coord found = $bmax $yunit"
echo "Minimum y coord found = $bmin $yunit"
endif


echo ""
echo "================="
echo "Template creation"
# MUST CHANGE THIS FROM MAXIMUM EXTENT OF ALL FILES TO
# MAX AREA OF ALL MOPRA FILES IN REGION (w/ slight padding)
# If filename includes 12co, 13co, or c18o, use mom0  
# If filename includes hco, use mom0 files containing 'main'
# If region files are disjoint, draw perimeter around all contained patches
# For multiple velo integrals, take max extent of all combined 

#X axis calculations: if angles, convert radians to degrees for regrid
if ($xunit == 'radians') then
    set newXrefval = `calc "(($lmax)+($lmin))*(180/pi)/2"`
    set newXrange  = `calc "abs(($lmax)-($lmin))*(180/pi)"`
    set newXdelta  = `calc "($delt1)*(180/pi)"`
    set xunit = degrees
else
    set newXrefval = `calc "(($lmax)+($lmin))/2"`		# This will be at *MEAN* pixel!!!
    #don't median-combine!
    set newXrange  = `calc "abs(($lmax)-($lmin))"`		# Equiv. to (pixels-1)*delta
    set newXdelta  = `calc "($delt1)"` #should all be the same?
endif
if ($yunit == 'radians') then
    # cosine(declination/latitude) correction
    set newXpixels = `calc -i "($newXrange/abs($newXdelta)*cos($rval2))+(1)"`
else
    set newXpixels = `calc -i "($newXrange/abs($newXdelta))+(1)"`
endif
set newXrefpix = `calc "($newXpixels+1)/2"`			# *Mean* pixel number (can be fractional)
echo "newXrefval = $newXrefval $xunit"
echo "newXrefpix = $newXrefpix"
echo "newXrange  = $newXrange $xunit"
echo "newXdelta  = $delt1 $xunit = last one found"
echo "newXpixels = $newXpixels"

#Y axis calculations
if ($yunit == 'radians') then
    set newYrefval = `calc "(($bmax)+($bmin))*(180/pi)/2"`
    set newYrange  = `calc "abs(($bmax)-($bmin))*(180/pi)"`
    set newYdelta  = `calc "($delt2)*(180/pi)"`
    set yunit = degrees
else
    set newYrefval = `calc "(($bmax)+($bmin))/2"`
    set newYrange  = `calc "abs(($bmax)-($bmin))"`
    set newYdelta  = `calc "($delt2)"`
endif
set newYpixels = `calc -i "($newYrange/abs($newYdelta))+(1)"`
set newYrefpix = `calc "($newYpixels+1)/2"`
echo "newYrefval = $newYrefval $yunit"
echo "newYrefpix = $newYrefpix"
echo "newYrange  = $newYrange $yunit"
echo "newYdelta  = $delt2 $yunit = last one found"
echo "newYpixels = $newYpixels"

echo ""
echo "====================="
echo "Regrid maps to mosaic"

foreach file ($files) #$files will need to be changed to 2 separate variables:
#$mfiles - list of mopra files to combine
#$ifiles - all other files in that folder   
#HOW TO KEEP ALL NON-MOPRA FILES AS SEPARATE EXTENSIONS???
    set filer = $file:r
    set filee = $file:e
    if ($starttype == 'fits') set filee = mir
    rm -rf $filer.regr
    regrid in=$filer.$filee out=$filer.regr axes=1,2 \
    desc=$newXrefval,$newXrefpix,$newXdelta,$newXpixels,$newYrefval,$newYrefpix,$newYdelta,$newYpixels
end
echo "Regrid Complete"

# something about the z axis is malformed as far as kvis is concerned
# z axis mostly ok with miriad, except when (eg) imstat tries to make rms spectrum plot, and z axis all screwed up

# imcomb will fail if $outs aren't properly contained in $filer.  
# then later, somehow the etaC scaling fails in PPlot for the mosaic, during the fits => mir conversion.  WTF???

echo ""
echo "============"
echo 'Combine maps'
#I think this is the step where I need to write to separate fits extensions
rm -rf $out1.mosaic$out2 $out1.mosaic$out2.fits
#set name = "$out1*$out2"
imcomb in="*.regr" out=$out1.mosaic$out2 options=relax
prthd in=$out1.mosaic$out2
fits in=$out1.mosaic$out2 out=$out1.mosaic$out2.fits op=xyout
echo ""
echo "Final mosaic map written to "$out1.mosaic$out2.fits

# Cleanup
foreach file ($files) #this looks unsafe...
    set filer = $file:r
    if ($delmir  == 'y') rm -rf $filer.mir
    if ($delregr == 'y') rm -rf $filer.regr
    if ($delcomb == 'y') rm -rf $out1.mosaic$out2
end
echo ""
echo 'All done'
end:
