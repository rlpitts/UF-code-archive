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
# 2015Oct19: Begin editing for use on CHaMP data -> changing variables -RLP

set datadir = $1 # Where the files are located; '.' =  cwd
# Will have to set datadir manually either at command line or here
# Set padding factor by region - would take me longer to code 

#if (! -e maplist) then		# If no list given, merge all fits files
#    ls *.fits >! maplist
#endif
set mfiles = `cat ./molist` #list of mopra files
# Note: cat command followed by file name prints contents of that file
set ifiles = `sed '/^ *#/d;s/#.*//' $datadir/imlist` #all other files in datadir not commented out
# Note: see http://www.grymoire.com/Unix/Sed.html for more info on Sed
set redomir = n
set iredomir = y
set redoreg = y

set delmir = n		#\
set delregr = y         # |
set deltrim = y		# |
set delmos = y		# > Cleanup options
set delcirc = y		# |
set delconv = y		# |
set delblr = y		# |
set delmskd = y		#/

#Notes 1:
# set assigns variables local to this shell
# set env assigns variables for this shell and all subshells

###OVERRIDE ABOVE PARAMETERS IF GIVEN ON COMMAND LINE
#***NOTE: this isn't working. All of my wat***
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
# **ALL FILENAMES ENTERED IN in= OR out= MUST BE IN ""
#   otherwise, miriad may misinterpret non-alphanumeric characters

foreach file ($mfiles)
    echo "now loading $file"
    set filee = $file:e #extension
    set filer = $file:r #file name root
    if ($filee == 'fits' || -d $file) then
    #if $file extension is fits OR $file is a directory
    	if (-e $filer.fits && ! -d $file && ! -e $filer.mir) then
	#if $file of given name exists AND is not a directory or miriad file
	    fits in="$filer.fits" out="$filer.mir" op=xyin #miriad fits conversion
	else if (-e $filer.fits && -e $filer.mir) then
	    if ($redomir == 'y') then
		echo "overwriting $filer.mir"
		rm -rf $filer.mir #remove existing version of $file if extant
		fits in="$filer.fits" out="$filer.mir" op=xyin
	    else
		echo "$filer.mir already exists, overwrite disabled"
	    endif
	else if (-d $file || $filee == 'mir') then
	    echo "$file is an existing directory, assuming miriad-style"
	    set delmir=n	# DO NOT DELETE!!!
	endif
    else
	echo "$file must be *.fits or miriad-style, skipping"
    endif
end

foreach file ($ifiles)
    set ifilee = $file:e #extension
    set ifiler = $file:r #file name root
    if ($ifilee == 'fits' || -d $file) then
        echo "now loading $file" 
	if (-e $ifiler.fits && ! -e $ifiler.mir && ! -d $file) then
	    fits in="$ifiler.fits" out="$ifiler.mir" op=xyin #miriad fits conversion
	    #op=xyin converts fits to mir; op=xyout converts mir to fits
	    echo "writing $ifiler.mir"
	else if (-e $ifiler.fits && -e $ifiler.mir) then
	    #contrary to last loop, default here is to replace .mir files
	    if ($iredomir == 'n') then
		echo "$ifiler.mir already exists, overwrite disabled"
	    else
		echo "overwriting $ifiler.mir"
		rm -rf $ifiler.mir #remove existing version of $file if extant
		fits in="$ifiler.fits" out="$ifiler.mir" op=xyin
	    endif
	else if (-d $file || $ifilee == 'mir') then
	    echo "$file is an existing directory, assuming miriad-style"
	    set delmir=n	# DO NOT DELETE!!!
	endif
    else
	echo "$file must be *.fits or miriad-style, skipping"
    endif
end
set fict = `ls -d R*/*.mir | wc -l`
echo "Found $fict usable image files"
echo ""
echo "==================================="
echo "Set up template area for Regridding"

foreach file ($mfiles) #mfiles should actually only have 1 file per folder, maybe 2
    if ($file:h == $datadir) then
	set filer = $file:r #get root of $file
    	set filet = $file:t #get tail of $file (file name + extension)
	if ($cwd:t == $datadir) then
	    set filer = $filet:r #get root of tail
	endif
	set t1 = $filer.mir #**file to regrid to make new template**
    	set filedim = `gethd in=$filer.mir/naxis` #get image dims from header
    	echo ""
    	echo "$filer"

    	#Accumulate X axis descriptors (I know all these terms)
    	set ctyp1 = `gethd in="$filer.mir/ctype1"`
    	set npix1 = `gethd in="$filer.mir/naxis1"`
    	set delt1 = `gethd in="$filer.mir/cdelt1"`
    	set rpix1 = `gethd in="$filer.mir/crpix1"`
    	set rval1 = `gethd in="$filer.mir/crval1"`
    	#echo $ctyp1 $npix1 $delt1 $rpix1 $rval1
	# Must pad lmin & lmax slightly or rotations will cut corners
    	set lmax = `calc "($rval1)+($delt1)*(($rpix1)-($npix1))"`
    	set lmin = `calc "($rval1)+($delt1)*(($rpix1)-(1))"`
    	echo "$lmax $lmin"

    	#Accumulate Y axis descriptors
    	set ctyp2 = `gethd in="$filer.mir/ctype2"`
    	set npix2 = `gethd in="$filer.mir/naxis2"`
    	set delt2 = `gethd in="$filer.mir/cdelt2"`
    	set rpix2 = `gethd in="$filer.mir/crpix2"`
    	set rval2 = `gethd in="$filer.mir/crval2"`
    	#echo $ctyp2 $npix2 $delt2 $rpix2 $rval2
    	set bmax = `calc "($rval2)-($delt2)*(($rpix2)-($npix2))"`
    	set bmin = `calc "($rval2)-($delt2)*(($rpix2)-(1))"`
    	echo "$bmax $bmin"
    else
	echo $file:h " not in $datadir"
    endif
end

#Use last map to define coord types - should not be any velocity axes by now
set cbase1 = `echo "$ctyp1" | cut -d - -f 1 -`
set cbase2 = `echo "$ctyp2" | cut -d - -f 1 -`
echo "ctypes: $cbase1, $cbase2"
#if ($filedim > 2) set cbase3 = `echo "$ctyp3" | cut -d - -f 1 -`
if ($cbase1 == 'RA' || $cbase1 == 'DEC' || $cbase1 == 'GLON' || $cbase1 == 'GLAT') then
    set xunit = radians
else if ($cbase1 == 'VELO') then
    set xunit = km/sec
endif
if ($cbase2 == 'RA' || $cbase2 == 'DEC' || $cbase2 == 'GLON' || $cbase2 == 'GLAT') then
    set yunit = radians
else if ($cbase2 == 'VELO') then
    set yunit = km/sec
endif
# Insert any other unit definitions here
echo "Maximum x coord = $lmax $xunit"
echo "Minimum x coord = $lmin $xunit"
echo "Maximum y coord = $bmax $yunit"
echo "Minimum y coord = $bmin $yunit"
echo "Minimum x resolution = $delt1 $xunit"
echo "Minimum y resolution = $delt2 $yunit"
endif

echo ""
echo "============================"
echo "Make template for 2nd regrid"
# MUST CHANGE THIS FROM MAXIMUM EXTENT OF ALL FILES TO
# MAX AREA OF ALL MOPRA FILES IN REGION (w/ slight padding)
# If filename includes 12co, 13co, or c18o, use mom0  
# If filename includes hco, use mom0 files containing 'main'
# If region files are disjoint, draw perimeter around all contained patches
# For multiple velo integrals, take max extent of all combined 

#X axis calculations: if angles, convert radians to degrees for regrid
set pixsiz = 12. #pixel size for 36" resolution
# Nyquist sampling with 2 pix per resolution element
# Some people like 3 or 4. Peter uses 3 so I'll stick with that.
set pad = `awk '/'$datadir'/ {print $4}' rlist-eq2galPAs.tbl`
# ^range-padding factor = sqrt(1+|sin(2a)|)
if ($xunit == 'radians') then
    set newXrefval = `calc "(($lmax)+($lmin))*(180/pi)/2"`
    set newXrange  = `calc "abs(($lmax)-($lmin))*(180/pi)*$pad"`
    #set newXdelta  = `calc "($delt1)*(180/pi)"`
    set xunit = degrees
else
    set newXrefval = `calc "(($lmax)+($lmin))/2"`	# This will be *MEAN* pixel!!! (don't median-combine)
    set newXrange  = `calc "abs(($lmax)-($lmin))*$pad"`	# = (pixels-1)*delta
    # pad newXrange to keep rotations from cutting corners (# of pixels not affected)
    # 4% is just a starting value
    #set newXdelta  = `calc "($delt1)"` 		# All the same?
endif
set newXdelta  = `calc "-($pixsiz)/3600"`
if ($yunit == 'radians') then
    # cosine(declination/latitude) correction
    set newXpixels = `calc -i "($newXrange/abs($newXdelta)*cos($rval2))+(1)"`
else
    set newXpixels = `calc -i "($newXrange/abs($newXdelta))+(1)"`
endif
set newXrefpix = `calc "($newXpixels+1)/2"`		# *Mean* pixel number (can be fractional)
echo "newXpixels = $newXpixels"
echo "newXrefval = $newXrefval $xunit"
echo "newXrefpix = $newXrefpix"
echo "newXrange  = $newXrange $xunit"
echo "newXdelta fixed at $newXdelta degrees"

#Y axis calculations
if ($yunit == 'radians') then
    set newYrefval = `calc "(($bmax)+($bmin))*(180/pi)/2"`
    set newYrange  = `calc "abs(($bmax)-($bmin))*(180/pi)*$pad"`
    set yunit = degrees
else
    set newYrefval = `calc "(($bmax)+($bmin))/2"`
    set newYrange  = `calc "abs(($bmax)-($bmin))*$pad"`
    # pad newYrange to keep rotations from cutting corners
endif
set newYdelta  = `calc "($pixsiz)/3600"`
set newYpixels = `calc -i "($newYrange/abs($newYdelta))+(1)"`
set newYrefpix = `calc "($newYpixels+1)/2"`
echo "newYpixels = $newYpixels"
echo "newYrefval = $newYrefval $yunit"
echo "newYrefpix = $newYrefpix"
echo "newYrange  = $newYrange $yunit"
echo "newYdelta fixed at $newYdelta degrees"

echo "Making template..."
set t2 = $t1:r.tin
if (! -e $t2) then
    regrid in="$t1" out="$t2" axes=1,2 \
    desc=$newXrefval,$newXrefpix,$newXdelta,$newXpixels,$newYrefval,$newYrefpix,$newYdelta,$newYpixels
else
    echo "$t2 exists"
endif
prthd in="$t2"

echo ""
echo "==============================="
echo "First Regrid: set-up for imcomb"

foreach file ($ifiles)
#Parameters shared for both regriddings: X/Yrange, X/Yrefval (crvals), x/yunit
#Parameters to be computed for each image/list of images: X/Ydelta, X/Ypixels, X/Yrefpix (crpix)
#crpix must be computed based on native resolution of each image, i.e. using cdelts
    set ifir = $file:r
    echo "Obtaining axis descriptors for $file..."
    #x-axis
    set cd1 = `gethd in="$ifir.mir/cdelt1"` #in radians - keep native resolution
    set Xdelta = `calc "($cd1)*(180/pi)"` #convert to degrees
    set Xpixels = `calc -i "($newXrange/abs($Xdelta)*cos($newYrefval*pi/180))+(1)"` #lat/lon correction (calc trig fxns take angles in rads)
    set Xrefpix = `calc "($Xpixels+1)/2"` # mean pixel
    #y-axis
    set cd2 = `gethd in="$ifir.mir/cdelt2"` #in radians - keep native resolution
    set Ydelta = `calc "($cd2)*(180/pi)"` #convert back to degrees
    set Ypixels = `calc -i "($newYrange/abs($Ydelta))+(1)"` #lat/lon correction
    set Yrefpix = `calc "($Ypixels+1)/2"` # mean pixel
    #echo "X-, Ydelta: $Xdelta $Ydelta"
    #echo "NAXIS1,-2: $Xpixels $Ypixels"
    #echo "CRVAL1,-2: $Xrefpix $Yrefpix"
    set input = "$ifir.mir" #default regrid input, changed if hpacs
    if (-e $ifir.regr && $redoreg == "n") then
	echo "$ifir.regr exists, redo disabled"
    else
	rm -rf $ifir.regr
	imblr in="$ifir.mir" out="$input" value=0 #convol wipes image if it contains nans
	#native resolution differs by band
	#cdelts should be the same for all bands of same instrument - same detector
	#PACS is a special snowflake that needs to be circularized
	# so it can be regridded & not f*** up everything later (beam dims & PA in hdr)*
	#*except they won't all have the same PA or BPA in hdr because they're spread across several regions
	# --> some were found in another region first & assigned a different PA & BPA accordingly.
	# BPA additionally depends on the scan speed. Good effing god.
    	if ($ifir =~ "*hpacs*") then
	    set bpai = `grep 'ADDBEAM.PY' $ifir.mir/history | grep $datadir | awk '{print substr($4,1)}'`
	    # for the beams I rounded up to the nearest integer larger than both bmaj & bmin (undefined if smaller)
	    puthd in="$ifir.mir/BPA" value=$bpai
	    # by this point bpa should be same for all files in $datadir
            rm -rf $ifir.blr
	    imblr in="$ifir.mir" out="$ifir.blr" value=0 #convol wipes image if it contains nans
	    if ($ifir =~ "*hpacs*B*") then
		set b = 13
	    else if ($ifir =~ "*hpacs*R*") then
		set b = 16
	    endif
	    rm -rf $ifir.circ
	    rm -rf $ifir.sharp
	    set input = $ifir.trim #<-- output of this subroutine, input of regrid
	    rm -rf $ifir.trim
	    convol map="$ifir.blr" fwhm=$b,$b out="$ifir.circ" options=final #circ for circularized beam
 	    imsharp in="$ifir.circ" out="$ifir.sharp" #strip badly-convolved edges (this and next line)
	    maths exp="<$ifir.circ>" mask="<$ifir.sharp>.gt.0" out="$input" 
	    # maths requires < > around all filenames containing +, -, /, or other operator-like chars
	    # .gt & .lt are FORTRAN syntax for > & < respectively 
        endif
    	regrid in="$input" out="$ifir.regr" axes=1,2 \
    	desc=$newXrefval,$Xrefpix,$Xdelta,$Xpixels,$newYrefval,$Yrefpix,$Ydelta,$Ypixels \
	rotate=0.0
    endif
    if (! -e $ifir.regr) then
	echo "Oh dear. Existence failure. Curse you $ifir.regr /dies"
	exit 1
    endif
    set mssg = `immask in="$ifir.regr"`
    set catch = `calc $mssg[5]/$mssg[8] | awk '{if ($1 < 0.3) print "bad"; else print "good"}'`
    # Miriad seems to have a problem with images that are >70% masked
    if ($catch == 'bad') then
	echo "$ifir.regr has too few unmasked pixels. Deleting to prevent further damage."
	rm -rf $ifir.regr
    else
	prthd in="$ifir.regr"
    endif
end
echo ""
echo "1st Regrid Complete. Decluttering..."
set ddir = "$datadir/"
rm -rf $ddir*.circ
rm -rf $ddir*.sharp
rm -rf $ddir*.trim
rm -rf $ddir*.blr #this ext is used again later for the 2nd round of imblr
echo "...done."
#exit 0 #needed for diagnosis
echo "================================================================"
echo "mask, imcomb, imblr, convol, & regrid to $t2 for each band in Region $datadir"
#Apparently tsch doesn't support lists of lists ...t(-_-t)
# hopefully this will do as a workaround

set bmd = 36
#Convolve down to 36" resolution, same as Mopra
# Initial params in $ifiles headers (even PACS & its special-snowflakyness)
#NOTE: convol expects BMAJ/BMIN and FWHM in arc sec but fits op=xyin expects BMAJ/BMIN in deg
# ^accounted for in fits header
set u = err
set bandarr = (8700atlasg19 5000hspire35 3500hspire24 2500hspire18 1700hpacsR15 0750hpacsB12 0240mipsga06 0220wisew412 0120wisew307 0046wisew206 0034wisew106 0213msxs3E18 0147msxs3D18 0121msxs3C18 0083msxs3A18 0080iracI412 0058iracI312 0045iracI212 0036iracI112 0080iracI406 0058iracI306 0045iracI206 0036iracI106) #file name bases & keys for setting input image lists
#Key:
# chars 1-4 = wavelength*10
# chars 5-10 = instrument (abbr)
# chars 11-12 = native resolution ['']*10

#Here's where shit gets complicated:
# imcomb requires a mask, but convol wipes the image if the mask is too extensive, so
# I have to mask the image so that imcomb knows not to take average of data & no data,
# then replace the masked values with 0s so convol doesn't choke
#For PACS, I have to do this twice since it has to be convolved twice. FTGE.

foreach band ($bandarr)
    set inputu = "" #make sure no residual inputu from previous iterations
    set input = ""
    set bname = ""
    echo "Making files of base $band..." #this is SUCH a pain. Holy crap.
    # I could do this in a fraction of the lines if I could work out how to invoke Python from here.
    if ($band == "8700atlasg19") then
	set input = `find $datadir/ATLASGAL*.regr -maxdepth 0`
	set bname = "ATLASGAL"

    else if ($band == "5000hspire35") then
	set input = `find $datadir/hspire*R.regr -maxdepth 0`
	set inputu = `find $datadir/hspire*Rerr.regr -maxdepth 0`
	set threshold = 0.00001
	set bname = "hspireR"

    else if ($band == "3500hspire24") then
	set input = `find $datadir/hspire*G.regr -maxdepth 0`
	set inputu = `find $datadir/hspire*Gerr.regr -maxdepth 0`
	set threshold = 0.00001
	set bname = "hspireG"

    else if ($band == "2500hspire18") then
	set input = `find $datadir/hspire*B.regr -maxdepth 0`
	set inputu = `find $datadir/hspire*Berr.regr -maxdepth 0`
	set threshold = 0.00001
	set bname = "hspireB"

    else if ($band == "1700hpacsR15") then
	set input = `find $datadir/hpacs*R.regr -maxdepth 0`
	set inputu = `find $datadir/hpacs*Rerr.regr -maxdepth 0`
	set threshold = 0.00001
	set bname = "hpacsR"

    else if ($band == "0750hpacsB12") then
	set input = `find $datadir/hpacs*B.regr -maxdepth 0`
	set inputu = `find $datadir/hpacs*Berr.regr -maxdepth 0`
	set threshold = 0.00001
	set bname = "hpacsB"

    else if ($band == "0240mipsga06") then
	set input = `find $datadir/MG????????_024.regr -maxdepth 0`
	set inputu = `find $datadir/MG*std*.regr -maxdepth 0`
	set threshold = 0.0001 #use 0.000001 for stdev image
	set bname = "MIPSGAL"

    else if ($band == "0220wisew412") then
	set input = `find $datadir/wise*w4-int*.regr -maxdepth 0`
	set inputu = `find $datadir/wise*w4-unc*.regr -maxdepth 0`
	set threshold = 0.0001
	set bname = "WISE band 4"

    else if ($band == "0120wisew307") then
	set input = `find $datadir/wise*w3-int*.regr -maxdepth 0`
	set inputu = `find $datadir/wise*w3-unc*.regr -maxdepth 0`
	set threshold = 0.00001
	set bname = "WISE band 3"

    else if ($band == "0046wisew206") then
	set input = `find $datadir/wise*w2-int*.regr -maxdepth 0`
	set inputu = `find $datadir/wise*w2-unc*.regr -maxdepth 0`
	set threshold = 0.000001
	set bname = "WISE band 2"

    else if ($band == "0034wisew106") then
	set input = `find $datadir/wise*w1-int*.regr -maxdepth 0`
	set inputu = `find $datadir/wise*w1-unc*.regr -maxdepth 0`
	set threshold = 0.000001
	set bname = "WISE band 1"

    else if ($band == "0213msxs3E18") then
	set input = `find $datadir/msx*E*.regr -maxdepth 0`
	set threshold = 0.00001
	set bname = "MSX band E"
    else if ($band == "0147msxs3D18") then
	set input = `find $datadir/msx*D*.regr -maxdepth 0`
	set threshold = 0.000001
	set bname = "MSX band D"

    else if ($band == "0121msxs3C18") then
	set input = `find $datadir/msx*C*.regr -maxdepth 0`
	set threshold = 0.000001
	set bname = "MSX band C"

    else if ($band == "0083msxs3A18") then
	set input = `find $datadir/msx*A*.regr -maxdepth 0`
	set threshold = 0.000001
	set bname = "MSX band A"

    else if ($band == "0080iracI412") then
	set invcar412 = `find $datadir/*VELACAR*I4*2.regr -maxdepth 0`
	set inglm412 = `find $datadir/GLM*I4*2.regr -maxdepth 0`
	set input = "$invcar412 $inglm412"
	set threshold = 0.000001
	set bname = "IRAC band 4 1.2''/px"

    else if ($band == "0058iracI312") then
	set invcar312 = `find $datadir/*VELACAR*I3*2.regr -maxdepth 0`
	set inglm312 = `find $datadir/GLM*I3*2.regr -maxdepth 0`
	set input = "$invcar312 $inglm312"
	set threshold = 0.000001
	set bname = "IRAC band 3 1.2''/px"

    else if ($band == "0045iracI212") then
	set invcar212 = `find $datadir/*VELACAR*I2*2.regr -maxdepth 0`
	set inglm212 = `find $datadir/GLM*I2*2.regr -maxdepth 0`
	set input = "$invcar212 $inglm212"
	set threshold = 0.000001
	set bname = "IRAC band 2 1.2''/px"

    else if ($band == "0036iracI112") then
	set invcar112 = `find $datadir/*VELACAR*I1*2.regr -maxdepth 0`
	set inglm112 = `find $datadir/GLM*I1*2.regr -maxdepth 0`
	set input = "$invcar112 $inglm112"
	set threshold = 0.000001
	set bname = "IRAC band 1 1.2''/px"

    else if ($band == "0080iracI406") then
	set invcar46 = `find $datadir/*VELACAR*I4*6.regr -maxdepth 0`
	set inglm46 = `find $datadir/GLM*I4*6.regr -maxdepth 0`
	set input = "$invcar46 $inglm46"
	set threshold = 0.000001
	set bname = "IRAC band 4 0.6''/px"

    else if ($band == "0058iracI306") then
	set invcar36 = `find $datadir/*VELACAR*I3*6.regr -maxdepth 0`
	set inglm36 = `find $datadir/GLM*I3*6.regr -maxdepth 0`
	set input = "$invcar36 $inglm36"
	set threshold = 0.000001
	set bname = "IRAC band 3 0.6''/px"

    else if ($band == "0045iracI206") then
	set invcar26 = `find $datadir/*VELACAR*I2*6.regr -maxdepth 0`
	set inglm26 = `find $datadir/GLM*I2*6.regr -maxdepth 0`
	set input = "$invcar26 $inglm26"
	set threshold = 0.000001
	set bname = "IRAC band 2 0.6'/px"

    else if ($band == "0036iracI106") then
	set invcar16 = `find $datadir/*VELACAR*I1*6.regr -maxdepth 0`
	set inglm16 = `find $datadir/GLM*I1*6.regr -maxdepth 0`
	set input = "$invcar16 $inglm16"
	set threshold = 0.000001
	set bname = "IRAC band 1 0.6''/px"

    else
	echo "Unmatched band $band - check your input list"
	exit 1
    endif
    set nf = `echo $input | wc -w`
    set nfu = `echo $inputu | wc -w`
    if ($nf == 0 && $nfu == 0) then
    	echo "$datadir contains no usable $bname images"
    else if ($nf == 0 && $nfu != 0) then
        echo "ERROR: Uncertainty maps for $bname should not exist without image maps."
	exit 1
    else
	echo "Input file(s) found"
	if ($band == 0080iracI412) echo $input
	set out0 = "$datadir/$band-$datadir.mos"
	set out1 = "$datadir/$band-$datadir.blr"
	set out2 = "$datadir/$band-$datadir.conv"
	set out3 = "$datadir/$band-$datadir.mr" #master/final file
    	rm -rf $out0
	rm -rf $out1
    	rm -rf $out2
    	rm -rf $out3
    	rm -rf $out3.fits
    	if ($nf > 1) then
	    echo "Stand by: combining $nf input files"
	    if ($band == "8700atlasg19" || $band =~ "*msxs3*") then
		echo "found ATLASGAL or MSX files --> no mask needed"
	    	imcomb in="$input" out="$out0"
		imblr in="$out0" out="$out1" value=0
    	        convol map="$out1" fwhm=$bmd out="$out2" options=final
	    else
		echo "Now masking for imcomb"
		set minput = ""
		foreach f ($input)
		    set outf = $f:r
		    rm -rf $outf.mskd
		    #IMPORTANT: mask is binary Boolean - where mask = True, data are NOT masked
		    maths exp="<$f>" mask="abs(<$f>).gt.$threshold" out="$outf.mskd"
		    set minput = "$minput $outf.mskd"
		end #by here, input --> minput via previous loop, hopefully
		echo "now convolving $nf input files"
    		imcomb in="$minput" out="$out0"
		imblr in="$out0" out="$out1" value=0
    		convol map="$out1" fwhm=$bmd out="$out2" options=final
	    endif
    	else if ($nf == 1) then
	    echo 'now convolving 1 input file'
	    imblr in="$input" out="$out1" value=0
    	    convol map="$out1" fwhm=$bmd out="$out2" options=final
    	endif
        regrid in="$out2" out="$out3" tin="$t2"
        fits in="$out3" out="$out3.fits" op=xyout
        echo "Now writing $out3.fits"
	if ($nfu != 0) then
	    set out0u = "$datadir/$band-$u-$datadir.mos"
	    set out1u = "$datadir/$band-$u-$datadir.blr"
	    set out2u = "$datadir/$band-$u-$datadir.conv"
	    set out3u = "$datadir/$band-$u-$datadir.mr" #master/final file
    	    rm -rf $out0u
    	    rm -rf $out1u
    	    rm -rf $out2u
    	    rm -rf $out3u
    	    rm -rf $out3u.fits	
    	    if ($nfu > 1) then #no ATLASGAL images here
		set minputu = ""
		foreach f ($inputu)
		    set outfu = $f:r
		    rm -rf $outfu.mskd
		    set thresholdu = 0.000001 #just leave this the same for all
		    maths exp="<$f>" mask="abs(<$f>).gt.$thresholdu" out="$outfu.mskd"
		    set minputu = "$minputu $outfu.mskd"
		end #by here, inputu --> minputu via previous loop, hopefully
    	        imcomb in="$minputu" out="$out0u" options=relax
		imblr in="$out0u" out="$out1u" value=0
    	        convol map="$out1u" fwhm=$bmd out="$out2u" options=final
    	    else if ($nfu == 1) then
		imblr in="$inputu" out="$out1u" value=0
    	        convol map="$inputu" fwhm=$bmd out="$out2u" options=final
    	    endif
            regrid in="$out2u" out="$out3u" tin="$t2"
            fits in="$out3u" out="$out3u.fits" op=xyout
            echo "Now writing $out3u.fits"
	endif
    endif
end

echo ""
echo "======="
echo "Cleanup"

# Cleanup
foreach file ($ifiles)
    set filer = $file:r
    if ($delmir  == 'y') rm -rf $filer.mir
end

# .blr part takes awhile - make option to detect & keep existing .blr files
if ($delmos == 'y') rm -rf $ddir*.mos
if ($delregr == 'y') rm -rf $ddir*.regr
if ($delconv == 'y') rm -rf $ddir*.conv
if ($delmskd == 'y') rm -rf $ddir*.mskd
if ($delblr == 'y') rm -rf $ddir*.blr
echo ""
echo 'All done'
end:
