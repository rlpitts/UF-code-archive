#!/bin/csh -f
# Requires molist, imlist to be premade
# external dependencies on imblrepair2.csh, strjoin.py

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
set delcirc = y		# > Cleanup options
set delmos = y		# |
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

    set mssg = `immask in="$ifiler.mir"`
    if ($mssg[5] < $mssg[8]) rm -rf $ifiler.mir/mask

    if ($iredomir == 'y') then
	if ($ifiler.mir =~ "*hpacs*" || $ifiler =~ "*hspire*") then
	    mv $ifiler.mir $ifiler.jnk
	    if ($ifiler.mir =~ "h*err*") then
	        maths exp="<$ifiler.jnk>" mask="<$ifiler.jnk>.lt.999" out="$ifiler.jnk2"
	    	./imblrepair2.sh "$ifiler.jnk2"
		maths exp="<$ifiler.jnk2>" mask="<$ifiler.jnk2>.ne.0" out="$ifiler.mir"
	    else
		#minm = `histo in="$ifiler.jnk" | grep Minimum | awk '{ if ($3 < 0) print $3; else print 0}'`
		maths exp="<$ifiler.jnk>" mask="<$ifiler.jnk>.ne.0" out="$ifiler.mir"
	    	./imblrepair2.sh "$ifiler.mir"
	    endif
	    rm -rf $ifiler.jnk*
	endif
    endif
    if (! -e $ifiler.mir) then
	echo "Error: $ifiler.mir was either deleted or not made"
	exit
    endif
end
#set fict = `ls -d R*/*.mir | wc -l`
#echo "Found $fict usable image files"
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
fits in="$t2" out="$t2.fits" op=xyout

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
    echo "X-, Ydelta: $Xdelta $Ydelta"
    echo "NAXIS1,-2: $Xpixels $Ypixels"
    echo "CRPIX1,-2: $Xrefpix $Yrefpix"
    set input = "$ifir.mir" #default regrid input, changed if hpacs
    if (-e $ifir.regr && $redoreg == "n") then
	echo "$ifir.regr exists, redo disabled"
    else
	rm -rf $ifir.regr
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
	    puthd in="$ifir.mir/bpa" value=$bpai
	    # by this point bpa should be same for all files in $datadir
	    #smooth or convol will take care of the rest
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
    set mcatch = `calc $mssg[5]/$mssg[8] | awk '{if ($1 < 0.25) print "bad"; else print "good"}'`
    # Miriad seems to have a problem with images that are too extensively masked
    set hcatch = `histo in="$ifir.regr" | grep All`
    #this will print 'All pixels are  0.0' (array length = 4) if region is all padding
    # for any variable, $# prints variable array length instead of variable itself
    if ($mcatch == 'bad' || $#hcatch == 4) then
	echo "$ifir.regr has too few unmasked pixels. Deleting to prevent further damage."
	rm -rf $ifir.regr
    else
	prthd in="$ifir.regr"
    endif
end
echo ""
echo "1st Regrid Complete." # Decluttering..."
echo "================================================================"
echo "mask, imcomb, imblr, convol, & regrid to $t2 for each band" | tee $datadir/imtile_log.txt
#don't add -a to tee here or else log file will grow with every run

#Apparently tsch doesn't support lists of lists ...t(-_-t)
# hopefully this will do as a workaround

set bmd = 37.
set bmrad=`calc "$bmd*pi/(180*3600)"`
set ddir = "$datadir/"
#Convolve down to 36" resolution, same as Mopra
# Initial params in $ifiles headers (even PACS & its special-snowflakyness)
#NOTE: convol expects BMAJ/BMIN and FWHM in arc sec but fits op=xyin expects BMAJ/BMIN in deg
# ^accounted for in fits header
set u = "err"
set bandarr = (8700atlasg19 5000hspire35 3500hspire24 2500hspire18 1600hpacsR15 0700hpacsB12 0240mipsga06 0220wisew412 0120wisew307) #0046wisew206 0034wisew106 0213msxs3E18 0147msxs3D18 0121msxs3C18 0083msxs3A18 0080iracI406 0058iracI306 0045iracI206 0036iracI106 #0080iracI412 0058iracI312 0045iracI212 0036iracI112 ) #file name bases & keys for setting input image lists
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
	set input = `ls $datadir/ATLASGAL_glon???.regr -d`
	set inputu = `ls $datadir/ATLASGAL_glon???_err.regr -d`
	set bname = "ATLASGAL"

    else if ($band == "5000hspire35") then
	set input = `ls $datadir/hspire*R.regr -d`
	set inputu = `ls $datadir/hspire*Rerr.regr -d`
	set bname = "hspireR"

    else if ($band == "3500hspire24") then
	set input = `ls $datadir/hspire*G.regr -d`
	set inputu = `ls $datadir/hspire*Gerr.regr -d`
	set bname = "hspireG"

    else if ($band == "2500hspire18") then
	set input = `ls $datadir/hspire*B.regr -d`
	set inputu = `ls $datadir/hspire*Berr.regr -d`
	set bname = "hspireB"

    else if ($band == "1600hpacsR15") then
	set input = `ls $datadir/hpacs*R.regr -d`
	set inputu = `ls $datadir/hpacs*Rerr.regr -d`
	set bname = "hpacsR"

    else if ($band == "0700hpacsB12") then
	set input = `ls $datadir/hpacs*B.regr -d`
	set inputu = `ls $datadir/hpacs*Berr.regr -d`
	set bname = "hpacsB"

    else if ($band == "0240mipsga06") then
	set input = `ls $datadir/MG????????_024.regr -d`
	set inputu = `ls $datadir/MG*err*.regr -d`
	set bname = "MIPSGAL"

    else if ($band == "0220wisew412") then
	set input = `ls $datadir/wise*w4-int*.regr -d`
	set inputu = `ls $datadir/wise*w4-err*.regr -d`
	set bname = "WISE band 4"

    else if ($band == "0120wisew307") then
	set input = `ls $datadir/wise*w3-int*.regr -d`
	set inputu = `ls $datadir/wise*w3-err*.regr -d`
	set bname = "WISE band 3"

#    else if ($band == "0046wisew206") then
#	set input = `ls $datadir/wise*w2-int*.regr -d`
#	set inputu = `ls $datadir/wise*w2-err*.regr -d`
#	set bname = "WISE band 2"

#    else if ($band == "0034wisew106") then
#	set input = `ls $datadir/wise*w1-int*.regr -d`
#	set inputu = `ls $datadir/wise*w1-err*.regr -d`
#	set bname = "WISE band 1"

#    else if ($band == "0213msxs3E18") then
#	set input = `ls $datadir/msx*E*.regr -d`
#	set bname = "MSX band E"

#    else if ($band == "0147msxs3D18") then
#	set input = `ls $datadir/msx*D*.regr -d`
#	set bname = "MSX band D"

#    else if ($band == "0121msxs3C18") then
#	set input = `ls $datadir/msx*C*.regr -d`
#	set bname = "MSX band C"

#    else if ($band == "0083msxs3A18") then
#	set input = `ls $datadir/msx*A*.regr -d`
#	set bname = "MSX band A"

#    else if ($band == "0080iracI406") then
#	set invcar46 = `ls $datadir/*VELACAR*I4*6.regr -d`
#	set inglm46 = `ls $datadir/GLM*I4*6.regr -d`
#	set input = "$invcar46 $inglm46"
#	set threshold = 0.000001
#	set bname = "IRAC band 4 0.6''/px"

#    else if ($band == "0058iracI306") then
#	set invcar36 = `ls $datadir/*VELACAR*I3*6.regr -d`
#	set inglm36 = `ls $datadir/GLM*I3*6.regr -d`
#	set input = "$invcar36 $inglm36"
#	set bname = "IRAC band 3 0.6''/px"

#    else if ($band == "0045iracI206") then
#	set invcar26 = `ls $datadir/*VELACAR*I2*6.regr -d`
#	set inglm26 = `ls $datadir/GLM*I2*6.regr -d`
#	set input = "$invcar26 $inglm26"
#	set bname = "IRAC band 2 0.6'/px"

#    else if ($band == "0036iracI106") then
#	set invcar16 = `ls $datadir/*VELACAR*I1*6.regr -d`
#	set inglm16 = `ls $datadir/GLM*I1*6.regr -d`
#	set input = "$invcar16 $inglm16"
#	set bname = "IRAC band 1 0.6''/px"

    else
	echo "Unmatched band $band - check your input list"
	exit 1
    endif
    #set nf = `echo $input | wc -w`
    #set nfu = `echo $inputu | wc -w`

    if ($#input == 0 && $#inputu == 0) then
    	echo "$datadir contains no usable $bname images"

    else if ($#input == 0 && $#inputu != 0) then
        echo "ERROR: Cannot compute uncertainty maps for $bname without image maps."
	exit 1

    else
	echo "$#input Input file(s) found"
	set out0 = "$datadir/$band-$datadir.mos"
	set out1 = "$datadir/$band-$datadir.blr"
	set out2 = "$datadir/$band-$datadir.conv"
	set out3 = "$datadir/$band-$datadir.mr" #master/final file
    	rm -rf $out0
	rm -rf $out1
    	rm -rf $out2
    	#rm -rf $out3
    	rm -rf $out3.fits

    	if ($#input > 1) then
	    echo "Now masking for imcomb"
	    set minput = ""
	    set max = ""
	    foreach f ($input)
		set outf = $f:r
		rm -rf $outf.mskd
		#IMPORTANT: mask is binary Boolean - where mask = True, data are NOT masked
		if ($band == "8700atlasg19") then
		    maths exp="<$f>" mask="abs(<$f>).ne.0" out="$outf.mskd"
		else
		    maths exp="<$f>" mask="abs(<$f>).gt.0" out="$outf.mskd"
		endif
		set minput = "$minput $outf.mskd"
		set newmax = `histo in="$f" | grep Maximum | awk '{ print $3 }' `
		set max = `echo "$newmax $max" | awk '{if ($1 > $2) print $1; else print $2}'`
		#^blank is replaced by newmax too.
	    end #by here, input --> minput via previous loop, hopefully

	    echo "combining the following $#input files: $minput" | tee -a $datadir/imtile_log.txt
    	    imcomb in="$minput" out="$out0" #do NOT delete $out0 yet!
	    set combmax = `histo in="$out0" | grep Maximum | awk '{ print $3 }' `
	    set checkmax = `echo "$combmax $max" | awk '{if ($1 >= 1.5*$2 || $1 <= 0.5*$2) print "bad"; else print "good"}'`
	    if ($checkmax == "bad") then
		echo "Error: imcomb failed to normalize overlap"
		exit 1
	    endif

	    if ($band == "8700atlasg19") then
		imblr in="$out0" out="$out1" value=0
		if (! -e $out1) then set out1 = $out0
	    else
		./imblrepair2.sh $out0 ofile=$out1 sweep=n axis=xy
	    endif
    	else if ($#input == 1) then
	    set mssg = `immask in="$input"`
	    if ($mssg[5] < $mssg[8]) then
		if ($band == "8700atlasg19") then
		    imblr in="$input" out="$out1"
		    if (! -e $out1) then set out1 = $input
	    	else
		    ./imblrepair2.sh $input ofile=$out1 sweep=$delregr axis=xy
	        endif
	    else
		set out1 = $input
	    endif		
    	endif

	echo "now convolving input files in band $band" | tee -a $datadir/imtile_log.txt
    	convol map="$out1" fwhm=$bmd pa=0.0 out="$out2" options=final | tee -a $datadir/imtile_log.txt
        #regrid in="$out2" out="$out3" tin="$t2" | tee -a $datadir/imtile_log.txt
        fits in="$out2" out="$out2.fits" op=xyout | tee -a $datadir/imtile_log.txt
	if ($delblr == 'y') rm -r $out1
	if ($delconv == 'y') rm -r $out2
        echo "Now writing $out3.fits"
	echo "regrid, '$out2.fits', '$t2.fits', '$out3.fits'" | idl
	#rm -r $out3
	#python confetti_sweep.py $out3.fits -t zero

	if ('$#inputu' != '0') then
	    set out00u = "$datadir/$band-$u-$datadir.sqr"
	    set out0u = "$datadir/$band-$u-$datadir.mos"
	    set out1u = "$datadir/$band-$u-$datadir.blr"
	    set out2u = "$datadir/$band-$u-$datadir.conv"
	    set out3u = "$datadir/$band-$u-$datadir.mr" #master/final file
    	    rm -rf $out00u
    	    rm -rf $out0u
    	    rm -rf $out1u
    	    rm -rf $out2u
    	    #rm -rf $out3u
    	    rm -rf $out3u.fits

    	    if ($#inputu > 1) then
		set minputu = ""
		foreach f ($inputu)
		    set outfu = $f:r
		    #set datf = $input[1]
		    rm -rf $outfu.mskd
		    set maxm = `histo in="$f" | grep Maximum | awk '{ print $3 }' ` 
		    set thresholdu = `calc $maxm/9 | awk '{if ($1 < 1) print 9; else print $1}'`
		    echo $thresholdu
		    #image maps are added, not multiplied; errors add in quadrature
		    #if $thresholdu is not defined, $#thresholdu will be empty & file should be deleted
		    if ($#thresholdu != 0) then
			maths exp="<$f>**2" mask="<$f>.lt.$thresholdu" out="$outfu.mskd" | tee -a $datadir/imtile_log.txt 
			set minputu = "$minputu $outfu.mskd"
		    else
			rm -rf $f 
		    endif
		end
		#set summ = `python strjoin.py $minputu -c '>+<' -p '<' -s '>'`
		echo "Adding error maps in quadrature..." | tee -a $datadir/imtile_log.txt
		imcomb in="$minputu" out="$out00u" options=nonormalise | tee -a $datadir/imtile_log.txt
		#this is where I use the image map to mask the error map
		maths exp="sqrt(<$out00u>)" mask="<$out0>.ne.0" out="$out0u" | tee -a $datadir/imtile_log.txt
		rm -rf $out00u

		set mssg = `immask in="$out0u"`
		if ($mssg[5] < $mssg[8]) then
		    set maxm = `histo in=$out0u | grep Maximum | awk '{ print $3 }' ` 
		    set thresholdu = `calc $maxm/9 | awk '{if ($1 < 1) print 9; else print $1}'`
		    imblr in="$out0u" out="$out1u" value=$thresholdu
		else
		    set out1u = $out0u
		endif

    	    else if ($#inputu == 1) then
		set mssg = `immask in="$inputu"`
		if ($mssg[5] < $mssg[8]) then
		    set maxm = `histo in="$inputu" | grep Maximum | awk '{ print $3 }' ` 
		    set thresholdu = `calc $maxm/9 | awk '{if ($1 < 1) print 9; else print $1}'`
		    imblr in="$inputu" out="$out1u" value=$thresholdu
		else
		    set out1u = $inputu
		endif
    	    endif

    	    convol map="$out1u" fwhm=$bmd pa=0.0 out="$out2u" options=final | tee -a $datadir/imtile_log.txt
	    set maxm = `histo in="$out2u" | grep Maximum | awk '{ print $3 }'`
	    fits in="$out2u" out="$out2u.fits" op=xyout | tee -a $datadir/imtile_log.txt
            #regrid in="$out2u" out="$out3u" tin="$t2" | tee -a $datadir/imtile_log.txt
	    if ($delregr == 'y') rm -r $out0u
	    if ($delblr == 'y') rm -r $out1u
	    if ($delconv == 'y') rm -r $out2u
            echo "regrid, '$out2u.fits', '$t2.fits', '$out3u.fits', MISSING=$maxm" | idl
            echo "Now writing $out3u.fits"
	    #rm -r $out3u
	    #python confetti_sweep.py $out3u.fits -t zero
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

if ($delmos == 'y') rm -rf $ddir*.mos
if ($delregr == 'y') rm -rf $ddir*.regr
if ($delconv == 'y') rm -rf $ddir*.conv.fits
if ($delmskd == 'y') rm -rf $ddir*.mskd
if ($delblr == 'y') rm -rf $ddir*.blr
echo ""
echo 'All done'
end:
