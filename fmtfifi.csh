#!/bin/csh -f
# Requires molist, imlist to be premade
# external dependencies on imblrepair2.csh, strjoin.py

#set datadir = $1 # Where the files are located; '.' =  cwd
# Will have to set datadir manually either at command line or here
# Set padding factor by region - would take me longer to code 

#if (! -e maplist) then		# If no list given, merge all fits files
#    ls *.fits >! maplist
#endif
set ifiles = `sed '/^ *#/d;s/#.*//' imlist` #all other files in datadir not commented out
# Note: see http://www.grymoire.com/Unix/Sed.html for more info on Sed
set iredomir = 'n'
set deljybm = 'y'
set delblr = 'y'
set delconv = 'y'
set deljypx = 'y'
set delregr = 'y'
###OVERRIDE ABOVE PARAMETERS IF GIVEN ON COMMAND LINE
#***NOTE: this isn't working. All of my wat***
foreach par ( $argv[*] )
    set check=`echo $par | awk -F= '{print NF}'`
    if ( "$check" >= 2 ) set $par
end

goto start
start:
echo ""
echo "===================================="
echo "Start with fits to miriad conversion"

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
	echo "$file must be *.fits or miriad-style; skipping"
    endif

    rm -rf $ifiler.mir/mask #because miriad snows on 64-bit data -_-

    echo "Scaling and convolving $ifiler.mir..."
    echo "First convert Jy/pixel to Jy/beam:"
#    set chan = `gethd in="$ifiler.mir/channel"`
#    echo $chan
    set wave = `gethd in="$ifiler.mir/crval3"`
    set switch = `calc -i $wave - 100`
    if ($switch < 0) then
	set cd1 = `gethd in="$ifiler.mir/cdelt1"`
	set cd2 = `gethd in="$ifiler.mir/cdelt2"`
	set bmaj = `gethd in="$ifiler.mir/bmaj"`
	set bmin = `gethd in="$ifiler.mir/bmin"`
	set ratio = `calc "pi*abs(bmaj*bmin)/abs(2.7725887*cd1*cd2)"` #beam area / pixel area
	#can't make ln work for anything so just use 4*ln(2) = 2.7725887 (better precision than calc's pi
	rm -rf "$ifiler.jybm"
	rm -rf "$ifiler.blr"
	rm -rf "$ifiler.conv"
	rm -r "$ifiler.jypx"
	maths exp="<$ifiler.mir>*$ratio" mask="<$ifiler.mir>.lt.0." out="$ifiler.jybm"
    	puthd in="$ifiler.jybm/bunit" value='Jy/beam' 
	imblr in="$ifiler.jybm" out="$ifiler.blr"
	convol map="$ifiler.jybm" fwhm=10. pa=0.0 out="$ifiler.conv" options=cube | tee -a fmtfifi_log.txt
	maths exp="<$ifiler.conv>/$ratio"  mask="<$ifiler.mir>.le.0." out="$ifiler.jypx"
    	puthd in="$ifiler.jypx/bunit" value='Jy/pixel'
    	puthd in="$ifiler.jypx/bmaj" value=`calc "10./3600."` type=double
    	puthd in="$ifiler.jypx/bmin" value=`calc "10./3600."` type=double
	if ($deljybm == 'y') rm -rf "$ifiler.jybm"	
	if ($delblr == 'y') rm -rf "$ifiler.blr"
	if ($delconv == 'y') rm -rf "$ifiler.conv"
    else
	cp -rf $ifiler.mir $ifiler.jypx
    endif
    echo "Done."
end

echo ""
echo "===================================="
echo "Begin regridding and output..."
set tmplt = "fifi_cii_158um_FLUX.mir"
set etmplt = "fifi_cii_158um_ERROR.mir"

foreach file ($ifiles)
    set ifiler = $file:r #file name root
    if ($ifiler != $tmplt:r && $ifiler != $etmplt:r) then
	prthd in=$ifiler.jypx
	echo "Regridding $ifiler.jypx to $tmplt..."
	rm -r $ifiler.regr
	regrid in="$ifiler.jypx" out="$ifiler.regr" tin=$tmplt axes=1,2 | tee -a fmtfifi_log.txt
	echo "Done. Writing to fits..."
	rm -r $ifiler.regr.fits
	fits in="$ifiler.regr" out="$ifiler.regr.fits" op=xyout
	echo "Cleanup."
	if ($deljypx == 'y') rm -r "$ifiler.jypx"
	if ($delregr == 'y') rm -r "$ifiler.regr"
    endif
end
echo 'Complete'
end:
