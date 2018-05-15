pro regrid, oldim, refim, newim, INTERP=interp, MISSING=missing

;+
; NAME:
;	regrid
; PURPOSE:
;	call HASTROM.PRO to up/down-sample, shift,
;	& rotate an image to match astrometry in
;	reference image header & write output to
;	a new file
; CALLING:
;	regrid, oldim, refim, newim, INTERP=2, MISSING=0
; INPUTS:
;	oldim = input fits file to be regridded
;	refim = fits file with astrometry to regrid to
;	newim = name to save regridded data & header to
;	INTERP = input for HASTROM to specify type of
;		interpolation: 0 for nearest-neighbor,
;		1 for bilinear, 2 for bicubic (default)
;	MISSING = input for HASTROM to specify fill-
;		value for areas where input image &
;		reference image don't overlap
; EXTERNAL CALLS:
;	hastrom (& all its dependencies),
;	sxaddpar, sxdelpar, writefits
;-
  ;;compile_opt strictarr
  ;;args = command_line_args(count=nargs)

  IF (N_elements(interp) EQ 0) THEN interp=2
  IF (N_elements(missing) EQ 0) THEN missing=0

  IF (isa(oldim,'string') EQ 1) THEN oldarr = readfits(oldim,oldhdr) ELSE $
    MESSAGE,"Error: OLDIM must be a file name"
  IF (isa(refim,'string') EQ 1) THEN refhdr = headfits(refim) ELSE $
    MESSAGE,"Error: REFIM must be a file name"
  IF (isa(newim,'string') EQ 0) THEN MESSAGE,"Error: NEWIM must be a file name"

  IF sxpar(refhdr,'naxis3') NE 0 THEN sxdelpar, refhdr, 'naxis3'

  HASTROM, oldarr, oldhdr, newarr, newhdr, refhdr, MISSING=missing, INTERP=interp

  writefits, newim, newarr, newhdr
  print, 'wrote regridded image & astrometry to '+newim
END
