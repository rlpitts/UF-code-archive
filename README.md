# Hello-World
misc. code dump

I did not write coordconv.py, but formatter.py needs it to work and it's open source. See coordconv.py file for license.

I designed formatter.py specifically to homogenize the spacecraft and other MIR/FIR/sub-mm data I downloaded from IRSA and the APEX website and organize it the way I found most convenient on my workstation. One of its chief tasks is to sort the IR-to-submm data by Regions in the CHaMP survey, using the footprints of the CO survey images to create links in each Region folder to properly formatted single-extension fits files in a "workfiles" subdirectory in the directory for the data from each instrument. The program's main formatter.mover() function will not work without CHaMP files (which are proprietary), a file containing a list thereof, or NantenTab.dat to give object coordinates. It also won't work without PACS_PAsByFile.tbl or resolns.tbl--those are needed for setbeam() and addbeam(). In short, formatter.py is unlikely to be useful as-is, but individual modules may be helpful to copy wholesale or in part. I'm keeping it here primarily as an online record of my work and as a back-up. If it helps anyone else in any way, that's just icing on the cake, but I take no responsibility for an

The c-shell script imtile8.csh, and any routines it calls, are also only here as back-ups and records of work, because they require the all of following to be True:
1. All IR/submm data files have been sorted and symlinked via formatter.py. (Note to self: Formatter.addbeam() was added later so make sure that was called, too.) Miriad will have a cow if any of the files still have multiple extensions.
2. All IR/sub-mm data must have been converted to Galactic (J2000) coordinates. I used imcctran. Miriad will throw an error message if you forget but it won't tell you very clearly that that's the problem.
3. There must be a list of the CHaMP survey images to use as regridding templates in the CWD.
4. Each Region folder must contain a list of its IR/sub-mm data files/symlinks (with any you don't intend to process commented out. This is important: Idk which task it is, but one of Miriad's tasks will absolutely destroy the result if you combine 2 images in a place where one has data and another is blanked.)
5. There must be a copy of regrid.pro, imblreplace2.sh, and molist in the CWD.

Even then there's no guarantee. I'd do it all in python but there's something wrong with Astropy's convolution program because even when I make sure to specify the FWHM in sigma, the fluxes are lower than expected by a factor of several. If you can find a way to convol and regrid your files in pure python or IDL, or really any software from this century, I'd highly recommend it and would love to know about it.

sofia_utils.py and fmtfifi.csh are just scripts I hacked together quickly to convolve, regrid, and analyze our FIFI-LS data. Even more so than any other code in this repository, they are unfinished and only serve a handful of purposes. A day may come when I sort these code snippets into folders, when I learn to use class definitions, and make them suitable for more general use...but it is not this day. This day is night, and I'm going to bed.
