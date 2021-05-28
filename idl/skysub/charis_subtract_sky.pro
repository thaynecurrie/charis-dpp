pro charis_subtract_sky,pfname,medbox=medbox,prefname=prefname,suffname=suffname,pick=pick,scalesky=scalesky,scaleindchan=scaleindchan,skyrad=skyrad,skylim=skylim,ladder=ladder,verbose=verbose,help=help

;****Performs sky subtraction on CHARIS data cubes.  Useful for removing detector noise (and sky noise at K band)
;09/15/19 - Second version.   Adds switch to read sky frames from .info file
;08/16/17 - First version of sky subtraction routine.  Asks user to identify sky frames, subtracts them.  Puts in suffix of 'skysub'.

;Procedure:
;*1. Chooses sky frames: From parameter file or asks user to identify sky frames --> creates a master sky cube from median at each pixel
;*2. Loops the data cubes, performs sky subtraction

;****
if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"charis_subtract_sky,pfname,prefname=prefname,suffname=suffname,pick=pick,scalesky=scalesky,scaleindchan=scaleindchan,skyrad=skyrad,verbose=verbose,help=help"
print," "
print,"Example: charis_subtract_sky,'HR8799_low.info',/scalesky"
print,""
print,"***Important Keywords***"
print,'*pfname - parameter file (e.g. HR8799_low.info)'
print,"pick - Pick the sky frames from GUI"
print,"scalesky - weight the sky frames by the relative median sky value of the red-most channel [outer regions]"
print,"scaleindchan-weight the sky frames by the relative median sky value for *each* channel [outer regions]"
print,"skyrad - radius beyond which to estimate sky contribution"
print,"skylim - limit in cts/s to define 'channels with substantial sky emission'"
print,"medbox [currently not included] - high-pass filter the data cube after subtracting sky"
print,"***"
goto,breakout
endif

; setupdir,reducdir=reducdir
reducdir='./reduc/'

;data directory
;datadir=reducdir+'expand/'
datadir=reducdir+'prep/'

;determine reduction subdirectory
;subdir='expand/'
subdir='prep/'
reducdir+=subdir

;create list of filenames

if ~keyword_set(ladder) then begin
param,'fnum_sat',flist,/get,pfname=pfname
endif else begin
param,'fnum_lad',flist,/get,pfname=pfname
endelse

param,'fwhm',fwhm,/get,pfname=pfname

;**** Prefix names for your files (Make sure to change these with different data!!!)**
if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='e'

if ~keyword_set(skyrad) then skyrad=60
if ~keyword_set(skylim) then skylim=25

filenum=nbrlist(flist)
filesout=filelist(filenum,prefix=prefname,suffix=suffname+'_skysub')
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)

;get dim from first image header
htest=headfits(datadir+files[0])
fimagetest=readfits(datadir+files[0],htest1,ext=1)

;assume square array
dim=(size(fimagetest,/dim))[0]
sz=size(fimagetest,/dim)
d0=sz[0]

;1. ****Select sky frames, Computes Sky Cube

if keyword_set(pick) then begin
skyframes=dialog_pickfile(Title="Select Sky Frames",/multiple_files)
nsky=n_elements(skyframes)
skycube=fltarr(sz[0],sz[1],sz[2],nsky)
skytime=fltarr(nsky)

;Loop to define the sky cube
for i=0L,nsky-1 do begin
a=readfits(skyframes[i],ext=1,h1sky)
h0sky=headfits(skyframes[i])
exp1time=sxpar(h0sky,'exp1time')
coadds=sxpar(h0sky,'coadds')
skytime[i]=exp1time*coadds
skycube[*,*,*,i]=a
endfor

endif else begin
;reading sky frames from the .info file 
param,'fnum_sky',flistsky,/get,pfname=pfname
filenumsky=nbrlist(flistsky)
filessky=filelist(filenumsky,nsky,prefix=prefname,suffix=suffname)

skycube=fltarr(sz[0],sz[1],sz[2],nsky)
skytime=fltarr(nsky)

;Loop to define the sky cube
for i=0L,nsky-1 do begin
a=readfits(datadir+filessky[i],ext=1,h1sky)
h0sky=headfits(datadir+filessky[i])
exp1time=sxpar(h0sky,'exp1time')
coadds=sxpar(h0sky,'coadds')
skytime[i]=exp1time*coadds
skycube[*,*,*,i]=a
endfor

endelse

;in case you have variable integration times for your sky frames
;not sure you do this!
;for i=0L,nsky-1 do skycube[*,*,*,i]*=(max(skytime)/skytime[i])

;now collapse to get a master sky
;if n_elements(skycube gt 1) then begin
if nsky gt 1 then begin
;stop
master_skycube=median(skycube,dimension=4,/even)
endif else begin
master_skycube=skycube
endelse

;2. ****Loop to perform sky subtraction on science data

for i=0L,nfiles-1 do begin

imcube=readfits(datadir+files[i],h1cube,ext=1)
sz=size(imcube,/dim)
h0cube=headfits(datadir+files[i])
exp1time=sxpar(h0cube,'exp1time')
coadds=sxpar(h0cube,'coadds')
exptime_sci=exp1time*coadds


;subtract sky.   For now do a simple subtraction since the PSF halo contaminates an estimate of the sky background
if (~keyword_set(scalesky) and ~keyword_set(scaleindchan)) then begin
imcube-=master_skycube*exptime_sci/max(skytime)
endif else begin
;take the last channel (the most contaminated one in low-res) scale by average signal
dist_circle,gg,201

if keyword_set(scalesky) then begin
targslice=imcube[*,*,sz[2]-1]
skyslice=master_skycube[*,*,sz[2]-1]
good=where(gg gt skyrad and targslice ne 0)
scalefactor=median(targslice[good])/median(skyslice[good],/even)
imcube-=master_skycube*scalefactor
endif else begin
for il=0L,sz[2]-1 do begin
targslice=imcube[*,*,il]
good=where(gg gt skyrad and targslice ne 0)
skyslice=master_skycube[*,*,il]
if median(targslice[good],/even) gt skylim  and median(skyslice[good]) gt skylim then begin
scalefactor=median(targslice[good])/median(skyslice[good],/even)
endif else begin
scalefactor = 1.0
endelse

imcube[*,*,il]-=master_skycube[*,*,il]*scalefactor
print,scalefactor,il,median(targslice[good]),median(skyslice[good],/even)
endfor

endelse

endelse



writefits,reducdir+filesout[i],0,h0cube
writefits,reducdir+filesout[i],float(imcube),h1cube,/append

if keyword_set(verbose) then writefits,'mastersky.fits',master_skycube
endfor

breakout:

end

