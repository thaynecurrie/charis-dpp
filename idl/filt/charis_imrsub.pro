pro charis_imrsub,pfname,$
mask=mask,rmax=rmax,psfsub=psfsub,prad=prad,pang=pang,sdi=sdi,$
medbox=medbox,prefname=prefname,suffname=suffname,fc=fc,nonorthup=nonorthup,help=help

;Spatially filters a CHARIS data cube.   Limited functionality right now
;***08/30/2018 - fixed issue causing SVD to crash in PSF Subtraction routines
;*** 5/25/2017 - quick and dirty fix to adapt to CHARIS data

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"charis_imrsub,pfname,prad=prad,psfsub=psfsub,pang=pang,sdi=sdi,medbox=medbox,fc=fc,nonorthup=nonorthup"
print,""
print,"Example: charis_imrsub,'HR8799_low.info',/prad"
print,""
print,"***Keywords***"
print,'*pfname - parameter file (e.g. HR8799_low.info)'
print,"*prad - subtracts a radial profile from each cube slice"
print,"*psfsub - subtracts a median-average cube from each cube"
print,"*mask - masks angles corresponding to spider"
print,"*medbox - spatially filters cube with boxwidth of some size"
print,"*pang - subtracts radial-average [untested]"
print,"*sdi - does simple SDI subtraction for each cube"
print,"*fc - switch thrown when charis_imrsub is called in forward-modeling/synthetic companions"
print,"*nonorthup - do not rotate cube to north-up after combining"
goto,skiptotheend
endif


;A. Spatial Filtering
;1. radial profile subtraction
;2. unsharp masking

;B. PSF subtraction
;1. If "psfsub" switch is thrown, does median combination of images --> reference PSF and subtracts
;2. If "sdi" switch is thrown, does simple classical SDI subtraction of data cube.



;box size for median filtering
;if ~keyword_set(medbox) then medbox=11

; setupdir,reducdir=reducdir
reducdir='./reduc/'

;data directory
datadir=reducdir+'reg/'

;determine reduction subdirectory
subdir='rsub/'
reducdir+=subdir

;create list of filenames
;param,'obsdate',date,/get,pfname=pfname & date=strtrim(date,2)
param,'fnum_sat',flist,/get,pfname=pfname
param,'fwhm',fwhm,/get,pfname=pfname


;**** Prefix names for your files (Make sure to change these with different data!!!)**
if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then begin
test=file_search(datadir+'*reg_cal.fits')
if (n_elements(test) gt 0) then begin
suffname='reg_cal'
endif else begin
suffname='reg'
endelse

endif

filenum=nbrlist(flist)

if ~keyword_set(fc) then begin
filesout=filelist(filenum,prefix=prefname,suffix='rsub')
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
endif else begin
filesout=filelist(filenum,prefix=prefname,suffix='rsub_fc')
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname+'_fc')
endelse

;Reads in parallactic angle
readcol,'reduc.log',filenum,allxc,allyc,allrsat,allha,allpa

paranglef=fltarr(nfiles)

param,'spang*',spang,/get,pfname=pfname

;defines width of spider mask ... in case you want to use it
param,'spmask',spmask,/get,pfname=pfname

;get dim from first image header
htest=headfits(datadir+files[0])
fimagetest=readfits(datadir+files[0],htest1,ext=1)

;assume square array
dim=(size(fimagetest,/dim))[0]
sz=size(fimagetest,/dim)
d0=sz[0]
if ~keyword_set(xc0) then xc0=dim/2
if ~keyword_set(yc0) then yc0=dim/2

;Maximum radial distance for radial profile subtraction
if ~keyword_set(rmax) then rmax=dim/2


;mask the sat spots?
if keyword_set(mask) then begin
spi=mkspider(d0,spmask,spang,/justind)
endif

;Loop to subtract

imcubeout=fltarr(sz[0],sz[1],sz[2])
imgiantcube=fltarr(sz[0],sz[1],sz[2],nfiles)

for i=0L,nfiles - 1 do begin
print,'reading in file number ',i+1,' with name ',files[i]
imcube=readfits(datadir+files[i],h1cube,ext=1)
h0cube=headfits(datadir+files[i])

paranglef[i]=allpa[i]
;loop on cube element

;*************FILTERING*****
;You have the choice of doing one or more of the following
;1. subtract reference PSF (not implemented)
;2. subtract radial profile
;3. Fourier filter (not yet implemented)
;4. unsharp masking


for ii=0L,sz[2]-1 do begin

 imslice=imcube[*,*,ii]
 bad=where(imslice eq 0 or finite(imslice) eq 0,complement=good)
 imslice2=imslice
 imslice2[bad]=!values.f_nan
 ;***Radial profile subtraction

if keyword_set(prad) then begin
profrad_tc,imslice2,2,allrsat[i],rmax,p2d=pr,p1d=p1d,rayon=r
imslice[good]-=pr[good]
endif


if keyword_set(pang) then begin
profang_tc,imslice,2,allrsat[i],rmax,p2d=pr,rayon=r
imslice-=pr
endif

;***filtering
if keyword_set(lowpass) then begin
imslice[d0/2-rmax:d0/2+rmax,d0/2-rmax:d0/2+rmax]=filter_image(imslice[d0/2-rmax:d0/2+rmax,d0/2-rmax:d0/2+rmax],median=0.5*fwhm)
endif


;****''unsharp masking''
if keyword_set(medbox) then begin
;NaN-mask the zero'd regions first
imslicemask=imslice
goodslice=where(imslicemask ne 0 and finite(imslicemask) eq 1,complement=badslice)
imslicemask[badslice]=!values.f_nan
imslicemask-=filter_image(imslicemask,median=medbox)
imslice=imslicemask
;imslice[d0/2-rmax:d0/2+rmax,d0/2-rmax:d0/2+rmax]=imslice[d0/2-rmax:d0/2+rmax,d0/2-rmax:d0/2+rmax] $
;- filter_image(imslice[d0/2-rmax:d0/2+rmax,d0/2-rmax:d0/2+rmax],median=medbox)

endif

;**** Fourier Filtering
;not well understood right now - 1/18/2016
;if keyword_set(fft) then begin
;imslice=fft_filt(imslice,/highpass,boxsize=fft)
;endif

if keyword_set(mask) then imslice[spi]=!values.f_nan
imcubeout[*,*,ii]=imslice
endfor

sxaddpar,h1cube,'rsub',1
writefits,reducdir+filesout[i],0,h0cube
writefits,reducdir+filesout[i],imcubeout,h1cube,/append
imgiantcube[*,*,*,i]=imcubeout
endfor

;Now, for combining and derotating
imgiantcuberot=imgiantcube
imgiantavgcube=imgiantcube

if keyword_set(psfsub) then begin
imcuberef=median(imgiantavgcube,dimension=4,/even)
endif

;SDI
if keyword_set(sdi) then begin 
imgiantsdicubeavg=imgiantcube
imgiantsdicuberot=imgiantcube
endif

for i=0L,nfiles-1 do begin
imt=readfits(datadir+files[i],ext=1,h1cubeold)
h0cubeold=headfits(datadir+files[i])

print,'Filtering/Subtracting File Number ',long(i)
;derotating

pardiff=paranglef[0]-paranglef[i]
;print,pardiff

imuse =imgiantcuberot[*,*,*,i]
if keyword_set(psfsub) then imuse-=imcuberef

;SDI
if keyword_set(sdi) then begin
imsdi=imuse
nslices=(size(imsdi,/dim))[2]
;charis_alignspeckle,imsdi,h0cubeold,h1cubeold,refslice=nslices-1,/nolocs
;writefits,'alignspeckle.fits',imsdi

dist_circle,roi,sz[1]
good=where(roi le 70)
get_charis_wvlh,h0cubeold,wvlhinput
;wvlhratio=(wvlhinput/wvlhinput[11])^(-2.)

for ij=0L,nslices-1 do begin 
imsditest=imsdi
charis_alignspeckle,imsditest,h0cubeold,h1cubeold,refslice=ij,/nolocs
imsdiavg=median(imsditest,dimension=3,/even)
weightpsf=total(imsditest[*,*,ij]*imsdiavg,/nan)/total(imsdiavg*imsdiavg,/nan)
imsdi[*,*,ij]-=weightpsf*imsdiavg
endfor
imuse=imsdi
writefits,'ralignsdi.fits',imuse
endif

imusef=rotat_cube(imuse,1*pardiff,missing=!values.f_nan)
imgiantcuberot[*,*,*,i]=imusef
endfor


imcubeavg=median(imgiantavgcube,dimension=4,/even)
writefits,'imcubefinal_avg.fits',0,htest
writefits,'imcubefinal_avg.fits',imcubeavg,htest1,/append
imfinalavg=median(imcubeavg,dimension=3,/even)
writefits,'imcubefinal_avg_collapsed.fits',0,htest
writefits,'imcubefinal_avg_collapsed.fits',imfinalavg,htest1,/append

imcubefinal=median(imgiantcuberot,dimension=4,/even)

if ~keyword_set(nonorthup) then $
angoffset=charis_get_constant(name='angoffset')
charis_northup,imcubefinal,htest,htest1,angoffset=angoffset
writefits,'imcubefinal_derot.fits',0,htest
writefits,'imcubefinal_derot.fits',imcubefinal,htest1,/append
imfinal=median(imcubefinal,dimension=3,/even)
writefits,'imcubefinal_derot_collapsed.fits',0,htest
writefits,'imcubefinal_derot_collapsed.fits',imfinal,htest1,/append

skiptotheend:

end



