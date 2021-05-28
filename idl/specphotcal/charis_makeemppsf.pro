pro charis_makeemppsf,pfname,pick=pick,unsat=unsat,filt=filt,fsize=fsize,prad=prad,psfsize=psfsize,nonorm=nonorm,nozero=nozero,ladder=ladder,$
suffname=suffname,outname=outname,help=help

;make an empirical PSF from a set of registered images
;12-26-2018 - switch to use unsaturated PSF instead of satellite spots
;02-05-2018 - now puts empirical PSF in its own directory for later use.
;10-11-2017 - added 'nonorm' switch to not normalize the PSF.  Useful for a quick and dirty flux calibration estimate from an already-processed data set.
;added 'zero' switch to ensure that the background is zero'd

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"charis_makeemppsf - create an empirical PSF model from the sat spots or unsaturated star"
print,""
print,"charis_makeemppsf,pfname,pick=pick,unsat=unsat,filt=filt,fsize=fsize,prad=prad,psfsize=psfsize,ladder=ladder,"
print,"nonorm=nonorm,nozero=nozero,suffname=suffname,outname=outname"
print,""
print,"Example: charis_makeemppsf,'HR8799_low.info'"
print,"*** Primary Keywords ***"
print,'*pfname - parameter file (e.g. HR8799_low.info)'
print,"ladder - construct a PSF model from ladder images (e.g. when your main science seq. lacks sat spots)"
print,"pick - just pick a single cube for a PSF model"
print,"unsat - create a PSF model from an unsaturated star, not sat spots"
print,"psfsize - size of the PSF model in x and y"
print,"filt - spatial filtering"
print,"fsize - boxsize for spatial filtering"

goto,skiptotheend
endif

if ~keyword_set(outname) then outname='psfcube.fits'

datadir='./reduc/'
datadir1=datadir+'reg/'
;create directory for PSF model if it is not there ...
reducdir='./psfmodel/'
file_mkdir,reducdir

if ~keyword_set(fsize) then fsize=21
dist_circle,subarray_dist,fsize

if ~keyword_set(pick) then begin

prefname='n'
if ~keyword_set(suffname) then begin 
test=file_search(datadir1+'*reg_cal.fits')
if n_elements(test) gt 0 then begin
suffname='reg_cal'
endif else begin
suffname='reg'
endelse
endif

if ~keyword_set(ladder) then begin
param,'fnum_sat',flist,/get,pfname=pfname
endif else begin
param,'fnum_lad',flist,/get,pfname=pfname
endelse

filenum=nbrlist(flist)
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)

test=readfits(datadir1+files[0],/exten,h1)
h0=headfits(datadir1+files[0],ext=0)

endif else begin

files=dialog_pickfile(Title="Pick Files",/multiple)
test=readfits(files[0],/exten,h1)
h0=headfits(files[0],ext=0)

endelse

;Now Get the Wavelength Vector
get_charis_wvlh,h0,wvlhs
lambda=wvlhs*1d-3

Dtel=charis_get_constant(name='Dtel')
;Dtel = 7.9d0 ;effective Subaru pupil diameter
fwhm=0.206*lambda/Dtel
pixscale=charis_get_constant(name='pixscale')
fwhm/=pixscale  ;charis pixel scale

nwvlh=n_elements(lambda)

nfiles=n_elements(files)

dim=sxpar(h1,'naxis1')
xc=dim/2 & yc=dim/2


dpsf=21
if keyword_set(psfsize) then dpsf=psfsize

nspot=4

psf_spot=fltarr(dpsf,dpsf,nspot)
psf_imagelam=fltarr(dpsf,dpsf,nwvlh,nfiles)


for il=0L,nwvlh-1 do begin

print,'analyzing wavelength slice ',il+1

for nf=0L,nfiles-1 do begin
if ~keyword_set(pick) then begin
h0=headfits(datadir1+files[nf],/silent,ext=0)
im=(readfits(datadir1+files[nf],h1,/exten,/silent))[*,*,il]
endif else begin
h0=headfits(files[nf],/silent,ext=0)
im=(readfits(files[nf],h1,/exten,/silent))[*,*,il]
endelse

if ~keyword_set(unsat) then begin
bad=where(im le 0 or finite(im) eq 0,complement=good)
im2=im
im2[bad]=!values.f_nan
profrad_tc,im2,p2d=pr
im[good]-=pr[good]

if keyword_set(filt) then begin
im-=filter_image(im,median=fsize)
endif


;sat spot
;psf_spot=fltarr(dim,dim,4)
for spot =0,3 do begin
 hdrval=double(strsplit(sxpar(h1,'SATS'+strtrim(il,2)+'_'+strtrim(spot,2)),' ',/extract))
 shifts=hdrval-floor(hdrval)

;break out of computation if NaN value for fits header

if ((size(hdrval))[1]) lt 2 then begin
print,"You need data cubes with sat spots!"
goto,breakoutspot
endif

if (hdrval[0] gt 1000. or hdrval[1] gt 1000. $
or finite(hdrval[0]) eq 0 or finite(hdrval[1]) eq 0) then begin
psf_spot[*,*,spot]=!values.f_nan
print,"You need data cubes with sat spots!"
goto,breakoutspot
endif

 s=subarr(im,dpsf,floor([hdrval[0],hdrval[1]]))

 psf_spot[*,*,spot]=shift_sub(s,-1*shifts[0],-1*shifts[1])


 if ~keyword_set(nozero) then begin
 dum=psf_spot[*,*,spot]
 bckgd=median(dum[where(subarray_dist gt 3*fwhm[il])],/even)
 dum-=bckgd
 psf_spot[*,*,spot]=dum
 endif

;get rid of zeroes, this usually messes things up
dum=psf_spot[*,*,spot]
bad=where(dum lt -0.01*max(dum))
dum[bad]=0
psf_spot[*,*,spot]=dum

if ~keyword_set(nonorm) then begin
psf_spot[*,*,spot]/=total(psf_spot[*,*,spot],/nan)
endif
breakoutspot:


endfor


psf_imagelam[*,*,il,nf]=median(psf_spot,dimension=3,/even)

endif else begin

;unsaturated, no spots
s=subarr(im,dpsf,[xc,yc])

if ~keyword_set(nozero) then begin
dum=s
bckgd=median(dum[where(subarray_dist gt 2*fwhm[il])],/even)
dum-=bckgd
s=dum
endif
psf_imagelam[*,*,il,nf]=s

endelse


endfor
endfor

if nfiles gt 1 then begin
psf_image=median(psf_imagelam,/even,dimension=4)
endif else begin
for il=0L,nwvlh-1 do begin
psf_imagelam[*,*,il]=fixpix_rs(psf_imagelam[*,*,il],iter=3)
endfor
psf_image=psf_imagelam
endelse

if ~keyword_set(nonorm) then begin
for il=0L,nwvlh-1 do psf_image[*,*,il]/=total(psf_image[*,*,il],/nan)
endif

if ~keyword_set(nonorm) then begin
writefits,reducdir+'psfcube_med'+'.fits',0
writefits,reducdir+'psfcube_med'+'.fits',psf_image,/append
endif else begin
writefits,reducdir+'psfcube_med_nonorm'+'.fits',0
writefits,reducdir+'psfcube_med_nonorm'+'.fits',psf_image,/append
endelse


cubecol=median(psf_image,dimension=3,/even)
if ~keyword_set(nonorm) then begin
writefits,reducdir+'psf_medcol.fits',cubecol
endif else begin
writefits,reducdir+'psf_medcolnonorm.fits',cubecol
endelse


skiptotheend:
end
