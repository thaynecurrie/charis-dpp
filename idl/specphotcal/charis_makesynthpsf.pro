pro charis_makesynthpsf,pfname,synthtype=synthtype,fwhm=fwhm,rmax=rmax,psfsize=psfsize,$
;fsize=fsize,
smoothval=smoothval,$
;scalefwhm=scalefwhm,$
suffname=suffname,outname=outname,help=help

;creates a synthetic PSF, a drop-in replacement for charis_makeemppsf 
;...in case you have a different type of source (e.g. a spatially extended protoplanet-y thing)

;****Critical Keywords ****
;synthtype - 'flat' [a constant intensity over some radius], 'gaussian' [a gaussian distribution with some fwhm]
;rmax - for 'flat' case, the radial extent of the flat intensity profile
;smoothval - for 'flat' case, the radius over which the edge of the profile is smoothed
;fwhm - for 'gaussian', the fwhm of the gaussian function

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"charis_makesynthpsf,pfname,synthtype=synthtype,fwhm=fwhm,rmax=rmax,psfsize=psfsize,"
print,"smoothval=smoothval,outname=outname"

print,""
print,"Example: charis_makesynthpsf,'HR8799_low.info',synthtype='gaussian',fwhm=10,outname='gauss3med.fits'"
print,"***Keywords***"
print,'*pfname - parameter file (e.g. HR8799_low.info)'
print,"psfsize - size of the synthetic PSF"
print,"synthtype - 'flat' [a constant intensity over some radius], 'gaussian' [a gaussian distribution with some fwhm]"
print,"rmax - for 'flat' case, the radial extent of the flat intensity profile"
print,"smoothval - for 'flat' case, the radius over which the edge of the profile is smoothed"
print,"fwhm - for 'gaussian', the fwhm of the gaussian function"
print,"outname - name of the output file"
goto,skiptotheend

endif

;if ~keyword_set(outname) then outname='psfcube.fits'

datadir='./reduc/'
datadir1=datadir+'reg/'
;create directory for PSF model if it is not there ...
reducdir='./psfmodel/'
file_mkdir,reducdir

;if ~keyword_set(fsize) then fsize=21

if ~keyword_set(psfsize) then psfsize=21
dist_circle,subarray_dist,psfsize

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

param,'fnum_sat',flist,/get,pfname=pfname
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
nwvlh=n_elements(lambda)
nfiles=n_elements(files)

dim=sxpar(h1,'naxis1')
xc=dim/2 & yc=dim/2

dpsf=21

if keyword_set(psfsize) then dpsf=psfsize

;dist_circle,subarray_dist,fsize
dist_circle,subarray_dist,psfsize

psf_image=fltarr(dpsf,dpsf,nwvlh)
cubecol=fltarr(dpsf,dpsf)

if ~keyword_set(synthtype) then synthtype='flat'
if ~keyword_set(smoothval) then smoothval=2

if ~keyword_set(fwhm) then begin
Dtel=charis_get_constant(name='Dtel')
fwhm=0.206*lambda/Dtel
goto,skipconstantfwhm
endif

;if ~keyword_set(scalefwhm) then begin
;if n_elements(fwhm) gt 1 then fwhm=median(fwhm,/even)
if n_elements(fwhm) eq 1 then fwhm=replicate(fwhm,nwvlh)
;endif

skipconstantfwhm:

pixscale=charis_get_constant(name='pixscale')

;assume that if you entered FWHM larger than 1, you 'meant' to do FWHM in pixels
if fwhm[0] lt 1 then fwhm/=pixscale  ;charis pixel scale

;endif
case synthtype of 
  'flat':begin 
           if ~keyword_set(rmax) then rmax=5  ;so size is ~10 lambda/D
           good=where(subarray_dist le rmax)
           cubecol[good]=1
           cubecol=smooth(cubecol,smoothval,/nan)
           cubecol/=total(cubecol)
           for il=0L,nwvlh-1 do psf_image[*,*,il]=cubecol
         end

  'gaussian':begin
           for il=0L,nwvlh-1 do begin
              psf_image[*,*,il]=psf_gaussian(npixel=dpsf,fwhm=fwhm[il],/double,/normalize)
             if keyword_set(rmax) then begin
             slice=psf_image[*,*,il]
             dist_circle,g,dpsf
             slice[where(g gt rmax)]=0
             slice/=total(slice)
             psf_image[*,*,il]=slice
             endif
           endfor
           cubecol=median(psf_image,dimension=3,/even)
           ;cubecol=psf_image[*,*,0]
             end

endcase

if ~keyword_set(outname) then begin
outname ='synthcube_med.fits'
outnamecol='synth_medcol.fits'
endif else begin
outnamecol0=strsplit(outname,'.',/extract)
outnamecol=outnamecol0[0]+'col.fits'
endelse

writefits,reducdir+outname,0
writefits,reducdir+outname,psf_image,/append

cube=median(psf_image,dimension=3,/even)
writefits,reducdir+outnamecol,cubecol

skiptotheend:
end
