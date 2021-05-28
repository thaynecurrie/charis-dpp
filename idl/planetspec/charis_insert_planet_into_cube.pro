pro charis_insert_planet_into_cube,pfname,cube,pri_header,ext_header,xpos,ypos,contrast_in,$
inputpsfcube,planetmodel,empspectrum=empspectrum,negplanet=negplanet,pick=pick,manual=manual,hardwired=hardwired,$
overheadfac=overheadfac,diagnostic=diagnostic,cube_out=cube_out,spec_out=spec_out,writecube=writecube

;subroutine to insert a ****MODEL*** planet spectrum into a CHARIS data cube.    
;--meant to be called by other programs (PSF subtraction, Forward-Modeling/Throughput Correction Codes)
;--but can work as a stand-alone program to insert a planet model into a cube
;-- 06/18/2018 - negative planet flag
;-- 06/13/2018 - added flag to put in an empirical spectrum: empspectrum

if N_PARAMS() eq 0 then begin

print,'charis_insert_planet_into_cube,pfname,cube,pri_header,ext_header,xpos,ypos,contrast_in,inputpsfcube,planetmodel,empspectrum=empspectrum,pick=pick,manual=manual'
print,'overheadfac=overheadfac,diagnostic=diagnostic,cube_out,writecube=writecube'
print,'...to perform program manually, select "pick" or "manual"'
goto,skiptoend
endif

;****hardwired version********
if keyword_set(hardwired) then begin
;skiphardwired
;****HARDWIRED STUFF ****
;---these are just dummy values.  
;hardwired test image for now ...
cube=readfits('./reduc/reg/n0004reg_cal.fits',ext=1,ext_header)
pri_header=headfits('./reduc/reg/n0004reg_cal.fits')
;hardwired contrast
contrast0=-3

;hardwired position
xp=73
yp=102
xp=double(xp)
yp=double(yp)



;hardwired model
;planetdir='~/idl_tools/ADI_dl/charisred/models/planetmodels/simple/'
planetdir=charis_path(pathname='planetdir')
;inputspec=planetdir+'whitelight.txt'
inputspec=planetdir+'simpplanet.txt'

;hardwired model PSF cube
psfcube=readfits('psfcube_med.fits',ext=1,hpsfint)
dpsf=sxpar(hpsfint,'naxis1')
goto,breakout
endif
skiphardwired:
;****end hardwired version*******

;*****pick or manual*************
if (keyword_set(pick) or keyword_set(manual)) then begin
;manually set the entries ...

inputfile=dialog_pickfile(Title="Select Datacube")
cube=readfits(inputfile,ext=1,ext_header)
pri_header=headfits(inputfile)

read,'Enter the log(contrast) at H band for the model planet: ',contrast0
read,'Enter the x position for the model planet: ',xp
read,'Enter the y position for the model planet: ',yp
xp=double(xp) 
yp=double(yp)

inputpsfcube=dialog_pickfile(Title="Select the input PSF cube")
psfcube=readfits(inputpsfcube,ext=1,hpsfint)
dpsf=sxpar(hpsfint,'naxis1')

modelpath=charis_path(pathname='planetdir')
;modelpath=charis_path(pathname='modeldir')
inputspec=dialog_pickfile(Title="Select Planet Model",PATH=modelpath)
goto,breakout
endif
;*****end pick or manual version****

xp=double(xpos)
yp=double(ypos)
;defining the PSF cube from command line.  If you are doing this, you are reading in an image array, not opening a fits file.
psfcube=inputpsfcube
;define dimension size from sz,/dim call
dpsf=(size(psfcube,/dim))[0]

;******Empirical Spectrum*****
;****this is used for throughput corrections: forward-modeling or iterative nulling of planet signal
;
;Assumes that: 1. the dimensions and wavelengths of the input spectrum match that of the data
;	       2. the units match
;	       3. the resolution is identical.   
;---> this is a prediction for exactly the measurement you should get in real CHARIS data, bypassing everything else.	 	

if keyword_set(empspectrum) then begin
;****assumes that the wavelength array of the empirical spectrum PERFECTLY matches that of the data.
readcol,empspectrum,dum,charis_model_spec_cal,/silent
;charis_model_spec_cal=empspectrum
goto,breakout
endif

;**********

;*********command line version******
;***based on data you have put in at command line, now figure out key variables***

;contrast
contrast0=contrast_in
;name of planet spectrum to input
inputspec=planetmodel

;*****end command line version******

breakout:

nsc=21   ;number of subcompanions, hardwired for now

;grab the CHARIS wavelength array
get_charis_wvlh,pri_header,lambda,filtname=filtname0
lambda*=1d-3


;normalize the PSF cube to the aperture radius in the fits header
;dpsf=15 ;in case you need to make a smaller sub-array

for ilambda=0L,n_elements(lambda)-1 do begin
 rap=sxpar(ext_header,'R_AP'+strtrim(ilambda,2))
 psfcube[*,*,ilambda]/=charis_myaper(psfcube[*,*,ilambda],dpsf/2,dpsf/2,rap)
endfor

if ~keyword_set(empspectrum) then begin
;determine how bright the spectrum should be in each filter to have the requisite H-band integrated contrast wrt star.
charis_model_spec_cal=charis_planet_photometric_calibration_calculation(pri_header,ext_header,lambda,$
modelspec=inputspec,contrast=contrast0,starspec=starspec,filtname=filtname0)
charis_star_spec=starspec
endif

;****negative planets: only use for negplanet code for now ...

if keyword_set(negplanet) then charis_model_spec_cal*=-1.

if keyword_set(diagnostic) then begin
plot,lambda,starspec,xrange=[1,2.5],xstyle=1
oplot,lambda,charis_model_spec_cal*10^(-1.*contrast0),linestyle=2,psym=-4
window,1
plot,lambda,charis_model_spec_cal,linestyle=2
print,'modelspec',charis_model_spec_cal
endif

cube_in=cube   ;the input image
cube_dim=(size(cube_in))[1]   ;image slice dimensions,assume a square array
lambda_dim=(size(cube_in))[3]  ;wavelength array length

;x,y coordinates of pixels of the PSF centered on pixel 0.0

xpsf=lindgen(dpsf)#replicate(1l,dpsf)-dpsf/2
ypsf=replicate(1l,dpsf)#lindgen(dpsf)-dpsf/2

;indices of pixels of the psf centered on pixel 0.0 in the big picture
ipsf=xpsf+ypsf*cube_dim


;now read in declination and latitude
param,'DEC',decdeg,/get,pfname=pfname
lat=float(sxpar(pri_header,'lat'))

;now read in the exposure time, coadds, hour angle, and parallactic angle
exptime=float(sxpar(pri_header,'exp1time'))
ncoadds=float(sxpar(pri_header,'coadds'))
ha=float(sxpar(pri_header,'HA'))
pa=float(sxpar(pri_header,'PA'))


xc=cube_dim/2  ;x & y center of image, assume square image slices
yc=cube_dim/2 

;**hardwired stuff**
dx=xp-xc
dy=yp-yc
rc=sqrt(dy^2.+dx^2.)

dpsf=(size(psfcube))[1]  ;again, assume square arrays for the PSF cube postage stamp

x=ha+(findgen(nsc)/(nsc-1.))*(exptime*ncoadds/3600.)   ;hour angle at each subcompanion time
pa0=parangle(x,decdeg,lat)
dpa=pa0-pa0[0]

if ~keyword_set(overheadfac) then begin
;if overheadfac is not computed, then use empirical estimate of 5.9s for readout: roughly fact 1.148 for 45s exposure
;overheadfact0=5.59
overheadfact0=0.  ;you expose and then read out the array. The PSF is smeared over the exposure, not the subsequent readout delay
overheadfac=(ncoadds*exptime+overheadfact0)/(float(ncoadds*exptime))

endif
;overheadfac=5 ;hardwired for now
dpa*=overheadfac
;center the dpa value, assume the PA quoted is the avg PA during the exposure.

if max(dpa) gt 0 then begin
dpa-=max(dpa)/2.
endif else begin
dpa-=min(dpa)/2.
endelse

;dpa[*]=0
;this is where you would loop on the planet[0] through [i_planet-1]

plan_coords=cv_coord(/double,from_rect=[dx,dy],/to_polar,/degrees)

;asc=((plan_coords[0]+mean(dpa))-(pa-allpa[0]))*!dtor
asc=(plan_coords[0]+dpa)*!dtor
xsc=rc*cos(asc)+xc
ysc=rc*sin(asc)+yc

for icube=0L,lambda_dim-1 do begin
 imslice=cube_in[*,*,icube]

 for isc=0L,nsc-1 do begin

  xe=floor(xsc[isc])  
  ye=floor(ysc[isc])
  dxf=xsc[isc]-xe & dyf=ysc[isc]-ye
  psfs=shift_sub(psfcube[*,*,icube],dxf,dyf)
  imslice[ipsf+xe+ye*cube_dim]+=psfs*charis_model_spec_cal[icube]/nsc
 endfor

 cube_in[*,*,icube]=imslice
endfor

;Diagnostic
;writefits,'cube.fits',cube_in
;writefits,'diff.fits',cube_in-cube

cube_out=cube_in
spec_out=charis_model_spec_cal

skiptoend:
end
