pro charis_specphot_cal_unsat,pfname,calmethod=calmethod,pick=pick,datacube=datacube,calcube=calcube,subskyannulus=subskyannulus,$
filtcal_slice=filtcal_slice,nopradcal=nopradcal,$
meancomb=meancomb,$
fixradius=fixradius,$
ap_factor=ap_factor,$
;attenfac=attenfac,ap_factor=ap_factor,$
filtername=fname,$
starlib=starlib,$
empspectrum=empspectrum,$
av=av,$
prefname=prefname,suffname=suffname,$
verbose=verbose,$
fluxunits=fluxunits,outfilename=outfilename,help=help

;Spectro-Photometric Calibration program (work in progress)
; - 2/7/2020 - Major overhaul.   Allows 0) cal on indiv cubes, 1) an average cube, 2) another selected cube
; - 11/09/2017 - still some possible issues in absolute calibration
; - uses images with satellite spots, info file with spectral type and brightness of star in H band
; - feeds this information into 'calculation' subroutine to get a flux-calibrated spectrum of your star
; - computes satellite spot attenuation based on the amplitude of modulation and the MJD date
; - flux calibrates data cubes based on model spectrum and on sat spots
;-... added information about the star's brightness/channel in ext=1 of the fits header, option to spatially filter cal slice/turn off prad removal - 11/09/2017
;-...'calculation.pro subroutine done -7/12/2017
;-assumes same exposure time right now -7/12/2017
;
;*to implement ...
; switches on 'attenfac', 'starname','filter','mag', etc ... to ...
; override the choice of star name, filter, brightness

if (N_PARAMS() eq 0 or keyword_set(help)) then begin

Print,'Calibrate Photometric Flux subroutine'
Print,'This program does flux calibration for a CHARIS datacube or sequence of data cubes'
Print,' using data of a known star taken with a similar setup
Print,'  '
Print,'****Calling Sequence ****'
Print,'  '
Print,'charis_specphot_cal_unsat,pfname,calmethod=calmethod,pick=pick,datacube=datacube,calcube=calcube,modamp=modamp,subskyannulus=subskyannulus,'
Print,'filtcal_slice=filtcal_slice,meancomb=meancomb,ap_factor=ap_factor,filtername=fname,starname=starname,mag=mag,starlib=starlib,'
Print,'av=av,prefname=prefname,suffname=suffname,test=test,verbose=verbose,fluxunits=fluxunits,outfilename=outfilename,help=help,guide=guide'
Print,' '
Print,"Example: charis_specphot_cal_unsat,'LkCa15_low.info',calmethod=2,empspectrum='LkCa15empirical.fits' [would use an empirical spec instead of library spec]"
Print,'Parameters ...'
Print,'pfname - input info file for applying correction to pre-determined sequence of registered images'
Print,'calmethod - spectrophotometric calibration method: 0) cal on indiv cubes, 1) an average cube, 2) another selected cube'
Print,'pick - manually select images to be calibrated'
Print,'datacube - Data to Calibrate (if you define manually)'
Print,'calcube - Datacube used for Flux Calibration (if you select manually)'
Print,'subskyannulus - subtract the local "sky" background for sat spot photometry (usually no)'
Print,'filtcal_slice - spatially filter (subtract median filter) before sat spot photometry (usually no)'
Print,'fixradius - fix the aperture radius to a constant instead of scaling with wavelength
Print,'meancomb - mean-combine the calibration cubes instead of median'
Print,'ap factor - rescale the aperture for photometry'
Print,'filtername - set manual value for the filter (low, JHK)'
Print,'star lib - 0,1,2 (Pickles, Kurucz or Manual)
Print,'pref/suffname - Manually override the prefix and suffix of file names used (useful for the pfname switch but applied to images other than "reg" images)'
print,'fluxunits - What Flux Density units? 0=mJy,1=Jy,2=W/m^2/um,3=ergs/s/cm^2/A,4=ergs/s/cm^2/Hz,5=contrast'
;print,''
;print,'Limitations ...'
goto,skiptoend
endif

;Aapplies a photometric calibration either to an entire sequence of files or to one or more identified image

;default is to apply to sequence of registered but not PSF subtracted images
;allow for override by choosing 'pick'

;***STEP 1: Preliminary Stuff***
;define data directory
reducdir='./reduc/'

;data directory
datadir=reducdir+'reg/'
subdir='reg/'
reducdir+=subdir

;information about the target properties.
param,'JMAG',jmag,/get,pfname=pfname
param,'EJMAG',ejmag,/get,pfname=pfname
param,'HMAG',hmag,/get,pfname=pfname
param,'EHMAG',ehmag,/get,pfname=pfname
param,'KMAG',kmag,/get,pfname=pfname
param,'EKMAG',ekmag,/get,pfname=pfname
param,'SPTYPE',spt,/get,pfname=pfname

;****flag in case you did not give the luminosity class
if strlen(spt) eq 2 then begin
spt=spt+'V'
print,'*****WARNING!!!!!! Assuming a Dwarf Luminosity Class!!!!!******'
endif

;******Which Stellar Library Do you Want to Use?*****
;0 - Pickles 1998 PASP library (empirical, incomplete data)
;1 - Kurucz 1993 (theoretical models)
;2 - empirical model (user-selected)

if ~keyword_set(starlib) then begin
read,"Select the Stellar Library Source (0 = Pickles, 1 = Kurucz, 2 = empirical model)",slib
endif else begin
slib=starlib
endelse

;if you want to keep the aperture equal to the image slice FWHM then do nothing.  If you want to change the scaling set the ap_factor keyword
if ~keyword_set(ap_factor) then ap_factor = 1

;********STEP 2: CALIBRATION CUBE *********
if ~keyword_set(calmethod) then begin
calmethod=0
endif else begin
calmethod = long(calmethod)
endelse

case calmethod of 
;the median average of registered cubes
 0: begin

if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='reg'

    param,'fnum_sat',flist,/get,pfname=pfname
    filenum=nbrlist(flist)
    calcube=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
imcaltest=readfits(datadir+calcube[0],h1test,/exten)
h0test=headfits(datadir+calcube[0])

    end

 1: begin

if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='reg'

    param,'fnum_sat',flist,/get,pfname=pfname
    filenum=nbrlist(flist)
    calcube=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
imcaltest=readfits(datadir+calcube[0],h1test,/exten)
h0test=headfits(datadir+calcube[0])

    end

 2: begin

    if ~keyword_set(calcube) then $
calcube=dialog_pickfile(Title="Select Your Datacube Used for Flux Calibration",/multiple_files)
;***for now just allow for one cube. change this later!!!
imcaltest=readfits(calcube[0],h1test,/exten)
h0test=headfits(calcube[0])

    end
endcase

ncal=n_elements(calcube)
dimcal=(size(imcaltest,/dim))

;*****astrogrid spot calibration*****

;if look for the ASTROGRID modulation amplitude in the fits header

;***star to satellite contrast ratio*********
;*Assumptions: the wavelength range for the cal file is the same as the science data you are wanting to calibrate.
;;             the calibration hasn't changed.  If it has you will have to enter the right value manually.
;;This has been determined for the 1550nm spot using the latest values for modulation. 
;;;the precision is good to ~2% over a factor of 2 mod size (~factor of 4 brightness)
;;;;this strictly speaking has only been determined for the low-res mode.  Need to test high res.

;First, get the CHARIS wavelengths
get_charis_wvlh,h0test,wvlh_test


;Now get the time of observation in modified Julien Days
mjd_in=sxpar(h0test,'MJD',count=mjdcount)
if mjdcount eq 0 then begin  ;if for some reason MJD is not printed, then enter this manually
Read,'Are these data taken after August 28 2017 (if yes, type 0)?',mjdanswer
if mjdanswer eq 0 then mjd_in=58000.0
endif


;Now from the modulation size, get the CHARIS attenuation factor/channel
attenfac=1.

;***** STEP 3: Select Datacubes to Flux-Calibrate ******

;*If you want to manually set it then do this

;if you are calling the to-be-flux-cal'd file from command line or calling this program from the 1d_spectrum program ...
if keyword_set(datacube) then begin
files_full_path=datacube
nfiles=n_elements(files_full_path)
files=files_full_path
filesbase=files_full_path

;pick out the directory under which your files reside.  WARNING: assumes that all files are in same subdirectory!!!!!
your_path=strpos(files_full_path[0],'/',/reverse_search)+1

;and put the new files exactly where they were before; when you open the file later, look in the same place.
reducdir=strmid((files_full_path[0]),0,your_path)
datadir=reducdir

;***trimming off path
for i=0L,nfiles-1 do begin
z=strsplit(files_full_path[i],'/',/extract)
files[i]=z[n_elements(z)-1]
suffposition=strpos(files[i],'.fits',/reverse_search)
filesbase[i]=strmid(files[i],0,suffposition)
endfor

suffname2='_cal'
goto,breakout
endif

if ~keyword_set(pick) then begin
;assume we are working in the 'reg' subdirectory
param,'fnum_sat',flist,/get,pfname=pfname
filenum=nbrlist(flist)

if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='reg'

;input files
suffname2=suffname+'_cal'
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)

;assume a simple overwrite
filesout=filelist(filenum,nfiles,prefix=prefname,suffix=suffname2)

goto,breakout
endif else begin

;******Hack-tastic inelegant solution until I figure out a cleaner way of parsing the IDL path.
files_full_path=dialog_pickfile(Title="Select Datacubes to Flux Calibrate",/multiple)
nfiles=n_elements(files_full_path)
files=files_full_path
filesbase=files_full_path

;pick out the directory under which your files reside.  WARNING: assumes that all files are in same subdirectory!!!!!
your_path=strpos(files_full_path[0],'/',/reverse_search)+1

;and put the new files exactly where they were before; when you open the file later, look in the same place.
reducdir=strmid((files_full_path[0]),0,your_path)
datadir=reducdir
;print,'REDUCDIR ',reducdir

;***trimming off path
for i=0L,nfiles-1 do begin
z=strsplit(files_full_path[i],'/',/extract)
files[i]=z[n_elements(z)-1]
suffposition=strpos(files[i],'.fits',/reverse_search)
filesbase[i]=strmid(files[i],0,suffposition)
endfor

suffname2='_cal'

goto,breakout

endelse

;****Selecting Data Cubes: Done****;

breakout:

;******STEP 4: Determine Flux Calibration Scaling ******;

;*****Now perform photometry on the calimage, record its integration time
;***open the file, pull the fits header from the file, measure the integrated signal, measure the noise, do snrcalc/wavelength slice.
ncal=n_elements(calcube)
;test calimage

print,calcube


;not doing the 'goodcode hex2bin stuff with GPI. assume all sat spots are good

imgiantcube=fltarr(dimcal[0],dimcal[1],dimcal[2],ncal)

texp_cal=sxpar(h0test,'exp1time')
ncoadd_cal=sxpar(h0test,'coadds')
;tint_cal=texp_cal*float(ncoadd_cal)
tint_cal=fltarr(ncal)

;Telescope Diameter for Subaru
Dtel=charis_get_constant(name='Dtel') ;7.9d0 ;visible pupil for SCExAO

;pixel scale 16.4 mas/pixel
pixscale=charis_get_constant(name='pixscale') ;0.0164 nominally
fwhm=1.0*(1.d-6*(wvlh_test*1d-3)/Dtel)*(180.*3600./!dpi)/pixscale

;Now nominally assume that the aperture radius is 1/2 the FWHM
aperrad=double(fwhm*0.5)
if keyword_set(verbose) then print,aperrad[0],fwhm[0]

if keyword_set(verbose) then print,'fwhm is',fwhm

for ical=0L,ncal-1 do begin

case calmethod of
2: begin
h0cal=headfits(calcube[ical])
imcal=readfits(calcube[ical],h1cal,ext=1)
   end
else: begin

h0cal=headfits(datadir+calcube[ical])
imcal=readfits(datadir+calcube[ical],h1cal,ext=1)
      end
endcase

texp=sxpar(h0cal,'exp1time')
ncoadd=sxpar(h0cal,'coadds')
tint_cal_indiv=texp*float(ncoadd)
tint_cal[ical]=tint_cal_indiv

;if the integration time for this particular cube is different than your test one ...
;time_ratio=tint_cal_indiv/tint_cal
;print,'time ratio is ',time_ratio
;if abs(time_ratio-1.0) gt 0.1 then $
;imcal/=time_ratio

for il=0L,dimcal[2]-1 do begin
imslice=imcal[*,*,il]

if keyword_set(filtcal_slice) then begin
imslice-=filter_image(imslice,median=10*fwhm[il])
endif

imgiantcube[*,*,il,ical]=imslice

endfor

endfor

;****Create master calibration cube if ther are more than one initial cal cubes

case calmethod of

 0: begin
    end

 1:begin
if ncal gt 1 then begin

if ~keyword_set(meancomb) then begin
imcalcube=median(imgiantcube,dimension=4,/even)
tint_calf=median(tint_cal,/even)
endif else begin
imcalcube=mean(imgiantcube,dimension=4,/nan,/double)
tint_calf=mean(tint_cal,/nan,/double)
endelse

endif else begin

imcalcube=reform(imgiantcube,dimcal[0],dimcal[1],dimcal[2])
tint_calf=tint_cal[0]
endelse

ncal=1
   end

 2:begin

if ncal gt 1 then begin

if ~keyword_set(meancomb) then begin
imcalcube=median(imgiantcube,dimension=4,/even)
tint_calf=median(tint_cal,/even)
endif else begin
imcalcube=mean(imgiantcube,dimension=4,/nan,/double)
tint_calf=mean(tint_cal,/nan,/double)
endelse

endif else begin

imcalcube=reform(imgiantcube,dimcal[0],dimcal[1],dimcal[2])
tint_calf=tint_cal[0]
endelse

;writefits,'imcalcube.fits',0,h0test
;writefits,'imcalcube.fits',imcalcube,h1test,/append

ncal=1
   end
endcase

satflux=fltarr(ncal,dimcal[2])
esatflux=fltarr(ncal,dimcal[2])


;***Get the CHARIS wavelengths***
if ~keyword_set(filtername) then begin
get_charis_wvlh,h0test,wvlh_charis,filtname=filtname0
endif else begin
filtname0=fname
endelse

;Now, decision tree: set the star's magnitude to be magnitude determined by filter 
;H for lowres/broadband and H, J for J, K for K 
 case filtname0 of 
  'lowres': begin
   starmagval=hmag
   end

   'broadband':begin
    starmagval=hmag
    end
   
   'J':begin
    starmagval=jmag
   end

   'H': begin
    starmagval=hmag
    end
   
   'K': begin
    starmagval=kmag
    end 
 endcase
    


for ical=0L,ncal-1 do begin

case calmethod of 
  0: begin
imcalcube=imgiantcube[*,*,*,ical] 
     end
  else:begin
       end
endcase


;****Initialize the sat spot fluxes
;**now you have the positions of the satellite spots.  time to extract the photometry.
;assume that you don't have any errors triggered like can happen with the GPI pipeline


for il=0L,dimcal[2]-1 do begin
imslice=imcalcube[*,*,il]

if ~keyword_set(fixradius) then begin
aperradf=aperrad[il]*ap_factor
endif else begin
aperrad[il]=fixradius
aperradf=fixradius
endelse

if keyword_set(verbose) then print,aperradf,fwhm[il],fwhm[il]/2.,ap_factor,aperrad[il]*ap_factor,aperrad[il]

;aperture photometry at each slice assuming a zero background
phpadu=1.0
;sat_skyrad=[2,6]*aperradf
sat_skyrad=[4,10]*aperradf ;set this to larger, since we are flux-cal'ing off of the primary

if keyword_set(verbose) then  print,cens[0,0,il],cens[1,0,il],imslice[cens[0,0,il],cens[1,0,il]]

if keyword_set(subskyannulus) then begin
if keyword_set(setskyannulus) then sat_skyrad=setskyannulus
aper, imslice,dimcal[0]/2,dimcal[1]/2,flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
 ;print,'flux and sky1 is',flux,sky
endif else begin

aper, imslice,dimcal[0]/2,dimcal[1]/2,flux,eflux,sky,skyerr,phpadu,aperradf,setskyval=0,[0,0],/flux,/exact,/nan,/silent
     ;sat1flux[il]=flux
 ;print,'flux and sky1 is',flux,sky

endelse

       satflux[ical,il]=flux
       esatflux[ical,il]=eflux


endfor
       print,satflux[ical,0],ical
       print,esatflux[ical,0],ical
endfor

if n_elements(ncal) gt 1 then begin

satfluxf=reform(median(satflux,/even,dimension=1))
esatfluxf=reform(median(esatflux,/even,dimension=1))
endif else begin
satfluxf=reform(satflux)
esatfluxf=reform(esatflux)
endelse

snratio_sats=satfluxf/esatfluxf

if keyword_set(verbose) then print,snratio_sats

if keyword_set(verbose) then begin
window,1
;ploterror,wvlh_charis*1.d-3,satflux,esatflux,xrange=[1,3]
ploterror,wvlh_charis*1.d-3,satflux*(wvlh_charis/wvlh_charis[0])^(-1*2.),esatflux,xrange=[1,3]
oplot,wvlh_charis*1.d-3,satflux[0]*(wvlh_charis/wvlh_charis[0])^(-1*4.),linestyle=1
endif

     ;*********Flux Units 
funitlist=[0,1,2,3,4,5]   ;which corresponds to ...
funitname=['mJy','Jy','W/m^2/um','ergs/s/cm^2/A','ergs/s/cm^2/Hz','contrast']

     ;*****Normalize over which bandpass???
;1. use a piece of get_charis_wvlh to determine which normalization to use


if ~keyword_set(fluxunits) then begin
fluxunits_cube= funitlist[0]
funitname_output=funitname[0]
endif
if keyword_set(fluxunits) then begin
fluxunits_cube=funitlist[fluxunits]
funitname_output=funitname[fluxunits]
endif

;***Flux-Calibrated Model Spectrum
if fluxunits_cube le 4 then begin
 ;Physical Units

speccal=charis_photometric_calibration_calculation(h0cal,h1cal,wvlh_charis*1d-3,spectype=spt,star_mag=starmagval,filtname=filtname0,units=fluxunits_cube,stellarlib=slib,empspectrum=empspectrum,av=av)
;speccal=charis_photometric_calibration_calculation(h0cal,h1cal,wvlh_charis*1d-3,spectype=spt,star_mag=starmagval,filtname=filtname0,units=fluxunits_cube,stellarlib=slib,av=av)
;speccal=charis_photometric_calibration_calculation(h0cal,h1cal,wvlh_charis*1d-3,spectype=spt,star_mag=hmag,units=fluxunits_cube,/diagnostic)

endif else begin
 ;Contrast Units
speccal=findgen(n_elements(wvlh_charis))*0+1
endelse

print,'speccal',speccal
print,'sats',satfluxf/attenfac
plot,wvlh_charis*1d-3,speccal
oplot,wvlh_charis*1d-3,satfluxf/attenfac


;*****Flux-Calibration in Physical Units
;;;conversion factor/channel: multiply the cube slices by this

conv_fact=fltarr(nfiles,dimcal[2])

case calmethod of
0: begin
satfluxf=satflux
   end
else: begin
satfluxf=replicate(1,nfiles)#satfluxf

      end
endcase

if keyword_set(verbose) then print,'satfluxf ', satfluxf
print,'speccal',speccal


;;;flux-cal-ing the data cubes

for icube=0L,nfiles-1 do begin

conv_fact[icube,*]=(1.0/satfluxf[icube,*])*speccal*attenfac

 a=readfits(datadir+files[icube],ext=1,h1)
 h0=headfits(datadir+files[icube])

texp=sxpar(h0,'exp1time')
ncoadd=sxpar(h0,'coadds')
tint_cal_indiv=texp*float(ncoadd)


 print,'the file name is',files[icube]
 print,'the data directory is ',datadir

 for ichannel=0L,n_elements(wvlh_charis)-1 do begin
 a[*,*,ichannel]*=conv_fact[icube,ichannel]
 endfor

time_ratio=1.0 ; I think CHARIS computes in counts/s
case calmethod of 
 1: begin
;if the integration time for this particular cube is different than your test one ...
;time_ratio=tint_calf/tint_cal_indiv
a[*,*,*]*=time_ratio
    end
 2: begin
;time_ratio=tint_calf/tint_cal_indiv
a[*,*,*]*=time_ratio
    end
 else: begin
     end
endcase


;print,time_ratio

 ;funitname_output=strtrim(funitname_output)

sxaddpar,h1,'FLUXUNIT',strtrim(funitname_output,1),"Flux Density Units in Cube"

;Flux Density of Star
for ichannel=0L,n_elements(wvlh_charis)-1 $
 do sxaddpar,h1,'FSTAR_'+strtrim(ichannel,2),speccal[ichannel]," Star Flux"


;aperture
for ichannel=0L,n_elements(wvlh_charis)-1 $
 do sxaddpar,h1,'R_AP'+strtrim(ichannel,2),aperrad[ichannel]*ap_factor,"Aperture Radius (Pixels)"

;sky annulus
;for ichannel=0L,n_elements(wvlh_charis)-1 $
; do sxaddpar,h1,'R_AP'+strtrim(ichannel,2),aperrad[ichannel]*ap_factor,"Aperture Radius (Pixels)"
sxaddpar,h1,'Sky_In',strtrim(strc(2)),"Inner Radius for Sky Annulus (in Units of Aperture Radius)"
sxaddpar,h1,'Sky_Out',strtrim(strc(6)),"Outer Radius for Sky Annulus (in Units of Aperture Radius)"

;scaling
for ichannel=0L,n_elements(wvlh_charis)-1 $
 do sxaddpar,h1,'FSCALE'+strc(ichannel),conv_fact[icube,ichannel],"scale to convert counts to "+strc(funitname_output)

;error in spectrophotometric calibration/slice
for ichannel=0L,n_elements(wvlh_charis)-1 $
 do sxaddpar,h1,'CERR'+strc(ichannel),1/snratio_sats[ichannel],"Cal Fractional error for slice "+strc(ichannel)

;temporary hack of file path 
 if (keyword_set(pick) or keyword_set(datacube)) then begin
 print,'file! ',reducdir+filesbase[icube]+suffname2+'.fits'
 writefits,reducdir+filesbase[icube]+suffname2+'.fits',0,h0
 writefits,reducdir+filesbase[icube]+suffname2+'.fits',a,h1,/append
 
 if keyword_set(datacube) then outfilename=reducdir+filesbase[icube]+suffname2+'.fits'
 
 endif else begin
 writefits,reducdir+filesout[icube],0,h0
 writefits,reducdir+filesout[icube],a,h1,/append
 endelse

endfor
skiptoend:
close,/all
 
end 
