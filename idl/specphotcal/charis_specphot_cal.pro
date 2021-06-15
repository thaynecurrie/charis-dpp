pro charis_specphot_cal,pfname,calmethod=calmethod,pick=pick,datacube=datacube,calcube=calcube,modamp=modamp,set_atten_fac=atten_fac0,subskyannulus=subskyannulus,$
leftspot=leftspot,rightspot=rightspot,lowerleftspot=lowerleftspot,lowerrightspot=lowerrightspot,upperleftspot=upperleftspot,upperrightspot=upperrightspot,$
fitbckgd=fitbckgd,$
;simpleerror=simpleerror,$
expandederror=expandederror,$
filtcal_slice=filtcal_slice,nopradcal=nopradcal,$
notimeratio=notimeratio,$
meancomb=meancomb,$
fixradius=fixradius,$
ap_factor=ap_factor,$
;attenfac=attenfac,ap_factor=ap_factor,$
filtername=fname,$
;starname=starname,mag=mag,$
starlib=starlib,$
empspectrum=empspectrum,$
av=av,$
prefname=prefname,suffname=suffname,$
verbose=verbose,$
fluxunits=fluxunits,outfilename=outfilename,help=help
;,guide=guide

;Spectro-Photometric Calibration program 
; - Note: a bit uncertain about whether we need to scale the integration time or not.  I.e. are counts "per read" or are they cumulative?
; - appears to be counts/read, so answer is no
; - 3/21/2021 - adds "fitbckgd" switch to do a plane-fit subtraction of the background surrounding the spots for spectrophotometry. and option to choose left or right spots
; - 1/19/2021 - sets lab-determined dispersion for sig-sat spot to default; fixes the ladder code, splitting cases up into two
; - 9/18/2020 - Changed the default spectral library to Kurucz since we almost never want to use Pickles.
; - 3/26/2020 - Allows 0) cal on indiv cubes, 1) an average cube, 2) another selected cube, 3) a "ladder" cube with sat spot contrast saved in ladder directory
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
Print,' using data of a known star taken with a similar setup'
Print,'  '
Print,'****Calling Sequence ****'
Print,'  '
Print,'charis_specphot_cal,pfname,calmethod=calmethod,pick=pick,datacube=datacube,calcube=calcube,modamp=modamp,set_atten_fac=atten_fac0,subskyannulus=subskyannulus,fitbckgd=fitbckgd,expandederror=expandederror,leftspot=leftspot,rightspot=rightspot,'
;Print,'charis_specphot_cal,pfname,calmethod=calmethod,pick=pick,datacube=datacube,calcube=calcube,modamp=modamp,subskyannulus=subskyannulus,simpleerror=simpleerror,'
Print,'filtcal_slice=filtcal_slice,nopradcal=nopradcal,meancomb=meancomb,ap_factor=ap_factor,filtername=fname,starname=starname,mag=mag,starlib=starlib,'
Print,'av=av,prefname=prefname,suffname=suffname,test=test,verbose=verbose,fluxunits=fluxunits,outfilename=outfilename,help=help,guide=guide'
Print,' '
Print,"Example: charis_specphot_cal,'HR8799_low.info',calmethod=0"
Print,' '
Print,'***Keywords** ...'
Print,'pfname - input info file for applying correction to pre-determined sequence of registered images'
Print,'calmethod - spectrophotometric calibration method: 0) cal on indiv cubes, 1) an average cube, 2) another selected cube, 3) a ladder sequence to adjust (0), 4) a ladder cube, 5) calibrated cube of another star'
Print,'pick - manually select images to be calibrated'
Print,'datacube - Data to Calibrate (if you define manually)'
Print,'calcube - Datacube used for Flux Calibration (if you select manually)'
Print,'modamp - to manually set the amplitude modulation for the satellite spots'
Print,'set_atten_fac - manually set the attenuation factor at 1.55 microns, overriding pipeline calculation'
Print,'fitbckgd - do a plane-fit subtraction of the background surrounding the spots for spectrophotometry.'
Print,'left/right spot - option to choose just either the left or right spots.'
Print,'subskyannulus - subtract the local "sky" background for sat spot photometry (usually no)'
Print,'expandederror - ignore internal source tests used to estimate systematic error, take the sigma of individual sat spot measurements'
;Print,'simpleerror - set systematic error/channel to 1% for median-comb or 2% for mean-comb [justified from internal source tests]'
Print,'filtcal_slice - spatially filter (subtract median filter) before sat spot photometry (usually no)'
Print,'nopradcal - do not subtract off the median-average signal as a function of r (not advised)'
Print,'fixradius - fix the aperture radius to a constant instead of scaling with wavelength
Print,'meancomb - mean-combine the calibration cubes instead of median (usually not recommended)'
Print,'ap_factor - rescale the aperture for photometry'
Print,'filtername - set manual value for the filter (low, JHK)'
Print,'starlib - 1,2,3 (Kurucz,Pickles, or Manual) [Kurucz is recommended]
Print,'pref/suffname - Manually override the prefix and suffix of file names used (useful for the pfname switch but applied to images other than "reg" images)'
print,'fluxunits - What Flux Density units? 0=mJy,1=Jy,2=W/m^2/um,3=ergs/s/cm^2/A,4=ergs/s/cm^2/Hz,5=contrast'
;print,''
;print,'Limitations ...'
goto,skiptoend
endif

;Applies a photometric calibration either to an entire sequence of files or to one or more identified image

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
;0 - Kurucz 1993 (theoretical models, pretty robust)
;1 - Pickles 1998 PASP library (empirical, incomplete data)
;2 - empirical model (user-selected)

if ~keyword_set(starlib) then begin
read,"Select the Stellar Library Source (1 = Kurucz, 2 = Pickles, 3 = empirical model) ",slib
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
;registered cubes individually
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
;an average cube
if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='reg'

    param,'fnum_sat',flist,/get,pfname=pfname
    filenum=nbrlist(flist)
    calcube=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
imcaltest=readfits(datadir+calcube[0],h1test,/exten)
h0test=headfits(datadir+calcube[0])

    end

;user-id'd cubes
 2: begin

    if ~keyword_set(calcube) then $
calcube=dialog_pickfile(Title="Select Your Datacube Used for Flux Calibration",/multiple_files)
;***for now just allow for one cube. change this later!!!
imcaltest=readfits(calcube[0],h1test,/exten)
h0test=headfits(calcube[0])

    end

;calibration determined from ladder sequence 
 3: begin

if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='reg'

    param,'fnum_sat',flist,/get,pfname=pfname
    filenum=nbrlist(flist)
    calcube=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
imcaltest=readfits(datadir+calcube[0],h1test,/exten)
h0test=headfits(datadir+calcube[0])

   end
 
;ladder/separate sat spot-on cubes
 4: begin
if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='reg'
    param,'fnum_lad',flist,/get,pfname=pfname
    filenum=nbrlist(flist)
    calcube=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
    imcaltest=readfits(datadir+calcube[0],h1test,/exten)
    h0test=headfits(datadir+calcube[0])

    end

;user-id'd cubes from another star
 5: begin

    if ~keyword_set(calcube) then $
calcube=dialog_pickfile(Title="Select Your Datacube Used for Flux Calibration (Another Star, Use FITS headers)")
;***for now just allow for one cube. change this later!!!
imcaltest=readfits(calcube[0],h1test,/exten)
h0test=headfits(calcube[0])

    end

endcase

ncal=n_elements(calcube)
dimcal=(size(imcaltest,/dim))

;*****astrogrid spot calibration*****

if calmethod eq 5 then begin

get_charis_wvlh,h0test,wvlh_test

goto,skipastrogridnumber
endif

if calmethod ne 3 then begin
;if look for the ASTROGRID modulation amplitude in the fits header

if keyword_set(modamp) then begin
modsize=modamp
print,'modsize is',modsize

endif else begin
modsize1=sxpar(h0test,'X_GRDAMP',count=modcount1)
modsize2=sxpar(h1test,'X_GRDAMP',count=modcount2)
if ((modcount1 eq 0 and modcount2 eq 0) or (modsize1 eq 'unavailable' and modsize2 eq 'unavailable')) then begin
read,"Manually Enter the ASTROGRID modulation amplitude in nanometers ",modsize
endif else begin
;go find the modulation size in the fits header
; if in the primary but not secondary use primary
if modcount1 eq 0 then  modsize=modsize2*1d3
; if in the secondary but not primary use secondary
if modcount2 eq 0 then  modsize=modsize1*1d3
; if in both, then use primary
if modcount1 ne 0 and modcount2 ne 0 then modsize=modsize1*1d3
endelse
print,'modsize is',modsize

endelse

endif else begin

ladderdir='./ladder/'
ladderfile='atten_fac_ladder'
readcol,ladderdir+ladderfile,atten_fac0
endelse

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
if (calmethod ne 3 and ~keyword_set(atten_fac0)) then begin
attenfac=get_charis_sat_contrast(modsize,wvlh_test*1d-3,mjd_in)
endif else begin
attenfac=get_charis_sat_contrast(25,wvlh_test*1d-3,mjd_in,manual=atten_fac0)
endelse

;***If you are just going to use the fits headers for another star, then skip all of this.
skipastrogridnumber:

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
5: begin
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

for il=0L,dimcal[2]-1 do begin
imslice=imcal[*,*,il]

;for now, do not mask sat spots
if ~keyword_set(nopradcal) then begin
bad=where(imslice eq 0 or finite(imslice) eq 0,complement=good)
imslice2=imslice
imslice2[bad]=!values.f_nan
profrad_tc,imslice2,1,1,dimcal[0]/2,p2d=pr
imslice[good]-=pr[good]
endif

if keyword_set(filtcal_slice) then begin
imslice-=filter_image(imslice,median=10*fwhm[il])
endif

imgiantcube[*,*,il,ical]=imslice

endfor

endfor

;*****IF YOU HAVE THE CONVERSION PRE-DEFINED, then do this, skip to end
if calmethod eq 5 then begin

conv_fact_tot=fltarr(ncal,dimcal[2])
cerr0=fltarr(ncal,dimcal[2])
fstar0=fltarr(ncal,dimcal[2])
for ical=0L,ncal-1 do begin
 for il=0L,dimcal[2]-1 do begin
  conv_fact_tot[ical,il]=sxpar(h1cal,'FSCALE'+strc(il))
  cerr0[ical,il]=sxpar(h1cal,'CERR'+strc(il))
  fstar0[ical,il]=sxpar(h1cal,'FSTAR_'+strtrim(il,2))
  print,fstar0[ical,il],cerr0[ical,il],ical,il
 endfor
 funitname_output=sxpar(h1cal,'FLUXUNIT')
endfor

if ncal gt 1 then begin
conv_fact0=median(conv_fact_tot,/even,dimension=1)
cerr=median(cerro0,/even,dimension=1)
speccal=median(fstar0,/even,dimension=1)
snratio_sats=1./cerr
endif else begin
conv_fact0=conv_fact_tot
cerr=reform(cerr0,dimcal[2])
;fstar=fstar0
speccal=fstar0
snratio_sats=1./cerr
endelse

conv_fact0=replicate(1,nfiles)#conv_fact0
print,cerr
goto,skiptofluxcalingcubes
endif

;;*****IF YOU DON"T HAVE THE CONVERSION PRE-DEFINED (calmethod /=5)
;****Create master calibration cube if there is more than one initial cal cube

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

ncal=1

   end

 3:begin

;assume that the integration time is constant!
; tint_calf=tint_cal[0] 

   end

4:begin

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
    

;****Get the Satellite Spot Positions
;*****NOTE: Assumes the Satellite Spot Positions are the Same in 
cens=fltarr(2,4,dimcal[2])

;loop on wavelength
 for s=0L,dimcal[2]-1 do begin
  for j=0,3 do begin  ;assume 4 sat spots
   tmp=fltarr(2)+!values.f_nan
    sf=string(s)
    tmp= sxpar(h1test,'SATS'+strtrim(long(sf),2)+'_'+strtrim(j,2))  
;**total, ugly hack-tastic workaround since I can't seem to get IDL to treat the string properly otherwise...
    tmpz=strsplit(tmp,' ',/extract)
    tmp=[tmpz[0],tmpz[1]]
    cens[*,j,s]=tmp
  endfor
 endfor

for ical=0L,ncal-1 do begin

case calmethod of 
  0: begin
imcalcube=imgiantcube[*,*,*,ical] 
     end
  3: begin
imcalcube=imgiantcube[*,*,*,ical] 
     end
  else:begin
       end
endcase


;****Initialize the sat spot fluxes
;**now you have the positions of the satellite spots.  time to extract the photometry.
;assume that you don't have any errors triggered like can happen with the GPI pipeline

; extract the flux of the satellite spots
  sat1flux = fltarr(n_elements(cens[0,0,*]))    ;;top left
  sat2flux = fltarr(n_elements(cens[0,0,*]))    ;;bottom left
  sat3flux = fltarr(n_elements(cens[0,0,*]))    ;;top right
  sat4flux = fltarr(n_elements(cens[0,0,*]))    ;;bottom right
  sat1noise = fltarr(n_elements(cens[0,0,*]))    ;;top left
  sat2noise = fltarr(n_elements(cens[0,0,*]))    ;;bottom left
  sat3noise = fltarr(n_elements(cens[0,0,*]))    ;;top right
  sat4noise = fltarr(n_elements(cens[0,0,*]))    ;;bottom right

for il=0L,dimcal[2]-1 do begin
imslice=imcalcube[*,*,il]

if ~keyword_set(fixradius) then begin
aperradf=aperrad[il]*ap_factor
endif else begin
aperrad[il]=fixradius
aperradf=fixradius
endelse

;if keyword_set(verbose) then print,aperradf,fwhm[il],fwhm[il]/2.,ap_factor,aperrad[il]*ap_factor,aperrad[il]

;aperture photometry at each slice assuming a zero background
phpadu=1.0
sat_skyrad=[2,6]*aperradf

;if keyword_set(verbose) then  print,cens[0,0,il],cens[1,0,il],imslice[cens[0,0,il],cens[1,0,il]]


;plane-fit the background? -- usually not needed but offered for better precision.
if keyword_set(fitbckgd) then begin
dist_circle,gdist,201
weights=0
for q=0L,3 do begin
sourceind=get_xycind(201,201,cens[0,q,il],cens[1,q,il],5*aperradf)
bckind0=get_xycind(201,201,cens[0,q,il],cens[1,q,il],15*aperradf)
sourcedist=sqrt((cens[0,q,il]-201/2)^2.+(cens[1,q,il]-201/2)^2.)
annulus=where(abs(gdist-sourcedist) lt 5*aperradf)
bckind0=intersect(bckind0,annulus)
background=intersect(bckind0,sourceind,/xor)
coeff=planefit(background mod 201, background / 201,imslice[background],weights,yfit)
;xinds1=sourceind mod 201 
;yinds1 = sourceind / 201
;xinds2=background mod 201
;yinds2=background / 201
xinds=bckind0 mod 201
yinds=bckind0 / 201
src_bkg_plane=coeff[0]+coeff[1]*xinds+coeff[2]*yinds
imslice[bckind0]-=src_bkg_plane
endfor
endif

;0-top left, 1-top-right, 2-bottom-right, 3-bottom-left

if keyword_set(subskyannulus) then begin
aper, imslice,cens[0,0,il],cens[1,0,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
     sat1flux[il]=flux
 ;print,'flux and sky1 is',flux,sky
aper, imslice,cens[0,1,il],cens[1,1,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
     sat2flux[il]=flux
 ;print,'flux and sky2 is',flux,sky
aper, imslice,cens[0,2,il],cens[1,2,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
     sat3flux[il]=flux
 ;print,'flux and sky3 is',flux,sky
aper, imslice,cens[0,3,il],cens[1,3,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
     sat4flux[il]=flux
endif else begin

aper, imslice,cens[0,0,il],cens[1,0,il],flux,eflux,sky,skyerr,phpadu,aperradf,setskyval=0,[0,0],/flux,/exact,/nan,/silent
     sat1flux[il]=flux
 ;print,'flux and sky1 is',flux,sky
aper, imslice,cens[0,1,il],cens[1,1,il],flux,eflux,sky,skyerr,phpadu,aperradf,setskyval=0,[0,0],/flux,/exact,/nan,/silent
     sat2flux[il]=flux
 ;print,'flux and sky2 is',flux,sky
aper, imslice,cens[0,2,il],cens[1,2,il],flux,eflux,sky,skyerr,phpadu,aperradf,setskyval=0,[0,0],/flux,/exact,/nan,/silent
     sat3flux[il]=flux
 ;print,'flux and sky3 is',flux,sky
aper, imslice,cens[0,3,il],cens[1,3,il],flux,eflux,sky,skyerr,phpadu,aperradf,setskyval=0,[0,0],/flux,/exact,/nan,/silent
     sat4flux[il]=flux


endelse

       imslicemask=imslice
;SNR calculation per slice to get the background rms within an aperture, perform this on images with the sat spots masked out.
       dist_circle,mask1,dimcal[0],[cens[0,0,il],cens[1,0,il]]
       dist_circle,mask2,dimcal[0],[cens[0,1,il],cens[1,1,il]]
       dist_circle,mask3,dimcal[0],[cens[0,2,il],cens[1,2,il]]
       dist_circle,mask4,dimcal[0],[cens[0,3,il],cens[1,3,il]]
       masksat=where(mask1 le aperradf or mask2 le aperradf or mask3 le aperradf or mask4 le aperradf,nmasksat)
       imslicemask[masksat]=!values.f_nan
       charis_snratio_sub,imslicemask,fwhm=aperradf*2,/finite,noisemap=noisemap_slice,/silent

       sat1noise[il]=noisemap_slice(cens[0,0,il],cens[1,0,il])
       sat2noise[il]=noisemap_slice(cens[0,1,il],cens[1,2,il])
       sat3noise[il]=noisemap_slice(cens[0,2,il],cens[1,2,il])
       sat4noise[il]=noisemap_slice(cens[0,3,il],cens[1,3,il])

;Satellite Spot Flux
;mean or median-combine sat fluxs to get scaling? 
       if keyword_set(meancomb) then begin
       ;satflux[ical,il]=mean([sat1flux[il],sat2flux[il],sat3flux[il],sat4flux[il]])

       resistant_mean,([sat1flux[il],sat2flux[il],sat3flux[il],sat4flux[il]]),3,bah
       satflux[ical,il]=bah

       endif else begin
       satflux[ical,il]=median(([sat1flux[il],sat2flux[il],sat3flux[il],sat4flux[il]]),/even)
       endelse

;Satellite Spot Flux Uncertainty

;Average/Median Empirical Noise Value [should be the same for all]
       if keyword_set(meancomb) then begin
       ;stdev_sat_in=mean([sat1noise[il],sat2noise[il],sat3noise[il],sat4noise[il]])

       resistant_mean,[sat1noise[il],sat2noise[il],sat3noise[il],sat4noise[il]],3,stdev_sat_in
       endif else begin
       stdev_sat_in=median([sat1noise[il],sat2noise[il],sat3noise[il],sat4noise[il]],/even)
       endelse

       ;stdev_sat_in=0.25*sqrt(sat1noise[il]^2.+sat2noise[il]^2+sat3noise[il]^2.+sat4noise[il]^2.)

;****override to select left or right spots
       if keyword_set(leftspot) then begin
       satflux[ical,il]=0.5*(sat1flux[il]+sat4flux[il])
       stdev_sat_in=stdev([sat1flux[il],sat4flux[il]])
       endif
       if keyword_set(rightspot) then begin
       satflux[ical,il]=0.5*(sat2flux[il]+sat3flux[il])
       stdev_sat_in=stdev([sat2flux[il],sat3flux[il]])
       endif

       if keyword_set(lowerleftspot) then begin
       satflux[ical,il]=sat4flux[il]
       stdev_sat_in=stdev([sat1flux[il],sat4flux[il]])
       endif
       if keyword_set(lowerrightspot) then begin
       satflux[ical,il]=sat3flux[il]
       stdev_sat_in=stdev([sat2flux[il],sat3flux[il]])
       endif

       if keyword_set(upperleftspot) then begin
       satflux[ical,il]=sat1flux[il]
       stdev_sat_in=stdev([sat1flux[il],sat4flux[il]])
       endif
       if keyword_set(upperrightspot) then begin
       satflux[ical,il]=sat2flux[il]
       stdev_sat_in=stdev([sat2flux[il],sat3flux[il]])
       endif

       if ~keyword_set(expandederror) then begin
        if keyword_set(meancomb) then begin 
            stdev_sat_sys=0.02*satflux[ical,il]
            endif else begin
           stdev_sat_sys=0.01*satflux[ical,il]
        endelse
       endif else begin
       stdev_sat_sys=robust_sigma([sat1flux[il],sat2flux[il],sat3flux[il],sat4flux[il]])
       endelse

       esatflux[ical,il]=sqrt(stdev_sat_in^2.+stdev_sat_sys^2.)


       if keyword_set(verbose) then print,'ical il SatFlux SatFlux/ESatFlux SatFlux/stdevSat',ical,il,satflux[ical,il],satflux[ical,il]/esatflux[ical,il],satflux[ical,il]/stdev_sat_in,format='(a,i,i,f,f,f)'
       if keyword_set(verbose) then print,'ical il Sat Fluxes are ',sat1flux[il],sat2flux[il],sat3flux[il],sat4flux[il],format='(a,f,f,f,f)'
       if keyword_set(verbose) then print,'ical il esats',esatflux[ical,il],stdev_sat_in,stdev_sat_sys,format='(a,f,f,f)'
       ;if keyword_set(verbose) then print,'ical il SatFlux SatFlux/ESatFlux ',ical,il,satflux[ical,il],esatflux[ical,il],satflux[ical,il]/esatflux[ical,il],format='(a,i,i,f,f,f)'
endfor
       print,'For Cube ',ical,' Flux and EFlux for Channel 1',satflux[ical,0],esatflux[ical,0]
endfor

if keyword_set(verbose) then print,satflux[*,19]/esatflux[*,19]

;if n_elements(ncal) gt 1 then begin
if ncal gt 1 then begin
satfluxf=reform(median(satflux,/even,dimension=1))
esatfluxf=reform(median(esatflux,/even,dimension=1))
endif else begin
satfluxf=reform(satflux)
esatfluxf=reform(esatflux)
endelse

;the median SNR of the sat spots over the entire sequence.  
;this is ok if you want to get the abs. cal. uncertainty for a sequence-combined cube
;make sure that you are doing PSF subtraction over the same cubes that you do for spec-phot calibration.
snratio_sats=satfluxf/esatfluxf

;print,''
;print,'SNR for sats'
if keyword_set(verbose) then print,snratio_sats
;help,esatflux,esatfluxf
;help,satflux,satfluxf
;help,snratio_sats
;stop

if keyword_set(verbose) then begin
window,1
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

if keyword_set(verbose) then print,'satfluxf ', satfluxf
print,'speccal',speccal

skiptofluxcalingcubes:
conv_fact=fltarr(nfiles,dimcal[2])

case calmethod of
0: begin
satfluxf=satflux
   end
3: begin
satfluxf=satflux
   end
5: begin
   end
else: begin
satfluxf=replicate(1,nfiles)#satfluxf

      end
endcase



;;;flux-cal-ing the data cubes

for icube=0L,nfiles-1 do begin

case calmethod of 
5: begin
conv_fact[icube,*]=conv_fact0[icube,*]
end
else: begin
conv_fact[icube,*]=(1.0/satfluxf[icube,*])*speccal*attenfac
end
endcase

 a=readfits(datadir+files[icube],ext=1,h1)
 h0=headfits(datadir+files[icube])

texp=sxpar(h0,'exp1time')
ncoadd=sxpar(h0,'coadds')
tint_cal_indiv=texp*float(ncoadd)


 print,'the file name is',files[icube]
 print,'the data directory is ',datadir

 for ichannel=0L,dimcal[2]-1 do begin
 ;for ichannel=0L,n_elements(wvlh_charis)-1 do begin
 a[*,*,ichannel]*=conv_fact[icube,ichannel]
 endfor

time_ratio=1.0 ; I think CHARIS computes in counts/s
case calmethod of 
 1: begin
;if the integration time for this particular cube is different than your test one ...
;time_ratio=tint_calf/tint_cal_indiv
if keyword_set(notimeratio) then time_ratio=1.0
a[*,*,*]*=time_ratio
    end
 2: begin
;time_ratio=tint_calf/tint_cal_indiv
if keyword_set(notimeratio) then time_ratio=1.0
a[*,*,*]*=time_ratio
    end
 3: begin
;time_ratio=tint_calf/tint_cal_indiv
if keyword_set(notimeratio) then time_ratio=1.0
a[*,*,*]*=time_ratio

    end
 else: begin
 time_ratio=1.0
     end
endcase


;print,time_ratio

 ;funitname_output=strtrim(funitname_output)

sxaddpar,h1,'FLUXUNIT',strtrim(funitname_output,1),"Flux Density Units in Cube"

;Flux Density of Star
for ichannel=0L,dimcal[2]-1 $
;for ichannel=0L,n_elements(wvlh_charis)-1 $
 do sxaddpar,h1,'FSTAR_'+strtrim(ichannel,2),speccal[ichannel]," Star Flux"

;aperture
for ichannel=0L,dimcal[2]-1 $
;for ichannel=0L,n_elements(wvlh_charis)-1 $
 do sxaddpar,h1,'R_AP'+strtrim(ichannel,2),aperrad[ichannel]*ap_factor,"Aperture Radius (Pixels)"

;sky annulus
; do sxaddpar,h1,'R_AP'+strtrim(ichannel,2),aperrad[ichannel]*ap_factor,"Aperture Radius (Pixels)"
sxaddpar,h1,'Sky_In',strtrim(strc(2)),"Inner Radius for Sky Annulus (in Units of Aperture Radius)"
sxaddpar,h1,'Sky_Out',strtrim(strc(6)),"Outer Radius for Sky Annulus (in Units of Aperture Radius)"


;scaling
for ichannel=0L,dimcal[2]-1 $
;for ichannel=0L,n_elements(wvlh_charis)-1 $
 do sxaddpar,h1,'FSCALE'+strc(ichannel),time_ratio*conv_fact[icube,ichannel],"scale to convert counts to "+strc(funitname_output)

;for ichannel=0L,n_elements(wvlh_charis)-1  do begin
;for ichannel=0L,n_elements(wvlh_charis)-1  do begin
;print,'convfact2',conv_fact0[icube,ichannel]
;endfor
;stop

;error in spectrophotometric calibration/slice
for ichannel=0L,dimcal[2]-1 $
;for ichannel=0L,n_elements(wvlh_charis)-1 $
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
