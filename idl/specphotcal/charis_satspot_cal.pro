pro charis_satspot_cal,pfname,calmethod=calmethod,pick=pick,calcube=calcube,subskyannulus=subskyannulus,$
siglim=siglim,$
nopradcal=nopradcal,$
meancomb=meancomb,$
ap_factor=ap_factor,$
filtername=fname,$
prefname=prefname,suffname=suffname,$
verbose=verbose,$
;outfilename=outfilename,$
help=help

;Calibrates the contrast between the sat spots and the PSF
; - 3/25/2020 - version 0.1

if (N_PARAMS() eq 0 or keyword_set(help)) then begin

Print,"charis_satspot_cal,pfname,calmethod=calmethod,pick=pick,calcube=calcube,subskyannulus=subskyannulus"
print,"siglim=siglim,nopradcal=nopradcal,meancomb=meancomb,ap_factor=ap_factor,filtername=fname,"
print,"verbose=verbose,outfilename=outfilename"

Print,'Calibrate Sat Spot Contrast subroutine (useful for a ladder sequence)'
Print,'This program calibrates the contrast of the satellite spots'
Print,'It requires a data cube with the spots and unsaturated PSF of the star'
print,""
print,"***Keywords***"
print,'*pfname - parameter file (e.g. HR8799_low.info)'
print,"(very similar to charis_specphot_cal: consult documentation on that one: [e.g. charis_specphot_cal,/help])"
goto,skiptoend
endif

;***STEP 1: Preliminary Stuff***
;define data directory
reducdir='./reduc/'

;data directory
datadir=reducdir+'reg/'
subdir='reg/'
reducdir+=subdir

newdir='./ladder/'
file_mkdir,newdir

if ~keyword_set(siglim) then siglim=5.

;if you want to keep the aperture equal to the image slice FWHM then do nothing.  If you want to change the scaling set the ap_factor keyword
if ~keyword_set(ap_factor) then ap_factor = 1

;********STEP 2: CALIBRATION CUBE *********
if ~keyword_set(calmethod) then begin
calmethod=0
endif else begin
calmethod = long(calmethod)
endelse
print,'calmethod',calmethod
help,calmethod

case calmethod of 
 0: begin

if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='reg'
    param,'fnum_lad',flist,/get,pfname=pfname
    filenum=nbrlist(flist)
    calcube=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
imcaltest=readfits(datadir+calcube[0],h1test,/exten)
h0test=headfits(datadir+calcube[0])

    end

;the median average of registered cubes
 1: begin

if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='reg'
    param,'fnum_lad',flist,/get,pfname=pfname
    filenum=nbrlist(flist)
    calcube=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
imcaltest=readfits(datadir+calcube[0],h1test,/exten)
h0test=headfits(datadir+calcube[0])

    end
endcase

ncal=n_elements(calcube)
dimcal=(size(imcaltest,/dim))

;First, get the CHARIS wavelengths
get_charis_wvlh,h0test,wvlh_test


;This (below) is what you are trying to determine)
;Now from the modulation size, get the CHARIS attenuation factor/channel
;attenfac=get_charis_sat_contrast(modsize,wvlh_test*1d-3,mjd_in)

;******STEP 3: Determine Flux Calibration Scaling ******;

;*****Now perform photometry on the calimage, record its integration time
;***open the file, pull the fits header from the file, measure the integrated signal, measure the noise, do snrcalc/wavelength slice.
ncal=n_elements(calcube)
;test calimage

print,calcube

;not doing the 'goodcode hex2bin stuff with GPI. assume all sat spots are good

imgiantcube=fltarr(dimcal[0],dimcal[1],dimcal[2],ncal)

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

h0cal=headfits(datadir+calcube[ical])
imcal=readfits(datadir+calcube[ical],h1cal,ext=1)

for il=0L,dimcal[2]-1 do begin
imslice=imcal[*,*,il]

;for now, do not mask sat spots
if ~keyword_set(nopradcal) then begin
bad=where(imslice eq 0 or finite(imslice) eq 0,complement=good)
imslice2=imslice
imslice2[bad]=!values.f_nan
profrad_tc,imslice2,1,2*fwhm[il],dimcal[0]/2,p2d=pr
imslice[good]-=pr[good]
endif

imgiantcube[*,*,il,ical]=imslice

endfor

endfor

;****Create master calibration cube if ther are more than one initial cal cubes
if keyword_set(verbose) then print,calmethod

case calmethod of
 1:begin
if ncal gt 1 then begin

if ~keyword_set(meancomb) then begin
imcalcube=median(imgiantcube,dimension=4,/even)
endif else begin
imcalcube=mean(imgiantcube,dimension=4,/nan,/double)
endelse

endif else begin

imcalcube=reform(imgiantcube,dimcal[0],dimcal[1],dimcal[2])
endelse

ncal=1
   end
 else:begin
      end

endcase

satflux=fltarr(ncal,dimcal[2])
esatflux=fltarr(ncal,dimcal[2])
starflux=fltarr(ncal,dimcal[2])


;***Get the CHARIS wavelengths***
if ~keyword_set(filtername) then begin
get_charis_wvlh,h0test,wvlh_charis,filtname=filtname0
endif else begin
filtname0=fname
endelse

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
aperradf=aperrad[il]*ap_factor

;aperrad=0.5*fwhm[il]*ap_factor
;aperture photometry at each slice assuming a zero background
phpadu=1.0
sat_skyrad=[2,6]*aperradf

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

aper,imslice,dimcal[0]/2,dimcal[1]/2,sflux,esflux,ssky,sskyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
     starflux[ical,il]=sflux

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

aper,imslice,dimcal[0]/2,dimcal[1]/2,sflux,esflux,ssky,sskyerr,phpadu,aperradf,setskyval=0,[0,0],/flux,/exact,/nan,/silent
     starflux[ical,il]=sflux

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


  ;      stddev_sat_flux=  fltarr(n_elements(cens[0,0,*]))

       stdev_sat_in=0.25*sqrt(sat1noise[il]^2.+sat2noise[il]^2+sat3noise[il]^2.+sat4noise[il]^2.)
       stdev_sat_sys=stdev([sat1flux[il],sat2flux[il],sat3flux[il],sat4flux[il]])

       esatflux[ical,il]=sqrt(stdev_sat_in^2.+stdev_sat_sys^2.)

       satflux[ical,il]=mean([sat1flux[il],sat2flux[il],sat3flux[il],sat4flux[il]])
       ;if keyword_set(verbose) then print,il,'SatFlux ESatFlux Ratio',satflux[il],esatflux[il],satflux[il]/esatflux[il]
endfor
endfor

if keyword_set(verbose) then print,ncal
if ncal gt 1 then begin

satfluxf=reform(median(satflux,/even,dimension=1))
esatfluxf=reform(median(esatflux,/even,dimension=1))
starfluxf=reform(median(starflux,/even,dimension=1))
endif else begin
satfluxf=reform(satflux)
esatfluxf=reform(esatflux)
starfluxf=reform(starflux)
endelse

if keyword_set(verbose) then begin
window,3
ploterror,wvlh_charis*1d-3,satfluxf/starfluxf,esatfluxf/starfluxf
oplot,wvlh_charis*1d-3,satflux[0,*]/starflux[0,*],linestyle=1
oplot,wvlh_charis*1d-3,satflux[1,*]/starflux[1,*],linestyle=2
oplot,wvlh_charis*1d-3,satflux[2,*]/starflux[2,*],linestyle=3
endif

snratio_sats=satfluxf/esatfluxf

if keyword_set(verbose) then begin
for i=0L,21 do print,i,esatfluxf[i],satfluxf[i],starfluxf[i],esatfluxf[i]/starfluxf[i]
endif

if keyword_set(verbose) then begin
window,1
ploterror,wvlh_charis*1.d-3,satfluxf/starfluxf,esatfluxf/starfluxf,xrange=[1,3],yrange=[0,max(satfluxf/starfluxf)]
oplot,wvlh_charis*1.d-3,satfluxf[10]/starfluxf[10]*(wvlh_charis/wvlh_charis[10])^(-2.),linestyle=2
endif

if keyword_set(verbose) then print,'satfluxf ', satfluxf

p0=(satfluxf/starfluxf)[9]  ;this is close to 1.55 microns
good=where(esatfluxf/satfluxf le 1./siglim) ;choose only n-sigma sat spots
output=mpfitexpr('P[0]*(X/1.55)^(-2.)',wvlh_charis[good]*1d-3,(satfluxf/starfluxf)[good],(esatfluxf/starfluxf)[good],p0)
atten_fac0=output[0]


if keyword_set(verbose) then begin
window,4
conpred=output[0]*(wvlh_charis*1d-3/1.55)^(-2.)
ploterror,wvlh_charis*1.d-3,satfluxf/starfluxf,esatfluxf/starfluxf,xrange=[1,2.5],yrange=[0,max(satfluxf/starfluxf)]
oplot,wvlh_charis*1.d-3,conpred,linestyle=1
endif
;skiptoend:
close,/all


;Now, write the output of the scaling to file called "atten_fac_ladder"
newfilename='atten_fac_ladder'
openw,1,newdir+newfilename
printf,1,';AttenFac From charis_satspot_cal (1.55 microns)'
printf,1,atten_fac0
close,1

skiptoend:
end 
