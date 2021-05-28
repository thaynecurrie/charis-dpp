pro charis_calc_raw_contrast,pfname,pick=pick,modsize=modsize,sub=sub,$
header=header,ap_factor=ap_factor,$
prefname=prefname,suffname=suffname,outputfile=outputfile,debug=debug

;****note, requires sat spots for now!!!!!!!

;v3.0 - Aug 19 2020 - Simplified and documented.   Still need to throw switch to pull sat flux information from fits headers
;v2.0 - Oct 4 2017 - Allows switch to compute the "contrast" as 5*halo brightness normalized by the satellite spots.  Sends output to diff directory
;v1.0 - Oct 1 2017 - reports contrast as rms of the noise (probably a bit dangerous since the residuals in a non-PSF subtracted, speckle filled image are non-gaussian)

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"charis_calc_raw_contrast:Raw contrast estimate (still under some development)"
print,"charis_calc_raw_contrast,pfname,pick=pick,r_ap=r_ap,datacube=datacube,calcube=calcube,modsize=modsize,nosub=nosub,"
print,"attenfac=attenfac,ap_factor=ap_factor,outputfile=outputfile"

print,""
print,"****Important Keywords"
print,""
print,"*pfname - info file name"
print,"*pick - select a single cube (otherwise, performs calculation for all images in info file"
print,"*modsize - manuallys set modulation amplitude of satellite spots"
print,"*sub - calculate rms vs. radius instead of raw intensity profile"
print,"*ap_factor - scale the aperture radius by some amount"
print,"*outputfile - outputfile for contrast curve"
print,"*header - use header value for the star flux [only available in combination with 'reg_cal' suffname]

goto,skiptotheend
endif

if keyword_set(sub) then begin
file_mkdir,'raw_contrasts_rms'
outputdir='raw_contrasts_rms'
endif else begin
file_mkdir,'raw_contrasts'
outputdir='raw_contrasts'
endelse

;applies a photometric calibration either to an entire sequence of files or to one or more identified image

;default is to apply to sequence of registered but not PSF subtracted images
;allow for override by choosing 'pick'

if keyword_set(sub) then nosubfac=1. else nosubfac=1/5.

;***Preliminary***
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


;aperture radius: use the FWHM of the image for now
;if you want to keep the aperture equal to the image slice FWHM then do nothing.  If you want to change the scaling set the ap_factor keyword

;*****Identify files to estimate contrasts****
if keyword_set(pick) then begin
files_full_path=dialog_pickfile(Title="Select Datacubes For Contrast Computation",/multiple)
nfiles=n_elements(files_full_path)
files=files_full_path
filesbase=files_full_path

;pick out the directory under which your files reside.  WARNING: assumes that all files are in same subdirectory!!!!!
your_path=strpos(files_full_path[0],'/',/reverse_search)+1

;and put the new files exactly where they were before; when you open the file later, look in the same place.
reducdir=strmid((files_full_path[0]),0,your_path)
;`datadir=reducdir
;print,'REDUCDIR ',reducdir

;if ~keyword_set(suffname) then suffname='reg'

;***trimming off path
for i=0L,nfiles-1 do begin
z=strsplit(files_full_path[i],'/',/extract)
files[i]=z[n_elements(z)-1]
suffposition=strpos(files[i],'.fits',/reverse_search)
filesbase[i]=strmid(files[i],0,suffposition)
endfor

suffname2='_cal'

goto,breakout


endif else begin

;assume we are working in the 'reg' subdirectory
param,'fnum_sat',flist,/get,pfname=pfname
filenum=nbrlist(flist)

if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='reg_cal'

;input files
suffname2=suffname+'_cal'
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)


;assume a simple overwrite
filesout=filelist(filenum,nfiles,prefix=prefname,suffix=suffname2)

endelse

breakout:

;print,files

if ~keyword_set(modsize) then begin
;if look for the ASTROGRID modulation amplitude in the fits header
test=readfits(reducdir+files[0],h1test,ext=1)
h0test=headfits(reducdir+files[0])

modsize1=sxpar(h0test,'X_GRDAMP',count=modcount1)
modsize2=sxpar(h1test,'X_GRDAMP',count=modcount2)

if (modcount1 eq 0 and modcount2 eq 0) then begin
read,"Manually Enter the ASTROGRID modulation amplitude in nanometers ",modsize
endif else begin
if modcount1 eq 0 then  modsize=modsize2*1d3
if modcount2 eq 0 then  modsize=modsize1*1d3
if (modcount1 ne 0 and modcount2 ne 0) then modsize=0.5*(modsize1+modsize2)*1d3
endelse

endif else begin
test=readfits(reducdir+files[0],h1test,ext=1)
h0test=headfits(reducdir+files[0])
modsize=modsize
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
attenfac=get_charis_sat_contrast(modsize,wvlh_test*1d-3,mjd_in)


;*********Contrast Measurements********
;*****Now perform photometry on the satspot_images
;***open the file, pull the fits header from the file, measure the integrated signal of the spots, do SNR map calc of image,$
;    get residual noise profile, divide by sat spots and multiply by 5, display, write contrast at 0.25", 0.5" and 0.9" to file


;test calimage
imcaltest=test
dimcal=(size(imcaltest,/dim))
profrad_tc,imcaltest[*,*,0],1,1,dimcal[0]/2,p1d=p1test,rayon=rayontest

;not doing the 'goodcode hex2bin stuff with GPI. assume all sat spots are good

satflux=fltarr(nfiles,dimcal[2])
esatflux=fltarr(nfiles,dimcal[2])

medsatflux=fltarr(nfiles)

;output contrasts
median_contrast=fltarr(nfiles,n_elements(p1test)) 
contrast_per_channel=fltarr(n_elements(p1test),dimcal[2],nfiles)
giantcube=fltarr(dimcal[0],dimcal[1],dimcal[2],nfiles)

noisemap_slicef=fltarr(n_elements(p1test),dimcal[2],nfiles)  ;hardwired for now, assume square arrays

for ical=0L,nfiles-1 do begin
print,'Calculating Raw Contrast for File ',files[ical]
h0cal=headfits(reducdir+files[ical])
imcal=readfits(reducdir+files[ical],h1cal,ext=1)

;***Get the CHARIS wavelengths***
get_charis_wvlh,h0cal,wvlh_test,filtname=filtname0
wvlh_charis=wvlh_test
lambda=1.d-3*wvlh_charis
Dtel=charis_get_constant(name='Dtel')
pixscale=charis_get_constant(name='pixscale')

fwhm=1.0*(1.d-6*(wvlh_charis*1d-3)/Dtel)*(180.*3600./!dpi)/pixscale

if keyword_set(verbose) then print,'fwhm is',fwhm

;Now nominally assume that the aperture radius is 1/2 the FWHM
aperrad=double(fwhm*0.5)
if keyword_set(verbose) then print,aperrad[0],fwhm[0]

if ~keyword_set(ap_factor) then ap_factor = 1.0

cens=fltarr(2,4,dimcal[2])

;loop on wavelength
 for s=0L,dimcal[2]-1 do begin

  for j=0,3 do begin  ;assume 4 sat spots
   tmp=fltarr(2)+!values.f_nan
    sf=string(s)
    tmp= sxpar(h1cal,'SATS'+strtrim(long(sf),2)+'_'+strtrim(j,2))
;**total, ugly hack-tastic workaround since I can't seem to get IDL to treat the string properly otherwise...
    tmpz=strsplit(tmp,' ',/extract)
    tmp=[tmpz[0],tmpz[1]]
    cens[*,j,s]=tmp
  endfor
 endfor

;if you have a calibrated cube, you can just use fits header for star photometry instead of calculating it all over again.


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
;do a radial profile subtraction of the image.  Should get close to zero'ing the background.
;assume for now that the satellite spot signal does not affect the median profile calculation.
imslice=imcal[*,*,il]
profrad_tc,imslice,1,1,dimcal[0]/2,p2d=pr
imslice-=pr

aperradf=aperrad[il]*ap_factor
if keyword_set(verbose) then print,aperradf,fwhm[il],fwhm[il]/2.,ap_factor,aperrad[il]*ap_factor,aperrad[il]

if (keyword_set(header) and suffname eq 'reg_cal') then begin
sf=string(il)
tmp=sxpar(h1cal,'FSTAR_'+strtrim(long(sf),2))
tmp2=sxpar(h1cal,'CERR'+strtrim(long(sf),2))*tmp
satflux[ical,il]=tmp
esatflux[ical,il]=tmp2

endif else begin
phpadu=1.0
sat_skyrad=[2.6]*aperradf

if keyword_set(verbose) then  print,cens[0,0,il],cens[1,0,il],imslice[cens[0,0,il],cens[1,0,il]]
  ;stop
aper, imslice,cens[0,0,il],cens[1,0,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent,setskyval=0
     sat1flux[il]=flux
 ;print,'flux and sky1 is',flux,sky
aper, imslice,cens[0,1,il],cens[1,1,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent,setskyval=0
     sat2flux[il]=flux
 ;print,'flux and sky2 is',flux,sky
aper, imslice,cens[0,2,il],cens[1,2,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent,setskyval=0
     sat3flux[il]=flux
 ;print,'flux and sky3 is',flux,sky
aper, imslice,cens[0,3,il],cens[1,3,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent,setskyval=0
     sat4flux[il]=flux
 ;print,'flux and sky4 is',flux,sky

       stdev_sat_sys=stdev([sat1flux[il],sat2flux[il],sat3flux[il],sat4flux[il]])
       stdev_sat_in=0.25*sqrt(sat1noise[il]^2.+sat2noise[il]^2+sat3noise[il]^2.+sat4noise[il]^2.)
       esatflux[ical,il]=sqrt(stdev_sat_in^2.+stdev_sat_sys^2.)

       satflux[ical,il]=mean([sat1flux[il],sat2flux[il],sat3flux[il],sat4flux[il]])/attenfac[il]
    
endelse
       imslicemask=imslice

;SNR calculation per slice to get the background rms within an aperture, perform this on images with the sat spots masked out.
       dist_circle,mask1,dimcal[0],[cens[0,0,il],cens[1,0,il]]
       dist_circle,mask2,dimcal[0],[cens[0,1,il],cens[1,1,il]]
       dist_circle,mask3,dimcal[0],[cens[0,2,il],cens[1,2,il]]
       dist_circle,mask4,dimcal[0],[cens[0,3,il],cens[1,3,il]]
       masksat=where(mask1 le aperrad or mask2 le aperrad or mask3 le aperrad or mask4 le aperrad,nmasksat)
       imslicemask[masksat]=!values.f_nan

       if keyword_set(sub) then begin
       charis_snratio_sub,imslicemask,fwhm=aperradf*2,/finite,noisemap=noisemap_slice,/silent,/zero,/filt,convimage=convimage
       profrad_tc,noisemap_slice,1,1,100,rayon=rayon,p1d=dum
       endif else begin
       charis_snratio_sub,abs(imcal[*,*,il]),fwhm=aperradf*2,/finite,noisemap=noisemap_slice,/silent,convimage=convimage
       profrad_tc,convimage,1,1,100,rayon=rayon,p1d=dum
       endelse

       ;writefits,'convimage.fits',convimage/(sat1flux[il]/attenfac[il])
       ;writefits,'imcal.fits',imcal[*,*,il]

       noisemap_slicef[*,il,ical]=dum

       sat1noise[il]=noisemap_slice(cens[0,0,il],cens[1,0,il])
       sat2noise[il]=noisemap_slice(cens[0,1,il],cens[1,2,il])
       sat3noise[il]=noisemap_slice(cens[0,2,il],cens[1,2,il])
       sat4noise[il]=noisemap_slice(cens[0,3,il],cens[1,3,il])




       if keyword_set(verbose) then print,satflux[ical,il],esatflux[ical,il],satflux[ical,il]/esatflux[ical,il]

giantcube[*,*,il,ical]=convimage

;print,sat1flux[il],sat2flux[il],sat3flux[il],sat4flux[il],satflux[ical,il],attenfac[il]
contrast_per_channel[*,il,ical]=nosubfac*5*noisemap_slicef[*,il,ical]/satflux[ical,il]

if keyword_set(debug) then begin

print,satflux[ical,il],sat1flux[il],sat2flux[il],sat3flux[il],sat4flux[il],attenfac[il],il,wvlh_charis[il]


;plot the contrast
set_plot,'x'
setcolors,/system_variables,/silent
 ;plot,rayon*pixscale,5*noisemap_slicef[*,il,ical],xrange=[0,1],yrange=[10,10000],/ylog,$
 plot,rayon*pixscale,5*contrast_per_channel[*,il,ical],xrange=[0,1],yrange=[5d-6,5d-3],/ylog,$
 xtitle='Separation (arc-sec)',ytitle=textoidl(' Contrast ( \Delta F)'),xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=1.5
 loadct,13,/silent

 ;for ilam=0L,dimcal[2]-1 do begin
 ;oplot,rayon*pixscale,nosubfac*5*noisemap_slicef[*,ilam,ical]/satflux[ical,ilam],color=(double(ilam+1)/dimcal[2])*256,thick=3
 ;endfor
 ;oplot,rayon*pixscale,nosubfac*5*noisemap_colf/medsatflux[ical],color=!cyan,thick=10

writefits,'noisemap.fits',noisemap_slicef[*,il,ical]
wait,1
endif

endfor

;calculate the contrast in the median_collapsed cube
;imcol=median(imcal,dimension=3,/even)

imcol=median(giantcube[*,*,*,ical],dimension=3,/even)
medsatflux[ical]=median(satflux[ical,*],/even)
print,medsatflux[ical]

writefits,'imcol.fits',imcol


if keyword_set(sub) then begin
;charis_snratio_sub,imcol,fwhm=median(aperrad*ap_factor,/even)*2.,/finite,noisemap=noisemap_col,/silent,/filt,/zero
profrad_tc,abs(imcol/0.6745),1,1,100,rayon=rayon,p1d=noisemap_colf
endif else begin
profrad_tc,imcol,1,1,100,rayon=rayon,p1d=noisemap_colf
endelse

;writefits,'imcol.fits',imcol
;writefits,'imcolf.fits',imcolf

;noise map of wavelength-collapsed image
;profrad_tc,noisemap_col,1,1,100,p1d=noisemap_colf,rayon=rayon

;debug
;writefits,'contrastmap1.fits',5*noisemap_col/medsatflux[ical]

;plot the contrast
if keyword_set(sub) then begin
plotrange=[5d-6,5d-3]
endif else begin
plotrange=[1d-5,2.5d-3]
endelse

;plot the contrast
set_plot,'x'
setcolors,/system_variables,/silent
 plot,rayon*pixscale,nosubfac*median(attenfac,/even)*5*noisemap_colf/medsatflux[ical],/nodata,xrange=[0,1],yrange=[plotrange[0],plotrange[1]],/ylog,$
 xtitle='Separation (arc-sec)',ytitle=textoidl(' Contrast ( \Delta F)'),xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=1.5
 loadct,13,/silent

 for ilam=0L,dimcal[2]-1 do begin
 oplot,rayon*pixscale,nosubfac*5*noisemap_slicef[*,ilam,ical]/satflux[ical,ilam],color=(double(ilam+1)/dimcal[2])*256,thick=3
 endfor
 oplot,rayon*pixscale,nosubfac*5*noisemap_slicef[*,4,ical]/satflux[ical,4],color=!blue,thick=5
oplot,rayon*pixscale,nosubfac*5*noisemap_slicef[*,11,ical]/satflux[ical,11],color=!green,thick=5
oplot,rayon*pixscale,nosubfac*5*noisemap_slicef[*,18,ical]/satflux[ical,18],color=!purple,thick=5
 oplot,rayon*pixscale,nosubfac*5*noisemap_colf/medsatflux[ical],color=!cyan,thick=10


set_plot,'ps'
device,filename='./'+outputdir+'/'+files[ical]+'_contrastcurve.eps',/encapsulated,/color,bits=8
setcolors,/system_variables,/silent
 plot,rayon*pixscale,nosubfac*median(attenfac,/even)*5*noisemap_colf/medsatflux[ical],/nodata,xrange=[0,1],yrange=[plotrange[0],plotrange[1]],/ylog,$
 xtitle='Separation (arc-sec)',ytitle=textoidl(' 5-\sigma Contrast ( \Delta F)'),xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=3
 loadct,13,/silent

;plot contrast in individual slices
 for ilam=0L,dimcal[2]-1,1 do begin
 if ilam eq 4 or ilam eq 11 or ilam eq 18 then begin
 oplot,rayon*pixscale,smooth(nosubfac*5*noisemap_slicef[*,ilam,ical]/satflux[ical,ilam],3),color=(double(ilam+1)/dimcal[2])*256,thick=8
 endif else begin
 oplot,rayon*pixscale,smooth(nosubfac*5*noisemap_slicef[*,ilam,ical]/satflux[ical,ilam],3),color=(double(ilam+1)/dimcal[2])*256,thick=3,$
linestyle=2
 endelse
 endfor

;plot wavelength-collapsed image contrast
setcolors,/system_variables,/silent
 oplot,rayon*pixscale,nosubfac*5*noisemap_colf/medsatflux[ical],color=!cyan,thick=10

loadct,13,/silent
setcolors,/system_variables,/silent
al_legend,color=[(5./22.)*256,(12/22.)*256,(19./22.)*256,!cyan],['J','H','Ks','Broadband'],/bottom,$
box=0,linestyle=0,thick=[3,3,3,10],charsize=1.25,charthick=3
device,/close

median_contrast[ical,*]=nosubfac*5*noisemap_colf/medsatflux[ical]

;debug
;print,'median contrast is ',median_contrast[ical,0:40]

;write to file: 1. r (arc-sec), 2. contrast in wavelength-collapsed cube, 3. J-slice contrast, 4. H-slice contrast, 5. K-slice contrast
writecol,'./'+outputdir+'/'+files[ical]+'_contrast.txt',rayon*pixscale,nosubfac*5*noisemap_colf/medsatflux[ical],$
nosubfac*5*noisemap_slicef[*,4,ical]/satflux[ical,4],nosubfac*5*noisemap_slicef[*,11,ical]/satflux[ical,11],nosubfac*5*noisemap_slicef[*,18,ical]/satflux[ical,18],fmt='(f,f,f,f,f)'

endfor

if nfiles gt 1 then begin
mediancube=median(giantcube,dimension=4,/even)  ;averaged cube
collapsedcube=median(mediancube,dimension=3,/even) ;averaged broadband image
endif else begin
mediancube=reform(giantcube,dimcal[0],dimcal[1],dimcal[2])  ;averaged cube
collapsedcube=median(mediancube,dimension=3,/even) ;averaged broadband image
endelse


medsatfluxcube=median(satflux,dimension=1,/even) ;averaged satellite flux/channel
medsatfluxall=median(medsatfluxcube,/even)  ;averaged satellite flux/channel averaged over channels


;now, loop one last time to get sequence-averaged contrast per channel and sequence-averaged, wavelength-collapsed contrast

;contrast_per_channel=fltarr(dimchan[0],dimchan[1],dimchan[2])
contrast_curve_per_channel=fltarr(dimcal[2],n_elements(p1test))

for i=0L,dimcal[2]-1 do begin
charis_snratio_sub,mediancube[*,*,i],fwhm=aperrad[i]*ap_factor*2,/finite,noisemap=noisemap_seq_per_channel,/silent,/filt,/zero
if ~keyword_set(nosub) then begin
profrad_tc,noisemap_seq_per_channel,1,1,100,p1d=dum,rayon=rdum
endif else begin
profrad_tc,mediancube[*,*,i],1,1,100,p1d=dum,rayon=rdum
endelse

contrast_curve_per_channel[i,*]=nosubfac*5*dum/medsatfluxcube[i]
endfor

;contrast curve per channel

charis_snratio_sub,collapsedcube,fwhm=median(aperrad*ap_factor,/even)*2.,/finite,noisemap=noisemap_seq_col,/zero,/filt
profrad_tc,noisemap_seq_col,1,1,100,p1d=dum,rayon=rdum
medcontrastcol=nosubfac*5*dum/medsatfluxall

;debug
;writefits,'colcube.fits',collapsedcube/medsatfluxall
;writefits,'noisemap2.fits',noisemap_seq_col
;print,'medcontrastcol',medcontrastcol[0:40]


;plot the contrast of sequence-combined image
set_plot,'ps'
if keyword_set(outputfile) then begin
prefplotname=outputfile
endif else begin
profplotname='combined_'
endelse
device,filename='./'+outputdir+'/'+prefplotname+'_contrastcurve.eps',/encapsulated,/color,bits=8
;device,filename='./raw_contrasts/'+prefplotname+'_contrastcurve.eps',/encapsulated,/color,bits=8
;device,filename='./raw_contrasts/n'+nbr2txt(ical,4)+'_contrastcurve.eps',/encapsulated,/color,bits=8
setcolors,/system_variables,/silent
 plot,rayon*pixscale,median(attenfac,/even)*medcontrastcol,/nodata,xrange=[0,1],yrange=[5d-6,5d-3],/ylog,$
 xtitle='Separation (arc-sec)',ytitle=textoidl('Sequence-Combined 5-\sigma Contrast ( \Delta F)'),xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=3
 loadct,13,/silent
 for ilam=0L,dimcal[2]-1,1 do begin
 if ilam eq 4 or ilam eq 11 or ilam eq 18 then begin
 oplot,rayon*pixscale,smooth(contrast_curve_per_channel[ilam,*],3),color=(double(ilam+1)/dimcal[2])*256,thick=8
 endif else begin
 oplot,rayon*pixscale,smooth(contrast_curve_per_channel[ilam,*],3),color=(double(ilam+1)/dimcal[2])*256,thick=3,linestyle=2
 endelse
 endfor
setcolors,/system_variables,/silent
 oplot,rayon*pixscale,medcontrastcol,color=!cyan,thick=10
loadct,13,/silent
setcolors,/system_variables,/silent
al_legend,color=[(5./22.)*256,(12/22.)*256,(19./22.)*256,!cyan],['J','H','Ks','Broadband'],/bottom,$
box=0,linestyle=0,thick=[3,3,3,10],charsize=1.25,charthick=3
device,/close
set_plot,'x'
plot,rayon*pixscale,medcontrastcol,/nodata,xrange=[0,1],yrange=[5d-6,5d-3],/ylog,$
xtitle='Separation (arc-sec)',ytitle='Contrast (deltaF)',xthick=5,ythick=5,xstyle=1,ystyle=1
setcolors,/system_variables,/silent
for ilam=0L,dimcal[2]-1,1 do begin
 if ilam eq 4 or ilam eq 11 or ilam eq 18 then begin
 oplot,rayon*pixscale,smooth(contrast_curve_per_channel[ilam,*],3),color=(double(ilam+1)/dimcal[2])*256,thick=8
 endif else begin
 oplot,rayon*pixscale,smooth(contrast_curve_per_channel[ilam,*],3),color=(double(ilam+1)/dimcal[2])*256,thick=3
 endelse
endfor

oplot,rayon*pixscale,medcontrastcol,linestyle=0,thick=5,color=!cyan

if ~keyword_set(outputfile) then outputfile='contrast_combined_broadband_jhk'
if keyword_set(nosub) then outputfile='nosub'+outputfile
writecol,'./'+outputdir+'/'+outputfile,rayon*pixscale,medcontrastcol,contrast_curve_per_channel[4,*],$
contrast_curve_per_channel[11,*],contrast_curve_per_channel[18,*]

skiptotheend:

end
