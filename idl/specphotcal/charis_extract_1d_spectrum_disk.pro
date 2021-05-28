pro charis_extract_1d_spectrum_disk,pfname,datacube=datacube,coords=coords,fitcollapse=fitcollapse,fitslice=fitslice,planefit=planefit,nozerosky=nozerosky,noprad=noprad,filt_slice=filt_slice,perpixel=perpixel,fixaprad=fixaprad,$
pick=pick,r_ap=r_ap,calcube=calcube,$
attenfac=attenfac,ap_factor=ap_factor,$
filter=fname,starname=starname,mag=mag,$
ndfilt=ndfilt,$
clipoutliers=clipoutliers,$
throughputcor=throughputcor,$
prefname=prefname,suffname=suffname,test=test,$
fluxunits=fluxunits,help=help,outputspec=outputspec

;Built off of GPI's spectral extraction program but probably more robust.
;***modified to do surface brightness calculations

;***Preliminary***
;define data directory
reducdir='./reduc/'

;data directory
datadir=reducdir+'reg/'
;subdir='reg/'
subdir='proc/'
reducdir+=subdir

;information about the target properties.
param,'HMAG',hmag,/get,pfname=pfname
param,'EHMAG',ehmag,/get,pfname=pfname
param,'SPTYPE',spt,/get,pfname=pfname

if keyword_set(datacube) then begin
extractcube=datacube
extractcube=reducdir+datacube
endif else begin
extractcube=dialog_pickfile(Title="Select the Data Cube From Which You Want to Extract a Spectrum")
endelse

;read in data cube, if the cube has been flux-calibrated then proceed, if not then do flux calibration

test=readfits(extractcube,h1,ext=1)

caltest=sxpar(h1,'FLUXUNIT',count=ct)

;If 'FLUXUNIT' has been found, then your data cube is flux-calibrated. proceed with extracting a spectrum.  If it has not been found, then do flux calibration first.


if ct eq 0 then begin    
;********porting charis_specphot_cal.pro here*****
;*************************************************

charis_specphot_cal,pfname,pick=pick,r_ap=r_ap,datacube=extractcube,calcube=calcube,$
filter=fname,starname=starname,mag=mag,$
prefname=prefname,suffname=suffname,test=test,$
fluxunits=fluxunits,help=help,outfilename=outfilename
print,'outfilename is',outfilename

extractube=outfilename
print,'extractcube is ..',extractcube
;the outfile is now the file from which you will extract a cube

endif 


;****Extract a Spectrum*********
;*******************************

;1. Define position
if ~keyword_set(coords) then begin
read,"Input X Coordinate ",xpos
read,"Input Y Coordinate ",ypos
endif else begin
xpos=coords[0]  ;x position
ypos=coords[1]  ;y position
endelse
xposf=xpos
yposf=ypos

;2. Read in the datacube
data=readfits(extractcube,h1,ext=1)
h0=headfits(extractcube)
get_charis_wvlh,h0,wvlh_charis ;pull wavelengths
;Telescope Diameter for Subaru
Dtel=charis_get_constant(name='Dtel') ;7.9d0 ;visible pupil for SCExAO
;pixel scale 16.4 mas/pixel
pixscale=charis_get_constant(name='pixscale') ;0.0164 nominally
fwhm=1.0*(1.d-6*(wvlh_charis*1d-3)/Dtel)*(180.*3600./!dpi)/pixscale

nlambda=(size(data,/dim))[2]
aperrad=fltarr(nlambda)

;3. Define extraction radius
if ~keyword_set(fixaprad) then begin
 for ir=0L,nlambda-1 do begin
  aperrad[ir]=sxpar(h1,'R_AP'+strtrim(string(ir),2))
 endfor
  ;print,aperrad
  ;stop
endif else begin
aperrad=0.5*fwhm
endelse

;4.If you want to fit the position of the feature in the collapsed cube, then throw the fitcollapse switch
;****warning: be ultra-careful with this for disk extraction. Will probably throw an error since extended emission is not going to be fit by a 
; gaussian PSF
if keyword_set(fitcollapse) then begin
data_col=median(data,dimension=3,/even)
med_wvlh=median(wvlh_charis,/even)
fwhm_col=1.0*(1.d-6*(med_wvlh*1d-3)/Dtel)*(180.*3600./!dpi)/pixscale
gcntrd,data_col,xpos,ypos,xposf,yposf,fwhm_col

endif


fluxdens_spectrum=fltarr(nlambda)
efluxdens_spectrum=fltarr(nlambda)
snrslice=fltarr(nlambda)

for il=0L,nlambda-1 do begin
fwhm_slice=fwhm[il]
aperradf=aperrad[il]
if keyword_set(ap_factor) then aperradf*=ap_factor
sat_skyrad=aperradf*[2,4]

if keyword_set(filt_slice) then begin
noisyslice=data[*,*,il]
noisyslice-=filter_image(noisyslice,median=10*fwhm_slice)
data[*,*,il]=noisyslice
endif

if ~keyword_set(noprad) then begin
profrad_tc,data[*,*,il],p2d=pr
data[*,*,il]-=pr
endif


;if you want to fit the companion position per wavelength slice then throw the fitslice switch
;*** warning: again, be careful about this with disk extraction!
if keyword_set(fitslice) then begin
gcntrd,data[*,*,il],xpos,ypos,xposf,yposf,fwhm_slice
print,'Fitted Centroid Position is ...',xposf,yposf,' ... for slice ',strtrim(il,2)
endif

;Aperture Photometry
phpadu=1.0 ;dummy number
if keyword_set(nozerosky) then begin
aper,data[*,*,il],xposf,yposf,flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
endif else begin
aper,data[*,*,il],xposf,yposf,flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent,setskyval=0.
endelse

;Now pull the spot calibration uncertainty from the fits header
spotcal_uncert=sxpar(h1,'CERR'+strtrim(long(il),2))

;now compute the SNR properly
;charis_snratio_sub,data[*,*,il],fwhm=2*aperrad[il],coord=[xposf,yposf],snrval=snrval,/finite,/zero,/silent,snrmap=snrmap,/fixpos
charis_snratio_sub,data[*,*,il],fwhm=2*aperradf,coord=[xposf,yposf],snrval=snrval,/finite,/zero,/silent,snrmap=snrmap,/fixpos
if il eq 11 then begin
writefits,'snrmap.fits',snrmap
endif

snrslice[il]=snrmap[xposf,yposf]

if ~keyword_set(perpixel) then begin
ff=filter_image(data[*,*,il],median=2*aperradf)
endif else begin
ff=data[*,*,il]
endelse
flux=ff[xposf,yposf]
;now convert to surface brightness: native units/square arcsec
fluxdens_spectrum[il]=flux/(pixscale^2.)
;fluxdens_spectrum[il]=flux/(!pi*aperradf^2.*pixscale^2.)
;efluxdens_spectrum[il]=(fluxdens_spectrum[il]/snrval)
efluxdens_spectrum[il]=sqrt((fluxdens_spectrum[il]/snrval)^2.+(flux*spotcal_uncert/pixscale)^2.)

print,'ff is ',fluxdens_spectrum[il],data[xposf,yposf,il],flux,il
data[*,*,il]=ff/pixscale^2.

endfor

if (keyword_set(throughputcor)) then begin
if (throughputcor ne 0) then begin
;if (keyword_set(throughputcor) or throughputcor eq 1) then begin
;1. search for synth_throughput*txt files
;2. if you find no files, stop program and tell user to run fwdmod first
;3. if you only find one file, then use this one
;4. if you find more than one, then ask user which one to use
test=file_search('synth_throughput*txt')
ntest=n_elements(test)
print,test
case ntest of
   0:begin
     print,'ERROR! You need to run forward-modeling to get a throughput correction estimate first!'
     stop
     end
   1:begin
     throughput_est=file_search('synth_throughput*txt')
     readcol,throughput_est,dumlam,throughput
     end
   else: begin
     throughput_est=dialog_pickfile(Title="Select Your Throughput Correction Map",filter='*.txt',/fix_filter)
     readcol,throughput_est,dumlam,throughput
     end
endcase

fluxdens_spectrum/=throughput
efluxdens_spectrum/=throughput
endif
endif


if ~keyword_set(outputspec) then outputspec='spectrum_disk.dat'
writecol,outputspec,wvlh_charis*1d-3,fluxdens_spectrum,efluxdens_spectrum,snrslice
writefits,'sbmap.fits',data,h1

if ~keyword_set(clipoutliers) then begin
set_plot,'ps'
device,filename='spectrum_disk.eps',bits=8,/encapsulated
plot,wvlh_charis*1d-3,fluxdens_spectrum,xrange=[1,2.5],/nodata,xthick=5,ythick=5,xtitle='Wavelength (Microns)',ytitle='Flux Density (mJy)/Square Arcsec'
oploterror,wvlh_charis*1d-3,fluxdens_spectrum,efluxdens_spectrum,errstyle=1,linestyle=0,thick=4
device,/close

set_plot,'x'
plot,wvlh_charis*1d-3,fluxdens_spectrum,xrange=[1,2.5],/nodata,xthick=5,ythick=5,xtitle='Wavelength (Microns)',ytitle='Flux Density (mJy)/Square Arcsec'
oploterror,wvlh_charis*1d-3,fluxdens_spectrum,efluxdens_spectrum,errstyle=1,linestyle=0,thick=4

endif else begin

outliers=where(abs(fluxdens_spectrum) gt 10*median(fluxdens_spectrum,/even) or throughput lt 0.03, complement=good)
set_plot,'ps'
device,filename='spectrum_disk.eps',bits=8,/encapsulated
plot,wvlh_charis[good]*1d-3,fluxdens_spectrum[good],xrange=[1,2.5],/nodata,xthick=5,ythick=5,xtitle='Wavelength (Microns)',ytitle='Flux Density (mJy)/Square Arcsec'
oploterror,wvlh_charis[good]*1d-3,fluxdens_spectrum[good],efluxdens_spectrum[good],errstyle=1,linestyle=0,thick=4
device,/close

set_plot,'x'
plot,wvlh_charis[good]*1d-3,fluxdens_spectrum[good],xrange=[1,2.5],/nodata,xthick=5,ythick=5,xtitle='Wavelength (Microns)',ytitle='Flux Density (mJy)/Square Arcsec'
oploterror,wvlh_charis[good]*1d-3,fluxdens_spectrum[good],efluxdens_spectrum[good],errstyle=1,linestyle=0,thick=4


endelse

end
