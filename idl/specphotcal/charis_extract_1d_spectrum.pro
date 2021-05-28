pro charis_extract_1d_spectrum,pfname,datacube=datacube,coords=coords,fitcollapse=fitcollapse,fitslice=fitslice,fitbckgd=fitbckgd,extended=extended,$
throughputcor=throughputcor,$
fitrange=fitrange,$
nozerosky=nozerosky,noprad=noprad,filt_slice=filt_slice,filtsnr=filtsnr,nocalerror=nocalerror,$
bklim=bklim,bklir=bklir,$
pick=pick,r_ap=r_ap,calcube=calcube,$
;attenfac=attenfac,ap_factor=ap_factor,$ currently not set
filter=fname,starname=starname,mag=mag,$
oplotstar=oplotstar,$
clipoutliers=clipoutliers,plotfromzero=plotfromzero,yclip=yclip,$
prefname=prefname,suffname=suffname,test=test,$
verbose=verbose,$
breakout=breakout,$
xposout=xposout,yposout=yposout,snrout=snrout,$
fluxunits=fluxunits,outputspec=outputspec,help=help

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"Spectral Extraction, Outputs: a) file with spectrum, b) SNR map of collapsed cube"
print,""
print,"charis_extract_1d_spectrum,pfname,datacube=datacube,coords=coords,fitcollapse=fitcollapse,fitslice=fitslice,fitbckgd=fitbckgd,extended=extended,fitrange=fitrange,throughputcor=throughputcor,"
print,"nozerosky=nozerosky,noprad=noprad,filt_slice=filt_slice,filtsnr=filtsnr,nocalerror=nocalerror"
print,"bklim=bklim,bklir=bklir,pick=pick,r_ap=r_ap,calcube=calcube,filter=fname,starname=starname,mag=mag,"
print,"oplotstar=oplotstar,clipoutliers=clipoutliers,plotfromzero=plotfromzero,yclip,yclip,breakout=breakout,xposout=xposout,yposout=yposout,snrout=snrout,fluxunits=fluxunits,outputspec=outputspec"
print,""
print,"*****Example: charis_extract_1d_spectrum,'HD1160_low.info',datacube='loci.fits',/throughputcor,/filtsnr,/fitcollapse"
print,""
print,"***Major Keywords***"
print,'*pfname - parameter file (e.g. HR8799_low.info)'
print,"*datacube - the data cube from which you want to extract a spectrum"
print,"*coords - the coordinates of the companion [can be approximate if you use 'fitcollapse']"
print,"*fitcollapse - do you fit for the centroid from the wavelength-collapsed image?"
print,"*fitslice - do you fit individual slices instead?"
print,"*fitrange - fit centroid over this range of wavelength channels [e.g. fitrange=[2:10] fits over the 3-11th channels]"
print,"*fitbckgd - do you fit and subtract the background [useful mainly if the source is sitting on disk emission]"
print,"*extended - is the source extended?"
print,"*throughputcor - Apply a throughput correction? [from forward-modeling, usually mandatory for a calibrated spectrum]"
print,"*nozerosky - subtract a sky annulus [you rarely want to do this as bckgd is usually flattened]"
print,"*noprad - do not subtract off a radial profile [typically doesn't matter, can matter in pathological reductions]"
print,"*filt_slice - spatially filter the wavelength slice before extracting spectrum [usually no]"
print,"*filtsnr - apply a spatial filter before computing SNR [usually advantageous]"
print,"*nocalerror - ignore the `errors' on the sat spot photometry [charis_specphot_cal can often overestimate them]"
print,"*bklim, bklir - the azimuthal and radial extents for fitbckgd"
print,"*pick - pick the cube from GUI"
print,""
print,"**Plotting/Output Keywords"
print,""
print,"*oplotstar - overplot a scaled stellar spectrum"
print,"*clipoutliers - clip outlying flux measurements when plotting spectrum [useful for noisy spectra]"
print,"*plotfromzero - self-explanatory"
print,"*yclip - clip out flux values   yclip-times the median"
print,"*breakout - through switch if the derived position is not determinable [avoids command line errors]"
print,"*xposout, yposout,snrout - x,y, and SNR output for companion"
print,"*outputspec - change the name of the output spectrum"

goto,skiptotheend
endif

;Built off of GPI's spectral extraction program but probably more robust.

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

if ~keyword_set(bklim) then bklim=4
if ~keyword_set(bklir) then bklir=2
bklim = bklim > 2.999999  ;ensure that the background region is greater than 4*r
bklir = bklir > .99999

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

goto,skipoverthisnow
if ct eq 0 then begin    
;********porting charis_specphot_cal.pro here*****
;*************************************************

charis_specphot_cal,pfname,pick=pick,r_ap=r_ap,datacube=extractcube,calcube=calcube,$
;attenfac=attenfac,ap_factor=ap_factor,$ ;not in use
filter=fname,starname=starname,mag=mag,$
prefname=prefname,suffname=suffname,test=test,$
fluxunits=fluxunits,help=help,outfilename=outfilename
print,'outfilename is',outfilename

extractube=outfilename
print,'extractcube is ..',extractcube

;the outfile is now the file from which you will extract a cube

endif 
skipoverthisnow:

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
 for ir=0L,nlambda-1 do begin
  aperrad[ir]=sxpar(h1,'R_AP'+strtrim(string(ir),2))
 endfor
  ;print,aperrad

;4.If you want to fit the position of the companion in the collapsed cube, then throw the fitcollapse switch

data_col=median(data,dimension=3,/even)
med_wvlh=median(wvlh_charis,/even)
fwhm_col=1.0*(1.d-6*(med_wvlh*1d-3)/Dtel)*(180.*3600./!dpi)/pixscale

if ~keyword_set(fitrange) then fitrange=[0:n_elements(aperrad)-1]

if keyword_set(fitcollapse) then begin
apmed=median(aperrad[fitrange],/even)

if keyword_set(fitbckgd) then begin

if ~keyword_set(fitrange) then begin
image=data_col
endif else begin
image=median(data[*,*,fitrange],dimension=3,/even)
data_col=image
endelse

writefits,'imagefits.fits',image

;to start, get indices of the ``source", defined as < 2*rad
sourceind=get_xycind(201,201,xposf,yposf,2*apmed)
bckind0=get_xycind(201,201,xposf,yposf,bklim*apmed)
sourcedist=sqrt((xposf-201/2)^2.+(yposf-201/2)^2.)
dist_circle,gdist,201
annulus=where(gdist gt sourcedist-bklir*apmed and gdist lt sourcedist+bklir*apmed)  ;an annulus +/- 1.5 aperture radii
bckind0=intersect(bckind0,annulus)
background=intersect(bckind0,sourceind,/xor)
weights=0
coeff=planefit(background mod 201, background / 201,image[background],weights,yfit)

;plane-fit's estimate of the background underneath the source
xinds=sourceind mod 201
yinds=sourceind / 201
src_bkg_plane=coeff[0]+coeff[1]*xinds+coeff[2]*yinds
image[sourceind]-=src_bkg_plane
gcntrd,image,xpos,ypos,xposf,yposf,fwhm_col
if (xposf lt 0 or yposf lt 0) then begin
cntrd,image,xpos,ypos,xposf,yposf,fwhm_col
endif


if keyword_set(extended) then begin
;if the source is extended
sourceind=get_xycind(201,201,xposf,yposf,2*apmed)
xinds=sourceind mod 201
yinds=sourceind / 201
;center of mass
xcen=total(1.*xinds*image[sourceind])/(total(image[sourceind])*1.)
ycen=total(1.*yinds*image[sourceind])/(total(image[sourceind])*1.)
print,'xposf',xposf,yposf
xposf=xcen
yposf=ycen
endif

endif else begin

cntrd,data_col,xpos,ypos,xposf,yposf,fwhm_col

if (xposf lt 0 or yposf lt 0) then begin
gcntrd,data_col,xpos,ypos,xposf,yposf,fwhm_col
endif

if keyword_set(extended) then begin
;if the source is extended
sourceind=get_xycind(201,201,xposf,yposf,2*apmed)
xinds=sourceind mod 201
yinds=sourceind / 201
;center of mass
xcen=total(1.*xinds*data_col[sourceind])/(total(data_col[sourceind])*1.)
ycen=total(1.*yinds*data_col[sourceind])/(total(data_col[sourceind])*1.)
xposf=xcen
yposf=ycen
endif


print,'centroid is ',xposf,yposf
endelse


;gcntrd,data_col,xpos,ypos,xposf,yposf,fwhm_col
print,'Fitted Centroid Position is ...',xposf,yposf
endif


if keyword_set(breakout) then begin
if finite(xposf) eq 0 or finite(yposf) eq 0 then begin
xposout=!values.f_nan
yposout=!values.f_nan
snrout=!values.f_nan
goto,skiptoend
endif

endif

fluxdens_spectrum=fltarr(nlambda)
efluxdens_spectrum=fltarr(nlambda)
snrslice=fltarr(nlambda)

for il=0L,nlambda-1 do begin
fwhm_slice=fwhm[il]
aperradf=aperrad[il]
sat_skyrad=aperradf*[2,bklim]

;if you want to fit the companion position per wavelength slice then throw the fitslice switch
if keyword_set(fitslice) then begin
gcntrd,data[*,*,il],xpos,ypos,xposf,yposf,fwhm_slice
print,'Fitted Centroid Position is ...',xposf,yposf,' ... for slice ',strtrim(il,2)
endif

;how to treat the background ...
;1.nominal: do nothing except subtract radial profile: assume background is zero'd by lstsquares psf sub (usually valid)
;2.fitbckgd: fit the background to an annular region excluding the companion, interpolate over companion's position
;3.apply a spatial filter of 10 psf footprints
;4.do absolutely nothing: usually valid, usually indistinguishable from #1

if ~keyword_set(noprad) then begin
profrad_tc,data[*,*,il],p2d=pr
data[*,*,il]-=pr
endif

if keyword_set(fitbckgd) then begin
image=data[*,*,il]
;to start, get indices of the ``source", defined as < 2*rad
sourceind=get_xycind(201,201,xposf,yposf,2*aperradf)
bckind0=get_xycind(201,201,xposf,yposf,sat_skyrad[1])
sourcedist=sqrt((xposf-201/2)^2.+(yposf-201/2)^2.)
dist_circle,gdist,201
annulus=where(gdist gt sourcedist-bklir*aperradf and gdist lt sourcedist+bklir*aperradf)  ;an annulus +/- 1.5 aperture radii
bckind0=intersect(bckind0,annulus)
background=intersect(bckind0,sourceind,/xor)
weights=0
coeff=planefit(background mod 201, background / 201,image[background],weights,yfit)

;plane-fit's estimate of the background underneath the source
xinds=sourceind mod 201
yinds=sourceind / 201
src_bkg_plane=coeff[0]+coeff[1]*xinds+coeff[2]*yinds

;plane-fit's estimate of the background  around the source
xwing=background mod 201
ywing=background / 201
wings=coeff[0]+coeff[1]*xwing+coeff[2]*ywing

bcksub=fltarr(201,201)
bcksub[sourceind]=src_bkg_plane

if keyword_set(verbose) then begin
blah=fltarr(201,201)
blah2=fltarr(201,201)
blah[sourceind]=src_bkg_plane
blah2[background]=wings
writefits,'blah.fits',blah
writefits,'blah2.fits',blah2
writefits,'total.fits',data[*,*,il]-(blah+blah2)

print,'background-subbing slice number ',il+1
if keyword_set(verbose) then $
writefits,'fitted_'+strtrim(il+1,2)+'.fits',data[*,*,il]-(blah+blah2)
endif

endif

if keyword_set(filt_slice) then begin
noisyslice=data[*,*,il]
noisyslice-=filter_image(noisyslice,median=10*fwhm_slice)
data[*,*,il]=noisyslice

if keyword_set(verbose) then $
writefits,'filtered_'+strtrim(il+1,2)+'.fits',noisyslice
endif




;Aperture Photometry
phpadu=1.0 ;dummy number

if keyword_set(nozerosky) then begin
aper,data[*,*,il],xposf,yposf,flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
endif else begin

if ~keyword_set(fitbckgd) then begin
aper,data[*,*,il],xposf,yposf,flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent,setskyval=0.
endif else begin
aper,(data[*,*,il]-bcksub),xposf,yposf,flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent,setskyval=0.

endelse


endelse

;Now pull the spot calibration uncertainty from the fits header
spotcal_uncert=sxpar(h1,'CERR'+strtrim(long(il),2))

;now compute the SNR per res element
;***major new assumption***
;1. if a point source, the commented out and new versions are the same
;2. if extended source, it covers more than one res element and you use a larger ap-radius
;   what you are really doing is summing over mult res elements.  
;   the calc then implicitly assumes that the SNR in the 'core' is a good proxy for the SNR 'out of' the core. 
;   e.g. that the intensity distribution is flat and these are independent samples.   
;   The alternative, summing over larger ap, overestimates finite element correction
if ~keyword_set(filtsnr) then begin
;charis_snratio_sub,data[*,*,il],fwhm=2*aperradf,coord=[xposf,yposf],snrval=snrval,/finite,zero=1,/silent,fixpos=1
charis_snratio_sub,data[*,*,il],fwhm=fwhm[il],coord=[xposf,yposf],snrval=snrval,/finite,zero=1,/silent,fixpos=1
endif else begin
;charis_snratio_sub,data[*,*,il],fwhm=2*aperradf,coord=[xposf,yposf],snrval=snrval,/finite,zero=1,/silent,/fixpos,filter=1,expfilt=3
charis_snratio_sub,data[*,*,il],fwhm=fwhm[il],coord=[xposf,yposf],snrval=snrval,/finite,zero=1,/silent,/fixpos,filter=1,expfilt=3
endelse

fluxdens_spectrum[il]=flux
if ~keyword_set(nocalerror) then begin
efluxdens_spectrum[il]=sqrt((fluxdens_spectrum[il]/snrval)^2.+(flux*spotcal_uncert)^2.)
endif else begin
efluxdens_spectrum[il]=sqrt((fluxdens_spectrum[il]/snrval)^2.+(0*flux*spotcal_uncert)^2.)
endelse
snrslice[il]=snrval

endfor

;Now get SNR of collapsed cube
if ~keyword_set(filtsnr) then begin
charis_snratio_sub,data_col,fwhm=fwhm_col,coord=[xposf,yposf],snrval=snrval,/finite,zero=1,/silent,fixpos=1,snrmap=snrmapcol
endif else begin
charis_snratio_sub,data_col,fwhm=fwhm_col,coord=[xposf,yposf],snrval=snrval,/finite,zero=1,/silent,fixpos=1,filter=1,expfilt=3,snrmap=snrmapcol
endelse

writefits,'snrmapcol.fits',snrmapcol


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
     readcol,throughput_est[0],dumlam,throughput ;to make this backwards-compatible with old readcol versions
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


maxval=max(fluxdens_spectrum) 
minval=min(fluxdens_spectrum)

if ~keyword_set(yclip) then yclip = 5
if keyword_set(clipoutliers) then begin
;print,maxval,yclip*median(fluxdens_spectrum,/even)
maxval = maxval < yclip*median(fluxdens_spectrum,/even)
goodclip=where(fluxdens_spectrum le yclip*median(fluxdens_spectrum,/even))
maxval=max(fluxdens_spectrum[goodclip])

endif

;decision tree to decide the minimum and maximum values
;broadband
if min(wvlh_charis) lt 1170 and max(wvlh_charis) gt 2300 then begin
minxrange=1.0
maxxrange=2.5
endif

;J
if min(wvlh_charis) lt 1170 and max(wvlh_charis) lt 1400 then begin
minxrange=1.1
maxxrange=1.4
endif

;H
if min(wvlh_charis) gt 1400 and min(wvlh_charis) lt 1600 then begin
minxrange=1.4
maxxrange=1.9
endif

;K
if min(wvlh_charis) gt 1850 then begin
minxrange=1.9
maxxrange=2.5
endif


;3. Define extraction radius
if keyword_set(oplotstar) then begin
starflux=fltarr(nlambda)
 for ir=0L,nlambda-1 do begin
  starflux[ir]=sxpar(h1,'FSTAR_'+strtrim(string(ir),2))
 endfor
starscalefac=1*median(starflux,/even)/median(fluxdens_spectrum,/even)
maxval=max(fluxdens_spectrum) > max(starflux)/starscalefac
minval=min(fluxdens_spectrum) < min(starflux)/starscalefac
endif

if keyword_set(plotfromzero) then minval=0
if ~keyword_set(outputspec) then outputspec='spectrum.dat'
writecol,outputspec,wvlh_charis*1d-3,fluxdens_spectrum,efluxdens_spectrum,snrslice

if ~keyword_set(yclip) then yclip =5

if ~keyword_set(clipoutliers) then begin
set_plot,'ps'
device,filename='spectrum.eps',bits=8,/encapsulated,/color
plot,wvlh_charis*1d-3,fluxdens_spectrum,xrange=[minxrange,maxxrange],yrange=[0.8*minval,1.1*maxval],/nodata,$
charsize=1.5,xthick=5,ythick=5,xtitle='Wavelength (Microns)',ytitle='Flux Density (mJy)',charthick=5
setcolors,/system_variables
oploterror,wvlh_charis*1d-3,fluxdens_spectrum,efluxdens_spectrum,errstyle=0,errcolor=!gray,linestyle=0,thick=10,errthick=5
if keyword_set(oplotstar) then oplot,wvlh_charis*1d-3,starflux/starscalefac,linestyle=2,color=!magenta,thick=10
device,/close

set_plot,'x'
;plot,wvlh_charis[good]*1d-3,fluxdens_spectrum[good],xrange=[1,2.5],yrange=[0.8*minval,1.1*maxval],/nodata,$
plot,wvlh_charis*1d-3,fluxdens_spectrum,xrange=[minxrange,maxxrange],yrange=[0.8*minval,1.1*maxval],/nodata,xthick=5,ythick=5,xtitle='Wavelength (Microns)',ytitle='Flux Density (mJy)',$
charsize=1.5,charthick=2,xstyle=1

oploterror,wvlh_charis*1d-3,fluxdens_spectrum,efluxdens_spectrum,errstyle=1,linestyle=0,thick=4
if keyword_set(oplotstar) then oplot,wvlh_charis*1d-3,starflux/starscalefac,linestyle=1

endif else begin

if ~keyword_set(throughput) then begin
outliers=where(abs(fluxdens_spectrum) gt yclip*median(fluxdens_spectrum,/even), complement=good)
endif else begin
outliers=where(abs(fluxdens_spectrum) gt yclip*median(fluxdens_spectrum,/even) or throughput lt 0.03, complement=good)
endelse

;print,abs(fluxdens_spectrum),5*median(fluxdens_spectrum,/even)

set_plot,'ps'
device,filename='spectrum.eps',bits=8,/encapsulated,/color
;plot,wvlh_charis[good]*1d-3,fluxdens_spectrum[good],xrange=[1,2.5],yrange=[0.8*minval,1.1*maxval],/nodata,$
plot,wvlh_charis[good]*1d-3,fluxdens_spectrum[good],xrange=[minxrange,maxxrange],yrange=[0.8*minval,1.1*maxval],/nodata,xthick=5,ythick=5,xtitle='Wavelength (Microns)',ytitle='Flux Density (mJy)',$
charsize=1.5,charthick=5,xstyle=1
oploterror,wvlh_charis[good]*1d-3,fluxdens_spectrum[good],efluxdens_spectrum[good],errstyle=1,linestyle=0,thick=10
setcolors,/system_variables
if keyword_set(oplotstar) then oplot,wvlh_charis*1d-3,starflux/starscalefac,linestyle=2,color=!magenta,thick=10
device,/close

set_plot,'x'
;plot,wvlh_charis[good]*1d-3,fluxdens_spectrum[good],xrange=[1,2.5],yrange=[0.8*minval,1.1*maxval],/nodata,$
plot,wvlh_charis[good]*1d-3,fluxdens_spectrum[good],xrange=[minxrange,maxxrange],yrange=[0.8*minval,1.1*maxval],/nodata,xthick=5,ythick=5,xtitle='Wavelength (Microns)',ytitle='Flux Density (mJy)',$
charsize=1.5,charthick=2,xstyle=1
oploterror,wvlh_charis[good]*1d-3,fluxdens_spectrum[good],efluxdens_spectrum[good],errstyle=1,linestyle=0,thick=4
if keyword_set(oplotstar) then oplot,wvlh_charis*1d-3,starflux/starscalefac,linestyle=1


endelse

if keyword_set(oplotstar) then print,'scale factor is ',starscalefac

errcen=(1/2.35)*fwhm_col/snrval

print,'Done !',' Spectrum Extracted at Position ',xposf,yposf,' [E,N] ',([xposf,yposf]-201/2)*0.0162,' with rough errors of ',errcen*0.0162
print,'the estimated SNR of the collapsed cube is ',snrval

xposout=(xposf-201/2)*0.0162
yposout=(yposf-201/2)*0.0162
snrout=snrval

skiptoend:

close,/all

skiptotheend:
end
