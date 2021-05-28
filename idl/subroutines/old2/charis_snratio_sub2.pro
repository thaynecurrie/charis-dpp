pro charis_snratio_sub1, image,coord=coord,fwhm=fwhm,rmax=rmax,exc=exc,filter=filter,fixpos=fixpos,zero=zero,finite=finite,expfilt=expfilt,outputpos=outputpos,$
snrval=snrval,snrmap=snrmap,noisemap=noisemap,imcol=imcol,silent=silent
;,left=left,right=right

;***subroutine version. returns values
;computes a SNR map of an image 
; -1. blindly
; -2. avoiding regions around a point source if "coord" switch is set
;note that not using 'alt' produces potentially skewed results with any non-zero value of exc

snrval=-99

if keyword_set(coord) then begin
xc=coord[0]
yc=coord[1]
endif

if ~keyword_set(fwhm) then fwhm=5.
if ~keyword_set(rmax) then rmax=100.
if ~keyword_set(exc) then exc=1.5

fwhm=double(fwhm)

im=image

filtfact=1
filtfact0=3
if keyword_set(expfilt) then filtfact=filtfact0
if keyword_set(filter) then im-=filter_image(im,median=5*fwhm*filtfact)

s=size(im)
dimx=s[1] & dimy=s[2]
xcen=dimx/2 & ycen=dimy/2

;if you want to turn on the finite sample correction from Mawet et al. (2014)
if keyword_set(finite) then begin

radpoint=findgen(dimx/2)+1
finitecorrgrid=fltarr(n_elements(radpoint))

for i=0L,n_elements(radpoint)-1 do begin
if 2*!pi*radpoint[i]/fwhm le 2 then begin
finitecorrgrid[i] = 999.
endif else begin

if keyword_set(coord) then begin
finitecorrgrid[i]=t_cvf(2.867d-7,2*!pi*radpoint[i]/fwhm-2)/5.*(sqrt(1+1/(2*!pi*radpoint[i]/fwhm-2)))
endif else begin
finitecorrgrid[i]=t_cvf(2.867d-7,2*!pi*radpoint[i]/fwhm-1)/5.*(sqrt(1+1/(2*!pi*radpoint[i]/fwhm-1)))
endelse
;print,'rads',radpoint[i]/fwhm,finitecorrgrid[i]

endelse
endfor

dist_circle,raddist,[dimx,dimy],[xcen,ycen]
finitecorr=interpol(finitecorrgrid,radpoint,raddist,/nan)
endif

;Summing within an Aperture
if keyword_set(silent) then begin
im1=sumaper_im(im,fwhm/2.,rmax,/nan,/silent)
endif else begin
im1=sumaper_im(im,fwhm/2.,rmax,/nan)
endelse
;;***no masking

;**map of the median-average deviation

if keyword_set(zero) then begin

;median summed pixel value
profrad_tc,im1,1,0,rmax,p1d=zsig,rayon=rsig2,p2d=immed
im1-=immed
profrad_tc,abs((im1)/0.6745),1,0,rmax,p1d=sig,rayon=rsig,p2d=imnoise

endif else begin

profrad_tc,abs((im1)/0.6745),1,0,rmax,p1d=sig,rayon=rsig,p2d=imnoise
endelse
;**

;SNR map
if keyword_set(finite) then imnoise*=finitecorr
snrmap=im1/imnoise
;if keyword_set(finite) then snrmap/=finitecorr

;writefits,'snrmap.fits',snrmap
writefits,'noise.fits',imnoise
writefits,'im1.fits',im1
;stop


if keyword_set(coord) then begin

;position of point source
if ~keyword_set(fixpos) then begin
gcntrd,im,xc,yc,xcf,ycf,fwhm,/silent
endif else begin
xcf=xc
ycf=yc
endelse

;***masking point source
dist_circle,radpt,[dimx,dimy],xcf,ycf,/double
good=where(radpt gt exc*fwhm/2.,complement=nbad)
im2=im & im2[nbad]=!values.f_nan

;Summing within an Aperture w/ masking
if keyword_set(silent) then begin
im2=sumaper_im(im2,fwhm/2.,rmax,/nan,/silent)
endif else begin
im2=sumaper_im(im2,fwhm/2.,rmax,/nan)
endelse

if keyword_set(zero) then begin

;median summed pixel value
profrad_tc,im2,1,0,rmax,p1d=zsig,rayon=rsig2,p2d=immed
im2-=immed
profrad_tc,abs((im2)/0.6745),1,0,rmax,p1d=sig,rayon=rsig,p2d=imnoise

endif else begin

profrad_tc,abs((im2)/0.6745),1,0,rmax,p1d=sig,rayon=rsig,p2d=imnoise
endelse

;SNR map
if keyword_set(finite) then imnoise*=finitecorr
snrmap=im1/imnoise

;writefits,'snrmap_obj.fits',snrmap
;writefits,'noise_objmasked.fits',imnoise

aper,im,xcf,ycf,flux,eflux,sky,skyerr,10,fwhm*0.5,/nan,/exact,setskyval=0,/flux,/silent
;print,'Estimated Position Is ',xcf,ycf
;print,'Integrated Signal Is ',flux
;print,'1-Sigma Noise Is ',imnoise[xcf,ycf]

if finite(flux) eq 0 then begin

endif else begin
snrval=flux/imnoise[xcf,ycf]
endelse


endif

if keyword_set(finite) then begin
;noisemap=imnoise*finitecorr
noisemap=imnoise
endif else begin
noisemap=imnoise
endelse
;snrmap=snrmap
imcol=im1

if keyword_set(coord) then begin
outputpos=[xcf,ycf]
endif

;end



end
