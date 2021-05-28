pro snratio, image,coord=coord,fwhm=fwhm,expfilt=expfilt,rmax=rmax,exc=exc,ext=ext,filter=filter,fixpos=fixpos,zero=zero,finite=finite,$
snrval=snrval,left=left,right=right
;computes a SNR map of an image 
; -1. blindly
; -2. avoiding regions around a point source if "coord" switch is set
;note that not using 'alt' produces potentially skewed results with any non-zero value of exc



if keyword_set(coord) then begin
xc=coord[0]
yc=coord[1]
endif

if ~keyword_set(fwhm) then fwhm=5.
if ~keyword_set(rmax) then rmax=100.
if ~keyword_set(exc) then exc=1.5

fwhm=double(fwhm)

if (~keyword_set(ext)) then begin
im=readfits(image,h1)
endif else begin
im=readfits(image,h1,/exten)
endelse

if keyword_set(left) then begin
generate_grids,xg,yg,(size(im,/dim))[0],/whole,/freqshift
xg-=min(xg)
yg-=min(yg)
bad=where(xg gt (size(im,/dim))[0]/2)
im[bad]=!values.f_nan

endif
if keyword_set(right) then begin
generate_grids,xg,yg,(size(im,/dim))[0],/whole,/freqshift
xg-=min(xg)
yg-=min(yg)
bad=where(xg lt (size(im,/dim))[0]/2)
im[bad]=!values.f_nan

endif

expfilt=1
expfilt0=3
if keyword_set(expfilt) then expfilt=expfilt0
if keyword_set(filter) then im-=filter_image(im,median=5*expfilt*fwhm)

s=size(im)
dimx=s[1] & dimy=s[2]
xcen=dimx/2 & ycen=dimy/2

;if you want to turn on the finite sample correction from Mawet et al. (2014)
if keyword_set(finite) then begin
radpoint=findgen(dimx/2)+1
finitecorrgrid=fltarr(n_elements(radpoint))

for i=0L,n_elements(radpoint)-1 do begin
if 2*!pi*radpoint[i]/fwhm le 2 then begin
finitecorrgrid[i] = 99.

endif else begin

;***fix the below (4/27/2019)
if ~(keyword_set(left) or keyword_set(right)) then begin
if keyword_set(coord) then begin
finitecorrgrid[i]=t_cvf(2.867d-7,2*!pi*radpoint[i]/fwhm-2)/5.*(sqrt(1+1/(2*!pi*radpoint[i]/fwhm-2)))
endif else begin
finitecorrgrid[i]=t_cvf(2.867d-7,2*!pi*radpoint[i]/fwhm-1)/5.*(sqrt(1+1/(2*!pi*radpoint[i]/fwhm-1)))
endelse

endif else begin

if keyword_set(coord) then begin
finitecorrgrid[i]=t_cvf(2.867d-7,!pi*radpoint[i]/fwhm-2)/5.*(sqrt(1+1/(!pi*radpoint[i]/fwhm-2)))
endif else begin
finitecorrgrid[i]=t_cvf(2.867d-7,!pi*radpoint[i]/fwhm-1)/5.*(sqrt(1+1/(!pi*radpoint[i]/fwhm-1)))
endelse

endelse

endelse
endfor

dist_circle,raddist,[dimx,dimy],[xcen,ycen]
finitecorr=interpol(finitecorrgrid,radpoint,raddist,/nan)
writefits,'finitecorr.fits',finitecorr
endif

;Summing within an Aperture
im1=sumaper_im(im,fwhm/2.,rmax,/nan)

if (keyword_set(left) or keyword_set(right)) then $
im1[bad]=!values.f_nan
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

if keyword_set(finite) then begin
unres=where(raddist/fwhm le 0.5)
snrmap[unres]=0
imnoise[unres]=0
endif

writefits,'snrmap.fits',snrmap
writefits,'noise.fits',imnoise
writefits,'im1.fits',im1

if keyword_set(coord) then begin

;position of point source
if ~keyword_set(fixpos) then begin
gcntrd,im,xc,yc,xcf,ycf,fwhm
endif else begin
xcf=xc
ycf=yc
endelse

;***masking point source
dist_circle,radpt,[dimx,dimy],xcf,ycf,/double
good=where(radpt gt exc*fwhm/2.,complement=nbad)
im2=im & im2[nbad]=!values.f_nan

;Summing within an Aperture w/ masking
im2=sumaper_im(im2,fwhm/2.,rmax,/nan)

if keyword_set(zero) then begin

;median summed pixel value
profrad_tc,im2,1,0,rmax,p1d=zsig,rayon=rsig2,p2d=immed
im2-=immed
profrad_tc,abs((im2)/0.6745),1,0,rmax,p1d=sig,rayon=rsig,p2d=imnoise

endif else begin

profrad_tc,abs((im2)/0.6745),1,0,rmax,p1d=sig,rayon=rsig,p2d=imnoise
endelse

if keyword_set(finitecorr) then imnoise*=finitecorr

;SNR map
snrmap=im1/imnoise

writefits,'snrmap_obj.fits',snrmap
writefits,'noise_objmasked.fits',imnoise

aper,im,xcf,ycf,flux,eflux,sky,skyerr,10,fwhm*0.5,/nan,/exact,setskyval=0,/flux,/silent
print,'Estimated Position Is ',xcf,ycf
print,'Integrated Signal Is ',flux
print,'1-Sigma Noise Is ',imnoise[xcf,ycf]

if finite(flux) eq 0 then begin
if ~keyword_set(finite) then begin
print,'approx SNR ratio is ',im1[xcf,ycf]/imnoise[xcf,ycf]
endif else begin
print,'approx SNR ratio is ',im1[xcf,ycf]/imnoise[xcf,ycf]
endelse
endif else begin
if ~keyword_set(finite) then begin
print,'SNR ratio is ',flux/imnoise[xcf,ycf]
snrval=flux/imnoise[xcf,ycf]
endif else begin
print,'SNR ratio is ',flux/imnoise[xcf,ycf]
snrval=flux/imnoise[xcf,ycf]
endelse
endelse

endif

end



end
