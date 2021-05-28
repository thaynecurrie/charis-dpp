pro charis_northup,image,header0,header1,angoffset=angoffset,outputrot=outputrot
;takes a CHARIS cube and rotates it north up.

;if ~keyword_set(angoffset) then angoffset=0
if ~keyword_set(angoffset) then angoffset=charis_get_constant(name='angoffset')

extast,header1,astr,res
getrot,header1,npa,cdelt,/silent



sz=size(image)
centery=sz[2]/2
centerx=sz[1]/2

nslice=1
if sz[0] gt 2 then $
nslice=sz[3]

d_PAR_ANG=-npa

d_PAR_ANG+=angoffset

outputrot=d_PAR_ANG

hrot2, image[*,*,0], header1, nearest , header1, -d_PAR_ANG, $
         centerx, centery, 2,  interp=0,  missing=!values.f_nan
sxaddpar, header1, 'ROTANG', -d_PAR_ANG, 'The rotation from Rotate North Up primitive'

;loop over cube
;print,sz

;print,npa

for i=0L,nslice-1 do begin
;for i=0L,sz[3]-1 do begin
 interpolated = rotat(image[*,*,i],-d_PAR_ANG,centerx,centery,missing=!values.f_nan)
 ;interpolated = rot(reform(image[*,*,i],sz[1],sz[2]), -d_PAR_ANG,1., cubic=-0.5, centerx,centery, missing=!values.f_nan)
 ;interpolated = rot(reform(image[*,*,i],sz[1],sz[2]), -d_PAR_ANG,1., cubic=-0.5, /pivot, missing=!values.f_nan)
 ;writefits,'interp.fits',interpolated
 ;nearest = rot(reform(image[*,*,i],sz[1],sz[2]), -d_PAR_ANG, 1., interp=0, /pivot, missing=!values.f_nan)
 ;nearest = rot(reform(image[*,*,i],sz[1],sz[2]), -d_PAR_ANG, 1., interp=0, centerx,centery, missing=!values.f_nan)
 ;nearest = rot(reform(image[*,*,i],sz[1],sz[2]), -d_PAR_ANG, 1., interp=0, /pivot, missing=!values.f_nan)
;wnan = where(~finite(interpolated), nanct)
;if nanct gt 0 then interpolated[wnan] = nearest[wnan]
image[*,*,i] = interpolated
endfor

;image0=rotat(image[*,*,0],-100.d0,hdr=header1)
image0=rotat(image[*,*,0],angoffset,hdr=header1)
;writefits,'image.fits',image
end
