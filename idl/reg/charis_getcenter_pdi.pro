function charis_getcenter_pdi,imt,xcg,ycg,roi,smallsteps=smallsteps

;very simple cross-correlation fitting for centroid position of unsatured CHARIS PSF 
;corlim=5
;corlim2=1.25
corlim=3
corlim2=1.

im=imt
writefits,'im.fits',im
s=size(im) & d0=s[1]
xc=d0/2 & yc=d0/2
generate_grids,xz,yz,201

xroi=roi mod d0
yroi=roi/d0

;**diagnostique
mask=fltarr(d0,d0)
mask[roi]=1.
;**
img=imt
img[roi]=!values.f_nan
writefits,'mask.fits',img
;stop

;"guess" de shift au centre
;guess the shift of the center
sx0=round(xc-xcg) & sy0=round(yc-ycg)

;print,sx0,sy0
;stop
;premiere iteration, precision de 2 pixels
ntry=21
step=2.
;bandaid
if keyword_set(smallsteps) then step=0.5
cor=fltarr(ntry,ntry)
sx=sx0+(findgen(ntry)-ntry/2)*step
sy=sy0+(findgen(ntry)-ntry/2)*step
for i=0,ntry-1 do begin
    for j=0,ntry-1 do begin
        ;im1=interpolate(im,xroi-sx[i],yroi-sy[j],cubic=-0.5,missing=!values.f_nan)
        ;im1=interpolate(im,xroi-sx[i],yroi-sy[j],cubic=-0.5,missing=0.)
        ;imz=interpolate(im,xz-0,yz-0,cubic=-0.5,missing=!values.f_nan)
        im1=interpolate(im,xz-sx[i],yz-sy[j],cubic=-0.5,missing=0.)
        im1[where(im1 lt 0)]=!values.f_nan
        ;good=where(im1 gt 0)
        ;writefits,'imz.fits',im1
        ;stop
        ;reverse(im1)=(image shiftee centrosymmetrique)[roi]

        ;and this is equally true that the region of interest is center-symmetric
        ; and that the region of interest indices are in ascending order
        ;ceci est vrai tant que la roi est centrosymmetrique
        ;et que les indices de roi sont en ordre croissant

        cor[i,j]=total(im1[roi]*reverse(im1[roi]),/nan)
        ;cor[i,j]=robust_corr(im1,reverse(im1))
        ;print,cor[i,j],sx[i],sy[j]
        if ( long(sx[i]) eq 23 and long(sy[i]) eq 24) then writefits,'imz.fits',imz
        ;stop
    endfor
endfor

m=max(cor,ind) & i=ind mod ntry & j=ind/ntry
sx0=sx[i] & sy0=sy[j]
sind=sort(-cor)
if cor[i,j]/robust_sigma(cor) lt corlim $
  or cor[sind[0]]/cor[sind[1]] lt corlim2 then begin
    ;intersection of maximum correlation of the two spiders
    ;print,'hi',cor[i,j]/robust_sigma(cor)
    ;print,'hihi',cor[sind[0]]/cor[sind[1]]
    inter=inter_spider(cor-median(cor,corlim),spang)
    ;retains the shift corresponding to this intersection
    sx0=interpolate(sx,inter[0],cubic=-0.5) & sx0=round(sx0*100.)/100.
    sy0=interpolate(sy,inter[1],cubic=-0.5) & sy0=round(sy0*100.)/100.
endif


;interate, precision up to 0.01
;try all shifts on a square grid

ntry=11
step=[0.5,0.2,0.1,0.05,0.01]
for k=0,n_elements(step)-1 do begin
    sx=sx0+(findgen(ntry)-ntry/2)*step[k]
    sy=sy0+(findgen(ntry)-ntry/2)*step[k]
    cor=fltarr(ntry,ntry)
    
    for i=0,ntry-1 do for j=0,ntry-1 do begin
        ;im1=interpolate(im,xroi-sx[i],yroi-sy[j],cubic=-0.5,missing=0.)
        im1=interpolate(im,xz-sx[i],yz-sy[j],cubic=-0.5,missing=0.)
        ;reverse(im1)=(image shiftee centrosymmetrique)[roi]
        ;cor[i,j]=total(im1*reverse(im1),/nan)
        cor[i,j]=total(im1[roi]*reverse(im1[roi]),/nan)
        ;if (k eq 1) then print,cor[i,j],sx[i],sy[j]
    endfor

    m=max(cor,ind) & i=ind mod ntry & j=ind/ntry
    sx0=sx[i] & sy0=sy[j]
endfor ;k

;stop  ;imm=translate(im,sx0,sy0)
      ;blinkdl,imm*mask,centrosymm(imm*mask),0.5,box=512
;print,xc,yc,sx0,sy0
return,[xc-sx0,yc-sy0]
end
