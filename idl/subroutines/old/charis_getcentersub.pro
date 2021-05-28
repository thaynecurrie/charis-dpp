function charis_getcentersub,imt,xcg,ycg,roi,smallsteps=smallsteps

;corlim=5
;corlim2=1.25
corlim=3
corlim2=1.

im=imt
imr=rotat(im,180)
s=size(im) & d0=s[1]
xc=d0/2 & yc=d0/2

;region of interest
xroi=roi mod d0
yroi=roi/d0

;print,xroi
;print,yroi

;**diagnostique
mask=fltarr(d0,d0)
mask[roi]=1.
;**
img=imt
img[roi]=!values.f_nan
writefits,'mask.fits',img

;"guess" de shift au centre
;guess the shift of the center
sx0=round(xc-xcg) & sy0=round(yc-ycg)
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
        im1=interpolate(im,xroi-sx[i],yroi-sy[j],cubic=-0.5,missing=0.)
        imr1=interpolate(imr,xroi-sx[i],yroi-sy[j],cubic=-0.5,missing=0.)
        ;im1=im
        ;reverse(im1)=(image shiftee centrosymmetrique)[roi]
        ;help,im1
        ;help,im1
        ;writefits,'im1.fits',im1
        ;stop
        ;and this is equally true that the region of interest is center-symmetric
        ; and that the region of interest indices are in ascending order
        ;ceci est vrai tant que la roi est centrosymmetrique
        ;et que les indices de roi sont en ordre croissant
        
        cor[i,j]=-total(im1*imr1,/nan)
        ;cor[i,j]=total(im1*reverse(im1),/nan)
        ;cor[i,j]=robust_corr(im1,reverse(im1))
        ;print,cor[i,j],'corrr',total(im1*im1,/nan)
        ;stop
        cor[i,j]=-total(abs(im1-imr1),/nan)
        ;cor[i,j]=sqrt(1./total(abs(im1-imr1),/nan)^2.)
    endfor
endfor

m=max(cor,ind) & i=ind mod ntry & j=ind/ntry
sx0=sx[i] & sy0=sy[j]
sind=sort(-cor)
print,sx0,sy0
stop
if cor[i,j]/robust_sigma(cor) lt corlim $
  or cor[sind[0]]/cor[sind[1]] lt corlim2 then begin
    ;intersection of maximum correlation of the two spiders
    print,'hi',cor[i,j]/robust_sigma(cor)
    print,'hihi',cor[sind[0]]/cor[sind[1]]
    inter=inter_spider(cor-median(cor,corlim),spang)
    ;retains the shift corresponding to this intersection
    sx0=interpolate(sx,inter[0],cubic=-0.5) & sx0=round(sx0*100.)/100.
    sy0=interpolate(sy,inter[1],cubic=-0.5) & sy0=round(sy0*100.)/100.
endif

;itere jusqu'a precision de 0.01
;essaie tous les shifts sur une grille carree

;interate, precision up to 0.01
;try all shifts on a square grid

ntry=11
step=[0.5,0.2,0.1,0.05,0.01]
for k=0,n_elements(step)-1 do begin
    sx=sx0+(findgen(ntry)-ntry/2)*step[k]
    sy=sy0+(findgen(ntry)-ntry/2)*step[k]
    cor=fltarr(ntry,ntry)
    
    for i=0,ntry-1 do for j=0,ntry-1 do begin
        im1=interpolate(im,xroi-sx[i],yroi-sy[j],cubic=-0.5,missing=0.)
        imr1=interpolate(imr,xroi-sx[i],yroi-sy[j],cubic=-0.5,missing=0.)
        ;reverse(im1)=(image shiftee centrosymmetrique)[roi]
        ;cor[i,j]=total(im1*imr1,/nan)
        ;cor[i,j]=total(im1*reverse(im1),/nan)
        cor[i,j]=-total(abs(im1-imr1),/nan)
        ;cor[i,j]=sqrt(1./total(abs(im1-imr1),/nan)^2.)
    endfor

    m=max(cor,ind) & i=ind mod ntry & j=ind/ntry
    sx0=sx[i] & sy0=sy[j]
endfor ;k

;stop  ;imm=translate(im,sx0,sy0)
      ;blinkdl,imm*mask,centrosymm(imm*mask),0.5,box=512
;print,xc,yc,sx0,sy0
return,[xc-sx0,yc-sy0]
end
