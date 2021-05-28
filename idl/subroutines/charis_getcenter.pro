function charis_getcenter,imt,xcg,ycg,roi,smallsteps=smallsteps

;not verified, probably not robust.  Use only if you have to.

;corlim=5
;corlim2=1.25
corlim=3
corlim2=1.

im=imt
imr=rotat(im,180)
;imr=reverse(im,2)
;writefits,'im.fits',imt
;writefits,'imr.fits',imr
;writefits,'diff.fits',imt-imr

s=size(im) & d0=s[1]
xc=d0/2 & yc=d0/2
xroi=roi mod d0
yroi=roi/d0

;print,xroi
;print,yroi

;stop

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
;print,sx0,sy0,xcg,ycg
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
        im1=interpolate(im,xroi-sx[i],yroi-sy[j],cubic=-0.5,missing=0.)
        imr1=interpolate(imr,xroi+sx[i],yroi+sy[j],cubic=-0.5,missing=0.)
        im1-=mean(im1,/nan)
        im1/=sqrt(total(im1^2.,/nan)) 
        imr1-=mean(imr1,/nan)
        imr1/=sqrt(total(imr1^2.,/nan))

        ;im1=im
        ;reverse(im1)=(image shiftee centrosymmetrique)[roi]
        ;writefits,'im1.fits',im1

        ;and this is equally true that the region of interest is center-symmetric
        ; and that the region of interest indices are in ascending order
        ;ceci est vrai tant que la roi est centrosymmetrique
        ;et que les indices de roi sont en ordre croissant
        ;help,im1,imr1
        ;cor[i,j]=correlate(im1,imr1)
        ;cor[i,j]=correl_images(im1,imr1)
        ;cor[i,j]=total(im1*imr1,/nan)
        ;cor[i,j]=total(im1*reverse(reverse(im1)),/nan)
        cor[i,j]=-total(abs(im1-imr1),/nan)
        ;cor[i,j]=-total(abs(im1-rotat(im1,180)),/nan)
        ;cor[i,j]=robust_corr(im1,reverse(reverse(im1,1),2))
        ;print,'corr ',cor[i,j],'corrr',sx[i],sy[j],sx0,sy0
        ;stop
;        print,mean(im1),mean(imr1),total(im1^2.),total(imr1^2.),sx[i],sy[j],cor[i,j]
    endfor
endfor
m=max(cor,ind) & i=ind mod ntry & j=ind/ntry
;print,xc-sx[i],yc-sy[j],sx[i],sy[j],m
;;plot,cor
;;stop
sx0=sx[i] & sy0=sy[j]
sind=sort(-cor)
;print,cor[i,j]/robust_sigma(cor)
;print,cor[sind[0]]/cor[sind[1]]
;stop
goto,skipoverme
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
skipoverme:
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
        imr1=interpolate(imr,xroi+sx[i],yroi+sy[j],cubic=-0.5,missing=0.)
        im1-=mean(im1,/nan)
        im1/=sqrt(total(im1^2.,/nan)) 
        imr1-=mean(imr1,/nan)
        imr1/=sqrt(total(imr1^2.,/nan))
        ;reverse(im1)=(image shiftee centrosymmetrique)[roi]
        ;cor[i,j]=total(im1*imr1,/nan)
        ;cor[i,j]=correlate(im1,imr1)
        ;cor[i,j]=total(im1*reverse(im1),/nan)
        cor[i,j]=-total(abs(im1-imr1),/nan)
        ;cor[i,j]=total(im1*imr1,/nan)
        ;cor[i,j]=total(im1*reverse(reverse(im1)),/nan)
        ;print,'corr ',cor[i,j],'corrr',sx[i],sy[j],sx0,sy0
        ;print,mean(im1),mean(imr1),total(im1^2.),total(imr1^2.),sx[i],sy[j],cor[i,j]
    endfor
    g=mpfit2dpeak(cor)
    m=max(cor,ind) & i=ind mod ntry & j=ind/ntry
    sx0=sx[i] & sy0=sy[j]
    ;plot,cor
    
endfor ;k
;stop
;stop  ;imm=translate(im,sx0,sy0)
      ;blinkdl,imm*mask,centrosymm(imm*mask),0.5,box=512
;print,xc,yc,sx0,sy0
;stop
return,[xc-sx0,yc-sy0]
end
