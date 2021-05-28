function charis_getrelshift,im1,im2,roi,sxg,syg,subpix=subpix,res=res

;im1=reference image
;im2=image to shift 

;shift initial
if n_params() lt 4 then sx0=0. else sx0=sxg
if n_params() lt 5 then sy0=0. else sy0=syg
if ~keyword_set(subpix) then sx0=round(sxg)
if ~keyword_set(subpix) then sy0=round(syg)
s=size(im1) & d0=s[1]
xc=d0/2 & yc=d0/2

;conserve seulement les pixels inclus dans roi pour im1
;retains only the pixels used in ROI of im1 
im1t=im1[roi]
im1t-=mean(im1t,/nan)
im1t/=sqrt(total(im1t^2,/nan))

xroi=roi mod d0 & yroi=roi/d0
if ~keyword_set(subpix) then begin
    ;premiere iteration, nombre entier de pixels
    ntry=21
    cor=fltarr(ntry,ntry)
    sx=sx0+(findgen(ntry)-ntry/2)
    sy=sy0+(findgen(ntry)-ntry/2)
    for i=0,ntry-1 do begin
        for j=0,ntry-1 do begin
            im2t=interpolate(im2,xroi-sx[i],yroi-sy[j],cubic=-0.5,missing=0.)
            im2t-=mean(im2t,/nan)
            im2t/=sqrt(total(im2t^2,/nan))
            cor[i,j]=total(im1t*im2t,/nan)
            if keyword_set(res) then cor[i,j]=-total(abs(im1t-im2t),/nan)
        endfor
    endfor
    m=max(cor,ind) & i=ind mod ntry & j=ind/ntry
    sx0=sx[i] & sy0=sy[j]
endif

;iterate until 10^k-pixel
;try all shifts on a square grid
ntry=21
for k=1,3 do begin
    sx=sx0+(findgen(ntry)-10.)*10.^(-k)
    sy=sy0+(findgen(ntry)-10.)*10.^(-k)
    cor=fltarr(ntry,ntry)
    
    for i=0,ntry-1 do begin
        for j=0,ntry-1 do begin
            im2t=interpolate(im2,xroi-sx[i],yroi-sy[j],cubic=-0.5,missing=0.)
            im2t-=mean(im2t,/nan)
            im2t/=sqrt(total(im2t^2,/nan))
            cor[i,j]=total(im1t*im2t,/nan)
            if keyword_set(res) then cor[i,j]=-total(abs(im1t-im2t),/nan)
            ;if(abs(sy[j] + 1.1) lt 0.1)then cor[i,j] =!values.f_nan
        endfor
    endfor
    g=mpfit2dpeak(cor)
    ;m=max(cor,ind) 
    m=max(g,ind)
    i=ind mod ntry & j=ind/ntry
    sx0=sx[i] & sy0=sy[j]
endfor ;k

;return,[sx0,sy0]
return,[xc-sx0,yc-sy0]
end
