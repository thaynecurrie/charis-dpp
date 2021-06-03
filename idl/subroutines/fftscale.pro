function fftscale,imin,scxin,scyin,rapprec,mtfpixcor=mtfpixcor,$
                  silent=silent,scxout=scxout,scyout=scyout,dimout=dimout,$
                  parms=parms
;+
; NAME: FFTSCALE
;   Change l'echelle spatiale d'une image par FFT. Determination des
;   echelles des images et scale des images
;   Change the spatial scale of an image using FFTs. Comments
;   translated for a small fee.
;
; INPUTS:
; 	imin		an image
; 	scxin		scale in X
; 	scyin		scale in Y  **MUST EQUAL SCXIN RIGHT NOW**
; 	rapprec		Fractional precision, typically 1e-7
;
; KEYWORDS:
;
; OUTPUTS:
;
; NOTES:
;       From Christian:
;       Call:
;	
;	im2=fftscale(im1,scalex,scaley,precision)
;	
;	where scalex and scaley are the scale parameters in X and Y. You 
;	currently needs to use the same scale in both axis since I haven't 
;	modified the soft to do a none equal spatial scale. The precision 
;	keyword is the scale precision that you want (Delta scale/scale, 
;	typically ~10^-7). You have the option to modify the pixel MTF, 
;	but that is not really useful for simulated data. The soft does a 
;	2 iterations image & Fourier planes zero padding to find the best 
;	possible combination of zero padding to do the image spatial scale. 
;	The worse case scenario is to do a ~1x scale. What you want to do in 
;	these cases is to do a combination of 2 scales to avoid a ~1x scale 
;	(like if you want to do a 1.01x scale, you do a 1.1x scale followed by 
;	a (1.01/1.1)x scale).
;
; HISTORY:
; 	By Christian Marois
; 	2008-03-07  Documentation in English added by Marshall Perrin
;       07.30.2012 - Rewrite to accomodate non-even number of pixels - ds
;-

scx=double(scxin)
scy=double(scyin)

im=double(imin)
s=size(im)
dimx=s[1]
dimy=s[2]
meilprec=1E10
found=0

if scx ne scy then begin
   message,'Currently, fftscale only supports scx = scy.',/continue
   return,-1
endif

odd = dimx mod 2.
if odd then imcent = (dimx-1)/2. else imcent = (dimx/2.-0.5)

if keyword_set(parms) then begin
   if parms.scx ne scx or parms.scy ne scy or parms.dimx ne dimx or parms.dimy ne dimy then begin
      message,'Parameter structure dimensions or scales do not match inputs or requested scales. Returning.',/continue
      return,-1
   endif
   padinit = parms.padinit
   padinit2 = parms.padinit2
   padfft = parms.padfft
   padfft2 = parms.padfft2
   scf = parms.scf
   meilprec = parms.meilprec

endif else begin
   if scx ge 1.d then begin
      for i=0,1023,2 do begin
         padinit=dimx+i
         for j=0,scx*padinit-1,2 do begin
            padfft=padinit+j
            scx1=double(padfft)/double(padinit)
            for k=0,(scx-scx1)*dimx,2 do begin
               padinit2=dimx+k
               padfft2=round((scx/scx1)*padinit2)
               scx2=double(padfft2)/double(padinit2) 
               
               prec=abs(scx1*scx2/scx-1.d)
               scf=scx1*scx2
               
               if abs(prec) lt meilprec and scx ge 1.d and scx1 ge 1.d and scx2 ge 1.d then begin
                  padinitF=padinit
                  padfftF=padfft
                  padinit2F=padinit2
                  padfft2F=padfft2
                  if keyword_set(scxout) then scxout=scf
                  if keyword_set(scyout) then scyout=scf
                  meilprec=abs(prec)
                  scfF=scf
                  if abs(prec) lt rapprec then found=1
                  if abs(prec) lt rapprec then break
               endif
               
            endfor
            if abs(prec) lt rapprec then break
         endfor
         if abs(prec) lt rapprec then break
      endfor
   endif else begin
      for i=0,1023,2 do begin
         padinit=dimx+i
         for j=0,padinit-scx*padinit-1,2 do begin
            padfft=padinit-j
            scx1=double(padfft)/double(padinit)
            for k=0,1023,2 do begin
               padinit2=dimx+k
               padfft2=round((scx/scx1)*padinit2)
               scx2=double(padfft2)/double(padinit2)
               
               prec=abs(scx1*scx2/scx-1.d)
               scf=scx1*scx2
               
               if abs(prec) lt meilprec and scx lt 1.d and padfft2 lt 1023 then begin
                  padinitF=padinit
                  padfftF=padfft
                  padinit2F=padinit2
                  padfft2F=padfft2
                  scfF=scf
                  if keyword_set(scxout) then scxout=scf
                  if keyword_set(scyout) then scyout=scf
                  meilprec=abs(prec)
                  if abs(prec) lt rapprec then found=1
                  if abs(prec) lt rapprec then break
               endif
            endfor
            if abs(prec) lt rapprec then break
         endfor
         if abs(prec) lt rapprec then break
      endfor
   endelse

   if found eq 0 then begin
      if not keyword_set(silent) then print,'Pas de scale trouve... Conservation du meilleur scale'
   endif

   padinit=padinitF
   padinit2=padinit2F
   padfft=padfftF
   padfft2=padfft2F
   scf=scfF
   parms = create_struct('padinit',padinit,'padinit2',padinit2,'padfft',padfft,'padfft2',padfft2,'scf',scf,$
                         'dimx',dimx,'dimy',dimy,'scx',scx,'scy',scy,'meilprec',meilprec)
endelse

if not keyword_set(silent) then begin
   print,'Dim init = ',dimx,' Dim init + pad = ',padinit
   print,'Dim FFT init = ',padinit,' Dim FFT init + pad = ',padfft
   
   print,'Dim init 2 = ',dimx,' Dim init + pad = ',padinit2
   print,'Dim FFT init 2 = ',padinit2,' Dim FFT init + pad = ',padfft2
   
   print,'Scale demande = ',scx,' Scale obtenu = ',scf
   print,'Precision = ',meilprec
endif

;ITERATION 1
if padinit ne padfft then begin
   padim=dblarr(padinit,padinit)
   if padinit mod 2. then initcent = (padinit-1.)/2. else initcent = (padinit/2.) - 0.5
   
   padim[initcent-imcent:initcent+imcent,initcent-imcent:initcent+imcent] = im
   padim=smartshift(padim,-padinit/2.,-padinit/2.,/nofftw) 
   
   imfft=fft(padim,1,/double)
   imfftr=smartshift(real_part(imfft),padinit/2.,padinit/2.,/nofftw)
   imffti=smartshift(imaginary(imfft),padinit/2.,padinit/2.,/nofftw)
   
   padfftr=dblarr(padfft,padfft)
   padffti=dblarr(padfft,padfft)
   if padfft mod 2. then fftcent = (padfft-1.)/2. else fftcent = (padfft/2.) - 0.5

   if padfft ge padinit then begin
      padfftr[fftcent-initcent:fftcent+initcent,fftcent-initcent:fftcent+initcent] = imfftr
      padffti[fftcent-initcent:fftcent+initcent,fftcent-initcent:fftcent+initcent] = imffti
   endif else begin
      padfftr = imfftr[initcent-fftcent:initcent+fftcent,initcent-fftcent:initcent+fftcent]
      padffti = imffti[initcent-fftcent:initcent+fftcent,initcent-fftcent:initcent+fftcent]
   endelse
   
   padfftim=dcomplex(smartshift(padfftr,-padfft/2.,-padfft/2.,/nofftw),smartshift(padffti,-padfft/2.,-padfft/2.,/nofftw))
   
   psf=fft(padfftim,-1,/double)   
   rpsf=smartshift(real_part(psf),padfft/2.,padfft/2.,/nofftw)
   
   if padfft ge dimx then begin
      im = rpsf[fftcent-imcent:fftcent+imcent,fftcent-imcent:fftcent+imcent]
   endif else begin
      im = dblarr(dimx,dimy)
      im[imcent-fftcent:imcent+fftcent,imcent-fftcent:imcent+fftcent] = rpsf
   endelse
endif

;ITERATION 2
if padinit2 ne padfft2 then begin
   padim=dblarr(padinit2,padinit2)
   if padinit2 mod 2. then initcent2 = (padinit2-1.)/2. else initcent2 = (padinit2/2.) - 0.5

   padim[initcent2-imcent:initcent2+imcent,initcent2-imcent:initcent2+imcent] = im
   padim=smartshift(padim,-padinit2/2.,-padinit2/2.,/nofftw) 
   
   imfft=fft(padim,1,/double)
   imfftr=smartshift(real_part(imfft),padinit2/2.,padinit2/2.,/nofftw)
   imffti=smartshift(imaginary(imfft),padinit2/2.,padinit2/2.,/nofftw)
   
   padfftr=dblarr(padfft2,padfft2)
   padffti=dblarr(padfft2,padfft2)
   if padfft2 mod 2. then fftcent2 = (padfft2-1.)/2. else fftcent2 = (padfft2/2.) - 0.5

   if padfft2 ge padinit2 then begin
      padfftr[fftcent2-initcent2:fftcent2+initcent2,fftcent2-initcent2:fftcent2+initcent2] = imfftr
      padffti[fftcent2-initcent2:fftcent2+initcent2,fftcent2-initcent2:fftcent2+initcent2] = imffti
   endif else begin
      padfftr = imfftr[initcent2-fftcent2:initcent2+fftcent2,initcent2-fftcent2:initcent2+fftcent2]
      padffti = imffti[initcent2-fftcent2:initcent2+fftcent2,initcent2-fftcent2:initcent2+fftcent2]
   endelse
   
   if keyword_set(mtfpixcor) then begin
      if not keyword_set(silent) then print,"Correction pour la difference des mtfs des pixels..."
      scalefact=double(double(scx))
      pixmtf1=double(mtfpix(3.6,1.6,padfft2,0.018,xc=0.,yc=0.,scfact=scalefact))
      pixmtf2=double(mtfpix(3.6,1.6,padfft2,0.018,xc=0.,yc=0.))
      padfftim=dcomplex(shift(padfftr*pixmtf2/pixmtf1,-floor(padfft2/2),-floor(padfft2/2)),shift(padffti*pixmtf2/pixmtf1,-floor(padfft2/2),-floor(padfft2/2)))
   endif else padfftim=dcomplex(smartshift(padfftr,-padfft2/2.,-padfft2/2.,/nofftw),smartshift(padffti,-padfft2/2.,-padfft2/2.,/nofftw))
   
   psf=fft(padfftim,-1,/double) 
   rpsf=smartshift(real_part(psf),padfft2/2.,padfft2/2.,/nofftw)   
   fpsf=dblarr(dimx,dimy)
   
   if keyword_set(dimout) then dimx=dimout
   if padfft2 ge dimx then $
      fpsf = rpsf[fftcent2-imcent:fftcent2+imcent,fftcent2-imcent:fftcent2+imcent] $
   else $
      fpsf[imcent-fftcent2:imcent+fftcent2,imcent-fftcent2:imcent+fftcent2] = rpsf
endif else fpsf = im

return,fpsf

end





