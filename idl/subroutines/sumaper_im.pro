function sumaper_im,im,rsum,rmax,silent=silent,nan=nan,reducdir=reducdir

;Retourne une image dans laquelle chaque pixel consiste en la somme de
;im sur une ouverture centree sur ce pixel.
;
;Parametres d'entree:
;-------------------
;im: image de la psf a analyser
;rsum: rayon de l'ouverture a prendre pour chaque pixel
;rmax: rayon maximum a considerer, par defaut=100
;
;=
;
;Returns an image in which each pixel consists of the sum of "im" over
;an aperture centered on that pixel. {Basically, it's convolution.}
;
;Input parameters:
;-----------------
;im: PSF image to be analyzed
;rsum: radius of the sampling aperture
;rmax: maximum radius to consider, by default = 100

im2=im*0.

s=size(im) & dimx=s[1] & dimy=s[2]
xc=dimx/2 & yc=dimy/2

;if (n_params() lt 3) then rmax=sqrt((dimx/2.)^2+(dimy/2.)^2)

;indices des pixels a analyser
;=
;indices of the pixels to be analyzed
; MWM: Edited code to produce get_cind.fits for all of our data

if ((n_params() lt 3) and keyword_set(reducdir)) then begin
   print, 'Beware: Using '+reducdir+'get_cind.fits in sumaper_im.pro.  MWM 04/29/09.'
   ind=readfits(reducdir+'get_cind.fits')
endif else begin
   ind=get_cind(dimx,dimy,rmax)
endelse

   ; if ~keyword_set(rmax) then rmax=715	;; RKmod
   ; ind=get_cind(dimx,dimy,rmax)		;; RKmod

npts=n_elements(ind)

;masque de l'ouverture et indices des pixels du masque (pour pixel central)
;=
;mask of the aperture, and indices of the mask pixels (for the central pixel)
ic=charis_myaper(im,xc,yc,rsum,mask,imask)
;ic=myaper(im,xc,yc,rsum,mask,imask)
;if (not keyword_set(silent)) then print,'Calcul de l''image:   0.0%',format='(a,$)'
if (not keyword_set(silent)) then print,'Image calculation:    0.0%',format='(a,$)'
for n=long(0),npts-1 do begin
    x=ind[n] mod dimx & y=ind[n]/dimx
    ;indices des pixels du masque (centre sur ce pixel)
    ;=
    ;indices of the mask pixels (centered on this pixel)
    i=imask+(x-xc)+(y-yc)*dimx
    ;intensite totale dans l'ouverture
    ;=
    ;total intensity in the aperture
    im2[x,y]=total(mask*im[i],nan=nan)

    if (not keyword_set(silent)) then $
      print,replicate(string(8b),6),(n+1)/float(npts)*100.,'%',format='(6a,f5.1,a,$)'
endfor ;n
;if (not keyword_set(silent)) then print
;stop
return,im2
end
