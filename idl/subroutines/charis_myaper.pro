function charis_myaper,im,xc,yc,r,mask,imask,nan=nan

;***** version of myaper.pro from the GPI Data Reduction Pipeline/David Lafreniere's pipeline specific for CHARIS (to avoid collisions with uses elsewhere)
;+
; NAME:
;      myaper
; PURPOSE:
;       Do simple aperture photometry, simple alternative to aper 
;
; EXPLANATION:
;       Using a mask 
;
; Calling SEQUENCE:
;      
;
; INPUT/OUTPUT:
;
;im: image, 2-d array
;xc, yc: center of aperture
;r: radius of aperture 
;mask: mask corresponding to aperture, centered on round(xc,yc)
;imask: indices of mask pixels, expressed wrt im
; OPTIONAL OUTPUT:
;       
;
; EXAMPLE:
;
;
; DEPENDENCIES:
; pixwt.pro
;
; NOTES: 
;      
;             
; REVISION HISTORY
;       Written  before 2008, Mathilde beaulieu/Jean-Francois Lavigne/David Lafreniere/Jerome Maire. 
;-



if (r lt 0.) then return,-1
if (r eq 0.) then return,im[round(xc),round(yc)]
s=size(im) & dimx=s[1] & dimy=s[2]

;half-dimension of subsection
re=ceil(r)
ixmin=round(xc)-re & ixmax=round(xc)+re
iymin=round(yc)-re & iymax=round(yc)+re
ixmin=ixmin>0      & iymin=iymin>0
ixmax=ixmax<(dimx-1) & iymax=iymax<(dimy-1)
nx=ixmax-ixmin+1     & ny=iymax-iymin+1

;x and y arrays for the image subsection
xm=findgen(nx)#replicate(1.,ny)+ixmin
ym=replicate(1.,nx)#findgen(1.,ny)+iymin
;pixel distances
distsq=(xm-xc)^2+(ym-yc)^2
;pixel index
imask=ym*dimx+xm

;constructing mask
mask=fltarr(nx,ny)
;indices of pixels within r+padding
ind=where(distsq le (r+0.70711)^2,nin) 
if (nin gt 0) then mask[ind]=1.
;fracitonal pixels
iedge=where(distsq le (r+0.70711)^2 and distsq ge ((r-0.70711)>0)^2,nedge)
for n=0,nedge-1 do mask[iedge[n]]=pixwt(xc,yc,r,xm[iedge[n]],ym[iedge[n]])
mask=mask>0.<1.
;--construit mask

return,total(im[imask]*mask,nan=nan)
end
