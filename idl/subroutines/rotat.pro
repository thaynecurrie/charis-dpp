function rotat,im,deg,x0,y0,hdr=hdr,missing=missing
;deg is clockwise angle of rotation
;x0,y0: pivot point, default is dimx/2,dimy/2

sz=size(im) & dimx=sz[1] & dimy=sz[2]
if n_params() lt 4 then y0=dimy/2
if n_params() lt 3 then x0=dimx/2
;sxaddpar,hdr,'crpix1',x0
;sxaddpar,hdr,'crpix2',y0

;return,rot(im,deg,1.0,x0,y0,missing=missing,cubic=-0.5,/pivot)

theta=deg*!dpi/180.
ct=cos(theta)
st=sin(theta)

;les lignes suivantes font la rotation sans interpoler la ou il y
;a des NAN, beaucoup plus rapide
imt=make_array(size=sz,value=!values.f_nan)
imtnoint=imt
igood=where(finite(im) eq 1,cgood)
if cgood eq 0 then return,imt

;coordonnees initiales des bons pixels
xi=igood mod dimx
yi=igood/dimx

;vecteur des coordonnees finales correspondant aux bons pixels
xf=round((x0-ct*x0-st*y0)+ct*xi+st*yi)
yf=round((y0+st*x0-ct*y0)-st*xi+ct*yi)

;comme les coordonnees initiales tournees peuvent ne pas couvrir
;entierement la grille de pixels, on inclut dans les coordonnees
;finales les 9 pixels les plus proches de chaque rotation de (xf,yf)

mask=bytarr(dimx,dimy)
mask[xf,yf]=1
mask=dilate(mask,bytarr(3,3)+1)
indf=where(mask eq 1)
xf=indf mod dimx
yf=indf/dimx

;coordonnees initiales correspondant a (xf,yf)
xi=(x0-ct*x0+st*y0)+ct*xf-st*yf
yi=(y0-st*x0-ct*y0)+st*xf+ct*yf

imt[indf]=interpolate(im,xi,yi,cubic=-0.5,missing=missing)
imtnoint[indf]=interpolate(im,xi,yi,cubic=0,missing=missing)
wnan=where(~finite(imt), nanct)
if nanct gt 0 then imt[wnan] = imtnoint[wnan]

;meme chose que rot
;theta=deg*!dpi/180.
;ct=cos(theta)
;st=sin(theta)
;p=[[x0-ct*x0+st*y0,-st],[ct,0.]]
;q=[[y0-st*x0-ct*y0,ct],[st,0.]]
;imt=poly_2d(im,p,q,cubic=-0.5)

;met astrometrie a jour, si elle est definie
if keyword_set(hdr) then begin
    extast, hdr, astr
    crpix=astr.crpix
    cd=astr.cd

    rot_mat=[ [ ct, st], [-st, ct] ] 

    ;nouvelles valeurs
    crpix=transpose(rot_mat)#(crpix-1-[x0,y0])+1+[x0,y0]
    cd=cd#rot_mat
    astr.crpix=crpix
    astr.cd=cd
    ;met dans header
    putast,hdr,astr
endif

return,imt
end
