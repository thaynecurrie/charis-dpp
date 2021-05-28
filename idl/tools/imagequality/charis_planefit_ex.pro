pro charis_planefit_ex

;an example of plane-fitting data

image=dialog_pickfile(Title="Select the image")
;image=readfits(image)
image=readfits(image,/ext)

sourcepos=[100,62]
rad=4

sourceind=get_xycind(201,201,sourcepos[0],sourcepos[1],2*rad)
bckind0=get_xycind(201,201,sourcepos[0],sourcepos[1],rad*4)
sourcedist=sqrt((sourcepos[0]-201/2)^2.+(sourcepos[1]-201/2)^2.)
dist_circle,g,201
annulus=where(g gt sourcedist-rad*3 and g lt sourcedist+rad*3)
bckind0=intersect(bckind0,annulus)
;bckind=setdifference(sourceind,bckind0)
;bckind=bckind0

;bckind=intersect(annulus,sourceind,/xor)
bckind=intersect(bckind0,sourceind,/xor)

background=bckind
;[where(finite(image) eq 1)]

weights=0
coeff=planefit(background mod 201, background / 201,image[background],weights,yfit)

xinds=sourceind mod 201
yinds=sourceind / 201

src_bkg_plane=coeff[0]+coeff[1]*xinds+coeff[2]*yinds

xwing=background mod 201
ywing=background / 201
wings=coeff[0]+coeff[1]*xwing+coeff[2]*ywing

for i=0L,n_elements(xwing)-1 do begin


print,'i ',i,' coords ',xwing[i],ywing[i]

endfor
blah=fltarr(201,201)

blah[sourceind]=src_bkg_plane
writefits,'blah.fits',blah

blah2=fltarr(201,201)
blah2[background]=wings

writefits,'blah2.fits',blah2

writefits,'total.fits',image-(blah+blah2)



end
