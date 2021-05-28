pro charis_calc_ee,pick=pick,position=position,il=il,prad=prad,dpsf=dpsf,bckgd=bckgd,help=help

print,"charis_calc_ee,pick=pick,position=position,il=il,prad=prad,dpsf=dpsf,bckgd=bckgd,help=help"
print,""
print,"quick script that calculates encircled energy from empirical PSF derived from sat spots"
print,"***Keywords"
print,"*pick - pick the file, otherwise perform on psfcube_med.fits"
print,"*position - position of psf [default: center of cube]"
print,"*il - wavelength channel"
print,"*prad - subtract off radial profile? [OK if the PSF is off-axis]"
print,"*dpsf - defines range of E.E. calculation"
print,"*bckgd - subtract a pre-defined background value"
print,"*


if ~keyword_set(rbk) then rbk=[10,15]

if ~keyword_set(il) then il=10


;quick script that calculates encircled energy of a given image



if keyword_set(pick) then  begin

file_in=dialog_pickfile(Title="Select Data Cube for Calculation")

cube=readfits(file_in,/ext,/silent)

endif else begin

datadir='./psfmodel/'
file_in='psfcube_med.fits'

cube=readfits(datadir+file_in,/ext,/silent)


endelse

;now calculate the encircled energe per passband

sz=size(cube,/dim)

if n_elements(sz) eq 2 then nwvlh=1 else nwvlh = sz[2]

if ~keyword_set(dpsf) then dpsf=41

dpsf= dpsf < sz[1]

ee=fltarr(dpsf,nwvlh)
fluxtot=fltarr(nwvlh)

rap=findgen(dpsf)/2.

if ~keyword_set(position) then position=[sz[0]/2,sz[1]/2]
print,position


for i=0L,nwvlh-1 do begin
;fluxtot[i]=total(cube[*,*,i])
if keyword_set(bckgd) then cube[*,*,i]-=bckgd
fluxtot[i]= charis_myaper(cube[*,*,i],position[0],position[1],dpsf/2.,/nan)


 for j=0L,dpsf-1 do begin
   fluxr=charis_myaper(cube[*,*,i],position[0],position[1],rap[j],/nan)
   ee[j,i]=fluxr/fluxtot[i]

   if i eq il then print,dpsf/2,rap[j],fluxr,fluxtot[i]
 endfor
endfor

print,max(rap),dpsf/2
print,fluxtot
;now, to plot

dum=cube[*,*,il]
if keyword_set(prad) then begin
profrad_tc,dum,1,1,101,p2d=pr

dum-=pr
endif

dist_circle,g,sz[1],position[0],position[1]

good=where(g le dpsf/4)
print,max(dum[good])
z=where(dum eq max(dum[good]))
profradxy,dum,position[0],position[1],1,dpsf/2,p1d=pr,rayon=rayon
pr/=max(dum[good])

  plot,rap,ee[*,il],xrange=[0,max(rap)+1]
  oplot,rayon,pr,linestyle=1

print,pr
print,max(dum[good])


end
