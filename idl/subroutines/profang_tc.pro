pro profang_tc,imt,drang,r1,r2,med4q=med4q,p2d=p2d,p1d=p1d,rayon=rayon,max=max

;Calculates a median angular profile.
; By default, calculates a angular profile 1d consisting of the median
; values of all pixels within an angular range set by res 
;certain distance from
;Then uses a cubic spline to extrapolate the profile of
;, length of each pixel to construct the 2D profile.
;
; res: resolution (in pixels) of the radial profile 1d, optional value
; Default is 1.
; p2d: 2d profile
; P1D: 1d Profile
; radius: radius corresponding P1D
; / med4q: 2d calculates the profile by taking the median of the 4 quadrants

s=size(imt) & dimx=s[1] & dimy=s[2]
sz=size(imt,/dim)
im=double(imt)

rmin=10
rmax=35

dist_circle,rrange,sz[0]
angr=angarr(sz[0],sz[1])
degr=(360/(2*!pi))*angr + 180 mod 360
ndistang=float(360/drang)
distang=findgen(ndistang)*drang

;calcul de profil 1d
;calculates 1d profile

angavg=fltarr(n_elements(distang))

for i=0L,n_elements(distang)-2 do begin
angrange=where(degr ge distang[i] and degr lt distang[i+1] and rrange ge rmin and rrange lt rmax,nang)
;print,distang[i],distang[i+1],nang,min(degr),max(degr),min(rrange),max(rrange),rmin,rmax
angavg[i]=median(im[angrange],/even)
endfor
;angrange=where(degr ge distang[n_elements(distang)-1] and degr gt distang[0],nangrange)
;bah=fltarr(sz[0],sz[1])

;print,nangrange,distang[n_elements(distang)-1],distang[0]
;angavg[n_elements(distang)-1]=median(im[angrange],/even)


;calculate 2d profile

p2d=dblarr(dimx,dimy)
sind=sort(degr)
p2d[sind]=spline(distang,angavg,degr[sind],/double)

p2d=float(p2d)
if arg_present(p2d) then p2d=float(p2d)
rayon=distang
end
