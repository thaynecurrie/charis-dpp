pro charis_est_lum_bd,d=d,k=k,spt=spt,bck=bck,setmbol=setmbol,help=help

if (keyword_set(help)) then begin
print,"estimates the luminosity of a BD comparison object given a K band magnitude and a distance"
print,""
print,"charis_est_lum_bd,d=d,k=k,spt=spt,bck=bck"
print,""
print,"***Keywords*"
print,""
print,"*d - distance in parsecs"
print,"*k - kband magnitude"
print,"*spt - numerical spectral type (10=L0, 15=L5, etc.)"
print,"*bck - manually enter the correction"
print,"*setmbol - manually enter the bolometric magnitude"

goto,skiptotheend
endif

;estimates the luminosity of a BD comparison object given a K band magnitude and a distance

;from Golimowski et al. 2004
;spt=10 --> L0
c0=3.9257
c1=-3.8338d-1
c2=5.3597d-2
c3=-2.655d-3
c4=4.0859d-5
mkcor=c0+c1*spt+c2*spt^2.+c3*spt^3.+c4*spt^4.   

if keyword_set(bck) then mkcor=bck
mbol=mkcor+k-5*alog10(d/10.)

if keyword_set(setmbol) then mbol=setmbol

;now take the difference with the sun mbol = +4.74

logl=alog10(10^(-0.4*(mbol-4.74)))
print,'logl',logl

skiptotheend:
end

