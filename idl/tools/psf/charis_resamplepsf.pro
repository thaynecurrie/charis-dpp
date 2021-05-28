pro charis_resamplepsf,inputpsf,outputpsf,scalefac=scalefac,help=help

;Resamples a PSF from some pixel scale to the CHARIS pixel scale
;useful for resampling disk models from another instrument to CHARIS

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"charis_resamplepsf,inputpsf,outputpsf,scalefac=scalefac"
print,""
print,"***Keywords"
print,"*scalefac= the scaling factor for resampling"

goto,skiptotheend
endif

if ~keyword_set(scalefac) then scalefac = 1.0

;assume square arrays
sz=(size(inputpsf,/dim))[1]

goodmap=where(finite(inputpsf))
goodinds=array_indices(inputpsf,goodmap)

c0=[(size(inputpsf,/dim))[0]/2,(size(inputpsf,/dim))[1]/2] # (fltarr(n_elements(goodmap))+1.)
coordsp = cv_coord(from_rect=goodinds - c0,/to_polar)

scl=scalefac

coordsj = cv_coord(from_polar=[coordsp[0,*],coordsp[1,*]*scl],/to_rect) + c0
;coordsj = cv_coord(from_polar=[coordsp[0,*],coordsp[1,*]],/to_rect) + c0

outputpsf=fltarr((size(inputpsf,/dim))[0],(size(inputpsf,/dim))[1])

tmp = make_array((size(inputpsf,/dim))[0:1],type=typ) + !values.f_nan

tmp[goodmap] =  interpolate(inputpsf,coordsj[0,*],coordsj[1,*],cubic=-0.5,missing=!values.f_nan)

outputpsf[goodmap]=tmp

writefits,'tmp.fits',tmp
outputpsf=subarr(outputpsf,sz)

skiptotheend:
end
