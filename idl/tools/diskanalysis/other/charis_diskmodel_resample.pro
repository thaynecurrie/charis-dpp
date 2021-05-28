function charis_diskmodel_resample,model,inputpixscale=inputpixscale,help=help

;takes a disk model of some pixel scale, returns resampled disk model in CHARIS units

if (N_PARAMS() eq 0 or keyword_set(help))) then begin
print,"charis_diskmodel_resample,model,inputpixscale=inputpixscale"
print,""
print,"***Keywords"
print,"*inputpixscale - input pixel scale of the image to be resampled"
goto,skiptotheend
endif

;read input disk model
inputdiskmodel=readfits(model)

modeldim=[(size(inputdiskmodel,/dim))[0]/2,(size(inputdiskmodel,/dim))[1]/2]

if modeldim[0] lt 201 then begin
inputdiskmodel=subarr(inputdiskmodel,400,/zeroout)
modeldim=[(size(inputdiskmodel,/dim))[0]/2,(size(inputdiskmodel,/dim))[1]/2]

endif


goodmap=where(finite(inputdiskmodel))
goodinds=array_indices(inputdiskmodel,goodmap)

c0=[(size(inputdiskmodel,/dim))[0]/2,(size(inputdiskmodel,/dim))[1]/2] # (fltarr(n_elements(goodmap))+1.)
coordsp = cv_coord(from_rect=goodinds - c0,/to_polar)

;if you don't set the input pixel scale, assume you are reading in
; a HiCIAO SEEDS disk model
if ~keyword_set(inputpixscale) then inputpixscale =0.0095


charis_size=[201,201]
pixscale=charis_get_constant(name='pixscale')
scl=inputpixscale/pixscale

coordsj = cv_coord(from_polar=[coordsp[0,*],coordsp[1,*]/scl],/to_rect) + c0

outputdiskmodel=fltarr((size(inputdiskmodel,/dim))[0],(size(inputdiskmodel,/dim))[1])

tmp = make_array((size(inputdiskmodel,/dim))[0:1],type=typ) + !values.f_nan
tmp[goodmap] =  interpolate(inputdiskmodel,coordsj[0,*],coordsj[1,*],cubic=-0.5,missing=!values.f_nan)

outputdiskmodel[goodmap]=tmp
outputdiskmodel=subarr(outputdiskmodel,charis_size)

return,outputdiskmodel

skiptotheend:

end
