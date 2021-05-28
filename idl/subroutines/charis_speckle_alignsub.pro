function charis_speckle_alignsub,cubein,refslice=refslice,lambda=lambda,$
locs=locs,fft=fft,silent=silent,reverse=reverse

;+
; NAME: 
;       charis_alignspecklesub
; PURPOSE:
;       Rescale image slices such that speckles (and sat spots) are in
;       the same location in each slice. 
;       


;;cube dimensions
sz = double(size(cubein,/dim))
typ = size(cubein,/type)
if n_elements(sz) ne 3 then begin
   message,'You must supply a 3D image cube.',/continue
   return,-1
endif

;;by default use 0th slice as reference and assume H band
if n_elements(refslice) eq 0 then refslice = 0

wvlhs=lambda
if keyword_set(reverse) then scl = lambda/lambda[refslice] else $
   scl = lambda[refslice]/lambda

;if the locations are set then use them to tune the central point from which to scale.
;otherwise, assume middle of the array
if keyword_set(locs) then begin
 ;locs=locs
 if n_elements(locs) gt 8 then locs=locs[*,*,refslice]
 ;;find coordinates of all finite points in the reference slice
 goodmap=where(finite(cubein[*,*,refslice]))
 goodinds=array_indices(cubein[*,*,refslice],goodmap)
 c0 = (total(locs,2)/4) # (fltarr(n_elements(goodmap))+1.)
endif else begin
 ;;find coordinates of all finite points in the reference slice
 goodmap=where(finite(cubein[*,*,refslice]))
 ;goodmap=where((cubein[*,*,refslice]))
; print,n_elements(goodmap)
 goodinds=array_indices(cubein[*,*,refslice],goodmap)

 c0=[100,100] # (fltarr(n_elements(goodmap))+1.)
;c0=[200,200] # (fltarr(n_elements(goodmap))+1.)
endelse

coordsp = cv_coord(from_rect=goodinds - c0,/to_polar)

if ~keyword_set(fft) then begin
out=make_array(sz,type=typ)
for j=0,sz[2]-1 do begin
      if j ne refslice then begin
         tmp = make_array(sz[0:1],type=typ) + !values.f_nan
         coordsj = cv_coord(from_polar=[coordsp[0,*],coordsp[1,*]/scl[j]],/to_rect) + c0
         tmp[goodmap] =  interpolate(cubein[*,*,j],coordsj[0,*],coordsj[1,*],cubic=-0.5,missing=!values.f_nan)
         ;tmp[goodmap] =  interpol(cubein[*,*,j],coordsj[0,*],coordsj[1,*],/lsquadratic)
      endif else tmp = cubein[*,*,j]
      out[*,*,j] = tmp
endfor

endif else begin

;goto,skipoverme

;FFT method
redoparms = 1
      parms = create_struct('padinit',0L,$
                            'padinit2',0L,$
                            'padfft',0L,$
                            'padfft2',0L,$
                            'scf',0d,$
                            'dimx',0L,$
                            'dimy',0L,$
                            'scx',0d,$
                            'scy',0d,$
                            'meilprec',0d)
      parmsarr = replicate(parms,sz[2])
;stop
;;remove NaNs from input
   in = make_array(sz[0] + (sz[0] mod 2.), sz[1] + (sz[1] mod 2.), sz[2])
   in[0:sz[0]-1,0:sz[1]-1,*] = cubein
   badmap = where(~FINITE(in))
   in[badmap] = 0.

   ;badmap = where(in eq 0)
   ;in[badmap] = !values.f_nan

  ;;apply scaling
   out = make_array(size(in,/dim),type=size(cubin,/type))
   for j = 0,sz[2]-1 do begin
      if j ne refslice then begin
         out[*,*,j] = fftscale(in[*,*,j],scl[j],scl[j],1e-6,silent=silent)
         ;if redoparms then begin
         ;   out[*,*,j] = fftscale(in[*,*,j],scl[j],scl[j],1e-7,parms=tmp,silent=silent)
         ;   parmsarr[j] = temporary(tmp)
         ;endif else out[*,*,j] = fftscale(in[*,*,j],scl[j],scl[j],1e-7,parms=parmsarr[j],silent=silent)
      endif else out[*,*,j] = in[*,*,j]
   endfor

   ;;put NaNs back
   out[badmap] = !VALUES.F_NAN
   out = out[0:sz[0]-1,0:sz[1]-1,*]
endelse
;skipoverme:
;goto,skipoverme

return,out
end

