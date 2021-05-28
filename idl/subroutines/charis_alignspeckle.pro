pro charis_alignspeckle,image,h0,h1,fft=fft,reverse=reverse,nolocs=nolocs,$
refslice=refslice
;built and adapted from gpi_specklealign
;Align speckles in an image prior to doing stuff



if ~keyword_set(refslice) then begin

refslice=0
k=refslice
endif else begin
;print,refslice
;stop
k=refslice
endelse

;filter=sxpar(h0,'FILTERNAME')
get_charis_wvlh,h0,wvlh_charis

;If you have the satellite locations determined, then use them
if ~keyword_set(nolocs) then begin
;sat=sxpar(h1,'SATMASK')
;goodcode=hex2bin(sat,(size(image,/dim))[2])
;good=long(where(goodcode eq 1))

;goodcode=ulong(make_array((size(image,/dim))[2]))
;goodcode[*]=1
;help,goodcode
;print,'goodcode',goodcode
;stop
;good=long(where(goodcode eq 1))
;print,good
cens=fltarr(2,4,(size(image,/dim))[2])
;print,'GOOD!',n_elements(good),(size(image,/dim))[2]

for s=0,(size(image,/dim))[2]-1 do begin
;for s=0,n_elements(good)-1 do begin
 for j=0,3 do begin
  tmp=fltarr(2)+!values.f_nan
   sf=string(s)
   tmp= sxpar(h1,'SATS'+strtrim(long(sf),2)+'_'+strtrim(j,2))
;**total, ugly hack-tastic workaround since I can't seem to get IDL to treat the string properly otherwise...
   tmpz=strsplit(tmp,' ',/extract)

   tmp=[tmpz[0],tmpz[1]]
;  reads,sxpar('SATS'+strtrim(long([s]),2)+'_'+strtrim(j,2),h1),
   cens[*,j,s]=tmp
 endfor
endfor

tmp=where(~finite(cens[*,*,k]),ct)
if ct ne 0 then begin
print,'stop!'
stop
endif


endif

badvalue=where(image eq 0)
;image[badvalue]=!values.f_nan

if ~keyword_set(nolocs) then begin
if ~keyword_set(reverse) then begin

if ~keyword_set(fft) then begin
im2=charis_speckle_alignsub(image,refslice=k,lambda=wvlh_charis,locs=cens[*,*,k])
endif else begin
im2=charis_speckle_alignsub(image,refslice=k,lambda=wvlh_charis,locs=cens[*,*,k],/fft)
endelse

endif else begin


if ~keyword_set(fft) then begin
im2=charis_speckle_alignsub(image,refslice=k,lambda=wvlh_charis,locs=cens[*,*,k],/reverse)
endif else begin
im2=charis_speckle_alignsub(image,refslice=k,lambda=wvlh_charis,locs=cens[*,*,k],/reverse,/fft)
endelse


endelse

endif else begin

if ~keyword_set(reverse) then begin


if ~keyword_set(fft) then begin
im2=charis_speckle_alignsub(image,refslice=k,lambda=wvlh_charis)
endif else begin
im2=charis_speckle_alignsub(image,refslice=k,lambda=wvlh_charis,/fft)
endelse

endif else begin

if ~keyword_set(fft) then begin
im2=charis_speckle_alignsub(image,refslice=k,lambda=wvlh_charis,/reverse)
endif else begin
im2=charis_speckle_alignsub(image,refslice=k,lambda=wvlh_charis,/reverse,/fft)
endelse

endelse

endelse

image=im2








end
