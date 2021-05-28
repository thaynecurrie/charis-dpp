function charis_combinefits,flist,mean=mean,dir=directory,cube=cube,nwght=nwght,varcube=varcube
;takes median of fits images or a robust mean


nfiles=n_elements(flist)

if ~keyword_set(cube) then begin
bah=headfits(directory+flist[0],/silent)
fsize=sxpar(bah,'naxis1')
image_in=fltarr(fsize,fsize,nfiles)

for i=0L,nfiles-1 do begin
im=readfits(directory+flist[i],h1,/silent)
image_in[*,*,i]=im
endfor

if ~keyword_set(mean) then begin
image=median(image_in,dimension=3,/even)
endif else begin
sigcut=3.
resistant_mean,image_in,sigcut,image,dimension=3
endelse

endif else begin
;if a datacube.  

test=readfits(directory+flist[0],h1,/silent,/exten)
sz=size(test)
fsize=sz[1]
n_lambda=sz[3]

image_in=fltarr(fsize,fsize,n_lambda,nfiles)

for i=0L,nfiles-1 do begin
im=readfits(directory+flist[i],h1,/exten,/silent)
;hyper-cube
image_in[*,*,*,i]=im
endfor

;noise weighting
if keyword_set(nwght) then begin

i1=fltarr(fsize,fsize,n_lambda)
i2=fltarr(fsize,fsize,n_lambda)
i1r=fltarr(fsize,fsize)
i2r=fltarr(fsize,fsize)
i1[*,*,*]=0
i2[*,*,*]=0

;loop file number

for i=0L,nfiles-1 do begin
 for il=0L,n_lambda-1 do begin
  i1r[*,*]=0
  i2r[*,*]=0
  duh1=image_in[*,*,il,i]
  duh2=varcube[*,*,il,i]
  good=where(finite(duh1) eq 1 and finite(duh2) eq 1)
  i1r[good]=1./duh2[good]
  i1[*,*,il]+=i1r
  i2r[good]=(duh1[good]/duh2[good])
  i2[*,*,il]+=i2r
 endfor

endfor

i1f=1./i1
image=i1f*i2

endif else begin

;no noise-weighting
if ~keyword_set(mean) then begin
image=median(image_in,dimension=4,/even)
endif else begin
sigcut=3.
resistant_mean,image_in,sigcut,image,dimension=4
endelse

endelse

endelse

return,image
end
