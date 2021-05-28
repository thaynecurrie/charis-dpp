pro charis_calc_final_contrast,pick=pick,datacube=datacube,databb=databb,noattenmap=noattenmap,$
maskcoord=maskcoord,$
meanadd=meanadd,$
nofinite=nofinite,nozero=nozero,nofilt=nofilt,$
help=help
;,$
;outfile=outfile

;v0.1 - Sept 3 2018 - Does contrast curve calculation for PSF-subbed data.

;Procedure ...
;1. will look in reduc/proc directory for your file.   If it cannot find file will prompt for file name; /pick allows you to find the file.
;2. will check that the file has been flux calibrated and that the star brightness is listed in one of the fits extension headers.
;3. if so, will compute the 'band-integrated' star signal across cube, will assume that this is ~signal in aperture ~averaged ap.
;4. will take each slice, divide by an attenuation map if provided, compute the rms, combine all slices and compute rms for broadband image.

if ((N_PARAMS() eq 0 and ~keyword_set(pick) and ~keyword_set(datacube)) or keyword_set(help)) then begin
print,'5-sigma contrast curve calculations ...'
print,"charis_calc_final_contrast,pick=pick,datacube=datacube,databb=databb,noattenmap=noattenmap,"
print,"maskcoord=maskcoord,meanadd=meanadd,nofinite=nofinite,nozero=nozero,nofilt=nofilt"
;print,"maskcoord=maskcoord,meanadd=meanadd,nofinite=nofinite,nozero=nozero,nofilt=nofilt,outfile=outfile"
print,""
print,"Example: charis_calc_final_contrast,datacube='asdi.fits'"
print,""
print,"***Keywords***"
print,"*datacube - the datacube in the ./reduc/proc directory for which you want to calculate a contrast curve"
print,"*databb - if you want to override the calculation of the wavelength-collapsed image to compute broadband contrast"
print,"*pick - pick the file for contrast curve from GUI"
print,"*noattenmap - do not include an signal loss map (due to annealing by PSF sub. method)"
print,"*meanadd - meancombine the channels for a collapsed cube instead of median"
print,"*nofinite,nozero,nofilt - turn off finite-element corrections, subtraction the zero value, and spatial filtering in the SNR map calculations"
;print,"*outfile - the output file name"

goto,skiptotheend
endif

reducdir='./reduc/proc/'
if (~keyword_set(pick) and keyword_set(datacube)) then begin
testname=file_search(reducdir+datacube,count=testcount)

;*****finding the reduced data cube

if testcount eq 1 then begin ;if you find exactly one file then use it
reduccube=testname
endif else begin ;if either zero files or more than one, prompt to select the cube you want.
reduccube=dialog_pickfile(Title="Select Your Reduced and Calibrated Data Cube")
endelse

endif else begin

reduccube=dialog_pickfile(Title="Select Your Reduced and Calibrated Data Cube")
endelse

;*****the attenuation map****
attenmapname='attenmap_lambda.fits'
attenmapcollapsedname='attenmap.fits'
if ~keyword_set(noattenmap) then begin
testname=file_search(attenmapname,count=counttest)
testname2=file_search(attenmapcollapsedname,count=counttest2)

if counttest eq 1 then begin
attenmap_slice=testname
endif else begin
attenmap_slice=dialog_pickfile(Title="Select the Attenuation Map (Cube) Yielding Pt Source Signal Loss/Wavelength")
endelse

if counttest2 eq 1 then begin
attenmap_col=testname2
endif else begin
attenmap_col=dialog_pickfile(Title="Select the Attenuation Map (Collapsed) Yielding Pt Source Signal Loss/Wavelength")
endelse

attenmap_slice=readfits(attenmap_slice)
attenmap_col=readfits(attenmap_col)
endif


file_mkdir,'final_contrast' ;where you store the output
outputdir='final_contrast/'

;information about the target properties.
;param,'HMAG',hmag,/get,pfname=pfname ;don't need to know this
;param,'EHMAG',ehmag,/get,pfname=pfname ;don't need to know this

;****Now read in cube and extract the information about the star flux

redcube=readfits(reduccube,h1,/exten)
h0=headfits(reduccube)

;loop to read out starflux
get_charis_wvlh,h0,wvlh_test,filtname=filtname0
lambda=1.d-3*wvlh_test
Dtel=charis_get_constant(name='Dtel')
;Dtel=7.9d0
pixscale=charis_get_constant(name='pixscale')
;pixscale=0.0164

dimcube=(size(redcube,/dim))

;nominally, you have a 201,201,21 cube [for low res].  for high-res the number of slices is different. alter accordingly
numchannels=dimcube[2]

;if you want the program to compute a broadband image for you, then do this.   
if ~keyword_set(databb) then begin

if ~keyword_set(meanadd) then begin
redcube_collapsed=median(redcube,dimension=3,/even)
redcube_j=median(redcube[*,*,0:4],dimension=3,/even)
redcube_h=median(redcube[*,*,7:13],dimension=3,/even)
redcube_k=median(redcube[*,*,16:numchannels-2],dimension=3,/even)
;endelse

endif else begin
resistant_mean,redcube,3.0,redcube_collapsed,dimension=3
resistant_mean,redcube[*,*,0:4],3.0,redcube_j,dimension=3
resistant_mean,redcube[*,*,7:13],3.0,redcube_h,dimension=3
resistant_mean,redcube[*,*,16:numchannnels-2],3.0,redcube_k,dimension=3
endelse

;endelse

endif else begin
testbb=file_search(bbimage,count=counttestbb)
if counttestbb eq 0 then testbb=dialog_pickfile(Title="Select Broadband Image")
redcube_collapsed=readfits(testbb,/exten)
redcube_j=median(redcube[*,*,0:4],dimension=3,/even)
redcube_h=median(redcube[*,*,7:13],dimension=3,/even)
redcube_k=median(redcube[*,*,16:numchannels-2],dimension=3,/even)
endelse

;endelse


starflux=fltarr(dimcube[2])
; the dimensions for starflux and fwhm should be the same
fwhm=1.0*(1.d-6*lambda/Dtel)*(180.*3600./!dpi)/pixscale
bbfwhm=median(fwhm,/even)
jfwhm=median(fwhm[0:4],/even)
hfwhm=median(fwhm[7:13],/even)
kfwhm=median(fwhm[16:numchannels-2],/even)
;endelse

for s=0L,dimcube[2]-1 do begin
 sf=string(s)
 tmp=sxpar(h1,'FSTAR_'+strtrim(long(sf),2))
 starflux[s]=tmp
 print,'starflux',starflux[s]
endfor

bbstarflux=median(starflux,/even)
print,"The Star's Median Flux Density Is ...",bbstarflux

jstarflux=median(starflux[0:4],/even)
hstarflux=median(starflux[7:13],/even)
kstarflux=median(starflux[16:numchannels-2],/even)

print,"The Star's Median Flux Density at J, H, and K Is ...",jstarflux,hstarflux,kstarflux


;for now, assume that the planet cal and signal is defined from a 0.0164"/pix scale

;****Loop for wavelength on the cube.   Assumes median combination of image
contrast_slice=fltarr(dimcube[0],dimcube[1],dimcube[2])
contrast_bb=fltarr(dimcube[0],dimcube[1])

;if you want to mask a bright source ...
nmask=0
if keyword_set(maskcoord) then begin
dist_circle,mask,dimcube[0],maskcoord[0],maskcoord[1]
maskobj=where(mask le 1.5*max(fwhm),nmask)

endif

;per wavelength channel
for il=0L,dimcube[2]-1 do begin
slice=redcube[*,*,il]
if keyword_set(maskcoords) then begin
if nmask gt 0 then slice[maskobj]=!values.f_nan
endif
;print,fwhm[il]
charis_snratio_sub,slice,fwhm=fwhm[il],/finite,noisemap=noisemap_slice,/silent,/zero,/filt
if ~keyword_set(noattenmap) then noisemap_slice/=attenmap_slice[*,*,il]
contrast_slice[*,*,il]=5*noisemap_slice/starflux[il]
endfor

;***J,H,K Images***

jimage=redcube_j
himage=redcube_h
kimage=redcube_k
if (keyword_set(maskcoords) and nmask gt 0) then begin
jimage[maskobj]=!values.f_nan
himage[maskobj]=!values.f_nan
kimage[maskobj]=!values.f_nan
endif

charis_snratio_sub,jimage,fwhm=jfwhm,/finite,noisemap=noisemap_j,/silent,/zero,/filt

;himage[maskobj]=!values.f_nan
charis_snratio_sub,himage,fwhm=hfwhm,/finite,noisemap=noisemap_h,/silent,/zero,/filt

;kimage[maskobj]=!values.f_nan
charis_snratio_sub,kimage,fwhm=kfwhm,/finite,noisemap=noisemap_k,/silent,/zero,/filt



if ~keyword_set(noattenmap) then begin


if ~keyword_set(meanadd) then begin
attenmap_j=median(attenmap_slice[*,*,0:4],dimension=3,/even)
attenmap_h=median(attenmap_slice[*,*,7:13],dimension=3,/even)
attenmap_k=median(attenmap_slice[*,*,16:numchannels-2],dimension=3,/even)
;endelse

endif else begin
resistant_mean,attenmap_slice[*,*,0:4],3.0,attenmap_j,dimension=3
resistant_mean,attenmap_slice[*,*,7:13],3.0,attenmap_h,dimension=3
resistant_mean,attenmap_slice[*,*,16:numchannnels-2],3.0,attenmap_k,dimension=3
endelse

noisemap_j/=attenmap_j
noisemap_h/=attenmap_h
noisemap_k/=attenmap_k
endif

contrast_j=5*noisemap_j/jstarflux
contrast_h=5*noisemap_h/hstarflux
contrast_k=5*noisemap_k/kstarflux


;***Broadband Image***
combimage=redcube_collapsed

if (keyword_set(maskcoords) and nmask gt 0) then combimage[maskobj]=!values.f_nan
charis_snratio_sub,combimage,fwhm=bbfwhm,/finite,noisemap=noisemap_collapsed,/silent,/zero,/filt
;writefits,'bah.fits',5*noisemap_collapsed/bbstarflux

if ~keyword_set(noattenmap) then noisemap_collapsed/=attenmap_col

contrast_bb=5*noisemap_collapsed/bbstarflux
writefits,outputdir+'contrast_slice.fits',contrast_slice,h1
writefits,outputdir+'contrast_bb.fits',contrast_bb,h1


;Now, write the 1-D profile out to a text file, a plot, and display

profrad_tc,contrast_bb,1,1,70,p1d=pr_bb,rayon=rbb
writecol,outputdir+'broadband_contrast.txt',rbb*pixscale,pr_bb

profrad_tc,contrast_j,1,1,70,p1d=pr_j,rayon=rj
writecol,outputdir+'j_contrast.txt',rj*pixscale,pr_j

profrad_tc,contrast_h,1,1,70,p1d=pr_h,rayon=rh
writecol,outputdir+'h_contrast.txt',rh*pixscale,pr_h

profrad_tc,contrast_k,1,1,70,p1d=pr_k,rayon=rk
writecol,outputdir+'k_contrast.txt',rk*pixscale,pr_k

prslice=fltarr(n_elements(rbb),dimcube[2])
for il=0L,dimcube[2]-1 do begin
profrad_tc,contrast_slice[*,*,il],1,1,70,p1d=pr,rayon=rbb
;help,pr,prslice
prslice[*,il]=pr
writecol,outputdir+'slice_'+strtrim(il,2)+'_contrast.txt',rbb*pixscale,prslice[*,il]
endfor

rarc=rbb*pixscale

;plot and display
set_plot,'x'
setcolors,/system_variables,/silent
plot,rarc,pr_bb,xrange=[0,1],yrange=[5e-7,1e-3],/ylog,xstyle=1,ystyle=1,$
xtitle='Separation (arc-sec)',ytitle=textoidl(' Contrast ( \Delta F)'),xthick=5,ythick=5,$
charsize=1.5,charthick=1.5,/nodata
loadct,13,/silent

 for ilam=0L,dimcube[2]-1 do begin
 oplot,rarc,prslice[*,ilam],color=(double(ilam+1)/dimcube[2])*256,thick=3
 endfor
 setcolors,/system_variables,/silent
 oplot,rarc,pr_bb,color=!cyan,thick=10



;plot and display
set_plot,'ps'
device,filename=outputdir+'contrast_slice.eps',/encapsulated,/color,bits=8
setcolors,/system_variables,/silent
plot,rarc,pr_bb,xrange=[0,1],yrange=[5e-7,5e-4],/ylog,xstyle=1,ystyle=1,$
xtitle='Separation (arc-sec)',ytitle=textoidl(' Contrast ( \Delta F)'),xthick=5,ythick=5,$
charsize=1.5,charthick=1.5,/nodata
loadct,13,/silent

 for ilam=0L,dimcube[2]-1 do begin
 if ilam eq 4 or ilam eq 11 or ilam eq 18 then begin
 ;oplot,rarc,prslice[*,ilam],color=(double(ilam+1)/dimcube[2])*256,thick=8
 endif else begin
 oplot,rarc,prslice[*,ilam],color=(double(ilam+1)/dimcube[2])*256,thick=3,linestyle=1
 endelse
 endfor
 setcolors,/system_variables,/silent
 oplot,rarc,pr_k,color=(double(19)/dimcube[2])*256,thick=10
 oplot,rarc,pr_h,color=(double(12)/dimcube[2])*256,thick=10
 oplot,rarc,pr_j,color=(double(5)/dimcube[2])*256,thick=10
 oplot,rarc,pr_bb,color=!magenta,thick=15

al_legend,color=[(5./22.)*256,(12/22.)*256,(19./22.)*256,!magenta],['J','H','Ks','Broadband'],/right,$
box=0,linestyle=0,thick=[3,3,3,10],charsize=1.25,charthick=3
device,/close



skiptotheend:
end
