pro charis_check_for_detection,datacube=datacube,threshold=threshold,excl=excl,rangeconstant=rangeconstant,colimage=colimage,snrmapout=snrmapout,outfile=outfile,reducpath=reducpath,verbose=verbose,numdet=numdet

numdet=0


;takes a PSF subtracted data cube, auto-checks the data cube for a statistically significant detection
if ~keyword_set(datacube) and ~keyword_set(colimage) then begin
print,"charis_check_for_detection: checks the data cube for a statistically significant detection"
print,''
print,"*Calling Sequence*"
print,'charis_check_for_detection,datacube=datacube,threshhold=threshhold,reducpath=reducpath,colimage=colimage,snrmap=snrmap,outfile=outfile
print,"Example: charis_check_for_detection,datacube='adi.fits'
print,''
print,"****Important Keywords****"
print,"datacube - name of the PSF-subtracted data cube in which you want to search for a detection"
print,"colimage - name of the PSF-subtracted, wavelength-collapsed image [either this or datacube must be set]"
print,"threshold - threshold in sigma to trigger a detection (default = 5.0)"
print,"excl - the exclusion threshold for pixel mask after a detection (to avoid double-counting a bright detection"
print,"snrmapout - output the SNR map to this file [default = not set]"
print,"reducpath - 0= ./reduc/proc/, 1=current path"
print,"outfile - output file to contain the position and rough contrast of detections [not throughput corrected]"

goto,skiptotheend

endif

if ~keyword_set(reducpath) then reducpath=0
case reducpath of 
  0: begin

     pathtocube='./reduc/proc/'
     end
  else: begin

     pathtocube='./'
     end
endcase

if keyword_set(colimage) then begin

psfsubimage=readfits(pathtocube+colimage,h1,ext=1)
h0=headfits(pathtocube+colimage)

endif else begin
psfsubname=(strsplit(datacube,'.',/extract))[0]+'_collapsed.fits'
psfsubimage=readfits(pathtocube+psfsubname,h1,ext=1)

endelse

;pull the wavelength array

;Information about Data
get_charis_wvlh,h1,wvlh_charis
;wvlh_charis*=1d-3
med_wvlh=median(wvlh_charis,/even)
Dtel=charis_get_constant(name='Dtel') ;7.9d0 ;visible pupil for SCExAO
;pixel scale 16.2 mas/pixel
pixscale=charis_get_constant(name='pixscale') ;0.0162 nominally
fwhm_col=1.0*(1.d-6*(med_wvlh*1d-3)/Dtel)*(180.*3600./!dpi)/pixscale


;Information about Reduction
;if this is PSF subtracted data, then set threshold nominally at 5-sigma
psfsubcheck1=sxpar(h1,'ZERO',count=checkcount1)
psfsubcheck2=sxpar(h1,'MEANADD',count=checkcount2)

;LOCI or KLIP?
rminloci=sxpar(h1,'LOCI_RMI',count=ctlociri)
rmaxloci=sxpar(h1,'LOCI_RMA',count=ctlocirout)
rminklip=sxpar(h1,'KLIP_RMI',count=ctklipri)
rmaxklip=sxpar(h1,'KLIP_RMA',count=ctkliprout)

if ctlociri ge 1 then begin
rmin=rminloci
rmax=rmaxloci
endif else begin
rmin=rminklip
rmax=rmaxklip
endelse

;if you
if ~keyword_set(threshold) then begin
if (checkcount1 ge 1 or checkcount2 ge 1) then begin
;set threshold at 5-sigma for PSF subbed data
threshold=5.0

endif else begin
;set threshold at 10-sigma for non-PSF subbed data
threshold=10.0

endelse
endif

if ~keyword_set(excl) then excl=2
if ~keyword_set(rangeconstant) then rangeconstant = 1.5

charis_snratio_sub,psfsubimage,fwhm=fwhm_col,/finite,zero=1,snrmap=snrmapcol
;charis_snratio_sub,psfsubimage,fwhm=fwhm_col,/finite,zero=1,snrmap=snrmapcol,filter=1,expfilt=3
snrmapout=snrmapcol

dist_circle,g,201

;Now, depending on rin and rout for the reduction, find peaks that are more than x-sigma significant
detectmap=fltarr(201,201)
good=where(snrmapcol ge threshold and g ge 1*rangeconstant*fwhm_col+rmin and g le -1*rangeconstant*fwhm_col+rmax,ngood) ;setting value to 1.5 typically is safest

if keyword_set(verbose) then writefits,'sn.fits',snrmapcol

;No Detections Found
if ngood eq 0 then begin
print,"No ",strtrim(float(threshold),2)+"-sigma sources detected between ",strtrim(float(rmin),2)," and ",strtrim(float(rmax),2)
numdet=0
goto,skiptotheend
endif

;Detections Found
detectmap[good]=1
;print,snrmapcol[good]
;print,good mod 201, good/201

;Order the signficant pixels by detection significance
good=reverse(good(sort(snrmapcol[good])))
;print,good,snrmapcol[good]

if keyword_set(verbose) then writefits,'detectmap.fits',detectmap

;Loop the significant pixels, set to zero any pixels within a 1.5xFWHM-wide aperture scaled by SNR of 'detection'
for i=0L,n_elements(good)-1 do begin

if detectmap[good[i]] eq 1 then begin
goodposition=[good[i] mod 201,good[i]/201]
dist_circle,g,201,goodposition[0],goodposition[1]
bad=where(g gt 1d-9 and g le (excl+1.5*alog10(snrmapcol[good[i]]/3.)*fwhm_col) and detectmap eq 1,nbad)
if nbad gt 0 then detectmap[bad]=0
endif

endfor

gooddet=where(detectmap eq 1,ngooddet)
;No Detections Found
if ngooddet eq 0 then begin
print,"No ",strtrim(float(threshold),2)+"-sigma sources detected between ",strtrim(float(rmin),2)," and ",strtrim(float(rmax),2)
goto,skiptotheend
endif

xdetposition=fltarr(ngooddet) ;xpos
ydetposition=fltarr(ngooddet) ;ypos
fluxdet=fltarr(ngooddet) ;flux

for i=0L,ngooddet-1 do begin
;initial positions
xdet=gooddet[i] mod 201
ydet=gooddet[i]/201
cntrd,psfsubimage,xdet,ydet,xout,yout,fwhm_col
xdetposition[i]=xout
ydetposition[i]=yout
aper,psfsubimage,xout,yout,flux,eflux,sky,skyerr,1.0,0.5*fwhm_col,[2,6]*fwhm_col,[0,0],/flux,/exact,/nan,/silent,setskyval=0
fluxdet[i]=flux
endfor

;now divide by star flux

;test for flux calibration
contrast=replicate(-9,ngooddet)
test=sxpar(h1,'FSTAR_0',count=starcount)
if starcount ge 1 then begin
starflux=fltarr(n_elements(wvlh_charis))
for ir=0L,n_elements(wvlh_charis)-1 do begin
starflux[ir]=sxpar(h1,'FSTAR_'+strtrim(string(ir),2))
endfor
medstarflux=median(starflux,/even)
contrast=alog10(fluxdet/medstarflux)
endif

for i=0L,ngooddet-1 do begin
print,"Detection Number ",strtrim(i+1,2)," At ",xdetposition[i],ydetposition[i]," with SNR ",snrmapcol[gooddet[i]]," and Contrast ",contrast[i]
endfor

good=where(xdetposition gt 0 and ydetposition gt 0)

if keyword_set(outfile) then $
writecol,outfile,xdetposition[good],ydetposition[good],contrast[good]

numdet=ngooddet
skiptotheend:
end
