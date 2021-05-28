pro charis_calcsptgravind,lambda=lambda,fnu=fnu,flam=flam,outputspt=outputspt,outputgrav=outputgrav,mode=mode,justind=justind

;Calculates a rough approximation of the Allers & Liu 2013 spectral type and gravity indices
;****Currently works only for Low-Resolution Mode!!!   Alpha-Version!

;If you do not give a wavelength array, assume you are using low-res data
if ~keyword_set(lambda) then begin
get_charis_wvlh,dum,wavelengths,manual='lowres'
lambda=1.d-3*wavelengths
endif

;if you do not set the mode, then assume low-resolution.   Make this only version for now
if ~keyword_set(mode) then begin

;if the data are in F-nu units then convert to pseudo-Flambda, if in F-lambda, leave alone
if (~keyword_set(flam) and keyword_set(fnu)) then begin
flux=fnu/lambda^2
endif 
if keyword_set(flam) then begin
flux=flam
endif
if (~keyword_set(flam) and ~keyword_set(fnu)) then begin
spectrum=dialog_pickfile(Title="Select Your Output Spectrum in Fnu Units")
readcol,spectrum,lambda,fnu,efnu
flux=fnu/lambda^2.
endif

;****Calculate indices

;SpT
;H20 ~avg (1.522-1.575)/avg (1.471-1.522)
;H20-D ~ avg (1.931-1.999)/2.067
;H20-1 ~ 1.329/avg(1.284-1.328)
;H20-2 ~ 2.067/2.139

;***index_h20

index=(flux[8]+flux[9])/(flux[7]+flux[8])
;index=1.1
c0=-83.5437 & c1=169.388 & c2=-104.424 & c3=24.0476
spth20=c0+c1*index+c2*index^2.+c3*index^3.
indexh20=index

;***index_h20d

index=0.5*(flux[15]+flux[16])/(flux[17])
c0=79.4477 & c1= -202.245 & c2=229.884 & c3=-97.230
spth20d=c0+c1*index+c2*index^2.+c3*index^3.
indexh20d=index

;***index_h201

index=flux[4]/(0.5*(flux[3]+flux[4]))
;c0=23.1927 & c1= 39.3513 & c2=-80.7404 & c3=28.5982
c0=12.1927 & c1= 39.3513 & c2=-80.7404 & c3=28.5982
spth201=c0+c1*index+c2*index^2.+c3*index^3.
indexh201=index

;***index_h202

;index=(0.75*flux[17]+0.25*flux[16])/(0.75*flux[18]+0.25*flux[19])
index=(flux[17])/(flux[18])
c0=10.8822 & c1=55.4580 & c2=-97.8144 & c3=37.5013
spth202=c0+c1*index+c2*index^2.+c3*index^3.
indexh202=index


;Gravity
;H2K ~ 2.14/2.21 ;or avg(2.14-2.213)/avg(2.213-2.29)
;H-cont ~ lam_c=1.575, lam1=1.471, lam2=1.686

index_h2k=(flux[18]+flux[19])/(flux[19]+flux[20])

lamc=lambda[9]
lam1=lambda[7]
lam2=lambda[11]
flux2=flux[11]
flux1=flux[7]
fluxc=flux[9]

index_hcont=((lamc-lam1)*flux2/(lam2-lam1)+(lam2-lamc)*flux1/(lam2-lam1))/fluxc


endif

;adjustments
;spth20-=1
;spth20d-=3
;spth202-=1
outputspt=[spth20,spth20d,spth201,spth202]
if keyword_set(justind) then begin
print,'indexs',indexh20,spth20
outputspt=[indexh20,indexh20d,indexh201,indexh202]
endif
outputgrav=[index_h2k,index_hcont]
print,'Spectral Type ',outputspt
print,'Gravity Index ',outputgrav

end


