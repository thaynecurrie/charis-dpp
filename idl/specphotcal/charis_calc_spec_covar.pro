pro charis_calc_spec_covar,data_cube,rho=rho,$
prad=prad,filt=filt,normnoise=normnoise,width=width,step=step,$
maxfit=maxfit,minfit=minfit,$
maskresspot=maskresspot,maskcoord=maskcoord,fixpos=fixpos,excval=excval,outparams=outparams,help=help

;07-21-2021 - code is now documented
;NOTE: CODE UNDER DEVELOPMENT/NOT SUPPORTED BY PIPELINE YET and NOT DOCUMENTED

;code to calculate the spectral covariance
;03-31-2018 - rough try

;Keywords
;rho = set the angular separation (in lambda/D units at 1.6 microns)
;prad = subtract a radial profile of the image slice first (can be helpful?)
;normnoise = normalize by the robust standard deviation of the noise first (recommended.)
;step = step size in rho*(lam[i]-lam[j])/lam_c to compute an averaged psi[i,j]
;maskcoord= mask a point source here.

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,'charis_calc_spec_covar: calculates the spectral covariance of a PSF subtracted data cube'
print,'written by T. Currie (2018)'
print,'**Calling Sequence**'
print,'charis_calc_spec_covar,data_cube,rho=rho,prad=prad,filt=filt,normnoise=normnoise,width=width,step=step,'
print,'maxfit=maxfit,minfit=minfit,maskresspot=maskresspot,maskcoord=maskcoord,fixpos=fixpos,excval=excval,outparams=outparams'
print,''
print,'Example:'
print,"charis_calc_spec_covar,'alocihd1160.fits',rho=18.9,maskcoord=[144,78],/outparams"
print,''
print,"***Important Keywords***"
print,'data_cube -- your PSF subtracted data cube.  (You have to give the full path)'
print,'*rho - angular separation in lambda/D units at 1.6 microns where the covariance is evaluated (e.g. rho=20 is ~0.86")'
print,'*prad - subtract a radial profile of the image slice first (sometimes helpful) '
print,'*normnoise - normalize by the robust standard deviation of the noise (usually helps)'
print,'*step - step size in rho*(lam[i]-lam[j])/lam_c to compute an averaged psi[i,j] (usually dont change this)'
print,'*maskcoord - x,y coordinates of a point source to mask (ALWAYS do this if you have a bright companion)'
print,''
print,'*outparams - do you want to save the covariance function to a file? (REQUIRED for charis_empbdplanspec.pro!!!!)'

goto,skiptoend
endif

;Procedure ...

;set the angle of interest
if ~keyword_set(rho) then rho = 10
if ~keyword_set(width) then width = 0.5
if ~keyword_set(step) then step=0.125
if ~keyword_set(excval) then excval=3.

;;*optional - normalize image by standard deviation in annulus

;***reading in file***
data_cube_in=readfits(data_cube,hcube1,ext=1)
hcube0=headfits(data_cube)

;***preliminaries
;** setting up, wavelength array values and size of array, FWHM, and image region in units of lambda/D
;lam_dim=(size(data_cube_in,/dim))[0]
get_charis_wvlh,hcube0,wavelengths
lambda=wavelengths*1.d-3

lambda_med=median(lambda,/even)


lam_dim=n_elements(lambda)

;Telescope Diameter for Subaru
Dtel=charis_get_constant(name='Dtel') ;7.9d0 ;visible pupil for SCExAO
pixscale=charis_get_constant(name='pixscale') ;0.0162 nominally

;FWHM of a point source
fwhm=1.*(1.d-6*lambda/Dtel)*(180.*3600./!dpi)/pixscale
med_fwhm=median(fwhm,/even) ; the median FWHM over the bandpass

;array of distances ... in pixel units
print,(size(data_cube_in,/dim))[0]
dist_circle,sep_lambd,(size(data_cube_in,/dim))[0]
sep_lambd/=med_fwhm


;an annulus of width 'width' lambda/D centered on rho lambda/D
roi=where(abs(sep_lambd - rho) le width)  

gah=fltarr((size(data_cube_in,/dim))[0],(size(data_cube_in,/dim))[0])
gah=data_cube_in[*,*,10]
gah[roi]=1

;mask point sources, subtract radial profile, normalize by stdev of pixels
if keyword_set(maskcoord) then begin
xpto=maskcoord[0]
ypto=maskcoord[1]
data_cube_col=median(data_cube_in,/even,dimension=3)

if ~keyword_set(fixpos) then begin
gcntrd,data_cube_col,xpto,ypto,xpt,ypt,med_fwhm
endif else begin
xpt=xpto
ypt=ypto
endelse

dist_circle,radpt,[(size(data_cube_in,/dim))[0],(size(data_cube_in,/dim))[0]],xpt,ypt
nopointsource=where(radpt gt excval*med_fwhm,complement=pointsource,nbad)

endif
for i=0L,lam_dim-1 do begin
gah=data_cube_in[*,*,i]

if keyword_set(maskcoord) then begin
gah[pointsource]=!values.f_nan
endif

;radial-profile subbing
if keyword_set(prad) then begin
profrad_tc,gah,1,1,60,p2d=pr
gah-=pr
endif

if keyword_set(filt) then begin
gah-=filter_image(gah,median=10*med_fwhm)
endif

if keyword_set(normnoise) then begin
profrad_tc,abs(gah)/.6745,1,1,60,p2d=prnoise
gah/=prnoise
endif

if keyword_set(maskresspot) then begin
resspot=where(gah gt 10.,nresspot)
if nresspot gt 1 then gah[resspot]=!values.f_nan
endif
data_cube_in[*,*,i]=gah
endfor

if keyword_set(maskcoord) then writefits,'masked_cube.fits',median(data_cube_in,dimension=3,/even)

;1. compute the average spectral correlation in an annulus

spec_cor=fltarr(lam_dim,lam_dim)
lam_ratio=fltarr(lam_dim,lam_dim)

for ii=0L,lam_dim-1 do begin
 for ij=0L,lam_dim-1 do begin
  slice_i=data_cube_in[*,*,ii]
  slice_j=data_cube_in[*,*,ij]
  slice_i=slice_i[roi]
  slice_j=slice_j[roi]

  lam_ratio[ii,ij]=(lambda[ii]-lambda[ij])/lambda_med
  num=mean(slice_i*slice_j,/nan)
  denom=sqrt(mean(slice_i^2.,/nan)*mean(slice_j^2.,/nan))
  spec_cor[ii,ij]=abs(double(num)/double(denom))
  ;spec_cor[ii,ij]=double(num)/double(denom)
  if (spec_cor[ii,ij] gt 0.3 and spec_cor[ii,ij] lt 0.9) then print,'high cor ',ii,ij
  if (ii ge 17 and ii gt ij) then print,'k cor',spec_cor[ii,ij],ii,ij
 endfor
endfor

print,max(abs(lam_ratio))

;*******Fiducial Model ******
;a mock-up model
arho=0.1
alam=0.3
adel=0.6
sigrho=0.38
siglam=0.75

;arho=0.38
;alam=0.18
;adel=0.44
;sigrho=0.82
;siglam=0.32

goto,skipoverse
arho=0.42
alam=0.19
adel=(1-arho-alam)
sigrho=0.68
siglam=0.35
skipoverse:

lam_ratio_sort=lam_ratio(sort(lam_ratio))
lam_ratio_sortf=lam_ratio_sort(where(lam_ratio_sort ne 0))
mod_speccor=arho*exp(-0.5*(rho*lam_ratio/sigrho)^2.)+alam*exp(-0.5*(lam_ratio/siglam)^2.)+adel*0
diag_el=findgen(lam_dim)+findgen(lam_dim)*lam_dim
mod_speccor[diag_el]+=adel
mod_speccorf=mod_speccor(where(lam_ratio_sort ne 0))
;*********



;******sample the distribution of rho*lam_ratio vs. spec_cor evenly
;try a sampled range in separations
lam_ratio_step=step
if ~keyword_set(maxfit) then begin
maxfit=max(round(lam_ratio*rho))
minfit=min(round(lam_ratio*rho))
endif 
lam_ratio_num=(maxfit-minfit)/float(lam_ratio_step)
;lam_ratio_num=(max(round(lam_ratio*rho))-min(round(lam_ratio*rho)))/float(lam_ratio_step)
;rholam_ratio_spline=findgen(lam_ratio_num+1)*lam_ratio_step+min(round(rho*lam_ratio))
rholam_ratio_spline=findgen(lam_ratio_num+1)*lam_ratio_step+minfit

spec_cor_avg=fltarr(lam_ratio_num)
espec_cor_avg=fltarr(lam_ratio_num)  ;not used for model fitting but for comparison purposes


for i=0L,lam_ratio_num-1 do begin
good=where(abs(lam_ratio*rho-rholam_ratio_spline[i]) le lam_ratio_step)

;straight average
spec_cor_avg[i]=mean(spec_cor[good])
;robust average
;resistant_mean,(spec_cor)[good],5,mean,sigma
;spec_cor_avg[i]=mean
;spec_cor_avg[i]=median(spec_cor[good],/even)

espec_cor_avg[i]=stddev((spec_cor)[good],/nan)

endfor
espec_cor_avg[*]=median(espec_cor_avg,/even)

;*****Setting up MP-FIT*****
;Starting values for L-M fit to speckle covariance model
p0=[arho,alam,sigrho,siglam,rho]

pi = replicate({value:0.D, fixed:0, limited:[0,0], $ 
 limits:[0.D,0]}, 5)
;help,pi
;stop
;set boundaries, initialize
pi(0:3).limited(0)=1
pi(0:3).limited(1)=1

;set the boundaries for A_rho and A_lambda between 0.1 and 0.8, adjust if necessary
pi(0).limits(0)=0.0005
pi(0).limits(1)=0.9
pi(1).limits(0)=0.0005
pi(1).limits(1)=0.9

;unfortunately cannot do this.
;pi(0).limits(1)=0.99-pi(1).limits(1)

;sigma_rho = the correlated noise part due to speckle noise
pi(2).limits(0)=0.01
pi(2).limits(1)=0.9
;pi(2).limits(1)=2.9
;sigma_lambda = noise due to spectral extraction
pi(3).limits(0)=0.025
pi(3).limits(1)=2.9

;fix the value for rho
pi(4).fixed=1
;set starting values
pi(*).value=[arho,alam,sigrho,siglam,rho]

;************

;*****MP-FIT Execution***
;do the fitting where i =/ j and where you have enough bins that the spec_cov estimate has some meaning.
;good=where(finite(espec_cor_avg) ne 0 and spec_cor_avg ne 1)

good=where(finite(espec_cor_avg) ne 0 and spec_cor_avg ne 1 and abs(rholam_ratio_spline) ge step)
good_zero=where(rholam_ratio_spline eq 0)
good_inclzero=where((finite(espec_cor_avg) ne 0 and spec_cor_avg ne 1 and abs(rholam_ratio_spline) ge step) or rholam_ratio_spline eq 0)
;zeroval=where(finite(espec_cor_avg) ne 0 and spec_cor_avg ne 1 and abs(rholam_ratio_spline) eq 0,)

;perform L-M minimization via mpfitfun.pro; 
;-define a function GBFUNC (residing in the same subdirectory), which expresses the functional form of spectral covariance from Greco & Brandt (2016)
;-fit the interpolated data to the function, assuming a weighting proportional to the spectral covariance
;-outputs A_rho, A_lambda, sigma_rho, sigma_lambda, and rho (which is held fixed)

;***
fitresult=mpfitfun('GBFUNC',rholam_ratio_spline[good]/rho,spec_cor_avg[good],espec_cor_avg[good],weights=spec_cor_avg[good],p0,parinfo=pi,/silent,/nan,/quiet)

;fitresult=mpfitfun('GBFUNC',rholam_ratio_avg[good]/rho,spec_cor_avg[good],espec_cor_avg[good],p0,parinfo=pi)
;fitresult=mpfitfun('GBFUNC',lam_ratio_sortf,mod_speccorf,mod_speccorf*0,p0,weights=1D/mod_speccorf,parinfo=pi)

;************************


;******The model-predicted covariances
model_speccov=fitresult[0]*exp(-0.5*(rholam_ratio_spline/fitresult[2])^2.)+fitresult[1]*exp(-0.5*(rholam_ratio_spline/(fitresult[3]*rho))^2.)
model_speccov[good_zero]=1-(fitresult[0]+fitresult[1])+fitresult[0]+fitresult[1]
;******

;for i=0L,n_elements(model_speccov)-1 do print,rho*lam_ratio_sortf[i],model_speccov[i]
;for i=0L,n_elements(model_speccov)-1 do print,rholam_ratio_spline[i],model_speccov[i]

print,'A_rho is ',fitresult[0]
print,'A_lambda is ',fitresult[1]
print,'A_delta is ',1-(fitresult[0]+fitresult[1])
print,'sigma_rho is ',fitresult[2]
print,'sigma_lambda is ',fitresult[3]


window,1
setcolors,/system_variables
plot,15*lam_ratio[*,10],spec_cor[*,10],/nodata,yrange=[-0.2,1.1],ystyle=1,xrange=[-3.5,3.5]

oplot,-1*step+findgen(10)*0,findgen(11)*0.1,linestyle=1
oplot,1*step+findgen(10)*0,findgen(11)*0.1,linestyle=1
;print,step
;print,findgen(11)*0.1
;stop

plotsym,0,/fill
;*the interpolated spec_covar values
oploterror,rholam_ratio_spline,spec_cor_avg,espec_cor_avg,psym=3,color=!orange,thick=5,errthick=3

;functional form for fiducial model fit
oplot,rho*lam_ratio(sort(lam_ratio)),mod_speccor(sort(lam_ratio)),color=!cyan,linestyle=2,thick=2
;functional form for model fit solution from MP-FIT
oplot,rholam_ratio_spline[good_inclzero],model_speccov[good_inclzero],color=!magenta,linestyle=0,thick=3
;print,rholam_ratio_spline[good_inclzero]
;stop
;oplot,rholam_ratio_spline,model_speccov,color=!magenta,linestyle=0,thick=3

;***raw covariance data
oplot,rho*lam_ratio,spec_cor,psym=8

plotsym,0,/fill
;****covariances for J band only***
oplot,rho*lam_ratio[0:8,0:8],spec_cor[0:8,0:8],psym=8,color=!blue

;****covariances for H band only***
oplot,rho*lam_ratio[9:14,9:14],spec_cor[9:14,9:14],psym=8,color=!green

;****covariances for K band only***
oplot,rho*lam_ratio[16:20,16:20],spec_cor[16:20,16:20],psym=8,color=!red

set_plot,'ps'
device,filename='spec_covar'+strtrim(rho,2)+'lam_D.eps',/encapsulated,/color,bits=8
!p.font=1
setcolors,/system_variables
plot,15*lam_ratio[*,10],spec_cor[*,10],/nodata,yrange=[-0.2,1.1],ystyle=1,xrange=[-4,4],xstyle=1,$
;plot,15*lam_ratio[*,10],spec_cor[*,10],/nodata,yrange=[-0.2,1.1],ystyle=1,xrange=[-5,5],xstyle=1,$
;xthick=5,ythick=5,xtitle=textoidl('\lambda_{i}')
xtitle=textoidl('\rho (\lambda_{i} - \lambda_{j})/\lambda_{c}'),ytitle=textoidl('\psi_{i,j}'),$
charsize=2,xthick=5,ythick=5,charthick=5

;oplot,-1*step+findgen(10)*0,findgen(11)*0.1,linestyle=1
;oplot,1*step+findgen(10)*0,findgen(11)*0.1,linestyle=1
;print,step
;print,findgen(11)*0.1
;stop

plotsym,0,/fill
;*the interpolated spec_covar values
oploterror,rholam_ratio_spline[good],spec_cor_avg[good],espec_cor_avg[good],psym=3,color=!orange,thick=5,errthick=3

;*the interpolated spec_covar values
;oploterror,rholam_ratio_spline,spec_cor_avg,espec_cor_avg,psym=3,color=!orange,thick=5,errthick=3


;***raw covariance data
oplot,rho*lam_ratio,spec_cor,psym=8,color=!gray,symsize=0.75

plotsym,0,/fill
;****covariances for J band only***
oplot,rho*lam_ratio[0:8,0:8],spec_cor[0:8,0:8],psym=8,color=!blue,symsize=0.75

;****covariances for H band only***
oplot,rho*lam_ratio[9:14,9:14],spec_cor[9:14,9:14],psym=8,color=!green,symsize=0.75

;****covariances for K band only***
oplot,rho*lam_ratio[16:20,16:20],spec_cor[16:20,16:20],psym=8,color=!red,symsize=0.75
;window,1


;plot,rho*lam_ratio_sortf,fitresult,color=!magenta,linestyle=1,thick=3
;functional form for fiducial model fit
;oplot,rho*lam_ratio(sort(lam_ratio)),mod_speccor(sort(lam_ratio)),color=!cyan,linestyle=2,thick=2
;functional form for model fit solution from MP-FIT
;oplot,rholam_ratio_spline[good],model_speccov[good],color=!magenta,linestyle=0,thick=3
oplot,rholam_ratio_spline[good_inclzero],model_speccov[good_inclzero],color=!magenta,linestyle=0,thick=3
;oplot,rholam_ratio_spline,model_speccov,color=!magenta,linestyle=0,thick=3


;***labels
al_legend,textoidl('\rho = ')+strmid(strtrim(rho,2),0,3),/top,box=0,charsize=1.5
al_legend,textoidl('A_{\rho} \sim ')+strmid(strtrim(fitresult[0],1),0,4),position=[2,0.95],box=0,charsize=1.5
al_legend,textoidl('A_{\lambda} \sim ')+strmid(strtrim(fitresult[1],1),0,4),position=[2,0.9],box=0,charsize=1.5
al_legend,textoidl('A_{\delta} \sim ')+strmid(strtrim(1-(fitresult[0]+fitresult[1]),1),0,4),position=[2,0.85],box=0,charsize=1.5
al_legend,textoidl('\sigma_{\rho} \sim ')+strmid(strtrim(fitresult[2],1),0,4),position=[2,0.8],box=0,charsize=1.5
al_legend,textoidl('\sigma_{\lambda} \sim ')+strmid(strtrim(fitresult[3],1),0,4),position=[2,0.75],box=0,charsize=1.5
device,/close

set_plot,'x'

;now write out the spec_covariance into a fits file
writefits,'spec_covar'+strtrim(rho,2)+'lam_D'+'.fits',spec_cor

goto,skipcontour
;contour plot
xrange=indgen(21)
yrange=indgen(21)
set_plot,'ps'
device,filename='contourcovar.eps'
contour,spec_cor
;,xrange,yrange
device,/close
skipcontour:

;if switch thrown, then return fit values and write to file...
if keyword_set(outparams) then begin


print,'A_rho is ',fitresult[0]
print,'A_lambda is ',fitresult[1]
print,'A_delta is ',1-(fitresult[0]+fitresult[1])
print,'sigma_rho is ',fitresult[2]
print,'sigma_lambda is ',fitresult[3]

;a-rho, a-lambda, a-delta, sigma-rho,sigma-lambda,rho
writecol,'covar_params',fitresult[0],fitresult[1],1-(fitresult[0]+fitresult[1]),$
fitresult[2],fitresult[3],rho
outparams=[fitresult[0],fitresult[1],1-(fitresult[0]+fitresult[1]),$
 fitresult[2],fitresult[3]]
endif

skiptoend:
end

 
;*****a second attempt with MPFITFUN
FUNCTION GBFUNC,X,P
;FUNCTION GBFUNC, lam_ratio_sortf,[arhof,alamf,sigrhof,siglamf]
YMOD=P[0]*exp(-0.5*(rho*X/P[2])^2.)+P[1]*exp(-0.5*(X/P[3])^2.)
;YMOD=arhof*exp(-0.5*(rho*lam_ratio_sortf/sigrho)^2.)+alamf*exp(-0.5*(lam_ratio_sortf/siglam)^2.)
return, YMOD
END
