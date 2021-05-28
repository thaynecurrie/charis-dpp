pro charis_recalc_disk_fwdmod_chisq,reducname=reducname,prefname=prefname,pickmodels=pickmodels,$
roirange=roirange,lrmin=lrmin,lrmax=lrmax,snoise=snoise,help=help

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"charis_recalc_disk_fwdmod_chisq,reducname=reducname,prefname=prefname,pickmodels=pickmodels"
print,"roirange=roirange,lrmin=lrmin,lrmax=lrmax,snoise=snoise"
print,"recalculates the chi-square of a forward-modeled disk image using a different region of interest"
print,""
print,"***Keywords*"
print,""
print,"*reducname - name of reduced cube (the actual data)"
print,"*pickmodels - select the forward-modeled disks for which you want to recalculate chisq"
print,"*roirang - x,y, and rotation angle of ROI for chisq calculation"
print,"*lrmin -inner radius for ROI"
print,"*lrmax -outer radius for ROI"
print,"*snoise - scaling of the noise"

goto,skiptotheend
endif

;recalculates the chi-square of a forward-modeled disk image using a different region of interest

;******
reducdir='./reduc/'
subdir='proc/'
reducdir1=reducdir+'reg/'
datadir=reducdir1
reducdir+=subdir

reducdirorg=reducdir

reducdir_model1=reducdir+'/model_in/'
reducdir_model2=reducdir+'/model_psfsub/'
reducdir_modelstat=reducdir+'/model_stat/'
reducdir_modelatten=reducdir+'/model_atten/'


openw,31,reducdir_modelstat+'NEW_fit_outcomes'+timestamp()+'.dat'
printf,31,'MODEL_NAME','G','KSI0','ALP_I','ALP_O','BETA','XDO','YDO','E','THETA0','CHISQ/dof',$
   format='(a20,a5,1x,8(a6,1x),a9)'


;****


;***stuff for CHARIS***
;Telescope Diameter for Subaru
Dtel=charis_get_constant(name='Dtel') ;7.9d0 ;visible pupil for SCExAO
;pixel scale 16.4 mas/pixel
pixscale=charis_get_constant(name='pixscale') ;0.0164 nominally


if ~keyword_set(reducname) then begin
reducname='final.fits'
endif

;****the data cubes ********
    hrefcube=headfits(reducdirorg+reducname,ext=0)
    refcube=readfits(reducdirorg+reducname,ext=1,h1)
    reducname_col=(strsplit(reducname,'.',/extract))[0]+'_collapsed.fits'
    hrefcol=headfits(reducdirorg+reducname_col,ext=0)
    refcol=readfits(reducdirorg+reducname_col,ext=1,h1col)
;stop
imx=sxpar(h1,'naxis1')
dimx=sxpar(h1,'naxis2')
dimy=dimx   ;assume square arrays
xc=dimx/2 & yc=dimy/2

;array of radial distances
dist_circle,rarray,[dimx,dimy],[xc,yc]

;Now Get the Wavelength Vector
get_charis_wvlh,hrefcube,lambda
lambda*=1d-3
nlambda=n_elements(lambda)
;determine array of FWHM
fwhm=1.0*(1.d-6*lambda/Dtel)*(180.*3600./!dpi)/pixscale
medfwhm=median(fwhm,/even)

;**********

;***SNR Map for Real Data***
;-self-consistently calculate the SNR map of the collapsed image
;- save SNR map, the aperture-summed image, and the uncertainty map
snratio_sub,refcol,fwhm=medfwhm,rmax=lrmax,/zero,/finite,snrmap=snrmap,noisemap=sigdata,imcol=refcol_con

if ~keyword_set(snoise) then snoise =1.
sigdata=sigdata/snoise

;*************************************
;generate the NEW region of interest

charis_generate_roi,output=roi_geom,roidim=[dimx,dimy],roirange=roirange
roitotal=roi_geom
roi_rad=where(rarray ge lrmin and rarray le lrmax)
roitotal=where(roitotal eq 1)
roitotal=intersect(roi_rad,roitotal)
roiregion=fltarr(dimx,dimy)
roiregion[roitotal]=1
writefits,'roiregion.fits',roiregion


;new region of interest defined
;*************************************

 if keyword_set(pickmodels) then begin
 modelname=dialog_pickfile(Title="Select Your Forward-Modeled Disk")
 endif else begin
 modelname=file_search(reducdir_model2+'*_psfsubcol.fits')
 endelse

;loop on the models

    for ij=0L,n_elements(modelname)-1 do begin

    imcol=readfits(modelname[ij],h_psf,ext=1,h1)
    writefits,'ii.fits',imcol

    ;now read in the important parameters
model_g=sxpar(h_psf,'G')
model_ksi=sxpar(h_psf,'KSI0')
model_alphain=sxpar(h_psf,'ALPHAIN')
model_alphaout=sxpar(h_psf,'ALPHAOUT')
model_beta=sxpar(h_psf,'BETA')
model_xdo=sxpar(h_psf,'XDO')
model_ydo=sxpar(h_psf,'YDO')
model_ecc=sxpar(h_psf,'E')
model_theta=sxpar(h_psf,'THETA0')
model_itilt=sxpar(h_psf,'ITILT')
model_r0=sxpar(h_psf,'r0')

;   define imcol as the forward-modeled disk image.

    ;define earlier, up top ;refcol_con=sumaper_im(refcol,raper,lrmax,/nan)
    imcol_con=sumaper_im(imcol,medfwhm*0.5,lrmax,/nan)


    diffcube=refcol_con-imcol_con
    diffcubeu=refcol-imcol

    chisqmap=fltarr(dimx,dimy)
    chisqmap[roitotal]=abs(diffcube[roitotal]/sigdata[roitotal])^2.
    charis_calc_chisqbin,chisqmap,medfwhm,chisq_bin_out,nbins,chisqmap_bin=chisqmap_bin
    chisq_col=chisq_bin_out
    print,chisq_col,nbins

    modelbasename=(strsplit(modelname[ij],'/',/extract,count=modcount))[modcount-1]
    modelbasename0=(strsplit(modelbasename,'.',/extract,count=modcount))[0]
   

    printf,31,modelbasename,model_g,model_ksi,model_alphain,model_alphaout,model_beta,$
              model_xdo,model_ydo,model_ecc,model_theta,chisq_col,format='(a20,8(f6.3,1x),f8.2,1x,f6.3)'


    endfor

skiptotheend:
end
