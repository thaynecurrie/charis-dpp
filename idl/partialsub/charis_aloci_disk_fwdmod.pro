pro charis_aloci_disk_fwdmod,pfname,reducname=reducname,suffname=suffname,$
;northup=northup,$
noconv=noconv,$
simplescale=simplescale,$
lrmin=lrmin,lrmax=lrmax,snrlim=snrlim,roirange=roirange,roimanual=roimanual,$
maskrad=maskrad,$
psf=psf,gausspsf=gausspsf,absflux=absflux,manscale=manscale,calcattenfact=calcattenfact,$
synthscale=synthscale,$
asdi=asdi,sdi_reducnamesdi=sdi_reducname,$
ffac=ffac,ntc=ntc,test=test


;Version 2.3 - updated documentation
;***2/7/2019**
;Version 2.2 - added absflux switch for physically motivated models (e.g. MCMax3D; not GRaTeR)
;               in this case, you just feed in a model with flux density as-is

;***8/11/2018**
;Version 2.1 - clean up documentation of code ...
;****2/5/2018**
;Version 2.0 - CHARIS/A-LOCI Forward-Modeling for Disks (complete).
;               code is cleaned up, new structure using PSF models (rewrite of charis_makeemppsf)
;               new output of attenuated model
;               new switch (calcattenfact) to only output the model cube and attenuation factor if you throw it
;****1/29/2018**
;Version 1.1 - CHARIS/A-LOCI Forward-Modeling for Disks

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"charis_aloci_disk_fwdmod: A-LOCI fwd-modeling of a disk source, Version 2.3 (June 2020)"
print,'Written by T. Currie (2014), adapted for CHARIS IFS Data (2019)'
print,''
print,"charis_aloci_disk_fwdmod,pfname,reducname=reducname,"
;print,"method=method,planetmethod=planetmethod,planetmodel=planetmodel,filt=filt,"
print,"rdi=rdi,refpath=refpath"
;print,"ntc=ntc,dlrstep=dlrstep,lrmin=lrmin,lrmax=lrmax"
print,"pickpsf=pickpsf,gauspsf=gauspsf,contrast=contrast"
print,"reducname=reducname,prefname=prefname"
print,""
print,'Example'
print,"charis_aloci_disk_fwdmod,'LkCa15_low.info',/rdi,refpath='../../V819Tau/charis_nearirsep0917',reducname='refaloci_final9cal.fits',method=1,/northup,/absflux,manscale=1"
print,""
print,"***Important Keywords****"
print,""
print,"*pfname - parameter file (e.g. LkCa15_low.info)"
print,"*reducname - name of ADI-reduced cube"
;print,"*filt - do you spatially filter the cube first?"
print,'*asdi - did you do ADI+SDI?'
print,'*sdi_reducnamesdi - name of SDI-reduced cube'
print,"*pickpsf, gauspsf - pick PSF from GUI instead of using emp. value; use gaussian PSF"
print,"*roirange - x,y, and rotation angle of ROI for chisq calculation"
print,"*roimanual - set the above one step at a time"
print,"*snrlim - select only regions above this S/N for evaluation"
print,"*lrmin,lrmax -inner,lrmax radius for ROI"
print,"*simplescale - do not do a noise-weighted minimization of disk model"
print,"*manscale - don't try to initially scale the disk model (before fwdmod) to minimize residuals"
print,"*synthscale - scale disk model by this amount initially"
print,"*absflux - select for physically motivated models (like MCMax3D)"
print,"*calcattenfact - calculate an attenuation cube to do throughput correction of real disk"
print,"*noconv - do not convolve disk model with instrumental PSF"
print,"savecube - save the output cube of the forward-model, not just the collapsed image"
goto,skiptotheend
endif

;if ~keyword_set(snrlim) then snrlim = 2.

if ~keyword_set(prefname) then begin
;***Prefixes and such***
prefname='n'
endif

if ~keyword_set(suffname) then suffname='reg_cal'
suffname2='regcal'

    
;*****************************

; calculates the attenuation from A-LOCI
; ntc = number of theta companion, default = 15
; / test to stop at the end of the loop to see
;And the residual flux for testing ntc

;setupdir,reducdir=reducdir

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
file_mkdir,reducdir_model1
file_mkdir,reducdir_model2
file_mkdir,reducdir_modelstat
if keyword_set(calcattenfact) then file_mkdir,reducdir_modelatten
openw,31,reducdir_modelstat+'fit_outcomes'+timestamp()+'.dat'
printf,31,'MODEL_NAME','G','KSI0','ALP_I','ALP_O','BETA','XDO','YDO','E','THETA0','CHISQ/dof',$
   format='(a20,a5,1x,8(a6,1x),a9)'
;******


;*****Reduction Type **** 0=ADI, 1=SDI, 2=ADIpSDI, 3=SDIpADI, 4=ASDI
if ~keyword_set(reductype) then reductype = 0  ;set ADI as default
if keyword_set(sdi) then reductype = 1 ;set to SDI
if keyword_set(adipsdi) then reductype = 2 ;set to ADIpSDI
if keyword_set(sdipadi) then reductype = 3 ;set to SDIpADI

;*******Reduced Files to Use*******
if ~keyword_set(reducname) then begin

;file name, with full path
reducname_full_path=dialog_pickfile(Title="Select Your Processed Image")

your_path=strpos(reducname_full_path,'/',/reverse_search)+1

;determining the long-form of the reduction directory
reducdirorg=strmid((reducname_full_path),0,your_path)

;now determining the base name of your reduced file
z=strsplit(reducname_full_path,'/',/extract)
reducname=z[n_elements(z)-1]
endif



if (reductype eq 2 or reductype eq 3) and ~keyword_set(sdi_reducname) then begin
;SDI

;file name, with full path
sdi_reducname_full_path=dialog_pickfile(Title="Select Your Processed Image (Assuming SDI Subtraction)")

sdi_your_path=strpos(sdi_reducname_full_path,'/',/reverse_search)+1

;determining the long-form of the reduction directory
sdi_reducdirorg=strmid((sdi_reducname_full_path),0,sdi_your_path)

;now determining the base name of your reduced file
z=strsplit(sdi_reducname_full_path,'/',/extract)
sdi_reducname=z[n_elements(z)-1]
endif


;***stuff for CHARIS***
;Telescope Diameter for Subaru
Dtel=charis_get_constant(name='Dtel') ;7.9d0 ;visible pupil for SCExAO
;pixel scale 16.4 mas/pixel
pixscale=charis_get_constant(name='pixscale') ;0.0164 nominally

;North PA offset
angoffset=charis_get_constant(name='angoffset')  ;nominally -2.2 east of north


if ~keyword_set(reducname) then begin
reducname='final.fits'
endif

;****the data cubes
    hrefcube=headfits(reducdirorg+reducname,ext=0)
    refcube=readfits(reducdirorg+reducname,ext=1,h1)
    reducname_col=(strsplit(reducname,'.',/extract))[0]+'_collapsed.fits'
    hrefcol=headfits(reducdirorg+reducname_col,ext=0)
    refcol=readfits(reducdirorg+reducname_col,ext=1,h1col)


if (keyword_set(sdi_reducname) or keyword_set(asdi)) then begin
    sdihrefcube=headfits(reducdirorg+sdi_reducname,ext=0)
    sdirefcube=readfits(reducdirorg+sdi_reducname,ext=1,sdi_h1)
    sdi_reducname_col=(strsplit(sdi_reducname,'.',/extract))[0]+'_collapsed.fits'
    sdihrefcol=headfits(reducdirorg+sdi_reducname_col,ext=0)
    sdirefcol=readfits(reducdirorg+sdi_reducname_col,ext=1,sdi_h1col)
endif


imx=sxpar(h1,'naxis1')
dimx=sxpar(h1,'naxis2')
dimy=dimx   ;assume square arrays
xc=dimx/2 & yc=dimy/2

;Now Get the Wavelength Vector
get_charis_wvlh,hrefcube,lambda
lambda*=1d-3
nlambda=n_elements(lambda)
;determine array of FWHM
fwhm=1.0*(1.d-6*lambda/Dtel)*(180.*3600./!dpi)/pixscale
medfwhm=median(fwhm,/even)

;nominally, we use medfwhm for the aperture in snrmap calculations and evaluating chi-sq.   
;option to change medfwhm to 2*aperture radius if need be but will require additional code changes where applicable.
;param,'raper',raper,/get,pfname=pfname


;since you're just a disk model at a single location.
nrc=1

;FOR DISK keep centered on central pixel!
rc=0.

;number of angles to position


header=headfits(reducdir+reducname,ext=1)

;***********
;ALOCI parameters from fits header in output file

;for ADI or SDI-only
znfwhm=sxpar(header,'loci_nfw')
zna=sxpar(header,'loci_na')
zgeom=sxpar(header,'loci_geo')
zdrsub=sxpar(header,'loci_drs')
zrmin=sxpar(header,'loci_rmi')
zrmax=sxpar(header,'loci_rma')

;zcorrlim=sxpar(header,'corr_lim',count=countcor)
zsvd=sxpar(header,'svd',count=countcutoff)
znref=sxpar(header,'loci_nre',count=countnref)
zpixmask=sxpar(header,'pixmaskv')
zzero=sxpar(header,'zero')
zmeanadd=sxpar(header,'meanadd')
zrsub=sxpar(header,'rsub')
zadiztype=sxpar(header,'adiztype')

if countcutoff eq 0 then zcutoff=1.d-99

;for ADI + SDI on residuals [2], SDI + ADI on residuals [3]
if (reductype eq 2 or reductype eq 3) then begin
    sdiheader0=headfits(reducdirorg+sdi_reducname,ext=0)
    sdirefcube=readfits(reducdirorg+sdi_reducname,ext=1,sdiheader)
    sdi_reducname_col=(strsplit(sdi_reducname,'.',/extract))[0]+'_collapsed.fits'
    sdihrefcol=headfits(reducdirorg+sdi_reducname_col,ext=0)
    sdirefcol=readfits(reducdirorg+sdi_reducname_col,ext=1,sdi_h1col)
    hrefcol=sdihrefcol
    refcol=sdirefcol
sdiznfwhm=sxpar(sdiheader,'loci_nfw')
sdizna=sxpar(sdiheader,'loci_na')
sdizgeom=sxpar(sdiheader,'loci_geo')
sdizdrsub=sxpar(sdiheader,'loci_drs')
sdizrmin=sxpar(sdiheader,'loci_rmi')
sdizrmax=sxpar(sdiheader,'loci_rma')
zsdiztype=sxpar(sdiheader,'sdiztype')

;sdizcorrlim=sxpar(sdiheader,'corr_lim',count=sdicountcor)
sdizsvd=sxpar(sdiheader,'svd',count=sdicountcutoff)
sdiznref=sxpar(sdiheader,'loci_nre',count=sdicountnref)
sdizpixmask=sxpar(sdiheader,'pixmaskv')
sdizzero=sxpar(sdiheader,'zero')
sdizmeanadd=sxpar(sdiheader,'meanadd')
if sdicountcutoff eq 0 then sdizcutoff=1.d-99
endif

if (reductype eq 3) then zrsub=sxpar(sdiheader,'rsub')


;define the temporary directory
tmpdir=reducdir+'tmp/'


;parameters of the sequence
param,'obsdate',date,/get,pfname=pfname & date=strtrim(date,2)
param,'fnum_sat',flist,/get,pfname=pfname

filenum=nbrlist(flist)
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
filesfc=filelist(filenum,nfiles,prefix=prefname,suffix=suffname+'_fc')


;header of the first image to determine dimension

h=headfits(datadir+files[0],ext=1)
dimx=sxpar(h,'naxis1')
dimy=sxpar(h,'naxis2')
xc=dimx/2 & yc=dimy/2
;print,'xc yc',xc,yc

intr_psf=fltarr(dimx,dimy,nlambda)

;****if you have no empirical PSF and are ok with just a gaussian ... (not optimal)
if keyword_set(gausspsf) then begin
for il=0L,nlambda -1 do begin
;print,il,il,lambda[il]

;determine FWHM
fwhm_channel=(1.d-6*lambda[il]/Dtel)*(180.*3600./!dpi)/pixscale

;approximate the PSF by a gaussian; should be okay for high Strehl
a=psf_gaussian(npixel=201,fwhm=fwhm_channel,centroid=[xc,yc],/normalize,/double)
intr_psf[*,*,il]=a
endfor
goto,breakoutpsf
endif 

;****if you want to use an empirical PSF (a better idea than a gaussian)
if keyword_set(psf) then begin

psffile=dialog_pickfile(Title="Choose Empirical PSF")
intr_psf=readfits(psffile,ext=1)
goto,breakoutpsf
endif

psfdir='./psfmodel/'
print,'Using model PSF in "./psfmodel/" directory: psfcube_med.fits'
intr_psf=readfits(psfdir+'psfcube_med.fits',ext=1)

breakoutpsf:

lat=1.*sxpar(h,'lat')
lng=1.*sxpar(h,'lng')

;reading hour angle and parallactic angle at the beginning of exposures

readcol,'reduc.log',filenum,allxc,allyc,allrsat,allha,allpa

;****radii for which to calculate the self subtraction
;****Note this should be EQUAL to what you actually do for PSF-subtracting the real data
if ~keyword_set(lrmax) then lrmax=dimx/2
if ~keyword_set(lrmin) then lrmin=5

;***SNR Map***
;-self-consistently calculate the SNR map of the collapsed image
;- save SNR map, the aperture-summed image, and the uncertainty map
if ~keyword_set(maskrad) then maskrad=0.1
charis_snratio_sub,refcol,fwhm=medfwhm,rmax=zrmax,/zero,/finite,snrmap=snrmap,noisemap=sigdata,imcol=refcol_con

;if you need to mask a region for disk emission
if keyword_set(maskrad) then begin
dist_circle,radmask,dimx
masknoise=where(radmask le maskrad,nmasknoise,complement=nomaskrad)
correct_noise=where(radmask gt maskrad and sigdata ne 0 and finite(sigdata) ne 0,ncorrectnoise)
maxcorrect_noise=max(sigdata[correct_noise])
sigdata[masknoise]=maxcorrect_noise
snrmap=refcol_con/sigdata
endif

;uncomment if you want this
;writefits,'refcol.fits',refcol_con
;writefits,'sigdata.fits',sigdata
;writefits,'snr.fits',snrmap

;*****defining region of interest for evaluating the disk model...
;-0. manually (based on geometry)
;-1. pre-defined region of interest
;-2. snrlimit
;-3. combination of 1 and 2

;array of radial distances
dist_circle,rarray,[dimx,dimy],[xc,yc]

if keyword_set(roimanual) or (~keyword_set(snrlim) and ~keyword_set(roirange)) then begin

read,'Set the Width of the Evaluation Region (set to 201 if no restriction)',gah1
read,'Set the Height of the Evaluation Region (set to 201 if no restriction)',gah2
read,'Set the Position Angle of the Evaluation Region (set to 0 if no restriction)',gah3
roirange=[gah1,gah2,gah3]
charis_generate_roi,output=roi_geom,roidim=[dimx,dimy],roirange=roirange
roitotal=roi_geom
goto,breakoutroi
endif

if keyword_set(snrlim) and ~keyword_set(roirange) then begin
charis_generate_roi,output=roi_snrmap,roidim=[dimx,dimy],snrmap=snrmap
roitotal=roi_snrmap
goto,breakoutroi
endif 

if keyword_set(roirange) and ~keyword_set(snrlim) then begin
charis_generate_roi,output=roi_geom,roidim=[dimx,dimy],roirange=roirange
roitotal=roi_geom
goto,breakoutroi
endif

if keyword_set(roirange) and keyword_set(snrlim) then begin
charis_generate_roi,output=roi_snrmap,roidim=[dimx,dimy],snrmap=snrmap
charis_generate_roi,output=roi_geom,roidim=[dimx,dimy],roirange=roirange
roitotal=long(roi_snrmap*roi_geom)
goto,breakoutroi
endif

breakoutroi:
roi_rad=where(rarray ge lrmin and rarray le lrmax)
roitotal=where(roitotal eq 1)
roitotal=intersect(roi_rad,roitotal)
roiregion=fltarr(dimx,dimy)
roiregion[roitotal]=1
;writefits,'roifit.fits',roiregion

;roitotal[roitotal_ind]=1
;writefits,'roitotal.fits',roitotal

;***PSF
;****modify for disk!
dpsf=dimx
;h=headfits(datadir+files[0])

;if you just want a mock disk, you can throw these two switches...
;if ~keyword_set(xext_fwhm) then xext_fwhm=120.
;if ~keyword_set(yext_fwhm) then yext_fwhm=10.

;Here is where you can put in a disk model from GRaTeR ...
;find all model files
modelname=dialog_pickfile(Title="Select Your Model Disks",/multiple)
ntc=n_elements(modelname)
print,ntc
print,modelname

;Now, loop on the model grid!

for itc=0,ntc-1 do begin ;loop sur les angles

;the input disk model
psf=readfits(modelname[itc],h_psf)

gg=sxpar(h_psf,'NAXIS')

if gg eq 0 then begin

psf=readfits(modelname[itc],h_psf,/ext)
gg=sxpar(h_psf,'NAXIS')
if gg eq 0 then begin
print,'ERROR!'
return
endif

endif

dimpsf=size(psf,/dim)

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

;scaled to match image (a guess)
if ~keyword_set(synthscale) then begin
synthcube_scale_guess=1.*total(psf[roitotal]*refcol[roitotal],/nan)/total(psf[roitotal]*psf[roitotal],/nan)
print,'scale is ',synthcube_scale_guess
endif else begin
synthcube_scale_guess=synthscale
endelse

;synthcube_scale_guess=1d-6
if ~keyword_set(absflux) then psf*=synthcube_scale_guess
psf+=1.d-10
psf_input=psf

;treat the PSF as a cube.
psf_chan=fltarr(dimx,dimy,nlambda)
psf_chanconv=fltarr(dimx,dimy,nlambda)

for il=0L,nlambda-1 do begin

if ~keyword_set(noconv) then begin
;convolve the input disk model with the PSF
if (n_elements(dimpsf) eq 2) then begin
psf_chan[*,*,il]=convolve(psf,intr_psf[*,*,il])
endif else begin
psf_chan[*,*,il]=convolve(psf[*,*,il],intr_psf[*,*,il])
endelse

endif else begin


if (n_elements(dimpsf) eq 2) then begin
psf_chan[*,*,il]=psf
endif else begin
psf_chan[*,*,il]=psf[*,*,il]
endelse

endelse

psf_chanconv[*,*,il]=psf_chan[*,*,il]

endfor

;the input model (convolved with PSF)
writefits,'psf_chan.fits',psf_chanconv
;... wavelength-collapsed
writefits,'psf.fits',median(psf_chanconv,dimension=3,/even)


;x,y coordinates of pixels of the PSF centered on pixel 0.0

xpsf=lindgen(dpsf)#replicate(1l,dpsf)-dpsf/2
ypsf=replicate(1l,dpsf)#lindgen(dpsf)-dpsf/2

;indices of pixels of the psf centered on pixel 0.0 in the big picture
ipsf=xpsf+ypsf*dimx

;stretched to the psf radius rc=

imt=fltarr(dpsf,dpsf,nrc)
nsc=21 ;number of sub-images for each image slice

;***the below through the smear_fac loop might be vestigal and might be delete-able
for n=0,nfiles-1 do begin

;*****Keywords for RA,DEC, individual exposure time (exp1time), and coadds (so exp1time*ncoadds = exptime)*****
    h=headfits(datadir+files[n],ext=0)
    radeg=sxpar(h,'ra')
    decdeg=sxpar(h,'dec')

    exptime=float(sxpar(h,'exp1time'))
    ncoadds=float(sxpar(h,'coadds'))
endfor

;*************************************************************
;***Change back to original version, which was the right way to calculate this!
;***smearing effect is negligible.  eliminated for a cleaner code.

smear_fac=fltarr(nrc)
smear_fac[*]=1.
imt=0
;-**


;loop calculates the angles and the self-subtraction
print,'ntc is ',ntc

print,'NEXT!',itc
;add in fake disks here
    ;constructed images false companions
    print,'itc='+strtrim(itc+1,2)+'/'+strtrim(ntc,2)

    print,' adding model disks ...'
;  "adding companions.."
    for n=0,nfiles-1 do begin

        ;im=readfits(reducdir+files[n],h,/silent)
        h0=headfits(datadir+files[n],ext=0,/silent)
        im=readfits(datadir+files[n],h,/silent,ext=1)
        iminit=im
;now get the north PA offset
        getrot,h,northpa,cdelt


;****** NOTE: In this convention, angoffset is SUBTRACTED.   CONSIDER CHANGING THIS
        northpa-=angoffset

        ;print,northpa,allpa[n]

        nwvlh=(size(im,/dim))[2]
;empty cube
        im[*]=0
        ;im[*]=1.d-15
        ;parallactic angle of sub companions

        exptime=float(sxpar(h0,'exp1time'))
        ncoadds=float(sxpar(h0,'coadds'))
        ha=allha[n]
        x=ha+(findgen(nsc)/(nsc-1.))*(exptime*ncoadds/3600.) ;SI HA DEBUT POSE
        pa0=parangle(x,decdeg,lat)
        dpa=pa0-pa0[0]
        ;dpa[*]=0
        pa=allpa[n]+dpa

;loop over wavelength channels
  for il=0L,nwvlh-1 do begin

            outsidefov=where(iminit[*,*,il] eq 0,noutsidefov)
        for irc=0,nrc-1 do begin
            ;determine angles des sous-compagnons

            ;asc=(replicate(tcomp[itc],nsc)-(pa-allpa[0]))*!dtor
            ;asc=(replicate(tcomp[itc]+(irc mod 10)*36.+(irc mod 2)*(180.-36.),nsc)+(pa-allpa[0]))*!dtor

            ;asc=(replicate(tcomp[itc]+(irc mod 10)*36.+(irc mod 2)*(180.-36.),nsc)-(pa-allpa[0]))*!dtor

            ;determine x,y positions of sub-companions

;            xsc=rc[irc]*cos(asc)+xc
;            ysc=rc[irc]*sin(asc)+yc

            ;xsc=rc[irc]*cos(-1*asc)+xc
            ;ysc=rc[irc]*sin(-1*asc)+yc
           
            for isc=0,nsc-1 do begin
                ;xe=floor(xsc[isc]) & ye=floor(ysc[isc])
                ;dxf=xsc[isc]-xe    & dyf=ysc[isc]-ye

; so basically don't shift the model disk at all.  very small effect.
                dxf=0 & dyf=0
                psfs=shift_sub(psf_chan[*,*,il],dxf,dyf)
              

;If you have an IFS data format written by normal human beings then the first line is right.
                im[*,*,il]+=(rotat(psfs,(-1*northpa))/smear_fac[irc]/nsc)
                scratch=im[*,*,il]
                bad=where(finite(scratch) eq 0)
                scratch[bad]=0
                ;if(noutsidefov gt 0) then scratch[outsidefov]=0
                im[*,*,il]=scratch
            endfor
        endfor
  endfor
        ;register image
        writefits,datadir+filesfc[n],0,h0
        writefits,datadir+filesfc[n],im,h,/append
        im[*]=0
    endfor

    ;*** ADI treatment
   
    if (zrsub gt 0) then begin
    charis_imrsub,pfname,prefname=prefname,suffname=suffname,rmax=80,/prad,/fc
    suffname='rsub_fc'
    endif


 print,'LOCI parameters are ',znfwhm,zdrsub,zna,zgeom,zrmin,zrmax


case reductype of

 0: begin   ;just ADI
;A-LOCI with Forward-Modeling

charis_adialoci,pfname,prefname=prefname,rsubval=zrsub,$
;charis_adialocif,pfname,prefname=prefname,rsubval=zrsub,zonetype=zadiztype,$
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,svd=zsvd,nref=znref,pixmask=zpixmask,zero=zzero,meanadd=zmeanadd,$
/fc,/fwdmod

    end

 1: begin    ;just SDI

charis_sdialoci,pfname,prefname=prefname,rsubval=zrsub,$
;charis_sdialocif,pfname,prefname=prefname,rsubval=zrsub,zonetype=zadiztype,$
;suffname='_alocisub',$
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,svd=zsvd,nref=znref,pixmask=zpixmask,zero=zzero,meanadd=zmeanadd,$
/fc,/fwdmod

    end

 2: begin    ;ADI, SDI on post-ADI residuals


charis_adialoci,pfname,prefname=prefname,rsubval=zrsub,$
;charis_adialocif,pfname,prefname=prefname,rsubval=zrsub,zonetype=zadiztype,$
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,svd=zsvd,nref=znref,pixmask=zpixmask,zero=zzero,meanadd=zmeanadd,$
/fc,/fwdmod,/norot

charis_sdialoci,/fwdmod,pfname,prefname=prefname,suffname='_alocisub',/postadi,$
;charis_sdialocif,pfname,prefname=prefname,suffname='_alocisub',/postadi,zonetype=zsdiztype,$
nfwhm=sdiznfwhm,drsub=sdizdrsub,na=sdizna,geom=sdizgeom,rmin=sdizrmin,rmax=sdizrmax,svd=sdizsvd,nref=sdiznref,pixmask=sdizpixmask,zero=sdizzero,meanadd=sdizmeanadd,$
/fc,/fwdmod

    end


 3: begin   ;SDI, ADI on post-SDI residuals


charis_sdialoci,/fwdmod,pfname,prefname=prefname,/norot,$
;charis_sdialocif,pfname,prefname=prefname,zonetype=zsdiztype,$
nfwhm=sdiznfwhm,drsub=sdizdrsub,na=sdizna,geom=sdizgeom,rmin=sdizrmin,rmax=sdizrmax,svd=sdizsvd,nref=sdiznref,pixmask=sdizpixmask,zero=sdizzero,meanadd=sdizmeanadd,$
/fc,/norot,/fwdmod

charis_adialoci,pfname,prefname=prefname,rsubval=zrsub,$
;charis_adialocif,pfname,prefname=prefname,rsubval=zrsub,zonetype=zadiztype,$
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,svd=zsvd,nref=znref,pixmask=zpixmask,zero=zzero,meanadd=zmeanadd,$
suffname='_sdialocisub',/fc,/postsdi,/fwdmod

    end


 else: begin

    print,'cannot find reduction method'
    print,'FWD-Mod Failure'
    goto,skiptotheend

       end

endcase


;Now read back in the files, rescale them to minimize chi-squared, add fits header information, and save them with unique input file name
    h0cube=headfits(reducdir+'final_fc.fits') 
    imcube=readfits(reducdir+'final_fc.fits',ext=1,h1cube)
    h0col=headfits(reducdir+'final_fc_collapsed.fits')
    imcol=readfits(reducdir+'final_fc_collapsed.fits',ext=1,h1col)

;model base name, extract string here and in line below to determine how to write out scaled input model and scaled PSF sub model
    modelbasename=(strsplit(modelname[itc],'/',/extract,count=modcount))[modcount-1]
    modelbasename0=(strsplit(modelbasename,'.',/extract,count=modcount))[0]

;and read back in the original collapsed PSF and PSF/vs. channel
    psf_input=readfits('psf.fits')
    psf_input_spec=readfits('psf_chan.fits')

   ;the collapsed cube 
;****
    ;define earlier, up top ;refcol_con=sumaper_im(refcol,raper,lrmax,/nan)
    imcol_con=sumaper_im(imcol,medfwhm*0.5,lrmax,/nan)

if keyword_set(simplescale) then begin
    synthcol_scale=total(refcol_con[roitotal]*imcol_con[roitotal],/nan)/total(imcol_con[roitotal]*imcol_con[roitotal],/nan)
endif else begin
    synthcol_scale=total((refcol_con*imcol_con/sigdata^2.)[roitotal],/nan)/total((imcol_con*imcol_con/sigdata^2.)[roitotal],/nan)
endelse

    ;synthcol_scale=total(refcol_con[roitotal]*imcol_con[roitotal],/nan)/total(imcol_con[roitotal]*imcol_con[roitotal],/nan)

;if you want to skip scaling in case it messes up
    if keyword_set(manscale) then synthcol_scale=1
    imcol_con*=synthcol_scale
    imcol*=synthcol_scale
    psf_input*=synthcol_scale


    diffcube=refcol_con-imcol_con
    diffcubeu=refcol-imcol

;now define the convolved real data and model flux for comparison
    ;chisq_col=total(abs(diffcube[roitotal]/sigdata[roitotal])^2.,/nan)/(n_elements(roitotal)-1)
    ;chisqmap=total(abs(diffcube[roitotal]/sigdata[roitotal])^2.,/nan)
    chisqmap=fltarr(dimx,dimy)
    chisqmap[roitotal]=abs(diffcube[roitotal]/sigdata[roitotal])^2.

    writefits,'imcol_con.fits',imcol_con
    writefits,'refcol_con.fits',refcol_con
    writefits,'diffcube.fits',diffcube
    writefits,'diffcubeu.fits',diffcubeu
    writefits,'sigdata.fits',sigdata
    writefits,'chisqmap.fits',chisqmap

;Now, bin the image to resolution, take the total chisquared value of the binned image, find number of bins, and 
; then compute the reduced chi-squared in the binned image.
     
    charis_calc_chisqbin,chisqmap,medfwhm,chisq_bin_out,nbins,chisqmap_bin=chisqmap_bin
    chisq_col=chisq_bin_out
    nres=nbins
    writefits,'chisqmapbin.fits',chisqmap_bin

    sxaddpar,h1col,'sfac',synthcol_scale
    sxaddpar,h1col,'csq_col (binned image)',chisq_col
    sxaddpar,h1col,'npts (binned image)',nres
    sxaddpar,h1col,'G',model_g
    sxaddpar,h1col,'ksi0',model_ksi
    sxaddpar,h1col,'alphain',model_alphain
    sxaddpar,h1col,'alphaout',model_alphaout
    sxaddpar,h1col,'beta',model_beta
    sxaddpar,h1col,'r0',model_r0
    sxaddpar,h1col,'xdo',model_xdo
    sxaddpar,h1col,'ydo',model_ydo
    sxaddpar,h1col,'e',model_ecc
    sxaddpar,h1col,'theta0',model_theta
    sxaddpar,h1col,'itilt',model_itilt

;PSF-subbed model
    writefits,reducdir_model2+modelbasename0+'_psfsubcol.fits',0,h0col
    writefits,reducdir_model2+modelbasename0+'_psfsubcol.fits',imcol,h1col,/append
;input model (with scaling applied)
    writefits,reducdir_model1+modelbasename0+'_inputcolsc.fits',0,h0col
    writefits,reducdir_model1+modelbasename0+'_inputcolsc.fits',psf_input,h1col,/append
;***  

if keyword_set(calcattenfact) then begin 
   ;the data cube
    synthcube_scale=fltarr((size(imcube,/dim))[2])
    for icube=0L,(size(imcube,/dim))[2]-1 do begin
    imcube_slice=imcube[*,*,icube]
    refcube_slice=refcube[*,*,icube]
    refcube_con=sumaper_im(refcube_slice,medfwhm/2.,lrmax,/nan)
    imcube_con=sumaper_im(imcube_slice,medfwhm/2.,lrmax,/nan)
   
;here, we aren't defining chisq in a robust way.  Leave that to the collapsed cube for now. 
    synthcube_scale[icube]=total(refcube_con[roitotal]*imcube_con[roitotal],/nan)/total(imcube_con[roitotal]*imcube_con[roitotal],/nan)
    ;synthcube_scale[icube]=total(refcube_slice[roitotal]*imcube_slice[roitotal],/nan)/total(imcube_slice[roitotal]*imcube_slice[roitotal],/nan)
    if keyword_set(manscale) then synthcube_scale[icube]=1
    imcube[*,*,icube]*=synthcube_scale[icube]
    psf_input_spec[*,*,icube]*=synthcube_scale[icube]
    diffslice=refcube_slice-imcube_slice
    chisqslice=total(abs(diffslice[roitotal])^2.,/nan)/(n_elements(roitotal)-1)


    sxaddpar,h1cube,'sfac_'+strtrim(string(icube),2),synthcube_scale[icube],after='sfac_'+strtrim(string(icube-1),2)
    sxaddpar,h1cube,'csq_'+strtrim(string(icube),2),chisqslice,after='csq_'+strtrim(string(icube-1),2)
    endfor
    
    sxaddpar,h1cube,'sfac',synthcol_scale
    sxaddpar,h1col,'csq_col (binned image)',chisq_col
    sxaddpar,h1col,'npts (binned image)',nres
    sxaddpar,h1cube,'G',model_g
    sxaddpar,h1cube,'ksi0',model_ksi
    sxaddpar,h1cube,'alphain',model_alphain
    sxaddpar,h1cube,'alphaout',model_alphaout
    sxaddpar,h1cube,'beta',model_beta
    sxaddpar,h1cube,'r0',model_r0
    sxaddpar,h1cube,'xdo',model_xdo
    sxaddpar,h1cube,'ydo',model_ydo
    sxaddpar,h1cube,'e',model_ecc
    sxaddpar,h1cube,'theta0',model_theta
    sxaddpar,h1cube,'itilt',model_itilt

    writefits,reducdir_model2+modelbasename0+'_psfsubcube.fits',0,h0cube
    writefits,reducdir_model2+modelbasename0+'_psfsubcube.fits',imcube,h1cube,/append
    writefits,reducdir_model1+modelbasename0+'_inputcubesc.fits',0,h0cube
    writefits,reducdir_model1+modelbasename0+'_inputcubesc.fits',psf_input_spec,h1cube,/append
    relatten=(psf_input_spec-imcube)/psf_input_spec
    notdisk=where(abs(relatten) gt 1.5,nnotdisk)
    if nnotdisk gt 0 then relatten[notdisk]=!values.f_nan
    writefits,reducdir_modelatten+modelbasename0+'_relattencube.fits',0,h0cube
    writefits,reducdir_modelatten+modelbasename0+'_relattencube.fits',relatten,h1cube,/append
endif
    

    print,modelbasename,' ','r-chisq is ',chisq_col

;just a diagnostic
    writefits,'diff.fits',0,hrefcol
    writefits,'diff.fits',diffcubeu,h1col,/append

;now that you have everything done, save the model parameters in a file with their chisq values
    printf,31,modelbasename,model_g,model_ksi,model_alphain,model_alphaout,model_beta,$
              model_xdo,model_ydo,model_ecc,model_theta,chisq_col,format='(a20,8(f6.3,1x),f8.0,1x,f6.3)'


endfor

skiptotheend:

end  ;end program
