pro charis_aloci_negspec_planet,pfname,reducname=reducname,sdi_reducname=sdi_reducname,reductype=reductype,prefname=prefname,gauspsf=gauspsf,method=method,$
adi=adi,sdi=sdi,adipsdi=adipsdi,sdipadi=sdipadi,$
nonorthup=nonorthup, planetmethod=planetmethod,planetmodel=planetmodel,pickpsf=pickpsf,$
magfile=magfile

;Version 1.1 - added documentation
;*****CAN ONLY DO ADI RIGHT NOW****
;*****CAN ONLY DO ONE PLANET AT A TIME****
;****06/18/2018**
;Version 1.0 - negative planet iteration to hone forward-modeling solution
;****

;*****************************

;setupdir,reducdir=reducdir

;determine reduction subdirectory
;pos1=strpos(pfname,'/')+1
;pos2=strpos(pfname,'.',/reverse_search)

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"charis_aloci_negspec_planet: A-LOCI negative planet subtraction, Version 1.1 (June 2020)"
print,'Written by T. Currie (2014), adapted for CHARIS IFS Data (2019)'
print,''
print,"charis_aloci_negspec_planet,pfname,reducname=reducname,sdi_reducname=sdi_reducname,"
print,"method=method,planetmethod=planetmethod,planetmodel=planetmodel,"
print,"rdi=rdi,refpath=refpath"
print,"ntc=ntc,dlrstep=dlrstep,lrmin=lrmin,lrmax=lrmax"
print,"pickpsf=pickpsf,gauspsf=gauspsf,contrast=contrast"
print,"reducname=reducname,prefname=prefname"
print,""
print,'Example (ADI):'
print,"charis_aloci_negspec_planet,'HR8799_low.info',reducname='aloci.fits',method=0"
print,""
print,"***Important Keywords****"
print,""
print,"*pfname - parameter file (e.g. HR8799_low.info)"
print,"*reducname - name of ADI-reduced cube"
print,"*method - 0 [manually enter x,y, and log(contrast) at prompt], 1 [select a file containing x,y,contrast from GUI]"
print,"*planetmethod - 0 [use .info file name], 1 [use this file in subdir, see below], 2 [pick one from GUI]"
print,"*planetmodel - if selected, will trigger planetmethod = 1"
;print,"*filt - do you spatially filter the cube first?"
print,"*adi - Did you do ADI? (sets reductype = 0, redundant for now)"
print,"*sdi - Did you do SDI? (sets reductype = 1)"
print,"*adipsdi - SDI on ADI residuals (sets reductype=2)"
print,"*sdipadi - ADI on SDI residuals (sets reductype=3)"
;print,"*sdi_reducname - name of post-ADI, SDI reduced cube [if /asdi is used]"
print,"*pickpsf, gauspsf - pick PSF from GUI instead of using emp. value; use gaussian PSF"
print,"*magfile- file with x,y, contrast"
goto,skiptotheend
endif

reducdir='./reduc/'
;subdir='procsub/'
subdir='proc/'
reducdir1=reducdir+'reg/'
datadir=reducdir1
reducdir+=subdir

reducdirorg=reducdir
;*****

if ~keyword_set(prefname) then begin
;***Prefixes and such***
prefname='n'
endif

suffname='reg_cal'

;*****Reduction Type **** 0=ADI, 1=SDI, 2=ADIpSDI, 3=SDIpADI, 4=ASDI
if ~keyword_set(reductype) then reductype = 0  ;set ADI as default
if keyword_set(adi) then reductype = 0 ;set to ADI
if keyword_set(sdi) then reductype = 1 ;set to SDI
;if keyword_set(asdi) then reductype = 4 ;set to ADIpSDI

;****to implement later
if keyword_set(adipsdi) then reductype = 2 ;set to ADIpSDI
if keyword_set(sdipadi) then reductype = 3 ;set to SDIpADI
;if keyword_set(asdi) then reductype = 3 ;set to ASDI

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


;***Selecting a Planet Model***
;-0/default: you select via .info file
;-1: you manually input the file name from the planetspec directory
;-2: you select it via gui

if ~keyword_set(planetmethod) then planetmethod='0'
if keyword_set(planetmodel) then planetmethod='1'

case planetmethod of
   '0':begin
     param,'planmod',inputmodel0,/get,pfname=pfname
       modelpath=charis_path(pathname='planetdir')
       inputmodel=modelpath+inputmodel0
       end

   '1':begin
       modelpath=charis_path(pathname='planetdir')
       inputmodel=modelpath+planetmodel
       end

   '2':begin
        inputmodel=dialog_pickfile(Title="Select Planet Spectrum")
       end

    else: begin

         modelpath=charis_path(pathname='modeldir')
         inputmodel=dialog_pickfile(Title="Select Planet Model",PATH=modelpath)
         end
endcase

;***Location of Simulated Planet***
;*****Method****

case method of
    0: begin
     ;read,'Enter the log(contrast) at H band for the model planet: ',contrast
     contrast=-4
     read,'Enter the x position for the model planet in the first image: ',xp0
     read,'Enter the y position for the model planet in the first image: ',yp0

       end

    1: begin


    if ~keyword_set(magfile) then begin
    magfile=dialog_pickfile(Title="Choose the input file (X,Y,H band Log(Contrast))")
    endif
readcol,magfile,xp0,yp0,contrast

       end
endcase

print,'The Planet Position is ',xp0,' ',yp0,' with a Log H-band Contrast of ',contrast

;***stuff for CHARIS***
;Telescope Diameter for Subaru
Dtel=charis_get_constant(name='Dtel') ;7.9d0 ;visible pupil for SCExAO
pixscale=charis_get_constant(name='pixscale') ;0.0164 nominally

;****the data cubes
    header0=headfits(reducdirorg+reducname,ext=0)
    refcube=readfits(reducdirorg+reducname,ext=1,header)
    reducname_col=(strsplit(reducname,'.',/extract))[0]+'_collapsed.fits'
    hrefcol=headfits(reducdirorg+reducname_col,ext=0)
    refcol=readfits(reducdirorg+reducname_col,ext=1,h1col)
    ;*** ADI treatment

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

 print,'LOCI parameters are ',znfwhm,zdrsub,zna,zgeom,zrmin,zrmax


imx=sxpar(header,'naxis1')
dimx=sxpar(header,'naxis2')
dimy=dimx   ;assume square arrays
xc=dimx/2 & yc=dimy/2

;Now Get the Wavelength Vector
get_charis_wvlh,header,lambda
lambda*=1d-3
nlambda=n_elements(lambda)
;determine array of FWHM
fwhm=1.0*(1.d-6*lambda/Dtel)*(180.*3600./!dpi)/pixscale

;****extraction radius*****
aperrad=fltarr(n_elements(lambda))
for ir=0L,n_elements(lambda)-1 do begin
aperrad[ir]=sxpar(header,'R_AP'+strtrim(string(ir),2))
endfor

;;***Now get basic image properties

;the below is actually not needed since the charis_northup function now modifies ROTANG to compensate for the north angle offset
;angoffset=charis_get_constant(name='angoffset')
case reductype of
   0: begin ;ADI
northpa=sxpar(header,'ROTANG')
      end
   1: begin  ;SDI
northpa=sxpar(header,'ROTANG')
      end
   2: begin ;ADI, p-SDI
northpa=sxpar(sdiheader,'ROTANG')
      end
   3: begin  ;SDI, p-ADI
northpa=sxpar(header,'ROTANG')
      end
   4: begin ;ASDI
northpa=sxpar(header,'ROTANG')
      end
endcase
;************


;define the temporary directory
tmpdir=reducdir+'tmp/'

;parameters of the sequence
param,'obsdate',date,/get,pfname=pfname & date=strtrim(date,2)
param,'fnum_sat',flist,/get,pfname=pfname

filenum=nbrlist(flist)
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
filesfc=filelist(filenum,nfiles,prefix=prefname,suffix=suffname+'_fc')


lat=1.*sxpar(header0,'lat')
lng=1.*sxpar(header0,'lng')

;reading hour angle and parallactic angle at the beginning of exposures

readcol,'reduc.log',filenum,allxc,allyc,allrsat,allha,allpa


;read in values from reduced data's header
;header=headfits(reducdirorg+reducname,ext=1)
;testfile=readfits(reducdirorg+reducname,ext=1)
;image_dim=(size(testfile))[1]

;***now define PSF size 
;***PSF
dpsf=21
intr_psf=fltarr(dpsf,dpsf,nlambda)
;***
;for now we are going to do just one planet
nplanet=1
print,'loop size is ',nplanet

;****Define Your Prior: The Intrinsic PSF****
;****

;1.
;use an empirical PSF

if ~keyword_set(gauspsf) then begin

if ~keyword_set(pickpsf) then begin
;-automatically
fname=findfile('./psfmodel/psfcube_med.fits',count=c)

if c eq 0 then begin
;-manually
fname=dialog_pickfile(Title="Pick Reference PSF Files",/MULTIPLE_FILES)
c=n_elements(fname)
endif

endif else begin
fname=dialog_pickfile(Title="Pick Reference PSF Files",/MULTIPLE_FILES)
c=n_elements(fname)

endelse

;loop in wavelength
for il=0,nlambda-1 do begin
;Since you have coronographic data and/or saturated data

; best to just use a model PSF, e.g. c=0
;if c gt 0 then begin
    fpeak=0.
    for n=0,c-1 do begin
        imcube=readfits(fname[n],ext=1)
        im=imcube[*,*,il]
        im=subarr(im,dpsf)
;make the PSF a unit PSF: signal of 1 within a FWHM-sized aperture
        im/=charis_myaper(im,dpsf/2,dpsf/2,0.5*fwhm[il])
        if n eq 0 then psf=im/c else psf+=im/c
    endfor
   intr_psf[*,*,il]=psf
endfor

;2.
;use a simple gaussian
endif else begin

for il=0L,nlambda -1 do begin

;approximate the PSF by a gaussian; should be okay for high Strehl and simple PSF structure (e.g. Subaru or VLT, not Keck)
a=psf_gaussian(npixel=dpsf,fwhm=fwhm[il],/double,/normalize)
a/=charis_myaper(a,dpsf/2,dpsf/2,0.5*fwhm[il])

intr_psf[*,*,il]=a
endfor
endelse

;11/25 - empirical PSF done


;now should have [xp0,yp0] in first image, contrast in H band, planet model, intrinsic PSF

;x,y coordinates of pixels of the PSF centered on pixel 0.0

xpsf=lindgen(dpsf)#replicate(1l,dpsf)-dpsf/2
ypsf=replicate(1l,dpsf)#lindgen(dpsf)-dpsf/2

;indices of pixels of the psf centered on pixel 0.0 in the big picture
ipsf=xpsf+ypsf*dimx

imt=fltarr(dpsf,dpsf,nplanet)
nsc=21 ;nombre de sous-compagnon


;*****Keywords for RA,DEC, individual exposure time (exp1time), and coadds (so exp1time*ncoadds = exptime)*****
    param,'RA',radeg,/get,pfname=pfname
    param,'DEC',decdeg,/get,pfname=pfname


res_flux=fltarr(nlambda,nplanet)  ;residual flux (1-attenuation)
;diff_r=fltarr(nlambda,nplanet)  ;astrometric bias in angular separation
;diff_ang=fltarr(nlambda,nplanet) ;astrometric bias in position angle

planet_input_spectrum=fltarr(nlambda,nfiles) ;the array of input planet spectra

;loop calculates the subtraction residuals for fake planets of a given brightness

;define the brightness steps


; xp and yp are given in first image, now calculate new xp and yp positions in subsequent images
dx0=xp0[0]-xc  ;dx in first image
dy0=yp0[0]-yc  ;dy in first image
rc0=sqrt(dx0^2.+dy0^2.)

;plan_coords0=cv_coord(/double,from_rect=[dx0,dy0],/to_polar,/degrees)
;asc=(plan_coords0[0]-(allpa-allpa[0]))*!dtor  ;array of angles

;if ~keyword_set(nonorthup) then asc+=northpa*!pi/180.

;xposarr=rc0*cos(asc)+xc  ;array of x positions
;yposarr=rc0*sin(asc)+yc  ;array of y positions

for n=0,nfiles-1 do begin

        h0sci=headfits(datadir+files[n],ext=0,/silent)
        im=readfits(datadir+files[n],h1sci,/silent,ext=1)

  plan_coords0=cv_coord(/double,from_rect=[dx0,dy0],/to_polar,/degrees)
asc=(plan_coords0[0]-(allpa[n]-allpa[0]))*!dtor  ;array of angles      
if ~keyword_set(nonorthup) then asc+=northpa*!pi/180.
xposarr=rc0*cos(asc)+xc  ;array of x positions
yposarr=rc0*sin(asc)+yc  ;array of y positions
    
;now, add planet to data cube

   charis_insert_planet_into_cube,pfname,im,h0sci,h1sci,xposarr,yposarr,contrast,intr_psf,inputmodel,cube_out=outputcube,spec_out=outputspec,$
    empspectrum=inputmodel,/negplanet

;writefits,'outputcube.fits',outputcube
planet_input_spectrum[*,n]=outputspec
;register image

        writefits,datadir+filesfc[n],0,h0sci
        writefits,datadir+filesfc[n],outputcube,h1sci,/append
        im[*]=0
 endfor

;stop
;average input planet spectrum
planet_avg_input_spectrum=median(planet_input_spectrum,dimension=2,/even)
plot,lambda,-1*planet_avg_input_spectrum,linestyle=0,psym=-4,xrange=[1,2.5]


    ;*** ADI treatment


    if (zrsub gt 0) then begin
    charis_imrsub,pfname,prefname=prefname,suffname=suffname,rmax=80,/prad,/fc
    suffname='rsub_fc'
    endif

 print,'LOCI parameters are ',znfwhm,zdrsub,zna,zgeom,zrmin,zrmax



itc=0

case reductype of

 0: begin   ;just ADI
;A-LOCI with Forward-Modeling

charis_adialoci,pfname,prefname=prefname,rsubval=zrsub,$
;charis_adialocif,pfname,prefname=prefname,rsubval=zrsub,zonetype=zadiztype,$
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,svd=zsvd,nref=znref,pixmask=zpixmask,zero=zzero,meanadd=zmeanadd,$
/fc,outfile='res'+string(itc,format='(i1)')+'_np_'+'.fits'

    end

 1: begin    ;just SDI

charis_sdialoci,pfname,prefname=prefname,rsubval=zrsub,$
;charis_sdialocif,pfname,prefname=prefname,rsubval=zrsub,zonetype=zadiztype,$
;suffname='_alocisub',$
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,svd=zsvd,nref=znref,pixmask=zpixmask,zero=zzero,meanadd=zmeanadd,$
/fc,outfile='res'+string(itc,format='(i1)')+'_np_'+'.fits'

    end

 2: begin    ;ADI, SDI on post-ADI residuals


charis_adialoci,pfname,prefname=prefname,rsubval=zrsub,$
;charis_adialocif,pfname,prefname=prefname,rsubval=zrsub,zonetype=zadiztype,$
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,svd=zsvd,nref=znref,pixmask=zpixmask,zero=zzero,meanadd=zmeanadd,$
/fc,/norot

charis_sdialoci,/fwdmod,pfname,prefname=prefname,suffname='_alocisub',/postadi,$
;charis_sdialocif,pfname,prefname=prefname,suffname='_alocisub',/postadi,zonetype=zsdiztype,$
nfwhm=sdiznfwhm,drsub=sdizdrsub,na=sdizna,geom=sdizgeom,rmin=sdizrmin,rmax=sdizrmax,svd=sdizsvd,nref=sdiznref,pixmask=sdizpixmask,zero=sdizzero,meanadd=sdizmeanadd,$
/fc,outfile='res'+string(itc,format='(i1)')+'_np_'+'.fits'

    end


 3: begin   ;SDI, ADI on post-SDI residuals


charis_sdialoci,/fwdmod,pfname,prefname=prefname,/norot,$
;charis_sdialocif,pfname,prefname=prefname,zonetype=zsdiztype,$
nfwhm=sdiznfwhm,drsub=sdizdrsub,na=sdizna,geom=sdizgeom,rmin=sdizrmin,rmax=sdizrmax,svd=sdizsvd,nref=sdiznref,pixmask=sdizpixmask,zero=sdizzero,meanadd=sdizmeanadd,$
/fc,/norot

charis_adialoci,pfname,prefname=prefname,rsubval=zrsub,$
;charis_adialocif,pfname,prefname=prefname,rsubval=zrsub,zonetype=zadiztype,$
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,svd=zsvd,nref=znref,pixmask=zpixmask,zero=zzero,meanadd=zmeanadd,$
suffname='_sdialocisub',/fc,/postsdi,outfile='res'+string(itc,format='(i1)')+'_np_'+'.fits'

    end


 else: begin

    print,'cannot find reduction method'
    print,'FWD-Mod Failure'
    goto,skiptotheend

       end

endcase


;Now read back in the files, rescale them to minimize chi-squared, add fits header information, and save them with unique input file name        
   
;your original file 
 ;  remember ...
    ;hrefcube=headfits(reducdir+reducname,ext=0)
    ;refcube=readfits(reducdir+reducname,ext=1,h1)
    ;reducname_col=(strsplit(reducname,'.',/extract))[0]+'_collapsed.fits'
    ;hrefcol=headfits(reducdir+reducname_col,ext=0)
    ;refcol=readfits(reducdir+reducname_col,ext=1,h1col)


;the empty cube, after processing
    h0cube=headfits(reducdir+'res'+string(itc,format='(i1)')+'.fits')
    imcube=readfits(reducdir+'res'+string(itc,format='(i1)')+'.fits',h1cube,ext=1)

    modelname='res'+string(itc,format='(i1)')+'.fits'
    modelbasename=(strsplit(modelname,'.',/extract,count=modcount))[0]

;wavelength-collapsed version.
    h0col=headfits(reducdir+'res'+string(itc,format='(i1)')+'_collapsed.fits')
    imcol=readfits(reducdir+'res'+string(itc,format='(i1)')+'_collapsed.fits',h1col,ext=1)

;for now do just one planet
nplanet=1
;***now Loop on Wavelength and then Planet to get Attenuation vs. Channel 
   for ilambda=0L,nlambda-1 do begin 
     imslice=imcube[*,*,ilambda]
    for iplanet=0L,nplanet-1 do begin 
    
     ;compute the flux at the original position 
     ;aper,imslice,xp,yp,flux,eflux,sky,skyerr,1,0.5*fwhm[ilambda],[2,6]*fwhm[ilambda],setskyval=0,/flux,/exact,/nan,/silent
     aper,imslice,xp0,yp0,flux,eflux,sky,skyerr,1,aperrad[ilambda],[2,6]*2*aperrad[ilambda],setskyval=0,/flux,/exact,/nan,/silent

     ;now compare to the original, input flux to get the attenuation.
     ;res_flux[ilambda]=flux/flux_planet[ilambda]
     res_flux[ilambda,iplanet]=flux/planet_avg_input_spectrum[ilambda]
     print,'Throughput for Planet ',string(iplanet+1),' at Wavelength ',string(lambda[ilambda]),' is ...',res_flux[ilambda,iplanet]
    endfor
   endfor

    ;writecol,reducdir+'synth_throughput.txt',lambda,res_flux
    ;writecol,'synth_throughput.txt',lambda,res_flux

;*******TO DO ****** 11/25****
;***Throughput for wavelength collapsed image
 ;   for iplanet=0L,nplanet-1 do begin
 ;    aper,imcol,xp[iplanet],yp[iplanet],flux,eflux,sky,skyerr,1,0.5*median(fwhm,/even),[2,4]*median(fwhm,/even),setskyval=0,/flux,/exact,/nan,/silent 
 ;   endfor
 
skiptotheend:

end
