pro charis_aloci_fakeplanet,pfname,reducname=reducname,sdi_reducnamesdi=sdi_reducname,$
nonorthup=nonorthup,$
prefname=prefname,gauspsf=gauspsf,method=method,$
planetmethod=planetmethod,planetmodel=planetmodel,pickpsf=pickpsf,$
reductype=reductype,$
adi=adi,sdi=sdi,adipsdi=adipsdi,sdipadi=sdipadi,$
subsky=subsky
;,filt=filt

;***inserts a planet into a cube, processes the data.

;Version 2.0 - Syntax changes.  Does multiple planets
;*****CAN ONLY DO ONE PLANET AT A TIME****
;***04/21/2018**
;Version 1.1 - Some syntax changes
;***02/08/2018**
;Version 1.0 - forward-modeling ...

;***12/4/2017**
;Version 0.1 -will work assuming that a) you saved the algorithm coefficients before and b) you used local masking
;- outputs a text file, synth_throughput.txt, which contains the estimated throughput vs. spectral channel at a given location.

;**11/22/2017**
;major re-write
;-puts the planet insert into a separate subroutine

;*****************************

;determine reduction subdirectory
;pos1=strpos(pfname,'/')+1
;pos2=strpos(pfname,'.',/reverse_search)

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"charis_aloci_fake_planet: inserts planet into cube, processes the data with A-LOCI"
print,'Written by T. Currie (2017), adapted for CHARIS IFS Data (2019)'
print,''
print,"charis_aloci_fake_planet,pfname,reducname=reducname,sdi_reducname=sdi_reducname,"
print,"method=method,planetmethod=planetmethod,planetmodel=planetmodel,adi=adi,sdi=sdi,adipsdi=adipsdi,sdipadi=sdipadi,"
print,"pickpsf=pickpsf,gauspsf=gauspsf"
print,"reducname=reducname,prefname=prefname"
print,""
print,'Example (ADI):'
print,"charis_aloci_fake_planet,'HR8799_low.info',reducname='aloci.fits',method=0"
print,""
print,'Example (SDI):'
print,"charis_aloci_fake_planet,'HR8799_low.info',reducname='aloci.fits',method=0,reductype=1"
print,""
print,'Example (ADI, p-SDI):'
print,"charis_aloci_fake_planet,'HR8799_low.info',reducname='aloci.fits',sdi_reducname='sdi.fits',reductype=2"
print,""
print,'Example (SDI,p-ADI):'
print,"charis_aloci_fake_planet,'HR8799_low.info',reducname='aloci.fits',sdi_reducname='sdi.fits',reductype=3"
print,""
print,"***Important Keywords****"
print,""
print,"*pfname - parameter file (e.g. HR8799_low.info)"
print,"*reducname - name of ADI-reduced cube"
print,"*reductype - 0 [ADI], 1 [SDI], 2 [ADI, p-SDI], 3 [SDI, p-ADI], 4 (to be added) [ASDI]"
print,"*method - 0 [manually enter x,y, and log(contrast) at prompt], 1 [select a file containing x,y,contrast from GUI]"
print,"*planetmethod - 0 [use .info file name], 1 [use this file in subdir, see below], 2 [pick one from GUI]"
print,"*planetmodel - if selected, will trigger planetmethod = 1"
;print,"*filt - do you spatially filter the cube first?"
print,"*adi - Did you do ADI? (sets reductype = 0, redundant for now)"
print,"*sdi - Did you do SDI? (sets reductype = 1)"
print,"*adipsdi - SDI on ADI residuals (sets reductype=2)"
print,"*sdipadi - ADI on SDI residuals (sets reductype=3)"
print,"*sdi_reducname - name of post-ADI, SDI reduced cube [if /asdi is used]"
print,"*pickpsf, gauspsf - pick PSF from GUI instead of using emp. value; use gaussian PSF"
;print,"*contrast - set the broadband contrast of the synthetic planet"
print,"*subsky - do you subtract the sky background radial profile before comparing b4 and after throughput?"
goto,skiptotheend
endif

reducdir='./reduc/'
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
if keyword_set(sdi) then reductype = 1 ;set to SDI

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

;SDI
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

    else: begin

         modelpath=charis_path(pathname='modeldir')
         inputmodel=dialog_pickfile(Title="Select Planet Model",PATH=modelpath)
         end
endcase

;***Location of Simulated Planet***
;*****Method****
;0. You manually input a log(contrast), x position, y position (Default)
;1. You read a file which has log(contrast), x position, y position
;reading in planet position, brightness ...

if ~keyword_set(method) then begin
;***default is a manually inputing value
     read,'Enter the log(contrast) at H band for the model planet: ',contrast
     read,'Enter the x position for the model planet in the first image: ',xp0
     read,'Enter the y position for the model planet in the first image: ',yp0

endif else begin

;get type name
;methodtype=typename(method)
;if methodtype eq 'INT' then method=strtrim(string(method))

case method of 
;*** manually entering some H-band contrast and x,y position
    'manual': begin
     read,'Enter the log(contrast) at H band for the model planet: ',contrast
     read,'Enter the x position for the model planet in the first image: ',xp0
     read,'Enter the y position for the model planet in the first image: ',yp0
              
              end   

    0: begin
     read,'Enter the log(contrast) at H band for the model planet: ',contrast
     read,'Enter the x position for the model planet in the first image: ',xp0
     read,'Enter the y position for the model planet in the first image: ',yp0
              
       end   
;*** pre-selecting a file at some H-band contrast and position
   'magfile': begin
    magfile=dialog_pickfile(Title="Choose the input file (X,Y,H band Log(Contrast))")
    readcol,magfile,xp0,yp0,contrast
    if n_elements(xp0) eq 1 then begin
    xp0=xp0[0]
    yp0=yp0[0]
    contrast=contrast[0]
    endif
              end
    1: begin
    magfile=dialog_pickfile(Title="Choose the input file (X,Y,H band Log(Contrast))")
    readcol,magfile,xp0,yp0,contrast
    if n_elements(xp0) eq 1 then begin
    xp0=xp0[0]
    yp0=yp0[0]
    contrast=contrast[0]
    endif
       end

    2: begin
    magfile=dialog_pickfile(Title="Choose the input file (X,Y,H band Log(Contrast))")
    readcol,magfile,xp0,yp0,contrast
    if n_elements(xp0) eq 1 then begin
    xp0=xp0[0]
    yp0=yp0[0]
    contrast=contrast[0]
    endif
       end
endcase

endelse

;number of planets
nplanet=n_elements(xp0)

print,'The Planet Position is ',xp0,' ',yp0,' with a Log H-band Contrast of ',contrast
;print,'xp0 stuff',xp0,yp0,contrast
;help,xp0,yp0,contrast
;stop

;***stuff for CHARIS***
;Telescope Diameter for Subaru
Dtel=charis_get_constant(name='Dtel') ;7.9d0 ;visible pupil for SCExAO
;pixel scale 16.4 mas/pixel
pixscale=charis_get_constant(name='pixscale') ;0.0164 nominally

;****the data cubes
    hrefcube=headfits(reducdirorg+reducname,ext=0)
    refcube=readfits(reducdirorg+reducname,ext=1,header)
    reducname_col=(strsplit(reducname,'.',/extract))[0]+'_collapsed.fits'
    adihrefcol=headfits(reducdirorg+reducname_col,ext=0)
    adirefcol=readfits(reducdirorg+reducname_col,ext=1,h1col)
    hrefcol=adihrefcol
    refcol=adirefcol

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


;lat=1.*sxpar(header0,'lat')
;lng=1.*sxpar(header0,'lng')

;reading hour angle and parallactic angle at the beginning of exposures

readcol,'reduc.log',filenum,allxc,allyc,allrsat,allha,allpa


;read in values from reduced data's header
header=headfits(reducdirorg+reducname,ext=1)
testfile=readfits(reducdirorg+reducname,ext=1)
image_dim=(size(testfile))[1]

;dx=fltarr(n_elements(xp))
;dy=fltarr(n_elements(yp))
;rc=fltarr(n_elements(xp))

;***now define PSF size 
;***PSF
dpsf=21
intr_psf=fltarr(dpsf,dpsf,nlambda)
;***
;for now we are going to do just one planet
;nplanet=1
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

planet_input_spectrum=fltarr(nlambda,nplanet,nfiles) ;the array of input planet spectra

;loop calculates the subtraction residuals for fake planets of a given brightness

;define the brightness steps

;flux=fltarr(nplanet)
;mag_flux=fltarr(nplanet)



; xp and yp are given in first image, now calculate new xp and yp positions in subsequent images
dx0=xp0-xc  ;dx in first image
dy0=yp0-yc  ;dy in first image
rc0=sqrt(dx0^2.+dy0^2.)

for n=0,nfiles-1 do begin

        h0sci=headfits(datadir+files[n],ext=0,/silent)
        im=readfits(datadir+files[n],h1sci,/silent,ext=1)
        outputcube=im
;set to empty data cube
        im[*]=1.d-15
;now, add planet to data cube
        outputcube0=im

   for iplanet=0,nplanet-1 do begin

 plan_coords0=cv_coord(/double,from_rect=[dx0[iplanet],dy0[iplanet]],/to_polar,/degrees)
 asc=(plan_coords0[0]-(allpa[n]-allpa[0]))*!dtor  ;array of angles

if ~keyword_set(nonorthup) then asc+=northpa*!pi/180.
xposarr=rc0[iplanet]*cos(asc)+xc  ;array of x positions
yposarr=rc0[iplanet]*sin(asc)+yc  ;array of y positions
;print,nfiles,n
   charis_insert_planet_into_cube,pfname,im,h0sci,h1sci,xposarr,yposarr,contrast[iplanet],intr_psf,inputmodel,cube_out=outputcube_indiv,spec_out=outputspec
   ;charis_insert_planet_into_cube,pfname,im,h0sci,h1sci,xposarr[n],yposarr[n],contrast[iplanet],intr_psf,inputmodel,cube_out=outputcube_indiv,spec_out=outputspec

planet_input_spectrum[*,iplanet,n]=outputspec
   outputcube0+=outputcube_indiv
   endfor
  
        outputcube+=outputcube0
;register image

        writefits,datadir+filesfc[n],0,h0sci
        writefits,datadir+filesfc[n],outputcube,h1sci,/append
;        im[*]=0
 endfor

 ;average input planet spectrum
 planet_avg_input_spectrum=median(planet_input_spectrum,dimension=3,/even)
 plot,lambda,planet_avg_input_spectrum[*,0],linestyle=0,psym=-4,xrange=[1,2.5]
 
 ;planet-to-star contrast
 starspectrum=fltarr(n_elements(lambda))
 ;a simple loop to pull out the star's spectrum
 ;
 for islice = 0L,n_elements(lambda)-1 do begin
 starspectrum[islice]=sxpar(h1sci,'FSTAR_'+strtrim(islice,2))
 endfor
 star_avg_input_spectrum=median(starspectrum,/even)
 
 for iplanet=0,nplanet-1 do begin
 planet_to_star_contrast = median(planet_avg_input_spectrum[*,iplanet],/even)/star_avg_input_spectrum
 print,'The Planet-To-Star Contrast in the Filter is ',planet_to_star_contrast,' for planet ',iplanet+1
 print,'Planet Avg. Flux Is ',median(planet_avg_input_spectrum[*,iplanet],/even),' for planet ',iplanet+1
 endfor
 print,'Star Avg. Flux Is ',star_avg_input_spectrum


    if (zrsub gt 0) then begin
    charis_imrsub,pfname,prefname=prefname,suffname=suffname,rmax=80,/prad,/fc

    suffname='rsub_fc'
    endif

print,'LOCI parameters are ',znfwhm,zdrsub,zna,zgeom,zrmin,zrmax,znref
    ;*** ADI treatment

itc=0

case reductype of

 0: begin   ;just ADI
;A-LOCI with Forward-Modeling

 charis_adialoci,pfname,prefname=prefname,rsubval=zrsub,$
;charis_adialocif,pfname,prefname=prefname,rsubval=zrsub,zonetype=zadiztype,$
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,svd=zsvd,nref=znref,pixmask=zpixmask,zero=zzero,meanadd=zmeanadd,$
/fc,outfile='res'+string(itc,format='(i1)')+'.fits'

    end

 1: begin    ;just SDI

   charis_sdialoci,pfname,prefname=prefname,rsubval=zrsub,$
;charis_sdialocif,pfname,prefname=prefname,rsubval=zrsub,zonetype=zadiztype,$
;suffname='_alocisub',$
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,svd=zsvd,nref=znref,pixmask=zpixmask,zero=zzero,meanadd=zmeanadd,$
/fc,outfile='res'+string(itc,format='(i1)')+'.fits'

    end

 2: begin    ;ADI, SDI on post-ADI residuals


charis_adialoci,pfname,prefname=prefname,rsubval=zrsub,$
;charis_adialocif,pfname,prefname=prefname,rsubval=zrsub,zonetype=zadiztype,$
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,svd=zsvd,nref=znref,pixmask=zpixmask,zero=zzero,meanadd=zmeanadd,$
/fc,/norot

charis_sdialoci,/fwdmod,pfname,prefname=prefname,suffname='_alocisub',/postadi,$
;charis_sdialocif,pfname,prefname=prefname,suffname='_alocisub',/postadi,zonetype=zsdiztype,$
nfwhm=sdiznfwhm,drsub=sdizdrsub,na=sdizna,geom=sdizgeom,rmin=sdizrmin,rmax=sdizrmax,svd=sdizsvd,nref=sdiznref,pixmask=sdizpixmask,zero=sdizzero,meanadd=sdizmeanadd,$
/fc,outfile='res'+string(itc,format='(i1)')+'.fits'

    end


 3: begin   ;SDI, ADI on post-SDI residuals


charis_sdialoci,/fwdmod,pfname,prefname=prefname,/norot,$
;charis_sdialocif,pfname,prefname=prefname,zonetype=zsdiztype,$
nfwhm=sdiznfwhm,drsub=sdizdrsub,na=sdizna,geom=sdizgeom,rmin=sdizrmin,rmax=sdizrmax,svd=sdizsvd,nref=sdiznref,pixmask=sdizpixmask,zero=sdizzero,meanadd=sdizmeanadd,$
/fc,/norot

charis_adialoci,pfname,prefname=prefname,rsubval=zrsub,$
;charis_adialocif,pfname,prefname=prefname,rsubval=zrsub,zonetype=zadiztype,$
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,svd=zsvd,nref=znref,pixmask=zpixmask,zero=zzero,meanadd=zmeanadd,$
suffname='_sdialocisub',/fc,/postsdi,outfile='res'+string(itc,format='(i1)')+'.fits'

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
;nplanet=1
;***now Loop on Wavelength and then Planet to get Attenuation vs. Channel
;   for ilambda=0L,nlambda-1 do begin
;     imslice=imcube[*,*,ilambda]
    for iplanet=0L,nplanet-1 do begin
     for ilambda=0L,nlambda-1 do begin
    
      imslice=imcube[*,*,ilambda]
 
      ;compute the flux at the original position
      ;aper,imslice,xp,yp,flux,eflux,sky,skyerr,1,0.5*fwhm[ilambda],[2,6]*fwhm[ilambda],setskyval=0,/flux,/exact,/nan,/silent

     if keyword_set(subsky) then begin
      profrad_tc,imslice,1,1,rmax,p2d=pr
      imslice-=pr
     endif

      aper,imslice,xp0[iplanet],yp0[iplanet],flux,eflux,sky,skyerr,1,aperrad[ilambda],[2,6]*2*aperrad[ilambda],setskyval=0,/flux,/exact,/nan,/silent

       print,'position and flux is ',xp0[iplanet],yp0[iplanet],flux
      ;now compare to the original, input flux to get the attenuation.
      res_flux[ilambda,iplanet]=flux/planet_avg_input_spectrum[ilambda,iplanet]
      print,'Throughput for Planet ',string(iplanet+1),' at Wavelength ',string(lambda[ilambda]),$
           ' with contrast of ',planet_avg_input_spectrum[ilambda,iplanet]/starspectrum[ilambda],' is ...',res_flux[ilambda,iplanet]
      ;print,'Throughput for Planet ',string(iplanet+1),' at Wavelength ',string(lambda[ilambda]),' is ...',res_flux[ilambda,iplanet]
     endfor
    endfor
 
     ;writecol,reducdir+'synth_throughput.txt',lambda,res_flux
     ;writecol,'synth_throughput.txt',lambda,res_flux

    for iplanet=0L,nplanet-1 do begin
    writecol,'synth_throughput'+strtrim(iplanet,2)+'.txt',lambda,res_flux[*,iplanet]
    endfor
  
  for iplanet=0L,nplanet-1 do begin
  planet_to_star_contrast = median(planet_avg_input_spectrum[*,iplanet],/even)/star_avg_input_spectrum
   print,'Planet Avg. Flux Is ',median(planet_avg_input_spectrum[*,iplanet],/even),' for planet ',iplanet+1
   print,'The Planet-to-Star Contrast in Broadband Is ',planet_to_star_contrast, ' for planet ',iplanet+1
  print,'The Planet-to-Star Contrast is JHK Is roughly ',median(planet_avg_input_spectrum[0:5,iplanet]/starspectrum[0:5],/even),median(planet_avg_input_spectrum[8:13,iplanet]/starspectrum[8:13],/even),$
     median(planet_avg_input_spectrum[15:20,iplanet]/starspectrum[15:20],/even)
  endfor
 print,'Star Avg. Flux Is ',star_avg_input_spectrum


;*******TO DO ****** 11/25****
;***Throughput for wavelength collapsed image
 ;   for iplanet=0L,nplanet-1 do begin
 ;    aper,imcol,xp[iplanet],yp[iplanet],flux,eflux,sky,skyerr,1,0.5*median(fwhm,/even),[2,4]*median(fwhm,/even),setskyval=0,/flux,/exact,/nan,/silent 
 ;   endfor
 

;this is the raw difference.
;    writefits,'diff.fits',refcol-imcol
;    writefits,'new.fits',imcol
;    ;stop
;    help,im2,im


skiptotheend:
end
