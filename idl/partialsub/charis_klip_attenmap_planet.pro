pro charis_klip_attenmap_planet,pfname,reducname=reducname,prefname=prefname,gauspsf=gauspsf,$
planetmethod=planetmethod,planetmodel=planetmodel,pickpsf=pickpsf,$
contrast=contrast,$
reductype=reductype,$
rdi=rdi,$
ntc=ntc,dlrstep=dlrstep,lrmin=lrmin,lrmax=lrmax,help=help

;Version 2.0 ** -- attenuation map calculation for KLIP.   Does ADI and RDI. ;DOESN'T have reductype syntax set yet
;****08/21/2020 -- documentation and help function included
;Version 1.0 
;****11/12/2019** - attenuation map calculation for KLIP.   Does ADI and RDI.   Does not yet do SDI.
;Version 0.1
;***09/05/2018** - attenuation map calculation.   Forward-models planet signal all over the place; interpolates results to get a full map.
;****for now does array sampling same as Lafreniere; spiral pattern 
;****
;*****************************

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"charis_klip_attenmap_planet: PCA/KLIP attenuation map, Version 2.0 (Aug 2020)"
print,'Written by T. Currie (2014), adapted for CHARIS IFS Data (2019)'
print,'Note: RDI,SDI not yet implemented!!!!!!'
print,''
print,"charis_klip_attenmap_planet,pfname,planetmethod=planetmethod,planetmodel=planetmodel,filt=filt"
print,"ntc=ntc,dlrstep=dlrstep,lrmin=lrmin,lrmax=lrmax"
;print,"filt=filt,rdi=rdi,asdi=asdi,"
print,"pickpsf=pickpsf,gauspsf=gauspsf,contrast=contrast,"
print,"reducname=reducname,prefname=prefname"
print,""
print,'Example (ADI):'
print,"charis_klip_attenmap_planet,'LkCa15_low.info',reducname='klip.fits',ntc=11"
print,""
print,"***Important Keywords****"
print,""
print,"*pfname - parameter file (e.g. LkCa15_low.info)"
print,"*reducname - name of ADI-reduced cube"
;print,"*reductype - 0 [ADI], 1 [SDI (to be added)], 2 [ADIpSDI (to be added), 3 [ASDI (to be added)], 4 RDI"
print,"*planetmethod - 0 [use .info file name], 1 [use this file in subdir, see below], 2 [pick one from GUI]"
print,"*planetmodel - if selected, will trigger planetmethod = 1"
print,"*filt - do you spatially filter the cube first?"
print,"*pickpsf, gauspsf - pick PSF from GUI instead of using emp. value; use gaussian PSF"
print,"*contrast - set the broadband contrast of the synthetic planet [shouldn't matter]"
print,"*ntc - number of angle ranges for which fwd-mod is repeated for attenmap calc"
print,"*dlrstep - radial step of synthetic planets"
print,"*lrmin, lrmax - minimum and maximum angular sep of synthetic planets"
goto,skiptotheend
endif


;determine reduction subdirectory
;pos1=strpos(pfname,'/')+1
;pos2=strpos(pfname,'.',/reverse_search)

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

;*****Reduction Type **** 0=ADI, 1=SDI, 2=ADIpSDI, 3=ASDI
if ~keyword_set(reductype) then reductype = 0  ;set ADI as default
;if keyword_set(sdi) then reductype=1 ; set to SDI
;if keyword_set(asdi) then reductype = 2 ;set to ADIpSDI
if keyword_set(rdi) then reductype = 4 ;set to RDI

;****to implement later
;if keyword_set(adipsdi) then reductype = 2 ;set to ADIpSDI
;if keyword_set(asdi) then reductype = 3 ;set to ASDI

;*******Reduced Files to Use*******
;ADI
if ~keyword_set(reducname) then begin

;file name, with full path
reducname_full_path=dialog_pickfile(Title="Select Your Processed Image (Assuming Just ADI Subtraction)")

your_path=strpos(reducname_full_path,'/',/reverse_search)+1

;determining the long-form of the reduction directory
reducdirorg=strmid((reducname_full_path),0,your_path)

;now determining the base name of your reduced file
z=strsplit(reducname_full_path,'/',/extract)
reducname=z[n_elements(z)-1]
endif




;*******Fake Planet Range******
; set the range in radius for fake planets
if ~keyword_set(dlrstep) then dlrstep = 2.5
;for now, hardwire the inner radius to be 10 pixels and outer radius to be 60 pixels
if ~keyword_set(lrmin) then lrmin = 10
if ~keyword_set(lrmax) then lrmax = 60
nrc=round((lrmax-lrmin)/dlrstep)+1
rc=findgen(nrc)*dlrstep+lrmin  
;rc = array of radial separations
;for now, we are going to do the same thing that Lafreniere does.
;

;ntc = number of iterations
if ~keyword_set(ntc) then ntc=5
;***********

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
;***Hardwire Contrast for now ...
;contrast = -6
if ~keyword_set(contrast) then contrast = -6

;***stuff for CHARIS***
;Telescope Diameter for Subaru
Dtel=charis_get_constant(name='Dtel') ;7.9d0 ;visible pupil for SCExAO
;pixel scale 16.4 mas/pixel
pixscale=charis_get_constant(name='pixscale') ;0.0164 nominally

;****the data cubes
    header0=headfits(reducdirorg+reducname,ext=0)
    refcube=readfits(reducdirorg+reducname,ext=1,header)
    reducname_col=(strsplit(reducname,'.',/extract))[0]+'_collapsed.fits'
    hrefcol=headfits(reducdirorg+reducname_col,ext=0)
    refcol=readfits(reducdirorg+reducname_col,ext=1,h1col)
    ;*** ADI treatment

;***********
;parameters from fits header in output file

znfwhm=sxpar(header,'klip_nfw')
zdrsub=sxpar(header,'klip_drs')
zrmin=sxpar(header,'klip_rmi')
zrmax=sxpar(header,'klip_rma')
znpca=sxpar(header,'klip_pca')
znzone=sxpar(header,'klip_nzo')
zzero=sxpar(header,'zero')
zmeansub=sxpar(header,'meansub')
zmeanadd=sxpar(header,'meanadd')
zrsub=sxpar(header,'rsub')

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
northpa=sxpar(header,'ROTANG')

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


;now should have the contrast in H band, planet model, intrinsic PSF

;x,y coordinates of pixels of the PSF centered on pixel 0.0

xpsf=lindgen(dpsf)#replicate(1l,dpsf)-dpsf/2
ypsf=replicate(1l,dpsf)#lindgen(dpsf)-dpsf/2

;indices of pixels of the psf centered on pixel 0.0 in the big picture
;ipsf=xpsf+ypsf*dimx

;imt=fltarr(dpsf,dpsf,nplanet)
nsc=21 ;nombre de sous-compagnon


;*****Keywords for RA,DEC, individual exposure time (exp1time), and coadds (so exp1time*ncoadds = exptime)*****
    param,'RA',radeg,/get,pfname=pfname
    param,'DEC',decdeg,/get,pfname=pfname

;*****throughput maps
throughput_map=fltarr(dimx,dimy) ;map of the throughput of the collapsed-cube
throughput_map_lambda=fltarr(dimx,dimy,nlambda) ;map of the throughput per channel
;initialize the maps to have NaN values, populate maps after outcome of calculations
throughput_map[*]=!values.f_nan
throughput_map_lambda[*]=!values.f_nan

;*****variables saving the throughput in each loop
res_flux=fltarr(nlambda,nrc,ntc) ;residual flux per wavelength per separation
res_flux_collapse=fltarr(nrc,ntc) ;residual collapsed-wavelength flux per separation

planet_input_spectrum=fltarr(nlambda,nfiles) ;the array of input planet spectra

;loop calculates the subtraction residuals for fake planets of a given brightness

;define the brightness steps

;flux=fltarr(nplanet)
;mag_flux=fltarr(nplanet)



tcomp=0.0
if ~keyword_set(ntc) then ntc=1

xplanet_position=fltarr(nrc,ntc)
yplanet_position=fltarr(nrc,ntc)
rplanet_position=fltarr(nrc,ntc)

;defining the planet position at image,0
gridfakes=fltarr(201,201)
generate_grids,xgridpos,ygridpos,201
for irc=0,nrc-1 do begin
for itc=0,ntc-1 do begin
tcmp=(long(itc)/float(ntc))*360
;ac=(tcmp+(irc mod 10)*36.+(irc mod 2)*(180.-36.))*!dtor
ac=(tcmp+irc*360./float(nrc) + 180*(irc mod 4))*!dtor
xplanet_position[irc,itc]=rc[irc]*cos(ac)+xc
yplanet_position[irc,itc]=rc[irc]*sin(ac)+yc 
dist_circle,dst,201,xplanet_position[irc,itc],yplanet_position[irc,itc]
occupied=where(dst le 2.5,nocc)
if nocc gt 0 then gridfakes[occupied]=1
endfor
endfor

;this is mostly for diagnostic purposes
;writefits,'gridoffakecomps.fits',gridfakes

for itc=0,ntc-1 do begin
suffname='reg_cal'
iternum=itc
for n=0,nfiles-1 do begin

        h0sci=headfits(datadir+files[n],ext=0,/silent)
        im=readfits(datadir+files[n],h1sci,/silent,ext=1)

        
    
;number of wavelength channels
;        nwvlh=(size(im,/dim))[2]

;set to empty data cube
        ;imff=im
        im[*]=1.d-15

;now, add planet to data cube
  numplanet=0
  outputcube=im

  print,'Adding Planets '+' to cube '+strtrim(n+1,2)+'/'+strtrim(nfiles,2)
  for irc= 0,nrc-1 do begin
   tcomp=(long(itc)/float(ntc))*360
   ;tcomp=itc*45.
    ;ac=(tcmp+(irc mod 10)*36.+(irc mod 2)*(180.-36.))*!dtor
   ;asc=(tcomp*90+(irc mod 2)*180.+0*(irc mod 2)*(180.)- allpa[n]+allpa[0])*!dtor
   ;asc=(tcomp+(irc mod 2)*180.+0*(irc mod 2)*(180.)- allpa[n]+allpa[0])*!dtor

   asc=(tcomp+irc*360./float(nrc) + 180*(irc mod 4) - allpa[n]+allpa[0])*!dtor
   if ~keyword_set(nonorthup) then asc+=northpa*!pi/180.

   ;asc=(tcomp+(irc mod itc/float(ntc))*360/float(itc+1)+(irc mod 2)*(180.-360/float(nrc))- allpa[n]+allpa[0])*!dtor
   ;asc=(tcomp+(irc mod 10)*36.+(irc mod 2)*(180.-36.)- allpa[n]+allpa[0])*!dtor
   ;asc=((tcomp mod 2)*45+(irc mod 4)*90.+0*(irc mod 2)*(180.-0.) - allpa[n]+allpa[0])*!dtor
   ;asc=((tcomp mod 2)*45+(irc mod 10)*36.+(irc mod 2)*(180.-90.) - allpa[n]+allpa[0])*!dtor
       xpos=rc[irc]*cos(asc)+xc  ;x position
       ypos=rc[irc]*sin(asc)+yc  ;y position
;   print,'xy is ',xpos,ypos,irc,itc
  ; print,'Adding Planet '+strtrim(numplanet+1,2)+'/'+strtrim(long(nrc*1L*ntc*1L),2)+' to cube '+strtrim(n,2)
   charis_insert_planet_into_cube,pfname,im,h0sci,h1sci,xpos,ypos,contrast,intr_psf,inputmodel,cube_out=outputcube0,spec_out=outputspec
   outputcube+=outputcube0
   numplanet+=1
  endfor

planet_input_spectrum[*,n]=outputspec
;register image

        writefits,datadir+filesfc[n],0,h0sci
        writefits,datadir+filesfc[n],outputcube,h1sci,/append
        im[*]=0
 endfor
;stop
;average input planet spectrum
planet_avg_input_spectrum=median(planet_input_spectrum,dimension=2,/even)
plot,lambda,planet_avg_input_spectrum,linestyle=0,psym=-4,xrange=[1,2.5]


    ;*** ADI treatment

if (zrsub gt 0) then begin
    charis_imrsub,pfname,prefname=prefname,suffname=suffname,rmax=lrmax,/prad,/fc

    suffname='rsub_fc'
endif



case reductype of 
 0: begin
charis_adiklip,pfname,prefname=prefname,rsubval=zrsub,$
nfwhm=znfwhm,drsub=zdrsub,npca=znpca,rmin=zrmin,rmax=zrmax,nzone=znzone,zero=zzero,meansub=zmeansub,meanadd=zmeanadd,$
/fc,/fwdmod,outfile='attenres'+strtrim(string(iternum,format='(i2)'),2)+'.fits'
    
    end
 4: begin
charis_rdiklip,pfname,prefname=prefname,rsubval=zrsub,$
refpath=refpath,$
drsub=zdrsub,npca=znpca,rmin=zrmin,rmax=zrmax,$
zero=zzero,meansub=zmeansub,meanadd=zmeanadd,$
/fc,/fwdmod,outfile='attenres'+strtrim(string(iternum,format='(i2)'),2)+'.fits'


    end

 else: begin

    print,'cannot find reduction method'
    print,'KLIP implementation is limited'
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
    h0cube=headfits(reducdir+'attenres'+strtrim(string(iternum,format='(i2)'),2)+'.fits')
    imcube=readfits(reducdir+'attenres'+strtrim(string(iternum,format='(i2)'),2)+'.fits',h1cube,ext=1)

    modelname='attenres'+strtrim(string(itc,format='(i2)'),2)+'.fits'
    modelbasename=(strsplit(modelname,'.',/extract,count=modcount))[0]

;wavelength-collapsed version.
    h0col=headfits(reducdir+'attenres'+strtrim(string(iternum,format='(i2)'),2)+'_collapsed.fits')
    imcol=readfits(reducdir+'attenres'+strtrim(string(iternum,format='(i2)'),2)+'_collapsed.fits',h1col,ext=1)
;endfor
;******stop at this point
;stop

;for now do just one planet
;nplanet=1
;***now Loop on Wavelength and then Planet to get Attenuation vs. Channel 
;res_flux=fltarr(nlambda,nrc)
   ;for irc=0,nrc-1 do begin
   for ilambda=0L,nlambda-1 do begin 
     imslice=imcube[*,*,ilambda]
     for irc=0,nrc-1 do begin
     openw,5,reducdir+'atten'+strtrim(itc,2)+'_'+strtrim(ilambda,2)+'.dat'
     xp0=xplanet_position[irc,itc]
     yp0=yplanet_position[irc,itc]
     rp0=sqrt((xp0-xc)^2.+(yp0-yc)^2.)

     ;compute the flux at the original position 
     aper,imslice,xp0,yp0,flux,eflux,sky,skyerr,1,aperrad[ilambda],[2,6]*aperrad[ilambda],setskyval=0,/flux,/exact,/nan,/silent

     ;now compare to the original, input flux to get the attenuation.
     res_flux[ilambda,irc,itc]=flux/planet_avg_input_spectrum[ilambda]
     throughput_map_lambda[xp0,yp0,ilambda]= res_flux[ilambda,irc]
     print,'Throughput for Planet ',strtrim(irc+1,2),'at radius ',strtrim(rc[irc],2),' at Wavelength ',string(lambda[ilambda]),' is ...',res_flux[ilambda,irc,itc]
     printf,5,rp0,xp0,yp0,res_flux[ilambda,irc,itc]
    close,5
    endfor
    ;close,5
   endfor

;res_flux_collapse=fltarr(nrc)
;****for collapsed cube***
;for iplanet=0L,nplanet-1 do begin
openw,5,reducdir+'atten_collapsed'+strtrim(itc,2)+'.dat'
for irc=0L,nrc-1 do begin
 xp0=xplanet_position[irc,itc]
 yp0=yplanet_position[irc,itc]
 rp0=sqrt((xp0-xc)^2.+(yp0-yc)^2.)
 rplanet_position[irc,itc]=rp0
aper,imcol,xp0,yp0,flux,eflux,sky,skyerr,1,0.5*median(fwhm,/even),[2,6]*median(fwhm,/even),setskyval=0,/flux,/exact,/nan,/silent
res_flux_collapse[irc,itc]=flux/median(planet_avg_input_spectrum,/even)
throughput_map[xp0,yp0]= res_flux_collapse[irc,itc]
printf,5,rp0,xp0,yp0,res_flux_collapse[irc,itc]
print,'Throughput for Planet ',strtrim(irc+1,2), 'at radius ',strtrim(rc[irc],2),' in Collapsed Cube is ...',res_flux_collapse[irc,itc]
endfor
close,5
endfor ;loop on ntc

;****Now you have populated the throughput maps, save them, and construct a radial profile***
writefits,'throughputmap.fits',throughput_map
writefits,'throughputmap_lambda.fits',throughput_map_lambda

dist_circle,rarray,dimx
goodarray=where(rarray ge lrmin-2.5 and rarray le lrmax+2.5,ngoodarray,complement=badarray)
attenmap=fltarr(dimx,dimy)
attenmap_lambda=fltarr(dimx,dimy,nlambda)
;radial profiles
;define a radial-profile of the throughput map
profrad_tc,throughput_map,median(fwhm,/even),lrmin,lrmax,p1d=pr,/nan,rayon=rsample
;attenmap=pr
attenmapo=interpol(pr,rsample,rarray)
attenmapo[badarray]=0
attenmap=reform(attenmapo,dimx,dimy)
for il=0L,nlambda-1 do begin
throughput_map_slice=throughput_map_lambda[*,*,il]
profrad_tc,throughput_map_slice,median(fwhm,/even),lrmin,lrmax,p1d=pr,/nan,rayon=rsample
attenmapslice=interpol(pr,rsample,rarray)
attenmapslice[badarray]=0
attenmap_lambda[*,*,il]=attenmapslice
endfor

;write output
writefits,'attenmap.fits',attenmap
writefits,'attenmap_lambda.fits',attenmap_lambda

skiptotheend:
end
