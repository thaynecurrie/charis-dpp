pro charis_aloci_negspec_disk,pfname,reducname=reducname,sdi_reducname=sdi_reducname,reductype=reductype,prefname=prefname,$
nonorthup=nonorthup, diskmodel=diskmodel

;Version 1.1 - added documentation
;****insert negative copy of disk into sequence, perform PSF subtraction on the disk-free cubes
;*****CAN ONLY DO ADI RIGHT NOW****
;****06/18/2018**
;Version 1.0 - negative planet iteration to hone forward-modeling solution
;****

;*****************************

;setupdir,reducdir=reducdir

;determine reduction subdirectory
;pos1=strpos(pfname,'/')+1
;pos2=strpos(pfname,'.',/reverse_search)

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"charis_aloci_negspec_disk: A-LOCI negative disk subtraction, Version 1.0 (May 2021)"
print,'Written by T. Currie (2014), adapted for CHARIS IFS Data (2019)'
print,''
print,"STILL IN DEVELOPMENT!!!!"
print,""
print,"charis_aloci_negspec_planet,pfname,reducname=reducname,sdi_reducname=sdi_reducname,"
print,"diskmodel=diskmodel"
print,""
print,'Example (ADI):'
print,"charis_aloci_negspec_disk,'HR8799_low.info',reducname='aloci.fits',method=0"
print,""
print,"***Important Keywords****"
print,""
print,"*pfname - parameter file (e.g. HR8799_low.info)"
print,"*reducname - name of ADI-reduced cube"
print,"*diskmodel - the name of the disk model (CHARIS cube)"
;print,"*filt - do you spatially filter the cube first?"
;print,"*asdi - Did you do SDI on the ADI residuals?"
;print,"*sdi_reducname - name of post-ADI, SDI reduced cube [if /asdi is used]"
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
if keyword_set(sdi) then reductype = 1 ;set to SDI
;if keyword_set(asdi) then reductype = 4 ;set to ADIpSDI

;****to implement later
;if keyword_set(adipsdi) then reductype = 2 ;set to ADIpSDI
;if keyword_set(sdipadi) then reductype = 3 ;set to SDIpADI
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


;The disk model

;this is a 3D cube
diskmodel=readfits(diskmodel,/ext)

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

nrc=1

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
angoffset=northpa

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


nsc=21 ;nombre de sous-compagnon

;*****Keywords for RA,DEC, individual exposure time (exp1time), and coadds (so exp1time*ncoadds = exptime)*****
    param,'RA',radeg,/get,pfname=pfname
    param,'DEC',decdeg,/get,pfname=pfname

;loop calculates the angles and the self-subtraction
;print,'ntc is ',ntc

;add in fake disks here
    ;constructed images false companions

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
        ;im[*]=0
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
            bad=where(im[*,*,il] eq 0)
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
                psfs=-1*shift_sub(diskmodel[*,*,il],dxf,dyf)
                ;psfs[bad]=0
                psfsrot=rotat(psfs,(-1*northpa))/nsc
                writefits,'psf.fits',psfsrot
                psfsrot[bad]=0
                psfsrot[where(finite(psfsrot) eq 0)]=0
;If you have an IFS data format written by normal human beings then the first line is right.
                im[*,*,il]+=psfsrot

;(rotat(psfs,(-1*northpa))/nsc)
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



itc=0

case reductype of

 0: begin   ;just ADI
;A-LOCI with Forward-Modeling

 charis_adialoci,pfname,prefname=prefname,rsubval=zrsub,$
;charis_adialocif,pfname,prefname=prefname,rsubval=zrsub,zonetype=zadiztype,$
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,svd=zsvd,nref=znref,pixmask=zpixmask,zero=zzero,meanadd=zmeanadd,$
/np,outfile='ressub'+string(itc,format='(i1)')+'_np_'+'.fits'

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

;SKIP OVER ALL THIS
goto,skipoverallthis

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

skipoverallthis:
 
skiptotheend:

end
