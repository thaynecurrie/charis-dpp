pro charis_sdialoci,pfname,prefname=prefname,suffname=suffname,nfwhm=nfwhm,drsub=drsub0,na=na,geom=geom,$
rsubval=rsubval,pixmask=pixmask,$
svd=svd,$
;nonorthup=nonorthup,$
meanadd=meanadd,$
zero=zero,$
nref=nref,$
savecoeff=savecoeff,$
usecoeff=usecoeff,$
channel=channel,$
debug=debug,$
postadi=postadi,$
rmin=rmin,rmax=rmax,$
fc=fc,fwdmod=fwdmod,$
norot=norot,angoffset=angoffset,help=help,outfile=outfile

;****Public, Pipeline-Release SDI/A-LOCI code****
;****03/05/2021***
;Version 2.0 -- more code cleanup, put in option to do postSDI ADI

;****12/2019***
;Version 1.1 - cleaned up code, fixed north-up treatment, help function added
;Version 1.0 - SDI/(A-)LOCI for CHARIS. Forward-modeling incorporated
;****09/02/2018***
;Version 0.2 - SDI/(A-)LOCI for CHARIS.   Code cleaned up.  Need to still incorporate forward-modeling.
;****08/31/2018***
;Version 0.1- SDI for CHARIS, quick and dirty to ensure it works.  do other stuff later

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,'charis_sdialoci.pro: ALOCI-SDI PSF subtraction method, Version 1.1'
print,'Written by T. Currie (2011-2014), adapted for CHARIS IFS Data (4/2017), updated 6/2020'
print,''
print,'**Calling Sequence**'
print,"charis_sdialoci,pfname,postadi=postadi,nfwhm=nfwhm,rmin=rmin,rmax=rmax,drsub=drsub0,na=na,geom=geom,svd=svd,nref=nref,pixmask=pixmask,meanadd=meanadd,zero=zero,rsubval=rsubval,"
print,"fwdmod=fwdmod,savecoeff=savecoeff,usecoeff=usecoeff,channel=channel,prefname=prefname,suffname=suffname,fc=fc"
print,"norot=norot,angoffset=angoffset,outfile=outfile"
print,''
print,'Example:'
print,"charis_sdialoci,'HR8799_low.info',/postadi,nfwhm=0.75,na=150,svd=2d-6,pixmask=1,nref=80,rmin=10,rmax=70,drsub=10,outfile='aloci.fits'"
print,''
print,"***Important Keywords***"
print,'*pfname - parameter file (e.g. HR8799_low.info)'
print,"*postadi - are you doing SDI on the post-ADI residuals [almost always 'yes', so set /postadi]"
print,"*nfwhm - rotation gap/exclusion zone (in lambda/D units)"
print,"*rmin - minimum radius of subtraction"
print,"*rmax - maximum radius of subtraction"
print,"*drsub - radial width of subtraction region"
print,"*na - optimization area (in units of PSF cores)"
print,"*geom - geometry of subtraction & optimization zones"
print,"*svd - SVD cutoff for covariance matrix inversion"
print,"*nref - construct a reference PSF from the 'nref'-most correlated frames"
print,"*pixmask - set to 1 to mask the subtraction zone (moving-pixel mask)"
print,"*meanadd - use a robust mean combination instead of median"
print,"*zero - do you subtract off the median of a region after PSF-subtracting?"
print,"*savecoeff - save the coefficients in a file [for forward-modeling]"
print,"*usecoeff - use saved coefficients [in forward-modeling]"
print,"*fwdmod - forward-modeling of a synthetic source [called in aloci_klip_fwdmod_planet/disk]"
print,"*channel - perform PSF subtraction only for this channel [integer value]"
print,"*outfile - the name of the output file"
print,"*norot - do NOT rotate images north-up [set this switch if you want to do ADI later]"

goto,endofprogram
endif

if ~keyword_set(angoffset) then angoffset=charis_get_constant(name='angoffset') ;nominally 2.2 deg


;Telescope Diameter for Subaru
Dtel=charis_get_constant(name='Dtel') ;7.9d0 ;visible pupil for SCExAO
;pixel scale 16.2 mas/pixel
pixscale=charis_get_constant(name='pixscale') ;0.0162 nominally

if ~keyword_set(outfile) then begin
if keyword_set(fc) then outfile='final_fc.fits'
if ~keyword_set(fc) then outfile='final.fits'
endif

;saving coefficients.
;coeff.dat

if keyword_set(savecoeff) then begin
;file_delete,'locicoeff.dat'
ff=file_search('sdilocicoeff.dat',count=filecount)
if filecount gt 0 then file_delete,'sdilocicoeff.dat'
openw,1,'sdilocicoeff.dat'
endif

if keyword_set(usecoeff) then readcol,'sdilocicoeff.dat',il_use,ir_use,it_use,nf_use,ck_use,format='i,i,i,i,f'

;drsub0=drsub
reducdir='./reduc/'

if ~keyword_set(rsubval) then begin
param,'rsub',rsubval,/get,pfname=pfname
endif

;determine reduction subdirectory
subdir='proc/'
;***edit: For now assume that you aren't doing radial profile subtraction

if (rsubval gt 0) then begin
reducdir1=reducdir+'rsub/'
endif else begin
reducdir1=reducdir+'reg/'
endelse

if keyword_set(postadi) then begin
reducdir1=reducdir+'proc/'
endif


datadir=reducdir1
reducdir+=subdir

;define a temporary directory
tmpdir=reducdir+'tmp/'
file_mkdir,tmpdir

;create list of filenames

param,'fnum_sat',flist,/get,pfname=pfname

;*** Prefixes***
if ~keyword_set(prefname) then prefname='n'
;***edit: again, assume no radial profile subtraction for now.

;****Suffixes****
;Overview: 
; a. normally just read in registered images (no keywords)
; b. read in spatially-filtered images (/rsub)
; c. read in ADI residual cubes (/postadi; note should throw flag of /norot in adi step)

;note: you need 
if (~keyword_set(suffname) and ~keyword_set(postadi)) then begin
if (rsubval gt 0) then begin
suffname='rsub'
endif else begin
test=file_search(datadir+'*reg_cal.fits')
if (n_elements(test) gt 0) then begin
suffname='reg_cal'  
endif else begin
suffname='reg'
endelse

endelse
endif

if (keyword_set(postadi) and ~keyword_set(suffname)) then suffname='_alocisub'

;*******


filenum=nbrlist(flist)

if ~keyword_set(fc) then begin
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
endif

if keyword_set(fc) then begin
if keyword_set(fwdmod) then begin
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
filesfwd=filelist(filenum,nfiles,prefix=prefname,suffix=suffname+'_fc')
endif else begin
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname+'_fc')
endelse

endif

filestmp=filelist(filenum,nfiles,prefix=prefname,suffix='_tmp',ext='.dat')

;**vestigal, keep in case you want to normalize by the radial profile (probably not)
filesprof=filelist(filenum,nfiles,prefix=prefname,suffix='_prof',ext='.fits')
;**

if ~keyword_set(fc) then begin
filesout=filelist(filenum,prefix=prefname,suffix='_sdialocisub')
endif

if keyword_set(fc) then begin
filesout=filelist(filenum,prefix=prefname,suffix='_sdialocisub_fc')
endif

nlist=indgen(nfiles)

;*************
;*************

;Define region for spider mask, the image FWHM, and the saturation radius
param,'fwhm',fwhm,/get,pfname=pfname ; I don't think this is needed.
param,'rsat',rsat,/get,pfname=pfname

;LOCI algorithm parameters
;nfwhm, drsub,na,and geom

if ~keyword_set(nfwhm) then begin
    param,'nfwhm',nfwhm,/get,pfname=pfname
    if nfwhm eq 0 then nfwhm=0.75
endif
if ~keyword_set(drsub0) then begin
    param,'dr',drsub0,/get,pfname=pfname
endif
if ~keyword_set(rmin) then begin
     param,'rmin',rmin,/get,pfname=pfname
endif
if ~keyword_set(rmax) then begin
     param,'rmax',rmax,/get,pfname=pfname
endif

if ~keyword_set(na) then begin
    param,'locina',na,/get,pfname=pfname
endif
if ~keyword_set(geom) then begin
    param,'locigeom',geom,/get,pfname=pfname
    if geom eq 0 then geom=1.
endif


;SVD
if ~keyword_set(svd) then begin
   param,'ALOCISVD',svd,/get,pfname=pfname
   svd=10^(1.0*svd)
endif

;Correlation-based Frame Selection
if ~keyword_set(nref) then begin
   param,'ALOCINRE',nref,/get,pfname=pfname
endif

if ~keyword_set(pixmask) then begin
   param,'pixmask',pixmask,/get,pfname=pfname
endif

if ~keyword_set(zero) then begin
   param,'zero',zero,/get,pfname=pfname
endif

if ~keyword_set(meanadd) then begin
  param,'meanadd',meanadd,/get,pfname=pfname
endif

;**get dim from first image header
print,'files',files[0]
test=readfits(reducdir1+files[0],/exten,h1)
h0=headfits(reducdir1+files[0],ext=0)

;get north-up value
northpa=sxpar(h0,'TOT_ROT',count=northcount)

if northcount eq 0 then begin
print,'no northup value'
stop
endif

dim=sxpar(h1,'naxis1')
xc=dim/2 & yc=dim/2
print,dim,xc,yc
print,size(test)

;Now Get the Wavelength Vector
filter=sxpar(h0,'FILTNAME')
get_charis_wvlh,h0,wavelengths
lambda=wavelengths*1.d-3

;**Which Telescope?? Latitude/Longitude
lat=double(sxpar(h0,'lat',count=latmatch))
lng=double(sxpar(h0,'lng',count=lngmatch))

;If you don't have an entry, assume you're at Maunakea
if latmatch eq 0 then lat = 19.825d0
if lngmatch eq 0 then lng = -155.4802d0


;**Exposure Time.  In imprep.pro, set all indiv. times to 'exp1time' and coadds to 'coadds'
exptime=sxpar(h0,'exp1time')
coadds=sxpar(h0,'coadds')

;if ~keyword_set(rmin) then rmin=(rsat-5)>5
;if ~keyword_set(rmax) then rmax=1.1*dim/2

;**For fits header keywords later
loci_nfwhm=nfwhm
loci_drsub=drsub0
loci_na=na
loci_geom=geom
loci_rmin=rmin
loci_rmax=rmax
loci_nref=nref

;load parallactic angles
readcol,'reduc.log',ffilenum,allxc,allyc,allrsat,allha,allpa,/silent

;Debugging
;print,total(ffilenum ne filenum),n_elements(filenum),n_elements(ffilenum)
;****in case there is a mismatch of reduc.log and the files in your .info file
if ((total(ffilenum ne filenum) gt 0) or (n_elements(filenum) ne n_elements(ffilenum)))  then begin
print,'ffilenum is',long(ffilenum),n_elements(ffilenum)
print,'filenum  is',filenum,n_elements(filenum)
stop
endif

nlambda=n_elements(lambda)
;Since lambda is defined, now you can define the temperorary files
filestmp=strarr(nlambda)
for il=0L,nlambda-1 do begin
filestmp[il]='lambda_'+nbr2txt(il,4)+'_tmp'+'.dat'
endfor

;determine radii
    nrsub=ceil((rmax-rmin)/drsub0)
    rsub=findgen(nrsub)*drsub0+rmin
    drsub=replicate(drsub0,nrsub)
drsub=drsub<(rmax-rsub)

;Array of FWHM
fwhm=1.0*(1.d-6*lambda/Dtel)*(180.*3600./!dpi)/pixscale

;array of distances and angles to determine indices in each section
distarr=shift(dist(dim),dim/2,dim/2)
ang=(angarr(dim)+2.*!pi) mod (2.*!pi)

;*****SDI Loop Structure *****
;Loop nfile (the data cube)
; Loop channel (channel il in cube)
;  - align speckles to channel
;  Loop radius
;   Loop sector
;    - subtract PSF
;    - save results in binary file
;   end Loop sector
;  end Loop radius
; end Loop channel
; rotate/align to true north
;end Loop nfile
;-open binary files, stitch together data cubes
;-combine data cubes together for master cube; combine master cube slices for master broadband image
;*******

for nf=0,nfiles-1 do begin   ;start nf

;Loop on images
print,'Reading in File '+strtrim(nf+1,2)+'/'+strtrim(nfiles,2)

;forward-modeling/saved coefficients

;use coeff
 if keyword_set(usecoeff) then begin
if (nf_use eq !null) eq 0 then begin
 coeffstouse=where(nf_use eq nf,ncoeffstouse)
 if ncoeffstouse gt 0 then begin
 c_use2=ck_use[coeffstouse]
 il_use2=il_use[coeffstouse]
 ir_use2=ir_use[coeffstouse]
 it_use2=it_use[coeffstouse]
 endif
 endif
 endif



;****probably should put this outside the n_file loop.  
;estimates the largest optimization radius needed
rimmax=0.
for ir=0,nrsub-1 do begin
    r=rsub[ir]
;for now, for simplicity, just assume that the area is the largest fwhm
    area=na*!pi*(max(fwhm)/2.)^2
    ;width of optimization radius desired
    dropt=sqrt(geom*area)
    nt=round((2*!pi*(r+dropt/2.)*dropt)/area)>1
    dropt=sqrt(r^2+(nt*area)/!pi)-r
    rimmax>=r+dropt
endfor ;end ir

;for SDI
rimmax<=1.1*dim/2
;****

;cut the image of the rings 5 pixels wide
;save in a file

drim=5.

nrim=ceil((rimmax-rmin)/drim)
rim=findgen(nrim)*drim+rmin
;determine indices of pixels included in each ring
;DRIM of pixels and save them to disk

rimmax_edge=65 ; the nominal edge of the image, where pixels interior are on detector

;***wavelength loop***

for il=0L,nlambda-1 do begin   ;start il

;***define an array that tells you maximum extent of scaled/magnified slices aligned to slice il (not needed?; save just in case)
rimmax_slice=rimmax_edge*lambda[il]/lambda

;use coeff
 if keyword_set(usecoeff) then begin
 if (il_use2 eq !null) eq 0 then begin
 coeffstouse=where(il_use2 eq il,ncoeffstouse)
 if ncoeffstouse gt 0 then begin
 c_use3=c_use2[coeffstouse]
 ir_use3=ir_use2[coeffstouse]
 it_use3=it_use2[coeffstouse]
 endif
 endif
 endif

for ir=0,nrim-1 do begin
    ri=rim[ir] & rf=ri+drim
    ia=where(distarr lt rf and distarr ge ri)
    openw,funit,tmpdir+'indices_a'+nbr2txt(ir,3)+'.dat',/get_lun
    writeu,funit,ia
    free_lun,funit
endfor ; end ir

;Cutting images into rings and place rings even nfile radius 
;in a single file

noise_im=fltarr(nrim,nlambda)


;Align Speckles

    h0=headfits(reducdir1+files[nf],/silent,ext=0)
    im_tmp=(readfits(reducdir1+files[nf],h1,/exten,/silent))

    charis_alignspeckle,im_tmp,h0,h1,refslice=il,/nolocs


    ;not needed? 
    ;bad=where(im_tmp eq 0 and distarr gt max(rimmax_slice),nbad)
    ;if nbad gt 0 then im_tmp[bad]=!values.f_nan
    for ig=0L,nlambda -1 do begin
    islice=im_tmp[*,*,ig]
    if ig ne il then begin
    bad=where(finite(islice) eq 0,nbad)
    if nbad gt 0 then islice[bad]=0
    im_tmp[*,*,ig]=islice
    endif
    endfor

    writefits,tmpdir+'align_'+strtrim(il,2)+files[nf],0,h0
    writefits,tmpdir+'align_'+strtrim(il,2)+files[nf],im_tmp,h1,/append


    if keyword_set(fwdmod) then begin
    im_fwd_tmp=(readfits(reducdir1+filesfwd[nf],/exten,horig,/silent))
    h0orig=headfits(reducdir1+filesfwd[nf])
    charis_alignspeckle,im_fwd_tmp,h0orig,horig,/nolocs,refslice=il
    
    writefits,tmpdir+'align_fwd_'+strtrim(il,2)+filesfwd[nf],0,h0orig
    writefits,tmpdir+'align_fwd_'+strtrim(il,2)+filesfwd[nf],im_fwd_tmp,horig,/append

    endif 

 for islice=0L,nlambda -1 do begin

    im=im_tmp[*,*,islice]
    if keyword_set(fwdmod) then im_fwd=im_fwd_tmp[*,*,islice]

    for ir=0,nrim-1 do begin
        ia=read_binary(tmpdir+'indices_a'+nbr2txt(ir,3)+'.dat',data_type=3)
        if ia[0] eq -1 then continue
        openw,funit,tmpdir+'values_a'+nbr2txt(ir,3)+'.dat',/get_lun,append=(islice gt 0)
        writeu,funit,float(im[ia])
        free_lun,funit


         if keyword_set(fwdmod) then begin
         openw,funit,tmpdir+'values_fwd'+nbr2txt(ir,3)+'.dat',/get_lun,append=(islice gt 0)
         writeu,funit,float(im_fwd[ia])
         free_lun,funit
         endif

    endfor

;**** For now, we don't care about field rotation.  We'll put that back in some day later.

endfor ;end islice

;DEFINE wavelength displacement for interior of each annulus
mvmts=(lambda[il]/lambda-1d0)#rsub

;quick diagnostic check (uncomment as needed)
;for ir=0L,n_elements(rsub)-1 do begin
; for ill=0L,n_elements(lambda)-1 do begin
;   print,mvmts[ill,ir],lambda[ill],rsub[ir],lambda[il]
; endfor
;endfor
;stop

;MAIN Radius LOOP      

iaim_loaded=-1
for ir=0,nrsub-1 do begin
    ri=rsub[ir] & dr=drsub[ir] & r=ri+dr/2. & rf=ri+dr

    print,'Image '+strtrim(nf+1,2)+'/'+strtrim(nfiles,2),' Wavelength ' + strtrim(il+1,2)+'/' + strtrim(nlambda,2) + ' Annulus '+strtrim(ir+1,2)+'/'+strtrim(nrsub,2)+' with radius '+$
        string(r,format='(f5.1)')+$
        ' [>='+string(ri,format='(f5.1)')+', <'+string(rf,format='(f5.1)')+']...'

    ;print,'Image '+strtrim(nf+1,2)+'/'+strtrim(nfiles,2),' Annulus '+strtrim(ir+1,2)+'/'+strtrim(nrsub,2)+' with radius '+$
        ;string(r,format='(f5.1)')+$
        ;' [>='+string(ri,format='(f5.1)')+', <'+string(rf,format='(f5.1)')+']...'


;forward-modeling/saved coefficients
  if keyword_set(usecoeff) then begin
     if (ir_use3 eq !null) eq 0 then begin
     coeffstouse=where(ir_use3 eq ir,ncoeffstouse)
     if ncoeffstouse gt 0 then begin
     c_use4=c_use3[coeffstouse]
     it_use4=it_use3[coeffstouse]
     endif
     endif
  endif

;****
;for now, for simplicity, set FWHM = max(FWHM).  So area grows with wavelength.
    area=na*!pi*(max(fwhm)/2.)^2.

    ;width of desired optimization annulus
    dropt=sqrt(geom*area)

    if dropt lt dr then begin
        print,'dropt < drsub !!!'
        print,'dropt: ',dropt
        print,'drsub: ',dr
        stop
    endif

    ;***determining the area optimization for this annulus removal
    ;for region removed in early reg_optimization

        r1opt=ri
        ;number of annulus section
        nt=round((2*!pi*(r1opt+dropt/2.)*dropt)/area)>1
        
        ;dropt for annulus with sections of exact area

        dropt=sqrt(r1opt^2+(nt*area)/!pi)-r1opt
        r2opt=r1opt+dropt

        if r2opt gt rim[nrim-1]+drim then begin
            r2opt=rim[nrim-1]+drim
            dropt=r2opt-r1opt
;check this
            nt=round((2*!pi*(r1opt+dropt/2.)*dropt)/area)>1
        endif

    ;determines what image to load into memory
    i1aim=floor((r1opt-rmin)/drim)
    i2aim=floor((r2opt-rmin)/drim)
        ;print,i2aim,i1aim,n_elements(rim),rim[i2aim],sqrt(geom*area)
    if i2aim eq nrim then i2aim-=1
    if rim[i2aim] eq r2opt then i2aim-=1
    iaim=indgen(i2aim-i1aim+1)+i1aim

    ;removes the annuli that aren't necessary

    if ir gt 0 then begin
        irm=where(distarr[ia] lt rim[i1aim] or distarr[ia] ge rim[i2aim]+drim,crm,complement=ikp)
        if crm gt 0 then remove,irm,ia
        if crm gt 0 then annuli=annuli[ikp,*]

         if keyword_set(fwdmod) then begin
         if crm gt 0 then annuli_fwd=annuli_fwd[ikp,*]
         endif

    endif

    ;instructs the missing annuli
    iaim_2load=intersect(iaim,intersect(iaim,iaim_loaded,/xor_flag))

    c2load2=where(iaim_2load ge 0,c2load)

    ;*debugging*
    ;if(c2load le 0)then begin
    ;if(ir eq 15)then begin
    ; c2load = 1
    ; iaim_2load=0
    ;endif

    for k=0,c2load-1 do begin
    
        ia_tmp=read_binary(tmpdir+'indices_a'+nbr2txt(iaim_2load[k],3)+'.dat',data_type=3)
        annuli_tmp=read_binary(tmpdir+'values_a'+nbr2txt(iaim_2load[k],3)+'.dat',data_type=4)

        if keyword_set(fwdmod) then $
         annuli_fwd_tmp=read_binary(tmpdir+'values_fwd'+nbr2txt(iaim_2load[k],3)+'.dat',data_type=4)

        annuli_tmp=reform(annuli_tmp,n_elements(ia_tmp),nlambda)
        if keyword_set(fwdmod) then $
         annuli_fwd_tmp=reform(annuli_fwd_tmp,n_elements(ia_tmp)*1L,nlambda*1L)
        
        if ir+k eq 0 then ia=ia_tmp else ia=[ia,ia_tmp]
        if ir+k eq 0 then annuli=annuli_tmp else annuli=[annuli,annuli_tmp]

        if keyword_set(fwdmod) then begin
        if ir+k eq 0 then annuli_fwd=annuli_fwd_tmp else annuli_fwd=[annuli_fwd,annuli_fwd_tmp]
        endif

    endfor

    ia_tmp=0 & annuli_tmp=0

    ;remembers the list of annuli changes
    iaim_loaded=iaim

    ;indices of pixels for optimization annulus
    iaopt=where(distarr[ia] ge r1opt and distarr[ia] lt r2opt)

    if keyword_set(debug) then begin
        ;print,
        testimagef=fltarr(dim,dim)
        testimagef[ia[iaopt]]=1
        writefits,'testimagef.fits',testimagef
    stop
    endif

    ;angle of annulus sections
    dt=2.*!pi/nt

    ;loop on angular sections

    for it=0,nt-1 do begin
        ;indices of pixels included in this section: i.e. the optmization region



;***NOTE: the nomenclature here works but is wonky:  'noiopt' determines 'opt' a few lines later. 'iopt' determines 'isub'.
        iopt=where(ang[ia[iaopt]] ge it*dt and ang[ia[iaopt]] lt (it+1)*dt,complement=noiopt)
        npix=n_elements(iopt)

        ;if npix lt 5 then continue
  
        iopt=iaopt[iopt]
        noiopt=iaopt[noiopt]

        ;if keyword_set(debug) then begin
        ;print,
        ;testimage=fltarr(dim,dim)
        ;testimage[ia[iopt]]=1
        ;testimage[ia[noiopt]]=2
        ;testimage[ia]+=4
        ;writefits,'testimagef2.fits',testimage
        ;endif

        ;indices of pixels to subtract

        isub=where(distarr[ia[iopt]] ge ri and distarr[ia[iopt]] lt rf)

        if (pixmask gt 0) then begin
        ;if keyword_set(pixmask) then begin
        ioptf=where(distarr[ia[noiopt]] ge ri and distarr[ia[noiopt]] lt rf)
        ioptf=noiopt[ioptf]
        endif else begin
        ioptf=where(distarr[ia] ge ri and distarr[ia] lt rf)
        ;ioptf=iopt[ioptf]
        ;ioptf=ia[ioptf]
        endelse

         if n_elements(isub) lt 2 then continue
        isub=iopt[isub]

        optreg=annuli[ioptf,*]
        if keyword_set(fwdmod) then begin 
         optreg_fwd=annuli_fwd[ioptf,*]
         optreg+=annuli_fwd[ioptf,*]
        endif
 
        ;nanmask=longarr(nlambda)
        nanmask=total(optreg,1)
        ;bad=where(finite(bah) eq 0)
        if keyword_set(debug) then begin
        testimage=fltarr(dim,dim)
        testimage[ia[isub]]=2
        testimage[ia[ioptf]]=1
        ;testimage[ia[iaopt[noiopt]]]=3

        writefits,'testimage3b.fits',testimage
        stop
        endif


;***** filtering of deviant pixels originally written for ADI version.   Not in use here.  Use NaN mask +SVD truncation instead.
        ;always keep isub defined here (before removing deviant pixels)
        ;otherwise bright sources (which are identified as deviant pixels)
        ;would be masked out in the result

        ;removed from the region opt pixels or there is a NAN or very 
        ;deviant point in at least one annulus

;        z=finite(optreg)

        ;the following three lines are equivalent to the following
        ;inoise=floor((distarr(ia[iopt])-rim[0])/drim)#replicate(1,90)
        ;tmp=abs(optreg/noise_im[inoise])
        ;z<=(tmp lt 15.)
;        z<=(abs(optreg/noise_im[floor((distarr(ia[iopt])-rim[0])/drim)#replicate(1,nfiles)]) lt 15.)

        z=optreg
        z/=(replicate(1,n_elements(iopt))#median(abs(z),dim=1))
        z/=(median(abs(z),dim=2)#replicate(1,nlambda))
        ;z/=(median(abs(z),dim=2)#replicate(1,nfiles))
        ;z=(abs(z) lt 7 and finite(z) eq 1)
        ;z=(finite(z) eq 1)

        zq=median(z,dimension=2,/even)

        ;for a given pixel [i,*] look to see whether an image pixel that is deviant

         igoodf=where(zq ne 0 and finite(zq) ne 0,cgood)

if keyword_set(postadi) then begin

        optreg=optreg[igoodf,*]
        iopt=iopt[igoodf,*]
        if keyword_set(fwdmod) then optreg_fwd=optreg_fwd[igoodf,*]
endif

;******
        ;if cgood lt 2 then continue
        ;optreg=optreg[igood,*]
        ;iopt=iopt[igood]

        ;there clues to avoid a build images later

        openw,lunit,tmpdir+'indices_images.dat',/get_lun,append=(ir+it gt 0)
        writeu,lunit,ia[isub]
        free_lun,lunit

       ;build large matrix of a linear system to solve
       ;aa = the covariance matrix without any frame selection/filtering of NaN'd frames

       aa=optreg##transpose(optreg)

;*******Define the movement - dpix_lam - of a point source at r in slice il from other slices due to magnification
            dpix_lam=abs(mvmts[*,ir]-mvmts[il,ir])

            ;OFFSET determined enough images for subtraction


            indim=where(dpix_lam gt (nfwhm*fwhm[il]) and finite(nanmask) ne 0,c1)
            ;indim=where(dpix_lam gt (nfwhm*fwhm[il]) and rimmax_slice gt 1.0*rf and finite(nanmask) ne 0,c1)

            ;correlation-based frame selection (A-LOCI) ; remove for public code
            ;if nref le c1 then begin
            ;index_corr=(reverse(sort(aa[indim,nf])))[0:nref-1 < c1]
            ;indim=indim[index_corr]
            c1=n_elements(indim)
            ;endif


            ;print,'nelements',c1,nfwhm,fwhm/ri*!radeg,max(dpa),dtpose[nf]
            igood=where(finite(annuli[isub,il]) eq 1,c2)
            if keyword_set(fwdmod) then igood_opt_fwd=where(finite(annuli[iopt,il]) eq 1,c2fwd)
            if c1 eq 0 or c2 lt 5 then begin
                ;*debug*
                ;diff=fltarr(n_elements(isub))*0
                ;if ~keyword_set(postadi) then begin
                diff=fltarr(n_elements(isub))+!values.f_nan
                difffwd=diff
                ;endif else begin
                ;diff=annuli[isub,il]
                ;endelse

            endif else begin

        ;forward-modeling/saved coefficients
  if keyword_set(usecoeff) then begin
 ;coeffstouse =where(il_use eq il and ir_use eq ir and it_use eq it and nf_use eq nf)
  if (it_use4 eq !null) eq 0 then begin
  coeffstouse=where(it_use4 eq it,ncoeffstouse)
  
  if ncoeffstouse gt 0 then begin
  c=c_use4[coeffstouse]
  endif
  endif
  if ~keyword_set(fwdmod) then goto,skipmatrixinversion
  endif

                if ~keyword_set(fwdmod) then begin
                ;matrix of linear system to solve
                a=(aa[indim,*])[*,indim]
                ;vector b of a linear system to solve
                b=aa[indim,il]
                
                ;solve the system
                ;c=invert(a,/double)#b

                ;*****pseudo-inverse of covariance matrix a
                if (keyword_set(svd) and n_elements(a) gt 2) then begin
                inv_a=svd_invert(a,svd,/double)
                c=inv_a#b
                endif else begin
                c=invert(a,/double)#b
                endelse

                ;construct the reference
                skipmatrixinversion:
                 
                ref=fltarr(n_elements(isub))
                for k=0,c1-1 do begin 
                  ref[igood]+=c[k]*annuli[isub[igood],indim[k]]
                  if keyword_set(savecoeff) then printf,1,long(il),long(ir),long(it),long(nf),c[k]
                endfor

                ;make the difference
                ;reftot=total(ref) 
                
                diff=annuli[isub,il]-ref
    
               if (zero gt 0) then diff-=median(diff,/even)
               ;if keyword_set(zero) then diff-=median(diff,/even)
                
                ;test
                ;diff=annuli[isub,il]
                ;diff=ref
                skipme:

              endif else begin
  
                 ;if you are doing forward modeling then ...
                 ;   - you already have the set of coefficients c_k
                 ;   - you need to solve for the perturbing coefficients beta
                  
                 ref=fltarr(n_elements(isub))
                 refpert=fltarr(n_elements(ioptf))
                 ;self-subtraction
                  ;coeffgood=where((it_use eq it) and (ir_use eq ir) and (il_use eq il) and (nf_use eq nf))
                  ;c=ck_use[coeffgood]
;                  print,c
;                  print,ri,nf,ir,it,il
                  ;stop
                 for k=0,c1-1 do begin
                  ;ref[igood]+=c[k]*annuli_fwd[isub[igood],indim[k]]
                  ;refpert[igood]+=c[k]*annuli_fwd[iopt[igood],indim[k]]
                  refpert[ioptf]+=c[k]*annuli_fwd[ioptf,indim[k]]
                 endfor
  
  ;just self-subtraction
                ; diff1=annuli_fwd[isub,il]-ref
                ; diff=diff1
;goto,skipoverpert
   ;the first part of the perturbed column vector b'
                 pert1=(annuli[ioptf,*])[*,indim]
                 ;pert1=(annuli[iopt,*])[*,indim]
  ;the second part of the perturbed column vector b'
                 pert2=annuli_fwd[ioptf,il]-refpert[ioptf]
                 ;pert2=annuli_fwd[iopt[igood_opt_fwd],il]-refpert[igood_opt_fwd]
  
  ;now based on the above self-subtraction, introduce as a perturbation on the linear system
                 ;previous A covariance matrix is the same as the original
                 ;vector A of linear system to solve
                 a=(aa[indim,*])[*,indim]
  
                 ;for every l in column, b,   b_l= sum(pixels)_ annuli_orig(pixels,set of ref images) mult by s     elf-subtracted planet signal(pixels).
                 ;part1##(transpose(part2) t    o mult and sum over pixels
                 b=pert2##transpose(pert1)
 
                 if (keyword_set(svd) and n_elements(a) gt 2) then begin
                  inv_a=svd_invert(a,svd,/double)
                  pert_c=inv_a#b
                  endif else begin
                  pert_c=invert(a,/double)#b
                 endelse

  
                 ref2=fltarr(n_elements(isub))
                 ;diff=ref
                 for k=0,c1-1 do begin
                   if (finite(pert_c[k]) ne 0) then begin
  
                   ;ref2[igood]+=pert_c[k]*annuli[isub[igood],indim[k]]
                   ;ref[igood]+=c[k]*annuli_fwd[isub[igood],indim[k]]+pert_c[k]*annuli[isub[igood],indim[k]]
                   ref[igood]+=c[k]*annuli_fwd[isub[igood],indim[k]]+pert_c[k]*annuli[isub[igood],indim[k]]+pert_c[k]*annuli_fwd[isub[igood],indim[k]]
                   endif
                 endfor
skipoverpert:
                diff=annuli_fwd[isub,il]-ref
                if (zero gt 0) then diff-=median(diff,/even)


              endelse   ; for fwdmod
              endelse   ;for enough pix to do psf sub.
        
            ;register the difference, add (append) the values of this annulus to 
            ;binary file image
            openw,lunit,tmpdir+filestmp[il],/get_lun,append=(ir+it gt 0)
            writeu,lunit,diff
            free_lun,lunit
    endfor  ;end loop it
endfor   ;end loop ir

;deletes files. .dat annuli
file_delete,file_search(tmpdir,'indices_a*.dat')
file_delete,file_search(tmpdir,'values_a*.dat')

endfor ; end loop il/wavelength?

;reading signs of pixels removed

ind=read_binary(tmpdir+'indices_images.dat',data_type=3)

;delete temporary indices
file_delete,tmpdir+'indices_images.dat'

;rebuild and turn images

 for il=0L,nlambda -1 do begin
    print,'Wavelength '+strtrim(il+1,2)+'/'+strtrim(nlambda,2)+': '+files[nf]+'...'

    ;reconstruct image
    print,' reconstruction...'

    im=make_array(dim,dim,type=4,value=!values.f_nan)
    im[ind]=read_binary(tmpdir+filestmp[il],data_type=4)

    ;delete temporary files
    ;file_delete,tmpdir+filestmp[il]

    ;get header
    h1=headfits(reducdir1+files[nf],/exten,/silent)
    h0=headfits(reducdir1+files[nf],ext=0,/silent)

;******astrometry*****

;******default is to rotate north-up using the tot-rot fits header keyword in the primary header
;throw /norot  switch if you do not want to derotate the image to a common value (e.g., if you are doing SDI after this step)

if ~keyword_set(norot) then begin
theta=-1*allpa[nf]

charis_northup,im,h0,h1,angoffset=angoffset
endif else begin

endelse

    if keyword_set(norot) then begin

    endif else begin
    sxaddhist,'A rotation of '+string(theta,format='(f8.3)')+$
      ' degrees was applied to the image.',h1
    endelse
; Add LOCI keywords
    sxaddpar,h1,'rsub',rsubval
    sxaddpar,h1,'loci_nfwhm ',loci_nfwhm
    sxaddpar,h1,'loci_na ',loci_na
    sxaddpar,h1,'loci_geom ',loci_geom
    sxaddpar,h1,'loci_drsub ',loci_drsub
    sxaddpar,h1,'loci_rmin ',loci_rmin
    sxaddpar,h1,'loci_rmax ',loci_rmax
    sxaddpar,h1,'svd ',svd
    sxaddpar,h1,'pixmaskv',pixmask
    sxaddpar,h1,'zero',zero
    sxaddpar,h1,'meanadd',meanadd

    suffix1='-loci'+strcompress(string(il),/REMOVE_ALL)
    writefits,tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix1+'.fits',0,h0
    writefits,tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix1+'.fits',im,h1,/append
endfor
endfor ;nfiles

;Okay, now combine the images together, construct datacubes, and then construct a combined datacube.
imt=dblarr(dim,dim,nlambda)
suffix0='-loci'

for nf=0,nfiles-1 do begin
 for il=0,nlambda-1 do begin
 h0=headfits(tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix0+strcompress(string(il),/REMOVE_ALL)+'.fits',/SILENT,ext=0)
 imt[*,*,il]=readfits(tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix0+strcompress(string(il),/REMOVE_ALL)+'.fits',/SILENT,/exten,h1)
 endfor

 writefits,reducdir+filesout[nf],0,h0
 writefits,reducdir+filesout[nf],imt,h1,/append
endfor

if ~keyword_set(mean) then begin
im=charis_combinefits(filesout,dir=reducdir,/cube)
endif else begin

im=charis_combinefits(filesout,/mean,dir=reducdir,/cube)
endelse


h0=headfits(reducdir+filesout[0],ext=0)
h1=headfits(reducdir+filesout[0],ext=1)

outname=reducdir+outfile
if keyword_set(fc) and ~keyword_set(outfile) then outname=reducdir+'final_fc.fits'
writefits,outname,0,h0
writefits,outname,im,h1,/append

;northup, regular cube and collapsed cube
im=readfits(outname,h1,/exten)
h0=headfits(outname,ext=0)

outpref=strsplit(outfile,'.',/extract)
outname_col=outpref[0]+'_collapsed'+'.fits'

h0_col=h0
h1_col=h1


if (meanadd gt 0) then begin
resistant_mean,im,3,im_collapsed,numrej,dimension=3
endif else begin
im_collapsed=median(im[*,*,*],/even,dimension=3)
endelse


writefits,outname,0,h0
writefits,outname,im,h1,/append

writefits,reducdir+outname_col,0,h0
writefits,reducdir+outname_col,im_collapsed,h1,/append

;file_delete,file_search(tmpdir,'*-loci*fits')

if (keyword_set(fwdmod) or keyword_set(savecoeff) or keyword_set(usecoeff)) then free_lun,1

endofprogram:
;close,/all
end
