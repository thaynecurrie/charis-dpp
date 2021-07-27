pro charis_adialoci,pfname,prefname=prefname,nfwhm=nfwhm,drsub=drsub0,na=na,geom=geom,$
svd=svd,$
pixmask=pixmask,$
meanadd=meanadd,$
zero=zero,$
nref=nref,$
savecoeff=savecoeff,$
usecoeff=usecoeff,$
channel=channel,$
rsubval=rsubval,$
suffname=suffname,$
postsdi=postsdi,$
rmin=rmin,rmax=rmax,$
fc=fc,fwdmod=fwdmod,$
norot=norot,angoffset=angoffset,$
outfile=outfile,$
 guide=guide,help=help

;****Public, Pipeline-Release ADI/A-LOCI code****

;***3/6/2021**
;Version 2.4 - now can do ADI on post-SDI cubes
;***10/6/2020**
;Version 2.3 - renamed code from charis_sublocirx to charis_adialoci
;***06/04/2020**
;Version 2.2 - cleaned up better.  
;***06/11/2018**
;Version 2.1 - now can do single channel reductions (useful for iterative nulling/neg. planets)

;***02/08/2018***
;Version 2.0 - incorporates forward-modeling using perturbed coefficients (T. Currie 2018, in prep.)

;***02/01/2018***
;Version 1.2 - some code cleanup, changed way we do north PA rotation.

;***10/16/2017***
;Version 1.1
;-cleaned up code syntax, removing unnecessary legacy LOCI keywords not used
;- to do: better version of correlation matrix-based frame selection
;         switch to change the optimization/subtraction zone geometry
;         Tikhonov Regularization/identity matrix damping
;         model PSF/lambda function damping
;***04/08/2017***
;Version 1.0
;adapted for CHARIS, for now

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,'charis_adialoci.pro: ALOCI PSF subtraction method, Version 2.2'
print,'Written by T. Currie (2011-2014), adapted for CHARIS IFS Data (4/2017), updated 6/2020'
print,''
print,'**Calling Sequence**'
print,"charis_adialoci,pfname,prefname=prefname,postsdi=postsdi,nfwhm=nfwhm,rmin=rmin,rmax=rmax,drsub=drsub0,na=na,geom=geom,svd=svd,nref=nref,pixmask=pixmask,meanadd=meanadd,zero=zero,rsubval=rsubval,"
print,"fwdmod=fwdmod,savecoeff=savecoeff,usecoeff=usecoeff,channel=channel,prefname=prefname,suffname=suffname,fc=fc"
print,"norot=norot,angoffset=angoffset,outfile=outfile"
print,''
print,'Example:'
print,"charis_adialoci,'HR8799_low.info',nfwhm=0.75,na=150,svd=2d-6,nref=80,rmin=10,rmax=70,drsub=10,outfile='aloci.fits'"
print,''
print,"***Important Keywords***"
print,'*pfname - parameter file (e.g. HR8799_low.info)'
print,"*postsdi - are you doing ADI on the post-SDI residuals [so set /postsdi]"
print,"*nfwhm - rotation gap/exclusion zone (in lambda/D units)"
print,"*rmin - minimum radius of subtraction"
print,"*rmax - maximum radius of subtraction"
print,"*drsub - radial width of subtraction region"
print,"*rsubval - set to 1 to use spatially filtered data [usually the right decision]"
print,"*na - optimization area (in units of PSF cores)"
print,"*geom - geometry of subtraction & optimization zones"
print,"*svd - SVD cutoff for covariance matrix inversion"
print,"*nref - construct a reference PSF from the 'nref'-most correlated frames"
print,"*pixmask - set to 1 to mask the subtraction zone (moving-pixel mask)"
print,"*meanadd - use a robust mean combination instead of median"
print,"*zero - do you subtract off the median of a region after PSF-subtracting?"
print,"*norot - do NOT rotate images north-up [set this switch if you want to do SDI later]"
print,"*savecoeff - save the coefficients in a file [for forward-modeling]"
print,"*usecoeff - use saved coefficients [in forward-modeling]"
print,"*fwdmod - forward-modeling of a synthetic source [called in aloci_klip_fwdmod_planet/disk]"
print,"*channel - perform PSF subtraction only for this channel [integer value]"
print,"*outfile - the name of the output file"

goto,endofprogram
endif

starttime=systime(/seconds)

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
ff=file_search('locicoeff.dat',count=filecount)
if filecount gt 0 then file_delete,'locicoeff.dat'
openw,1,'locicoeff.dat'
endif

if keyword_set(usecoeff) then readcol,'locicoeff.dat',il_use,ir_use,it_use,nf_use,ck_use,format='i,i,i,i,f'

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

if keyword_set(postsdi) then begin
reducdir1=reducdir+'proc/'
endif

datadir=reducdir1
reducdir+=subdir

;define a temporary directory
tmpdir=reducdir+'tmp/'
file_mkdir,tmpdir

;create list of filenames

;param,'obsdate',date,/get,pfname=pfname & date=strtrim(date,2)

param,'fnum_sat',flist,/get,pfname=pfname

;*** Prefixes***
if ~keyword_set(prefname) then prefname='n'
;***edit: again, assume no radial profile subtraction for now.

if (~keyword_set(suffname) and ~keyword_set(postsdi)) then begin
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

if (keyword_set(postsdi) and ~keyword_set(suffname)) then suffname='_sdialocisub'

;*****

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
filesout=filelist(filenum,prefix=prefname,suffix='_alocisub')
endif

if keyword_set(fc) then begin
filesout=filelist(filenum,prefix=prefname,suffix='_alocisub_fc')
endif

nlist=indgen(nfiles)

;*************
;print,files
;stop
;*************

;Define region for spider mask, the image FWHM, and the saturation radius
param,'spang*',spang,/get,pfname=pfname
;param,'spmask',spmask,/get,pfname=pfname
param,'fwhm',fwhm,/get,pfname=pfname
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
    if size(na,/type) eq 7 then na=float(strsplit(na,', ',/extract))
    if na[0] eq 0 then na=200.
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
test=readfits(reducdir1+files[0],/exten,h1)
h0=headfits(reducdir1+files[0],ext=0)

;get north-up value
northpa=sxpar(h0,'TOT_ROT',count=northcount)

if northcount eq 0 then begin
print,'no northup value'
print,'something wrong with fits headers'
print,'check them'
goto,endofprogram
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
;print,drsub,na,geom,rmin,rmax

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


dtmean=mean((abs(allpa-shift(allpa,-1)))[0:nfiles-2])*!dtor

;determine radii
if n_elements(drsub0 eq 1) then begin
    nrsub=ceil((rmax-rmin)/drsub0)
    rsub=findgen(nrsub)*drsub0+rmin
    drsub=replicate(drsub0,nrsub)
;above line is why code sometimes crashes in loop!
endif else begin
    nrsub=0
    r=rmin
    rsub=fltarr(1000)
    drsub=fltarr(1000)
    while r lt rmax do begin
        dr=((0.5+atan((r-drsub0[2])/drsub0[3])/!pi)*(drsub0[1]-drsub0[0])+drsub0[0])
        rsub[nrsub]=r
        drsub[nrsub]=dr
        r+=dr
        nrsub+=1
    endwhile
    rsub=rsub[0:nrsub-1] & drsub=drsub[0:nrsub-1]
endelse
drsub=drsub<(rmax-rsub)

;array of distances and angles to determine indices in each section
distarr=shift(dist(dim),dim/2,dim/2)
ang=(angarr(dim)+2.*!pi) mod (2.*!pi)

;Wavelength Loop, ADI per Wavelength
;we want ADI for datacubes, i.e. several spectral channels but also for
      ;other type of data: collapsed datacubes, single spectral channel ADI, ADI after SDI,etc...
      ; so we have to verify the dimension of ADI inputs hereafter:

for il=0,n_elements(lambda)-1 do begin

if keyword_set(channel) then begin
if il ne long(channel) then goto,skipthischannel
endif


;use coeff
if keyword_set(usecoeff) then begin
if (il_use eq !null) eq 0 then begin
coeffstouse=where(il_use eq il,ncoeffstouse)
if ncoeffstouse gt 0 then begin
c_use2=ck_use[coeffstouse]
ir_use2=ir_use[coeffstouse]
it_use2=it_use[coeffstouse]
nf_use2=nf_use[coeffstouse]
endif
endif
endif


print,'LOCI Wavelength '+strtrim(il+1,2)+'/'+strtrim(n_elements(lambda),2)

;I think the GPI pipeline has this wrong.  Redo.
;Put this outside of the loop so you save time.
fwhm=1.0*(1.d-6*lambda[il]/Dtel)*(180.*3600./!dpi)/pixscale

print,fwhm,lambda[il],Dtel,pixscale

;estimates the largest optimization radius needed
rimmax=0.
for ir=0,nrsub-1 do begin
    r=rsub[ir]
    if n_elements(na) eq 1 then area=na*!pi*(fwhm/2.)^2 $
      else area=((0.5+atan((r-na[2])/na[3])/!pi)*(na[1]-na[0])+na[0])*!pi*(fwhm/2.)^2
    ;width of optimization radius desired
    dropt=sqrt(geom*area)
    nt=round((2*!pi*(r+dropt/2.)*dropt)/area)>1
    dropt=sqrt(r^2+(nt*area)/!pi)-r
    rimmax>=r+dropt
endfor
rimmax<=1.1*dim/2 ;for CHARIS

;cut the image of the rings 5 pixels wide
;save in a file

drim=5.

nrim=ceil((rimmax-rmin)/drim)
rim=findgen(nrim)*drim+rmin
print,'stuff',nrim,drim,rmin
;determine indices of pixels included in each ring
;DRIM of pixels and save them to disk

for ir=0,nrim-1 do begin
    ri=rim[ir] & rf=ri+drim
    ia=where(distarr lt rf and distarr ge ri)
    openw,funit,tmpdir+'indices_a'+nbr2txt(ir,3)+'.dat',/get_lun
    writeu,funit,ia
    free_lun,funit
endfor

;Cutting images into rings and place rings even nfile radius 
;in a single file

el=dblarr(nfiles) & az=dblarr(nfiles)
dec=dblarr(nfiles) & decdeg=dblarr(nfiles)
dtpose=dblarr(nfiles)
noise_im=fltarr(nrim,nfiles)

for nf=0,nfiles-1 do begin
    h0=headfits(reducdir1+files[nf],/silent,ext=0)
    im=(readfits(reducdir1+files[nf],h1,/exten,/silent))[*,*,il]
    if keyword_set(fwdmod) then begin
    im_fwd=(readfits(reducdir1+filesfwd[nf],/exten,horig,/silent))[*,*,il]
    endif

    ;print,'reading in file 'files[nf]

    norm=0
    if norm then begin
        ;normalize image by radial profile noise
        profrad,abs(im),2.,0.,rimmax,p2d=pr
        im/=pr
        writefits,tmpdir+filesprof[nf],pr
    endif

    for ir=0,nrim-1 do begin
        ia=read_binary(tmpdir+'indices_a'+nbr2txt(ir,3)+'.dat',data_type=3)
        if ia[0] eq -1 then continue
        openw,funit,tmpdir+'values_a'+nbr2txt(ir,3)+'.dat',/get_lun,append=(nf gt 0)
        writeu,funit,float(im[ia])
        free_lun,funit
        if keyword_set(fwdmod) then begin
        openw,funit,tmpdir+'values_fwd'+nbr2txt(ir,3)+'.dat',/get_lun,append=(nf gt 0)
        writeu,funit,float(im_fwd[ia])
        free_lun,funit
        endif

        ;calcule le bruit dans cet anneau
        noise_im[ir,nf]=median(abs(im[ia]-median(im[ia])))/0.6745
    endfor

;**note: the GPI pipeline just reads in 'DEC' from the keyword.  Here, we read it from the fits header.

    dec[nf]=sxpar_charis(h0,'DEC',/justfirst,/silent)*!dtor
    decdeg[nf]=dec[nf]*!radeg
    dtpose[nf]=abs(rot_ratef(allha[nf],decdeg[nf],lat))*exptime*coadds*!radeg
endfor

;MAIN LOOP      
;on all rings, determined annulus ref, removed

iaim_loaded=-1
for ir=0,nrsub-1 do begin

;coeffs
if keyword_set(usecoeff) then begin
if (ir_use2 eq !null) eq 0 then begin
coeffstouse=where(ir_use2 eq ir,ncoeffstouse)
if ncoeffstouse gt 0 then begin
c_use3=c_use2[coeffstouse]
it_use3=it_use2[coeffstouse]
nf_use3=nf_use2[coeffstouse]
endif
endif
endif

    ri=rsub[ir] & dr=drsub[ir] & r=ri+dr/2. & rf=ri+dr
    print,'ALOCI Wavelength '+strtrim(il+1,2)+'/'+strtrim(n_elements(lambda),2),' Annulus '+strtrim(ir+1,2)+'/'+strtrim(nrsub,2)+' with radius '+$
        string(r,format='(f5.1)')+$
        ' [>='+string(ri,format='(f5.1)')+', <'+string(rf,format='(f5.1)')+']...'

    if n_elements(na) eq 1 then area=na*!pi*(fwhm/2.)^2 $
      else area=((0.5+atan((r-na[2])/na[3])/!pi)*(na[1]-na[0])+na[0])*!pi*(fwhm/2.)^2

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

    if 1 then begin
        r1opt=ri
        ;number of annulus section
        nt=round((2*!pi*(r1opt+dropt/2.)*dropt)/area)>1
        ;print,'nt is',nt,r1opt,dropt,area

        ;dropt for annulus with sections of exact area

        dropt=sqrt(r1opt^2+(nt*area)/!pi)-r1opt
        r2opt=r1opt+dropt

        if r2opt gt rim[nrim-1]+drim then begin
            r2opt=rim[nrim-1]+drim
            dropt=r2opt-r1opt
            nt=round((2*!pi*(r1opt+dropt/2.)*dropt)/area)>1
        endif
    endif

    ;subtracted for region in central reg_optimization
    if 0 then begin
        ;number of annulus section

        nt=round((2*!pi*r*dropt)/area)>1

         ;and for optimization the center annulus on r
         ;dr_opt for sections with exact area

        dropt=(area*nt)/(2.*!pi*r)
        r1opt=r-dropt/2.
        r2opt=r+dropt/2.
    
        if r1opt lt rmin then begin
            r1opt=rmin
            dropt=sqrt(geom*area)
            nt=round((2*!pi*(r1opt+dropt/2.)*dropt)/area)>1
            r2opt=sqrt((area*nt)/!pi+r1opt^2)
            dropt=r2opt-r1opt
        endif
        if r2opt gt rim[nrim-1]+drim then begin
            r2opt=rim[nrim-1]+drim
            dropt=sqrt(geom*area)
            nt=round((2*!pi*(r2opt-dropt/2.)*dropt)/area)>1
            r1opt=sqrt(r2opt^2-(area*nt)/!pi)
            dropt=r2opt-r1opt
        endif
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
    
    ;*debugging*
    ;    print,'k is ',k,' ir is ',ir,c2load2,iaim_2load,c2load
    ;    print,iaim_2load[k]
    ;    print,'  ','indices_a'+nbr2txt(iaim_2load[k],3)+'.dat'

        ia_tmp=read_binary(tmpdir+'indices_a'+nbr2txt(iaim_2load[k],3)+'.dat',data_type=3)
        annuli_tmp=read_binary(tmpdir+'values_a'+nbr2txt(iaim_2load[k],3)+'.dat',data_type=4)
 
        if keyword_set(fwdmod) then $
        annuli_fwd_tmp=read_binary(tmpdir+'values_fwd'+nbr2txt(iaim_2load[k],3)+'.dat',data_type=4)
        annuli_tmp=reform(annuli_tmp,n_elements(ia_tmp),nfiles)
        if keyword_set(fwdmod) then $
        annuli_fwd_tmp=reform(annuli_fwd_tmp,n_elements(ia_tmp)*1L,nfiles*1L)

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
    ;print,'hihihihi',r1opt,r2opt,ir,nrim-1
    ;if (ir eq nrim-1) then print,'hihi',r1opt,r2opt,ir,nrim-1 
    ;if (ir eq nrim-1) then stop
    
     if (pixmask gt 0) then iaopt2=where(distarr[ia] ge r1opt +dr and distarr[ia] lt r2opt)

    ;angle of annulus sections
    dt=2.*!pi/nt
    ;print,'sqrt ',sqrt(area*geom),'r1opt is ',r1opt,' r2opt is',r2opt,' dropt is',r2opt-r1opt

    ;loop on angular sections

    for it=0,nt-1 do begin
        ;indices of pixels included in this section: i.e. the optmization region

    if keyword_set(usecoeff) then begin
    if (it_use3 eq !null) eq 0 then begin
    coeffstouse=where(it_use3 eq it,ncoeffstouse)
    if ncoeffstouse gt 0 then begin
    c_use4=c_use3[coeffstouse]
    nf_use4=nf_use3[coeffstouse]
    endif
    endif
    endif

        iopt=where(ang[ia[iaopt]] ge it*dt and ang[ia[iaopt]] lt (it+1)*dt) 
        if (pixmask gt 0) then $
        iopt2=where(ang[ia[iaopt2]] ge it*dt and ang[ia[iaopt2]] lt (it+1)*dt)

        npix=n_elements(iopt)

        if npix lt 5 then continue
 
        iopt=iaopt[iopt]
         
        ;instructs the region of optimization in memory

        optreg=annuli[iopt,*]
        if keyword_set(fwdmod) then begin
          optreg_fwd=annuli_fwd[iopt,*]
          optreg+=annuli_fwd[iopt,*]
        endif
        ;indices of pixels to subtract

        isub=where(distarr[ia[iopt]] ge ri and distarr[ia[iopt]] lt rf)
        if n_elements(isub) lt 2 then continue
        isub=iopt[isub]

        if (pixmask gt 0) then begin
        ;if keyword_set(pixmask) then begin
         iopt=iaopt2[iopt2]
         optreg=annuli[iopt,*]
         if keyword_set(fwdmod) then begin 
           optreg_fwd=annuli_fwd[iopt,*]
           optreg+=annuli_fwd[iopt,*]
         endif
        endif

;        ****masking deviants****

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
        z/=(replicate(1,n_elements(iopt))#median(abs(z),dim=1,/even))
        z/=(median(abs(z),dim=2,/even)#replicate(1,nfiles))
        z=(abs(z) lt 7 and finite(z) eq 1)
        ; z=(finite(z) eq 1)

        zq=median(z,dimension=2,/even)

        ;for a given pixel [i,*] look to see whether an image pixel that is deviant
;a patch for now: SVD screws up for NaNs/zeroes in outer image slice regions when using post-SDI
if ~keyword_set(postsdi) then begin
        igoodf=where(min(z,dim=2),cgood)
endif else begin
        igoodf=where(zq ne 0 and finite(zq) ne 0,cgood)
        optreg=optreg[igoodf,*]
        iopt=iopt[igoodf,*]
        if keyword_set(fwdmod) then optreg_fwd=optreg_fwd[igoodf,*]
endelse

        ;if cgood lt 5 then continue
        ;optreg=optreg[igoodf,*]
        ;iopt=iopt[igoodf]
        ;if keyword_set(fwdmod) then optreg_fwd=optreg_fwd[igoodf,*]

        ;there clues to avoid a build images later
        
        openw,lunit,tmpdir+'indices_images.dat',/get_lun,append=(ir+it gt 0)
        writeu,lunit,ia[isub]
        free_lun,lunit
        ;build large matrix of a linear system to solve
       
       if ~keyword_set(fwdmod) then begin   
       aa=optreg##transpose(optreg)
       endif else begin
       optreg=optreg
       ;optreg=optreg+optreg_fwd
       aa=optreg##transpose(optreg)
       endelse


        ;loop on all images and made the last
        for nf=0,nfiles-1 do begin
            ;separation angulaire de toutes les images par rapport a image n
            ; angular separation of all pictures from a photo n

            dpa=abs(allpa-allpa[nf])

            ;OFFSET determined enough images for subtraction

            indim=where(dpa gt (nfwhm*fwhm/ri*!radeg+dtpose[nf]),c1)

;correlation-based frame selection

            if nref le c1 then begin
            index_corr=(reverse(sort(aa[indim,nf])))[0:nref-1 < c1]
            indim=indim[index_corr]
            c1=n_elements(indim)
            endif


            ;print,'nelements',c1,nfwhm,fwhm/ri*!radeg,max(dpa),dtpose[nf]
            igood=where(finite(annuli[isub,nf]) eq 1,c2)
            ;if c1 eq 0 then begin
            if c1 eq 0 or c2 lt 2 then begin
                ;*debug*
                ;diff=fltarr(n_elements(isub))*0
                diff=fltarr(n_elements(isub))+!values.f_nan
                difffwd=diff
            endif else begin
;***edit: For now, let's not do damped-LOCI.  Leave that for a future iteration

;now include option to use saved coefficients.
if keyword_set(usecoeff) then begin
;coeffstouse =where(il_use eq il and ir_use eq ir and it_use eq it and nf_use eq nf)
if (nf_use4 eq !null) eq 0 then begin
coeffstouse=where(nf_use4 eq nf,ncoeffstouse)

if ncoeffstouse gt 0 then begin
c=c_use4[coeffstouse]
endif
endif

;if ~keyword_set(fwdmod) then goto,skipmatrixinversion


endif
                ;if ~keyword_set(fwdmod) then begin
                ;matrix of linear system to solve
                a=(aa[indim,*])[*,indim]
                ;vector b of a linear system to solve
                b=aa[indim,nf]
                
                ;solve the system
                ;c=invert(a,/double)#b

                if (keyword_set(svd) and n_elements(a) gt 2) then begin
                inv_a=svd_invert(a,svd,/double)
                c=inv_a#b
                endif else begin
                c=invert(a,/double)#b
                endelse

                if ~keyword_set(fwdmod) then begin
                ;construct the reference
                skipmatrixinversion:

                ref=fltarr(n_elements(isub))
                for k=0,c1-1 do begin 
                 ref[igood]+=c[k]*annuli[isub[igood],indim[k]]
                 if keyword_set(savecoeff) then printf,1,long(il),long(ir),long(it),long(nf),c[k]
                endfor

                ;make the difference
        ;        reftot=total(ref) 
                
                ;*debug*
                ;if(finite(reftot) eq 0)then begin
                ;diff=0*(annuli[isub,nf]-ref)
                ;goto,skipme
                ;endif

                diff=annuli[isub,nf]-ref
    
               if (zero gt 0) then diff-=median(diff,/even)
                skipme:

;fwdmod
               endif else begin 

               ;if you are doing forward modeling then ...
               ;   - you already have the set of coefficients c_k
               ;   - you need to solve for the perturbing coefficients beta

               ref=fltarr(n_elements(isub))
               refpert=fltarr(n_elements(iopt))
               ;self-subtraction
               for k=0,c1-1 do begin
                ;ref[igood]+=c[k]*annuli_fwd[isub[igood],indim[k]]
                refpert[igood]+=c[k]*annuli_fwd[iopt[igood],indim[k]]
               endfor

;just self-subtraction
               ;diff1=annuli_fwd[isub,nf]-refpert

;the first part of the perturbed column vector b'
               pert1=(annuli[iopt,*])[*,indim]
               ;pert1+=(annuli_fwd[iopt,*])[*,indim]
;the second part of the perturbed column vector b'
               pert2=annuli_fwd[iopt,nf]-refpert

;now based on the above self-subtraction, introduce as a perturbation on the linear system
               ;previous A covariance matrix is the same as the original
               ;vector A of linear system to solve
               a=(aa[indim,*])[*,indim]

               ;for every l in column, b,   b_l= sum(pixels)_ annuli_orig(pixels,set of ref images) mult by self-subtracted planet signal(pixels).
               ;part1##(transpose(part2) t    o mult and sum over pixels
               b=pert2##transpose(pert1)

               if (keyword_set(svd) and n_elements(a) gt 2) then begin
                inv_a=svd_invert(a,svd,/double)
                pert_c=inv_a#b
                endif else begin
                pert_c=invert(a,/double)#b
               endelse

               ref2=fltarr(n_elements(isub))

               for k=0,c1-1 do begin
                 if (finite(pert_c[k]) ne 0) then begin

                 ;ref2[igood]+=pert_c[k]*annuli[isub[igood],indim[k]]
                 ;ref[igood]+=c[k]*annuli_fwd[isub[igood],indim[k]]+pert_c[k]*annuli[isub[igood],indim[k]]

                 ref[igood]+=c[k]*annuli_fwd[isub[igood],indim[k]]+pert_c[k]*annuli[isub[igood],indim[k]]+pert_c[k]*annuli_fwd[isub[igood],indim[k]]

                 endif
               endfor
               diff=annuli_fwd[isub,nf]-ref
               if (zero gt 0) then diff-=median(diff,/even)
               ;diff=annuli_fwd[isub,nf]
             endelse
            endelse
        
            ;register the difference, add (append) the values of this annulus to 
            ;binary file image
            openw,lunit,tmpdir+filestmp[nf],/get_lun,append=(ir+it gt 0)
            writeu,lunit,diff
            free_lun,lunit
        endfor
    endfor
endfor

;deletes files. .dat annuli
;file_delete,file_search(tmpdir,'indices_a*.dat')
;file_delete,file_search(tmpdir,'values_a*.dat')

;reading signs of pixels removed

ind=read_binary(tmpdir+'indices_images.dat',data_type=3)

;delete temporary indices
file_delete,tmpdir+'indices_images.dat'

;rebuild and turn images


for nf=0,nfiles-1 do begin
    print,'Image '+strtrim(nf+1,2)+'/'+strtrim(nfiles,2)+': '+files[nf]+'...'

    ;reconstruct image
    print,' reconstruction...'

    im=make_array(dim,dim,type=4,value=!values.f_nan)
    im[ind]=read_binary(tmpdir+filestmp[nf],data_type=4)
    ;delete temporary files
    file_delete,tmpdir+filestmp[nf]

    if norm eq 1 then begin

        ;multiplied by the radial profile of noise
        pr=readfits(tmpdir+filesprof[nf],/silent)
        im*=pr
        file_delete,tmpdir+filesprof[nf]
    endif

    ;get header
    h1=headfits(reducdir1+files[nf],/exten,/silent)
    h0=headfits(reducdir1+files[nf],ext=0,/silent)

; ;***removed a lot of legacy code here useful for Gemini/NIRI but otherwise not.
;;**** see older code versions
   
    ;rotate to bring the first image

    print,' rotation...'


;******astrometry*****

;****** default is to rotate north-up given fits header information and the north PA offset

;throw the /norot switch if you do not want to derotate the image to a common value (e.g., if you are doing SDI on residuals after this step)

if ~keyword_set(norot) then begin
theta=-1*allpa[nf]

charis_northup,im,h0,h1,angoffset=angoffset
endif else begin

endelse

    ;rotskip:


    ;register the difference

 ;register the difference
    if keyword_set(norot) then begin

    endif else begin
    sxaddhist,'A rotation of '+string(theta,format='(f8.3)')+$
      ' degrees was applied to the image.',h1
    endelse
; Add ALOCI keywords
    sxaddpar,h1,'rsub',rsubval
    sxaddpar,h1,'loci_nfwhm ',loci_nfwhm
    sxaddpar,h1,'loci_na ',loci_na
    sxaddpar,h1,'loci_geom ',loci_geom
    sxaddpar,h1,'loci_drsub ',loci_drsub
    sxaddpar,h1,'loci_rmin ',loci_rmin
    sxaddpar,h1,'loci_rmax ',loci_rmax
    sxaddpar,h1,'loci_nref ',loci_nref
    sxaddpar,h1,'svd ',svd
    
    sxaddpar,h1,'pixmaskv',pixmask
    sxaddpar,h1,'zero',zero
    sxaddpar,h1,'meanadd',meanadd
    suffix1='-loci'+strcompress(string(il),/REMOVE_ALL)

    writefits,tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix1+'.fits',0,h0
    writefits,tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix1+'.fits',im,h1,/append
endfor


skipthischannel:
endfor ;wav

breakoutwavelengthloop:

;endtime=systime(/seconds)
;Okay, now combine the images together, construct datacubes, and then construct a combined datacube.
imt=dblarr(dim,dim,n_elements(lambda))
;imtot=dlbarr(dim,dim,n_elements(lambda),nfiles)
suffix0='-loci'

for nf=0,nfiles-1 do begin
 for il=0,n_elements(lambda)-1 do begin

 if ~keyword_set(channel) then begin
 h0=headfits(tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix0+strcompress(string(il),/REMOVE_ALL)+'.fits',/SILENT,ext=0)
 imt[*,*,il]=readfits(tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix0+strcompress(string(il),/REMOVE_ALL)+'.fits',/SILENT,/exten,h1)
 endif else begin
 h0=headfits(tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix0+strcompress(string(channel),/REMOVE_ALL)+'.fits',/SILENT,ext=0)
 imt[*,*,il]=readfits(tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix0+strcompress(string(channel),/REMOVE_ALL)+'.fits',/SILENT,/exten,h1)
 endelse

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

file_delete,file_search(tmpdir,'*-loci*fits')

if (keyword_set(fwdmod) or keyword_set(savecoeff) or keyword_set(usecoeff)) then free_lun,1

;close,/all

print,'ANGOFFSET WAS',angoffset

endtime=systime(/seconds)
print,"elapsed time is ",endtime-starttime," seconds"
endofprogram:

end
