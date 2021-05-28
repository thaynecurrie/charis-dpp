pro charis_rdiklip,pfname,refpath=refpath,prefname=prefname,suffname=suffname,nfwhm=nfwhm,drsub=drsub0,$
npca=npca,meansub=meansub,$
rsubval=rsubval,$
meanadd=meanadd,$
zero=zero,$
fwdmod=fwdmod,$
savecoeff=savecoeff,$
usecoeff=usecoeff,$
rmin=rmin,rmax=rmax,fc=fc,$
norot=norot,nonorthup=nonorthup,angoffset=angoffset,$
outfile=outfile,$
 guide=guide,help=help

;***10/15/2019** - version 2 of PCA/KLIP-RDI - still does not do zonal (local) subtraction, but annular sub

;***12/13/2018** - version 1 of PCA/KLIP-RDI code - does reference star differential imaging using a PSF library
;***01/30/2018**
;Version 1.2 of PCA/smart-KLIP for CHARIS
;-added KLIP forward-modeling capabilities a la Pueyo (2016).   Cleaned up code.


if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"KLIP-RDI PSF subtraction, Version 2.0 (January 2019)"
print,""
print,"charis_rdiklip,pfname,npca=npca,refpath=refpath,prefname=prefname,suffname=suffname,nfwhm=nfwhm,drsub=drsub0,meansub=meansub"
print,"rsubval=rsubval,zero=zero,fwdmod=fwdmod,savecoeff=savecoeff,usecoeff=usecoeff"
print,"rmin=rmin,rmax=rmax,fc=fc"
print,"norot=norot,nonorthup=nonorthup,angoffset=angoffset,outfile=outfile"
print,''
print,'Example:'
print,"charis_rdiklip,'LkCa15_low.info',refpath='../V819Tau_low/',nfwhm=0.75,npca=10,rmin=2,rmax=70,drsub=68,outfile='rdiklip.fits'"
print,''
print,"***Important Keywords***"
print,'*pfname - parameter file (e.g. LkCa15_low.info)'
print,"*refpath - (local) path of reference library folder"
print,"*npca - number of KL modes retained"
print,"*rmin - minimum radius of subtraction"
print,"*rmax - maximum radius of subtraction"
print,"*drsub - radial width of subtraction region"
print,"*meansub - do you subtract off the mean first? [usually yes]"
print,"*rsubval - do you PSF-subtract spatially-filtered images instead of non-filtered ones? [usually no]"
print,"*zero - do you subtract off the median of a region after PSF-subtracting?"
print,"*savecoeff - save the eigenvectors/eigenvalues in a file [for forward-modeling]"
print,"*usecoeff - use eigenvectors/values [in forward-modeling]"
print,"*fwdmod - forward-modeling of a synthetic source [called in charis_klip_fwdmod_planet/disk]"
print,"*outfile - the name of the output file"

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

if keyword_set(savecoeff) then openw,1,'rdiklipcoeff.dat'
if keyword_set(usecoeff) then readcol,'rdiklipcoeff.dat',il_use,ir_use,nf_use,eigenvect_use,ck_use,format='i,i,i,f20.15,d20.10'

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

if ~keyword_set(suffname) then begin
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

filenum=nbrlist(flist)
;files=filelist(filenum,nfiles,prefix=prefname,suffix='_f',/gz)

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
filesout=filelist(filenum,prefix=prefname,suffix='_adisub')
endif

if keyword_set(fc) then begin
filesout=filelist(filenum,prefix=prefname,suffix='_adisub_fc')
endif

nlist=indgen(nfiles)

;*************


;********Reference Library*****
if keyword_set(refpath) then begin
if ~keyword_set(rsubval) then begin
;reffilenames=file_search(refpath+'/reduc/reg/'+'*.fits')
reffilenames=file_search(refpath+'/reduc/reg/'+'*'+suffname+'.fits')
endif else begin
;reffilenames=file_search(refpath+'/reduc/rsub/'+'*.fits')
reffilenames=file_search(refpath+'/reduc/rsub/'+'*'+suffname+'.fits')
;print,reffilenames
endelse
;reffilenames=file_search(refpath+'*.fits')

endif else begin
reffilenames=dialog_pickfile(title="Select Reference PSF library Files",/multiple)
endelse

nreffiles=n_elements(reffilenames)
;print,nreffiles
;print,reffilenames
;print,refpath
;*************

;Define region for spider mask, the image FWHM, and the saturation radius
param,'spang*',spang,/get,pfname=pfname
param,'spmask',spmask,/get,pfname=pfname
param,'fwhm',fwhm,/get,pfname=pfname
param,'rsat',rsat,/get,pfname=pfname

;ADI algorithm parameters

;**get dim from first image header
if ~keyword_set(nfwhm) then begin
    param,'nfwhm',nfwhm,/get,pfname=pfname
    if nfwhm eq 0 then nfwhm=0.75
endif
if ~keyword_set(drsub0) then begin
    param,'dr',drsub0,/get,pfname=pfname
    ;if size(drsub0,/type) eq 7 then drsub0=float(strsplit(drsub0,', ',/extract))
    ;if drsub0[0] eq 0 then drsub0=5.
    ;if drsub0[0] eq 0 then drsub0=3.
endif
if ~keyword_set(rmin) then begin
     param,'rmin',rmin,/get,pfname=pfname
endif
if ~keyword_set(rmax) then begin
     param,'rmax',rmax,/get,pfname=pfname
endif
if ~keyword_set(npca) then begin
    param,'KLIPNPCA',npca,/get,pfname=pfname
endif

if ~keyword_set(nzone) then begin
    param,'KLIPNZON',nzone,/get,pfname=pfname
endif

if ~keyword_set(meansub) then begin
    param,'meansub',meansub,/get,pfname=pfname
endif

if ~keyword_set(meanadd) then begin
    param,'meanadd',meanadd,/get,pfname=pfname
endif

if ~keyword_set(zero) then begin
   param,'zero',zero,/get,pfname=pfname
endif

;vestigal, remove later
if ~keyword_set(na) then begin
    param,'na',na,/get,pfname=pfname
    if size(na,/type) eq 7 then na=float(strsplit(na,', ',/extract))
    if na[0] eq 0 then na=200.
endif
print,'files',files[0]
;h=headfits(reducdir1+files[0])
test=readfits(reducdir1+files[0],/exten,h1)
h0=headfits(reducdir1+files[0],ext=0)
;test=readfits(reducdir1+files[0],/exten)
;get north-up value
northpa=sxpar(h0,'TOT_ROT',count=northcount)

dim=sxpar(h1,'naxis1')
xc=dim/2 & yc=dim/2
print,dim,xc,yc
print,size(test)

;Now Get the Wavelength Vector
filter=sxpar(h0,'FILTNAME',count=filtcount)
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

;if ~keyword_set(rmin) then rmin=(rsat-5)>2
if ~keyword_set(rmax) then rmax=1.1*dim/2

;**For fits header keywords later
;**For fits header keywords later
adi_nfwhm=nfwhm
adi_drsub=drsub0
adi_npca=npca
adi_rmin=rmin
adi_rmax=rmax

;load parallactic angles
readcol,'reduc.log',ffilenum,allxc,allyc,allrsat,allha,allpa,/silent

;Debugging
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
;we want ADI for datacubes, i.e. several specral channels but also for
      ;other type of data: collapsed datacubes, single spectral channel ADI, ADI after SDI,etc...
      ; so we have to verify the dimension of ADI inputs hereafter:

for il=0,n_elements(lambda)-1 do begin
;use coeff
if keyword_set(usecoeff) then begin
coeffstouse=where(il_use eq il,ncoeffstouse)

if ncoeffstouse gt 0 then begin
c_use2=ck_use[coeffstouse]
;il_use2=il_use[coeffstouse]
eigenvect_use2=eigenvect_use[coeffstouse]
ir_use2=ir_use[coeffstouse]
nf_use2=nf_use[coeffstouse]
endif

endif


print,'PCA/KLIP Wavelength '+strtrim(il+1,2)+'/'+strtrim(n_elements(lambda),2)

;I think the GPI pipeline has this wrong.  Redo.
;Put this outside of the loop so you save time.
fwhm=1.*(1.d-6*lambda[il]/Dtel)*(180.*3600./!dpi)/pixscale

print,fwhm,lambda[il],Dtel,pixscale

;cut the image of the rings 5 pixels wide
;save in a file

drim=5.

nrim=ceil((rmax-rmin)/drim)
rim=findgen(nrim)*drim+rmin

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
    ;print,'reading in file ',files[nf]

    ;help,im
    ;writefits,'im.fits',im
    ;stop

    for ir=0,nrim-1 do begin
        ia=read_binary(tmpdir+'indices_a'+nbr2txt(ir,3)+'.dat',data_type=3)
        if ia[0] eq -1 then continue
        openw,funit,tmpdir+'values_a'+nbr2txt(ir,3)+'.dat',/get_lun,append=(nf gt 0)
        writeu,funit,im[ia]
        free_lun,funit

        if keyword_set(fwdmod) then begin
        openw,funit,tmpdir+'values_fwd'+nbr2txt(ir,3)+'.dat',/get_lun,append=(nf gt 0)
        writeu,funit,im_fwd[ia]
        free_lun,funit
        endif

        ;calcule le bruit dans cet anneau
        noise_im[ir,nf]=median(abs(im[ia]-median(im[ia])))/0.6745
    endfor
    el[nf]=sxpar(h0,'altitude')*!dtor
     
;**note: the GPI pipeline just reads in 'DEC' from the keyword.  Here, 
; keep original method of using DEC as determined from the elevation and azimuth

    dec[nf]=sxpar_charis(h0,'DEC',/justfirst,/silent)*!dtor
    decdeg[nf]=dec[nf]*!radeg
    dtpose[nf]=abs(rot_ratef(allha[nf],decdeg[nf],lat))*exptime*coadds*!radeg
endfor


for nfref=0L,nreffiles-1 do begin
    im=(readfits(reffilenames[nfref],h,/silent,/exten))[*,*,il]
    
    for ir=0,nrim-1 do begin
        ia=read_binary(tmpdir+'indices_a'+nbr2txt(ir,3)+'.dat',data_type=3)
        if ia[0] eq -1 then continue
;**** the reference data
        openw,funit,tmpdir+'values_aref'+nbr2txt(ir,3)+'.dat',/get_lun,append=(nfref gt 0)
        writeu,funit,float(im[ia])
        free_lun,funit
    endfor
endfor


;MAIN LOOP      
;on all rings, determined annulus ref, removed

iaim_loaded=-1
for ir=0,nrsub-1 do begin

;coeffs
if keyword_set(usecoeff) then begin
coeffstouse=where(ir_use2 eq ir,ncoeffstouse)

if ncoeffstouse gt 0 then begin
c_use3=c_use2[coeffstouse]
eigenvect_use3=eigenvect_use2[coeffstouse]
nf_use3=nf_use2[coeffstouse]
endif

endif

    ri=rsub[ir] & dr=drsub[ir] & r=ri+dr/2. & rf=ri+dr
    print,'PCA Wavelength '+strtrim(il+1,2)+'/'+strtrim(n_elements(lambda),2),' Annulus '+strtrim(ir+1,2)+'/'+strtrim(nrsub,2)+' with radius '+$
        string(r,format='(f5.1)')+$
        ' [>='+string(ri,format='(f5.1)')+', <'+string(rf,format='(f5.1)')+']...'

    ;print,' Annulus '+strtrim(ir+1,2)+'/'+strtrim(nrsub,2)+' at radius '+$
    ;  string(r,format='(f5.1)')+$
    ;  ' [>='+string(ri,format='(f5.1)')+', <'+string(rf,format='(f5.1)')+']...'

    ;aire de region a ce rayon
    ;area of this region has radius

    if n_elements(na) eq 1 then area=na*!pi*(fwhm/2.)^2 $
      else area=((0.5+atan((r-na[2])/na[3])/!pi)*(na[1]-na[0])+na[0])*!pi*(fwhm/2.)^2

    ;largeur de l'anneau d'optimisation desiree
    ;width of desired optimization annulus
    ;dropt=sqrt(geom*area)

    dropt=dr

    ;if dropt lt dr then begin
    ;    print,'dropt < drsub !!!'
    ;    print,'dropt: ',dropt
    ;    print,'drsub: ',dr
    ;    stop
    ;endif

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

    ;if ir gt 0 then begin
   ;     irm=where(distarr[ia] lt rim[i1aim] or distarr[ia] ge rim[i2aim]+drim,crm,complement=ikp)
   ;     if crm gt 0 then remove,irm,ia
   ;     if crm gt 0 then annuli=annuli[ikp,*]
   ;     if keyword_set(fwdmod) then begin
   ;     if crm gt 0 then annuli_fwd=annuli_fwd[ikp,*]
   ;     endif
   ; endif
    
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

        ia_tmp=read_binary(tmpdir+'indices_a'+nbr2txt(iaim_2load[k],3)+'.dat',data_type=3)
        annuli_tmp=read_binary(tmpdir+'values_a'+nbr2txt(iaim_2load[k],3)+'.dat',data_type=4)

        ia_tmpref=ia_tmp
        annuli_tmpref=read_binary(tmpdir+'values_aref'+nbr2txt(iaim_2load[k],3)+'.dat',data_type=4)

        if keyword_set(fwdmod) then $
        annuli_fwd_tmp=read_binary(tmpdir+'values_fwd'+nbr2txt(iaim_2load[k],3)+'.dat',data_type=4)

        annuli_tmp=reform(annuli_tmp,n_elements(ia_tmp),nfiles)
        annuli_tmpref=reform(annuli_tmpref,n_elements(ia_tmpref),nreffiles)

        if keyword_set(fwdmod) then $
        annuli_fwd_tmp=reform(annuli_fwd_tmp,n_elements(ia_tmp)*1L,nfiles*1L)

        if ir+k eq 0 then ia=ia_tmp else ia=[ia,ia_tmp]
        if ir+k eq 0 then annuli=annuli_tmp else annuli=[annuli,annuli_tmp]
        if ir+k eq 0 then iaref=ia_tmpref else iaref=[iaref,ia_tmpref]
        if ir+k eq 0 then annuliref=annuli_tmpref else annuliref=[annuliref,annuli_tmpref]

        if keyword_set(fwdmod) then begin
         if ir+k eq 0 then annuli_fwd=annuli_fwd_tmp else annuli_fwd=[annuli_fwd,annuli_fwd_tmp]
        endif

    endfor

 ;not sure if the below is needed

     ;removes the annuli that aren't necessary
    if ir gt 0 then begin
        irm=where(distarr[ia] lt rim[i1aim] or distarr[ia] ge rim[i2aim]+drim,crm,complement=ikp)
        if crm gt 0 then remove,irm,ia
        if crm gt 0 then remove,irm,iaref
        if crm gt 0 then annuli=annuli[ikp,*]
        if crm gt 0 then annuliref=annuliref[ikp,*]

        if keyword_set(fwdmod) then begin
        if crm gt 0 then annuli_fwd=annuli_fwd[ikp,*]
        endif

    endif


    ia_tmp=0 & annuli_tmp=0 & ia_tmpref=0 & annuli_tmpref=0

    ;remembers the list of annuli changes
    iaim_loaded=iaim


   ;indices of pixels for optimization annulus

;****FOR NOW, assuming you are doing ADI PSF sub in one annulus at a time!

        ;indices of pixels included in this section: i.e. the optmization region

        ;instructs the region of optimization in memory


        ;indices of pixels to subtract

        isub=where(distarr[ia] ge ri and distarr[ia] lt rf)

;region of 'optimization', really just the array of n_pixels_annulus by n_images

        ;there clues to avoid a build images later
        openw,lunit,tmpdir+'indices_images.dat',/get_lun,append=(ir gt 0)
        writeu,lunit,ia[isub]
        free_lun,lunit

;Subtract mean/median from each annulus
;you need to do this to do KLIP properly.  Usually, I do this in the imrsub.pro step
       if (meansub gt 0) then begin
       meanopt=fltarr(nfiles)
       meanoptref=fltarr(nreffiles)
       if keyword_set(fwdmod) then meanopt_fwd=fltarr(nfiles)

       for nf=0L,nfiles-1 do begin
       ;meanopt[nf]=median(annuli[isub,nf])
       meanopt[nf]=mean(annuli[isub,nf],/nan)
       annuli[isub,nf]-=meanopt[nf]
       ;optreg[*,nf]-=mean(optreg[isub,nf],/nan)
 
       if keyword_set(fwdmod) then begin
         meanopt_fwd[nf]=mean(annuli_fwd[isub,nf],/nan)
         annuli_fwd[isub,nf]-=meanopt_fwd[nf]
       endif
       endfor

       for nfref=0,nreffiles-1 do begin
        meanoptref[nfref]=mean(annuliref[isub,nfref],/nan)
        annuliref[isub,nfref]-=meanoptref[nfref]
       endfor


       endif

        ;optreg=annuli[isub,*]
        optreg=annuliref[isub,*]
        optreg_sci=annuli[isub,*]
        if keyword_set(fwdmod) then optreg_fwd=annuli_fwd[isub,*]

 ;the following three lines are equivalent to the following
        ;inoise=floor((distarr(ia[iopt])-rim[0])/drim)#replicate(1,90)
        ;tmp=abs(optreg/noise_im[inoise])
        ;z<=(tmp lt 15.)
;        z<=(abs(optreg/noise_im[floor((distarr(ia[iopt])-rim[0])/drim)#replicate(1,nfiles)]) lt 15.)

        ;other way
        z=optreg_sci
        z/=(replicate(1,n_elements(isub))#median(abs(z),dim=1))
        z/=(median(abs(z),dim=2)#replicate(1,nfiles))
        z=(abs(z) lt 7 and finite(z) eq 1)


        ;for a given pixel [i,*] look to see whether an image pixel that is deviant

        igoodf=where(min(z,dim=2),cgood)
        if cgood lt 1 then continue
        ;optreg=optreg[igoodf,*]
        ;optreg_sci=optreg_sci[igoodf,*]
        ;isub=isub[igoodf]
        ;if keyword_set(fwdmod) then optreg_fwd=optreg_fwd[igoodf,*]

        ;loop on all images and made the last
        for nf=0,nfiles-1 do begin
            ;separation angulaire de toutes les images par rapport a image n
            ; angular separation of all pictures from a photo n

            ;dpa=abs(allpa-allpa[nf])

            ;OFFSET determined enough images for subtraction

            ;indim=where(dpa gt (nfwhm*fwhm/ri*!radeg+dtpose[nf]),c1)
            c1=nreffiles

;comment the igood stuff out for now
            igood=where(finite(annuli[isub,nf]) eq 1,c2)


;;Setting up the PCA

           if c1 lt 2 or c2 lt 5 then begin
                diff=fltarr(n_elements(isub))+!values.f_nan
                ;if keyword_set(fwdmod) then diffwwd=fltarr(n_elements(isub))+!values.f_nan
                difffwd=diff
            endif else begin

           npca = npca < c1
           data_mat=optreg[*,*]
           if keyword_set(fwdmod) then data_mat_fwd=optreg_fwd[*,*]
           ;data_mat=transpose(ref_library)

           covmatrix=matrix_multiply(data_mat,data_mat,/atranspose)
;;***forward-modeling
if keyword_set(usecoeff) then begin
coeffstouse=where(nf_use3 eq nf,ncoeffstouse)

if ncoeffstouse gt 0 then begin
;c=ck_use[coeffstouse]
c=c_use3[coeffstouse]
eigenvect=eigenvect_use3[coeffstouse]

;now, reform the eigenvectors and eigenvalues
;help,eigenvect,c1,npca,data_mat_fwd
eigenvect=reform(eigenvect,c1,npca)
eigenval=(reform(c,c1,npca))[0,*]
endif

goto,skipmatrixdecomp
endif



           eigenval=la_eigenql(covmatrix,eigenvectors=eigenvect,$
            range=[c1-npca,c1-1],/double)

           eigenval=reverse(eigenval)
           if npca eq 1 then eigenvect = transpose (eigenvect) else eigenvect=reverse(eigenvect,2)

           skipmatrixdecomp:
           pc=matrix_multiply(eigenvect,data_mat,/atranspose,/btranspose)
     ;construct the reference

;****forward-modeling, calculate the perturbed KL-modes
           ;comment out for RDI
           ;if keyword_set(fwdmod) then dpc=matrix_multiply(eigenvect,data_mat_fwd,/atranspose,/btranspose)

            ref=fltarr(n_elements(isub))

             if keyword_set(fwdmod) then ref_fwd=fltarr(n_elements(isub))
             if keyword_set(savecoeff) then begin

             for k=0L,npca-1 do begin
             for jk=0,c1-1 do begin

             if (npca gt 1) then begin
             printf,1,long(il),long(ir),long(nf),eigenvect[jk,k],eigenval[k],jk,k,format='(i3,i3,i3,f20.10,1x,e20.10,1x,i3,i3)'
             endif else begin
             printf,1,long(il),long(ir),long(nf),(transpose(eigenvect))[jk,k],eigenval[k],jk,k,format='(i3,i3,i3,f20.10,1x,e20.10,1x,i3,i3)'
             endelse
;diagnostic tests             ;print,eigenvect[jk,k],eigenval[k],coeff*reform(pc[k,jk]),npca,c1
             endfor
             endfor
             endif

            
           
            for k=0L,npca-1 do begin

             pc[k,*]=pc[k,*]/sqrt(eigenval[k])
             coeff=total(reform(optreg_sci[*,nf]*reform(pc[k,*])))

;diagnostic tests
             ;ref=ref+coeff*reform(pc[k,*])

             ref=ref+coeff*reform(pc[k,*])
             ;ref[igoodf]=ref[igoodf]+coeff*reform(pc[k,*])

              ;print,coeff

             if keyword_set(fwdmod) then begin
             ;*****forward-model****
             ;perturbed KL models,k: not needed for RDI
             ;dpc[k,*]=dpc[k,*]/sqrt(eigenval[k])
             coeff_oversub=total(reform(optreg_fwd[*,nf]*reform(pc[k,*]))) ;<S|KL>
;             coeff_selfsub1=total(reform(optreg[*,nf]*reform(dpc[k,*])))   ;<N|DKL>KL
;             coeff_selfsub2=total(reform(optreg[*,nf]*reform(pc[k,*])))    ;<N|KL>DKL
             ;coeff_selfsub1=total(reform(optreg_fwd[*,nf]*reform(dpc[k,*])))   ;<S|DKL>KL
             ;coeff_selfsub2=total(reform(optreg_fwd[*,nf]*reform(pc[k,*])))    ;<S|KL>DKL
             ref_fwd=ref_fwd+coeff_oversub*reform(pc[k,*])

             endif
            endfor

            reftot=0.

                diff=annuli[isub,nf]-ref
                ;diff=annuli[isub,nf]
                ;diff=ref

;if you are forward-modeling, then output the annealed section of the synthetic signal, not the subtracted image section
                if keyword_set(fwdmod) then difffwd=annuli_fwd[isub,nf]-ref_fwd
                ;difffwd=annuli_fwd[isub,nf]

                ;diff=annuli[isub,nf]
                diff=float(diff)
                if keyword_set(fwdmod) then difffwd=float(difffwd)

                if (zero gt 0) then diff-=median(diff,/even)
                if (zero gt 0) and keyword_set(fwdmod) then difffwd=float(difffwd - median(difffwd,/even))

                ;if keyword_set(fwdmod) then difffwd=float(ref_fwd)
;test images - debugging
goto,skipme
if keyword_set(verbose) then begin
                testimage=fltarr(dim,dim)
                testimage[ia[isub]]=annuli[isub,nf]

                writefits,'test.fits',testimage
                testimage2=testimage*0
                testimage2[ia[isub]]=ref
                writefits,'test2.fits',testimage2

                ;testc=testimage*0
                ;testc[ia[isub]]=2
                ;writefits,'testc.fits',testc

                ;test
                ;diff=annuli[isub,nf]
endif
                skipme:

            endelse


            ;register the difference, add (append) the values of this annulus to 
            ;binary file image
            openw,lunit,tmpdir+filestmp[nf],/get_lun,append=(ir gt 0)
            if ~keyword_set(fwdmod) then begin
            writeu,lunit,diff
            endif else begin
            ;writeu,lunit,diff-difffwd
            writeu,lunit,difffwd
            ;writeu,lunit,diff
            endelse
            free_lun,lunit
        ;endfor
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

    ;diagnostic ;bah=read_binary(tmpdir+filestmp[nf],data_type=4)

    im[ind]=read_binary(tmpdir+filestmp[nf],data_type=4)

    ;delete temporary files
    file_delete,tmpdir+filestmp[nf]

    ;get header
    h1=headfits(reducdir1+files[nf],/exten,/silent)
    h0=headfits(reducdir1+files[nf],ext=0,/silent)
    ;h=headfits(reducdir+files[nf])

; ;***removed a lot of legacy code here useful for Gemini/NIRI but otherwise not.
;;**** see older code versions   

    ;rotate to bring the first image

    print,' rotation...'
    ;theta=-(allpa[nf]-allpa[0])+angoffset
    ;theta=-(allpa[nf]-allpa[0])

    ;if nf ne 0 then begin 

    ;im=rotat(im,theta,missing=!values.f_nan)
    ;endif

;******astrometry*****

;******default is to rotate north-up using the tot-rot fits header keyword in the primary header
;********add option to ...
;******** - add an additional north PA offset (useful for future precise astro cals)
;******** - completely ignore the north PA rotation (potentially useful for backwards compatibility or a quick reduction)

;****if nothing then just rotate north-up*****
;throw this switch if you do not want to derotate the image to a common value (e.g., if you are doing SDI after this step)

if ~keyword_set(norot) then begin
;if ~keyword_set(nonorthup) then begin
;northpa= sxpar(h0,'TOT_ROT',count=northcount)

theta=-(allpa[nf]-allpa[0])

charis_northup,im,h0,h1,angoffset=angoffset
endif else begin

;theta=-(allpa[nf]-allpa[0])

endelse
    ;rotskip:

    ;if nf ne 0 then im=rotat(im,theta,missing=!values.f_nan,hdr=h)

;register the difference
    if keyword_set(norot) then begin
    ;if keyword_set(nonorthup) then begin
    ;sxaddhist,'A rotation of '+string(theta,format='(f8.3)')+$
    ;  ' degrees was applied to the image.',h1
    endif else begin
    sxaddhist,'A rotation of '+string(theta,format='(f8.3)')+$
      ' degrees was applied to the image.',h1
    endelse

; Add PCA/KLIP keywords
    sxaddpar,h1,'rsub',rsubval
    sxaddpar,h1,'klip_pca ',adi_npca
    sxaddpar,h1,'klip_nfwhm ',adi_nfwhm
    sxaddpar,h1,'klip_drsub ',adi_drsub
    sxaddpar,h1,'klip_rmin ',adi_rmin
    sxaddpar,h1,'klip_rmax ',adi_rmax
    sxaddpar,h1,'meansub',meansub
    sxaddpar,h1,'zero',zero
    sxaddpar,h1,'meanadd',meanadd

    suffix1='-klip'+strcompress(string(il),/REMOVE_ALL)
    ;fname=tmpdir+prefix+nbr2txt(nlist[nf],4)+suffix1+'.fits'

    writefits,tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix1+'.fits',0,h0
    writefits,tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix1+'.fits',im,h1,/append

;,/compress
;    writefits,fname,subarr(im,400),h,/compress
endfor
endfor ;wav

;Okay, now combine the images together, construct datacubes, and then construct a combined datacube.

imt=dblarr(dim,dim,n_elements(lambda))
;imtot=dlbarr(dim,dim,n_elements(lambda),nfiles)
suffix0='-klip'

for nf=0,nfiles-1 do begin

 for il=0,n_elements(lambda)-1 do begin
 h0=headfits(tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix0+strcompress(string(il),/REMOVE_ALL)+'.fits',/SILENT,ext=0)
 imt[*,*,il]=readfits(tmpdir+outfile+'_'+nbr2txt(nlist[nf],4)+suffix0+strcompress(string(il),/REMOVE_ALL)+'.fits',/SILENT,/exten,h1)
 endfor

;suffix=suffix0
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
;im_collapsed=median(im[*,*,4:35],/even,dimension=3)
endelse

;rotate northup
;if ~keyword_set(nonorthup) then begin
;charis_northup,im,h0,h1
;charis_northup,im_collapsed,h0_col,h1_col
;endif

writefits,outname,0,h0
writefits,outname,im,h1,/append

writefits,reducdir+outname_col,0,h0
writefits,reducdir+outname_col,im_collapsed,h1,/append

file_delete,file_search(tmpdir,'*-klip*fits')

if (keyword_set(fwdmod) or keyword_set(savecoeff) or keyword_set(usecoeff)) then free_lun,1

;print,'ANGOFFSET WAS',angoffset

endofprogram:
end


