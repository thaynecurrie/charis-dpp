pro charis_adiklip,pfname,prefname=prefname,suffname=suffname,nfwhm=nfwhm,drsub=drsub0,nzone=nzone,$
npca=npca,meansub=meansub,$
rsubval=rsubval,$
meanadd=meanadd,$
zero=zero,$
fwdmod=fwdmod,$
savecoeff=savecoeff,$
usecoeff=usecoeff,$
rmin=rmin,rmax=rmax,fc=fc,$
nonorthup=nonorthup,angoffset=angoffset,$
outfile=outfile,$
verbose=verbose,$
 help=help

;****02/04/2020**
;Version 2.1
;- added indirect self-subtraction terms a la Pueyo (2016).

;***09/26/2019**
;Version 2.0.1 of PCA/smart-KLIP for CHARIS
;- adding zonal geometry, so you can do locally-opt KLIP instead of full annulus

;***01/30/2018**
;Version 1.2 of PCA/smart-KLIP for CHARIS
;-added KLIP forward-modeling capabilities a la Pueyo (2016).   Cleaned up code.

;***10/16/2017***
;Version 1.1 for CHARIS
;-cleaned up code syntax, removing legacy LOCI keywords, improve help function.
;- to do: implement coefficient saving/importing for forward-modeling.
;         implement switch to use fake companions for contrast curves (fc)

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,'charis_adiklip.pro: PCA/KLIP PSF subtraction method, Version 1.2 (January 2018)'
print,'Written by T. Currie (2014), adapted for CHARIS IFS Data (4/2017)'
print,''
print,'charis_adiklip,pfname,nfwhm=nfwhm,drsub=drsub0,npca=npca,nzone=nzone,meansub=meansub,rsubval=rsubval,'
print,'meanadd=meanadd,zero=zero,fwdmod=fwdmod,savecoeff=savecoeff,usecoff=usecoeff,'
print,'rmin=rmin,rmax=rmax,fc=fc,nonorthup=nonorthup,angoffset=angoffset,outfile=outfile,prefname=prefname,suffname=suffname,help=help'
print,''
print,'Example:'
print,"charis_adiklip,'HR8799l.info',nfwhm=0.75,npca=10,meansub=1,rmin=10,rmax=70,drsub=10,nzone=3,outfile='klip.fits'"
print,""
print,"***Important Keywords***"
print,'*pfname - parameter file (e.g. LkCa15_low.info)'
print,"*nfwhm - rotation gap/exclusion zone (in lambda/D units)"
print,"*npca - number of KL modes retained"
print,"*rmin - minimum radius of subtraction"
print,"*rmax - maximum radius of subtraction"
print,"*drsub - radial width of subtraction region"
print,"*nzone - number of subtraction zones per annulus"
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

if keyword_set(savecoeff) then openw,1,'klipcoeff.dat'
if keyword_set(usecoeff) then readcol,'klipcoeff.dat',il_use,ir_use,it_use,nf_use,eigenvect_use,ck_use,format='i,i,i,i,f20.15,d20.10',/quick

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
;print,files
;stop
;*************

;Define region for spider mask, the image FWHM, and the saturation radius
param,'spang*',spang,/get,pfname=pfname
param,'spmask',spmask,/get,pfname=pfname
param,'fwhm',fwhm,/get,pfname=pfname
param,'rsat',rsat,/get,pfname=pfname

;ADI algorithm parameters

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

;**get dim from first image header
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


;**For fits header keywords later
adi_nfwhm=nfwhm
adi_drsub=drsub0
adi_npca=npca
adi_rmin=rmin
adi_rmax=rmax

;load parallactic angles
readcol,'reduc.log',ffilenum,allxc,allyc,allrsat,allha,allpa,/silent

;Debugging
;print,total(ffilenum ne filenum),n_elements(filenum),n_elements(ffilenum)
;stop
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
;vestigal: probably do not need this
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
if (il_use eq !null) eq 0 then begin
coeffstouse=where(il_use eq il,ncoeffstouse)

if ncoeffstouse gt 0 then begin
c_use2=ck_use[coeffstouse]
;il_use2=il_use[coeffstouse]
eigenvect_use2=eigenvect_use[coeffstouse]
ir_use2=ir_use[coeffstouse]
it_use2=it_use[coeffstouse]
nf_use2=nf_use[coeffstouse]
endif
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
    ;print,'reading in file ',files[nf]


    norm=0
    if norm then begin
        print,'normalized'
        ;++normalise l'image par son profil radial de bruit
        ;normalize image by radial profile noise
        profrad,abs(im),2.,0.,rimmax,p2d=pr
        im/=pr
        writefits,tmpdir+filesprof[nf],pr
    endif

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

        ;calcule the noise in the annulus
        noise_im[ir,nf]=median(abs(im[ia]-median(im[ia])))/0.6745
    endfor

    el[nf]=sxpar(h0,'altitude')*!dtor

;**note: the GPI pipeline just reads in 'DEC' from the keyword.  Here, 
; keep original method of using DEC as determined from the elevation and azimuth

    dec[nf]=sxpar_charis(h0,'DEC',/justfirst,/silent)*!dtor
    decdeg[nf]=dec[nf]*!radeg

    ;print,'dec',decdeg[nf]
    ;comp=sxpar(h0,'DEC')
    ;print,comp,3600*(comp-decdeg[nf])

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
eigenvect_use3=eigenvect_use2[coeffstouse]
nf_use3=nf_use2[coeffstouse]
endif
endif
endif

    ri=rsub[ir] & dr=drsub[ir] & r=ri+dr/2. & rf=ri+dr
    print,'PCA Wavelength '+strtrim(il+1,2)+'/'+strtrim(n_elements(lambda),2),' Annulus '+strtrim(ir+1,2)+'/'+strtrim(nrsub,2)+' with radius '+$
        string(r,format='(f5.1)')+$
        ' [>='+string(ri,format='(f5.1)')+', <'+string(rf,format='(f5.1)')+']...'

    ;area of this region has radius of ...

    dropt=dr
    r1opt=ri
    r2opt=r1opt+dropt
    ;if dropt lt dr then begin
    ;    print,'dropt < drsub !!!'
    ;    print,'dropt: ',dropt
    ;    print,'drsub: ',dr
    ;    stop
    ;endif

    ;***determining the area optimization for this annulus removal

    r2opt = r2opt < rmax

    ;determines what image to load into memory
    i1aim=floor((r1opt-rmin)/drim)
    i2aim=floor((r2opt-rmin)/drim)
        ;print,i2aim,i1aim,n_elements(rim),rim[i2aim],sqrt(geom*area)
    if i2aim eq nrim then i2aim-=1
    if rim[i2aim] eq r2opt then i2aim-=1
    iaim=indgen(i2aim-i1aim+1)+i1aim


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

    ;removes the annuli that aren't necessary

    if ir gt 0 then begin
        irm=where(distarr[ia] lt rim[i1aim] or distarr[ia] ge rim[i2aim]+drim,crm,complement=ikp)
        if crm gt 0 then remove,irm,ia
        if crm gt 0 then annuli=annuli[ikp,*]
        if keyword_set(fwdmod) then begin
        if crm gt 0 then annuli_fwd=annuli_fwd[ikp,*]
        endif
    endif

    ia_tmp=0 & annuli_tmp=0

    ;remembers the list of annuli changes
    iaim_loaded=iaim

         iaopt=where(distarr[ia] ge r1opt and distarr[ia] lt r2opt)
   ;indices of pixels for optimization annulus

;Select annular zone

        dt=2.*!pi/float(nzone)
        
        for it=0,nzone-1 do begin
 ;indices of pixels included in this section: i.e. the optmization region

    if keyword_set(usecoeff) then begin
    if (it_use3 eq !null) eq 0 then begin
    coeffstouse=where(it_use3 eq it,ncoeffstouse)
    if ncoeffstouse gt 0 then begin
    c_use4=c_use3[coeffstouse]
    nf_use4=nf_use3[coeffstouse]
    eigenvect_use4=eigenvect_use3[coeffstouse]
    endif
    endif
    endif
       
        ;iopt=where(ang[ia] ge it*dt and ang[ia] lt (it+1)*dt) 

        iopt=where(ang[ia[iaopt]] ge it*dt and ang[ia[iaopt]] lt (it+1)*dt)
        iopt=iaopt[iopt]
        isub=where(distarr[ia[iopt]] ge ri and distarr[ia[iopt]] lt rf)
        isub=iopt[isub]

        ;indices of pixels included in this section: i.e. the optmization region

        ;instructs the region of optimization in memory

        ;indices of pixels to subtract


;region of 'optimization', really just the array of n_pixels_annulus by n_images
        
        ;there clues to avoid a build images later
        openw,lunit,tmpdir+'indices_images.dat',/get_lun,append=(ir +it gt 0)
        writeu,lunit,ia[isub]
        free_lun,lunit

;Subtract mean/median from each annulus
;you need to do this to do KLIP properly.  Can do this in the imrsub.pro step instead

       if (meansub gt 0) then begin
       meanopt=fltarr(nfiles)
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

       endif

        optreg=annuli[isub,*]
        if keyword_set(fwdmod) then optreg_fwd=annuli_fwd[isub,*]

;Retain David Lafreniere's outlier masking from original LOCI code; usually gains ~10% in SNR
 ;the following three lines are equivalent to the following
        ;inoise=floor((distarr(ia[iopt])-rim[0])/drim)#replicate(1,90)
        ;tmp=abs(optreg/noise_im[inoise])
        ;z<=(tmp lt 15.)
;        z<=(abs(optreg/noise_im[floor((distarr(ia[iopt])-rim[0])/drim)#replicate(1,nfiles)]) lt 15.)

        ;other way
        z=optreg
        z/=(replicate(1,n_elements(isub))#median(abs(z),dim=1))
        z/=(median(abs(z),dim=2)#replicate(1,nfiles))
        z=(abs(z) lt 7 and finite(z) eq 1)


        ;for a given pixel [i,*] look to see whether an image pixel that is deviant

        igoodf=where(min(z,dim=2),cgood)
        if cgood lt 1 then continue
        optreg=optreg[igoodf,*]
        ;isub=isub[igoodf]
        if keyword_set(fwdmod) then optreg_fwd=optreg_fwd[igoodf,*]

        ;loop on all images and made the last
        for nf=0,nfiles-1 do begin
            ;separation angulaire de toutes les images par rapport a image n
            ; angular separation of all pictures from a photo n

            dpa=abs(allpa-allpa[nf])

            ;OFFSET determined enough images for subtraction

            indim=where(dpa gt (nfwhm*fwhm/ri*!radeg+dtpose[nf]),c1)

;comment the igood stuff out for now
            igood=where(finite(annuli[isub,nf]) eq 1,c2)


;;Setting up the PCA

           if c1 lt 2 or c2 lt 5 then begin
                ;*debug*
                ;diff=fltarr(n_elements(isub))*0
                diff=fltarr(n_elements(isub))+!values.f_nan
                ;if keyword_set(fwdmod) then diffwwd=fltarr(n_elements(isub))+!values.f_nan
                difffwd=diff
            endif else begin
;print,il,nf,ir,it,c1
           npcaf = npca < c1
           data_mat=optreg[*,indim]
           if keyword_set(fwdmod) then data_mat_fwd=optreg_fwd[*,indim]
           ;data_mat=transpose(ref_library)

           covmatrix=matrix_multiply(data_mat,data_mat,/atranspose)
;;***forward-modeling
if keyword_set(usecoeff) then begin
;help,nf_use4,nf
;print,nf_use4
;print,nf
if (nf_use4 eq !null) eq 0 then begin
coeffstouse=where(nf_use4 eq nf,ncoeffstouse)
;print,ncoeffstouse
;print,c_use4
if ncoeffstouse gt 0 then begin
;c=ck_use[coeffstouse]
c=c_use4[coeffstouse]
eigenvect=eigenvect_use4[coeffstouse]

;now, reform the eigenvectors and eigenvalues
eigenvect=reform(eigenvect,c1,npcaf)
eigenval=(reform(c,c1,npcaf))[0,*]
;help,eigenvect,c,c1,npcaf
;help,eigenvect
endif
endif
goto,skipmatrixdecomp
endif


;usual eigenvalues/eigenvectors
           eigenval=la_eigenql(covmatrix,eigenvectors=eigenvect,$
            range=[c1-npcaf,c1-1],/double)

           eigenval=reverse(eigenval)
           ;help,eigenvect
           ;if npcaf eq 1 then eigenvect = transpose (eigenvect) else eigenvect=reverse(eigenvect,2)
           if npcaf gt 1 then eigenvect=reverse(eigenvect,2)

           skipmatrixdecomp:

;the full eigenvalues/eigenvectors
if keyword_set(fwdmod) then begin
           eigenval_full=la_eigenql(covmatrix,eigenvectors=eigenvect_full,$
            range=[0,c1-1],/double)
           eigenval_full=reverse(eigenval_full)
            ;if npcaf eq 1 then eigenvect_full = transpose (eigenvect_full) else eigenvect_full=reverse(eigenvect_full,2)
            if npcaf gt 1 then eigenvect_full=reverse(eigenvect_full,2)

          pc_full=matrix_multiply(eigenvect_full,data_mat,/atranspose,/btranspose)

endif

;print,'after ',il,nf,ir,it,c1
;help,eigenvect,data_mat
;if nf eq 100 then stop
           pc=matrix_multiply(eigenvect,data_mat,/atranspose,/btranspose)
     ;construct the reference

;****forward-modeling, calculate the self-subtraction component of the perturbed KL-modes
           if keyword_set(fwdmod) then dpc=matrix_multiply(eigenvect,data_mat_fwd,/atranspose,/btranspose)

;and the cross-term
           if keyword_set(fwdmod) then CadeltaS=matrix_multiply(data_mat_fwd,data_mat,/atranspose)+matrix_multiply(data_mat,data_mat_fwd,/atranspose)
           

            ref=fltarr(n_elements(isub))

             if keyword_set(fwdmod) then ref_fwd=fltarr(n_elements(isub))
             if keyword_set(savecoeff) then begin

             for k=0L,npcaf-1 do begin
             for jk=0,c1-1 do begin

             if (npcaf gt 1) then begin
             printf,1,long(il),long(ir),long(it),long(nf),eigenvect[jk,k],eigenval[k],jk,k,format='(i3,1x,i3,1x,i3,1x,i3,1x,f20.10,1x,e20.10,1x,i3,i3)'
             endif else begin
             ;printf,1,long(il),long(ir),long(it),long(nf),(transpose(eigenvect))[jk,k],eigenval[k],jk,k,format='(i3,i3,i3,i3,f20.10,1x,e20.10,1x,i3,i3)'
             printf,1,long(il),long(ir),long(it),long(nf),eigenvect[jk,k],eigenval[k],jk,k,format='(i3,1x,i3,1x,i3,1x,i3,1x,f20.10,1x,e20.10,1x,i3,i3)'
             endelse
;diagnostic tests             ;print,eigenvect[jk,k],eigenval[k],coeff*reform(pc[k,jk]),npcaf,c1
             endfor
             endfor
             endif



            for k=0L,npcaf-1 do begin

             pc[k,*]=pc[k,*]/sqrt(eigenval[k])
             coeff=total(reform(optreg[*,nf]*reform(pc[k,*])))


;diagnostic tests

             ref[igoodf]=ref[igoodf]+coeff*reform(pc[k,*])


             if keyword_set(fwdmod) then begin

             ;*****forward-model****
             ;perturbed KL models,k

;****1. direct self-subtraction
             direct_selfsub=dpc[k,*]/sqrt(eigenval[k])

;****2.indirect self-subtraction terms

;*******part 1
             ;help,transpose(eigenvect[*,k]),CadeltaS,eigenvect[*,k],pc[k,*],transpose(eigenvect_full[*,k]),eigenvect_full[*,k]
             first_indirect_term=(-0.5/eigenval[k])*transpose(eigenvect[*,k])#CadeltaS#eigenvect[*,k]#pc[k,*]

             second_indirect_term=first_indirect_term*0
;*******part 2
;loop for 2nd
             for jk=0,c1-1 do begin
             ;for jk=0,npcaf-1 do begin

             ;if (jk gt k) then begin
             if (jk ne k) then begin
              eigenval_con=sqrt(eigenval_full[jk]/eigenval[k])/(eigenval[k]-eigenval_full[jk])
              second_indirect_term+= eigenval_con*transpose(eigenvect_full[*,jk])#CadeltaS#eigenvect_full[*,k]#pc_full[jk,*]/sqrt(eigenval_full[jk])
             endif

             endfor

;****diagnostic tests
             ;dpc[k,*]=direct_selfsub
             ;dpc[k,*]=first_indirect_term
             ;dpc[k,*]=second_indirect_term
             ;dpc[k,*]=second_indirect_term+first_indirect_term


             ;dpc[k,*]=direct_selfsub+first_indirect_term
             ;dpc[k,*]=direct_selfsub+second_indirect_term

;*****deltaZ_k = V_k^TA_fwd/sqrt(lambda_k) - 0.5*V_k^T#C_adels#V_k#Z_k + sum(NR, j /=k) sqrt(lambda_j/lambda_k)*V_j^T#C_adels#V_k#Z_j/(lambda_k-lambda_j)
             dpc[k,*]=direct_selfsub+first_indirect_term+second_indirect_term

             coeff_oversub=total(reform(optreg_fwd[*,nf]*reform(pc[k,*]))) ;<S|KL>
             coeff_selfsub1=total(reform(optreg[*,nf]*reform(dpc[k,*])))   ;<N|DKL>KL
             coeff_selfsub2=total(reform(optreg[*,nf]*reform(pc[k,*])))    ;<N|KL>DKL

             ref_fwd[igoodf]=ref_fwd[igoodf]+coeff_oversub*reform(pc[k,*])+$
                1.*coeff_selfsub1*reform(pc[k,*]) + $
                1.*coeff_selfsub2*reform(dpc[k,*])
                ;0*coeff_selfsub1*reform(pc[k,*]) + $
                ;0*coeff_selfsub2*reform(dpc[k,*])


             endif
            endfor

            reftot=0.

                diff=annuli[isub,nf]-ref

;diagnostic tests
                ;diff=annuli[isub,nf]
                ;diff=-1*ref

;if you are forward-modeling, then output the annealed section of the synthetic signal, not the subtracted image section

                if keyword_set(fwdmod) then difffwd=annuli_fwd[isub,nf]-ref_fwd

                ;if keyword_set(fwdmod) then difffwd=annuli_fwd[isub,nf]
                ;if keyword_set(fwdmod) then difffwd=-1*ref_fwd

                diff=float(diff)
                if keyword_set(fwdmod) then difffwd=float(difffwd)

                if (zero gt 0) then diff-=median(diff,/even)
                if (zero gt 0) and keyword_set(fwdmod) then difffwd=float(difffwd - median(difffwd,/even))

;test images - debugging
if keyword_set(verbose) then begin
                testimage=fltarr(dim,dim)
                testimage[ia[isub]]=annuli[isub,nf]

                writefits,'test.fits',testimage
                testimage2=testimage*0
                testimage2[ia[isub]]=ref
                writefits,'test2.fits',testimage2

                ;test
                ;diff=annuli[isub,nf]
endif
                skipme:

            endelse


            ;register the difference, add (append) the values of this annulus to 
            ;binary file image
            openw,lunit,tmpdir+filestmp[nf],/get_lun,append=(ir +it gt 0)
            if ~keyword_set(fwdmod) then begin
            writeu,lunit,diff
            endif else begin
            writeu,lunit,difffwd
            endelse
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

    ;im=make_array(dim,dim,type=4,value=0)
    im=make_array(dim,dim,type=4,value=!values.f_nan)

    ;diagnostic ;bah=read_binary(tmpdir+filestmp[nf],data_type=4)
    ;bah=read_binary(tmpdir+filestmp[nf],data_type=4)

    im[ind]=read_binary(tmpdir+filestmp[nf],data_type=4)
    ;delete temporary files
    ;file_delete,tmpdir+filestmp[nf]

    if norm eq 1 then begin

        ;multiplied by the radial profile of noise
        pr=readfits(tmpdir+filesprof[nf],/silent)
        im*=pr
        file_delete,tmpdir+filesprof[nf]
    endif

    ;get header
    h1=headfits(reducdir1+files[nf],/exten,/silent)
    h0=headfits(reducdir1+files[nf],ext=0,/silent)
    ;h=headfits(reducdir+files[nf])

; ;***removed a lot of legacy code here useful for Gemini/NIRI but otherwise not.
;;**** see older code versions   

    ;rotate to bring the first image

    print,' rotation...'

;******astrometry*****

;******default is to rotate north-up using the tot-rot fits header keyword in the primary header
;********add option to ...
;******** - add an additional north PA offset (useful for future precise astro cals)
;******** - completely ignore the north PA rotation (potentially useful for backwards compatibility or a quick reduction)

;****if nothing then just rotate north-up*****

if ~keyword_set(nonorthup) then begin

theta=-1*allpa[nf]

charis_northup,im,h0,h1,angoffset=angoffset
endif else begin

theta=-(allpa[nf]-allpa[0])

if nf ne 0 then begin
im=rotat(im,theta,hdr=h1,missing=!values.f_nan)
endif

endelse

    ;rotskip:

    ;if nf ne 0 then im=rotat(im,theta,missing=!values.f_nan,hdr=h)

    ;register the difference
    if keyword_set(nonorthup) then begin
    sxaddhist,'A rotation of '+string(theta,format='(f8.3)')+$
      ' degrees was applied to the image.',h1
    endif else begin
    sxaddhist,'A rotation of '+string(northpa,format='(f8.3)')+$
      ' degrees was applied to the image.',h1
    endelse
; Add PCA/KLIP keywords
    sxaddpar,h1,'rsub',rsubval
    sxaddpar,h1,'klip_pca ',adi_npca
    sxaddpar,h1,'klip_nfwhm ',adi_nfwhm
    sxaddpar,h1,'klip_drsub ',adi_drsub
    sxaddpar,h1,'klip_rmin ',adi_rmin
    sxaddpar,h1,'klip_rmax ',adi_rmax
    sxaddpar,h1,'klip_nzone ',nzone
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


writefits,outname,0,h0
writefits,outname,im,h1,/append

writefits,reducdir+outname_col,0,h0
writefits,reducdir+outname_col,im_collapsed,h1,/append

file_delete,file_search(tmpdir,'*-klip*fits')
endofprogram:

if (keyword_set(fwdmod) or keyword_set(savecoeff) or keyword_set(usecoeff)) then free_lun,1

end
