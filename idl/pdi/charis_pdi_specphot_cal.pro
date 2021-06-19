pro charis_pdi_specphot_cal,pfname,$
calmethod=calmethod,pick=pick,datacube=datacube,calcube=calcube,modamp=modamp,subskyannulus=subskyannulus,$
filtcal_slice=filtcal_slice,nopradcal=nopradcal,$
notimeratio=notimeratio,$
meancomb=meancomb,$
fixradius=fixradius,$
ap_factor=ap_factor,$
filtername=filtername,$
starlib=starlib,$
empspectrum=empspectrum,$
av=av,$
prefname=prefname,$
Lsuffname=Lsuffname,$
Rsuffname=Rsuffname,$
verbose=verbose,$
fluxunits=fluxunits,$
outfilename=outfilename,$
fov_mask_pars=fov_mask_pars,$
unsat=unsat,$
help=help

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
    print,'charis_pdi_specphot_cal.pro: Flux calibrates PDI data using charis_specphot_cal to compute scaling factors for the single sum sequence.'
    print,"Typically executed after 'charis_pdi_register_cube' and before 'charis_pdi_hwp_match'."
    print,'Written by K. Lawson (2020)'
    print,''
    print,'**Calling Sequence**'
    print,'charis_pdi_specphot_cal,pfname,prefname=prefname,Lsuffname=Lsuffname,Rsuffname=Rsuffname,fov_mask_pars=fov_mask_pars,'
    print,'calmethod=calmethod,pick=pick,datacube=datacube,calcube=calcube,modamp=modamp,subskyannulus=subskyannulus,'
    print,'filtcal_slice=filtcal_slice,nopradcal=nopradcal,meancomb=meancomb,ap_factor=ap_factor,filtername=fname,starname=starname,mag=mag,starlib=starlib,'
    print,'av=av,prefname=prefname,suffname=suffname,test=test,verbose=verbose,fluxunits=fluxunits,outfilename=outfilename,help=help'
    print,''
    print,'Example:'
    print,"charis_pdi_specphot_cal, 'abaur_low.info', fov_mask_pars=[54,126,-26.9166], modamp=25, ap_factor=5, starlib=3, empspectrum='abaur_prism_071210.fits'"
    print,''
    print,"***Keywords***"
    print,'pfname=STRING - The parameter file (e.g. HR8799_low.info)'
    print,"prefname=STRING - The prefix of files to use as input. Defaults to 'n'."
    print,"Lsuffname=STRING - The suffix of left pol files to use as input. Defaults to 'leftreg'."
    print,"Rsuffname=STRING - The suffix of right pol files to use as input. Defaults to 'rightreg'."
    print,"/unsat - Use charis_specphotcal_unsat instead of charis_specphot_cal; appropriate for unsaturated, unocculted data"
    print,"fov_mask_pars=ARRAY - Three element array providing the x-axis width (pixels), y-axis width (pixels), and angle (degrees) to use in creating a mask of the valid field of view. Defaults to the full frame."
    print,"Other keywords as accepted by 'charis_specphot_cal'."
    goto,endofprogram
endif

reducdir='./reduc/'

;data directory
datadir=reducdir+'reg/'
subdir='reg/'
reducdir+=subdir

if ~keyword_set(prefname) then prefname = 'n'
if ~keyword_set(Lsuffname) then Lsuffname = 'leftreg'
if ~keyword_set(Rsuffname) then Rsuffname = 'rightreg'

param,'fnum_sat',flist,/get,pfname=pfname

filenum_pol = nbrlist(flist)
files_polleft = reducdir+filelist(filenum_pol, nfiles, prefix=prefname, suffix=Lsuffname)
files_polright = reducdir+filelist(filenum_pol, nfiles, prefix=prefname, suffix=Rsuffname)
files_tmp = reducdir+filelist(filenum_pol, nfiles, prefix='n', suffix='_tmp')

imtest = readfits(files_polleft[0],h1test,/exten,/silent)
sz = size(imtest, /dim)
dimx = sz[0]
dimy = sz[1]

if ~keyword_set(fov_mask_pars) then fov_mask_pars=[dimx,dimy,0]
charis_generate_roi,output=fov_mask,roirange=fov_mask_pars,roidim=[dimx,dimy]
outside_fov = where(fov_mask eq 0)

for i=0,nfiles-1 do begin
    Lcube = readfits(files_polleft[i],h1L,/exten,/silent) ; Left pol
    Rcube = readfits(files_polright[i],h1R,/exten,/silent) ; Right pol

    badl=where(Lcube eq 0)
    badr=where(Rcube eq 0)
    Lcube[badl] = !values.f_nan
    Rcube[badr] = !values.f_nan
    for Li=0,sz[2]-1 do begin
        left_slice = Lcube[*,*,Li]
        left_slice[outside_fov] = !values.f_nan
        Lcube[*,*,Li] = left_slice
        right_slice = Rcube[*,*,Li]
        right_slice[outside_fov] = !values.f_nan
        Rcube[*,*,Li] = right_slice
    endfor

    Icube = Lcube+Rcube ; Single sum cube on which to perform specphotcal
    bad = where(~finite(Icube))
    Icube[bad] = 0.
    
    h0 = headfits(files_polleft[i], ext=0, /silent)
    writefits,files_tmp[i],0,h0
    writefits,files_tmp[i],Icube,h1L,/append
endfor 

; Run charis_specphot_cal on our intensity images.
if keyword_set(unsat) then begin
    charis_specphot_cal_unsat, pfname, calmethod=calmethod, pick=pick, datacube=datacube, calcube=calcube, subskyannulus=subskyannulus, filtcal_slice=filtcal_slice, nopradcal=nopradcal,$
    meancomb=meancomb, fixradius=fixradius, ap_factor=ap_factor, filtername=filtername, starlib=starlib, empspectrum=empspectrum, av=av,$
    prefname='n', suffname='_tmp', verbose=verbose, fluxunits=fluxunits, outfilename=outfilename, help=help
endif else begin
    charis_specphot_cal, pfname, calmethod=calmethod, pick=pick, datacube=datacube, calcube=calcube, modamp=modamp, subskyannulus=subskyannulus, filtcal_slice=filtcal_slice, nopradcal=nopradcal, notimeratio=notimeratio,$
    meancomb=meancomb, fixradius=fixradius, ap_factor=ap_factor, filtername=filtername, starlib=starlib, empspectrum=empspectrum, av=av,$
    prefname='n', suffname='_tmp', verbose=verbose, fluxunits=fluxunits, outfilename=outfilename, help=help
endelse

files_tmp_cal = reducdir+filelist(filenum_pol, nfiles, prefix='n', suffix='_tmp_cal')
files_polleft_cal = reducdir+filelist(filenum_pol, nfiles, prefix=prefname, suffix=Lsuffname+'_cal')
files_polright_cal = reducdir+filelist(filenum_pol, nfiles, prefix=prefname, suffix=Rsuffname+'_cal')

for i=0,nfiles-1 do begin
    h1cal = headfits(files_tmp_cal[i], ext=1, /silent)

    fluxunit = sxpar(h1cal,'fluxunit',comment=fluxunit_comment)
    sky_in = sxpar(h1cal,'sky_in',comment=sky_in_comment)
    sky_out = sxpar(h1cal,'sky_out',comment=sky_out_comment)

    Lcube = readfits(files_polleft[i],h1L,/exten,/silent)
    h0L = headfits(files_polleft[i],ext=0,/silent)
    sxaddpar,h1L,'fluxunit',fluxunit,fluxunit_comment
    sxaddpar,h1L,'sky_in',sky_in,sky_in_comment
    sxaddpar,h1L,'sky_out',sky_out,sky_out_comment

    Rcube = readfits(files_polright[i],h1R,/exten,/silent)
    h0R = headfits(files_polright[i],ext=0,/silent)
    sxaddpar,h1R,'fluxunit',fluxunit,fluxunit_comment
    sxaddpar,h1R,'sky_in',sky_in,sky_in_comment
    sxaddpar,h1R,'sky_out',sky_out,sky_out_comment

    for Li=0,sz[2]-1 do begin ;Li == lambda index
        sLi = string(Li)
        sLi = sLi.trim()
        fstar_Li = sxpar(h1cal,'fstar_'+sLi,comment=fstar_Li_comment)
        r_ap_Li = sxpar(h1cal,'r_ap'+sLi,comment=r_ap_Li_comment)
        fscale_Li = sxpar(h1cal,'fscale'+sLi,comment=fscale_Li_comment)
        cerr_Li = sxpar(h1cal,'cerr'+sLi,comment=cerr_Li_comment)
        
        Lcube[*,*,Li] *= fscale_Li
        sxaddpar, h1L, 'fstar_'+sLi, fstar_Li, fstar_Li_comment
        sxaddpar, h1L, 'r_ap'+sLi, r_ap_Li, r_ap_Li_comment
        sxaddpar, h1L, 'fscale'+sLi, fscale_Li, fscale_Li_comment
        sxaddpar, h1L, 'cerr'+sLi, cerr_Li, cerr_Li_comment

        Rcube[*,*,Li] *= fscale_Li
        sxaddpar, h1R, 'fstar_'+sLi, fstar_Li, fstar_Li_comment
        sxaddpar, h1R, 'r_ap'+sLi, r_ap_Li, r_ap_Li_comment
        sxaddpar, h1R, 'fscale'+sLi, fscale_Li, fscale_Li_comment
        sxaddpar, h1R, 'cerr'+sLi, cerr_Li, cerr_Li_comment
    endfor

    writefits,files_polleft_cal[i],0,h0L
    writefits,files_polleft_cal[i],Lcube,h1L,/append

    writefits,files_polright_cal[i],0,h0R
    writefits,files_polright_cal[i],Rcube,h1R,/append
    spawn,'rm '+files_tmp[i]
    spawn,'rm '+files_tmp_cal[i]
endfor
endofprogram:
end           