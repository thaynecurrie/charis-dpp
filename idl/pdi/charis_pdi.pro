pro charis_pdi, pfname, prefname=prefname, Lsuffname=Lsuffname, Rsuffname=Rsuffname, hwpinfo_name=hwpinfo_name, outfile=outfile,$
fov_mask_pars=fov_mask_pars, meanadd=meanadd, angoffset=angoffset, first_order_ip_correct=first_order_ip_correct, fo_ip_rlims=fo_ip_rlims,$
azimuthal_stokes=azimuthal_stokes, phi_offset=phi_offset, single_sum_zeros=single_sum_zeros, help=help

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
    print,'charis_pdi.pro: Carries out single+double summing/differencing of PDI sequences, outputting Stokes parameter image cubes and wavelength collapsed images.'
    print,"Typically executed after 'charis_pdi_hwp_match'."
    print,'Written by K. Lawson (2020)'
    print,''
    print,'**Calling Sequence**'
    print,'charis_pdi, pfname, prefname=prefname, Lsuffname=Lsuffname, Rsuffname=Rsuffname, hwpinfo_name=hwpinfo_name, outfile=outfile,'
    print,'fov_mask_pars=fov_mask_pars, meanadd=meanadd, angoffset=angoffset, first_order_ip_correct=first_order_ip_correct, fo_ip_rlims=fo_ip_rlims,'
    print,'azimuthal_stokes=azimuthal_stokes, phi_offset=phi_offset, help=help'
    print,''
    print,'Example:'
    print,"charis_pdi, 'ab_aur_low.info', fov_mask_pars=[54,126,-26.9166], /first_order_ip_correct, /mean, fo_ip_rlims=[10,30], phi_offset=-0.271"
    print,''
    print,"***Keywords***"
    print,'pfname=STRING - The parameter file (e.g. HR8799_low.info)'
    print,"prefname=STRING - The prefix of files to use as input. Defaults to 'n'."
    print,"Lsuffname=STRING - The suffix of left pol files to use as input. Defaults to 'leftreg_cal', or 'leftreg' if no images are found for 'leftreg_cal'."
    print,"Rsuffname=STRING - The suffix of right pol files to use as input. Defaults to 'rightreg_cal', or 'rightreg' if no images are found for 'rightreg_cal'."
    print,"hwpinfo_name=STRING - The path for the half-wave plate info file produced by charis_pdi_hwp_match. Defaults to 'hwp_info.txt'."
    print,"outfile=STRING - prefix to use in naming final sequence-averaged PDI images and image cubes. Defaults to 'pdi'."
    print,"fov_mask_pars=ARRAY - Three element array providing the x-axis width (pixels), y-axis width (pixels), and angle (degrees) to use in creating a mask of the valid field of view. Defaults to the full frame."
    print,"/meanadd - (as in charis_adialoci) Use an outlier resistant mean rather than the median to combine exposures and to wavelength collapse image cubes."
    print,"angoffset=FLOAT - Offset of the nominal image PA from true north. Defaults to the value for 'angoffset' from charis_get_constant."
    print,"/first_order_ip_correct - Perform simple first order instrumental polarization correction. Should generally be used until Mueller Matrix corrections available."
    print,"fo_ip_rlims=ARRAY - Two element array indicating the inner and outer radii of the annulus to use in first order IP correction; no effect if /first_order_ip_correct is not set. Defaults to [10,30]."
    print,"/azimuthal_stokes - Calculate azimuthal stokes parameters (Qphi, Uphi) in addition to typical linear Stokes parameters."
    print,"phi_offset - Offset for the angle of linear polarization (AOLP) in radians to assume when computing azimuthal Stokes parameters. No effect if /azimuthal_stokes is not set. Defaults to 0."
    print,"/single_sum_zeros - Apply zeros instead of nans to pixels beyond the defined FOV for the output single sum image cubes."
    goto,endofprogram
endif

reducdir='./reduc/'
datadir=reducdir+'reg/'
procdir=reducdir+'proc/'

if ~keyword_set(angoffset) then angoffset=charis_get_constant(name='angoffset') ;nominally 2.2 deg

if ~keyword_set(hwpinfo_name) then hwpinfo_name = 'hwp_info.txt'

readcol,hwpinfo_name,name,hwp_pos,hwp_ang,hwp_cycle,/silent,skipline=1,format="(A,I,F,I)"

if ~keyword_set(prefname) then prefname = 'n'

if ~keyword_set(Lsuffname) then begin
    files=FILE_SEARCH(datadir+prefname+'*leftreg_cal.fits', count=cal_count)
    if cal_count eq 0 then Lsuffname='leftreg' else Lsuffname='leftreg_cal'
endif

if ~keyword_set(Rsuffname) then begin
    files=FILE_SEARCH(datadir+prefname+'*rightreg_cal.fits', count=cal_count)
    if cal_count eq 0 then Rsuffname='rightreg' else Rsuffname='rightreg_cal'
endif

if ~keyword_set(fo_ip_rlims) then fo_ip_rlims=[10,30]

if ~keyword_set(outfile) then outfile='pdi'

if ~keyword_set(phi_offset) then phi_offset=0.

readcol,'reduc_left.log',ffilenum,allxc,allyc,allrsat,allha,allpa,/silent

param,'fnum_sat',flist,/get,pfname=pfname

filenum_pol= nbrlist(flist)
files_polleft = datadir+filelist(filenum_pol,nfiles_polpair,prefix=prefname,suffix=Lsuffname)
files_polright = datadir+filelist(filenum_pol,nfiles_polpair,prefix=prefname,suffix=Rsuffname)

imtest = readfits(files_polleft[0],h1test,/exten,/silent)
sz = size(imtest, /dim)
dimx = sz[0]
dimy = sz[1]

if ~keyword_set(fov_mask_pars) then fov_mask_pars=[dimx,dimy,0]
charis_generate_roi,output=fov_mask,roirange=fov_mask_pars,roidim=[dimx,dimy]
outside_fov = where(fov_mask eq 0)

valid_cycles = where(hwp_cycle ne -1)
ncycles= n_elements(uniq(hwp_cycle[valid_cycles], sort(hwp_cycle[valid_cycles]))) ; Exclude the unused cycle=-1 frames (if any)
cyclenums='1-'+(string(ncycles)).trim()
cycle_filenums = nbrlist(cyclenums)

IQout= datadir+filelist(cycle_filenums, ncyc, prefix='c', suffix='_iq'); Saving tot intensity products to 'reduc/reg/' so they can be used with PSFsub programs without alteration
IUout= datadir+filelist(cycle_filenums, ncyc, prefix='c', suffix='_iu')
Iout= datadir+filelist(cycle_filenums, ncyc, prefix='c', suffix='_i')

Qout= procdir+filelist(cycle_filenums, ncyc, prefix='c', suffix='_q')
Uout= procdir+filelist(cycle_filenums, ncyc, prefix='c', suffix='_u')
PIout= procdir+filelist(cycle_filenums, ncyc, prefix='c', suffix='_pi')

q_dpas = fltarr(ncycles)
u_dpas = fltarr(ncycles)

; Separately computing all single sums and writing to ./reduc/reg with standard total intensity filenames:
for Ti=0,nfiles_polpair-1 do begin
    h0L = headfits(files_polleft[Ti], ext=0, /silent)
    leftcube= readfits(files_polleft[Ti],h1L,/exten,/silent)
    rightcube= readfits(files_polright[Ti],h1R,/exten,/silent)
    imcube = leftcube+rightcube
    for Li=0,sz[2]-1 do begin
        imslice = imcube[*,*,Li]
        if keyword_set(single_sum_zeros) then imslice[outside_fov] = 0. else imslice[outside_fov] = !values.f_nan
        imcube[*,*,Li] = imslice
    endfor
    writefits,(files_polleft[Ti].replace(Lsuffname, 'reg_cal')),0,h0L
    writefits,(files_polleft[Ti].replace(Lsuffname, 'reg_cal')),imcube,h1L,/append
endfor

for ci=0,ncycles-1 do begin
    cycle_inds = where(hwp_cycle eq ci, Ncyclefiles)
    if Ncyclefiles ne 4 then begin
        s = string(ci)
        print,'For cycle ',s.trim(),' more/fewer than 4 images are indicated in ',hwpinfo_name
        print,'Please ensure that each cycle index (except -1) appears exactly 4 times!'
        break
    endif

    Qp_ind= cycle_inds[where(hwp_pos[cycle_inds] eq 0)]
    Up_ind= cycle_inds[where(hwp_pos[cycle_inds] eq 1)]
    Qm_ind= cycle_inds[where(hwp_pos[cycle_inds] eq 2)]
    Um_ind= cycle_inds[where(hwp_pos[cycle_inds] eq 3)]

    cycle_pas = [allpa[Qp_ind], allpa[Up_ind], allpa[Qm_ind], allpa[Um_ind]]
 
    left_cycle = fltarr(sz[0],sz[1],sz[2],4)
    right_cycle = fltarr(sz[0],sz[1],sz[2],4)

    left_cycle[*,*,*,0]= readfits(files_polleft[Qp_ind],h1_qpL,/exten,/silent) ; IQ+L
    left_cycle[*,*,*,1]= readfits(files_polleft[Up_ind],h1_upL,/exten,/silent) ; IU+L
    left_cycle[*,*,*,2]= readfits(files_polleft[Qm_ind],h1_qmL,/exten,/silent) ; IQ-L
    left_cycle[*,*,*,3]= readfits(files_polleft[Um_ind],h1_umL,/exten,/silent) ; IU-L

    right_cycle[*,*,*,0]= readfits(files_polright[Qp_ind],h1_qpR,/exten,/silent) ; IQ+R
    right_cycle[*,*,*,1]= readfits(files_polright[Up_ind],h1_upR,/exten,/silent) ; IU+R
    right_cycle[*,*,*,2]= readfits(files_polright[Qm_ind],h1_qmR,/exten,/silent) ; IQ-R
    right_cycle[*,*,*,3]= readfits(files_polright[Um_ind],h1_umR,/exten,/silent) ; IU-R

    badl=where(left_cycle eq 0)
    badr=where(right_cycle eq 0)
    left_cycle[badl] = !values.f_nan
    right_cycle[badr] = !values.f_nan
    csz = size(left_cycle,/dim)
    for Li=0,csz[2]-1 do begin
        for i=0,csz[3]-1 do begin
            left_slice = left_cycle[*,*,Li,i]
            left_slice[outside_fov] = !values.f_nan
            left_cycle[*,*,Li,i] = left_slice
            right_slice = right_cycle[*,*,Li,i]
            right_slice[outside_fov] = !values.f_nan
            right_cycle[*,*,Li,i] = right_slice
        endfor
    endfor

    ; Single sums
    IQp = left_cycle[*,*,*,0] + right_cycle[*,*,*,0]
    IUp = left_cycle[*,*,*,1] + right_cycle[*,*,*,1]
    IQm = left_cycle[*,*,*,2] + right_cycle[*,*,*,2]
    IUm = left_cycle[*,*,*,3] + right_cycle[*,*,*,3]

    ; Double sums
    IQ = 0.5*(IQp + IQm)
    IU = 0.5*(IUp + IUm)

    ; Total intensity
    total_intensity = 0.5*(IQ + IU)

    ; Single diffs
    Qp = left_cycle[*,*,*,0] - right_cycle[*,*,*,0]
    Up = left_cycle[*,*,*,1] - right_cycle[*,*,*,1]
    Qm = left_cycle[*,*,*,2] - right_cycle[*,*,*,2]
    Um = left_cycle[*,*,*,3] - right_cycle[*,*,*,3]

    ; Double diffs
    PQ = 0.5*(Qp-Qm)
    PU = 0.5*(Up-Um)

    if keyword_set(first_order_ip_correct) then begin
        dim = sz[0]
        distarr = shift(dist(dim),dim/2,dim/2)
        annulus = where((distarr ge fo_ip_rlims[0]) and (distarr le fo_ip_rlims[1]))
        for Li=0,sz[2]-1 do begin
            Qslice = PQ[*,*,Li]
            IQslice = IQ[*,*,Li]
            Uslice = PU[*,*,Li]
            IUslice = IU[*,*,Li]

            Qratio = Qslice/IQslice
            Uratio = Uslice/IUslice
            cq = median(Qratio[annulus], /even)
            cu = median(Uratio[annulus], /even)

            Qips = Qslice-cq*IQslice
            Uips = Uslice-cu*IUslice

            PQ[*,*,Li] = Qips
            PU[*,*,Li] = Uips
        endfor
    endif 

    ; Polarized intensity
    polarized_intensity = sqrt(PQ^2 + PU^2)

    qpa = 0.5*(cycle_pas[0] + cycle_pas[2])
    upa = 0.5*(cycle_pas[1] + cycle_pas[3])

    h0 = headfits(files_polleft[Up_ind], ext=0, /silent)
    h1 = h1_upL

    q_dpas[ci] = qpa-allpa[Up_ind]
    u_dpas[ci] = upa-allpa[Up_ind]

    ; Saving IQ, Q products...
    orignames=sxpar(h1_qpL,'origname')+', '+sxpar(h1_qmL,'origname')
    sxaddpar,h1,'cycle',orignames,'origname entries for frames in this cycle'
    writefits,IQout[ci],0,h0
    writefits,IQout[ci],IQ,h1,/append
    writefits,Qout[ci],0,h0
    writefits,Qout[ci],PQ,h1,/append

    ; Saving IU, U products...
    orignames=sxpar(h1_upL,'origname')+', '+sxpar(h1_umL,'origname')
    sxaddpar,h1,'cycle',orignames,'origname entries for frames in this cycle'
    writefits,IUout[ci],0,h0
    writefits,IUout[ci],IU,h1,/append
    writefits,Uout[ci],0,h0
    writefits,Uout[ci],PU,h1,/append

    ; Saving I, PI products...
    orignames=sxpar(h1_qpL,'origname')+', '+sxpar(h1_upL,'origname')+', '+sxpar(h1_qmL,'origname')+', '+sxpar(h1_umL,'origname')
    sxaddpar,h1,'cycle',orignames,'origname entries for frames in this cycle'
    writefits,Iout[ci],0,h0
    writefits,Iout[ci],total_intensity,h1,/append
    writefits,PIout[ci],0,h0
    writefits,PIout[ci],polarized_intensity,h1,/append
endfor

qhcube = fltarr(sz[0],sz[1],sz[2],ncycles)
uhcube = fltarr(sz[0],sz[1],sz[2],ncycles)
for ci=0,ncycles-1 do begin
    qh0 = headfits(Qout[ci], ext=0, /silent)
    qcube = readfits(Qout[ci], qh1, /exten, /silent)
    qangoffset = angoffset + q_dpas[ci]
    charis_northup,qcube,qh0,qh1,angoffset=qangoffset,outputrot=qoutputrot

    Qout_derot = Qout[ci].replace('.fits', '_derot.fits')
    writefits,Qout_derot,0,qh0
    writefits,Qout_derot,qcube,qh1,/append

    uh0 = headfits(Uout[ci], ext=0, /silent)
    ucube = readfits(Uout[ci], uh1, /exten, /silent)
    uangoffset = angoffset + u_dpas[ci]
    charis_northup,ucube,uh0,uh1,angoffset=uangoffset,outputrot=uoutputrot

    Uout_derot = Uout[ci].replace('.fits', '_derot.fits')
    writefits,Uout_derot,0,uh0
    writefits,Uout_derot,ucube,uh1,/append

    qhcube[*,*,*,ci] = qcube
    uhcube[*,*,*,ci] = ucube

    if ci eq 0 then begin
        h0_main = qh0
        h1_main = qh1
        h1q_main = qh1
        h1u_main = uh1
    endif

    cycle_str = 'cyc_'+((string(ci+1)).trim())
    qorignames = sxpar(qh1,'cycle')
    sxaddpar,h1q_main,cycle_str,qorignames
    uorignames = sxpar(uh1,'cycle')
    sxaddpar,h1u_main,cycle_str,uorignames
    orignames = qorignames+', '+uorignames
    sxaddpar,h1_main,cycle_str,orignames
endfor

h1s = list(h1_main, h1q_main, h1u_main)
for i=0,2 do begin
    h1 = h1s[i]
    sxdelpar,h1,'cycle'; Removing "cycle" entry from original combo
    sxaddpar,h1,'prefname',prefname
    sxaddpar,h1,'Lsuffname',Lsuffname
    sxaddpar,h1,'Rsuffname',Rsuffname
    sxaddpar,h1,'hwpinfo_name',hwpinfo_name
    sxaddpar,h1,'fov_mask_pars',strjoin((string(fov_mask_pars)).trim(),', ')
    sxaddpar,h1,'meanadd',meanadd
    sxaddpar,h1,'angoffset',angoffset
    sxaddpar,h1,'first_order_ip_correct',first_order_ip_correct
    sxaddpar,h1,'fo_ip_rlims',strjoin((string(fo_ip_rlims)).trim(),', ')
    sxaddpar,h1,'phi_offset',phi_offset
    h1s[i] = h1
endfor
h1_main = h1s[0]
h1q_main = h1s[1]
h1u_main = h1s[2]

if ~keyword_set(meanadd) then begin
    qcube=median(qhcube, dimension=4, /even)
    ucube=median(uhcube, dimension=4, /even)
endif else begin
    resistant_mean,qhcube,3.,qcube,dimension=4
    resistant_mean,uhcube,3.,ucube,dimension=4
endelse

writefits,procdir+outfile+'_q.fits',0,h0_main
writefits,procdir+outfile+'_q.fits',qcube,h1q_main,/append

writefits,procdir+outfile+'_u.fits',0,h0_main
writefits,procdir+outfile+'_u.fits',ucube,h1u_main,/append

pi_cube = sqrt(qcube^2 + ucube^2)
writefits,procdir+outfile+'_pi.fits',0,h0_main
writefits,procdir+outfile+'_pi.fits',pi_cube,h1_main,/append

; Collapse Q and U wavelength cubes to a single broadband image, calculate products again
if ~keyword_set(meanadd) then begin
    qcol=median(qcube, dimension=3, /even)
    ucol=median(ucube, dimension=3, /even)
endif else begin
    resistant_mean,qcube,3.,qcol,dimension=3
    resistant_mean,ucube,3.,ucol,dimension=3
endelse

writefits,procdir+outfile+'_q_collapsed.fits',0,h0_main
writefits,procdir+outfile+'_q_collapsed.fits',qcol,h1q_main,/append
writefits,procdir+outfile+'_u_collapsed.fits',0,h0_main
writefits,procdir+outfile+'_u_collapsed.fits',ucol,h1u_main,/append

pi_col = sqrt(qcol^2 + ucol^2)
writefits,procdir+outfile+'_pi_collapsed.fits',0,h0_main
writefits,procdir+outfile+'_pi_collapsed.fits',pi_col,h1_main,/append

if keyword_set(azimuthal_stokes) then begin
    phi = angarr(dimx)
    phi = transpose(phi)
    phi += phi_offset

    lam_offsets = [0.0819,  0.1237,  0.1562,  0.1893,  0.2035,  0.2062,  0.1885, 0.1431,  0.0918,  0.048 , -0.0088, -0.0441, -0.0746, -0.092 , -0.1026, -0.1178, -0.1329, -0.1422, -0.1464, -0.1419, -0.1346, -0.1192]

    phi_cube = fltarr(sz[0],sz[1],sz[2])
    for Li=0,sz[2]-1 do phi_cube[*,*,Li] = phi + lam_offsets[Li]

    qphi_cube = (-qcube*cos(2.*phi_cube))-(ucube*sin(2.*phi_cube))
    uphi_cube = (qcube*sin(2.*phi_cube))-(ucube*cos(2.*phi_cube))

    writefits,procdir+outfile+'_qphi.fits',0,h0_main
    writefits,procdir+outfile+'_qphi.fits',qphi_cube,h1_main,/append
    writefits,procdir+outfile+'_uphi.fits',0,h0_main
    writefits,procdir+outfile+'_uphi.fits',uphi_cube,h1_main,/append

    qphi_col = (-qcol*cos(2.*phi))-(ucol*sin(2.*phi))
    uphi_col = (qcol*sin(2.*phi))-(ucol*cos(2.*phi))
    writefits,procdir+outfile+'_qphi_collapsed.fits',0,h0_main
    writefits,procdir+outfile+'_qphi_collapsed.fits',qphi_col,h1_main,/append
    writefits,procdir+outfile+'_uphi_collapsed.fits',0,h0_main
    writefits,procdir+outfile+'_uphi_collapsed.fits',uphi_col,h1_main,/append
endif
endofprogram:
end
