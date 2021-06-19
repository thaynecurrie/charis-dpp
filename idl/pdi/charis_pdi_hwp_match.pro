pro charis_pdi_hwp_match, pfname, hwpfilein=hwpfilein, prefname=prefname, suffname=suffname, hwpfileout=hwpfileout, help=help

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
    print,'charis_pdi_hwp_match.pro: Uses HWP log to assign a HWP angle to each frame,'
    print,'then match to appropriate frames for PDI functionality.'
    print,"Typically executed after 'charis_pdi_specphot_cal' and before 'charis_pdi'."
    print,'Written by K. Lawson (2021)'
    print,''
    print,'**Calling Sequence**'
    print,"charis_pdi_hwp_match,pfname,hwpfile=hwpfile, prefname=prefname, suffname=suffname, hwpfileout=hwpfileout"
    print,''
    print,'Example:'
    print,"charis_pdi_hwp_match,'abaur_low.info'"
    print,''
    print,"***Keywords***"
    print,'pfname=STRING - The parameter file (e.g. HR8799_low.info)'
    print,"hwpfilein=STRING - path to HWP log from Subaru. Default: searches for filename like 'vampHWPLog_*.csv'"
    print,"prefname=STRING - The prefix of files to use as input. Defaults to 'n'."
    print,"suffname=STRING - The suffix of files to use as input. Defaults to 'leftreg_cal' if present, or 'leftreg' otherwise."
    print,"hwpfileout=STRING - File into which to save HWP angles and matches. Defaults to 'hwp_info.txt'."
    goto,endofprogram
endif

reducdir='./reduc/'
datadir=reducdir+'reg/'

if ~keyword_set(hwpfilein) then begin
    hwp_search = file_search('./', 'vampHWPLog_*.csv', count=nfiles)
    if nfiles eq 1 then begin
        hwpfilein = hwp_search[0]
        print,'Using HWP log file:'
        print,hwpfilein
    endif else begin
        if nfiles eq 0 then begin
            print,'No HWP log file found (e.g. "vampHWPLog_*.csv).'
            print,'Please add a HWP log file to your working directory,'
            print,'or explicitly indicate the file to use by setting'
            print,'the keyword "hwpfilein".'
        endif else begin
            print,'Multiple HWP log files found:'
            print, hwp_search
            print,'Please explicitly indicate the file to use by setting'
            print,'the keyword "hwpfilein" to the file you wish to use.'
        endelse
        goto,endofprogram
    endelse
endif

param,'fnum_sat',flist,/get,pfname=pfname
if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then begin
    test=file_search(datadir+'*leftreg_cal.fits')
    if n_elements(test) gt 1 then begin
        suffname='leftreg_cal'
    endif else begin
        test=file_search(datadir+'*leftreg.fits')
        if n_elements(test) gt 1 then begin
            suffname='leftreg'
        endif else begin
            datadir=reducdir+'prep/'
            suffname='left'
        endelse
    endelse
endif

print,'Using prefname: ',prefname
print,'Using suffname: ',suffname

if ~keyword_set(hwpfileout) then hwpfileout='hwp_info.txt'

filenum=nbrlist(flist)
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)

hwp_ang=fltarr(nfiles)
hwp_pos=intarr(nfiles)
utcexp=fltarr(nfiles)
parangs=fltarr(nfiles)

readcol,hwpfilein,angs,utstart_string,utend_string,target,format='(f,a,a,a)'
uni_angs = angs[uniq(angs, sort(angs))]

; There are at least two utilized formats for PDI logs: 
; 1) start and end time columns are formatted as: hh:mm:ss.sss
; 2) yyyymmddThhmmss (where the T is literally the letter T)
; Current implementation will cause problems if dates are needed.
utchwp=fltarr(n_elements(angs))
for i=0L,n_elements(angs)-1 do begin
    pos = STRPOS(utend_string[i], ':')
    if pos eq -1 then begin
        datetime = strsplit(utend_string[i],/extract,'T')
        hour = strmid(datetime[1],0,2)
        minute = strmid(datetime[1],2,2)
        second = strmid(datetime[1],4,2)
        utchwp[i]=ten(hour,minute,second)
        time_format = 'fits'
    endif else begin
        time=strsplit(utend_string[i],/extract,':')   
        utchwp[i]=ten(time[0],time[1],time[2])
        time_format = 'utc'
    endelse
endfor

; There is some ambiguity with fits headers here. 
; 'UTC-TIME' claims to be average exposure time (call this case 1), but this would mean HWP moves DURING many exposures.
; I suspect it really gives the start time of the exposure (case 2), but could also be the end time (case 3)
; Here, we assume case 2 is correct. For this, we don't need the exposure time and should just match each 'UTC-TIME' 
; with the immediately preceding hwp movement end time. 
; i.e. the last HWP pos we finished moving to before starting the observation must be the HWP pos of that observation

; ; This version assumes that 'UTC-TIME' gives the start time of the observation.
; Note: this appears to give the incorrect solution for the TW Hya data. However, other versions give incorrect soln for AB Aur...
; Compare single diffs for TW Hya exposures: n0026 and n0027. Matching this way indicates the same HWP position for them
; But the HWP positions clearly differ (average cubes to H-band to make difference obvious)
if time_format eq 'utc' then begin
    for i=0L,nfiles-1 do begin
        h=headfits(datadir+files[i],ext=1,/silent)
        timestring=sxpar(h,'UTC-TIME')
        time=strsplit(timestring,/extract,':')
        exptime = float(sxpar(h, 'EXPTIME'))
        utcexp[i]=ten(time[0],time[1],time[2])
        parangs[i]=sxpar(h,'PARANG')
        dts = abs(utchwp-utcexp[i])
        min_dt = min(dts, argmin) ;Get index of nearest hwp move end (in time) to utcexp[i]
        if utchwp[argmin] lt utcexp[i] then hwpind = argmin ;If corresponding time is before utcexp[i], then that index gives the matching hwp pos
        if utchwp[argmin] gt utcexp[i] then hwpind = argmin-1 ;If corresponding time is after utcexp[i], then we want entry for one index prior!
        hwp_ang[i] = angs[hwpind]
        min_diff = min(abs(hwp_ang[i]-uni_angs), angpos)
        hwp_pos[i] = angpos
    endfor
endif

; ; This version assumes 'UTC-TIME' gives the mean time for the observation.
; for i=0L,nfiles-1 do begin
;     h=headfits(datadir+files[i],ext=1,/silent)
;     timestring=sxpar(h,'UTC-TIME')
;     time=strsplit(timestring,/extract,':')
;     exptime = float(sxpar(h, 'EXPTIME'))
;     utcexp[i]=ten(time[0],time[1],time[2])-((exptime/3600.)/2.)
;     parangs[i]=sxpar(h,'PARANG')
;     dts = abs(utchwp-utcexp[i])
;     min_dt = min(dts, argmin) ;Get index of nearest hwp move end (in time) to utcexp[i]
;     if utchwp[argmin] lt utcexp[i] then hwpind = argmin ;If corresponding time is before utcexp[i], then that index gives the matching hwp pos
;     if utchwp[argmin] gt utcexp[i] then hwpind = argmin-1 ;If corresponding time is after utcexp[i], then we want entry for one index prior!
;     hwp_ang[i] = angs[hwpind]
;     min_diff = min(abs(hwp_ang[i]-uni_angs), angpos)
;     hwp_pos[i] = angpos
; endfor

; ; This version assumes 'UTC-TIME' gives the end time for the observation.
if time_format eq 'fits' then begin
    for i=0L,nfiles-1 do begin
        h=headfits(datadir+files[i],ext=1,/silent)
        timestring=sxpar(h,'UTC-TIME')
        time=strsplit(timestring,/extract,':')
        exptime = float(sxpar(h, 'EXPTIME'))
        utcexp[i]=ten(time[0],time[1],time[2])-(exptime/3600.)
        parangs[i]=sxpar(h,'PARANG')
        dts = abs(utchwp-utcexp[i])
        min_dt = min(dts, argmin) ;Get index of nearest hwp move end (in time) to utcexp[i]
        if utchwp[argmin] lt utcexp[i] then hwpind = argmin ;If corresponding time is before utcexp[i], then that index gives the matching hwp pos
        if utchwp[argmin] gt utcexp[i] then hwpind = argmin-1 ;If corresponding time is after utcexp[i], then we want entry for one index prior!
        hwp_ang[i] = angs[hwpind]
        min_diff = min(abs(hwp_ang[i]-uni_angs), angpos)
        hwp_pos[i] = angpos
    endfor
endif

nangs = n_elements(uni_angs)
ang_counts = intarr(nangs)
for i=0,nangs-1 do ang_counts[i] = total((hwp_pos eq i))

ncycles = min(ang_counts, min_occ_pos)
min_occ_msk = where(hwp_pos eq min_occ_pos)
inds = indgen(nfiles)

hwp_cyc1=intarr(nfiles)
hwp_cyc1[*]=-1
used_inds=[]
for i=0L,ncycles-1 do begin
    ind_min_occ = inds[min_occ_msk[i]]
    used_inds=[used_inds, ind_min_occ]
    hwp_cyc1[ind_min_occ] = i
    dpas = abs(sin((parangs[ind_min_occ] - parangs)/180.*!PI))
    for j=0,nangs-1 do begin
        if j ne min_occ_pos then begin
            ang_inds = inds[where(hwp_pos eq j, ang_count)]
            ang_inds_notused = []; I hope there's a better way to do this in IDL...
            for k=0,ang_count-1 do begin
                used = where(used_inds eq ang_inds[k], used_count)
                if used_count eq 0 then ang_inds_notused=[ang_inds_notused, ang_inds[k]]
            endfor
            if n_elements(ang_inds_notused) gt 0 then begin
                min_dpa = min(dpas[ang_inds_notused], dpa_argmin)
                nearest_ind = ang_inds_notused[dpa_argmin]
                used_inds=[used_inds, nearest_ind]
                hwp_cyc1[nearest_ind] = i
            endif
        endif
    endfor
endfor 

hwp_cyc2=intarr(nfiles)
hwp_cyc2[*]=-1
used_inds = []
for ind=0L,ncycles-1 do begin
    i = ncycles-1-ind ; reverse order...
    ind_min_occ = inds[min_occ_msk[i]]
    used_inds=[used_inds, ind_min_occ]
    hwp_cyc2[ind_min_occ] = i
    dpas = abs(sin((parangs[ind_min_occ] - parangs)/180.*!PI))
    for j=0,nangs-1 do begin
        if j ne min_occ_pos then begin
            ang_inds = inds[where(hwp_pos eq j, ang_count)]
            ang_inds_notused = []; I hope there's a better way to do this in IDL...
            for k=0,ang_count-1 do begin
                used = where(used_inds eq ang_inds[k], used_count)
                if used_count eq 0 then ang_inds_notused=[ang_inds_notused, ang_inds[k]]
            endfor
            if n_elements(ang_inds_notused) gt 0 then begin
                min_dpa = min(dpas[ang_inds_notused], dpa_argmin)
                nearest_ind = ang_inds_notused[dpa_argmin]
                used_inds=[used_inds, nearest_ind]
                hwp_cyc2[nearest_ind] = i
            endif
        endif
    endfor
endfor 

paranges1 = []
paranges2 = []
for i=0L,ncycles-1 do begin
    cycleinds1 = where(hwp_cyc1 eq i)
    first1 = where(utcexp[cycleinds1] eq min(utcexp[cycleinds1]), complement=notfirst1)
    dpa1 = abs(((360 + parangs[cycleinds1[first1]]-parangs[cycleinds1[notfirst1]]) - 360) mod 360)
    parange1 = max(dpa1)
    paranges1 = [paranges1, parange1]

    cycleinds2 = where(hwp_cyc2 eq i)
    first2 = where(utcexp[cycleinds2] eq min(utcexp[cycleinds2]), complement=notfirst2)
    dpa2 = abs(((360 + parangs[cycleinds2[first2]]-parangs[cycleinds2[notfirst2]])-360) mod 360)
    parange2 = max(dpa2)
    paranges2 = [paranges2, parange2]
endfor

if max(paranges1) lt max(paranges2) then hwp_cyc = hwp_cyc1 else hwp_cyc = hwp_cyc2

openw,lunit,hwpfileout,/get_lun
printf, lunit, 'FILE, ', 'HWP_POS, ', 'HWP_ANG, ', 'CYCLE, '
for i=0,nfiles-1 do begin
    printf, lunit, filenum[i], hwp_pos[i], hwp_ang[i], hwp_cyc[i]
endfor
free_lun,lunit
endofprogram:
end