pro charis_pdi_flatfield, pfname, prefname=prefname, suffname=suffname, darkflatdir=darkflatdir, brightflatdir=brightflatdir, help=help

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
    print,'charis_pdi_flatfield.pro: Computation and application of CHARIS PDI master flatfield'
    print,"If the working directory already contains 'master_pdiflat.fits', it will be used."
    print,"Otherwise, 'master_pdiflat.fits' will be computed from bright and dark PDI flats."
    print,"Typically executed after 'charis_imprep' and before 'charis_subtract_sky'."
    print,'Written by K. Lawson (2021)'
    print,''
    print,'**Calling Sequence**'
    print,"charis_pdi_flatfield, pfname, prefname=prefname, suffname=suffname, darkflatdir=darkflatdir, brightflatdir=brightflatdir, help=help"
    print,''
    print,'Example:'
    print,"charis_pdi_flatfield, 'abaur_low.info'"
    print,''
    print,"***Keywords***"
    print,'pfname=STRING - The parameter file (e.g. HR8799_low.info)'
    print,"prefname=STRING - The prefix of files to use as input. Defaults to 'n'."
    print,"suffname=STRING - The suffix of files to use as input. Defaults to 'e'."
    print,"darkflatdir=STRING - The path to the directory containing dark flats. Defaults to './darkflats/'."
    print,"brightflatdir=STRING - The path to the directory containing bright flats. Defaults to './brightflats/'."
    goto,endofprogram
endif

;********Directory Structure*****
reducdir='./reduc/'

;data directory
datadir=reducdir+'prep/'

;determine reduction subdirectory
subdir='prep/'
reducdir+=subdir
;*******end Directory Structure

;******Keywords and Filenames*****
;create list of filenames
param,'fnum_sat',flist,/get,pfname=pfname
param,'fnum_sky',flistsky,/get,pfname=pfname

if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='e'

test=file_search('./master_pdiflat.fits',count=count)
if count eq 0 then begin
    ; If the master flat doesn't already exist, create one.
    if ~keyword_set(darkflatdir) then darkflatdir='./darkflats/'
    if ~keyword_set(brightflatdir) then brightflatdir='./brightflats/'

    darkflat_files = file_search(darkflatdir+'*.fits', count=ndark)
    brightflat_files = file_search(brightflatdir+'*.fits', count=nbright)

    h0 = headfits(darkflat_files[0], ext=0, /silent)
    imtest = readfits(darkflat_files[0],h1,/exten,/silent)
    sz = size(imtest, /dim)

    dark_array = fltarr(sz[0],sz[1],sz[2],ndark)
    for i=0,ndark-1 do dark_array[*,*,*,i] = readfits(darkflat_files[i],/exten,/silent)
    darkcube = median(dark_array, dimension=4, /even)

    bright_array = fltarr(sz[0],sz[1],sz[2],nbright)
    for i=0,nbright-1 do bright_array[*,*,*,i] = readfits(brightflat_files[i],/exten,/silent)
    brightcube = median(bright_array, dimension=4, /even)

    charis_generate_roi,output=fov_mask,roirange=[125,125,-26.9166],roidim=[sz[0],sz[1]]
    fov_mask = where(fov_mask eq 1)
    flatcube = brightcube-darkcube
    for i=0,sz[2]-1 do begin
        flatslice = flatcube[*,*,i]
        flatslice = flatslice/(median(flatslice[fov_mask], /even))
        flatslice[where(flatslice eq 0)] = 1
        flatcube[*,*,i] = flatslice
    endfor
    writefits,'master_pdiflat.fits',0,h0
    writefits,'master_pdiflat.fits',flatcube,h1,/append

endif else begin
    ; If it DOES exist already, then load it.
    flatcube=readfits('master_pdiflat.fits',h1,/exten,/silent)
endelse

filenum=nbrlist(flist)
if flistsky eq '' then skyfilenum = [] else skyfilenum = nbrlist(flistsky)
filenum = [filenum,skyfilenum]
files=datadir+filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
filesout=reducdir+filelist(filenum,prefix=prefname,suffix='e_flat')

for i=0,nfiles-1 do begin
    h0=headfits(files[i],ext=0,/silent)
    imcube=readfits(files[i],h1,/exten,/silent)
    imcube_flat = imcube/flatcube
    writefits,filesout[i],0,h0
    writefits,filesout[i],imcube_flat,h1,/append
endfor

endofprogram:
end