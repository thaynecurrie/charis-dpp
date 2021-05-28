pro charis_pdi_split_pols,pfname,prefname=prefname,suffname=suffname,polcent_offset=polcent_offset,shift_to_center=shift_to_center,help=help

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
    print,'charis_pdi_split_pols.pro: Splits CHARIS PDI cubes containing slices of two orthogonally polarized'
    print,'images to individual left and right polarized images.'
    print,"Typically executed after 'charis_subtract_sky' and before 'charis_pdi_register_cube'."
    print,'Written by K. Lawson (2020)'
    print,''
    print,'**Calling Sequence**'
    print,"charis_pdi_split_pols,pfname,polcent_offset=polcent_offset,prefname=prefname,suffname=suffname,shift_to_center=shift_to_center"
    print,''
    print,'Example:'
    print,"charis_pdi_split_pols,'abaur_low.info'"
    print,''
    print,"***Keywords***"
    print,'pfname=STRING - The parameter file (e.g. HR8799_low.info)'
    print,"prefname=STRING - The prefix of files to use as input. Defaults to 'n'."
    print,"suffname=STRING - The suffix of files to use as input. Defaults to 'e'."
    print,"polcent_offset=ARRAY - Two element array providing the pixel coordinate offset to the center of the right pol from the nominal center. Defaults to [31.758, 16.542]."
    print,"/shift_to_center - Shifts the left and right polarized images to approximate center. In general, should NOT be used."
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

;**** Prefix names for your files (Make sure to change these with different data!!!)**
if ~keyword_set(prefname) then prefname='n'
; This could be used on any of three products: 
; a) PDI flatfielded + skysubbed data
; b) Just PDI flatfielded
; c) Just sky subbed
if ~keyword_set(suffname) then begin
    test=file_search(datadir+'*flat_skysub.fits')
    if n_elements(test) gt 1 then begin
        suffname='e_flat_skysub'
    endif else begin
        test=file_search(datadir+'*flat.fits')
        if n_elements(test) gt 1 then begin
            suffname='e_flat'
        endif else begin
            test=file_search(datadir+'*skysub.fits')
            if n_elements(test) gt 1 then begin
                suffname='e_skysub'
            endif else begin
                suffname='e'
            endelse
        endelse
    endelse
endif
print,'Using prefname: ',prefname
print,'Using suffname: ',suffname

;**** Intializing unsigned offsets for polarimetric centers from nominal image center
if ~keyword_set(polcent_offset) then polcent_offset=[31.758, 16.542]
dx0_pol=polcent_offset[0]
dy0_pol=polcent_offset[1]

filenum=nbrlist(flist)
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
filesoutleft=filelist(filenum,prefix=prefname,suffix='left')
filesoutright=filelist(filenum,prefix=prefname,suffix='right')

h=headfits(datadir+files[0],/exten)
h0=headfits(datadir+files[0])

nx=sxpar(h,'naxis1')
ny=sxpar(h,'naxis2')
nL=sxpar(h,'naxis3')

x0=(nx-1.)/2.
y0=(ny-1.)/2.

;; Construct new dimensions, being careful about scalars
xvals=findgen(nx)
yvals=findgen(ny)

xg=fltarr(nx,ny)
for i=0,nx-1 do xg[*,i]=xvals
yg=fltarr(nx,ny)
for i=0,ny-1 do yg[i,*]=yvals

polmap=(dx0_pol/dy0_pol)*(xg-x0)+(yg-y0)

Lidx=where(polmap lt 0)
Lmsk=fltarr(nx,ny)
Lmsk[Lidx]=1.

Ridx=where(polmap gt 0)
Rmsk=fltarr(nx,ny)
Rmsk[Ridx]=1.

for iT=0L,nfiles-1 do begin
    ;read in file
    a=readfits(datadir+files[iT],h1,/exten,/silent)
    h0=headfits(datadir+files[iT],ext=0,/silent)
    sz=size(a)

    aleft=fltarr(nx,ny,nL)
    aright=fltarr(nx,ny,nL)
    for iL=0L,nL-1 do begin
        aleft[*,*,iL]=a[*,*,iL]*Lmsk
        aright[*,*,iL]=a[*,*,iL]*Rmsk
        if keyword_set(shift_to_center) then begin
            aleft[*,*,iL]=shift_sub(aleft[*,*,iL], dx0_pol, dy0_pol)
            aright[*,*,iL]=shift_sub(aright[*,*,iL], -1.*dx0_pol, -1.*dy0_pol)
        endif
    endfor

    writefits,reducdir+filesoutright[iT],0,h0
    writefits,reducdir+filesoutright[iT],aright,h1,/append

    writefits,reducdir+filesoutleft[iT],0,h0
    writefits,reducdir+filesoutleft[iT],aleft,h1,/append  
endfor
endofprogram:
end