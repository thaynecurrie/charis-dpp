pro charis_register_cube,pfname,prefname=prefname,suffname=suffname,$
method=method,$
guessoffsets=guessoffsets,$
rsub=rsub,$
refcube=refcube,$
medbox=medbox,$
ladder=ladder,$
xcorr=xcorr,$
revise=revise,$
splitpsf=splitpsf,$
keepprev=keepprev,$
nosmooth=nosmooth,$
fwhmlim=fwhmlim,$
smallsteps=smallsteps,$
;,bandaid=bandaid,$
;refslice=refslice,$
smask=smask,$
checkquality=checkquality,$
astrogrid=astrogrid,$
verbose=verbose,help=help

;05/24/2021 - Edited to add 'astrogrid' keyword: if using PDI mode, set to 11.2
;03/27/2020 - Edited to handle cal spot frame used to register images without spots
;03/25/2020 - Includes switch to register ladder sequences that select different files and do not produce a reduc.log file
;07/26/2019 - Changed the Reference Slice to Channel 11; ***CAUTION*** ONLY GOOD FOR BROADBAND RIGHT NOW
;02/16/2018 - Heavily revised and simplified.   Registration organized into three methods: sat spots, psf centroiding, halo centroiding (to do!!!)
;09/27/2017 - Switch to estimate spot positions in registered images where spot positions cannot be calculated
;09/07/2017 - Switch to escape the spot placements to -9999 in case of failure.  Allows to you to see whether some frames are just off
;	or whether there was a centroid adjustment later.
;08/13/2017 - Option to find satellite spot positions from x-corr of region of interests instead of gaussian centroiding.
;07/10/2017 - Renamed to charis_imregister.pro
;07/10/2017 - Now at least for the sat-spot = true case, saves the newsatspot positions in fits header.
;04/13/2017 - Added switch to do x-corr fitting of the halo instead of satellite spots
;04/08/2017 - Redone for CHARIS to find the satellite positions
;03/09/2014 - Quick hack-version for GPI image registration across cube
; components.  

;Limitations ... 
;- cross-correlation currently does not work
;- gaussian centroiding ONLY for method = 'psf'
;- does not yet do method = 'halo'

;*****Nominal Program (no switches)*****
;1. Uses an initial guess for the header positions from the reference slice that is hardwired for now.
;2. Computes centroid position for this reference slice.  
;3. Goes on to compute satellite spot positions at other wavelengths.
;4. Computes the centroid from the spot positions in each slice
;5. Fits a 2rd-order polynominal to the centroid vs. slice
;6. Shifts image slice to put the centroid at dimx/2 dimy/2

;Result: puts the star at the dimx/2 dimy/2 position. 

;*****Method****
;method = [sats, psf, calsat]
;****notes: halo [matching just on the basis of the PSF halo] is TBD but lower priority

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"charis_register_cube,pfname,method=method,guessoffsets=guessoffsets,astrogrid=astrogrid,rsub=rsub,medbox=medbox"
print,",refcube=refcube,ladder=ladder,xcorr=xcorr,splitpsf=splitpsf,keepprev=keepprev,nosmooth=nosmooth,revise=revise,fwhmlim=fwhmlim"
print,"smallsteps=smallsteps,smask=smask,checkquality=checkquality,verbose=verbose,help=help"
print,""
print,"***Major Keywords***"
print,""
print,'*pfname - parameter file (e.g. HR8799_low.info)'
print,""
print,"*method - [which method used for registering cubes?]"
print," - choices: "
print," 0/'sats' [use sat spots], Examples: charis_register_cube,'HR8799_low.info',method=0"
print,"OR","                                  charis_register_cube,'HR8799_low.info',method='sats'"
print,""
print," 1/'psf' [use unsat. PSF], Examples: charis_register_cube,'LkCa15_low.info',method=1"
print,"OR","                                  charis_register_cube,'LkCa15_low.info',method='psf'"
print,""
print," 2/'calsat' [use sat spots in cal. cube to register cubes w/o spots]"
print,"                          Examples: charis_register_cube,'ABAur_low.info',method=2"
print,"OR","                                  charis_register_cube,'ABAur_low.info',method='calsat'"
print,""
print,"*guessoffsets - guess the starting centroid offset from array center [default=0]"
print,"*rsub - subtract radial profile b4 centroiding? [default=no]; *medbox - box-width for unsharp masking [default=15]"
print,""
print,"***other keywords"
print,"*refcube [reference cube for centroiding, method=1],*ladder [is this a ladder sequence?],*xcorr [are we doing cross-correlating?, method=0]"
print,"*revise [in case of bad centroiding],*splitpsf [fix for a splitting PSF core], *keepprev [keeps orig. centroid estimate incaseof failure]"
print,"*nosmooth [do not do sub-PSF core smoothing b4 centroiding],*fwhmlim[boxsize for centroid search]"
print,"*smallsteps[for method=0+xcorr flag thrown,*smask[mask pixels surrounding bright companion]"
print,"*checkquality[calculate core-to-halo ratio, plot results],*verbose[for debugging]"
print,"*astrogrid [Alternate astrogrid spacing in lambda / D, to use for spacings other than 15.9 lambda/D; 15.9 lambda/D for total intensity, usually 11.2 lambda/D for PDI]"
goto,skiptotheend
endif

if ~keyword_set(method) then begin 
method = 0

goto,skipmethodlist
endif

if method eq 'sats' then method = 0
if method eq 'psf' then method = 1
if method eq 'calsats' then method =2
skipmethodlist:

;if ~keyword_set(refslice) then refslice=10  ;hardwire this for now

refslice=10 ;hardwire this for now

if ~keyword_set(ladder) then begin
;*****Output Log File
;***Make file saving the centroid positions (just for accounting purposes), the inner saturation radius (don't worry),
;the hour angle (need it for LOCI/ADI), and the Parallactic Angle (*Definitely* need for LOCI/ADI)
;Writing out the reduc.log file header and first entry
openw,1,'reduc.log'
printf,1,'File_no,','XC,','YC,','Rsat','HA','Par.Angle'
;*****
endif

;********Directory Structure*****
; setupdir,reducdir=reducdir
reducdir='./reduc/'

;data directory
datadir=reducdir+'prep/'

;determine reduction subdirectory
subdir='reg/'
reducdir+=subdir

;define a temporary directory
;tmpdir=reducdir+'tmp/'
;file_mkdir,tmpdir
;*******end Directory Structure


;******Keywords and Filenames*****
;create list of filenames
param,'obsdate',date,/get,pfname=pfname & date=strtrim(date,2)
if ~keyword_set(ladder) then begin
param,'fnum_sat',flist,/get,pfname=pfname
endif else begin
param,'fnum_lad',flist,/get,pfname=pfname
endelse
param,'rsat',rsat,/get,pfname=pfname
param,'DEC',decdeg,/get,pfname=pfname
lat = 19.825d0 ;mauna kea

if ~keyword_set(fwhmlim) then fwhmlim=7

;**** Prefix names for your files (Make sure to change these with different data!!!)**
if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then begin
;if you can find sky subtracted frames, then use these.  If not, then use non-subbed ones
test=file_search(datadir+'*skysub.fits')
if n_elements(test) gt 1 then begin
suffname='e_skysub'
endif else begin
suffname='e'
endelse
endif

filenum=nbrlist(flist)
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
filesout=filelist(filenum,prefix=prefname,suffix='reg')

;Use spot position from the 'ladder' cube which has satellite spots present; cross-correlate the halo to register all cubes
if (method eq 2) then begin
param,'fnum_lad',flistcalsat,/get,pfname=pfname
filenumcal=nbrlist(flistcalsat)
filescalsat=filelist(filenumcal,nfilescal,prefix=prefname,suffix=suffname)
if nfilescal gt 1 then begin
reffile=filescalsat[0]
endif else begin
reffile=filescalsat
endelse

endif

;*************

paranglef=fltarr(nfiles)

;***To be considered in a future release
param,'spang*',spang,/get,pfname=pfname
param,'cenmask',cenmask,/get,pfname=pfname
param,'spmask',spmask,/get,pfname=pfname
param,'RA',ra2000,/get,pfname=pfname
cenmask=strsplit(cenmask,',',/extract)
cenmasks=cenmask
;***

if ~keyword_set(medbox) then medbox=15
;******end Keywords and Filenames*****

;***source mask file****
if keyword_set(smask) then begin
readcol,'smask.dat',xsource,ysource
endif

;get dim from first image header
h=headfits(datadir+files[0],/exten)
h0=headfits(datadir+files[0])
print,files[0]
dim=sxpar(h,'naxis1')
if ~keyword_set(xc0) then xc0=dim/2
if ~keyword_set(yc0) then yc0=dim/2


;You probably want to set for array element 0 for now.
if ~keyword_set(refcube) then refcube=0

;******Image Registration Methods******
;1. Use Satellite Spots  ; Limitation: x-correlation currently doesn't work
;2. Use the Unsaturated Star PSF ; Limitation: x-correlation currently doesn't work
;3. Use the Halo  ; ---- TO DO!!!

case method of 
 0: begin
 ;'sats': begin

;******Method = 'sats'

;take initial frame and grab the filter name
filtname=strtrim(sxpar(h0,'CAL_BAND'))
;filtnamelist=['
case filtname of 
'J': begin
       sat_xguess=charis_get_constant(name='jspots_x')
       sat_yguess=charis_get_constant(name='jspots_y')
     end
'H': begin
       sat_xguess=charis_get_constant(name='hspots_x')
       sat_yguess=charis_get_constant(name='hspots_y')
       ;sat_xguess=[64,112,137,89]
       ;sat_yguess=[114,138,89,64]
       ;refslice=0; until you fix above
     end
'K': begin
       sat_xguess=charis_get_constant(name='kspots_x')
       sat_yguess=charis_get_constant(name='kspots_y')
       ;sat_xguess=[50,116,150,84]
       ;sat_yguess=[117,150,83,49]
       ;refslice=0; until you fix above
     end
'lowres': begin
       sat_xguess=charis_get_constant(name='lowspots_x')
       sat_yguess=charis_get_constant(name='lowspots_y')
     end
  
 else: begin

;endif else begin
;if not highres then default to low res spots for now
       sat_xguess=charis_get_constant(name='lowspots_x')
       sat_yguess=charis_get_constant(name='lowspots_y')
       ;  sat_xguess=[59,113,141,87]
       ;  sat_yguess=[114,141,86,59]
       end
endcase

if keyword_set(astrogrid) then begin
  sat_xguess = (sat_xguess-xc0)*(astrogrid/15.9)+xc0
  sat_yguess = (sat_yguess-yc0)*(astrogrid/15.9)+yc0
endif

if keyword_set(guessoffsets) then begin
sat_xguess+=guessoffsets[0]
sat_yguess+=guessoffsets[1]
endif

;print,sat_xguess
;print,sat_yguess
;stop


for i=0L,nfiles-1 do begin

print,'Registering File Number ',i+1,' ',files[i]
;read in file
a=readfits(datadir+files[i],h1,/exten,/silent)
h0=headfits(datadir+files[i],ext=0,/silent)
a_in=a
sz=size(a)


parangle0=float(sxpar(h0,'PA'))
ha0=float(sxpar(h0,'HA'))

;***This *should* work.  Manually edit if it results in PA being discontinuous in reduc.log

;center the pa value on zero
;need to fix this...
;if parangle0 gt 180 then parangle0-=360
if parangle0 gt 180 and decdeg lt lat then parangle0-=360

;***

;*****Image Registration Loop.
;**- If you have satellite spots then you have a choice of getting the centroid from gaussian-fitting or from joint cross-correlation

get_charis_wvlh,h0,wvlhs
scl=wvlhs[refslice]/wvlhs


;****do an initial centroid measurement
;***assumes that the reference channel sat spots are visible.

a_firstchan=a[*,*,refslice]
a_firstchan-=filter_image(a_firstchan,median=medbox)
if keyword_set(rsub) then begin
profrad_tc,a_firstchan,2,1,dim/2,p2d=pr
;profrad_tc,imslice2,2,allrsat[i],rmax,p2d=pr,p1d=p1d,rayon=r
a_firstchan-=pr
endif
 gcntrd,a_firstchan,sat_xguess[0],sat_yguess[0],satcx0,satcy0,fwhmlim
 gcntrd,a_firstchan,sat_xguess[1],sat_yguess[1],satcx1,satcy1,fwhmlim
 gcntrd,a_firstchan,sat_xguess[2],sat_yguess[2],satcx2,satcy2,fwhmlim
 gcntrd,a_firstchan,sat_xguess[3],sat_yguess[3],satcx3,satcy3,fwhmlim

 
if (satcx0 lt 0 or satcy0 lt 0) then begin 
cntrd,a_firstchan,sat_xguess[0],sat_yguess[0],satcx0,satcy0,fwhmlim,/keepcenter
endif
if (satcx1 lt 0 or satcy1 lt 0) then begin 
cntrd,a_firstchan,sat_xguess[1],sat_yguess[1],satcx1,satcy1,fwhmlim,/keepcenter
endif
if (satcx2 eq -1 or satcy2 eq -1) then begin 
cntrd,a_firstchan,sat_xguess[2],sat_yguess[2],satcx2,satcy2,fwhmlim,/keepcenter
endif
if (satcx3 eq -1 or satcy3 eq -1) then begin 
cntrd,a_firstchan,sat_xguess[3],sat_yguess[3],satcx3,satcy3,fwhmlim,/keepcenter
endif

if (satcx0 lt 0 or satcy0 lt 0) then begin
cntrd,a_firstchan,sat_xguess[0],sat_yguess[0],satcx0,satcy0,fwhmlim,extendbox=2
endif
if (satcx1 lt 0 or satcy1 lt 0) then begin
cntrd,a_firstchan,sat_xguess[1],sat_yguess[1],satcx1,satcy1,fwhmlim,extendbox=2
endif
if (satcx2 eq -1 or satcy2 eq -1) then begin
cntrd,a_firstchan,sat_xguess[2],sat_yguess[2],satcx2,satcy2,fwhmlim,extendbox=2
endif
if (satcx3 eq -1 or satcy3 eq -1) then begin
cntrd,a_firstchan,sat_xguess[3],sat_yguess[3],satcx3,satcy3,fwhmlim,extendbox=2
endif

 cens0=[[satcx0,satcy0],[satcx1,satcy1],[satcx2,satcy2],[satcx3,satcy3]]
 
 cens=dblarr(2,4,sz[3])
 cens[*,*,0]=cens0
 c0=(total(cens0,2)/4) # (fltarr(4)+1.)
 cens0p=cv_coord(from_rect=cens0-c0,/to_polar)
initial_centroid=median(cens0,dimension=2,/even)

;dist_circle,g,[sz[1],sz[2]]

PSFcens=fltarr(2,sz[3])
deltax=fltarr(sz[3])
deltay=fltarr(sz[3])

for ii=0L,sz[3]-1 do begin
 cens[*,*,ii]=cv_coord(from_polar=[cens0p[0,*],cens0p[1,*]/scl[ii]],/to_rect)+c0

 a_test=a[*,*,ii]
;medbox=15
 a_test-=filter_image(a_test,median=medbox)  ; this is roughly 5 lambda/D = safe
if keyword_set(rsub) then begin
profrad_tc,a_test,2,1,dim/2,p2d=pr
a_test-=pr
endif
if ~keyword_set(nosmooth) then a_test=smooth(a_test,3)

;****if you have a very bright pt source that is skewing the sat spot determination, apply an NAN mask here.
if keyword_set(smask) then begin
for sm=0L,n_elements(x)-1 do begin
gcntrd,a_test,xsource[sm],ysource[sm],xsourceout,ysourceout
dist_circle,sloc,[sz[1],sz[2]],xsourceout,ysourceout
maskme=where(sloc le 2)
a_test[maskme]=!values.f_nan
endfor
endif

;***
;****1. Gaussian Centroiding, 2. Cross-Correlation Fitting to find the star center.
;do gaussian centroid estimate for satellite spot positions.  

if ~keyword_set(xcorr) then begin
;*Gaussian Centroiding

 gcntrd,a_test,cens[0,0,ii],cens[1,0,ii],satcx0,satcy0,fwhmlim
 gcntrd,a_test,cens[0,1,ii],cens[1,1,ii],satcx1,satcy1,fwhmlim
 gcntrd,a_test,cens[0,2,ii],cens[1,2,ii],satcx2,satcy2,fwhmlim
 gcntrd,a_test,cens[0,3,ii],cens[1,3,ii],satcx3,satcy3,fwhmlim

if (satcx0 lt 0 or satcy0 lt 0) then begin 
cntrd,a_test,cens[0,0,ii],cens[1,0,ii],satcx0,satcy0,fwhmlim,/keepcenter
endif
if (satcx1 lt 0 or satcy1 lt 0) then begin 
;print,'v1',satcx1,satcy1
cntrd,a_test,cens[0,1,ii],cens[1,1,ii],satcx1,satcy1,fwhmlim,/keepcenter
;print,'v2',satcx1,satcy1
endif
if (satcx2 lt 0 or satcy2 lt 0) then begin 
cntrd,a_test,cens[0,2,ii],cens[1,2,ii],satcx2,satcy2,fwhmlim,/keepcenter
endif
if (satcx3 lt 0 or satcy3 lt 0) then begin 
cntrd,a_test,cens[0,3,ii],cens[1,3,ii],satcx3,satcy3,fwhmlim,/keepcenter
endif

if keyword_set(revise) then begin

if (satcx0 lt 0 or satcy0 lt 0) then begin
cntrd,a_test,cens[0,0,ii],cens[1,0,ii],satcx0,satcy0,fwhmlim,extendbox=2
endif
if (satcx1 lt 0 or satcy1 lt 0) then begin
;print,'v1',satcx1,satcy1
cntrd,a_test,cens[0,1,ii],cens[1,1,ii],satcx1,satcy1,fwhmlim,extendbox=2
print,'v3',satcx1,satcy1
;stop
endif
if (satcx2 lt 0 or satcy2 lt 0) then begin
cntrd,a_test,cens[0,2,ii],cens[1,2,ii],satcx2,satcy2,fwhmlim,extendbox=2
endif
if (satcx3 lt 0 or satcy3 lt 0) then begin
cntrd,a_test,cens[0,3,ii],cens[1,3,ii],satcx3,satcy3,fwhmlim,extendbox=2
endif

endif
 cens[0,0,ii]=satcx0
 cens[1,0,ii]=satcy0
 cens[0,1,ii]=satcx1
 cens[1,1,ii]=satcy1
 cens[0,2,ii]=satcx2
 cens[1,2,ii]=satcy2
 cens[0,3,ii]=satcx3
 cens[1,3,ii]=satcy3

;kill switch
if (satcx0 lt 0 or satcy0 lt 0 or satcx1 lt 0 or satcy1 lt 0 or satcx2 lt 0 or satcy2 lt 0 or satcx3 lt 0 or satcy3 lt 0) then begin
cens[*,*,ii]=!values.f_nan
endif

PSFcens[*,ii]=[(median(cens[0,*,ii],/even)),(median(cens[1,*,ii],/even))]
deltax[ii]=xc0-PSFcens[0,ii]
deltay[ii]=yc0-PSFcens[1,ii]

if keyword_set(verbose) then print,'Center from Gaussian Fitting is ',PSFcens[*,ii],' at slice',ii

endif else begin
;***Cross-Correlation Fitting


;calculating satellite spot locations
;for s=0L,sz[3]-1 do begin

;*Define region of interest
;sat 1
dist_circle,sat1,[sz[1],sz[2]],cens[0,0,ii],cens[1,0,ii]

;sat 2
dist_circle,sat2,[sz[1],sz[2]],cens[0,1,ii],cens[1,1,ii]

;sat 3
dist_circle,sat3,[sz[1],sz[2]],cens[0,2,ii],cens[1,2,ii]

;sat 4
dist_circle,sat4,[sz[1],sz[2]],cens[0,3,ii],cens[1,3,ii]

search_radius=2.
roi=where((sat1 le search_radius) or (sat2 le search_radius) or (sat3 le search_radius) or (sat4 le search_radius),nroi)

print,initial_centroid

if ~keyword_set(smallsteps) then begin
fidcenter=charis_getcenter(a_test,initial_centroid[0],initial_centroid[1],roi)

endif else begin
fidcenter=charis_getcenter(a_test,initial_centroid[0],initial_centroid[1],roi,/smallsteps)
endelse

PSFcens[0,ii]=fidcenter[0]
PSFcens[1,ii]=fidcenter[1]

deltax[ii]=xc0-PSFcens[0,ii]
deltay[ii]=yc0-PSFcens[1,ii]
if keyword_set(verbose) then print,'Center from X-Corr Fitting is ',PSFcens[*,ii],' at slice',ii


endelse

endfor

;****Polynomial fit to individual slice centroid measurements****
slices=findgen(sz[3])
;-refslice
;good=where(finite(PSFcens[0,*]) ne 0 and finite(PSFcens[1,*]) ne 0,ngood)
good=where(finite(PSFcens[0,2:sz[3]-1-2]) ne 0 and finite(PSFcens[1,2:sz[3]-1-2]) ne 0,ngood)
print,'number of slices for fit ', ngood

;****for now put in break to exit in case the program cannot find any centroids
if ngood gt 0 then begin
fitx=robust_poly_fit((slices[2:sz[3]-1-2])[good]-refslice,(PSFcens[0,2:sz[3]-1-2])[good]-PSFcens[0,refslice],2,numit=15)
fity=robust_poly_fit((slices[2:sz[3]-1-2])[good]-refslice,(PSFcens[1,2:sz[3]-1-2])[good]-PSFcens[1,refslice],2,numit=15)

refoffsetx=PSFcens[0,refslice]-xc0
refoffsety=PSFcens[1,refslice]-yc0

deltax=fitx[0]+fitx[1]*(slices-refslice)+fitx[2]*(slices-refslice)^2.+refoffsetx
deltay=fity[0]+fity[1]*(slices-refslice)+fity[2]*(slices-refslice)^2.+refoffsety

print,refoffsetx,xc0,PSFcens[0,refslice]

if keyword_set(guessoffsets) then begin
plot,slices,PSFcens[0,*],yrange=[sz[1]/2-7+long(guessoffsets[0]),sz[1]/2+7+long(guessoffsets[1])],ystyle=1,/nodata
endif else begin
plot,slices,PSFcens[0,*],yrange=[sz[1]/2-7,sz[1]/2+7],ystyle=1,/nodata
endelse
oplot,slices,PSFcens[0,*],psym=4,thick=5
oplot,slices,deltax+xc0,linestyle=1
oplot,slices,PSFcens[1,*]+3,psym=4,thick=5
oplot,slices,deltay+yc0+3,linestyle=1

endif else begin
deltax[*]=-9999
deltay[*]=-9999
refoffsetx=-9999
refoffsety=-9999
endelse


;if you cannot find the sat spot position in some slice then use reference channel sat position and extrapolate
; cens[*,*,ii]=cv_coord(from_polar=[cens0p[0,*],cens0p[1,*]/scl[ii]],/to_rect)+c0
cenxref=cens[0,*,0]-deltax[0]
cenyref=cens[1,*,0]-deltay[0]

for s=0L,sz[3]-1 do begin
a_in[*,*,s]=shift_sub(a_in[*,*,s],-1*deltax[s],-1*deltay[s])
;***new stuff
cens[0,*,s]-=deltax[s]
cens[1,*,s]-=deltay[s]
PSFcens[*,s]=[(median(cens[0,*,s],/even)),(median(cens[1,*,s],/even))]

;switch if sat spot position is NaN ...
centotal=total(cens[*,*,s])
cenref=(total(cens[*,*,0],2)/4) # (fltarr(4)+1.)

if finite(centotal) eq 0 then begin
cens0=cens[*,*,0]
 cens0p=cv_coord(from_rect=cens0-cenref,/to_polar)
 cens[*,*,s]=cv_coord(from_polar=[cens0p[0,*],cens0p[1,*]/scl[s]],/to_rect)+cenref
PSFcens[*,s]=[(median(cens[0,*,s],/even)),(median(cens[1,*,s],/even))]
endif

 for j=0L,3 do begin

  sxaddpar,h1,'SATS'+strtrim(s,2)+'_'+strtrim(j,2),$
   string(strtrim(cens[*,j,s],2),format='(F7.3," ",F7.3)'),$
    'Location of sat. spot '+strtrim(j,2)+' of slice '+strtrim(s,2)
 endfor

sxaddpar,h1,'PSFC_'+strtrim(s,2),$
 string(strtrim(PSFcens[*,s],2),format='(F7.3," ",F7.3)'),$
  'PSF Center of slice '+strtrim(s,2)

sxaddpar,h1,'shift_'+strtrim(s,2),$
 string(strtrim([-1*deltax[s],-1*deltay[s]],2),format='(F7.3," ",F7.3)'),$
  'Image shift for slice '+strtrim(s,2)

if keyword_set(verbose) then print,'PSFCENS is ',i,PSFcens[*,s],' shifts are ',deltax[s],deltay[s]

endfor



psfcenx=median(PSFcens[0,*],/even)
psfceny=median(PSFcens[1,*],/even)
sxaddpar,h1,"PSFCENTX",psfcenx,'Mean PSF center X' 
sxaddpar,h1,"PSFCENTY",psfceny,'Mean PSF center Y' 
if keyword_set(verbose) then print,'Median PSF Center is ',psfcenx,psfceny
a=a_in


writefits,reducdir+filesout[i],0,h0
writefits,reducdir+filesout[i],a,h1,/append
if ~keyword_set(ladder) then printf,1,filenum[i],psfcenx,psfceny,rsat,ha0,parangle0

endfor
        end
;******end method - sats

;******Method = 'psf'
;here, take the first cube in the sequence, copy to a dummy file, get its relative offset

     1: begin
 ;'psf': begin

;*****ref cube: not used for registration right now!*****
;reference data cube
        refimage=readfits(datadir+files[refcube],/exten,h1ref,/silent)
        h0ref=headfits(datadir+files[refcube],/silent)

        ;figure out the registration for this cube using cross-correlation
        refsz=size(refimage)
;********

         ;wavelength array
         get_charis_wvlh,h0ref,wvlhs
         scl=wvlhs[refslice]/wvlhs

         ;circular region centered on middle of image
         if ~keyword_set(guessoffsets) then begin
           guessoffsets=fltarr(2)
           guessoffsets[0]=0 
           guessoffsets[1]=0
         endif
         xoo=refsz[1]/2+guessoffsets[0] & yoo = refsz[2]/2+guessoffsets[1]

         ;dist_circle,g,[refsz[1],refsz[2]],xoo,yoo

;Telescope Diameter for Subaru
;Dtel=7.9d0 ;visible pupil for SCExAO
Dtel=charis_get_constant(name='Dtel')
;pixel scale 16.4 mas/pixel
pixscale=charis_get_constant(name='pixscale')
;pixscale=0.0164
fwhm=1.*(1.d-6*wvlhs*1d-3/Dtel)*(180.*3600./!dpi)/pixscale
         
         ;Define variables
         PSFcensref=fltarr(2,refsz[3])
         deltax=fltarr(refsz[3])
         deltay=fltarr(refsz[3])

;now Register Images

       
;slices=findgen(refsz[3]-refslice)
slices=findgen(refsz[3])
for i=0L,nfiles-1 do begin


print,'Registering File Number ',i+1,' ',files[i]
;read in file
a=readfits(datadir+files[i],h1,/exten,/silent)
h0=headfits(datadir+files[i],ext=0,/silent)
a_in=a
sz=size(a)


         PSFcens=fltarr(2,sz[3])
         deltax=fltarr(sz[3])
         deltay=fltarr(sz[3])

parangle0=float(sxpar(h0,'PA'))
ha0=float(sxpar(h0,'HA'))

;center the pa value on zero
if parangle0 gt 180 and decdeg lt lat then parangle0-=360


for il=0L,refsz[3]-1 do begin

gcntrd,a[*,*,il],xc0+guessoffsets[0],yc0+guessoffsets[1],xc1,yc1,fwhmlim

if (xc1 lt 0 or yc1 lt 0) then begin
cntrd,a[*,*,il],xc0+guessoffsets[0],yc0+guessoffsets[1],xc1,yc1,fwhmlim,/keepcenter
endif

if keyword_set(revise) then begin

if (xc1 lt 0 or yc1 lt 0) then begin
cntrd,a[*,*,il],xc0+guessoffsets[0],yc0+guessoffsets[1],xc1,yc1,fwhmlim,extendbox=3
endif

if keyword_set(splitpsf) then begin

az=sumaper_im(smooth(a[*,*,il],2),fwhmlim,30,/nan)
;cntrd,az,xc0+guessoffsets[0],yc0+guessoffsets[1],xc1,yc1,fwhm[il],extendbox=3
cntrd,az,xc0+guessoffsets[0],yc0+guessoffsets[1],xc1,yc1,fwhmlim
endif

if keyword_set(keepprev) then begin
;xc1=!values.f_nan
;yc1=!values.f_nan
xc1=PSFcens[0,il]
yc1=PSFcens[1,il]
endif

endif

PSFcens[0,il]=xc1
PSFcens[1,il]=yc1

if (xc1 lt 0 or yc1 lt 0) then begin
PSFcens[*,il]=!values.f_nan
endif

deltax[il]=xc0-PSFcens[0,il]
deltay[il]=yc0-PSFcens[1,il]

print,'Center from PSF Gaussian Fitting is ',PSFcens[*,il],' at slice',il
endfor

;*****modeling the centriods vs. wavelength as 2nd order polynomial
refoffsetx=PSFcens[0,refslice]-xc0
refoffsety=PSFcens[1,refslice]-yc0
;slices=findgen(refsz[3]-refslice)
slices=findgen(refsz[3])

good=where(finite(PSFcens[0,2:refsz[3]-1-2]) ne 0 and finite(PSFcens[1,2:refsz[3]-1-2]) ne 0,ngood)

if ngood gt 0 then begin
fitx=robust_poly_fit((slices[2:refsz[3]-1-2])[good]-refslice,(PSFcens[0,2:refsz[3]-1-2])[good]-PSFcens[0,refslice],2,numit=15)
fity=robust_poly_fit((slices[2:refsz[3]-1-2])[good]-refslice,(PSFcens[1,2:refsz[3]-1-2])[good]-PSFcens[1,refslice],2,numit=15)
 
deltax=fitx[0]+fitx[1]*(slices-refslice)+fitx[2]*(slices-refslice)^2.+refoffsetx
deltay=fity[0]+fity[1]*(slices-refslice)+fity[2]*(slices-refslice)^2.+refoffsety

if keyword_set(guessoffsets) then begin
plot,slices,PSFcens[0,*],yrange=[sz[1]/2-7+long(guessoffsets[0]),sz[1]/2+7+long(guessoffsets[1])],ystyle=1,/nodata
endif else begin
plot,slices,PSFcens[0,*],yrange=[sz[1]/2-7,sz[1]/2+7],ystyle=1,/nodata
endelse
oplot,slices,PSFcens[0,*],psym=4,thick=5
oplot,slices,deltax+xc0,linestyle=1
oplot,slices,PSFcens[1,*]+3,psym=4,thick=5
oplot,slices,deltay+yc0+3,linestyle=1


;wait,2.5
endif else begin
deltax[*]=-9999
deltay[*]=-9999
refoffsetx=-9999
refoffsety=-9999
endelse
;*********

;loop to realign
for ii=0L,sz[3]-1 do begin
a_in[*,*,ii]=shift_sub(a_in[*,*,ii],-1*deltax[ii],-1*deltay[ii])
PSFcens[0,ii]-=deltax[ii]
PSFcens[1,ii]-=deltay[ii]
sxaddpar,h1,'PSFC_'+strtrim(ii,2),$
 string(strtrim(PSFcens[*,ii],2),format='(F7.3," ",F7.3)'),$
  'PSF Center of slice '+strtrim(ii,2)
endfor


psfcenx=median(PSFcens[0,*],/even)
psfceny=median(PSFcens[1,*],/even)
sxaddpar,h1,"PSFCENTX",psfcenx,'Mean PSF center X'
sxaddpar,h1,"PSFCENTY",psfceny,'Mean PSF center Y'
a=a_in

writefits,reducdir+filesout[i],0,h0
writefits,reducdir+filesout[i],a,h1,/append
if ~keyword_set(ladder) then printf,1,filenum[i],psfcenx,psfceny,rsat,ha0,parangle0

endfor


        end

;******end method - psf


;******Method = 'calsat'
;here, take the first cube in the sequence, copy to a dummy file, get its relative offset
  2: begin
 ;'calsat': begin
 
        filtname=strtrim(sxpar(h0,'CAL_BAND'))
;filtnamelist=['
case filtname of
'J': begin
       sat_xguess=charis_get_constant(name='jspots_x')
       sat_yguess=charis_get_constant(name='jspots_y')
     end
'H': begin
       sat_xguess=charis_get_constant(name='hspots_x')
       sat_yguess=charis_get_constant(name='hspots_y')
       ;refslice=0; until you fix above
     end
'K': begin
       sat_xguess=charis_get_constant(name='kspots_x')
       sat_yguess=charis_get_constant(name='kspots_y')
       ;refslice=0; until you fix above
     end
'lowres': begin
       sat_xguess=charis_get_constant(name='lowspots_x')
       sat_yguess=charis_get_constant(name='lowspots_y')
     end

 else: begin

;endif else begin
;if not highres then default to low res spots for now
       sat_xguess=charis_get_constant(name='lowspots_x')
       sat_yguess=charis_get_constant(name='lowspots_y')
       end
endcase

if keyword_set(guessoffsets) then begin
sat_xguess+=guessoffsets[0]
sat_yguess+=guessoffsets[1]
endif


;*****ref cube: not used for registration right now!*****
;reference data cube
        refimage=readfits(datadir+reffile,/exten,h1ref,/silent)
        h0ref=headfits(datadir+reffile,/silent)

        ;figure out the registration for this cube using cross-correlation
        sz=size(refimage)
;********

         ;wavelength array
         get_charis_wvlh,h0ref,wvlhs
         scl=wvlhs[refslice]/wvlhs

;register the sat spots for the refcal

;****do an initial centroid measurement from channel 1
;***assumes that the reference channel sat spots are visible.

ref_firstchan=refimage[*,*,refslice]
ref_firstchan-=filter_image(ref_firstchan,median=medbox)
if keyword_set(rsub) then begin
profrad_tc,ref_firstchan,2,1,dim/2,p2d=pr
ref_firstchan-=pr
endif
 gcntrd,ref_firstchan,sat_xguess[0],sat_yguess[0],satcx0,satcy0,fwhmlim
 gcntrd,ref_firstchan,sat_xguess[1],sat_yguess[1],satcx1,satcy1,fwhmlim
 gcntrd,ref_firstchan,sat_xguess[2],sat_yguess[2],satcx2,satcy2,fwhmlim
 gcntrd,ref_firstchan,sat_xguess[3],sat_yguess[3],satcx3,satcy3,fwhmlim

 cens0=[[satcx0,satcy0],[satcx1,satcy1],[satcx2,satcy2],[satcx3,satcy3]]

 cens=dblarr(2,4,sz[3])
 cens[*,*,0]=cens0
 c0=(total(cens0,2)/4) # (fltarr(4)+1.)
 cens0p=cv_coord(from_rect=cens0-c0,/to_polar)
initial_centroid=median(cens0,dimension=2,/even)

dist_circle,g,[sz[1],sz[2]]

PSFcensref=fltarr(2,sz[3])
PSFcens=fltarr(2,sz[3])
deltax=fltarr(sz[3])
deltay=fltarr(sz[3])

for ii=0L,sz[3]-1 do begin
 cens[*,*,ii]=cv_coord(from_polar=[cens0p[0,*],cens0p[1,*]/scl[ii]],/to_rect)+c0

a_test=refimage[*,*,ii]

 gcntrd,a_test,cens[0,0,ii],cens[1,0,ii],satcx0,satcy0,fwhmlim
 gcntrd,a_test,cens[0,1,ii],cens[1,1,ii],satcx1,satcy1,fwhmlim
 gcntrd,a_test,cens[0,2,ii],cens[1,2,ii],satcx2,satcy2,fwhmlim
 gcntrd,a_test,cens[0,3,ii],cens[1,3,ii],satcx3,satcy3,fwhmlim

  cens[0,0,ii]=satcx0
 cens[1,0,ii]=satcy0
 cens[0,1,ii]=satcx1
 cens[1,1,ii]=satcy1
 cens[0,2,ii]=satcx2
 cens[1,2,ii]=satcy2
 cens[0,3,ii]=satcx3
 cens[1,3,ii]=satcy3

PSFcensref[*,ii]=[(median(cens[0,*,ii],/even)),(median(cens[1,*,ii],/even))]
deltax[ii]=xc0-PSFcensref[0,ii]
deltay[ii]=yc0-PSFcensref[1,ii]

endfor

;****Polynomial fit to individual slice centroid measurements****
slices=findgen(sz[3])
;-refslice
good=where(finite(PSFcensref[0,2:sz[3]-1-2]) ne 0 and finite(PSFcensref[1,2:sz[3]-1-2]) ne 0,ngood)
print,'good is', ngood

if ngood eq 0 then stop

fitx=robust_poly_fit((slices[2:sz[3]-1-2])[good]-refslice,(PSFcensref[0,2:sz[3]-1-2])[good]-PSFcensref[0,refslice],2,numit=15)
fity=robust_poly_fit((slices[2:sz[3]-1-2])[good]-refslice,(PSFcensref[1,2:sz[3]-1-2])[good]-PSFcensref[1,refslice],2,numit=15)

refoffsetx=PSFcensref[0,refslice]-xc0
refoffsety=PSFcensref[1,refslice]-yc0

deltaxref=fitx[0]+fitx[1]*(slices-refslice)+fitx[2]*(slices-refslice)^2.+refoffsetx
deltayref=fity[0]+fity[1]*(slices-refslice)+fity[2]*(slices-refslice)^2.+refoffsety

print,'Offsets of Reference Image ',deltaxref,deltayref


if keyword_set(guessoffsets) then begin
plot,slices,PSFcensref[0,*],yrange=[sz[1]/2-7+long(guessoffsets[0]),sz[1]/2+7+long(guessoffsets[1])],ystyle=1,/nodata
endif else begin
plot,slices,PSFcensref[0,*],yrange=[sz[1]/2-7,sz[1]/2+7],ystyle=1,/nodata
endelse
oplot,slices,PSFcensref[0,*],psym=4,thick=5
oplot,slices,deltaxref+xc0,linestyle=1
oplot,slices,PSFcensref[1,*]+3,psym=4,thick=5
;oplot,slices,PSFcens[1,*],linestyle=0
oplot,slices,deltayref+yc0+3,linestyle=1

PSFcensrefcorr=fltarr(2,sz[3])
for s=0L,sz[3]-1 do begin
refimage[*,*,s]=shift_sub(refimage[*,*,s],-1*deltaxref[s],-1*deltayref[s])
cens[0,*,s]-=deltaxref[s]
cens[1,*,s]-=deltayref[s]
PSFcensrefcorr[*,s]=[(median(cens[0,*,s],/even)),(median(cens[1,*,s],/even))]

;switch if sat spot position is NaN ...
centotal=total(cens[*,*,s])
cenref=(total(cens[*,*,0],2)/4) # (fltarr(4)+1.)

if finite(centotal) eq 0 then begin
cens0=cens[*,*,0]
 cens0p=cv_coord(from_rect=cens0-cenref,/to_polar)
 cens[*,*,s]=cv_coord(from_polar=[cens0p[0,*],cens0p[1,*]/scl[s]],/to_rect)+cenref
PSFcensrefcorr[*,s]=[(median(cens[0,*,s],/even)),(median(cens[1,*,s],/even))]
endif

for j=0L,3 do begin

  sxaddpar,h1ref,'SATS'+strtrim(s,2)+'_'+strtrim(j,2),$
   string(strtrim(cens[*,j,s],2),format='(F7.3," ",F7.3)'),$
    'Location of sat. spot '+strtrim(j,2)+' of slice '+strtrim(s,2)
 endfor

sxaddpar,h1ref,'PSFC_'+strtrim(s,2),$
 string(strtrim(PSFcensrefcorr[*,s],2),format='(F7.3," ",F7.3)'),$
  'PSF Center of slice '+strtrim(s,2)

endfor
writefits,'refimage.fits',0,h0ref
writefits,'refimage.fits',refimage,h1ref,/append

;Now, loop through the science images and cross-correlate the halo to find offsets.

pixscale=charis_get_constant(name='pixscale')
for i=0L,nfiles-1 do begin
 
print,'Registering File Number ',i+1,' ',files[i]
 ;read in file
 a=readfits(datadir+files[i],h1,/exten,/silent)
 h0=headfits(datadir+files[i],ext=0,/silent)
 a_in=a
 sz=size(a)

 parangle0=float(sxpar(h0,'PA'))
 ha0=float(sxpar(h0,'HA'))
 
 ;center the pa value on zero
 ;need to fix this...
 if parangle0 gt 180 and decdeg lt lat then parangle0-=360
 
 ;*****Image Registration Loop.

;1. Define an annulus between the location of the occulting spot and the satellite spots
;2. cross correlate between reference and science images with the reference's initial offset as your starting guess

   print,PSFcensrefcorr[0,10],PSFcensrefcorr[1,10]
   dist_circle,cendistance,201,PSFcensrefcorr[0,10],PSFcensrefcorr[1,10]  
   roi=where(cendistance gt 0.2/pixscale and cendistance le 0.5/pixscale)
   cendistance[roi]=10

 ;**- If you have satellite spots then you have a choice of getting the centroid from gaussian-fitting or from joint cross-correlation

   for s=0L,sz[3]-1 do begin
    
    fidcenter=charis_getrelshift(refimage[*,*,s],a_in[*,*,s],roi,xc0-PSFcensref[0,s],yc0-PSFcensref[0,s])

    PSFcens[0,s]=fidcenter[0]
    PSFcens[1,s]=fidcenter[1]

   endfor

slices=findgen(sz[3])

good=where(finite(PSFcens[0,2:sz[3]-1-2]) ne 0 and finite(PSFcens[1,2:sz[3]-1-2]) ne 0,ngood)

fitx=robust_poly_fit((slices[2:sz[3]-1-2])[good]-refslice,(PSFcens[0,2:sz[3]-1-2])[good]-PSFcens[0,refslice],2,numit=15)
fity=robust_poly_fit((slices[2:sz[3]-1-2])[good]-refslice,(PSFcens[1,2:sz[3]-1-2])[good]-PSFcens[1,refslice],2,numit=15)

refoffsetx=PSFcens[0,refslice]-xc0
refoffsety=PSFcens[1,refslice]-yc0
deltax=fitx[0]+fitx[1]*(slices-refslice)+fitx[2]*(slices-refslice)^2.+refoffsetx
deltay=fity[0]+fity[1]*(slices-refslice)+fity[2]*(slices-refslice)^2.+refoffsety


if keyword_set(guessoffsets) then begin
plot,slices,PSFcens[0,*],yrange=[sz[1]/2-7+long(guessoffsets[0]),sz[1]/2+7+long(guessoffsets[1])],ystyle=1,/nodata
endif else begin
plot,slices,PSFcens[0,*],yrange=[sz[1]/2-7,sz[1]/2+7],ystyle=1,/nodata
endelse
oplot,slices,PSFcens[0,*],psym=4,thick=5
oplot,slices,deltax+xc0,linestyle=1
oplot,slices,PSFcens[1,*]+3,psym=4,thick=5
oplot,slices,deltay+yc0+3,linestyle=1

cenxref=cens[0,*,0]-deltax[0]
cenyref=cens[1,*,0]-deltay[0]

for s=0L,sz[3]-1 do begin
a_in[*,*,s]=shift_sub(a_in[*,*,s],-1*deltax[s],-1*deltay[s])
;***new stuff
PSFcens[0,s]-=deltax[s]
PSFcens[1,s]-=deltay[s]
sxaddpar,h1,'PSFC_'+strtrim(s,2),$
 string(strtrim(PSFcens[*,s],2),format='(F7.3," ",F7.3)'),$
  'PSF Center of slice '+strtrim(s,2)
endfor

a=a_in
psfcenx=median(PSFcens[0,*],/even)
psfceny=median(PSFcens[1,*],/even)
sxaddpar,h1,"PSFCENTX",psfcenx,'Mean PSF center X'
sxaddpar,h1,"PSFCENTY",psfceny,'Mean PSF center Y'

writefits,reducdir+filesout[i],0,h0
writefits,reducdir+filesout[i],a,h1,/append

if ~keyword_set(ladder) then printf,1,filenum[i],psfcenx,psfceny,rsat,ha0,parangle0
endfor

        end

;******end method - calsat

endcase

;****Check PSF quality ****
if (keyword_set(checkquality) and method eq 0) then begin
;if (keyword_set(checkquality)) then begin

test=headfits(reducdir+filesout[0])

get_charis_wvlh,h0,wvlhs
;Telescope Diameter for Subaru
;Dtel=7.9d0 ;visible pupil for SCExAO
Dtel=charis_get_constant(name='Dtel')
;pixel scale 16.4 mas/pixel
pixscale=charis_get_constant(name='pixscale')
;pixscale=0.0164
fwhm=1.*(1.d-6*wvlhs*1d-3/Dtel)*(180.*3600./!dpi)/pixscale

openw,2,'PSFquality'
sattohalo=fltarr(nfiles)
for ifile=0L,nfiles-1 do begin

inputcube=readfits(reducdir+filesout[ifile],/ext,h1input)

cens=fltarr(2,4,sz[3])

;for il=0L,sz[3]-1 do begin
for il=0L,0 do begin
;input slice
imslice=inputcube[*,*,il]

for j=0,3 do begin  ;assume 4 sat spots
   tmp=fltarr(2)+!values.f_nan
    sf=string(il)
    tmp= sxpar(h1input,'SATS'+strtrim(long(sf),2)+'_'+strtrim(j,2))
;**total, ugly hack-tastic workaround since I can't seem to get IDL to treat the string properly otherwise...
    tmpz=strsplit(tmp,' ',/extract)
    tmp=[tmpz[0],tmpz[1]]
    cens[*,j,il]=tmp
  endfor

aperradf=0.5*fwhm[il]
sat_skyrad=[2,6]*aperradf
phpadu=1.0

;***Sat flux
aper, imslice,cens[0,0,il],cens[1,0,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
     sat1flux=flux
 ;print,'flux and sky1 is',flux,sky
aper, imslice,cens[0,1,il],cens[1,1,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
     sat2flux=flux
 ;print,'flux and sky2 is',flux,sky
aper, imslice,cens[0,2,il],cens[1,2,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
     sat3flux=flux
 ;print,'flux and sky3 is',flux,sky
aper, imslice,cens[0,3,il],cens[1,3,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
     sat4flux=flux

satfluxf=median([sat1flux,sat2flux,sat3flux,sat4flux])


;***Halo flux
profrad_tc,imslice,1,1,60,p2d=prhalo,rayon=rayon
dist_circle,g,201
haloroi=where(g ge 0.2/pixscale and g le 0.4/pixscale,nhaloroi)
haloflux=median(prhalo[haloroi],/even)

sattohalo[ifile]=satfluxf/haloflux

endfor ;il

printf,2,long(ifile),' ',filesout[ifile],sattohalo[ifile]
endfor ;ifile

set_plot,'ps'
device,filename='halototcore.eps',/encapsulated
plot,long(findgen(nfiles)),sattohalo,psym=4,/nodata,xthick=5,ythick=5,charsize=1.25,xtitle='File Index',ytitle='Sat/Halo'
oplot,long(findgen(nfiles)),sattohalo,psym=4
device,/close

if ~keyword_set(bsize) then bsize=0.5
device,filename='halotocoredist.eps',/encapsulated
plothist,sattohalo,bin=bsize,peak=1,/nodata,xthick=5,ythick=5,xtitle='Sat/Halo'
plothist,sattohalo,bin=bsize,peak=1,thick=5,/overplot
device,/close

endif

close,1
close,2

skiptotheend:

end
