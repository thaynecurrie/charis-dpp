; Procedure written by David Lafreniere
; modified by Markus Janson
; modified for simplicity and speed by Thayne Currie for CHARIS, added additional keywords to reduce user interactions
;v 3.0 - 08112020- added new target search parameters to fix SIMBAD errors
;v 2.0 - redone 

pro charis_newobs, date=date,object=object,filter=filter,help=help,outfilename=outfilename,cfa=cfa

if (keyword_set(help)) then begin
print,"charis_newobs,date=[date],object=[objectname],filter=[filter]('low','J','H','K')"
print,""
print,"***Keywords***"
print,"[keywords are self-explanatory]"
print,"***"
print,"Example: charis_newobs,date='20190929',object='HR 8799',filter='low'"
goto,skiptoend
endif

if (N_PARAMS() eq 0 and ~keyword_set(object) and ~keyword_set(date) and ~keyword_set(filter)) then begin
print,"charis_newobs,date=[date],object=[objectname],filter=[filter]('low','J','H','K')"
print,""
print,"***Keywords***"
print,"[keywords are self-explanatory]"
print,"***"
print,"Example: charis_newobs,date='20190929',object='HR 8799',filter='low'"
goto,skiptoend
endif

;***Set up directory structure, move files to ./data/raw/

charis_setupdir
;SPAWN,'mv CRSA*fits ./data/raw/'
if (!version.os_family EQ 'Windows') then begin
SPAWN,'move CRSA*cube.fits ./data/raw/'
endif else begin
SPAWN,'mv CRSA*cube.fits ./data/raw/'
endelse

;***

;***INCOMPLETE as of now

if ~keyword_set(date) then begin
date=''
read,'Enter date: ',date
endif
if ~keyword_set(object) then begin
object=''
read,'Enter object name: ',object
endif
if ~keyword_set(filter) then begin
filter=''
read,'Enter filter name: (low, J, H, or K)',filter
endif

suffix='_'+string(filter)

obj=strcompress(object,/remove_all)

obj_main = obj
pfname=obj+suffix+'.info'
outfilename=pfname
param,'name',object,/set,pfname=pfname

;;******OBJECT INFORMATION *******
print,'Searching for ',object
   charis_querysimbadf,obj,ra,dec,id,found=found
   if found eq 0 then begin
    print,'Trying ','* ',obj
    charis_querysimbadf,'* '+obj,ra,dec,id,found=found
   endif

   if found eq 0 then begin
    print,'Trying ','* ',object
    charis_querysimbadf,'* '+object,ra,dec,id,found=found
   endif

   if found eq 0 then begin
    print,'Trying ','V* ',obj
    charis_querysimbadf,'V* '+obj,ra,dec,id,found=found
   endif
   if found eq 0 then begin
    print,'Trying ','V* ',object
    charis_querysimbadf,'V* '+object,ra,dec,id,found=found
   endif

if found eq 0 then begin
  Print,'Cannot Find Star, Check the name'
  stop
endif

   print,ra,dec
;2mass database for JHK
   flag_s2mass=0
   s2mass=charis_queryvizierf('2MASS-PSC',[ra,dec],20./60.,/allcolumns)

;HIP database for parallax, vmag, stype etc.

   flag_ship =0
   ship=charis_queryvizierf('I/239/hip_main',[ra,dec],20./60.,/allcolumns)
   if size(ship,/type) ne 8 then flag_ship=1
   ;print,ship.e_plx


;Put parameters into parameter file

;from charis_querysimbadf
param,'ra',ra,'RA, equinox J2000, epoch J2000',/set,pfname=pfname
param,'dec',dec,'DEC, equinox J2000, epoch J2000',/set,pfname=pfname
param,'ra_hms',deg2hms(ra),'RA, equinox J2000, epoch J2000',/set,pfname=pfname
param,'dec_dms',deg2dms(dec),'DEC, equinox J2000, epoch J2000',/set,pfname=pfname


if flag_ship eq 1 then begin
plx=0
eplx=0
pmra=0
epmra=0
pmde=0
epmde=0
sptype='-999'
distance=plx
endif else begin

;Parallax and error on parallax might not have been found
plxset = where(strmatch( tag_names(ship), 'plx', /fold_case) eq 1)
eplxset = where(strmatch( tag_names(ship), 'e_plx', /fold_case) eq 1)
pmraset = where(strmatch( tag_names(ship), 'pmra', /fold_case) eq 1)
epmraset = where(strmatch( tag_names(ship), 'e_pmra', /fold_case) eq 1)
pmdecset = where(strmatch( tag_names(ship), 'pmde', /fold_case) eq 1)
epmdecset = where(strmatch( tag_names(ship), 'e_pmde', /fold_case) eq 1)
sptypeset=where(strmatch( tag_names(ship),'sptype', /fold_case) eq 1)


if (plxset gt 0) then plx=ship.plx else begin
    plx=0
    print, 'Parallax not found. Setting parallax=0.'
endelse
if (eplxset gt 0) then eplx=ship.e_plx else begin 
    eplx=0
    print, 'Error on parallax not found. Setting error on parallax=0.'
endelse

;sometimes the HIP catalog triggers a stupid error ...
if n_elements(plx) gt 1 then plx=plx[0]
if n_elements(eplx) gt 1 then eplx=eplx[0]

if plx gt 0 then distance = 1.d3/plx else distance =-999

if (pmraset gt 0) then begin
pmra=ship.pmra
epmra=ship.e_pmra
endif else begin
;if (pmraset gt 0) then (pmra=ship.pmra) and (epmra=ship.epmra) else begin
    pmra=0
    epmra=0
    print, 'PM RA not found. Setting (error on) PMRA=0.'
endelse
if (pmdecset gt 0) then begin
pmde=ship.pmde
epmde=ship.e_pmde
endif else begin 
    pmde=0
    epmde=0
    print, 'PM DEC not found. Setting (error on) PMDEC=0.'
endelse

if (sptypeset gt 0) then sptype=ship.sptype else begin
sptype = '-999'
Print, 'Spectral type not found, setting to dummy variable.'
endelse

endelse

if n_elements(pmra) gt 1 then pmra=pmra[0]
if n_elements(epmra) gt 1 then epmra=epmra[0]
if n_elements(pmde) gt 1 then pmde=pmde[0]
if n_elements(epmde) gt 1 then epmde=epmde[0]
if n_elements(sptype) gt 1 then sptype=sptype[0]

;Write parallax and proper motion + associated errors
param,'plx',plx,' Parallax, mas',/set,pfname=pfname
param,'distance',distance,' Distance (pc)',/set,pfname=pfname
param,'eplx',eplx,' Error on parallax, mas',/set,pfname=pfname
param,'pmra',pmra,' Proper motion RA [*cos(DEC)], mas/yr',/set,pfname=pfname
param,'epmra',epmra,' Error on PMRA',/set,pfname=pfname
param,'pmdec',pmde,' Proper motion DEC, mas/yr',/set,pfname=pfname
param,'epmdec',epmde,' Error on PMDEC',/set,pfname=pfname
param,'sptype',sptype,' Spectral Type',/set,pfname=pfname


;Presumably, your source is in the 2MASS point source catalog!
  ;print,'radec',ra,dec
  s=charis_queryvizierf('II/246/out',[ra,dec],20./60.,/allcolumns)
  ind=sort(s._r)
  s=s[ind[0]]

;Put in spectral type if you have it

param,'jmag',s.jmag,'J magnitude, 2MASS',/set,pfname=pfname
param,'ejmag',s.e_jmag,'Error on JMAG',/set,pfname=pfname
param,'hmag',s.hmag,'H magnitude, 2MASS',/set,pfname=pfname
param,'ehmag',s.e_hmag,'Error on HMAG',/set,pfname=pfname
param,'kmag',s.kmag,'Ks magnitude, 2MASS',/set,pfname=pfname
param,'ekmag',s.e_kmag,'Error on KMAG',/set,pfname=pfname


param,'obsdate',date,'Observation date',/set,pfname=pfname
param,'fnum_sat','1-100','File #s, sat, main seq',/set,pfname=pfname
param,'fnum_sky','101-103','File #s, sky frames',/set,pfname=pfname
param,'fnum_lad','1-5','File #s, ladder',/set,pfname=pfname

;default, canned values for other parameters
;other algorithm parameters
;unclear whether the below is needed
param,'procdir','./proc/','Subdirectory to put processed files',/set,pfname=pfname
param,'rsat',-1.,'PSF saturation radius',/set,pfname=pfname
param,'fwhm',2.63,'Broadband PSF FWHM in pixels',/set,pfname=pfname


;Algorithm Parameters

;LOCI and KLIP
param,'****','******PSF Sub Range',/set,pfname=pfname

param,'nfwhm',0.75,'Min displacement for sub, in # of PSF FWHM',/set,pfname=pfname
param,'dr',5.,'Width of subtraction annuli',/set,pfname=pfname
;param,'dr','5.,50.,120.,1.','Width of subtraction annuli, LOCI',/set,pfname=pfname
param,'rmin',5,'Min radius for PSF sub',/set,pfname=pfname
param,'rmax',65,'Max radius for PSF sub',/set,pfname=pfname

param,'*****','******Switches for Sub/Comb',/set,pfname=pfname
param,'meanadd',0,'median-comb(=0), mean-comb(=1)',/set,pfname=pfname
param,'zero',0,'no sub median annulus(=0), sub med. annulus(=1)',/set,pfname=pfname
param,'rsub',0,'no rad-prof sub.(=0), rad-prof sub. b4 PSF sub(=1)',/set,pfname=pfname

;(A-)LOCI parameters
param,'******','******(A-)LOCI Params',/set,pfname=pfname
;param,'***Comment***','(A-)LOCI Params',/set,pfname=pfname
param,'locina',100.,'Area of optimization regions, LOCI',/set,pfname=pfname
param,'locigeom',1.,'Geometry of optimization regions, LOCI',/set,pfname=pfname
param,'alocisvd',-6,'log of the SVD cutoff, A-LOCI',/set,pfname=pfname
param,'alocinref',50,'frame selection for A-LOCI',/set,pfname=pfname
param,'pixmask',0,'no pixmask(=0),pixmask(=1)',/set,pfname=pfname
;param,'zonetype',0,'A-LOCI optzone type',/set,pfname=pfname

param,'*******','******KLIP Params',/set,pfname=pfname
param,'klipnpca',5,'Number of KL Modes, KLIP',/set,pfname=pfname
param,'klipnzone',1,'Number of azimuthal zones per annulus, KLIP',/set,pfname=pfname
param,'meansub',1,'no subtract mean(=0), subtract mean(=1)',/set,pfname=pfname

;***Planet Model Parameters
param,'********','PlntFwdMod',/set,pfname=pfname
param,'planmod','L0_g.fits','Model Plan Spec',/set,pfname=pfname

skiptoend:


end
