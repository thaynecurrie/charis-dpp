pro charis_imprep,pfname,$
gz=gz,$
fixpix=fixpix,$
charismanual=charismanual,help=help
;Establishes uniform fits header information for all science frames you want to reduce
;individual integration time, coadds (so int.*coadds = total int), altitude, azimuth, LST, PA, beam (if dither), latitude, longitude, HA, WCS

;08/11/2020 -- added ability to edit list of files [procedure copied from HCI]
;04/18/2019 -- streamlined for CHARIS.  removed commands for other instruments
;03/06/2018 -- copied from imprep.pro to package with CHARIS post-processing pipeline

if ~keyword_set(pfname) then begin
possiblenames=file_search('*.info')
print,'possible ',possiblenames

if strlen(possiblenames) gt 4 then begin
pfname=possiblenames
print,'Only one info file found ',pfname
print,'...using this one'
endif else begin

pfname0=''
read,pfname0,Prompt="Select Info File (Next time, Just enter it at command line) "
pfname=file_search(pfname0)
pfname=pfname0
endelse

endif

if (keyword_set(help)) then begin
print,"charis_imprep,pfname,gz=gz,fixpix=fixpix,charismanual=charismanual"
print,"usually, do not need to enter any keywords *except* for pfname"
print,""
print,"***Keywords***"
print,"pfname - parameter file name"
print,"fixpix - do bad pixel repair in the cube"
print,"*****"
print,"Example: charis_imprep,'HR8799_low.info'"
goto,breakout
endif

if ~keyword_set(prefname) then prefname='n'
datadir='./data/raw/'

reducdir1='./data/'
if ~keyword_set(reducdir) then begin
subdir='prep/'
reducdir=reducdir1+subdir
endif else begin
reducdir=reducdir1
endelse

;if keyword_set(cal) then begin
;reducdir='./data/cal/'
;endif

file_mkdir,reducdir

print,'reducdir is',reducdir

;First, if you just have the CHARIS data cubes in your current directory, move them to datadir
a=file_test('CRSA*fits')
if a eq 1 then begin
com='mv '+'CRSA*fits ' +datadir
SPAWN,com
endif


;Second, take all the files in datadir and save their names as a string array 'files'
files=file_search(datadir+'*fit*')

;****This is not needed for CHARIS right now.  Keep some of the source code in case problem pops up.
goto,skipsort
if keyword_set(sort) then begin
times=dblarr(n_elements(files))
for i=0L,n_elements(files)-1 do begin
h1=headfits(files[i])

;****THIS IS A HACK RIGHT NOW!!!!!

;times[i]=double(sxpar(h1,'mjd'))
times[i]=ten(sxpar(h1,'lst'))
;times[i]=ten(sxpar(h1,'p_hststr'))

if times[i] lt 20 then times[i]+=24
print,'times',times[i],' ',files[i]
endfor

files=files(sort(times))
times=times(sort(times))

endif
skipsort:
;*****



nfiles=n_elements(files)

print,'nfiles is',nfiles

;Instrument Class: SCExAO/CHARIS
;we are at Maunakea
  lat=19.825504d0
  lng=-155.47602d0

;test the first image, see if it has the fits header info you need.

bbtest=readfits(files[0],h1test,ext=1)
htest=headfits(files[0])

;if you want to manually set the exposure time and coadds then input them here...
if keyword_set(charismanual) then begin

read,'**Enter exposure time (seconds) ',exp1time_manual

read,'**Enter number of coadds ',coadds_manual
endif

;Obs log
obslogname='obslog.txt'
openw,1,obslogname

;openw,1,'reduc.log'
;printf,1,'File_no,','XC,','YC,','Rsat','HA','Par.Angle'
;printf,1,'File_in','File_out','ExpTime','ParAngle','Health_Flag'
printf,1,'File_in','','File_out','','ExpTime','','ParAngle','','Health_Flag'

for i=0L,nfiles-1 do begin
bb=readfits(files[i],h1,/exten)
h0=headfits(files[i],ext=0)

sz=size(bb)

;if you don't manually set the exposure time and coadds, then use the number of reads from the file to determine exptime
  if keyword_set(charismanual) then begin
  exp1time=exp1time_manual
  endif else begin
  exp1time=1.475
  endelse
   
  if keyword_set(charismanual) then begin
  coadds=coadds_manual
  endif else begin
  firstread=sxpar(h1,'firstrd')
  lastread=sxpar(h1,'lastrd')
  coadds=long((lastread-firstread)+1)
  endelse


  lst1=sxpar(h0,'st',count=lst1count)
  lst2=sxpar(h0,'lst',count=lst2count)
  if lst1count gt 0 then lst=lst1
  if lst2count gt 0 then lst=lst2

 ;Hour Angle, probably is not entered in here.
  ha=sxpar(h0,'HA',count=hacount)

  if hacount gt 0 then begin
  ha=ten(ha)
  sxdelpar,h0,'HA'
  sxaddpar,h0,'HA',ha,after='lst',format="f9.5"
  endif

 ;;get the date, then calculate the epoch 
  dateobs=sxpar(h0,'UTC-DATE')
  utcdate=ten(sxpar(h0,'UTC-TIME'))
  hms0=sixty(utcdate)


  ymd0 = double(strsplit(dateobs,'-',/extract))
  jd0=double(julday(ymd0[1], ymd0[2], ymd0[0],hms0[0],hms0[1],hms0[2]))

  epoch0 = (jd0 - 2451545d0)/365.25d0 + 2000d0
   ;print,jd0,epoch0

  ;***if the LST is not entered, then you need to compute it based on the MJD and longitude of Maunakea
  if lst1count eq 0 and lst2count eq 0 then begin
   jd=2400000.5+double(sxpar(h0,'mjd'))
   ct2lst,lst,lng,blah,jd

  ;now pull the right ascension and compute the difference between LST and RA
   ;put ra and dec in degrees for now.

  ;ra_string=strsplit(sxpar(h0,'ra'),':',/extract)
  ;ra=ten(ra_string[0],ra_string[1],ra_string[2])

  ra=15*ten(sxpar_charis(h0,'ra',count=racount,/justfirst))
  dec=ten(sxpar_charis(h0,'dec',count=racount,/justfirst))

  endif

  ;precess coordinates to current epoch
  rap=ra
  decp=dec
  precess,rap,decp,2000d0,epoch0

  ;compute the hour angle for current epoch, ra2000, and lst
  ha=lst-rap/15.


  if hacount eq 0 then begin
  sxaddpar,h0,'HA',ha,after='SHUTTER',format="f9.5"
  endif
  
  sxaddpar,h0,'lst',lst,after='HA',format="f9.5"


;now put in integration time
  sxaddpar,h0,'exp1time',exp1time,after='altitude',format="f9.5"
  sxaddpar,h0,'coadds',coadds,after='exp1time',format="f9.5"
  
  sxaddpar,h0,'telescop','Subaru',after='maskivar',format="a10"
  sxaddpar,h0,'lat',lat,after='telescop',format="f9.5"
  sxaddpar,h0,'lng',lng,after='lat',format="f12.7"

;And, parallactic angle
  pa0=sxpar(h0,'PARANG')
  sxaddpar,h0,'PA',pa0,after='PARANG',format="f12.6"


;****wavelength information****
;...will add later
sxaddpar,h1,'CRPIX3',1,after='CRPIX2',format="i4"
sxaddpar,h1,'CUNIT3','microns',after='CRPIX3',format="a8"
;sxaddpar,h1,'CRVAL3',crvalue3


;***Obs Log ***
;if keyword_set(obslog) then begin
;input file name: files[i]
;output file name: fileoutname
;exposure time: exp1time*coadds
;PA: pa0
;Health Flag: science, sky, bad frame

filtname=sxpar(h0,'filtname')
filebasename=files[i]
filebasename=strsplit(filebasename,/extract,'/')
fileinname=filebasename[n_elements(filebasename)-1]
fileoutname=prefname+nbr2txt(i+1,4)+'e.fits'

cubehealthflag=charis_check_health(bb)

printf,1,fileinname,fileoutname,filtname,exp1time*coadds,pa0,cubehealthflag,format='(a,1x,a,1x,a,1x,f6.2,1x,f7.2,1x,a)'

sz=size(bb,/dim) 
if keyword_set(fixpix) then begin
for j=0L,sz[2]-1 do begin

slice=bb[*,*,j]
slice=fixpix_rs(slice,iter=10)
bb[*,*,j]=slice
endfor
endif

;For now, just assume CHARIS will create the data cubes

fileoutname=prefname+nbr2txt(i+1,4)+'e.fits'

reducdir='./reduc/prep/'
;if for some reason you want to g-zip the data cubes
if keyword_set(gz) then begin
writefits,reducdir+fileoutname+'.gz',0,h0,/compress
writefits,reducdir+fileoutname+'.gz',bb,h1,/compress,/append
endif else begin
writefits,reducdir+fileoutname,0,h0
writefits,reducdir+fileoutname,bb,h1,/append
endelse
endfor

close,1

charis_set_files,pfname,'FNUM_SAT','science',obslogname
charis_set_files,pfname,'FNUM_SKY','sky',obslogname






breakout:

end
