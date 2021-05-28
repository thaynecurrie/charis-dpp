pro charis_check_psfquality,pfname,bsize=bsize,suffname=suffname,ladder=ladder

param,'fnum_sat',flist,/get,pfname=pfname
if keyword_set(ladder) then $
param,'fnum_lad',fladder,/get,pfname=pfname

;********Directory Structure*****
; setupdir,reducdir=reducdir
reducdir='./reduc/'

;determine reduction subdirectory
subdir='reg/'
reducdir+=subdir

prefname='n'
if ~keyword_set(suffname) then suffname='reg'
filenum=nbrlist(flist)
filesout=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)

if keyword_set(ladder) then begin 
filenumlad=nbrlist(fladder)
filesladder=filelist(filenumlad,nfileslad,prefix=prefname,suffix=suffname)
endif


h0=headfits(reducdir+filesout[0])
test=readfits(reducdir+filesout[0],ext=1)

sz=size(test)
print,sz
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

;sat1flux=fltarr(sz[3])
;sat2flux=fltarr(sz[3])
;sat3flux=fltarr(sz[3])
;sat4flux=fltarr(sz[3])
;satfluxf=fltarr(sz[3])
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




end
