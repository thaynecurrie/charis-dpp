pro charis_ssdi,pfname,rsub=rsub,weight=weight,filt=filt,suffname=suffname,aloci=aloci,klip=klip,northup=northup,help=help

;05/25/2021 - updated!

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"Simple SSDI of spatially-filtered cubes or PSF subtracted cubes..."
print,"charis_ssdi,pfname,rsub=rsub,weight=weight,filt=filt,suffname=suffname,aloci=aloci,pca=pca,help=help"
print,""
print,"Example: charis_ssdi,'HR8799_low.info',/weight,/aloci"
print,""
print,"***Keywords***"
print,"*rsub - if thrown, SSDI the spatially-filtered (not PSF subtracted) cubes"
print,"*aloci - do SSDI on ADI/A-LOCI residuals"
print,"*klip - do SSDI on ADI/KLIP residuals"
print,"*weight - residual-minimized weight applied to each slice for SDI subtraction"
print,"*filt - spatially filter after SSDI subtraction"
print,"*northup - rotate northup"

goto,skiptotheend
endif

;dumb,simple classical SDI of post-ADI residuals

prefname='n'

param,'fnum_sat',flist,/get,pfname=pfname
filenum=nbrlist(flist)

reducdir='./reduc/'
if ~keyword_set(rsub) then begin
subdir='proc/'
endif else begin
subdir='rsub/'
endelse

if keyword_set(rsub) then begin
readcol,'reduc.log',a1,a2,a3,a4,a5,pa
theta=-1*(pa-pa[0])
endif

datadir=reducdir+subdir
reducdir=datadir

sufftype='_alocisub'
if keyword_set(aloci) then sufftype='_alocisub'
if keyword_set(klip) then sufftype='_adisub'

if ~keyword_set(suffname) then begin
if ~keyword_set(rsub) then begin
files=filelist(filenum,nfiles,prefix=prefname,suffix=sufftype)
endif else begin
files=filelist(filenum,nfiles,prefix=prefname,suffix='rsub')
endelse
endif else begin
files=filelist(filenum,nfiles,prefix=prefname,suffix='_'+suffname)

endelse

filesout=filelist(filenum,prefix=prefname,suffix='_ssdisub')
filesref=filelist(filenum,prefix=prefname,suffix='_refssdi')

;test image
testimage=readfits(datadir+files[0],ext=1,h1test)
h0test=headfits(datadir+files[0])

sz=size(testimage,/dim)

get_charis_wvlh,h0test,wavelengths
lambda=wavelengths*1.d-3
nlambda=n_elements(lambda)
;print,lambda

;loop on image
 ;loop on wavelength

 for i=0L,nfiles-1 do begin
 print,"SDI'ing Data Cube Number ",long(i+1)
 imcube=readfits(datadir+files[i],h1,ext=1)
 h0=headfits(datadir+files[i])
 refimage=fltarr(sz[0],sz[1],sz[2])

  imcubebefore=imcube

  for il=0L,nlambda -1 do begin
    imcuberef=imcubebefore
    ;imcuberef=imcube
    ;print,h0,h1
    charis_alignspeckle,imcuberef,h0,h1,/nolocs,refslice=il
    imsliceref=median(imcuberef,/even,dimension=3)
    
    refimage[*,*,il]=imsliceref

    if ~keyword_set(weight) then begin
    imcube[*,*,il]-=imsliceref
    endif else begin
    weightpsf=total(imcube[*,*,il]*imsliceref,/nan)/total(imsliceref*imsliceref,/nan)
    imcube[*,*,il]-=weightpsf*imsliceref
    endelse

    if keyword_set(filt) then imcube[*,*,il]-=filter_image(imcube[*,*,il],median=15)

    ;if keyword_set(rsub) then begin
    ;rotate cube
    
    ;imf=rotat(imcube[*,*,il],theta[i],missing=!values.f_nan)
    ;imcube[*,*,il]=imf
    ;endif 
    ;stop

  endfor
    if keyword_set(rsub) then begin
    charis_northup,imcube,h0,h1
    endif

    writefits,reducdir+filesout[i],0,h0
    writefits,reducdir+filesout[i],imcube,h1,/append
 endfor
    imcombine=medfitsme(filesout,dir=reducdir,/cube)

    if keyword_set(northup) then begin
     imcomb0=imcombine
     charis_northup,imcomb0,h0test,h1test
   endif
    writefits,'ssdi.fits',0,h0test
    writefits,'ssdi.fits',imcombine,h1test,/append

 ;wavelength collapsed
    ;added filtering
    if keyword_set(filt) then for i=0L,21 do imcombine[*,*,i]-=filter_image(imcombine[*,*,i],median=9)
    writefits,'ssdicol.fits',0,h0test
    writefits,'ssdicol.fits',median(imcombine,/even,dimension=3),h1test,/append

skiptotheend:
end
