pro charis_pdi_register_cube,pfname,prefname=prefname,$
Lsuffname=Lsuffname,$
Rsuffname=Rsuffname,$
polcent_offset=polcent_offset,$
astrogrid=astrogrid,$
method=method,$
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
smask=smask,$
checkquality=checkquality,$
verbose=verbose,$
help=help

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
    print,'charis_pdi_register_cube.pro: Registers left and right pol cubes using charis_register_cube.'
    print,"Typically executed after 'charis_pdi_split_pols' and before 'charis_pdi_specphot_cal'."
    print,'Written by K. Lawson (2020)'
    print,''
    print,'**Calling Sequence**'
    print,"charis_pdi_register_cube,pfname,prefname=prefname,Lsuffname=Lsuffname,Rsuffname=Rsuffname,polcent_offset=polcent_offset,"
    print,"astrogrid=astrogrid,method=method,guessoffsets=guessoffsets,rsub=rsub,medbox=medbox,refcube=refcube,ladder=ladder,"
    print,"xcorr=xcorr,splitpsf=splitpsf,keepprev=keepprev,nosmooth=nosmooth,revise=revise,fwhmlim=fwhmlim,"
    print,"smallsteps=smallsteps,smask=smask,checkquality=checkquality,verbose=verbose,help=help"
    print,''
    print,'Example:'
    print,"charis_pdi_register_cube,'abaur_low.info',astrogrid=11.2"
    print,''
    print,"***Keywords***"
    print,'pfname=STRING - The parameter file (e.g. HR8799_low.info)'
    print,"prefname=STRING - The prefix of files to use as input. Defaults to 'n'."
    print,"Lsuffname=STRING - The suffix of left pol files to use as input. Defaults to 'left'."
    print,"Rsuffname=STRING - The suffix of right pol files to use as input. Defaults to 'right'."
    print,"polcent_offset=ARRAY - Two element array providing the pixel coordinate offset to the center of the right pol from the nominal center. Defaults to [31.758, 16.542]."
    print,"Other keywords as accepted by 'charis_register_cube'."
    goto,endofprogram
endif

if ~keyword_set(polcent_offset) then polcent_offset=[31.758, 16.542]
if ~keyword_set(Lsuffname) then Lsuffname='left'
if ~keyword_set(Rsuffname) then Rsuffname='right'

;* Figuring out what charis_register_cube is going to name our products:
datadir='./reduc/prep/'
reducdir='./reduc/reg/'

if ~keyword_set(ladder) then param,'fnum_sat',flist,/get,pfname=pfname else param,'fnum_lad',flist,/get,pfname=pfname
if ~keyword_set(prefname) then prefname='n'

filenum=nbrlist(flist)
filesout=filelist(filenum,nfiles,prefix=prefname,suffix='reg')
;* Names to change defaults to for left and right products:
Lfilesout=filelist(filenum,prefix=prefname,suffix=Lsuffname+'reg')
Rfilesout=filelist(filenum,prefix=prefname,suffix=Rsuffname+'reg')

; Registration for left pol frames:
charis_register_cube,pfname,prefname=prefname,suffname=Lsuffname,method=method,guessoffsets=-1*polcent_offset,$
  rsub=rsub, refcube=refcube, medbox=medbox, ladder=ladder, xcorr=xcorr, revise=revise, splitpsf=splitpsf,$
  keepprev=keepprev, nosmooth=nosmooth, fwhmlim=fwhmlim, smallsteps=smallsteps, smask=smask,$
  checkquality=checkquality, astrogrid=astrogrid, verbose=verbose, help=help

spawn,'mv reduc.log reduc_left.log'
for i=0L,nfiles-1 do spawn,'mv '+reducdir+filesout[i]+' '+reducdir+Lfilesout[i]

; Registration for right pol frames:
charis_register_cube,pfname,prefname=prefname,suffname=Rsuffname,method=method,guessoffsets=polcent_offset,$
  rsub=rsub, refcube=refcube, medbox=medbox, ladder=ladder, xcorr=xcorr, revise=revise, splitpsf=splitpsf,$
  keepprev=keepprev, nosmooth=nosmooth, fwhmlim=fwhmlim, smallsteps=smallsteps, smask=smask,$
  checkquality=checkquality, astrogrid=astrogrid, verbose=verbose, help=help

spawn,'cp reduc.log reduc_right.log'; Just copy instead of moving this one to leave reduc.log for total intensity PSF sub later.
for i=0L,nfiles-1 do spawn,'mv '+reducdir+filesout[i]+' '+reducdir+Rfilesout[i]

endofprogram:
end
