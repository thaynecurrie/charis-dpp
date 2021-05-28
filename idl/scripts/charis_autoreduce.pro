pro charis_autoreduce,quicklook=quicklook,help=help,target=target,liststeps=liststeps

if keyword_set(liststeps) then begin
print,'1.charis_newobs'
print,'2.charis_imprep'
print,'3.charis_subtract_sky'
print,'4.charis_register_cube'
print,'5.charis_specphot_cal'
print,'6.charis_makeemppsf'
print,'7.charis_imrsub'
print,''
print,'8.charis_adialoci'
print,'AND/OR'
print,'charis_sdialoci'
print,'OR'
print,'charis_adiklip/charis_rdiklip'
print,''
print,'9.charis_aloci_fwdmod_planet'
print,'OR'
print,'charis_klip_fwdmod_planet'
print,''
print,'10.charis_extract_1d_spectrum'
print,'11.charis_aloci_attenmap_planet'
print,'12.charis_calc_final_contrast'
goto,skiptotheend

endif

if keyword_set(help) then begin

print,"CHARIS DPP Auto-Reduce Program"
print,"Does automatic reduction of CHARIS IFS data "
print,""
print,"Example: charis_autoreduce
print,""
print,"****Keywords****"
print,"quicklook - Is this a quicklook reduction [no calibration] or an end-to-end reduction?"
print,"target - manually sets the target name"
print,"liststeps - do not run autoreduce but just list the reduction steps: helpful for novice users"
goto,skiptotheend
endif

Print,"CHARIS DPP Auto-Reduce Program"
Print,"Use at your own risk! Not suitable for Publication-Grade Results!  You Must Check Output!"

;Step 0. 
;Find data cubes in home directory or in ./data/raw
;Take out the middle one of the list [to avoid header errors]
;Pull the target name from the headers

;find cubes
datacubes=file_search('*cube.fits',count=ndatacubes)
if ndatacubes lt 1 then begin

datacubes=file_search('./data/raw/*cube.fits',count=ndatacubes)
if ndatacubes lt 1 then begin
print,"cannot find any data cubes.  Make sure you are in the correct home directory!"
goto,skiptotheend
endif

endif

;take the middle
excube=long(float(ndatacubes)/2.)
scratch=headfits(datacubes[excube])
if ~keyword_set(target) then begin
targname=sxpar(scratch,'OBJECT',count=ntargcount)
if ntargcount eq 0 then begin
targname=''
read,'Target name not found in header, select your target name: ',targname
endif

;if the target name is not given, search on SIMBAD for the real target name


endif else begin

targname=target
endelse

;print,ntargcount,ra,dec
filts=strtrim(sxpar(scratch,'FILTNAME'),2)
if filts eq 'Broadband' then filts = 'low'
obsdate=strjoin(strsplit(sxpar(scratch,'UTC-Date'),/extract,'-'))
print,"Your Target Is ",targname
print,'date is',obsdate
print,'your filter is ',filts

step1:
close,/all
charis_newobs,object=targname,date=obsdate,filter=filts,outfilename=pfname

step2:
close,/all
charis_imprep,pfname

step3:
close,/all

;Now, check to see if there are any sky frames
readcol,'obs*txt',n1,n2,n3,n4,n5,n6,format='(a,a,a,a,a,a)'
skyframe=where(n6 eq 'sky',numbersky) 

if numbersky ge 1 then begin
charis_subtract_sky,pfname

step4:
close,/all
charis_register_cube,pfname
endif else begin
charis_register_cube,pfname,suffname='e'
endelse

step5:
close,/all
charis_specphot_cal,pfname,modamp=25,starlib=1

step6:
close,/all
charis_makeemppsf,pfname

step7:
close,/all
charis_imrsub,pfname,/prad,/psfsub

step8:
close,/all
charis_adialoci,pfname,rsubval=1,outfile='aloci.fits'
;now go find sources in the image
;charis_adiklip,pfname,rsubval=1,outfile='klip.fits'

charis_check_for_detection,datacube='aloci.fits',outfile='aloci_det.targ',numdet=numdet

if numdet eq 0 then goto,step11


step9:
close,/all
charis_aloci_fwdmod_planet,pfname,reducname='aloci.fits',method=1,targfile='aloci_det.targ'


step10:
close,/all
readcol,'aloci_det.targ',xpos,ypos,contrast
for i=0L,n_elements(xpos)-1 do begin
charis_extract_1d_spectrum,pfname,coords=[xpos[i],ypos[i]],datacube='aloci.fits',/throughputcor,/fitcollapse,outputspec='aloci_auto_spec_'+strtrim(i+1,2)+'.dat',/filtsnr
endfor

step11:
close,/all
charis_aloci_attenmap_planet,pfname,ntc=5,reducname='aloci.fits'

step12:
close,/all
charis_calc_final_contrast,datacube='aloci.fits'


print,'all done!'


;Step 1. 



skiptotheend:
end
