pro charis_empbdplanspec,planetspec=planetspec,planetcovar=planetcovar,library=library,pick=pick,band=band,delay=delay,saveplot=saveplot,calcind=calcind,justind=justind

;compares a CHARIS planet spectrum to spectral libraries w/ option to include covariance
;*** Currently ONLY WORKS with Low-Res!!!! ****

;if you include the covariance, calculate it here once, so you do not have to recalculate it
if keyword_set(planetcovar) then begin

readcol,planetcovar,arho,alam,adel,sigrho,siglam,rho

;now populate the matrix C to get errors
endif


if keyword_set(saveplot) then begin
file_mkdir,'empbdplanspec_out'
plotdir='./empbdplanspec_out/'
endif

;initiate array of CHARIS data
get_charis_wvlh,dum,wavelengths,manual='lowres'
lambda=1d-3*wavelengths

limit1=1d-1
limit2=2.d-1

;now down-select the wavelengths to identify those that are not covered over the water bands.
;lmbad=[1400.,1900.17]*1d-3 ;wavelengths covering the water bands
;lmbad=[1374.62,1867.17]*1d-3 ;wavelengths covering the water bands
;bad=where(abs(lambda-lmbad[0]) lt limit1 or abs(lambda-lmbad[1]) lt limit2,nbad,complement=good)

bad=where(lambda eq lambda[5] or lambda eq lambda[6] or lambda eq lambda[7] $
or lambda eq lambda[14] or lambda eq lambda[15] or lambda eq lambda[16],nbad,complement=goodinit)

if ~keyword_set(band) then begin 
good=goodinit
endif else begin
;further down-select by wavelength
case band of
     'j': begin
        bad2=where(lambda gt lambda[4],badband,complement=goodband)
          end
     'h': begin
        bad2=where(lambda lt lambda[8] or lambda gt lambda[13],badband,complement=goodband)
          end
     'k': begin
        bad2=where(lambda lt lambda[17],badband,complement=goodband)
          end
     'hk': begin
        bad2=where(lambda lt lambda[8],badband,complement=goodband)
           end
     'jhk': begin
        bad2=where(lambda lt lambda[1],badband,complement=goodband)
            end
endcase
good=intersect(goodband,goodinit)
;ngood=n_elements(good)
endelse

ngood=n_elements(lambda[good])
print,'Wavelengths Considered',lambda[good],' out of ',lambda
print,'Number of Wavelengths',n_elements(lambda[good])


;goto,skipoverme
;*****Read in Planet Spectrum*****
;***assumes F_nu (mJy) units
if ~keyword_set(planetspec) then begin
planetspec=dialog_pickfile(Title="Select Your CHARIS Planet Spectrum")
endif
readcol,planetspec,wavelength_planet,planet_flux,eplanet_flux
skipoverme:

;****Now, decision tree to determine whether this is lowres, J, H, or K...
minwvlhp=1d3*min(wavelength_planet)
maxwvlhp=1d3*max(wavelength_planet)
if minwvlhp lt 1200 and maxwvlhp gt 2300 then filtername0='lowres'
if minwvlhp lt 1200 and maxwvlhp lt 1400 then filtername0='J'
if minwvlhp lt 1500 and minwvlhp gt 1400 and maxwvlhp lt 1850 then filtername0='H'
if minwvlhp gt 1850 then filtername0='K'


;***Calculate Spectral Covariance in each combination of i and j channels

;if you include the covariance, calculate it here once, so you do not have to recalculate it
if keyword_set(planetcovar) then begin

readcol,planetcovar,arho,alam,adel,sigrho,siglam,rho
nwvlh=n_elements(wavelengths)
speccovar=fltarr(nwvlh,nwvlh)
psi=fltarr(nwvlh,nwvlh)
lamarr=lambda
lam_med=median(lamarr,/even)

;now populate the matrix C to get errors

for i=0L,21 do begin
 for j=0L,21 do begin
  lamratio=(lamarr[i]-lamarr[j])/lam_med
  if i ne j then begin
  psi[i,j]=arho*exp(-0.5*(rho*lamratio/sigrho)^2.) +$
           alam*exp(-0.5*(lamratio/siglam)^2.) +$
           0*adel
  endif else begin

  psi[i,j]=arho+alam+adel
  ;psi[i,j]=1
  endelse

  speccovar[i,j]=(psi[i,j]*sqrt(eplanet_flux[i]^2.*eplanet_flux[j]^2.))  ; spectral covariance
 endfor
endfor

endif

;***

if ~keyword_set(library) then begin
library='montreal'
dirname='Montreal_Spectral_library'
endif

montrealdir=charis_path(pathname='montrealdir')
if library eq 'montreal' then begin 
librarydir=montrealdir
dirname='Montreal_Spectral_library'
endif

bonnefoydir=charis_path(pathname='bonnefoydir')
if library eq 'bonnefoy' then begin 
librarydir=bonnefoydir
dirname='bonnefoy_library'
endif

file_mkdir,dirname


if ~keyword_set(pick) then begin
openw,2,dirname+'/'+'fitoutcome.txt'
endif else begin
openw,2,dirname+'/'+'fitoutcomepick.txt'
endelse

if ~keyword_set(calcind) then begin
printf,2,'File Name ','Spectral Type ','Gravity Class', 'Chi-Sq'
endif else begin
printf,2,'File Name ','Spectral Type ','Gravity Class', 'Chi-Sq',' SpT Ind(4) ',' Grav Ind (2) '
endelse


;*****Read in String of Spectral Library Collection*****

;check
if ~keyword_set(pick) then begin
montreal_spec=file_search(librarydir+'*fits')
endif else begin
montreal_spec=dialog_pickfile(Title="Select the Model Spectra",/multiple_files,path=librarydir)
endelse

;read in file that lists the Montreal spectrum name, the spectral type and the gravity
readcol,montrealdir+'montreal_library_full.tsv',delimiter=string(9b),name,ra,dec,wave,yn,spt,format='a,a,a,a,a,a,a'

;a test: loop the montreal_spec names, match them to the names in the .tsv library, seems to work

mont_spt=strarr(n_elements(montreal_spec))
mont_gravclass=intarr(n_elements(montreal_spec))
mont_name=strarr(n_elements(montreal_spec))

;gravity class
;0 field
;1 intermediate
;2 low
;3 very low

name=strcompress(name,/remove_all)

print,n_elements(montreal_spec)

for ispec=0L,n_elements(montreal_spec) -1 do begin

;extract the root file name
your_path=strpos(montreal_spec[ispec],'/',/reverse_search)+1
fname=strmid(montreal_spec[ispec],your_path,strlen(montreal_spec[ispec])-1)
fname=(strsplit(fname,'.',/extract))[0]

for itxt=0L,n_elements(name) -1 do begin
;name[g]=strtrim(strcompress(name[g]))

mm=fname.contains(name[itxt])

if mm eq 1 then begin
mont_spt[ispec]=spt[itxt]
mont_name[ispec]=name[itxt]

;intermediate gravity
intgrav=mont_spt[ispec].contains('int')
if intgrav eq 1 then mont_gravclass[ispec] = 1

;low gravity
lowgrav=mont_spt[ispec].contains('low')
if lowgrav eq 1 then mont_gravclass[ispec]=2

;very low gravity
very=mont_spt[ispec].contains('very')
if very eq 1 then mont_gravclass[ispec]=3
endif

endfor



;now clean up the naming convention for the spectral type and gravity class
;...if the spectral type has an '<' or ':' in it, then remove ...

bb=strsplit(mont_spt[ispec],'<',/extract,count=ltcount)
if ltcount gt 0 then $
mont_spt[ispec]=strjoin(bb)

bb=strsplit(mont_spt[ispec],':',/extract,count=semicount)
if semicount gt 0 then mont_spt[ispec]=bb[0]
;mont_spt[ispec]=strjoin(bb)

;space
bb=strsplit(mont_spt[ispec],' ',/extract,count=spacecount)
if spacecount gt 0 then begin
;space_count+=1
;print,'spacecount',space_count
mont_spt[ispec]=bb[0]
endif


;now smooth spectrum to CHARIS resolution ...

if library eq 'montreal' or ~keyword_set(library) then begin
print,'fname is ',fname
print,montreal_spec[ispec]
inputmodel=readfits(montreal_spec[ispec],h1)
inputwvlh=inputmodel[*,0]
inputspec=inputmodel[*,1]
inputespec=inputmodel[*,2]

charis_convert_planetspec,models=montreal_spec[ispec],/filerot,modeldir='/',outputspec=outputspec_pred,outputerrspec=outputerrspec_pred,filtname=filtername0

endif else begin
;bonnefoy library

inputmodel=readfits(montreal_spec[ispec],h1)
inputwvlh=inputmodel[0,*]
inputspec=inputmodel[1,*]
;inputespec=inputmodel[2,*]
inputespec=0.001*inputspec
charis_convert_planetspec,models=fname+'.fits',modeldir=librarydir,outputspec=outputspec_pred,outputerrspec=outputerrspec_pred,filtname=filtername0

endelse

;convert from Flam to Fnu for model spectrum
;****for charis-converted model
clight=2.99792458d14  
conv_fact=((wavelength_planet[*]^2)/clight) ; this is in um*s
                                conv_fact*=1e4 ; now in A*s
                                conv_fact*=(1.0/10.0^(-23.0))
                                conv_fact*=1d3  ;from Jy to mJy
outputspec_pred*=conv_fact
outputerrspec_pred*=conv_fact

;****for raw model
conv_factinput=((inputwvlh[*]^2)/clight) ; this is in um*s
                                conv_factinput*=1e4 ; now in A*s
                                conv_factinput*=(1.0/10.0^(-23.0))
                                conv_factinput*=1d3  ;from Jy to mJy
inputspec*=conv_factinput
inputespec*=conv_factinput
;****

goto,skipoverthisplot
;window,3
set_plot,'x'
setcolors,/system_variables,/silent

plot,wavelength_planet,outputspec_pred,xrange=[1,2.5],linestyle=0,/nodata,yrange=[min(outputspec_pred),max(outputspec_pred)],ystyle=1,xstyle=1
oplot,wavelength_planet,outputspec_pred,linestyle=0,color=!blue,psym=-4,thick=3
al_legend,fname,/bottom,charsize=1.25,box=0
skipoverthisplot:


estmed=total(outputspec_pred[good]*planet_flux[good]/(eplanet_flux[good]^2.+outputerrspec_pred[good]^2.))/total(outputspec_pred[good]^2./(eplanet_flux[good]^2.+outputerrspec_pred[good]^2.))
;estmed=total(outputspec_pred[good]*planet_flux[good]/eplanet_flux[good]^2.)/total(outputspec_pred[good]^2./eplanet_flux[good]^2.)

;print,'planet flux',planet_flux[good]
;print,'eplan',eplanet_flux[good]
;print,'model',outputspec_pred[good]
;print,'emodel',outputerrspec_pred[good]
;print,estmed
;plot,wavelength_planet,outputspec_pred*estmed
;stop

;outputerrspec_pred[*]=outputspec_pred*0.01
baderror=where(outputerrspec_pred gt 0.5*outputspec_pred,nbaderror)

if ~keyword_set(planetcovar) then begin
chisqcalc,obs=planet_flux[good],eobs=eplanet_flux[good],model_in=outputspec_pred[good],err_model_in=outputerrspec_pred[good],alpha_in=estmed,chisqout=chisqout,ndof=ndof,alpha_out=alpha_out
endif else begin

chisqcalcwcovar,obs=planet_flux,eobs=eplanet_flux,model_in=outputspec_pred,err_model_in=outputerrspec_pred,alpha_in=estmed,chisqout=chisqout,ndof=ndof,alpha_out=alpha_out,speccovar=speccovar,good=good
endelse

print,'Fit Output is ...',fname,' ',mont_spt[ispec],' ',mont_gravclass[ispec],' ',chisqout/float(ngood-1)


;do you do the index calculations or not?

if ~keyword_set(calcind) then begin
if nbaderror gt 0 then begin
printf,2,fname,mont_spt[ispec],mont_gravclass[ispec],-999,format='(a,1x,a,1x,i3,1x,f8.3)'
endif else begin
printf,2,fname,mont_spt[ispec],mont_gravclass[ispec],chisqout/float(ngood-1),format='(a,1x,a,1x,i3,1x,f8.3)'
endelse

endif else begin

if ~keyword_set(justind) then begin
charis_calcsptgravind,fnu=outputspec_pred,outputspt=outputspt,outputgrav=outputgrav
endif else begin
charis_calcsptgravind,fnu=outputspec_pred,outputspt=outputspt,outputgrav=outputgrav,justind=justind
endelse

if nbaderror gt 0 then begin
printf,2,fname,mont_spt[ispec],mont_gravclass[ispec],-999,outputspt[0],outputspt[1],outputspt[2],outputspt[3],outputgrav[0],outputgrav[1],format='(a,1x,a,1x,i3,1x,f8.3,1x,6(f8.3,1x))'
endif else begin
printf,2,fname,mont_spt[ispec],mont_gravclass[ispec],chisqout/float(ngood-1),outputspt[0],outputspt[1],outputspt[2],outputspt[3],outputgrav[0],outputgrav[1],format='(a,1x,a,1x,i3,1x,f8.3,1x,6(f8.3,1x))'
endelse

endelse

outputspec_pred*=alpha_out
outputerrspec_pred*=alpha_out
inputspec*=alpha_out
inputespec*=alpha_out

;****Now report and plot the output

;window,4
set_plot,'x'
setcolors,/system_variables,/silent
plot,wavelength_planet,outputspec_pred,xrange=[1,2.5],linestyle=0,/nodata,yrange=[0.9*min(outputspec_pred),1.1*max(outputspec_pred)],ystyle=1,xstyle=1
;plot,wavelength_planet,outputspec_pred,xrange=[1,2.5],linestyle=0,/nodata,yrange=[min(outputspec_pred),max(outputspec_pred)],ystyle=1,xstyle=1
;oplot,wavelength_planet,alpha_out*outputspec_pred,linestyle=0,color=!blue,psym=-4,thick=3
;****the raw spectrum
oploterror,inputwvlh,inputspec,inputespec,linestyle=1,thick=2,/nohat,color=!yellow,errcolor=!gray,errthick=0.25
;****the charis-converted spectrum
oploterror,wavelength_planet,outputspec_pred,outputerrspec_pred,linestyle=0,color=!blue,psym=-4,thick=3
;**** the planet spectrum
oploterror,wavelength_planet,planet_flux,eplanet_flux,color=!green,psym=-3,thick=2
al_legend,fname,/bottom,charsize=1.25,box=0

if nbaderror gt 0 then begin
print,'BAD!!!
;stop
;wait,5
endif

if keyword_set(saveplot) then begin
set_plot,'ps'
device,filename=plotdir+fname+'.eps',/encapsulated
setcolors,/system_variables,/silent
plot,wavelength_planet,outputspec_pred,xrange=[1,2.5],linestyle=0,/nodata,yrange=[0.9*min(outputspec_pred),1.1*max(outputspec_pred)],ystyle=1,xstyle=1,$
xthick=5,ythick=5,charsize=1.25,xtitle='Wavelength (Microns)',ytitle='Flux Density (mJy)'
;plot,wavelength_planet,outputspec_pred,xrange=[1,2.5],linestyle=0,/nodata,yrange=[min(outputspec_pred),max(outputspec_pred)],ystyle=1,xstyle=1
;oplot,wavelength_planet,alpha_out*outputspec_pred,linestyle=0,color=!blue,psym=-4,thick=3
;****the raw spectrum
oploterror,inputwvlh,inputspec,inputespec,linestyle=1,thick=1.5,/nohat,color=!gray,errcolor=!gray,errthick=0.25
;oplot,inputwvlh,inputspec,linestyle=1,thick=3
;****the charis-converted spectrum
oploterror,wavelength_planet,outputspec_pred,outputerrspec_pred,linestyle=0,color=!green,psym=-4,thick=10,errthick=5
;oplot,wavelength_planet,outputspec_pred,linestyle=0,color=!green,psym=-4,thick=3
;**** the planet spectrum
oploterror,wavelength_planet,planet_flux,eplanet_flux,color=!blue,psym=-3,thick=10,errthick=5
al_legend,fname,charsize=0.75,charthick=2,box=0,position=[1.7,min(outputspec_pred)]
al_legend,textoidl('\chi^2 =')+string(sigfig(chisqout/(ngood-1),4)),/left,box=0,charsize=0.9
device,/close
endif

if keyword_set(delay) then wait,2

endfor
close,/all

;now for your target
charis_calcsptgravind,fnu=planet_flux,outputspt=outputsptplan,outputgrav=outputgravplan
print,outputsptplan
print,outputgravplan


end
