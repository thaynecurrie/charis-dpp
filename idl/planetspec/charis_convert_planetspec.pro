pro charis_convert_planetspec,models=models,mode=mode,filerot=filerot,spex=spex,usetext=usetext,modeldir=modeldir,plotql=plotql,outputspec=outputspec,outputerrspec=outputerrspec,filtname=filtname0

;models - name of models
;mode - lowres, J, H, or K

if ~keyword_set(filtname) then filtname0 = 'lowres'

;defines wavelength array
get_charis_wvlh,dum,wavelengths,manual=filtname0
lambda=1d-3*wavelengths

if ~keyword_set(modeldir) then begin
;modeldir='~/idl_tools/ADI_dl/charisred/models/'
modeldir=charis_path(pathname='modeldir')
;if keyword_set(spex) then modeldir='~/Research/Planets/DI/spex/'
if keyword_set(spex) then modeldir=charis_path(pathname='spexdir')
endif

if ~keyword_set(models) then begin
models=dialog_pickfile(Title="Select Planet Model Spectra",PATH=modeldir)

endif else begin
models=modeldir+models
endelse

nmodels=n_elements(models)

ftype=0
if keyword_set(usetext) then ftype=1

for i=0L,nmodels - 1 do begin


if ftype lt 1 then begin
aaa=readfits(models[i],haaa)

if ~keyword_set(filerot) then begin
mwvlh=aaa[0,*]
mflux=aaa[1,*]
endif else begin
mwvlh=aaa[*,0]
mflux=aaa[*,1]
endelse
endif else begin
readcol,models[i],mwvlh,mflux

endelse

;plot,mwvlh,mflux
if ~keyword_set(filerot) then begin
model_lowres=charis_convert_planetspec_to_lowres_calculation(lambda,modelspec=models[i],ftype=ftype,model_errest=model_errest,filtname=filtname0)
;model_lowres=charis_convert_planetspec_to_lowres_calculation(lambda,modelspec=models[i],ftype=ftype,model_errest=model_errest)
endif else begin
model_lowres=charis_convert_planetspec_to_lowres_calculation(lambda,modelspec=models[i],ftype=ftype,/filerot,model_errest=model_errest,filtname=filtname0)
;model_lowres=charis_convert_planetspec_to_lowres_calculation(lambda,modelspec=models[i],ftype=ftype,/filerot,model_errest=model_errest)
endelse

;ploterror,lambda,model_lowres,model_errest


;name of the spectrum
subdir=strpos(models[i],'/',/reverse_search)+1
lenpath=strlen(models[i])
modelname=strmid(models[i],subdir,lenpath)

if keyword_set(plotql) then begin
window,0
set_plot,'x'
setcolors,/system_variables
plot,mwvlh,mflux,yrange=[min(model_lowres),max(model_lowres)],xrange=[1,2.5],linestyle=0,/nodata
oplot,mwvlh,mflux,linestyle=0,color=!green
oplot,lambda,model_lowres,linestyle=0,psym=-4,color=!blue,thick=3
al_legend,modelname,/bottom,charsize=1.25,box=0

window,1
setcolors,/system_variables
plot,mwvlh,mflux*mwvlh^2.,yrange=[min(model_lowres*lambda^2.),max(model_lowres*lambda^2.)],xrange=[1,2.5],linestyle=0,/nodata
oplot,mwvlh,mflux*mwvlh^2.,linestyle=0,color=!green
oplot,lambda,model_lowres*lambda^2.,linestyle=0,color=!blue,psym=-4,thick=3
al_legend,modelname,/bottom,charsize=1.25,box=0

endif

;wait,1
endfor

;if you are using this program as a subroutine and want an output spectrum for analysis, then do the following

outputspec=model_lowres
;if keyword_set(outputerrspec) then outputerrspec=model_errest
outputerrspec=model_errest

end




