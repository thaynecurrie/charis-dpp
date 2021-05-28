function  charis_convert_planetphot_to_lowres_calculation,model_flux,model_wavelengths_init,lambda_out,filtname=filtname

;function: charis_convert_planetphot_to_lowres_calculation
; call this function to take a model or empirical input spectrum and get a photometric prediction

;use this in code for atmosphere model comparisons where you already have everything written out to a file. 

; model input as ergs/s/cm^2/A (or, more generally, F_lambda units

;***Spectrophotometric calibration subroutine, following GPI DRP v1.4 methods
;-takes in planet/model spectrum
;-bins and interpolates at resolution of CHARIS
;-return the binned/interpolated CHARIS spectrum.

;********

model_wavelengths=model_wavelengths_init*1d3

filtresponsedir=charis_path(pathname='filtresponsedir')

;find the filter file

case filtname of 

   'J': begin
    ;charis_std_filter='j2m'
    charis_std_filter='Jband'
    end
   'H': begin
     charis_std_filter='Hband'
    end
   'K': begin
     charis_std_filter='Ksband'
    end

;others
   'Y' : begin
     charis_std_filter='Ybandest'
    end
   
   'Lp': begin
     charis_std_filter='Lpband'
    end
   
   'Ms': begin
     charis_std_filter='Msband' 
    end
   
    else: begin
    print,'Cannot find filter'
    stop
    end 
endcase


readcol,filtresponsedir+charis_std_filter+'.dat',charis_std_filter_wvlh,charis_std_filter_response

;plot,charis_std_filter_wvlh,charis_std_filter_response

good=where(charis_std_filter_response gt 10)

charis_std_filter_wvlh=charis_std_filter_wvlh[good]
charis_std_filter_response=charis_std_filter_response[good]

charis_std_filter_response=charis_std_filter_response(sort(charis_std_filter_wvlh))
charis_std_filter_wvlh=charis_std_filter_wvlh(sort(charis_std_filter_wvlh))

nfilterwvlh=n_elements(charis_std_filter_wvlh)
modelwvlhrange=where(model_wavelengths ge charis_std_filter_wvlh[0] and model_wavelengths le charis_std_filter_wvlh[nfilterwvlh-1],nmodelwvlhrange)

;model_wavelengths=model_wavelengths[modelwvlhrange]
;model_flux=model_flux[modelwvlhrange]



;print,charis_std_filter,' compare numbers',nmodelwvlhrange,nfilterwvlh 
;print,min(model_wavelengths),max(model_wavelengths),charis_std_filter_wvlh[0],charis_std_filter_wvlh[nfilterwvlh-1]

if nmodelwvlhrange gt nfilterwvlh then begin

model_flux_f=interpol(model_flux,model_wavelengths,charis_std_filter_wvlh)

;if charis_std_filter eq 'Jband' then begin
;plot,model_wavelengths,model_flux_f
;stop
;endif

num_model=int_tabulated(charis_std_filter_wvlh,charis_std_filter_wvlh*double(charis_std_filter_response)*model_flux_f,/double,/sort)
denom_model=int_tabulated(charis_std_filter_wvlh,charis_std_filter_wvlh*double(charis_std_filter_response),/double)

;set_plot,'x'
;window,3
;plot,charis_std_filter_wvlh,charis_std_filter_response

;print,'OPTION 1'
;wait,1

endif else begin

;print,'OPTION 2'
;wait,1

;we get an interpolation error otherwise
model_wavelengthsf=model_wavelengths[modelwvlhrange]
model_fluxf=model_flux[modelwvlhrange]
charis_std_filter_response_f=interpol(charis_std_filter_response,charis_std_filter_wvlh,model_wavelengthsf,/spline)
num_model=int_tabulated(model_wavelengthsf,model_wavelengthsf*double(charis_std_filter_response_f)*model_fluxf,/double,/sort)
denom_model=int_tabulated(model_wavelengthsf,model_wavelengthsf*double(charis_std_filter_response_f),/double)

;if charis_std_filter eq 'Jband' then begin
;set_plot,'x'
;window,3
;plot,charis_std_filter_wvlh,charis_std_filter_response
;plot,model_wavelengths,charis_std_filter_response_f
;oplot,charis_std_filter_wvlh,charis_std_filter_response,linestyle=1
;stop
;endif

endelse

model_flux_out=num_model/denom_model

;plot,model_wavelengths*1d-3,smooth(model_flux,20),xrange=[0.8,1.2]
;setcolors,/system_variables
;oplot,[1.03],[model_flux_out],psym=4,symsize=10,color=!green,thick=10

;set_plot,'x'
;window,0
;plot,model_wavelengths*1d-3,model_flux,xrange=[0.8*median(charis_std_filter_wvlh)*1d-3,2]
;plot,model_wavelengths*1d-3,model_flux,xrange=[0.8*median(charis_std_filter_wvlh)*1d-3,1d-3*1.2*median(charis_std_filter_wvlh)]
;setcolors,/system_variables
;oplot,[median(charis_std_filter_wvlh)*1d-3],[model_flux_out],psym=4,symsize=10,color=!green,thick=10

;if charis_std_filter eq 'Jband' then stop
  return,model_flux_out

close,/all
end
