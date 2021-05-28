function  charis_convert_planetspec_to_lowres_calculation2,model_flux,model_wavelengths,lambda_out,filtname=filtname

;function: charis_convert_planetspec_to_lowres_calculation2
; call this function to take a model or empirical input spectrum and degrade it to CHARIS resolution.  


;use this in code for atmosphere model comparisons where you already have everything written out to a file. 

; model input as ergs/s/cm^2/A (or, more generally, F_lambda units

;***Spectrophotometric calibration subroutine, following GPI DRP v1.4 methods
;-takes in planet/model spectrum
;-bins and interpolates at resolution of CHARIS
;-return the binned/interpolated CHARIS spectrum.

;********
;Bin the spectrum from full resolution to matching CHARIS's resolution

;Hard-wire for now
;**** resolution
;filter = 'Broadband'
  ;average spectral resolution per filter
  case filtname of 
   'lowres': specresolution=20.
   'J': specresolution=75.
   'H': specresolution=75.
   'K': specresolution=75.
  endcase

;****array of wavelengths
get_charis_wvlh,dum,lambda,manual=filtname
lambda*=1d-3   ;nm to microns

;*******
  ;plot,model_flux
  ;take wavelength spacing of model spectrum, compare it to the CHARIS wavelength array.  Degrade spectrum resolution accordingly.
  fwhmloc=VALUE_LOCATE(model_wavelengths,[(lambda[0]),(lambda[1])])
  if max(fwhmloc) lt 0 then fwhmloc=VALUE_LOCATE(model_wavelengths,[(lambda[n_elements(lambda)-2]),(lambda[n_elements(lambda)-1])])
  fwhm_fr=float(fwhmloc[1]-fwhmloc[0])
  gaus=PSF_GAUSSIAN(Npixel=3.*fwhm_fr, FWHM=fwhm_fr, NDIMEN =1, /NORMAL )
;smooth the model spectrum to charis resolution
  good=where(model_flux gt 0)
  charis_model_flux_smooth=CONVOL( reform(model_flux[good]), gaus , /EDGE_TRUNCATE )
 

;smooth the transmission function to charis resolution
;-since the filter is usually a subset of the full JHK CHARIS low res coverage, use different approach

  ;now interpolate the smoothed model flux onto the CHARIS grid
  charis_model_flux=interpol(charis_model_flux_smooth,model_wavelengths[good],lambda)   

  ;return the model spectrum: binned/smoothed, over the charis wavelengths, and scaled to the desired contrast
  return,charis_model_flux

  end
