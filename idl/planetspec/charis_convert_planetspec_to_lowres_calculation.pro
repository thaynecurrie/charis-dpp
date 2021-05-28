function  charis_convert_planetspec_to_lowres_calculation,lambda,modelspec=modelspec,model_errest=model_errest,$
usefits=usefits,usetxt=usetxt,ftype=ftype,filerot=filerot,filtname=filtname

;function: charis_convert_planetspec_to_lowres_calculation
; call this function to take a model or empirical input spectrum and degrade it to CHARIS resolution.   Currently works only for broadband
;keywords: 
;-pri_header & ext_header: the primary and first extension of the charis fits header.  Use to pull information about the star's brightness
;-modelspec - the model spectrum you wish to use.  for now, must be in units of ergs/s/cm^2/A
;-contrast - the H band contrast between the star and the companion

;*******CAUTION******* FOR NOW, THIS IS HARDWIRED FOR ...
; data cube input as mJy
; model input as ergs/s/cm^2/A (or, more generally, F_lambda units

;***Spectrophotometric calibration subroutine, following GPI DRP v1.4 methods
;-takes in planet/model spectrum
;-bins and interpolates at resolution of CHARIS
;-return the binned/interpolated CHARIS spectrum.

;*****Spectrum of the Planet
;ftype
;1 - reading in an ascii text file
;2 - reading in a fits file

if ftype ge 1 then begin
;read in the planet/model spectrum
readcol,modelspec,model_wavelengths,model_flux,model_eflux   ;assume units are ergs/s/cm^2/A for now and wavelength is in angstroms
endif else begin
aaa=readfits(modelspec,haaa)
sz=size(aaa,/dim)
if ~keyword_set(filerot) then begin
model_wavelengths=aaa[0,*]
model_flux=aaa[1,*]
if sz[1] eq 3 then begin
model_eflux=aaa[2,*]
endif else begin
model_eflux=0.1*model_flux
endelse
endif else begin

;if you have a normal file then do this ...
model_wavelengths=aaa[*,0]
model_flux=aaa[*,1]
model_eflux=aaa[*,2]

endelse
endelse

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

  charis_model_err_smooth=CONVOL( reform(model_eflux[good]), gaus, /EDGE_TRUNCATE )

  ;assume very coarsely that S/N ~ sqrt(S) and that locally the S ~ constant with channel
  charis_model_err_smooth=charis_model_err_smooth/sqrt(fwhm_fr)
  charis_model_eflux=interpol(charis_model_err_smooth,model_wavelengths[good],lambda)   
  model_errest=charis_model_eflux
  
  ;--convert to F_nu units (mJy)

  ;****assume star is in mJy and model in ergs/s/cm^2/A for now***
                          c=2.99792458d14                     ; um / s
                                conv_fact=((lambda[*]^2)/c) ; this is in um*s
                                conv_fact*=1e4 ; now in A*s
                                conv_fact*=(1.0/10.0^(-23.0))
                                conv_fact*=1d3  ;from Jy to mJy
  goto,skipoverme
  ;the model planet flux: smoothed/binned over charis wavelengths and converted to mJy
   charis_model_flux_fnu0=conv_fact*charis_model_flux

  ;- integrate model over response function for 'model photometry'
   charis_model_raw_photflux_pseudo=int_tabulated(lambda,charis_model_flux_fnu0*charis_std_filter_response_bin,/sort)/int_tabulated(lambda,charis_std_filter_response_bin,/sort)

  ;integrate starflux over filter response for 'star photometry'
   charis_star_photflux_pseudo_fnu=int_tabulated(lambda,charis_star_spec*charis_std_filter_response_bin,/sort)/int_tabulated(lambda,charis_std_filter_response_bin,/sort)

  ;Now, use this information to get the additional mult. factor to make model/star = contrast
   contrast_conv_fac=10^(double(contrast))*charis_star_photflux_pseudo_fnu/charis_model_raw_photflux_pseudo

  ;And finally, multiply the model spectrum to make it the right contrast
   charis_model_flux_cal_fnu=charis_model_flux_fnu0*contrast_conv_fac
;*10^(contrast)
  
  ;return the stellar spectrum (as a diagnostic) 
  starspec=charis_star_spec
  skipoverme:

  ;return the model spectrum: binned/smoothed, over the charis wavelengths, and scaled to the desired contrast
  return,charis_model_flux

  end
