function charis_planet_photometric_calibration_calculation,pri_header,ext_header,lambda,$
modelspec=modelspec,filtname=filtname,contrast=contrast,starspec=starspec
;mag=star_mag,units=units,$
;diagnostic=diagnostic

;function: charis_planet_photometric_calibration_calculation
; call this function to figure out how bright a model planet is in data cube units for a given requested contrast (as integrated over a bandpass)
;keywords: 
;-pri_header & ext_header: the primary and first extension of the charis fits header.  Use to pull information about the star's brightness
;-modelspec - the model spectrum you wish to use.  for now, must be in units of ergs/s/cm^2/A
;-contrast - the H band contrast between the star and the companion

;*******CAUTION******* FOR NOW, THIS IS HARDWIRED FOR ...
; data cube input as mJy
; model input as ergs/s/cm^2/A (or, more generally, F_lambda units
; + only tested for broadband mode

;***Spectrophotometric calibration subroutine, following GPI DRP v1.4 methods
;-reads in header, wavelength array, spectrum file, and requested log(contrast) between star and companion
;-takes in planet/model spectrum
;-bins and interpolates at resolution of CHARIS
;-bin and interpolate filter transmission to CHARIS lambda-scale
;-convolve with filter function to get pseudo-photometry in H band.
;-grab flux from fits header for star
;-determine multiplicative factor that gives pseudo-photometry*factor=star*contrast
;-multiply the binned/interpolated CHARIS spectrum by this factor
;-return the calibrated, binned/interpolated CHARIS spectrum.

;****Filter Responses****
;standard filter response functions subdirectory
;we are going to use the MKO H band filter for calibration

filtresponsedir=charis_path(pathname='filtresponsedir')
;filtresponsedir='~/idl_tools/ADI_dl/charisred/tools/filters/filter_response/'

case filtname of
   'J': begin
    charis_std_filter='Jband'
    end
   'H': begin
     charis_std_filter='Hband'
    end
   'K': begin
     charis_std_filter='Ksband'
    end
    'lowres':begin
     charis_std_filter='Hband'
     end
    'broadband':begin
     charis_std_filter='Hband'
     end
endcase

;charis_std_filter='Hband'
readcol,filtresponsedir+charis_std_filter+'.dat',charis_std_filter_wvlh,charis_std_filter_response,/silent
;readcol,filtresponsedir+charis_std_filter+'.dat',charis_std_filter_wvlh,charis_std_filter_response
;convert to microns
charis_std_filter_wvlh*=1d-3

good=where(charis_std_filter_response gt 10)

charis_std_filter_wvlh=charis_std_filter_wvlh[good]
charis_std_filter_response=charis_std_filter_response[good]

charis_std_filter_response=charis_std_filter_response(sort(charis_std_filter_wvlh))
charis_std_filter_wvlh=charis_std_filter_wvlh(sort(charis_std_filter_wvlh))

;*****Spectrum of the Star: Brightness of Star in each CHARIS channel
;-first, identify the number of wavelength channels.  It should be 22 for the broadband mode
nlambda=sxpar(ext_header,'naxis3')
charis_star_spec=fltarr(nlambda)
for ilambda=0L,nlambda-1 do begin
 charis_star_spec[ilambda]=sxpar(ext_header,'FSTAR_'+strtrim(string(ilambda),2))
; print,charis_star_spec[ilambda]
endfor

;*****Spectrum of the Planet
;1. If it's a fits file, then assume [*,0] is wavelength and [*,1] is flux density in F-lambda units
;read in the planet/model spectrum
modelspecname=strsplit(modelspec,'.',/extract)
;if its a fits file then do this
if n_elements(modelspecname) eq 1 then begin
;assume its a .txt file that can be read by readcol
readcol,modelspec,model_wavelengths,model_flux,/silent   ;assume units are ergs/s/cm^2/A for now and wavelength is in angstroms
endif
if modelspecname[1] eq 'fits' then begin
a=readfits(modelspec,/silent)
model_wavelengths=a[*,0]
model_flux=a[*,1]
;bad=where(model_flux lt 0)
;model_flux[bad]=0
endif else begin
;assume its a .txt file that can be read by readcol
readcol,modelspec,model_wavelengths,model_flux,/silent   ;assume units are ergs/s/cm^2/A for now and wavelength is in angstroms
endelse

;*****CAUTION***** For now, only tested with broadband.
;Bin the spectrum from full resolution to matching CHARIS's resolution
 filter=sxpar(pri_header,'FILTNAME')
 dloglam=sxpar(pri_header,'DLOGLAM')
 ;specresolution=median(lambda/(lambda-shift(lambda,1)))
 ;print,filter

  ;average spectral resolution per filter
  case filter of 
   'Broadband': specresolution=20.
   'ND': specresolution=20.
   'J': specresolution=75.
   'H': specresolution=75.
   'K': specresolution=75.
    else:specresolution=20.
  endcase

  ;take wavelength spacing of model spectrum, compare it to the CHARIS wavelength array.  Degrade spectrum resolution accordingly.
  fwhmloc=VALUE_LOCATE(model_wavelengths,[(lambda[0]),(lambda[1])])
  fwhm_fr=float(fwhmloc[1]-fwhmloc[0])
  gaus=PSF_GAUSSIAN(Npixel=3.*fwhm_fr, FWHM=fwhm_fr, NDIMEN =1, /NORMAL )
;smooth the model spectrum to charis resolution
  charis_model_flux_smooth=CONVOL( reform(model_flux), gaus , /EDGE_TRUNCATE,/nan)

;smooth the transmission function to charis resolution
;-since the filter is usually a subset of the full JHK CHARIS low res coverage, use different approach

  charis_avg_pseudores=median(lambda/(lambda-shift(lambda,1)),/even)  
  filter_avg_pseudores=median(charis_std_filter_wvlh/(charis_std_filter_wvlh-shift(charis_std_filter_wvlh,1)),/even)
  ;assume that the filter spacing is finer so charis_avg_pseudores will always be smaller.  code will crash otherwise
  if filter_avg_pseudores lt charis_avg_pseudores then begin
    print,'error!'
  stop
  endif
  if filter_avg_pseudores eq charis_avg_pseudores then fwhm_fr=1
  if filter_avg_pseudores gt charis_avg_pseudores then fwhm_fr=float(filter_avg_pseudores)/float(charis_avg_pseudores)

  gaus=PSF_GAUSSIAN(Npixel=3.*fwhm_fr, FWHM=fwhm_fr, NDIMEN =1, /NORMAL )
  charis_std_filter_response_smooth=CONVOL( reform(charis_std_filter_response),gaus,/EDGE_TRUNCATE)

  ;interpolate filter response and model flux to charis wavelength array
  charis_std_filter_response_bin=interpol(charis_std_filter_response_smooth,charis_std_filter_wvlh,lambda)

  ;set the bad values (those outside the wavelength range) to zero
  badresponse=where(charis_std_filter_response_bin lt 0)
  charis_std_filter_response_bin[badresponse]=0
  
  ;now interpolate the smoothed model flux onto the CHARIS grid
  charis_model_flux=interpol(charis_model_flux_smooth,model_wavelengths,lambda)   

  ;integrate model over filter response
  ;-we are going to first switch the model/input spectrum to fnu units and integrate over lambda (int(fnu*phi_lambda*dlambda)).   Technically, this does not give us a flux density 
  ;-- since it should be int(flambda*phi_lambda*dlambda) = int(fnu*phi_nu*dnu)
  ;-- but since we are doing int(fnu*phi_nu*dnu) for both the star and the model, this should give us a good estimate of the mult. factor for the model to input a planet of a given 
  ;-- contrast (over a filter)

  ;****assume star is in mJy and model in ergs/s/cm^2/A for now***
                          c=2.99792458d14                     ; um / s
                                conv_fact=((lambda[*]^2)/c) ; this is in um*s
                                conv_fact*=1e4 ; now in A*s
                                conv_fact*=(1.0/10.0^(-23.0))
                                conv_fact*=1d3  ;from Jy to mJy

  ;the model planet flux: smoothed/binned over charis wavelengths and converted to mJy
   charis_model_flux_fnu0=conv_fact*charis_model_flux
;   print,lambda,charis_model_flux_fnu0
;   for ill=0L,n_elements(lambda)-1 do begin
;   print,lambda[ill],charis_model_flux_fnu0[ill]
;   endfor

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
  ;return the model spectrum: binned/smoothed, over the charis wavelengths, and scaled to the desired contrast
  return,charis_model_flux_cal_fnu

  end
