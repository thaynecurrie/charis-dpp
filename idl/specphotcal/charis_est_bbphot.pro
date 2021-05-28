pro charis_est_bbphot,inputspectrumfile,inputunits

;Need a better accounting of the error propagation.   You are probably overestimating errors slightly.

if N_PARAMS() eq 0 then  begin
print, "charis_est_bb_phot,inputspectrumfile,inputunits"
print, "Takes a spectrum extracted from a CHARIS data cube and then estimates JHK photometry"
print, "inputunits are 0-mJy, 1-Jy, 2-W/m^2/um, 3-ergs/s/cm^2/A, 4-ergs/s/cm^2/Hz"
print, "default is 0 (mJy)"
print, "CAUTION: right now the units are hardwired to mJy!!! (will add additional unit conversion later)"
goto,skiptoend
endif 

if ~keyword_set(inputunits) then inputunits = 0

;read in the extracted spectrum
readcol,inputspectrumfile,lambda,flux,eflux

filtresponsedir=charis_path(pathname='filtresponsedir')

;Decision Tree for filter...
jflag=0 & hflag=0 & kflag = 0

if min(lambda) lt 1.17 and max(lambda) gt 2.3 then begin
jflag=1 & hflag=1 & kflag = 1 
endif

if min(lambda) lt 1.17 and max(lambda) lt 1.4 then begin
jflag=1 & hflag=0 & kflag = 0
endif

if min(lambda) gt 1.4 and min(lambda) lt 1.6 then begin
jflag=0 & hflag=1 & kflag = 0
endif 

if min(lambda) gt 1.85 then begin
jflag=0 & hflag=0 & kflag = 1
endif

filt_j='Jband.dat'
filt_h='Hband.dat'
filt_k='Ksband.dat'

if jflag gt 0 then begin
;Jband
readcol,filtresponsedir+filt_j,charis_std_filter_wvlh_j,charis_std_filter_response_j
good_j=where(charis_std_filter_response_j gt 10)
charis_std_filter_wvlh_j=charis_std_filter_wvlh_j[good_j]
charis_std_filter_response_j=charis_std_filter_response_j[good_j]
charis_std_filter_response_j=charis_std_filter_response_j(sort(charis_std_filter_wvlh_j))
charis_std_filter_wvlh_j=charis_std_filter_wvlh_j(sort(charis_std_filter_wvlh_j))

;rebin
fwhmloc=VALUE_LOCATE(charis_std_filter_wvlh_j*1e-3,[(lambda[0]),(lambda[1])])
fwhm=float(fwhmloc[1]-fwhmloc[0])
gaus=PSF_GAUSSIAN(Npixel=3.*fwhm, FWHM=fwhm, NDIMEN =1, /NORMAL )
charis_std_filter_response_jf=CONVOL( reform(charis_std_filter_response_j), gaus , /EDGE_TRUNCATE )
charis_std_filter_response_jf=interpol(charis_std_filter_response_jf,charis_std_filter_wvlh_j*1e-3,lambda)

endif

if hflag gt 0 then begin
;Hband
readcol,filtresponsedir+filt_h,charis_std_filter_wvlh_h,charis_std_filter_response_h
good_h=where(charis_std_filter_response_h gt 10)
charis_std_filter_wvlh_h=charis_std_filter_wvlh_h[good_h]
charis_std_filter_response_h=charis_std_filter_response_h[good_h]
charis_std_filter_response_h=charis_std_filter_response_h(sort(charis_std_filter_wvlh_h))
charis_std_filter_wvlh_h=charis_std_filter_wvlh_h(sort(charis_std_filter_wvlh_h))
;rebin
fwhmloc=VALUE_LOCATE(charis_std_filter_wvlh_h*1e-3,[(lambda[7]),(lambda[8])])
fwhm=float(fwhmloc[1]-fwhmloc[0])
gaus=PSF_GAUSSIAN(Npixel=3.*fwhm, FWHM=fwhm, NDIMEN =1, /NORMAL )
charis_std_filter_response_hf=CONVOL( reform(charis_std_filter_response_h), gaus , /EDGE_TRUNCATE )
charis_std_filter_response_hf=interpol(charis_std_filter_response_hf,charis_std_filter_wvlh_h*1e-3,lambda)

endif

if kflag gt 0 then begin
;Ksband
readcol,filtresponsedir+filt_k,charis_std_filter_wvlh_k,charis_std_filter_response_k
good_k=where(charis_std_filter_response_k gt 10)
charis_std_filter_wvlh_k=charis_std_filter_wvlh_k[good_k]
charis_std_filter_response_k=charis_std_filter_response_k[good_k]
charis_std_filter_response_k=charis_std_filter_response_k(sort(charis_std_filter_wvlh_k))
charis_std_filter_wvlh_k=charis_std_filter_wvlh_k(sort(charis_std_filter_wvlh_k))
;rebin
if (jflag gt 0 and hflag gt 0) then begin
fwhmloc=VALUE_LOCATE(charis_std_filter_wvlh_k*1e-3,[(lambda[16]),(lambda[17])])
endif else begin
fwhmloc=VALUE_LOCATE(charis_std_filter_wvlh_k*1e-3,[(lambda[8]),(lambda[9])])
endelse
fwhm=float(fwhmloc[1]-fwhmloc[0])
gaus=PSF_GAUSSIAN(Npixel=3.*fwhm, FWHM=fwhm, NDIMEN =1, /NORMAL )
charis_std_filter_response_kf=CONVOL( reform(charis_std_filter_response_k), gaus , /EDGE_TRUNCATE )
charis_std_filter_response_kf=interpol(charis_std_filter_response_kf,charis_std_filter_wvlh_k*1e-3,lambda)
endif

if keyword_set(diagnostic) then begin
plot,charis_std_filter_wvlh_j*1e-3,charis_std_filter_response_j,xrange=[1,2]
oplot,charis_std_filter_wvlh_h*1e-3,charis_std_filter_response_h
oplot,lambda,charis_std_filter_response_hf,linestyle=1,psym=-4
oplot,lambda,charis_std_filter_response_jf,linestyle=1,psym=-4
oplot,lambda,charis_std_filter_response_kf,linestyle=1,psym=-4

endif
;stop

case inputunits of 
 0: begin
    c=2.99792458d14 
    conv_fact=1d-3 ;from mJy to Jy
    conv_fact*=(1.0/10.0^(23.0)) ;from Jy to ...
    conv_fact*=1e-4 ; 
    conv_fact*=(c/(lambda[*]^2.))  ; now to ergs/s/cm^2/A

    flux_orig=flux
    eflux_orig=eflux

    flux*=conv_fact
    eflux*=conv_fact   
    end
  3: begin
    c=2.99792458d14 
    conv_fact=1.0 
    flux_orig=flux
    eflux_orig=eflux
    flux*=conv_fact
    eflux*=conv_fact
    end
endcase
;now convolve/sum
;    help,charis_std_filter_response_jf,flux
    print,flux

    if jflag gt 0 then begin
    good_j=where(charis_std_filter_response_jf gt 10)

    num_j=int_tabulated(lambda[good_j],lambda[good_j]*double(charis_std_filter_response_jf[good_j])*flux[good_j],/double,/sort)
    denom_j=int_tabulated(lambda[good_j],lambda[good_j]*charis_std_filter_response_jf[good_j],/double,/sort)
    num_ej=int_tabulated(lambda[good_j],lambda[good_j]*double(charis_std_filter_response_jf[good_j])*eflux[good_j],/double,/sort)
    num_jl=int_tabulated(lambda[good_j],double(charis_std_filter_response_jf[good_j])*lambda[good_j],/double,/sort)
    denom_jl=int_tabulated(lambda[good_j],charis_std_filter_response_jf[good_j],/double,/sort)
    jflux=num_j/denom_j
    ejflux=num_ej/denom_j
    effj=num_jl/denom_jl

    endif else begin
    jflux=1d20
    ejflux=1d20
    effj=1.25
    endelse

    if hflag gt 0 then begin
    good_h=where(charis_std_filter_response_hf gt 10)
    num_h=int_tabulated(lambda[good_h],lambda[good_h]*charis_std_filter_response_hf[good_h]*flux[good_h],/double)
    denom_h=int_tabulated(lambda[good_h],lambda[good_h]*charis_std_filter_response_hf[good_h],/double)
    num_eh=int_tabulated(lambda[good_h],lambda[good_h]*charis_std_filter_response_hf[good_h]*eflux[good_h],/double)
    num_hl=int_tabulated(lambda[good_h],charis_std_filter_response_hf[good_h]*lambda[good_h],/double)
    denom_hl=int_tabulated(lambda[good_h],charis_std_filter_response_hf[good_h],/double)
    hflux=num_h/denom_h
    ehflux=num_eh/denom_h
    effh=num_hl/denom_hl
    endif else begin
    hflux=1d20
    ehflux=1d20
    effh=1.65
    endelse
   
    if kflag gt 0 then begin 
    good_k=where(charis_std_filter_response_kf gt 10)
    num_k=int_tabulated(lambda[good_k],lambda[good_k]*charis_std_filter_response_kf[good_k]*flux[good_k],/double)
    denom_k=int_tabulated(lambda[good_k],lambda[good_k]*charis_std_filter_response_kf[good_k],/double)
    num_ek=int_tabulated(lambda[good_k],lambda[good_k]*charis_std_filter_response_kf[good_k]*eflux[good_k],/double)
    num_kl=int_tabulated(lambda[good_k],charis_std_filter_response_kf[good_k]*lambda[good_k],/double)
    denom_kl=int_tabulated(lambda[good_k],charis_std_filter_response_kf[good_k],/double)
    kflux=num_k/denom_k
    ekflux=num_ek/denom_k
    effk=num_kl/denom_kl
    endif else begin
    kflux=1d20
    ekflux=1d20
    effk=2.19
    endelse

    ;print,'jhkflux is ',num_j/denom_j,num_h/denom_h,num_k/denom_k
    if jflag gt 0 then print,'J = ',jflux,' +/- ',ejflux
    if hflag gt 0 then print,'H = ',hflux,' +/- ',ehflux
    if kflag gt 0 then print,'K = ',kflux,' +/- ',ekflux
    

    flags=[jflag,hflag,kflag]
    fluxes=[jflux,hflux,kflux]
    efluxes=[ejflux,ehflux,ekflux]
    waves=[1.25,1.65,2.19]
    effwaves=[effj,effh,effk]
    goodflag=where(flags gt 0)
    plot,lambda,flux
    ;oplot,[1.25,1.65,2.19],[num_j/denom_j,num_h/denom_h,num_k/denom_k],psym=4
    oplot,waves[goodflag],fluxes[goodflag]

;now convert back to 
  
    print,'eff wvlhs',effj,effh,effk
    fluxdensity_filt_flam=fluxes
;[num_j/denom_j,num_h/denom_h,num_k/denom_k]
    efluxdensity_filt_flam=efluxes
;[num_ej/denom_j,num_eh/denom_h,num_ek/denom_k]
    effwvlhs=effwaves
;[num_jl/denom_jl,num_hl/denom_hl,num_kl/denom_kl]

    ;now convert back to mJy
    conv_fact=((effwvlhs[*]^2)/c)
     conv_fact*=1e4 ; now in A*s
    conv_fact*=(1.0/10.0^(-23.0))
    conv_fact*=1d3  ;from Jy to mJy
    fluxdensity_filt_mjy=fluxdensity_filt_flam*conv_fact
    efluxdensity_filt_mjy=efluxdensity_filt_flam*conv_fact

    plot,lambda,flux_orig
    oplot,effwvlhs[goodflag],fluxdensity_filt_mjy[goodflag],psym=4
    print,'Flux Density ',fluxdensity_filt_mjy[goodflag]
    print,'eFlux Density',efluxdensity_filt_mjy[goodflag]
    fluxzero=[1560,1040,670]*1d3
    print,-2.5*alog10(fluxdensity_filt_mjy[goodflag]*1./fluxzero[goodflag]),1.0857*(efluxdensity_filt_flam[goodflag]/fluxdensity_filt_flam[goodflag])

skiptoend:
end
