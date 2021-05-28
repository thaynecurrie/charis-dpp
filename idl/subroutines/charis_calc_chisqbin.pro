pro charis_calc_chisqbin,chisq_map,fwhm,chisq_bin_out,nbins,chisqmap_bin=chisqmap_bin

;read in the chi-squared map.
;rebin by FWHM
;calc total chisquared
; divide by n_bin

sz=size(chisq_map,/dim)
print,sz[0]

chisq_bin_map=congrid(chisq_map,sz[0]/fwhm,sz[1]/fwhm)

roi_bin=where(chisq_bin_map ne 0 and finite(chisq_bin_map) ne 0,nroi_bin)

;reduced chi squared

chisq_bin_out=total(chisq_bin_map,/nan)/double(nroi_bin-1)
nbins=double(nroi_bin)

;writefits,'chisq_bin_map.fits',chisq_bin_map
chisqmap_bin=chisq_bin_map

end
