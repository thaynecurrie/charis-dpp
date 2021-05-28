pro chisqcalcwbin,obs=obs,eobs=eobs,widthobs=widthobs,model_in=model_in,err_model_in=err_model_in,alpha_in=alpha_in,$
chisqout=chisqout,ndof=ndof,alpha_out=alpha_out,phot=phot,ephot=ephot,widthphot=widthphot,modelphot=modelphot

;Feb 11, 2021 - chi-squared without covariance calculations but with bin widths

;July 7, 2018 - simple chi-squared without covariance calculations

;number of flux measurements
n_measure=n_elements(obs) 
ndof=n_measure

ncovar=fltarr(n_measure,n_measure)

chisqtot=1d30

;if the bin width is not set, then assume it is 1
if ~keyword_set(widthobs) then begin
widthobs=obs*0+1.
endif

sfac=(findgen(1d5)*1d-3+1d-3)*alpha_in
nsfac=n_elements(sfac)
sfac_tot=-99

for i=0L,nsfac-1 do begin
diff=obs-sfac[i]*model_in
if ~keyword_set(err_model_in) then begin
gval=total(widthobs(diff/eobs)^2.)

if keyword_set(phot) then begin 
diffphot=phot-sfac[i]*modelphot
gval+=total((diffphot/ephot)^2.*widthphot)
endif else begin

gval/=total(widthphot)
endelse

endif else begin

;****these two are equivalent
;1. 
;cq=diag_matrix(eobs^2.+sfac[i]^2.*err_model_in^2.)
;gval=(transpose(diff)#invert(cq,/double)#diff)*widthobs/widthobs
;2.

gval=total(widthobs*diff^2./(eobs^2.+sfac[i]^2.*err_model_in^2.))

endelse
if gval lt chisqtot then begin
chisqtot=gval
sfac_tot=sfac[i]
endif
endfor

alpha_out=sfac_tot
chisqout=chisqtot

end
