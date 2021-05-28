pro grater_grid
;generates a grid of grater files

;set these
itilt_in=[56,60]
;r0_in=46.5
beta_in=0

r0_range=[25,36.9,39.1,41.3]
gscat_range=[0,0.1,0.2]
ksi_range=[3,5]
alphain_range=[1,2.5,5,10]
alphaout_range=[-2.5]
;alphaout_range=[-5,-10]
pa_range=[-1.2]
;ydo_range=[4.7,6.1,7.4]
;xdo_range=[2.4,2.7,3]
ydo_range=[4.3,6.1,7.9]
xdo_range=[0.5,1.6,2.7]
;e_range=[0,0.05,0.1]
;e_range=[0.05,0.1,0.15,0.2,0.3]
;e_range=[0,0.1,0.2,0.3]
;e_in=0.
e_range=[0.]
;theta_range=[0,45,90,135,180]
theta_range=[0.]
openw,18,'modellist'
openw,19,'model_params'
printf,18,'MODEL_LIST'
printf,19,'MODEL_NAME','G','KSI0','ALP_I','ALP_O','BETA','XDO','YDO','INC','THETA0',$
   format='(a25,a6,1x,8(a7,1x))'
;printf,19,'MODEL_NAME','G','KSI0','ALP_I','ALP_O','BETA','XDO','YDO','E','THETA0',$
   ;format='(a25,a6,1x,8(a7,1x))'

;print,'gscat',n_elements(gscat_range)

print,n_elements(gscat_range)*n_elements(ksi_range)*n_elements(alphain_range)*n_elements(alphaout_range)*n_elements(pa_range)*n_elements(r0_range)*n_elements(xdo_range)*n_elements(ydo_range)*n_elements(itilt_in)
;print,n_elements(gscat_range)*n_elements(ksi_range)*n_elements(alphain_range)*n_elements(alphaout_range)*n_elements(pa_range)*n_elements(r0_range)*n_elements(xdo_range)*n_elements(ydo_range)
;stop

for q=0L,n_elements(r0_range) - 1 do begin
for i=0L,n_elements(gscat_range)-1 do begin
 for ii=0L,n_elements(ksi_range) -1 do begin
  for iii=0L,n_elements(alphain_range)-1 do begin
   for iv=0L,n_elements(alphaout_range)-1 do begin
    for j=0L,n_elements(pa_range)-1 do begin
     for ji=0L,n_elements(xdo_range)-1 do begin  
      for jii=0L,n_elements(ydo_range)-1 do begin
       for jiii=0L,n_elements(e_range)-1 do begin
        for jiv=0L,n_elements(itilt_in)-1 do begin
        ;for jiv=0L,n_elements(theta_range)-1 do begin
;outfilename='grater'+strtrim(i,1)+'.fits'
;outfilename='grater'+string(i,ii,iii,format='i3')+'.fits'
;outfilename='grater'+strtrim(i,1)+strtrim(ii,1)+strtrim(iii,1)+strtrim(iv,1)+$
outfilename='egrater'+strtrim(i,1)+strtrim(ii,1)+strtrim(iii,1)+strtrim(iv,1)+$
   strtrim(j,1)+strtrim(ji,1)+strtrim(jii,1)+strtrim(jiii,1)+strtrim(jiv,1)+strtrim(q,1)+'_model.fits'
print,outfilename
;stop
;;stop
;goto,skipoverme
charis_call_grater,g=gscat_range[i],ksi0=ksi_range[ii],alphain=alphain_range[iii],$
  alphaout=alphaout_range[iv],pa=pa_range[j],xdo=xdo_range[ji],ydo=ydo_range[jii],$
   e=e_in,theta0=theta_range[0],outfile=outfilename,$
   itilt=itilt_in[jiv],r0=r0_range[q],$
dstar=48.2,pfov=0.0162,nx=201,ny=201
;and now the hardwired stuff
  
skipoverme:
printf,18,outfilename,format='(a20)'
printf,19,outfilename,gscat_range[i],ksi_range[ii],alphain_range[iii],alphaout_range[iv],pa_range[j],$
              xdo_range[ji],ydo_range[jii],itilt_in[jiv],theta_range[0],format='(a25,9(f7.3,1x))'
;printf,19,outfilename,gscat_range[i],ksi_range[ii],alphain_range[iii],alphaout_range[iv],pa_range[j],$
              ;xdo_range[ji],ydo_range[jii],e_range[jiii],theta_range[jiv],format='(a25,9(f7.3,1x))'
         endfor
        endfor
       endfor
      endfor
     endfor
    endfor
   endfor
  endfor
 endfor
 endfor

end
