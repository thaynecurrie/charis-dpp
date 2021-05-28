pro grater_demo,outfilename=outfilename

beta_in=1
;r0_range=96
r0_range=50
;gscat_range=0.25
gscat_range=0.5
ksi_range=3.
alphain_range=3
alphaout_range=-2
pa_range=98.8+180
;theta_range=112.5-90-60
;theta_range=112.5-90-120
theta_range=112.4
;itilt_in=86.2
itilt_in=80.
;itilt_in=0
;xdo=-10

e_in=0.0

charis_call_grater,g=gscat_range,ksi0=ksi_range,alphain=alphain_range,beta=beta_in,$
  alphaout=alphaout_range,pa=pa_range,xdo=xdo,$
   e=e_in,theta0=theta_range,outfile=outfilename,$
   itilt=itilt_in,r0=r0_range,$
;dstar=48.2
dstar=48.2,pfov=0.0164,nx=601,ny=601

end

