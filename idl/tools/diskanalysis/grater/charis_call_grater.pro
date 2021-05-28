;pro call_grater,g=g,dstar=dstar,pfov=pfov,nx=nx,ny=ny,r0=r0,alphain=alphain,alphaout=alphaout,ksi0=ksi0,beta=beta,xdo=xdo,ydo=ydo,e=e,itilt=itilt0,pa=pa,theta0=theta0,outfile=outfile
pro charis_call_grater,g=g,dstar=dstar,pfov=pfov,nx=nx,ny=ny,$
 r0=r0,itilt=itilt,pa=pa,alphain=alphain,alphaout=alphaout,$
  ksi0=ksi0,beta=beta,xdo=xdo,ydo=ydo,e=e,theta0=theta0,outfile=outfile,$
   conv=conv


  ;nx = 281                    ; number of pixels along the x-axis
  ;ny = 281                      ; number of pixels along the y-axis
  ;pfov = 0.014                  ; pixel field of view, arcsec

;goto,skipoverme

if ~keyword_set(g) then g = 0.
if ~keyword_set(r0) then r0=0.
if ~keyword_set(itilt) then itilt=0.
if ~keyword_set(pa) then pa = 0.
if ~keyword_set(nx) then nx=201
if ~keyword_set(ny) then ny=201
if ~keyword_set(pfov) then pfov=charis_get_constant(name='pixscale')
if ~keyword_set(dstar) then dstar=110.5
  if ~keyword_set(xdo) then xdo = 0.                      ; disk offset along the x-axis in the disk frame, AU
  if ~keyword_set(ydo) then ydo = 0.                     ; disk offset along the y-axis in the disk frame, AU
  if ~keyword_set(e) then e = 0.0                        ; disk eccentricity

  if ~keyword_set(alphain) then alphain = 10.d0               ; power law index of the radial volume density profile for r<r0
  if ~keyword_set(alphaout) then alphaout = -5.d0              ; power law index of the radial volume density profile for r>r0
  if ~keyword_set(ksi0) then ksi0 = 5                   ; disk vertical scale height at r0, AU
  if ~keyword_set(beta) then beta = 1.d0                   ; flaring index
  if ~keyword_set(theta0) then theta0 = 0.                   ; argument of pericenter, degrees

skipoverme:
print,g,dstar,pfov,nx,ny,r0,itilt,pa
print,alphain,alphaout,ksi0,beta,theta0
print,xdo,ydo,e
print,'theta is',theta0
;stop
;stop

;itilt,r0,pa
;stop
  ;loadct,3
;if N_PARAMS() eq 0 then begin
;
  ;-- Parameters:
  ;denstype = 'gaussian'
  denstype = '2powerlaw'
  
goto,skipoverparams
  if ~keyword_set(g) then g =  0.2                    ; anisotropic scattering factor, -1<g<1

 if ~keyword_set(r0) then r0 = 47.                      ; reference radius, AU (usually close to the peak density position)

skipoverparams:



  gamma = 2.d0                  ; shape of the vertical profile (eg: 1.=exponential, 2.=gaussian)

  ;dens = {type:'gaussian', density0:1.d0, $
  ;        r0:r0, alphain:alphain, alphaout:alphaout, $ ; radial structure
  ;        beta:beta, gamma:gamma, ksi0:ksi0, $         ; vertical structure
  ;        rmin:0.d0, rmax:0.d0}                        ; rmin and/or rmax = 0 means: let the code compute it for you
  dens = {type:'2powerlaw', density0:1.d0, $
          r0:r0, alphain:alphain, alphaout:alphaout, $ ; radial structure
          beta:beta, gamma:gamma, ksi0:ksi0, $         ; vertical structure
          rmin:0.d0, rmax:0.d0}                        ; rmin and/or rmax = 0 means: let the code compute it for you


  ;-- Compute the scattered light image:
  map = scimage( $
        g=g, $                  ; anisotropic scattering factor, -1<g<1
        xdo=xdo,ydo=ydo, $      ; disk offsets in the disk frame, in AU
        e=e,theta0=theta0, $    ; eccentricity and argument of pericenter
        itilt=itilt, $          ; disk inclination, in degrees (0=pole-on)
        PA = pa, $              ; position angle of the disk major axis (East of North), in degrees
        nx=nx,ny=ny, $          ; numbers of pixels along the x- and y-axis
        pfov=pfov, $            ; pixel field of view, in arcsec
        dens = dens, $          ; density structure
        dstar=dstar)            ; stellar distance, in parsec
  

  ;-- Write the results in a fits file:
;goto,skipheader
sxaddpar,h1,'simple','T'
sxaddpar,h1,'naxis',long(2)
sxaddpar,h1,'naxis1',nx
sxaddpar,h1,'naxis2',ny
sxaddpar,h1,'g',g
sxaddpar,h1,'ksi0',ksi0
sxaddpar,h1,'alphain',alphain
sxaddpar,h1,'alphaout',alphaout
sxaddpar,h1,'beta',beta
sxaddpar,h1,'r0',r0
sxaddpar,h1,'xdo',xdo
sxaddpar,h1,'ydo',ydo
sxaddpar,h1,'e',e
sxaddpar,h1,'pa',pa
sxaddpar,h1,'theta0',theta0
sxaddpar,h1,'itilt',itilt
skipheader:

if keyword_set(conv) then begin
psfg=psf_gaussian(npix=nx,fwhm=3)
map=convolve(map,psfg)
endif

if ~keyword_set(outfile) then begin
  ;writefits,'myimage.fits',0,h1
  ;writefits,'myimage.fits',map,/append
  writefits,'myimage.fits',map,h1
endif else begin
;writefits,outfile,0,h1
;writefits,outfile,map,/append
writefits,outfile,map,h1
endelse
  ;tvim,map,/sc,title='Scattered light image'
  ;stop
end
