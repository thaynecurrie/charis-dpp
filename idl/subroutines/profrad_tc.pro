Pro profrad_tc,imt,res,r1,r2,med4q=med4q,p2d=p2d,p1d=p1d,rayon=rayon,max=max,nan=nan

;Calculates a median radial profile.
; By default, calculates a radial profile 1d consisting of the median
; values of all pixels within a certain distance from
; center of the image. This interval distance is determined by
; res. Then uses a cubic spline to extrapolate the profile of
;, length of each pixel to construct the 2D profile.
;
; res: resolution (in pixels) of the radial profile 1d, optional value
; Default is 1.
; p2d: 2d profile
; P1D: 1d Profile
; radius: radius corresponding P1D
; / med4q: 2d calculates the profile by taking the median of the 4 quadrants

s=size(imt) & dimx=s[1] & dimy=s[2]
im=double(imt)

;calcul de profil 1d
;calculates 1d profile

if (arg_present(p1d) or (not keyword_set(med4q) and arg_present(p2d))) then begin
    if (n_params() eq 1) then res=1. 
    distarr=shift(dist(dimx,dimy),dimx/2,dimy/2)
    if n_elements(r1) eq 0 then rmin=0. else rmin=r1
    if n_elements(r2) eq 0 then rmax=max(distarr) else rmax=r2
    if res ne 0 then rdistarr=round(distarr/res)*res else rdistarr=distarr
    rayon=rdistarr[0:dimx/2,0:dimy/2]
    rayon=rayon[where(rayon ge rmin and rayon le rmax)]
    sind=sort(rayon) & rayon=rayon[uniq(rayon,sind)]
    p1d=dblarr(n_elements(rayon))

    i=where(distarr ge rmin and distarr le rmax)
    sind=i[sort(rdistarr[i])]
    imax=n_elements(sind)-1
    i1=0l
    for r=0,n_elements(rayon)-1 do begin
        for i2=i1,imax do if rdistarr[sind[i2]] ne rayon[r] then break
        p1d[r]=median(im[sind[i1:i2-1]],/even)
        if keyword_set(max) then p1d[r]=max((im[sind[i1:i2-1]]))
        i1=i2
    endfor
   if keyword_set(nan) then begin
     p1d=interpol(p1d,rayon,rayon,/spline,/nan)
   endif
;    for r=0,n_elements(rayon)-1 do $
;      p1d[r]=median(im[where(rdistarr eq rayon[r])],/even)
endif

;calcul du profil 2d
;calculate 2d profile

if (not keyword_set(med4q) and arg_present(p2d)) then begin
    ;interpole le profil sur l'image au complet d'un coup
    ;interpolates the profile on the entire image at once

    ;p2d=dblarr(dimx,dimy)
    ;sind=sort(distarr)
    ;p2d[sind]=spline(rayon,p1d,distarr[sind])

    ;fait un seul quadrant et copie sur les 3 autres quadrants
    ;cette methode est 3X plus rapide pour l'interpolation

    ;is one quadrant and copy the other 3 quadrants
    ;and this method is 3x faster interpolation

    distarr=distarr[0:dimx/2,0:dimy/2]
    rdistarr=rdistarr[0:dimx/2,0:dimy/2]
    q1=dblarr(dimx/2+1,dimy/2+1)       
    p1dt=p1d
    rayont=rayon
    ;traite les NAN
    i=where(finite(p1d) eq 0,c)
    if c gt 0 then begin
        for n=0l,c-1 do begin
            ii=where(rdistarr eq rayon[i[n]],cc)
            if cc gt 0 then q1[ii]=!values.f_nan
        endfor
        remove,i,rayont,p1dt
    endif
    sind=sort(distarr)
    for nmin=0l,n_elements(q1)-1 do if finite(q1[sind[nmin]]) eq 1 and distarr[sind[nmin]] ge rmin then break
    for nmax=nmin+1,n_elements(q1)-1 do if finite(q1[sind[nmax]]) eq 0 or distarr[sind[nmax]] ge rmax then break
    if ~keyword_set(nan) then begin
    q1[sind[nmin:nmax-1]]=spline(rayont,p1dt,distarr[sind[nmin:nmax-1]])
    endif else begin
    q1[sind[nmin:nmax-1]]=interpol(rayont,p1dt,distarr[sind[nmin:nmax-1]],/spline,/nan)
    endelse

    if (dimx mod 2 eq 0) then i1x=1 else i1x=0
    if (dimy mod 2 eq 0) then i1y=1 else i1y=0

    ;pas besoin de prendre de transposes ici car les profils sont
    ;symmetriques aux transposes
    ;not need to take to implement here because the profiles are 
    ;symmetrical to transposes

    p2d=dblarr(dimx,dimy)
    p2d[0:dimx/2,0:dimy/2]=q1
    p2d[0:dimx/2,dimy/2:dimy-1]=reverse(q1[*,i1y:dimy/2],2)
    p2d[dimx/2:dimx-1,dimy/2:dimy-1]=$
      reverse(reverse(q1[i1x:dimx/2,i1y:dimy/2],2),1)
    p2d[dimx/2:dimx-1,0:dimy/2]=reverse(q1[i1x:dimx/2,*],1)
endif

if (keyword_set(med4q) and arg_present(p2d)) then begin
    ;place les 4 quadrants dans un cube, fait la mediane, et replace cette
    ;mediane dans les quadrants
    ;(equivalent a faire un cube avec 4 rotation et prendre la mediane)
    
    ;up to 4 quadrants in a cube, is the median, and replace this
    ; medial quadrants
    ;(equivalent to a cube with 4 rotation and take the median)
    q=dblarr(dimx/2+1,dimy/2+1,4)
    if (dimx mod 2 eq 0) then i1x=1 else i1x=0
    if (dimy mod 2 eq 0) then i1y=1 else i1y=0
    q[*,*,0]=im[0:dimx/2,0:dimy/2]
    q[*,i1y:dimy/2,1]=reverse(im[0:dimx/2,dimy/2:dimy-1],2)
    q[*,*,1]=transpose(q[*,*,1])
    q[i1x:dimx/2,i1y:dimy/2,2]=reverse(reverse(im[dimx/2:dimx-1,dimy/2:dimy-1],2),1)
    q[i1x:dimx/2,*,3]=reverse(im[dimx/2:dimx-1,0:dimy/2],1)
    q[*,*,3]=transpose(q[*,*,3])
    q=median(q,dimension=3,/even)
    ;q=dmedian(q,dimension=3,/even)
;    q=min(q,dimension=3)
    qt=transpose(q)
    p2d=dblarr(dimx,dimy)
    p2d[0:dimx/2,0:dimy/2]=q
    p2d[0:dimx/2,dimy/2:dimy-1]=reverse(qt[*,i1y:dimy/2],2)
    p2d[dimx/2:dimx-1,dimy/2:dimy-1]=$
      reverse(reverse(q[i1x:dimx/2,i1y:dimy/2],2),1)
    p2d[dimx/2:dimx-1,0:dimy/2]=reverse(qt[i1x:dimx/2,*],1)

endif

if arg_present(p2d) then p2d=float(p2d)
end
