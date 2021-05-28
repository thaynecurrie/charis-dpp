pro charis_est_astrom_grid,pfname,datacube=datacube,xrange=xrange,yrange=yrange,deltapos=deltapos,fitrange=fitrange,bklim=bklim,bklir=bklir,noexten=noexten,help=help,nothroughputcor=nothroughputcor

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
print,"a simple wrapper to estimate position for a source of unknown intensity distribution [not a point source]"
print,""
print,"charis_est_astrom_grid,pfname,datacube=datacube,xrange=xrange,yrange=yrange,deltapos=deltapos,fitrange=fitrange,bklim=bklim,bklir=bklir,noexten=noexten"
print,""
print,"Example: charis_est_astrom_grid,'LkCa15_low.info',pfname,datacube='pca.fits',xrange=[99,102],yrange=[61,64]"
print,""
print,"***Keywords***"
print,""
print,"*datacube - the datacube"
print,"*xrange,yrange - x and y ranges of initial positions"
print,"*deltapos - delta x,y step size"
print,"*fitrange - channels over which to estimate the centroid"
print,"*bklim,bklir - azimuthal and radial regions to compute background"
print,"*noexten - assume a point source distribution"

goto,skiptotheend
endif


;a simple wrapper to estimate position for a source of unknown intensity distribution

if ~keyword_set(xrange) then xrange=[99,102]
if ~keyword_set(yrange) then yrange=[61,64]
if ~keyword_set(deltapos) then deltapos=1.0
if ~keyword_set(fitrange) then fitrange=[0:21]

if ~keyword_set(bklim) then bklim=3
if ~keyword_set(bklir) then bklir=1

nx=long((xrange[1]-xrange[0])/deltapos)+1
ny=long((yrange[1]-yrange[0])/deltapos)+1

xcguess=findgen(nx)*deltapos+xrange[0]
ycguess=findgen(ny)*deltapos+yrange[0]

snrange=fltarr(nx,ny)
xoutrange=fltarr(nx,ny)
youtrange=fltarr(nx,ny)

xin=fltarr(nx,ny)
yin=fltarr(nx,ny)


for i=0L,nx-1 do begin


 for j=0L,ny-1 do begin

 xin[i,j]=xcguess[i]
 yin[i,j]=ycguess[j]
if ~keyword_set(noexten) then begin
 if ~keyword_set(nothroughputcor) then begin
 charis_extract_1d_spectrum,pfname,coords=[xcguess[i],ycguess[j]],datacube=datacube,/nocalerror,/plotfromzero,/fitbckgd,/filtsnr,/fitcollapse,/throughputcor,bklim=bklim,bklir=bklir,/extended,fitrange=fitrange,$
xposout=xposout,yposout=yposout,snrout=snrout,/breakout
 endif else begin
 charis_extract_1d_spectrum,pfname,coords=[xcguess[i],ycguess[j]],datacube=datacube,/nocalerror,/plotfromzero,/fitbckgd,/filtsnr,/fitcollapse,bklim=bklim,bklir=bklir,/extended,fitrange=fitrange,$
xposout=xposout,yposout=yposout,snrout=snrout,/breakout

 endelse
endif else begin
 charis_extract_1d_spectrum,pfname,coords=[xcguess[i],ycguess[j]],datacube=datacube,/nocalerror,/plotfromzero,/fitbckgd,/filtsnr,/fitcollapse,/throughputcor,bklim=bklim,bklir=bklir,fitrange=fitrange,$
xposout=xposout,yposout=yposout,snrout=snrout,/breakout
endelse
 snrange[i,j]=snrout
 xoutrange[i,j]=xposout
 youtrange[i,j]=yposout
 print,'i j ',i,j,' x y and snr are ',xposout,yposout,snrout

 endfor

endfor

;calculate the median position and standard deviation
xoutrange=reform(xoutrange,nx*ny*1L)
youtrange=reform(youtrange,nx*ny*1L)
snrange=reform(snrange,nx*ny*1L)

xin=reform(xin,nx*ny*1L)
yin=reform(yin,nx*ny*1L)

writecol,'astrom_output.txt',xin,yin,xoutrange,youtrange,snrange

good=where(finite(xoutrange) ne 0 and finite(youtrange) ne 0)
xavg=median(xoutrange[good],/even)
yavg=median(youtrange[good],/even)
delxavg=stdev(xoutrange[good])
delyavg=stdev(youtrange[good])

;weighted value

xweight=total(xoutrange[good]*snrange[good])/total(snrange[good],/nan)
yweight=total(youtrange[good]*snrange[good])/total(snrange[good],/nan)

print,'averaged astrometry is ',xavg,' +/- ',delxavg,' ',yavg,' +/- ',delyavg


writecol,'avgastrom.txt',xavg,delxavg,yavg,delyavg

print,'weighted astrometry is ',xweight,' +/- ',delxavg,' ',yweight,' +/- ',delyavg
writecol,'weighted_astrom.txt',xweight,delxavg,yweight,delyavg

skiptotheend:
end
