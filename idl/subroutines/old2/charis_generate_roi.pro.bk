pro charis_generate_roi,output=output,snrmap=snrmap,snrlim=snrlim,roirange=roirange,roidim=roidim

if (~keyword_set(snrmap) and ~keyword_set(snrlim) and ~keyword_set(roirange)) then begin
print,'pro charis_generate_roi,output=output,snrmap=snrmap,snrlim=snrlim,roidim=roidim,roirange=roirange'
print,'***required input***'
print,'1a. snrmap=snrmap    ; SNR map of some square dimension (probably 201x201)...or '
print,'1b. roirange=roirange     ;dx,dy,pa of ROI '
print,'2. output=output ; the output array: otherwise nothing gets saved'
print,'...'
print,'**optional inputs**'
print,'snrlim - limiting SNR value in map to define ROI (default is snrlim=2)'
print,'...'
print,'For now, defaults to CHARIS dimensions (i.e. roidim=[201,201])
print,'Code will crash if the SNRMAP does not have dimensions of [201,201]***'
goto,skiptoend
endif

;do this to define a "region of interest", defined 

;Inputs ...
;if snrmap is set ...
;...use a signal to noise map and define the limits based on some SNR cutoff
;if roirange is set ...
;... define a geometric area

;Outputs ...
;output - an array of [201,201] that has 1 for the roi and 0 for outside it

;***dimensions****
;for simplicity, assume the dimensions are 201x201 by default

if ~keyword_set(roidim) then begin 
dimx=201 & dimy=201
endif else begin
dimx=roidim[0] & dimy=roidim[1]
endelse

;***roi snrmap calc is started
if keyword_set(snrmap) then begin

;read in the SNRmap array

;dimensions are implicit

if ~keyword_set(snrlim) then snrlim=2.

outputsnrmap=snrmap*0
roisnrmap=where(snrmap gt snrlim,nroi)
if nroi gt 0 then outputsnrmap[roisnrmap]=1

endif else begin

outputsnrmap=fltarr(dimx,dimy)
outputsnrmap[*]=1
endelse
;
;***snrmap roi calc is done

;***roi geometry is started

if ~keyword_set(roirange) then begin
;assume the entire image is your evaluation area ...

outputroirange=fltarr(dimx,dimy)
outputroirange[*]=1
endif else begin
dx=roirange[0]
dy=roirange[1]
rotang=roirange[2]
;endelse

outputroirange=fltarr(dimx,dimy)
outputcen=[dimx/2,dimy/2]
print,-dx/2+outputcen[0],dx/2+outputcen[0]
outputroirange[-dx/2+outputcen[0]:dx/2+outputcen[0],-dy/2+outputcen[1]:dy/2+outputcen[1]]=1

if rotang ne 0 then begin
outputroirange=rotat(outputroirange,rotang)
good=where(outputroirange gt 0 and finite(outputroirange) ne 0,ngood,complement = bad)
if ngood gt 0 then outputroirange[good]=1 & outputroirange[bad]=0
endif

endelse

;now you have roi from snrmap and from geometry set.   

;if keyword_set(snrmap) then begin
output=fltarr(dimx,dimy)
goodsnrmap=where(outputsnrmap eq 1,nsnr)
goodroimap=where(outputroirange eq 1,nmap)
goodoutput=intersect(goodsnrmap,goodroimap)
;print,goodoutput,goodsnrmap
;print,nsnr,nmap
;stop
output[goodoutput]=1
output=long(output)
;writefits,'outputsnrmap.fits',outputsnrmap
;writefits,'outputroi.fits',outputroirange
skiptoend:

end
