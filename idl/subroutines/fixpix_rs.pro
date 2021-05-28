;;;
;;; function names modified for REDSPEC by SSKIM (Jan 2001)
;;;


FUNCTION noise_rs,im

COMPILE_OPT idl2, hidden

; function to find the std dev of a block of pixels

; sort the array

imsort = im[sort(im)]

; find the number of pixels

npix = (size(imsort))[1]

; pick the integer that is closest to 2/3 of the number of pixels
; this is an integer

num = 2*npix/3 + 1

; find the noise as a the minimum difference between
; the averages of two pairs of points approximately 2/3 of the
; entire array length apart

noise = min(imsort[num-1:npix-1]+imsort[num-2:npix-2]-imsort[0:npix-num]-imsort[1:npix-num+1])/2.0

; now correct the noise for the fact that it's not exactly over 2/3 of the
; elements

noise = noise / (1.9348+(float(num-2)/float(npix) - 0.66667)/.249891)

; now correct for the finite number of pixels

if npix gt 100 then $
        noise = noise/(1.0-0.6141*exp(-0.5471*alog(npix))) $
        else $
        noise = noise/(1.0-0.2223*exp(-0.3294*alog(npix)))

return, noise

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION fixpix_rs, im,maskarr,iter=iternum,quiet=quiet

; sigma size problem fixed by SSKIM (Sep 2000)

COMPILE_OPT idl2, hidden

; INITIALIZATIONS

; in case of an error return to the calling prog

;on_error, 2

; determine the median value of the entire image

medimg = median(im, /even)

; create a new image to be the cleaned up output

imnew = im

; set the maskarray to zero

maskarr = im*0.0

; if the number of iterations is not explicitly declared set it to 3

if (not keyword_set(iternum)) then iternum =3

; set mparm to 1 as a default


mparm = 1

qf = (not keyword_set(quiet))


;; READ THE HEADER FILE
;
;if (not sxpar_rs(h,"D80P")) then begin
;    if qf then print, "Warning: This image has not been run through D80P"
;    mparm = sxpar_rs(h,"MPARM201")
;    if (mparm eq 0) then begin
;	 if qf then print, "Warning: No MPARM found in the header file. Setting to default."
;	 mparm = 1
;    endif
;    if qf then print, format='("mparm = ",i0)',mparm
;endif	
;

; DETERMINE THE AVERAGE NOISE IN THE ENTIRE IMAGE

; scan through 5x5 boxes to determine some average estimate of the
; image noise

xs = (size(im))[1]
ys = (size(im))[2]
;sigma = fltarr((xs-1)/5*(ys-1)/5)
sigma = fltarr(xs/5*ys/5)  ; Corrected by SSKIM

n=0L
for i=2,xs-3,5 do begin

    for j=2,ys-3,5 do begin

        tmp=im[i-2:i+2,j-2:j+2]

        srt=tmp[sort(tmp)]

        sigma[n] = (min(srt[16:24]+srt[15:23]-srt[0:8]-srt[1:9])/2.0)/1.540

        n = n+1
    endfor

endfor

; define the median value of the sigmas, and the sigma of the sigmas

medsig = median(sigma,/even)
;sigsig = noise(sigma)
sigsig = noise_rs(sigma)  ; Modified by SSKIM


; BEGIN SCANNING FOR HOT & COLD PIXELS

; start at (4,4) (so that any pixel in the box can have a 5x5 
; square centered on it.  find the hottest and coldest pixels

; loop through several iterations of replactments

for iter=1,iternum do begin

    hotcnt = 0                  ; variables to count the replacements of hot and cold pixels
    coldcnt = 0
    addon = ((iter+1) mod 2) *2

    for i=4 + addon,xs-5,5 do begin

        for j=4 + addon,ys-5, 5 do begin

            box = imnew[i-2:i+2,j-2:j+2]
            
            hotval = max(box,hoti)
            hotx = hoti mod 5 + i-2 ; coords in original image
            hoty = hoti/5 + j-2

            coldval = min(box,coldi)
            coldx = coldi mod 5 + i-2
            coldy = coldi/5 + j-2

                                ; begin the decision process for the hottest pixel
            
            hot = imnew[hotx-2:hotx+2,hoty-2:hoty+2]
            med8 = median([hot[1:3,1],hot[1:3,3],hot[1,2],hot[3,2]],/even)
            med16 = median([hot[0:4,0],hot[0:4,4],transpose(hot[0,1:3]),transpose(hot[4,1:3])],/even)
            srt = hot[sort(hot)]
            sig = (min(srt[16:24]+srt[15:23]-srt[0:8]-srt[1:9])/2.0)/1.540

                                ; decide from the noise in the box if we 
                                ; are on a feature or a gaussian background

            if sig gt (medsig + 2.0*sigsig) then $
              sig = max([medsig+2.0*sigsig,sqrt(5.0+.210*abs(med8)/mparm)*mparm])

                                ; decide whether to replace pixel

            if ((med8+2.0*med16)/3.0 -medimg) gt 2.0*sig then begin
                if (imnew[hotx,hoty] gt (2.0*med8-medimg+3.0*sig)) then begin
                    imnew[hotx,hoty] = med8
                    maskarr[hotx,hoty]=1
                    hotcnt=hotcnt+1
                endif
            endif else begin
                if ((imnew[hotx,hoty]-(med8+2.0*med16)/3) gt 5.0*sig) then begin
                    imnew[hotx,hoty] = med8
                    maskarr[hotx,hoty]=1
                    hotcnt=hotcnt+1
                endif
            endelse

                                ; begin the decision process for the coldest pixel
            
            cld = imnew[coldx-2:coldx+2,coldy-2:coldy+2]
            med8 = median([cld[1:3,1],cld[1:3,3],cld[1,2],cld[3,2]],/even)
            med16 = median([cld[0:4,0],cld[0:4,4],transpose(cld[0,1:3]),transpose(cld[4,1:3])],/even)
            srt = cld[sort(cld)]
            sig = (min(srt[16:24]+srt[15:23]-srt[0:8]-srt[1:9])/2.0)/1.540

                                ; decide from the noise in the box if we 
                                ; are on a feature or a gaussian background

            if sig gt (medsig + 2.0*sigsig) then $
              sig = max([medsig+2.0*sigsig,sqrt(5.0+.210*abs(med8))])

                                ; decide whether to replace pixel

            if ((med8+2.0*med16)/3.0 -medimg) lt -2.0*sig then begin
                if (imnew[coldx,coldy] lt (2.0*med8-medimg-3.0*sig)) then begin
                    imnew[coldx,coldy] = med8
                    maskarr[coldx,coldy]=1
                    coldcnt=coldcnt+1
                endif
            endif else begin
                if ((imnew[coldx,coldy]-(med8+2.0*med16)/3) lt -5.0*sig) then begin
                    imnew[coldx,coldy] = med8
                    maskarr[coldx,coldy]=-1
                    coldcnt = coldcnt+1
                endif
            endelse
            

        endfor		
        
    endfor

    if qf then print, format='(i5,i5,"  hot and cold pixels in iter. ",i0)',hotcnt,coldcnt,iter
endfor

;; add keyword to header file
;sxaddpar,h,'FIX_ITER',iternum

return, imnew

end

