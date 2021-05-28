function mkspider,dim,st,angles,justind=justind,cen=cen

;dim: dimension de l'image a retourner
;st: epaisseur d'un bras du spider en pixels
;angles: angles des bras du spider, sens horaire a partir des y positifs
;/justind: pour retourner seulement les indices des pixels inclus
;          dans le spider

; dim: dimension of the image return
; st: thickness of an arm of the spider in pixels
; angles: angles of the spider arms, clockwise from positive y
; / justind: to return only indices of the pixels included
;And in the spider


;array du spider
;array of spider

spider=dblarr(dim,dim)
;centre du spider
;center of spider

if keyword_set(cen) then begin
    xc=cen[0] & yc=cen[1]
endif else begin
    xc=dim/2 & yc=dim/2
endelse

for nang=0,n_elements(angles)-1 do begin
    t=angles[nang]
    t=(t+360.) mod 360.
    if (t mod 90. eq 0.) then begin
        case t of
            0.:   begin & x1=xc-ceil(st/2.) & x2=xc+ceil(st/2.) & y1=yc & y2=dim-1 & end
            90.:  begin & x1=xc & x2=dim-1  & y1=yc-ceil(st/2.) & y2=yc+ceil(st/2.) & end
            180.: begin & x1=xc-ceil(st/2.) & x2=xc+ceil(st/2.) & y1=0 & y2=yc & end
            270.: begin & x1=0 & x2=xc      & y1=yc-ceil(st/2.) & y2=yc+ceil(st/2.) & end
        endcase

        if keyword_set(justind) then begin
            dx=x2-x1+1 & dy=y2-y1+1
            ind=(replicate(1,dx)#lindgen(dy)+y1)*dim+(lindgen(dx)#replicate(1,dy)+x1)
            ind=reform(ind,n_elements(ind))
            if (n_elements(indall) eq 0) then indall=ind else indall=[indall,ind]
        endif else begin
            if (t eq 0. or t eq 180.) then begin
                for x=x1,x2 do begin
                    frac=((x+0.5) < (xc+st/2.d)) - ((x-0.5) > (xc-st/2.d))
                    spider[x,y1:y2]=frac>spider[x,y1:y2]
                endfor
            endif else begin
                for y=y1,y2 do begin
                    frac=((y+0.5) < (yc+st/2.d)) - ((y-0.5) > (yc-st/2.d))
                    spider[x1:x2,y]+=frac
                    spider[x1:x2,y]+=frac
                endfor
            endelse
        endelse            

    endif else begin
        ;equation de la droite du spider
        ;equation of the right of the spider

        a=1/tan(t*!pi/180.) & b=yc
        ;indices des pixels qui interceptent le spider
        ;indices of pixels that intercept the spider

        xarr=lindgen(dim,dim) mod dim & yarr=lindgen(dim,dim)/dim

        ;on ramene l'angle entre 0 et 90 deg
        ; we reduce the angle between 0 and 90 deg

        if (t gt 0 and t lt 90) then begin
            t=t*!pi/180.
            ;x1=xc & x2=dim-1
            x1=xc-st/2.*cos(t) & x2=dim-1
            ;ymin=yc & ymax=dim-1
            sign=1.
            quad=1
        endif
        if (t ge 90 and t lt 180) then begin
            t=(180.-t)*!pi/180.
            x1=xc-st/2.*cos(t) & x2=dim-1
            sign=-1.
            quad=2
        endif
        if (t gt 180 and t le 270) then begin
            t=(t-180.)*!pi/180.
            x1=0 & x2=xc+st/2.*cos(t)
            sign=1.
            quad=3
        endif
        if (t gt 270 and t lt 360) then begin
            t=(360.-t)*!pi/180.
            x1=0 & x2=xc+st/2.*cos(t)
            sign=-1.
            quad=4
       endif
        ;interval en y des pixels qui interceptent les droites limites
        ;interval in pixels that will intercept the straight boundaries

        dy=st/sin(t)+1./tan(t)+1.
        ;interval entre les droites limites
        ;interval between straight boundaries
        db=st/sin(t)
        ind=-1l
        x1=floor(x1) & x2=ceil(x2)
        for x=x1,x2 do begin
            ;y central du spider
            y=yc+sign*(x-xc)/tan(t)
            ;ymin et ymax des pixels du spider
            y1=floor(y-dy/2.)>0
            y2=ceil(y+dy/2.)<(dim-1)
            ;indices de ces pixels
            if (y2-y1 ge 0) then $
              newi=(lindgen(y2-y1+1)+y1)*dim+replicate(x,y2-y1+1) else newi=-1l
            ;newi=where(xarr eq x and yarr ge ymin and yarr le ymax)
            if (newi[0] ne -1) then $
              if (ind[0] eq -1) then ind=newi else ind=[ind,newi]

        endfor

        if keyword_set(justind) then begin
            if (n_elements(indall) eq 0) then indall=ind else indall=[indall,ind]
        endif else begin
            ;integrale sur la surface des pixels
            b1=b-db/2. & b2=b+db/2.
            for i=long(0),n_elements(ind)-1 do begin
                xi=(ind[i] mod dim)-0.5 & yi=(ind[i]/dim)-0.5
                npas=100 ;nombre de sous-division de pixel pour integration
                x=xi+findgen(npas)/float(npas)
                ymin=yi & ymax=yi+1.
                if quad eq 1 or quad eq 4 then ymin=((xc-x)/a+b)>ymin
                if quad eq 2 or quad eq 3 then ymax=((xc-x)/a+b)<ymax
                yf=ymax<(a*(x-xc)+b2)
                yi=(ymin>(a*(x-xc)+b1))<yf
                aire=total(yf-yi)/float(npas)
                spider[ind[i]]+=aire
            endfor
        endelse ;keyword_set(justind)
    endelse ;(t eq 0 or t eq 180)
endfor
spider=spider<1.

if (keyword_set(justind)) then begin
    indall=indall[sort(indall)]
    indall=indall[uniq(indall)]
    return,indall
endif

return,1.d0-spider
end
