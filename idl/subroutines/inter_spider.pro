function inter_spider,im,spang

;trouve le point d'intersection des 2 spiders dans im
;find the intersection of the 2 spiders in image

sz=size(im) & dimx=sz[1] & dimy=sz[2]

;premier spider
;first spider
t1=spang[0]
a1=1./tan(t1*!dtor)
if (a1 ge 0) then begin
    bmax=dimy-1
    bmin=floor(-a1*(dimx-1))
endif else begin
    bmax=ceil((dimy-1)-a1*(dimx-1))
    bmin=0
endelse
;if(finite(bmin) eq 0) then bmin=0
print,bmax,bmin,sz,dimx,dimy,a1,t1
print,spang
b1=findgen(bmax-bmin+1)+bmin
l1=lineint(im,a1,b1)
m1=max(l1,i1) & b1=b1[i1]

;stop

;deuxieme spider
;second spider
t2=spang[1]
a2=1./tan(t2*!dtor)
if (a2 ge 0) then begin
    bmax=dimy-1
    bmin=floor(-a2*(dimx-1))
endif else begin
    bmax=ceil((dimy-1)-a2*(dimx-1))
    bmin=0
endelse
b2=findgen(bmax-bmin+1)+bmin
l2=lineint(im,a2,b2)
m2=max(l2,i2) & b2=b2[i2]

;stop

;trouve l'intersection des deux droites
;is the intersection of two straight lines
xi=(b2-b1)/(a1-a2)
yi=a1*xi+b1

return,[xi,yi]
end
