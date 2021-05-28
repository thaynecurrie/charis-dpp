function angarr,dim1,dim2
;returns an array where each element has a value
; its angle to the center of the image
; of pi-pi relative to the axis of x positive

if (n_params() eq 1) then dim2=dim1

x=findgen(dim1)#replicate(1.,dim2)-dim1/2
y=replicate(1.,dim1)#findgen(dim2)-dim2/2
ang=atan(y,x)
i=where((x eq 0.) and (y eq 0.))
if (i[0] ne -1) then ang[i]=0.

return,ang
end
