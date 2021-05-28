function charis_check_health,inputcube

;this function takes the input cube, does some basic statistics and determines if it is:
;Goal: 1. a good science cube, 2. a good sky cube, 3. a bad science/sky cube (usually caused by bad extraction)
;For now: just check whether it is a good cube or a bad cube

;Is the cube good or bad?
avgimage=mean(inputcube,dimension=3)
dist_circle,g,(size(avgimage,/dim))[0]
bad=where(g le 60 and finite(avgimage) eq 0,nbad)

if nbad gt 5 then begin

healthflag='bad'
return,healthflag

endif else begin

healthflag='good'
;return,healthflag

endelse

;Is the cube a science frame or a sky frame?
;TO DO!!!

inner=where(g lt 20 and finite(avgimage) ne 0 and avgimage ne 0,ninner)
outer=where(g gt 60 and finite(avgimage) ne 0 and avgimage ne 0,nouter)

if mean(avgimage[inner]) gt 3*mean(avgimage[outer]) then begin
healthflag='science'
endif else begin
healthflag='sky'
endelse

return, healthflag

;2. take median of region beyond 60 pixels, take median of region within 20 pixels
;3. if median within 20 pixels is > 10x higher then this is a science frame


end
