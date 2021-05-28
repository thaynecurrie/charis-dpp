function nbrlist,tlist,rns=rns

;Retourne un array de type integer contenant tous les nombres
;specifies par une liste du format '1-5,10,11,18-20'
;
;tlist: string, specification de la liste de nombres
;               "A-B": de A jusqu'a B inclusivement  
;               ",": pour separer les series de nombre consecutifs
;               ex. '1-10,12,15-18'
;
;keyword optionel:
; rns=[r,s]: "Read'N Skip", pour repeter un patron dans une sequence
;            a partir du debut de la liste, inclut les r premiers
;            nombres et saute les s suivants, et recommence
;            jusqu'a la fin de la liste

;Returns an integer array containing all numbers
; specifies a list of -5,10,11,18-20 format '1 '
;
; tlist: string, specification of the list of numbers
; "AB" from A to B inclusive
; "," To separate the series of consecutive numbers
; Ex. '1 -10,12,15-18 '
;
; optional keyword:
; Rns = [r, s]: "Read'N Skip" to repeat a pattern in a sequence
; Has from the beginning of the list includes the first r
; Numbers and skips the following s, and again
; Until the end of the list

txtlist=tlist
nchar=strlen(txtlist)

if (keyword_set(rns)) then period=rns[0]+rns[1]

while (nchar gt 0) do begin
    comma=strpos(txtlist,',')
    if (comma eq -1) then comma=nchar+1
    dash=strpos(txtlist,'-')
    if (dash gt comma) then dash=-1
    ;if there is a dash in this part of the string
    if (dash ne -1) then begin
        ni=fix(strmid(txtlist,0,dash))
        nf=fix(strmid(txtlist,dash+1,comma-dash-1))

        if (keyword_set(rns)) then begin
            for n=ni,nf do begin
              if ((n-ni) mod period lt rns[0]) then $
                  if (n_elements(list) eq 0) then list=n else list=[list,n]
            endfor
        endif else begin
            for n=ni,nf do $
                if (n_elements(list) eq 0) then list=n else list=[list,n]
        endelse
    ;if no dash in this part of the string
    endif else begin
        n=fix(strmid(txtlist,0,comma))
        if (n_elements(list) eq 0) then list=n else list=[list,n]
    endelse

    ;number of char not read yet
    nchar=nchar-comma-1
    nchar=nchar > 0
    ;trim txtlist to remove characters read
    if (nchar gt 0) then txtlist=strmid(txtlist,comma+1,nchar)
endwhile

return,list
end
