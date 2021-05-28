pro param,parname,parval,parcom,set=set,get=get,del=del,count=count,$
          pfname=pfname,before=before,after=after

; parname: string, name of the parameter
; ParVal: value of parameter
; PARCOM: string, comment the parameter, optional
; / set: to establish a parameter and its value
; / get: to return the value of the parameter in ParVal
; count: To return the number of keyword match found
; Has parname
; pfname: string, if another file name that is desired param.dat

if ~keyword_set(pfname) then pfname='param.dat'

if file_test(pfname) then begin
    ;fichier_exist, verifie si on peut le lire
    r=file_test(pfname,/read)
    if (r eq 0) then begin
        print,'Impossible de lire '+pfname
        return
    endif
    ;on peut le lire, ouvre et lit array des parametres
    openr,funit,pfname,/get_lun
    line=''
    while (not eof(funit)) do begin
        readf,funit,line
        if (n_elements(p) eq 0) then p=[line] else p=[p,line]
    endwhile
    free_lun,funit
endif else begin
    ;cree array des parametres
    p=['END     '] ;les espaces sont importants
endelse

if keyword_set(get) then parval=sxpar_param(p,parname,comment=parcom,count=count)

if n_params() lt 3 then parcom=''

if keyword_set(set) then sxaddpar_dl,p,parname,parval,parcom,before=before,after=after

if keyword_set(del) then sxdelpar,p,parname

if (file_test(pfname) eq 1 and file_test(pfname,/write) eq 0) then begin
    print,'Impossible d''ecrire dans '+pfname
    return
endif
openw,funit,pfname,/get_lun
printf,funit,p
free_lun,funit

end
