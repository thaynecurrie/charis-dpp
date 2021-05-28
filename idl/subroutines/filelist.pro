function filelist,flist,nfiles,prefix=prefix,suffix=suffix,dir=dir,gz=gz,$
                  nzero=nzero,rns=rns,ext=ext

if n_elements(flist) eq 1 or size(flist,/type) ne 7 then begin
    if size(flist,/type) eq 7 then nlist=nbrlist(flist,rns=rns) else nlist=flist
    ;build the list of file names
    ;construit la liste des noms des fichiers
    if (not keyword_set(suffix)) then suffix=''
    if (not keyword_set(nzero)) then nzero=4
    if ~keyword_set(ext) then ext='.fits'
    if keyword_set(gz) then ext+='.gz'
    fnames=prefix+nbr2txt(nlist,nzero)+suffix+ext
endif else fnames=flist
if ~keyword_set(dir) then dir=''
fnames=dir+fnames
nfiles=n_elements(fnames)

return,fnames
end
