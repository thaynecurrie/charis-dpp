pro charis_set_files,pfname,parvariable,keyword,obslog,inputarray=inputarray

;pfname - parameter file
;keyword - the keyword
;obslog - observational log

param,parvariable,parval,parcom,/get,pfname=pfname

if ~keyword_set(inputarray) then begin
readcol,obslog,fnameorig,fname,filter,exptime,panew,filetype,format='(a,a,a,f,f,a)'
;now pick out the frames that match the keyword
framesel=where(filetype eq keyword,nframes)


if nframes eq 0 then begin 
print,'No matching files of ',keyword,' found!'
goto,endoffile
endif

fnumber=fltarr(nframes)
for i=0L,nframes-1 do fnumber[i]=long((strsplit(strtrim(fname(framesel[i]),2),'n',/extract))[0])
endif else begin
fnumber=inputarray
nframes=n_elements(fnumber)
endelse

framestring=strtrim(long(fnumber[0]),2)

for i=1L,nframes-2 do begin

if long(fnumber[i]) eq long(fnumber[i-1])+1 then begin

if long(fnumber[i+1]) gt long(fnumber[i])+1 then begin
framestring+='-'+strtrim(long(fnumber[i]),2)
endif 
goto,endofloop
endif else begin

framestring+=','+strtrim(long(fnumber[i]),2)

endelse


endofloop:
endfor
if long(fnumber[nframes-1]) eq long(fnumber[nframes-2])+1 then begin
framestring+='-'+strtrim(long(fnumber[nframes-1]),2)
endif else begin
framestring+=','+strtrim(long(fnumber[nframes-1]),2)
endelse

param,parvariable,framestring,parcom,/set,pfname=pfname

endoffile:
end
