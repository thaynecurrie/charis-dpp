pro charis_del_files,pfname,parvariable,filestodel,stringinput=stringinput

;pfname - parameter file
;parvariable - the variable
;[filestoadd] - the file numbers to add
;stringinput - input the new files as a string array then 

param,parvariable,parval,parcom,/get,pfname=pfname

origarray=nbrlist(parval)
if keyword_set(stringinput) then filestodel=nbrlist(filetoadd)
;todelete=uniq(intersect(origarray,filestodel))
;todelete=where(origarray eq filestodel)
match,origarray,filestodel,todelete
;print,origarray[todelete]
;print,origarray
remove,todelete,origarray
updatedarray=origarray
;print,' '
;print,updatedarray

obslog='none'

charis_set_files,pfname,parvariable,keyword,obslog,inputarray=updatedarray

end
