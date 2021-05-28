pro charis_add_files,pfname,parvariable,filetoadd,stringinput=stringinput

;pfname - parameter file
;parvariable - the variable
;[filestoadd] - the file numbers to add
;stringinput - input the new files as a string array then 

param,parvariable,parval,parcom,/get,pfname=pfname

origarray=nbrlist(parval)
if keyword_set(stringinput) then filetoadd=nbrlist(filetoadd)
updatedarray=[origarray,filetoadd]

obslog='none'

charis_set_files,pfname,parvariable,keyword,obslog,inputarray=updatedarray

end
