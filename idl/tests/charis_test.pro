pro charis_test

;very simple program to test whether you have your path installed correctly

print,'Checking for correct CHARIS-DPP Path'

charis_which,'charis_get_constant',/quiet,output=setupout

if strlen(setupout) gt 0 then begin 
print,'Setup test program found: ',setupout 
print,'CHARIS-DPP Installation likely successful'
endif else begin 
print,'Setup test program not found!'
print,'CHARIS-DPP Installation likely failed'
endelse


end 
