pro charis_get_phot_zeropoint,filtername=fname,zeropoint=zpt,lambda=lambdao,$
verbose=verbose,help=help

if keyword_set(verbose) then Print,'The Filter Name string is case sensitive and stringent.  Type carefully'
if keyword_set(verbose) then Print,'**Warning: default is MKO filters, not 2MASS'

if keyword_set(help) then begin
Print,'Function to help determine the photometric zero point'
Print,''
Print,'Filters supported: YJHKLpMp (MKO), JHKs (2MASS), H20 and PAH (3.05 and 3.3)'
Print,'The Filter Name string is case sensitive and stringent.  Type carefully'
goto,skiptotheend
endif

;case (filtername eq 'Y' or  filtername eq 'Y MKO' or filtername eq 'YMKO')
case fname of
;zpt=2026

;**Y band**
'Y': begin 
 zpt=2026 
 lambdao=1.03
 end
'YMKO': begin 
zpt=2026
lambdao=1.03
end
'Y MKO': begin
zpt=2026
lambdao=1.03
end

;**J band**

'J': begin
zpt=1560
lambdao=1.25
end
'JMKO': begin
zpt=1560
lambdao=1.25
end
'J MKO': begin
zpt=1560
lambdao=1.25
end

'J 2MASS': begin
zpt=1594
lambdao=1.26
end
'J2MASS': begin
zpt=1594
lambdao=1.26
end

;**H band**

'H': begin
zpt=1050
lambdao=1.65
end
'HMKO': begin
zpt=1050
lambdao=1.65
end
'H MKO': begin
zpt=1050
lambdao=1.65
end

'H 2MASS': begin
zpt=1024
lambdao=1.66
end
'H2MASS': begin
zpt=1024
lambdao=1.66
end


;**K band**

'K': begin
zpt=667
lambdao=2.15
end
'KMKO': begin
zpt=667
lambdao=2.15
end
'K MKO': begin
zpt=667
lambdao=2.15
end

'K 2MASS': begin
zpt=666.7
lambdao=2.16
end
'K2MASS': begin
zpt=666.7
lambdao=2.16
end
'Ks': begin
zpt=666.7
lambdao=2.16
end


;**H20**
'H20': begin
zpt=356.
lambdao=3.05
end

'[3.05]': begin 
zpt=356.
lambdao=3.05
end


;**PAH**
'PAH': begin 
zpt=313.
lambdao=3.29
end
'[3.3]': begin
zpt=313.
lambdao=3.29
end


;**Lp**
'L': begin
zpt=248.
lambdao=3.78
end
'Lp': begin 
zpt=248.
lambdao=3.78
end

;**[4.05]
'Bralph': begin
zpt=207.
lambdao=4.05
end
'[4.05]': begin 
zpt=207.
lambdao=4.05
end

;'**Ms'
'Mp': begin 
zpt=163.
lambdao=4.68
end
'Ms': begin 
zpt=163.
lambdao=4.68
end

'Barr M': begin 
zpt=154.
lambdao=4.78
end
'Ml': begin 
zpt=154.
lambdao=4.78
end


endcase

skiptotheend:

end
