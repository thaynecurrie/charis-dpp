pro charis_get_numspt,sptarrayin=sptarrayin,sptarrayout=sptarrayout,name=name

;takes an array of Spectral Types, Returns Numerical Spectral Types

spt=sptarrayin

 numspectype=intarr(n_elements(spt))
 numspectype[*]=-999999
 subspectype=fltarr(n_elements(spt))
 numerical_spt=fltarr(n_elements(spt))



for i=0L,n_elements(spt)-1 do begin
 sg=(spt[i])
 ;K,M,L,T=[-10,0,10,20]

 ;K
 mtype=strmatch(sg,'*K*')
 ;print,mtype
if mtype ne 0 then begin
numspectype[i]=-10
sspectype='K'
goto,breakout
 endif
 ;M
 mtype=strmatch(sg,'*M*')
 if mtype ne 0 then begin
 numspectype[i]=0
sspectype='M'
 goto,breakout
 endif
 ;L
 mtype=strmatch(sg,'*L*')
 if mtype ne 0 then begin
 numspectype[i]=10
 sspectype='L'
 goto,breakout
 endif
 ;T
 mtype=strmatch(sg,'*T*')
 if mtype ne 0 then begin
 numspectype[i]=20
sspectype='T'
 goto,breakout
 endif
 breakout:

 subs=(strsplit(sg,sspectype,/extract))
;[n_elements(sg)-1]
 subspectype[i]=subs[n_elements(subs)-1]

numerical_spt[i]=numspectype[i]+subspectype[i]
if keyword_set(name) then print,name[i],numerical_spt[i],i
if ~keyword_set(name) then print,numerical_spt[i],i
endfor

sptarrayout=numerical_spt

END


