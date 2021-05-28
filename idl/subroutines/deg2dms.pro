function deg2dms,deg

nn=n_elements(deg)
if nn eq 1 then return,(['-','+'])[deg ge 0]+string(sixty(abs(deg)),format='(i02,1h:,i02,1h:,f04.1)')

out=strarr(nn)
for n=0,nn-1 do out[n]=(['-','+'])[deg[n] ge 0]+string(sixty(abs(deg[n])),format='(i02,1h:,i02,1h:,f04.1)')
return,out

end
