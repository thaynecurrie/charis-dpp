function deg2hms,deg

nn=n_elements(deg)
if nn eq 1 then return,(['-',' '])[deg ge 0]+string(sixty(abs(deg)/15.),format='(i02,1h:,i02,1h:,f05.2)')

out=strarr(nn)
for n=0,nn-1 do out[n]=(['-',' '])[deg[n] ge 0]+string(sixty(abs(deg[n])/15.),format='(i02,1h:,i02,1h:,f05.2)')
return,out

end
