pro plotchisqdist_charis,sample=sample,secondsample=secondsample,filein=filein,delay=delay,justind=justind,hk=hk,$
maskbad=maskbad
;plots chisq distribution for empirical comparisons

;*****UNSUPPORTED**** Use at your own risk and understand the source code and how to change it!
;CAUTION!!!! - Currently, only the Montreal Library is supported.  MUST EDIT SOURCE CODE TO GET BONNEFOY LIBRARY TO WORK!


;*****Default = Montreal Library *****
prefdir='./'
if ~keyword_set(filein) then filein='fitoutcomepick.txt'
if ~keyword_set(sample) then sample = 'montreal'

case sample of

  'montreal': begin
   if keyword_set(hk) then outputdir=prefdir+'Montreal_spectral_library_hk/' else outputdir=prefdir+'Montreal_spectral_library/'
   readcol,outputdir+filein,name,spt,grav,chisq,h20,h20d,h201,h202,h2k,hcont,format='(a,a,i,f,f,f,f,f,f,f)'
   
	      end

  'bonnefoy': begin
   if keyword_set(hk) then outputdir=prefdir+'bonnefoy_library_hk/' else outputdir=prefdir+'bonnefoy_library/'
   readcol,outputdir+filein,name,spt,grav,chisq,h20,h20d,h201,h202,h2k,hcont,format='(a,a,i,f,f,f,f,f,f,f)'
 
              end

  ;'cruztemplate': begin
  ;                end
endcase

;print,spt
;stop
charis_get_numspt,sptarrayin=spt,sptarrayout=numerical_spt

;sample='bonnefoy'

if keyword_set(secondsample) then begin

if ~keyword_set(filein) then filein='fitoutcomepick.txt'
;filein='fitoutcomepick.txt'

case sample of

  'montreal': begin
   if keyword_set(hk) then outputdir=prefdir+'Montreal_spectral_library_hk/' else outputdir=prefdir+'Montreal_spectral_library/'
   readcol,outputdir+filein,name2,spt2,grav2,chisq2,h202,h20d2,h2012,h2022,h2k2,hcont2,format='(a,a,i,f,f,f,f,f,f,f)'
   
	      end

  'bonnefoy': begin
   if keyword_set(hk) then outputdir=prefdir+'bonnefoy_library_hk/' else outputdir=prefdir+'bonnefoy_library/'
   readcol,outputdir+filein,name2,spt2,grav2,chisq2,h202,h20d2,h2012,h2022,h2k2,hcont2,format='(a,a,i,f,f,f,f,f,f,f)'
             print,outputdir+filein 
              end

  ;'cruztemplate': begin
  ;                end
endcase


charis_get_numspt,sptarrayin=spt2,sptarrayout=numerical_spt2
endif


;H20 plot

;goto,skipplothist
set_plot,'x'

;H20
good=where(abs(h20) lt 100 and numerical_spt ge 5 and numerical_spt le 14)
offseth20=median(numerical_spt[good]-h20[good],/even)
print,median(numerical_spt[good]-h20[good]-offseth20,/even)
print,stdev(numerical_spt[good]-h20[good]-offseth20)
print,'offseth20 is',offseth20,' c0 is ',-83.5437+offseth20

plothist,numerical_spt[good]-(h20[good]+offseth20),bin=0.5,xrange=[-10,10],peak=1
setcolors,/system_variables
plot,numerical_spt,h20+offseth20,xrange=[-10,30],yrange=[-10,30],xthick=5,ythick=5,psym=4
oplot,findgen(200)*0.1,findgen(200)*0.1,linestyle=0,color=!blue
;stop
if keyword_set(delay) then wait,5

;H20-D
good=where(abs(h20d) lt 100 and numerical_spt ge 10 and numerical_spt le 18)
offseth20d=median(numerical_spt[good]-h20d[good],/even)
print,median(numerical_spt[good]-h20d[good]-offseth20d,/even)
print,stdev(numerical_spt[good]-h20d[good]-offseth20d)
print,'offseth20d is',offseth20d,' c0 is ',79.4477+offseth20d

plothist,numerical_spt[good]-(h20d[good]+offseth20d),bin=0.5,xrange=[-10,10],peak=1
if keyword_set(delay) then wait,5


;H20-1
good=where(abs(h201) lt 100 and numerical_spt ge 4 and numerical_spt le 15)
offseth201=median(numerical_spt[good]-h201[good],/even)
print,median(numerical_spt[good]-h201[good]-offseth201,/even)
print,stdev(numerical_spt[good]-h201[good]-offseth201)
print,'offseth201 is',offseth201,' c0 is ',12.1927+offseth201

plothist,numerical_spt[good]-(h201[good]+offseth201),bin=0.5,xrange=[-10,10],peak=1
if keyword_set(delay) then wait,5


;H20-2
good=where(abs(h202) lt 100 and numerical_spt ge 4 and numerical_spt le 12)
offseth202=median(numerical_spt[good]-h202[good],/even)
print,median(numerical_spt[good]-h202[good]-offseth202,/even)
print,stdev(numerical_spt[good]-h202[good]-offseth202)
print,'offseth202 is',offseth202,' c0 is ',10.8822+offseth202

plothist,numerical_spt[good]-(h202[good]+offseth202),bin=0.5,xrange=[-10,10],peak=1

if keyword_set(delay) then wait,5



;stop
;plothist,numerical_spt-h20d,bin=0.5,xrange=[-10,10],peak=1
;plothist,numerical_spt-h201,bin=0.5,xrange=[-10,10],peak=1
;plothist,numerical_spt-h202,bin=0.5,xrange=[-10,10],peak=1
;stop

;*****Now, take the robust_mean of the stdev
skipplothist:


set_plot,'ps'
device,filename='h20range.eps',/encapsulated,bits=8,/color
setcolors,/system_variables
plot,numerical_spt,h20,xrange=[-10,30],yrange=[-10,30],xthick=5,ythick=5,psym=4,/nodata
oplot,numerical_spt,h20+offseth20,psym=4
oplot,findgen(200)*0.1,findgen(200)*0.1,linestyle=0,color=!blue
device,/close
;stop

device,filename='h20drange.eps',bits=8,/color
plot,numerical_spt,h20,xrange=[-10,30],yrange=[-10,30],xthick=5,ythick=5,psym=4,/nodata
setcolors,/system_variables
oplot,numerical_spt,h20d+offseth20d,psym=4
oplot,findgen(200)*0.1,findgen(200)*0.1,linestyle=0,color=!blue
device,/close

device,filename='h201range.eps',bits=8,/color
plot,numerical_spt,h20,xrange=[-10,30],yrange=[-10,30],xthick=5,ythick=5,psym=4,/nodata
setcolors,/system_variables
oplot,numerical_spt,h201+offseth201,psym=4
oplot,findgen(200)*0.1,findgen(200)*0.1,linestyle=0,color=!blue
device,/close

device,filename='h202range.eps',bits=8,/color
plot,numerical_spt,h20,xrange=[-10,30],yrange=[-10,30],xthick=5,ythick=5,psym=4,/nodata
setcolors,/system_variables
oplot,numerical_spt,h202+offseth202,psym=4
oplot,findgen(200)*0.1,findgen(200)*0.1,linestyle=0,color=!blue
device,/close

;****gravity indices

if keyword_set(maskbad) then begin
;badarray=[126,127,141,242,259,267,273,277,285,286,293,300,305,308,329,332,351]
;badarray=[141,242,259,267,273,277,285,286,293,300,305,308,329,332,351]
;badarray=[141,165,242,259,267,273,277,285,286,293,300,305,308,329,332,351]
;badarray=[141,165,242,259,267,273,277,285,286,291,293,300,305,308,329,332,351]
;badarray=[141,165,242,253,259,267,273,277,285,286,291,293,300,305,308,329,332,351]
;forhk
;if keyword_set(hk) then badarray[*]+=1
print,badarray
;stop
chisq[badarray]=!values.f_nan
h2k[badarray]=!values.f_nan
hcont[badarray]=!values.f_nan

bad=where(hcont gt 1.04 and hcont lt 1.4 and grav eq 0 and numerical_spt gt 9)
print,name[bad]
print,bad
;bandaid fix
if ~keyword_set(hk) then  hcont[bad]=!values.f_nan
;stop

endif

device,filename='hcont.eps',bits=8,/color
plot,numerical_spt,hcont,xrange=[0,20],yrange=[0.75,1.15],xthick=5,ythick=5,psym=4,/nodata
setcolors,/system_variables
field=where(grav eq 0)
int=where(grav eq 1)
low=where(grav eq 2)
vlow=where(grav eq 3)
plotsym,0,/fill
oplot,numerical_spt[field],hcont[field],psym=8,color=!gray,symsize=1.25
oplot,numerical_spt[int],hcont[int],psym=8,color=!blue,symsize=1.25
oplot,numerical_spt[low],hcont[low],psym=8,color=!green,symsize=1.25
device,/close


device,filename='h2k.eps',bits=8,/color
plot,numerical_spt,h2k,xrange=[0,20],yrange=[0.95,1.2],xthick=5,ythick=5,psym=4,/nodata,ystyle=1
setcolors,/system_variables
field=where(grav eq 0)
int=where(grav eq 1)
low=where(grav eq 2)
vlow=where(grav eq 3)
plotsym,0,/fill
oplot,numerical_spt[field],h2k[field],psym=8,color=!gray,symsize=1.25
oplot,numerical_spt[int],h2k[int],psym=8,color=!blue,symsize=1.25
oplot,numerical_spt[low],h2k[low],psym=8,color=!green,symsize=1.25
device,/close

comb=hcont+1*(1-h2k)

device,filename='comb.eps'
plot,numerical_spt,hcont+1*(1-h2k),xrange=[0,20],yrange=[0.5,1.25],xthick=5,ythick=5,psym=4,/nodata,ystyle=1
setcolors,/system_variables
field=where(grav eq 0)
int=where(grav eq 1)
low=where(grav eq 2)
vlow=where(grav eq 3)
plotsym,0,/fill
oplot,numerical_spt[field],hcont[field]+1*(1-h2k[field]),psym=8,color=!gray,symsize=1.25
oplot,numerical_spt[int],hcont[int]+1*(1-h2k[int]),psym=8,color=!blue,symsize=1.25
oplot,numerical_spt[low],hcont[low]+1*(1-h2k[low]),psym=8,color=!green,symsize=1.25
device,/close


;bad=where(grav eq 0 and comb gt 1.05)
;bad=where(grav eq 1 and comb gt 0.95 and numerical_spt gt 9)
;bad=where(grav eq 1 and h2k gt 1.1)
;print,name[bad]


;sptrange=['M0','M5','L0','L5','T0','T5']
;nsptrange=[0,5,10,15,20,25]

;sptrange=['K0','K5','M0','M5','L0','L5','T0','T5']
;nsptrange=[-10,-5,0,5,10,15,20,25]
sptrange=['G5','K0','K5','M0','M5','L0','L5','T0','T5']
nsptrange=[-15,-10,-5,0,5,10,15,20,25]

;sptrange=['M5','L0','L5','T0']
;nsptrange=[5,10,15,20]

;bad=where(grav eq 0 and abs(chisq) le 3.)
;print,name[bad]
;for i=0L,n_elements(bad)-1 do print,name[bad[i]],grav[bad[i]],chisq[bad[i]],bad[i]
;stop


;remove the one beta object with sizeable errors
if keyword_set(maskbad) then begin
badarray2=[126,127]
chisq[badarray2]=!values.f_nan
endif

if ~keyword_set(hk) then begin
good=where(abs(chisq) le 1.8,ngood)
endif else begin
good=where(abs(chisq) le 1.85,ngood)
endelse
print,'we are good!',spt[good],chisq[good],name[good]
;print,good
;stop
print,'we are good2!',numerical_spt[good]
;print,name[178],chisq[178],chisq[67],chisq[126]
print,ngood
if ngood gt 1 then begin
if ~keyword_set(hk) and n_elements(good) gt 0 then numerical_spt[good[1]]=10.25
endif




;bad=where(numerical_spt ge 8 and numerical_spt le 15 and grav eq 0 and abs(chisq) le 3)


;for i=0L,n_elements(bad)-1 do print,name[bad[i]],grav[bad[i]],' ',spt[bad[i]],chisq[bad[i]],bad[i],hcont[bad[i]],h2k[bad[i]]


comb=hcont+1*(1-h2k)

bad=where((comb gt 1.05 and comb lt 2 and grav eq 0) or comb gt 1.3)
hcont[bad]=!values.f_nan
h2k[bad]=!values.f_nan
comb[bad]=!values.f_nan

;print,''
;print,'mynameis',name[where(grav eq 0 and abs(chisq) le 1.5)]
;print,where(grav eq 0 and abs(chisq) le 1.5)
;print,'mynameis',name[where(grav eq 1 and abs(chisq) le 1.5)]
;stop
;stop

device,filename='chisqdist.eps'
!p.font=1
if ~keyword_set(hk) then begin
plot,numerical_spt,chisq,xrange=[min(nsptrange),max(nsptrange)],yrange=[0.2,400],/ylog,ystyle=1,/nodata,$
;plot,numerical_spt,chisq,xrange=[0,20],yrange=[0.75,100],/ylog,ystyle=1,/nodata,$
;plot,numerical_spt,chisq,xrange=[0,20],yrange=[0,10],ystyle=1,/nodata,$
;xtickname=nsptrange,xticks=5
xtickname=sptrange,xtickv=nsptrange,charthick=2.5,xthick=5,ythick=5,charsize=1.5,$
xtitle='Spectral Type',ytitle=textoidl('\chi^{2}_{\nu}'),$
xticks=8
endif else begin
plot,numerical_spt,chisq,xrange=[min(nsptrange),max(nsptrange)],yrange=[0.5,200],/ylog,ystyle=1,/nodata,$
;plot,numerical_spt,chisq,xrange=[0,20],yrange=[0.75,100],/ylog,ystyle=1,/nodata,$
;plot,numerical_spt,chisq,xrange=[0,20],yrange=[0,10],ystyle=1,/nodata,$
;xtickname=nsptrange,xticks=5
xtickname=sptrange,xtickv=nsptrange,charthick=2.5,xthick=5,ythick=5,charsize=1.5,$
xtitle='Spectral Type',ytitle=textoidl('\chi^{2}_{\nu} (HK)')
endelse
setcolors,/system_variables
oplot,numerical_spt[field],chisq[field],psym=8,color=!gray
oplot,numerical_spt[int],chisq[int],psym=8,color=!blue,symsize=1.25
oplot,numerical_spt[low],chisq[low],psym=8,color=!green,symsize=1.25
oplot,numerical_spt[vlow],chisq[vlow],psym=8,color=!orange,symsize=1.25
if keyword_set(secondsample) then begin
oplot,numerical_spt2,chisq2,psym=8,color=!magenta,symsize=1.25
endif

;oplot,findgen(100)*0.25,findgen(100)*0+1.721,linestyle=1,thick=3

;if ~keyword_set(hk) then begin
;oplot,findgen(100)*0.25,findgen(100)*0+1.67,linestyle=1,thick=3
;endif else begin
;oplot,findgen(100)*0.25,findgen(100)*0+1.86,linestyle=1,thick=3
;endelse

setcolors,/system_variables
al_legend,[textoidl('MTL Field'),textoidl('MTL \beta'),textoidl('MTL \gamma'),textoidl('MTL \delta')],color=[!gray,!blue,!green,!orange],box=0,psym=8,charthick=2,/right

;al_legend,[textoidl('MTL Field'),textoidl('MTL \beta'),textoidl('MTL \gamma'),textoidl('MTL \delta'),$
;textoidl('Bonnefoy Lib.')],color=[!gray,!blue,!green,!orange,!magenta],box=0,psym=8,charthick=2,/right

setcolors,/system_variables
plotsym,3,/fill
;al_legend,[textoidl('\kappa And b')],color=!cyan,psym=8,symsize=1.75,position=[19.6,42],box=0,charthick=2,/right
device,/close

;****gravity indices

device,filename='hcont.eps',bits=8,/color
;plot,numerical_spt,hcont,xrange=[5,20],yrange=[0.75,1.15],psym=4,/nodata,$
plot,numerical_spt,hcont,xrange=[0,20],yrange=[0.75,1.15],psym=4,/nodata,ystyle=1,$
xtickname=sptrange,xtickv=nsptrange,charthick=2.5,xthick=5,ythick=5,charsize=1.5,$
xtitle='Spectral Type',ytitle=textoidl('H_{cont.}')
;xtickname=sptrange,xtickv=nsptrange,charthick=2.5,xthick=5,ythick=5,charsize=1.25
setcolors,/system_variables
field=where(grav eq 0)
int=where(grav eq 1)
low=where(grav eq 2)
vlow=where(grav eq 3)
plotsym,0,/fill
oplot,numerical_spt[field],hcont[field],psym=8,color=!gray
oplot,numerical_spt[int],hcont[int],psym=8,color=!blue,symsize=1.25
oplot,numerical_spt[low],hcont[low],psym=8,color=!green,symsize=1.25
oplot,numerical_spt[vlow],hcont[vlow],psym=8,color=!orange,symsize=1.25
if keyword_set(secondsample) then begin
oplot,numerical_spt2,hcont2,psym=8,color=!magenta,symsize=1.25
endif
setcolors,/system_variables
al_legend,[textoidl('MTL Field'),textoidl('MTL \beta'),textoidl('MTL \gamma'),textoidl('MTL \delta'),$
textoidl('Bonnefoy Lib.')],color=[!gray,!blue,!green,!orange,!magenta],box=0,psym=8,charthick=2,/right
setcolors,/system_variables
plotsym,3,/fill
;kap and: 1.070 +/- 0.039
oploterror,[10.5],[1.070],[0.039],psym=8,color=!cyan,symsize=2.5,errthick=5
;oplot,[10.5],[1.070],psym=8,color=!cyan,symsize=2.
al_legend,[textoidl('\kappa And b')],color=!cyan,psym=8,symsize=1.75,position=[19.6,1.0425],box=0,charthick=2,/right
device,/close


device,filename='h2k.eps',bits=8,/color
;plot,numerical_spt,h2k,xrange=[5,20],yrange=[0.9,1.2],psym=4,/nodata,$
plot,numerical_spt,h2k,xrange=[0,20],yrange=[0.95,1.2],psym=4,/nodata,$
xtickname=sptrange,xtickv=nsptrange,charthick=2.5,xthick=5,ythick=5,charsize=1.5,$
xtitle='Spectral Type',ytitle=textoidl('H_{2}K')
setcolors,/system_variables
field=where(grav eq 0)
int=where(grav eq 1)
low=where(grav eq 2)
vlow=where(grav eq 3)
plotsym,0,/fill
oplot,numerical_spt[field],h2k[field],psym=8,color=!gray
oplot,numerical_spt[int],h2k[int],psym=8,color=!blue,symsize=1.25
oplot,numerical_spt[low],h2k[low],psym=8,color=!green,symsize=1.25
oplot,numerical_spt[vlow],h2k[vlow],psym=8,color=!orange,symsize=1.25
if keyword_set(secondsample) then begin
oplot,numerical_spt2,h2k2,psym=8,color=!magenta,symsize=1.25
endif
setcolors,/system_variables
al_legend,[textoidl('MTL Field'),textoidl('MTL \beta'),textoidl('MTL \gamma'),textoidl('MTL \delta'),$
textoidl('Bonnefoy Lib.')],color=[!gray,!blue,!green,!orange,!magenta],box=0,psym=8,charthick=2,/left
setcolors,/system_variables
plotsym,3,/fill
oploterror,[10.5],[1.055],[0.041],psym=8,color=!cyan,symsize=2.5,errthick=5
;oplot,[10.5],[1.055],psym=8,color=!cyan,symsize=2.
al_legend,[textoidl('\kappa And b')],color=!cyan,psym=8,symsize=1.75,position=[0.36,1.1325],box=0,charthick=2,/left
device,/close



device,filename='comb.eps'
;plot,numerical_spt,hcont+1*(1-h2k),xrange=[5,20],yrange=[0.5,1.25],xthick=5,ythick=5,psym=4,/nodata
plot,numerical_spt,hcont+1*(1-h2k),xrange=[0,20],yrange=[0.6,1.15],ystyle=1,psym=4,/nodata,$
xtickname=sptrange,xtickv=nsptrange,charthick=2.5,xthick=5,ythick=5,charsize=1.5,$
xtitle='Spectral Type',ytitle=textoidl('H_{cont}+(1-H_{2}K)')
setcolors,/system_variables
field=where(grav eq 0)
int=where(grav eq 1)
low=where(grav eq 2)
vlow=where(grav eq 3)
plotsym,0,/fill
oplot,numerical_spt[field],hcont[field]+1*(1-h2k[field]),psym=8,color=!gray
oplot,numerical_spt[int],hcont[int]+1*(1-h2k[int]),psym=8,color=!blue,symsize=1.25
oplot,numerical_spt[low],hcont[low]+1*(1-h2k[low]),psym=8,color=!green,symsize=1.25
oplot,numerical_spt[vlow],hcont[vlow]+1*(1-h2k[vlow]),psym=8,color=!orange,symsize=1.25
if keyword_set(secondsample) then begin
oplot,numerical_spt2,hcont2+1*(1-h2k2),psym=8,color=!magenta,symsize=1.25
endif
setcolors,/system_variables
al_legend,[textoidl('MTL Field'),textoidl('MTL \beta'),textoidl('MTL \gamma'),textoidl('MTL \delta'),$
textoidl('Bonnefoy Lib.')],color=[!gray,!blue,!green,!orange,!magenta],box=0,psym=8,charthick=2,/right
setcolors,/system_variables
plotsym,3,/fill

;kap And b
;oploterror,[10.5],[1.07+(1-1.055)],[0.057],psym=8,color=!cyan,symsize=2.5,errthick=5
;oplot,[10.5],[1.07+(1-1.055)],psym=8,color=!cyan,symsize=2.
;al_legend,[textoidl('\kappa And b')],color=!cyan,psym=8,symsize=1.75,position=[19.6,1.0],box=0,charthick=2,/right
device,/close

;bah=where(numerical_spt gt 5 and numerical_spt lt 7 and h2k lt 0.97 and h2k gt 0.95)

;limit=1.67*1.25
;limit=1.86*1.25
;limit=2.5*1.25

;limit=2.3125
;limit=2.69
;limit=2.435
limit=1.67

;good=where(abs(chisq) le 1.721)
good=where(abs(chisq) le limit)
for i=0L,n_elements(good)-1 do print,name[good[i]],numerical_spt[good[i]],grav[good[i]],chisq[good[i]]

if keyword_set(secondsample) then begin
good=where(abs(chisq2) le limit)
for i=0L,n_elements(good)-1 do print,name2[good[i]],numerical_spt2[good[i]],chisq2[good[i]]
endif
;,grav[good[i]]

;Now print top 10 sources
namer=name(sort(abs(chisq)))
sptr=spt(sort(abs(chisq)))
chisqr=chisq(sort(abs(chisq)))

for i=0L,n_elements(good)-1 do print,i,namer[i],' ',sptr[i],chisqr[i]

end
