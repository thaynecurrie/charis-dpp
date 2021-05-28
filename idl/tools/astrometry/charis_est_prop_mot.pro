pro charis_est_prop_mot,object=object,startdate=startdate,nyears=nyears,ndays=ndays,position=position,inffac=inffac,dpc=dpc,aberration=aberration,flipxaxis=flipxaxis,posfile=posfile,ql=ql,tozero=tozero,verbose=verbose,plotranges=plotranges


;*****Alpha Version!!!!

if ~keyword_set(object) then object='HD 1160'
;sandbox to do a proper motion estimation

if ~keyword_set(startdate) then startdate='10-18-2018'


hci_querysimbadf,object,ra,dec,id,found=found

ship=hci_queryvizierf('I/239/hip_main',[ra,dec],20./60.,/allcolumns,found=found)
;ra and dec in degrees

;parallax and proper motion pull
;Parallax and error on parallax might not have been found
plxset = where(strmatch( tag_names(ship), 'plx', /fold_case) eq 1)
pmraset = where(strmatch( tag_names(ship), 'pmra', /fold_case) eq 1)
pmdecset = where(strmatch( tag_names(ship), 'pmde', /fold_case) eq 1)


if (plxset gt 0) then plx=ship.plx else begin
    plx=0
    print, 'Parallax not found. Setting parallax=0.'
endelse

plx*=1d-3 ;to arc-sec

if keyword_set(dpc) then plx=1./dpc

;sometimes the HIP catalog triggers a stupid error ...
if n_elements(plx) gt 1 then plx=plx[0]

if (pmraset gt 0) then begin
pmra=ship.pmra
endif else begin
pmra=0
print,'NO RA PM!'
endelse
print,pmra
if (pmdecset gt 0) then begin
pmdec=ship.pmde
endif else begin
pmdec=0
print,'NO DEC PM!'
endelse

ymd0=double(strsplit(startdate,'-',/extract))
jd0=julday(ymd0[0], ymd0[1], ymd0[2])
epoch0 = (jd0 - 2451545d0)/365.256363d0 + 2000d0


ra0=ra
dec0=dec

ra0precess=ra0
dec0precess=dec0
ra0precess=ra0+pmra/3.6d6/cos(dec0*!dtor)*(epoch0-2000d0)
dec0precess=dec0+pmdec/3.6d6*(epoch0-2000d0)
precess,ra0precess,dec0precess,2000d0,epoch0

;now you have the real ra and dec: 

rastar=ra0precess
decstar=dec0precess
rastar=ra0
decstar=dec0

print,rastar,decstar
;stop
;convert form celestial coordinates to ecliptic coordinates
euler,rastar,decstar,lambda,beta,3

;**** Loop****

;do date loop over 2 calendar years

if keyword_set(nyears) then begin

ndays=365.25*nyears
endif else begin
if ~keyword_set(ndays) then ndays=365
endelse

parday=fltarr(ndays)
fracday=fltarr(ndays)

jdnew=jd0

del_lambda=fltarr(ndays)
del_beta=fltarr(ndays)

pmraf=pmra*1d-3/365.256363d0 ;proper motion per day
pmdecf=pmdec*1d-3/365.256363d0 ;proper motion per day

;pmraf/=cos(decstar*!pi/180.)

if ~keyword_set(inffac) then inffac=1.0

for i=0L,ndays-1 do begin
;for i=1L,ndays-1 do begin
;julian date
parday[i]=jdnew

;compute Sun's position on JD
sunpos,jdnew,rasun,decsun
euler,rasun,decsun,lamsun,betsun,3

;compute the fractional date
fracday[i]=(jdnew - 2451545d0)/365.256363d0 + 2000d0

;the calendar year
caldat,double(jdnew),month,day,year
;if keyword_set(verbose) then print,'month day year',i,month,day,year,fracday[i]


;****Now, get the offset

;offset in RA & DEC
;del_lambda[i]=1*plx*sin((lamsun-lambda)*!pi/180.)/cos(beta*!pi/180.) + pmraf*double(i)*inffac
del_lambda[i]=1*plx*sin((lamsun-lambda)*!pi/180.)/cos(0*beta*!pi/180.) + pmraf*double(i)*inffac

;del_lambda[i]=1*plx*sin((lamsun-lambda)*!pi/180.) + pmraf*double(i)*inffac
;del_lambda[i]=1*plx*sin((lamsun-lambda)*!pi/180.)/cos(1*beta*!pi/180.) + pmraf*double(i)*inffac/cos(1*beta*!pi/180.)

del_beta[i]=-1*plx*cos((lamsun-lambda)*!pi/180.)*sin(beta*!Pi/180.) + pmdecf*double(i)*inffac

;del_lambda[i]=plx*sin((rastar-rasun)*!pi/180.)/cos(decstar*!pi/180.) + pmraf*double(i)*inffac
;del_beta[i]=plx*cos((rastar-rasun)*!pi/180.) + pmdecf*double(i)*inffac


dellamo=plx*sin((lamsun-lambda)*!pi/180.)/cos(beta*!pi/180.)
delbeto=plx*cos((lamsun-lambda)*!pi/180.)*sin(beta*!pi/180.) + pmdecf*double(i)*inffac

;to convert to the motion of a background star

;del_lambda[i]*=-1.0
;del_beta[i]*=-1.0

;print,'month day year',i,month,day,year,' jd ',jdnew,' delta ',del_lambda[i],del_beta[i],dellamo,delbeto,rasun,betsun,lambda

if keyword_set(verbose) then print,'stuff!',i,month,day,year,del_lambda[i],del_beta[i]


;index the date by one day
jdnew+=1


endfor

del_lambda[*]-=del_lambda[0]
del_beta[*]-=del_beta[0]

offsetfromstar=fltarr(ndays,ndays)
xoffset=fltarr(ndays)
yoffset=fltarr(ndays)

if ~keyword_set(position) then begin
position=[0,0]
endif

xoffset=position[0]-del_lambda
yoffset=position[1]-del_beta


if keyword_set(ql) then begin
;xoffset=-1*xoffset
window,0
set_plot,'x'
plot,del_lambda,del_beta,xrange=[max(del_lambda)+0.05,min(del_lambda)-0.05],yrange=[min(del_beta)-0.05,max(del_beta)+0.05]
;plot,del_lambda,del_beta,xrange=[max(del_lambda)+0.05,min(del_lambda)-0.05],yrange=[min(del_beta)-0.05,max(del_beta)+0.05]
;plot,del_lambda,del_beta,xrange=[-0.6,-0.95],yrange=[4.63,4.89]
plot,xoffset,yoffset,xrange=[max(xoffset),min(xoffset)],yrange=[min(yoffset),max(yoffset)]

window,1
plot,del_lambda,del_beta,xrange=[max(del_lambda)*1.1,min(del_lambda)*0.9],yrange=[min(del_beta)*0.9,max(del_beta)*1.1]

endif

set_plot,'ps'
device,filename='prop_mot.eps',/encapsulated,/color,bits=8
!p.font=1
!p.thick=2
!x.thick=5
!y.thick=5

minyoffset=min(yoffset)
maxyoffset=max(yoffset)
minxoffset=min(xoffset)
maxxoffset=max(xoffset)
if keyword_set(tozero) then begin
minxoffset = minxoffset < 0
maxoffset =maxxoffset > 0
minyoffset=minyoffset < 0
maxyoffset = maxyoffset > 0
endif

if keyword_set(plotranges) then begin
maxxoffset=plotranges[0]
minxoffset=plotranges[1]
maxyoffset=plotranges[2]
minyoffset=plotranges[3]

endif

if keyword_set(posfile) then begin
setcolors,/system_variables
readcol,posfile,epos,npos,sigepos,signpos,dateobs,format='(f,f,f,f,a)'
help,dateobs

;now recompute all the del-lambda stuff
;first, convert back to original
del_lambda[*]+=(del_lambda[0])
del_beta[*]+=del_beta[0]

print,n_elements(epos)

;now convert from UT observing date to JD to fractional date
fracdayobs=fltarr(n_elements(epos))
predastrom=fltarr(2,2,n_elements(epos))
nepos=n_elements(epos)


ymd0=double(strsplit(dateobs[0],'-',/extract))
jd0=julday(long(ymd0[1]), long(ymd0[2]), long(ymd0[0]))
fracdaystart=(jd0 - 2451545d0)/365.256363d0 + 2000d0
diff=abs(fracdaystart-fracday)
good=where(diff eq min(diff))

;help,del_lambda
;print,del_lambda[good],del_beta[good]
del_lambda[*]-=(del_lambda[good])[0]
del_beta[*]-=(del_beta[good])[0]

;print,del_lambda
;help,del_lambda
;stop

xoffset=position[0]-del_lambda
yoffset=position[1]-del_beta

;xoffset-=xoffset[good]
;yoffset-=yoffset[good]

for i=0L,n_elements(epos)-1 do begin
ymd0=double(strsplit(dateobs[i],'-',/extract))
jd0=julday(long(ymd0[1]), long(ymd0[2]), long(ymd0[0]))
fracdayobs[i]=(jd0 - 2451545d0)/365.256363d0 + 2000d0
diff=abs(fracdayobs[i]-fracday)
good=where(diff eq min(diff))
predastrom(0,0,i)=xoffset[good]
predastrom(1,0,i)=yoffset[good]
predastrom(0,1,i)=epos[i]
predastrom(1,1,i)=npos[i]
print,'pred obs ',i,' ',dateobs[i],fracdayobs[i],fracday[good],predastrom(0,0,i),predastrom(1,0,i),predastrom(0,1,i),predastrom(1,1,i)
endfor
print,'starting position is ',xoffset[0],yoffset[0]


endif

if ~keyword_set(flipxaxis) then begin
plot,xoffset,yoffset,xrange=[maxxoffset,minxoffset],yrange=[minyoffset,maxyoffset],$
/nodata,xtitle=textoidl('E offset (arc-sec)'),ytitle=textoidl('N offset (arc-sec)'),xminor=4,yminor=4,xthick=5,ythick=5,charsize=1.5,xstyle=1,ystyle=1
oplot,xoffset,yoffset,linestyle=0,thick=5
endif else begin
plot,xoffset,yoffset,xrange=[-1*minxoffset,-1*maxxoffset],yrange=[minyoffset,maxyoffset],$
/nodata,xtitle=textoidl('W offset (arc-sec)'),ytitle=textoidl('N offset (arc-sec)'),xminor=4,yminor=4,xthick=5,ythick=5,charsize=1.5
oplot,-1*xoffset,yoffset,linestyle=0,thick=5
endelse

setcolors,/system_variables
if keyword_set(posfile) then begin
plotsym,0,/fill

for i=1L,nepos-1 do begin
oplot,predastrom[0,*,i],predastrom[1,*,i],linestyle=2,thick=2,color=!black
endfor

oploterror,epos[0],npos[0],sigepos[0],signpos[0],color=!magenta,psym=8,thick=5,symsize=1.5

oploterror,epos[1:nepos-1],npos[1:nepos-1],sigepos[1:nepos-1],signpos[1:nepos-1],color=!green,psym=8,thick=5,symsize=1.5
endif

;legend
setcolors,/system_variables
al_legend,['first epoch position','subsequent positions'],color=[!magenta,!green],box=0,/right,psym=8,symsize=1.5,charsize=1.5
;al_legend,['first epoch position','subsequent positions'],color=[!magenta,!green],box=0,/left,psym=8,symsize=1.5,charsize=1.5
al_legend,['background object path'],position=[0.,-0.05],box=0,/right,charsize=1.5,linestyle=0,thick=5

oplot,[1.615],[0.252],psym=4,symsize=2
al_legend,['Predicted Position (Oct. 2015)'],box=0,position=[1.7,.1],charsize=1.5,psym=4,symsize=2

nfac=1

;now, estimate the chi-squared 
diffx=predastrom[0,0,1]-predastrom[0,1,1]
diffy=predastrom[1,0,1]-predastrom[1,1,1]

chisq=diffx^2/(nfac*sigepos[1])^2.+diffy^2./(nfac*signpos[1])^2.

print,chisq
print,diffx,diffy
print,predastrom(0,0,1),predastrom(1,0,1)
device,/close

end
