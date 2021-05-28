function rot_ratef,h,dec,phideg,r=r,dx=dx
;print,'phi is ',phi,phi*!dtor
phi=phideg*!dtor
;usage r=rot_rate(h,dec)
;returns the field of view rotation rate in radian/sec for an object
;of declination dec at an hour angle h, observed from Mauna Kea
;
;h=hour angle in hours
;dec=declination in decimal degrees
;
;or
;usage: t=rot_rate(h,dec,r=r,dx=dx)
;returns the time (in sec) required for the field of view rotation to induce
;a displacement of dx at an angular separation of r when observing a target
;of declination dec at an hour angle h
;r and dx must be in same units

;retourne le taux de rotation du champ en rad/sec pour
;un objet de declinaison dec, a un angle horaire h
;au mauna kea
;
;h en heures
;dec en degres
;
;pour avoir le temp requis pour un deplacement de dx a un rayon donne, specifier:
;r=r en asec
;dx=deplacement en asec

;latitude
;phi=19.823801447d0*!dtor
;phi=19.825d0*!dtor
;phi=ten(31,41,19.6)*1.d0*!dtor

;angle au zenith
z=acos( (sin(dec*!dtor)*sin(phi)+cos(dec*!dtor)*cos(phi)*cos(h*15.*!dtor))<1.>(-1.) )

;angle d'azimuth
a=acos( ((cos((90.-dec)*!dtor)-sin(phi)*cos(z)) / (cos(phi)*sin(z)))<1.>(-1.) )

;angle de rotation en rad/sec
w=7.2925e-5*cos(a)*cos(phi)/sin(z)
if ~keyword_set(dx) or ~keyword_set(r) then return,w

;angle qu'on doit tourner en asec
theta=(float(dx)/r)
;temp requis pour tourner de ce temp la
t=abs(theta/w)
return,t

end
