function density,r,theta,z,par
;+
; PURPOSE : compute the density of particles in the unit of the normalized density value
;    (generaly in : number of particles / cm^3).
;
; INPUT:
;       r : radius in AU
;       theta : polar angle in radians between -Pi and Pi
;       z : distance to the midplane in AU
;       par : structure containing all relevant informations on the density profile
;
; OUTPUT : density
;-

type = par.type ; type of model of density profile

case type of
    
    ;=========================================================
    ;-- Two power law density profile:
    '2powerlaw': begin
        ;print,'-> Two power law density profile'
        x = r / par.r0
        rho = par.density0 * $
          sqrt(2.d0) / sqrt(x^(-2.d0*par.alphain)+x^(-2.d0*par.alphaout)) * $
          exp(-(abs(z)/par.ksi0/x^par.beta)^par.gamma)
        id = where(float(r) lt float(par.rmin))
        if (size(id))(0) ne 0 then rho(id) = 0.d0
        id = where(float(r) gt float(par.rmax))
        if (size(id))(0) ne 0 then rho(id) = 0.d0
    end

    ;=========================================================
    ;-- Gaussian density profile:
    'gaussian': begin
        ;print,'-> Gaussian density profile'
        rho = par.density0 * $
          exp(-(r-par.r0)^2/2.d0/par.dx^2) * $
          exp(-(abs(z)/par.ksi0/(r/par.r0)^par.beta)^par.gamma)
    end


endcase


return,double(rho)
end
