;*****a second attempt with MPFITFUN
FUNCTION GBFUNC,X,P
;FUNCTION GBFUNC, lam_ratio_sortf,[arhof,alamf,sigrhof,siglamf]
YMOD=P[0]*exp(-0.5*(P[4]*X/P[2])^2.)+P[1]*exp(-0.5*(X/P[3])^2.)
;YMOD=arhof*exp(-0.5*(rho*lam_ratio_sortf/sigrho)^2.)+alamf*exp(-0.5*(lam_ratio_sortf/siglam)^2.)
return, YMOD
END
