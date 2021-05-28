;---------------------------------------------------------------------------
; Invert a matrix matA by Singular Value Decomposition (SVD)
; If m is the maximal singular values, only the singular values larger
; than m*eps are kept during the inversion
; singVal contains the vector of singular values ordered in decreasing order
; dimKernel is the dimension of the kernel (number of singular values
; excluded).
; matKernel is the coordinate transformation matrix from the original basis
; into the kernel basis (coordinate of the kernel vectors)
; vecWKernel is the singular values of the kernel vectors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function svd_invert, matA, eps, dimKernel=dimKernel, matKernel=matKernel,vecWKernel = vecWKernel, singVal=singVal,double=double

; Step 1

svdc, matA, vecW, matU, matV, /column, double=double ; /column is the trick !

; Step 2

n = n_elements(vecW)
matWInv = dblarr(n,n)
m = max(vecW)
;print,'m',m
for i=0,n-1 do begin
    if vecW(i)/m ge eps then matWInv(i,i) = 1 / vecW(i)
endfor

; Step 3

matAInv = matV # matWInv # transpose(matU)
; checked !!!

; Step 4

singVal=vecW(reverse(sort(vecW)))
w = where( vecW/m lt eps)
if w(0) eq -1 then begin
    dimKernel = 0
    matKernel = 0d
endif else begin
    dimKernel = n_elements(w)
    matKernel = matV(*,w)
    vecWKernel = vecW(w)
endelse

return, matAInv

end
