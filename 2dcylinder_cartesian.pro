;------------------------------------------------------
;-- 2D Parallel Plate Capacitor -----------------------
;------------------------------------------------------
;-- This code solves Laplace's Equation for a ---------
;-- parallel plate capacitor on a 2D grid -------------
;-- with the boundary grounded (V=0). -----------------
;------------------------------------------------------
;-- written by M.R. Stoneking for ---------------------
;-- Physics 23, Spring Term 2000 ----------------------
;------------------------------------------------------

;-- user sets parameters in this section --------------
imax = 100		;-- horizontal grid index maximum (must be even)
jmax = 100		;-- vertical grid index maximum (must be even)
Vo = 10.		;-- potential difference between cap plates.
			;-- top plate will be at +Vo/2, bottom at -Vo/2.
maxits = 1000	;-- maximum number of passes through the grid
eps = 1.e-7	;-- tolerance in total error
omega = 1.8	;-- over-relaxation parameter
d = 8			;-- plate separation (must be an even number of grid points).
w = 40		;-- width of cap plates (must be an even number of grid points)
			;-- note that the plates are always two grid points thick.
;------------------------------------------------------


;-- dimension arrays  ---------------------------------
V = fltarr(imax+1,jmax+1)		;-- potential array
jflag = intarr(imax+1,jmax+1)	;-- array of flags to indicate if
						;-- grid site is to be updated (jflag=0)
						;-- or fixed (jflag=1).  All boundary sites
						;-- must be fixed as well as cap plates.
err_tot = fltarr(maxits)		;-- array of total error values... permits tracking
						;-- of convergence toward best solution.
Ex = fltarr(imax+1,jmax+1)		;-- x-component of the electric field.
Ey = fltarr(imax+1,jmax+1)		;-- y-component of the electric field.
Etot = fltarr(imax+1,jmax+1)	;-- magnitude of the electric field.
;-------------------------------------------------------


;-- set boundary conditions ----------------------------
for i=0,imax do begin

	V(i,jmax) = 0.		;-- ground top boundary
	V(i,0) = 0.		;-- ground bottom boundary
	V(0,i) = 0.		;-- ground left boundary
	V(imax,i) = 0.		;-- ground right boundary
	
	V(i,jmax/2) = 0.	;-- ground the midplane. The potential must be
					;-- antisymmetric about the midplane.  We will 
					;-- relax toward solution on top half plane only
					;-- and make use of symmetry to reduce copmutation time.							
endfor

ileft = (imax-w)/2		;-- location of left end of plate
iright = (imax+w)/2		;-- location of right end of plate
jup = (jmax+d)/2		;-- location of bottom face of top plate
jdown = (jmax-d)/2		;-- location of top face of bottome plate

for i=ileft,iright do begin

	V(i,jup) = Vo/2.		;-- top plate is at half the potential difference.
	V(i,jup+1) = Vo/2.
	jflag(i,jup) = 1		;-- do not want to update the potential at these
	jflag(i,jup+1) = 1		;-- grid sites.
	V(i,jdown) = -Vo/2.		;-- bottom plate is at minus half the potential diff.
	V(i,jdown-1) = -Vo/2.	

endfor	
;-------------------------------------------------------

;-- relaxation section ---------------------------------
for l=0,(maxits-1) do begin

	err_it = 0.			;-- reset accumulated error for this iteration.
	
	for i=1,(imax-1) do begin	;-- do not include boundaries in relaxation (i=0 & imax)
		
		for j=(jmax/2),(jmax-1) do begin ;--only pass over the top half plane
		
		if (jflag(i,j) eq 1) then goto, nextone
			;-- check to see if grid site is fixed.
			
		delta = 0.25*(V(i-1,j)+V(i+1,j)+V(i,j-1)+V(i,j+1))-V(i,j)
			;-- difference between current value and current best value
						
		err_it = err_it + abs(delta)
			;-- accumulate error for this iteration
			
		V(i,j) = V(i,j) + omega*delta	;-- update grid site
		V(i,jmax-j) = - V(i,j)			;-- make use of symmetry to set
								;-- values in bottome half plane

		nextone:		
				
		endfor
		
	endfor
	
	err_tot(l) = err_it

	;-- check convergence -------	
	if (l gt 0) then begin
	if abs((err_tot(l-1)-err_tot(l))/err_tot(0)) lt eps then goto, done
	endif
	
endfor
print,'maximum iterations exceeded'
done:
print,'convergence criterion met'
;------------------------------------------------------------------

;-- calculate the electric field.----------------------------------
for i=1,(imax-1) do begin

	for j=1,(jmax-1) do begin
	
	Ex(i,j) = .5*(V(i+1,j)-V(i-1,j))
	Ey(i,j) = .5*(V(i,j+1)-V(i,j-1))
	Etot(i,j) = sqrt(Ex(i,j)^2+Ey(i,j)^2)

	endfor
endfor
;------------------------------------------------------------------


;-- plotting section ----------------------------------------------


;------------------------------------------------------------------
	

end

