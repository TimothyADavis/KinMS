;##############################################################################
;
; Copyright (C) 2016, Timothy A. Davis
; E-mail: DavisT -at- cardiff.ac.uk
;
; Updated versions of the software are available through github:
; https://github.com/TimothyADavis/KinMS
; 
; If you have found this software useful for your research,
; I would appreciate an acknowledgment to the use of the
; "KINematic Molecular Simulation (KinMS) routines of Davis et al., (2013)".
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;##############################################################################
;+
; NAME:
;    rad_transfer
;
; PURPOSE:
;    Add (very simplistic) radiative transfer functionality into KinMS
;
; CALLING SEQUENCE:
;   out= rad_transfer(x2,y2,z2,los_vel,cloudsize,[intsigma=intsigma])
;
; REQUIRED INPUTS:
;   x2			A vector containing the x co-ordinate for each cloud
;   y2			A vector containing the x co-ordinate for each cloud	
;   z2			A vector containing the x co-ordinate for each cloud
;   los_vel 	A vector containing the line of sight velocity for 
;				each cloud in km/s
;   cloudsize	The assumed size for a typical cloud. All clouds within
;               cloudsize of one another with a velocity difference  
;				smaller than intsigma block radiation from one another
;
; OPTIONAL INPUTS:
;  intsigma		A scalar value for the internal velocity dispersion of 
;				each cloud (km/s). All clouds within cloudsize of one another
;               with a velocity difference smaller than intsigma block 
;				radiation from one another. [Default=3 km/s]
;
; OUTPUT:
; This function returns 1 for each cloud that is unobscured, and 0 for every cloud whose light is blocked.
;
; RESTRICTIONS:
;    Requires IDL 5.0 or higher (square bracket array
;    syntax). Requires the latest IDL astrolib. 
;    Requires J.D Smiths hist_nd procedure (included in this package)
;
; USAGE EXAMPLE:
;    As set up currently can be run simply by supplying the /radtransfer flag to KinMS
;    
;##############################################################################

function rad_transfer, x2,y2,z2,los_vel,cloudsize,intsigma=intsigma

;;; setup defaults ;;;
if not keyword_set(intsigma) then intsigma=4. ;; The value given by Larsons Laws for a 30pc cloud
flux=fltarr(n_elements(x2))+1.0 ;;; create the output flux vector initially assuming the clouds are opticially thin.
;;;

;;; work out which areas of the source have potentially overlapping clouds ;;;
hist=HIST_ND(transpose([[x2],[y2],[los_vel]]),[cloudsize*2.,cloudsize*2.,intsigma*2],REVERSE_INDICES=ri) ;; use a histogram with cells twice as large, to ensure all cases caught
s=size(hist)
lookat=where(hist gt 1) ;; these cells need checking
;;;

;;; loop over cells that need checking ;;;
for i=0,n_elements(lookat)-1 do begin
	
	;;; find the clouds that lie in this cell ;;;
	index=ARRAY_INDICES(hist, lookat[i])
    ind=[index[0]+s[1]*(index[1]+s[2]*index[2])]
    winbin=ri[ri[ind]:ri[ind+1]-1]
	;;;

	;;; for each of these clouds check if it overlaps or not, and deal with it
	for j=0,n_elements(winbin)-1 do begin
		rad=sqrt(((x2[winbin]-x2[winbin[j]])^2) + ((y2[winbin]-y2[winbin[j]])^2)) ;; spatial distance between this cloud and the others
		los_dif=abs(los_vel[winbin]-los_vel[winbin[j]]) ;; velocity difference between this cloud and the others
		w=where(rad lt cloudsize AND los_dif lt intsigma,count) ;; select clouds that truely overlap
		if count gt 1 then begin ;; if this cloud overlaps with more than itself
			order=w[sort(z2[winbin[w]])]	;; sort clouds 
			flux[winbin[order[0:-2]]]=0.0   ;; assign flux only to the nearest
		endif
	endfor
	;;;
endfor
;;;

return,flux
end