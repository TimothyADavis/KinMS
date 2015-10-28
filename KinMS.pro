;##############################################################################
;
; Copyright (C) 2012, Timothy A. Davis
; E-mail: tdavis -at- eso.org
;
; Updated versions of the software are available from my web page
; http://www.eso.org/~tdavis/
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
;    KinMS
;
; PURPOSE:
;    Simulate rotating molecular gas distributions, and produce
;    datacubes in FITS format.
;
; CALLING SEQUENCE:
;    KinMS,xs,ys,vs,dx,dy,dv,beamsize,inc,sbprof=,sbrad=,velrad=
;      ,velprof=,[posang=posang,gassigma=,noiselev=,galname=,diskthick=,cleanout=
;      ,ra=,dec=,nsamps=,cubeout=,flux_clouds=,inclouds=,vlos_clouds=]
;
; REQUIRED INPUTS:
;    XS        X-axis size for resultant cube (in arcseconds)
;    YS        Y-axis size for resultant cube (in arcseconds)
;    VS        Velocity axis size for resultant cube (in km/s)
;    DX        Pixel size in x-direction (arcsec/pixel)
;    DX        Pixel size in y-direction (arcsec/pixel)
;    DX        Channel size in velocity  (km/s/chan)
;    BEAMSIZE  Scalar or three element vector for size of convolving
;              beam (in arcseconds).  If a scalar then beam is assumed
;              to be circular. If a vector then denotes beam major
;              axis size in element zero, and the beam minor axis in
;              element one. The beam position angle should be given in
;              element two (default=zero). beamsize=[bmaj,bmin,bpa]
;    INC       Inclination angle of the gas disc on the sky
;              (degrees). Can input a constant or a vector,
;              giving the inclination as a function of the
;              radius vector SBRAD (in order to model warps etc)
;
;
;
; OPTIONAL INPUTS:
;    SBPROF    Surface brightness profile you wish to use to crate the
;              gas model. The output gas distribution will have this
;              distribution convolved with the beam. If INTFLUX is not
;              set then model will be normalised to the integral
;              of this SBPROF function over radius and galaxy position
;              angle.
;    SBRAD     Radius vector corresponding to SBPROF. Units of
;              arcseconds. If you have small structures in your SBPROF
;              make sure they are well sampled in radius.
;    VELPROF   Vector containing the velocity profile the gas follows
;              as a function of radius.
;    VELRAD    Radius vector corresponding to VELPROF. Units of
;              arcseconds.
;    INCLOUDS  If your required gas distribution is not symmetric you
;              may input vectors containing the position of the
;              clouds you wish to simulate. This 3-vector should
;              contain the X, Y and Z positions, in units of arcseconds
;              from the phase centre. If this variable is used, then
;              DISKTHICK is ignored. 
;              Example: INCLOUDS=[[0,0,0],[10,-10,2],...,[xpos,ypos,zpos]]
;    VLOS_CLOUDS This vector should contain the LOS velocity for
;                each point defined in INCLOUDS, in units of km/s. If
;                not supplied then INCLOUDS is assumed to be the -face
;                on- distribution and that VELPROF/VELRAD should be
;                used, and the distribution projected. If this
;                variable is used then GASSIGMA/INC are ignored.
;    FLUX_CLOUDS This vector can be used to supply the flux of each
;                point in INCLOUDS. If used alone then total flux in the model is equal
;                to total(FLUX_INCLOUDS). If INTFLUX used then this vector denotes
;                the relative brightness of the points in
;                INCLOUDS. Note- This option reduces execution speed
;    POSANG    Position angle of the gas distribution on the sky
;              (degrees CCW from north to the receding side of the
;              gas distribution). Default is 0 deg.  Can input a constant or
;              a vector, giving the position angle as a function of the
;              radius vector SBRAD (in order to model warps etc)
;    VPOSANG   Position angle of the line of nodes in the rotation
;              distribution of the gas (degrees CCW from north to the receding side 
;              of the gas distribution). Default is the same as POSANG.  Can input a constant or
;              a vector, giving the position angle as a function of the
;              radius vector SBRAD (in order to model warps etc)  
;    GASSIGMA  Gas velocity dispersion (in km/s). Either a constant or
;              a vector, giving the velocity dispersion as a function of the
;              radius vector SBRAD. Default is no dispersion.
;    DISKTHICK Give the molecular gas disk a non zero thickness.
;              Currently gas clouds are uniformly distributed in the
;              vertical direction. Either a constant or
;              a vector, giving the disk thickness as a function of the
;              radius vector SBRAD. Default is zero (thin disk).
;    NSAMPS    How many "GMC type" point sources to use when creating
;              the model. Default is 5e5. Lower numbers can be used to
;              simulate resolved GMCs/flocculent gas, and larger
;              numbers may be required at high angular resolutions.
;    INTFLUX   Total integrated flux you want the output gas to
;              have. (In Jy/km/s). If not set then model will be
;              normalised to 1) the integral of the SBPROF function over
;              radius and galaxy position angle, if used. 2)
;              unnormalised, if the flux vector input FLUX_CLOUDS is set
;    CLEANOUT  If set then do not convolve with the beam, and output the 
;              "clean components". Useful to create input for other
;              simulation tools (e.g sim_observe in CASA).
;    RA        RA to use in the header of the output cube (in degrees).
;    DEC       DEC to use in the header of the output cube (in degrees).
;    FILENAME  String to append to the filename for the output
;              cube. If not set then no file is written (leave out if you 
;              only wish to use the CUBEOUT option, for instance). 
;    VSYS      Systemic velocity (km/s). Default=0 km/s.
;    RESTFREQ  Rest-frequency of spectral line of choice (in Hz). Only
;              matters if you are outputting a FITS file (with
;              FILENAME keyword set). Default: 115.271e9 - 12CO(1-0)
;    PHASECEN  Offset the simulated object from the pointing centre of
;              the cube by this amount. Two element vector specifying
;              offset in ARCSECONDS in x and y. Default [0,0]
;    VPHASECEN Offset the centre of the velocity field from PHASECEN. Units of
;              ARCSECONDS. Default is no offset from PHASECEN.
;    VOFFSET   Offset the systemic velocity of the object from the
;              central velocity of the cube (set in the VSYS
;              keyword). Units of km/s. Default=0 km/s.
;    FIXSEED   Fix the seeds supplied to the random number generators
;              so that the output will always be the same for given input
;              parameters
;
;  OPTIONAL OUTPUTS:
;
;    CUBEOUT       Returns the created cube as a 3 dimensional vector
;
;
; SIDE EFFECTS:
;    The output file is overwritten if it exists.
;
; RESTRICTIONS:
;    Requires IDL 5.0 or higher (square bracket array
;    syntax). Requires the latest IDL astrolib, along with the
;    hist_nd routine (a version of which is distributed with the
;    package for ease).
;    
;
; USAGE EXAMPLE:
;    Simple usage examples are given in the KinMS_testsuite.pro file
;    
;
; MODIFICATION HISTORY:
; tdavis@eso.org
; http://www.eso.org/~tdavis/
; 05/11/2012. v1.00 - Cleaned up code and documentation ready for public release.  
; 07/03/2013. v1.10 - Added VSYS, PHASECEN AND RESTFREQ commands. Changed output fits headers to match
;                    WCS standards. 
; 05/04/2013. v1.11 - Changes to speed up the code, with thanks to
;                     Daniel Bramich at ESO.
; 30/04/2013. v1.12 - Bug fixes for INCLOUDS mechanism, thanks to John
;                     Lopez, UNSW.
; 16/10/2013. v1.20 - Changed rotation code, output normalisation,
;                     method for inputting flux
; 22/12/2013  v1.30 - Changed seeds, adding FIXSEED to allow
;                     repeatable results  
; 11/03/2014  v1.40 - Changes in implementation to allow speed
;                     increases in the convolution, and modelling of
;                     warps.
; 09/05/2014  v1.45 - Added the vradial mechanism to allow modelling
;                     of inflow and outflow
; 02/09/2015  v1.5 - Update to handling of position angles, added
;                    seperate velocity position angle/phasecentre
;
;##############################################################################


function kinms_samplefromarbdist_onesided,sbrad,sbprof,nsamps,seed,r_flat=r_flat,diskthick=diskthick

  ;;;; This function takes the input radial distribution and draws
  ;;;; n samples from under it. It also accounts for disk thickness
  
     px=dblarr(n_elements(sbprof))
     sbprofboost=sbprof*(2*!dpi*abs(sbrad))                                                        ;; boost flux in outer parts to compensate for spreading over 2pi radians
     for i=1,n_elements(sbrad)-1 do px[i] = INT_TABULATED((sbrad[0:i]),(sbprofboost[0:i]),/double) ;; integrate surf brightness w.r.t radius   
     px/=max(px)                                                                                   ;; normalize
     pick=randomu(seed[0],nsamps)                                                                  ;; pick random y to transform
     phi=randomu(seed[1],nsamps)*2*!dpi                                                            ;; random angle
     r_flat=interpol(sbrad,px,pick)                                                                ;; invert to get correct distribution
     if n_elements(diskthick) gt 1 then diskthick_here=interpol(diskthick,sbrad,r_flat) else diskthick_here=diskthick
     zpos=(diskthick_here)*((randomu(seed[2],nsamps)*2)-1)                                         ;; do disk thickness
     r_3d = sqrt((r_flat^2)+zpos^2)                                                                ;; 3d distance to each
     theta=acos(zpos/r_3d)                                                                         ;; angle out of plane
     xpos=((r_3d*cos(phi)*sin(theta)))                                                             ;; go to cartesian
     ypos=((r_3d*sin(phi)*sin(theta)))                                                             ;; go to cartesian
     return,TRANSPOSE([[xpos],[ypos],[zpos]])                                                      ;; return INCLOUDS format
end

function kinms_create_velfield_onesided,velrad,velprof,r_flat,inc,posang,gassigma,seed,xpos,ypos,vphasecent,vposang=vposang,vradial=vradial,inc_rad=inc_rad,posang_rad=posang_rad


  ;;;; This function takes the input circular velocity distribution
  ;;;; and creates the velocity field taking into account warps,
  ;;;; inflow/outflow etc as required by the keywords passed.
  
  vrad=interpol(velprof,velrad,r_flat)
  los_vel=dblarr(n_elements(vrad))

  ;;; in case of warps/inflow/outflow then setup vectors to deal with this
  if n_elements(vradial) gt 1 then vradial_rad=interpol(vradial,velrad,r_flat) else vradial_rad=fltarr(n_elements(r_flat))+vradial
  if n_elements(vposang) gt 1 then vposang_rad=interpol(vposang,velrad,r_flat) else vposang_rad=vposang
  ;;;;

   ;;;;; include random motions
  if keyword_set(gassigma) then begin
     veldisp=randomn(seed[3],n_elements(xpos)) 
     if n_elements(gassigma) gt 1 then veldisp=temporary(veldisp)*interpol(gassigma,velrad,r_flat) else veldisp=temporary(veldisp)*gassigma
  endif else  veldisp=fltarr(n_elements(xpos))
  ;;;;;

  ;;;;; project velocities taking into account inclination, which may change with radius
  w=where(xpos+vphasecent[0] ne 0.0 and ypos+vphasecent[1] ne 0.0)
  ang2rot=(90-(posang_rad-vposang_rad))
  los_vel[w]=veldisp[w]                                                                                                                              ;; random motions
  los_vel[w]+=vrad[w]*(cos(atan(float(ypos[w]+vphasecent[1])/float(xpos[w]+vphasecent[0]))+(!dtor*(90-ang2rot[w])))*sin(inc_rad[w]*!dtor))           ;; circular motions
  los_vel[w]+=vradial_rad[w]*(sin(atan(float(ypos[w]+vphasecent[1])/float(xpos[w]+vphasecent[0]))+(!dtor*(90-ang2rot[w])))*sin(inc_rad[w]*!dtor))    ;; inflow/outflow
  w=where(xpos+vphasecent[0] eq 0.0 and ypos+vphasecent[1] eq 0.0)
  if w[0] ne -1 then los_vel[w]=veldisp[w]+(vrad[w]*sin(inc_rad[w]*!dtor))+ vradial_rad[w]*(sin(atan(float(ypos[w])/float(xpos[w]))+(!dtor*(90-ang2rot[w])))*sin(inc_rad[w]*!dtor))
  w=where(xpos+vphasecent[0] gt 0.0)
  los_vel[w]=temporary(los_vel[w])*(-1)
  ;;;;;;
  
  return,los_vel
  end


pro KinMS,xs,ys,vs,dx,dy,dv,beamsize,inc,gassigma=gassigma,sbprof=sbprof,sbrad=sbrad,velrad=velrad,velprof=velprof,filename=galname,diskthick=diskthick,cleanout=cleanout,ra=ra,dec=dec,nsamps=nsamps,cubeout=cubeout,posang=posang,intflux=intflux,inclouds=inclouds,vlos_clouds=vlos_clouds,flux_clouds=flux_clouds,vsys=vsys,restfreq=restfreq,phasecen=phasecen,voffset=voffset,fixseed=fixseed,vradial=vradial,vphasecen=vphasecen,vposang=vposang,sb_sampfunc=sb_sampfunc,vel_func=vel_func
;!EXCEPT = 2


;;;; set defaults ;;;;
  if not keyword_set(nsamps) then nsamps=10000
  if size(beamsize,/N_DIMENSIONS) eq 0 then beamsize=[beamsize,beamsize,0.0] else if n_elements(beamsize) eq 2 then beamsize=[beamsize[0],beamsize[1],0.0]
  if not keyword_set(diskthick) then diskthick=0.
  if not keyword_set(gassigma) then gassigma=0.
  if not keyword_set(posang) then posang=0.
  if not keyword_set(vposang) then vposang=posang
  if not keyword_set(vradial) then vradial=0.
  if not keyword_set(nsamps) then nsamps=10000
  if not keyword_set(vsys) then vsys=0
  if not keyword_set(voffset) then voffset=0
  if not keyword_set(ra) then ra=0
  if not keyword_set(dec) then dec=0
  if not keyword_set(restfreq) then restfreq=115.271e9 ;; 12CO(1-0)
  if not keyword_set(velprof) and not keyword_set(INCLOUDS) then stop,"Please define either a SB profile/vel curve or individual cloud positions with the SBPROF/INCLOUDS mechanisms."
  if not keyword_set(velrad) and not keyword_set(INCLOUDS) then stop,"Please define either a SB profile/vel curve or individual cloud positions with the SBPROF/INCLOUDS mechanisms."
  if not keyword_set(sbprof) and not keyword_set(INCLOUDS) then stop,"Please define either a SB profile/vel curve or individual cloud positions with the SBPROF/INCLOUDS mechanisms."
  if not keyword_set(sbrad) and not keyword_set(INCLOUDS) then stop,"Please define either a SB profile/vel curve or individual cloud positions with the SBPROF/INCLOUDS mechanisms."
  if n_elements(gassigma) gt 1 and n_elements(gassigma) ne n_elements(sbrad) then stop,"If gassigma is a vector must be the same size as SBRAD"
  if n_elements(diskthick) gt 1 and n_elements(diskthick) ne n_elements(sbrad) then stop,"If diskthick is a vector must be the same size as SBRAD"
  if keyword_set(inclouds) then if size(inclouds,/N_DIMENSIONS) ne 2 then stop,"INCLOUDS must have 3 dimensions"
  if keyword_set(inclouds) then if not keyword_set(vlos_clouds) then if not keyword_set(velprof) then stop,"If using INCLOUDS and not using VLOS_CLOUDS then VELPROF/VELRAD must be set"
  if keyword_set(vlos_clouds) then if size(vlos_clouds,/N_DIMENSIONS) ne 1 then stop,"VLOS_CLOUDS must be a vector"
  if keyword_set(phasecen) then if (size(phasecen,/N_DIMENSIONS) ne 1 OR n_elements(phasecen) ne 2) then stop,"PHASECEN must be a 2 element vector with X and Y offsets in arcseconds."
  if keyword_set(vphasecen) then if (size(vphasecen,/N_DIMENSIONS) ne 1 OR n_elements(vphasecen) ne 2) then stop,"VPHASECEN must be a 2 element vector with X and Y offsets in arcseconds."
  if not keyword_set(ra) then ra=0.
  if not keyword_set(dec) then dec=0.
  if not keyword_set(phasecen) then phasecen=[0,0]
  if not keyword_set(vphasecen) then vphasecen=phasecen
  if keyword_set(fixseed) then begin
     seed=indgen(4)+100 ;; fix seeds 
  endif else seed=randomu(seedit,4)*100. ;; random seeds
  if not keyword_set(sb_sampfunc) then sb_sampfunc="kinms_samplefromarbdist_onesided"
  if not keyword_set(vel_func) then vel_func="kinms_create_velfield_onesided"
;;;;
  
;;;; work out images sizes ;;;;
  xsize = round(xs/dx)
  ysize = round(ys/dy)
  vsize = round(vs/dv)
  cent=[(xsize/2.)+(phasecen[0]/dx),(ysize/2.)+(phasecen[1]/dy),(vsize/2.)+(voffset/dv)]
  velcube=(findgen(vsize)-cent[2])*dv
  vphasecent=(vphasecen-phasecen)/[dx,dy]
;;;;

  if not keyword_set(inclouds) then begin
     inclouds=call_FUNCTION(sb_sampfunc,sbrad,sbprof,nsamps,seed,r_flat=r_flat,diskthick=diskthick) ;; use the function in sb_sampfunc to setup inclouds, if its not already set
     r_flat/=dx
  endif
  
;;;; set up INCLOUDS for use
  xpos=reform(INCLOUDS[0,*]/dx)
  ypos=reform(INCLOUDS[1,*]/dy)
  zpos=reform(INCLOUDS[2,*]/dx)
  if n_elements(r_flat) eq 0 then r_flat=sqrt(((xpos-(phasecen[0]/dx))^2) + ((ypos-(phasecen[1]/dy))^2))
  if keyword_set(VLOS_CLOUDS) then los_vel=VLOS_CLOUDS
  



;;;; create velocity structure ;;;;
  if not keyword_set(VLOS_CLOUDS) then begin

     if n_elements(inc) gt 1 then inc_rad=interpol(inc,velrad,r_flat) else inc_rad=fltarr(n_elements(r_flat))+inc
     if n_elements(posang) gt 1 then posang_rad=interpol(posang,velrad,r_flat) else posang_rad=posang
     
     los_vel=call_FUNCTION(vel_func,velrad/dx,velprof,r_flat,inc,posang,gassigma,seed,xpos,ypos,vphasecent,vposang=vposang,vradial=vradial,inc_rad=inc_rad,posang_rad=posang_rad)
     ;;;; project face on clouds to desired inclination ;;;;
     c = cos(inc_rad/!radeg)
     s = sin(inc_rad/!radeg)
     x2 =  xpos
     y2 =  c*ypos + s*zpos
     z2 = -s*ypos + c*zpos       
     if n_elements(posang) ge 1 then begin  
        ;; rotate gas on sky to required angle
        ang=90-posang_rad
        c = cos(ang/!radeg)
        s = sin(ang/!radeg)
        x3 =  c*x2 + s*y2
        y3 = -s*x2 + c*y2
        x2=x3
        y2=y3
     endif
  endif else begin
     x2=xpos
     y2=ypos
     z2=zpos
  endelse
  if not keyword_set(vlos_clouds) then  vlos_clouds=los_vel

  
  
;;;; now add the flux into the cube
  los_vel_dv_cent2=round((los_vel/dv)+cent[2])
  x2_cent0=round(x2+cent[0])
  y2_cent1=round(y2+cent[1])
  ;Get good subs
  subs = where(((x2_cent0 ge 0L) AND (x2_cent0 lt xsize) AND (y2_cent1 ge 0L) AND (y2_cent1 lt ysize) AND (los_vel_dv_cent2 ge 0L) AND (los_vel_dv_cent2 lt vsize)) EQ 1B, nsubs)

 ;;;For robustness should really test nsubs GT 0 before continuing
  if nsubs gt 0L then begin
     if not keyword_set(FLUX_CLOUDS) then begin
       cube=double(hist_nd(transpose([[x2_cent0[subs]],[y2_cent1[subs]],[los_vel_dv_cent2[subs]]]),1,min=[0,0,0],max=[xsize,ysize,vsize]-1))
     endif else begin
        cube = dblarr(xsize,ysize,vsize)
        ;Get subs in 3D array (repeats are ok)
        subs3d = x2_cent0[subs] + (y2_cent1[subs]*xsize) + (los_vel_dv_cent2[subs]*(xsize*ysize))
        for i = 0L,(nsubs-1L) do begin
           const = FLUX_CLOUDS[subs[i]]
           csub = subs3d[i]
           cube[csub] = cube[csub] + const
        endfor
     endelse
  end


;;;;


;;;; do beamsize convolution ;;;;
  if not keyword_set(cleanout) then begin 
      psf=rot(psf_gaussian(Npixel=[xsize,ysize],fwhm=[beamsize[0]/dx,beamsize[1]/dy],NDIMEN=2,/double) >1e-7,90-beamsize[2],/interp,missing=0.0)
      doconvolve=where(total(total(cube,1),1) gt 0.0) 
     ;; only convolve the planes with flux in them for speed
     for nplane=0,n_elements(doconvolve)-1 do begin
        ;; if you want a convolved cube, then convolve it
          cube[*,*,doconvolve[nplane]] = convol_fft(cube[*,*,doconvolve[nplane]],psf,KERNEL_FFT=psf_FT)
       endfor
  endif
;;;;
;
;;; normalize fluxes, if requested
  if keyword_set(intflux) then begin
      if not keyword_set(cleanout) then cube *= ((intflux*total(psf))/((total(cube)*dv))) else cube *= ((intflux)/((total(cube)*dv))) ;; make total flux equal to influx
  endif else begin
     if keyword_set(flux_clouds) then begin
        cube*=(total(flux_clouds)/total(cube)) ;; make total flux equal to total(flux_clouds)
     endif else begin
        cube*=((INT_TABULATED( sbrad, sbprof))/((total(cube)*dv)/(sqrt(beamsize[0]*beamsize[1])^2))) ;; normalize to get same flux as input sb profile variable
     endelse
  endelse

;;;; write cube ;;;;

  if keyword_set(galname) then begin
     mkhdr, h, cube
     sxdelpar,h,'DATE'
     sxdelpar,h,'COMMENT'
     sxaddpar,h,'NAXIS4',1
     sxaddpar,h,'NAXIS',4
     sxaddpar,h,'CDELT1',(dx)/(-3600d0)
     sxaddpar,h,'CDELT2',(dy)/3600d0
     sxaddpar,h,'CDELT3',(dv)*1000d,"m/s"
     sxaddpar,h,'CDELT4',1d
     sxaddpar,h,'CRPIX1',double(cent[0]-1)
     sxaddpar,h,'CRPIX2',double(cent[1])
     sxaddpar,h,'CRPIX3',double(cent[2])
     sxaddpar,h,'CRPIX4',1d
     sxaddpar,h,'CRVAL1',double(ra)
     sxaddpar,h,'CRVAL2',double(dec)
     sxaddpar,h,'CRVAL3',double(vsys*1000d),"m/s"
     sxaddpar,h,'CRVAL4',1d
     sxaddpar,h,'CUNIT1','deg'
     sxaddpar,h,'CUNIT2','deg'
     sxaddpar,h,'CUNIT3','m/s'
     sxaddpar,h,'CUNIT4',''
     sxaddpar,h,'BSCALE',1d,'PHYSICAL = PIXEL*BSCALE + BZERO'              
     sxaddpar,h,'BZERO',0d                                                                                     
     sxaddpar,h,'BMIN',min(beamsize[0:1]/3600d0)
     sxaddpar,h,'BMAJ',max(beamsize[0:1]/3600d0)
     sxaddpar,h,'BTYPE','Intensity'  
     sxaddpar,h,'BPA',beamsize[2]
     sxaddpar,h,'CTYPE1','RA---SIN' 
     sxaddpar,h,'CTYPE2','DEC--SIN'
     sxaddpar,h,'CTYPE3','VRAD'
     sxaddpar,h,'CTYPE4','STOKES'  
     sxaddpar,h,'EQUINOX',2000d
     sxaddpar,h,'RADESYS','FK5'
     sxaddpar,h,'BUNIT','Jy/beam'
     sxaddpar,h,'RESTFRQ',double(restfreq),'[Hz]'
     sxaddpar,h,'SPECSYS','BARYCENT'
     writefits,galname+"_simcube.fits",cube,h
  endif
  
  cubeout=cube
  
end

