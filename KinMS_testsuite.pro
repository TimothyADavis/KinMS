pro KinMStest_ngc4324,_extra=_extra
;;;;;;;;;;;
;
; A test procedure to demonstrate the KinMS code, and check if it
; works on your system. This procedure makes a basic simulation of the
; molecular gas ring in NGC4324, and plots the simulation moment zero,
; one and PVD against the observed ones from the CARMA observations
; of Alatalo et al., 2012. Any default paramaters can be changed by 
; specifying them at the command line (see KinMS.pro for the full
; details of all the avaliable parameters).
;
;;;;;;;;;;;

;;; Define the simulated observation parameters ;;;
   xsize=128 ;; arcseconds
   ysize=128 ;; arcseconds
   vsize=800 ;; km/s
   cellsize=1 ;; arcseconds/pixel
   dv=10 ;; km/s/channel
   beamsize=4. ;; arcseconds
;;;

;;; Define the gas distribution required ;;;
   diskthick=1. ; arcseconds
   inc=65. ; degrees
   posang=230. ; degrees
   x=(findgen(64)*10.)/10.
   x1=x
   fx = gaussian(x,[0.1,21,3])
   vel=interpol([0,50,100,175,175,175,175],[0.01,1,3,5,7,10,200],x)
   phasecen=[-1,-1]
   voffset=-10
;;;

;;; Run KinMS
   KinMS,xsize,ysize,vsize,cellsize,dv,beamsize,inc,sbprof=fx,sbrad=x,velrad=x1,velprof=vel,diskthick=dickthick,nsamps=nsamps,cubeout=f,posang=posang,intflux=27.2,_extra=_extra,phasecen=phasecen,voffset=voffset
;;;

;;; Create rough moment zero and one maps and a major axis PVD from
;;; the simulation 
   mom0rot=total(f,3)
    x1=(findgen((xsize/cellsize))-((xsize/cellsize)/2.))*cellsize
   v1=(findgen(vsize/dv)-((vsize/dv)/2.))*dv
   mom1=mom0rot*0.0
   for i=0,127 do for j=0,127 do mom1[i,j]=mean(v1[where(f[i,j,*] eq max(f[i,j,*]))])
   mom1[where(mom0rot lt 0.3*max(mom0rot))] =-10000
   pvdcube=f*0.0
   for i=0,79 do pvdcube[*,*,i]=rot(f[*,*,i],-230+180,/interp)
   pvd=reform(pvdcube[*,64,*])

;;;
   
;;; Get the observed data from the IDL save file supplied with the
;;; test suite. 
   restore,"test_suite/KinMStest_NGC4324.sav"
;;;
   

;;; Create the plot to your x-device
   device,decomposed=0,RETAIN=2
   loadct,39,/silent
   !p.multi=[0,2,2]
   contour,mom0rot,x1,x1,levels=max(mom0rot)*((findgen(5)+5)/10.),xrange=[-25,25],yrange=[-25,25],xtitle="RA (arcsec)",ytitle="DEC (arcsec)",/xstyle,/ystyle,background=255,color=0
   contour,mom0gal,x2,x2,/overplot,color=250,nlevels=10,levels=max(mom0gal,/nan)*((findgen(8)+2)/10.)
   al_legend,["Moment Zero"],/top,/right,box=0,textcolor=0
   al_legend,["Data","Model"],/bottom,/left,box=0,textcolor=0,color=[250,0],linestyle=[0,0] 
   contour,pvd,x1,v1,levels=max(pvd)*((findgen(9)+3)/10.),xrange=[-30,30],yrange=[-200,200],/xstyle,/ystyle,xtitle="Position (arcsec)",ytitle="Velocity (km/s)",background=255,color=0
   contour,pvdgal,x3,v3-5,levels=max(pvdgal)*((findgen(9)+3)/10.),/overplot,color=250
   al_legend,["PVD (PA=50 deg)"],/top,/left,box=0,textcolor=0
   contour,mom1gal-1662.,x2,x2,/fill,levels=v3[1:43],xrange=[-25,25],yrange=[-25,25],/cell,color=0,xtitle="RA (arcsec)",ytitle="DEC (arcsec)",xmargin=[15,0],/xstyle,/ystyle
   al_legend,["Moment One (data)"],/top,/right,box=0,textcolor=0
   contour,mom1,x1,x1,/fill,levels=v3[1:43],xrange=[-25,25],yrange=[-25,25],/cell,color=0,xtitle="RA (arcsec)",ytickname=replicate(' ',10),xmargin=[0,15],/xstyle,/ystyle
   al_legend,["Moment One (model)"],/top,/right,box=0,textcolor=0
  ; COLORBAR,/vert,/right,range=[min(v3[1:43],/nan),max(v3[1:43],/nan)],DIVISIONS=6,title="Velocity (km/s)",position=[!x.window[1]+0.01,!y.window[0],!x.window[1]+0.05,!y.window[1]],ncolors=255,charsize=1,color=0



   !p.multi=[0,1,1]
;;;
end
pro KinMStest_expdisk,scalerad,inc,_extra=_extra
;;;;;;;;;;;
;
; A test procedure to demonstrate the KinMS code, and check if it
; works on your system. This procedure demonstrates how to create a simulation of an
; exponential disk of molecular gas. The user should input values for the
; scalerad and inc variables, and the procedure will the create the simulation
; and display it to screen. Any default paramaters can be changed by specifying 
; them at the command line (see KinMS.pro for the full
; details of all the avaliable parameters).
;
;  INPUTS:
;       Scalerad -    Scale radius for the exponential disk (arcseconds) 
;       Inc      -    Inclination to project the disk (degrees).
;
;;;;;;;;;;;

;;;; Validate input ;;;;
Fcent=1.
  if keyword_set(Fcent) + keyword_set(scalerad) ne 2 then stop,"Please pass as the first argument the scale radius for the molecular disk (and the inclination)" 
  if not keyword_set(inc) then begin
     inc=0.0
     print,"Inc not set: will assume 0 degrees (face on)"
  endif

;;;; Setup cube parameters ;;;;
   xsize=128
   ysize=128
   vsize=1400
   cellsize=1
   dv=10
   beamsize=4.
   nsamps=5e5
;;;;

;;;; Set up exponential disk SB profile/velocity ;;;;
x=findgen(1000)/10.
x1=x
fx = Fcent*exp(-x/scalerad)
vel=interpol([0,50,100,210,210],[0.01,0.5,1,3,500],x)
;;;;


;;;; Simulate ;;;;
KinMS,xsize,ysize,vsize,cellsize,dv,beamsize,inc,sbprof=fx,sbrad=x,velrad=x1,velprof=vel,cubeout=f,nsamps=nsamps,_extra=_extra,intflux=30.,posang=270
;;;;


;;;; Create plot data from cube ;;;;
   mom0rot=total(f,3)
   x1=(findgen((xsize/cellsize))-((xsize/cellsize)/2.))*cellsize
   v1=(findgen(vsize/dv)-((vsize/dv)/2.))*dv
   mom1=(mom0rot*0.0)-10000.0
   for i=0,(xsize/cellsize)-1 do for j=0,(ysize/cellsize)-1 do begin
      if mom0rot[i,j] gt 0.1*max(mom0rot) then begin
         mom1[i,j]=mean(v1[where(f[i,j,*] eq max(f[i,j,*]))])
      endif
   endfor
   pvdcube=f
   if keyword_set(_extra) then begin
      if TAG_EXIST(_extra, posang,/quiet) then for i=0,n_elements(v1)-1 do pvdcube[*,*,i]=rot(f[*,*,i],_extra.posang-90,/interp)
   endif
   pvd=reform(pvdcube[*,64,*])
   spec=(total(total(f,1),1))/(mean(beamsize)^2)
;;;;

;;;; Plot the results ;;;;
   device,decomposed=0,RETAIN=2
   loadct,39,/silent
   levs=v1[where(abs(v1) le max(vel*sin(!dtor*inc))+20.)]
   !p.multi=[0,2,2]
   contour,mom0rot,x1,x1,/cell_fill,levels=(max(mom0rot)*((findgen(9)+1)/10.)),color=0,background=255,/xstyle,/ystyle,xtitle="RA (arcsec)",ytitle="DEC (arcsec)",/iso
   al_legend,["Moment 0"],box=0,/top,/right,textcolor=0
   contour,mom1,x1,x1,/cell_fill,levels=levs,color=0,background=255,/xstyle,/ystyle,xtitle="RA (arcsec)",ytitle="DEC (arcsec)" ,/iso
   al_legend,["Moment 1"],box=0,/top,/right,textcolor=0
   contour,pvd,x1,v1,levels=max(pvd)*((findgen(9)+1)/10.),/xstyle,/ystyle,xtitle="Position (arcsec)",ytitle="Velocity (km/s)",background=255,color=0,yrange=[min(levs)-30,max(levs)+30]
   al_legend,["Major Axis PVD"],box=0,/top,/right,textcolor=0
   plot,v1,spec,color=0,psym=10,xtitle="Velocity (km/s)",ytitle="Flux",/xstyle,/ystyle
   al_legend,["Spectrum"],box=0,/top,/right,textcolor=0
   !p.multi=[0,1,1]
;;;;
end
pro KinMStest_inclouds,_extra=_extra
;;;;;;;;;;;
;
; A test procedure to demonstrate the KinMS code, and check if it
; works on your system. This procedure demonstrates how to use the
; INCLOUDS parameter set to create simulations, in this case of a very
; unrealistic object. Once you understand this example then see the
; INFITS and INCLOUDS_SPIRAL test for more realistic examples.
;
;;;;;;;;;;;

;;;; Setup cube parameters ;;;;
   xsize=128
   ysize=128
   vsize=1400
   cellsize=1
   dv=10
   beamsize=4.
   inc=0.
;;;;

;;;; define where clouds are in each dimension (x,y,z) ;;;;
inclouds=transpose([[ 40.0000, 0.00000, 0.00000],[ 39.5075, 6.25738, 0.00000],[ 38.0423, 12.3607, 0.00000],[ 35.6403, 18.1596, 0.00000],[ 32.3607, 23.5114, 0.00000],[ 28.2843, 28.2843, 0.00000],[ 23.5114, 32.3607, 0.00000],[ 18.1596, 35.6403, 0.00000],[ 12.3607, 38.0423, 0.00000],[ 6.25737, 39.5075, 0.00000],[ 0.00000, 40.0000, 0.00000],[-6.25738, 39.5075, 0.00000],[-12.3607, 38.0423, 0.00000],[-18.1596, 35.6403, 0.00000],[-23.5114, 32.3607, 0.00000],[-28.2843, 28.2843, 0.00000],[-32.3607, 23.5114, 0.00000],[-35.6403, 18.1596, 0.00000],[-38.0423, 12.3607, 0.00000],[-39.5075, 6.25738, 0.00000],[-40.0000, 0.00000, 0.00000],[-39.5075,-6.25738, 0.00000],[-38.0423,-12.3607, 0.00000],[-35.6403,-18.1596, 0.00000],[-32.3607,-23.5114, 0.00000],[-28.2843,-28.2843, 0.00000],[-23.5114,-32.3607, 0.00000],[-18.1596,-35.6403, 0.00000],[-12.3607,-38.0423, 0.00000],[-6.25738,-39.5075, 0.00000],[ 0.00000,-40.0000, 0.00000],[ 6.25738,-39.5075, 0.00000],[ 12.3607,-38.0423, 0.00000],[ 18.1596,-35.6403, 0.00000],[ 23.5114,-32.3607, 0.00000],[ 28.2843,-28.2843, 0.00000],[ 32.3607,-23.5114, 0.00000],[ 35.6403,-18.1596, 0.00000],[ 38.0423,-12.3607, 0.00000],[ 39.5075,-6.25737, 0.00000],[ 15.0000, 15.0000, 0.00000],[-15.0000, 15.0000, 0.00000],[-19.8504,-2.44189, 0.00000],[-18.0194,-8.67768, 0.00000],[-14.2856,-13.9972, 0.00000],[-9.04344,-17.8386, 0.00000],[-2.84630,-19.7964, 0.00000],[ 3.65139,-19.6639, 0.00000],[ 9.76353,-17.4549, 0.00000],[ 14.8447,-13.4028, 0.00000],[ 18.3583,-7.93546, 0.00000],[ 19.9335,-1.63019, 0.00000]])
;;;;

vlos_clouds=FLTARR(52) ;; in this unrealistic case each cloud is at the systemic velocity


;;;; run the simulation with VLOS_CLOUDS ;;;;
   KinMS,xsize,ysize,vsize,cellsize,dv,beamsize,inc,intflux=30.,inclouds=inclouds,vlos_clouds=vlos_clouds,posang=90.,cubeout=f,_extra=_extra,filename="Testsmile"
;;;;
;;;; run the simulation with a velocity curve ;;;;
   x=findgen(1000)/10.
   vel=interpol([0,50,100,210,210],[0.01,0.5,1,3,500],x)
   inc=55.
   KinMS,xsize,ysize,vsize,cellsize,dv,beamsize,inc,intflux=30.,inclouds=inclouds,velprof=vel,velrad=x,cubeout=f2,posang=90.,_extra=_extra
;;;;

;;;; Plot the results ;;;;
   device,decomposed=0,RETAIN=2
   loadct,0,/silent
   mom0=total(f,3)
   mom0_2=total(f2,3)
   v1=(findgen(vsize/dv)-((vsize/dv)/2.))*dv
   levs=v1[where(abs(v1) le max(vel*sin(!dtor*inc))+20.)]
   mom1=mom0*0.0
   s=size(mom1)
   for i=0,s[1]-1 do for j=0,s[2]-1 do mom1[i,j]=mean(v1[where(f2[i,j,*] eq max(f2[i,j,*]))])
   mom1[where(mom0_2 lt 0.1*max(mom0_2))] =-10000
   !p.multi=[0,2,1]
   contour,mom0,/fill,levels=(max(mom0,/nan)*((findgen(9)+1)/10.)),/xstyle,/ystyle,/iso,background=255,color=0,xtitle="RA (arcsec)",ytitle="DEC (arcsec)"
   loadct,39,/silent
   contour,mom1,/fill,levels=levs,/xstyle,/ystyle,/iso,background=255,color=0,xtitle="RA (arcsec)",ytitle="DEC (arcsec)"
   !p.multi=[0,1,1]
   
;;;;

end

pro KinMStest_inclouds_spiral,_extra=_extra
;;;;;;;;;;;
;
; A test procedure to demonstrate the KinMS code, and check if it
; works on your system. This procedure demonstrates how to use the
; INCLOUDS parameter set to create simulations, in this case of
; molecular gas in a two armed spiral pattern. Any default paramaters
; can be changed by specifying them at the command line (see KinMS.pro 
; for the full details of all the avaliable parameters).
;
;
;;;;;;;;;;;

;;;; Setup cube parameters ;;;;
   xsize=128*2
   ysize=128*2
   vsize=1400
   cellsize=1
   dv=10
   beamsize=4.
   inc=55.
;;;;

;;;; define where clouds are in each dimension (x,y,z) using a logarithmic spiral ;;;;

t=((findgen(4000))/100.)-20.
a=0.002*60
b=0.5
x1=a*exp(b*t)*cos(t)
y1=a*exp(b*t)*sin(t)
x2=a*(-1)*exp(b*t)*cos(t)
y2=a*(-1)*exp(b*t)*sin(t)
inclouds=[[x1,x2],[y1,y2],[y1*0.0,y2*0.0]]

;;;;


;;;; define velocity curve ;;;;
   x=findgen(1000)/10.
   vel=interpol([0,50,100,210,210],[0.01,0.5,1,3,500],x)
;;;;

;;;; run the simulation ;;;;
   KinMS,xsize,ysize,vsize,cellsize,dv,beamsize,inc,intflux=30.,inclouds=inclouds,velprof=vel,velrad=x,cubeout=f2,flux_clouds=fluxclouds,_extra=_extra,filename="spiraltest"
;;;;

;;;; Plot the results ;;;;
   device,decomposed=0,RETAIN=2
   loadct,0,/silent
   mom0_2=total(f2,3,/nan)
   v1=(findgen(vsize/dv)-((vsize/dv)/2.))*dv
   levs=v1[where(abs(v1) le max(vel*sin(!dtor*inc))+20.)]
   mom1=mom0_2*0.0
   s=size(mom1)
   for i=0,s[1]-1 do for j=0,s[2]-1 do mom1[i,j]=mean(v1[where(f2[i,j,*] eq max(f2[i,j,*]))])
   mom0_2[where(mom0_2 lt 0.0005)] =0
   mom1[where(mom0_2 eq 0)]=-100000
   mom0_2=alog10(mom0_2)
   mom0_2-=min(mom0_2,/nan)
   !p.multi=[0,2,1]
   contour,mom0_2,/fill,levels=(max(mom0_2,/nan)*((findgen(59)+1)/60.)),/xstyle,/ystyle,/iso,background=255,color=0,xtitle="RA (arcsec)",ytitle="DEC (arcsec)"
   loadct,39,/silent
   contour,mom1,/fill,levels=levs,/xstyle,/ystyle,/iso,background=255,color=0,xtitle="RA (arcsec)",ytitle="DEC (arcsec)"
   !p.multi=[0,1,1]
;;;;

end



pro KinMStest_infits,_extra=_extra
;;;;;;;;;;;
;
; A test procedure to demonstrate the KinMS code, and check if it
; works on your system. This procedure demonstrates how to use the
; an input FITS image to create a simulation of what the molecular gas
; may look like, with a given instrument (in this case CARMA). We use a
; GALEX (Morrissey et al., 2007) FUV image of NGC1437A, and scale it
; assuming the FUV emission comes from star-formation and thus
; molecular gas, and that the galaxy has a total integrated CO flux of
; 30 Jy km/s. We use the FITS image to set the surface-brightness, and
; impose a flat velocity gradient.
;
;;;;;;;;;;;

;;;; Setup cube parameters ;;;;
   xsize=128
   ysize=128
   vsize=1400
   cellsize=1
   dv=10
   beamsize=5.
   inc=0.
;;;;


;;;; Read in the FITS file and create the INCLOUDS variables based on it ;;;;
   phasecent=[88,61] ;; point we wish to correspond to the phase centre in the simulation
   fin=readfits("test_suite/NGC1437A_FUV.fits",hdr)
   s=size(fin)
   xvec=(findgen(s[1])-phasecent[0])*(sxpar(hdr,'cdelt1')*3600.)
   yvec=(findgen(s[2])-phasecent[1])*(sxpar(hdr,'cdelt2')*3600.)
   w=where(fin gt 0.002) ;; clip the image to avoid FUV noise entering the simulation
   flux_clouds=fin[w]
   index=array_indices(fin,w)
   x=reform(xvec[index[0,*]])
   y=reform(yvec[index[1,*]])
   inclouds=[[x],[y],[y*0.0]]
   ang=80*!dtor
   vlos_clouds=interpol([-200,0,200],[-60,0,60],y*sin(ang)+x*cos(ang)) ;; impose a flat velocity profile
;;;;

;;;; run the simulation ;;;;
   KinMS,xsize,ysize,vsize,cellsize,dv,beamsize,inc,intflux=30.,inclouds=inclouds,vlos_clouds=vlos_clouds,cubeout=f,posang=90.,_extra=_extra,flux_cloud=flux_clouds
;;;;

;;;; Display the output ;;;;
   device,decomposed=0,RETAIN=2
   loadct,39,/silent
   !p.multi=[0,2,2]
   x1=(findgen((xsize/cellsize))-((xsize/cellsize)/2.))*cellsize
   v1=(findgen(vsize/dv)-((vsize/dv)/2.))*dv
   contour,fin,xvec,yvec,/fill,levels=(max(fin)*((findgen(19)+1)/20.)),/xstyle,/ystyle,xrange=[min(x1),max(x1)],yrange=[min(x1),max(x1)],xtitle="RA (arcsec)",ytitle="DEC (arcsec)",color=0,background=255,/iso
   al_legend,["Input"],/top,/right,box=0,textcolor=0
   mom0=total(f,3)
   contour,mom0,x1,x1,/fill,levels=(max(mom0)*((findgen(19)+1)/20.)),/xstyle,/ystyle,xtitle="RA (arcsec)",ytitle="DEC (arcsec)",color=0,background=255,/iso
   al_legend,["Moment 0"],/top,/right,box=0,textcolor=0
   mom1=mom0*0.0
   for i=0,xsize-1 do for j=0,ysize-1 do mom1[i,j]=mean(v1[where(f[i,j,*] eq max(f[i,j,*]))])
   mom1[where(mom0 lt 0.05*max(mom0))] =-10000
   contour,mom1,x1,x1,/cell_fill,levels=v1[70-10:70+10],/xstyle,/ystyle,xtitle="RA (arcsec)",ytitle="DEC (arcsec)",color=0,background=255,/iso
   al_legend,["Moment 1"],/top,/right,box=0,textcolor=0
   plot,v1,total(total(f,1),1),psym=10,xtitle="Velocity (km/s)",ytitle="Flux",xrange=[-300,300],color=0,background=255
   al_legend,["Spectrum"],/top,/right,box=0,textcolor=0
   !p.multi=[0,1,1]
;;;;
end
pro KinMStest_veldisp,_extra=_extra
;;;;;;;;;;;
;
; A test procedure to demonstrate the KinMS code, and check if it
; works on your system. This procedure demonstrates how to create a simulation of an
; exponential disk of molecular gas with a velocity dispersion that
; varies with radius. Any default paramaters
; can be changed by specifying them at the command line (see KinMS.pro 
; for the full details of all the avaliable parameters).
;
;
;;;;;;;;;;;



;;;; Setup cube parameters ;;;;
   xsize=128
   ysize=128
   vsize=1400
   cellsize=1
   dv=10
   beamsize=2.
   nsamps=5e5
;;;;

;;;; Set up exponential disk SB profile/velocity ;;;;
fcent=10.
scalerad=20.
inc=30.
x=findgen(1000)/10.
x1=x
fx = Fcent*exp(-x/scalerad)
gassigma=interpol([50,8,8],[0,20,500],x)
vel=interpol([0,50,100,210,210],[0.01,0.5,1,3,500],x)
;;;;


;;;; Simulate ;;;;
KinMS,xsize,ysize,vsize,cellsize,dv,beamsize,inc,sbprof=fx,sbrad=x,velrad=x1,velprof=vel,cubeout=f,nsamps=nsamps,_extra=_extra,intflux=30.,gassigma=gassigma,posang=270
;;;;


;;;; Create plot data from cube ;;;;
   mom0rot=total(f,3)
   x1=(findgen((xsize/cellsize))-((xsize/cellsize)/2.))*cellsize
   v1=(findgen(vsize/dv)-((vsize/dv)/2.))*dv
   mom1=mom0rot*0.0
   for i=0,xsize-1 do for j=0,ysize-1 do mom1[i,j]=mean(v1[where(f[i,j,*] eq max(f[i,j,*]))])
   mom1[where(mom0rot lt 0.1*max(mom0rot))] =-10000
   pvdcube=f
   pvd=reform(pvdcube[*,64,*])
   spec=(total(total(f,1),1))/(mean(beamsize)^2)
;;;;

;;;; Plot the results ;;;;
   device,decomposed=0,RETAIN=2
   loadct,39,/silent
   levs=v1[where(abs(v1) le max(vel*sin(!dtor*inc))+50.)]
   !p.multi=[0,2,2]
   contour,mom0rot,x1,x1,/cell_fill,levels=(max(mom0rot)*((findgen(9)+1)/10.)),color=0,background=255,/xstyle,/ystyle,xtitle="RA (arcsec)",ytitle="DEC (arcsec)"
   al_legend,["Moment 0"],box=0,/top,/right,textcolor=0
   contour,mom1,x1,x1,/cell_fill,levels=levs,color=0,background=255,/xstyle,/ystyle,xtitle="RA (arcsec)",ytitle="DEC (arcsec)" 
   al_legend,["Moment 1"],box=0,/top,/right,textcolor=0
   contour,pvd,x1,v1,levels=max(pvd)*((findgen(9)+1)/10.),/xstyle,/ystyle,xtitle="Position (arcsec)",ytitle="Velocity (km/s)",background=255,color=0,yrange=[min(levs)-30,max(levs)+30]
   al_legend,["Major Axis PVD"],box=0,/top,/right,textcolor=0
   plot,v1,spec,color=0,psym=10,xtitle="Velocity (km/s)",ytitle="Flux",/xstyle,/ystyle
   al_legend,["Spectrum"],box=0,/top,/right,textcolor=0
   !p.multi=[0,1,1]
;;;;
end
pro KinMStest_diskthick,_extra=_extra
;;;;;;;;;;;
;
; A test procedure to demonstrate the KinMS code, and check if it
; works on your system. This procedure demonstrates how to create a simulation of an
; exponential disk of molecular gas with a thickness that
; varies with radius. Any default paramaters
; can be changed by specifying them at the command line (see KinMS.pro 
; for the full details of all the avaliable parameters).
;
;
;;;;;;;;;;;



;;;; Setup cube parameters ;;;;
   xsize=128
   ysize=128
   vsize=1400
   cellsize=1
   dv=10
   beamsize=2.
   nsamps=5e5
;;;;

;;;; Set up exponential disk SB profile/velocity ;;;;
fcent=110.
scalerad=2000.
inc=90.
x=dindgen(1000)/10.
x1=x
fx = Fcent*exp(-x/scalerad)
diskthick=interpol([1,1,5,15],[0,10,15,20],x)
vel=interpol([0,50,100,210,210],[0.01,0.5,1,3,500],x)
;;;;


;;;; Simulate ;;;;
KinMS,xsize,ysize,vsize,cellsize,dv,beamsize,inc,sbprof=fx,sbrad=x,velrad=x1,velprof=vel,cubeout=f,nsamps=nsamps,_extra=_extra,intflux=30.,diskthick=diskthick,posang=270
;;;;


;;;; Create plot data from cube ;;;;
   mom0rot=total(f,3)
   x1=(findgen((xsize/cellsize))-((xsize/cellsize)/2.))*cellsize
   v1=(findgen(vsize/dv)-((vsize/dv)/2.))*dv
   mom1=mom0rot*0.0
   for i=0,xsize-1 do for j=0,ysize-1 do mom1[i,j]=mean(v1[where(f[i,j,*] eq max(f[i,j,*]))])
   mom1[where(mom0rot lt 0.1*max(mom0rot))] =-10000
   pvdcube=f
   pvd=reform(total(pvdcube[*,60:68,*],2))
   spec=(total(total(f,1),1))/(mean(beamsize)^2)
;;;;

;;;; Plot the results ;;;;
   device,decomposed=0,RETAIN=2
   loadct,39,/silent
   levs=v1[where(abs(v1) le max(vel*sin(!dtor*inc))+20.)]
   !p.multi=[0,2,2]
   contour,mom0rot,x1,x1,/cell_fill,levels=(max(mom0rot)*((findgen(9)+1)/10.)),color=0,background=255,/xstyle,/ystyle,xtitle="RA (arcsec)",ytitle="DEC (arcsec)"
   al_legend,["Moment 0"],box=0,/top,/right,textcolor=0
   contour,mom1,x1,x1,/cell_fill,levels=levs,color=0,background=255,/xstyle,/ystyle,xtitle="RA (arcsec)",ytitle="DEC (arcsec)" 
   al_legend,["Moment 1"],box=0,/top,/right,textcolor=0
   contour,pvd,x1,v1,levels=max(pvd)*((findgen(9)+1)/10.),/xstyle,/ystyle,xtitle="Position (arcsec)",ytitle="Velocity (km/s)",background=255,color=0,yrange=[min(levs)-30,max(levs)+30]
   al_legend,["Major Axis PVD"],box=0,/top,/right,textcolor=0
   plot,v1,spec,color=0,psym=10,xtitle="Velocity (km/s)",ytitle="Flux",/xstyle,/ystyle
   al_legend,["Spectrum"],box=0,/top,/right,textcolor=0
   !p.multi=[0,1,1]
;;;;
end

pro KinMStest_warp
;;;;;;;;;;;
;
; A test procedure to demonstrate the KinMS code, and check if it
; works on your system. This procedure demonstrates how to create a simulation of a
; warped exponential disk of molecular gas. Any default paramaters
; can be changed by specifying them at the command line (see KinMS.pro 
; for the full details of all the avaliable parameters).
;
;
;;;;;;;;;;;


;;;; Setup cube parameters ;;;;
   xsize=128
   ysize=128
   vsize=1400
   cellsize=1
   dv=10
   beamsize=4.
   nsamps=5e5
;;;;

;;;; Set up exponential disk SB profile/velocity ;;;;
Fcent=1.
scalerad=20.

x=findgen(400)/0.25
x1=x
fx = Fcent*exp(-x/scalerad)
inc=60.
vel=interpol([0,50,100,210,210],[0.01,0.5,1,3,500],x)
;;;;

;;;; Set up warp ;;;;
posang=interpol([270,270,300,300],[0.0,15,50,500],x)
;;;;

;;;; Simulate ;;;;
KinMS,xsize,ysize,vsize,cellsize,dv,beamsize,inc,sbprof=fx,sbrad=x,velrad=x1,velprof=vel,cubeout=f,nsamps=nsamps,_extra=_extra,intflux=30.,posang=posang
;;;;


;;;; Create plot data from cube ;;;;
   mom0rot=total(f,3)
   x1=(findgen((xsize/cellsize))-((xsize/cellsize)/2.))*cellsize
   v1=(findgen(vsize/dv)-((vsize/dv)/2.))*dv
   mom1=(mom0rot*0.0)-10000.0
   for i=0,(xsize/cellsize)-1 do for j=0,(ysize/cellsize)-1 do begin
      if mom0rot[i,j] gt 0.1*max(mom0rot) then begin
         mom1[i,j]=mean(v1[where(f[i,j,*] eq max(f[i,j,*]))])
      endif
   endfor
   pvdcube=f
   if keyword_set(_extra) then begin
      if TAG_EXIST(_extra, posang,/quiet) then for i=0,n_elements(v1)-1 do pvdcube[*,*,i]=rot(f[*,*,i],_extra.posang-90,/interp)
   endif
   pvd=reform(pvdcube[*,64,*])
   spec=(total(total(f,1),1))/(mean(beamsize)^2)
;;;;

;;;; Plot the results ;;;;
   device,decomposed=0,RETAIN=2
   loadct,39,/silent
   levs=v1[where(abs(v1) le max(vel*sin(!dtor*inc))+20.)]
   !p.multi=[0,2,2]
   contour,mom0rot,x1,x1,/cell_fill,levels=(max(mom0rot)*((findgen(9)+1)/10.)),color=0,background=255,/xstyle,/ystyle,xtitle="RA (arcsec)",ytitle="DEC (arcsec)",/iso
   al_legend,["Moment 0"],box=0,/top,/right,textcolor=0
   contour,mom1,x1,x1,/cell_fill,levels=levs,color=0,background=255,/xstyle,/ystyle,xtitle="RA (arcsec)",ytitle="DEC (arcsec)" ,/iso
   al_legend,["Moment 1"],box=0,/top,/right,textcolor=0
   contour,pvd,x1,v1,levels=max(pvd)*((findgen(9)+1)/10.),/xstyle,/ystyle,xtitle="Position (arcsec)",ytitle="Velocity (km/s)",background=255,color=0,yrange=[min(levs)-30,max(levs)+30]
   al_legend,["Major Axis PVD"],box=0,/top,/right,textcolor=0
   plot,v1,spec,color=0,psym=10,xtitle="Velocity (km/s)",ytitle="Flux",/xstyle,/ystyle
   al_legend,["Spectrum"],box=0,/top,/right,textcolor=0
   !p.multi=[0,1,1]
;;;;
end



pro KinMS_testsuite

;;;;;;;;;;;
;
;A master test procedure to demonstrate the KinMS code, running each
;individual test one after another. 
;
;;;;;;;;;;;
meep=""

print,"Test - simulate the gas ring in NGC4324"
KINMSTEST_NGC4324
print,"[Enter to continue]"
read,meep
print,"Test - simulate an exponential disk"
KINMSTEST_EXPDISK,10,50
print,"[Enter to continue]"
read,meep
print,"Test - using the INCLOUDS mechanism - unrealistic"
KINMSTEST_INCLOUDS
print,"[Enter to continue]"
read,meep
print,"Test - using the INCLOUDS mechanism - realistic"
KINMSTEST_INCLOUDS_SPIRAL
print,"[Enter to continue]"
read,meep
print,"Test - using a FITS file as input"
KINMSTEST_INFITS
print,"[Enter to continue]"
read,meep
print,"Test - using variable velocity dispersion"
KINMSTEST_VELDISP
print,"[Enter to continue]"
print,"Test - using variable disk thickness"
read,meep
KINMSTEST_DISKTHICK
print,"[Enter to continue]"
read,meep
print,"Test - simulate a warped exponential disk"
KINMSTEST_WARP

end
