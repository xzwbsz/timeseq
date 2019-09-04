	program timeseq_test
c
c Program to calculate switching
c times for deceleration of molecules in a Stark decelerator. Input file
c contains dimensions of acceleration array, selected part of that array,
c dimensions decelerator and machine, velocity, phase and mass of incoming
c molecule. (for instance timeseqOH.dat). Output generated: T2jump.out: 
c burst file that can be sent to the machine, and T2jump.dat: file 
c containing velocity, energy, etc for the different stages.
c
c
c Last modification: June 2011 by Janneke Blokland.
c
c
	IMPLICIT none
c
	character*128 dummychar
c
c        
        character*128 filename,arg        
        integer iarg,iargc,icm,indx   
        integer IXRED,IYRED,IZRED
	PARAMETER(IXRED=200,IYRED=200,IZRED=200)

	real*8 accx(IXRED,IYRED,IZRED),accy(IXRED,IYRED,IZRED)
	real*8 accz(IXRED,IYRED,IZRED)
	common /basis/ accx,accy,accz
c
	real*8 LtoDEC
	common /basis/  LtoDEC
c
	real*8 gu		
	common /basis/ 	gu
c
	integer ni,nj,nk
	common /basis/ ni,nj,nk
c
	integer odd
	common /basis/ odd
c
	integer nt,nbunch,ntotal
	common /basis/ nt,nbunch,ntotal
c	
	integer nb,ne,nd
	common /basis/ nb,ne,nd
c	
	integer hit
	common /basis/ hit
c
c
        real*8 phase,phase_deg,phase2,phase_deg2
        real*8 trigBU,v_buncher
	real*8 x,y,z,vx,vy,vz
	real*8 x_old,y_old,z_old,vx_old,vy_old,vz_old
	real*8 vx_end_dec,vx_begin,m
	real*8 finalx,finaly,finalz,finaltof
	real*8 r,theta,vr,vtheta
c
	real*8 mint,maxt
	real*8 minv,maxv
	real*8 mass,J2wavnr
c
	integer n,inunit,outunit,out2unit,out3unit,nseg,nswitch
c
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c       INTEGRATION PACKAGE
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        INTEGER NEQ,LRW,LIW
        PARAMETER (NEQ=6, LRW=33+7*NEQ, LIW=34)
        REAL*8 YSTART(NEQ)
        REAL*8 RWORK(LRW)
        INTEGER IWORK(LIW)
        REAL*8  TSTART,T,TEND, TWANT
        INTEGER INFO(15)
        INTEGER IDID
        REAL*8 RTOL, ATOL
        REAL*8 RPAR(10)
        INTEGER IPAR(10)
        EXTERNAL DDERKF
        EXTERNAL F, TRAP
c
	REAL*8 A,B,EPS,ETA
	REAL*8 V
c
c
	integer i,j,k,ii
	integer dummyi,dummyj,dummyk
	real*8 dummy
	integer overtone
	real*8 s
c
	real*8 t2jump(300),tof,tof_old,TOFFSET,t_step
c
	real*8 LA,LB,LC,L1,L2,L3,d,L4,L5,L6,L7,L8,L9,L10,Ladd
	real*8 Ltot,LtoLAST,LtoEndDec
c
	real*8 kf,kfoverm
	common /beter/ kfoverm
c
	logical log_switchBU1(0:3)
	logical log_switchBU2(0:3)
	logical log_switchBU(0:3)
	integer int_switchBU
	character*128 ptu
        common /ptu/ ptu
	character*4 hex_switchBU
	character*4 hex
	character*8 pulse
        integer int_switchBU1,int_switchBU2
	real*8 ns
c	
	real*8 last_T2j	
c
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		SETUP
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		INPUT FILES 
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
	if (iargc().eq.0) then
          write(*,*) 'Usage: timeseq.exe -i "input file molecule" '             
   	stop
c
	else
          icm=0
          iarg = 0
          indx = iargc()
            do while (iarg .lt. indx)
	      iarg = iarg + 1
	      call getarg(iarg,arg)
              if (arg .eq. '-i') then
                iarg = iarg + 1
                call getarg(iarg,filename)
              end if
	   end do
      	end if
c
        write(*,*) 'READING: ',filename
c
	inunit=34
        open(unit=inunit,file=filename)
c
        read(inunit,*) dummychar
        read(inunit,*) dummychar
        read(inunit,*) dummychar
c       Simion array: 
        read(inunit,*) ni			!Dimensions acceleration array 
	read(inunit,*) nj
	read(inunit,*) nk
	read(inunit,*) nb			!omit first data points
	read(inunit,*) ne			!omit last data points
	    nd = ne-nb
	read(inunit,*) gu			!#gridunits/meter
c
        read(inunit,*) dummychar
        read(inunit,*) dummychar
        read(inunit,*) dummychar	
c       Hexapole
        read(inunit,*) dummy
        read(inunit,*) dummy
        read(inunit,*) dummy
        read(inunit,*) dummy
        read(inunit,*) dummy
c
        read(inunit,*) dummychar
        read(inunit,*) dummychar
        read(inunit,*) dummychar
c       Physical size set-up
        read(inunit,*) dummy
	read(inunit,*) LA			!nozzle-skimmer
	read(inunit,*) LB			!skimmer-hexapole
	LB=LA+LB
        read(inunit,*) hex
	read(inunit,*) LC			!hex not used
	read(inunit,*) L1			!hex. used
	read(inunit,*) L2			!hex not used
	read(inunit,*) L3			!hex.-first stage
             if(hex.eq.'yes') then
	         LtoDEC  =  LB+LC+L1+L2+L3
	     else
	         LtoDEC = LB
	     endif	 	 	
	read(inunit,*) L4			!end decel-detection
	read(inunit,*) nt	        	!Total number of stages in dec.
	d  = dble(nd)/dble(gu)
        Ltot=LtoDEC+d*nt+L4
        read(inunit,*) dummy			!number of switching times dec. 
c
	read(inunit,*) s                        !overtone number
c            
        read(inunit,*) dummychar
        read(inunit,*) dummychar
        read(inunit,*) dummychar	
c       Parameters used in calculation:
        read(inunit,*) vx			 !velocity synchronous particle
	vx_begin=vx
	read(inunit,*) t_step			 !time stepsize
	     t_step=t_step*1.0D-9
	read(inunit,*) phase_deg 		 !phase	
           phase =  (d/2) * ( (phase_deg/90.) -1.) !converted to distance
 	read(inunit,*) mass			 !Mass of molecule used in AMU
	     mass=mass*1.6605402e-27			!Convert AMU to kg
c
        read(inunit,*) dummychar
        read(inunit,*) dummychar
        read(inunit,*) dummychar
c       Parameters used in fly
        read(inunit,*) dummy
        read(inunit,*) dummy
        read(inunit,*) pulse
        read(inunit,*) dummy
        read(inunit,*) dummy
        read(inunit,*) dummy
        read(inunit,*) dummy
        read(inunit,*) dummy
	read(inunit,*) dummy
        read(inunit,*) dummy
        read(inunit,*) dummy
        read(inunit,*) dummy
c
        close(inunit)
c
	write(*,*) ''
	write(*,*) 'Distance to detector: ', Ltot
	write(*,*) 'Velocity synchr. mol  ', vx
	write(*,*) 'Phase (deg.)          ',phase_deg
	write(*,*) 'Phase (mm):           ',phase
	write(*,*) 'time resolution (ns): ', t_step*1e9
	write(*,*) 
	write(*,*) 'ni            :',ni
	write(*,*) 'nj            :',nj
	write(*,*) 'nk            :',nk
	write(*,*) 'nb            :',nb
	write(*,*) 'ne            :',ne
	write(*,*) 'nd            :',nd
	write(*,*)
	write(*,*) 'LA+LB         :',LB*1000.
	write(*,*) 'LC            :',LC*1000.
	write(*,*) 'L1            :',L1*1000.
	write(*,*) 'L2            :',L2*1000.
	write(*,*) 'L3            :',L3*1000.
	write(*,*) 'nt            :',nt
        write(*,*) 'nbunch        :',nbunch
	write(*,*) 'd             :',d*1000.
	write(*,*) 'L5            :',L5*1000.
	write(*,*) 'L6            :',L6*1000.
	write(*,*) 'LtoDEC        :',LtoDEC*1000.
	write(*,*) 'Ltot          :',Ltot*1000.
	write(*,*)
	write(*,*) 'vx            :',vx
	write(*,*) 't_step        :',t_step
	write(*,*) 'phase (mm)    :',phase*1000.
        write(*,*) '# electrodes -1:' ,nt
	write(*,*) 'mass          :',mass
c
	ns=1e-9
	J2wavnr=1/1.986447e-23			!conversion of J to cm-1
c				
c	
	filename='../output/outax.dat'
	write(*,*)
	write(*,*) 'READING ACCELERATION FILE: ',filename
c
	open(unit=inunit,file=filename)
	read(inunit,*) ni,nj,nk
	write(*,*)'Dimensions a_x array:'
	write(*,*) ni,nj,nk
c
	do i=1,ni
	  do j=1,nj
	    do k=1,nk
	      read(inunit,*) dummyi,dummyj,dummyk,accx(i,j,k)
              accy(i,j,k)=0.   !neglect transverse forces
	      accz(i,j,k)=0.
            end do
	  end do
	end do 
c
	close(inunit)
c
	write(*,*) 'DAT FILES READ'
c
	outunit  = 35
	out2unit = 36
	out3unit = 37
	open(unit = outunit,file  = 
     >	'../output/T2jump.out')
	write(*,*) 'OPENING T2jump.dat'
c
	open(unit = out2unit,file = 
     >	'../output/T2jump.dat')
c
	open(unit = out3unit,file = 
     >	'../output/P2jump.dat')
c
        write(*,*) 'output written to: ../output/T2jump.dat'
        write(*,*) ''
	write(*,*) 'STARTING CALCULATION'
c	
        write(outunit,2008) vx_begin,vx_end_dec,nt,
     .                      phase_deg,mass/1.6605402e-27
        write(outunit,2007)
        write(outunit,2007)
        write(outunit,2009)
c
        write(*,*) 'phase decelerator (deg.):', phase_deg
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		INITIALISATION
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
	maxv = -1.0d60
	minv =  1.0d60
c
	maxt = -1.0d60
	mint =  1.0d60
c
	A = -1.0d59
	B =  1.0d59
	EPS = 1.0e-6
	ETA = 0.0e0
c
c       XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		START CALCULATION
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
	  tof = 0.
c
	  x   = 0.
	  y   = 0.
	  z   = 0.
c
	  vx  = vx_begin
	  vy  = 0.
	  vz  = 0.
c
	  hit = 0
	  last_T2j = 0.
c
	  ptu='free flight'
	  call ptu_2_hex(hex_switchBU)
          write(out2unit,2002) tof*1e6,(tof-last_T2j)*1e6,
     .          vx,(0.5*mass*vx**2*J2wavnr),x*1e3,
     .          hex_switchBU
c
          write(out3unit,*) x-0.09836
c
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		I EXCITATION AND HEXAPOLE
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
          if(hex.eq.'no') then	
              x = x+LB
	      tof = tof + LB/vx
	      TOFFSET  =   tof-1.01*1e-6
	  else
	 	x   =        x+LB+LC	
		tof =        tof + (LB+LC)/vx	
	 	TOFFSET  =   tof-1.01*1e-6
	  	last_T2j =   tof	
	  	ptu='hexapole'
	  	call ptu_2_hex(hex_switchBU)
	  	write(outunit,1001) int((tof-TOFFSET)/ns),hex_switchBU
	  	write(out2unit,2001) tof*1e6,(tof-last_T2j)*1e6,
     .    	      vx,(0.5*mass*vx**2*J2wavnr),x*1e3,
     .          hex_switchBU
                write(out3unit,*) x-0.09836
c
c		XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c			II HEXAPOLE
c		XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
	  	last_T2j  = tof
	  	x   =  x + L1	
	  	tof = tof + (L1)/vx	
	  	ptu='free flight'
	  	call ptu_2_hex(hex_switchBU)
	  	write(outunit,1002) int((tof-TOFFSET)/ns),hex_switchBU
	  	write(out2unit,2002) tof*1e6,(tof-last_T2j)*1e6,
     .          	vx,(0.5*mass*vx**2*J2wavnr),x*1e3,
     .          	hex_switchBU
                write(out3unit,*) x-0.09836
c
c		XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c			III HEXAPOLE TO DECELERATOR
c		XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c	
	  	last_T2j  = tof
	  	x   =  x + L2 + L3	
	  	tof =  tof + (L2+L3)/vx
	  endif	  	 	
c
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		IV STARK DECELERATOR
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c	 
c        INITIALIZATION OF THE INTEGRATION
          do ii=1,15
             INFO(ii)=0
          enddo
          ATOL = 1D-7
          RTOL = 1D-7
c
	  odd = 1				! 90 degr. rotation of stages
c
         if(s.eq.1) then
	    nswitch = nt
	 else
	    nswitch = int((nt-1)/s) + 1   
         endif
c
         write(*,*) 'Overtone number s:', s
         write(*,*) 'Number of switch times:', nswitch          
c
	  do  i=1,nswitch      ! begin loop over stages 	  	
c
	    if(i.lt.(nswitch+1)) then       
              if(odd.eq.1) then
		  ptu = 'decelerator_vertical'
	       else
		  ptu = 'decelerator_horizontal'
	       endif	
            endif
c
c
	    call ptu_2_hex(hex_switchBU)	    
c	
            if(i.eq.1) then
	        write(outunit,1005) int((tof-TOFFSET)/ns),hex_switchBU
	        write(out2unit,2011) tof*1e6,(tof-last_T2j)*1e6,
     .              vx,(0.5*mass*vx**2*J2wavnr),x*1e3,
     .              hex_switchBU,i-1
	        write(out3unit,*) x-0.09836
            else
	          write(outunit,1003) int((tof-TOFFSET)/ns),hex_switchBU,
     .                i-1
  	          write(out2unit,2003) tof*1e6,(tof-last_T2j)*1e6,
     .                vx,(0.5*mass*vx**2*J2wavnr),x*1e3,
     .                hex_switchBU,i-1
                  write(out3unit,*) x-0.09836
            endif
c
	    last_T2j  = tof
c
            if(i.eq.1) then  !(first stage, use s=1 always!!!)
	       overtone=1
	    else
	       overtone=s
	    endif      
c
	    do while((x-(d+overtone*(dfloat(i)-1)*d)-LtoDEC).lt.phase) !integrate to next stage
c
c                 (first stage, use s=1 always!!)
c
	      TSTART = tof			!Start time integration
	      TEND   = tof + t_step 		!increase with min. time step burst!!!
	      TWANT  = TEND
c	      IFAIL  = 0
c	   
	      if(odd.eq.1) then
	        YSTART(1) = x			!Initial conditions
	        YSTART(2) = y
	        YSTART(3) = z
c	
	        YSTART(4) = vx
	        YSTART(5) = vy
	        YSTART(6) = vz
	      else
	        YSTART(1) = x 	!switching equals 90 degr. rotation (change y
	        YSTART(2) = z	!and z) + shift over length of single stage.
	        YSTART(3) = y   !Interchanging  is done here, shifting in F(..)
c	
	        YSTART(4) = vx
	        YSTART(5) = vz
	        YSTART(6) = vy
	      endif	
c
c
	      if(hit.eq.0) then
		x_old=x		! save data previous step
		y_old=y
		z_old=z
		vx_old=vx
		vy_old=vy
		vz_old=vz
		tof_old = tof
c
        call DDERKF(F, NEQ, TSTART, YSTART, TWANT, INFO,
     .   RTOL, ATOL, IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR)
c
	tof = TWANT
c
        INFO(1)=0
        if(IDID.ne.2) then
c          IDID=2 iff the integration went OK to the TEND
           write(*,*) 'ERROR T'
        endif
c
c
	        if(odd.eq.1) then
	          x  = YSTART(1)
	          y  = YSTART(2)
	          z  = YSTART(3)
c
	          vx = YSTART(4)
	          vy = YSTART(5)
	          vz = YSTART(6)
	        else 
	          x  = YSTART(1)
	          y  = YSTART(3)
	          z  = YSTART(2)
c
    	          vx = YSTART(4)
	          vy = YSTART(6)
	          vz = YSTART(5)
  	        end if 
c
                if (vx.le.0.0) then
		   stop 'velocity less than zero in Stark dec.'
		endif
c
c	        HSTART = HNEXT	!Use last step size for next entry
c
	      else
c	      
	        write(*,*) 'ERROR HIT'
c
	      end if
c
	    end do 		!end while (integration) loop
c	  
	    odd = odd*-1
c
	    if (abs((x-(overtone*dfloat(i)*d)-LtoDEC)-phase).gt.
     .          abs((x_old-(overtone*dfloat(i)*d)-LtoDEC)-phase)) then	   
c
	        tof =  tof_old		!check for timing closest to desired
		x =  x_old		!phase.		
		y  = y_old
		z  = z_old
		vx = vx_old
		vy = vy_old
		vz = vz_old
	    endif
c
	  end do 	!End number of stages loop
c
	  vx_end_dec=vx
c
          ptu='free flight'
          call ptu_2_hex(hex_switchBU)
	  write(outunit,1003) int((tof-TOFFSET)/ns),hex_switchBU,
     . 	        i-1
	  write(out2unit,2003) tof*1e6,(tof-last_T2j)*1e6,
     .          vx,(0.5*mass*vx**2*J2wavnr),x*1e3,
     .          hex_switchBU,i-1
          write(out3unit,*) x-0.09836
c
c       XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c                     VII FLY TO DETECTOR
c       XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
c	  last_T2j  = tof
c          x = x - phase + L4    ! fly to center last electrodes
c          tof = tof + (- phase + L4)/vx   ! phase is negative!!

          last_T2j = tof
	  LtoEndDec = (Ltot-L4)-x
	  x = Ltot
	  tof = tof + (LtoEndDec+L4)/vx

c  
          ptu='free flight'
          call ptu_2_hex(hex_switchBU)
	  write(out2unit,2013) tof*1e6,(tof-last_T2j)*1e6,
     .          vx,(0.5*mass*vx**2*J2wavnr),x*1e3,
     .          hex_switchBU
          write(out3unit,*) x-0.09836
c
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		Write last lines of files 
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c       
	write(outunit,2007)
	write(outunit,2007)
	write(outunit,2008) vx_begin,vx_end_dec,nt,
     .                      phase_deg,mass/1.6605402e-27
c
	close(out2unit)
	close(outunit)
	close(out3unit)
	
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		Formats output files 
c	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
1001	FORMAT(i12,5x,'0x',A,5x,'!  HEXAPOLE')
1002    FORMAT(i12,5x,'0x',A,5x,'!  ')
1003    FORMAT(i12,5x,'0x',A,5x,'!  ',i3)
1004    FORMAT(i12,5x,'0x',A,5x,'!  TRAP')
1005    FORMAT(i12,5x,'0x',A,5x,'!  ',i3,'   DEC')
1006    FORMAT(i12,5x,'0x',A,5x,'!  ',i3,'   BUNCHER')
2001	FORMAT(2(f8.1,1x),f7.1,1x,2(f7.2,1x),A,' !  HEXAPOLE')
2002    FORMAT(2(f8.1,1x),f7.1,1x,2(f7.2,1x),A,' !  ')
2003    FORMAT(2(f8.1,1x),f7.1,1x,2(f7.2,1x),A,' !  ',i3)
2011    FORMAT(2(f8.1,1x),f7.1,1x,2(f7.2,1x),A,' !  ',i3,'   DEC')
2012    FORMAT(2(f8.1,1x),f7.1,1x,2(f7.2,1x),A,' !  ',i3,'   BUNCHER')
2013    FORMAT(2(f8.1,1x),f7.1,1x,2(f7.2,1x),A,' !  DETECTION')
2007    FORMAT('#')
2008	FORMAT('# vx_i:',f7.2,1x,'vx_f:',f7.2,1x,'stages:',i3,1x,
     .          'phase:',f7.2,1x,'mass:',f7.2)
2009	FORMAT('[1]') 
c
c
	stop 'normal termination'
	end
c	
c
c
	SUBROUTINE F(T,Y,YP,RPAR,IPAR)
c	
        real*8 RPAR
        integer IPAR
c
	REAL*8 	T
	REAL*8  Y(*),YP(*)
c	
	REAL*8  YP1,YPN
c
	real*8 accx(200,200,200),accy(200,200,200),accz(200,200,200)
	common /basis/ accx,accy,accz
c
	real*8 LtoDEC
	common /basis/ LtoDEC
c
	real*8 gu
	common /basis/ gu
c
	integer  ni,nj,nk
	common /basis/ ni,nj,nk	
c		
	integer odd
	common /basis/ odd
c
	integer nt,nbunch,ntotal
	common /basis/ nt,nbunch,ntotal
c	
	integer nb,ne,nd
	common /basis/ nb,ne,nd
c
	integer hit
	common /basis/ hit
c
c
	integer i,j,k,xi,yj,zk
c	  
	g0i = -LtoDEC	
	g0j = dfloat(nj+1)/(2.*gu)
	g0k = dfloat(nk+1)/(2.*gu)
c
	xi = dint((Y(1)+g0i)*gu)	!position in decelerator in gridunits,
	yj = dint((Y(2)+g0j)*gu)	!center-offset taken into account by
	zk = dint((Y(3)+g0k)*gu)	!g0i-g0k
c
c
	if(((xi.ge.0).and.(xi.le.nt*nd)).and. 	!check if inside decelerator
     .	   ((yj.ge.1).and.(yj.le.nj)).and.
     .	   ((zk.ge.1).and.(zk.le.nk)))   then
c
c
	  xx = (Y(1)+g0i)*gu - dfloat(xi)	!between 0 and 1
	  yy = (Y(2)+g0j)*gu - dfloat(yj)
	  zz = (Y(3)+g0k)*gu - dfloat(zk)
c
c	
	if(odd.eq.-1) then			 !check odd or even (vertical or
	  xi=xi+nd			        !horizontal) stages switched
c
	endif
c
c
c	write(*,*) 'xi',xi
c 
	  xi = mod(xi,(2*nd))
c
c	
c 
	  if(xi.lt.nd) then  
c
	    xi = xi+nb
c
c
	YP(4)=((1.0 - xx)*(1.0 - yy)*(1.0 - zz)*accx(xi,yj,zk)
     .         +(xx)      *(1.0 - yy)*(1.0 - zz)*accx(xi+1,yj,zk)
     .         +(1.0 - xx)*(yy)*      (1.0 - zz)*accx(xi,yj+1,zk)
     .         +(xx)      *(yy)*      (1.0 - zz)*accx(xi+1,yj+1,zk)
     .         +(1.0 - xx)*(1.0 - yy)*(zz)      *accx(xi,yj,zk+1) 
     .         +(xx)      *(1.0 - yy)*(zz)      *accx(xi+1,yj,zk+1) 
     .         +(1.0 - xx)*(yy)      *(zz)      *accx(xi,yj+1,zk+1) 
     .         +(xx)      *(yy)      *(zz)      *accx(xi+1,yj+1,zk+1))
c
c	
c
	YP(5)=((1.0 - xx)*(1.0 - yy)*(1.0 - zz)*accy(xi,yj,zk)
     .         +(xx)      *(1.0 - yy)*(1.0 - zz)*accy(xi+1,yj,zk)
     .         +(1.0 - xx)*(yy)*      (1.0 - zz)*accy(xi,yj+1,zk)
     .         +(xx)      *(yy)*      (1.0 - zz)*accy(xi+1,yj+1,zk)
     .         +(1.0 - xx)*(1.0 - yy)*(zz)      *accy(xi,yj,zk+1) 
     .         +(xx)      *(1.0 - yy)*(zz)      *accy(xi+1,yj,zk+1) 
     .         +(1.0 - xx)*(yy)      *(zz)      *accy(xi,yj+1,zk+1) 
     .         +(xx)      *(yy)      *(zz)      *accy(xi+1,yj+1,zk+1))
c
c	
c
	YP(6)=((1.0 - xx)*(1.0 - yy)*(1.0 - zz)*accz(xi,yj,zk)
     .         +(xx)      *(1.0 - yy)*(1.0 - zz)*accz(xi+1,yj,zk)
     .         +(1.0 - xx)*(yy)*      (1.0 - zz)*accz(xi,yj+1,zk)
     .         +(xx)      *(yy)*      (1.0 - zz)*accz(xi+1,yj+1,zk)
     .         +(1.0 - xx)*(1.0 - yy)*(zz)      *accz(xi,yj,zk+1) 
     .         +(xx)      *(1.0 - yy)*(zz)      *accz(xi+1,yj,zk+1) 
     .         +(1.0 - xx)*(yy)      *(zz)      *accz(xi,yj+1,zk+1) 
     .         +(xx)      *(yy)      *(zz)      *accz(xi+1,yj+1,zk+1))
c	
	else 
c
	  xi = nb + (nd - (xi-nd))
c
c Note: Fx * -1
c Stage is anti-symmetric with respect to rods at HV
c
	 YP(4)=-((1.0 - xx)*(1.0 - yy)*(1.0 - zz)*accx(xi,yj,zk)
     .         +(xx)      *(1.0 - yy)*(1.0 - zz)*accx(xi-1,yj,zk)
     .         +(1.0 - xx)*(yy)*      (1.0 - zz)*accx(xi,yj+1,zk)
     .         +(xx)      *(yy)*      (1.0 - zz)*accx(xi-1,yj+1,zk)
     .         +(1.0 - xx)*(1.0 - yy)*(zz)      *accx(xi,yj,zk+1) 
     .         +(xx)      *(1.0 - yy)*(zz)      *accx(xi-1,yj,zk+1) 
     .         +(1.0 - xx)*(yy)      *(zz)      *accx(xi,yj+1,zk+1) 
     .         +(xx)      *(yy)      *(zz)      *accx(xi-1,yj+1,zk+1))
c
c	
c
	YP(5)=((1.0 - xx)*(1.0 - yy)*(1.0 - zz)*accy(xi,yj,zk)
     .         +(xx)      *(1.0 - yy)*(1.0 - zz)*accy(xi-1,yj,zk)
     .         +(1.0 - xx)*(yy)*      (1.0 - zz)*accy(xi,yj+1,zk)
     .         +(xx)      *(yy)*      (1.0 - zz)*accy(xi-1,yj+1,zk)
     .         +(1.0 - xx)*(1.0 - yy)*(zz)      *accy(xi,yj,zk+1) 
     .         +(xx)      *(1.0 - yy)*(zz)      *accy(xi-1,yj,zk+1) 
     .         +(1.0 - xx)*(yy)      *(zz)      *accy(xi,yj+1,zk+1) 
     .         +(xx)      *(yy)      *(zz)      *accy(xi-1,yj+1,zk+1))
c
c	
c
	YP(6)=((1.0 - xx)*(1.0 - yy)*(1.0 - zz)*accz(xi,yj,zk)
     .         +(xx)      *(1.0 - yy)*(1.0 - zz)*accz(xi-1,yj,zk)
     .         +(1.0 - xx)*(yy)*      (1.0 - zz)*accz(xi,yj+1,zk)
     .         +(xx)      *(yy)*      (1.0 - zz)*accz(xi-1,yj+1,zk)
     .         +(1.0 - xx)*(1.0 - yy)*(zz)      *accz(xi,yj,zk+1) 
     .         +(xx)      *(1.0 - yy)*(zz)      *accz(xi-1,yj,zk+1) 
     .         +(1.0 - xx)*(yy)      *(zz)      *accz(xi,yj+1,zk+1) 
     .         +(xx)      *(yy)      *(zz)      *accz(xi-1,yj+1,zk+1))
	  endif
c
c 
	else 
c	  
	  YP(4) = 0
	  YP(5) = 0
	  YP(6) = 0
c
     	  if( (((yj.lt.1).or.(yj.gt.nj)).or.
     .	      ((zk.lt.1).or.(zk.gt.nk))).and.
     .         (xi.le.(nt+2)*nd.and.xi.ge.0))   then
c
	    hit = 1
c
	  endif
c
c	
	endif  
c		
	YP(1) = Y(4)
	YP(2) = Y(5)
	YP(3) = Y(6)
c	
	RETURN
	END
c
c
        SUBROUTINE ptu_2_hex(hex)
	character*128 ptu
	common /ptu/ ptu
	character*4 hex
	hex='0000'
	if(ptu.eq.'free flight') then
	   hex = '0000'
	endif 
	if(ptu.eq.'hexapole') then
	   hex = '0001'
	endif
        if(ptu.eq.'decelerator_horizontal') then
           hex = '0020'
        endif
        if(ptu.eq.'decelerator_vertical') then
           hex = '0010'
        endif
	ptu='free flight'
        RETURN
	END
c
