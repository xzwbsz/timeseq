	program timesequence
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
c Last modification: Thu 07 Jun 2012 02:19:36 PM CEST  by Dongdong Zhang.
c
c
	IMPLICIT none
c
	character*128 dummychar
c
c        
	character*128 filename 
	character*128 filename1
	character*128 filename2
	character*128 arg 
        integer iarg,iargc,icm,indx  
        integer IXRED,IYRED,IZRED
	PARAMETER(IXRED=200,IYRED=200,IZRED=200)

	real*8 accx(IXRED,IYRED,IZRED)
	real*8 accy(IXRED,IYRED,IZRED)
	real*8 accz(IXRED,IYRED,IZRED)
	common /basis/ accx,accy,accz
c
	real*8 l_exc_dec
	common /basis/  l_exc_dec
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
	integer nt
	common /basis/ nt
c	
	integer nb,ne,nd
	common /basis/ nb,ne,nd
c	
	integer hit
	common /basis/ hit
c
c
        real*8 phase
        real*8 trigBU
	real*8 x,y,z,vx,vy,vz
	real*8 x_old,y_old,z_old,vx_old,vy_old,vz_old
	real*8 vx_end_dec,vx0
	integer m
c
	real*8 l_stage,l_exc_det,l_sk_hex
	real*8 l_dec_det,l_sk_dec,r_sk,l_x0_sk,l_nozzle_hex
        real*8 l_exc_hex, l_hex
	real*8 LA,LB,LC,L1,L2,L3,d,L4,L5,L6,L7,L8,L9,L10,Ladd
	real*8 l_exc_enddec
	common /dimensions/ l_exc_hex,l_hex
c
	real*8 mint,maxt
	real*8 minv,maxv
	real*8 mass,J2wavnr
	real*8 ns
	PARAMETER(ns=1e-9,J2wavnr=1/1.986447e-23)
c
	integer n,inunit,outunit,out2unit,nbswitch
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
        EXTERNAL F, F2, TRAP
c
	REAL*8 A,B,EPS,ETA
c
c
	integer i,j,k,ii
	integer dummyi,dummyj,dummyk
	real*8 dummy	
c
	real*8 t2jump(600),tof,tof_old,t_step,TOFFSET
	real*8 syncpositions(600)
c
c
c
	character*128 ptu
        common /ptu/ ptu
	character*4 hex_switchBU
        character*8 hexapole
	character*8 pulse
c	
	real*8 last_T2j	
c
c
	integer pline
	character*8 packetfilename
c
c
	integer load_integer
	real*8 load_real
	character*8 load_char
	character*128 load_longchar
	character*128 label
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		SETUP
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		INPUT FILES 
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
	if (iargc().eq.0) then
          write(*,*) 'Usage: -i1 "input file molecule" '
	  write(*,*) '       -i2 "input z positions to switch" '               
   	stop
c
	else
          icm=0
          iarg = 0
          indx = iargc()
            do while (iarg .lt. indx)
              iarg = iarg + 1
              call getarg(iarg,arg)
              if (arg .eq. '-i1') then
                iarg = iarg + 1
                call getarg(iarg,filename1)
              elseif (arg.eq. '-i2') then
                iarg = iarg + 1
                call getarg(iarg,filename2)
              end if
           end do
      	end if
c
        write(*,*) 'READING: ',filename1
c
c
	label='ni'              !Dimensions acceleration array
        ni=load_integer(label,filename1)
        label='nj'
        nj=load_integer(label,filename1)
        label='nk'
        nk=load_integer(label,filename1)
c Omit first and last data points; used only for fitting purposes:
	label='nbegin'
        nb=load_integer(label,filename1)
	label='neind'
        ne=load_integer(label,filename1)
	nd = ne-nb	
c Read dimensions
	label='gu'              !#gridunits/meter
        gu=load_integer(label,filename1)
c
c
	label='LA'		!Distance center package to skimmer
	l_x0_sk=load_real(label,filename1)
	label='LB'		!Dist. skimmer-hexapole
        l_sk_hex=load_real(label,filename1)
	l_nozzle_hex=l_x0_sk+l_sk_hex
        l_exc_hex=l_x0_sk+l_sk_hex
        l_sk_hex = l_x0_sk + l_sk_hex
c
	label='Hex'		!hexapole installed?
        hexapole=load_char(label,filename1)
	label='LC'		!hexapole (not used)
        LC=load_real(label,filename1)
	label='L1'		!hexapole (used)
        L1=load_real(label,filename1)
	l_hex = L1
	label='L2'               !hexapole (not used)
        L2=load_real(label,filename1)
	label='L3'		!hexapole-first stage
        L3=load_real(label,filename1)
c
	if(hexapole.eq.'yes') then
	   l_exc_dec  =  l_sk_hex+LC+L1+L2+L3
	else
	   l_exc_dec = l_sk_hex
	endif	
c
	l_stage = (nd/gu)	!Length of single stage
	label='L4'		!Distance decelerator to detection
        l_dec_det=load_real(label,filename1)
	label='nt'		!Number of stages (physical size)
        nt=load_integer(label,filename1)
	l_exc_det=l_exc_dec+(nt*l_stage)+l_dec_det
c
c
c Parameters used in calculation
	label='t_step'		!mass
        t_step=load_real(label,filename1)
	t_step=t_step*1.0D-9
	label='mass'		!mass
        mass=load_real(label,filename1)
        mass = mass*1.6605402D-27 ! mass in kg
c	
	label='trigBU'		!trigger time difference diss. and burst unit
	trigBU=load_real(label,filename1)
        trigBU=trigBU/1.0e6
	label='packet'		!beam or block or packet pulse
        pulse=load_char(label,filename1)
c
	if (pulse.eq.'beam') then
	   label='mean_v'	!average velocity incoming package
	   vx0=load_real(label,filename1)
	endif
	if (pulse.eq.'packet') then
	   label='pfile'	!file to load the packet
	   packetfilename=load_char(label,filename1)
	   write(*,*) packetfilename
	   label='pline'        !line to load 
           pline=load_integer(label,filename1)
	endif
c
c
	write(*,*) 'l_exc_hex             : ',l_exc_hex*1e3
	write(*,*) 'l_hex                 : ',l_hex*1e3
c
c
	write(*,*) ''
	write(*,*) 'Distance to detector: ', l_exc_det*1e3
	write(*,*) 'Velocity synchr. mol  ', vx0
	write(*,*) 'time resolution (ns): ', t_step*1e9
	write(*,*) 
	write(*,*) 'ni            :',ni
	write(*,*) 'nj            :',nj
	write(*,*) 'nk            :',nk
	write(*,*) 'nb            :',nb
	write(*,*) 'ne            :',ne
	write(*,*) 'nd            :',nd
	write(*,*)
	write(*,*) 'LA+LB         :',l_sk_hex*1000.
	write(*,*) 'LC            :',LC*1000.
	write(*,*) 'L1            :',L1*1000.
	write(*,*) 'L2            :',L2*1000.
	write(*,*) 'L3            :',L3*1000.
	write(*,*) 'l_stage       :',l_stage*1000.
	write(*,*) 'l_exc_dec     :',l_exc_dec*1000.
	write(*,*)
        write(*,*) '# electrodes -1:' ,nt
	write(*,*) 'mass          :',mass
c       x0      = 0
c	x0      = 0.10922		!x=0 at excitation region
c
c End reading parameters
c
        write(*,*) 'END OF READING PARAMETERS: ',filename1
c
c Begin to read the switching positions
c
	write(*,*) 'BEGIN READING POSITION: ',filename2
c
	inunit=33
        open(unit=inunit,file=filename2)
        read(inunit,*) nbswitch ! number of switch times
	do ii=1,nbswitch
	   read(inunit,*) syncpositions(ii)
	enddo
	close(inunit)
	write(*,*) 'END READING POSITION: ',filename2
c
c end reading switching positions
c	ns=1e-9
c	J2wavnr=1/1.986447e-23			!conversion of J to cm-1
c
	filename='../output/outax.dat'
	write(*,*)
	write(*,*) 'READING ACCELERATION FILE: ',filename
	inunit=32
        open(unit=inunit,file=filename)
        read(inunit,*) ni,nj,nk
        write(*,*)'Dimensions a_x array:'
        write(*,*) ni,nj,nk
        do i=1,ni
          do j=1,nj
            do k=1,nk
              read(inunit,*) dummyi,dummyj,dummyk,accx(i,j,k)
              accy(i,j,k)=0.   !neglect transverse forces
              accz(i,j,k)=0.
            end do
          end do
        end do
	close(inunit)
c
	write(*,*) 'DAT FILES READ'
c
	outunit  = 35
	out2unit = 36
	open(unit = outunit,file  = '../output/T2jump.out')
	write(*,*) 'OPENING T2jump.dat'
c
	open(unit = out2unit,file ='../output/T2jump.dat')
c 
        write(*,*) 'output written to: output/T2jump.dat'
        write(*,*) ''
	write(*,*) 'STARTING CALCULATION'
c	
        write(outunit,2008) vx0,vx_end_dec,nt,	! first time to write in T2jump.out file
     .                      phase,mass/1.6605402e-27
        write(outunit,2007)
        write(outunit,2007)
        write(outunit,2009)
c
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		INITIALISATION
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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
c	IFAIL = 0
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		START CALCULATION
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
	  tof = 0.
c
	  x   = 0.
	  y   = 0.
	  z   = 0.
c
	  vx  = vx0
	  vy  = 0.
	  vz  = 0.
c
	  hit = 0
	  last_T2j = 0.
c
	  ptu='free flight'
c	  ***************************
	  call ptu_2_hex(hex_switchBU)
c	  ***************************
          write(out2unit,2002) tof*1e6,(tof-last_T2j)*1e6,
     .          vx,(0.5*mass*vx**2*J2wavnr),x*1e3,
     .          hex_switchBU
c	  write(*,*) 'check initial values:',
c     .		x, vx, tof, hit,last_T2j 
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		I EXCITATION AND HEXAPOLE
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
          if(hexapole.eq.'no') then	
              x = x+l_sk_hex
	      tof = tof + l_sk_hex/vx
	      TOFFSET  =   tof-1.01*1e-6
	  else
	 	x   =        x+l_sk_hex+LC	
		tof =        tof + (l_sk_hex+LC)/vx	
	 	TOFFSET  =   tof-1.01*1e-6
	  	last_T2j =   tof	
	  	ptu='hexapole'
c	       ****************************
	  	call ptu_2_hex(hex_switchBU)
c	       ****************************
	  	write(outunit,1001) int((tof-TOFFSET)/ns),hex_switchBU
	  	write(out2unit,2001) tof*1e6,(tof-last_T2j)*1e6,
     .    	      vx,(0.5*mass*vx**2*J2wavnr),x*1e3,
     .          hex_switchBU
c	  write(*,*) 'check from nozzle to hexapole:',
c     .		x, vx, tof, hit,last_T2j 

c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		II HEXAPOLE
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
	  	last_T2j  = tof
	  	x   =  x + L1	
	  	tof = tof + (L1)/vx	
	  	ptu='free flight'
c	       ****************************
	  	call ptu_2_hex(hex_switchBU)
c	       ****************************
	  	write(outunit,1002) int((tof-TOFFSET)/ns),hex_switchBU
	  	write(out2unit,2002) tof*1e6,(tof-last_T2j)*1e6,
     .          	vx,(0.5*mass*vx**2*J2wavnr),x*1e3,
     .          	hex_switchBU
c	  write(*,*) 'check end of hexapole:',
c     .		x, vx, tof, hit,last_T2j 

c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		III HEXAPOLE TO DECELERATOR
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c	
	  	last_T2j  = tof
	  	x   =  x + L2 + L3	
	  	tof =  tof + (L2+L3)/vx
	  endif	 
c	  write(*,*) 'check from hexapole to decelerator:',
c     .		x, vx, tof, hit,last_T2j  	 	
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		IV STARK DECELERATOR
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c	 
c INITIALIZATION OF THE INTEGRATION
          do ii=1,15
             INFO(ii)=0
          enddo
          ATOL = 1D-7
          RTOL = 1D-7
c
	  odd = 1				! 90 degr. rotation of stages
c
c
         write(*,*) 'Number of switch times:', nbswitch          
c
	  do  i=1,nbswitch      ! begin loop over stages nbswitch=281	  	
c
	    if(i.lt.(nbswitch+1)) then       
              if(odd.eq.1) then
		  ptu = 'decelerator_vertical'
	       else
		  ptu = 'decelerator_horizontal'
	       endif	
            endif
c
c	    ****************************
	    call ptu_2_hex(hex_switchBU)
c	    ****************************	    
c	
            if(i.eq.1) then
	        write(outunit,1005) int((tof-TOFFSET)/ns),hex_switchBU
	        write(out2unit,2011) tof*1e6,(tof-last_T2j)*1e6,
     .              vx,(0.5*mass*vx**2*J2wavnr),x*1e3,
     .              hex_switchBU,i-1
            else
	          write(outunit,1003) int((tof-TOFFSET)/ns),hex_switchBU,
     .                i-1
  	          write(out2unit,2003) tof*1e6,(tof-last_T2j)*1e6,
     .                vx,(0.5*mass*vx**2*J2wavnr),x*1e3,
     .                hex_switchBU,i-1
            endif
c
	    last_T2j  = tof
c      
c
	    do while(x.lt.l_exc_dec+syncpositions(i))  !integrate to next stage
c
c (first stage, use s=1 always!!)
c
	      TSTART = tof			!Start time integration
	      TEND   = tof + t_step 		!increase with min. time step burst!!!
	      TWANT  = TEND
c IFAIL  = 0
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
c	******************************************************
        call DDERKF(F, NEQ, TSTART, YSTART, TWANT, INFO,
     .   RTOL, ATOL, IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR)
c	******************************************************
c
	tof = TWANT
c
        INFO(1)=0
        if(IDID.ne.2) then
c IDID=2 iff the integration went OK to the TEND
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
c HSTART = HNEXT	!Use last step size for next entry
c
	      else
c	      
	        write(*,*) 'ERROR HIT'
c
	      end if
c
	    end do 		!end while (integration) loop
c	  
	    odd = odd*(-1)
c  here there is a bug in the overtone version 24.05.2012
c	   
	    if (abs(x-(l_exc_dec+syncpositions(i))).gt.
     .	         abs(x_old-(l_exc_dec+syncpositions(i)))) then
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
c	  write(*,*) 'check every switching time:',
c     .		x, vx, tof, hit,last_T2j 
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
c
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c                     VII FLY TO DETECTOR
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
c last_T2j  = tof
c x = x - phase + L4    ! fly to center last electrodes
c tof = tof + (- phase + L4)/vx   ! phase is negative!!

          last_T2j = tof
	  l_exc_enddec = (l_exc_det-L4)-x
	  x = l_exc_det
	  tof = tof + (l_exc_enddec+L4)/vx

c  
          ptu='free flight'
c	  ****************************
          call ptu_2_hex(hex_switchBU)
c	  ****************************
	  write(out2unit,2013) tof*1e6,(tof-last_T2j)*1e6,
     .          vx,(0.5*mass*vx**2*J2wavnr),x*1e3,
     .          hex_switchBU
c	  write(*,*) 'check to the end of detector:',
c     .		x, vx, tof, hit,last_T2j 
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		Write last lines of files 
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c       
	write(outunit,2007)
	write(outunit,2007)
	write(outunit,2008) vx0,vx_end_dec,nt,
     .                      phase,mass/1.6605402e-27
c
	close(out2unit)
	close(outunit)
	
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		Formats output files 
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c		Force 
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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
	real*8 l_exc_dec
	common /basis/ l_exc_dec
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
	integer nt
	common /basis/ nt
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
	g0i = -l_exc_dec	
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
c write(*,*) 'xi',xi
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
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c			PTU
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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
           hex = '1000'
        endif
        if(ptu.eq.'decelerator_vertical') then
           hex = '2000'
        endif
        if(ptu.eq.'decelerator_MW') then
           hex = '0002'
        endif
	ptu='free flight'
        RETURN
	END
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c              load integer
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
        integer function load_integer(label,filename)
        implicit none
c
	character*128 label
	character*128 filename
c
	character*128 dummych1,dummych2,dummych3
	integer dummyint
	integer inunit,dummycount,il
c
	inunit=34
        open(unit=inunit,file=filename)
	dummych2='aaaaa'
	dummycount=0
        do while((dummych2.ne.label).and.(dummych1.ne.'end-of-file'))
	   read(inunit,*) dummych1,dummych2,dummych3
	   dummycount=dummycount+1
	end do
	close(inunit)
        if (dummych2.ne.label) then 
	   write(*,*) 'not found ',label
	endif
	open(unit=inunit,file=filename)
	do il=1,dummycount-1
	   read(inunit,*) dummych1,dummych2,dummych3
	end do
	read(inunit,*) dummyint
	close(inunit)
	load_integer=dummyint
        return
        end
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c              load real
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
        real*8 function load_real(label,filename)
        implicit none
c
        character*128 label
        character*128 filename
c
        character*128 dummych1,dummych2,dummych3
	real*8 dummyreal
        integer inunit,dummycount,il
c
        inunit=34
        open(unit=inunit,file=filename)
        dummych2='aaaaa'
        dummycount=0
        do while((dummych2.ne.label).and.(dummych1.ne.'end-of-file'))
           read(inunit,*) dummych1,dummych2,dummych3
           dummycount=dummycount+1
        end do
        close(inunit)
        if (dummych2.ne.label) then 
	   write(*,*) 'not found ',label
	endif
        open(unit=inunit,file=filename)
        do il=1,dummycount-1
           read(inunit,*) dummych1,dummych2,dummych3
        end do
        read(inunit,*) dummyreal
        close(inunit)
        load_real=dummyreal
        return
        end
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c              load char
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
        character*8 function load_char(label,filename)
        implicit none
c
        character*128 label
        character*128 filename
c
        character*128 dummych1,dummych2,dummych3
	character*8 dummych
        integer inunit,dummycount,il	
c
        inunit=34
        open(unit=inunit,file=filename)
        dummych2='aaaaa'
        dummycount=0
        do while((dummych2.ne.label).and.(dummych1.ne.'end-of-file'))
           read(inunit,*) dummych1,dummych2,dummych3
           dummycount=dummycount+1
        end do
        close(inunit)
        if (dummych2.ne.label) then 
	   write(*,*) 'not found ',label
	endif
        open(unit=inunit,file=filename)
        do il=1,dummycount-1
           read(inunit,*) dummych1,dummych2,dummych3
        end do
        read(inunit,*) dummych
        close(inunit)
        load_char=dummych
        return
        end
c
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c              load char
c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c
        character*8 function load_longchar(label,filename)
        implicit none
c
        character*128 label
        character*128 filename
c
        character*128 dummych1,dummych2,dummych3
        character*128 dummych
        integer inunit,dummycount,il
c
        inunit=34
        open(unit=inunit,file=filename)
        dummych2='aaaaa'
        dummycount=0
        do while((dummych2.ne.label).and.(dummych1.ne.'end-of-file'))
           read(inunit,*) dummych1,dummych2,dummych3
           dummycount=dummycount+1
        end do
        close(inunit)
        if (dummych2.ne.label) then 
	   write(*,*) 'not found ',label
	endif
        open(unit=inunit,file=filename)
        do il=1,dummycount-1
           read(inunit,*) dummych1,dummych2,dummych3
        end do
        read(inunit,*) dummych
        close(inunit)
        load_longchar=dummych
        return
        end
c
