      MODULE share_dim
         IMPLICIT NONE
!   System Length/Cells  ==>   nx=2**6=64 
         INTEGER, PARAMETER :: lpx=9,nx=2**lpx, nxp1=nx+1, nxh=nx/2
         INTEGER, PARAMETER :: nsp=2,npg1=512,npg2=512
         INTEGER, PARAMETER :: npg=npg1+npg2,np=nx*npg
         INTEGER, PARAMETER :: ihparm=18,nheadr=ihparm+8*nsp
      END MODULE share_dim
!
      MODULE share_data
         USE share_dim
!
         CHARACTER(50) :: home
         REAL, DIMENSION(nheadr) :: header
         INTEGER, DIMENSION(nsp) :: nions
!         REAL, DIMENSION(nsp) :: umx,umy,umz,vx2,vy2,vz2,vken,tpal,tper
         REAL, DIMENSION(nsp) :: rdn,betain,vdr,anis,qi,ai,qm,gg,vthpal,vthper
         REAL, DIMENSION(0:nxp1) :: by,bz
         REAL, DIMENSION(0:nxp1) :: ex,ey,ez
         REAL, DIMENSION(np) :: xp0,xp1,vxp,vyp,vzp,qmr
         REAL, DIMENSION(3) :: bf0
!         REAL :: bfldx,bfldy,bfldz,efldx,efldy,efldz
!         REAL :: bfldenrgy,efldenrgy,vkinenrgy,totenrgy
         REAL :: xlen,dx,dt,tmax,betaen,wpiwci,dus,bnamp,slx,eta
         REAL :: bx,game
         INTEGER :: itmax,lfld,ifield,iparticles,ienergy,ifilter,itrestart
         INTEGER :: it1,it2,it3,it4,iseed
      END MODULE share_data
!
      PROGRAM hyb1d_rrkfft
!
!***********************************************************************
!      Hybrid simulation code with RRK time-advancing for particles    *
!      and fields for experiment with multi-ion-system. The code       *
!      uses the Rational Runge-Kutta (RRK) for time evolution and the  *
!      FFT algorithm for spatial derivative.                           *
!                                                                      *
!         Created by Adolfo F. Vinas at NASA-GSFC                      *
!                            (05/23/1993) for <1D-x, 3D-v, 3-comp.>    *
!***********************************************************************
!
! Program to Solve the Electromagnetic Particle-In-Cell (PIC) hybrid simulation
! in 1D-x and 3D-v explicitly using pseudo-spectral Fourier Methods. 
! Ions are kinetics, electrons are a massless fluid (me=0). Periodic boundary 
! conditions are assumed.
!
! Spatial
!  Grid:  The setup for the grid is as follows: The length of the physical grid is XL.
!         and the physical grid lies between 1 to Nx. There are two ghost cells: one at
!         location 0 and the other at location Nx+1. Thus we have:
!
!                    |------------------- L -------------------|
!
!         -----|-----|-----|-----|-----...... -----|-----|-----|-----|-----|
!              0     1     2     3               Nx-2  Nx-1   Nx   Nx+1
!
! NORMALIZATION: (We use CGS-units)
!  e = Electron charge
!  mp = proton mass
!  c = speed of light
!  betae0 = 8*pi*No*Tp0/(Bo*Bo) ==>based on the total electron density & temperature
!  Wthp = Proton thermal velocity Wthp=Sqrt(2*Tp0/mp)
!
!  t = Time is normalized to the proton cyclotron frequency : Wcp=e*Bo/(mp*c)
!        t = t_units * Wcp
!  x = position is normalized to the proton inertial length Lp=c/Wpp = Va/Wcp
!        x = x_units/Lde
!  v = All velocities are normalized to the Alfven velocity Va=Bo/Sqrt(4*pi*No*mp)
!        v = v_units/Va
!  B = magnetic field is normalized to the mean field Bo
!        B = (B_units/Bo)
!  E = since the electric & magnetic fields in CGS have the same units,
!      the electric field is normalized as the magnetic field.
!       E = c*E_units/(Va*Bo)
!
! VARIABLES:
!  xp(n) = particle position for all particles, n=1,npt
!  vxp(n),vyp(n),vzp(n) = particle velicity components, n=1,npt
!  bx =constant & by(i), bz(i) = magnetic field components, i=1,nx
!  ex(i), ey(i), ez(i) = electric field components, i=1,nx
!  rho(i) = total ion moment density, i=1,nx
!  ux(i), uy(i), uz(i) = velocity moment components, i=1,nx
!
! INPUTS: The input is done via 2- Namelists: inputdat & ions. This is
!         in a file usually named as "hyb1d_yourlabel.dat". See this last file for
!         a description of the input variables.
!
! OUTPUTS: The output is generated in three files: 
!          fields.d10 contains magnetic & electric fields, moment density &
!                    velocities
!          partcls.d11 contains the particle position & velocities 
!                    (i.e.distribution functions)
!          energy.d12  contains kinetic, magnetic, electric and total energies
!          restart.d13 data to restart & continue the simulation (Not yet
!                     implememented)
!
!            
      USE share_data
      REAL :: t,dts, time
      INTEGER :: i, j, it, ist, itime
!
      call cpu_time(t_begin)
      ist=0
      it=0
      t=0.
!
! Initialization
!
      CALL init     
      call energy(it,t)
      call outdat(it,t)
!
! Restarting Simulation 
!
      IF (itrestart /= 0) THEN
!        call restart(t)
      END IF

!
! Skipping Initialization for Restart 
!
!      IF (itrestart /= 0) GO TO 1
!
! Initial Condition
!
      it=0
      itime=it
      t=it*dt
      time=t
!
      call pushx
      call moment(it)
      call rrkfld(it)
!
!***********************************************************************
!     looping part
!***********************************************************************
!
      it1=1
      it2=1
      it3=1
      it4=0
      if (ifield.eq.1)    it1=0
      if (iparticles.eq.1) it2=0
      if (ienergy.eq.1)    it3=0
      if (mod(it,ifilter).eq.it4) call filter2
!
      do it=1,itmax
	  time=it*dt
!
        if (mod(it,ienergy).eq.it3) call energy(it,time)
!
        if (mod(it,ifield).eq.it1) call outdat(it,time)
!
        call pushv
!
        call repx
!
        call pushx
!
        call moment(it)
!
        call rrkfld(it)
!
        if (mod(it,ifilter).eq.it4) call filter2
!
      enddo
      it=itmax+1
!
        if (mod(it,ienergy).eq.it3) call energy(it,time)
!
        if (mod(it,ifield).eq.it1) call outdat(it,time)
!
      close(10)
	  close(11)
	  close(12)
      print *,' *** Program Terminated ***'
	  call cpu_time(t_end)
	  tcpu=t_end-t_begin
	  print *,' Total CPU-Time =',tcpu
	  PAUSE ' *** Press Carriage Return (ENTER) to TERMINATE ***'
      END PROGRAM hyb1d_rrkfft
!
      SUBROUTINE init
      USE share_data
!
      CHARACTER(72) :: version,ptitle
      character(90) :: outfil0, outfil1, outfil2, outfil3, outfil4, outfil5
      REAL :: hold,dth,dthdx,dtdx,dxdt,betan
      INTEGER :: indx,is,n1,n2,n,ix,indh
!
      common /indxy/ indx
      NAMELIST /inputdat/version,xlen,dt,tmax,betaen,wpiwci,game,dus,bnamp,bf0,  &
              itrestart,ifield,iparticles,ienergy,ifilter,eta
      NAMELIST /ions/ptitle,nisp,qi,ai,betain,vdr,rdn,anis
!
      home='./'   !for WSL Linux
      outfil0 = trim(home)//'hyb1drrkfft.dat'
      open(2,file=outfil0)
      read(2, inputdat)
      read(2, ions)
      close(2)
      write(*, inputdat)
      write(*, ions)
!
      outfil1= trim(home)//'avgflds.d9'    
      outfil2= trim(home)//'fields.d10'
      outfil3= trim(home)//'partcls.d11'
      outfil4= trim(home)//'energy.d12'
!      outfil5= trim(home)//'restart.dat'

      open(9,file=trim(home)//'avgflds.d9',form='formatted')
      open(10,file=trim(home)//'fields.d10',form='unformatted')
      open(11,file=trim(home)//'partcls.d11',form='unformatted')
      open(12,file=trim(home)//'energy.d12',form='unformatted')
!	  open(13,file=outfil5,form='unformatted')
!
      if(nisp /= nsp) then
        print *,'<<<<< SIMULATION STOPPED >>>>'
        print *,'---- Particle Array Mismatch (np < npt) -----'
        print *, 'nisp, nsp = ',nisp,nsp
        stop
      endif
!
      indx =1
      lfld =1
      time=0.0
      iseed=-53241
      hold=ran2(iseed)
!
!***********************************************************************
!     initialization
!***********************************************************************
!
!
! Check initial constants for normalization
!
!      itrestart=0
      dx=xlen/float(nx)
      slx=float(nx)*dx
      dth=0.5*dt
      dthdx=dth/dx
      dtdx=dt/dx
      dxdt=1./dtdx
      itmax=INT(tmax/dt)+1
	  print *,'*** tmax & itmax ***',tmax,itmax       
!
! Check CFL-Condition
!
!      IF(dxdt <= cvs) then
!        print *,'<<<<< CFL Condition is Not satisfied >>>>'
!        print *,'----  -----'
!        print *, 'dxdt, c/Wthp = ',dxdt,cvs
!        stop
!      ENDIF
!
!
      betan=0.
      DO is=1,nsp
        nions(is)=np/nsp
        gg(is)=(rdn(is)/qi(is))*float(nx)/nions(is)
        betan=betan+betain(is)
      END DO
      betan=betan+betaen
!
!
!   delta=charge/mass ratio-vector array
!
      n2=0
      DO is=1,nsp
        n1=n2+1
        n2=n2+nions(is)
        qm(is)=qi(is)/ai(is)
        DO n=n1,n2
          qmr(n)=qm(is)
        END DO
      END DO
!
! Initialize Particle position and velocities
!
      call initp
!
! Initialize the E-B Fields
!
      CALL initf
!
! Initialize Wave vectors
!
      CALL kvec(dx)
!
!
!      IF(itime == 0) then
        header(1)=dt
        header(2)=dx
        header(3)=float(nx)
        header(4)=float(itmax)
        header(5)=float(lfld)
        header(6)=float(nsp)
        header(7)=float(np)
        header(8)=float(npg)
        header(9)=betaen
        header(10)=wpiwci
        header(11)=bf0(1)
        header(12)=bf0(2)
        header(13)=bf0(3)
        header(14)=float(ifield)
        header(15)=float(iparticles)
        header(16)=float(ienergy)
        header(17)=float(ifilter)
        header(18)=float(itrestart)
!
        DO is=1,nsp
          indh=ihparm+8*(is-1)
          header(indh+1)=rdn(is)
          header(indh+2)=vdr(is)
          header(indh+3)=betain(is)
          header(indh+4)=anis(is)
          header(indh+5)=qi(is)
          header(indh+6)=ai(is)
          header(indh+7)=float(nions(is))
          header(indh+8)=gg(is)
        END DO
        write(10) header
        write(11) header
        write(12) header
!      ENDIF
!
!
      print *,' *** Pause To Check INITIAL Parameters ***'     
      print *,'****** nx, nsp, npg *****',nx, nsp, npg
      print *,'****** Lx, dx, dt *****',slx, dx, dt
      PAUSE ' *** Press Carriage Return (ENTER) to Continue ***'

!
      END SUBROUTINE init
!
      SUBROUTINE initf
      USE share_data
!
      REAL :: sumy,sumz,ran2
	  INTEGER :: ix
!
!***********************************************************************
!
!      Initialize field with  random noise only (abs(bnamp))
!
!***********************************************************************
!
      ex=0.
      ey=0.
      ez=0.
!
      bx=bf0(1)
      DO ix=1,nx
        by(ix)=bnamp*(ran2(iseed)-0.5)*2.+bf0(2)
        bz(ix)=bnamp*(ran2(iseed)-0.5)*2.+bf0(3)
      END DO
!
      sumy=SUM(by(1:nx)-bf0(2))/float(nx)
      sumz=SUM(bz(1:nx)-bf0(3))/float(nx)
!
      DO ix=1,nx
        by(ix)=by(ix)-sumy
        bz(ix)=bz(ix)-sumz
      END DO
!
!  Periodic boundary conditions
!
      by(  0)=by(nx)
      by(nxp1)=by( 1)
      bz(  0)=bz(nx)
      bz(nxp1)=bz( 1)
      END SUBROUTINE initf
!
      SUBROUTINE initp
      USE share_data
!
       REAL, DIMENSION(np) :: smm,smo,smp
       INTEGER, DIMENSION(np) :: imx,iox,ipx
       REAL, DIMENSION(nx,nsp) :: dnsp0,dnsp1
       REAL, DIMENSION(nx) :: denp1,tex,pe
       REAL :: pi,twopi,xl,rangle,cosa,sina,sumd,ran2,gasdev,gam,hold
       REAL :: vdrpal,vdrper
       INTEGER :: n,n1,n2,ix,is
       common /igridx/imx,iox,ipx
       common /shpfun/smm,smo,smp
!
      common /momnt1/dnsp0,dnsp1,denp1
!
      pi=4.*atan(1.)
      twopi=2.*pi
	  xl=float(nx)*dx
	  rangle=atan2(bf0(2),bf0(1))
	  cosa=cos(rangle)
	  sina=sin(rangle)
!
!
!*****************************************************************
!      subroutine meshx will initialize the mesh and assign the
!      indices to the mesh-grids including the periodic boundary
!      conditions on the ghost-grids
!*****************************************************************
!
      CALL meshx(dx)
!
!
!***********************************************************************
!      initialize particle position and velocities
!           shifted maxwell distribution with thermal spread
!***********************************************************************
!
      n2=0
      DO is=1,nsp
        n1=n2+1
        n2=n2+nions(is)
        DO n=n1,n2
           xp0(n)=xl*ran2(iseed)
        END DO
      END DO
!
!
! //////////////////// shape fun ///////////////////////////////////////
!
      call shapefun(dx,xp0)
!
      n2=0
      DO is=1,nsp
        n1=n2+1
        n2=n2+nions(is)
        vthpal(is)=SQRT(betain(is)/ai(is))
        vthper(is)=SQRT(betain(is)*anis(is)/ai(is))
        vdrpal= vdr(is)
        vdrper=0.0     
        DO n=n1,n2
!        vmag=sqrt(-alog(1.0-0.9999999*ran2(iseed)))
!		th=2.*pi*ran2(iseed)
!		vxt=vthpal(is)*vmag*cos(th)+vdrpal
!        vmag=sqrt(-alog(1.0-0.9999999*ran2(iseed)))
!		th=2.*pi*ran2(iseed)
!		vyt=vthper(is)*vmag*cos(th)+vdrper
!		vzt=vthper(is)*vmag*sin(th)+vdrper
!		vxp(n)=vxt*cosa-vyt*sina
!		vyp(n)=vyt*cosa+vxt*sina
!		vzp(n)=vzt
!
          vmag=gasdev(iseed)
          vxp(n)=vthpal(is)*vmag+vdrpal
          vmag=gasdev(iseed)
          vyp(n)=vthper(is)*vmag+vdrper
          vmag=gasdev(iseed)
          vzp(n)=vthper(is)*vmag+vdrper
!
!         vxt=vthpal(is)*vmag*cos(th)+vdrpal
!         vmag=gasdev(iseed)
!         th=2.*pi*ran2(iseed)
!         vyt=vthper(is)*vmag*cos(th)+vdrper
!         vzt=vthper(is)*vmag*sin(th)+vdrper
!         vxp(n)=vxt*cosa-vyt*sina
!         vyp(n)=vyt*cosa+vxt*sina
!         vzp(n)=vzt
        END DO
      END DO
!
!***********************************************************************
!       Take the moment of particles over each specie for density
!       dnsp1(ix.is) : real density of is-specie (ion1, ion2,...etc
!                      at the position whose mesh number is ix
!***********************************************************************
!
      DO is=1,nsp
        DO ix=1,nx
          dnsp1(ix,is)=0.0
        END DO
      END DO
!
      n2=0
      DO is=1,nsp
        n1=n2+1
        n2=n2+nions(is)
        DO n=n1,n2
          dnsp1(imx(n),is)=dnsp1(imx(n),is) + smm(n)
          dnsp1(iox(n),is)=dnsp1(iox(n),is) + smo(n)
          dnsp1(ipx(n),is)=dnsp1(ipx(n),is) + smp(n)
        END DO
      END DO
!
      DO is=1,nsp
        DO ix=1,nx
          dnsp1(ix,is)=dnsp1(ix,is)*gg(is)
        END DO
      END DO
!
      gam=game-1.0
      DO ix=1,nx
        sumd=0.
        DO is=1,nsp
          sumd=sumd + dnsp1(ix,is)*ai(is)
        END DO
        denp1(ix)=sumd
        tex(ix)=betaen*denp1(ix)**gam
        pe(ix)=denp1(ix)*tex(ix)
      END DO
      END SUBROUTINE initp
!
      SUBROUTINE meshx(dx)
      USE share_dim
!
      REAL, DIMENSION(0:nxp1) :: xi
      INTEGER, DIMENSION(0:nxp1) :: igmx,igox,igpx
      REAL :: dx
      INTEGER :: ix
      common /xmesh/xi,igmx,igox,igpx
!
!  xi= the mesh location
!  igox= the mesh index number
!  igmx= the previous mesh index
!  igpx= the next mesh index
!
      DO ix=0,nxp1
        xi(ix)=dx*float(ix)
        igmx(ix)=ix-1
        igox(ix)=ix
        igpx(ix)=ix+1
      END DO
!
!  Periodic boundary conditions
!
      igmx(0)=nx-1
      igmx(1)=nx
      igox(0)=nx
      igox(nxp1)=1
      igpx(nx)=1
      igpx(nxp1)=2
!
      END SUBROUTINE meshx
!
      SUBROUTINE shapefun(dx,xp)
        USE share_dim
!
	REAL, DIMENSION(np) :: xp,smm,smo,smp
	INTEGER, DIMENSION(np) :: imx,iox,ipx
        REAL, DIMENSION(0:nxp1) :: xi
	INTEGER, DIMENSION(0:nxp1) :: igmx,igox,igpx
        REAL :: dx,gdx,dlx
        INTEGER :: n,ix
        common /xmesh/xi,igmx,igox,igpx
        common /igridx/imx,iox,ipx
        common /shpfun/smm,smo,smp
!
!  xp= particle position
!  ix= nearest mesh number
!  dlx= distance between particle and the nearest mesh
!
      gdx=1.0/dx
      DO n=1,np
        ix=xp(n)*gdx+0.5
        dlx=xp(n)*gdx-float(ix)
        imx(n)=igmx(ix)
        iox(n)=igox(ix)
        ipx(n)=igpx(ix)
        smm(n)=(0.5*dlx**2 - 0.5*dlx + 0.125)
        smo(n)=(-dlx**2 + 0.750)
        smp(n)=(0.5*dlx**2 + 0.5*dlx + 0.125)
      END DO
      END SUBROUTINE shapefun
!
      SUBROUTINE pushx
      USE share_data
!
      REAL, DIMENSION(0:nxp1) :: xi
      INTEGER, DIMENSION(0:nxp1) :: igmx,igox,igpx
      REAL :: xmax,xmin
      INTEGER :: n
      common /xmesh/xi,igmx,igox,igpx
!
!***********************************************************************
!     Time integration of particle position
!
!              dx
!              -- = vx
!              dt
!
!***********************************************************************
!     xmin,xmax : boundary
!***********************************************************************
!
      xmax=xi(nx)
      xmin=0.
      DO n=1,np
        xp1(n)=xp0(n)+dt*vxp(n)
      END DO
!
!***********************************************************************
!     Periodic boundary condition
!***********************************************************************
!
      DO n=1,np
        IF(xp1(n) > xmax) then
           xp1(n)=xp1(n)-xmax
        ELSE IF(xp1(n) < xmin) then
           xp1(n)=xp1(n)+xmax
        ENDIF
      END DO
      END SUBROUTINE pushx
!
      SUBROUTINE repx
      USE share_data
!
      INTEGER :: n
!
      DO n=1,np
        xp0(n)=xp1(n)
      END DO
      END SUBROUTINE repx
!
      SUBROUTINE moment(it)
      USE share_data
!
	  REAL, DIMENSION(np) :: smm,smo,smp
	  INTEGER, DIMENSION(np) :: imx,iox,ipx
      REAL, DIMENSION(0:nxp1) :: xi,denp1,denh,uxh,uyh,uzh
      INTEGER, DIMENSION(0:nxp1) :: igmx,igox,igpx
      REAL, DIMENSION(nx,nsp) :: dnsp0,dnsp1,dnsh,uxsh,uysh,uzsh
      REAL :: sumdh,sumd1,sumvx,sumvy,sumvz
      REAL :: dnhavg,dnp1avg,uxhavg,uyhavg,uzhavg
      INTEGER :: i,is,ix,n,n1,n2,it
!
      common /xmesh/xi,igmx,igox,igpx
      common /igridx/imx,iox,ipx
      common /shpfun/smm,smo,smp
      common /momnt1/dnsp0,dnsp1,denp1
      common /momnt2/dnsh,uxsh,uysh,uzsh
      common /momnth/denh,uxh,uyh,uzh
!
!***********************************************************************
!     Take the moment of particles - 2/2 method (2nd order )
!
!     ix: nearest mesh number
!     dl: length between particle and the nearest mesh
!     dnsp0(ix,is) : real density of is-specie(ion1 or ion2)
!                    at the position whose mesh number is ix
!     uxsh(ix,is): dnsp0(ix,is)*(velocity of is-specie
!                  at the position whose mesh number is ix)
!
!                <charge neutrality>
!     dens(elec)=dens(ion1)*qi(ion1)+qi(ion2)*dens(ion2)+...
!
!     uxh: mean velocity
!
!          dens(ion1)*uxsh(ion1)*qi(ion1)+qi(ion2)*dens(ion2)*uxsh(ion2)+...
!     uxh=-----------------------------------------------------------
!                                dens(elec)
!***********************************************************************
!
      DO is=1,nsp
        DO ix=1,nx
          dnsp0(ix,is)=dnsp1(ix,is)
          dnsp1(ix,is)=0.
          uxsh (ix,is)=0.
          uysh (ix,is)=0.
          uzsh (ix,is)=0.
        END DO
      END DO
!
      n2=0
      DO is=1,nsp
           n1=n2+1
           n2=n2+nions(is)
!
        DO n=n1,n2
          uxsh(imx(n),is)=uxsh(imx(n),is)+ smm(n)*vxp(n)
          uxsh(iox(n),is)=uxsh(iox(n),is)+ smo(n)*vxp(n)
          uxsh(ipx(n),is)=uxsh(ipx(n),is)+ smp(n)*vxp(n)
!
          uysh(imx(n),is)=uysh(imx(n),is)+ smm(n)*vyp(n)
          uysh(iox(n),is)=uysh(iox(n),is)+ smo(n)*vyp(n)
          uysh(ipx(n),is)=uysh(ipx(n),is)+ smp(n)*vyp(n)
!
          uzsh(imx(n),is)=uzsh(imx(n),is)+ smm(n)*vzp(n)
          uzsh(iox(n),is)=uzsh(iox(n),is)+ smo(n)*vzp(n)
          uzsh(ipx(n),is)=uzsh(ipx(n),is)+ smp(n)*vzp(n)
        END DO
      END DO
!
! /////////////////// Shape fun ////////////////////////////////////////
!
      CALL shapefun(dx,xp1)
!
      n2=0
      DO is=1,nsp
        n1=n2+1
        n2=n2+nions(is)
!
        DO n=n1,n2
          dnsp1(imx(n),is)=dnsp1(imx(n),is) + smm(n)
          dnsp1(iox(n),is)=dnsp1(iox(n),is) + smo(n)
          dnsp1(ipx(n),is)=dnsp1(ipx(n),is) + smp(n)
!
          uxsh(imx(n),is)=uxsh(imx(n),is)+ smm(n)*vxp(n)
          uxsh(iox(n),is)=uxsh(iox(n),is)+ smo(n)*vxp(n)
          uxsh(ipx(n),is)=uxsh(ipx(n),is)+ smp(n)*vxp(n)
!
          uysh(imx(n),is)=uysh(imx(n),is)+ smm(n)*vyp(n)
          uysh(iox(n),is)=uysh(iox(n),is)+ smo(n)*vyp(n)
          uysh(ipx(n),is)=uysh(ipx(n),is)+ smp(n)*vyp(n)
!
          uzsh(imx(n),is)=uzsh(imx(n),is)+ smm(n)*vzp(n)
          uzsh(iox(n),is)=uzsh(iox(n),is)+ smo(n)*vzp(n)
          uzsh(ipx(n),is)=uzsh(ipx(n),is)+ smp(n)*vzp(n)
        END DO
      END DO
!
      DO is=1,nsp
        DO ix=1,nx
          uxsh (ix,is)=uxsh (ix,is)*gg(is)*0.5
          uysh (ix,is)=uysh (ix,is)*gg(is)*0.5
          uzsh (ix,is)=uzsh (ix,is)*gg(is)*0.5
          dnsp1(ix,is)=dnsp1(ix,is)*gg(is)
        END DO
      END DO
!
      DO is=1,nsp
        DO ix=1,nx
          dnsh(ix,is)=(dnsp1(ix,is)+dnsp0(ix,is))*0.5
        END DO
      END DO
!
!
      DO ix=1,nx
        sumdh=0.0
        sumd1=0.0
        sumvx=0.0
        sumvy=0.0
        sumvz=0.0
        DO is=1,nsp
          sumdh=sumdh + dnsh(ix,is)*qi(is)
          sumd1=sumd1 + dnsp1(ix,is)*qi(is)
          sumvx=sumvx+uxsh(ix,is)*qi(is)
          sumvy=sumvy+uysh(ix,is)*qi(is)
          sumvz=sumvz+uzsh(ix,is)*qi(is)
        END DO
        denh(ix)=sumdh
        denp1(ix)=sumd1
        if(sumdh > 0.) then 
          uxh(ix)=sumvx/denh(ix)
          uyh(ix)=sumvy/denh(ix)
          uzh(ix)=sumvz/denh(ix)
        else
          uxh(ix)=0.
          uyh(ix)=0.
          uzh(ix)=0.
        endif
      END DO
!
      CALL filter(denh)
      CALL filter(denp1)
      CALL filter(uxh)
      CALL filter(uyh)
      CALL filter(uzh)
!
!***********************************************************************
!        Periodic Boundary Conditions
!***********************************************************************
!
      denh(0)   =denh(nx)
      denp1(0)  =denp1(nx)
      denh(nxp1) =denh(1)
      denp1(nxp1)=denp1(1)
      uxh(0)=uxh(nx)
      uyh(0)=uyh(nx)
      uzh(0)=uzh(nx)
      uxh(nxp1)=uxh(1)
      uyh(nxp1)=uyh(1)
      uzh(nxp1)=uzh(1)
!
!
!  Write Average Moments
!
!      IF((mod(it,ifield) == it1)) then
!      dnhavg=SUM(denh(1:nx))/float(nx)
!      dnp1avg=SUM(denp1(1:nx))/float(nx)
!      uxhavg=SUM(uxh(1:nx))/float(nx)
!      uyhavg=SUM(uyh(1:nx))/float(nx)
!      uzhavg=SUM(uzh(1:nx))/float(nx)
!      print *,' Avg Moments =',it,dnhavg,dnp1avg,uxhavg,uyhavg,uzhavg
! 19   format(I6,3x,4(e11.3,3x))
!      ENDIF
!
      END SUBROUTINE moment
!
      SUBROUTINE rrkfld(it)
      USE share_data
!
      REAL, DIMENSION(nx) :: gy1,gz1,gy2,gz2,byd,bzd
      REAL :: ddt,g11,g12,g22,gg22,bxd,etan
	  REAL :: byavg,bzavg,bperpavg
      INTEGER :: l,ix,it
!
      common /params/etan
!
      etan=eta
      IF (lfld <= 0) RETURN
      ddt=dt/float(lfld)
!
!***********************************************************************
!     Time integration of magnetic field
!            RRK method
!             dB
!             -- = - rot(E)
!             dt
!
!                                  betae    grad(ne)     (rot(B) x B)
!             E = - v(mean) x B - ------ * ---------  + -------------
!                                    2        ne              ne
!
!***********************************************************************
!
!//////////// RRK looping //////////////////////////////////////////////
!
      bxd=bx
      DO l=1,lfld
!
!***********************************************************************
!     RRK advancing : first step
!***********************************************************************
!
        DO ix=1,nx
          byd(ix)=by(ix)
          bzd(ix)=bz(ix)
        END DO
!   ***************************
        CALL rrksub(bxd,byd,bzd,gy1,gz1)
!   ***************************
        DO ix=1,nx
          byd(ix)=by(ix)+ddt*0.5*gy1(ix)
          bzd(ix)=bz(ix)+ddt*0.5*gz1(ix)
        END DO
!
!**********************************************************
!     RRK advancing : second step
!**********************************************************
!
!     ***************************
        CALL rrksub(bxd,byd,bzd,gy2,gz2)
!     ***************************
        g11=0.
        g12=0.
        g22=0.
        DO ix=1,nx
          gy2(ix)=2.*gy1(ix)-gy2(ix)
          gz2(ix)=2.*gz1(ix)-gz2(ix)
!
          g11=g11+gy1(ix)**2+gz1(ix)**2
          g22=g22+gy2(ix)**2+gz2(ix)**2
          g12=g12+gy1(ix)*gy2(ix)+gz1(ix)*gz2(ix)
        END DO
!
        gg22=1./g22 *ddt
!
        DO ix=1,nx
          gyy=gg22*(2.*g12*gy1(ix)-g11*gy2(ix))
          gzz=gg22*(2.*g12*gz1(ix)-g11*gz2(ix))
          by(ix)=by(ix)+gyy
          bz(ix)=bz(ix)+gzz
        END DO
      END DO
!
!***************************************************
!       periodic boundary condition
!***************************************************
!
      by(0  )=by(nx)
      bz(0  )=bz(nx)
      by(nxp1)=by(1 )
      bz(nxp1)=bz(1 )
!
!
!  Write Average Magnetic Field
!
!      IF((mod(it,ifield) == it1)) then
!      byavg=SUM(by(1:nx)-bf0(2))/float(nx)
!      bzavg=SUM(bz(1:nx)-bf0(3))/float(nx)
!      bperpavg=SUM(SQRT((by(1:nx)-bf0(2))**2+(bz(1:nx)-bf0(3))**2))/float(nx)
!      write(9,19) it,byavg,bzavg,bperpavg
! 19   format(I6,3x,3(e11.3,3x))
!      ENDIF
      END SUBROUTINE rrkfld
!
      SUBROUTINE rrksub(bxd,byd,bzd,cury,curz)
      USE share_dim
!
      REAL, DIMENSION(0:nxp1) :: denh,uxh,uyh,uzh
      REAL, DIMENSION(nx) :: cury,curz,byd,bzd
      REAL, DIMENSION(nx) :: wrk1,wrk2,wrk3
      REAL :: bxd,etan
      INTEGER :: ix
      common /momnth/denh,uxh,uyh,uzh
      common /params/etan
!
!***********************************************************
!     calculation of -rot(E)
!***********************************************************
!
      DO ix=1,nx
        wrk2(ix)=byd(ix)
        wrk3(ix)=bzd(ix)
      END DO
!      CALL curl1(wrk1,wrk2,wrk3,cury,curz)
      CALL curl2(wrk1,wrk2,wrk3,cury,curz)
!
!
!***********************************************************
!      0
!     wrk2 : -E
!     wrk3
!***********************************************************
!
      DO ix=1,nx
        if(denh(ix) > 0.) then
          wrk2(ix)=uzh(ix)*bxd-uxh(ix)*bzd(ix)   &
                 -( curz(ix)*bxd)/denh(ix)-etan*cury(ix)
          wrk3(ix)=uxh(ix)*byd(ix)-uyh(ix)*bxd   &
                 -(-cury(ix)*bxd)/denh(ix)-etan*curz(ix)
        else
          wrk2(ix) =0.
          wrk3(ix) =0.
        endif
      END DO
!
!      CALL curl1(wrk1,wrk2,wrk3,cury,curz)
      CALL curl2(wrk1,wrk2,wrk3,cury,curz)
      END SUBROUTINE rrksub
!
      SUBROUTINE pushv
      USE share_data
!
      REAL, DIMENSION(np) :: smm,smo,smp
      INTEGER, DIMENSION(np) :: imx,iox,ipx
      REAL, DIMENSION(0:nxp1) :: denp1,denh,uxh,uyh,uzh
      REAL, DIMENSION(nx,nsp) :: dnsp0,dnsp1,dnsh,uxsh,uysh,uzsh
      REAL, DIMENSION(nx) :: cury,curz,grdnx,wrk1,wrk2,wrk3
      REAL, DIMENSION(nx,3) :: fnp1,un
      REAL, DIMENSION(nx,3,nsp) :: run
      REAL, DIMENSION(np,3) :: g1,g2,vpp
      REAL :: gdx,g11,g12,g22
      INTEGER :: is,ix,n,n1,n2,k
!
      common /igridx/imx,iox,ipx
      common /shpfun/smm,smo,smp
      common /momnt1/dnsp0,dnsp1,denp1
      common /momnt2/dnsh,uxsh,uysh,uzsh
      common /momnth/denh,uxh,uyh,uzh
!
!***********************************************************************
!     Time integration of particle velocity
!            RRK method
!
!             dv
!             -- = qmr * ( E + v * B )
!             dt
!                        (qmr : ion1 = 1,ion2 = delta)
!
!              E = - U(mean) * B                     <=   unknown term
!
!                      grad(P)      (rot(B) * B
!                  - ------------ + -----------      <=   known terms
!                    2.dens(elec)    dens(elec)
!***********************************************************************
!
      gdx =1./dx
!
      DO is=1,nsp
        DO k=1,3
          DO ix=1,nx
            run(ix,k,is)=0.
          ENDDO   !Close ix-loop
        ENDDO   !Close k-loop
      ENDDO   !Close is-loop
!
!**********************************************************************
!
!    grdnx = gradx(dens(elec))
!
!     0    }
!          }
!    cury  }  = rot(B)
!          }
!    curz  }
!
!**********************************************************************
!
      DO ix=1,nx
        wrk1(ix)=denp1(ix)
      ENDDO
!      call grad1(wrk1,grdnx)
      call grad2(wrk1,grdnx)
!
      DO ix=1,nx
        wrk2(ix)=by(ix)
        wrk3(ix)=bz(ix)
      ENDDO
!      call curl1(wrk1,wrk2,wrk3,cury,curz)
      call curl2(wrk1,wrk2,wrk3,cury,curz)
!
!**********************************************************************
!    fnp1(i,1)
!    fnp1(i,2): known term of acceleration
!    fnp1(i,3)
!**********************************************************************
!
      DO ix=1,nx
        if(denp1(ix) > 0.) then
           gdnnp1=1./denp1(ix)
        else
           gdnnp1=0.
        endif
!
        fnp1(ix,1)= (-betaen* grdnx(ix)*0.5   &
                    + cury(ix)*bz(ix)- curz(ix)*by(ix))* gdnnp1
        fnp1(ix,2)= ( curz(ix)*bx )* gdnnp1+eta*cury(ix)
        fnp1(ix,3)= (-cury(ix)*bx )* gdnnp1+eta*curz(ix)
      ENDDO   !Close ix-loop
!
!**********************************************************************
!    RRK advancing : first step
!**********************************************************************
!
!
      DO k=1,3
        DO n=1,np
          g2(n,k)= (fnp1(imx(n),k)*smm(n)   &
                   + fnp1(iox(n),k)*smo(n)   &
                   + fnp1(ipx(n),k)*smp(n))*qmr(n)
        ENDDO   !Close n-loop
      ENDDO    !Close k-loop
!
      DO n=1,np
        g1(n,1)=g2(n,1)+qmr(n)*   &
              (((vyp(n)-uyh(imx(n)))*bz(imx(n))  &
                 -(vzp(n)-uzh(imx(n)))*by(imx(n)))*smm(n)  &
              +((vyp(n)-uyh(iox(n)))*bz(iox(n))  &
                 -(vzp(n)-uzh(iox(n)))*by(iox(n)))*smo(n)  &
              +((vyp(n)-uyh(ipx(n)))*bz(ipx(n))  &
                 -(vzp(n)-uzh(ipx(n)))*by(ipx(n)))*smp(n))
!
        g1(n,2)=g2(n,2)+qmr(n)*   &
              (((vzp(n)-uzh(imx(n)))*bx  &
                 -(vxp(n)-uxh(imx(n)))*bz(imx(n)))*smm(n)  &
              +((vzp(n)-uzh(iox(n)))*bx  &
                 -(vxp(n)-uxh(iox(n)))*bz(iox(n)))*smo(n)  &
              +((vzp(n)-uzh(ipx(n)))*bx  &
                 -(vxp(n)-uxh(ipx(n)))*bz(ipx(n)))*smp(n))
!
        g1(n,3)=g2(n,3)+qmr(n)*   &
              (((vxp(n)-uxh(imx(n)))*by(imx(n))  &
                 -(vyp(n)-uyh(imx(n)))*bx)*smm(n)  &
              +((vxp(n)-uxh(iox(n)))*by(iox(n))  &
                 -(vyp(n)-uyh(iox(n)))*bx)*smo(n)  &
              +((vxp(n)-uxh(ipx(n)))*by(ipx(n))  &
                 -(vyp(n)-uyh(ipx(n)))*bx)*smp(n))
      ENDDO   !Close n-loop
!
!
      DO n=1,np
        vpp(n,1)=vxp(n)+0.5*dt*g1(n,1)
        vpp(n,2)=vyp(n)+0.5*dt*g1(n,2)
        vpp(n,3)=vzp(n)+0.5*dt*g1(n,3)
      ENDDO   !Close n-loop
!
      n2=0
      DO is=1,nsp
        n1=n2+1
        n2=n2+nions(is)
        q =qi(is)
        DO k=1,3
          DO n=n1,n2
            run(imx(n),k,is)=run(imx(n),k,is)+q*vpp(n,k)*smm(n)
            run(iox(n),k,is)=run(iox(n),k,is)+q*vpp(n,k)*smo(n)
            run(ipx(n),k,is)=run(ipx(n),k,is)+q*vpp(n,k)*smp(n)
          ENDDO   !Close n-loop
        ENDDO   !Close k-loop
      ENDDO   !Close is-loop
!
      DO ix=1,nx
        if(denp1(ix) > 0.) then
          gdnnp1=1./denp1(ix)
        else
          gdnnp1=0.
        endif
!        gdnnp1=1./denp1(ix)
        rvsumx=0.
        rvsumy=0.
        rvsumz=0.
        DO is=1,nsp
          rvsumx=rvsumx+run(ix,1,is)*gg(is)
          rvsumy=rvsumy+run(ix,2,is)*gg(is)
          rvsumz=rvsumz+run(ix,3,is)*gg(is)
        ENDDO   !Close is-loop
        un(ix,1)=rvsumx*gdnnp1
        un(ix,2)=rvsumy*gdnnp1
        un(ix,3)=rvsumz*gdnnp1
      ENDDO   !Close ix-loop
!
      DO ix=1,nx
        ex(ix)=-un(ix,2)*bz(ix)+un(ix,3)*by(ix)+fnp1(ix,1)
        ey(ix)=-un(ix,3)*bx+un(ix,1)*bz(ix)+fnp1(ix,2)
        ez(ix)=-un(ix,1)*by(ix)+un(ix,2)*bx+fnp1(ix,3)
      ENDDO   !Close ix-loop
!
!**********************************************************************
!    RRK advancing : second step
!           vn(i,1): (x-component of mean velocity)
!**********************************************************************
!
      DO n=1,np
        g2(n,1)=g2(n,1)+qmr(n)*   &
              (((vpp(n,2)-un(imx(n),2))*bz(imx(n))  &
                 -(vpp(n,3)-un(imx(n),3))*by(imx(n)))*smm(n)  &
              +((vpp(n,2)-un(iox(n),2))*bz(iox(n))  &
                 -(vpp(n,3)-un(iox(n),3))*by(iox(n)))*smo(n)  &
              +((vpp(n,2)-un(ipx(n),2))*bz(ipx(n))  &
                 -(vpp(n,3)-un(ipx(n),3))*by(ipx(n)))*smp(n))
!
        g2(n,2)=g2(n,2)+qmr(n)*   &
              (((vpp(n,3)-un(imx(n),3))*bx  &
                 -(vpp(n,1)-un(imx(n),1))*bz(imx(n)))*smm(n)  &
              +((vpp(n,3)-un(iox(n),3))*bx  &
                 -(vpp(n,1)-un(iox(n),1))*bz(iox(n)))*smo(n)  &
              +((vpp(n,3)-un(ipx(n),3))*bx  &
                 -(vpp(n,1)-un(ipx(n),1))*bz(ipx(n)))*smp(n))
!
        g2(n,3)=g2(n,3)+qmr(n)*   &
              (((vpp(n,1)-un(imx(n),1))*by(imx(n))  &
                 -(vpp(n,2)-un(imx(n),2))*bx)*smm(n)  &
              +((vpp(n,1)-un(iox(n),1))*by(iox(n))  &
                 -(vpp(n,2)-un(iox(n),2))*bx)*smo(n)  &
              +((vpp(n,1)-un(ipx(n),1))*by(ipx(n))  &
                 -(vpp(n,2)-un(ipx(n),2))*bx)*smp(n))
      ENDDO   !Close n-loop
!
!
       g11=0.
       g22=0.
       g12=0.
       DO k=1,3
         DO n=1,np
           g2(n,k)=2.*g1(n,k)-g2(n,k)
           g11=g11+g1(n,k)*g1(n,k)
           g22=g22+g2(n,k)*g2(n,k)
           g12=g12+g1(n,k)*g2(n,k)
         ENDDO   !Close n-loop
       ENDDO   !Close k-loop
!
       g22=dt/g22
       DO n=1,np
         vxp(n)=vxp(n)+g22*(2.*g12*g1(n,1)-g11*g2(n,1))
         vyp(n)=vyp(n)+g22*(2.*g12*g1(n,2)-g11*g2(n,2))
         vzp(n)=vzp(n)+g22*(2.*g12*g1(n,3)-g11*g2(n,3))
       ENDDO   !Close n-loop
       END SUBROUTINE pushv
!
      SUBROUTINE curl1(x,y,z,ry,rz)
      USE share_dim
!
      REAL, DIMENSION(nx) :: x,y,z,rx,ry,rz
	  REAL, DIMENSION(nx) :: ek,attenu
      REAL, DIMENSION(nx) :: axr,bxr,axi,bxi,c1r,c1i,c2r,c2i
      REAL, DIMENSION(nxh*lpx,2) :: txr,txi
      INTEGER, DIMENSION(nxh*lpx) :: listx
      INTEGER :: i
!
      common/k_vec/ek
      common/attenuation/attenu
!
!  y  Components of Curl
!
      do i=1,nx
        axr(i)=z(i)
        axi(i)=0.0
      enddo
!
!  Forward FFT (transform into Fourier k-space)
!
      call fftx1(axr,axi,bxr,bxi,txr,txi,listx,nxh,1,lpx)
      do i=1,nx
        axr(i)=-ek(i)*bxi(i)
        axi(i)= ek(i)*bxr(i)
      enddo
!
!  Inverse Fourier transform (transform into configuration space)
!
      call fftx1(axr,axi,bxr,bxi,txr,txi,listx,nxh,2,lpx)
      do i=1,nx
        ry(i)=-bxr(i)
      enddo
!
!
!  z Components of Curl
!
      do i=1,nx
        axr(i)=y(i)
        axi(i)=0.0
      enddo
!
!  Forward FFT (transform into Fourier k-space)
!
      call fftx1(axr,axi,bxr,bxi,txr,txi,listx,nxh,1,lpx)
      do i=1,nx
        axr(i)=-ek(i)*bxi(i)
        axi(i)= ek(i)*bxr(i)
      enddo
!
!  Inverse Fourier transform (transform into configuration space)
!
      call fftx1(axr,axi,bxr,bxi,txr,txi,listx,nxh,2,lpx)
      do i=1,nx
        rz(i)=bxr(i)
      enddo
      return
!
! **********************************************************************
      entry grad1(x,rx)
!
! ---  x-component
!
      do i=1,nx
        axr(i)=x(i)
        axi(i)=0.0
      enddo
!
!  forward fft (transform into fourier k-space)
!
      call fftx1(axr,axi,bxr,bxi,txr,txi,listx,nxh,1,lpx)
      do i=1,nx
        axr(i)=-ek(i)*bxi(i)
        axi(i)= ek(i)*bxr(i)
      enddo
!
!  inverse fourier transform (transform into configuration space)
!
      call fftx1(axr,axi,bxr,bxi,txr,txi,listx,nxh,2,lpx)
      do i=1,nx
        rx(i)=bxr(i)
      enddo
      END SUBROUTINE curl1
!
      SUBROUTINE curl2(x,y,z,ry,rz)
      USE share_dim
!
      REAL, DIMENSION(nx) :: x,y,z,rx,ry,rz
	  REAL, DIMENSION(nx) :: ek,attenu
      REAL, DIMENSION(nx) :: axr,bxr,axi,bxi,c1r,c1i,c2r,c2i
      REAL, DIMENSION(nxh*lpx,2) :: txr,txi
      INTEGER, DIMENSION(nxh*lpx) :: listx
      INTEGER :: i
!
      common/k_vec/ek
      common/attenuation/attenu
!
!  y & z Components of Curl
!
      do i=1,nx
        axr(i)=z(i)
        axi(i)=y(i)
      enddo
!
!  Forward FFT (transform into Fourier k-space)
!
      call fftx1(axr,axi,bxr,bxi,txr,txi,listx,nxh,1,lpx)
      call tran1(bxr,bxi,c1r,c1i,c2r,c2i,nx)
      do i=1,nx
        axr(i)=-ek(i)*c1i(i)
        axi(i)= ek(i)*c1r(i)
        bxr(i)=-ek(i)*c2i(i)
        bxi(i)= ek(i)*c2r(i)
      enddo
!
!  Inverse Fourier transform (transform into configuration space)
!
      call tran2(c1r,c1i,axr,axi,bxr,bxi,nx)
      call fftx1(c1r,c1i,bxr,bxi,txr,txi,listx,nxh,2,lpx)
      do i=1,nx
        ry(i)=-bxr(i)
        rz(i)=+bxi(i)
      enddo
!
      return
!
! **********************************************************************
      entry grad2(x,rx)
!
! ---  x-component
!
      do i=1,nx
        axr(i)=x(i)
        axi(i)=0.0
      enddo
!
!  forward fft (transform into fourier k-space)
!
      call fftx1(axr,axi,bxr,bxi,txr,txi,listx,nxh,1,lpx)
      do i=1,nx
        axr(i)=-ek(i)*bxi(i)
        axi(i)= ek(i)*bxr(i)
      enddo
!
!  inverse fourier transform (transform into configuration space)
!
      call fftx1(axr,axi,bxr,bxi,txr,txi,listx,nxh,2,lpx)
      do i=1,nx
        rx(i)=bxr(i)
      enddo
      END SUBROUTINE curl2
!
      SUBROUTINE kvec(dx)
      USE share_dim
!
      REAL, DIMENSION(nx) :: ek,attenu
      REAL :: pi,xl,hal,hk
      INTEGER :: m,ia,ib,nx3,k
      common/k_vec/ek
      common/attenuation/attenu
!
!  xl= maximum length of the box in units of ion-inertial length
!  ek= wavevector in units of the ion-inertial length
!
      pi=4.0*atan(1.0)
      xl=float(nx)*dx
      hal=2.0*pi/xl
      ek(1)=0.0
      attenu(1)=1.0
      do m=1,nxh-1
        ia=m+1
        ib=nx-m+1
        hk=hal*float(m)
        ek(ia)=hk
        ek(ib)=-hk
        attenu(ia)=1.0
        attenu(ib)=1.0
      enddo
      ek(nxh+1)=pi/dx
      attenu(nxh+1)=0.0
      return
!
      entry da_kvec
!
! Dealiasing
!
      nx3=nx/3
      attenu(1)=1.0
      do k=1,nxh-1
        ia=k+1
        ib=nx-k+1
        if(k.le.nx3) then
        attenu(ia)=1.0
        attenu(ib)=1.0
        else
        attenu(ia)=0.0
        attenu(ib)=0.0
        endif
      enddo
      attenu(nxh+1)=0.0
!
      END SUBROUTINE kvec
!
      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
!      USES ran2
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran2
      SAVE iset,gset
      DATA iset/0/
      if (iset == 0) then
1       v1=2.*ran2(idum)-1.
        v2=2.*ran2(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.) goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      END FUNCTION gasdev
!
      FUNCTION ran2(idum)
      INTEGER :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL :: ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,  &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,   &
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER :: idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum <= 0) then
        idum=max(-idum,1)
        idum2=idum
        do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum < 0) idum=idum+IM1
          if (j <= NTAB) iv(j)=idum
        END DO
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum < 0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2 < 0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy < 1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      END FUNCTION ran2
!
      SUBROUTINE energy(it,time)
      USE share_data
!
      REAL, DIMENSION(nsp) :: umx,umy,umz,vx2,vy2,vz2,vkenrgy,tpal,tper,tpery,tperz
      REAL, DIMENSION(nsp) :: wkpal,wkper
      REAL, DIMENSION(3) :: pmom
      REAL, DIMENSION(6) :: enrgys
      REAL :: bfldx,bfldy,bfldz,efldx,efldy,efldz
      REAL :: bfldenrgy,efldenrgy,wktotenrgy,totenrgy
      REAL :: ggs,fac1,fac2,time,ksum,wksum1,wksum2,wkpalenrgy,wkperenrgy
      INTEGER :: n1,n2,n,i,j,ix,is,it
!
      common /tempi/tpal,tper,umx,enrgys
      common /wkenergy/wkpal,wkper
!
!-----------------------------------------------------------------------
!  #   time table of energy & momentum consevation
!-----------------------------------------------------------------------
!
!
! ---- energy (particles)
!
      enrgys(:) = 0.0
      wksum = 0.0
      wksum1 = 0.0
      wksum2 = 0.0
      n2=0       
      DO is=1,nsp
        n1=n2+1
        n2=n2+nions(is)
        fac1=1.0/float(nions(is))   !
        fac2=0.5*ai(is)*gg(is)
!
        umx(is)=0.0
        umy(is)=0.0
        umz(is)=0.0
        DO n=n1,n2
          umx(is)=umx(is)+vxp(n)*fac1
          umy(is)=umy(is)+vyp(n)*fac1
          umz(is)=umz(is)+vzp(n)*fac1
!          umx(is)=fac1*SUM(vxp(n1:n2))
!          umy(is)=fac1*SUM(vyp(n1:n2))
!          umz(is)=fac1*SUM(vzp(n1:n2))
        ENDDO   !Close n-loop
!
        wkpal(is) = 0.0
        wkper(is) = 0.0
        vx2(is)=0.0
        vy2(is)=0.0
        vz2(is)=0.0
        tpal(is)=0.0
        tpery(is)=0.0
        tperz(is)=0.0
        DO n=n1,n2
          vx2(is)=vx2(is)+vxp(n)**2 *fac2
          vy2(is)=vy2(is)+vyp(n)**2 *fac2
          vz2(is)=vz2(is)+vzp(n)**2 *fac2
          tpal(is)=tpal(is)+(vxp(n)-umx(is))**2*fac2
          tpery(is)=tpery(is)+(vyp(n)-umy(is))**2*fac2
          tperz(is)=tperz(is)+(vzp(n)-umz(is))**2*fac2
        ENDDO   !Close n-loop
        tper(is)=0.5*(tpery(is)+tperz(is))
        wkpal(is) = wkpal(is) + vx2(is)
        wkper(is) = wkper(is) + 0.5*(vy2(is)+vz2(is))
!
!          vx2(is)=fac2*SUM(vxp(n1:n2)**2)
!          vy2(is)=fac2*SUM(vyp(n1:n2)**2)
!          vz2(is)=fac2*SUM(vzp(n1:n2)**2)
!          tpal(is)=fac2*SUM((vxp(n1:n2)-umx(is))**2)
!          tpery(is)=fac2*SUM((vyp(n1:n2)-umy(is))**2)
!          tperz(is)=fac2*SUM((vzp(n1:n2)-umz(is))**2)
!
        vkenrgy(is)=(vx2(is)+vy2(is)+vz2(is))
!        vkenrgy(is)=tpal(is)+tper(is)
        wksum1 = wksum1 + wkpal(is)
        wksum2 = wksum2 + wkper(is)
        wksum = wksum + vkenrgy(is)
      ENDDO   !Close is-loop
      wkpalenrgy = wksum1
      wkperenrgy = wksum2
      wktotenrgy = wksum
!
! ---- energy (fields)
!
      bfldy =0.0
      bfldz =0.0
      efldx =0.0
      efldy =0.0
      efldz =0.0
!        DO i=1,nx
!          bfldy=bfldy+(by(i)-bf0(2))**2*0.5
!          bfldz=bfldz+(bz(i)-bf0(3))**2*0.5
!          efldx=efldx+ex(i)**2*0.5
!          efldy=efldy+ey(i)**2*0.5
!          efldz=efldz+ez(i)**2*0.5
!        ENDDO
!
!
          bfldy=0.5*SUM((by(1:nx)-bf0(2))**2)
          bfldz=0.5*SUM((bz(1:nx)-bf0(3))**2)
          efldx=0.5*SUM(ex(1:nx)**2)
          efldy=0.5*SUM(ey(1:nx)**2)
          efldz=0.5*SUM(ez(1:nx)**2)
!
      bfldenrgy=bfldy+bfldz
      efldenrgy=efldx+efldy+efldz
!      print *,' Time, Efld & Bfld Energies: ',time,efldenrgy,bfldenrgy
!      print *,' Time= ',time
!      DO ix=1,nx
!       print *,'  &, Ex, Ey, Ez : ',ex(ix),ey(ix),ez(ix)
!      ENDDO
!
! ----
!
      totenrgy=bfldenrgy + efldenrgy
      DO is=1,nsp
        totenrgy=vkenrgy(is)+totenrgy
      ENDDO
!      print *,' Time, B, E, K and Tot. energies = ', time,bfldenrgy,efldenrgy,vkinenrgy,totenrgy
!      print *,'tpal and tper.  =',tpar(1),tper(1)
!
      print *,'wkpal and wkper eenrgies =',wkpalenrgy,wkperenrgy
      print *,'wktot and tot eenrgies =',wktotenrgy,totenrgy
       enrgys(1) = bfldenrgy
       enrgys(2) = efldenrgy
       enrgys(3) = wkpalenrgy
       enrgys(4) = wkperenrgy
       enrgys(5) = wktotenrgy
       enrgys(6) = totenrgy
!
! --- energy monitor ---
!
!      if (itime == 1) then
!        smemo=total
!      else if(itime /= 1)
!        iflge=0
!        psen=abs(total-smemo)/smemo*100.0
!        if (psen >= 5.) then
!        if (psen >= 10.) then
!           call enepri
!           write(6,55) itime
! 55        format(1h ,'*** poor energy conservation ***, itime=',i5)
!          iflge=1
!       endif
!     endif
!
! ---- momentum
!
      DO ic=1,3
        pmom(ic)=0.0
      ENDDO
!
      n2=0
      DO is=1,nsp
        n1=n2+1
        n2=n2+nions(is)
        ggs=ai(is)*gg(is)
!        DO n=n1,n2
!          pmom(1)=pmom(1)+vxp(n)*ggs
!          pmom(2)=pmom(2)+vyp(n)*ggs
!          pmom(3)=pmom(3)+vzp(n)*ggs
!
          pmom(1)=pmom(1)+ggs*SUM(vxp(n1:n2))
          pmom(2)=pmom(2)+ggs*SUM(vyp(n1:n2))
          pmom(3)=pmom(3)+ggs*SUM(vzp(n1:n2))
!        ENDDO
      ENDDO
      END SUBROUTINE energy
!
!
      SUBROUTINE outdat(itime,time)
      USE share_data
!
      REAL, DIMENSION(nx,nsp) :: dnsp0,dnsp1,dnsh,uxsh,uysh,uzsh
      REAL, DIMENSION(0:nxp1) :: denp1,denh,uxh,uyh,uzh
      REAL, DIMENSION(2) :: tim
      REAL, DIMENSION(6) :: enrgys
      REAL, DIMENSION(nsp) :: tpal,tper,umx
      REAL, DIMENSION(nsp) :: wkpal, wkper
      REAL :: bfldenrgy, efldenrgy,vkinenrgy,totenrgy
      REAL :: time
      INTEGER :: n1,n2,n,i,ip,ix,iy,k,itime
!
      common /momnt1/dnsp0,dnsp1,denp1
      common /momnt2/dnsh,uxsh,uysh,uzsh
      common /momnth/denh,uxh,uyh,uzh
      common /tempi/tpal,tper,umx,enrgys
      common /wkenergy/wkpal,wkper
!
!-----------------------------------------------------------------------
!  #    save field-data in file on unit=10
!  #    save particle-data in file on unit=11
!  #    save energy/diagnostic-data in file on unit=12
!  #    save last time step in unit=13 for restarting process
!  
!  #    In all data sets we save the same header data
!
!-----------------------------------------------------------------------
!
!
!
! --- Time & field data
!
      IF((mod(itime,ifield) == it1)) then
          tim(1)=float(itime)
          tim(2)=time
!
! --- magnetic, velocity and density fields
!
!
	  print *,' *** Writing Field Data ***',itime
      write(10) tim,by(1:nx),bz(1:nx),ex,ey,ez,  &
                uxsh,uysh,uzsh,dnsh,tpal,tper
      call flush(10)
!           
      ENDIF
!
! --- particle data
!
      IF((mod(itime,iparticles) == it2)) then
          tim(1)=float(itime)
          tim(2)=time
	      print *,' *** Writing Particle Data ***',itime
          write(11) tim,xp1,vxp,vyp,vzp
          call flush(11)
      ENDIF
!
! --- Energy data
!
      IF((mod(itime,ienergy) == it3)) then
          tim(1)=float(itime)
          tim(2)=time
!          print *,' *** Writing Energy Data ***',itime,time,umx,vkinenrgy,totenrgy
          write(12) tim,tpal,tper,wkpal,wkper,umx,enrgys
          call flush(12)
      ENDIF
!
! Print Last Time step data for Restarting Code
!
!
! --- Writing restarting file
!
      IF(time >= tmax) then
         open(13,file=trim(home)//'restart.d13',form='binary')	     
!
! ---     Writing Header
!
	     write(13) header
!
! --- Writing Field data
!
         tim(1)=float(itime)
         tim(2)=time
         write(13) tim,by(1:nx),bz(1:nx),ex,ey,ez,  &
               uxsh,uysh,uzsh,dnsh,tpal,tper
         call flush(13)
!
! --- Writing Particle data
!
         write(13) tim,xp1,vxp,vyp,vzp
         call flush(13)
         close(13)
      ENDIF
      END SUBROUTINE outdat
!
      SUBROUTINE filter(a)
      USE share_dim
!
      REAL, DIMENSION(0:nxp1) :: a,g
      REAL, DIMENSION(0:nxp1) :: xi
      INTEGER, DIMENSION(0:nxp1) :: igmx,igox,igpx
      INTEGER :: i
!
      common /xmesh/xi,igmx,igox,igpx
!
! -- binomial filter
!
      do i=0,nxp1
        g(i)=(2.0*a(i)+a(igpx(i))+a(igmx(i)))*0.25
      enddo
!
! -- compensation for binomial filter
!
      do i=0,nxp1
        a(i)=(6.0*g(i)-g(igpx(i))-g(igmx(i)))*0.25
      enddo
      END SUBROUTINE filter
!
      SUBROUTINE filter2
      USE share_data
!
      REAL, DIMENSION(0:nxp1) :: gy,gz,xi
      INTEGER, DIMENSION(0:nxp1) :: igmx,igox,igpx
      INTEGER :: i
!
      common /xmesh/xi,igmx,igox,igpx
!
! -- binomial filter
!
      do i=0,nxp1
        gy(i)=(2.0*by(i)+by(igpx(i))+by(igmx(i)))*0.25
        gz(i)=(2.0*bz(i)+bz(igpx(i))+bz(igmx(i)))*0.25
      enddo
!
! -- compensation for binomial filter
!
      do i=0,nxp1
        by(i)=(6.0*gy(i)-gy(igpx(i))-gy(igmx(i)))*0.25
        bz(i)=(6.0*gz(i)-gz(igpx(i))-gz(igmx(i)))*0.25
      enddo
      END SUBROUTINE filter2
!
      SUBROUTINE smooth(a)
      USE share_dim
!
      REAL, DIMENSION(0:nxp1) :: a,g
      REAL, DIMENSION(0:nxp1) :: xi
      INTEGER, DIMENSION(0:nxp1) :: igmx,igox,igpx
      INTEGER :: i
!
      common /xmesh/xi,igmx,igox,igpx
!
! -- binomial filter
!
      do i=0,nxp1
        g(i)=(2.0*a(i)+a(igpx(i))+a(igmx(i)))*0.25
      enddo
!
      do i=0,nxp1
        a(i)=g(i)
      enddo
      END SUBROUTINE smooth
!
      SUBROUTINE inpdat(itime,time)
      USE share_data
!
!
      REAL, DIMENSION(nheadr) :: hedr
      REAL, DIMENSION(0:nxp1) :: denp1,denh,uxh,uyh,uzh
      REAL, DIMENSION(nx,nsp) :: dnsp0,dnsp1,dnsh,uxsh,uysh,uzsh
      REAL, DIMENSION(2) :: tim
      REAL, DIMENSION(6) :: enrgys
      REAL, DIMENSION(nsp) :: tpal,tper,umx
      REAL :: sumdh,sumd1,sumvx,sumvy,sumvz,time
!      REAL :: bfldenrgy,efldenrgy,vkinenrgy,totenrgy
      INTEGER :: n1,n2,n,i,ip,ix,iy,k,itime0
!
      common /momnt1/dnsp0,dnsp1,denp1
      common /momnt2/dnsh,uxsh,uysh,uzsh
      common /momnth/denh,uxh,uyh,uzh
      common /tempi/tpal,tper,umx,enrgys
!
!-----------------------------------------------------------------------
!  This program reads the last time step in restart.dat (unit 13) 
!  #    read data from datafile in unit=13
!
!-----------------------------------------------------------------------
!
!
      write(6,10)
 10   format(1h,'* ---- input header from unit 16 *')
      read(13) hedr
!
!	  WHERE(hedr /= header)
!		    print *,' Header(i)= ',hedr
!      ENDWHERE
!
!      IF((nx /= jnx).or.(ny /= jny).or.(jnp /= np).or.(jsp /= nis)) then
!         write(6,30) jnx,jny,jnp,jsp
! 30      format(1h,//,'* * * parameter mistmatch * * *',4(i6,1x))
!         stop
!      ENDIF
!
      it1=1
      it2=1
      IF(ifield == 1) it1=0
      IF(iparticles == 1) it2=0
!
! --- input time & magnetic field data & velocity and density per specie 
!
      write(6,40)
 40   format(1h ,'--- input field data from unit 10 ---')
      read(13) tim,by(1:nx),bz(1:nx),ex,ey,ez,  &
               uxsh,uysh,uzsh,dnsh,tpal,tper
!
      itime0=tim(1)
      time0=tim(2)
!
! --- particle data
!
          write(6,50)
 50       format(1h ,'--- input particle data from unit 16---')
          read(13) tim,xp0,vxp,vyp,vzp
!
! periodic boundary conditions for b-field
!
        by(0)=by(nx)
        by(nxp1)=by(1)
        bz(0)=bz(nx)
        bz(nxp1)=bz(1)
!
        END SUBROUTINE inpdat
!
      SUBROUTINE reinitp
      USE share_data
!
! Reinitializes the simulation after reading input data from input
!
!
      INTEGER, DIMENSION(np) :: imx,iox,ipx,imy,ioy,ipy
      REAL, DIMENSION(np) :: smm,smo,smp,som,soo,sop,spm,spo,spp
      REAL, DIMENSION(0:nxp1) :: denp1,denh,uxh,uyh,uzh
      REAL, DIMENSION(nx,nsp) :: dnsp0,dnsp1,dnsh,uxsh,uysh,uzsh
      REAL, DIMENSION(2) :: tim
      REAL, DIMENSION(6) :: enrgys
      REAL, DIMENSION(nsp) :: tpal,tper,umx
      REAL :: sumnp1,sumnh,sumvx,sumvy,sumvz
      REAL :: bfldenrgy,efldenrgy,vkinenrgy,totenrgy
      INTEGER :: n1,n2,n,i,ip,ix,iy,k
!
      common /igridxy/imx,iox,ipx,imy,ioy,ipy
      common /shpfun/smm,smo,smp,som,soo,sop,spm,spo,spp
      common /momnt1/dnsp0,dnsp1,denp1
      common /momnt2/dnsh,uxsh,uysh,uzsh
      common /momnth/denh,uxh,uyh,uzh
      common /tempi/tpal,tper,umx,enrgys
!
!*****************************************************************
!      subroutine meshx will initialize the mesh and assign the
!      indices to the mesh-grids including the periodic boundary
!      conditions on the ghost-grids
!*****************************************************************
!
      call meshx(dx)
!
!
!***********************************************************************
!       take the moment of particles
!         2/2 method (2nd order )
!           ix: nearest mesh number
!           dl: length between particle and the nearest mesh
!          rwn(k,i) : real density of k particle(ion1 or ion2)
!                       at the position whose mesh number is i
!***********************************************************************
!
      do is=1,nsp
	    qm(is)=qi(is)/ai(is)
        do ix=1,nx
          dnsp1(ix,is)=0.0
        enddo
      enddo
!
! //////////////////// shape fun ///////////////////////////////////////
!
      call shapefun(dx,xp1)
!
!
      n2=0
      do is=1,nsp
          n1=n2+1
          n2=n2+nions(is)
!
        do n=n1,n2
          dnsp1(imx(n),is)=dnsp1(imx(n),is) + smm(n)
          dnsp1(iox(n),is)=dnsp1(iox(n),is) + smo(n)
          dnsp1(ipx(n),is)=dnsp1(ipx(n),is) + smp(n)
        enddo
      enddo
!
      do is=1,nsp
        do ix=1,nx
          dnsp1(ix,is)=dnsp1(ix,is)*gg(is)
        enddo
      enddo
!
      do ix=1,nx
        sumnp1=0.
        sumnh=0.
        do is=1,nsp
          sumnp1=sumnp1+dnsp1(ix,is)*qi(is)
          sumnh=sumnh+dnsh(ix,is)*qi(is)
        enddo
        denp1(ix)=sumnp1
        denh(ix)=sumnh
      enddo
!
!
!***********************************************************************
!        periodic boundary conditions
!***********************************************************************
!
           denh  (0)=denh  (nx)
           denp1(0)=denp1(nx)
           denh  (nxp1)=denh  (1)
           denp1(nxp1)=denp1(1)
!
!      call filter(denh)
!      call filter(denp1)
!
        do ix=1,nx
          sumvxh=0.
          sumvyh=0.
          sumvzh=0.
          gdn=1./denh(ix)
          do is=1,nsp
            sumvxh=sumvxh+uxsh(ix,is)*qi(is)
            sumvyh=sumvyh+uysh(ix,is)*qi(is)
            sumvzh=sumvzh+uzsh(ix,is)*qi(is)
          enddo
          uxh(ix)=sumvxh*gdn
          uyh(ix)=sumvyh*gdn
          uzh(ix)=sumvzh*gdn
        enddo
!
!***********************************************************************
!        Periodic boundary conditions
!***********************************************************************
!
           uxh  (0)=uxh  (nx)
           uyh  (0)=uyh  (nx)
           uzh  (0)=uzh  (nx)
           uxh  (nxp1)=uxh  (1)
           uyh  (nxp1)=uyh  (1)
           uzh  (nxp1)=uzh  (1)
!
      END SUBROUTINE reinitp
!
      SUBROUTINE fftx1(ar,ai,br,bi,tr,ti,list,n2,key,lp)
      IMPLICIT NONE
!
! fft1 (fftsub,ffttbl,fftlst)
! transferred from fft2 in 'super-computer' by murata et al.(maruzen)
!
      INTEGER :: indx,indy,n2,key,lp,i,k,l,m
      REAL, DIMENSION(n2*2) :: ar,br,ai,bi
      REAL, DIMENSION(n2*lp,2) :: tr,ti
      INTEGER, DIMENSION(n2*lp) :: list
      REAL :: div
      common /indxy/ indx,indy
!
!.....initialize table & list
!
      IF(indx == 1) then
        call ffttblx1(n2,tr,ti,br,bi,lp)
        call fftlstx1(list,n2,lp)
        indx=0
      ENDIF
!
      k=1
      l=n2
      DO i=1,lp
        m=(i-1)*n2+1
          IF(k == 1) then
            call fftsubx1(ar,ai,br,bi,tr(m,key),ti(m,key),list(m),l,n2)
          ELSE
            call fftsubx1(br,bi,ar,ai,tr(m,key),ti(m,key),list(m),l,n2)
          ENDIF
        k=k*(-1)
        l=l/2
      ENDDO
!
!  if key=2 then inverse fft
!           else normal  fft
!
      IF(key == 2) then
        div=1.e0/(n2*2)
          IF(k == 1) then
            DO i=1,n2*2
              br(i)=ar(i)*div
              bi(i)=ai(i)*div
            ENDDO
          ELSE
            DO i=1,n2*2
              br(i)=br(i)*div
              bi(i)=bi(i)*div
            ENDDO
          ENDIF
      ELSE IF (k == 1) then
        DO i=1,n2*2
          br(i)=ar(i)
          bi(i)=ai(i)
        ENDDO
      ENDIF
      END SUBROUTINE fftx1
!
      SUBROUTINE fftsubx1(ar,ai,br,bi,tr,ti,list,l,n2)
      IMPLICIT NONE
      INTEGER :: i,n2,l
      REAL, DIMENSION(n2*2) :: ar,ai
      REAL, DIMENSION(n2,2) :: br,bi
      REAL, DIMENSION(n2) :: tr,ti
      INTEGER, DIMENSION(n2) :: list
!
      DO i=1,n2
        br(i,1)=ar(list(i))+ar(list(i)+l)*tr(i)-ai(list(i)+l)*ti(i)
        bi(i,1)=ai(list(i))+ai(list(i)+l)*tr(i)+ar(list(i)+l)*ti(i)
      ENDDO
      DO i=1,n2
        br(i,2)=ar(list(i))-ar(list(i)+l)*tr(i)+ai(list(i)+l)*ti(i)
        bi(i,2)=ai(list(i))-ai(list(i)+l)*tr(i)-ar(list(i)+l)*ti(i)
      ENDDO
      END SUBROUTINE fftsubx1
!
      SUBROUTINE ffttblx1(n2,tr,ti,br,bi,lp)
      IMPLICIT NONE
      INTEGER :: n2,lp,i,j,k,l,nn
      REAL, DIMENSION(n2,2) :: br,bi
      REAL, DIMENSION(n2,lp,2) :: tr,ti
      REAL :: pi,trr,tii
!
      pi=atan(1.e0)*4.e0
      DO i=1,n2
        trr=cos(2.e0*pi*(i-1)/(n2*2))
        tii=sin(2.e0*pi*(i-1)/(n2*2))
        br(i,1)= trr
        bi(i,1)=-tii
        br(i,2)= trr
        bi(i,2)= tii
      ENDDO
!
      k=1
      nn=n2
!
      DO l=1,lp
        DO j=0,k-1
          DO i=1,nn
            tr(i+j*nn,l,1)=br(1+nn*j,1)
            ti(i+j*nn,l,1)=bi(1+nn*j,1)
            tr(i+j*nn,l,2)=br(1+nn*j,2)
            ti(i+j*nn,l,2)=bi(1+nn*j,2)
          ENDDO
        ENDDO
        k=k*2
        nn=nn/2
      ENDDO
      END SUBROUTINE ffttblx1
!
      SUBROUTINE fftlstx1(list,n2,lp)
      IMPLICIT NONE
      INTEGER :: n1,n2,lp,nn,i,j,k,m
      INTEGER, DIMENSION(n2,lp) :: list
!
      n1=n2
      nn=1
      DO k=1,lp
        m=0
        DO j=1,nn
          DO i=1,n1
            m=m+1
            list(m,k)=i+(j-1)*2*n1
          ENDDO
        ENDDO
        n1=n1/2
        nn=nn*2
      ENDDO
      END SUBROUTINE fftlstx1
!
      SUBROUTINE tran1(br,bi,ar1,ai1,ar2,ai2,nx)
      IMPLICIT NONE
      INTEGER :: nx,nx2,nxh1,j
      REAL, DIMENSION(nx) :: br,bi,ar1,ai1,ar2,ai2
!
      ar1(1)=br(1)
      ai1(1)=0.0
      ar2(1)=bi(1)
      ai2(1)=0.0
      nx2=nx+2
      nxh1=nx/2+1
      DO j=2,nxh1
        ar1(j)=+0.5*(br(j)+br(nx2-j))
        ai1(j)=+0.5*(bi(j)-bi(nx2-j))
        ar2(j)=+0.5*(bi(j)+bi(nx2-j))
        ai2(j)=-0.5*(br(j)-br(nx2-j))
      ENDDO
!
      DO j=2,nxh1
        ar1(nx2-j)=+0.5*(br(j)+br(nx2-j))
        ai1(nx2-j)=-0.5*(bi(j)-bi(nx2-j))
        ar2(nx2-j)=+0.5*(bi(j)+bi(nx2-j))
        ai2(nx2-j)=+0.5*(br(j)-br(nx2-j))
      ENDDO
	  return
!
      ENTRY tran2(br,bi,ar1,ai1,ar2,ai2,nx)
      DO j=1,nx
        br(j)=ar1(j)-ai2(j)
        bi(j)=ai1(j)+ar2(j)
      ENDDO
      END SUBROUTINE tran1
! &INPUTDAT
! VERSION        = 'Testing Akimotos Proton-Core/Proton-beam Instabilities'
! XLEN           = 256.0
! DT           = 0.05
! TMAX         = 200.0
! BETAEN       = 1.0
! WPIWCI       = 1.e4
! GAME         = 1.0
! DUS          = 10.0
! BNAMP        = 1.e-6
! BF0          = 1.0,0.0,0.0
! ITRESTART    = 0
! IFIELD       = 1
!! IPARTICLES   = 8
! IENERGY      = 1
! IFILTER      = 4
! ETA          = 0.0
! &END
 


