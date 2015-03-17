      subroutine fdist(NZ,TEMP,TFERMI,WR,WZ,XMK,PAR,CPF,ISEED,JCOUNT,IP
     & ,SETMAX,DIVE,IDNUM,IDN)
      use constant
      implicit real*8(a-h,o-z)
      parameter(NDP=40,IDIM=3)
      dimension PAR(NZ,7)
      dimension XPP(NZ),PXP(NZ),YPP(NZ),PYP(NZ),ZPP(NZ),PZP(NZ)
      dimension NNX(NZ),NNY(NZ),NNZ(NZ)
      dimension DIVE(2,IDN),IDNUM(IDN),SETMAX(6)
c
c --- check 
c      write(*,*)'NZ',NZ
c      write(*,*)'TFERMI',TFERMI
c      write(*,*)'TEMP',TEMP,TEMP/TFERMI
c      write(*,*)SQRT(8d0*XMK*FERE),'PM'
c      write(*,*)SQRT(8d0*FERE/XMK)/OMEGA,'XM'
c      write(*,*)BoltzK*TEMP,'KT'
c      write(*,*)OMEGA,'wkav'
c      write(*,*)CPF,'CPF'
c      write(*,*)CPF/(TFermi*BoltzK),'CPF/FE dist'
c --- idle running
c      if (itimes .eq. 0)then
c         itimes=100
c         ISEED=321695
c      endif
c ---
!      if (IP .EQ. 1)write(*,*)'Fermion'
!      if (IP .EQ. 0)write(*,*)'Boson Ground State'
!      if (IP .EQ. -1)write(*,*)'Boson Excited State'
c      write(*,*)'IP',IP
c      do j=1,1000
c         dummy=prand(ISEED)
c      enddo

c --- test
c      write(*,*)TFERMI,'TFERMI'
c      write(*,*)TEXP,'TEMP',TEXP/TFERMI
c      write(*,*)'CPF/EF',CPF/(TFermi*BoltzK)
c ---
c --- OPT
c      call RPOPT(NZ,TEMP,TFERMI,OMEGA,CPF,IP,EOPT)      
c ---
!      write(*,*)'ISEED IN BEFORE',ISEED
      OMEGA=(WR**2*WZ)**(1d0/3d0)
c --- start MONTE CALRO
         call MKPTD(XPP,PXP,YPP,PYP,ZPP,PZP,NZ,WR,WZ,CPF,TEMP,TFERMI,
     &     ISEED,XMK,R,P,JCOUNT,IP,SETMAX,DIVE,IDNUM,IDN)
c         write(*,*)XPP(NZ-120),PXP(NZ-25)
c --- conversion random number to SI
         call RSI(XPP,PXP,YPP,PYP,ZPP,PZP,
     &        XMK,OMEGA,TEMP,TFERMI,NZ,PAR,IP)
c         do JJ=1,NZ
c         write(*,*)JJ,PAR(JJ,2)
c         enddo
c
!         write(*,*)'ISEED IN AFTER',ISEED

c ---
c      do jj=1,NZ
c      write(101,100)XPP(jj),' ',PXP(jj)
c      write(102,100)YPP(jj),' ',PYP(jj)
c      write(103,100)ZPP(jj),' ',PZP(jj)
c      write(104,100)R(jj),' ',P(jj)
c      enddo
c ---      
c
c      write(*,*)JCOUNT
c ---
c         do jj=1,NZ
c         if(IP .EQ. -1)then
c            write(*,*)jj,'soto'
c            do kk=1,7
c            write(*,*)PAR(jj,kk)
c            enddo
c         endif
c         enddo
c ---

 100  format(f0.15,a1,f0.15)
      return
      end
c===============================================
c     fermi distribution at 3D
c================================================
      function FERPTD(X,PX,Y,PY,Z,PZ,XMK,WR,WZ,CPF,TEMP,TFERMI,FNM)
      use constant
      implicit real*8(a-h,o-z)
!c
      BETA=BoltzK*TEMP
      E=(PX**2+PY**2+PZ**2)/(2d0*XMK)
      EZ=XMK*WZ**2*(Z**2)/2d0
      E=E+EZ+XMK*WR**2*(X**2+Y**2)/2d0-CPF
c      E=(PXT**2+PYT**2+PZT**2)/(2d0*XMK)
c      E=E+XMK*OMEGA**2*(XT**2+YT**2+ZT**2)/2d0-CPF
c
      FERPTD=1d0/(DEXP(E/BETA)+1d0)
c      F1=FNM(CPF,TEMP)
c --- without normalize
      FERPTD=FERPTD/FNM
c ---
c
      return
      end
c==============================================================
c     BOSE distribution at 3d
c==============================================================
      function BOSETD(X,PX,Y,PY,Z,PZ,XMA,WR,WZ,CPB,TEMP,TCRIT,IP)
      USE CONSTANT
      implicit real*8(a-h,o-z)
c
!      parameter(pi=3.141592654d0)
!      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)
c      data itimes/0/
c
      BETA=BoltzK*TEMP
      if(IP .EQ. 0)then
!c --- new write wigner distribution
        ER1=(PX**2+PY**2)/XMA/WR/hbar
        ER2=(X**2+Y**2)*XMA*WR/hbar
        EZ1=(PZ**2)/XMA/WZ/hbar
        EZ2=(Z**2)*XMA*WZ/hbar
        EP=ER1+EZ1
        ER=ER2+EZ2
c --- old
c        BOSETD=DEXP(-EP/XMA/OMEGA/hbar)*DEXP(-XMA*OMEGA*ER/hbar/2d0)
c --- new
        BOSETD=DEXP(-EP)*DEXP(-ER)
c ---
         else   
!c ---
      E=(PX**2+PY**2+PZ**2)/(2d0*XMA)
      EZ=XMA*WZ**2*Z**2/2d0
      E=E+EZ+XMA*WR**2*(X**2+Y**2)/2d0-CPB
c
      BOSETD=1d0/(DEXP(E/BETA)-1d0)
c      ZAA=DEXP(CPB/BETA)
c      ZAA=ZAA/(1d0-ZAA)
c      BH=1d0/(DEXP((hbar*OMEGA-CPB)/BETA)-1d0)
c      if (TEMP/TCRIT .le. 1d0)then
c --- without normalize
!      BH=1d0/(DEXP(-CPB/BETA)-1d0)
!         if(BH .le. 1d0)BOSETD=BOSETD/BH
c ---
!      else
!         BH=1d0/(DEXP(-CPB/BETA)-1d0)
!      endif
c --- not need
!c      if(IP .EQ. 0)then
!c      BH=(ZAA+BH)/2d0
!c      endif
c ---
c      write(*,*)'BB,BH,IP',BOSETD,BH,IP
!      BOSETD=BOSETD/BH
      endif
c
! 220  continue
c      BOSETD=BOSETD
c      write(*,*)'BOSETD',BOSETD,BH,IP
c
c --- new write thomas fermi distribution
c      if(IP .EQ. 0)then
c         EZERO=XMK*OMEGA**2*(XT**2+YT**2+ZT**2)/2d0
c         BOSETD=(CPB-EZERO)/CPB
c      endif
c ---
c      write(*,*)BOSETD,'CHECK'
c      write(*,*)ZAA,BH,'check'
c ---
c
      return
      end
c===============================================
      subroutine CONXP(X,PX,Y,PY,Z,PZ,XMAX,PXMAX,YMAX,PYMAX,ZMAX,PZMAX
     &,IP)
!     IP=-1:random to SI         IP=1:SI to random
c!     IP=-2:random to SI for BEC IP=2:SI to random for BEC
c      use constant
      implicit real*8(a-h,o-z)
c --- 
      if (IP .EQ. -1)then
      PXT=(PX-0.5d0)
      XT=(X-0.5d0)
      PYT=(PY-0.5d0)
      YT=(Y-0.5d0)
      PZT=(PZ-0.5d0)
      ZT=(Z-0.5d0)
c ---
      X=PXT*XMAX
      PX=XT*PXMAX
      Y=PYT*YMAX
      PY=YT*PYMAX
      Z=PZT*ZMAX
      PZ=ZT*PZMAX
      endif
c ---
      if (IP .EQ. 1)then
      X=(X/XMAX+0.5d0)
      PX=(PX/PXMAX+0.5d0)
      Y=(Y/YMAX+0.5d0)
      PY=(PY/PYMAX+0.5d0)
      Z=(Z/ZMAX+0.5d0)
      PZ=(PZ/PZMAX+0.5d0)
      endif
c ---    
      return
      end
c===============================================
c     make distribution (return values are SI UNIT)
c===============================================
      subroutine MKPTD(XPP,PXP,YPP,PYP,ZPP,PZP,NZ,WR,WZ,
     & CPS,TEMP,TSP,ISEED,XMA,R,P,JCOUNT,IP,SETMAX,DIVE,IDNUM,IDN)
      use constant
      implicit real*8(a-h,o-z)
c      parameter (IDN=300)
      dimension XPP(NZ),PXP(NZ),YPP(NZ),PYP(NZ),ZPP(NZ),PZP(NZ)
      dimension R(NZ),P(NZ),NNX(NZ),NNY(NZ),NNZ(NZ)
      dimension DIVE(2,IDN),IDNUM(IDN),INUM(IDN),SETMAX(6)
      external FERPTD
      external BOSETD
      external PRAND
      external FNM
      external FEN
c      data itimes/0/
c
c ---
c      call RPOPT(NZ,TEMP,TFERMI,OMEGA,CPF,IP,EOPT)      
c --- initialize
      OMEGA=(WR**2*WZ)**(1d0/3d0)
      NNX=0
      NNY=0
      NNZ=0
      INUM=0
c ---
      JCOUNT=0
      JLIM=NZ*20
c --- initialize
      XPP=0d0
      PXP=0d0
      YPP=0d0
      PYP=0d0
      ZPP=0d0
      PZP=0d0
c ---      
      jp=0
      NCOUNT=0
c --- max setting
c      write(*,*)'start max set'
      XMAX=SETMAX(1)
      PXMAX=SETMAX(2)
      YMAX=SETMAX(3)
      PYMAX=SETMAX(4)
      ZMAX=SETMAX(5)
      PZMAX=SETMAX(6)   
!      if (IP .NE. 0)then
!         call SETCM(XMAX,PXMAX,YMAX,PYMAX,ZMAX,PZMAX,
!     &        XMA,WR,WZ,TEMP,TSP,BETA,FE,IP,CPS,ENMAX)
!         call divide(IP,CPS,TEMP,OMEGA,DIVE,IDNUM,ENMAX,IDN)
!c         write(*,*)XMAX,PXMAX,'X,PX'
!      else
!         call SETCZ(XMAX,PXMAX,YMAX,PYMAX,ZMAX,PZMAX,
!     &        XMA,WR,WZ,TEMP,TSP,BETA,FE)
!      endif
c      if (IP .EQ. 1)F1=FNM(CPS,TEMP)
      if (IP .EQ. 1)F1=1d0
c ---
c      write(*,*)'X-ZMAX',XMAX,ZMAX
c      write(*,*)'PX-PZMAX',PXMAX,PZMAX

      do while(NCOUNT .LT. NZ)
      jp=jp+1
      XPC=prand(ISEED)
      PXC=prand(ISEED)
      YPC=prand(ISEED)
      PYC=prand(ISEED)
      ZPC=prand(ISEED)
      PZC=prand(ISEED)
      CPP=prand(ISEED)
c      if (MOD(NCOUNT,1000) .EQ. 0)write(*,*)XPC,PXC,'BEFORE'
      call CONXP(XPC,PXC,YPC,PYC,ZPC,PZC,XMAX,PXMAX,YMAX,PYMAX,ZMAX,
     & PZMAX,-1)
c      if (MOD(NCOUNT,1000) .EQ. 0)write(*,*)XPC,PXC,'AFTER'
c
      if (IP .EQ. 1)then
      CHP=FERPTD(XPC,PXC,YPC,PYC,ZPC,PZC,XMA,WR,WZ,CPS,TEMP,TSP,F1)
      else
      CHP=BOSETD(XPC,PXC,YPC,PYC,ZPC,PZC,XMA,WR,WZ,CPS,TEMP,TSP,IP)
      endif
c      if(IP .eq. -1)then
c      write(211,*)CHP,CPP
c      endif
c
      CHKEN=fen(XMA,WR,WZ,XPC,YPC,ZPC,PXC,PYC,PZC)
      if (IP .NE. 0)then
      do is=1,IDN
         if (CHKEN .le. DIVE(1,is))then !if1
            if (INUM(is) .lt. IDNUM(is))then !if2
               CHP=CHP/DIVE(2,is)
               if (CPP .LE. CHP)then !if3
                  NCOUNT=NCOUNT+1
                  XPP(NCOUNT)=XPC
                  PXP(NCOUNT)=PXC
                  YPP(NCOUNT)=YPC
                  PYP(NCOUNT)=PYC
                  ZPP(NCOUNT)=ZPC
                  PZP(NCOUNT)=PZC

c --- 
c         if (IP. EQ. -1)then
c            write(*,*)XPC,PXC,YPC,PYC,ZPC,PZC
c            write(*,*)NCOUNT
c         endif
c ---
c         if(IP .EQ. -1)write(*,*)NCOUNT,'genarated'
c         if(IP .EQ. 0)write(*,*)NCOUNT,'genarated'
c============================================================
!c --- pauli principle new
         if (IP .EQ. 1)then!1        !pauli
         call SETN(PXP(NCOUNT),XPP(NCOUNT),PYP(NCOUNT),YPP(NCOUNT), !pauli
     &    PZP(NCOUNT),ZPP(NCOUNT),WR,WZ,XMA,NNX(NCOUNT),NNY(NCOUNT), !pauli
     &    NNZ(NCOUNT),TEMP,TSP)                                       !pauli
         !write(*,*)NNX(NCOUNT),NNY(NCOUNT),NNZ(NCOUNT) !pauli
         if (NCOUNT .eq. 1)then !pauli
            ISIGN=1 !pauli
            else    !pauli
               call PAULIN(NNX,NNY,NNZ,NCOUNT,NZ,ISIGN) !pauli
         endif         !pauli
         if (ISIGN .EQ. 0)NCOUNT=NCOUNT-1     !pauli
         if (ISIGN .EQ. 0)JCOUNT=JCOUNT+1     !pauli
         if (ISIGN .EQ. 1)INUM(is)=INUM(is)+1 !pauli

         endif !1    !pauli
!c --- pauli new end
!        !!check
!c         !if (ISIGN .EQ. 1)write(*,*)'generated',NCOUNT
!c ---
!c
!c --- pauli principle old
!c         if (itimes .eq. 0)then
!c            itimes=100
!c            ISIGN=1
!c         else
!c            call PAULI(XPP,PXP,YPP,PYP,ZPP,PZP,
!c     &           XMK,OMEGA,TEMP,TSP,NCOUNT,ISIGN)  
!c         endif
!c         if (ISIGN .EQ. 0)NCOUNT=NCOUNT-1
!c         if (ISIGN .EQ. 0)JCOUNT=JCOUNT+1
!c ---
!c         if (ISIGN .EQ. 1)write(*,*)'generated',NCOUNT !CHECK
!c ---
         if (IP .NE. 1)then              !pauli
                  INUM(is)=INUM(is)+1
               endif                     !pauli
                  exit
c ---
               else
                  exit
               endif !if3
            else
               exit
            endif !if2
         endif !if1
      enddo
!c
               if(JCOUNT .GT. JLIM)then !can't create distribution
                  write(*,*)'stopped'
                  exit
               endif
c --- for BEC
         else
          if (CPP .LE. CHP)then
             NCOUNT=NCOUNT+1
             XPP(NCOUNT)=XPC
             PXP(NCOUNT)=PXC
             YPP(NCOUNT)=YPC
             PYP(NCOUNT)=PYC
             ZPP(NCOUNT)=ZPC
             PZP(NCOUNT)=PZC
          endif
c ---
      endif



c============================================================
      enddo
c
      if (NCOUNT .GE. NZ)write(*,*)'distribution created'
      if (IP .EQ. 1)then
      write(*,*)'excluded',JCOUNT
      write(*,*)'fermion count',NCOUNT,jp
      write(1231,*)TEMP,jp
      else
c      write(*,*)'no excluded particle'
      write(*,*)'boson count',NCOUNT,jp
      if (IP .EQ. 0)write(1232,*)TEMP,jp
      if (IP .EQ. -1)write(1233,*)TEMP,jp
      endif
c      write(*,*)TEMP/TFERMI,'TEMP/TFERMI'
c      write(1001,*)TEMP/TFERMI,NCOUNT-1,JCOUNT
c
c      write(*,*)'I',I
c      write(*,*)NZ
c      write(*,*)'NCOUNT',NCOUNT
 100  format(f0.15,a1,f0.15)
      return
      end
c 
c==================================================
c==================================================
c     set XMAX,PXMAX...
c==================================================
      subroutine SETCM(XMAX,PXMAX,YMAX,PYMAX,ZMAX,PZMAX,
     & XMK,WR,WZ,TEMP,TFERMI,BETA,FE,IP,CPS,ENMAX)
      use constant
      implicit real*8(a-h,o-z)
c
      FE=TFERMI*BoltzK
      BETA=TEMP*BoltzK
c      AA=EOPT*2d0
c --- vs kbT
c      CPE=CPS/BETA
c      ENMAX=25d0*BETA+CPS
c      AA=((25d0*BETA+CPS))*2d0      
      AA=ENMAX*2d0
c --- variable
c      if (IP .EQ. 1)then
c         AA=50d0*TEMP/TFERMI
c      endif
c      if (IP .EQ. -1)then
c         AA=40d0*TEMP/TFERMI
c      endif      
c
c --- default
c      AA=40d0*TEMP/TFERMI
c ---
c      AA=50d0*TEMP/TFERMI
c ---     
c      if (AA .le. 8d0)then
c         AA=8d0
c      endif
!      write(*,*)'AA',AA
c
c      XMAX=(AA*FE/(XMK*OMEGA**2))**(0.5d0)
c      PXMAX=(AA*FE*XMK)**(0.5d0)
c      YMAX=(AA*FE/(XMK*OMEGA**2))**(0.5d0)
c      PYMAX=(AA*FE*XMK)**(0.5d0)
c      ZMAX=(AA*FE/(XMK*OMEGA**2))**(0.5d0)
c      PZMAX=(AA*FE*XMK)**(0.5d0)
      XMAX=(AA/(XMK*WR**2))**(0.5d0)
      PXMAX=(AA*XMK)**(0.5d0)
      YMAX=(AA/(XMK*WR**2))**(0.5d0)
      PYMAX=(AA*XMK)**(0.5d0)
      ZMAX=(AA/(XMK*WZ**2))**(0.5d0)
      PZMAX=(AA*XMK)**(0.5d0)
c
      return
      end
c
c=================================================
c     set XMAX,PXMAX... under TC boson
c==================================================
      subroutine SETCZ(XMAX,PXMAX,YMAX,PYMAX,ZMAX,PZMAX,
     & XMK,WR,WZ,TEMP,TCRIT,BETA,FE)
      use constant
      implicit real*8(a-h,o-z)
c
!      parameter(pi=3.141592654d0)
!      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)
c
c      FE=TCRIT*BoltzK
c      FE=CPB
c      BETA=TEMP*BoltzK
c --- default
c      FE=hbar*WR*60d0
c      FZ=hbar*WZ*60d0
c ---
      FE=hbar*WR*40d0
      FZ=hbar*WZ*40d0
      AA=1d0
c      AA=hbar*OMEGA/FE*20.0d0
c      write(*,*)AA,'AA'

c      XMAX=(AA*FE/(XMK*OMEGA**2))**(0.5)
c      PXMAX=(AA*FE*XMK)**(0.5)
c      YMAX=(AA*FE/(XMK*OMEGA**2))**(0.5)
c      PYMAX=(AA*FE*XMK)**(0.5)
c      ZMAX=(AA*FE/(XMK*OMEGA**2))**(0.5)
c      PZMAX=(AA*FE*XMK)**(0.5)
      XMAX=(AA*FE/(XMK*WR**2))**(0.5d0)
      PXMAX=(AA*FE*XMK)**(0.5d0)
      YMAX=(AA*FE/(XMK*WR**2))**(0.5d0)
      PYMAX=(AA*FE*XMK)**(0.5d0)
      ZMAX=(AA*FZ/(XMK*WZ**2))**(0.5d0)
      PZMAX=(AA*FZ*XMK)**(0.5d0)
c
c      write(*,*)XMAX,'XMAX'
c
      return
      end

c================================================
c     E=n*hbar*OMEGA > n=E/(hbar*OMEGA)
c================================================
      subroutine SETN(PXP,XPP,PYP,YPP,PZP,ZPP,WR,WZ,
     & XMK,NNX,NNY,NNZ,TEMP,TFERMI)
      use constant
      implicit real*8(a-h,o-z)
!   
c      EX=(PXP**2)/(2d0*XMK)+XMK*(OMEGA*XPP)**2/2d0
c      EY=(PYP**2)/(2d0*XMK)+XMK*(OMEGA*YPP)**2/2d0
c      EZ=(PZP**2)/(2d0*XMK)+XMK*(OMEGA*ZPP)**2/2d0
      EX=(PXP**2)/(2d0*XMK)+XMK*(WR*XPP)**2/2d0
      EY=(PYP**2)/(2d0*XMK)+XMK*(WR*YPP)**2/2d0
      EZ=(PZP**2)/(2d0*XMK)+XMK*(WZ*ZPP)**2/2d0
c --- old one
c      NNX=IDINT(EX/(hbar*OMEGA))
c      if (EX/(hbar*OMEGA)-NNX .GE. 0.5d0)then
c         NNX=NNX+1
c      endif
c ---
      NNX=IDINT(EX/(hbar*WR))
      if (EX/(hbar*WR)-NNX .GE. 0.5d0)then
         NNX=NNX+1
      endif
      NNY=IDINT(EY/(hbar*WR))
      if (EY/(hbar*WR)-NNY .GE. 0.5d0)then
         NNY=NNY+1
      endif
      NNZ=IDINT(EZ/(hbar*WZ))
      if (EZ/(hbar*WZ)-NNZ .GE. 0.5d0)then
         NNZ=NNZ+1
      endif
c
      return
      end
c================================================
c     pauli principle NX,NY,NZ
c================================================
      subroutine PAULIN(NNX,NNY,NNZ,NCOUNT,NZ,ISIGN)
      use constant
      implicit real*8(a-h,o-z)
      dimension NNX(NZ),NNY(NZ),NNZ(NZ)
c
      ISIGN=1
      do i=1,NCOUNT-1
         NA=ABS(NNX(NCOUNT)-NNX(i))
         NB=ABS(NNY(NCOUNT)-NNY(i))
         NC=ABS(NNZ(NCOUNT)-NNZ(i))
         if (NA+NB+NC .EQ. 0)then
         ISIGN=0
         exit
         endif
      enddo
c
      return
      end
c
c================================================
c     pauli principle
c=================================================
      subroutine PAULI(XPC,PXC,YPC,PYC,ZPC,PZC,
     & XMK,OMEGA,TEMP,TFERMI,NCOUNT,ISIGN)
      use constant
      implicit real*8(a-h,o-z)
c     
!      parameter(pi=3.141592654d0)
!      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)
      dimension XPC(NCOUNT),PXC(NCOUNT),YPC(NCOUNT),PYC(NCOUNT)
      dimension ZPC(NCOUNT),PZC(NCOUNT)
      external prand
c
      IC=0
c      h=2d0*pi*hbar
      CA=(1d0/(2d0*XMK))**(0.5d0)
      CB=(XMK/2d0)**(0.5d0)*OMEGA
      AA=1.0d0
c      AA=2.0d0
      call SETCM(XMAX,PXMAX,YMAX,PYMAX,ZMAX,PZMAX,
     & XMK,OMEGA,TEMP,TFERMI,BETA,FE)
      PXN=(PXC(NCOUNT)-0.5d0)*AA*PXMAX
      XPN=(XPC(NCOUNT)-0.5d0)*AA*XMAX
      PYN=(PYC(NCOUNT)-0.5d0)*AA*PYMAX
      YPN=(YPC(NCOUNT)-0.5d0)*AA*YMAX
      PZN=(PZC(NCOUNT)-0.5d0)*AA*PZMAX
      ZPN=(ZPC(NCOUNT)-0.5d0)*AA*ZMAX
c
      do ii=1,10000 !1

      do jj=1,NCOUNT-1 !2
         PXT=(PXC(jj)-0.5d0)*AA*PXMAX
         XPT=(XPC(jj)-0.5d0)*AA*XMAX
         PYT=(PYC(jj)-0.5d0)*AA*PYMAX
         YPT=(YPC(jj)-0.5d0)*AA*YMAX
         PZT=(PZC(jj)-0.5d0)*AA*PZMAX
         ZPT=(ZPC(jj)-0.5d0)*AA*ZMAX
c
         DPX=ABS(PXN-PXT)
         DXP=ABS(XPN-XPT)
         DPY=ABS(PYN-PYT)
         DYP=ABS(YPN-YPT)
         DPZ=ABS(PZN-PZT)
         DZP=ABS(ZPN-ZPT)
c
         DALL=DPX*DXP*DPY*DYP*DPZ*DZP
         DR=SQRT(DXP**2+DYP**2+DZP**2)
         DP=SQRT(DPX**2+DPY**2+DPZ**2)
c
c --- DR*DP .LT. h
         if(DR*DP .LE. h)then
            R=DSQRT(XPN**2+YPN**2+ZPN**2)
            P=DSQRT(PXN**2+PYN**2+PZN**2)
            theta=2d0*pi*prand(ISEED)
            phi=2d0*pi*prand(ISEED)
            xi=2d0*pi*prand(ISEED)
            eta=2d0*pi*prand(ISEED)
            XPN=R*sin(theta)*cos(phi)            
            YPN=R*sin(theta)*sin(phi)
            ZPN=R*cos(theta)
            PXN=P*sin(xi)*cos(eta) 
            PYN=P*sin(xi)*sin(eta)
            PZN=P*cos(xi)
c ---     
            ISIGN=0            
            exit
c ---
            else
            ISIGN=1
         endif
c
c
c --- DX*DPX*PY*DPY*DZ*DPZ .LT. h**3
c         if(DALL .LE. h**3)then
c            ISIGN=0
c            exit
c         else
c            ISIGN=1
c         endif
c
c --- DX*DPX .LT. h .and. DY*DPY .LT. h .and. DZ*DPZ .LT. h 
c         if(DPX*DXP .LT. h)then
c            if(DPY*DYP .LT. h)then
c               if(DPZ*DZP .LT. h)then
c            ISIGN=0
c            exit
c         else
c            ISIGN=1c
c         endif
c         else
c            ISIGN=1
c         endif
c         else
c            ISIGN=1
c         endif
c ---
c
      enddo !2
      if (ISIGN .EQ. 1)then
         PXC(NCOUNT)=0.5d0+PXN/PXMAX
         XPC(NCOUNT)=0.5d0+XPN/XMAX
         PYC(NCOUNT)=0.5d0+PYN/PYMAX
         YPC(NCOUNT)=0.5d0+YPN/YMAX
         PZC(NCOUNT)=0.5d0+PZN/PZMAX
         ZPC(NCOUNT)=0.5d0+ZPN/ZMAX
         exit
      endif

      enddo !1
c      
      return
      end
c
c==============================================
c     convert random to SI 
c==============================================
      subroutine RSI(XPC,PXC,YPC,PYC,ZPC,PZC,
     & XMK,OMEGA,TEMP,TFERMI,NZ,PAR,IP)
      USE CONSTANT
      implicit real*8(a-h,o-z)
      dimension XPC(NZ),PXC(NZ),YPC(NZ),PYC(NZ)
      dimension ZPC(NZ),PZC(NZ),PAR(NZ,7)
c
c===================================================
!      write(*,*)'check RSI'
!      if (IP .EQ. 0)then
!         call SETCZ(XMAX,PXMAX,YMAX,PYMAX,ZMAX,PZMAX,
!     & XMK,OMEGA,TEMP,TFERMI,BETA,BE)
!c         write(*,*)XMAX,YMAX,ZMAX
!c         write(*,*)PXMAX,PYMAX,PZMAX
!         else
!      call SETCM(XMAX,PXMAX,YMAX,PYMAX,ZMAX,PZMAX,
!     & XMK,OMEGA,TEMP,TFERMI,BETA,FE,EOPT)
!      endif
!c
!      AA=1d0
!c      AA=2d0
      do jj=1,NZ
!c ---
!c         if(IP .EQ. -1)then
c            write(*,*)jj,'mae'
c            write(*,*)PXC(jj),XPC(jj),PYC(jj),YPC(jj),PZC(jj),ZPC(jj)
!c         endif
!c ---
         X=XPC(jj)
         Y=YPC(jj)
         Z=ZPC(jj)
         PX=PXC(jj)
         PY=PYC(jj)
         PZ=PZC(jj)
         PR=X*Y*Z*PX*PY*PZ
c         if(PXC(jj)*XPC(jj)*PYC(jj)*YPC(jj)*PZC(jj)*ZPC(jj) .eq. 0)
c     &    then
         if(PR .eq. 0)then
         PAR(jj,7)=1
         else
c         PAR(jj,1)=XPC(jj)
c         PAR(jj,2)=YPC(jj)
c         PAR(jj,3)=ZPC(jj)
c         PAR(jj,4)=PXC(jj)
c         PAR(jj,5)=PYC(jj)
c         PAR(jj,6)=PZC(jj)
         PAR(jj,1)=X
         PAR(jj,2)=Y
         PAR(jj,3)=Z
         PAR(jj,4)=PX
         PAR(jj,5)=PY
         PAR(jj,6)=PZ
c ---
         PAR(jj,7)=0
         endif
c ---
c         if(IP .EQ. -1)then
c            write(*,*)jj,'ato'
c            do kk=1,7
c            write(*,*)PAR(jj,kk)
c            enddo
c         endif
c ---
      enddo
c
      return
      end
c
c=====================================================
c     decide chemical potential
c
c     Fermion:IFB=1,Boson:IFB=-1
c======================================================
      subroutine CPBF(NZ,TEMP,TSP,OMEGA,CPP,IFB)
      use constant
      implicit real*8(a-h,o-z)
c
!      parameter(pi=3.141592654d0)
!      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)
!      parameter(h=2d0*pi*hbar)
      External GE,DE        !integrand
c    
c ---
c      write(*,*)'NZ',NZ
c      write(*,*)'TSP',TSP
c      write(*,*)'OMEGA',OMEGA
c      write(*,*)'IFB',IFB
c ---
c --- episilon,delta
        EPS=TSP*BOLTZK*10d-10
        DEL=EPS*0.1d0
c ---         
c      write(*,*)'TEMP',TEMP/TSP
c      write(*,*)'EPS',EPS
c      write(*,*)'DEL',DEL
c
c ---     chemical potential
c
         CPPOINT=10000
c
         if(IFB .EQ. 1)then
!            write(*,*)'Fermion Temp',TEMP/TSP
         CPMAX=1d0*TSP*BoltzK
         else
!            write(*,*)'Boson Temp',TEMP/TSP
         CPMAX=-10d-45
         endif
c
         CPMIN=-15d0*TSP*BoltzK
         CPSTEP=(CPMAX-CPMIN)/(CPPOINT-1)
c     
         VAL2=dble(NZ)
         CP=CPMIN
         do INUM=1,CPPOINT
            CP1=CP
            CP=CPMIN+(CPSTEP*(INUM-1))
            CP2=CP
c
c ---     integration
c
            emin=0d0
c --- default max
c            emax=10d0*TSP*BoltzK
c --- not default
            emax=100d0*TSP*BoltzK
c ---
            INTE=10000
c     
            VAL1=VAL2
           call trap(emin,emax,INTE,VAL,NZ,CP,TEMP,OMEGA,IFB)
            VAL2=VAL
c            !write(*,*)VAL
c
            JCHK=0
            if (VAL1*VAL2 .LT. 0d0)then
               JCHK=1
c --- check 
c               !write(*,*)CP1/(TSP*BOLTZK),'CP1'
c               !write(*,*)CP2/(TSP*BOLTZK),'CP2'
c --- start bisection
               do intc=1,100000 !start bisection
                  CPT=(CP1+CP2)/2d0
        call trap(emin,emax,INTE,fx,NZ,CP2,TEMP,OMEGA,IFB)
        call trap(emin,emax,INTE,fy,NZ,CPT,TEMP,OMEGA,IFB)
               if(fx*fy .ge. 0)then
                  CP2=CPT
                  else
                  CP1=CPT
               endif
               if (abs(CP1-CP2) .le. EPS*100d0)then
                  CP2=(CP1+CP2)/2d0
                  CPP=CP2
                  VAL2=0d0
                  exit
               endif
               if (intc .GE. 10000)then
                  write(*,*)'failed at bisection'
                  VAL2=0d0
                  exit
               endif
               enddo !end bisection
c              ! write(*,*)CP2,'CPT'
c ---
c --- start newton's method
               do intb=1,10000 !start newton
                  CP1=CP2
        call trap(emin,emax,INTE,fa,NZ,CP1+del,TEMP,OMEGA,IFB)
        call trap(emin,emax,INTE,fb,NZ,CP1-del,TEMP,OMEGA,IFB)
        call trap(emin,emax,INTE,fc,NZ,CP1,TEMP,OMEGA,IFB)
                  fd=(fa-fb)/(2d0*del)
                  CP2=CP1-fc/fd
c
                  if (abs(CP2-CP1) .lt. EPS)then                     
                     !write(50,*)TEMP/TCrit,CP2/(TCrit*BoltzK)
                     CPP=CP2
                     VAL2=0d0
                     exit
                  endif
                     if (intb .EQ. 10000)write(*,*)'failed at newtons'
                enddo  !end newton
c --------------------------     
                exit
               endif
c     
            enddo
c --- if can't find chemical potential
            if(JCHK .EQ. 0)then
            write(*,*)'can not find chemical potential'
            CPP=0d0
            endif
c ---
c --- check chemical potential
          call trap(emin,emax,INTE,AVAL,NZ,CPP,TEMP,OMEGA,IFB)
          if (IFB .EQ. 1)then
!          write(*,*)NZ-AVAL,'ferimion count CP'
          else
          call FU(CPP,TEMP,BVAL)
!          write(*,*)NZ-AVAL-BVAL,'excite'
!          write(*,*)BVAL,'ground'
!          write(*,*)NZ-AVAL,'boson count CP'
          endif
c          write(50,*)TEMP/TCRIT,(NZ-AVAL-BVAL)/NZ*100d0
c          write(60,*)TEMP/TCRIT,BVAL/NZ*100d0
c
c      close(50)
 
c      write(*,*)TEMP/TSP,CPP/(TSP*BoltzK)
      return
      end
c
c============================================================
c     subroutine trapezium integration
c=============================================================
      subroutine trap(EN1,ENN,INTE,VAL,N,CP,TEMP,OMEGA,IFB)
      use constant
      implicit real*8(A-H,O-Z)
!      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)
      External DE,GE
c      dimension XZ(9999),WZ(9999),WK(9999)
      dimension XZ(99999),WZ(99999),WK(99999)
c      data itimes/0/
c
c --- simpson
         if(mod(INTE,2).eq.0)then
            INTE=INTE-1
         endif
         Xstep=1.0d0/(INTE-1)
         Simp=4.0d0/3.0d0*Xstep
         do nn=2,INTE-1,2
            XZ(nn)=Xstep*(nn-1)
            WZ(nn)=Simp
         enddo
         Simp=2.0d0/3.0d0*Xstep
         do nn=1,INTE-2,2
            XZ(nn)=Xstep*(nn-1)
            WZ(nn)=Simp
         enddo
         XZ(1)=0.0d0
         XZ(INTE)=1.0d0
         WZ(1)=Xstep/3.0d0
         WZ(INTE)=Xstep/3.0d0
c
         AA=ENN-EN1
         BB=EN1
         do nn=1,INTE
            XZ(nn)=AA*XZ(nn)+BB
            WZ(nn)=AA*WZ(nn)
         enddo
c
      VAL=0d0
      if(IFB .EQ. 1)then
         do nn=1,INTE
            VAL=VAL+DE(XZ(nn),CP,TEMP,OMEGA)*WZ(nn)
         enddo
         VAL=N-VAL
      else
         do nn=1,INTE
            VAL=VAL+GE(XZ(nn),CP,TEMP,OMEGA)*WZ(nn)
         enddo
         call FU(CP,TEMP,BVAL)
c     write(*,*)BVAL,TEMP/TCRIT
         VAL=VAL+BVAL
         VAL=N-VAL
      endif
c --- trapezium integration 
c      CONST=1d0/(2d0*(hbar*OMEGA)**3)
c      ESTEP=(ENN-EN1)/(INTE-1)
c      VAL=DE(EN1,CP,TEMP)+DE(ENN,CP,TEMP)
c      VAL=0.5d0*VAL
c
c      do NUM=2,INTE-1
c         EN=EN1+(NUM-1)*ESTEP
c         VAL=VAL+DE(EN,CP,TEMP)
c      ENDDO
c
c      VAL=CONST*VAL*ESTEP
c      VAL=N-VAL
c ---
      return
      end
c========================================
c=========================================
c     function D(E)F(E)*1/CONST fermion
c=========================================
      function DE(EN,CP,TEMP,OMEGA)
      use constant
      implicit real*8(A-H,O-Z)
!      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)
c
      alpha=BoltzK*TEMP
      EE=DEXP(EN/alpha)
      EM=DEXP(-CP/alpha)
      DE=1d0/(EE*EM+1d0)
      CONST=EN**2/(2d0*(hbar*OMEGA)**3)
      DE=DE*CONST
c     
      return
      end
c===========================================c
c=========================================
c     function D(E)G(E)*1/CONST boson
c=========================================
      function GE(EN,CP,TEMP,OMEGA)
      use constant
      implicit real*8(A-H,O-Z)
!      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)
c
      alpha=BoltzK*TEMP
      EE=DEXP(EN/alpha)
      EM=DEXP(-CP/alpha)
      GE=1d0/(EE*EM-1d0)
      CONST=EN**2/(2d0*(hbar*OMEGA)**3)
c
      GE=GE*CONST
c     
      return
      end
c============================================    
c     FUGACITY Z
c============================================      
      subroutine FU(CP,TEMP,BVAL)
      use constant
      implicit real*8(A-H,O-Z)
!      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)
c
      alpha=BoltzK*TEMP
      Z=DEXP(CP/alpha)
      BVAL=Z/(1d0-Z)

      return
      end
c      
c======================================
c     making random number
c======================================
      function prand(idum)
c --- old; cycle=~10^9
c      function prand(iseed)
c      implicit real*8(a-h,o-z)
c      data maxint/2147483647/
c      data mprod/48828125/
c      data rinvm/4.65661287524d-10/
c
c      if (iseed.le.1) then
c         iseed=487219165
c      endif
c      iseed=iseed*mprod
c      iseed=iand(iseed,maxint)
c      prand=iseed*rinvm
c      end
c --- new; from NR ran2.f
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      DOUBLE PRECISION prand,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      prand=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software oV)#S.




c=====================================
c      twobody 
c=====================================
      subroutine TWOB(XMK,XMRB,PARUX,PARUY,PARUZ,PARUPX,PARUPY,
     & PARUPZ,PARDX,PARDY,PARDZ,PARDPX,PARDPY,PARDPZ,TWX,TWY,TWZ,
     & TPX,TPY,TPZ,OMEGAK,OMEGARB,OEMGAKZ,OMEGARBZ,NX,NY,NZ)
      use constant
      implicit real*8(a-h,o-z)
c      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)

      TPX=PARUPX+PARDPX
      TPY=PARUPY+PARDPY
      TPZ=PARUPZ+PARDPZ
c      
      RMA=XMK+XMRB
c      
      TWX=(XMK*PARUX+XMRB*PARDX)/RMA
      TWY=(XMK*PARUY+XMRB*PARDY)/RMA
      TWZ=(XMK*PARUZ+XMRB*PARDZ)/RMA
c
c      OMEGA=(OMEGARB+OMEGAK)/2d0
c      OMEGAZ=(OMEGARBZ+OMEGAKZ)/2d0
      OMEGA=174d0*2d0*pi
      OMEGAZ=2.93d0*2d0*pi
c
      EX=TPX**2/2d0/RMA+RMA*(OMEGA*TWX)**2/2d0
      EY=TPY**2/2d0/RMA+RMA*(OMEGA*TWY)**2/2d0
      EZ=TPZ**2/2d0/RMA+RMA*(OMEGAZ*TWZ)**2/2d0

      NX=IDINT(EX/(hbar*OMEGA))
      if (EX/(hbar*OMEGA)-NX .GE. 0.5d0)then
         NX=NX+1
      endif
      NY=IDINT(EY/(hbar*OMEGA))
      if (EY/(hbar*OMEGA)-NY .GE. 0.5d0)then
         NY=NY+1
      endif
      NZ=IDINT(EZ/(hbar*OMEGAZ))
      if (EZ/(hbar*OMEGAZ)-NZ .GE. 0.5d0)then
         NZ=NZ+1
      endif

!      NX=EX/(hbar*OMEGA)
!      if (EX/(hbar*OMEGA)-NX .GE. 0.5d0)then
!         NX=NX+1
!     endif
!      NY=EY/(hbar*OMEGA)
!      if (EY/(hbar*OMEGA)-NY .GE. 0.5d0)then
!         NY=NY+1
!      endif 
c     
c      NZ=EZ/(hbar*OMEGA)
c      if (EZ/(hbar*OMEGA)-NZ .GE. 0.5d0)then
c         NZ=NZ+1
c      endif
c              
!      NZ=EZ/(hbar*OMEGAZ)
!      if (EZ/(hbar*OMEGAZ)-NZ .GE. 0.5d0)then
!         NZ=NZ+1
!      endif         
c
      return
      end
c
c=================================================
c
c================================================
c     Pauli blocking (this calculate is suspcious
c                    because we don't know real frequency of KRB)
c===============================================
      subroutine PAULIB(NX,NY,NZ,LCOUNT,ICHK)
      use constant
      implicit real*8(a-h,o-z)
c      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)
      dimension NX(LCOUNT),NY(LCOUNT),NZ(LCOUNT)
c      
      ICHK=1
      do i=1,LCOUNT-1
         NA=ABS(NX(LCOUNT)-NX(i))
         NB=ABS(NY(LCOUNT)-NY(i))
         NC=ABS(NZ(LCOUNT)-NZ(i))
         if (NA+NB+NC .EQ. 0)then
         ICHK=0
         exit
         endif
      enddo      
      
      return
      end
c====================================================

c============================================================
c     check atoms inside the BEC ;Thomas-Fermi radius
c===========================================================
      subroutine INTFR(X,Y,Z,PX,PY,PZ,XMK,XMRB,WKR,WKZ,WRBR,WRBZ,CPB,
     &  dely,IC)
      use constant
      implicit real*8(a-h,o-z)
c
!      parameter(pi=3.141592654d0)   
!      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)
!      parameter(h=2d0*pi*hbar)
c
      AA=1d0
      OMEGARB=(WRBR**2*WRBZ)**(1d0/3d0)
c --- Thomas Fermi
c      FE=CPB
c      FE=(PX**2+PY**2+PZ**2)/(2d0*XMRB)
c      XMAX=(AA*FE/(XMRB*WRBR**2))**(0.5d0)
c      PXMAX=(AA*FE*XMRB)**(0.5d0)
c      YMAX=(AA*FE/(XMRB*WRBR**2))**(0.5d0)
c      PYMAX=(AA*FE*XMRB)**(0.5d0)
c      ZMAX=(AA*FE/(XMRB*WRBZ**2))**(0.5d0)
c      PZMAX=(AA*FE*XMRB)**(0.5d0)
c --- wigner
      XMAX=(hbar/XMRB/OMEGARB)**(0.5d0)
      PXMAX=(xmrb*OMEGARB*hbar)**(0.5d0)
      YMAX=XMAX
      PYMAX=PXMAX
      ZMAX=XMAX
      PZMAX=PXMAX      
c
      AX=1d0/XMAX
      AY=1d0/YMAX
      AZ=1d0/ZMAX
      APX=1d0/PXMAX
      APY=1d0/PYMAX
      APZ=1d0/PZMAX
c     
c      YMAX=YMAX+dely
c      write(*,*)'intfr max',XMAX,YMAX,ZMAX,PXMAX,PYMAX,PZMAX
c      
      TX=X*AX
      TY=Y*AY
      TZ=Z*AZ
      TPX=PX*APX
      TPY=PY*APY
      TPZ=PZ*APZ
c
      TDY=dble(dely*AY)
c      write(*,*)'intrf unnor',X,Y,Z,PX,PY,PZ
c      write(*,*)'intfr norm',TX,TY,TZ,TPX,TPY,TPZ,TDY
c      
c      E=(PX**2+PY**2+PZ**2)/(2d0*XMK)
c      E=E+XMK*(OMEGA**2*(X**2+Y**2)+OMEGAZ**2*Z**2)/2d0
c
c      D=(-TX)**2+(TDY-TY)**2+(-TZ)**2+
c     &     (-TPX)**2+(-TPY)**2+(-TPZ)**2
c
      D=(-TX)**2+(TDY-TY)**2+(-TZ)**2   
c      write(*,*)'D',D-1d0
      if(D-1d0 .le. 0d0)then
         IC=1
      else
         IC=0
      endif
      
      return
      end
c==========================================================
c
c==============================================================
c     make load file
c==============================================================
      subroutine makeload(TEMP,NTLOOP)
      implicit real*8(a-h,o-z)
c
      ITEMP=TEMP*100
      write(*,*)'TEMP',ITEMP
      if (NTLOOP .EQ. 1)then
         write(7100,*)'set xr[0:2.0]'
         write(7100,*)'set yr[0:30]'
         write(7100,*)'set term postsc eps enhan col'
         write(7100,*)'set output "TDCR.EPS"'
         write(7100,*)'p "CR" wi lines'
      endif
c
      return
      end
c====================================================
c  new write check distribution at Z axis
c ND(N,N,1=x or 2=y or 3=z) 
c check R and P distribution
c===================================================
      subroutine chkzdis(NK,PARU,NDISU,XXI,PPI,RRI,RPK)
      use constant
      implicit real*8(a-h,o-z)
      dimension PARU(NK,7),NDISU(100,100,3)
      dimension XXI(100),PPI(100),RRI(100,2),RPK(100,100)
!c
      do i=1,NK
c --- x direction
         if (PARU(i,7) .eq. 1)then
         else
         do lnum=1,3
            mnum=lnum+3
            do j=1,100
               if (PARU(i,lnum) .LE. XXI(j))then
                  do k=1,100
                     if (PARU(i,mnum) .LE. PPI(k))then
                        NDISU(j,k,lnum)=NDISU(j,k,lnum)+1
                        exit
                     endif
                  enddo
                  exit
               endif            
            enddo 
         enddo
c --- end
      endif
      enddo

c      do i=1,NK
c  --- rp direction
c         if (PARU(i,7) .eq. 1)then
c         else
c            RA=DSQRT(PARU(i,1)**2+PARU(i,2)**2+PARU(i,3)**2)
c            PA=DSQRT(PARU(i,4)**2+PARU(i,5)**2+PARU(i,6)**2)
c            do j=1,100
c               if (RA .le. RRI(j,1))then
c                  do k=1,100
c                     if (PA .le. RRI(k,2))then
c                        RPK(j,k)=RPK(j,k)+1
c                        exit
c                     endif
c                  enddo
c                  exit
c               endif
c            enddo
c  --- end
c         endif
c      enddo
c ---
      do i=1,NK
c  --- x-y direction (integrate p,z)
         if (PARU(i,7) .eq. 1)then
         else
            RA=PARU(i,1)
            PA=PARU(i,2)
            do j=1,100
               if (RA .le. RRI(j,1))then
                  do k=1,100
                     if (PA .le. RRI(k,2))then
                        RPK(j,k)=RPK(j,k)+1
                        exit
                     endif
                  enddo
                  exit
               endif
            enddo
c  --- end
         endif
      enddo
c
      return
      end         
c=====================================================

c=====================================================
      subroutine SETA(NA,XM,WR,WZ,WA,IP,ISO,TSP,PAR)
      use constant
      implicit real*8(a-h,o-z)
      dimension PAR(5)
      
      XM=XM*1.67262158d-27
      WR=WR*2d0*pi
      WZ=WZ*2d0*pi
      WA=(WR**2d0*WZ)**(1.0d0/3.0d0) !avarage of omega

      if (IP .EQ. 1)then
         TSP=(6d0*NA)**(1.0d0/3.0d0)*hbar*WA/BoltzK !Fermi Temp
      else
         TSP=0.94d0*hbar*WA*NA**(1d0/3d0)/BoltzK !Critical Temp
      endif
      
      PAR(1)=XM
      if (ISO .eq. 1)then
         PAR(2)=WA
      else
         PAR(2)=WR
      endif
      if (ISO .eq. 1)then
         PAR(3)=WA
      else
         PAR(3)=WZ
      endif
      PAR(4)=WA
      PAR(5)=TSP

      return
      end
c=============================================
      subroutine SETC(NA,NZE,TSP,TEMP,OMEGA,IP,CP)
      use constant
      implicit real*8(a-h,o-z)
      
      if (IP .eq. 1)then
      call CPBF(NA,TEMP,TSP,OMEGA,CP,1)
      write(*,*)CP/(TSP*BoltzK),'CPF/FE'
      NZE=0
      else
      call CPBF(NA,TEMP,TSP,OMEGA,CP,-1)
      write(*,*)CP/(TSP*BoltzK),'CPB/CE'
      call FU(CP,TEMP,BVAL)      
      if(BVAL-INT(BVAL) .GE. 0.5d0)then
         NZE=IDINT(BVAL)+1
         else
         NZE=IDINT(BVAL)
      endif
      NEX=NA-NZE
      write(*,*)'N ground',NZE
      write(*,*)'N excite',NEX
      endif
      
      return
      end
c==============================================
c==============================================
      subroutine CRP(RM,PM,RRI,N)
      implicit real*8(a-h,o-z)
      dimension RRI(N,2)
c ---
      RA=RM/dble(N-1)
      PA=PM/dble(N-1)
      do i=1,N
      RRI(i,1)=RA*dble(i-1)     
      RRI(i,2)=PA*dble(i-1) 
      enddo
c      write(*,*)'CRP',RA,PA

      return
      end
c========================================
c==================================================== 
c check R and P distribution for molecule
c===================================================
      subroutine moldis(PMOL,RRI,PRDEL,NK)
      use constant
      implicit real*8(a-h,o-z)
      parameter (NA=500)
      dimension RRI(NA,2),PMOL(NK,3),PRDEL(NA,NA)
!c
      do i=1,NK
c --- rp direction
         if (PMOL(i,3) .eq. 0)then
c            write(*,*)'OUT',i
            exit
         endif
         do j=1,NA
            if (PMOL(i,1) .le. RRI(j,1))then
               do k=1,NA
                  if (PMOL(i,2) .le. RRI(k,2))then
                     PRDEL(j,k)=PRDEL(j,k)+1
                     exit
                  endif
               enddo
               exit
            endif
         enddo
c ---
      enddo

      return
      end
c================================================
