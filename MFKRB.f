      program molfraction
      use constant
      implicit real*8(a-h,o-z)
      dimension PAR1(5),PAR2(5),CM(4),IMODE(2),RPMEAN(4)
      NAMELIST /INF/NA,NB,IFA,IFB,XMA,XMB,WAR,WAZ,WBR,WBZ,GANMA,
     &     NTLOOP,NTDIV,NTMAX,TMAX,TMIN,
     &     MODE1,MODE2,TEMP,ISEED,ISEED2,ISO,ISAG
c
      open(1,file='input')
      read(1,INF)
      write(*,INF)
      close(1)

c --- setting
      call SETA(NA,XMA,WAR,WAZ,WAA,IFA,ISO,TSP1,PAR1)
      call SETA(NB,XMB,WBR,WBZ,WBA,IFB,ISO,TSP2,PAR2)
      IMODE(1)=MODE1
      IMODE(2)=MODE2
c --- 
      AMASS=(XMA*XMB)/(XMA+XMB) !reduced mass
c      AMASS=(XMK+XMRB) !molecule mass
c --- gravitational sag
      if (ISAG .eq. 1)dely=grav*(1d0/WBR**2-1d0/WAR**2) !use default
      if (ISAG .eq. 2)dely=grav*(1d0/WAA**2-1d0/WBA**2)
      if (ISAG .eq. 0)dely=0d0
c --- check
 !     write(*,*)'NK',NA
 !     write(*,*)'NRB',NB
 !     write(*,*)'K MASS',PAR1(1)
 !     write(*,*)'RB MASS',PAR2(1)      
 !     write(*,*)'TFERMI',PAR1(5)
 !     write(*,*)'TCRITICAL',PAR2(5)
 !     write(*,*)'WAR',PAR1(2)
 !     write(*,*)'WAZ',PAR1(3)
 !    write(*,*)'WAA',PAR1(4)
 !     write(*,*)'WBR',PAR2(2)
 !     write(*,*)'WBZ',PAR2(3)
 !     write(*,*)'WBA',PAR2(4)      
 !     write(*,*)'delta y',dely
 !     write(*,*)'AMASS',AMASS
c ---
      if (NTLOOP .EQ. 1)write(*,*)'Temperature Loop'
      if (NTLOOP .EQ. 0)write(*,*)'z distribution'
      if (NTLOOP .EQ. 2)write(*,*)'number loop'
      if (NTLOOP .EQ. 3)write(*,*)'gamma loop'
c    
c --- time evaluation 
      call cpu_time(altime1)
      TEMPS=TEMP*TSP2
c --- start temparature cycle
      if (NTLOOP .EQ. 1)then
      IMODE(1)=1
      TSTEP=(TMIN-TMAX)/(NTDIV-1)
      do ITEMP=1,NTDIV !1
         TEMPS=(TMAX+TSTEP*(ITEMP-1))*TSP2
         if ((ITEMP .eq. 1).and.(NTDIV .eq. 1))TEMPS=TEMP*TSP2

         call CRKRB(NA,NB,PAR1,PAR2,IFA,IFB,IMODE,
     &       TEMPS,DELY,NTLOOP,AMASS,NTMAX,GANMA,CM,
     &       PSPDK,PSPDRB,PSPDBE,PSPDTH,ISEED,ISEED2,DRSUM,DPSUM,RPMEAN)
c
         TT=TEMPS/TSP2
         write(600,*)TT,CM(1)
         write(601,*)TT,PSPDK
         write(602,*)TT,PSPDRB
         write(603,*)TT,PSPDBE
         write(604,*)TT,PSPDTH
         write(605,*)TT,CM(2)
         write(606,*)TT,CM(3)
         write(607,*)TT,CM(4)
         write(608,*)TT,DRSUM/abo
         write(609,*)TT,DPSUM/hbar
         write(112,*)TT,RPMEAN(1),RPMEAN(3)
         write(113,*)TT,RPMEAN(2),RPMEAN(4)
         write(*,*)'check times',(DRSUM*DPSUM)/(GANMA*h)
c
      enddo !1
      endif      
c ---
      if (NTLOOP .EQ. 0)then
c
         call CRKRB(NA,NB,PAR1,PAR2,IFA,IFB,IMODE,
     &       TEMPS,DELY,NTLOOP,AMASS,NTMAX,GANMA,CM,
     &       PSPDK,PSPDRB,PSPDBE,PSPDTH,ISEED,ISEED2,DRSUM,DPSUM,RPMEAN)
c         
         write(700,*)TEMPS/TSP2,CM(1)
         write(608,*)TEMPS/TSP2,DRSUM/abo
         write(609,*)TEMPS/TSP2,DPSUM/hbar
         write(112,*)TEMPS/TSP2,RPMEAN(1),RPMEAN(3)
         write(113,*)TEMPS/TSP2,RPMEAN(2),RPMEAN(4)
         
         write(*,*)'check times',(DRSUM*DPSUM)/(GANMA*h)         
c     
      endif
c ---
      
      if (NTLOOP .EQ. 2)then
         write(720,*)
         open(600,position='append')
         open(607,position='append')
         open(608,position='append')         
         open(605,position='append')
         open(606,position='append')         
c
         call CRKRB(NA,NB,PAR1,PAR2,IFA,IFB,IMODE,
     &       TEMPS,DELY,NTLOOP,AMASS,NTMAX,GANMA,CM,
     &       PSPDK,PSPDRB,PSPDBE,PSPDTH,ISEED,ISEED2,DRSUM,DPSUM,RPMEAN)
c        
         TT=TEMPS/TSP2
         c1=dble(NA)/dble(NA+NB)
         c2=dble(NA)/100d0
         write(600,*)NA,NB,CM(1)
         write(607,*)NA,NB,CM(1)*c1
         write(608,*)NA,NB,CM(1)*c2
         write(605,*)NA,NB,CM(2)*c2
         write(606,*)NA,NB,CM(2)*c2
         close(600)
         close(605)
         close(606)
         close(607)
         close(608)
c     
      endif
c ---
!     loop for gamma
      if (NTLOOP .EQ. 3)then
      open (609,file='GCR',position='append')
      IMODE(1)=1
      TSTEP=(TMIN-TMAX)/(NTDIV-1)
      do ITEMP=1,NTDIV !1
         TEMPS=(TMAX+TSTEP*(ITEMP-1))*TSP2
         if ((ITEMP .eq. 1).and.(NTDIV .eq. 1))TEMPS=TEMP*TSP2

         call CRKRB(NA,NB,PAR1,PAR2,IFA,IFB,IMODE,
     &       TEMPS,DELY,NTLOOP,AMASS,NTMAX,GANMA,CM,
     &       PSPDK,PSPDRB,PSPDBE,PSPDTH,ISEED,ISEED2,DRSUM,DPSUM,RPMEAN)
c
         TT=TEMPS/TSP2
         write(609,*)GANMA,TT,CM(1)

      enddo !1
      write(609,*)
      close(609)
      endif       
c
      call makeload(TEMPS/TSP2,NTLOOP)
c --- time evaluation
      call cpu_time(altime2)
      altime=altime2-altime1
      write(*,*)'ALL TIME',altime
c ---
      stop
      end
c
c=============================================================
      subroutine CRKRB(NK,NRB,PAR1,PAR2,IFA,IFB,IMO,
     & TEMP,DELY,NTLOOP,AMASS,NTMAX,GANMA,CMA,
     & PSPDK,PSPDRB,PSPDBE,PSPDTH,ISEED,ISEED2,DRSUM,DPSUM,RPMEAN)
      use constant
      implicit real*8(a-h,o-z)
      real*8, allocatable::PARZU(:,:),PAREU(:,:)
      real*8, allocatable::PARZD(:,:),PARED(:,:)      
      real*8, allocatable::PARU(:,:),PARD(:,:)
      real*8, allocatable::PMOL(:,:),MOLN(:,:)
      real*8, allocatable::ZMOL(:,:)
      real*8, allocatable::DIVE1(:,:),DIVE2(:,:),DIVE3(:,:)
      real*8, allocatable::DIVE4(:,:)
      integer*4, allocatable::IDNUM1(:),IDNUM2(:)
      integer*4, allocatable::IDNUM3(:),IDNUM4(:)      
      parameter (IDN=300)
      parameter(IDIM=3)!NTMAX:number of trial !IDIM:dimension
      dimension PAR1(5),PAR2(5),CMA(4),PSPD(4),IMO(2)
      dimension NDISU(100,100,3),NDISD(100,100,3),MDIS(100,100,3)
      dimension NDAU(100,100,3),NDAD(100,100,3)
      dimension LTH(100,100,3),LBE(100,100,3)
      dimension XXL(100),PPL(100),ENDIV(100)
      dimension ENDU(100),ENDD(100)
      dimension SUMNU(3),SUMND(3),SUMBE(3),SUMTH(3),SUMOL(3)
      dimension SUMAU(3),SUMAD(3)      
      dimension TWX(NK),TWY(NK),TWZ(NK),
     &  TPX(NK),TPY(NK),TPZ(NK),
     &  NTX(NK),NTY(NK),NTZ(NK),AVG(6),VAR(6)
      dimension RRI(500,2),RRS(100,2),RPK(100,100),RPRB(100,100)
      dimension ARPK(100,100),ARPRB(100,100)
      dimension PRDEL(500,500),PRDUM(100,100)
      dimension SETMAX1(6),SETMAX2(6),SETMAX3(6),SETMAX4(6)
      dimension EDIM1(1000),EDIM2(1000)
      dimension RPMEAN(4)
c      dimension DIVE1(2,IDN),DIVE2(2,IDN),DIVE3(2,IDN),DIVE4(2,IDN)
c      dimension IDNUM1(IDN),IDNUM2(IDN),IDNUM3(IDN),IDNUM4(IDN)
c=============================================================
c --- substitute
      XMK=PAR1(1)
      XMRB=PAR2(1)
      WKR=PAR1(2)
      WKZ=PAR1(3)
      OMEGAK=PAR1(4)
      WRBR=PAR2(2)
      WRBZ=PAR2(3)
      OMEGARB=PAR2(4)
      TFERMI=PAR1(5)
      TCRIT=PAR2(5)
c --- check
!      write(*,*)'NA',NK
!      write(*,*)'NB',NRB
!      write(*,*)'MASS 1',XMK
!      write(*,*)'MASS 2',XMRB 
!      write(*,*)'IFA',IFA     
!      write(*,*)'TSP1',TFERMI
!      write(*,*)'IFB',IFB
!      write(*,*)'TSP2',TCRIT
!      write(*,*)'WAR',WKR
!      write(*,*)'WAZ',WKZ
!      write(*,*)'WAA',OMEGAK
!      write(*,*)'WBR',WRBR
!      write(*,*)'WBZ',WRBZ
!      write(*,*)'WBA',OMEGARB      
!      write(*,*)'delta y',dely
!      write(*,*)'REDUCED MASS',AMASS
!      write(*,*)'ISEED initial',ISEED
!      write(*,*)'TEMP',TEMP
c --- allocate
      allocate(PARU(NK,7))
      allocate(PARD(NRB,7))
      allocate(PMOL(NK,3))
      allocate(MOLN(NK,4))
      allocate(ZMOL(NK,7))
c ---      
      NDISU=0
      NDISD=0
      LBE=0
      LTH=0
      JBE=0
      JTH=0
      AAA=100d0/dble(NK)
      BBB=100d0/dble(NRB)
      DRSUM=0d0
      DPSUM=0d0
c --- volume of particles 
      RMEAN1=0d0
      PMEAN1=0d0
      RMEAN2=0d0
      PMEAN2=0d0
c --- lattice
      if(NTLOOP .EQ. 0)then
c      XK=3d-5
c      PK=10d-28
c --- twice  
      XK=6d-5
      PK=50d-28
      RM=6d-5
      PM=10d-28
      RM2=6d-5
      PM2=25d-28
c      ENMAX=1d-29      
c --- four times
c      XK=12d-5
c      PK=20d-28
      ENMAX=5d-29
c ---
      AAA=100d0/dble(NK)
      BBB=100d0/dble(NRB)
      call CXP(XK,PK,XXL,PPL,100) 
      call CEN(ENMAX,ENDIV,100)
      call CRP(RM,PM,RRI,500)
      call CRP2(RM2,PM2,RRS,100)
      endif
c ---
      CMAVG=0d0
      PBAVG=0d0
      NPESUM=0d0
c --- FIX TEMP
      DUMTEMP=TEMP/TFERMI
      TEMPA=TEMP
c      if (DUMTEMP .le. 2.2d0)TEMPA=2.2d0*TFERMI !limit temperature
      if (DUMTEMP .le. 0d0)TEMPA=2.2d0*TFERMI 
c --- chemical potential
      call  SETC(NK,NZE1,TFERMI,TEMPA,OMEGAK,IFA,CPF)
      NEX1=NK-NZE1
      call  SETC(NRB,NZE2,TCRIT,TEMP,OMEGARB,IFB,CPB)
      NEX2=NRB-NZE2
c      NZE2=NZE2
c ---
      if (NZE2 .GE. 1)then
         NBCHK=1
      else
         NBCHK=0
      endif      
c --- allocate
      allocate(PARZU(NZE1,7))
      allocate(PAREU(NEX1,7))
      allocate(PARZD(NZE2,7))
      allocate(PARED(NEX2,7))
c ---
      NDCOUNT=0
      ENMAX1=35d0
      ENMAX2=35d0
      IDE1=220
      IDE2=315
      write(*,*)'before',ENMAX1,ENMAX2
      call setdiv(IFA,ENMAX1,CPF,TEMPA,OMEGAK,IDN1,EDIM1,IDE1)
      allocate(DIVE1(2,IDN1),DIVE2(2,IDN1),IDNUM1(IDN1),IDNUM2(IDN1))
      call setdiv(IFB,ENMAX2,CPB,TEMP,OMEGARB,IDN2,EDIM2,IDE2)
      allocate(DIVE3(2,IDN2),DIVE4(2,IDN2),IDNUM3(IDN2),IDNUM4(IDN2)) 
      write(*,*)'after',ENMAX1,ENMAX2
c --- set max status
      call MAXSET(SETMAX1,XMK,WKR,WKZ,TEMPA,TFERMI,IFA,CPF,
     &       DIVE1,IDNUM1,IDN1,EDIM1,ENMAX1)
      if (IFA .eq. -1)then
         call MAXSET(SETMAX2,XMK,WKR,WKZ,TEMPA,TFERMI,0,CPF,
     &        DIVE2,IDNUM2,IDN1,EDIM1,ENMAX1)
      endif
      call MAXSET(SETMAX3,XMRB,WRBR,WRBZ,TEMP,TCRIT,IFB,CPB,
     &       DIVE3,IDNUM3,IDN2,EDIM2,ENMAX2)
      if (IFB .eq. -1)then
         call MAXSET(SETMAX4,XMRB,WRBR,WRBZ,TEMP,TCRIT,0,CPB,
     &        DIVE4,IDNUM4,IDN2,EDIM2,ENMAX2)  
      endif      
c      do ich=1,IDN2
c         write(*,*)DIVE3(:,ich),IDNUM3(ich)
c      enddo
c ---
c
      do NT=1,NTMAX !2
c ---
         write(*,*)'ISEED before',ISEED
         PARU=0
         PARD=0
         PAREU=0
         PARZU=0
         PARED=0
         PARZD=0
         PMOL=0
         ZMOL=0
         ZMOL(:,7)=1
         ISB=NT*192
         do ISA=1,ISB
            dum=prand(ISEED)
            dum=prand(ISEED2)
         enddo

         write(*,*)'Times',NT
!         write(*,*)'ISEED before',ISEED
c --- evaluating time 1
         call cpu_time(time1)     !  Startup
c --- particle 1
      if (IFA .eq. -1)then
      call fdist(NZE1,TEMPA,TFERMI,WKR,WKZ,XMK,PARZU,CPF,ISEED,NPE1,0
     & ,SETMAX2,DIVE2,IDNUM2,IDN1)
      call fdist(NEX1,TEMPA,TFERMI,WKR,WKZ,XMK,PAREU,CPF,ISEED,NPE1,-1
     & ,SETMAX1,DIVE1,IDNUM1,IDN1)
      do i=1,NEX1
         PARU(i,1)=PAREU(i,1)
         PARU(i,2)=PAREU(i,2)
         PARU(i,3)=PAREU(i,3)
         PARU(i,4)=PAREU(i,4)
         PARU(i,5)=PAREU(i,5)
         PARU(i,6)=PAREU(i,6)
         PARU(i,7)=PAREU(i,7)        
      enddo
      if (NZE1 .ge. 1)then
      do i=1,NZE1
         mm=i+NEX1
         PARU(mm,1)=PARZU(i,1)
         PARU(mm,2)=PARZU(i,2)
         PARU(mm,3)=PARZU(i,3)
         PARU(mm,4)=PARZU(i,4)
         PARU(mm,5)=PARZU(i,5)
         PARU(mm,6)=PARZU(i,6)
         PARU(mm,7)=PARZU(i,7)           
      enddo
      endif
      else
c         write(*,*)'coco'
      call fdist(NEX1,TEMPA,TFERMI,WKR,WKZ,XMK,PARU,CPF,ISEED,NPE1,1
     & ,SETMAX1,DIVE1,IDNUM1,IDN1)
      endif
      write(*,*)
      call cpu_time(FETIME)
c --- particle 2
      if (IFB .EQ. -1)then
      call fdist(NZE2,TEMP,TCRIT,WRBR,WRBZ,XMRB,PARZD,CPB,ISEED2,NPE2,0
     & ,SETMAX4,DIVE4,IDNUM4,IDN2)
      call fdist(NEX2,TEMP,TCRIT,WRBR,WRBZ,XMRB,PARED,CPB,ISEED2,NPE2,-1
     & ,SETMAX3,DIVE3,IDNUM3,IDN2)
      do i=1,NEX2
         PARD(i,1)=PARED(i,1)
         PARD(i,2)=PARED(i,2)+dely
         PARD(i,3)=PARED(i,3)
         PARD(i,4)=PARED(i,4)
         PARD(i,5)=PARED(i,5)
         PARD(i,6)=PARED(i,6)
         PARD(i,7)=PARED(i,7)        
      enddo
      if (NZE2 .GE. 1)then
      do i=1,NZE2
         mm=i+NEX2
         PARD(mm,1)=PARZD(i,1)
         PARD(mm,2)=PARZD(i,2)+dely
         PARD(mm,3)=PARZD(i,3)
         PARD(mm,4)=PARZD(i,4)
         PARD(mm,5)=PARZD(i,5)
         PARD(mm,6)=PARZD(i,6)
         PARD(mm,7)=PARZD(i,7)           
      enddo
      endif
      else
      call fdist(NEX2,TEMP,TCRIT,WRBR,WRBZ,XMRB,PARD,CPB,ISEED2,NPE2,1
     & ,SETMAX3,DIVE3,IDNUM3,IDN2)
      do i=1,NEX2
         PARD(i,2)=PARD(i,2)+dely
      enddo
      endif
      call cpu_time(BOTIME)
c --- raw data of RB distribution
      write(*,*)''
c --- evaluating time 2
         call cpu_time(time2)
         EPTFE=FETIME-time1
         EPTBO=BOTIME-FETIME
         elapt=time2-time1
         write(*,*)'elapsed:',elapt,'PAR1',EPTFE,'PAR2',EPTBO
c ---
         write(*,*)'ISEED after',ISEED
      write(*,*)''

      NPESUM=NPESUM+NPE1
c --- K & RB z-pz distribution
      if (NTLOOP .EQ. 0)then
      call chkzdis(NK,PARU,NDISU,XXL,PPL,RRS,RPK)
      call chkzdis(NRB,PARD,NDISD,XXL,PPL,RRS,RPRB)
      call chkENdis(NK,PARU,WKR,WKZ,XMK,TFERMI,ENDU,ENDIV,dely,0)
      call chkENdis(NRB,PARD,WRBR,WRBZ,XMRB,TCRIT,ENDD,ENDIV,dely,1)
      endif

c --- volume
      call varia(PARU,NK,AVG,VAR)
      RMDUM=DSQRT(AVG(1)**2+AVG(2)**2+AVG(3)**2)
      PMDUM=DSQRT(AVG(4)**2+AVG(5)**2+AVG(6)**2)
      RMEAN1=RMEAN1+RMDUM
      PMEAN1=PMEAN1+PMDUM
      call varia(PARD,NRB,AVG,VAR)
      RMDUM=DSQRT(AVG(1)**2+AVG(2)**2+AVG(3)**2)
      PMDUM=DSQRT(AVG(4)**2+AVG(5)**2+AVG(6)**2)
      RMEAN2=RMEAN2+RMDUM
      PMEAN2=PMEAN2+PMDUM
c ---

!      write(*,*)'sum ENDIS',sum(ENDU),sum(ENDD)
!      write(*,*)'sum repack',sum(NDISD)
!      write(*,*)'sum raw',sum(LBE)+sum(LTH)
c      do i=1,100
c         do j=1,100
c            write(2222,*)NDISD(i,j)
c            write(2223,*)LBE(i,j)+LTH(i,j)
c         enddo
c      enddo
c ---
c --- start calculate conversion rate
c     when Distribution Check on do not calculate conversion rate
c 
c --- cputime
      call cpu_time(time3)
c --- 
c      if(NTLOOP .EQ. 1)then  ! conversion rate caluculation
      if(IMO(1) .NE. 0)then  ! conversion rate caluculation
!         write(*,*)'in'
c --
      LEXCOUNT=0
      LCOUNT=0
      do i=1,NK !3
         if (PARU(i,7) .eq. 1)then
            write(*,*)i,'stopped count'
            exit
         endif
!c
c
c --- substitute
            P1=PARU(i,1)
            P2=PARU(i,2)
            P3=PARU(i,3)
            P4=PARU(i,4)/XMK
            P5=PARU(i,5)/XMK
            P6=PARU(i,6)/XMK
c ---
         do j=1,NRB !4 !do j=1,NRB
c --- inside BEC
!            if (TEMP .LE. TCRIT)then
!c            if(NBCHK .EQ. 1)then
!c               write(*,*)'enter intfr'
!               call INTFR(PARU(i,1),PARU(i,2),PARU(i,3),PARU(i,4),
!     & PARU(i,5),PARU(i,6),XMK,XMRB,WKR,WKZ,WRBR,WRBZ,CPB,dely,ICC)
!               if(ICC .EQ. 1)then
c --- count K inside BEC
!                  JBE=JBE+1
c ---                  
!c                  write(*,*)'inside BEC'
!c                  PARU(i,7)=0
!c                  LCOUNT=LCOUNT+1
!c                  write(4121,*)PARU(i,1),PARU(i,2),PARU(i,3)
!c                  write(4122,*)PARU(i,4),PARU(i,5),PARU(i,6)
!c                  
!c                  exit
!               endif
!            endif
            
            if (PARD(j,7) .eq. 0)then !partcile already form molecule
               Q4=PARD(j,4)/XMRB
               Q5=PARD(j,5)/XMRB
               Q6=PARD(j,6)/XMRB
               DXP=(P1-PARD(j,1))
               DYP=(P2-PARD(j,2))
               DZP=(P3-PARD(j,3))
               DPX=(P4-Q4)*AMASS
               DPY=(P5-Q5)*AMASS
               DPZ=(P6-Q6)*AMASS
               DR=DSQRT(DXP**2+DYP**2+DZP**2)
c               DP=DSQRT(DPX**2+DPY**2+DPZ**2)*AMASS
c --- angular momentum
               DP=DSQRT(DPX**2+DPY**2+DPZ**2)
c               DA=(DYP*DPZ-DZP*DPY)**2
c               DB=(DZP*DPX-DXP*DPZ)**2
c               DC=(DXP*DPY-DYP*DPX)**2
c               DD=DSQRT(DA+DB+DC)
c ---
c               if(DABS(DR*DP) .LT. ganma*h)then !close in phase space
c --- according to olsen's thesis
               if(DABS(DR*DP) .LT. ganma*h)then !close in phase space
c --- angluar momentum criterion
c               if(DD .LT. ganma*h)then   
                  PARU(i,7)=1
                  PARD(j,7)=1
!c                  MCOUNT=MCOUNT+1
                  LCOUNT=LCOUNT+1
                  PMOL(LCOUNT,1)=DR
                  PMOL(LCOUNT,2)=DP
                  PMOL(LCOUNT,3)=1
c                   
c --- prepare to check Pauli Blocking (P.B.) 
c --- frequency difference include
       call TWOB(XMK,XMRB,PARU(i,1),PARU(i,2),PARU(i,3),PARU(i,4),
     &  PARU(i,5),PARU(i,6),PARD(j,1),PARD(j,2),PARD(j,3),PARD(j,4)
     &  ,PARD(j,5),PARD(j,6),TWX(LCOUNT),TWY(LCOUNT),TWZ(LCOUNT),
     &  TPX(LCOUNT),TPY(LCOUNT),TPZ(LCOUNT),
     &  wkr,wrbr,wkz,wrbz,NTX(LCOUNT),NTY(LCOUNT),NTZ(LCOUNT))
c --- substitute
       ZMOL(LCOUNT,1)=TWX(LCOUNT)
       ZMOL(LCOUNT,2)=TWY(LCOUNT)
       ZMOL(LCOUNT,3)=TWZ(LCOUNT)
       ZMOL(LCOUNT,4)=TPX(LCOUNT)
       ZMOL(LCOUNT,5)=TPY(LCOUNT)
       ZMOL(LCOUNT,6)=TPZ(LCOUNT)
       ZMOL(LCOUNT,7)=0
       MOLN(LCOUNT,1)=NTX(LCOUNT)
       MOLN(LCOUNT,2)=NTY(LCOUNT)
       MOLN(LCOUNT,3)=NTZ(LCOUNT)
       MOLN(LCOUNT,4)=1
c --- check P.B.
!       if(LCOUNT .NE. 1)then !not equal to 1
!                call PAULIB(NTX,NTY,NTZ,LCOUNT,ICHK)
!                if (ICHK .EQ. 1)then !check
c                   NDCOUNT=NDCOUNT+1
c                  
c --- check combined Rb themal or BEC
                  if (j .LT. NEX2)then
                     JTH=JTH+1
                  else
                     JBE=JBE+1
                  endif
c ---
                  exit
!                else
!                   write(*,*)'PB',NTX(LCOUNT),NTY(LCOUNT),NTZ(LCOUNT)
!                PARU(i,7)=0
!                PARD(j,7)=0
!                LCOUNT=LCOUNT-1
!                LEXCOUNT=LEXCOUNT+1
!                endif !check
!             endif ! not equal to 1
c ---
!!           exit !if not check P.B.,'exit' needs here
c --- end check P.B.  
          endif !close in phase space
!c               
            endif !partcile already form molecule
         enddo !4
      enddo !3
c
      CM=dble(LCOUNT)/dble(i-1)*100d0
      CMAVG=CMAVG+CM
      write(*,*)'conversion',CM,'percent'
      write(111,*)CM
      PB=dble(LEXCOUNT)/dble(i-1)*100d0
      PBAVG=PBAVG+PB
      write(*,*)'pauli blocking',PB,'percent'      
      
c      write(*,*)dble(LCOUNT)/dble(i-1)*100d0,'percent'
c
c --- check molecular RP
      do imdum=1,LCOUNT
         DRDUM=PMOL(imdum,1)/dble(LCOUNT)
         DPDUM=PMOL(imdum,2)/dble(LCOUNT)
         DRSUM=DRSUM+DRDUM
         DPSUM=DPSUM+DPDUM
      enddo
c      write(*,*)DRDUM,DPDUM
      call moldis(PMOL,RRI,PRDEL,NK)
      call chkzdis(NK,ZMOL,MDIS,XXL,PPL,RRS,PRDUM)
      endif ! conversion rate calculation
c --- check distribution after molecule creation
      call chkzdis(NK,PARU,NDAU,XXL,PPL,RRS,ARPK)
      call chkzdis(NRB,PARD,NDAD,XXL,PPL,RRS,ARPRB)
c
c ---
      write(*,*)''
      enddo !2
c --- cputime
      call cpu_time(time4)
      elapt=time4-time3
      write(*,*)'elapsed2',elapt
      write(*,*)''
c --- MOLECULAR RPCHECK
c      write(*,*)DRSUM,DPSUM
      DRSUM=DRSUM/dble(NTMAX)
      DPSUM=DPSUM/dble(NTMAX)
c ---
      do ix=1,500
         RRA=RRI(ix,1)/abo
         SUMR=0d0
         do jp=1,500
            RRB=RRI(jp,2)/hbar
            write(9940,*)RRA,RRB,dble(PRDEL(ix,jp))/dble(NTMAX)
            SUMR=SUMR+dble(PRDEL(ix,jp))/dble(NTMAX)
         enddo
         write(9940,*)
         write(9941,*)RRA,SUMR
      enddo

      do ix=1,500
         RRA=RRI(ix,2)/hbar
         SUMP=0d0
         do jp=1,500
            SUMP=SUMP+dble(PRDEL(jp,ix))/dble(NTMAX)
         enddo
         write(9942,*)RRA,SUMP
      enddo

c ---
c
c      write(*,*)''
c      exavg=dble(NPESUM)/dble(NTMAX)/dble(NK)*100d0
!c      write(*,*)exavg,'exavg'
      if (NTLOOP .EQ. 1)then !write conversion rate
      write(*,*)'TEMP',TEMP/TCRIT
c      write(*,101)'excluded avg  ',exavg,' percent'
 101  format(a14,f8.3,a8)
!c      write(700,*)TEMP/TFERMI,dble((NPESUM/(NTMAX*2d0))/NK*100)
c
      write(*,101)'conv rate avg ',CMAVG/dble(NTMAX),' percent'
!c      write(600,*)TEMP/TCRIT,CMAVG/NTMAX
c
      write(*,101)'P.B. rate avg ',PBAVG/dble(NTMAX),' percent'
cc      write(900,*)TEMP/TCRIT,PBAVG/NTMAX
c      write(*,*)''
c
      endif !write conversion rate
c 
      CMAVG=CMAVG/dble(NTMAX)
      CMA(1)=CMAVG
      PBAVG=PBAVG/dble(NTMAX)
      CMA(4)=PBAVG
c
c --- end calcurate conversion rate
c============================================
c  zdist check
c============================================
      if(NTLOOP .EQ. 0)then
         DELX=(XXL(2)-XXL(1))/2d0
         DELP=(PPL(2)-PPL(1))/2d0
         DNT=dble(NTMAX)
         SUM1=0d0
         SUM2=0d0
         SUM3=0d0
c ---         
      do ix=1,100
         XXI=XXL(ix)-DELX
c         SUMUP=0
c         SUMDN=0
         SUMNU=0d0
         SUMND=0d0
         SUMOL=0d0
         SUMAU=0d0
         SUMAD=0d0
         do jp=1,100
c            PPI=PK/25d0*dble(jp)-2d0*PK
            PPI=PPL(jp)-DELP
c --- K
            SUMNU(1)=SUMNU(1)+dble(NDISU(ix,jp,1)/DNT)
            SUMNU(2)=SUMNU(2)+dble(NDISU(ix,jp,2)/DNT)
            SUMNU(3)=SUMNU(3)+dble(NDISU(ix,jp,3)/DNT)
c --- RB
            SUMND(1)=SUMND(1)+dble(NDISD(ix,jp,1)/DNT)
            SUMND(2)=SUMND(2)+dble(NDISD(ix,jp,2)/DNT)
            SUMND(3)=SUMND(3)+dble(NDISD(ix,jp,3)/DNT)
c --- MOLCULE
            SUMOL(1)=SUMOL(1)+dble(MDIS(ix,jp,1)/DNT)
            SUMOL(2)=SUMOL(2)+dble(MDIS(ix,jp,2)/DNT)
            SUMOL(3)=SUMOL(3)+dble(MDIS(ix,jp,3)/DNT)
c --- K AFTER
            SUMAU(1)=SUMAU(1)+dble(NDAU(ix,jp,1)/DNT)
            SUMAU(2)=SUMAU(2)+dble(NDAU(ix,jp,2)/DNT)
            SUMAU(3)=SUMAU(3)+dble(NDAU(ix,jp,3)/DNT)
c --- Rb AFTER
            SUMAD(1)=SUMAD(1)+dble(NDAD(ix,jp,1)/DNT)
            SUMAD(2)=SUMAD(2)+dble(NDAD(ix,jp,2)/DNT)
            SUMAD(3)=SUMAD(3)+dble(NDAD(ix,jp,3)/DNT)
            

         enddo

         write(4001,*)XXI,SUMNU(1)*AAA
         write(4002,*)XXI,SUMNU(2)*AAA
         write(4003,*)XXI,SUMNU(3)*AAA
         write(4101,*)XXI,SUMAU(1)*AAA
         write(4102,*)XXI,SUMAU(2)*AAA
         write(4103,*)XXI,SUMAU(3)*AAA
         write(5001,*)XXI,SUMND(1)*BBB
         write(5002,*)XXI,SUMND(2)*BBB
         write(5003,*)XXI,SUMND(3)*BBB
         write(5101,*)XXI,SUMAD(1)*BBB
         write(5102,*)XXI,SUMAD(2)*BBB
         write(5103,*)XXI,SUMAD(3)*BBB
         write(6001,*)XXI,SUMOL(1)
         write(6002,*)XXI,SUMOL(2)
         write(6003,*)XXI,SUMOL(3)         

      enddo
c --- write data
!      CKK=100d0/SUM2
!      CBO=100d0/SUM3
!      CTH=100d0/SUM1
!      do ix=1,100
!         XXI=XXL(ix)
!         do jp=1,100
!c            PPI=PK/25d0*dble(jp)-2d0*PK
!            PPI=PPL(jp)
!            call ENXP(XXI-dely,PPI,XMRB,OMEGARB,TEMP,CPB,dum,VAL1,-1)
!            call ENXP(XXI,PPI,XMK,OMEGAK,TEMP,CPF,dum,VAL2,1)
!            call ENXP(XXI-dely,PPI,XMRB,OMEGARB,TEMP,CPB,dum,VAL3,0)
!            write(9960,*)XXI,PPI,VAL1*CTH
!            write(9961,*)XXI,PPI,VAL2*CKK
!            write(9962,*)XXI,PPI,VAL3*CBO
!         enddo
!         write(9960,*)
!         write(9961,*)
!         write(9962,*)
!      enddo

c --- pz direction
c
      do jp=1,100
c         PPI=PK/25d0*dble(jp)-2d0*PK
         PPI=PPL(jp)-DELP
c         SUMUP=0
c         SUMDN=0
         SUMNU=0d0
         SUMND=0d0
         SUMOL=0d0
         SUMAU=0d0
         SUMAD=0d0
         do ix=1,100
c --- K
            SUMNU(1)=SUMNU(1)+dble(NDISU(ix,jp,1)/DNT)
            SUMNU(2)=SUMNU(2)+dble(NDISU(ix,jp,2)/DNT)
            SUMNU(3)=SUMNU(3)+dble(NDISU(ix,jp,3)/DNT)
c --- RB
            SUMND(1)=SUMND(1)+dble(NDISD(ix,jp,1)/DNT)
            SUMND(2)=SUMND(2)+dble(NDISD(ix,jp,2)/DNT)
            SUMND(3)=SUMND(3)+dble(NDISD(ix,jp,3)/DNT)
c --- MOLECULE
            SUMOL(1)=SUMOL(1)+dble(MDIS(ix,jp,1)/DNT)
            SUMOL(2)=SUMOL(2)+dble(MDIS(ix,jp,2)/DNT)
            SUMOL(3)=SUMOL(3)+dble(MDIS(ix,jp,3)/DNT)
c --- K AFTER
            SUMAU(1)=SUMAU(1)+dble(NDAU(ix,jp,1)/DNT)
            SUMAU(2)=SUMAU(2)+dble(NDAU(ix,jp,2)/DNT)
            SUMAU(3)=SUMAU(3)+dble(NDAU(ix,jp,3)/DNT)
c --- Rb AFTER
            SUMAD(1)=SUMAD(1)+dble(NDAD(ix,jp,1)/DNT)
            SUMAD(2)=SUMAD(2)+dble(NDAD(ix,jp,2)/DNT)
            SUMAD(3)=SUMAD(3)+dble(NDAD(ix,jp,3)/DNT)

         enddo
         write(4004,*)PPI,SUMNU(1)*AAA
         write(4005,*)PPI,SUMNU(2)*AAA
         write(4006,*)PPI,SUMNU(3)*AAA
         write(4104,*)PPI,SUMAU(1)*AAA
         write(4105,*)PPI,SUMAU(2)*AAA
         write(4106,*)PPI,SUMAU(3)*AAA
         write(5004,*)PPI,SUMND(1)*BBB
         write(5005,*)PPI,SUMND(2)*BBB
         write(5006,*)PPI,SUMND(3)*BBB
         write(5104,*)PPI,SUMAD(1)*BBB
         write(5105,*)PPI,SUMAD(2)*BBB
         write(5106,*)PPI,SUMAD(3)*BBB
         write(6004,*)PPI,SUMOL(1)
         write(6005,*)PPI,SUMOL(2)
         write(6006,*)PPI,SUMOL(3) 
c
      enddo

c --- RPSPACE
         DELAX=(RRS(2,1)-RRS(1,1))/2d0
         DELAY=(RRS(2,2)-RRS(1,2))/2d0
      do ira=1,100
         RIX=RRS(ira,1)-DELAX
         RSUM1=0d0
         RSUM2=0d0
         RSUM3=0d0
         RSUM4=0d0
         do ipa=1,100
            PIX=RRS(ipa,2)-DELAY
            write(3001,*)RIX,PIX,RPK(ira,ipa)
            write(3002,*)RIX,PIX,RPRB(ira,ipa)
            write(3003,*)RIX,PIX,ARPK(ira,ipa)
            write(3004,*)RIX,PIX,ARPRB(ira,ipa)
            RSUM1=RSUM1+RPK(ira,ipa)
            RSUM2=RSUM2+RPRB(ira,ipa)
            RSUM3=RSUM3+ARPK(ira,ipa)
            RSUM4=RSUM4+ARPRB(ira,ipa)
         enddo
         write(3001,*)
         write(3002,*)
         write(3003,*)
         write(3004,*)
         write(3005,*)RIX,RSUM1
         write(3006,*)RIX,RSUM2
         write(3007,*)RIX,RSUM3
         write(3008,*)RIX,RSUM4
      enddo
c --- RPSPACE end
         
c      
      endif 
c
      write(*,*)'check'
      call ENDOMAIN(NK,XMK,OMEGAK,TFERMI,TEMP,ENDIV,ENDU,CPF,IFA,
     & 1,NTMAX)
      call ENDOMAIN(NRB,XMRB,OMEGARB,TCRIT,TEMP,ENDIV,ENDD,CPB,IFB,
     & 2,NTMAX)
c  --- phase space peak density
      PSPDK=maxval(NDISU)/dble(NTMAX)*AAA
      PSPDRB=maxval(NDISD)/dble(NTMAX)*BBB
      PSPDBE=maxval(LBE)/dble(NTMAX)*BBB
      PSPDTH=maxval(LTH)/dble(NTMAX)*BBB
c ---
c ---  volume of particles
      RPMEAN(1)=RMEAN1/dble(NTMAX)/abo
      RPMEAN(3)=RMEAN2/dble(NTMAX)/abo
      RPMEAN(2)=PMEAN1/dble(NTMAX)/hbar
      RPMEAN(4)=PMEAN2/dble(NTMAX)/hbar
c ---
      ABE=dble(JBE/NTMAX)*AAA
      ATH=dble(JTH/NTMAX)*AAA
      CMA(2)=ABE
      CMA(3)=ATH
c
c=============================================
c --- deallocate       
      deallocate(PARU)
      deallocate(PARD)
      deallocate(PARZU,PAREU)
      deallocate(PARZD)
      deallocate(PARED)
      deallocate(DIVE1,DIVE2,DIVE3,DIVE4,IDNUM1,IDNUM2,IDNUM3,IDNUM4)
c ---
      return
      end
c===============================================

c==============================================
      subroutine CXP(XK,PK,XXI,PPI,N)
      implicit real*8(a-h,o-z)
      dimension XXI(N),PPI(N)
!     N must EVEN INTEGER
      if (MOD(N,2) .NE. 0)N=N+1
      DIV=dble(N)/2d0
      XA=XK/DIV
      PA=PK/DIV
      do i=1,N
      XXI(i)=XA*dble(i)-XK     
      PPI(i)=PA*dble(i)-PK 
      enddo

      return
      end
c==============================================
      subroutine CEN(ENMAX,ENDIV,N)
      implicit real*8(a-h,o-z)
      dimension ENDIV(N)
!     N must EVEN INTEGER
      if (MOD(N,2) .NE. 0)N=N+1
      DIV=dble(N)/2d0
      ENS=ENMAX/DIV
      do i=1,N
c      XXI(i)=XA*dble(i)-XK     
c      PPI(i)=PA*dble(i)-PK 
      ENDIV(i)=ENS*dble(i)              
c      XXI(i)=XK/25d0*dble(i)-2d0*XK     
c      PPI(i)=PK/25d0*dble(i)-2d0*PK
      enddo

      return
      end
c=============================================
      subroutine CRP2(RM2,PM2,RRS,N)
c --- for atomic distribution
      implicit real*8(a-h,o-z)
      dimension RRS(N,2)
c --- rp      
!      RSTEP=RM2/dble(N)
!      PSTEP=PM2/dble(N)
!      do i=1,N
!         RRS(i,1)=RSTEP*dble(i)
!         RRS(i,2)=PSTEP*dble(i)
!      enddo
c --- x-y
!     N must EVEN INTEGER
      if (MOD(N,2) .NE. 0)N=N+1
      DIV=dble(N)/2d0
      XA=RM2/DIV
      do i=1,N
      RRS(i,1)=XA*dble(i)-RM2
      RRS(i,2)=XA*dble(i)-RM2
      enddo

      return
      end
c==============================================
c     function G(E):IP=-1 F(E):IP=1 BEC:IP=0
c==============================================
      subroutine ENXP(X,P,XM,OMEGA,TEMP,CP,SUM,VAL,IP)
      use constant
      implicit real*8 (a-h,o-z)
      external GEF
      
      if (IP .EQ. 1)then
         EN=P**2/(2d0*XM)+XM*OMEGA**2*X**2/(2d0)
         EA=GEF(0d0,CP,TEMP,OMEGA,IP)
c         ENXP=GEF(EN,CP,TEMP,OMEGA,IP)*TOP/EA
         VAL=GEF(EN,CP,TEMP,OMEGA,IP)/EA
      else
         if (IP .EQ. -1)then
            EN=P**2/(2d0*XM)+XM*OMEGA**2*X**2/(2d0)
            BH=1d0/(DEXP((hbar*OMEGA-CP)/(BoltzK*TEMP))-1d0)
c            ENXP=GEF(EN,CP,TEMP,OMEGA,IP)*TOP/BH
            VAL=GEF(EN,CP,TEMP,OMEGA,IP)/BH
         else
            EP=P**2
            ER=X**2
            VAL=DEXP(-EP/XM/OMEGA/hbar)*DEXP(-XM*OMEGA*ER/hbar)
         endif
      endif
      SUM=SUM+VAL
      
      return
      end
c==============================================
c     function G(E):IP=-1 or F(E):IP=1 
c=========================================
      function GEF(EN,CP,TEMP,OMEGA,IP)
      use constant
      implicit real*8(A-H,O-Z)
!      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)
c
      alpha=BoltzK*TEMP
      EE=DEXP(EN/alpha)
      EM=DEXP(-CP/alpha)
      if(IP .EQ. -1)GEF=1d0/(EE*EM-1d0)
      if(IP .EQ. 1)GEF=1d0/(EE*EM+1d0)
c      CONST=EN**2/(2d0*(hbar*OMEGA)**3)
c
c      GE=GE*CONST
c     
      return
      end
c==============================================
c     converge?
c===================================================
      subroutine ENCONV(TEMP,CP,TP,OMEGA)
      use constant
      implicit real*8(a-h,o-z)
!      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)
     
      EN=80d0*BoltzK*TP
      VAL=GEF(EN,CP,TEMP,OMEGA)
      write(*,*)'CONVERGENCE',VAL
      
      return
      end
c=====================================================
c     optimize RMAX,PMAX
c=====================================================
      subroutine RPOPT(NZ,TEMP,TSP,OMEGA,CPP,IFB,EOPT)
      use constant
      implicit real*8(a-h,o-z)
c
!      parameter(pi=3.141592654d0)
!      parameter(hbar=1.05457148d-34,BoltzK=1.380d-23)
!      parameter(h=2d0*pi*hbar)
      External GE,DE        !integrand
c    
         if(IFB .EQ. 1)then
            write(*,*)'Fermion Temp',TEMP/TSP
c         CPMAX=1d0*TSP*BoltzK
         else
            write(*,*)'Boson Temp',TEMP/TSP
c         CPMAX=-10d-45
         endif
c
c         CPMIN=-15d0*TSP*BoltzK
c         CPSTEP=(CPMAX-CPMIN)/(CPPOINT-1)
c     
         CP=CPP
         emin=0d0
c         EN=emin
c         ISTEP=10
c         ESTEP=0.1d0*BoltzK*TSP
         EN=0.02d0*BoltzK*TSP
         ISTEP=2
         ESTEP=0.02d0*BoltzK*TSP
         
c         emax=100d0
c         EN=ENMIN
c         ESTEP=(EMAX-EMIN)/(EPOINT-1)
         VAL=dble(NZ)
         NTIMES=0
c         INTE=3000
c         INTE=0
         INTE=2
         do while(VAL .GE. 0d0)
            EN=EN+ESTEP
            INTE=INTE+ISTEP

c            INTE=4000
c            INTE=1200
!            INTE=2000
!            INTE=1000
c     
c            VAL1=VAL2
           call trap(emin,EN,INTE,VAL,NZ,CP,TEMP,OMEGA,IFB)
c           EN=EN+ESTEP
           NTIMES=NTIMES+1
c           write(*,*)EN/(BoltzK*TSP),VAL
           enddo
           write(*,*)'EN in optimization',EN/(BoltzK*TSP),NTIMES,VAL
           EOPT=EN
          call trap(emin,EOPT,INTE,AVAL,NZ,CPP,TEMP,OMEGA,IFB)
          call FU(CPP,TEMP,BVAL)
          write(*,*)NZ-AVAL-BVAL,'excite in opt'
          call trap(emin,EOPT/2d0,INTE,AVAL,NZ,CPP,TEMP,OMEGA,IFB)
          call FU(CPP,TEMP,BVAL)
          write(*,*)NZ-AVAL-BVAL,'excite in opt/2'
          EOPT=EOPT/(BoltzK*TSP)

      return
      end
c
c=====================================================-
      function FNM(CPF,TEMP)
      use constant
      implicit real*8(a-h,o-z)

      BETA=BOLTZK*TEMP
      d1=dexp(-CPF/BETA)
      FNM=1d0/(d1+1d0)
      
      return
      end
c====================================================
c     check distribution at ENERGY axis
c===================================================
      subroutine chkENdis(NK,PARU,WKR,WKZ,XMK,TFERMI,EDIS,ENDIV,dely,IP)
      use constant
      implicit real*8(a-h,o-z)
      dimension PARU(NK,7),EDIS(100)
      dimension ENDIV(100)
c
      ds=0d0
      if (IP .EQ. 1)ds=dely

      CR1R=XMK*WKR**2/2d0
      CR1Z=XMK*WKZ**2/2d0
      CP1=1d0/(2d0*XMK)
      do i=1,NK
         R=PARU(i,1)**2+(PARU(i,2)-ds)**2
         P=PARU(i,4)**2+PARU(i,5)**2+PARU(i,6)**2 
         EN=R*CR1R+PARU(i,3)**2*CR1Z+P*CP1
         do j=1,100
            if (EN .LE. ENDIV(j))then
               EDIS(j)=EDIS(j)+1
               exit
            endif
         enddo 
      enddo


      return
      end

c=====================================================
      subroutine ENDOMAIN(N,XMK,OMEGA,TSP,TEMP,EN,ENDIV,CPS,IP,IS,NTMAX)
      use constant 
      implicit real*8(a-h,o-z)
      dimension EN(100),ENDIV(100)
 
      CR=XMK*OMEGA**2/2d0
      CP=1d0/(2d0*XMK)
      AA=100d0/dble(N)
      INTE=10000
      const=BOLTZK*TSP
      do i=1,100
         if (i .eq. 1)E1=0d0
         E2=EN(i)
         call trap2(E1,E2,INTE,VAL,N,CPS,TEMP,OMEGA,IP)
         E1=E2
         CE2=E2/const
         if(IS .EQ. 1)then
         write(7001,*)CE2,VAL
         write(7002,*)CE2,ENDIV(i)/dble(NTMAX)
         else
         write(7003,*)CE2,VAL
         write(7004,*)CE2,ENDIV(i)/dble(NTMAX)
         endif
      enddo
c
      return
      end

c============================================================
c     subroutine trapezium integration
c=============================================================
      subroutine trap2(EN1,ENN,INTE,VAL,N,CP,TEMP,OMEGA,IFB)
      use constant
      implicit real*8(A-H,O-Z)
      External DE,GE
      dimension XZ(99999),WZ(99999),WK(99999)
      data itimes/0/
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
c         VAL=N-VAL
      else
         do nn=1,INTE
            VAL=VAL+GE(XZ(nn),CP,TEMP,OMEGA)*WZ(nn)
         enddo
      endif
c
      return
      end
c========================================      
c========================================
      subroutine varia(PAR,N,AVG,VAR)
c     obtain avarage of (X,P) and variance
c=======================================
      implicit real*8(a-h,o-z)
      dimension PAR(N,7),AVG(6),VAR(6)

      AVG=0d0
      VAR=0d0
      do i=1,N
         do j=1,6
            AVG(j)=AVG(j)+DABS(PAR(i,j))/dble(N)
         enddo
      enddo

      do i=1,N
         do j=1,6
            VAR(j)=(AVG(j)-DABS(PAR(i,j)))**2/dble(N)
         enddo
      enddo

      return
      end
c========================================      
! old 
!      subroutine divide(IFA,CPF,TEMP,OMEGA,DIVE,IDNUM,ENMAX,N)
!      use constant
!      implicit real*8(a-h,o-z)
!c      parameter (INTE=10000,ENOV=2d-28)
!      dimension DIVE(2,N),IDNUM(N)
!
!      INTE=10000
!      ENOV=2d-28
!      EN1=ENMAX/dble(N)
!      EMIN=0d0
!      EMAX=EN1
!      BETA=BOLTZK*TEMP
!      VAL=0d0
!      M=0
!c      write(*,*)ENMAX,ENOV,OMEGA,CPF,IFA,TEMP,VAL,INTE
!
!c      write(*,*)'2?'
!c      call trap(ENMAX,ENOV,INTE,VAL,M,CPF,TEMP,OMEGA,IFA)
!      call trap2(ENMAX,ENOV,INTE,VAL,M,CPF,TEMP,OMEGA,IFA)
!      VAL=DABS(VAL)
!      AM=VAL/dble(N)
!      if (AM .lt. 1d0)AM=1d0
!      write(*,*)'AM',AM
!
!c      write(*,*)'1?'
!      do i=1,N
!c        call trap(emin,emax,INTE,VAL,M,CPF,TEMP,OMEGA,IFA)
!        call trap2(EMIN,EMAX,INTE,VAL,M,CPF,TEMP,OMEGA,IFA)
!        VAL=DABS(VAL)+AM
!        IDNUM(i)=DINT(VAL)
!        DIVE(1,i)=EMAX
!        if (IFA .eq. 1)DIVE(2,i)=1d0/(dexp((EMIN-CPF)/BETA)+1d0)
!        if (IFA .eq. -1)DIVE(2,i)=1d0/(dexp((EMIN-CPF)/BETA)-1d0)
!        EMIN=EMIN+EN1
!        EMAX=EMAX+EN1
!      enddo
!        
!c ---  output
!c      if (IFA .eq. -1)then
!c      do i=1,N
!c         write(2011,*)DIVE(1,i),IDNUM(i)
!c         write(2012,*)DIVE(1,i),DIVE(2,i)         
!c      enddo
!c      endif
!
!      return
!      end
c===========================================
      function fen(XMA,WR,WZ,X,Y,Z,PX,PY,PZ)
      use constant
      implicit real*8(a-h,o-z)

      E=(PX**2+PY**2+PZ**2)/(2d0*XMA)
      EZ=XMA*WZ**2*Z**2/2d0
      FEN=E+EZ+XMA*WR**2*(X**2+Y**2)/2d0
      
      return
      end
c============================================
      subroutine MAXSET(SETMAX,XMA,WR,WZ,TEMP,TSP,IP,CPS,
     & DIVE,IDNUM,IDN,EDIM,ENMAX)
      use constant
      implicit real*8(a-h,o-z)
      dimension DIVE(2,IDN),IDNUM(IDN),SETMAX(6),EDIM(1000)

      if (IP .NE. 0)then
c --- fix check
            BETA=BOLTZK*TEMP
         if (IP .EQ. 1)then
c            BETA=BOLTZK*TEMP
            ECH=ENMAX
            ECH=ECH/(TSP*BOLTZK)
            if (ECH .le. 5)then
               ENMAX=TSP*BOLTZK*5d0
               write(*,*)'ENMAX CHANGED'
            endif
         endif
         write(*,*)IP,ENMAX/BETA
c ---
         call SETCM(XMAX,PXMAX,YMAX,PYMAX,ZMAX,PZMAX,
     &        XMA,WR,WZ,TEMP,TSP,BETA,FE,IP,CPS,ENMAX)
         OMEGA=(WR**2*WZ)**(1d0/3d0)
         call divide(IP,CPS,TEMP,OMEGA,DIVE,IDNUM,ENMAX,IDN,EDIM)
         SETMAX(1)=XMAX
         SETMAX(2)=PXMAX
         SETMAX(3)=YMAX
         SETMAX(4)=PYMAX
         SETMAX(5)=ZMAX
         SETMAX(6)=PZMAX        
      else
         call SETCZ(XMAX,PXMAX,YMAX,PYMAX,ZMAX,PZMAX,
     &        XMA,WR,WZ,TEMP,TSP,BETA,FE)
         SETMAX(1)=XMAX
         SETMAX(2)=PXMAX
         SETMAX(3)=YMAX
         SETMAX(4)=PYMAX
         SETMAX(5)=ZMAX
         SETMAX(6)=PZMAX     
      endif       

      return
      end
c =======================================================
      subroutine setdiv(IP,ENEND,CPS,TEMP,OMEGA,IDN,EDIM,IDE)
      use constant
      implicit real*8(a-h,o-z)
      dimension EDIM(1000)
      external GEF
      
      BETA=BOLTZK*TEMP
      INTE=1000
      M=0
      EN=0d0
      ENA=0d0
      ENB=EN
      ENMAX=ENEND*BETA+CPS
      ESTEP=OMEGA*hbar
c      ESTEP=ENMAX/100d0
      BA=0d0
      BCHK=0d0
      VAL=1d0
      do while(BCHK .ge. 0)
         ENA=ENA+ESTEP
         call trap2(ENB,ENA,INTE,VAL,M,CPS,TEMP,OMEGA,IP)
c         write(*,*)VAL
         BCHK=VAL-BA
         BA=VAL
         ENB=ENA
      enddo
c      ENMAX=ENEND*BETA+CPS
      write(*,*)'ENA',ENA,ENA/ENMAX
c      ECHK1=ENA
      
      do while(VAL .ge. 1)
         ENA=ENA+ESTEP
         call trap2(ENB,ENA,INTE,VAL,M,CPS,TEMP,OMEGA,IP)
c         write(*,*)VAL
         ENB=ENA
      enddo      
      write(*,*)'ENA2',ENA,ENA/ENMAX
      ECHK1=ENA

c      IEC1=300
      IEC1=IDE
c      ECHK1=ENA/2d0
c      ECHK2=ENA*4d0
      ESTEP=ECHK1/dble(IEC1)
c      ESTEP2=(ENMAX-ECHK1)/200d0
c      ESTEP3=(ENMAX-ECHK2)/20d0

      EDIM(1)=EN
c      BETA=BOLTZK*TEMP
c      FMIN=GEF(EN,CPS,TEMP,OMEGA,IP)
      IDN=0
c      ESTEP=hbar*OMEGA
c      ENMAX=ENEND*BETA+CPS
c      ESTEP=ENMAX/100d0
c      ESTEP2=ENMAX/300d0
c      ESTEP3=ENMAX/50d0
c      ESTEP3=ENMAX
c      ECHK1=ENMAX/4d0
c      ECHK2=ENMAX/3d0*2d0
c      ENMAX=ENEND*BETA+CPS
c      write(*,*)ENMAX,ESTEP,ECH,FMIN
c      VALSUM=0d0
      do i=1,IEC1
         EN=EN+ESTEP
         IDN=IDN+1
         IA=IDN+1
         EDIM(IA)=EN
      enddo
      write(*,*)'ECHECK',EDIM(IA),ECHK1,EN
c      do while(EN .lt. ECHK2)
c         EN=EN+ESTEP2
c         IDN=IDN+1
c         IA=IDN+1
c         EDIM(IA)=EN
c      enddo
      
      call trap2(ENA,ENMAX,INTE,VAL,M,CPS,TEMP,OMEGA,IP)
      VAL=VAL+1
      IEC2=DINT(VAL)
      ESTEP3=(ENMAX-ENA)/dble(IEC2)

      do i=1,IEC2
         EN=EN+ESTEP3
         IDN=IDN+1
         IA=IDN+1
         EDIM(IA)=EN
      enddo
      write(*,*)'ECHECK',EDIM(IA),ENMAX,EN      
      
c      IDA=IDN+1
c      EDIM(IDA)=EN
      
c      IDN=IDN+1
c      IX=IDN+1
c      EDIM(IX)=ENMAX
      ENEND=ENMAX
      
c      write(*,*)IP,EDIM
      return
      end
c==============================================================
      subroutine divide(IFA,CPF,TEMP,OMEGA,DIVE,IDNUM,ENMAX,N,EDIM)
      use constant
      implicit real*8(a-h,o-z)
      dimension DIVE(2,N),IDNUM(N),EDIM(1000)

      INTE=1000
      ENOV=100d-28
      EN1=ENMAX/dble(N)
      EMIN=EDIM(1)
c      EMAX=EN1
      BETA=BOLTZK*TEMP
      VAL=0d0
      M=0
c
      call trap2(ENMAX,ENOV,INTE,VAL,M,CPF,TEMP,OMEGA,IFA)
      VAL=DABS(VAL)
      AM=VAL/dble(N)
      if (AM .lt. 1d0)AM=1d0
c      
      write(*,*)'AM',AM
      do i=1,N
         il=i+1
         EMAX=EDIM(il)
c        call trap(emin,emax,INTE,VAL,M,CPF,TEMP,OMEGA,IFA)
        call trap2(EMIN,EMAX,INTE,VAL,M,CPF,TEMP,OMEGA,IFA)
        VAL=DABS(VAL)+AM
        IDNUM(i)=DINT(VAL)
        DIVE(1,i)=EMAX
        if (IFA .eq. 1)DIVE(2,i)=1d0/(dexp((EMIN-CPF)/BETA)+1d0)
        if (IFA .eq. -1)DIVE(2,i)=1d0/(dexp((EMIN-CPF)/BETA)-1d0)
c        write(*,*)IFA,DIVE(1,i)/ENMAX,DIVE(2,i),IDNUM(i)
        EMIN=EMAX
      enddo
      write(*,*)i,SUM(IDNUM)
        
c ---  output
c      if (IFA .eq. -1)then
c      do i=1,N
c         write(2011,*)DIVE(1,i),IDNUM(i)
c         write(2012,*)DIVE(1,i),DIVE(2,i)         
c      enddo
c      endif

      return
      end      
