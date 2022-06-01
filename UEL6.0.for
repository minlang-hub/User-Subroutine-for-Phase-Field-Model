! ==========================// UEL //=================================== 
! User Subroutine UEL for Abaqus:
! Element type U1: displacement and phase-field 2D element.
! DOF 1-2 : Displacement
! DOF 11  : Phase field
!     
! Unified phase field model of Jian-Ying Wu(2017):
!     Geometry function: a(d) = 2d - d*d
!     Cornelissen softening cohesive model is adopted:
!     p=2, a2=1.3868, a3=0.6567
!
! Material properties to be given through the input file (*.inp), are:
! PROPS(1) = Length scale parameter (lc)
! PROPS(2) = Critical energy release rate (G)
! PROPS(3) = Young's modulus (E)
! PROPS(4) = Poisson's ratio (nu)
! PROPS(5) = thickness of the plate (t)
!
! ----SVARS updated every increment
! SVARS for each integral point:
!     1:  ENG, the elastic energy of increment n
!     2:  HIST, the history variable of increment n
!________________________________________________________________________    
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
      INCLUDE 'ABA_PARAM.INC'
!________________________________________________________________________
! ************************************************************
! *  DATA DICTIONARY: Define variable types, definitions.    *
! ************************************************************
      REAL(8)::ZERO=0.D0,ONE=1.D0,TWO=2.D0,HALF=0.5D0
      INTEGER NSVARS,NRHS,NDOFEL,MDLOAD,MLVARX,NPREDF,MCRD,NPROPS,NNODE,
     1 JTYPE,KINC,NDLOAD,NJPROP
      REAL(8) SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),
     1 COORDS(MCRD,NNODE),U(NDOFEL),DU(MLVARX,1),V(NDOFEL),A(NDOFEL),
     2 TIME(2),PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     3 DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*),
     4 AMATRX(NDOFEL,NDOFEL),RHS(MLVARX,1),DTIME,PNEWDT,PERIOD
      INTEGER JELEM
      INTEGER i,j,l,k,ipt                 ! DO loop variables
      REAL(8) xi, eta                     ! Local coordinates
      REAL(8) GAUSS(NNODE,2),WEIGHT(NNODE)! Local coordinates and weights of ipts
      REAL(8) det_J                       ! Jacobian determinant
      REAL(8) Np(1,4),Np_T(4,1)           ! [N1 N2 N3 N4]
      REAL(8) Bp(2,4),Bp_T(4,2)           ! B matrix in phase field
      REAL(8) Bu(3,2*NNODE),Bu_T(8,3)     ! Strain-displacement matrix
      REAL(8) dP(2,1)                     ! Gradient of phase field
      REAL(8) CMAT(3,3)                   ! Elastic matrix
      REAL(8) UU(8)                       ! Dof of displacement
      REAL(8) PHI(NNODE,1)                ! Dof of phase-field
      REAL(8) STRAIN(3,1),STRESS(3,1)     ! Array of strain & stress
      REAL(8) PHASE                       ! Phase field value of current ipt
      REAL(8) Ru(2*NNODE)                 ! Residual of displacment
      REAL(8) Rp(NNODE)                   ! Residual of phase field
      REAL(8) Kuu(8,8)                    ! Element stiffness matrix of displacement
      REAL(8) Kpp(4,4)                    ! Element stiffness matrix of phase-field
      REAL(8) ENG,HIST
      REAL(8) lb, Gf, EMOD, ENU, THCK
      REAL(8) Kuu11(8,3),Kuu12(8,8),Kuu1(8,8)
      REAL(8) Kpp1(4,4),Kpp11(4,4),Kpp21(4,2),Kpp22(4,4)
      REAL(8) Rp11(4,1)
      REAL(8) omega, domega, ddomega, dalpha, ddalpha, c0
      REAL(8) omega0,domega0,ddomega0,dalpha0,ddalpha0,PHASE0
      
      ! Initialization
      dNdx=0.d0; det_J=0.d0; UU=0.d0; PHI=0.d0; ENG=0.d0
      Np=0.d0; Np_T=0.d0; Bp=0.d0; Bp_T=0.d0; Bu=0.d0; Bu_T=0.d0
      dP=0.d0; CMAT=0.d0; CMAT1=0.d0; STRESS=0.d0; STRAIN=0.d0;
      PHASE=0.d0; Ru=0.d0; Rp=0.d0;AMATRX=0.d0; Kuu=0.d0; Kpp=0.d0;
      RHS=0.d0; STRAIN_=0.d0; STRESS_=0.d0; A2=0.d0;Kuu11=0.d0; 
      Kuu12=0.d0; Kuu1=0.d0; omega=0.d0; domega=0.d0; ddomega=0.d0
      dalpha=0.d0; ddalpha=0.d0; omega0=0.d0; domega0=0.d0;
      ddomega0=0.d0; dalpha0=0.d0; ddalpha0=0.d0; PHASE0=0.d0
      Kpp1=0.d0;Kpp11=0.d0; Kpp21=0.d0; Kpp22=0.d0; Rp11=0.d0
      T1=0.d0; T2=0.d0; T2_T=0.d0; HIST=0.d0; theta=0.d0
! ______________________________________________________________________
! **************************************************************
! *  Constructing elemet TYPE U1( Displacement & Phase-field ) *
! **************************************************************
      ! < 1. Material parameters >
      lb =PROPS(1)
      Gf =PROPS(2)
      EMOD  =PROPS(3)
      ENU   =PROPS(4)
      THCK  =PROPS(5)

      ! < 2. Materials stiffness matrix (plane stress) >                
      CMAT(1,1)=EMOD/(1-ENU**TWO)
      CMAT(2,2)=EMOD/(1-ENU**TWO)
      CMAT(3,3)=EMOD/(TWO*(ONE+ENU))
      CMAT(1,2)=EMOD/(1-ENU**TWO)*ENU
      CMAT(2,1)=EMOD/(1-ENU**TWO)*ENU
      
      ! < 3. Initialization >
      DO i =1,NNODE
          ! DOF of displacement
          UU(2*i-1)=U(3*i-2); UU(2*i)=U(3*i-1)
          ! DOF of phase-field
          PHI(i,1)=U(3*i)
      END DO
      GAUSS(1,1)=-0.577350269189626; GAUSS(1,2)=-0.577350269189626
      GAUSS(2,1)= 0.577350269189626; GAUSS(2,2)=-0.577350269189626
      GAUSS(3,1)= 0.577350269189626; GAUSS(3,2)= 0.577350269189626
      GAUSS(4,1)=-0.577350269189626; GAUSS(4,2)= 0.577350269189626
      WEIGHT=ONE
      c0=3.1415926535897932384626433832795028841971693993751d0

      ! < 4. Calculating properties at each integral point¡ý >
      LOOP_IPT :DO ipt=1,4
      ! ---------------------------< LOOP:Begin >-------------------------------
          ! < Shape functions and its local derivatives >
          xi  = GAUSS(ipt,1); eta = GAUSS(ipt,2)  ! Local coordinates of current ipt
          CALL SHAPEFUN(Np,Bu,Bp,det_J,xi,eta,COORDS)
          
          ! < phase-field value and its gradient at ipt >
          PHASE=0.d0; dP=0.d0
          DO i=1,NNODE
              PHASE=PHASE+Np(1,i)*PHI(i,1)
              DO j=1,2
                  dP(j,1)=dP(j,1)+Bp(j,i)*PHI(i,1)
              END DO
          END DO
          
          ! < Strains & Stresses of ipt(GCS) >
          DO i=1,3
              STRAIN(i,1)=0.d0
              DO j=1,2*NNODE
                  STRAIN(i,1)=STRAIN(i,1)+Bu(i,j)*UU(j)
              END DO
          END DO
          CALL UMATMUL(CMAT,STRAIN,STRESS,3,1,3)
          
          ! <  Hitory variables >
          CALL CAL_ENERGY_PLUS(ENG,STRAIN,EMOD,ENU)
          SVARS(2*(ipt-1)+1)=ENG
          IF (SVARS(2*(ipt-1)+1).GT.SVARS(2*(ipt-1)+2)) THEN
              HIST=SVARS(2*(ipt-1)+1)
          ELSE
              HIST=SVARS(2*(ipt-1)+2)
          END IF
          dalpha0  =  2.D0
          ddalpha0 =  -2.D0
          PHASE0   = ZERO
          CALL DEGRADFUNC(omega0,domega0,ddomega0,PHASE0,PROPS,c0)
          IF (HIST<-2.d0*Gf/c0/lb/dOMEGA0) THEN
              HIST=-2.d0*Gf/c0/lb/dOMEGA0     ! Critical value
          END IF
          SVARS(2*(ipt-1)+2)=HIST
          
          ! < Derivatives of geometry function & degradation function  >
          dalpha  =   2.D0-2.D0*PHASE
          ddalpha =  -2.D0
          CALL DEGRADFUNC(omega,domega,ddomega,PHASE,PROPS,c0)
          
          ! < Element stiffness matrix >
          ! --------------------- Kuu ----------------------
          CALL UTRANSPOSE(Bu,Bu_T,3,8)
          CALL UMATMUL(Bu_T,CMAT,Kuu11,8,3,3)
          CALL UMATMUL(Kuu11,Bu,Kuu12,8,8,3)
          Kuu1 = omega*Kuu12*det_J*WEIGHT(ipt)*THCK
          Kuu = Kuu+Kuu1
          ! --------------------- Kpp -----------------------
          CALL UTRANSPOSE(Np,Np_T,1,4)
          CALL UTRANSPOSE(Bp,Bp_T,2,4)
          CALL UMATMUL(Np_T,Np,Kpp11,4,4,1)
          CALL UMATMUL(Bp_T, Bp,Kpp22,4,4,2)
          Kpp1 = ( (Gf/c0*ddalpha/lb+ddomega*HIST)*Kpp11 + 
     1        Gf/c0*2*lb*Kpp22 )*det_J*WEIGHT(ipt)*THCK
          Kpp = Kpp+Kpp1
          
          ! < Internal forces >
          ! --------------------- Ru -----------------------
          DO i=1,2*NNODE
              DO j=1,2*NNODE
                  Ru(i)=Ru(i)-Kuu1(i,j)*UU(j)
              END DO
          END DO
          ! --------------------- Rp -----------------------
          CALL UMATMUL(Bp_T,dP,Rp11,4,1,2)
          DO i=1,NNODE
              Rp(i)=Rp(i)-( Gf/c0*2*lb*Rp11(i,1)+
     1            Np_T(i,1)*(Gf/c0*dalpha/lb+domega*HIST) )*
     1                                WEIGHT(ipt)*det_J*THCK
          END DO
 
      ! ---------------------------< LOOP:End >--------------------------------- 
      END DO LOOP_IPT
      
      ! < 5. Updatding AMATRX & RHS >
      DO i=1,4
          DO j=1,4
              AMATRX(3*i-2,3*j-2) =Kuu(2*i-1,2*j-1)
              AMATRX(3*i-1,3*j-2) =Kuu(2*i,2*j-1)
              AMATRX(3*i-2,3*j-1) =Kuu(2*i-1,2*j)
              AMATRX(3*i-1,3*j-1) =Kuu(2*i,2*j)
              AMATRX(3*i,3*j)     =Kpp(i,j)
          END DO
      END DO
      DO i=1,NNODE
          RHS(3*i-2,1)=Ru(2*i-1)
          RHS(3*i-1,1)=Ru(2*i)
          RHS(3*i,1)  =Rp(i)
      END DO
      
      END SUBROUTINE UEL
! =============================// END of UEL//======================================= 
      
! < 1 >
!______________________________________________________________________
      SUBROUTINE SHAPEFUN(N,Bu,Bp,det_J,xi,eta,COORDS)                !|
!----------------------------------------------------------------------|
!     SUBROUTINE SHAPEFUN                                              |
!         Calculate the shape function and its derivatives             |
!     IN : N, Bu, Bp, det_J, xi, eta                                   |
!     OUT: N, Bu, Bp, det_J                                            |
!----------------------------------------------------------------------|
          IMPLICIT NONE                                               !|
          INTEGER i,j                                                 !|
          REAL(8) N(1,4),dNdx(2,4),dNdxi(2,4),COORDS(2,4)             !|
          REAL(8) JACOBI(2,2),INV_JACOBI(2,2), det_J                  !|
          REAL(8) Bu(3,8),Bp(2,4)                                     !|
          REAL(8) xi,eta                                              !|
          REAL,PARAMETER::ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0     !|
                                                                      !|
          ! Initialization                                             |
          N=0.d0; dNdx=0.d0; dNdxi=0.d0; JACOBI=0.d0                  !|
          INV_JACOBI=0.d0; det_J=0.d0; Bu=0.d0; Bp=0.d0               !|
          ! Values of shape functions as a function of local coord.   !|
          N(1,1) = ONE/FOUR*(ONE-xi)*(ONE-eta)                        !|
          N(1,2) = ONE/FOUR*(ONE+xi)*(ONE-eta)                        !|
          N(1,3) = ONE/FOUR*(ONE+xi)*(ONE+eta)                        !|
          N(1,4) = ONE/FOUR*(ONE-xi)*(ONE+eta)                        !|
                                                                      !|
          ! Derivatives of shape functions respect to local coordinates|
          dNdxi(1,1) = ONE/FOUR*(MONE+eta)                            !|
          dNdxi(1,2) = ONE/FOUR*(ONE-eta)                             !|
          dNdxi(1,3) = ONE/FOUR*(ONE+eta)                             !|
          dNdxi(1,4) = ONE/FOUR*(MONE*ONE-eta)                        !|
          dNdxi(2,1) = ONE/FOUR*(MONE*ONE+xi)                         !|
          dNdxi(2,2) = ONE/FOUR*(MONE*ONE-xi)                         !|
          dNdxi(2,3) = ONE/FOUR*(ONE+xi)                              !|
          dNdxi(2,4) = ONE/FOUR*(ONE-xi)                              !|
                                                                      !|
          ! Jacobi matrix                                             !|
          DO j=1,4                                                    !|
              JACOBI(1,1)=JACOBI(1,1)+dNdxi(1,j)*COORDS(1,j)          !|
              JACOBI(1,2)=JACOBI(1,2)+dNdxi(1,j)*COORDS(2,j)          !|
              JACOBI(2,1)=JACOBI(2,1)+dNdxi(2,j)*COORDS(1,j)          !|
              JACOBI(2,2)=JACOBI(2,2)+dNdxi(2,j)*COORDS(2,j)          !|
          END DO                                                      !|
                                                                      !|
          ! Jacobi determinant                                         |
          det_J=JACOBI(1,1)*JACOBI(2,2)-JACOBI(1,2)*JACOBI(2,1)       !|
                                                                      !|
          ! Inverse of JACOBI                                         !|
          INV_JACOBI(1,1)= JACOBI(2,2)/det_J                          !|
          INV_JACOBI(1,2)=-JACOBI(1,2)/det_J                          !|
          INV_JACOBI(2,1)=-JACOBI(2,1)/det_J                          !|
          INV_JACOBI(2,2)= JACOBI(1,1)/det_J                          !|
                                                                      !|
          ! dNdx                                                       |
          CALL UMATMUL(INV_JACOBI,dNdxi,dNdx,2,4,2)                   !|
                                                                      !|
          ! Bp & Bu matrix( Bp=dNdx )                                 !|
          DO i=1,4                                                    !|
              Bp(1,i)=dNdx(1,i);      Bp(2,i)=dNdx(2,i)               !|
              Bu(1,i*2-1)=dNdx(1,i);  Bu(2,i*2)  =dNdx(2,i)           !|
              Bu(3,i*2-1)=dNdx(2,i);  Bu(3,i*2)  =dNdx(1,i)           !|
          END DO                                                      !|
      END SUBROUTINE SHAPEFUN                                         !|
!______________________________________________________________________|     
          
! < 3 >
!______________________________________________________________________
      SUBROUTINE CAL_ENERGY_PLUS(PSI_PLUS,STRAIN,E,NU)                !|
!---------------------------------------- -----------------------------|
!     SUBROUTINE CAL_ENERGY_PLUS                                       |
!         Calculate  the dagadation energy                             |
!     IN : STRAIN, E, NU                                               |
!     OUT: PSI_PLUS                                                    |
!----------------------------------------------------------------------|
      IMPLICIT NONE                                                   !|
      REAL(8) PSI_PLUS                                                !|
      REAL(8) STRAIN(3,1), epsilon_x, epsilon_y, gamma ,p ,r          !|
      REAL(8) E, NU, LAMBDA, MU                                       !|
      REAL(8) epsilon1, epsilon2, trace                               !|
      REAL(8) epsilon1_plus,epsilon2_plus, trace_plus                 !|
                                                                      !|
      ! Lame constant                                                 !|
      LAMBDA = E*NU/(1-NU**2)             ! Plane stress              !|
      !LAMBDA = E*(1-NU)/((1+NU)*(1-2*NU)) ! Plane strain             !|
      MU=E/(2.d0*(1.d0+NU))                                           !|
                                                                      !|
      epsilon_x=STRAIN(1,1);  epsilon_y=STRAIN(2,1)                   !|
      gamma=STRAIN(3,1)                                               !|
      p=(epsilon_x+epsilon_y)/2.d0                                    !|
      r=sqrt(((epsilon_x-epsilon_y)/2.d0)**2.d0+(gamma/2.d0)**2.d0)   !|
                                                                      !|
      ! Calculating principle stresses                                !|
      epsilon1=p+r;  epsilon2=p-r                                     !|
      trace = epsilon1+epsilon2                                       !|
                                                                      !|
      ! Calculating PSI_PLUS                                          !|
      CALL Macaulay(epsilon1,epsilon1_plus)                           !|
      CALL Macaulay(epsilon2,epsilon2_plus)                           !|
      CALL Macaulay(trace,trace_plus)                                 !|
      PSI_PLUS=0.5*LAMBDA*trace_plus**2.d0                            !|
     1    +MU*(epsilon1_plus**2.d0 + epsilon2_plus**2.d0)             !|
                                                                      !|
      END SUBROUTINE CAL_ENERGY_PLUS                                  !|
!______________________________________________________________________|


! < 5 >
!______________________________________________________________________
      SUBROUTINE DEGRADFUNC(OMEGA,dOMEGA,ddOMEGA,pf,PROPS,c0)         !|
!---------------------------------------- -----------------------------|
!     SUBROUTINE DEGRADFUNC                                            |
!         Calculate  the dagadation function and its derivatives       |
!     IN : pf, PROPS, c0                                               |
!     OUT: ENERGY_PLUS                                                 |
!----------------------------------------------------------------------|
      IMPLICIT NONE                                                   !|
      REAL(8):: OMEGA, dOMEGA, ddOMEGA, pf, PROPS(7), c0              !|
!     local varibles                                                   |
      REAL(8):: g, dg, ddg, fai, dfai, ddfai, E, Gf, ft, lb           !|
      REAL(8):: p, a1, a2, a3                                         !|
                                                                      !|
      E=PROPS(3); Gf=PROPS(2); ft=PROPS(7); lb=PROPS(1)               !|
      a1=4.d0/(c0*lb)*E*Gf/(ft*ft)                                    !|
      a2=1.3868d0;  p=2.d0                                            !|
      a3=0.6567d0;                                                    !|
                                                                      !|
      g=(1.d0-pf)**p; dg=-p*(1.d0-pf)**(p-1.d0);                      !|
      ddg=p*(p-1.d0)*(1.d0-pf)**(p-2.d0)                              !|
                                                                      !|
      fai=g+a1*pf+a1*a2*pf**2.d0+a1*a2*a3*pf**3.d0                    !|
      dfai=dg+a1+2.d0*a1*a2*pf+3.d0*a1*a2*a3*pf**2.d0                 !|
      ddfai=ddg+2.d0*a1*a2+6.d0*a1*a2*a3*pf                           !|
                                                                      !|
      omega=g/fai                                                     !|
                                                                      !|
      domega=(dg*fai-g*dfai)/(fai**2.d0)                              !|
      ddomega=((ddg*fai-g*ddfai)*fai-2.d0*(dg*fai-g*dfai)*dfai)       !|
     &        /(fai**3.d0)                                            !|
                                                                      !|
      END SUBROUTINE degradfunc                                       !|
!______________________________________________________________________|
                                                                      
                                                                      
! < 3 >
!______________________________________________________________________
      SUBROUTINE Macaulay(x,y)                                        !|
!----------------------------------------------------------------------|
!     SUBROUTINE Macaulay                                              |
!         Fcuntion of Macaulay brackets                                |
!     IN : x                                                           |
!     OUT: y                                                           |
!----------------------------------------------------------------------|
      IMPLICIT NONE                                                   !|
      REAL(8) x,y                                                     !|
      y=0.5*(DABS(x)+x)                                               !|
      END SUBROUTINE Macaulay                                         !|
!______________________________________________________________________|

      
! < 4 >
!______________________________________________________________________      
      SUBROUTINE UMATMUL(A,B,C,m,n,l)                                 !|
!----------------------------------------------------------------------|
!     SUBROUTINE Umatmul                                               |
!         Matrix mutiplication                                         |
!     IN : A,B,C,m,n,l                                                 |
!     OUT: C                                                           |
!----------------------------------------------------------------------|
      IMPLICIT NONE                                                   !|
      INTEGER i,j,k,m,n,l                                             !|
      REAL(8) A(m,l),B(l,n),C(m,n)                                    !|
      DO i=1,m                                                        !|
        DO j=1,n                                                      !|
            C(i,j)=0.d0                                               !|
            DO k=1,l                                                  !|
                C(i,j)=C(i,j)+A(i,k)*B(k,j)                           !|
            END DO                                                    !|
        END DO                                                        !|
      END DO                                                          !|
      END SUBROUTINE UMATMUL                                          !|
!______________________________________________________________________|

! < 5 >
!______________________________________________________________________
      SUBROUTINE UTRANSPOSE(A,A_T,m,n)                                !|
!----------------------------------------------------------------------|
!     SUBROUTINE UTRANSPOSE                                            |
!         Matrix transpose                                             |
!     IN : A,A_T,m,n                                                   |
!     OUT: A_T                                                         |
!----------------------------------------------------------------------|
      IMPLICIT NONE                                                   !|
      INTEGER i,j,m,n                                                 !|
      REAL(8) A(m,n),A_T(n,m)                                         !|
      DO i=1,m                                                        !|
          DO j=1,n                                                    !|
              A_T(j,i)=A(i,j)                                         !|
          END DO                                                      !|
      END DO                                                          !|
      END SUBROUTINE UTRANSPOSE                                       !|
!______________________________________________________________________|