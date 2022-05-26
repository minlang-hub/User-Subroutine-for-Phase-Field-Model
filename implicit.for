! ==========================// UEL //=================================== 
! User Subroutine UEL for Abaqus:
! Element type U1: displacement and phase-field element.
! DOF 1-2 : Displacement
! DOF 11  : Phase field    
!   
! Material properties to be given through the input file (*.inp), are:
! PROPS(1) = lb    --- length scale parameter
! PROPS(2) = GcI   --- critical energy release rate
! PROPS(3) = E     --- Young's modulus
! PROPS(4) = nu    --- Possion's ratio
! PROPS(5) = thck  --- thickness
!
! ----SVARS updated every increment
! SVARS for each integral point:
!     1:  hist,   the history variable
! Total number of sate-dependent variables: 4
!________________________________________________________________________    
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
      INCLUDE 'ABA_PARAM.INC'
      INTEGER NSVARS,NRHS,NDOFEL,MDLOAD,MLVARX,NPREDF,MCRD,NPROPS,NNODE,
     1 JTYPE,KINC,NDLOAD,NJPROP
      REAL(8) SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),
     1 COORDS(MCRD,NNODE),U(NDOFEL),DU(MLVARX,1),V(NDOFEL),A(NDOFEL),
     2 TIME(2),PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(16,1),
     3 DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*),
     4 AMATRX(NDOFEL,NDOFEL),RHS(MLVARX,1),DTIME,PNEWDT,PERIOD
      
      CALL ISOTROPIC(SVARS,PROPS,COORDS,U,AMATRX,RHS, NSVARS,
     1      NPROPS,MCRD,NNODE,NDOFEL,MLVARX,KINC,JELEM)
      
      END SUBROUTINE UEL
! ////////////////////////////////////////////////////////////////////////    
      
! Calculating AMATRX and RHS
      SUBROUTINE ISOTROPIC(SVARS,PROPS,COORDS,U,AMATRX,RHS, NSVARS,
     1      NPROPS,MCRD,NNODE,NDOFEL,MLVARX,KINC,JELEM)
      IMPLICIT NONE
      INTEGER NSVARS,NPROPS,MCRD,NNODE,NDOFEL,MLVARX,KINC,i,j,k,l,ipt,
     1 JELEM
      REAL(8) SVARS(NSVARS),PROPS(NPROPS),COORDS(MCRD,NNODE),U(NDOFEL),
     1 AMATRX(NDOFEL,NDOFEL),RHS(16,1),
     2 xi, eta, GAUSS(NNODE,2),det_J,WEIGHT(NNODE),Np(1,4),Np_T(4,1),
     3 Bp(2,4),Bp_T(4,2),Bu(3,2*NNODE),Bu_T(2*NNODE,3),dP(2,1),D(3,3),
     4 UU(2*NNODE,1),PHI(NNODE,1),STRAIN(3,1), STRESS(3,1), 
     5 phase,Ru(2*NNODE),Rp(NNODE),Kuu(8,8), Kpp(4,4),hist, eng, 
     6 omega, domega, ddomega, lb, GcI, E, nu, thck,
     7 Kuu11(8,3),Kuu12(8,8),Kuu1(8,8),Kpp1(4,4),Kpp11(4,4),Kpp21(4,2),
     8 Kpp22(4,4), Rp11(4,1), Ru1(8,1)
  
      ! Initialization
      AMATRX=0.d0; RHS=0.d0; det_J=0.d0; Np=0.d0; Np_T=0.d0; Bp=0.d0
      Bp_T=0.d0; Bu=0.d0; Bu_T=0.d0; dP=0.d0; D=0.d0; UU=0.d0
      PHI=0.d0; STRAIN=0.d0; STRESS=0.d0; phase=0.d0; Ru=0.d0; Rp=0.d0
      Kuu=0.d0; Kpp=0.d0; Kuu11=0.d0; hist=0.d0; eng=0.d0
      Kuu12=0.d0; Kuu1=0.d0; omega=0.d0; domega=0.d0; ddomega=0.d0
      Kpp1=0.d0; Kpp11=0.d0; Kpp21=0.d0; Kpp22=0.d0; Rp11=0.d0
! ______________________________________________________________________
! **************************************************************
! *  Constructing elemet TYPE U1( Displacement & Phase-field ) *
! **************************************************************
      ! < 1. Material parameters >
      lb=PROPS(1); GcI=PROPS(2); E=PROPS(3); nu=PROPS(4); thck=PROPS(5)
      ! ----------------------------------------------------
      ! < 2. Materials stiffness matrix (plane strain) >                
      CALL CONSTITUTIVE_MATRX(E,nu,D)
      ! --------------------------------------------------
      ! < 3. Initialization >
      DO i =1,NNODE
          ! DOF of displacement
          UU(2*i-1,1)=U(3*i-2); UU(2*i,1)=U(3*i-1)
          ! DOF of phase-field
          PHI(i,1)=U(3*i)
          END DO
      CALL GAUSS_IPT(GAUSS,WEIGHT)
      ! --------------------------------------------------
      ! 4. Calculating properties at each integral point¡ý
      LOOP_IPT :DO ipt=1,4
      ! ------------------------------< LOOP:Begin >----------------------------------
          xi  = GAUSS(ipt,1); eta = GAUSS(ipt,2)  ! Local coordinates of ipt
          ! < Shape functions and its local derivatives >
          CALL BASIC_MATRX(Np,Bu,Bp,det_J,xi,eta,COORDS)
          ! <  phase-field value and its gradient at ipt >
          phase=0.d0; dP=0.d0
          DO i=1,NNODE
              phase=phase+Np(1,i)*PHI(i,1)
              DO j=1,2
                  dP(j,1)=dP(j,1)+Bp(j,i)*PHI(i,1)
              END DO
          END DO
          
          ! < Degradation function >
          CALL DEGRADATION(phase, omega, domega, ddomega)
          
          ! < Strains & Stresses of ipt(GCS) >
          CALL UMATMUL(Bu,UU,STRAIN,3,1,8)
          CALL UMATMUL(D,STRAIN,STRESS,3,1,3)
          
          ! <  Hitory variable >
          CALL ENERGY_SPLIT(eng,STRAIN,E,nu)
          IF(eng.GT.SVARS(1*(ipt-1)+1)) THEN
              hist = eng
          ELSE
              hist = SVARS(1*(ipt-1)+1)
          END IF
          SVARS(1*(ipt-1)+1) = hist
          hist = hist/GcI
          
          
          ! < Element stiffness matrix >
          ! --------------------- Kuu ----------------------
          CALL UTRANSPOSE(Bu,Bu_T,3,8)
          CALL UMATMUL(Bu_T,D,Kuu11,8,3,3)
          CALL UMATMUL(Kuu11,Bu,Kuu12,8,8,3)
          Kuu1 = omega*Kuu12*det_J*WEIGHT(ipt)*thck
          Kuu = Kuu+Kuu1
          ! --------------------- Kpp -----------------------
          CALL UTRANSPOSE(Np,Np_T,1,4)
          CALL UTRANSPOSE(Bp,Bp_T,2,4)
          CALL UMATMUL(Np_T,Np,Kpp11,4,4,1)
          CALL UMATMUL(Bp_T,Bp,Kpp22,4,4,2)

          Kpp1 = ( (1/lb+ddomega*hist)*Kpp11 + lb*Kpp22 )*
     1        det_J*WEIGHT(ipt)*thck
          Kpp = Kpp+Kpp1
          ! < Internal forces >
          ! --------------------- Ru -----------------------
          CALL UMATMUL(Kuu1,UU,Ru1,8,1,8)
          DO i=1,8
              Ru(i)=Ru(i)-Ru1(i,1)
          END DO
          ! --------------------- Rp -----------------------
          CALL UMATMUL(Bp_T,dP,Rp11,4,1,2)
          DO i=1,NNODE
              Rp(i)=Rp(i)-( lb*Rp11(i,1)+
     1                        Np_T(i,1)*(1/lb*phase+domega*hist) )*
     1                                            WEIGHT(ipt)*det_J*thck
          END DO
          
      ! ------------------------------< LOOP:End >------------------------------------ 
      END DO LOOP_IPT
      ! 5. Updatding AMATRX & RHS
      DO i=1,4
        DO j=1,4
          AMATRX(3*i-2,3*j-2)=Kuu(2*i-1,2*j-1);
          AMATRX(3*i-1,3*j-2)=Kuu(2*i,2*j-1)
          AMATRX(3*i-2,3*j-1)=Kuu(2*i-1,2*j)
          AMATRX(3*i-1,3*j-1)=Kuu(2*i,2*j)
          AMATRX(3*i,3*j)=Kpp(i,j)
        END DO
      END DO
      
      DO i=1,NNODE
          RHS(3*i-2,1)=Ru(2*i-1); RHS(3*i-1,1)=Ru(2*i)
          RHS(3*i,1)  =Rp(i)
      END DO
      
      END SUBROUTINE ISOTROPIC

      
      
      
      
! Other subroutines might be used:
! < 1 >
!______________________________________________________________________
      SUBROUTINE GAUSS_IPT(GAUSS,WEIGHT)
      IMPLICIT NONE
      REAL(8) GAUSS(4,2),WEIGHT(4), x
      x = 0.577350269189626
      GAUSS(1,1)=-x; GAUSS(1,2)=-x
      GAUSS(2,1)= x; GAUSS(2,2)=-x
      GAUSS(3,1)= x; GAUSS(3,2)= x
      GAUSS(4,1)=-x; GAUSS(4,2)= x
      WEIGHT=1.d0
      END SUBROUTINE GAUSS_IPT
!______________________________________________________________________
      
! < 2 >
!______________________________________________________________________
      SUBROUTINE BASIC_MATRX(N,Bu,Bp,det_J,xi,eta,COORDS)
! Calculate the shape function and its derivatives
      IMPLICIT NONE
      INTEGER i,j
      REAL(8) N(1,4),dNdx(2,4),dNdxi(2,4),COORDS(2,4)
      REAL(8) JACOBI(2,2),INV_JACOBI(2,2), det_J
      REAL(8) Bu(3,8),Bp(2,4)
      REAL(8) xi,eta
      ! Initialization
      N=0.d0; dNdx=0.d0; dNdxi=0.d0; JACOBI=0.d0
      INV_JACOBI=0.d0; det_J=0.d0; Bu=0.d0; Bp=0.d0
      ! Values of shape functions as a function of local coord.
      N(1,1) = 1.d0/4.d0*(1.d0-xi)*(1.d0-eta)
      N(1,2) = 1.d0/4.d0*(1.d0+xi)*(1.d0-eta)
      N(1,3) = 1.d0/4.d0*(1.d0+xi)*(1.d0+eta)
      N(1,4) = 1.d0/4.d0*(1.d0-xi)*(1.d0+eta)
      ! Derivatives of shape functions respect to local coordinates
      dNdxi(1,1) = 1.d0/4.d0*(-1.d0+eta)
      dNdxi(1,2) = 1.d0/4.d0*(1.d0-eta)
      dNdxi(1,3) = 1.d0/4.d0*(1.d0+eta)
      dNdxi(1,4) = 1.d0/4.d0*(-1.d0-eta)
      dNdxi(2,1) = 1.d0/4.d0*(-1.d0+xi)
      dNdxi(2,2) = 1.d0/4.d0*(-1.d0-xi)
      dNdxi(2,3) = 1.d0/4.d0*(1.d0+xi)
      dNdxi(2,4) = 1.d0/4.d0*(1.d0-xi)
      ! Jacobi matrix
      DO j=1,4
          JACOBI(1,1)=JACOBI(1,1)+dNdxi(1,j)*COORDS(1,j)
          JACOBI(1,2)=JACOBI(1,2)+dNdxi(1,j)*COORDS(2,j)
          JACOBI(2,1)=JACOBI(2,1)+dNdxi(2,j)*COORDS(1,j)
          JACOBI(2,2)=JACOBI(2,2)+dNdxi(2,j)*COORDS(2,j)
      END DO
      ! Jacobi determinant
      det_J=JACOBI(1,1)*JACOBI(2,2)-JACOBI(1,2)*JACOBI(2,1)
      ! Inverse of JACOBI
      INV_JACOBI(1,1)= JACOBI(2,2)/det_J
      INV_JACOBI(1,2)=-JACOBI(1,2)/det_J
      INV_JACOBI(2,1)=-JACOBI(2,1)/det_J
      INV_JACOBI(2,2)= JACOBI(1,1)/det_J
      ! dNdx
      CALL UMATMUL(INV_JACOBI,dNdxi,dNdx,2,4,2)
      ! Bp & Bu matrix( Bp=dNdx )
      DO i=1,4
          Bp(1,i)=dNdx(1,i);      Bp(2,i)=dNdx(2,i)
          Bu(1,i*2-1)=dNdx(1,i);  Bu(2,i*2)  =dNdx(2,i)
          Bu(3,i*2-1)=dNdx(2,i);  Bu(3,i*2)  =dNdx(1,i)
      END DO
      END SUBROUTINE BASIC_MATRX
!______________________________________________________________________
  
! < 3 >
!______________________________________________________________________
      SUBROUTINE CONSTITUTIVE_MATRX(E,nu,D)
      IMPLICIT NONE
      REAL(8) E, nu, D(3,3)
      D=0.d0
      ! Planae stress
      !D(1,1)=E/(1.d0-nu**2)
      !D(2,2)=E/(1.d0-nu**2)
      !D(3,3)=E/(2*(1.d0+nu))
      !D(1,2)=E*nu/(1.d0-nu**2)
      !D(2,1)=E*nu/(1.d0-nu**2)
      
      ! Plane strain
      D(1,1)=E/((1.d0+nu)*(1.d0-2.d0*nu))*(1.d0-nu)
      D(2,2)=E/((1.d0+nu)*(1.d0-2.d0*nu))*(1.d0-nu)
      D(3,3)=E/((1.d0+nu)*(1.d0-2.d0*nu))*(0.5-nu)
      D(1,2)=E/((1.d0+nu)*(1.d0-2.d0*nu))*nu
      D(2,1)=E/((1.d0+nu)*(1.d0-2.d0*nu))*nu
      END SUBROUTINE CONSTITUTIVE_MATRX
!______________________________________________________________________ 
      
! < 4 >
!______________________________________________________________________
      SUBROUTINE DEGRADATION(phase, omega, domega, ddomega)
      IMPLICIT NONE
      REAL(8) phase, omega, domega, ddomega
      omega = (1-phase)**2
      domega= -2.d0*(1-phase)
      ddomega= 2.d0
      IF (omega.GT.1) THEN
          omega=1.d0
      ELSE IF(omega.LT.0.d0) THEN
          omega=0.d0
      END IF
      END SUBROUTINE DEGRADATION
!______________________________________________________________________ 
      
! < 5 >
!______________________________________________________________________
      SUBROUTINE ENERGY_SPLIT(ENERGY_PLUS,STRAIN,E,nu)
      IMPLICIT NONE
      REAL(8) ENERGY_PLUS, STRAIN(3,1), epsx, epsy, gamma ,p ,r, E, nu, 
     1 lambda, mu, ep1, ep2, ep, ep1_plus, ep2_plus, ep_plus

      ! Lame constants
      lambda=E*nu/(1.d0+nu)/(1.d0-2.d0*nu);  mu=E/(2.d0*(1.d0+nu))
      epsx=STRAIN(1,1);  epsy=STRAIN(2,1); gamma=STRAIN(3,1)
      p=(epsx+epsy)/2.d0
      r=sqrt(((epsx-epsy)/2.d0)**2.d0+(gamma/2.d0)**2.d0)
      ep1=p+r; ep2=p-r; ep=ep1+ep2
      CALL MACAULAY(ep1, ep1_plus)
      CALL MACAULAY(ep2, ep2_plus)
      CALL MACAULAY(ep, ep_plus)

      ENERGY_PLUS=0.5*lambda*ep_plus**2.d0+
     1                            mu*(ep1_plus**2.d0+ep2_plus**2.d0)
      END SUBROUTINE ENERGY_SPLIT
!______________________________________________________________________

! < 6 >
!______________________________________________________________________
      SUBROUTINE MACAULAY(x,y)
      IMPLICIT NONE
      REAL(8) x,y
      y=0.5*(DABS(x)+x)
      END SUBROUTINE MACAULAY
!______________________________________________________________________

      
! < 7 >
!______________________________________________________________________      
      SUBROUTINE UMATMUL(A,B,C,m,n,l)
! User defined matrix multiplication
      IMPLICIT NONE
      INTEGER i,j,k,m,n,l
      REAL(8) A(m,l),B(l,n),C(m,n)
      DO i=1,m
        DO j=1,n
            C(i,j)=0.d0
            DO k=1,l
                C(i,j)=C(i,j)+A(i,k)*B(k,j)
            END DO
        END DO
      END DO
      END SUBROUTINE UMATMUL
!______________________________________________________________________

! < 8 >
!______________________________________________________________________
      SUBROUTINE UTRANSPOSE(A,A_T,m,n)                                
! User defined matrix transpose
      IMPLICIT NONE
      INTEGER i,j,m,n
      REAL(8) A(m,n),A_T(n,m)
      DO i=1,m
          DO j=1,n
              A_T(j,i)=A(i,j)
          END DO
      END DO
      END SUBROUTINE UTRANSPOSE
!______________________________________________________________________