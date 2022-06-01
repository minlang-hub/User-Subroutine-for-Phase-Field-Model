      SUBROUTINE VUEL(nblock,rhs,amass,dtimeStable,svars,nsvars,
     &                energy,
     &                nnode,ndofel,props,nprops,jprops,njprops,
     &                coords,mcrd,u,du,v,a,
     &                jtype,jElem,
     &                time,period,dtimeCur,dtimePrev,kstep,kinc,
     &                lflags,
     &                dMassScaleFactor,
     &                predef,npredef,
     &                jdltyp, adlmag)
C
      include 'vaba_param.inc'

C     operational code keys
      parameter ( jMassCalc            = 1,
     &            jIntForceAndDtStable = 2,
     &            jExternForce         = 3)

C     flag indices
      parameter (iProcedure = 1,
     &           iNlgeom    = 2,
     &           iOpCode    = 3,
     &           nFlags     = 3)

C     energy array indices
      parameter ( iElPd = 1,
     &            iElCd = 2,
     &            iElIe = 3,
     &            iElTs = 4,
     &            iElDd = 5,
     &            iElBv = 6,
     &            iElDe = 7,
     &            iElHe = 8,
     &            iElKe = 9,
     &            iElTh = 10,
     &            iElDmd = 11,
     &            iElDc = 12,
     &            nElEnergy = 12)

C     predefined variables indices
      parameter ( iPredValueNew = 1,
     &            iPredValueOld = 2,
     &            nPred         = 2)

C     time indices
      parameter (iStepTime  = 1,
     &           iTotalTime = 2,
     &           nTime      = 2)

		   dimension rhs(nblock,ndofel),amass(nblock,ndofel,ndofel),
     &          dtimeStable(nblock),
     &          svars(nblock,nsvars),energy(nblock,nElEnergy),
     &          props(nprops),jprops(njprops),
     &          jElem(nblock),time(nTime),lflags(nFlags),
     &          coords(nblock,nnode,mcrd),
     &          u(nblock,ndofel), du(nblock,ndofel),
     &          v(nblock,ndofel), a(nblock, ndofel),
     &          dMassScaleFactor(nblock),  
     &          predef(nblock,nnode,npredef,nPred),
     &          adlmag(nblock)
      
      ! 4 node quadrilateral 8 DODs
      ! local subroutines variables
      ! displacement field
      real(8):: dis_mass(8,8),dd_mass(4,4)

      real(8):: dis_rhs(8,1),dd_rhs(4,1)

      real(8):: coordinate(2,4), copy_props(8)
      real(8):: sub_u(8,1)
      real(8):: sub_d(4,1)
      real(8):: val_svars(1)

      integer:: i, j, k
 
      ! properties
      do i=1,8
          copy_props(i)=props(i)
      enddo
      
      if ( lflags(iOpCode).eq.jMassCalc ) then
          ! define mass matrice
          do kblock = 1,nblock
              ! coordinate 1、2、3、4
                          !x1 x2 x3 x4
                          !y1 y2 y3 y4
              do i=1,2
                  do j=1,4                  !单元号  节点号  x/y坐标
                      coordinate(i,j)=coords(kblock,j,i)
                  enddo
              enddo

              call gen_mass_matrix(dis_mass,dd_mass,
     & copy_props,coordinate)

              do i=1,4
                  amass(kblock,3*i-2,3*i-2)=dis_mass(2*i-1,2*i-1)   !对角阵
                  amass(kblock,3*i-1,3*i-1)=dis_mass(2*i,2*i)
                  amass(kblock,3*i,3*i)    =dd_mass(i,i)
              enddo
          enddo          
      elseif ( lflags(iOpCode) .eq. jIntForceAndDtStable) then
          !defined RHS
          do kblock = 1,nblock
              ! coordinate
              do i=1,2
                  do j=1,4
                      coordinate(i,j)=coords(kblock,j,i)
                  enddo
              enddo
              
              do i=1,4
                  sub_u(2*i-1,1)=u(kblock,3*i-2)
                  sub_u(2*i,1)  =u(kblock,3*i-1)
                  sub_d(i,1)    =u(kblock,3*i)
              enddo
              
              !====== must be refined when considering multiple integral points =========
              do i=1,1
                  val_svars(i)=svars(kblock,i)    !val_svars四个状态变量（四个积分点）
              enddo
              
 
              call gen_rhs_vector(dis_rhs,dd_rhs,copy_props,
     &                                coordinate,sub_u,sub_d,val_svars)
              
              do i=1,1
                  svars(kblock,i)=val_svars(i)    !更新状态变量
              enddo
              
              
              do i=1,4
                  rhs(kblock,3*i-2)=dis_rhs(2*i-1,1)    !将dis_rhs、dd_rhs按要求位置放入rhs中
                  rhs(kblock,3*i-1)=dis_rhs(2*i,1)
                  rhs(kblock,3*i)  =dd_rhs(i,1)
              enddo
              
              ! stable time step
              dtimeStable(kblock)=copy_props(8)
              energy(kblock, iElKe) = 0.d0    !动能为零
              energy(kblock, iElIe) = 0.d0    !Internal energy为零

          enddo
      endif
      
      END subroutine VUEL
!====================== subroutines ==============================================================
      !====================== generate mass matrix ===============================================
      subroutine gen_mass_matrix(dis_mass,dd_mass,props,coords)
      ! props(8)  [E mu thickness Gf lb v_eta ro time_step ]
      !            1  2     3      4  5   6    7     8        
                                                !密度
      implicit none
      real(8):: dis_mass(8,8), dd_mass(4,4)
      real(8):: props(8),coords(2,4)
      
      real(8):: dis_mass1(8,8), dd_mass1(4,4)
      
      ! local variables
      real(8):: guassp(1), gwight(1), ro, thickness
      integer:: numgus
      real(8):: jacb(2,2), inv_jacb(2,2), det_jacb
      real(8):: nu(2,8), b(3,8)
      real(8):: nd(1,4),bd(2,4)
      real(8):: v_eta
      integer:: i,j,k

      real(8):: Gf, lb

      ! initialization
      dis_mass=0.d0; dd_mass=0.d0
      
      dis_mass1=0.d0; dd_mass1=0.d0
      
      jacb=0.d0; inv_jacb=0.d0; det_jacb=0.d0; b=0.d0; nu=0.d0
      nd=0.d0; bd=0.d0

      ! guass integration
      guassp(1)=0.d0; gwight(1)=2.d0
      numgus=1
      !numgus=2
      !guassp(1)=-1.d0/dsqrt(3.d0); guassp(2)=1.d0/dsqrt(3.d0)
      !gwight(1)=1.d0;  gwight(2)=1.d0
      
      ro=props(7); thickness=props(3); v_eta=props(6)
      Gf=props(4)
      lb=props(5)

      
      
      do i=1,numgus
          do j=1,numgus
              call jacob_matrix(jacb,inv_jacb,det_jacb,guassp(i),
     &                            guassp(j),coords)
              call b_matrix(nu,nd,bd,b,inv_jacb,guassp(i),guassp(j))

              dis_mass1=dis_mass1+gwight(i)*gwight(j)*
     &               (ro*matmul(transpose(nu),nu))*det_jacb*thickness
              dd_mass1=dd_mass1+gwight(i)*gwight(j)*
     &             (v_eta*matmul(transpose(nd),nd))*det_jacb*thickness
              
          enddo
      enddo
      
      do i=1,8
          do j=1,8
              dis_mass(i,i)=dis_mass(i,i)+dis_mass1(i,j)
          enddo
      enddo
      
      do i=1,4
          do j=1,4
              dd_mass(i,i)=dd_mass(i,i)+dd_mass1(i,j)
          enddo
      enddo
      end subroutine gen_mass_matrix
      !========================= generate RHS ============================================
      subroutine gen_rhs_vector(dis_rhs,dd_rhs,props,coords,u,sub_d,
     & val_svars)
      ! props(8)  [E mu thickness Gf lb v_eta ro time_step ]
      !            1  2     3      4  5   6    7     8    
      implicit none
      real(8):: dis_rhs(8,1),dd_rhs(4,1),props(8),coords(2,4),u(8,1)
      real(8):: sub_d(4,1)
      real(8):: val_svars(1)
      
      ! local variables
      real(8):: guassp(1), gwight(1), thickness
      integer:: numgus
      real(8):: d(3,3), jacb(2,2), inv_jacb(2,2), det_jacb 
      real(8):: b(3,8), nu(2,8)
      real(8):: nd(1,4), bd(2,4)
      real(8):: uu(8,1), dd(4,1)
      integer:: i,j,k
      
      ! phase field problem
      !==================================================
      real(8):: ppf(1,1)
      real(8):: pf, energy_plus
      real(8):: dalpha,ddalpha,c0
      real(8):: omega,domega,ddomega
              
      real(8):: Gf, lb
              
      real(8):: H_wave
      real(8):: eff_stress(3,1)
      real(8):: mul_bbd(4,1)
      
      real(8):: dd_mass1(4,4)
      real(8):: dd_mass(4,4)
      real(8):: d_c(4,1), dd1(4,1), r_c(4,1)
      real(8):: d_c0(4,1), dd0(4,1), r_c0(4,1)
      real(8):: v_eta
      real(8):: stress_plus(3,1), stress_minus(3,1)
      
      dd_mass1=0.d0
      dd_mass=0.d0
      c0=2.d0        
      !======================================

      ! initialization
      dis_rhs=0.d0; dd_rhs=0.d0
      d=0.d0; jacb=0.d0; inv_jacb=0.d0; det_jacb=0.d0
      nu=0.d0; b=0.d0
      nd=0.d0; bd=0.d0
      
      ! guass integration
      guassp(1)=0.d0; gwight(1)=2.d0
      numgus=1
      !numgus=2
      !guassp(1)=-1.d0/dsqrt(3.d0); guassp(2)=1.d0/dsqrt(3.d0)
      !gwight(1)=1.d0;  gwight(2)=1.d0
      
      thickness=props(3)
      Gf=props(4); lb=props(5); v_eta=props(6)
      
      
      do i=1,8
          uu(i,1)=u(i,1)        !!!!!!!!!!!!!!!!!
      enddo  
      
      do i=1,4
          dd(i,1)=sub_d(i,1)
          if (dd(i,1)>1.d0) then
              dd(i,1)=1.d0
          elseif (dd(i,1)<0.d0) then
              dd(i,1)=0.d0
          endif
      enddo
      
      do i=1,numgus
          do j=1,numgus
              call jacob_matrix(jacb,inv_jacb,det_jacb,guassp(i),
     &                            guassp(j),coords)
              call b_matrix(nu,nd,bd,b,inv_jacb,guassp(i),guassp(j))
              call d_matrix(d,props)
 
              ppf=matmul(nd,dd)
              pf=ppf(1,1)
              call gener_engeryp(stress_plus,stress_minus,
     &             energy_plus,b,uu,props)
              
              if (val_svars(2*(i-1)+j)>energy_plus) then
                  energy_plus=val_svars(2*(i-1)+j)
              else
                  val_svars(2*(i-1)+j)=energy_plus    !H
              endif
              
              call alpha_dd(dalpha,pf,ddalpha)
              call degradfunc(omega,domega,ddomega,pf)
              
              H_wave=domega*energy_plus+Gf/c0/lb*dalpha
              
              mul_bbd=matmul(transpose(bd),matmul(bd,dd))
              
              dd_rhs=dd_rhs+gwight(i)*gwight(j)*( 2.d0*Gf/c0*lb*mul_bbd
     &         +H_wave*transpose(nd) )*det_jacb*thickness
              
     
              if (omega>1.d0) then
                  omega=1.d0
              elseif (omega<0.d0) then
                  omega=0.d0
              endif
              
              eff_stress=omega*matmul(d,matmul(b,uu))      !有效应力
              !eff_stress=omega*stress_plus+stress_minus
              dis_rhs=dis_rhs+gwight(i)*gwight(j)*
     &        (matmul(transpose(b),eff_stress))*det_jacb*thickness
     
              dd_mass1=dd_mass1+gwight(i)*gwight(j)*
     &             (v_eta*matmul(transpose(nd),nd))*det_jacb*thickness
          enddo
      enddo
      
      do i=1,4
          do j=1,4
              dd_mass(i,i)=dd_mass(i,i)+dd_mass1(i,j)
          enddo
      enddo
      
      do i=1,4
          dd1(i,1)=1.d0
      enddo

      d_c=dd1-dd;                      !!!!!!!!
      d_c=d_c/props(8)
      
      r_c=matmul(dd_mass,d_c)
       
      do i=1,4
          if (-dd_rhs(i,1)>r_c(i,1)) then
              dd_rhs(i,1)=-r_c(i,1)
          endif
      enddo
      
      do i=1,4
          dd0(i,1)=0.d0
      enddo

      d_c0=dd0-dd;                      !!!!!!!!
      d_c0=d_c0/props(8)
      
      r_c0=matmul(dd_mass,d_c0)
       
      do i=1,4
          if (-dd_rhs(i,1)<r_c0(i,1)) then
              dd_rhs(i,1)=-r_c0(i,1)
          endif
      enddo      

      end subroutine gen_rhs_vector
      !================traditional material propoty matrix================================
      subroutine d_matrix(d,props)
      ! props(8)  [E mu thickness Gf lb v_eta ro time_step ]
      !            1  2     3      4  5   6    7     8    
      implicit none
      real(8):: d(3,3), props(8)
      real(8):: e, mu, a1, a2
   
      d=0.d0
          
      e=props(1);  mu=props(2)
      a1=e/(1.d0-mu*mu)
      a2=(1.d0-mu)/2.d0
     
      d(1,1)=1.d0;     d(1,2)=mu
      d(2,1)=d(1,2);   d(2,2)=1.d0
                                      d(3,3)=a2
      d=a1*d
      end subroutine d_matrix      
      !==================shape function and its derivative with xi and eta======================    
      subroutine shapefuc(n,dn_xieta,xi,eta)      
      implicit none      
      real(8):: n(4), dn_xieta(2,4),xi, eta

      n(1)=1.d0/4.d0*(1.d0-xi)*(1.d0-eta)
      n(2)=1.d0/4.d0*(1.d0+xi)*(1.d0-eta)
      n(3)=1.d0/4.d0*(1.d0+xi)*(1.d0+eta)
      n(4)=1.d0/4.d0*(1.d0-xi)*(1.d0+eta)
           
!     dn_xieta=[ n1_xi   n2_xi   n3_xi   n4_xi
!                n1_eta  n2_eta  n3_eta  n4_eta]
      dn_xieta(1,1)=-1.d0/4.d0*(1.d0-eta)
      dn_xieta(1,2)= 1.d0/4.d0*(1.d0-eta)
      dn_xieta(1,3)= 1.d0/4.d0*(1.d0+eta)
      dn_xieta(1,4)=-1.d0/4.d0*(1.d0+eta)
      
      dn_xieta(2,1)=-1.d0/4.d0*(1.d0-xi)
      dn_xieta(2,2)=-1.d0/4.d0*(1.d0+xi)
      dn_xieta(2,3)= 1.d0/4.d0*(1.d0+xi)
      dn_xieta(2,4)= 1.d0/4.d0*(1.d0-xi)
      end subroutine shapefuc
      !======================jacob matrix=========================================================     
      subroutine jacob_matrix(jacb,inv_jacb,det_jacb,xi,eta,coords)
      implicit none
      real(8):: jacb(2,2), inv_jacb(2,2), det_jacb, xi, eta, coords(2,4)
      
!     local varibles
      real(8):: n(4), dn_xieta(2,4)
      integer:: i
      
      jacb=0.d0; inv_jacb=0.d0; det_jacb=0.d0
      
      call shapefuc(n,dn_xieta,xi,eta)
      
      do i=1,4
          jacb(1,1)=jacb(1,1)+dn_xieta(1,i)*coords(1,i)
          jacb(1,2)=jacb(1,2)+dn_xieta(1,i)*coords(2,i)
          jacb(2,1)=jacb(2,1)+dn_xieta(2,i)*coords(1,i)
          jacb(2,2)=jacb(2,2)+dn_xieta(2,i)*coords(2,i)
      enddo
      
      det_jacb=jacb(1,1)*jacb(2,2)-jacb(1,2)*jacb(2,1)
      
      inv_jacb(1,1)=jacb(2,2);  inv_jacb(1,2)=-jacb(1,2)
      inv_jacb(2,1)=-jacb(2,1); inv_jacb(2,2)=jacb(1,1)
      
      inv_jacb=1.d0/det_jacb*inv_jacb
      end subroutine jacob_matrix
!===============traditional b matrix==============================================      
      subroutine b_matrix(nu,nd,bd,b,inv_jacb,xi,eta)
      implicit none
      real(8):: nd(1,4), bd(2,4), b(3,8), inv_jacb(2,2), xi, eta
      real(8):: nu(2,8)
      
!     local varibles
      real(8):: n(4), dn_xieta(2,4), dn_x(4), dn_y(4)
      integer:: i, j
      
!     initialize varibles
      b=0.d0; nu=0.d0
      
      call shapefuc(n,dn_xieta,xi,eta)
      
      do i=1,4
          nd(1,i)=n(i)
      enddo
      
      do i=1,4
          dn_x(i)=inv_jacb(1,1)*dn_xieta(1,i)
     &             +inv_jacb(1,2)*dn_xieta(2,i)
          dn_y(i)=inv_jacb(2,1)*dn_xieta(1,i)
     &             +inv_jacb(2,2)*dn_xieta(2,i)
      enddo
      
      do j=1,4
          b(1,2*(j-1)+1)=dn_x(j)
                                   b(2,2*(j-1)+2)=dn_y(j)
          b(3,2*(j-1)+1)=dn_y(j);  b(3,2*(j-1)+2)=dn_x(j)
      enddo
      
      do j=1,4
          bd(1,j)=dn_x(j)
          bd(2,j)=dn_y(j)
      enddo
      
      do i=1,4
          nu(1,2*i-1)=n(i)
          nu(2,2*i)  =n(i)
      enddo

      end subroutine b_matrix
!===============degradation function and its derivations========================
      subroutine degradfunc(omega,domega,ddomega,pf)
      implicit none
      real(8):: omega, domega, ddomega, pf
      omega=(1-pf)**2.d0    
      domega=-2.d0*(1-pf)
      ddomega=2.d0
!      if (omega.gt.1) then
!          omega=1.d0
!      else if(omega.lt.0.d0) then
!          omega=0.d0
!      end if     
      end subroutine degradfunc
!=====================r alpha===================================================
      subroutine alpha_dd(dalpha,pf,ddalpha)
      implicit none
      real(8):: dalpha, pf, ddalpha     
      dalpha=2.d0*pf
      ddalpha=2.d0
      end subroutine alpha_dd     
!==================energy_plus==================================================
      subroutine gener_engeryp(stress_plus,stress_minus,energy_plus,
     &                                                 b,uu,props)
      ! props(8)  [E mu thickness Gf lb v_eta ro time_step ]
      !            1  2     3      4  5   6    7     8   
      implicit none
      real(8):: stress_plus(3,1), stress_minus(3,1)
      real(8):: energy_plus, b(3,8), uu(8,1), props(8)
      real(8):: epsi(3,1), strain(2,2)
      real(8):: E, mu, lamd, G
      real(8):: valu(2), eigv(2,2), en1, pn1, pn2
      real(8):: fn1, jn1, gn1, gn2
      real(8):: n1(2), n2(2)
      real(8):: n11(2,2), n22(2,2)
      real(8):: stress_p(2,2), stress_m(2,2)
      integer:: i, j
      !initialize varibles
      stress_plus=0.d0; stress_minus=0.d0
      stress_p=0.d0; stress_m=0.d0
      
      E=props(1);   mu=props(2)
      !平面应变
      !lamd=E*mu/((1.d0+mu)*(1.d0-2.d0*mu)); G=E/(2.d0*(1.d0+mu))
      !平面应力
      lamd=E*mu/(1.d0-mu**2.d0); G=E/(2.d0*(1.d0+mu))
      
      epsi=matmul(b,uu)
      strain(1,1)=epsi(1,1);   strain(1,2)=0.5d0*epsi(3,1)
      strain(2,1)=strain(1,2); strain(2,2)=epsi(2,1);
      
      call eig(strain,valu,eigv)
      
      do i=1,2
          n1(i)=eigv(i,1)
          n2(i)=eigv(i,2)
      enddo
      
      do i=1,2
          do j=1,2
              n11(i,j)=n1(i)*n1(j)
              n22(i,j)=n2(i)*n2(j)
          enddo
      enddo
      
      
      en1=valu(1)+valu(2)
      fn1=0.5d0*(en1+dabs(en1))
      jn1=0.5d0*(en1-dabs(en1))
      
      pn1=0.5d0*(dabs(valu(1))+valu(1))
      gn1=0.5d0*(valu(1)-dabs(valu(1)))
      
      pn2=0.5d0*(dabs(valu(2))+valu(2))
      gn2=0.5d0*(valu(2)-dabs(valu(2))) 
      
      stress_p(1,1)=lamd*fn1
      stress_p(2,2)=lamd*fn1
       
      stress_m(1,1)=lamd*jn1
      stress_m(2,2)=lamd*jn1
       
      do i=1,2
          do j=1,2
              stress_p(i,j)=stress_p(i,j)+(2*G*pn1)*n11(i,j)+
     &        (2*G*pn2)*n22(i,j)
              stress_m(i,j)=stress_m(i,j)+(2*G*gn1)*n11(i,j)+
     &        (2*G*gn2)*n22(i,j)       
          enddo
      enddo
           
      
      stress_plus(1,1)=stress_p(1,1)
      stress_plus(2,1)=stress_p(2,2)
      stress_plus(3,1)=stress_p(1,2)
      
      stress_minus(1,1)=stress_m(1,1)
      stress_minus(2,1)=stress_m(2,2)
      stress_minus(3,1)=stress_m(1,2)
      
      energy_plus=lamd/2.d0*(0.5d0*(dabs(en1)+en1))**2.d0
     &                +G*(pn1*pn1+pn2*pn2)      
      
      end subroutine gener_engeryp    
!==================phase decomposite============================================
      subroutine eig(a,valu,eigv)
      implicit none
      real(8):: a(2,2),ca(2,2),valu(2),eigv(2,2)
      real(8):: r(2,2),eigg(2,2)
      real(8):: tol,t,coss,sins,eps
      integer:: i,j
      
      eps=1.d-10
      
      eigv=0.d0
      eigv(1,1)=1.d0; eigv(2,2)=1.d0
      do i=1,2
          do j=1,2
              ca(i,j)=a(i,j)
          enddo
      enddo
      
      do while (dabs(ca(1,2))>eps .or. dabs(ca(2,1))>eps)
          tol=(ca(1,1)-ca(2,2))/(2.d0*ca(1,2))
          t=tol/dabs(tol)/(dabs(tol)+dsqrt(1.d0+tol**2.d0))
          coss=(1.d0+t**2.d0)**(-0.5d0)
          sins=t*coss
          r(1,1)=coss; r(1,2)=-sins
          r(2,1)=sins; r(2,2)=coss
          ca=matmul(matmul(transpose(r),ca),r)
          eigv=matmul(eigv,r)
      enddo
      
      if (ca(1,1)>ca(2,2)) then
          valu(1)=ca(1,1)
          valu(2)=ca(2,2)
      else
          valu(1)=ca(2,2)
          valu(2)=ca(1,1)
          eigg(1,1)=eigv(1,2);  eigg(1,2)=eigv(1,1)
          eigg(2,1)=eigv(2,2);  eigg(2,2)=eigv(2,1)
          do i=1,2
              do j=1,2
                  eigv(i,j)=eigg(i,j)
              enddo
          enddo
      endif

      end subroutine eig
!=============================== end of the subroutines ===========================================