      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
!
      INCLUDE 'ABA_PARAM.INC'
!
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      call element_isotropic(rhs,amatrx,props,coords,u,svars)  
      
      end subroutine UEL   
!================subroutine element_isotropic========================================================     
!!!           in: props,coords,u,svars
!!!           out: rhs,amatrx
      subroutine element_isotropic(rhs,amatrx,props,coords,u,svars)
      implicit none
!     props=[ E  mu  t  Gc ft lb] 
!     coordinate 1、2、3、4
!               x1 x2 x3 x4
!               y1 y2 y3 y4      
      real(8):: rhs(12,1), amatrx(12,12)
      real(8):: props(6), coords(2,4), u(12), svars(4)
!     local varibles
      real(8):: ku(8,8), kd(4,4)
      real(8):: ru(8,1), rd(4,1)      
      real(8):: uu(8,1)
      real(8):: dd(4,1)      
      integer:: numgus, i, j
      real(8):: Gf, lb
      real(8):: guassp(2), gwight(2), jacb(2,2), inv_jacb(2,2), det_jacb
      real(8):: nd(1,4), bd(2,4), b(3,8) 
      real(8):: ppf(1,1), pf, d(3,3), energy_plus
      real(8):: dalpha, ddalpha, c0      
      real(8):: omega, domega, ddomega 
      real(8):: dalpha0, pf0, ddalpha0, c00
      real(8):: omega0, domega0, ddomega0
      real(8):: ph_energy, d_ph_energy            
      real(8):: kd1(4,4)
      real(8):: rd1(4,1)
      real(8):: nd_t(4,1), bd_t(4,2), b_t(8,3)
      real(8):: bd_t_mul_bd(4,4), bd_t_mul_bd_mul_dd(4,1)
      real(8):: nd_t_mul_nd(4,4), b_t_mul_d(8,3)
      real(8):: b_t_mul_d_mul_b(8,8)
!     initialize varibles
      rhs=0.d0; amatrx=0.d0; ku=0.d0; kd=0.d0; ru=0.d0; rd=0.d0
      uu=0.d0; dd=0.d0; numgus=0
      guassp=0.d0; gwight=0.d0; jacb=0.d0
      inv_jacb=0.d0; det_jacb=0.d0; nd=0.d0; bd=0.d0; b=0.d0
      ppf=0.d0; pf=0.d0; d=0.d0; energy_plus=0.d0
      dalpha=0.d0; ddalpha=0.d0; c0=0.d0 
      omega=0.d0; domega=0.d0; ddomega=0.d0
      dalpha0=0.d0; pf0=0.d0; ddalpha0=0.d0; c00=0.d0
      omega0=0.d0; domega0=0.d0; ddomega0=0.d0
      ph_energy=0.d0; d_ph_energy=0.d0
      kd1=0.d0; rd1=0.d0
      nd_t=0.d0; bd_t=0.d0; b_t=0.d0
      bd_t_mul_bd=0.d0; bd_t_mul_bd_mul_dd=0.d0
      nd_t_mul_nd=0.d0; b_t_mul_d=0.d0
      b_t_mul_d_mul_b=0.d0     
      
      do i=1,4
          uu(2*i-1,1)=u(3*i-2)
          uu(2*i,1)=u(3*i-1)
          dd(i,1)=u(3*i)
      enddo
      
      Gf=props(4);  lb=props(6)        
      
      numgus=2
      guassp(1)=-0.577350269189626; guassp(2)=0.577350269189626        
      gwight(1)=1.d0;  gwight(2)=1.d0
          
      do i=1,numgus
          do j=1,numgus 
              !jacob、inv_jacob、det_jacob
              call jacob_matrix(jacb,inv_jacb,det_jacb,guassp(i),
     &                            guassp(j),coords)    
              !nd,bd,b
              call b_matrix(nd,bd,b,inv_jacb,guassp(i),guassp(j))     
              
              call matmul_uel(nd,dd,ppf,1,1,4)
              pf=ppf(1,1)        
              
              call d_matrix(d,props)        
          
              call gener_engeryp(energy_plus,b,uu,props)         
              if (svars(2*(i-1)+j)>energy_plus) then
                  energy_plus=svars(2*(i-1)+j)
              else
                  svars(2*(i-1)+j)=energy_plus
              endif

              call alpha_dd(dalpha,pf,ddalpha,c0)
              call degradfunc(omega,domega,ddomega,pf,props,c0)
                      
              pf0=0.d0
              call alpha_dd(dalpha0,pf0,ddalpha0,c00)
              call degradfunc(omega0,domega0,ddomega0,pf0,props,c00)                                               
              if (energy_plus<-2.d0*Gf/c0/lb/domega0) then
                  energy_plus=-2.d0*Gf/c0/lb/domega0
              endif
     
              ph_energy=domega*energy_plus+Gf/c0*1.d0/lb*dalpha
              d_ph_energy=ddomega*energy_plus+Gf/(c0*lb)*ddalpha

              call transpose_uel(nd,nd_t,1,4)      
              call transpose_uel(bd,bd_t,2,4)  
              call transpose_uel(b,b_t,3,8)  
              call matmul_uel(bd_t,bd,bd_t_mul_bd,4,4,2)
              call matmul_uel(bd_t_mul_bd,dd,bd_t_mul_bd_mul_dd,4,1,4)
              call matmul_uel(nd_t,nd,nd_t_mul_nd,4,4,1)
              call matmul_uel(b_t,d,b_t_mul_d,8,3,3)
              call matmul_uel(b_t_mul_d,b,b_t_mul_d_mul_b,8,8,3)  
                          
              rd1=-gwight(i)*gwight(j)*
     &            (ph_energy*nd_t+Gf/c0*(2.d0*lb*
     &        bd_t_mul_bd_mul_dd))*det_jacb*props(3)
              rd=rd+rd1
              
              kd1=gwight(i)*gwight(j)*
     &        ((d_ph_energy)*
     &        nd_t_mul_nd+2.d0*lb*Gf/c0*
     &        bd_t_mul_bd)*det_jacb*props(3)
              kd=kd+kd1
     
     
              ku=ku+gwight(i)*gwight(j)*
     &        (omega*b_t_mul_d_mul_b)*det_jacb*props(3)
       
          enddo
      enddo
      call matmul_uel(-ku,uu,ru,8,1,8)
      
      do i=1,4
        do j=1,4
          amatrx(3*i-2,3*j-2)=ku(2*i-1,2*j-1);
          amatrx(3*i-1,3*j-2)=ku(2*i,2*j-1)
          amatrx(3*i-2,3*j-1)=ku(2*i-1,2*j)
          amatrx(3*i-1,3*j-1)=ku(2*i,2*j)
          amatrx(3*i,3*j)=kd(i,j)
        end do
      end do
      
      do i=1,4
          rhs(3*i-2,1)=ru(2*i-1,1)
          rhs(3*i-1,1)=ru(2*i,1)
          rhs(3*i,1)  =rd(i,1)
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !write(*,*)'debug',rhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      end subroutine element_isotropic
!================subroutines========================================================
!================traditional material propoty matrix================================
      subroutine d_matrix(d,props)
      implicit none
      real(8):: d(3,3), props(6)
      real(8):: e, mu, a1, a2
      
      d=0.d0
      e=props(1);  mu=props(2)
      a1=e/(1.d0-mu*mu)
      a2=(1.d0-mu)/2.d0
      
      d(1,1)=1.d0; d(1,2)=mu
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
!!!                 in: xi,eta,coords
!!!            out: jacb,inv_jacb,det_jacb
      subroutine jacob_matrix(jacb,inv_jacb,det_jacb,xi,eta,coords)
      implicit none
      real(8):: jacb(2,2), inv_jacb(2,2), det_jacb, xi, eta, coords(2,4)
      
!     local varibles
      real(8):: n(4), dn_xieta(2,4)
      integer:: i
      
      jacb=0.d0
      
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
      subroutine b_matrix(nd,bd,b,inv_jacb,xi,eta)
      implicit none
      real(8):: nd(1,4), bd(2,4), b(3,8), inv_jacb(2,2), xi, eta
      
!     local varibles
      real(8):: n(4), dn_xieta(2,4), dn_x(4), dn_y(4)
      integer:: i, j
      
!     initialize varibles
      b=0.d0
      
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
      
      end subroutine b_matrix    
!===============degradation function and its derivations========================
      subroutine degradfunc(omega,domega,ddomega,pf,props,c0)
      implicit none
      real(8):: omega, domega, ddomega, pf, props(6), c0
!     local varibles
      real(8):: g, dg, ddg, fai, dfai, ddfai, E, Gf, ft, lb
      real(8):: p, a1, a2, a3

      E=props(1); Gf=props(4); ft=props(5); lb=props(6)
      a1=4.d0/(c0*lb)*E*Gf/(ft*ft)
      a2=1.3868d0;  p=2.d0
      a3=0.6567d0;
      !a2=-0.5d0;  p=2.d0
      !a3=0.d0;
      
      g=(1.d0-pf)**p; dg=-p*(1.d0-pf)**(p-1.d0); 
      ddg=p*(p-1.d0)*(1.d0-pf)**(p-2.d0)
      
      fai=g+a1*pf+a1*a2*pf**2.d0+a1*a2*a3*pf**3.d0
      dfai=dg+a1+2.d0*a1*a2*pf+3.d0*a1*a2*a3*pf**2.d0
      ddfai=ddg+2.d0*a1*a2+6.d0*a1*a2*a3*pf
      
      omega=g/fai
      
      domega=(dg*fai-g*dfai)/(fai**2.d0)
      ddomega=((ddg*fai-g*ddfai)*fai-2.d0*(dg*fai-g*dfai)*dfai)
     &        /(fai**3.d0)

      if (omega.gt.1) then
          omega=1.d0
      else if(omega.lt.0.d0) then
          omega=0.d0
      end if          
      end subroutine degradfunc
!=====================r alpha===================================================
      subroutine alpha_dd(dalpha,pf,ddalpha,c0)
      implicit none
      real(8):: dalpha, pf, ddalpha, c0
      
      c0=3.1415926535897932384626433832795028841971693993751d0
      dalpha=2.d0-2.d0*pf
      ddalpha=-2.d0
      end subroutine alpha_dd     
!==================energy_plus==================================================
      subroutine gener_engeryp(energy_plus,b,uu,props)
      implicit none
      real(8):: energy_plus, b(3,8), uu(8,1), props(6)
      real(8):: epsi(3,1)   
      real(8):: E, mu, lamd, G
      real(8):: epsx, epsy, gamma, p, r
      real(8):: valu(2), en1, pn1, pn2      
      E=props(1);   mu=props(2)
!!!!plane stress
      lamd=E*mu/(1.d0-mu**2.d0); G=E/(2.d0*(1.d0+mu))    
!!!!plane strain
      !lamd=E*mu/((1.d0+mu)*(1.d0-2.d0*mu)); G=E/(2.d0*(1.d0+mu)) 
      epsi=matmul(b,uu)
      epsx=epsi(1,1); epsy=epsi(2,1); gamma=epsi(3,1)
      p=(epsx+epsy)/2.d0
      r=sqrt(((epsx-epsy)/2.d0)**2.d0+(gamma/2.d0)**2.d0)
      valu(1)=p+r; valu(2)=p-r
      en1=valu(1)+valu(2)
      
      pn1=0.5d0*(dabs(valu(1))+valu(1))
      pn2=0.5d0*(dabs(valu(2))+valu(2))     
      energy_plus=lamd/2.d0*(0.5d0*(dabs(en1)+en1))**2.d0
     &                +G*(pn1*pn1+pn2*pn2)
      end subroutine gener_engeryp        
!==================matrix multiplication=======================================
      subroutine matmul_uel(a,b,c,m,n,l)
      implicit none
      integer i, j, k, m, n, l
      real(8):: a(m,l),b(l,n),c(m,n)
      do i=1,m
        do j=1,n
            c(i,j)=0.d0
            do k=1,l
                c(i,j)=c(i,j)+a(i,k)*b(k,j)
            end do
        end do
      end do 
      end subroutine matmul_uel
!==================matrix transpose=======================================
      subroutine transpose_uel(a,a_t,m,n)           
      implicit none
      integer i,j,m,n
      real(8):: a(m,n),a_t(n,m)
      do i=1,m
          do j=1,n
              a_t(j,i)=a(i,j)
          end do
      end do
      end subroutine transpose_uel
!==================end of the subroutines=======================================
