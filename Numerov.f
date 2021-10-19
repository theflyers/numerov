cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module mesh 
c     用来储存全局变量中关于格点（mesh）的信息
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       implicit none 
       real*8 :: hcm ! Numerov方法中的步长
       integer :: irmatch !Numerov方法中最大的格点数
       integer :: nr ! 高斯积分需要用的变量
       real*8,allocatable,dimension(:) :: rr,rrw ! 高斯积分需要用的变量
      end module 
c-----------------------------------------------------------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module systems 
c      两体系统的变量
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       implicit none 
       real*8 :: z1,z2 ! 两体系统的电荷数
       real*8 :: mass1, mass2 ! 两体系统的质量数
       real*8 :: be ! 结合能 be < 0 
      end module 
c-----------------------------------------------------------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module constants
c      常用物理量
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       implicit none
       real*8,parameter :: hbarc=197.3269718d0         ! NIST Ref 02.12.2014   ! MeV.fm           
       real*8,parameter :: finec=137.03599d0
       real*8,parameter :: amu=931.49432d0  !MeV
       real*8,parameter :: e2=1.43997d0         ! MeV.fm
       real*8,parameter :: zero=0.0d0 
       real*8,parameter :: convert=1.0d0    ! convert integer to real
      end module constants
c-----------------------------------------------------------------------

      
      program np_bound_r_space
      use mesh 
      use systems
      implicit none 
      
      ! 定义格点信息
      hcm=0.01 ! fm 
      irmatch=2000 
      nr=40
      allocate(rr(1:nr),rrw(1:nr))
      call gauleg(nr,0.0d0,hcm*irmatch,rr,rrw) ! 初始化高斯积分
      
      !定义系统信息
      z1=1.0
      z2=0.0
      mass1=1.0078 
      mass2=1.0086 
      be=-2.224 ! MeV
      
      call npbound()


      
      end program 
      
      
      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine npbound()
c      求解薛定谔方程确定np相互作用势的深度
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh
       use systems
       use constants
       implicit none
       integer :: ir,i,irmid ! 格点信息
       real*8 :: mu,ecm ! 折合质量， 与质心系能量
       real*8 :: norm !  归一化系数
       real*8 :: eta ! 索末菲参数
       integer :: IE ! Whittaker function系数
       real*8,dimension(0:0) :: WK,WKD
       real*8 :: knp ! wave number 
       real*8,dimension(0:irmatch) :: vpot,vpot1  ! np相互作用势
       real*8,dimension(0:irmatch) ::  phi ! np波函数
       real*8 :: r 
       real*8 :: V0 ! 初始势阱深度
       real*8 :: gausspot ! 高斯势函数
       real*8 :: fpin,fpout
       real*8 :: delta
       real*8 :: n !势阱深度的系数
       real*8 :: const
       real*8 :: fl !  matching point wave function
       real*8 :: FFR4 ! interpolation function
       real*8,dimension(0:irmatch) :: fin
       real*8,dimension(1:irmatch) :: kl
       real*8 :: f0,phivphi
       integer :: l


       phi=0.0d0
       vpot=0.0d0
       
       ! 定义势能 
       ! 请给定给定任意初始势能深度，注意V0为负数
       V0=  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! 请给定给定任意初始势能深度，注意V0为负数
       
       do ir=0, irmatch
         r=ir*hcm 
         vpot1(ir) = gausspot(r,v0,0.0d0,1.484d0)  ! 注意此处没有加库伦力
         write(8,*) r, vpot1(ir)
       end do 
       
       call flush(8)
       
       ! 计算索末菲参数
       ecm=be
       mu=amu*(mass2*mass1)/(mass1+mass2)
       knp=sqrt(-2*mu*ecm/(hbarc**2))
       eta=z2*z1*e2*mu/hbarc/hbarc/knp

       ! 定义相遇点的位置
       irmid=nint(2.0d0/hcm) ! arbitrary radius

       ! 初始化delta的数值
       delta=0.0d0
       n=1
       l=0 ! s-wave 
       
       do !寻找合适的势阱深度
        n=n*(1+delta)
        vpot(0:irmatch)=n * vpot1(0:irmatch)
c****
c      Numerov method to solve the differential radial equation
c      from zero
       phi(0)=0      ! boundary condition
       phi(1)=hcm    ! arbitrary value

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccc请完善Numerov方法cccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
       fl=phi(irmid)
! 求解phi(irmid) 点出的微分， 需要phi(irmid-2)， phi(irmid-1)， phi(irmid+1)， phi(irmid+2) 的数值
       fpin=(-phi(irmid+2)+8.*phi(irmid+1)-8.*phi(irmid-1)
     &                    +phi(irmid-2))/12./hcm
       ! 储存numerov的结果从0到irmid区间段
       fin(0:irmid)=phi(0:irmid)



c      from infinity
        IE=0
        call WHIT(eta,hcm*(irmatch),knp,ecm,l,WK,WKD,IE)
        phi(irmatch)=WK(l)

        call WHIT(eta,hcm*(irmatch-1),knp,ecm,l,WK,WKD,IE)
        phi(irmatch-1)=WK(l)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccc请完善Numerov方法cccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! 求解phi(irmid) 点出的微分， 需要phi(irmid-2)， phi(irmid-1)， phi(irmid+1)， phi(irmid+2) 的数值
       	fpout=(-phi(irmid+2)+8.*phi(irmid+1)-8.*phi(irmid-1)
     &       	+phi(irmid-2))/12./hcm

        ! 计算0到irmid 与 irmid到无穷远 的系数
        norm=phi(irmid)/fl
        
        !前半段波函数乘以上面的系数 
        phi(0:irmid)=fin(0:irmid)*norm
        fpin=fpin*norm

        ! 计算\int phi(r) V(r) phi(r) dr 
        phivphi=0.0d0
        do i=1,nr
         f0=FFR4(rr(i)/hcm,vpot,irmatch+1)*
     &              abs(FFR4(rr(i)/hcm,phi(0:irmatch),irmatch+1))**2
         phivphi=phivphi+f0*rrw(i)
        end do ! 
        
        ! 计算\delta的数值
        delta=real(phi(irmid)*(fpout-fpin)/phivphi)
        
        ! 当\delta的数值足够小的时候我们可以认为该势能支持束缚态的存在
        if (abs(delta)<1e-6) exit

      end do ! for potential



c     归一化波函数

      norm=0.0d0
      do i=1,nr
      norm=norm+abs(FFR4(rr(i)/hcm,phi(0:irmatch),irmatch+1))**2*rrw(i)
      end do
      norm=1.0d0/norm
      phi=phi*sqrt(norm)
c****
      write(*,99)V0,n*V0,n
99    format(2x,"Adjust potential depth from ",F8.3," to ",F8.3,
     +                  " with scaling factor " F7.3)


        !输出波函数
       	do ir=0,irmatch
       		write (7,*) hcm*ir, phi(ir)
       	end do
        write(7,*)"&"

       write(*,*) ">>>>>!!!!!!!"

      end subroutine npbound
c-----------------------------------------------------------------------



c     interpolation function for uniform grids 
      FUNCTION FFR4(Y,F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 F(N),P,P1,P2,Q,X,FFR4
      REAL*8 Y
      PARAMETER(X=.16666666666667)
      P=Y
      I=P
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=P-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFR4=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFR4=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFR4=F(N)
      RETURN
      END function



C>     this subroutine calculates standard gauss Legendre points
C>     between x1 and x2 (usually -1.0d0 and  1.0d0)
C>     N is the number of mesh points required.
C>     The grid and the weights are stored in the arrays X and W
C>     @param[in] x1 lower boundary
C>     @param[in] x2 upper boundary
C>     @param[in] N number of grid points
C>     @param[out] X grid points
C>     @param[out] W integration weights
       SUBROUTINE gauleg(N,x1,x2,X,W)
        IMPLICIT NONE
        INTEGER N
        REAL*8 x1,x2,X(N),W(N)
        REAL*8 z1,z,xm,xl,pp,p3,p2,p1,pi,tol
        INTEGER m,i,j

        pi=acos(-1.0)
        tol=1.E-12

        m=(n+1)/2
        xm=0.5*(x2+x1)
        xl=0.5*(x2-x1)

        DO 10 i=1,m
         z=cos(pi*(i-0.25)/(N+0.5))

 20      CONTINUE
         p1=1.0E0
         p2=0.0E0
         DO 30 j=1,N
          p3=p2
          p2=p1
          p1=((2*j-1)*z*p2-(j-1)*p3)/j
 30      CONTINUE
         pp=N*(z*p1-p2)/(z*z-1.0E0)
         z1=z
         z=z1-p1/pp
         IF( abs(z1-z) .GT. tol) GOTO 20 ! Scheifenende

         X(i) = xm - xl*z
         X(n+1-i) = xm + xl*z
         W(i) = 2.E0*xl/((1.0-z*z)*pp*pp)
         W(n+1-i) = W(i)
 10     CONTINUE
       END SUBROUTINE gauleg



      SUBROUTINE WHIT(HETA,R,XK,E,LL,F,FD,IE)
C
C     CALCULATES  WHITTAKER  FUNCTION  WL(K,R)  WITH
C     ASYMPTOTIC  FORM  EXP(-(KR + ETA(LOG(2KR)))
C     E  IS  NEGATIVE
C     If IE = 0, allowed to return result e**IE larger than Whittaker,
C                for the IE value returned.
C     If IE > 0, must scale results by that amount.
C
C   input : 
C           HETA : Sommerfeld parameter
C           R : radius 
C           XK: module of wavenumber in fm^{-1}
c           E :  C.M. energy in MeV 
c           LL :  partial wave 
C           IE :  normally set to 0 
c   output:
c           F(LL+1) : WHITTAKER  FUNCTION
C           FD(LL+1) : derivative WHITTAKER  FUNCTION


!	use drier !  AMoro
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(LL+1),FD(LL+1) ,T(12),S(7)
!! AMoro: to replace drier module
      REAL*8 FPMAX
!	acc8 = epsilon(acc8);
      fpmax = huge(acc8)**0.8d0 
!! ------------------------------

      L = LL+1
C              NOW L = NO. OF VALUES TO FIND
      EE=-1.0
      AK=XK
      ETA=HETA
      LP1=L+1
      RHO=AK*R
	S(:) = 0
      IF(L-50)1,1,2
    1 LM=60
      GO TO 3
    2 LM=L+10
    3 LMP1=LM+1
      IS=7
      PJE=30.0*RHO+1.0
      H=max(INT(PJE),4)
      H=RHO/H
      RHOA=10.0*(ETA+1.0)
      IF(RHOA-RHO)13,13,14
   13 IFEQL=1
      RHOA=RHO
      GO TO 15
   14 IFEQL=0
   15 PJE=RHOA/H+0.5
      RHOA=H*INT(PJE)
      IF(IFEQL)16,16,18
   16 IF(RHOA-RHO-1.5*H)17,18,18
   17 RHOA=RHO+2.0*H
   18 IF(EE)55,55,19
   19 STOP 'WHIT'
   27 A=2.0-10.0/12.0*H*H*EE
      B=1.0/6.0*H*ETA
      C=1.0+1.0/12.0*H*H*EE
      M1=INT(RHOA/H-0.5)
      M2=INT(RHO/H-1.5)
      T(2)=B/FLOAT(M1+1)
      T(3)=B/FLOAT(M1)
      JS=M1
      DO 29 IS=M2,M1
      DO 28 I=1,6
      S(I)=S(I+1)
   28 CONTINUE
      T(1)=T(2)
      T(2)=T(3)
      T(3)=B/FLOAT(JS-1)
      S(7)=((A+10.0*T(2))*S(6)-(C-T(1))*S(5))/(C-T(3))
      JS=JS-1
      IF(ABS(S(7)).LE.FPMAX) GO TO 29
       DO 285 I=2,7
  285   S(I) = S(I) / FPMAX
   29 CONTINUE
      T(1)=S(4)
      T(2)=(1.0/60.0*(S(1)-S(7))+0.15*(S(6)-S(2))+0.75*(S(3)-S(5)))/H
      GO TO 60
   55 C=1.0/RHOA
      A=1.0
      B=1.0-C*ETA
      F(1)=A
      FD(1)=B
      DO 56 M=1,26
      D=0.5*(ETA+FLOAT(M-1))*(ETA+FLOAT(M))*C/FLOAT(M)
      A=-A*D
      B=-B*D-A*C
      F(1)=F(1)+A
      FD(1)=FD(1)+B
   56 CONTINUE
      A=-ETA*LOG(2.0*RHOA)-RHOA
      FPMINL = -LOG(FPMAX)
      if(IE.eq.0.and.A.LT.FPMINL) IE = INT(FPMINL-A)
      A=EXP(A+IE)
      F(1)=A*F(1)
c      FD(1)=A*FD(1)
      FD(1)=A*FD(1) * (-1d0 - 2*ETA/(RHOA))
      IF(IFEQL)57,57,61
   57 S(IS)=F(1)
      IF(IS-7)27,58,27
   58 IS=6
      RHOA=RHOA+H
      GO TO 55
   60 F(1)=T(1)
      FD(1)=T(2)
   61 C=1.0/RHO
      DO 63 M=1,L-1
      A=ETA/FLOAT(M)
      B=A+C*FLOAT(M)
      F(M+1)=(B*F(M)-FD(M))/(A+1.0)
      FD(M+1)=(A-1.0)*F(M)-B*F(M+1)
   63 CONTINUE
      DO 65 M=1,L
      FD(M)=AK*FD(M)
   65 CONTINUE
      RETURN
      END SUBROUTINE
      
      
c *** Gaussian
      function gausspot(r,v0,r0,a)
      implicit none
       real*8 r,v0,r0,gausspot,a
         if (a.gt.1e-6) then
           gausspot=V0*exp(-(r-r0)**2/a**2)
             else
               write(*,*)'a too small in gausspot!'
               stop
         endif
         return
      end function
c Coulomb potential
      FUNCTION VCOUL(R,z12,Rc)
          use constants
          implicit none
          real*8 r,rc,rc2,aux,vcoul,z12

          RC2=RC*2d0
          aux=e2*Z12
          vcoul=0
          if (z12.lt.1e-4) return
          if (rc.lt.1e-6) rc=1e-6

          IF(R.GT.RC)GO TO 1
          VCOUL=AUX*(3.-(R/RC)**2)/RC2
          RETURN
1         VCOUL=AUX/R
          RETURN
        END function
