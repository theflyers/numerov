
      program test
      implicit none
      ! 定义格点信息
      integer :: nr ! 高斯积分需要用的变量
      integer :: i
      real*8,allocatable,dimension(:) :: rr,rrw ! 高斯积分需要用的变量
      real*8 :: f,f1
      nr=40  !高斯格点数
      allocate(rr(1:nr),rrw(1:nr))
      call gauleg(nr,0.0d0,40.0d0,rr,rrw) ! 初始化高斯积分，输出权重与高斯点
      
      !计算从0到40的x^2的积分
      f=0.0d0
      do i=1,nr
          f=rr(i)**2*rrw(i)+f         !使用输出的格点与权重进行积分
      enddo
      f1=40.d0**3/3.d0             !计算理论值
      end program
      
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
