
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
