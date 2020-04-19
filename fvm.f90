MODULE fvm
    USE pardiso_sparse_solver_v3
    IMPLICIT NONE
    SAVE
    REAL(8):: p
    INTEGER(8):: maxiters,d,point1,point2,pvaries
    INTEGER(8), dimension(:), allocatable:: DAT3,DAT4
CONTAINS

    SUBROUTINE EV_F_CHIP(N, X, NEW_X, F, DAT, IDAT, IERR)
        IMPLICIT NONE
        INTEGER(8) N, NEW_X
        REAL(8) F, X(N)
        REAL(8) DAT(3*d**2-2*d,2), IDAT(d+1,d+1)
        INTEGER(8) IERR
        
        IF (NEW_X) THEN
            CALL fvm_simulate(X,d,DAT,IDAT)
        ENDIF

        F = 0.5*dot_product(DAT(:,2),DAT(:,2))/(1d0*INT(N))
        IERR = 0
        RETURN
    END

    SUBROUTINE fvm_simulate(v, d,DAT,IDAT)
        INTEGER(8), intent(in)::d
        INTEGER(8)	:: l,b,a,i, j, ii, jj, m,tot,el ! number of design variables
        REAL(8)     ::  lengthx, lengthy, lengthz, Q, Ts
        REAL(8)     ::  n, dx, dy, hx, hy
        REAL(8)     ::  bdir,  Sdir
        REAL(8)     :: ke, kw, kn, ks !c
        REAL(8)		:: dx_east, dx_west, dy_north, dy_south
        REAL(8), DIMENSION(5)  :: xi
        REAL(8), DIMENSION((d-1)**2), INTENT(in)   :: v
        REAL(8), DIMENSION(d*d) :: f
        REAL(8), DIMENSION(d+1,d+1), INTENT(inout) :: IDAT
        REAL(8), DIMENSION(3*d*d-2*d ,2),INTENT(inout) :: DAT
        
        ! Assign variables
        lengthx = 0.005d0
        lengthy = 0.005d0
        lengthz = 0.001d0
        Q = 0.5d0/(lengthz*lengthx*lengthy)
        Ts = 293d0
        DAT3(2)=49532 
        a = d
        b = a
        m = (a-1)**2
        n = a*b !number of state grid cells.
        dx = lengthx/(a-1) 
        dy = lengthy/(b-1)  
        DAT=0d0 

        ! Create discretized conductivity (design) field
        DO ii = 2,d
            DO jj = 2,d
                IDAT(ii,jj) = 0.2d0+(65d0-0.2d0)*v(ii-1+(a-1)*(jj-2))**p
            ENDDO
        ENDDO

        IDAT(1,2:d)= IDAT(2,2:d)
        IDAT(d+1,2:d)=IDAT(d,2:d)
        IDAT(:,1)=IDAT(:,2)
        IDAT(:,d+1)=IDAT(:,d)

        f = 0
        tot =1

        DO l = 1,n
            dat3(l)=tot
            el=0
            i= 1+ mod(l-1,a)
            j =  ceiling(1d0*l/b)
    
            ! Check cell dimensions
            IF ( i==a .OR. i==1)  THEN
             hx = dx/2
            ELSE
            hx = dx
            ENDIF

            IF (j==1 .OR. j==b) THEN
            hy = dy/2
            ELSE
            hy = dy
            ENDIF

            ! Initialize xi & interpolate .OR.conductivities of each side
            kw =  2* (IDAT(i,j)**(-1)+IDAT(i,j+1)**(-1))**(-1)
            ke =  2*(1/IDAT(i+1,j)+1/IDAT(i+1,j+1))**(-1)
            kn = 2*(1/IDAT(i,j+1) + 1/IDAT(i+1,j+1))**(-1)
            ks = 2*(1/IDAT(i,j) + 1/IDAT(i+1,j))**(-1)
            xi = 0
            bdir = 0
            Sdir = 0
            
            ! Start building stiffness matrix
            IF (.NOT. (i==1)) THEN
                IF (i==a .OR. i==2) THEN
                    dx_west = 0.75d0*dx
                ELSE
                    dx_west = dx
                ENDIF
                xi(1) =  - (kw*hy)/dx_west
            ENDIF

            IF (.NOT. (i==a)) THEN
                IF (i==a-1 .OR. i==1) THEN
                    dx_east = 0.75d0*dx
                ELSE
                    dx_east = dx
                ENDIF
                xi(2) =  -(ke*hy)/dx_east
                tot=tot+1
                el=el+1
                DAT(DAT3(l)+el,1) = xi(2)
                DAT4(DAT3(l)+el)=l+1
            ENDIF

            IF (.NOT. (j==1)) THEN
                IF (j==2 .OR. j==b) THEN
                    dy_south = 0.75d0*dy
                ELSE
                    dy_south = dy
                ENDIF
                xi(3) =  - (ks*hx)/dy_south
            ENDIF

            IF (.NOT. (j==b)) THEN
                IF ((j==b-1) .OR. (j==1)) THEN
                    dy_north = 0.75d0*dy
                ELSE
                    dy_north = dy
                ENDIF
                xi(4) = -(kn*hx)/dy_north
                el= el+1
                tot=tot+1
                DAT(DAT3(l)+el,1) = xi(4)
                DAT4(DAT3(l)+el)=l+a
            ENDIF

            ! Define source terms .OR.Dirichlet Boundary
            IF (( ((j-1)*dy >= 0.003d0) .and. ((j-1)*dy <= 0.005d0))) THEN
                IF (i == 1) THEN
                    bdir = 2* Ts*kw*hy/hx
                    Sdir = 2*kw*hy/hx
                ENDIF
            ENDIF


            tot = tot+1
            DAT(DAT3(l),1)= - (sum(xi) ) + Sdir
            f(l) = Q*hx*hy + bdir
            DAT4(DAT3(l))=l
        ENDDO
        DAT3(d**2+1)=tot
        CALL pardiso_sym_solver2(DAT(:,1),DAT3,DAT4, f, d, DAT(1:d**2,2))
    END SUBROUTINE



    SUBROUTINE harmonic_drv(x,c,dharm)
        REAL(8), intent(out) :: dharm
        REAL(8), intent(in) :: x,c
        dharm= 2*(c**2/(x+c)**2)
    END SUBROUTINE



    SUBROUTINE getdir(d,ii,jj,dy,dx,Ts,dharm,theta_k_v,dbdir,dSdir)
        INTEGER(8), intent(in)  :: d
        INTEGER(8)			    :: ii, jj
        REAL(8)                 :: dx, dy, Ts, dbdir, dSdir, dharm
        REAL(8), DIMENSION(d+1,d+1), intent(inout) :: theta_k_v
        dbdir =  4*Ts*dharm*theta_k_v(ii+1,jj+1)*dy/dx
        dSdir = 4*dharm*theta_k_v(ii+1,jj+1)*dy/dx
    END SUBROUTINE



    SUBROUTINE EV_GRAD_F_CHIP(N, X, NEW_X, GRAD, DAT, IDAT, IERR)
        IMPLICIT NONE
        INTEGER(8) NEW_X, IERR
        INTEGER N
        REAL(8)  X(N), GRAD(N)
        REAL(8) DAT( 3*d**2-2*d, 2),IDAT(d+1,d+1)
     
        IF (NEW_X) THEN
            CALL fvm_simulate(X,d,DAT,IDAT)
        ENDIF
            CALL ADJOINT(X,d,GRAD,DAT,IDAT)
            GRAD=GRAD/(1d0*N)
            IERR = 0
        RETURN
    END SUBROUTINE


    SUBROUTINE ADJOINT(v,d,GRAD,DAT,IDAT)
        IMPLICIT NONE
        INTEGER(8) :: d,n
        REAL(8):: GRAD((d-1)**2), v((d-1)**2)
        REAL(8):: DAT(3*d*d-2*d,2)
        REAL(8):: IDAT(d+1,d+1)
        INTEGER(8)   :: l,b,a,i, j, ii, jj, m !number of design variables
        REAL(8)      :: lengthx, lengthy, lengthz, Q, Ts
        REAL(8)      ::   dx, dy, hx, hy, dkdk
        REAL(8)      ::   dbdir, dSdir, dharm
        REAL(8)		 :: dx_east, dx_west, dy_north, dy_south
        REAL(8), DIMENSION(3)  :: dxi_dv, vec
        REAL(8), DIMENSION(d*d) :: PSI, theta_R_v
        REAL(8), DIMENSION(d+1,d+1) :: theta_k_v
        REAL(8), DIMENSION(d*d,5) :: dxi_dk
               
        CALL pardiso_sym_solver3(DAT(:,1),DAT3,DAT4, -DAT(1:d**2,2), d, PSI)

        lengthx = 0.005d0
        lengthy = 0.005d0
        lengthz = 0.001d0
        Q = 0.5d0/(lengthz*lengthx*lengthy)
        Ts = 293d0

        a = d
        b = a
        m = (a-1)**2
        n = a*b !number of state grid cells.
        dx = lengthx/(a-1)
        dy = lengthy/(b-1)  !dy should p

        theta_k_v=0
        dxi_dk=0
        DO ii = 2,d
            DO jj = 2,d
                theta_k_v(ii,jj)= p*(65d0-0.2d0)*v(ii-1+(a-1)*(jj-2))**(p-1d0)
            ENDDO
        ENDDO

        theta_k_v(1,2:d)= theta_k_v(2,2:d)
        theta_k_v(d+1,2:d)=theta_k_v(d,2:d)
        theta_k_v(:,1)=theta_k_v(:,2)
        theta_k_v(:,d+1)=theta_k_v(:,d)

        DO l = 1,d**2
            i= 1+ mod(l-1,a)
            j =  ceiling(1d0*l/b)
    
            IF ( i==a .OR. i==1)  THEN
                hx = dx/2
            ELSE
                hx = dx
            ENDIF
                
            IF (j==1 .OR. j==b) THEN
                hy = dy/2
            ELSE
                hy = dy
            ENDIF
   
            IF (.NOT. (i==1)) THEN
                IF (i==a .OR. i==2) THEN
                    dx_west = 0.75d0*dx
                ELSE
                    dx_west = dx
                ENDIF
                dxi_dk(l,2) = -hy/dx_west
            ENDIF

            IF (.NOT. (i==a)) THEN
                IF (i==a-1 .OR. i==1) THEN
                    dx_east = 0.75d0*dx
                ELSE
                    dx_east = dx
                ENDIF
                dxi_dk(l,4) =  -hy/dx_east
            ENDIF

            IF (.NOT. (j==1)) THEN
                IF (j==2 .OR. j==b) THEN
                    dy_south = 0.75d0*dy
                ELSE
                    dy_south = dy
                ENDIF
                dxi_dk(l,1) = -hx/dy_south
            ENDIF

            IF (.NOT. (j==b)) THEN
                IF ((j==b-1) .OR. (j==1)) THEN
                    dy_north = 0.75d0*dy
                ELSE
                    dy_north = dy
                ENDIF
                dxi_dk(l,5) = -hx/dy_north
            ENDIF
            dxi_dk(l,3) = dxi_dk(l,1)+dxi_dk(l,2)+dxi_dk(l,4)+dxi_dk(l,5)
        ENDDO


        DO l = 1, (d-1)**2
            theta_R_v = 0
            ii  =  1+ mod(l-1,a-1)
            jj =  ceiling((1.d0*l)/(b-1))
            i =  1+ mod(l-1,a)
            j =  ceiling((1.d0*l)/(b))
            m = ii+(a)*(jj-1)

            dxi_dv = 0d0

            !! 1st ELEMENT
            dbdir = 0d0
            dSdir = 0d0

            IF (jj==1) THEN
                dkdk = 1d0
            ELSE
                CALL harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj),dkdk)
            ENDIF

            dxi_dv(1) = dxi_dk(m,4)*dkdk*theta_k_v(ii+1,jj+1)!east

            IF (ii==1) THEN
                dkdk = 1
                IF (( ((jj-1)*dy >= 0.003d0)  .AND. ((jj-1)*dy <= 0.005d0))) THEN
                    CALL harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj), dharm)
                    CALL getdir(d,ii,jj,dy,dx,Ts,dharm,theta_k_v,dbdir,dSdir)
                ELSE
                    dbdir = 0d0
                    dSdir = 0d0
                ENDIF
            ELSE
                CALL harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii,jj+1),dkdk)
            ENDIF

            dxi_dv(2) = dxi_dk(m,5)*dkdk*theta_k_v(ii+1,jj+1) ! North
            dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2)) + dSdir

            vec = (/DAT(m+1,2),DAT(m+a,2),DAT(m,2)/)
            theta_R_v(m) =  dot_product(dxi_dv, vec) - dbdir
             
            !! 2nd ELEMENT
            dbdir = 0d0
            dSdir = 0d0

            IF (jj==1) THEN
                dkdk = 1d0
            ELSE
                CALL harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj),dkdk)
            ENDIF
    
            dxi_dv(1) = dxi_dk(m+1,2)*dkdk*theta_k_v(ii+1,jj+1) ! West

            IF (ii==a-1) THEN
                dkdk = 1d0
                dbdir =  0d0
                dSdir = 0d0
            ELSE
                CALL harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+2,jj+1),dkdk)
                dbdir =  0d0
                dSdir = 0d0
            ENDIF
            dxi_dv(2) = dxi_dk(m+1,5)*dkdk*theta_k_v(ii+1,jj+1) ! North
            dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2)) + dSdir
            vec = (/DAT(m,2),DAT(m+1+a,2),DAT(m+1,2)/)
            theta_R_v(m+1) =  dot_product(dxi_dv, vec) - dbdir
           
            dbdir =  0d0
            dSdir = 0d0
            !! 3rd ELEMENT
            IF (jj==b-1) THEN
                dkdk = 1d0
            ELSE
                CALL harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj+2),dkdk)
            ENDIF

            dxi_dv(1) = dxi_dk(m+a,4)*dkdk*theta_k_v(ii+1,jj+1)!east
            IF (ii==1) THEN
                dkdk = 1d0
                IF (( ((jj)*dy >= 0.003d0)  .and. ((jj)*dy <= 0.005d0)))  THEN
                    CALL harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj+2), dharm)
                    CALL getdir(d,ii,jj,dy,dx,Ts,dharm,theta_k_v,dbdir,dSdir)
                ELSE
                    dbdir =  0d0
                    dSdir = 0d0
                ENDIF
            ELSE
                CALL harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii,jj+1),dkdk)
                dbdir = 0d0
                dSdir = 0d0
            ENDIF
                dxi_dv(2) = dxi_dk(m+a,1)*dkdk*theta_k_v(ii+1,jj+1)!south
                dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2))+dSdir
                vec = (/DAT(m+a+1,2),DAT(m,2),DAT(m+a,2)/)
                theta_R_v(m+a) =  dot_product(dxi_dv, vec) - dbdir
            
            !! 4th ELEMENT
            dbdir =  0d0
            dSdir = 0d0
            IF (jj==b-1) THEN
                dkdk = 1d0
            ELSE
                CALL harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj+2),dkdk)
            ENDIF
    
            dxi_dv(1) = dxi_dk(m+a+1,2)*dkdk*theta_k_v(ii+1,jj+1) ! West

            IF (ii==a-1) THEN
                dkdk = 1d0
                dbdir =  0d0
                dSdir = 0d0
            ELSE
                CALL harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+2,jj+1),dkdk)
                dbdir = 0d0
                dSdir = 0d0
            ENDIF
            dxi_dv(2) = dxi_dk(m+a+1,1)*dkdk*theta_k_v(ii+1,jj+1)!south
            dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2)) + dSdir
            vec = (/DAT(m+a,2),DAT(m+1,2),DAT(m+a+1,2)/)
            theta_R_v(m+a+1) =  dot_product(dxi_dv, vec) - dbdir
            GRAD(l) = dot_product(PSI,theta_R_v)
        ENDDO
    END SUBROUTINE


    SUBROUTINE EV_G(N, X, NEW_X, M, G, DAT, IDAT, IERR)
        IMPLICIT NONE
        INTEGER N, NEW_X, M
        DOUBLE PRECISION G(M), X(N), ONES(N), lengthx
        REAL(8) DAT( 3*d**2-2*d, 2),IDAT(d+1,d+1)
        INTEGER(8) IERR

        lengthx=0.005d0
        IF (NEW_X) THEN
            CALL fvm_simulate(X,d,DAT,IDAT)
        ENDIF

        ONES=(lengthx)**2/(1d0*N)
        G(1)=dot_product(ONES,X)
        IERR = 0
        RETURN
    END SUBROUTINE

    SUBROUTINE EV_JAC_G(TASK, N, X,NEW_X, M, NZ, ACON, AVAR, A,DAT,IDAT, IERR)
        INTEGER TASK, N,  M, NZ,IERR,NEW_X
        REAL(8) DAT( 3*d**2-2*d, 2),IDAT(d+1,d+1)
        DOUBLE PRECISION  A(NZ),X(N)
        INTEGER ACON(NZ), AVAR(NZ), I
        INTEGER AVAR1(N), ACON1(N)

      
        ! Structure of Jacobian:
        DO I=1,N
            AVAR1(I)=I
        ENDDO
        ACON1= 1

        IF( TASK.eq.0 ) THEN
            DO I = 1, N
                AVAR(I) = AVAR1(I)
                ACON(I) = ACON1(I)
            ENDDO
        ELSE
            DO I=1,N
                A(I) = (0.01d0)**2/N
            ENDDO
        ENDIF
        IERR = 0
        RETURN
    END SUBROUTINE
	
	
	  subroutine ITER_CB(ALG_MODE, ITER_COUNT,OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT, DAT, ISTOP)
      implicit none
      integer ALG_MODE, ITER_COUNT, LS_TRIAL
      double precision OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
      double precision ALPHA_DU, ALPHA_PR
      
       REAL(8) DAT( 3*d**2-2*d, 2),IDAT(d+1,d+1)
      integer ISTOP
    
	 
	  if (ITER_COUNT .eq. 1) then
	  write(*,*) 'Test with p varying from 1 to 3'
	  endif
	  if (ITER_COUNT .lt. point1) p = 1
	  
	  if (ITER_COUNT .eq. point1) then
	  write(*,*) 'p increase begins'
	  endif
	  if ((ITER_COUNT .gt. point1-1 ) .AND. (ITER_COUNT .lt. point2)) then
				p=1 + 2*(ITER_COUNT-point1)/(point2-point1)
	endif
	
	  if (ITER_COUNT .eq. point2) then
	  write(*,*) 'p increase finished'
	  endif
                 if (ITER_COUNT .gt. point2-1) p = 3
	
	  !write(*,*) 'p=', p,'iter =', ITER_COUNT
	  ! example stopping condition
     ! if (INF_PR.le.1D-04) ISTOP = 1
	 
	 
	 if (pvaries .eq. 0) then
	   p =3
	   endif
	   

      return
      end
	
	
	
	
	
	
END MODULE
