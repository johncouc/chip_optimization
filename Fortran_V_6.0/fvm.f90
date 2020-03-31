module fvm

     ! USE sparse_system_solvers
      USE pardiso_sparse_solver  	

IMPLICIT NONE
SAVE
real(8):: p
integer(8):: maxiters,d
integer(8), dimension(:), allocatable:: DAT3,DAT4  
CONTAINS




  subroutine EV_F_CHIP(N, X, NEW_X, F, DAT, IDAT, IERR)
      implicit none
      integer(8) N, NEW_X
       
      REAL(8) F, X(N)
      !REAL(8) DAT( (floor(sqrt(N*1d0))+1)**2, (floor(sqrt(N*1d0))+1)**2+1),IDAT(floor(sqrt(N*1d0)) + 2,(floor(sqrt(N*1d0)) + 2))
      REAL(8) DAT(3*d**2-2*d,2), IDAT(d+1,d+1)
      !REAL(8) DAT(*),IDAT(*)  
      integer(8) IERR
      !integer(8) d
      !d = floor(sqrt(N*1d0))+1
      if (NEW_X) then
        call fvm_simulate(X,d,DAT,IDAT)
              end if
!      write(*,*) "T = ", DAT(:,2)

      F = 0.5*dot_product(DAT(:,2),DAT(:,2))/(1d0*INT(N))
!        F=dot_product(X,X)
      IERR = 0
      return
      end

  SUBROUTINE fvm_simulate(v, d,DAT,IDAT)
        integer(8), intent(in)::d
        integer(8)				       :: l,b,a,i, j, ii, jj, m,tot,el !number of design variables
	REAL(8)                                   ::  lengthx, lengthy, lengthz, Q, Ts
	REAL(8)                                   ::  n, dx, dy, hx, hy 
	REAL(8)                                   ::  bdir,  Sdir
	REAL(8)                                   :: ke, kw, kn, ks !c
	REAL(8)				       :: dx_east, dx_west, dy_north, dy_south
	REAL(8), DIMENSION(5)                     :: xi
	REAL(8), DIMENSION((d-1)**2), INTENT(in)      :: v
	REAL(8), DIMENSION(d*d) :: f 
        REAL(8), DIMENSION(d+1,d+1), INTENT(inout) :: IDAT
        REAL(8), DIMENSION(3*d*d-2*d ,2),INTENT(inout) :: DAT
        
        ! assign variables
        lengthx = 0.005d0
        lengthy = 0.005d0
        lengthz = 0.001d0
        Q = 0.5d0/(lengthz*lengthx*lengthy)
        Ts = 293d0
        DAT3(2)=49532 
        a = d
        b = a
        m = (a-1)**2
        n = a*b      !number of state grid cells. 
        dx = lengthx/(a-1) 
        dy = lengthy/(b-1)  
        DAT=0d0 
       ! IDAT=0d0        

    ! create discretized conductivity (design) field 
       do ii = 2,d
          do jj = 2,d
          IDAT(ii,jj) = 0.2d0+(65d0-0.2d0)*v(ii-1+(a-1)*(jj-2))**p  
        enddo
    enddo


         IDAT(1,2:d)= IDAT(2,2:d)
         IDAT(d+1,2:d)=IDAT(d,2:d)
         IDAT(:,1)=IDAT(:,2)
         IDAT(:,d+1)=IDAT(:,d)        


        f = 0
       tot =1

        do l = 1,n
           dat3(l)=tot
           el=0     
           i= 1+ mod(l-1,a) 
           j =  ceiling(1d0*l/b) 
    
! check cell dimensions 

        if ( i==a .or. i==1)  then
         hx = dx/2 
        else 
        hx = dx 
        endif

        if (j==1 .or. j==b) then
        hy = dy/2 
        else
        hy = dy 
        endif

! initialize xi & interpolate for conductivities of each side

    kw =  2* (IDAT(i,j)**(-1)+IDAT(i,j+1)**(-1))**(-1) 
    ke =  2*(1/IDAT(i+1,j)+1/IDAT(i+1,j+1))**(-1) 
    kn = 2*(1/IDAT(i,j+1) + 1/IDAT(i+1,j+1))**(-1) 
    ks = 2*(1/IDAT(i,j) + 1/IDAT(i+1,j))**(-1)   
    xi = 0
    bdir = 0
    Sdir = 0
! start building stiffness matrix       

        if (.not. (i==1)) then
         if (i==a .or. i==2) then
                dx_west = 0.75d0*dx 
         else
                 dx_west = dx 
        end if
        xi(1) =  - (kw*hy)/dx_west 
!        DAT(l,l-1) =  xi(1) 
     
        end if

        if (.not. (i==a)) then
                if (i==a-1 .or. i==1) then
                  dx_east = 0.75d0*dx 
                else
                 dx_east = dx 
                endif
        xi(2) =  -(ke*hy)/dx_east 
 !       DAT(l,l+1) =  xi(2) 
        tot=tot+1
        el=el+1
        DAT(DAT3(l)+el,1) = xi(2)
        DAT4(DAT3(l)+el)=l+1
       



        end if

        if (.not. (j==1)) then
                if (j==2 .or. j==b) then
                 dy_south = 0.75d0*dy 
        else
                dy_south = dy 
    end if
        xi(3) =  - (ks*hx)/dy_south 
!        DAT(l,l-a) =  xi(3) 
    !dxi_dk(l,l-a) = -hx/dy_south 
        endif

    if (.not. (j==b)) then
        if ((j==b-1) .or. (j==1)) then
                dy_north = 0.75d0*dy 
    else
    dy_north = dy 
    end if
         xi(4) = -(kn*hx)/dy_north 
 !        DAT(l,l+a) = xi(4) 
        el= el+1
        tot=tot+1
        DAT(DAT3(l)+el,1) = xi(4)
       DAT4(DAT3(l)+el)=l+a
       

        end if

! define source terms for Dirichlet Boundary

        if (( ((j-1)*dy >= 0.003d0) .and. ((j-1)*dy <= 0.005d0))) then
                if (i == 1) then
 
                bdir = 2* Ts*kw*hy/hx
                Sdir = 2*kw*hy/hx 
        end if
    !    if (i== a) then
     
     !           bdir =  2*Ts*ke*hy/hx
     !           Sdir = 2*ke*hy/hx 
     !   end if
        endif

    !    DAT(l,l) =  -(sum(xi))+(Sdir)

       tot = tot+1
       DAT(DAT3(l),1)= - (sum(xi) ) + Sdir
        f(l) = Q*hx*hy + bdir
       DAT4(DAT3(l))=l
        enddo
         DAT3(d**2+1)=tot
!        CALL solver_v1(DAT(1:(d)**2,1:(d)**2), f, d*d, DAT(:,(d)**2+1))
       ! CALL solver_v2(DAT(1:(d)**2,1:(d)**2), f, d*d, DAT(:,(d)**2+1))
        call pardiso_sym_solver2(DAT(:,1),DAT3,DAT4, f, d, DAT(1:d**2,2)) 
       ! open (unit = 4, file = "temp")
       ! do l=1,N
       ! write(4,*)DAT(l,d**2+1)
       ! enddo          
       ! close(4)
        end subroutine

        subroutine harmonic_drv(x,c,dharm) 
        real(8), intent(out) :: dharm
        real(8), intent(in) :: x,c
        dharm= 2*(c**2/(x+c)**2)
        end subroutine 



        subroutine getdir(d,ii,jj,dy,dx,Ts,dharm,theta_k_v,dbdir,dSdir)
        integer(8)    , intent(in)             :: d
	integer(8)			:: ii, jj
	REAL(8)                  :: dx, dy, Ts, dbdir, dSdir, dharm
	REAL(8), DIMENSION(d+1,d+1), intent(inout) :: theta_k_v

            dbdir =  4*Ts*dharm*theta_k_v(ii+1,jj+1)*dy/dx
            dSdir = 4*dharm*theta_k_v(ii+1,jj+1)*dy/dx  
        end subroutine



    subroutine EV_GRAD_F_CHIP(N, X, NEW_X, GRAD, DAT, IDAT, IERR)
      implicit none

        integer(8) N, NEW_X
      
        REAL(8)  X(N), GRAD(N)
        REAL(8) DAT( 3*d**2-2*d, 2),IDAT(d+1,d+1)

        integer(8) IERR
       ! integer(8) d
       ! d = floor(sqrt(N*1d0))+1
 
               ! CALLS = CALLS + 1
               ! if (CALLS.eq.maxiters/3) then
               ! write(*,*) CALLS,"evaluations, changing p : p=2"
               ! p=2d0
               ! call fvm_simulate(X,floor(sqrt(N*1d0) + 1),DAT,IDAT)
                
               ! endif
               ! if (CALLS.eq.maxiters/2) then
               ! write(*,*) CALLS,"evaluations, changing p : p=3"                       
               ! p=2.2d0
               ! call fvm_simulate(X,floor(sqrt(N*1d0) + 1),DAT,IDAT)
                
               ! endif

        if (NEW_X) then
        call fvm_simulate(X,d,DAT,IDAT)
        end if
        call ADJOINT(X,d,GRAD,DAT,IDAT)
        GRAD=GRAD/(1d0*N)
        IERR = 0
        return
        end subroutine


        subroutine ADJOINT(v,d,GRAD,DAT,IDAT)
              implicit none
              integer(8) :: d,n
              REAL(8):: GRAD((d-1)**2), v((d-1)**2)
              REAL(8):: DAT(3*d*d-2*d,2)
              REAL(8):: IDAT(d+1,d+1)
        
                integer(8)       :: l,b,a,i, j, ii, jj, m !number of design variables
        	REAL(8)                                   :: lengthx, lengthy, lengthz, Q, Ts
        	REAL(8)                                   ::   dx, dy, hx, hy, dkdk
        	REAL(8)                                   ::   dbdir, dSdir, dharm
!        	REAL(8)                                   :: ke, kw, kn, ks, dharm !,c
        	REAL(8)				       :: dx_east, dx_west, dy_north, dy_south
        	REAL(8), DIMENSION(3)                     :: dxi_dv, vec
!        	REAL(8), DIMENSION(5)                     :: xi
        	REAL(8), DIMENSION(d*d) :: PSI, theta_R_v
	        REAL(8), DIMENSION(d+1,d+1) :: theta_k_v 
        	REAL(8), DIMENSION(d*d,5) :: dxi_dk  
!FIELD ADJOINT EQ:
               !read(*,*) 
               
              ! CALL solver_v2(DAT(1:d*d,1:d*d),-DAT(1:d*d,d*d+1),d*d,PSI)
               call pardiso_sym_solver3(DAT(:,1),DAT3,DAT4, -DAT(1:d**2,2), d, PSI) 

                lengthx = 0.005d0
                lengthy = 0.005d0
                lengthz = 0.001d0
                Q = 0.5d0/(lengthz*lengthx*lengthy)
                Ts = 293d0
                
                a = d
                b = a
                m = (a-1)**2
                n = a*b      !number of state grid cells. 
                dx = lengthx/(a-1) 
                dy = lengthy/(b-1)  !dy should p

                theta_k_v=0
                dxi_dk=0
                do ii = 2,d
                        do jj = 2,d

                          theta_k_v(ii,jj)= p*(65d0-0.2d0)*v(ii-1+(a-1)*(jj-2))**(p-1d0) 
                        enddo
                enddo


            theta_k_v(1,2:d)= theta_k_v(2,2:d)
            theta_k_v(d+1,2:d)=theta_k_v(d,2:d)
            theta_k_v(:,1)=theta_k_v(:,2)
            theta_k_v(:,d+1)=theta_k_v(:,d)        




            do l = 1,d**2

 
                   i= 1+ mod(l-1,a) 
                   j =  ceiling(1d0*l/b) 
    
                if ( i==a .or. i==1)  then
                   hx = dx/2 
                else 
                   hx = dx 
                endif
                
                if (j==1 .or. j==b) then
                   hy = dy/2 
                else
                   hy = dy 
                endif
   
                if (.not. (i==1)) then
                    if (i==a .or. i==2) then
                    dx_west = 0.75d0*dx 
                    else
                    dx_west = dx 
                    end if
                dxi_dk(l,2) = -hy/dx_west 
   
     
                end if

                if (.not. (i==a)) then
                    if (i==a-1 .or. i==1) then
                    dx_east = 0.75d0*dx 
                            else
                    dx_east = dx 
                    endif
                    dxi_dk(l,4) =  -hy/dx_east 
                end if

                if (.not. (j==1)) then
                    if (j==2 .or. j==b) then
                    dy_south = 0.75d0*dy 
                    else
                    dy_south = dy 
                    end if
                    dxi_dk(l,1) = -hx/dy_south 
                endif

                if (.not. (j==b)) then
                    if ((j==b-1) .or. (j==1)) then
                    dy_north = 0.75d0*dy 
                    else
                    dy_north = dy 
                    end if
                 dxi_dk(l,5) = -hx/dy_north
                 end if

                dxi_dk(l,3) = dxi_dk(l,1)+dxi_dk(l,2)+dxi_dk(l,4)+dxi_dk(l,5)

                enddo


                do l = 1, (d-1)**2

 
                    theta_R_v = 0
                    ii  =  1+ mod(l-1,a-1)
                    jj =  ceiling((1.d0*l)/(b-1)) 
                    i =  1+ mod(l-1,a)
                    j =  ceiling((1.d0*l)/(b)) 
                    m = ii+(a)*(jj-1)

   
                dxi_dv = 0d0

! 1st ELEMENT

                dbdir = 0d0
                dSdir = 0d0

                if (jj==1) then
                    dkdk = 1d0
                else
                    call harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj),dkdk)
                endif

                dxi_dv(1) = dxi_dk(m,4)*dkdk*theta_k_v(ii+1,jj+1)!east

                if (ii==1) then
                    dkdk = 1
                    if (( ((jj-1)*dy >= 0.003d0)  .and. ((jj-1)*dy <= 0.005d0))) then
                	call harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj), dharm)
                    	call getdir(d,ii,jj,dy,dx,Ts,dharm,theta_k_v,dbdir,dSdir)
                    else
                 	dbdir = 0d0
                   	dSdir = 0d0
                    end if
                else
                    call harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii,jj+1),dkdk) 
                endif

                dxi_dv(2) = dxi_dk(m,5)*dkdk*theta_k_v(ii+1,jj+1)!north
                dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2)) + dSdir


                vec = (/DAT(m+1,2),DAT(m+a,2),DAT(m,2)/)
                theta_R_v(m) =  dot_product(dxi_dv, vec) - dbdir
             
!! 2nd ELEMENT

                dbdir = 0d0
                dSdir = 0d0

                if (jj==1) then
                    dkdk = 1d0
                else
                    call harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj),dkdk)
                end if
        
                dxi_dv(1) = dxi_dk(m+1,2)*dkdk*theta_k_v(ii+1,jj+1)!west

                if (ii==a-1) then
                    dkdk = 1d0
                   ! if (( ((jj-1)*dy >= 0.003d0)  .and. ((jj-1)*dy <= 0.007d0))) then
                   ! call harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj), dharm)
                   ! call getdir(d,ii,jj,dy,dx,Ts,dharm,theta_k_v,dbdir,dSdir) 
                   ! else
                    dbdir =  0d0
                    dSdir = 0d0
                   ! end if
                else
                    call harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+2,jj+1),dkdk)
                    dbdir =  0d0
                    dSdir = 0d0
                end if
                dxi_dv(2) = dxi_dk(m+1,5)*dkdk*theta_k_v(ii+1,jj+1)!north
                dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2)) + dSdir
                vec = (/DAT(m,2),DAT(m+1+a,2),DAT(m+1,2)/)
                theta_R_v(m+1) =  dot_product(dxi_dv, vec) - dbdir
               
                dbdir =  0d0
                dSdir = 0d0
!! 3rd ELEMENT
                if (jj==b-1) then
                    dkdk = 1d0
                else
                    call harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj+2),dkdk)
                endif

                dxi_dv(1) = dxi_dk(m+a,4)*dkdk*theta_k_v(ii+1,jj+1)!east
                if (ii==1) then
                    dkdk = 1d0
                    if (( ((jj)*dy >= 0.003d0)  .and. ((jj)*dy <= 0.005d0)))  then
                    call harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj+2), dharm)
                    call getdir(d,ii,jj,dy,dx,Ts,dharm,theta_k_v,dbdir,dSdir) 	
                     else
                    dbdir =  0d0
                    dSdir = 0d0
                    endif
                else
                    call harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii,jj+1),dkdk)
                    dbdir = 0d0
                    dSdir = 0d0
                endif
                dxi_dv(2) = dxi_dk(m+a,1)*dkdk*theta_k_v(ii+1,jj+1)!south
                dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2))+dSdir
                vec = (/DAT(m+a+1,2),DAT(m,2),DAT(m+a,2)/)
                theta_R_v(m+a) =  dot_product(dxi_dv, vec) - dbdir
!! 4th ELEMENT

                dbdir =  0d0
                dSdir = 0d0
                if (jj==b-1) then
                    dkdk = 1d0
                else
                    call harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj+2),dkdk)
                endif
        
                dxi_dv(1) = dxi_dk(m+a+1,2)*dkdk*theta_k_v(ii+1,jj+1)!west

                if (ii==a-1) then
                    dkdk = 1d0
                    !if (( ((jj)*dy >= 0.003d0)  .and. ((jj)*dy <= 0.007d0))) then
                !	call harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+1,jj+2), dharm)
                !       call getdir(d,ii,jj,dy,dx,Ts,dharm,theta_k_v,dbdir,dSdir) 
                !   else
                    dbdir =  0d0
                    dSdir = 0d0
                !   endif
                else
                    call harmonic_drv(IDAT(ii+1,jj+1),IDAT(ii+2,jj+1),dkdk)
                    dbdir = 0d0
                    dSdir = 0d0
                endif
                dxi_dv(2) = dxi_dk(m+a+1,1)*dkdk*theta_k_v(ii+1,jj+1)!south
                dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2)) + dSdir
                vec = (/DAT(m+a,2),DAT(m+1,2),DAT(m+a+1,2)/)
                theta_R_v(m+a+1) =  dot_product(dxi_dv, vec) - dbdir
                GRAD(l) = dot_product(PSI,theta_R_v)
        enddo
       ! open (unit = 3, file = "symmetric_grad")
       ! do l=1,N
       ! write(3,*)GRAD(l)
       ! enddo
       !  open (unit = 4, file = "PSI")
       ! do l=1,N
       ! write(4,*)PSI(l)
       ! enddo
       ! close(4)
        end subroutine


 subroutine EV_G(N, X, NEW_X, M, G, DAT, IDAT, IERR)
      implicit none
      integer N, NEW_X, M
      double precision G(M), X(N), ONES(N), lengthx
            
      REAL(8) DAT( 3*d**2-2*d, 2),IDAT(d+1,d+1)
 
       integer(8) IERR

        lengthx=0.005d0
      if (NEW_X) then
                    call fvm_simulate(X,d,DAT,IDAT)
         !   write(*,*) "gamw to mouni tis manas sas"
      end if


      
            
      ONES=(lengthx)**2/(1d0*N) 
      G(1)=dot_product(ONES,X)
      IERR = 0
      !write(*,*) "ev_g exit, cnst",G(1)
      return
      end

      subroutine EV_JAC_G(TASK, N, X,NEW_X, M, NZ, ACON, AVAR, A,DAT,IDAT, IERR)
      integer TASK, N,  M, NZ,IERR,NEW_X
       REAL(8) DAT( 3*d**2-2*d, 2),IDAT(d+1,d+1)

      double precision  A(NZ),X(N)
      integer ACON(NZ), AVAR(NZ), I
 integer AVAR1(N), ACON1(N)
 !     save  AVAR1, ACON1

      
!     structure of Jacobian:

     
      do I=1,N
         AVAR1(I)=I
      enddo
     ! data  AVAR1 /1, 2, 3, 4, 1, 2, 3, 4/
        ACON1= 1
      !data  ACON1 /1, 1, 1, 1, 2, 2, 2, 2/

      if( TASK.eq.0 ) then
        do I = 1, N
          AVAR(I) = AVAR1(I)
          ACON(I) = ACON1(I)
        enddo
      else
       do I=1,N
        A(I) = (0.01d0)**2/N
        enddo 
      endif
        IERR = 0
      return
      !write (*,*) "JAC G"
           end

END MODULE




