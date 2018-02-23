Program Berreman4x4

  Use io

  Implicit None
  Double Precision, Parameter :: Pi=3.141592653589793D0, eta0=376.98_8
  Complex*16, Parameter ::  I=(0.0D0,1.0D0)
  Double Precision, parameter :: theta=0.0D0, ng=1.5D0 
  Double Precision :: lz, alpha=0.0D0
  Double Precision :: eperpa,epara
  Double Precision :: no=1.62D0, ne=1.72D0
  Double Precision :: n_ito
  Double Precision :: lamb, lamb_i,lamb_f,dlamb
  Double Precision :: kx, k0, Xi
  Double Precision :: n(3) ,p0,q0,RWORK(8)
  integer :: ii,jj,kk,IPIV(4),INFO,Nz, tt,LWORK=16
  Complex*16, dimension(4,4)  :: Bij, Aij,Qij,Pij, CH_matrix, q_diag_matrix=(0.D0,0D0)
  Complex*16, dimension(4,4)  :: Atij, Arij, Iij
  Complex*16, dimension(4)    :: q_i
  Complex*16, dimension(4,1)  :: gamma_i, psi_i, psi_tr
  Complex*16, dimension(2)    :: E_i
  Complex*16 :: WORK(12),Vl(4,4),Vr(4,4)
  Double Precision,Dimension(3,3) :: eij
  Double Precision :: deltaE,dz,dz_ito, dz_polimide
  Double Precision :: transmitance, reflectance, incidence  
  Double Precision, Allocatable :: phi(:)
  Character (Len=32) :: passed_value
  Integer :: number_of_passed_arguments


  print*, "Welcome to  Berreman4x4 software V0.1 (23/02/18)"
  number_of_passed_arguments = COMMAND_ARGUMENT_COUNT()


  !Check if the user passed the right number of arguments:
  if ( number_of_passed_arguments < 2) then

     print*, "Error: Wrong number of arguments. You need to pass at least 2 arguments to this software. The remaining ones will be ignored."
     print*, "Usage: berreman4x4 number_of_lc_layers  cell_length"

       call EXIT(0)

  end if
  

  !Parsing input line arguments:  
  call get_command_argument(1,passed_value)
  read (passed_value,*) lz

  
  call get_command_argument(2,passed_value)
  read (passed_value,'(I10)') Nz


  !Parse data imput file:
  call parse_input_file(Nz,lz)
  call open_data_file()
  
  
  alpha=alpha*Pi/180.0D0
  
  lamb_i=0.35D0
  lamb_f=0.75D0
  dlamb=0.002D0

  
  p0=0.250D0
  q0=2*Pi/p0

  
  Allocate( phi(Nz) )
  !dz_polimide=0.98
  !dz_ito=0.025
  dz=lz/(Nz)
  E_i(1)=(1.0D0,0.0D0)
  E_i(2)=(0.0D0,0.0D0)


  incidence=Real( E_i(1)*conjg(E_i(1))/(cos(alpha)**2)+E_i(2)*conjg(E_i(2)) ) 


  eperpa=no**2
  epara=ne**2
  deltaE=epara-eperpa
  Xi=dsin(alpha)


  Iij=(0.0D0,0.0D0)
  forall(ii=1:4)   Iij(ii,ii)=(1.0D0,0.0D0)
  
  !Initiating input light bean
  psi_i(1,1)=E_i(1)
  psi_i(2,1)=ng*E_i(1)/cos(alpha)
  psi_i(3,1)=E_i(2)
  psi_i(4,1)=ng*cos(alpha)*E_i(2)


  

  Atij=(0.0D0,0.0D0)
  Atij(1,1)=(1.0D0,0.0D0)
  Atij(2,1)=ng/cos(alpha)
  Atij(3,2)=(1.0D0,0.0D0)
  Atij(4,2)=ng*cos(alpha)

  Arij=(0.0D0,0.0D0)
  Arij(1,3)=(1.0D0,0.0D0)
  Arij(2,3)=-ng/cos(alpha)
  Arij(3,4)=(1.0D0,0.0D0)
  Arij(4,4)=-ng*cos(alpha)


  
  time_DO: Do tt=1,1

     lamb=lamb_i
     call read_next_snapshot(phi,Nz)     
     
     omega_DO : Do while(lamb .le. lamb_f)

        Bij=Iij


        !Initiating electrical parameters:        
        k0=2.0D0*Pi/lamb
        kx=k0*dsin( alpha )
        

        
        !Filing the ito Berreman matrix:
        n_ito=2.525-1.271*lamb
        !call Fill_isotropic_berreman_matrix(Pij,n_ito,kz,dz_ito,k0)
        !Bij=matmul(Pij,Bij)


        !Filling the polimide Berreman matrix:
        !call Fill_isotropic_berreman_matrix(Pij,np,kz,dz_polimide,k0)
        !Bij=matmul(Pij,Bij)


        !Initiating the the system input light:

        lc: Do jj=1,Nz
           !Filling Berreman matrix Qij (dont confuse with the Lc order parameter):
           
           n(1)=cos(theta)*cos(phi(jj))
           n(2)=cos(theta)*sin(phi(jj))
           n(3)=0.0D0
           
           !n(1)=dcos( q0*dz*(jj-1) )
           !n(2)=dsin( q0*dz*(jj-1) )
           !n(3)=0.0D0
           
           eij(1,1)=eperpa+deltaE*n(1)**2
           eij(1,2)=deltaE*n(1)*n(2)
           eij(1,3)=deltaE*n(1)*n(3)

           eij(2,1)=eij(1,2)
           eij(2,2)=eperpa+deltaE*n(2)**2
           eij(2,3)=deltaE*n(2)*n(3)
           
           eij(3,1)=eij(1,3)
           eij(3,2)=eij(2,3)
           eij(3,3)=eperpa+deltaE*n(3)**2


           
           Qij(1,1)= -Xi*eij(3,1)/eij(3,3)
           Qij(1,2)=  -Xi**2/eij(3,3)+1D0           
           Qij(1,3)= -Xi*eij(3,2)/eij(3,3) 
           Qij(1,4)= (0.0D0,0.0D0 )

           Qij(2,1)=-eij(1,3)*eij(3,1)/eij(3,3) +eij(1,1)
           Qij(2,2)=-Xi*eij(1,3)/eij(3,3)
           Qij(2,3)=-eij(1,3)*eij(3,2)/eij(3,3)+eij(1,2)
           Qij(2,4)=(0.0D0,0.0D0)

           Qij(3,1)=(0.0D0,0.0D0)
           Qij(3,2)=(0.0D0,0.0D0)
           Qij(3,3)=(0.0D0,0.0D0)
           Qij(3,4)=(1D0,0D0)

           Qij(4,1)=-eij(2,3)*eij(3,1)/eij(3,3)+eij(2,1)
           Qij(4,2)=-Xi*eij(2,3)/eij(3,3)
           Qij(4,3)=-Xi**2-eij(2,3)*eij(3,2)/eij(3,3)+eij(2,2)
           Qij(4,4)=(0.0D0,0.0D0)

           !Calculating the Berraman matrix eigenvalues and eigenvectors:           
           
           call zgeev('V','V',4,Qij,4,q_i,VL,4,VR,4,WORk,LWORK,RWORK,INFO)
           if( INFO .ne. 0) then
           
           
              print*,  "convergence failure at Nz=", Nz,"info=",INFO, "lamb=",lamb 
              stop
           
           end if
           
           Vl=transpose(Vl)
           
           call normalize_left_and_right_eigenvectors(Vl,VR,4,1)

           q_diag_matrix=(0D0,0D0)
           forall(kk=1:4) q_diag_matrix(kk,kk)=zexp(-dz*k0*q_i(kk)*I)
           
           
              
           Pij=matmul( VR , matmul(q_diag_matrix,Vl) )

          

           
           Bij=matmul(Pij,Bij)

           
        end Do lc



        Aij=Atij-matmul(Bij,Arij)
        psi_tr=matmul(Bij,psi_i)

        

        call ZGESV( 4, 1, Aij, 4, IPIV, psi_tr, 4, INFO )
        if( INFO .ne. 0) print*, "info=",INFO, "lamb=" ,lamb , "deu ruim (berreman vector)!!!"
        
        transmitance=Real( psi_tr(1,1)*conjg(psi_tr(1,1))/(cos(alpha)**2)+psi_tr(2,1)*conjg(psi_tr(2,1)) )/incidence
        reflectance=Real( psi_tr(3,1)*conjg(psi_tr(3,1))/(cos(alpha)**2)+psi_tr(4,1)*conjg(psi_tr(4,1)) )/incidence


        write(70,*) lamb, (transmitance), (reflectance), (transmitance)+(reflectance)


        lamb=lamb+dlamb

     end Do omega_DO

     write(70,*)
     write(70,*)
     
  end Do time_DO


  call Close_data_files()


contains

  Subroutine normalize_left_and_right_eigenvectors (left_vectors,right_vectors,N,M)

    Implicit None
    Integer, Intent(IN) :: M,N
    Complex*16, Intent(INOUT) :: left_vectors(N*M,N*M), right_vectors(N*M,N*M)
    Integer :: ii
    Double Precision :: Norm

    Do ii=1,N*M


       norm=real(sum(left_vectors(ii,:)*conjg(right_vectors(:,ii))))

       
       if(Real(norm) < 0D0) then

          norm=sqrt(-norm)
          left_vectors(ii,:)=left_vectors(ii,:)/norm
          right_vectors(:,ii)=-right_vectors(:,ii)/norm

       else

          norm=sqrt(norm)
          left_vectors(ii,:)=left_vectors(ii,:)/norm
          right_vectors(:,ii)=right_vectors(:,ii)/norm

       end if

    end do
    
  end Subroutine normalize_left_and_right_eigenvectors
  

!  Subroutine Fill_isotropic_berreman_matrix(Pij,ng,alpha,dz,k0)
!
!    Implicit None
!    Complex*16, Intent(out) :: Pij(4,4)
!    Double Precision, intent(in) :: ng,dz, k0, alpha
!    Double Precision :: k, kz
!    kz=k0*dcos( alpha )
!    k=ng*k0
!
!    Pij=cmplx(0.0_8,0.0_8,kind=8)
!
!    Pij(1,1)=cmplx(cos(kz*dz),0.0_8,Kind=8)
!    Pij(1,2)=cmplx(0.0_8,-sin(kz*dz)*kz/(k*ng),Kind=8)
!
!    Pij(2,1)=cmplx(0.0_8,-sin(kz*dz)*(k*ng)/kz,Kind=8)
!    Pij(2,2)=cmplx(cos(kz*dz),0.0_8,Kind=8)
!
!    Pij(3,3)=cmplx(cos(kz*dz),0.0_8,Kind=8)
!    Pij(3,4)=cmplx(0.0_8,-sin(kz*dz)*k/(kz*ng),Kind=8)
!
!    Pij(4,3)=cmplx(0.0_8,-sin(kz*dz)*(kz*ng)/k,Kind=8)
!    Pij(4,4)=cmplx(cos(kz*dz_ito),0.0_8,Kind=8)
!
!  end Subroutine Fill_isotropic_berreman_matrix



end Program Berreman4x4



