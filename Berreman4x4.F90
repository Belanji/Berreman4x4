Program Berreman4x4

  Use io

  Implicit None
  Double Precision, Parameter :: Pi=3.141592653589793D0, eta0=376.98_8
  Complex*16, Parameter ::  I=(0.0D0,1.0D0)
  Double Precision, parameter :: theta=0.0D0, ng=1.5D0 
  Double Precision :: lz, alpha=0.0D0, beta=0D0
  Double Precision :: eperpa,epara
  Double Precision :: no=1.52D0, ne=1.78D0
  Double Precision :: n_ito
  Double Precision :: lamb, lamb_i,lamb_f,dlamb
  Double Precision :: kx, k0, Xi, ky
  Double Precision :: n(3) ,p0,q0,RWORK(8)
  integer :: ii,jj,IPIV(4),INFO,Nz, tt
  Complex*16, dimension(4,4)  :: Bij, Aij,Pij, CH_matrix
  Complex*16, dimension(4,4)  :: Atij, Arij, Iij
  Complex*16, dimension(4,1)  :: gamma_i, psi_i, psi_tr
  Complex*16, dimension(2)    :: E_i
  Double Precision, dimension(4)   :: q_i
  Double Precision, dimension(4,4) :: Qij
  Complex*16 :: WORK(12),Vl(4,4),Vr(4,4)
  Double Precision,Dimension(3,3) :: eij
  Double Precision :: deltaE,dz,dz_ito, dz_polimide
  Double Precision :: transmitance, reflectance, incidence  
  Double Precision, Allocatable :: phi(:)



  call parse_input_file(Nz,lz)
  call open_data_file()
  Allocate( phi(Nz) )
  
  alpha=alpha*Pi/180.0D0
  beta=beta*Pi/180.0D0
  
  lamb_i=0.3D0
  lamb_f=0.7D0
  dlamb=0.001D0


  !lamb_i=0.52D0
  !lamb_f=0.700D0
  !dlamb=0.5D0

  
  p0=0.338D0
  Nz=5000
  lz=50.7D0
  q0=2*Pi/p0


  dz_polimide=0.98
  dz_ito=0.025
  dz=lz/(Nz)
  E_i(1)=(1.0D0,0.0D0)
  E_i(2)=(0.0D0,0.0D0)


  incidence=Real( E_i(1)*conjg(E_i(1))/(cos(alpha)**2)+E_i(2)*conjg(E_i(2)) ) 

  !print*, incidence

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
     !call read_next_snapshot(phi,Nz)     
     
     omega_DO : Do while(lamb .le. lamb_f)

        Bij=Iij


        !Initiating electrical parameters:        
        k0=2.0D0*Pi/lamb
        kx=1.5D0*dsin( alpha )*dcos( beta )
        ky=1.5D0*dsin( alpha )*dsin( beta )

        
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
           
           !n(1)=cos(theta)*cos(phi(jj))
           !n(2)=cos(theta)*sin(phi(jj))
           !n(3)=sin( theta )
           n(3)=0.0D0
           n(1)=dcos( q0*dz*(jj-1) )
           n(2)=dsin( q0*dz*(jj-1) )

           eij(1,1)=eperpa+deltaE*n(1)**2
           eij(1,2)=deltaE*n(1)*n(2)
           eij(1,3)=deltaE*n(1)*n(3)

           eij(2,1)=eij(1,2)
           eij(2,2)=eperpa+deltaE*n(2)**2
           eij(2,3)=deltaE*n(2)*n(3)
           
           eij(3,1)=eij(1,3)
           eij(3,2)=eij(2,3)
           eij(3,3)=eperpa+deltaE*n(3)**2


           
           Qij(1,1)= -kx*eij(3,1)/eij(3,3)
           Qij(1,2)= 1D0-kx**2/eij(3,3)
           Qij(1,3)= -kx*eij(3,2)/eij(3,3) 
           Qij(1,4)= -kx*ky/eij(3,3)

           Qij(2,1)=eij(1,1)-eij(1,3)*eij(3,1)/eij(3,3)-ky*ky
           Qij(2,2)=Qij(1,1)
           Qij(2,3)=eij(1,2)+kx*ky-eij(1,3)*eij(3,2)/eij(3,3)
           Qij(2,4)=-ky*eij(1,3)/eij(3,3)

           Qij(3,1)=-ky*eij(1,3)/eij(3,3)
           Qij(3,2)=-kx*ky/eij(3,3)
           Qij(3,3)=-ky*eij(2,3)/eij(3,3)
           Qij(3,4)=1D0-kx*ky/eij(3,3)

           Qij(4,1)=Qij(2,3)
           Qij(4,2)=Qij(1,3)
           Qij(4,3)=eij(2,2)-kx**2-eij(2,3)*eij(3,2)/eij(3,3)
           Qij(4,4)=-ky*eij(2,3)/eij(3,3)

           !Calculating the Berraman matrix eigenvalues:           

           !Aij=Qij
           !call zgeev('N','N',4,Aij,4,q_i,VL,4,VR,4,WORk,12,RWORK,INFO)
           !if( INFO .ne. 0) then
           !
           !
           !   print*,  "info=",INFO, "Nz=",Nz," lamb=",lamb
           !   stop
           !
           !end if
                       
           q_i(1)= sqrt(eperpa-Xi**2)
           q_i(2)=-sqrt(eperpa-Xi**2)

           q_i(3)=(-Xi*deltaE*n(1)*n(3) + sqrt(eperpa*epara)*sqrt(  (eperpa+deltaE*n(3)**2)-Xi**2*( 1.0_8-deltaE*n(2)**2/epara )  ))/(eperpa+deltaE*n(3)**2)
           q_i(4)=(-Xi*deltaE*n(1)*n(3) - sqrt(eperpa*epara)*sqrt(  (eperpa+deltaE*n(3)**2)-Xi**2*( 1.0_8-deltaE*n(2)**2/epara )  ))/(eperpa+deltaE*n(3)**2)
           !q_i(3)=Qij(1,1) + sqrt( Qij(1,2)*Qij(2,1)+Qij(1,3)*Qij(2,3)/Qij(1,1) )
           !q_i(4)=Qij(1,1) - sqrt( Qij(1,2)*Qij(2,1)+Qij(1,3)*Qij(2,3)/Qij(1,1) )
           

           
           !Writing the sistem of equations for the Caley-Hamilton parameters:
           Do ii=1,4

              CH_matrix(ii,1)=  ( 1.0D0,0.0D0 )
              CH_matrix(ii,2)=  ( -I*k0*dz*q_i(ii) )
              CH_matrix(ii,3)=  ( -I*k0*dz*q_i(ii) )**2
              CH_matrix(ii,4)=  ( -I*k0*dz*q_i(ii) )**3
              gamma_i(ii,1)=cos(dz*k0*q_i(ii) )-I*sin( dz*k0*q_i(ii) )

           end Do

           !print*, gamma_i
           !print*, "   "
           
           
           call ZGESV( 4, 1, CH_matrix, 4, IPIV, gamma_i, 4, INFO )

           


           if( INFO .ne. 0) then


              print*, "info=",INFO, "Nz=",Nz, "deu ruim!!!"
              stop

           end if
           
              
           Pij=gamma_i(1,1)*Iij + gamma_i(2,1)*( -I*dz*k0 )*Qij + gamma_i(3,1)*( ( -I*dz*k0 )**2)*matmul(Qij,Qij) &
                & +gamma_i(4,1)*( ( -I*dz*k0 )**3)*matmul(Qij,matmul(Qij,Qij))

          

           
           Bij=matmul(Pij,Bij)

           !print*, "jj=", jj
           !print*, Bij
           !print*, " "
           !print*, " "
           
        end Do lc

        !Filling the polimide Berreman matrix:
        !call Fill_isotropic_berreman_matrix(Pij,np,kz,dz_polimide,k0)
        !Bij=matmul(Pij,Bij)



        !Filing the ito Berreman matrix:
        !n_ito=2.525-1.271*lamb
        !call Fill_isotropic_berreman_matrix(Pij,n_ito,kz,dz_ito,k0)
        !Bij=matmul(Pij,Bij)


        Aij=Atij-matmul(Bij,Arij)
        psi_tr=matmul(Bij,psi_i)



        !Abij=Atij+matmul(Bij,Arij)
        !call matinv4(Abij,Abinv_ij)

        call ZGESV( 4, 1, Aij, 4, IPIV, psi_tr, 4, INFO )

        if( INFO .ne. 0) print*, "info=",INFO, "lamb=" ,lamb , "deu ruim!!!"
        
        transmitance=Real( psi_tr(1,1)*conjg(psi_tr(1,1))/(cos(alpha)**2)+psi_tr(2,1)*conjg(psi_tr(2,1)) )/incidence
        reflectance=Real( psi_tr(3,1)*conjg(psi_tr(3,1))/(cos(alpha)**2)+psi_tr(4,1)*conjg(psi_tr(4,1)) )/incidence


        write(70,*) lamb, (transmitance), (reflectance), (transmitance)+(reflectance)


        lamb=lamb+dlamb

     end Do omega_DO

     write(70,*)
     write(70,*)
     
  end Do time_DO


  call Close_data_files()


!contains


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



