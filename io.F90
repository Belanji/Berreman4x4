module io

    !This module contains the routines used to read the data form the archive "constants.dat"
  
  Implicit none
  Character(len=10), Private, Parameter :: input_file_name="data.input"
  character(len=8), Parameter, public ::  run_info="run.log"
  Integer, private :: error_handling 

  
Contains

  
  Subroutine open_data_file()

    Implicit None


    Open(unit=50, file="phi_time.dat",status="old",iostat=error_handling,position="rewind",action="read")

    If( error_handling  /= 0)  Stop "Failed to open the data file (in the open_data_file subroutine)."

    
    Open(unit=70, file="ref_trans_time.dat",status="replace",iostat=error_handling,action="write")
    
    If( error_handling  /= 0)  Stop "Failed to create the output file ( int the open_data_file subroutine)."



    !Open(unit=700, file="testiando.dat",status="replace",iostat=error_handling,action="write")

  end Subroutine open_data_file


  function read_next_snapshot(phi,Nz) result(read_status)

    Implicit None
    Integer, Intent(in):: Nz
    Double precision   :: phi(:)
    Character(len=40) :: time
    Double Precision :: z
    Integer :: ii, read_status

    
    read(50,*,Iostat=read_status) time

    read_if: if(read_status == -1 ) then

       print*, "  "
       print*, "Berreman4x4 sucessfully executed."
       print*, "Finishing the program."

       
    else
   

       write(70,*) time
       
       Write(*,*) "Reading snapshot for ", time(2:30)
    
       Do ii=1,Nz

          read(50,*,Iostat=read_status)    z , phi(ii)

          !Check if data was successfully read:
          if(read_status /= 0 ) then

             Print*, "Failure reading the snapshot. at line", ii
             Print*, "Aborting the pfrogram."
             exit
             
          end if


       end Do

       
       read(50,*)
       read(50,*)

       
    end if read_if
    
    
       
    
  end function read_next_snapshot
  
  

  
  Subroutine Close_data_files()

    Close(50)
    close(70)
    
  end Subroutine Close_data_files
  
!  Subroutine parse_input_file(Nz,lz)
!
!    Implicit None
!    Integer, Intent(out):: Nz
!    Double precision, Intent(out) ::  lz
!    Character(len=40) :: character_holder
!    Double precision ::  first_argument !, second_argument
!    Integer :: error_handling
!    
!    Open(unit=90, file="data.input",status="old",iostat=error_handling,position="rewind",action="read")
!
!    !If( error_handling  /= 0)  Stop "Failed to open the data input file (in the open_data_file subroutine)."
!
!    !read(90,*) character_holder, first_argument
!    !Nz=int(first_argument)
!
!    
!    !call read_next_word(90,character_holder)
!    !read(90,FMT=*,ADVANCE='yes' ) Nz
!
!    !read(90,*) character_holder, first_argument
!    !lz=first_argument
!
!    
!    close(90)
!    
!    
!  end Subroutine parse_input_file


  !Subroutine read_next_word(ipFIle, word_holder)
  !
  !  Integer :: ipFile
  !  character(len=40) :: word_holder
  !
  !  
  !  !read(90,FMT=*,ADVANCE='no' ) character_holder
  !
  !end Subroutine read_next_word
  
    
END module io

  
