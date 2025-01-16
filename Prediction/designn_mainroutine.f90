        program designnmain
        implicit none
        real, allocatable,dimension(:) :: xa,ya,za,x,y,z
        character(len=2),allocatable,dimension(:) :: aatms,atms
        character(len=50) :: filename
        character(len=50) :: full_file,shell_file
        integer anAtoms,nAtoms
        
        shell_file='Shell'
        full_file='Full'

        call line_counter(full_file,anAtoms)
        allocate(xa(anAtoms),ya(anAtoms),za(anAtoms),aatms(anAtoms))
        call file_reading(aatms,xa,ya,za,anAtoms,full_file)

        call line_counter(shell_file,nAtoms)
        allocate(x(nAtoms),y(nAtoms),z(nAtoms),atms(nAtoms))
        call file_reading(atms,x,y,z,nAtoms,shell_file)

        call mapping(anAtoms,nAtoms,xa,ya,za,x,y,z,atms,aatms,& 
     & shell_file,full_file)

        stop
        end 

        


