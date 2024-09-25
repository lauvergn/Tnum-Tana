
     integer, parameter :: ncart_act = 9
     real(8) :: Hcc(ncart_act,ncart_act)
     integer :: i,j

     read(5,*) ((Hcc(i,j),i=1,j),j=1,ncart_act)
     DO j=1,ncart_act
     DO i=1,j-1
       Hcc(j,i) = Hcc(i,j)
     END DO
     END DO

     DO j=1,ncart_act
       write(6,*) j,Hcc(:,j)
     END DO

END
