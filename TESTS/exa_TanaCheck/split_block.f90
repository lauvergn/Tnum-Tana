SUBROUTINE read_arg()
  USE rec_vec_m
  IMPLICIT NONE
  integer :: nv


  character(len=:), allocatable :: arg,arg2
  integer :: long,i,err_read

  namelist / rec_vec / type_vec,max_layer,max_v,max_b,nrho,cart,   &
                       cos_th,frame0vec,type100,zmat_order,SpheConv_xzy


  DO i=1, COMMAND_ARGUMENT_COUNT(),2
    CALL GET_COMMAND_ARGUMENT( NUMBER=i, LENGTH=long )
    allocate( character(len=long) :: arg )
    CALL GET_COMMAND_ARGUMENT( NUMBER=i, VALUE=arg )

    CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, LENGTH=long )
    allocate( character(len=long) :: arg2 )
    CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, VALUE=arg2 )

    SELECT CASE(arg)
    CASE("-type_vec","-tv")
      read(arg2,*) type_vec
    CASE("-v","-vec")
      read(arg2,*) max_v
    CASE("-b","-block")
      read(arg2,*) max_b
    CASE("-l","-layer")
      read(arg2,*) max_layer
    CASE("-cart")
      read(arg2,*) cart
    CASE("-cos_th")
      read(arg2,*) cos_th
    CASE("-frame0vec")
      read(arg2,*) frame0vec
    CASE("-nrho")
      read(arg2,*) nrho
    CASE("-SpheConv_xzy","-SCxzy")
      read(arg2,*) SpheConv_xzy
    END SELECT

    print *,"Argument de rang ", i, " ==> ", arg , ' arg_val: ',arg2

    deallocate(arg)
    deallocate(arg2)
  END DO



  read(5,rec_vec,IOSTAT=err_read)
  IF (err_read /= 0) THEN
    write(6,*) ' WARNING: no namelist or errors while reading the namelist "rec_vec"'
  END IF
  !write(6,rec_vec)
  write(6,*) 'type_vec  ',type_vec
  write(6,*) 'max_v     ',max_v
  write(6,*) 'max_b     ',max_b
  write(6,*) 'max_layer ',max_layer
  write(6,*) 'cart,cos_th,frame0vec,SpheConv_xzy ',cart,cos_th,frame0vec,SpheConv_xzy
  write(6,*) 'nrho      ',nrho



  nv = max_v
  write(6,*) 'split ',nv ; flush(6)

  CALL Split_blocks_new(nv)

  write(6,*) '=================================='
  write(6,*) '=================================='

END SUBROUTINE read_arg

PROGRAM test
  USE rec_vec_m
  IMPLICIT NONE
  integer :: nv


  character(len=:), allocatable :: arg,arg2
  integer :: long,i

  write(6,*) '=================================='
  write(6,*) 'AUTOMATIC GENERATION of Tana data'
  write(6,*) '=================================='
  write(6,*) 'Code written by David Lauvergnat [1]'
  write(6,*) '[1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France'
  write(6,*) '=================================='

  CALL read_arg()

  nv = max_v
  write(6,*) 'split ',nv ; flush(6)

  CALL Split_blocks_new(nv)

  write(6,*) '=================================='
  write(6,*) '=================================='

END PROGRAM test
