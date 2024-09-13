!===============================================================================
!===============================================================================
!This file is part of FOR_EVRT library.
!
!===============================================================================
! MIT License
!
! Copyright (c) 2022 David Lauvergnat
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!===============================================================================
!===============================================================================
!=======================================================================================
!> @brief 
!> module for implementating MPI
!> to use long integer, use 64-bit compiled MPI
!>
!> include:
!> basic MPI variables,
!> asisstant MPI variables,
!> varilables for recording time usage,
!> varilables for process control (to be added to certain type later),
!>  
!> @note there is a bug in openMPI for working on 64-bit fortran USE MPI_F08
!>        which will possiblely be solved in future release MPI_v4.x or MPI_v5.x.
!>        The current strategy apply USE MPI or 
!>        keep using 32-bit with extra subroutine for large array, need further test
!
!> @note further implemention: using 32-bit MPI in 64-bit Fortran
!=======================================================================================
MODULE mod_MPI
#if(run_MPI)
  USE MPI
  !USE MPI_F08
  IMPLICIT NONE

  !> note the difference in "USE MPI" and "USE MPI_F08"
  !-USE MPI-----------------------------------------------------------------------------
  Integer(kind=MPI_INTEGER_KIND)                 :: MPI_stat(MPI_STATUS_SIZE) !< status of MPI process
  Integer(kind=MPI_INTEGER_KIND)                 :: Int_MPI                   !< integer type of default MPI
  Integer(kind=MPI_INTEGER_KIND)                 :: Real_MPI                  !< real type for fortran (Rkind)
  Integer(kind=MPI_INTEGER_KIND)                 :: Cplx_MPI                  !< complex type for fortran (Rkind)
  Integer(kind=MPI_INTEGER_KIND)                 :: Int_fortran               !< integer type of default fortran
  Integer(kind=MPI_INTEGER_KIND)                 :: MPI_NODE_0_WORLD
  Integer(kind=MPI_INTEGER_KIND)                 :: MPI_NODE_0_GROUP
  Integer(kind=MPI_INTEGER_KIND)                 :: MPI_NODE_0_COMM           !< communcator for masters on each nodes
  Integer(kind=MPI_INTEGER_KIND)                 :: MPI_COMM_current
  !-USE MPI_F08-------------------------------------------------------------------------
!  TYPE(MPI_Status)                              :: MPI_stat                  !< status of MPI process
!  TYPE(MPI_Datatype)                            :: Int_MPI                   !< integer type of default MPI
!  TYPE(MPI_Datatype)                            :: Real_MPI                  !< real type for fortran (Rkind)
!  TYPE(MPI_Datatype)                            :: Cplx_MPI                  !< complex type for fortran (Rkind)
!  TYPE(MPI_Datatype)                            :: Int_fortran               !< integer type of default fortran
!  TYPE(MPI_Group)                               :: MPI_NODE_0_WORLD
!  TYPE(MPI_Group)                               :: MPI_NODE_0_GROUP
!  TYPE(MPI_Comm)                                :: MPI_NODE_0_COMM           !< communcator for masters on each nodes
!  TYPE(MPI_Comm)                                :: MPI_COMM_current

  !-------------------------------------------------------------------------------------
  Integer(kind=MPI_INTEGER_KIND)                 :: MPI_NODE_0_ERR

  Integer(kind=MPI_INTEGER_KIND)                 :: MPI_err                   !< error flag for MPI
  Integer(kind=MPI_INTEGER_KIND)                 :: MPI_id                    !< rocess ID, 0~MPI_np-1
  Integer(kind=MPI_INTEGER_KIND)                 :: MPI_id_node=-1            !< rocess ID in node comm
  Integer(kind=MPI_INTEGER_KIND)                 :: MPI_np                    !< total number of MPI threads
  Integer(kind=MPI_INTEGER_KIND)                 :: MPI_np_node=-1            !< total number of MPI threads
  Integer(kind=MPI_INTEGER_KIND)                 :: MPI_tag1                  !< tag for MPI send and receive
  Integer(kind=MPI_INTEGER_KIND)                 :: MPI_tag2                  !< tag for MPI send and receive
  
  Integer(kind=MPI_INTEGER_KIND)                 :: source_MPI                !< for MPI_RECV, not used
  Integer                                        :: integer_MPI               !< get integer type for MPI

  Integer(kind=MPI_INTEGER_KIND)                 :: root_MPI=0                !< used for MPI functions, master ID 
  Integer(kind=MPI_INTEGER_KIND)                 :: size1_MPI=1               !< used for MPI functions, size 1
  Integer                                        :: ONE_1=1
  Integer(kind=MPI_INTEGER_KIND)                 :: i_MPI                     !< fake MPI thread id
  Integer                                        :: bound1_MPI                !< up boundary of works for threads
  Integer                                        :: bound2_MPI                !< dn boundary of works for threads
  Integer                                        :: nb_per_MPI                !< number of works per thread
  Integer                                        :: nb_rem_MPI                !< remainder of distribed works
  Integer,allocatable                            :: iGs_MPI(:,:)              !< iG boundary of theards in action
  Integer,allocatable                            :: iGs_MPI0(:,:)
  Integer,allocatable                            :: iGs_MPI_mc(:,:,:)         !< iG boundary of theards in action for memory control
  Integer,allocatable                            :: iGs_MPI_mc0(:,:,:)
  Integer,allocatable                            :: bounds_MPI(:,:)           !< boundary of each theards

  Integer                                        :: MPI_scheme=0              !< MPI working on MPI scheme #.; 0 auto
  Logical                                        :: keep_MPI=.FALSE.          !< if keep valables on current thread 
  Logical                                        :: Srep_MPI=.FALSE.          !< MPI working on full Smolyak Rep. 
  Logical                                        :: MPI_S2_L2=.FALSE.         !< enable level2 distribution of S2
  Integer                                        :: Num_L2=6                  !<1/Num_L2 threads used for level2 distributor
  Integer                                        :: n_level2                  !< number of level2 distributor 

  ! use type later to collect these parameters
  Integer                                        :: MPI_mc=1                  !< increase to reduce memory used in action
  Integer,allocatable                            :: MPI_nodes_np(:)           !< available np on each nodes
  Integer                                        :: MPI_nodes_num=0
  Character(LEN=MPI_MAX_PROCESSOR_NAME),allocatable :: MPI_nodes_name(:)      !< nodes name
  Character(LEN=MPI_MAX_PROCESSOR_NAME)          :: MPI_node_name
  Logical                                        :: MPI_nodes_p0=.FALSE.      !< the qusei-master p on each node
  Integer,allocatable                            :: MPI_nodes_p00(:)          !< information of MPI_nodes_p0 on master
  Logical,allocatable                            :: MPI_nodes_00(:)           !< if processor 0 on each node, known by MPI_id 0
  Integer                                        :: MPI_node_id               !< nodes rank
  Integer                                        :: MPI_node_p0_id            !< quasi-master id
  Integer                                        :: MPI_sub_id(2)             !< processors affiliated
  Integer                                        :: MPI_fake_nodes=0          !< fake nodes for S3
  Logical                                        :: MPI_iGs_auto=.FALSE.      !< if auto-adjust the distribution of Smolyak terms
  Integer                                        :: MPI_mem_node=0            !< in unit of GB

  TYPE MPI_S_type 
    Logical                                      :: davidson=.FALSE.
    Logical                                      :: prop=.FALSE.
    Integer                                      :: num_resetH=1
    Integer                                      :: nb_channels=1
    Integer                                      :: nb_psi=1                  ! nb of levels to converge
  END TYPE 
  TYPE(MPI_S_type)                               :: MPI_S

  !-------------------------------------------------------------------------------------
  !> varilables for recording time usage
  Integer                                        :: time_MPI_action=0         !< total time used in action 
  Integer                                        :: time_comm=0               !< total time used in comm. in action
  Integer,allocatable                            :: time_MPI_calcu(:)         !< time used in each action for sharable claculation
  Integer,allocatable                            :: time_MPI_local(:)         !< local time used in each action
  Integer,allocatable                            :: time_MPI_act_all(:)       !< time used in each action
  Integer                                        :: time_point1               !< tags for save current time
  Integer                                        :: time_point2
  Integer                                        :: time_temp1
  Integer                                        :: time_temp2
  Integer                                        :: time_rate                 !< for function system_clock()
  Integer                                        :: time_max                  !< for function system_clock()
  ! time_rate in kind=4: COUNT in system_clock represents milliseconds
  ! time_rate in kind>4: COUNT in system_clock represents micro- or nanoseconds
  
  !-------------------------------------------------------------------------------------
  ! varilables for process control, add to certain type later
  ! Logical                        :: Grid_allco=.True.

!=======================================================================================
  Contains  

  !-------------------------------------------------------------------------------------
  !> @brief MPI initialization 
  SUBROUTINE ini_MPI()
    USE mod_NumParameters
    IMPLICIT NONE

    CALL MPI_Init(MPI_err)
    CALL MPI_Comm_rank(MPI_COMM_WORLD,MPI_id,MPI_err)
    CALL MPI_Comm_size(MPI_COMM_WORLD,MPI_np,MPI_err)

    !> get Fortran default bit
    IF(sizeof(integer_MPI)==4) THEN
      Int_fortran=MPI_integer4
    ELSEIF(sizeof(integer_MPI)==8) THEN
      Int_fortran=MPI_integer8
    ELSE
      STOP 'integer neither 4 or 8 in default fortran'
    ENDIF
    
    !> define Fortran Real type according to Rkind
    IF(Rkind==8) THEN
      Real_MPI=MPI_Real8
      Cplx_MPI=MPI_Complex8
    ELSEIF(Rkind==16) THEN
      Real_MPI=MPI_Real16
      Cplx_MPI=MPI_Complex16
    ELSE
      STOP 'Rkind neither 64 or 128 bit, define Real_MPI in sub_module_MPI.f90'
    ENDIF
    
    !> get MPI default bit
    IF(MPI_INTEGER_KIND==4) THEN
      Int_MPI=MPI_integer4
    ELSEIF(MPI_INTEGER_KIND==8) THEN
      Int_MPI=MPI_integer8
    ELSE
      STOP 'integer neither 4 or 8 in default MPI'
    ENDIF

    IF(sizeof(integer_MPI)/=MPI_INTEGER_KIND) THEN
      write(out_unit,*) 'Please check the Integer type of Compiler and MPI used,'
      STOP 'Integer type of default Fortran Compiler and MPI not agree.'
    ENDIF

    write(out_unit,*) 'Initiaize MPI with ', MPI_np, 'cores.'
    write(out_unit,*) 'NOTE: MPI in progress. If get memory error, check if           &
                                      the variables are just allocated on root threads.'
    write(out_unit,*)
    write(out_unit,*) 'Integer type of default Fortran Compiler:',                    &
                                         sizeof(integer_MPI),', MPI: ',MPI_INTEGER_KIND
    IF(sizeof(integer_MPI)/=MPI_INTEGER_KIND)                                          &
                                 STOP 'Please use same integer type for Fortran and MPI'

  END SUBROUTINE ini_MPI

  !-------------------------------------------------------------------------------------
  !> @brief MPI_Finalize
  SUBROUTINE end_MPI()
    USE mod_NumParameters
    IMPLICIT NONE

    IF(MPI_id==0) THEN
      write(out_unit,*) 'time used for action: ',                                     &
             real(time_MPI_action,kind=Rkind)/real(time_rate,kind=Rkind),' from ',MPI_id
      write(out_unit,*) 'time used for MPI communication: ',                          &
                  real(time_comm,kind=Rkind)/real(time_rate,kind=Rkind),' from ', MPI_id
    ENDIF

    CALL MPI_Finalize(MPI_err)
  END SUBROUTINE end_MPI

  !-------------------------------------------------------------------------------------
  !> @brief get nodes information
  ! MPI_nodes_p0 are initialized, MPI_nodes_np ready
  SUBROUTINE get_nodes_info_MPI()
    USE mod_NumParameters
    IMPLICIT NONE

    Integer                                      :: act_len
    Integer                                      :: ii
    Integer                                      :: jj
    Logical                                      :: match
    Logical                                      :: nodes_p0_temp

    Integer                                      :: nb_per_node
    Integer                                      :: nb_rem_node
    Integer,ALLOCATABLE                          :: node_d1(:)
    Integer,ALLOCATABLE                          :: node_d2(:)

    allocate(MPI_nodes_name(0:MPI_np-1))
    allocate(MPI_nodes_np(0:MPI_np-1))
    IF(MPI_id==0) THEN
      allocate(MPI_nodes_p00(0:MPI_np-1))
      allocate(MPI_nodes_00(0:MPI_np-1))
      MPI_nodes_00=.FALSE.
    ENDIF
    MPI_nodes_np=0
    MPI_nodes_name=''

    IF(MPI_fake_nodes>0) THEN
      !-working with fake nodes---------------------------------------------------------
      ! checked
      MPI_nodes_num=MIN(MPI_fake_nodes,MPI_np)
      allocate(node_d1(0:MPI_nodes_num-1))
      allocate(node_d2(0:MPI_nodes_num-1))
      MPI_nodes_p0=.FALSE.

      nb_per_node=MPI_np/MPI_nodes_num
      nb_rem_node=mod(MPI_np,MPI_nodes_num)

      DO ii=0,MPI_nodes_num-1
        node_d1(ii)=ii*nb_per_node+1+MIN(ii,nb_rem_node)
        node_d2(ii)=(ii+1)*nb_per_node+MIN(ii,nb_rem_node)+merge(1,0,nb_rem_node>ii)
        MPI_nodes_np(ii)=node_d2(ii)-node_d1(ii)+1
        IF(MPI_id>=node_d1(ii)-1 .AND. MPI_id<=node_d2(ii)-1) THEN
          MPI_node_id=ii
          MPI_node_p0_id=node_d1(ii)-1
        ENDIF
        IF(MPI_id==node_d1(ii)-1) MPI_nodes_p0=.TRUE.
        IF(MPI_id==0) THEN
          MPI_nodes_p00(ii)=node_d1(ii)-1
          MPI_nodes_00(node_d1(ii)-1)=.TRUE.
        ENDIF
      ENDDO

      deallocate(node_d1)
      deallocate(node_d2)

    ELSE
      !-working on real multi-nodes-----------------------------------------------------
      CALL MPI_Get_processor_name(MPI_node_name,act_len,MPI_err)

      IF(MPI_id==0) THEN
        MPI_nodes_name(0)=MPI_node_name
        MPI_nodes_np(0)=MPI_nodes_np(0)+1
        MPI_node_id=0
        MPI_nodes_p00(0)=0
        MPI_nodes_p0=.TRUE.
        MPI_nodes_00(MPI_id)=.TRUE.
      ENDIF

      ii=0
      IF(MPI_id/=0) THEN
        CALL MPI_send(MPI_node_name,MPI_MAX_PROCESSOR_NAME,MPI_Character,root_MPI,     &
                      MPI_id,MPI_COMM_WORLD,MPI_err)
        CALL MPI_Recv(MPI_nodes_p0,size1_MPI,MPI_Logical,root_MPI,MPI_id,              &
                      MPI_COMM_WORLD,MPI_stat,MPI_err)
      ELSE
        DO i_MPI=1,MPI_np-1
          CALL MPI_Recv(MPI_node_name,MPI_MAX_PROCESSOR_NAME,MPI_Character,i_MPI,      &
                        i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
          
          nodes_p0_temp=.FALSE.

          ! check match
          match=.FALSE.
          DO jj=0,ii
            IF(MPI_node_name==MPI_nodes_name(jj)) THEN
              MPI_nodes_np(jj)=MPI_nodes_np(jj)+1
              match=.TRUE.
            ENDIF
          ENDDO

          IF(.NOT. match) THEN
            ii=ii+1
            MPI_nodes_name(ii)=MPI_node_name
            MPI_nodes_np(ii)=MPI_nodes_np(ii)+1
            nodes_p0_temp=.TRUE.
            MPI_nodes_p00(ii)=i_MPI
            MPI_nodes_00(i_MPI)=.TRUE.
          ENDIF

          CALL MPI_send(nodes_p0_temp,size1_MPI,MPI_Logical,i_MPI,                      &
                        i_MPI,MPI_COMM_WORLD,MPI_err)
        ENDDO ! i_MPI

        MPI_nodes_num=ii+1
        IF(.NOT. ALL(MPI_nodes_np(0:MPI_nodes_num-1)>=2))                              &
                                     write(out_unit,*) 'warning, only one core in node'

        write(out_unit,*) 'MPI working on',MPI_nodes_num,'nodes:'
        DO ii=0,MPI_nodes_num-1
          write(out_unit,11) MPI_nodes_name(ii),MPI_nodes_np(ii)
        ENDDO
        write(out_unit,*) ' '
11      format('   ',a10,' ',i3,' cores')

      ENDIF ! MPI_id

      CALL MPI_Bcast(MPI_nodes_num,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
      CALL MPI_Bcast(MPI_nodes_np(0:MPI_np-1),MPI_np,Int_MPI,root_MPI,                 &
                    MPI_COMM_WORLD,MPI_err)

      ! record node info
      jj=0
      DO ii=0,MPI_nodes_num-1
        IF(MPI_id<=MPI_nodes_np(0)-1) THEN
          MPI_node_id=0
          MPI_node_p0_id=0
        ELSEIF(MPI_id>Sum(MPI_nodes_np(0:ii))-1                                        &
              .AND. MPI_id<=Sum(MPI_nodes_np(0:ii+1))-1) THEN
          MPI_node_id=ii+1
          MPI_node_p0_id=Sum(MPI_nodes_np(0:ii))
        ENDIF
      ENDDO

    ENDIF ! MPI_fake_nodes

    IF(MPI_id==0 .OR. MPI_nodes_p0) THEN
      MPI_sub_id(1)=MPI_id+1
      MPI_sub_id(2)=MPI_id+MPI_nodes_np(MPI_node_id)-1
      write(out_unit,*) 'MPI_nodes_id check:',MPI_id,MPI_node_id,MPI_node_p0_id,MPI_sub_id(1),MPI_sub_id(2),MPI_nodes_p0
    ENDIF

  ENDSUBROUTINE get_nodes_info_MPI

!=======================================================================================

#else
  !< fake variables, for convenience
  Integer                                        :: MPI_err=0
  Integer                                        :: MPI_id=0
  Integer                                        :: MPI_np=1
  Integer                                        :: MPI_id_node=0
  Integer                                        :: MPI_np_node=1

#if(int8)
  Integer,parameter                              :: MPI_INTEGER_KIND=8
  Integer,parameter                              :: MPI_ADDRESS_KIND=8
#else
  Integer,parameter                              :: MPI_INTEGER_KIND=4
  Integer,parameter                              :: MPI_ADDRESS_KIND=4
#endif
  Integer                                        :: MPI_scheme=0
  Integer                                        :: iG1_MPI      
  Integer                                        :: iG2_MPI  
  Integer(kind=MPI_INTEGER_KIND)                 :: root_MPI=0      
  Integer(kind=MPI_INTEGER_KIND)                 :: size1_MPI=1
  Integer(kind=MPI_INTEGER_KIND)                 :: i_MPI    
  Logical                                        :: Srep_MPI=.FALSE.  
  Integer,allocatable                            :: iGs_MPI(:,:)  
  Integer,allocatable                            :: bounds_MPI(:,:) 
  Integer                                        :: Num_L2
  Integer                                        :: n_level2
  Integer                                        :: MPI_mc=1
  Integer                                        :: MPI_fake_nodes=0
  Logical                                        :: MPI_iGs_auto=.True.
  Logical                                        :: keep_MPI=.True.
  Integer                                        :: MPI_mem_node=0
  Logical                                        :: MPI_nodes_p0=.TRUE.
  Integer                                        :: MPI_sub_id(2)


  TYPE MPI_S_type 
    Logical                                      :: davidson=.FALSE.
    Logical                                      :: prop=.FALSE.
    Integer                                      :: num_resetH=1
    Integer                                      :: nb_channels=1
    Integer                                      :: nb_psi=1                  ! nb of levels to converge
  END TYPE 
  TYPE(MPI_S_type)                               :: MPI_S


  Integer                                        :: time_rate       !< for function system_clock()
  Integer                                        :: time_max        !< for function system_clock()

  Contains
  SUBROUTINE ini_MPI()
  ENDSUBROUTINE ini_MPI

  SUBROUTINE end_MPI()
  ENDSUBROUTINE end_MPI
#endif

END MODULE mod_MPI


