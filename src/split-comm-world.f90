subroutine split_comm_world(COMM_UNIVERSE4, COMM_WORLD4, MY_UNIVERSE_RANK, MY_WORLD_RANK, INSTANCE)

  implicit none

  include 'mpif.h'

  integer(4) :: MY_UNIVERSE_RANK, MY_WORLD_RANK ! ranks of the different communicators
  integer(4) :: COMM_UNIVERSE4, COMM_WORLD4     ! MPI communicators, UNIVERSE is used for inter-code communications
                                                ! WORLD is used for intra-code communications

  integer(4) :: UNIVERSE_SIZE4, WORLD_SIZE4 ! size of the different communicators
  integer(4) :: mapps                       ! number of different codes being executed
  integer(4) :: rank_master2=-1             ! rank of the process that will perform inter-code communications
  integer(4) :: istat4, ierro4
  integer(4) :: INSTANCE

  integer :: ii
  integer :: path_f, ext_i ! variables to obtain the code name

  character(128)          :: app_type         ! name of the code
  character(128), pointer :: app_type_arra(:) ! array with the name of the codes

  integer(4)          :: app_type_id         ! code identifier to perform the split
  integer(4), pointer :: app_type_id_arra(:) ! array with the identifiers of the codes

  !
  ! Nullify pointers & initializations
  !
  nullify(app_type_id_arra)
  nullify(app_type_arra)
  app_type_id = 0
  mapps       = 0

  !
  ! Get the names of the different applications launched
  ! Remove the path of the executable and the extension
  !
  !call GETARG(0_4,app_type)
  write(app_type, "(A8,I1)") "hemepure", INSTANCE

  !path_f = 0
  !ext_i  = len_trim(app_type)+1
  !do ii = 1, len_trim(app_type)
  !   if( app_type(ii:ii) == '/' )then
  !      path_f = ii
  !   end if
  !end do
  !do ii = path_f, len_trim(app_type)
  !   if( app_type(ii:ii) == '.' )then
  !      ext_i = ii
  !      exit
  !   end if
  !end do
  !
  !app_type = trim(app_type(path_f+1:ext_i-1))

  !
  ! Get the size and your rank in MPI_COMM_WORLD
  !
  COMM_UNIVERSE4 = MPI_COMM_WORLD

  call MPI_Comm_size(COMM_UNIVERSE4,UNIVERSE_SIZE4,istat4)
  call MPI_Comm_rank(COMM_UNIVERSE4,MY_UNIVERSE_RANK,istat4)

  if( UNIVERSE_SIZE4 > 1_4 ) then
     !
     ! Rank 0 gathers the names, compares them and assigns a code identifier (app_type_id) to each process
     !
     if( MY_UNIVERSE_RANK == 0 ) then
        allocate(app_type_arra(UNIVERSE_SIZE4))
        allocate(app_type_id_arra(UNIVERSE_SIZE4))
     else
        allocate(app_type_arra(1))
        allocate(app_type_id_arra(1))
     end if

     call MPI_Gather(app_type,128,MPI_CHARACTER,app_type_arra,128,MPI_CHARACTER,0_4,COMM_UNIVERSE4,istat4)

     if( MY_UNIVERSE_RANK == 0 ) then ! I am de MASTER OF THE UNIVERSE

        app_type_id      =  1
        app_type_id_arra = -1
        do ii = 1_4,UNIVERSE_SIZE4-1_4

           if( trim(app_type_arra(ii)) /= trim(app_type_arra(ii+1)) ) then
              app_type_id_arra(ii)   = app_type_id
              app_type_id_arra(ii+1) = app_type_id + 1
              app_type_id            = app_type_id + 1
              !
              ! This will be the process that will perform intercode communications
              !
              rank_master2           = ii
           else
              app_type_id_arra(ii)   = app_type_id
              app_type_id_arra(ii+1) = app_type_id
           endif
        end do

        mapps = app_type_id
     end if

     !
     ! Rank 0 in COMM_UNIVERSE sends the number of codes being executed to other ranks
     !
     call MPI_Bcast(mapps,1,MPI_INTEGER,0,COMM_UNIVERSE4,istat4)

     !
     ! Rank 0 scatters the code identifier (type_id_arra), and the split of COMM_UNIVERSE is performed
     !
     call MPI_Scatter(app_type_id_arra,1,MPI_INTEGER,app_type_id,1,MPI_INTEGER,0,COMM_UNIVERSE4,istat4)

     if( mapps > 1 ) then
        !
        ! Split communicator. This defines the communicator COMM_WORLD
        !
        call MPI_Comm_split(COMM_UNIVERSE4,app_type_id,MY_UNIVERSE_RANK,COMM_WORLD4,istat4)
        call MPI_Comm_size(COMM_WORLD4,WORLD_SIZE4,istat4)
        call MPI_Comm_rank(COMM_WORLD4,MY_WORLD_RANK,istat4)
     else
        COMM_WORLD4   = MPI_COMM_WORLD
        MY_WORLD_RANK = MY_UNIVERSE_RANK
     end if

     !
     ! Deallocate arrays
     !
     deallocate(app_type_id_arra)
     deallocate(app_type_arra)

  else
     COMM_WORLD4   = MPI_COMM_WORLD
     MY_WORLD_RANK = MY_UNIVERSE_RANK
  end if

  !print*, TRIM(app_type),' ',mapps, MY_UNIVERSE_RANK, MY_WORLD_RANK
  call MPI_Barrier(COMM_UNIVERSE4,istat4)
  !if( MY_UNIVERSE_RANK == 0 ) print*, "MASTER2 is ", rank_master2
end subroutine split_comm_world
