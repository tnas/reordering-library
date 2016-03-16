! COPYRIGHT (c) 2003 Council for the Central Laboratory
!         of the Research Councils
!
! Version 2.4.0
!
! History: See ChangeLog
!

   module hsl_mc73_double
       use hsl_mc65_double
       use hsl_zd11_double

! This module is for computing the
! spectral ordering a matrix A with a symmetric sparsity pattern.
! A multilevel algorithm is used. An approx. Fiedler vector of the
! (weighted) Laplacian associated with A is computed.
! The Fiedler vector is found by finding
! that of all the components and agglomerate.

      implicit none
      private
      public mc73_fiedler,mc73_order
      public mc73_print_message,mc73_initialize,mc73_control

      integer, parameter :: wpr = kind(0.0d0)

! Error/warning parameters
! memory alloc failure
      integer, public, parameter :: MC73_ERR_MEMORY_ALLOC = -1
! memory deallocate failure
      integer, public, parameter :: MC73_ERR_MEMORY_DEALLOC = -2
! N < 0
      integer, public, parameter :: MC73_ERR_N_NONPOSITIVE = -3
! LIRN < IP(N+1)-1
      integer, public, parameter :: MC73_ERR_LIRN_TOOSMALL = -4
! IP is not monotonically increasing
      integer, public, parameter :: MC73_ERR_RANGE_IP = -5
! JOB not valid
      integer, public, parameter :: MC73_ERR_JOB_WRONG = -6
! Array A not long enough
      integer, public, parameter :: MC73_ERR_A_TOOSMALL = -7

! IRN is not within the  range of [1,n]. Such entries are excluded
      integer, public, parameter :: MC73_WARN_RANGE_IRN = 1
! duplicate entries were found
      integer, public, parameter :: MC73_WARN_DUP_ENTRY = 2
! max number of Rayleigh Quotient iterations reached (not job = 1)
      integer, public, parameter :: MC73_WARN_MAXIT = 4
! MC67 did not reduce the profile
      integer, public, parameter :: MC73_WARN_MC67 = 8
! *****************************************************************
      type mc73_control
      integer :: lp ! Unit for error messages
      integer :: wp ! Unit for warning messages
      integer :: mp ! Unit for monitor output
      integer :: print_level ! Level of diagnostic printing
      integer :: mglevel ! Number of levels in the multilevel grid
      integer :: mlancz ! maximum number of Lanczos vectors used
      integer :: coarsest_size
! Stop coarsening when coarse graph has at most coarsest_size nodes
      integer :: maxit
! max. number of Rayleigh Quotient iterations (not job = 1)
      integer :: hager_exchange
! controls Hager exchanges
!            <=0 = no exchanges
!             >0 = down/up applied max. of hager_exchange times
      real (kind = wpr) :: tol  ! eigenvalue tolerance
      real (kind = wpr) :: tol1 ! RQI tolerance
      real (kind = wpr) :: rtol ! SYMMLQ tolerance
! minimum and maximum grid reduction factors
! that must be achieved during coarsening (not job = 1)

! If cgrid%size is greater than max_reduction*grid%size
! or cgrid%size is less than min_reduction*grid%size
! then carry on coarsening

      real (kind = wpr) :: min_reduction ! (default 0.1)
      real (kind = wpr) :: max_reduction ! (default 0.8)
      end type mc73_control
! *****************************************************************

      type mc73_multigrid
! size of this level (number of rows)
      integer :: size
! this level of matrix
      type (zd11_type), pointer :: graph
! where each row of this level of matrix will go (ie ordering for
! this level)
      integer, pointer, dimension (:) :: where
! row weight: the number of vertices this vertex of the coarse graph
! matrix represents
      integer, pointer, dimension (:) :: row_wgt
! the level
      integer :: level
! pointers to the fine and coarse grid (matrix)
      type (mc73_multigrid), pointer :: coarse, fine
! the prolongation operator
      type (zd11_type), pointer :: p

      end type mc73_multigrid
! *****************************************************************

      type list_node_type
      integer :: id
      type (list_node_type), pointer :: next,prev
      end type list_node_type
! *****************************************************************

! this is the first node in the bucket, any new node will be added
! before this node and become the first node
! thus this is a last-in-first-out scheme
      type list_node_pointer
      type (list_node_type), pointer :: current_node
      end type list_node_pointer
! *****************************************************************

      type queue_type
      integer :: low_bound,up_bound
! maxnodes = max no. of nodes allowed (as dimension of mynode)
! maxgain = max gain so far in the buckets (gives the highest non-empty bucket)
! num_nodes = total number of nodes in the buckets
      integer :: maxnodes,maxgain,num_nodes
! the starting node of each bucket (point to NULL if empty)
      type (list_node_pointer), pointer, dimension(:) :: buckets
! the location of each node in the buckets
      type (list_node_type), pointer, dimension (:) :: mynode
      end type queue_type
! *****************************************************************
      type mc73_multigrid_eig
! size of this level (number of rows)
      integer :: size
! this level of matrix
      type (zd11_type), pointer :: graph
! the eigenvector and eigenvalue
      real (kind = wpr), pointer, dimension (:) :: eigvector
      real (kind = wpr) :: lambda
! row weight: the number of vertices this vertex of the coarse graph
! matrix represents
      integer, pointer, dimension (:) :: row_wgt
! the level
      integer :: level
! pointers to the fine and coarse grid (matrix)
      type (mc73_multigrid_eig), pointer :: coarse, fine
! the prolongation operator
      type (zd11_type), pointer :: p

      end type mc73_multigrid_eig
! *****************************************************************
      interface mc73_fiedler
        module procedure mc73_fiedler
      end interface

      interface mc73_order
        module procedure mc73_order
      end interface

      interface mc73_print_message
        module procedure print_message
      end interface

      contains

      subroutine mc73_initialize(control)
      type(mc73_control), intent(out) :: control
      control%lp = 6
          control%wp = 6
          control%mp = 6
          control%print_level = 0
          control%mglevel = 10
          control%mlancz = 300
          control%coarsest_size = 200
          control%maxit = 10
          control%hager_exchange = 0

          control%tol = 0.001_wpr
          control%tol1 = 0.001_wpr
          control%rtol = 0.01_wpr
          control%min_reduction = 0.1_wpr
          control%max_reduction = 0.8_wpr
      end subroutine mc73_initialize
! *****************************************************************
      subroutine mc73_fiedler(n,lirn,irn,ip,list,fvector,info,a) BIND(C, NAME='mc73_fiedler')

      real (kind = wpr), parameter :: zero = 0.0_wpr
! the matrix
! n = order of A
! lirn = length of irn (at least number of entries in lower triangular
!        part of A)
! irn holds row indices
! ip holds column pointer indices
      integer,  intent (in) :: n
      integer,  intent (in) :: lirn
      integer,  intent (in) :: irn(lirn)
      integer,  intent (inout) :: ip(n+1)
! added following optional arguement.
! If present must hold matrix weights (otherwise, standard
! Laplacian used)
      real (kind = wpr),  intent (in), dimension(:), optional :: a

! information flag and control arrays
      integer, intent (out) :: info(1:10)
      !type(mc73_control),intent(in) :: control

! On exit, the number of the component within the
! adjacency graph of A that variable i belongs to is list(i)
      integer, intent (out) :: list(n)

! Fiedler vector to be calculated
      real (kind = wpr), intent (out) :: fvector(n)

! the matrix in its internal form
      type (zd11_type) :: matrix

! llap is set to true if a is present (weighted Laplacian)
      logical :: llap
! local error flag
      integer   :: ierr,i
! local control and info. arrays
      real (kind = wpr)  :: cntl(10)
      integer :: icntl(10)
! print level and error/warning streams
      integer :: err,wrn,mp,print_level
      logical :: lerr,lwrn

      icntl(1) = 6!control%lp
      icntl(2) = 6!control%wp
      icntl(3) = 6!control%mp
      icntl(4) = 0!control%print_level
      icntl(5) = max(1,10)!max(1,control%mglevel)
      icntl(6) = 200!control%coarsest_size
      icntl(7) = 10!control%maxit
      icntl(8) = 300!control%mlancz

      cntl(1) = 0.001_wpr!control%tol
      cntl(2) = max(0.1_wpr,0.1_wpr)!max(0.1_wpr,control%min_reduction)
      cntl(3) = 0.8_wpr!control%max_reduction
      if (cntl(2) > cntl(3)) cntl(2) = 0.1_wpr
      cntl(4) = 0.001_wpr!control%tol1
      cntl(5) = min(0.01_wpr*0.001_wpr,10.0_wpr*0.001_wpr)!min(control%rtol*0.001_wpr,10.0_wpr*control%tol)

      if (icntl(4) < 0) print_level = 0
! The default is icntl(4) = 0
      if (icntl(4) == 0) print_level = 1
      if (icntl(4) == 1) print_level = 2
      if (icntl(4) > 1) print_level = 3
      mp = icntl(3)
      if (mp < 0) print_level = 0
! Set error controls
      lerr = icntl(1) .ge. 0 .and. print_level.gt.0
      err  = icntl(1)
! Set warning controls
      lwrn = icntl(2).ge.0 .and. print_level.gt.0
      wrn = icntl(2)
      
! alteracao do vetor ip
      do i=1,n+1
        ip(i) = ip(i) + 1
      end do
      
      info(1:10) = 0
! initial printing
      if (print_level >= 2) then
         write (mp,'(/a)')    ' MC73_FIEDLER called'
         write (mp,'(a,i8,a,i10)') ' On entry N = ',N,'  LIRN = ',LIRN
         if (icntl(5) >1) then
           write (mp,'(a,i8)') ' Max. number of grid levels = ',icntl(5)
           write (mp,'(a,i8)') ' Max. size of coarsest grid = ',icntl(6)
           write (mp,'(a,i8)')&
                ' Max. no. of Rayleigh Quotient iterations = ',icntl(7)
           write (mp,'(a,i8)')&
                ' Max. no. of Lanczos vectors used = ',icntl(8)
           write (mp,'(a,es12.3)') &
         ' Eigenvalue tolerance          = ',cntl(1)
           write (mp,'(a,es12.3)') &
         ' RQI convergence tolerance     = ',cntl(4)
           write (mp,'(a,es12.3)') &
         ' SYMMLQ convergence tolerance  = ',cntl(5)
           write (mp,'(a,2es12.3)') &
         ' Grid reduction parameters     = ',cntl(2:3)
          else
           write (mp,'(a,i8)') ' Single level used'
           write (mp,'(a,es12.3)') &
         ' Eigenvalue tolerance          = ',cntl(1)
           write (mp,'(a,es12.3)') &
         ' RQI convergence tolerance     = ',cntl(4)
           write (mp,'(a,es12.3)') &
         ' SYMMLQ convergence tolerance  = ',cntl(5)
          end if
          if (present(a)) &
            write (mp,'(a)') ' Weighted Laplacian used'
      end if

! for graph of single vertex we can return immediately
     if (n == 1) then
        list(1) = 0
        fvector(1) = zero
        return
     end if

! Check input parameters
      if (n <= 0) then
         info(1) = MC73_ERR_N_NONPOSITIVE
      else if (lirn < ip(n+1)-1) then
         info(1) = MC73_ERR_LIRN_TOOSMALL
      else if (present(a)) then
         if (size(a) < ip(n+1)-1) then
            info(1) = MC73_ERR_A_TOOSMALL
         end if
      end if
      if (info(1) < 0) then
         if (lerr) call print_message(info(1),err,' MC73_FIEDLER')
         return
      end if
   
     do i = 1,n
        if (ip(i) > ip(i+1)) then
           ierr = -3
           go to 998
        end if
     end do

! convert to zd11 format and check data
      if (present(a)) then
        call mc65_matrix_construct(matrix,n,n,ip,irn,ierr,&
                                   val=a,checking=1)
! Ensure we have positive weights
        do i = 1,ip(n+1)-1
           matrix%val(i) = abs(matrix%val(i))
        end do
      else
        call mc65_matrix_construct(matrix,n,n,ip,irn,ierr,&
                                   type='general',checking=1)
      end if

! check for error return
      if (ierr < 0) then
         if (ierr == -12) ierr = -3
         go to 998
      end if
      if (ierr == 1) then
         info(1) = MC73_WARN_RANGE_IRN
         if (lwrn) call print_message(info(1),wrn,' MC73_FIEDLER')
      else if (ierr == 5) then
         info(1) = MC73_WARN_DUP_ENTRY
         if (lwrn) call print_message(info(1),wrn,' MC73_FIEDLER')
      end if

! make the matrix symmetric and remove diagonal
! (if a not present, graph of matrix A+A^T is returned)
      if (present(a)) then
         call mc65_matrix_symmetrize(matrix,ierr)
         call mc65_matrix_remove_diagonal(matrix)
         llap = .true.
      else
         call mc65_matrix_symmetrize(matrix,ierr,graph=.true.)
         llap = .false.
      end if
      if (ierr /= 0 ) go to 998
! Ordering depends how graph is coarsened. we can get consistent
! results by preordering (sorting)
!      call mc65_matrix_sort(matrix,ierr)
      if (ierr /= 0 ) go to 998

! Ready to compute Fielder vector

      call fiedler_graph(.false.,n,matrix,list,fvector,info,icntl,  &
                          cntl,llap)
      if (info(1) < 0) go to 999
      go to 999

! error returns
998   if (ierr == -1) then
         info(1) = MC73_ERR_MEMORY_ALLOC
      else if (ierr == -2) then
         info(1) = MC73_ERR_MEMORY_DEALLOC
      else if (ierr == -3) then
         info(1) = MC73_ERR_RANGE_IP
      else
         info(1) = ierr
      end if
      if (info(1) < 0 .and. lerr) &
            call print_message(info(1),err,' MC73_FIEDLER')
      call mc65_matrix_destruct(matrix,ierr)
! do not test error message since matrix may not have been allocated
      return
      
999   if (print_level >= 2 .and. info(1) >= 0)   &
         write (mp,fmt='(/a,i3,4i8)') ' On exit INFO(1:2) = ',info(1:2)
! alteracao do vetor ip
      do i=1,n+1
        ip(i) = ip(i) - 1
      end do
      
      call mc65_matrix_destruct(matrix,ierr)      
      return

      end subroutine mc73_fiedler
! *****************************************************************
      subroutine mc73_order(job,n,lirn,irn,ip,perm,control,info,rinfo,a)

      real (kind = wpr), parameter :: zero = 0.0_wpr
! job controls which algorithm is employed
      integer,  intent (in) :: job
! the matrix
! n = order of A
! lirn = length of irn (at least number of entries in lower triangular
!        part of A)
! irn holds row indices
! ip holds column pointer indices
      integer,  intent (in) :: n
      integer,  intent (in) :: lirn
      integer,  intent (in) :: irn(lirn)
      integer,  intent (in) :: ip(n+1)
! If a is present, must hold the matrix weights
      real (kind = wpr),  intent (in), dimension(:), optional :: a

! information flag and control arrays
      real (kind = wpr), intent (out)  :: rinfo(1:20)
      integer, intent (out) :: info(1:10)
      type(mc73_control),intent(in) :: control

! permutation
      integer, intent (out) :: perm(n)
! Fiedler vector to be calculated

! the matrix in its internal form
      type (zd11_type) :: matrix

! local work arrays
      real (kind = wpr), dimension (:), allocatable :: fvector
      real (kind = wpr), dimension (:), allocatable :: w
      integer, dimension (:), allocatable :: iw
! llap is set to true if a is present (weighted Laplacian)
      logical :: llap
! local error flag
      integer   :: ierr,i,hager
      integer   :: liw,nz
      integer   :: pair,temp_perm,temp_perm1
! local control and info. arrays
      real (kind = wpr)  :: cntl(10),cntl61(5),rinfo61(15),rinfo67(10)
      real (kind = wpr)  :: weight(2),rtemp(4),rtemp1(4)
      integer :: icntl(10),icntl61(10),icntl67(10),jcntl(2)
      integer :: info60(4),info61(10),info67(10)
! print level and error/warning streams
      integer :: err,wrn,mp,print_level
      logical :: lerr,lwrn

      external mc60cd,mc60fd,mc61id,mc61ad,mc67id,mc67ad

      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%print_level
      icntl(5) = max(1,control%mglevel)
      icntl(6) = control%coarsest_size
      icntl(7) = control%maxit
      icntl(8) = control%mlancz

      cntl(1) = control%tol
      cntl(2) = max(0.1_wpr,control%min_reduction)
      cntl(3) = control%max_reduction
      if (cntl(2) > cntl(3)) cntl(2) = 0.1_wpr
      cntl(4) = control%tol1
      cntl(5) = min(control%rtol,10.0_wpr*control%tol)

      if (icntl(4) < 0) print_level = 0
! The default is icntl(4) = 0
      if (icntl(4) == 0) print_level = 1
      if (icntl(4) == 1) print_level = 2
      if (icntl(4) > 1) print_level = 3
      mp = icntl(3)
      if (mp < 0) print_level = 0
! Set error controls
      lerr = icntl(1).ge.0 .and. print_level.gt.0
      err  = icntl(1)
! Set warning controls
      lwrn = icntl(2).ge.0 .and. print_level.gt.0
      wrn = icntl(2)

      info(1:10) = 0
      rinfo(1:20) = zero
      rinfo67(1:10) = zero

      hager = control%hager_exchange
      if (hager .le. 0) hager = 0
      if (present(a)) hager = 0

! initial printing (depends on job)
      if (print_level >= 2) then
         write (mp,'(/a,i4)')  ' MC73_ORDER called with JOB = ',JOB
         write (mp,'(a,i8,a,i10)') ' On entry N = ',N,'  LIRN = ',LIRN
         if (job == 1) then
          if (icntl(5) >1) then
           write (mp,'(a,i8)') ' Max. number of grid levels = ',icntl(5)
           write (mp,'(a,i8)') ' Max. size of coarsest grid = ',icntl(6)
           write (mp,'(a,2es12.3)')' Grid reduction parameters  = ',cntl(2:3)
          else
           write (mp,'(a,i8)') ' Single level. MC61 is used.'
          end if
         else if (job <= 3) then
          if (icntl(5) >1) then
           write (mp,'(a,i8)') ' Max. number of grid levels = ',icntl(5)
           write (mp,'(a,i8)') ' Max. size of coarsest grid = ',icntl(6)
           write (mp,'(a,i8)')&
                ' Max. no. of Rayleigh Quotient iterations = ',icntl(7)
           write (mp,'(a,i8)')&
                ' Max. no. of Lanczos vectors used = ',icntl(8)
           write (mp,'(a,es12.3)') &
         ' Eigenvalue tolerance          = ',cntl(1)
           write (mp,'(a,es12.3)') &
         ' RQI convergence tolerance     = ',cntl(4)
           write (mp,'(a,es12.3)') &
         ' SYMMLQ convergence tolerance  = ',cntl(5)
           write (mp,'(a,2es12.3)') &
         ' Grid reduction parameters     = ',cntl(2:3)
          else
           write (mp,'(a,i8)') ' Single level used'
           write (mp,'(a,es12.3)') &
         ' Eigenvalue tolerance          = ',cntl(1)
           write (mp,'(a,es12.3)') &
         ' RQI convergence tolerance     = ',cntl(4)
           write (mp,'(a,es12.3)') &
         ' SYMMLQ convergence tolerance  = ',cntl(5)
          end if
         end if
         if (present(a)) &
            write (mp,'(a)') ' Weighted Laplacian used'
         if (hager > 0) &
            write (mp,'(a)') ' MC67 (Hager) used to refine ordering.'
      end if

! for graph of single vertex we can return immediately
     if (n == 1) then
        perm(1) = 1
        return
     end if

! Check input parameters
      if (n <= 0) then
         info(1) = MC73_ERR_N_NONPOSITIVE
      else if (lirn < ip(n+1)-1) then
         info(1) = MC73_ERR_LIRN_TOOSMALL
      else if (job < 1 .or. job > 3) then
         info(1) = MC73_ERR_JOB_WRONG
      else if (present(a) .and. job /= 1) then
         if (size(a) < ip(n+1)-1) then
            info(1) = MC73_ERR_A_TOOSMALL
         end if
      end if
      if (info(1) < 0) then
         if (lerr) call print_message(info(1),err,' MC73_ORDER')
         return
      end if

     do i = 1,n
        if (ip(i) > ip(i+1)) then
           ierr = -3
           go to 998
        end if
     end do

! convert to zd11 format
! Include checking of data
      if (present(a) .and. job /= 1) then
        call mc65_matrix_construct(matrix,n,n,ip,irn,ierr,&
                                   val=a,checking=1)
! Ensure we have positive weights
        do i = 1,ip(n+1)-1
           matrix%val(i) = abs(matrix%val(i))
        end do
      else
        call mc65_matrix_construct(matrix,n,n,ip,irn,ierr,&
                                   type='general',checking=1)
      end if
! check for error return
      if (ierr < 0) then
         if (ierr == -12) ierr = -3
         go to 998
      end if
      if (ierr == 1) then
         info(1) = MC73_WARN_RANGE_IRN
         if (lwrn) call print_message(info(1),wrn,' MC73_ORDER')
      else if (ierr == 5) then
         info(1) = MC73_WARN_DUP_ENTRY
         if (lwrn) call print_message(info(1),wrn,' MC73_ORDER')
      end if

! make the matrix symmetric and remove diagonal
      if (present(a) .and. job /= 1) then
         call mc65_matrix_symmetrize(matrix,ierr)
         call mc65_matrix_remove_diagonal(matrix)
         llap = .true.
      else
         call mc65_matrix_symmetrize(matrix,ierr,graph=.true.)
         llap = .false.
      end if
      if (ierr /= 0 ) go to 998
! Ordering depends how graph is coarsened. we can get consistent
! results by preordering (sorting)
!      call mc65_matrix_sort(matrix,ierr)
      if (ierr /= 0 ) go to 998


     if (job == 1 .and. icntl(5) == 1) then
! Single level ==> multilevel reduces to Sloan algorithm
! Call MC61 and return
        call mc61id(icntl61,cntl61)
! Switch off printing
        icntl61(1) = -1
        icntl61(2) = -1
! Continue if duplicates/out-of-range entries found
        icntl61(3) = 1
! we are not using supervariables (as we want irn,ip to be for
! original matrix, not condensed matrix)
        icntl61(4) = 1
        liw = 8*n + 2
        allocate(iw(liw),w(n),stat=ierr)
        if (ierr /= 0) then
           ierr = -1
           go to 998
        end if
        nz = matrix%ptr(n+1) - 1
        call mc61ad(1,n,nz,matrix%col,matrix%ptr,perm,liw,iw,w,&
                    icntl61,cntl61,info61,rinfo61)
        info(2:5) = info61(4:7)
        rinfo(1:8) = rinfo61(1:8)

! Optionally call Hager to refine order
        if (hager > 0) then
          liw = size(iw)
          call mc67id(icntl67)
          icntl67(1) = -1
          icntl67(2) = -1
          icntl67(3) = 0
          icntl67(4) = hager
! down/up exchanges
          icntl67(5) = 0
          icntl67(6) = 0
          call mc67ad(n,nz,matrix%col,matrix%ptr,perm,  &
                      liw,iw,icntl67,info67,rinfo67)
          if (rinfo67(1) > zero) then
! mc67 has reduced profile
             iw(1:n) = 1
             call mc60fd(n,n,nz,matrix%col,matrix%ptr,iw,perm,   &
                         iw(n+1),rinfo(5))
          end if
        end if

        deallocate(iw,w,stat=ierr)
        if (ierr /= 0) then
           ierr = -2
           go to 998
        end if

        if (print_level >= 2) then
           write (mp,fmt='(/a)')  ' Original ordering:'
           call write_rinfo(n,mp,rinfo)
           if (hager == 0) write (mp,fmt='(/a)')  ' MC61 ordering:'
           if (hager > 0) write (mp,fmt='(/a)') &
                        ' MC61 + Hager ordering:'
           call write_rinfo(n,mp,rinfo(5:8))
        end if
        if (rinfo67(1) == zero .and. hager > 0) then
! mc67 has not reduced profile
          info(1) = info(1) + MC73_WARN_MC67
          if (lwrn) call print_message(MC73_WARN_MC67,wrn,' MC73')
        end if

        go to 999

      end if

! Ready to compute ordering  ... algorithm depends on parameter job
      if (job == 1) then
! Multilevel algorithm
         call front_matrix(n,matrix,perm,icntl,cntl,info)
         if (info(1) < 0) go to 999
! Compute statistics for initial and final order
         liw = 4*n+1
         if (hager > 0) liw = 7*n+2
         allocate(iw(liw),stat=ierr)
         if (ierr /= 0) then
            ierr = -1
            go to 998
         end if
         iw(1:n) = 1
! statistics for initial ordering
         do i = 1,n
            iw(n+i) = i
         end do
         nz = matrix%ptr(n+1)-1
         call mc60fd(n,n,nz,matrix%col,matrix%ptr,  &
                     iw(1),iw(n+1),iw(2*n+1),rinfo)
         if (print_level >= 2) then
           write (mp,fmt='(/a)')  ' Original ordering:'
           call write_rinfo(n,mp,rinfo)
         end if

! Optionally call Hager to refine order (default setting used)
         if (hager > 0) then
           liw = size(iw)
           call mc67id(icntl67)
           icntl67(1) = -1
           icntl67(2) = -1
           icntl67(3) = 0
           icntl67(4) = hager
           icntl67(5) = 0
           icntl67(6) = 0
           call mc67ad(n,nz,matrix%col,matrix%ptr,perm,liw,iw,   &
                       icntl67,info67,rinfo67)
         end if

! statistics for final ordering
         iw(1:n) = 1
         call mc60fd(n,n,nz,matrix%col,matrix%ptr,iw(1),perm,   &
                     iw(n+1),rinfo(5))
         if (print_level >= 2) then
           if (hager == 0) write (mp,fmt='(/a)')  ' Multilevel ordering:'
           if (hager > 0) write (mp,fmt='(/a)') &
                        ' Multilevel + Hager ordering:'
            call write_rinfo(n,mp,rinfo(5:8))
         end if
         if (rinfo67(1) == zero .and. hager > 0) then
! mc67 has not reduced profile
           info(1) = info(1) + MC73_WARN_MC67
           if (lwrn) call print_message(MC73_WARN_MC67,wrn,' MC73')
         end if


         deallocate(iw,stat=ierr)
         if (ierr /= 0) then
            ierr = -2
            go to 998
         end if


      else

! job = 2 or 3. spectral ordering (possibly refined with Sloan)
! we need to allocate an array of length n
         allocate(fvector(n),stat=ierr)
         if (ierr /= 0) then
            ierr = -1
            go to 998
         end if

         call fiedler_graph(.true.,n,matrix,perm,fvector,info,icntl,  &
                             cntl,llap)

         deallocate(fvector,stat=ierr)
         if (info(1) < 0) go to 999
         if (ierr /= 0) then
            ierr = -2
            go to 998
         end if

         liw = 7*n+2
         if (job == 2 .and. hager == 0) liw = 4*n+1
         if (job == 3) liw = liw + n
         allocate(iw(liw),stat=ierr)
         if (ierr /= 0) then
            ierr = -1
            go to 998
         end if
         iw(1:n) = 1
! statistics for initial ordering
         do i = 1,n
            iw(n+i) = i
         end do
         nz = matrix%ptr(n+1)-1
         call mc60fd(n,n,nz,matrix%col,matrix%ptr,   &
                     iw(1),iw(n+1),iw(2*n+1),rinfo)
         if (print_level >= 2) then
           write (mp,fmt='(/a)')  ' Original ordering:'
           call write_rinfo(n,mp,rinfo)
         end if
! statistics for spectral ordering
         call mc60fd(n,n,nz,matrix%col,matrix%ptr,   &
                     iw(1),perm,iw(n+1),rinfo(5))

! if asked for by user, compute hybrid ordering
         if (job == 3) then
            rinfo(9:12) = rinfo(5:8)
            pair = 1 + n + 3*n + 1
            temp_perm = pair + 2*n
            temp_perm1 = temp_perm + n
            allocate(w(n),stat=ierr)
            if (ierr /= 0) then
               ierr = -1
               go to 998
            end if

            if (rinfo(5) > rinfo(1)) then
! spectral ordering worse than original ordering and so use
! original ordering in hybrid algorithm
! (ordering chosen on basis of profile)
               do i = 1,n
                  perm(i) = i
               end do
            end if

            jcntl(1) = 0
            jcntl(2) = 2
! Do 2 runs (different weights)
            weight(1) = 1.0_wpr
            weight(2) = 2.0_wpr
            iw(temp_perm:temp_perm+n-1) = perm(1:n)
            iw(temp_perm1:temp_perm1+n-1) = perm(1:n)
            call mc60cd(n,n,nz,matrix%col,matrix%ptr,iw(1),jcntl, &
                        iw(temp_perm),weight,iw(pair),info60,iw(n+1),w)
            call mc60fd(n,n,nz,matrix%col,matrix%ptr,   &
                        iw(1),iw(temp_perm),iw(n+1),rtemp)

            weight(1) = 16.0_wpr
            weight(2) = 1.0_wpr
            call mc60cd(n,n,nz,matrix%col,matrix%ptr,iw(1),jcntl, &
                        iw(temp_perm1),weight,iw(pair),info60,iw(n+1),w)
            call mc60fd(n,n,nz,matrix%col,matrix%ptr,iw(1), &
                        iw(temp_perm1),iw(n+1),rtemp1)
! If this is better than first set of weights, keep this ordering
            if (rtemp1(1) < rtemp(1)) then
               if (rtemp1(1) < rinfo(5)) then
! we have improved on the spectral ordering
                  rinfo(5:8) = rtemp1(1:4)
                  perm(1:n) = iw(temp_perm1:temp_perm1+n-1)
               end if
            else
               if (rtemp(1) < rinfo(5)) then
! we have improved on the spectral ordering
                  rinfo(5:8) = rtemp(1:4)
                  perm(1:n) = iw(temp_perm:temp_perm+n-1)
               end if
            end if
            deallocate(w,stat=ierr)
            if (ierr /= 0) then
               ierr = -2
               go to 998
            end if

         end if

! Optionally call Hager to refine order (default setting used)
         if (hager > 0) then
           liw = size(iw)
           call mc67id(icntl67)
           icntl67(1) = -1
           icntl67(2) = -1
           icntl67(3) = 0
           icntl67(4) = hager
           icntl67(5) = 0
           icntl67(6) = 0
           call mc67ad(n,nz,matrix%col,matrix%ptr,   &
                       perm,liw,iw,icntl67,info67,rinfo67)

           if (rinfo67(1) > zero) then
! mc67 has reduced profile
! Compute statistics for final ordering
             iw(1:n) = 1
             call mc60fd(n,n,nz,matrix%col,matrix%ptr,  &
                         iw(1),perm,iw(n+1),rinfo(5))
           end if
         end if

         deallocate(iw,stat=ierr)
         if (ierr /= 0) then
            ierr = -2
            go to 998
         end if

         if (print_level >= 2) then
           if (hager == 0) then
             if (job == 2) then
               write (mp,fmt='(/a)')  ' Spectral ordering:'
             else if (job == 3) then
               write (mp,fmt='(/a)')  ' Spectral ordering:'
               call write_rinfo(n,mp,rinfo(9:12))
               write (mp,fmt='(/a)')  ' Hybrid ordering:'
             end if
           else
             if (job == 2) then
               write (mp,fmt='(/a)')  ' Spectral + Hager ordering:'
             else if (job == 3) then
               write (mp,fmt='(/a)')  ' Hybrid + Hager ordering:'
             end if
           end if
           call write_rinfo(n,mp,rinfo(5:8))
         end if
         if (rinfo67(1) == zero .and. hager > 0) then
! mc67 has not reduced profile
           info(1) = info(1) + MC73_WARN_MC67
           if (lwrn) call print_message(MC73_WARN_MC67,wrn,' MC73')
         end if

      end if
      go to 999

! error returns
998   if (ierr == -1) then
         info(1) = MC73_ERR_MEMORY_ALLOC
      else if (ierr == -2) then
         info(1) = MC73_ERR_MEMORY_DEALLOC
      else if (ierr == -3) then
         info(1) = MC73_ERR_RANGE_IP
      else
         info(1) = ierr
      end if
      if (info(1) < 0 .and. lerr) &
            call print_message(info(1),err,' MC73_ORDER')
      call mc65_matrix_destruct(matrix,ierr)
! do not test error message since matrix may not have been allocated
      return

999   if (print_level >= 2 .and. info(1) >= 0)   &
         write (mp,fmt='(/a,i3,4i8)') ' On exit INFO(1:2) = ',info(1:2)
      call mc65_matrix_destruct(matrix,ierr)
      return


      end subroutine mc73_order
! *****************************************************************
      subroutine write_rinfo(n,mp,rinfo)
      integer n,mp
      real (kind=wpr) rinfo(4)

      write (mp,fmt='(a,i12/a,i12/a,i12/a,i12/a,es12.4)')  &
           ' maximum wavefront = ',int(rinfo(2)),  &
           ' bandwidth         = ',int(rinfo(3)),  &
           ' profile           = ',int(rinfo(1)),  &
           ' envelope size     = ',int(rinfo(1)) - n,  &
           ' r.m.s. wavefront  = ',rinfo(4)
      end subroutine write_rinfo
! *****************************************************************

      subroutine front_matrix(n,matrix,perm,icntl,cntl,info)
! matrix: the matrix to be reordered.
! n: row dimension of the matrix
      integer,  intent (in) :: n
      type (zd11_type), target :: matrix
      integer,  intent (inout) :: info(10)
      integer,  intent(in) :: icntl(10)
      real (kind = wpr), intent (in) :: cntl(10)

! ordering to be calculated
      integer :: perm(:)

! the multilevel of graphs (matrices)
      type (mc73_multigrid) :: grid

! the number of levels, default 10
      integer :: i

! ============== used for tracking components of the graph
! root: the root of the current component
      integer :: root

! comp_nvtx: number of vertices in this component
! comp_nz: number of edges*2 in the component
! ncomp: number of components
      integer :: comp_nvtx,comp_nz,ncomp

! mask: mask that is initialised to zero. After a
! vertex is visited, its mask is set to the new index
! of this vertex in the component
! comp_vtx_list: a list of the vertices in the current component
! mask and comp_vtx_list satisfies mask(comp_vtx(i)) = i,
! i = 1,...,comp_nvtx
      integer, dimension (:), pointer :: mask, comp_vtx_list

! submatrix: the subgraph of this component
      type (zd11_type), pointer :: submatrix

! number of levels in the multilevel grid (default 10)
      integer :: mglevel
! print level
      integer :: err,wrn,mp,print_level
      logical :: lerr,lwrn

      integer degree,k,ierr,st
! number of nodes that have been ordered
      integer :: lstnum
      integer :: mglevel_cur

      if (icntl(4) < 0) print_level = 0
! The default is icntl(4) = 0
      if (icntl(4) == 0) print_level = 1
      if (icntl(4) == 1) print_level = 2
      if (icntl(4) > 1) print_level = 3
      mp = icntl(3)
      if (mp < 0) print_level = 0
! Set error controls
      lerr = icntl(1).ge.0 .and. print_level.gt.0
      err  = icntl(1)
! Set warning controls
      lwrn = icntl(2).ge.0 .and. print_level.gt.0
      wrn = icntl(2)

      mglevel = icntl(5)

! work on each connected component
      allocate(mask(n),stat=st)
      if (st /= 0) info(1) = MC73_ERR_MEMORY_ALLOC
      allocate(comp_vtx_list(n),stat=st)
      if (st /= 0) info(1) = MC73_ERR_MEMORY_ALLOC
      if (info(1) < 0) then
         if (lerr) call print_message(info(1),err,' MC73_ORDER')
         return
      end if

! ncomp is number of non-trivial components
      ncomp = 0
      mask = 0
! Order all nodes of degree zero
      lstnum = 0
      do i = 1,n
         k = matrix%ptr(i)
         degree = matrix%ptr(i+1) - k
         if (degree == 0) then
! single node in component ... number it first
           lstnum = lstnum + 1
           perm(i) = lstnum
           mask(i) = 1
         end if
      end do
! If matrix is diagonal, then we have now finished since all nodes
! have degree zero.
      if (lstnum == n) then
         deallocate(mask,comp_vtx_list,stat=st)
         if (st /= 0) then
           info(1) = MC73_ERR_MEMORY_DEALLOC
           if (lerr) call print_message(info(1),err,' MC73_ORDER')
         end if
         if (print_level == 3) write (mp,'(a//)') 'finished'
         return
      end if

      do i = 1,n
         if (mask(i) == 0) then
            root = i
            call mx_component(matrix,root,mask,comp_nvtx,comp_nz,&
                              comp_vtx_list,ierr)
! ierr holds stat parameter
            if (ierr /= 0 ) then
               info(1) = ierr
               if (lerr) call print_message(info(1),err,' MC73_ORDER')
               return
            end if

            ncomp = ncomp + 1
! store number of non-trivial components
            info(2) = ncomp
            if (print_level == 3) then
              write (mp,'(a,i4,a,i10,a,i10)') 'Component ',ncomp,&
                  ' number of nodes = ',comp_nvtx
            end if

            if (comp_nvtx == n) then
! if this is a connected graph anyway
                submatrix => matrix
            else
               allocate(submatrix)
               call mx_crop(matrix,comp_nvtx,comp_vtx_list,mask,submatrix,&
                            ierr,comp_nz)
               if (ierr < 0) then
                  info(1) = ierr
                  if (lerr) call print_message(info(1),err,' MC73_ORDER')
                  return
               end if
           end if

! construct the multrid at this level
             grid%size = submatrix%n
             grid%level = 1
             grid%graph => submatrix
             nullify(grid%p)

             allocate(grid%where(grid%size),grid%row_wgt(grid%size),stat=st)
             if (st /= 0) info(1) = MC73_ERR_MEMORY_ALLOC
             if (info(1) < 0) then
                if (lerr) call print_message(info(1),err,' MC73_ORDER')
                return
             end if

             grid%row_wgt = 1
! maximum level of grids allowed for this bisection
             mglevel_cur = mglevel

             call multilevel(grid,icntl,cntl,mglevel_cur,ierr)
             if (ierr /= 0) then
               if (ierr == -1) info(1) = MC73_ERR_MEMORY_ALLOC
               if (ierr == -2) info(1) = MC73_ERR_MEMORY_DEALLOC
               if (lwrn) call print_message(info(1),wrn,' MC73ORDER')
               return
             end if

! the ordering for this components will start from lstnum+1
             if (comp_nvtx == n) then
                perm = grid%where
             else
                perm(comp_vtx_list(1:comp_nvtx)) = grid%where + lstnum
             end if
! deallocate the finest level
             call multigrid_deallocate_first(comp_nvtx,n,grid,ierr)
             if (ierr /= 0) then
               info(1) = ierr
               if (lerr) call print_message(info(1),err,' MC73_ORDER')
               return
             end if
             lstnum = lstnum + comp_nvtx

          end if
       end do

      deallocate(mask,comp_vtx_list,stat=st)
      if (st /= 0) info(1) = MC73_ERR_MEMORY_DEALLOC
      if (info(1) < 0) then
         if (lerr) call print_message(info(1),err,' MC73ORDER')
         return
      end if

      end subroutine front_matrix

! ********************************************************
      recursive subroutine multilevel(grid,icntl,cntl,mglevel_cur,info)
      real (kind = wpr), parameter :: half = 0.5_wpr
      real (kind = wpr), parameter :: one = 1.0_wpr
! this level of matrix (grid)
      type (mc73_multigrid), intent (inout), target :: grid
! the coarse level grid
      type (mc73_multigrid), pointer :: cgrid
      integer,  intent(in) :: icntl(10)
      real (kind = wpr), intent (in) :: cntl(10)
      integer :: mglevel_cur
! Error flag
      integer :: info

! the coarse grid prolongator
      type (zd11_type), pointer :: p
! the partition on fine and coarse grid
      integer, dimension (:), pointer :: fwhere, cwhere
! the coarse graph
      type (zd11_type), pointer :: cgraph
! the fine graph
      type (zd11_type), pointer :: graph
! fine and coarse graph vertex weight
      integer, dimension (:), pointer :: crow_wgt

! cnvtx: number of vertex (rows) in the coarse matrix
      integer :: cnvtx

      integer :: st
      integer :: lpermsv

! Grid reduction factor
      real (kind = wpr) :: grid_rdc_fac_min,grid_rdc_fac_max
! smallest problem size to stop coarsening (default 100)
      integer :: coarsest_size
      integer   :: mp, print_level

! priority value of the fine and coarse grid vertices
      real (kind = wpr), dimension (:), pointer :: cpriority,fpriority

      info = 0
      if (icntl(4) < 0) print_level = 0
! The default is icntl(4) = 0
      if (icntl(4) == 0) print_level = 1
      if (icntl(4) == 1) print_level = 2
      if (icntl(4) > 1) print_level = 3
      mp = icntl(3)
      if (mp < 0) print_level = 0

      coarsest_size = max(2,icntl(6))
      if (print_level == 3) &
         call level_print(mp,'size of grid on level ',grid%level,&
         ' is ',real(grid%size,wpr))

      grid_rdc_fac_min = cntl(2)
! max grid reduction factor must be at least half and at most one
      grid_rdc_fac_max = max(half,cntl(3))
      grid_rdc_fac_max = min(one,grid_rdc_fac_max)

! operator complexity and grid complexity
!!      oc = oc + grid%graph%ptr(grid%size+1)-1
!!      gc = gc + grid%size

! Test to see if this is either the last level or
! if the matrix size too small
      if (grid%level >= mglevel_cur .or. grid%size <= coarsest_size) then
          nullify(grid%coarse)
          if (print_level == 3) &
              call level_print(mp,'end of level ',grid%level)

! coarsest level in multilevel so order with Sloan algorithm
          call mc60_wrapper(grid%graph,grid%where,grid%row_wgt,info)
          if (info < 0) return

          if (print_level == 3) &
             call level_print(mp,'after presmoothing ',grid%level)
          return
      end if

! Coarsest level not yet reached so carry on coarsening
      call coarsen(mp,print_level,grid,info)
      if (info < 0) return
      cgrid => grid%coarse
      cnvtx = cgrid%size
! allocate coarse grid quantities
      allocate(cgrid%where(cnvtx),cgrid%row_wgt(cnvtx),stat=st)
      if (st /= 0) then
         info = -1
         return
      end if

! see if the grid reduction is achieved, if not, set the allowed
! maximum level to current level and partition this level
! deallocate the coarse grid quantities that haves been allocated so far
       if (cgrid%size/real(grid%size) >  grid_rdc_fac_max .or. &
           cgrid%size/real(grid%size) <  grid_rdc_fac_min) then

         if (print_level == 3) then
            write (mp,'(a,i10,a,f12.4)') 'at level ',grid%level,&
               ' further coarsening gives reduction factor', &
                 cgrid%size/real(grid%size)
            write (mp,'(a,i10)') 'current size = ',grid%size
         end if

         mglevel_cur = grid%level
         call multilevel(grid,icntl,cntl,mglevel_cur,info)
         if (info < 0) return
         call multigrid_deallocate_last(cgrid,info)
         return
      end if

! restriction ================

! form the coarse grid graph and matrix
! cmatrix = P^T*matrix = R*matrix
      p => cgrid%p
      graph => grid%graph
      cgraph => cgrid%graph

! get the coarse matrix
      call galerkin_graph(graph,p,cgraph,info)
      if (info < 0) return

! row weight cw = R*w
      crow_wgt => cgrid%row_wgt
! numerical experiments found that row_wgt of one
! is best when calling mc60, i.e., do not utilise the row weight.
      crow_wgt = 1

! next level
      call multilevel(cgrid,icntl,cntl,mglevel_cur,info)

! prolongation ================

! injection of the priority from coarse grid to the
! fine grid, since cwhere(i) is the index of the
! i-th vertex in the new ordering, the priority
! of this vertex should be where(i)
! grid%where = P*priority_on_coarse_grid
! here P is a special matrix with only one non-zero entry per row
      fwhere => grid%where
      cwhere => cgrid%where
! maximal ind. vertex set based coarsening
! needs to be rounded to integer since P*cwhere may not be integer
      allocate(fpriority(size(fwhere)),cpriority(cgrid%size),stat=st)
      if (st /= 0) then
         info = -1
         return
      end if
      cpriority = cwhere
      call mc65_matrix_multiply_vector(p,cpriority,fpriority,info)
      fwhere = int(fpriority)
      deallocate(fpriority,cpriority,stat=st)
      if (st /= 0) then
         info = -2
         return
      end if

! post smoothing (refinement) using MC60
     lpermsv = 1
     call mc60_wrapper(grid%graph,grid%where,grid%row_wgt,info,lpermsv)
     if (info < 0) return

     if (print_level == 3) &
           call level_print(mp,' after post smoothing ',grid%level)
! deallocate the previous level
     call multigrid_deallocate(cgrid,info)

      end subroutine multilevel
!***************************************************************
      subroutine coarsen(mp,print_level,grid,ierr)
! coarsen the grid and set up the
! coarse grid equation, the prolongator and restrictor

      integer,  intent(in) :: mp,print_level
      integer,  intent(inout) :: ierr

      type (mc73_multigrid), intent (inout), target :: grid
      type (mc73_multigrid), pointer :: cgrid

      allocate(cgrid)

      cgrid%fine => grid
      grid%coarse => cgrid

      if (print_level == 3) &
         call level_print(mp,'before coarsening',grid%level)

! find the prolongator
      call prolng_max_indset(grid,ierr)

      if (print_level == 3) &
         call level_print(mp,'after coarsening ',grid%level)

      cgrid%level = grid%level + 1

      end subroutine coarsen
!***************************************************************

      subroutine prolng_max_indset(grid,ierr)

! input fine grid
      type (mc73_multigrid), intent (inout) :: grid
      integer :: ierr
! coarse grid based on the fine grid
      type (mc73_multigrid), pointer :: cgrid

! the fine grid matrix
      type (zd11_type), pointer :: matrix
! the coarse grid prolongator
      type (zd11_type), pointer :: p

! the number of coarse grid vertices
      integer :: cnvtx

! n: size of the matrix
! i,j: loop index
      integer :: st

! color: the color of each vertex
!  = C_NODE if this is the i-th C ndoe
!  = F_NODE if colored F

      integer, dimension (:), pointer :: color

      cgrid => grid%coarse
      allocate(cgrid%p)
      allocate(cgrid%graph)

      matrix => grid%graph
      p => cgrid%p

      allocate(color(grid%size),stat=st)
      if (st /= 0) then
         ierr = MC73_ERR_MEMORY_ALLOC
         return
      end if
      call color_init(matrix,color,ierr)
      if (ierr /= 0) then
         if (ierr == -1) ierr = MC73_ERR_MEMORY_ALLOC
         if (ierr == -2) ierr = MC73_ERR_MEMORY_ALLOC
         return
      end if

      call color_getweight(matrix,color,p,cnvtx,ierr)
      if (ierr /= 0) then
         if (ierr == -1) ierr = MC73_ERR_MEMORY_ALLOC
         if (ierr == -2) ierr = MC73_ERR_MEMORY_ALLOC
         return
      end if

! size of coarse grid
      cgrid%size = cnvtx

      end subroutine prolng_max_indset
!***************************************************************
      subroutine color_init(matrix,color,ierr)

!      integer  :: UNCOLORED, C_NODE, F_NODE, FF_NODE
!      parameter (UNCOLORED=10, C_NODE=0, F_NODE=-2, FF_NODE=-1)
      integer  :: UNCOLORED, C_NODE, F_NODE
      parameter (UNCOLORED=10, C_NODE=0, F_NODE=-2)

! matrix: the matrix for which a coarse matrix is to be
! generated
    type (zd11_type), intent (in) :: matrix
    integer :: ierr

! color: the color of each vertex
! = UNCOLORED if not colored
! = C_NODE if colored C
! = FF_NODE if a forced F node (node with no strong connection at all)
! = F_NODE if colored F
      integer, dimension (:), intent (out) :: color

! the queue that holds the uncolored nodes
      type (queue_type) :: queue

! gain: the gain of each vertex, which is defined
! roughly as the number of vertices a vertex controls (is depended by),
! and more accurately is defined as
! (number_of_uncoloured_neighbors +  2* number_of_F_neighbors)
! during the initial colouring stage.
      integer, allocatable, dimension(:) :: gain

! maxgain: maximum gain that could be possibly achieved
      integer :: maxgain

! thegain: gain of the node of concern
! slave: the node that is controlled by a node of concern
! best: the uncolored node with the highest gain
      integer :: thegain,slave,best

! n: size of the matrix
      integer :: n

! stat parameter
      integer :: st

! i,j,k: loop index
      integer :: i,j,k

      n = matrix%n
      allocate (gain(n),stat=st)
      if (st /= 0) then
        ierr = -1
        return
      end if
! initialise the gain
      do i = 1,n
         gain(i) = matrix%ptr(i+1) - matrix%ptr(i)
      end do
      maxgain = maxval(gain)*2

! set up the queue for the gain buckets
      call queue_init(queue,n,maxgain,0,ierr)
      if (ierr /= 0) return

! initially all vertices to be uncolored
      color = UNCOLORED

      do i = 1,n
         if (gain(i) == 0) then
! always color vertices that no one depends on as F node
            color(i) = F_NODE
         else
            call queue_add(queue,i,gain(i))
         end if
      end do

      COLORING: do
! find the remaining uncolored nodes with the highest gain
         call queue_getmax(queue,best,thegain)
! if no more uncoloured vertex we are done
         if (best == -1) exit COLORING
! color this node as C node and delete from the queue
         color(best) = C_NODE
         call queue_delete(queue,best,thegain)
! update the gain of neighbor
         NEIGHB: do j = matrix%ptr(best), matrix%ptr(best+1)-1
            slave = matrix%col(j)
! if slave neighbor of this C node is uncolored, color it as F node,
            if (color(slave) /= UNCOLORED) cycle NEIGHB
            thegain = gain(slave)
            call queue_delete(queue,slave,thegain)
            color(slave) = F_NODE
! up the gain of the slave neighbors of neighbor F nodes
            NEIGHBNEIGHB: do k = matrix%ptr(slave),matrix%ptr(slave+1)-1
               slave = matrix%col(k)
               if (color(slave) /= UNCOLORED) cycle NEIGHBNEIGHB
               thegain = gain(slave)
               call queue_update(queue,slave,thegain,thegain+1)
               gain(slave) = thegain + 1
            end do NEIGHBNEIGHB
         end do NEIGHB
      end do COLORING
      call queue_deallocate(queue,ierr)

      end subroutine color_init
!********************************************
      subroutine color_getweight(matrix,color,p,cnvtx,ierr)
      integer :: C_NODE, F_NODE
      parameter (C_NODE=0, F_NODE=-2)
      real (kind = wpr) :: one
      parameter (one=1.0_wpr)
! color: the color of each vertex
! = C_NODE for C node
! = F_NODE if colored F
! on exit:
! = i if this is the i-th C ndoe

! matrix: the matrix for which a coarse matrix is to be
! generated
      type (zd11_type), intent (in) :: matrix
! p: the prolongation matrix
      type (zd11_type), intent (out) :: p

      integer ierr

      integer, dimension (:), intent (inout) :: color
! cnvtx: the number of coarse grid vertices.
      integer, intent (out) :: cnvtx

! nz: the total number of nonzeros in the prolongation matrix
      integer :: nz

! n: size of the matrix
      integer :: n

! i,j: loop index
      integer :: i,j,k1,k2

! mask: initialised to zero, this
! tells which F node is this C node last strongly connected to
!      integer, dimension (matrix%n) :: mask

! master: nodes which control a node of concern
      integer :: master

! neighbpt: temp integer to hold the current entry in the row i
      integer :: neighbpt

      n = matrix%n
! work out the number of nonzeros in the matrix
      nz = 0
      cnvtx = 0
      ALL_NODES0: do i = 1,n
! if this is a C point, then interpolation weight is just one;
! if it is a normal F point, interpolation is needed.
         if (color(i) >= C_NODE) then
            nz = nz + 1
            cnvtx = cnvtx + 1
            color(i) = cnvtx
            cycle
         end if

         C_i0: do j = matrix%ptr(i),matrix%ptr(i+1) - 1
            master = matrix%col(j)
! ignore F neighbors first
            if (color(master) >= C_NODE) then
! mask C neighbor
               nz = nz + 1
            end if
         end do C_i0
      end do ALL_NODES0

      call mc65_matrix_construct(p,n,nz,ierr,cnvtx)
      if (ierr < 0) return

      p%ptr(1) = 1
!      mask = 0
      ALL_NODES: do i = 1,n
! if this is a C point, then interpolation weight is just one;
! if it is a normal F point, interpolation is needed.
         if (color(i) >= C_NODE) then
            neighbpt = p%ptr(i)
! add coarse node color(i) at row i
            p%col(neighbpt) = color(i)
            p%val(neighbpt) = one
            p%ptr(i+1) = neighbpt + 1
            cycle ALL_NODES
         end if

! mask all the C points that are connected to i
         neighbpt = p%ptr(i)
         C_i: do j = matrix%ptr(i),matrix%ptr(i+1) - 1
            master = matrix%col(j)
! ignore F neighbors first
            if (color(master) <= F_NODE) cycle C_i
! mask C neighbor
!            mask(master) = neighbpt
! this F node i is connected with C node master with a weight of 1/deg(i)
            p%col(neighbpt) = color(master)
            neighbpt = neighbpt + 1
         end do C_i
         p%ptr(i+1) = neighbpt
         k1 = p%ptr(i)
         k2 = p%ptr(i+1)
         p%val(k1:k2-1) = one/(k2 - k1)

      end do ALL_NODES

      end subroutine color_getweight

!*************************************************
      subroutine queue_init(queue,maxnodes,maxgain,mingain,ierr)

! initialise the queue
      type (queue_type), intent (out) :: queue
! maxnodes is max. no. entries allowed in queue
      integer, intent (in) :: maxnodes,maxgain,mingain
      integer :: ierr,i,st

! allocate buckets
      allocate(queue%buckets(mingain:maxgain),stat=st)
      if (st /= 0) then
         ierr = -1
         return
      end if
      queue%low_bound = mingain
      queue%up_bound = maxgain

! no maxgain yet (all buckets empty)
      queue%maxgain = queue%low_bound - 1

! no nodes
      queue%num_nodes = 0
! nullify the buckets
      do i = queue%low_bound, queue%up_bound
         nullify(queue%buckets(i)%current_node)
      end do

! allocate node pointers which say where in the buckets it locates
      allocate(queue%mynode(maxnodes),stat=st)
      if (st /= 0) then
         ierr = -1
         return
      end if
      queue%maxnodes = maxnodes

      end subroutine queue_init
!*************************************************

      subroutine queue_add(queue,id,gain)
! add into the queue a node with gain 'gain'

      type (queue_type), intent (inout) :: queue
      integer, intent (in) :: id,gain

      type (list_node_type), pointer :: new_node

!    if (gain > queue%up_bound.or. gain < queue%low_bound) then
!       write (*,*) 'GAIN = ',GAIN, ' BOUND = ',queue%low_bound,queue%up_bound
!       stop 'in queue_add, gain out of range'
!    end if

!    if (id < 0 .or. id > queue%maxnodes) then
!       stop 'in queue_add, id out of range'
!    end if

! node add one
      queue%num_nodes = queue%num_nodes + 1

! update the current maxgain
      if (gain > queue%maxgain) queue%maxgain = gain

      new_node => queue%mynode(id)

! fill its content
      new_node%id = id

! add this node at the front of the queue
      new_node%next => queue%buckets(gain)%current_node
      nullify(new_node%prev)

! creat the link provided that new_node is not the only one in the list
      if (associated(new_node%next)) new_node%next%prev => new_node

! the bucket now starts from this node
      queue%buckets(gain)%current_node => new_node

      end subroutine queue_add
!*************************************************
      subroutine queue_deallocate(queue,ierr)
! deallocate a queue
      type (queue_type), intent (inout) :: queue
      integer :: ierr

! deallocate buckets array
      deallocate(queue%buckets,stat=ierr)
      if (ierr /= 0) then
         ierr = -2
         return
      end if

! deallocate node pointers which say where in the buckets
! it locates
      deallocate(queue%mynode,stat=ierr)
      if (ierr /= 0) ierr = -2

      end subroutine queue_deallocate

!*************************************************
      subroutine queue_update(queue,id,oldgain,newgain)
! add into the queue a node with gain 'newgain'
! delete into the queue a node with gain 'oldgain'
      type (queue_type), intent (inout) :: queue
      integer, intent (in) :: id,oldgain,newgain

      type (list_node_type), pointer :: new_node

      new_node => queue%mynode(id)

! ============== delete the node from oldgain bucket ==============
! if the node to be deleted is not the last in the bucket
      if (associated(new_node%next)) then
         new_node%next%prev => new_node%prev
      end if

! if the node to be deleted is not the first in the bucket
      if (associated(new_node%prev)) then
         new_node%prev%next => new_node%next
! if it is the first then bucket will start from the next
      else
         queue%buckets(oldgain)%current_node => new_node%next
      end if

! =============== add newgain to bucket ===============
! add this node at the front of the queue
      new_node%next => queue%buckets(newgain)%current_node
      nullify(new_node%prev)

! create the link provided that new_node is not the only one in the list
      if (associated(new_node%next)) then
         new_node%next%prev => new_node
      end if

! the bucket now start from this node
      queue%buckets(newgain)%current_node => new_node

! =============== update the current maxgain =============
      if (newgain >= queue%maxgain) then
         queue%maxgain = newgain
      else if (oldgain == queue%maxgain) then
! find the highest none empty buckets
        do while (.not.associated(queue%buckets(queue%maxgain)%current_node))
           queue%maxgain = queue%maxgain - 1
        end do
      end if

      end subroutine queue_update
!*************************************************

      subroutine queue_getmax(queue,id,gain)
! give the node with the largest gain in the queue
! last in, first out
      type(queue_type), intent(in) :: queue
      integer, intent (out) :: id,gain

      if (queue%num_nodes == 0) then
         gain = 0
         id = -1
         return
      end if
      gain = queue%maxgain
      id = queue%buckets(gain)%current_node%id

      end subroutine queue_getmax
!*************************************************

      subroutine queue_delete(queue,id,gain)
! delete the node id from bucket gain

      type (queue_type), intent (inout) :: queue
      integer, intent (in) :: id,gain

      type (list_node_type), pointer :: new_node

!      if (gain > queue%up_bound.or. gain < queue%low_bound) then
!         stop 'in queue_delete, gain out of range'
!      end if

!      if (id < 0 .or. id > queue%maxnodes) then
!         stop 'in queue_delete, id out of range'
!      end if

! node less one
      queue%num_nodes = queue%num_nodes - 1

      new_node => queue%mynode(id)

! if the node to be deleted is not the last in the bucket
    if (associated(new_node%next)) then
       new_node%next%prev => new_node%prev
    end if

! if the node to be deleted is not the first in the bucket
      if (associated(new_node%prev)) then
         new_node%prev%next => new_node%next
! if it is the first then bucket will start from the next
      else
         queue%buckets(gain)%current_node => new_node%next
      end if

! update the max gain
      if (gain == queue%maxgain) then
! find the highest none empty buckets
         if (queue%num_nodes == 0) then
            queue%maxgain = queue%low_bound - 1
            return
         end if
         do while (.not.associated(queue%buckets(queue%maxgain)%current_node))
            queue%maxgain = queue%maxgain - 1
         end do
      end if

      end subroutine queue_delete
! *****************************************************************
      subroutine multigrid_deallocate(grid,info)
! deallocate a grid (at given level between last and first)
      type (mc73_multigrid), pointer :: grid
      integer :: st,info,info65

      call mc65_matrix_destruct(grid%graph,info65)
      if (info65 /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      call mc65_matrix_destruct(grid%p,info65)
      if (info65 /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      deallocate(grid%graph, grid%p, grid%where, grid%row_wgt, stat=st)
      if (st /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      deallocate(grid, stat=st)
      if (st /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      end subroutine multigrid_deallocate

! *****************************************************************
      subroutine multigrid_deallocate_last(grid,info)
! deallocate a grid (at the last level). In this case the matrix grid%graph
! has not been formed yet
      type (mc73_multigrid), pointer :: grid
      integer, intent(inout)  :: info
      integer :: st,info65

      call mc65_matrix_destruct(grid%p,info65)
      if (info65 /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      deallocate(grid%graph, grid%p, grid%where, grid%row_wgt, stat=st)
      if (st /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      deallocate(grid,stat=st)
      if (st /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      end subroutine multigrid_deallocate_last
! *****************************************************************
      subroutine multigrid_deallocate_first(nvtx,n,grid,info)
! deallocate a grid (at the first level). In this case the matrix grid%p
! does not exist
      type (mc73_multigrid) :: grid
      integer, intent(inout)  :: info
      integer :: n,nvtx
      integer :: ierr

      info = 0
      if (nvtx /= n) then
         call mc65_matrix_destruct(grid%graph,ierr)
         if (ierr /= 0) then
            info = MC73_ERR_MEMORY_DEALLOC
            return
         end if
      end if
!    deallocate(grid%graph, grid%where, grid%row_wgt, stat=ierr)
! in subroutine front, grid%graph is not allocated but is pointed to the
! finest level graph. so no need to deallocate
      nullify(grid%graph)
      deallocate(grid%where,grid%row_wgt,stat=ierr)
      if (ierr /= 0) info = MC73_ERR_MEMORY_DEALLOC

      end subroutine multigrid_deallocate_first
! ***************************************************
      subroutine mc60_wrapper(matrix,permsv,row_wgt,ierr,lpermsv)

! wrapper for mc60
      type (zd11_type), intent (in) :: matrix
! the vertex weight, used as the super-variable weight
      integer, dimension (:) :: row_wgt
! Error flag
! -1 for allocation error
! -2 for deallocation error
      integer :: ierr

! ordering of the rows,
! on entry: the variable with lower value of permsv
! is likely to be chosen first
! on exit: new index for variable I is given by permsv(i)
      integer, intent (inout), dimension (matrix%n) :: permsv

! if lpermsv is present, then permsv must contain input global
! priority function
      integer, intent(in), optional :: lpermsv

! MC60 related variables
! n is the order of the finest matrix
      integer :: n

! the number of super variables (number of coarse graph vertices)
      integer :: nsup

! integer variable set to the length of irn, not altered
      integer :: lirn

! integer array of length nsup. vars(is) is the number of
! variables in super variable (coarse vertex) is. not altered
      integer, dimension (:), allocatable :: vars

! controls:
! jcntl(1) = 0  (sloan's algorithm for reducing profile and wavefront)
!            1  (reverse cuthill-mckee for reducing bandwidth)
! jcntl(2) = 0 (automatic choice of pseudoperipheral pairs
!            1 (pseudoperipheral pairs specified in pair)
!            2 (global priority vector given in permsv
!                (sloan's alg. only))
      integer jcntl(2)


! real of length two. for sloan (jcntl(1) = 0) it must be set to
! w1 and w2 of the priority function (w1*deg(s) + w2*v*glob(s))
! that is minimized when chosing the next supervariable.
      real (kind = wpr) :: weight(2)

! the root pairs. of jcntl(2) = 0, it need not set on entry and on return
! pair(1,ic), pair(2,ic) holds the pseudoperipheral pair for nontrivial
! component ic, ic = 1, 2, ... info(1). the first component is largest.
! if jcntl(2) = 1, it must be set on entry
! to the pseudoperipheral pairs of the components and is not
! altered.  the first
! component need not be the largest. if jcntl(2) = 2
! it is not used.
      integer, dimension (:,:),allocatable :: pair

! integer of length 4, not set on entry
! info(1) number of nontrival components (two or more nodes)
! info(2) numbert of variables in the largest component of the graph
! info(3) number of level sets in the level-set structure of the largest
!         component
! info(4) width of the level-set structure of the largest component
      integer :: info(4)

! integer of length 3*nsup+1 as work space
      integer, dimension (:), allocatable :: iw
! real work array
      real (kind = wpr), dimension (:), allocatable :: w

! real of length 4 which holds profile, maximum wavefront,
! semi-bandwidth and root-mean-square
      real (kind = wpr), dimension (4) :: rinfo
      real (kind = wpr), dimension (4) :: rinfo_temp
      integer, dimension (:), allocatable :: temp
! stat parameter
      integer st

      external mc60cd,mc60fd

      ierr = 0
      nsup = matrix%n
      n = sum(row_wgt)

      lirn = matrix%ptr(nsup+1)-1

      allocate(pair(2,nsup/2),vars(nsup),iw(3*nsup+1),w(nsup),temp(nsup),&
               stat=st)
      if (st /= 0) then
         ierr = -1
         return
      end if
      vars = row_wgt

! use sloan algorithm and global priority vector (if present)
      jcntl(1) = 0
      if (present(lpermsv)) then
         jcntl(2) = 2
         weight(1) = 1.0_wpr
         weight(2) = 2.0_wpr
         temp = permsv
      else
! when no global priority function is present, do Sloan with (2,1)
         weight(1) = 2.0_wpr
         weight(2) = 1.0_wpr
         jcntl(2) = 0
      end if

      call mc60cd(n,nsup,lirn,matrix%col,matrix%ptr,vars,jcntl,&
                  permsv,weight,pair,info,iw,w)
      call mc60fd(n,nsup,lirn,matrix%col,matrix%ptr,vars,permsv,iw,rinfo)

!      write (*,'('prof = ',g16.8,' maxw = ',&
!           &g16.8,' sband = ',g16.8,' rms = ',g16.8)') rinfo(1:4)

! second set of weights for the priority function
      weight(1) = 16.0_wpr
      weight(2) = 1.0_wpr

      call mc60cd(n,nsup,lirn,matrix%col,matrix%ptr,vars,jcntl,&
                  temp,weight,pair,info,iw,w)
      call mc60fd(n,nsup,lirn,matrix%col,matrix%ptr,vars,temp,iw,rinfo_temp)

      if (rinfo_temp(1) < rinfo(1)) then
! selecting ordering on basis of profile
         permsv = temp
         rinfo = rinfo_temp
      end if
!      write (*,'('prof = ',g16.8,' maxw = ',&
!           &g16.8,' sband = ',g16.8,' rms = ',g16.8)') rinfo(1:4)

      deallocate(pair,vars,iw,w,temp,stat=st)
      if (st /= 0) then
         ierr = -2
         return
      end if

      end subroutine mc60_wrapper

! *****************************************************************
      subroutine fiedler_graph(lspec,n,matrix,perm,fvector,info,icntl,  &
                               cntl,llap)

! Subroutine to compute Fiedler vector or spectral ordering
! Fiedler vector is returned in fvector
! Spectral ordering is returned in perm

! If lspec = .false., on return fvector holds the computed
! Fiedler vector.
! If lspec = .true., on return perm holds the
! spectral ordering and fvector holds the computed
! sorted Fiedler vector.

! if llap = .true., then we work with weighted Laplacian (this
! is the case if user supplied values in array a)

! The position of variable I in the spectral ordering is perm(I),
! I = 1, 2, ..., N.  If perm(I) = J then
! fvector(J) is the value of variable I in the Fiedler vector.
! If the graph has more than one component,
! a change to a new component is indicated by a switch from a
! positive value to a negative value.  If a component consists of a
! single vertex, the corresponding variable  is ordered first
! and is assigned a value 0 in fvector.



      real (kind = wpr), parameter :: zero = 0.0_wpr
! control parameters
      integer,  intent (in) :: icntl(10)
      real (kind = wpr), intent (in) :: cntl(10)
! lspec = .true. if spectral permutuation wanted
! lspec = .false. if fiedler is wanted
      logical,  intent (in) :: lspec
      logical,  intent (in) :: llap

! Order of matrix
      integer,  intent (in) :: n
! Spectral ordering if lspec = .true.; otherwise
! used to identify the component to which variable belongs
      integer,  intent (out) :: perm(n)
! information array
      integer,  intent (inout) :: info(10)
! Fiedler vector
      real (kind = wpr), intent (out) :: fvector(n)

! the matrix whose Fiedler eigenvector is to be calculated
      type (zd11_type), target :: matrix

! the tolerance for Lanczos, rayleigh iteration and symmlq
! typically tol = tol1 and rtol = 10*tol or 100*tol
      real (kind = wpr)  :: tol
      real (kind = wpr)  :: tol1
      real (kind = wpr)  :: rtol

! the multilevel of graphs (matrices)
      type (mc73_multigrid_eig) :: grid

! root: the root of the current component
      integer :: root
! comp_nvtx: number of vertices in current component
! comp_nz: number of edges*2 in current component
! ncomp: number of components
      integer :: comp_nvtx,comp_nz,ncomp

! mask: mask that is initialised to zero. After a
! vertex is visited, its mask is set to the new index
! of this vertex in the component
! comp_vtx_list: a list of the vertices in the current component
! mask and comp_vtx_list satisfies mask(comp_vtx(i)) = i,
! i = 1,...,comp_nvtx
      integer,  dimension (:), pointer :: mask, comp_vtx_list
      integer,  dimension (:), allocatable :: comp_perm,iw

! submatrix: the subgraph of this component
      type (zd11_type), pointer :: submatrix

! number of levels in the multilevel grid (default 10)
      integer :: mglevel

! print level
      integer :: err,wrn,mp,print_level
      logical :: lerr,lwrn

! error flag
      integer :: ierr
! stat parameter
      integer :: st

! i: loop index
      integer :: i,ii,j,k
! number of nodes that have been ordered
      integer :: lstnum
      integer :: mglevel_cur
      integer :: mlancz
! Degree of node
      integer :: degree

! for verification only
!      real (kind = wpr) :: lambda
!      real (kind = wpr), dimension (n) :: y
!      real (kind = wpr) :: dnrm2

      external kb07ad

      if (icntl(4) < 0) print_level = 0
! The default is icntl(4) = 0
      if (icntl(4) == 0) print_level = 1
      if (icntl(4) == 1) print_level = 2
      if (icntl(4) > 1) print_level = 3
      mp = icntl(3)
      if (mp < 0) print_level = 0
! Set error controls
      lerr = icntl(1).ge.0 .and. print_level.gt.0
      err  = icntl(1)
! Set warning controls
      lwrn = icntl(2).ge.0 .and. print_level.gt.0
      wrn = icntl(2)

      mglevel = icntl(5)
      mlancz = icntl(8)

! tolerance for Lanczos, rayleigh iter. and symmlq
      tol = cntl(1)
      tol1 = cntl(4)
      rtol = cntl(5)

! work on each connected component
      allocate(mask(n),stat=st)
      if (st /= 0) info(1) = MC73_ERR_MEMORY_ALLOC
      allocate(iw(n),stat=st)
      if (st /= 0) info(1) = MC73_ERR_MEMORY_ALLOC
      allocate(comp_vtx_list(n),stat=st)
      if (st /= 0) info(1) = MC73_ERR_MEMORY_ALLOC
      if (info(1) < 0) then
         if (lerr) call print_message(info(1),err,' MC73')
         return
       end if

! ncomp is number of non-trivial components
      ncomp = 0
      mask = 0
! Order all nodes of degree zero
      lstnum = 0
      if (lspec) then
         do i = 1,n
            k = matrix%ptr(i)
            degree = matrix%ptr(i+1) - k
            if (degree == 0) then
! single node in component ... number it first
              lstnum = lstnum + 1
              perm(i) = lstnum
              mask(i) = 1
              fvector(lstnum) = zero
            end if
         end do
      else
         do i = 1,n
            k = matrix%ptr(i)
            degree = matrix%ptr(i+1) - k
            if (degree == 0) then
! single node in component
              lstnum = lstnum + 1
              perm(i) = 0
              mask(i) = 1
              fvector(i) = zero
            end if
         end do
      end if

! If matrix is diagonal, then we have now finished since all nodes
! have degree zero.
      if (lstnum == n) then
         deallocate(mask,comp_vtx_list,stat=st)
         if (st /= 0) then
           info(1) = MC73_ERR_MEMORY_DEALLOC
           if (lerr) call print_message(info(1),err,' MC73')
         end if
         if (print_level == 3) write (mp,'(a//)') 'finished'
         return
      end if

      do i = 1,n
        if (lstnum == n) exit
! Find first node not yet in a component
         if (mask(i) == 0) then
           root = i
! Find component of graph rooted at node i
           call mx_component(matrix,root,mask,comp_nvtx,comp_nz,&
                             comp_vtx_list,ierr)
! ierr holds stat parameter
           if (ierr /= 0 ) then
              info(1) = ierr
              if (lerr) call print_message(info(1),err,' MC73')
              return
           end if

           ncomp = ncomp + 1
! store number of non-trivial components
           info(2) = ncomp
           if (print_level == 3) then
             write (mp,'(a,i4,a,i10,a,i10)') 'Component ',ncomp,&
                  ' number of nodes = ',comp_nvtx
           end if

           if (comp_nvtx == n) then
! this is a connected graph
             submatrix => matrix
           else
! graph has more than one component
             allocate(submatrix)
! generate submatrix
             call mx_crop(matrix,comp_nvtx,comp_vtx_list,mask,submatrix,&
                          ierr,comp_nz)
             if (ierr < 0) then
                info(1) = ierr
                if (lerr) call print_message(info(1),err,' MC73')
                return
             end if
           end if

! construct the components of grid (which is of derived datatype multigrid)
           grid%size = submatrix%n
           grid%level = 1
           grid%graph => submatrix
           nullify (grid%p)

           allocate(grid%eigvector(grid%size),stat=st)
           if (st /= 0) then
              info(1) = MC73_ERR_MEMORY_ALLOC
              if (lerr) call print_message(info(1),err,' MC73')
              return
           end if

           allocate(grid%row_wgt(grid%size),stat=st)
           if (st /= 0) then
              info(1) = MC73_ERR_MEMORY_ALLOC
              if (lerr) call print_message(info(1),err,' MC73')
              return
           end if

          grid%row_wgt = 1
! maximum number of grid levels allowed for this bisection
           mglevel_cur = mglevel
           ierr = 0
           call multilevel_eig(grid,tol,tol1,rtol,icntl,cntl,  &
                               mglevel_cur,mlancz,ierr,llap)
           if (ierr < 0) then
              if (ierr == -1) info(1) = MC73_ERR_MEMORY_ALLOC
              if (ierr == -2) info(1) = MC73_ERR_MEMORY_DEALLOC
              return
           else if (ierr > 0) then
              if (ierr == 4 .and. info(1) < 4) &
                  info(1) = info(1) + MC73_WARN_MAXIT
              if (lwrn) call print_message(MC73_WARN_MAXIT,wrn,' MC73')
           end if


           if (comp_nvtx == n) then
! connected graph so we are done (this is the 'usual' case)
             fvector = grid%eigvector
!!               write (10,*) fvector

! verify (ie check || L*fvector - lambda*fvector||)
! this is only any use in case of single component
!    lambda = grid%lambda
!    call lxv(fvector,y,grid%graph,llap)
!    write (*,*) 'lambda  = ',lambda
!    write (*,*) 'residual  = ',dnrm2(n,(fvector*lambda+y),1)/dnrm2(n,y,1)

              if (lspec) then
                 do ii = 1,n
                    perm(ii) = ii
                 end do
                 call kb07ad(fvector,n,perm)

! need inverse
! so that new index number for node j is in perm(j)
! (node that is likely to be picked first when
! refined using Sloan approach is i such that perm(i)=1)
                 iw(1:n) = perm(1:n)
                 do ii = 1,n
                    j = iw(ii)
                    perm(j) = ii
                 end do
              else
                 perm(1:n) = 1
              end if

           else

! more than one component
              if (lspec) then
! allocate a permutation vector for this component
                 deallocate(comp_perm,stat=st)
                 allocate(comp_perm(comp_nvtx),stat=st)
                 if (st /= 0) then
                    info(1) = MC73_ERR_MEMORY_ALLOC
                    if (lerr) &
                    call print_message(info(1),err,' MC73')
                    return
                  end if

! Get spectral permutation on component
                 do ii = 1,comp_nvtx
                    comp_perm(ii) = ii
                 end do
                 call kb07ad(grid%eigvector,comp_nvtx,comp_perm)

! Take inverse
                iw(1:comp_nvtx) = comp_perm(1:comp_nvtx)
                do ii = 1,comp_nvtx
                   j = iw(ii)
                   comp_perm(j) = ii
                end do

                 do ii = 1,comp_nvtx
                   j = comp_vtx_list(ii)
                   fvector(ii+lstnum) = grid%eigvector(ii)
                   perm(j) = comp_perm(ii) + lstnum
                 end do

              else
! User requires the Fiedler vector
                 do ii = 1,comp_nvtx
                   j = comp_vtx_list(ii)
                   fvector(j) = grid%eigvector(ii)
                   perm(j) = ncomp
                 end do
              end if

           end if

           lstnum = lstnum + comp_nvtx

! deallocate the finest level
! if the graph has a single component, we do not want to destruct matrix
           call multigrid_eig_deallocate_first(comp_nvtx,n,grid,ierr)
           if (ierr /= 0) then
              info(1) = ierr
              if (lerr) call print_message(info(1),err,' MC73')
              return
           end if

         end if

      end do

      deallocate(mask,comp_vtx_list,iw,stat=st)
      if (st /= 0) info(1) = MC73_ERR_MEMORY_DEALLOC
      if (lspec .and. comp_nvtx /= n) then
         deallocate(comp_perm,stat=st)
         if (st /= 0) info(1) = MC73_ERR_MEMORY_DEALLOC
      end if
      if (info(1) < 0) then
         if (lerr) call print_message(info(1),err,' MC73')
         return
      end if
      if (print_level == 3) write (mp,'(a//)') 'finished'

      end subroutine fiedler_graph

! ********************************************************
      subroutine mx_crop(matrix,comp_nvtx,comp_vtx_list,mask,submatrix,&
                         ierr,comp_nz)
! mx_crop: generate a submatrix describing a component of
! the graph (matrix).

! matrix: the graph to be inspected. Must be a symmetric
! matrix with no diagonal.
     type (zd11_type), intent (in) :: matrix

! comp_nvtx: number of vertices in this component
! comp_nz: number of edges*2 in the component (optional)
      integer,  intent (in) :: comp_nvtx
      integer,  intent (in), optional :: comp_nz
      integer,  intent (out) :: ierr
!  ierr = 0 if successful
!   = MC73_ERR_MEMORY_ALLOC if memory allocation failed
!   = MC73_ERR_MEMORY_DEALLOC if memory deallocation failed

! mask: has the dimension as the number of vertices.
! The mask of a vertex is set to the new index
! of this vertex in the component to which it belongs
! comp_vtx_list: a list of the vertices in the current component
      integer,  dimension (:), intent (in) :: mask,comp_vtx_list

! submatrix: the submatrix describing the component
      type (zd11_type), intent (out) :: submatrix

! i: loop index
! nz: number of nonzeros in the submatrix
! n: number of rows and cols in the submatrix
! j: a vertex in its original index
! l1: starting index
! l2: stopping index
! sl1: starting index
! sl2: stopping index
      integer :: i,nz,n,j,l1,l2,sl1,sl2

      ierr = 0

      if (.not.present(comp_nz)) then
         nz = 0
         do i = 1,comp_nvtx
            j = comp_vtx_list(i)
            nz = nz + matrix%ptr(j+1) - matrix%ptr(j)
         end do
      else
         nz = comp_nz
      end if

      n = comp_nvtx
! storage allocation for col. indices and values
      call mc65_matrix_construct(submatrix,n,nz,ierr,n,   &
                                 type=zd11_get(matrix%type))
      if (ierr < 0) then
         if (ierr == -1) ierr = MC73_ERR_MEMORY_ALLOC
         if (ierr == -2) ierr = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      submatrix%ptr(1) = 1
      do i = 1,comp_nvtx
         j = comp_vtx_list(i)
         l1 = matrix%ptr(j)
         l2 = matrix%ptr(j+1)-1
         submatrix%ptr(i+1) = submatrix%ptr(i) + l2 - l1 + 1
         sl1 = submatrix%ptr(i)
         sl2 = submatrix%ptr(i+1)-1
 ! convert original index to new index in the component
         submatrix%col(sl1:sl2) = mask(matrix%col(l1:l2))
         if (zd11_get(matrix%type) /= 'pattern' ) &
            submatrix%val(sl1:sl2) = matrix%val(l1:l2)
      end do

      end subroutine mx_crop
! ************************************************************
      subroutine mx_component(matrix,root,mask,comp_nvtx,    &
                              comp_nz,comp_vtx_list,ierr)

! matrix_component: this subroutine finds the component
! of the graph starting from the root. This component
! is stored in the list of vertices 'comp_vtx_list',
! with 'comp_nvtx' vertices and 'comp_nz'/2 edges
! (thus there are comp_nz nonzeros in the matrix describing
! the component).

! matrix: the graph to be inspected. Must be a symmetric
! matrix with no diagonal.
        type (zd11_type), intent (in) :: matrix

! root: the root of the current component
        integer,  intent (in) :: root

! comp_nvtx: number of vertices in this component
! comp_nz: number of edges*2 in the component
        integer,  intent (out) :: comp_nvtx,comp_nz

! ierr: error tag
        integer,  intent (out) :: ierr

! mask: has the dimension as the number of vertices.
! It is zero if the vertex is not yet visited,
! nonzero if visited. When a vertex is visited the first time,
! its mask is set to the new index of this vertex in the component
! comp_vtx_list: a list of the vertices in the current component.
! mask and comp_vtx_list satisfies mask(comp_vtx(i)) = i,
! i = 1,...,comp_nvtx
        integer,  dimension (:), intent (inout) :: mask, &
             comp_vtx_list

! front_list: this array hold the list of vertices that
!         form the front in the breadth first search
!         starting from the root
        integer,  dimension (:), allocatable :: front_list

! front_sta: the starting position in the front_list of the
!        current front.
! front_sto: the ending position in the front_list of the current
! front.
        integer :: front_sta, front_sto

! n: number of vertices in the graph
        integer :: n
! v: a vertex
! u: neighbor of v
! i: loop index
! j: loop index
        integer :: i,j,u,v

        ierr = 0
        n = matrix%m

! the first front
        allocate(front_list(n),stat=ierr)
        if (ierr /= 0) then
           ierr = MC73_ERR_MEMORY_ALLOC
           return
        end if
        front_sta = 1
        front_sto = 1
        front_list(1) = root

        comp_nvtx = 1
        comp_nz = 0
! one vertex in the component so far
        comp_vtx_list(1) = root
! mask the root
        mask(root) = comp_nvtx

        do while (front_sto-front_sta >= 0)
           do i = front_sta, front_sto
! pick a vertex from the front
              v = front_list(i)
! count link to all neighbors as edges
              comp_nz = comp_nz + matrix%ptr(v+1) - matrix%ptr(v)
              do j = matrix%ptr(v),matrix%ptr(v+1)-1
! pick its neighbor
                 u = matrix%col(j)
                 if (mask(u) /= 0) cycle
! found a unmasked vertex
                 comp_nvtx = comp_nvtx + 1
! mask this vertex
                 mask(u) = comp_nvtx
! add this vertex to the component
                 comp_vtx_list(comp_nvtx) = u
! also add it to the front
                 front_list(comp_nvtx) = u
              end do
           end do
           front_sta = front_sto + 1
           front_sto = comp_nvtx
        end do

        deallocate(front_list,stat = ierr)
        if (ierr /= 0) ierr = MC73_ERR_MEMORY_DEALLOC

      end subroutine mx_component
!***************************************************************
      recursive subroutine multilevel_eig(grid,tol,tol1,rtol,icntl, &
                                     cntl,mglevel_cur,mlancz,info,llap)

      real (kind = wpr), parameter :: half = 0.5_wpr
      real (kind = wpr), parameter :: one = 1.0_wpr

! this level of matrix (grid)
      type (mc73_multigrid_eig), intent (inout), target :: grid

! the tolerance for Lanczos, Rayleigh quotient iteration and symmlq
! typically tol = tol1 and rtol = 10*tol or 100*tol
      real (kind = wpr), intent (in) :: tol,tol1,rtol

      integer,  intent(in) :: icntl(10)
      real (kind = wpr), intent (in) :: cntl(10)
      integer :: mglevel_cur
      integer :: mlancz
      logical :: llap
! Error flag
      integer :: info
! the coarse level grid
      type (mc73_multigrid_eig), pointer :: cgrid

! the  prolongator and restrictor
      type (zd11_type), pointer :: p

! the eigvector in fine and coarse grid
      real (kind = wpr), dimension (:), pointer :: feigvector, ceigvector

! the coarse graph
      type (zd11_type), pointer :: cgraph

! the fine graph
      type (zd11_type), pointer :: graph

! fine and coarse graph vertex weight
      integer,  dimension (:), pointer :: row_wgt, crow_wgt
      real (kind = wpr) :: grid_rdc_fac_min,grid_rdc_fac_max

! smallest problem size to stop coarsening (default 100)
      integer :: coarsest_size
      integer   :: mp, print_level
      integer   :: maxit


      info = 0
      if (icntl(4) < 0) print_level = 0
! The default is icntl(4) = 0
      if (icntl(4) == 0) print_level = 1
      if (icntl(4) == 1) print_level = 2
      if (icntl(4) > 1) print_level = 3
      mp = icntl(3)
      if (mp < 0) print_level = 0

      coarsest_size = max(2,icntl(6))
      maxit = max(1,icntl(7))

      grid_rdc_fac_min = cntl(2)
! max grid reduction factor must be at least half and at most one
      grid_rdc_fac_max = max(half,cntl(3))
      grid_rdc_fac_max = min(one,grid_rdc_fac_max)

      if (print_level == 3) &
         call level_print(mp,'size of grid on level ',grid%level,&
         ' is ',real(grid%size,wpr))

! if this is either the last level or grid size too small
! then compute eigenvector
      if (grid%level >= mglevel_cur .or. grid%size <= coarsest_size) then
         nullify (grid%coarse)
         if (print_level == 3)    &
         call level_print(mp,'end of level ',grid%level)
         call lanczos(grid%graph,grid%graph%n,grid%eigvector,grid%lambda,  &
                      mlancz,tol,6,1,info,llap)
         return
      end if

! Otherwise, continue coarsening
      call coarsen_eig(mp,print_level,grid,info)
      if (info < 0) return

      cgrid => grid%coarse

!      if (print_level == 3) &
!         call level_print(mp,'size of cgrid on level ',grid%level+1,&
!         ' is ',real(cgrid%size,wpr))

! see if the grid reduction is within range. If not, set the allowed
! maximum level to current level and partition this level
! deallocate the coarse grid quantities that has been allocated so far
      if (cgrid%size/real(grid%size) >  grid_rdc_fac_max .or. &
          cgrid%size/real(grid%size) <  grid_rdc_fac_min) then

         if (print_level == 3) then
            write (mp,'(a,i10,a,f12.4)') 'at level ',grid%level,&
               ' further coarsening gives reduction factor', &
                 cgrid%size/real(grid%size)
            write (mp,'(a,i10)') 'current size = ',grid%size
         end if

         mglevel_cur = grid%level
         call multilevel_eig(grid,tol,tol1,rtol,icntl,cntl,  &
                             mglevel_cur,mlancz,info,llap)
         if (info < 0) return
         call multigrid_eig_deallocate_last(cgrid,info)
         return
      end if

! restriction
! form the coarse grid graph and matrix
! cmatrix = P^T*matrix = R*matrix
      p => cgrid%p
      graph => grid%graph
      cgraph => cgrid%graph

! get the coarse matrix
      call galerkin_graph(graph,p,cgraph,info)
      if (info < 0) return

! row weight cw = R*w
      row_wgt => grid%row_wgt
      crow_wgt => cgrid%row_wgt
      call mc65_matrix_multiply_vector(p,row_wgt,crow_wgt,info,trans=.true.)

! next level
       call multilevel_eig(cgrid,tol,tol1,rtol,icntl,cntl,  &
                           mglevel_cur,mlancz,info,llap)
       if (info < 0) return

        grid%lambda = cgrid%lambda
! prolongation
! injection of the eigvector from coarse grid to the fine grid
! feigvector = P*ceigvector
       feigvector => grid%eigvector
       ceigvector => cgrid%eigvector
       call mc65_matrix_multiply_vector(p,ceigvector,feigvector,info)

! post smoothing (refinement) using Rayleigh Quotient Iteration
! note: grid%lambda is reset by rayleigh to the shift
! maxit is max. number of RQI iterations
       call rayleigh(grid%graph,feigvector,grid%lambda,tol1,rtol, &
                     maxit,info,print_level,mp,llap)
       if (info < 0) return

       if (print_level == 3) &
           call level_print(mp,' after post smoothing ',grid%level)

! deallocate the previous level
      call multigrid_eig_deallocate(cgrid,info)

      end subroutine multilevel_eig

!***************************************************************
      subroutine coarsen_eig(mp,print_level,grid,ierr)
! coarsen the grid and set up the
! coarse grid equation, the prolongator and restrictor

      integer,  intent(in) :: mp,print_level
      integer,  intent(inout) :: ierr
      type (mc73_multigrid_eig), intent (inout), target :: grid
      type (mc73_multigrid_eig), pointer :: cgrid

      allocate(cgrid)

      cgrid%fine => grid
      grid%coarse => cgrid

      if (print_level == 3)  &
         call level_print(mp,'before coarsening',grid%level)

! find the prolongator
      call prolng_heavy_edge_eig(grid,ierr)

      if (print_level == 3)  &
         call level_print(mp,'after coarsening ',grid%level)

      cgrid%level = grid%level + 1

      end subroutine coarsen_eig

!********************************************
      subroutine prolng_heavy_edge_eig(grid,ierr)

!    calculate the prolongator:
!    match the vertices of the heaviest edges

      integer,  intent(inout) :: ierr
! input fine grid
      type (mc73_multigrid_eig), intent (inout) :: grid

! coarse grid based on the fine grid
      type (mc73_multigrid_eig), pointer :: cgrid

! the fine grid row connectivity graph
      type (zd11_type), pointer :: graph

! the coarse grid prolongator
      type (zd11_type), pointer :: p

! the number of fine and coarse grid vertices
      integer :: nvtx,cnvtx

! working variables
      integer :: v,u,j,i,k,st

      integer, pointer, dimension (:) :: the_row
      real (kind = wpr), pointer, dimension (:) :: the_row_val

! whether a vertex is matched already
      integer, parameter :: unmatched = -1

! matching status of each vertex
      integer, pointer, dimension (:) :: match

! maximum weight and index of edges connected to the current vertex
      real (kind = wpr) :: maxwgt
      integer :: maxind

! order of which vertex is visited for matching
      integer, allocatable, dimension (:) :: order
      integer, allocatable, dimension (:) :: nrow

      integer :: nz
      integer :: info

! allocate the prolongation matrix pointers
      cgrid => grid%coarse
      graph => grid%graph
      allocate(cgrid%p)

! allocate the graph and matrix pointer and the mincut pointer
! so that everything is defined
      allocate(cgrid%graph)

      p => cgrid%p
      nvtx = graph%n

! prolongator start here ================================

! initialise the matching status
      allocate(match(nvtx),stat = st)
      if (st /= 0) then
         ierr = MC73_ERR_MEMORY_ALLOC
         return
      end if
      match = unmatched

! randomly permute the vertex order
      allocate (order(nvtx),stat=st)
      if (st /= 0) then
         ierr = MC73_ERR_MEMORY_ALLOC
         return
      end if
      do i = 1,nvtx
         order(i) = i
      end do

      allocate (nrow(grid%size),stat=st)
      if (st /= 0) then
         ierr = MC73_ERR_MEMORY_ALLOC
         return
      end if

! loop over each vertex and match along the heaviest edge
      cnvtx = 0
      nz = 0
      nrow = 0
      do i = 1, nvtx
        v = order(i)
! If already matched, next vertex please
        if (match(v) /= unmatched) cycle
! access the col. indices of row v
        call mc65_matrix_getrow(graph,v,the_row)
! access the entry values of row v
        call mc65_matrix_getrowval(graph,v,the_row_val)
        maxwgt = -huge(0.0_wpr)
! in the case no match is found then match itself
        maxind = v
! Loop over entries in row v
        do j = 1, size(the_row)
! u is col index of en entry in row v (so u is neighbor of v)
           u = the_row(j)
! heavy edge matching
! if u is unmatched and value of the entry in col. u is greater
! than maxwgt, select u as the matching.
           if (match(u)==unmatched .and. maxwgt < abs(the_row_val(j))) then
              maxwgt = abs(the_row_val(j))
              maxind = u
           end if
        end do
! the neighbor with heaviest weight
        match(v) = maxind
! mark maxind as having been matched
        match(maxind) = v
! increase number of vertices in coarse graph by 1
        cnvtx = cnvtx + 1
! construct the prolongation matrix: find vertex v and maxind is linked
! with the coarse grid vertex cnvtx
        nz = nz + 1
        nrow(v) = nrow(v) + 1
        if (maxind /= v) then
           nz = nz + 1
           nrow(maxind) = nrow(maxind) + 1
        end if
      end do

! storage allocation for col. indices and values of prolongation
! matrix P (order nvtx * cnvtx)
      call mc65_matrix_construct(p,nvtx,nz,info,n=cnvtx,type = 'general')
      p%val = 0.0_wpr

! Set column pointers for P
      p%ptr(1) = 1
      do i = 1,nvtx
         p%ptr(i+1) = p%ptr(i) + nrow(i)
      end do

! the row counter: nothing in matrix P filled yet
      nrow = 0

! Now fill in entries of matrix P
      match = unmatched
! loop over each vertex and match along the heaviest edge
      cnvtx = 0
      do i = 1,nvtx
        v = order(i)
! if already matched, next vertex please
        if (match(v) /= unmatched) cycle
        call mc65_matrix_getrow(graph,v,the_row)
        call mc65_matrix_getrowval(graph,v,the_row_val)
        maxwgt = -huge(0.0_wpr)
! in the case no match is found then match itself
        maxind = v
        do j = 1, size(the_row)
           u = the_row(j)
! heavy edge matching
           if (match(u)==unmatched .and. maxwgt < abs(the_row_val(j))) then
              maxwgt = abs(the_row_val(j))
              maxind = u
           end if
        end do
! the neighbor with heaviest weight
        match(v) = maxind
        match(maxind) = v
        cnvtx = cnvtx + 1
! construct the prolongation matrix: vertex v and maxind are linked
! with the coarse grid vertex cnvtx
! k points to start of row v and now(v) holds number of
! entries in row v that have been filled (so k + nrow(v) is
! first free entry in row v)
        k = p%ptr(v)
        p%col(k+nrow(v)) = cnvtx
        p%val(k+nrow(v)) = 1.0_wpr
        nrow(v) = nrow(v) + 1

        if (maxind /= v) then
           k = p%ptr(maxind)
           p%col(k+nrow(maxind)) = cnvtx
           p%val(k+nrow(maxind)) = 1.0_wpr
           nrow(maxind) = nrow(maxind) + 1
        end if

      end do

! deallocate the match pointer for vertices
      deallocate(match,stat=st)
      if (st /= 0) then
         ierr = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      deallocate(order,nrow,stat=st)
      if (st /= 0) then
         ierr = MC73_ERR_MEMORY_DEALLOC
         return
      end if

! prolongator end
! size of coarse grid
      cgrid%size = cnvtx

! allocate coarse grid quantities
      allocate(cgrid%eigvector(cnvtx),cgrid%row_wgt(cnvtx),stat=st)
      if (st /= 0) ierr = MC73_ERR_MEMORY_ALLOC

      end subroutine prolng_heavy_edge_eig

!*******************************************************************
      subroutine level_print(mp,title1,level,title2,res)

! print the title in a indented fasion depend on which level.
! will give output like this
!   ===== title on level   1 =====
!  ===== title on level   2 =====
!     ===== title on level   3 =====
!        ===== title on level   4 =====
!           ===== title on level   5 =====
!              ===== title on level   6 =====
!                  .......
!
    character (len = *), intent(in) :: title1
    integer,  intent(in) :: mp,level
    real (kind = wpr), optional, intent (in) :: res
    character (len = *), optional, intent(in) :: title2
    character (len=80) fmt
    integer :: char_len1,char_len2

    char_len1=len_trim(title1)

      if (present(res) .and. present(title2)) then
       char_len2=len_trim(title2)
       write (fmt, &
         "('(',i4,'X,''===== '',a',i4,',i4,a',i4,',g14.3,'' ====='')')") &
            level*3, char_len1,char_len2
       write (mp,fmt) title1,level,title2,res
      else if (present(res)) then

       write (fmt,  &
         "('(',i4,'X,''===== '',a',i4,',i4,'' is'',g14.3,'' ====='')')") &
            level*3, char_len1
       write (mp,fmt) title1,level,res
     else
       write (fmt,"('(',i4,'X,''===== '',a',i4,',i4,'' ====='')')") &
            level*3, char_len1
       write (mp,fmt) title1,level
      end if

      end subroutine level_print
! *******************************************************************
      subroutine multigrid_eig_deallocate(grid,info)
! deallocate a grid (at given level between last and first)
      type (mc73_multigrid_eig), pointer :: grid
      integer :: st,info,info65

      call mc65_matrix_destruct(grid%graph,info65)
      if (info65 /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      call mc65_matrix_destruct(grid%p,info65)
      if (info65 /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      deallocate(grid%graph, grid%p, grid%eigvector, grid%row_wgt, stat=st)
      if (st /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      deallocate(grid, stat=st)
      if (st /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      end subroutine multigrid_eig_deallocate

! *****************************************************************
      subroutine multigrid_eig_deallocate_last(grid,info)
! deallocate a grid (at the last level). In this case, the matrix grid%graph
! has not been formed yet
      type (mc73_multigrid_eig), pointer :: grid
      integer, intent(inout)  :: info
      integer :: st,info65

      call mc65_matrix_destruct(grid%p,info65)
      if (info65 /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      deallocate(grid%graph, grid%p, grid%eigvector, grid%row_wgt, stat=st)
      if (st /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      deallocate(grid,stat=st)
      if (st /= 0) then
         info = MC73_ERR_MEMORY_DEALLOC
         return
      end if

      end subroutine multigrid_eig_deallocate_last

! *********************************************************
      subroutine multigrid_eig_deallocate_first(nvtx,n,grid,info)
! deallocate a grid (at the first level). In this case the matrix grid%p
! does not exist
      type (mc73_multigrid_eig) :: grid
      integer, intent(out)  :: info
      integer :: n,nvtx
      integer :: ierr,st

      info = 0
      if (nvtx /= n) then
         call mc65_matrix_destruct(grid%graph,ierr)
         if (ierr /= 0) then
            info = MC73_ERR_MEMORY_DEALLOC
            return
         end if
      end if
! in subroutine fiedler, grid%graph is not allocated but is pointed to the
! finest level graph. so no need to deallocate
      nullify (grid%graph)
      deallocate(grid%eigvector, grid%row_wgt, stat=st)
      if (st /= 0) info = MC73_ERR_MEMORY_DEALLOC

      end subroutine multigrid_eig_deallocate_first

! *********************************************************
      subroutine rayleigh(graph,eigvector,lambda,tol1,rtol,maxit,info,  &
                          plevel,unit,llap)

! the rayleigh quotient iteration algorithm for the
! second smallest eigenvector (held in eigvector)
      type (zd11_type) :: graph
      real (kind = wpr), dimension (:) :: eigvector
! eigenvalue
      real (kind = wpr) :: lambda
! tolerance of rayleigh and symmlq
      real (kind = wpr), intent (in) :: tol1,rtol
      integer, intent(in) :: maxit
      integer, intent(inout) :: info,plevel,unit
      logical :: llap
! work arrays
      real (kind = wpr), dimension (:), allocatable :: w,v1
      integer :: st

      deallocate (w,v1,stat=st)
      allocate (w(7*graph%m),v1(graph%m),stat=st)
      if (st /= 0) info = MC73_ERR_MEMORY_ALLOC
      if (info < 0 )return

      call ray(graph,graph%m,unit,plevel,tol1,rtol,eigvector,v1,w,  &
               maxit,lambda,info,llap)
      deallocate (w,v1,stat=st)
      if (st /= 0) info = MC73_ERR_MEMORY_DEALLOC

      end subroutine rayleigh

! *********************************************************
      subroutine ray(graph,n,unit,level,tol1,rtol,v,v1,w,maxit,theta,  &
                     info,llap)

! this is the Rayleigh Quotient Iteration algorithm
! It takes the graph and an approximation
! of the eigenvector as the input,
! and outputs an even more accurate eigenvector
! written in 1993, converted to F90 12 Aug 1998. Yifan Hu
!
! DO WHILE (rho < tol)
! theta = x^T Lx
! (L - theta*I)*x = x
! x = x/||x^{(k+1)}
! rho = sqrt(||Lx||**2 - x^T Lx)/||x||
! END DO
! (L-theta*I)*x = x solved using symmlq algorithm

      real (kind = wpr), parameter :: zero = 0.0_wpr
      type (zd11_type) :: graph
! the dimension of the matrix
      integer, intent (in) :: n
! output channel
      integer, intent (in) :: unit
! max. number of iterations
      integer, intent (in) :: maxit
! warning flag
      integer, intent (inout) :: info
      logical :: llap
! the eigenvector
      real (kind = wpr) :: v(n)
! Work arrays
      real (kind = wpr) :: v1(n),w(7*n)

! eigenvalue
      real (kind = wpr) :: theta,rho,tol1,vnorm
      real (kind = wpr) :: dnrm2

! machine precision and termination tolerance of
! symmlq and residual norm

      real (kind = wpr), intent (in) :: rtol
      real (kind = wpr) :: eps,rnorm,w1

      integer :: i,ir,i1,i2,i3,i4,i5,i6,itnlim,itn,istop
      integer :: level,nout

      info = 0
! deflate v: v = v - (eTv/eTe)e
      w1 = sum(v) / real(n,wpr)
      v = v - w1

! v=x/||x||
      vnorm = dnrm2(n,v,1)
      v = v / vnorm

! Compute L*eigenvector ie v1 = L*v
      call lxv(v,v1,graph,llap)
      vnorm = dnrm2(n,v1,1)

! theta = (v^T)*L*v
      theta = dot_product(v1,v)

! Initialise iteration count
      ir = 0

   10 continue
      ir = ir+1
      if (unit >=0 .and. level > 2) write (unit,'(/a,i5,a)') &
         '------- start the ',ir, '  Rayleigh iteration -----'

! solve (L-theta*I)v1 = v using symmlq algorithm

! divide up workspace
      i1 = 1
      i2 = i1 + n
      i3 = i2 + n
      i4 = i3 + n
      i5 = i4 + n
      i6 = i5 + n

! output channel (if nout=-1 then switched off)
      nout = -1

! machine precision (use Fortan 90 numerical inquiry function epsilon)
      eps = epsilon(zero)

! iteration limit
      itnlim = n

! solve (L - theta*I)v1 = v with starting point 0
      call symmlq(n,v,w(i1),w(i2),w(i3),w(i4),w(i5),w(i6), &
                  theta,nout,itnlim,eps,rtol, &
                  istop,itn,rnorm,graph,llap)

      if (unit >=0 .and. level > 2) then
!        write (unit,*) 'after ',itn,' iterations of symmlq'
!        write (unit,*) 'residual norm = ',rnorm,' istop = ',istop
!        write (unit,*) 'shift = ',theta
        write (unit,'(a,i6)') ' No. iterations of symmlq = ',itn
      end if
! if istop > 3 then symmlq has not worked and we are in trouble.
! Can this happen?
      if (istop > 3) return

! the solution of (L-theta*I)x = v
      v1(1:n) = v(1:n)
      v(1:n)  = w(i5:i5+n-1)

! v = x/||x||
! In the unlikely event that v = 0 so that vnorm = 0, we will have to
! return with v unchanged
      vnorm = dnrm2(n,v,1)
      if (vnorm == zero) then
         v = v1
         return
      end if
      v = v/vnorm

! Lv
! call lxv(v,v1,graph)
! to save one matrix multiply vector
! calculate v1 using v1 = theta*v + oldv/||x||
      do i = 1,n
         v1(i) = theta*v(i) + v1(i) / vnorm
      end do

! theta = (v^T)*L*v
      theta= dot_product(v1,v)

! rho = sqrt{(Lv)*(Lv)-theta*theta}
! (This is from Barnard and Simon 1994)
      rho = dot_product(v1,v1)
      rho = sqrt(abs(rho - theta*theta))

! Not sure why Yifan has the following ... we comment it out and use
! Barnard and Simon
!      rho = 1 / vnorm

!  if (unit >= 0 .and. level > 2) then
!      write (unit, &
!     '(' Rayleigh tol===== ',g12.4,' 1/xnorm = ',g12.4)') &
!      rho,1/vnorm
!  end if

! see if we have got a good approximation to the eigenvalue
! or if we have exceeded max. no. of Rayleigh quotient iterations
      if (rho > tol1 .and. ir < maxit) go to 10

!  if (unit >= 0 .and. level > 2) then
!     write (unit,*) 'ir = ',ir, &
!      ' Rayleigh tol===== ',rho,' 1/xnorm= ',1/vnorm, &
!      ' lambda = ',theta
!  end if
      if (rho > tol1 .and. ir >= maxit) info = 4

      end subroutine ray

! *********************************************************
      subroutine symmlq( n, b, r1, r2, v, w, x, y, shift, &
                         nout, itnlim, eps, rtol, &
                         istop, itn, rnorm, graph, llap )
!  a routine from Michael A. Saunders. Converted to F90 12 Aug 1998.
      type (zd11_type) :: graph
      real(wpr), parameter :: zero = 0.0_wpr
      real(wpr), parameter :: one = 1.0_wpr
      integer :: n, nout, itnlim, istop, itn
      real (kind = wpr) :: shift, eps, rtol, anorm, acond, rnorm
      real (kind = wpr) :: b(n), r1(n), r2(n), v(n), w(n), x(n), y(n)
      logical :: llap
! symmlq  is designed to solve the system of linear equations
!
!            A*x = b
!
! where  A  is an  n x n  symmetric matrix and  b  is a given vector.
! The matrix  A  is not required to be positive definite.
! (If  A  is known to be definite, the method of conjugate gradients
! may be used -- it will require about the same number of iterations
! as  symmlq  but slightly less work per iteration.)
!
! The matrix  A  is intended to be large and sparse.  It is accessed
! by means of a subroutine call of the form
!
!            call lxv(x,y)
!
! which must return the product  y = A*x  for any given vector  x.
!
!
! More generally,  symmlq  is designed to solve the system
!
!            (A - shift*I)*x = b
!
! where   shift  is a specified scalar value.  If  shift  and  b
! are suitably chosen, the computed vector  x  may approximate an
! (unnormalized) eigenvector of  A,  as in the methods of
! inverse iteration and/or Rayleigh quotient iteration.
! Again, the matrix  (A - shift*I)  need not be positive definite.
! The work per iteration is very slightly less if  shift = 0.
!
!
! Parameters
! ----------
!
! N       input      n, the dimension of the matrix  A.
!
! B(N)    input      The rhs vector  b.
!
! R1(N)   output     Returns the final residual vector,
!                       r = b - (A - shift*I)*x.
!
! R2(N)   workspace
! V(N)    workspace
! W(N)    workspace
!
! X(N)    output     Returns the computed solution  x.
!
! Y(N)    workspace
!
! LXV     external   A subroutine defining the matrix  A.
!                    For a given vector  x,  the statement
!
!                          call lxv( x,y )
!
!                    must return the product  y = A*x
!                    without altering the vector  x.
!
! SHIFT   input      Should be zero if the system  A*x = b  is to
!                    be solved.  Otherwise, it could be an
!                    approximation to an eigenvalue of  A,  such as
!                    the Rayleigh quotient  b(t)*A*b/(b(t)*b)
!                    corresponding to the vector  b.
!                    If  b  is sufficiently like an eigenvector
!                    corresponding to an eigenvalue near SHIFT,
!                    then the computed x may have very large
!                    components.  When normalized,  x may be
!                    closer to an eigenvector than b.
!
! NOUT    input      The output file number.  If positive,
!                    a summary will be printed on nout  NOUT.
!
! ITNLIM  input      An upper limit on the number of iterations.
!
! EPS     input      The machine precision.
!
! RTOL    input      A user-specified tolerance.  symmlq terminates
!                    if it appears that  norm(Rbar)  is smaller than
!                    RTOL * norm(Abar) * norm(y),  where
!                    Abar  is the transformed matrix operator
!
!                      Abar = P * (A - shift*I) * P
!
!                    and  Rbar  is the transformed residual vector
!
!                      Rbar = P * ( b - (A - shift*I)*x ).
!
!                    If  shift = 0  and PRECON = .false., symmlq
!                    terminates if  norm(b - A*x)  is smaller than
!                    RTOL * norm(A) * norm(x).
!
! ISTOP   output     An integer giving the reason for termination...
!
!           0        b = 0,  so the exact solution is  x = 0.
!                    No iterations were performed.
!
!           1        Norm(Rbar)  appears to be less than
!                    the value  RTOL * norm(Abar) * norm(y).
!                    The solution in  X  should be acceptable.
!
!           2        Norm(Rbar)  appears to be less than
!                    the value   EPS * norm(Abar) * norm(y).
!                    This means that the residual is as small as
!                    seems reasonable on this machine.
!
!           3        Norm(Abar)*norm(y)  exceeds  norm(b)/EPS,
!                    which should indicate that  x  has essentially
!                    converged to an eigenvector of  A
!                    corresponding to the eigenvalue  shift.
!
!           4        ACOND  (see below)  has exceeded  0.1/EPS,  so
!                    the matrix  Abar  must be very ill-conditioned.
!                    X  may not contain an acceptable solution.
!
!           5        The iteration limit was reached before any of
!                    the previous criteria were satisfied.
!
!           6        An inner product of the form  x(t)*M**(-1)*x
!                    was not positive, so the preconditioning matrix
!                    M  does not appear to be positive definite.
!                    X  will not contain an acceptable solution.
!!!! we have M = I so should not get this message
!
! ITN     output     The number of iterations performed.
!
! ANORM   output     An estimate of the norm of the matrix operator
!                    Abar = P*(A - shift*I)*P,  where P = M**(-1/2).
!
! ACOND   output     An estimate of the condition of  Abar  above.
!                    This will usually be a substantial
!                    under-estimate of the true condition.
!
! RNORM   output     The norm of the final residual vector,
!                       b - (A - shift*I)*x.
!
! XNORM   output     The norm of the final solution vector  x.
!
!
!
! To change precision
! -------------------
!
! Alter the words
!        DOUBLE PRECISION,
!        daxpy, dcopy, ddot, dnrm2
! throughout.
!
! Also make sure  EPS  is set correctly in the calling program.
! ------------------------------------------------------------------
!
!
! This routine is an implementation of the algorithm described in
! the following references:
!
! C.C. Paige and M.A. Saunders,  Solution of Sparse Indefinite
!      Systems of Linear Equations,
!      SIAM J. Numer. Anal. 12, 4, September 1975, pp. 617-629.
!
! J.G. Lewis,  Algorithms for Sparse Matrix Eigenvalue Problems,
!      Report STAN-CS-77-595, Computer Science Department,
!      Stanford University, Stanford, California, March 1977.
!
! Applications of symmlq and the theory of preconditioning
! are described in the following references:
!
! D.B. Szyld and O.B. Widlund,  Applications of Conjugate Gradient
!      Type Methods to Eigenvalue Calculations,
!      in R. Vichnevetsky and R.S. Steplman (editors),
!      Advances in Computer Methods for Partial Differential
!      Equations -- III, IMACS, 1979, 167-173.
!
! D.B. Szyld,  A Two-level Iterative Method for Large Sparse
!      Generalized Eigenvalue Calculations,
!      Ph. D. dissertation, Department of Mathematics,
!      New York University, New York, October 1983.
! ------------------------------------------------------------------
!
!            Michael A. Saunders
!            Department of Operations Research
!            Stanford University
! ------------------------------------------------------------------
!
!
! Subroutines and functions
!
! USER       lxv
! BLAS       daxpy, dcopy, ddot , dnrm2
!

! Functions and local variables

      real (kind = wpr) :: ddot, dnrm2
      integer  :: i
      real (kind = wpr) :: s,t,z,b1,cs,sn,alfa,beta,dbar,diag,    &
                          epsr,epsx,gbar,gmax,gmin,         &
                          oldb,rhs1,rhs2,x1cg,zbar,         &
                          beta1,bstep,delta,denom,epsln,gamma,   &
                          tnorm,ynorm,cgnorm,qrnorm,snprod,ynorm2
      intrinsic abs,max,min,mod,sqrt

      external daxpy,dcopy,ddot,dnrm2

!     ------------------------------------------------------------------
!
!
!     Print heading and initialize.
!
      IF (nout >= 0) write (NOUT, 1000) N,  SHIFT, NOUT, ITNLIM, EPS, RTOL
      istop  = 0
      itn    = 0
      anorm  = zero
      acond  = zero
      rnorm  = zero
!      xnorm  = zero
      ynorm  = zero

      x = zero
      w = zero

! Set up  v,  the first vector in the Lanczos sequence.

     call dcopy ( n,b,1,y,1 )
     call dcopy ( n,b,1,r1,1 )
     b1     = y(1)
     beta1  = ddot( n,r1,1,y,1 )
     if (beta1 < zero) istop = 6
     if (beta1 <= zero) go to 900
     beta1  = sqrt( beta1 )
     s  = one/beta1

     v  = s*y

! Set up  y  for the second Lanczos vector.

     call lxv(v,y,graph,llap)
     call daxpy ( n,(-shift),v,1,y,1 )
     alfa   = ddot( n,v,1,y,1 )
     call daxpy ( n,(-alfa/beta1),r1,1,y,1 )

!  Make sure  r2  will be orthogonal to the first  v.

     z      = ddot( n,v,1,y,1 )
     s      = ddot( n,v,1,v,1 )
     call daxpy ( n,(-z/s),v,1,y,1 )
     call dcopy ( n,y,1,r2,1 )
     oldb   = beta1
     beta   = ddot( n,r2,1,y,1 )
      if (beta < zero) then
         istop = 6
         go to 900
      end if

      beta = sqrt( beta )
! Cause termination (later) if beta is essentially zero.
      if (beta <= eps) then
         istop = -1
      end if


! See if the local reorthogonalization achieved anything.

     denom  = sqrt(s) * dnrm2( n,r2,1 ) + eps
     s      = z/denom
     t      = ddot( n,v,1,r2,1 )/denom
     if (nout > 0) write (nout, 1100) beta1,alfa
     if (nout > 0) write (nout, 1120) s,t

! initialize other quantities.

     cgnorm = beta1
     gbar   = alfa
     dbar   = beta
     rhs1   = beta1
     rhs2   = zero
     bstep  = zero
     snprod = one
     tnorm  = alfa**2
!!!     tnorm  = alfa**2 + beta**2
     ynorm2 = zero
     gmax   = zero
     gmin   = one/eps

!     ------------------------------------------------------------------
!     Main iteration loop.
!     ------------------------------------------------------------------

! Estimate various norms and test for convergence.

300  anorm  = sqrt(tnorm)
     ynorm  = sqrt(ynorm2)
     epsx   = anorm*ynorm*eps
!
! Yifan Hu: do not use such loose tolerance
     epsr   = anorm*ynorm*rtol
!      epsr   = anorm*rtol
     diag   = gbar
! Yifan Hu On some examples EPSA = 0    IF (DIAG == ZERO) DIAG = EPSA
     if (diag == zero) diag = eps

     qrnorm = snprod*beta1
     cgnorm = qrnorm*beta / abs(diag)

!     Estimate  COND(A).
!     In this version we look at the diagonals of  L  in the
!     factorization of the tridiagonal matrix,  T = L*Q.

     denom  = min( gmin, abs(diag) )
     acond  = gmax/denom


! See if any of the stopping criteria are satisfied.
! In rare cases, istop is already -1 from above (Abar = const * I).

     if (istop == 0) then
         if (itn    >= itnlim ) istop = 5
         if (acond  >= 0.1/eps) istop = 4
         if (epsx   >= beta1  ) istop = 3
         if (cgnorm <= epsx   ) istop = 2
         if (cgnorm <= epsr   ) istop = 1
     end if
!     ==================================================================

! See if it is time to print something.

     if (nout <=  0)          go to 600
     if (n    <= 40)          go to 400
     if (itn  <= 10)          go to 400
     if (itn  >= itnlim - 10) go to 400
     if (mod(itn,10)  ==   0) go to 400
     if (cgnorm <= 10.0*epsx) go to 400
     if (cgnorm <= 10.0*epsr) go to 400
     if (acond  >= 0.01/eps ) go to 400
     if (istop  /= 0)         go to 400
     go to 600

! Print a line for this iteration.

 400 zbar   = rhs1/diag
     z      = (snprod*zbar + bstep)/beta1
     x1cg   = x(1) + w(1)*zbar + b1*z

     if (itn == 0) write (nout, 1200)
     write (nout, 1300) itn,x1cg,cgnorm,bstep,anorm,acond
     if (mod(itn,10) == 0) write (nout, 1500)
!     ==================================================================

!     Obtain the current Lanczos vector  V = (1/BETA)*Y
!     and set up  Y  for the next iteration.

 600 if (istop /= 0) go to 800
     s = one/beta

     v = s*y

     call lxv(v,y,graph,llap)
     call daxpy ( n,(-shift),v,1,y,1 )
     call daxpy ( n,(-beta/oldb),r1,1,y,1 )
     alfa   = ddot( n,v,1,y,1 )
     tnorm  = tnorm + (alfa**2) + 2.0*(beta**2)
     call daxpy ( n,(-alfa/beta),r2,1,y,1 )
     call dcopy ( n,r2,1,r1,1 )
     call dcopy ( n,y,1,r2,1 )
     oldb   = beta
     beta   = ddot( n,r2,1,y,1 )
     if (beta <= zero) then
        istop = 6
        go to 800
     end if
     beta   = sqrt( beta )
!!!     tnorm  = tnorm  +  alfa**2  +  oldb**2  +  beta**2

!     Compute the next plane rotation for  Q.

     gamma  = sqrt(gbar**2 + oldb**2)
     cs     = gbar/gamma
     sn     = oldb/gamma
     delta  = cs*dbar + sn*alfa
     gbar   = sn*dbar - cs*alfa
     epsln  = sn*beta
     dbar   =         - cs*beta

!     Update  X.

     z      = rhs1/gamma
     s      = z*cs
     t      = z*sn

     do i = 1, n
        x(i)  = (w(i)*s  +  v(i)*t)  +  x(i)
        w(i)  =  w(i)*sn -  v(i)*cs
     end do

!     Accumulate the step along the direction  B,
!     and go round again.

     bstep  = snprod*cs*z + bstep
     snprod = snprod*sn
     gmax   = max( gmax, gamma )
     gmin   = min( gmin, gamma )
     ynorm2 = z**2 + ynorm2
     rhs1   = rhs2 - delta*z
     rhs2   =      - epsln*z
     itn    = itn  + 1
     go to 300

!     ------------------------------------------------------------------
!     End of main iteration loop.
!     ------------------------------------------------------------------

! Move to the CG point.
! (In this version of SYMMLQ, we never stop at an LQ point.)

800  zbar   = rhs1/diag
     bstep  = snprod*zbar + bstep
     ynorm  = sqrt(ynorm2 + zbar**2)
     rnorm  = cgnorm
     call daxpy ( n,zbar,w,1,x,1 )

! Add the step along  B.
     bstep  = bstep/beta1
     call dcopy ( n,b,1,y,1 )
     call daxpy ( n,bstep,y,1,x,1 )

! Compute the final residual,  r1 = b - (a - shift*i)*x.
!     call lxv(x,y,graph)
!     call daxpy ( n,(-shift),x,1,y,1 )
!     r1 = b - y
!     rnorm  = dnrm2( n,r1,1 )
!     write (6,*) 'rnorm,cgnorm',rnorm,cgnorm
!     xnorm  = dnrm2( n,x,1 )
!
!     ==================================================================
!     Display final status.
!     ==================================================================
 900 if (nout  <= 0) go to 950
!     write (nout, 2000) istop, itn, anorm, acond, rnorm, xnorm
     write (nout, 2000) istop, itn, anorm, acond, rnorm
     if (istop == 0) write (nout, 3000)
     if (istop == 1) write (nout, 3100)
     if (istop == 2) write (nout, 3200)
     if (istop == 3) write (nout, 3300)
     if (istop == 4) write (nout, 3400)
     if (istop == 5) write (nout, 3500)
     if (istop == 6) write (nout, 3600)
 950 return

!     ------------------------------------------------------------------
1000 format( ' symmlq.            solution of symmetric  ax = b'&
       & , ' n    =', i7, 6x,  'shift =', es23.14&
       & , ' nout =', i7, 6x, 'itnlim =', i6, 6x, 'eps   =', es11.2, 5x,&
       &   ' rtol =', es11.2)
1100 format( ' beta1  =', es12.2 / ' alpha1 =', es12.2)
1120 format( ' (v1,v2) before and after ', es15.2&
       &       , ' local reorthogonalization', es15.2 )
1200 format( 5x, 'itn', 7x, 'x1(cg)',&
       &    10x, 'norm(r)', 5x, 'bstep', 7x, 'norm(a)', 3x, 'cond(a)')
1300 format(i8, es19.10, es11.2, es14.5, 2(es10.2))
1500 format(1x)
2000 format( ' exit symmlq.',    14x, 'istop =', i3, 18x, 'itn =', i8&
       &   , ' anorm =', es13.5, 6x, 'acond =', es13.5, 5x,&
       &     ' rnorm =', es13.5)
!       &     ' rnorm =', es13.5, 6x, 'xnorm =', es13.5)
3000 format( ' The exact solution is  x = 0.')
3100 format( ' Requested accuracy achieved, as determined by  rtol.')
3200 format( ' Reasonable accuracy achieved, given  eps.')
3300 format( ' x  has converged to an eigenvector.')
3400 format( ' acond  has exceeded  0.1/eps.')
3500 format( ' The iteration limit was reached.')
3600 format( ' istop = 6')
!     ------------------------------------------------------------------
      end subroutine symmlq

! ***************************************************

      subroutine lxv(x,y,graph,llap)
! subroutine to multiply the Laplacian by a vector
! y = L*x
      type (zd11_type) :: graph
      real (kind = wpr), intent (in), dimension (*) :: x
      real (kind = wpr), intent (out), dimension (*) :: y
      logical :: llap
      integer :: i,j
      integer,  dimension (:), pointer :: neighbor
      real (kind = wpr),  dimension (:), pointer :: val

      if (llap) then
         do i = 1,graph%m
            call mc65_matrix_getrow(graph,i,neighbor)
            call mc65_matrix_getrowval(graph,i,val)
            y(i) = 0.0_wpr
            do j = 1,size(neighbor)
              y(i) = y(i) + val(j)*x(neighbor(j))
              y(i) = y(i) - val(j)*x(i)
            end do
         end do
      else
         do i = 1,graph%m
            call mc65_matrix_getrow(graph,i,neighbor)
            y(i) = 0.0_wpr
            do j = 1,size(neighbor)
              y(i) = y(i) + x(neighbor(j))
            end do
            y(i) = y(i) - size(neighbor)*x(i)
         end do
      end if

    end subroutine lxv
! ***************************************************
      subroutine lanczos(graph,n,eigv,lambda,mlancz,tol,unit,level,ierr, &
                         llap)

! this is the Lanczos algorithm that calculates the
! eigenvector corresponding to the smallest positive eigenvalue of
! the Laplacian matrix of a graph. The product
! of the matrix with a vector is given by the
! external subroutine lxv. To save memory,
! intermediary vectors are not saved, the Lanczos algorithm
! is run twice.
!
! First written 1993
! converted into F90 on 12/8/1998.
! Yifan Hu, Daresbury Laboratory. Y.F.Hu@dl.ac.uk

! accuracy for the eigenvector calculation
  real (kind = wpr) :: tol

      type (zd11_type) :: graph
! dimension of matrix (number of vertices of the graph)
      integer :: n
! error flag
      integer :: ierr
      logical :: llap

! output channel
      integer :: unit

! print level
      integer :: level

! the eigenvector and eigenvalue
      real (kind = wpr), intent (out) :: eigv(n),lambda

! working variable
      real (kind = wpr) :: w1

! mlancz: maximum number of Lanczos vectors used
      integer :: mlancz,mlancz1
! maximum number of lanczos vector, not
! exceeding the size of the problem
      mlancz1 = mlancz
      if (mlancz1 > n) mlancz1 = n

      call egvector(mlancz1,n,unit,level,tol,eigv,lambda,graph,ierr,llap)
      if (ierr < 0) return

!  normalize
      w1 = dot_product(eigv,eigv)
      w1 = sqrt(w1)
      eigv(1:n) = eigv/w1

      end subroutine lanczos

! ***************************************************
      subroutine egvector(mlancz,n,unit,level,tol,eigv,lambda,graph, &
                          ierr,llap)

! mlancz: maximum number of lanczos vectors used
      real (kind = wpr), parameter :: one = 1.0_wpr
      real (kind = wpr), parameter :: two = 2.0_wpr
      real (kind = wpr), parameter :: twelve = 12.0_wpr
      real (kind = wpr), parameter :: zero = 0.0_wpr

! unit: print unit
! level : print level
! n size of the system
      integer, intent (in) :: unit,level,mlancz,n
! tol:  tolerance for the Lanczos process
      real (kind = wpr), intent (in) :: tol
! eigv: the eigenvector correspomds to the seconds largest eigenvalue of -L
      real (kind = wpr), intent (out) :: eigv(n)
! lambda: the second largest eigenvalue of -L, the Laplacian
      real (kind = wpr), intent (out) :: lambda
      type (zd11_type) :: graph
! error flag
      integer,  intent (inout) :: ierr
      logical :: llap

! working variables
      real (kind = wpr) :: dwork,dwork1,hs1,hs2
! vectors used in the Lanczos process
      real (kind = wpr), allocatable, dimension (:) :: v,u,r
! the coefficient and eggenvectors of the tridigonal systems
      real (kind = wpr), allocatable, dimension (:) :: alpha,beta,eig

      integer st
! DONE: where this is the first or the second cycle of the Lanczos
      logical done

! working variables
      integer :: nwork,i,j,k,jjj,seed
      real (kind = wpr) :: w
! random number generator
      real (kind = wpr) :: fa14ad

      external fa14id,fa14ad

      call fa14id(seed)

      allocate(v(n),u(n),r(n),stat=st)
      if (st /= 0) then
        ierr = MC73_ERR_MEMORY_ALLOC
        return
      end if

      allocate(alpha(mlancz),beta(mlancz),eig(mlancz),stat=st)
      if (st /= 0) then
        ierr = MC73_ERR_MEMORY_ALLOC
        return
      end if

      done = .false.

! calculate the elements for the Householder vector
   1  continue
      nwork = n
      dwork = real(nwork,wpr)
      dwork1 = sqrt(dwork + sqrt(dwork))
      hs2 = one/dwork1
      hs1 = hs2*(one + sqrt(dwork))

      dwork1 = sqrt(dwork*(dwork + one)*(dwork - one)/twelve)
      dwork = (dwork + one)/two

! initial eigenvector
       eigv = zero

! initial Lanczos vector
      do i = 1,n
         r(i) = (i - dwork)/dwork1
         v(i) = r(i)
      end do

   7  continue
! Compute u = L*r
      call lxv(r,u,graph,llap)

      do j = 1,n
        if (.not.done)  alpha(j) = dot_product(u,v)
        r = u - alpha(j)*v

! orthogonalize the lanzcos vector r
        dwork = zero
        do k = 2,n
           dwork = dwork + r(k)
        end do
        dwork = dwork*hs2 + r(1)*hs1
        dwork = dwork*hs2
        dwork1 = zero
        do k = 2,n
           r(k) = r(k) - dwork
           dwork1 = dwork1 + r(k)
        end do
        r(1) = -hs1*hs2*dwork1
        dwork1 = dwork1*hs2*hs2
        if (.not.done)  beta(j) = r(1)*r(1)

        do k = 2,n
           r(k) = r(k) - dwork1
           if (.not.done)  beta(j) = beta(j) + r(k)*r(k)
         end do

        if (.not.done)  beta(j) = sqrt(beta(j))

! if beta is already very small in the first iteration
! then eigvector found, but better perturb it in case
! it is not the one corresponding to the second largest
! eigenvalue

!!    if (j == 1 .and. beta(j) < tiny(1.0_wpr)) then
        if (j == 1 .and. beta(j) < epsilon(zero)) then
           if (unit >= 0 .and. level > 1) then
              write (unit,'(a)') &
                ' Eigenvector found on first iteration'
           end if
           do k = 1,n
              w = fa14ad(seed,-1)
              r(k) = 1 - w*w
           end do

! normalise
           dwork = dot_product(r,r)
           dwork = sqrt(dwork)
           v = r/dwork
           r = v
           go to 7
        end if

        if (j >= 2 .and. (.not.done)) then
           call tridiagonal(j,mlancz,alpha,beta,eig,lambda)
!  normalize the eig
           dwork = dot_product(eig(1:j),eig(1:j))
           dwork = sqrt(dwork)
           if (unit >= 0 .and. level > 1) write (unit,*)' j = ',j, &
             ' beta*eig = ',beta(j)*eig(j)/dwork
           if ((abs(beta(j)*eig(j)/dwork) < tol) &
             .or. (j >= mlancz) .or. (j >= n)) then
! converged or too many iters.
              if (unit >= 0 .and. level > 1) then
                 write (unit,*)' j = ',j,                  &
               ' beta*eig = ',beta(j)*eig(j)/dwork,        &
               ' lambda = ',lambda
              end if
! eigenvalue and eigenvector of tridiagonal found
              done = .true.
! jjj is the number of lanczos vectors needed
              jjj = j
! recalculate lanczos vector.
              go to 1
           end if
        end if

! r is the new lanczos vector.
! Be careful: on some small test problems, found beta(j) = 0
        if (beta(j) < epsilon(zero)) beta(j) = epsilon(zero)
        r = r/beta(j)

! Compute u = L*r
        call lxv(r,u,graph,llap)
        u = u - beta(j)*v

! if done then get the eigenvector of the original system
        if (done) eigv = eigv + v*eig(j)

! update lanczos vector
        v = r

! get out if reached the same number of lanczos vectors
        if (done) then
          if (j == jjj) exit
        end if

! Compute u = L*r
!   call lxv(r,u,graph)
!   u = u - beta(j)*v

     end do
     deallocate(u,v,r,alpha,beta,eig,stat = st)
     if (st /= 0) ierr = MC73_ERR_MEMORY_DEALLOC

     end subroutine egvector

! *******************************************************

      subroutine tridiagonal(j,m,alpha,beta,eig,lambda)

! subroutine to get the eigenvalue of the tridiagonal system

      real (kind = wpr), parameter :: one = 1.0_wpr
      real (kind = wpr), parameter :: two = 2.0_wpr
      real (kind = wpr), parameter :: four = 4.0_wpr
      real (kind = wpr), parameter :: zero = 0.0_wpr


! the order of the tridiagonal
       integer :: j,m
! the eigenvector
       real (kind = wpr) :: eig(m)
! coefficients of the tridigonal system
       real (kind = wpr), intent (in) :: alpha(m),beta(m)
! the eigenvalue (on entry it is that of the j-1 th order system,
! on exit it is that of the jth order system)
      real (kind = wpr) :: lambda
! up and low bound within which eigenvalue is to be found,
! the eigenvalue and a working variable
      real (kind = wpr) :: upbound,lowbound,root, work

      integer :: i

      if (j == 2) then
         work = alpha(1) - alpha(2)
         work = work*work + beta(1)*beta(1)*four
         work = sqrt(work)
         lambda = (alpha(1) + alpha(2) + work)/two
         root = lambda
         go to 10
      end if
      upbound = zero
      lowbound = lambda
      call solution(j,m,upbound,lowbound,alpha,beta,root)
      lambda = root
   10 continue
! calculate the eigenvector of the tridigonal system
      eig(1) = one
      eig(2) = eig(1)*(root - alpha(1))/beta(1)
      do i = 3,j
         eig(i) = -(beta(i-2)*eig(i-2) + &
                    (alpha(i-1) - root)*eig(i-1))/beta(i-1)
      end do

      end subroutine tridiagonal
! *******************************************************

      subroutine solution(j,m,upbound,lowbound,alpha,beta,root)
! use the bisection algorithm to find the eigenvalue
! of the j-th order tridiagonal system, using
! the upper bound 0 and lower bound the second largest
! eigenvalue of the (j-1)-th tridiagonal system.
! (given that the eigenvalues of the two systems inter-leave)
! the coefficient of the tridiagonal system
      real (kind = wpr), parameter :: zero = 0.0_wpr
      integer :: j,m
      real (kind = wpr), intent (in) :: alpha(m),beta(m)
! initial up and low bound, the eigenvalue to be found and the on going
! up and low bound
      real (kind = wpr) :: upbound,lowbound,root,up,low,root1
! the ratio of the jth order determinant and the (j-1)th order
! determinant for the tridiagonal system (T-root*I)
      real (kind = wpr) :: f
      integer :: flag,nflag
      real (kind = wpr) :: eps

      eps = epsilon(zero)
      nflag = 0
      root1 = zero
      up = upbound
      low = lowbound
10    root = (up + low)/2.0_wpr
!  Hence we have added this check since we found that code
!  could get 'stuck' in a loop with up-low < eps*100.0_wpr never satisfied
! (only found to get stuck in single precision code)
      if (abs(root-root1) < eps) go to 20
      root1 = root
30    flag = 0
      call delta(j,m,alpha,beta,root,f,flag)
      if (flag == 1) then
         root = (root + up)/2.0_wpr
         nflag = nflag + 1
         if (nflag > 100) then
            go to 15
         end if
         go to 30
      end if
15    continue
      if ((abs(f) < 1.0d-5) .or. (up-low < eps*100.0_wpr)) go to 20
      if (f > 0) then
         low = root
      else
         up = root
      end if
      go to 10
20    continue

      end subroutine solution

! ******************************************************************
      subroutine delta(j,m,alpha,beta,root,f,flag)
      real (kind = wpr), parameter :: small = 1.0e-10_wpr
! calculate the  ratio of the jth order determinant and the (j-1)th order
! determinant for the tridiagonal system (T-root*I)
! Order of the system
      integer :: j,m
! the coefficient of the tridiagonal system
      real (kind = wpr), intent (in) :: alpha(m),beta(m)
! the eigenvalue and the ratio of determinants
      real (kind = wpr) :: root,f
! flag = 0 if normal, flag = 1 if one of the intermediate
! determinant is zero (very small)
      integer :: flag,i

      f = alpha(1) - root
      if (abs(f) < small) then
         flag = 1
         return
      end if
      do  i = 2,j
         if (abs(f) < small) then
            flag = 1
            exit
         end if
         f = alpha(i) - root - (beta(i-1)*beta(i-1))/f
      end do

     end subroutine delta

!*************************************************
      subroutine galerkin_graph(matrix,p,cmatrix,info)

! Given matrix on fine grid and a prolongation operator p,
! find the coarse matrix R*A*P

! matrix: fine grid matrix
      type (zd11_type), intent (in) :: matrix
! p: prolongation operator
      type (zd11_type), intent (in) :: p
! cmatrix: coarse grid matrix
      type (zd11_type), intent (out) :: cmatrix
! r: restriction operator
      type (zd11_type) :: r

! nvtx,cnvtx: size of fine and coarse grid
      integer :: nvtx,cnvtx
      integer :: nz

      integer, intent (inout) :: info
      integer :: info65

      call mc65_matrix_transpose(p,r,info65)
      if (info65 < 0) then
         info = info65
         return
      end if

      nvtx = matrix%n
      cnvtx = p%n

! get the size of the coarse matrix first
      call galerkin_graph_rap_size(nvtx,cnvtx,nz,              &
                 p%ptr(nvtx+1)-1,p%col,p%ptr,                  &
                 matrix%ptr(nvtx+1)-1,matrix%col,matrix%ptr,   &
                 r%ptr(cnvtx+1)-1,r%col,r%ptr,info)
      if (info < 0) return

      call mc65_matrix_construct(cmatrix,cnvtx,nz,info65)
      if (info65 < 0) then
         info = info65
         return
      end if

      call galerkin_graph_rap(nvtx,cnvtx,                              &
                p%ptr(nvtx+1)-1,p%val,p%col,p%ptr,                     &
                matrix%ptr(nvtx+1)-1,matrix%val,matrix%col,matrix%ptr, &
                r%ptr(cnvtx+1)-1,r%val,r%col,r%ptr,                    &
                nz,cmatrix%val,cmatrix%col,cmatrix%ptr,info)
      if (info < 0) return

      call mc65_matrix_destruct(r,info65)
      if (info65 < 0) then
         info = info65
         return
      end if

      end subroutine galerkin_graph

!*************************************************

      subroutine galerkin_graph_rap_size(nvtx,cnvtx,nz,&
                   nzp,pcol,pptr,nzaa,acol,aptr,nzr,rcol,rptr,info)
! get the number of nonzeros in R*A*P
! nvtx: size of aa matrix
! cnvtx: size of ca matrix
      integer,  intent (in) :: nvtx,cnvtx
! nz: number of nonzeros in R*A*P
      integer,  intent (out) :: nz

! P: matrix
      integer,  intent (in) :: nzp
      integer,  intent (in), dimension (nzp) :: pcol
      integer,  intent (in), dimension (nvtx+1) :: pptr
! aa: matrix
      integer,  intent (in) :: nzaa
      integer,  intent (in), dimension (nzaa) :: acol
      integer,  intent (in), dimension (nvtx+1) :: aptr
! R: matrix
      integer,  intent (in) :: nzr
      integer,  intent (in), dimension (nzr) :: rcol
      integer,  intent (in), dimension (cnvtx+1) :: rptr

! mask: masking array to see if an entry has been seen before
      integer,  dimension (:), allocatable :: mask
! i,j,k: loop index
      integer :: i,j,k
! nz: number of nonzeros so far in ca
      integer :: nz1
! various neighbors
      integer :: neigh,neighneigh

! col: column index of a row of r*matrix
      integer,  dimension (:), allocatable :: col

      integer, intent (inout) :: info
      integer :: st

      allocate(mask(nvtx),col(nvtx),stat=st)
      if (st /= 0) then
        info = MC73_ERR_MEMORY_ALLOC
        return
      end if
      mask = 0
      nz = 0
! loop over coarse grid points
      do i = 1,cnvtx
! first form row i of (r*matrix)
       nz1 = 0
! for each vertex D that restricts to C (including itself).
       do j = rptr(i),rptr(i+1)-1
          neigh = rcol(j)
! find D's neighbor
          do k = aptr(neigh),aptr(neigh+1)-1
             neighneigh = acol(k)
             if (mask(neighneigh) /= i) then
                nz1 = nz1 + 1
                col(nz1) = neighneigh
                mask(neighneigh) = i
             end if
          end do
       end do
! form row i of (r*matrix)*p
       do j = 1,nz1
          neigh = col(j)
          do k = pptr(neigh),pptr(neigh+1)-1
             neighneigh = pcol(k)
             if (mask(neighneigh) /= -i .and. neighneigh /= i) then
                nz = nz + 1
                mask(neighneigh) = -i
             end if
          end do
       end do
     end do
     deallocate(mask,col,stat=st)
     if (st /= 0) info = MC73_ERR_MEMORY_DEALLOC

     end subroutine galerkin_graph_rap_size
! ******************************************************
  subroutine galerkin_graph_rap(nvtx,cnvtx,&
       nzp,pa,pcol,pptr,nzaa,aa,acol,aptr,&
       nzr,ra,rcol,rptr,nzca,ca,ccol,cptr,info)
! multiply R*A*P to get CA
! nvtx: size of aa matrix
! cnvtx: size of ca matrix
    integer,  intent (in) :: nvtx,cnvtx

! p: matrix
    integer,  intent (in) :: nzp
    real (kind = wpr), intent (in), dimension (nzp) :: pa
    integer,  intent (in), dimension (nzp) :: pcol
    integer,  intent (in), dimension (nvtx+1) :: pptr
! aa: matrix
    integer,  intent (in) :: nzaa
    real (kind = wpr), intent (in), dimension (:) :: aa
    integer,  intent (in), dimension (nzaa) :: acol
    integer,  intent (in), dimension (:) :: aptr
! r: matrix
    integer,  intent (in) :: nzr
    real (kind = wpr), intent (in), dimension (nzr) :: ra
    integer,  intent (in), dimension (nzr) :: rcol
    integer,  intent (in), dimension (cnvtx+1) :: rptr
! ca: matrix
    integer,  intent (in) :: nzca
    real (kind = wpr), intent (inout), dimension (nzca) :: ca
    integer,  intent (inout), dimension (nzca) :: ccol
    integer,  intent (inout), dimension (cnvtx+1) :: cptr


! mask: masking array to see if an entry has been seen before
      integer,  dimension (:), allocatable :: mask
! i,j,k,l: loop index
      integer :: i,j,k
! nz: number of nonzeros so far in ca
      integer :: nz,nzz,nz1
! various neighbors
      integer :: neigh,neighneigh
! r_ij: (i,j) element of r
      real (kind = wpr) :: r_ij
! col: column index of a row of r*matrix
! a: values of a row of r*matrix
      integer, dimension (:), allocatable :: col
      real (kind = wpr), dimension (:), allocatable :: a

      integer, intent (inout) :: info
      integer st

      allocate(col(nvtx),a(nvtx),mask(nvtx),stat=st)
      if (st /= 0) then
        info = MC73_ERR_MEMORY_ALLOC
        return
      end if
! now get the entries of the coarse matrix
      cptr(1) = 1
      mask = 0
      nz = 0
! loop over every coarse grid point
      do i = 1,cnvtx
! fist form row i of (r*matrix)
        nz1 = 0
! foreach each vertex D that restricts to C (including itself).
        do j = rptr(i),rptr(i+1)-1
          neigh = rcol(j)
          r_ij = ra(j)
! find D's neighbor
          do k = aptr(neigh),aptr(neigh+1)-1
             neighneigh = acol(k)
             nzz = mask(neighneigh)
             if (nzz == 0) then
                nz1 = nz1 + 1
                col(nz1) = neighneigh
                a(nz1) = r_ij*aa(k)
                mask(neighneigh) = nz1
             else
                a(nzz) = a(nzz) + r_ij*aa(k)
             end if
          end do
        end do
        mask(col(1:nz1)) = 0

! form row i of (r*matrix)*p
        do j = 1,nz1
          neigh = col(j)
          r_ij = a(j)
          do k = pptr(neigh),pptr(neigh+1)-1
             neighneigh = pcol(k)
             if (neighneigh == i) cycle
             nzz = mask(neighneigh)
             if (nzz == 0) then
                nz = nz + 1
                mask(neighneigh) = nz
                ca(nz) = r_ij*pa(k)
                ccol(nz) = neighneigh
             else
                ca(nzz) = ca(nzz) + r_ij*pa(k)
             end if
          end do
        end do

        mask(ccol(cptr(i):nz)) = 0
        cptr(i+1) = nz+1
      end do

      deallocate(col,a,mask,stat=st)
      if (st /= 0) info = MC73_ERR_MEMORY_DEALLOC

      end subroutine galerkin_graph_rap
! ***************************************************
      subroutine print_message(info,unit,message)
! printing error message
! info: is an integer scaler of INTENT (IN).
! It is the information flag
!   whose corresponding error/warning message is to be printed.
! unit: is the unit number
!     the user wish to print the message to.
!     If this number
!     is negative, printing is supressed.
! message: is an OPTIONAL assumed size CHARACTER array
!      of INTENT (IN).
!      It must be set by the user to hold the message
!      to be printed ahead of the error or warning message.
      integer, intent (in) :: info
      integer, intent (in) :: unit
      character (len = *), optional, intent (in) :: message
      integer :: length

      if (unit <= 0) return

      if (info > 0) then
         write (unit,advance = "yes", fmt = "(' Warning: ')")
      else if (info < 0) then
         write (unit,advance = "yes", fmt = "(' Error return: ')")
      end if

      if (present(message)) then
         length = len_trim(message)
         write (unit,advance = "no", fmt = "(a,' : ')") message(1:length)
      end if

      select case (info)
      case (0)
         write (unit,'(a)') ' successful completion'
      case (MC73_ERR_MEMORY_ALLOC)
         write (unit,'(a,i3)') ' memory allocation failure. INFO(1) = ', &
         MC73_ERR_MEMORY_ALLOC
      case (MC73_ERR_MEMORY_DEALLOC)
         write (unit,'(a,i3)') ' memory deallocate failure. INFO(1) = ', &
         MC73_ERR_MEMORY_DEALLOC
      case (MC73_ERR_N_NONPOSITIVE)
         write (unit,'(a,i3)') ' The order N is less than 1. INFO(1) = ', &
         MC73_ERR_N_NONPOSITIVE
      case (MC73_ERR_LIRN_TOOSMALL)
         write (unit,'(a,i3)') ' LIRN is too small. INFO(1) = ', &
         MC73_ERR_LIRN_TOOSMALL
      case (MC73_ERR_A_TOOSMALL)
         write (unit,'(a,i3)') ' Size of A is too small. INFO(1) = ', &
         MC73_ERR_A_TOOSMALL
      case (MC73_ERR_RANGE_IP)
         write (unit,'(a,i3)') &
       ' IP is not monotonically increasing. INFO(1) = ',MC73_ERR_RANGE_IP
      case (MC73_ERR_JOB_WRONG)
         write (unit,'(a,i3)') ' JOB is invalid. INFO(1) = ', &
         MC73_ERR_JOB_WRONG

      case (MC73_WARN_RANGE_IRN)
         write (unit,'(a,i3)') &
       ' Out of range entries found in IRN. INFO(1) = ',MC73_WARN_RANGE_IRN
      case (MC73_WARN_DUP_ENTRY)
         write (unit,'(a,i3)') &
       ' Duplicate entries were found in IRN. INFO(1) = ',MC73_WARN_DUP_ENTRY
      case (MC73_WARN_MAXIT)
         write (unit,'(a,i3)') &
        ' Max. no. Rayleigh Quotient Iterations reached. INFO(1) = ',&
         MC73_WARN_MAXIT
      case (MC73_WARN_MC67)
         write (unit,'(a,i3)') &
        ' Hager did not reduce the profile. INFO(1) = ', MC73_WARN_MC67
      case default
         write (unit,'(a,i10,a)') ' INFO(1) flag has value ', info, &
              ' This is not a recognized value'
      end select
      end subroutine print_message

   end module hsl_mc73_double



