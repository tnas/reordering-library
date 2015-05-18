! COPYRIGHT (c) 2006 Council for the Central Laboratory
!                    of the Research Councils
! This package may be copied and used in any application, provided no
! changes are made to these or any other lines.
! Original date 21 February 2006. Version 1.0.0.
! 6 March 2007 Version 1.1.0. Argument stat made non-optional

MODULE HSL_ZD11_single

!  ==========================
!  Sparse matrix derived type
!  ==========================

  TYPE, PUBLIC :: ZD11_type
    INTEGER :: m, n, ne
    CHARACTER, ALLOCATABLE, DIMENSION(:) :: id, type
    INTEGER, ALLOCATABLE, DIMENSION(:) :: row, col, ptr
    REAL ( KIND( 1.0E+0 ) ), ALLOCATABLE, DIMENSION(:) :: val
  END TYPE

CONTAINS

   SUBROUTINE ZD11_put(array,string,stat)
     CHARACTER, allocatable :: array(:)
     CHARACTER(*), intent(in) ::  string
     INTEGER, intent(OUT) ::  stat

     INTEGER :: i,l

     l = len_trim(string)
     if (allocated(array)) then
        deallocate(array,stat=stat)
        if (stat/=0) return
     end if
     allocate(array(l),stat=stat)
     if (stat/=0) return
     do i = 1, l
       array(i) = string(i:i)
     end do

   END SUBROUTINE ZD11_put

   FUNCTION ZD11_get(array)
     CHARACTER, intent(in):: array(:)
     CHARACTER(size(array)) ::  ZD11_get
! Give the value of array to string.

     integer :: i
     do i = 1, size(array)
        ZD11_get(i:i) = array(i)
     end do

   END FUNCTION ZD11_get

END MODULE HSL_ZD11_single


! COPYRIGHT (c) 1995 Council for the Central Laboratory
!                    of the Research Councils
! Original date 23 March 2001
!  March 2001: threadsafe version of HSL_FA04
!  December 2003. Long lines shortened.

!-*-*-*-*-*-*-*-*  H S L _ F A 1 4   M O D U L E  *-*-*-*-*-*-*-*-*-*-*-*-*-

! 12th July 2004 Version 1.0.0. Version numbering added.

      MODULE HSL_FA14_SINGLE

!  Portable random number generator by Linus Schrange, TOMS 5, 1979, pp 132-138
!  Fortran 95 version by Nick Gould and John Reid, RAL, September 1995

         IMPLICIT NONE

         PRIVATE
         PUBLIC :: FA14_RANDOM_REAL, FA14_RANDOM_INTEGER, FA14_GET_SEED,      &
                   FA14_SET_SEED, FA14_INITIALIZE

!  Define the working precision to be single

         INTEGER, PARAMETER :: working = KIND( 1.0E+0 )

         TYPE, PUBLIC :: FA14_seed
            PRIVATE
            INTEGER :: ix = 65535
         END TYPE
         INTEGER, PARAMETER :: a = 16807, b15 = 32768
         INTEGER, PARAMETER :: b16 = 65536, p = 2147483647
         INTEGER, PARAMETER :: b30 = 1073741824, q = 1073741823

      CONTAINS

!-*-*-*-*-*- F A 1 4 _ R A N D O M _ R E A L  S U B R O U T I N E  *-*-*-*-*-

         SUBROUTINE FA14_RANDOM_REAL ( seed, positive, random )

!  Real random number in the range [0, 1] ( if positive is .TRUE. )
!  or [-1, 1] ( if positive is .FALSE. )

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

         TYPE (FA14_seed), INTENT( INOUT ) :: seed
         LOGICAL, INTENT( IN ) :: positive
         REAL ( KIND = working ), INTENT( OUT ) :: random

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

         INTEGER :: fhi, k, leftlo, xalo, xhi
         REAL ( KIND = working ) :: x
         REAL ( KIND = working ), PARAMETER ::                                &
                one = 1.0_working, two = 2.0_working, rb16 = two ** 16,       &
                big  = one / ( two ** 31 - one ), big2 = two * big

         xhi = seed%ix / b16

! Get 16 lo bits of seed%ix and form lo product

         xalo = ( seed%ix - xhi * b16 ) * a

!  Get 15 hi order bits of lo product

         leftlo = xalo / b16

!  Form the 31 highest bits of full product

         fhi = xhi * a + leftlo

!  Get overflopast 31st bit of full product

         k = fhi / b15

!  Assemble all the parts and presubtract P. The parentheses are essential

         seed%ix = ( ( (xalo-leftlo*b16) - p ) + (fhi-k*b15)*b16 ) + k

!  Add p back in if neccessary

         IF ( seed%ix < 0 ) seed%ix = seed%ix + p

!  Multiply by 1/(2**31-1)

         xhi = seed%ix / b16
         x = FLOAT( xhi ) * rb16 + FLOAT( seed%ix - xhi * b16 )
         IF ( positive ) THEN
            random = x * big
         ELSE
            random = x * big2 - one
         END IF

         END SUBROUTINE FA14_RANDOM_REAL

!-*-*-*-*  F A 1 4 _ R A N D O M _ I N T E G E R   S U B R O U T I N E  *-*-

         SUBROUTINE FA14_RANDOM_INTEGER ( seed, n, random )

!  Integer random number in the range [1,n] if n > 1.
!  Otherwise, the value n is returned

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

         TYPE (FA14_seed), INTENT( INOUT ) :: seed
         INTEGER, INTENT( IN ) :: n
         INTEGER, INTENT( OUT ) :: random

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

         INTEGER :: be1, be2, c, d, f, fhi, g, k, leftlo
         INTEGER :: mhi, mlo, mu, nu, xalo, xhi, xlo

         IF ( n > 1 ) THEN

            xhi = seed%ix / b16

!  Get 16 lo bits of seed%ix and form lo product

            xalo = ( seed%ix - xhi * b16 ) * a

!  Get 15 hi order bits of lo product

            leftlo = xalo / b16

!  Form the 31 highest bits of full product

            fhi = xhi * a + leftlo

!  Get overflopast 31st bit of full product

            k = fhi / b15

!  Assemble all the parts and presubtract P. The parentheses are essential

            seed%ix = (((xalo-leftlo*b16) - p) + (fhi-k*b15) * b16) +  k

!  Add p back in if neccessary

            IF ( seed%ix < 0 ) seed%ix = seed%ix + p

!  Multiply by n and divide by 2**31-1 in integer arithmetic.
!  Split seed%ix and n into hi and lo parts

            xhi = seed%ix / b15 ; xlo = seed%ix - b15 * xhi
            mhi = n / b15 ; mlo = n - b15 * mhi

!  Calculate intermediate product and split into hi and lo parts.
!  Presubtract p

            f = ( xhi * mlo - p ) + xlo * mhi

!  f is > 0 if intermediate product would have overflowed

            IF ( f <= 0 ) THEN
               f = f + p ; be1 = f / b15 ; be2 = f - be1 * b15
            ELSE
               f = f - 1 ; be1 = f / b15 ; be2 = f - be1 * b15; be1 = be1 + b16
            ENDIF

!  Form product of lo parts and add in lo part of intermediate product
!  to get lo part of complete product

            g = b15 * be2 + xlo * mlo

!  Represent lo part of full product in base 2**30

            d = g / b30 ; c = xhi / 2

!  Calculate full product divided by 2**30

            f = (( 2 * ( c * mhi - q ) - 1) + mhi * ( xhi - 2 * c )) + d + be1

!  Get full product divided in base 2**31

            IF ( f <= 0 ) THEN
               f = f + p ; nu = f / 2 ; mu = f - nu * 2
            ELSE
               f = f - 1 ; nu = f / 2 ; mu = f - 2 * nu ; nu = nu + b30
            ENDIF

!  Calculate remainder of product divided by 2**31

            f = ( b30 * mu - p ) + nu + ( g - b30 * d )
            random = nu + 1

!  Add one if remainder is not < 2**31-1

            IF ( f >= 0 ) random = random + 1
         ELSE

!  If n is less than or equal to 1, set random to n.

            random = n
         END IF

         END SUBROUTINE FA14_RANDOM_INTEGER

!-*-*-*-*-*-*-  F A 1 4 _ G E T _ S E E D  S U B R O U T I N E  *-*-*-*-*-*-*-

         SUBROUTINE FA14_GET_SEED ( seed, value )

!  Determine the current word generator.

         TYPE (FA14_seed), INTENT( IN ) :: seed
         INTEGER, INTENT( OUT ) :: value

         value = seed%ix

         END SUBROUTINE FA14_GET_SEED

!-*-*-*-*-*-*-  F A 1 4 _ S E T _ S E E D   S U B R O U T I N E  *-*-*-*-*-*-

         SUBROUTINE FA14_SET_SEED ( seed, value )

!  Reset the word generator to value if value lies in the
!  interval [1, 2**31 - 1]. More generally, value is set
!  to ( value - 1 ) mod (2**31 -1) + 1

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

         TYPE (FA14_seed), INTENT( OUT ) :: seed
         INTEGER, INTENT( IN ) :: value

         seed%ix = MOD( value - 1, p ) + 1

         END SUBROUTINE FA14_SET_SEED

!-*-*-*-*-*-*-  F A 1 4 _ G E T _ S E E D  S U B R O U T I N E  *-*-*-*-*-*-*-

         SUBROUTINE FA14_INITIALIZE ( seed )

!  Set the word generator to its default value.

         TYPE (FA14_seed), INTENT( OUT ) :: seed

         seed%ix = 65535

         END SUBROUTINE FA14_INITIALIZE

      END MODULE HSL_FA14_SINGLE

! COPYRIGHT (c) 2001 Council for the Central Laboratory
!                    of the Research Councils
!
! Version 2.2.0
! See ChangeLog for version history.
!
 module HSL_MC65_single

  use hsl_zd11_single
  implicit none
  private

  INTEGER, PARAMETER :: myreal = KIND( 1.0 )
  INTEGER, PARAMETER :: myint = KIND( 1)

  integer (kind = myint), public, parameter :: &
       MC65_ERR_MEMORY_ALLOC = -1, &  ! memory alloc failure
       MC65_ERR_MEMORY_DEALLOC = -2, & ! memory deallocate failure
       ! dimension mismatch in matrix_sum
       MC65_ERR_SUM_DIM_MISMATCH = -3,&
       ! dimension mismatch in matrix_multiply
       MC65_ERR_MATMUL_DIM_MISMATCH = -4, &
       ! dimension mismatch in matrix_multiply_graph
       MC65_ERR_MATMULG_DIM_MISMATCH = -5, &
       ! y = Ax with A pattern only and x real.
       MC65_ERR_MATVEC_NOVALUE = -6, &
       ! no vacant unit has been found MC65_MATRIX_WRITE
       MC65_ERR_NO_VACANT_UNIT = -7,&
       ! if in MC65_MATRIX_READ, the file <file_name> does not exist.
       MC65_ERR_READ_FILE_MISS = -8,&
       ! if in MC65_MATRIX_READ, opening of
       ! the file <file_name> returns iostat /= 0
       MC65_ERR_READ_OPEN = -9, &
       ! in MC65_matrix_read, string length of matrix type
       ! exceeds maximum allowed (2000)
       MC65_ERR_READ_MAXLEN = -10, &
       ! in MC65_matrix_read, the particular FORM
       ! is not supported
       MC65_ERR_READ_WRONGFORM = -11, &
       ! in MC65_MATRIX_CONDENSE, try to condense
       ! a matrix that is not of type "pattern"
       MC65_ERR_CONDENSE_NOPAT = -12,&
       ! if PTR(1:M+1) is not monotonically increasing
       ! MATRIX_CONSTRUCT
       MC65_ERR_RANGE_PTR = -13

  integer (kind = myint), public, parameter :: &
       !  COL is not within the range of [1,n] MATRIX_CONSTRUCT,
       ! such entries are excluded in the construct of the matrix
       MC65_WARN_RANGE_COL = 1, &
       ! IRN is not within the  range of {\tt [1,M]}. In
       ! matrix_construct using coordinate formatted input data.
       ! related entries excluded
       MC65_WARN_RANGE_IRN = 2, &
       ! JCN is not within the  range of {\tt [1,N]}. In
       ! matrix_construct using coordinate formatted input data.
       ! related entries excluded
       MC65_WARN_RANGE_JCN = 3,&
       ! IN MC65_MATRIX_DIAGONAL_FIRST, some diagonal elements
       ! are missing!
       MC65_WARN_MV_DIAG_MISSING = 4,&
       ! if duplicate entries were
       !  found when cleaning the matrix
       MC65_WARN_DUP_ENTRY = 5,&
       ! if both out-of-range and duplicate entries are found. In
       ! matrix construct
       MC65_WARN_ENTRIES = 6,&
       ! if both out-of-range IRN and JCN entries found, such entries
       ! are excluded in the construct of the matrix
       MC65_WARN_RANGE_BOTH = 7
  public &
       MC65_print_message, &
       MC65_matrix_construct, &
       MC65_matrix_destruct, &
       MC65_matrix_reallocate, &
       MC65_matrix_transpose, &
       MC65_matrix_copy, &
       MC65_matrix_clean, &
       MC65_matrix_sort, &
       MC65_matrix_sum, &
       MC65_matrix_symmetrize, &
       MC65_matrix_getrow, &
       MC65_matrix_getrowval, &
       MC65_matrix_is_symmetric, &
       MC65_matrix_is_pattern,&
       MC65_matrix_is_different, &
       MC65_matrix_multiply, &
       MC65_matrix_multiply_graph, &
       MC65_matrix_multiply_vector, &
       MC65_matrix_to_coo, &
       MC65_matrix_remove_diagonal,&
       MC65_matrix_diagonal_first, &
       MC65_matrix_write, &
       MC65_matrix_read,&
       MC65_matrix_condense

!       MC65_matrix_fill, &
!       MC65_matrix_component,&
!       MC65_matrix_crop_unsym, &
!       MC65_matrix_crop



  interface MC65_print_message
     module procedure csr_print_message
  end interface

  interface MC65_matrix_construct
     module procedure csr_matrix_construct
     module procedure csr_to_csr_matrix
     module procedure coo_to_csr_format
  end interface
  interface MC65_matrix_is_symmetric
     module procedure csr_matrix_is_symmetric
  end interface
  interface MC65_matrix_is_pattern
     module procedure csr_matrix_is_pattern
  end interface
  interface MC65_matrix_sort
     module procedure csr_matrix_sort
  end interface
  interface MC65_matrix_clean
     module procedure csr_matrix_clean
  end interface
  interface MC65_matrix_multiply
     module procedure csr_matrix_multiply
  end interface
  interface MC65_matrix_multiply_graph
     module procedure csr_matrix_multiply_graph
  end interface
  interface MC65_matrix_to_coo
     module procedure csr_matrix_to_coo
  end interface
  interface MC65_matrix_write
     module procedure csr_matrix_write
  end interface
  interface MC65_matrix_read
     module procedure csr_matrix_read
  end interface
  interface MC65_matrix_reallocate
     module procedure csr_matrix_reallocate
  end interface

  interface MC65_matrix_destruct
     module procedure csr_matrix_destruct
  end interface
  interface MC65_matrix_transpose
     module procedure csr_matrix_transpose
  end interface
  interface MC65_matrix_symmetrize
     module procedure csr_matrix_symmetrize
  end interface
  interface MC65_matrix_getrow
     module procedure csr_matrix_getrow
  end interface
  interface MC65_matrix_getrowval
     module procedure csr_matrix_getrowval
  end interface
  interface MC65_matrix_sum
     module procedure csr_matrix_sum
  end interface
  interface MC65_matrix_copy
     module procedure csr_matrix_copy
  end interface
  interface MC65_matrix_multiply_vector
     module procedure csr_matrix_multiply_rvector
     module procedure csr_matrix_multiply_ivector
  end interface
  interface MC65_matrix_condense
     module procedure csr_matrix_condense
  end interface
  interface MC65_matrix_diagonal_first
     module procedure csr_matrix_diagonal_first
  end interface
  interface MC65_matrix_remove_diagonal
     module procedure csr_matrix_remove_diagonal
  end interface
!!$  interface MC65_matrix_fill
!!$     module procedure csr_matrix_fill
!!$  end interface
!!$  interface MC65_matrix_component
!!$     module procedure csr_matrix_component
!!$  end interface
!!$  interface MC65_matrix_crop
!!$     module procedure csr_matrix_crop
!!$  end interface
!!$  interface MC65_matrix_crop_unsym
!!$     module procedure csr_matrix_crop_unsym
!!$  end interface
  interface MC65_matrix_is_different
     module procedure csr_matrix_diff
  end interface
  interface expand
     module procedure expand1
     module procedure iexpand1
  end interface


contains

  subroutine csr_print_message(info,stream,context)
    ! printing error message
    ! info: is an integer scaler of INTENT (IN).
    ! It is the information flag
    !       whose corresponding error message is to be printed.
    ! stream: is an OPTIONAL integer scaler of INTENT (IN).
    !         It is the unit number
    !         the user wish to print the error message to.
    !         If this number
    !         is negative, printing is supressed. If not supplied,
    !         unit 6 is the default unit to print the error message.
    ! context: is an OPTIONAL assumed size CHARACTER array
    !          of INTENT (IN).
    !          It describes the context under which the error occurred.
    integer (kind = myint), intent (in) :: info
    integer (kind = myint), intent (in), optional :: stream
    character (len = *), optional, intent (in) :: context
    integer (kind = myint) :: unit,length



    if (present(stream)) then
       unit = stream
    else
       unit = 6
    end if
    if (unit <= 0) return

    if (info > 0) then
       write(unit,advance = "yes", fmt = "(' WARNING: ')")
    else if (info < 0) then
       write(unit,advance = "yes", fmt = "(' ERROR: ')")
    end if

    if (present(context)) then
       length = len_trim(context)
       write(unit,advance = "no", fmt = "(a,' : ')") context(1:length)
    end if

    select case (info)
    case (0)
       write(unit,"(a)") "successful completion"
    case (MC65_ERR_MEMORY_ALLOC)
       write(unit,"(A)") "memory allocation failure failure"
    case(MC65_ERR_MEMORY_DEALLOC)
       write(unit,"(A)") "memory deallocate failure"
    case(MC65_ERR_SUM_DIM_MISMATCH)
       write(unit,"(A)") &
            "dimension mismatch of matrices in MC65_matrix_sum"
    case(MC65_ERR_MATMUL_DIM_MISMATCH)
       write(unit,"(A)") "dimension mismatch in matrix_multiply"
    case(MC65_ERR_MATMULG_DIM_MISMATCH)
       write(unit,"(A)") "dimension mismatch in matrix_multiply_graph"
    case(MC65_ERR_MATVEC_NOVALUE)
       write(unit,"(A)") "try to compute y = Ax or y = A^Tx "
       write(unit,"(A)") "with A pattern only and x real"
    case(MC65_ERR_RANGE_PTR)
       write(unit,"(A)") "PTR(1:M+1) is not monotonically &
            &increasing  in MATRIX_CONSTRUCT"
    case(MC65_ERR_NO_VACANT_UNIT)
       write(unit,"(A)") "no vacant I/O unit has been &
            &found MC65_MATRIX_WRITE"
    case(MC65_ERR_READ_FILE_MISS)
       write(unit,"(A)") "in MC65_MATRIX_READ, &
            &the file to be read does not exist"
    case(MC65_ERR_READ_OPEN)
       write(unit,"(A)") "error opening file in MC65_MATRIX_READ"
    case(MC65_ERR_READ_MAXLEN)
       write(unit,"(A)") "in MC65_matrix_read, &
            &string length of matrix type"
       write(unit,"(A)") "exceeds maximum allowed (2000)"
    case(MC65_ERR_READ_WRONGFORM)
       write(unit,"(A)") "in MC65_matrix_read, &
            &the supplied input format is not supported"
    case(MC65_ERR_CONDENSE_NOPAT)
       write(unit,"(A)") "in MC65_MATRIX_CONDENSE, try to condense"
       write(unit,"(A)") "a matrix that is not of type 'pattern'"
       ! warnings ===========
    case(MC65_WARN_MV_DIAG_MISSING)
       write(unit,"(A)") "IN MC65_MATRIX_DIAGONAL_FIRST, &
            &some diagonal elements"
       write(unit,"(A)") "are missing"
    case(MC65_WARN_RANGE_COL)
       write(unit,"(A)") "Some column indices in COL array are not "
       write(unit,"(A)") "within the range of [1,N] &
            &in MATRIX_CONSTRUCT"
       write(unit,"(A)") "such entries are excluded &
            &in the constructed matrix"
    case(MC65_WARN_RANGE_IRN)
       write(unit,"(A)") "IRN is not within the &
            &range of [1,M] in MATRIX_CONSTRUCT"
       write(unit,"(A)") "using coordinate formatted input data. &
            &Such entries are excluded"
    case(MC65_WARN_RANGE_JCN)
       write(unit,"(A)") "JCN is not within the &
            &range of [1,N] in matrix_construct"
       write(unit,"(A)") " using coordinate formatted &
            &input data. Such entries are excluded"
    case (MC65_WARN_DUP_ENTRY)
       write(unit,"(A)") "duplicate entries were found&
            &and had been merged (summed)"
    case (MC65_WARN_ENTRIES)
       write(unit,"(A)") "both duplicate entries and out-of-range entries &
            &were found &
            &and had been treated"
    case (MC65_WARN_RANGE_BOTH)
       write(unit,"(A)") "both out-of-range IRN and JCN entries found. &
            &Such entries are excluded in the constructed matrix"
    case default
       write(*,"(a,i10,a)") "you have supplied an info flag of ",&
            info," this is not a recognized info flag"
    end select
  end subroutine csr_print_message


  subroutine csr_matrix_construct(matrix,m,nz,info,n,type,stat)
    ! subroutine csr_matrix_construct(matrix,m,nz,info[,n,type])
    ! 1) initialise the matrix row/column dimensions and
    ! matrix type, allocate row pointers
    ! by default n = m and type = "general"
    ! 2) allocate a space of nz for the column indices matrix%col
    ! and when matrix%type /= "pattern",
    ! also allocate space for entry values matrix%val.

    ! matrix: of the derived type ZD11_type, INTENT (inOUT),
    !         the sparse matrix in compressed sparse row format
    type (zd11_type), intent (inout) :: matrix

    ! m: is an integer scaler of INTENT (IN), hold the row dimension
    integer (kind = myint), intent (in) :: m

    ! nz: integer scaler of INTENT (IN).
    !     It contains the storage needs to be allocated for the
    !     entries of matrix.
    !     including the storage needs to be allocated for
    !     the array that holds
    !     the column indices; and when type /= "pattern",
    !     also the space allocated to hold the values of the matrix.
    !     $nz$ must be greater than or equal to the number of
    !     entries in the matrix
    !     if $nz < 0$, a storage  of 0 entry is allocated.
    integer (kind = myint), intent (in) :: nz

    ! info: is an integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory allocation failed
    integer (kind = myint), intent (out) :: info

    ! n: is an optional integer scaler of INTENT (IN). If supplied
    !    it contains the column dimension of the matrix
    integer (kind = myint), optional, intent (in) :: n

    ! type: is an optional character array of unspecified length.
    !       INTENT (IN).
    !       If supplied gives the type of the matrix.
    !       If not supplied,
    !       the type of the matrix will be set to "general"
    character (len=*), optional, intent (in) :: type

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! ==================== locals =========================
    ! ierr: error tag for deallocation
    integer (kind = myint) :: ierr
    ! nz1: space allocated
    integer (kind = myint) :: nz1

    if (present(stat)) stat = 0
    info = 0
    ierr = 0

    matrix%m = m
    if (present(n)) then
       matrix%n = n
    else
       matrix%n = m
    end if

    ! CHECK THAT storage is not already existing and if so
    ! see if it is already large enough
    if (allocated(matrix%ptr)) then
       if (size(matrix%ptr) < m+1) then
          deallocate(matrix%ptr,stat = ierr)
          if (present(stat)) stat = ierr
          if (ierr /= 0) then
             info = MC65_ERR_MEMORY_DEALLOC
             return
          end if
          allocate(matrix%ptr(m+1),stat = ierr)
          if (present(stat)) stat = ierr
       end if
    else
       allocate(matrix%ptr(m+1),stat = ierr)
       if (present(stat)) stat = ierr
    end if
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

    ! get rid of the storage of type anyway without
    ! checking to see if
    ! it has large enough storage.
    if (allocated(matrix%type)) then
       deallocate(matrix%type, stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          info = MC65_ERR_MEMORY_DEALLOC
          return
       end if
    end if

    ! by default the matrix has entry values.
    if (present(type)) then
       call zd11_put(matrix%type,type,stat = ierr)
    else
       call zd11_put(matrix%type,"general",stat = ierr)
    end if
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

    ! allocate space for the column indices and possibly the values
    nz1 = nz
    if (nz1 < 0) then
       nz1 = 0
    end if

    ! check to see if col is already allocated and
    ! if so if the storage
    ! is large enough
    if (allocated(matrix%col)) then
       if (size(matrix%col) < nz1) then
          deallocate(matrix%col,stat = ierr)
          if (present(stat)) stat = ierr
          if (ierr /= 0) then
             info = MC65_ERR_MEMORY_DEALLOC
             return
          end if
          allocate(matrix%col(nz1),stat = ierr)
          if (present(stat)) stat = ierr
       end if
    else
       allocate(matrix%col(nz1),stat = ierr)
       if (present(stat)) stat = ierr
    end if
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

    if (ZD11_get(matrix%type) == "pattern")  then
       if (allocated(matrix%val)) then
          deallocate(matrix%val,stat = ierr)
          if (present(stat)) stat = ierr
          if (ierr /= 0) then
             info = MC65_ERR_MEMORY_DEALLOC
          end if
       end if
       return
    end if

    ! check if matrix%val is allocated already, if so
    ! cehck if the size is large enough
    if (allocated(matrix%val)) then
       if (size(matrix%val) < nz1) then
          deallocate(matrix%val,stat = ierr)
          if (present(stat)) stat = ierr
          if (ierr /= 0) then
             info = MC65_ERR_MEMORY_DEALLOC
             return
          end if
          allocate(matrix%val(nz1),stat = ierr)
          if (present(stat)) stat = ierr
       end if
    else
       allocate(matrix%val(nz1),stat = ierr)
       if (present(stat)) stat = ierr
    end if
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

  end subroutine csr_matrix_construct

  subroutine csr_matrix_destruct(matrix,info,stat)
    ! subroutine csr_matrix_destruct(matrix,info):
    !
    ! destruct the matrix object by deallocating all
    ! space occupied by
    ! matrix. including matrix%ptr, matrix%col and matrix%type.
    ! when matrix%type /= "pattern", also
    ! deallocate matrix%val

    ! matrix: is of the derived type ZD11_type,
    !         with INTENT (INOUT). It
    !         the sparse matrix object to be destroyed.
    type (zd11_type), intent (inout) :: matrix

    ! info: is an integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
    integer (kind = myint), intent (out) :: info

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! ===================== local variables =============
    ! ierr: error tag for deallocation
    integer (kind = myint) :: ierr

    info = 0
    if (present(stat)) stat = 0

    deallocate(matrix%col,matrix%ptr, stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if


    if (ZD11_get(matrix%type) == "pattern")  then
       deallocate(matrix%type, stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          info = MC65_ERR_MEMORY_DEALLOC
       end if
       return
    end if

    deallocate(matrix%type, matrix%val, stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if

  end subroutine csr_matrix_destruct

  subroutine csr_matrix_reallocate(matrix,nz,info,stat)
    ! subroutine csr_matrix_reallocate(matrix,nz,info)
    !
    ! reallocate a space of nz for the column index array
    ! and when matrix%type /= "pattern", for the value array
    ! of the sparse matrix in the compact sparse row format.
    ! if nz < 0, a space of 0 will be allocated.
    ! The original context of the array(s),
    ! will be copied to the beginning of the reallocated
    ! array(s). When nz is smaller than the size of
    ! the said array(s), only part of the original
    ! array(s) up to the nz-th element is copied.

    ! matrix: is of the derived type ZD11_type with INTENT (INOUT),
    !         this is the sparse matrix to be reallocated.
    type (zd11_type), intent (inout) :: matrix

    ! nz:  is an integer scaler of INTENT (IN). It holds the
    !     space to be reallocated for the array of
    !     the column indices; and when the matrix
    !     is of type "pattern",
    !     also holds the space to be reallocated
    !     for the array of the entry
    !     values of the matrix.
    !     If nz < 0, a space of 0 will be allocated.

    integer (kind = myint), intent (in) :: nz

    ! info: is an integer scaler of INTENT (OUT).
    !       = 0  if the subroutine returns successfully;
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed.
    integer (kind = myint), intent (out) :: info

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! ================= local variables ================
    ! nz1: actual space allocated
    integer (kind = myint) :: nz1, ierr

    info = 0
    if (present(stat)) stat = 0

    if (nz < 0) then
       nz1 = 0
    else
       nz1 = nz
    end if

    if (size(matrix%col) /= nz1) then
       call expand(matrix%col,nz1,info,ierr)
       if (present(stat)) stat = ierr
       if (info < 0) then
          return
       end if
    end if

    if (ZD11_get(matrix%type) == "pattern")  return

    if (size(matrix%val) /= nz1) then
       call expand(matrix%val,nz1,info,ierr)
       if (present(stat)) stat = ierr
    end if


  end subroutine csr_matrix_reallocate

  subroutine csr_matrix_transpose(MATRIX1,MATRIX2,info,merge,pattern,stat)
    ! subroutine csr_matrix_transpose(MATRIX1,MATRIX2[,merge])
    !
    ! find the transpose of matrix MATRIX1 and
    ! put it in matrix MATRIX2. MATRIX2 = MATRIX1^T
    ! If merge is present and merge = .true.,
    ! repeated entries in MATRIX1 will be summed
    ! in MATRIX2
    ! If pattern is present and pattern = .true.,
    ! the new matrix will be of type "pattern" and
    ! will have no values

    ! MATRIX1: of the derived type ZD11_type, INTENT (INOUT),
    !    the matrix to be transposed
    type (zd11_type), intent (in) :: MATRIX1
    ! MATRIX2: of the derived type ZD11_type, INTENT (INOUT),
    !    MATRIX2 = MATRIX1^T. If merge = .true.,
    !    repeated entries of MATRIX1 are
    !    summed in MATRIX2.
    type (zd11_type), intent (inout) :: MATRIX2

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
    integer (kind = myint), intent (out) :: info

    ! merge: optional logical of INTENT (IN).
    !        if merge = .true, repeated entries of
    !        MATRIX1 will be merged in MATRIX2
    !        by summing the entry values. By default,
    !        no merge is performed.
    logical, intent (in), optional :: merge

    ! pattern: optional logical of INTENT (IN).
    !          if pattern = .true., MATRIX2 will have only
    !          the pattern of MATRIX1^T
    logical, intent (in), optional :: pattern

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! ================== local variables ===============
    !
    ! rowsz: number of entries in each row of MATRIX2
    integer (kind = myint), allocatable, dimension (:) :: rowsz


    ! mask: a masking array to see if a particular entry in a row
    !    has already appeared before.
    integer (kind = myint), dimension (:), allocatable :: mask

    integer (kind = myint) :: n,m,ierr,nz

    ! to_merge: whether repeated entries should be
    ! merged by summing the entry values
    logical :: to_merge

    ! pattern_only: whether the transpose matrix will be
    !    pattern only.
    logical :: pattern_only

    info = 0
    if (present(stat)) stat = 0

    to_merge = .false.
    pattern_only = .false.
    if (present(pattern)) then
       pattern_only = pattern
    end if
    if (ZD11_get(MATRIX1%type) == "pattern") pattern_only = .true.

    if (present(merge)) then
       to_merge = merge
    end if

    m = MATRIX1%n
    n = MATRIX1%m

    ! number of entries in each row of MATRIX2
    allocate(rowsz(m),stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if
    rowsz = 0

    if (to_merge) then
       allocate(mask(m),stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          info = MC65_ERR_MEMORY_ALLOC
          return
       end if
       mask = 0
       CALL csr_matrix_transpose_rowsz(n,MATRIX1%ptr,MATRIX1%col,rowsz,mask)
    else
       CALL csr_matrix_transpose_rowsz(n,MATRIX1%ptr,MATRIX1%col,rowsz)
    end if

    nz = sum(rowsz)

    if (pattern_only) then
       call csr_matrix_construct(MATRIX2,m,nz,info,n = n,&
            type = "pattern",stat=ierr)
    else
       call csr_matrix_construct(MATRIX2,m,nz,info,n = n,&
            type = ZD11_get(MATRIX1%type),stat=ierr)
    end if
    if (present(stat)) stat = ierr

    if (info < 0 ) return

    if (pattern_only) then
       if (to_merge) then
          call csr_matrix_transpose_pattern(m,n,MATRIX1%ptr,MATRIX1%col, &
                                  MATRIX2%ptr,MATRIX2%col,rowsz,mask)
       else
          call csr_matrix_transpose_pattern(m,n,MATRIX1%ptr,MATRIX1%col, &
                                  MATRIX2%ptr,MATRIX2%col,rowsz)
       end if
    else
       if (to_merge) then
          call csr_matrix_transpose_values(m,n,MATRIX1%ptr,MATRIX1%col, &
                   MATRIX1%val,MATRIX2%ptr,MATRIX2%col,MATRIX2%val,rowsz,mask)
       else
          call csr_matrix_transpose_values(m,n,MATRIX1%ptr,MATRIX1%col, &
                   MATRIX1%val,MATRIX2%ptr,MATRIX2%col,MATRIX2%val,rowsz)
       end if
    end if

    if (to_merge) then
       deallocate(mask,stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          info = MC65_ERR_MEMORY_DEALLOC
          return
       end if
    end if

    deallocate(rowsz,stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
    end if

  end subroutine csr_matrix_transpose


  subroutine csr_matrix_transpose_rowsz(n,ia,ja,rowsz,mask)
    ! get the number of elements in each row
    ! this is an internal private routine
    ! called by csr_matrix_transpose

    ! n: row dimension of A
    ! ia: row pointer of A
    ! ja: column indices of A
    ! rowsz: the number of elements in each row
    ! mask: zero on entry. On exit mask <= 0. If mask(k) = -i,
    !       column index
    !       k has been seen in row i of A.
    !       this argument is optional, and only present
    !       if repeated entries in A is to be merged when forming A^T
    integer (kind = myint), intent (in) :: n,ia(*),ja(*)
    integer (kind = myint), intent (inout) :: rowsz(*)
    integer (kind = myint), intent (inout), optional :: mask(*)
    integer (kind = myint) :: i,j,k

    if (present(mask)) then
       do i=1,n
          do j = ia(i),ia(i+1)-1
             k = ja(j)
             if (mask(k) /= -i) then
                mask(k) = -i
                rowsz(k) = rowsz(k) + 1
             end if
          end do
       end do
    else
       do i=1,n
          do j = ia(i),ia(i+1)-1
             k = ja(j)
             rowsz(k) = rowsz(k) + 1
          end do
       end do
    end if
  end subroutine csr_matrix_transpose_rowsz



  subroutine csr_matrix_transpose_values(m,n,ia,ja,aa,ib,jb,&
       bb,rowsz,mask)
    ! transpose the matrix A=(ia,ja,aa) that is not pattern only
    ! to give B=(ib,jb,bb). If mask is present,
    ! repeated entries in A will be merged
    ! this is an internal private routine
    ! called by csr_matrix_transpose

    ! n: row dimension of A
    ! m: column dimension of A
    ! ia: row pointer of A
    ! ja: column indices of A
    ! aa: values of A
    ! ib: row pointer of A
    ! jb: column indices of A
    ! bb: values of B
    ! rowsz: the number of elements in each row
    ! mask: zero on entry. On exit mask <= 0. If mask(k) = -i,
    !       column index
    !       k has been seen in row i of A.
    !       this argument is optional, and only present
    !       if repeated entries in A is to be merged when forming A^T
    integer (kind = myint), intent (in) :: m,n,ia(*),ja(*)
    integer (kind = myint), intent (out) :: ib(*),jb(*)
    real (kind = myreal), intent(in) :: aa(*)
    real (kind = myreal), intent(out) :: bb(*)
    integer (kind = myint), intent (inout) :: rowsz(*)
    integer (kind = myint), intent (inout), optional :: mask(*)

    integer (kind = myint) :: i,j,k,l1,l2

    ib(1) = 1
    do i = 1, m
       ib(i+1) = ib(i) + rowsz(i)
    end do
    rowsz(1:m) = ib(1:m)

    if (present(mask)) then
       do i=1,n
          l1 = ia(i); l2 = ia(i+1)-1
          do j = l1,l2
             k = ja(j)
             if (mask(k) <= 0) then
                ! record the position of (k,i) in the
                ! transpose matrix
                mask(k) = rowsz(k)
                jb(rowsz(k)) = i
                bb(rowsz(k)) = aa(j)
                rowsz(k) = rowsz(k) + 1
             else
                bb(mask(k)) = bb(mask(k))+aa(j)
             end if
          end do
          ! Note: Can't use following array statement instead of the next loop
          ! as there may be repeated entries in a column.
          ! mask(ja(l1:l2)) = 0
          do j = l1, l2
             mask(ja(j)) = 0
          end do
       end do

    else
       do i=1,n
          do j = ia(i),ia(i+1)-1
             k = ja(j)
             jb(rowsz(k)) = i
             bb(rowsz(k)) = aa(j)
             rowsz(k) = rowsz(k) + 1
          end do
       end do
    end if

  end subroutine csr_matrix_transpose_values




  subroutine csr_matrix_transpose_pattern(m,n,ia,ja,ib,jb,rowsz,mask)
    ! transpose the matrix A=(ia,ja) that is  pattern only
    ! to give B=(ib,jb). If mask is present,
    ! repeated entries in A will be merged
    ! this is an internal private routine
    ! called by csr_matrix_transpose

    ! n: row dimension of A
    ! m: column dimension of A
    ! ia: row pointer of A
    ! ja: column indices of A
    ! aa: values of A
    ! ib: row pointer of A
    ! jb: column indices of A
    ! bb: values of B
    ! rowsz: the number of elements in each row
    ! mask: zero on entry. On exit mask <= 0. If mask(k) = -i,
    !       column index
    !       k has been seen in row i of A.
    !       this argument is optional, and only present
    !       if repeated entries in A is to be merged when forming A^T
    integer (kind = myint), intent (in) :: m,n,ia(*),ja(*)
    integer (kind = myint), intent (out) :: ib(*),jb(*)
    integer (kind = myint), intent (inout) :: rowsz(*)
    integer (kind = myint), intent (inout), optional :: mask(*)

    integer (kind = myint) :: i,j,k,l1,l2
    ib(1) = 1
    do i = 1, m
       ib(i+1) = ib(i) + rowsz(i)
    end do
    rowsz(1:m) = ib(1:m)

    if (present(mask)) then
       do i=1,n
          l1 = ia(i); l2 = ia(i+1)-1
          do j = l1,l2
             k = ja(j)
             if (mask(k) <= 0) then
                ! record the position of (k,i) in the
                ! transpose matrix
                mask(k) = rowsz(k)
                jb(rowsz(k)) = i
                rowsz(k) = rowsz(k) + 1
             end if
          end do
          ! Note: Can't use following array statement instead of the next loop
          ! as there may be repeated entries in a column.
          ! mask(ja(l1:l2)) = 0
          do j = l1, l2
             mask(ja(j)) = 0
          end do
       end do

    else
       do i=1,n
          do j = ia(i),ia(i+1)-1
             k = ja(j)
             jb(rowsz(k)) = i
             rowsz(k) = rowsz(k) + 1
          end do
       end do
    end if

  end subroutine csr_matrix_transpose_pattern


  subroutine csr_matrix_copy(MATRIX1,MATRIX2,info,stat)
    ! subroutine csr_matrix_copy(MATRIX1,MATRIX2,info)
    !
    ! Subroutine {\tt MC65\_MATRIX\_COPY} creates a new sparse matrix
    ! {\tt B} which is a copy of the existing sparse matrix
    ! {\tt MATRIX1}.
    ! Storage is allocated for {\tt MATRIX2} and this storage should
    ! be deallocated when {\tt MATRIX2} is no longer used by calling
    ! subroutine {\tt MC65_MATRIX_DESTRUCT}.


    ! MATRIX1: of the derived type ZD11_type, INTENT (IN),
    !    the matrix to be copied from
    type (zd11_type), intent (in) :: MATRIX1

    ! MATRIX2: of the derived type ZD11_type, INTENT (IN),
    !    the matrix to be copied to
    type (zd11_type), intent (inout) :: MATRIX2

    ! info: integer scaler of INTENT (INOUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
    integer (kind = myint), intent (out) :: info

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! n: column size
    ! m: row size
    ! nz: number of entries in a.
    integer (kind = myint) :: nz,n,m

    if (present(stat)) stat = 0

    nz = MATRIX1%ptr(MATRIX1%m+1)-1
    n = MATRIX1%n
    m = MATRIX1%m
    if (present(stat)) then
       call csr_matrix_construct(MATRIX2,m,nz,info, n = n,&
            type = zd11_get(MATRIX1%type),stat=stat)
    else
       call csr_matrix_construct(MATRIX2,m,nz,info, n = n,&
            type = zd11_get(MATRIX1%type))
    end if
    if (info < 0) return
    MATRIX2%ptr(1:m+1) = MATRIX1%ptr(1:m+1)
    MATRIX2%col(1:nz) = MATRIX1%col(1:nz)
    if (ZD11_get(MATRIX1%type) == "pattern")  return
    MATRIX2%val(1:nz) = MATRIX1%val(1:nz)
  end subroutine csr_matrix_copy



  subroutine csr_matrix_clean(matrix,info,realloc,stat)
    ! subroutine csr_matrix_clean(matrix,info[,realloc])
    !
    ! cleaning the matrix by removing redundant entries
    ! for pattern only matrix, or
    ! by summing it with existing ones for non-pattern only matrix

    ! matrix:  is of the derived type {\tt ZD11\_type} with
    !     INTENT (INOUT),
    !     it is the  the sparse matrix to be cleaned
    type (zd11_type), intent (inout) :: matrix

    ! info: is a integer scaler of INTENT (OUT).
    ! {\tt INFO = 0} if the subroutine completes successfully;
    ! {\tt INFO = MC65\_ERR\_MEMORY\_ALLOC} if memory allocation
    !     failed; and
    ! {\tt INFO = MC65\_ERR\_MEMORY\_DEALLOC} if memory
    !     deallocation failed;
    ! INFO = MC65_WARN_DUP_ENTRY if duplicate entries were
    !  found when cleaning the matrix
    ! Error message can be printed by calling subroutine
    !    {\tt MC65\_ERROR\_MESSAGE}
    integer (kind = myint), intent (out) :: info

    ! realloc:  is an optional real scaler of INTENT (IN).
    !     It is used to control the reallocation of storage.
    ! \begin{itemize}
    !  \item{} if {\tt REALLOC < 0}, no reallocation.
    !  \item{} if {\tt REALLOC == 0}, reallocation is carried out
    !  \item{} if {\tt REALLOC > 0}, reallocation iff
    !       memory saving is greater than
    !      {\tt REALLOC*100}\%. For example,
    !      if {\tt REALLOC = 0.5}, reallocation will be carried out
    !      if saving in storage is greater than 50\%
    !  \end{itemize}
    ! If {\tt REALLOC} is not present, no reallocation
    !    is carried out.
    real (kind = myreal), INTENT(IN), optional :: realloc

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! m: row dimension of matrix
    ! n: column dimension of the matrix
    ! nz: number of entries after the matrix is cleaned
    integer (kind = myint) :: n,m,nz

    ! ia: row pointer of the original matrix
    ! ja: column entries
    ! ib:  row pointer of the cleaned matrix
    integer (kind = myint), allocatable, dimension(:) :: ib

    integer (kind = myint) :: ierr


    ! mask: mask array which tells if an entry has been seen before
    integer (kind = myint), allocatable, dimension (:) :: mask

    ! nz_old: space occupied before cleaning
    ! nentry_old: number of entries in the matrix before cleaning
    integer (kind = myint) :: nz_old, nentry_old

    info = 0
    if (present(stat)) stat = 0
    nz_old = size(matrix%col)
    m = matrix%m
    n = matrix%n
    nentry_old = matrix%ptr(m+1)-1

    allocate(mask(n), ib(m+1), stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if
    mask = 0

    if (ZD11_get(matrix%type) /= "pattern" ) then
       call csr_matrix_clean_private(m,nz,ib,mask,matrix%ptr,matrix%col, &
                                     a = matrix%val)
    else
       call csr_matrix_clean_private(m,nz,ib,mask,matrix%ptr,matrix%col)
    end if

    if (present(realloc)) then
       if (realloc == 0.0_myreal.or.(realloc > 0.and.&
            &nz_old > (1.0_myreal + realloc)*nz)) then
          call csr_matrix_reallocate(matrix,nz,info,stat=ierr)
          if (present(stat)) stat = ierr
          if (info < 0) return
       end if
    end if

    deallocate(mask,ib,stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
    end if

    if (nentry_old /= matrix%ptr(m+1)-1) then
       INFO = MC65_WARN_DUP_ENTRY
    END if

  end subroutine csr_matrix_clean


  subroutine csr_matrix_clean_private(m,nz,ib,mask,ia,ja,a)
    ! clean the matrix a = (ia,ja,a) by summing repeated entries
    ! or clean the matrix a = (ia,ja) by removing repeated entries
    ! this is a private routine called by csr_matrix_clean.
    !
    ! m: row dimension
    ! nz: number of nonzeros in the cleaned matrix
    ! ia: row pointer
    ! ja: column indices
    ! a: values. Optional
    ! ib: new row pointer for cleaned matrix
    ! mask: naming array set to zero on entry
    integer (kind = myint), intent (in) :: m
    integer (kind = myint), intent(inout) :: ia(*),ja(*),mask(*)
    integer (kind = myint), intent(out) :: nz,ib(*)
    real (kind = myreal), intent (inout), optional :: a(*)

    integer (kind = myint) :: i,j,jj

    ib(1) = 1
    nz = 1
    if (present(a)) then
       do i = 1, m
          do j = ia(i),ia(i+1)-1
             jj = ja(j)
             if (mask(jj) == 0) then
                mask(jj) = nz
                ja(nz) = jj
                a(nz) = a(j)
                nz = nz + 1
             else
                a(mask(jj)) = a(mask(jj)) + a(j)
             end if
          end do
          ib(i+1) = nz
          do j = ib(i),nz-1
             mask(ja(j)) = 0
          end do
       end do
    else
       do i = 1, m
          do j = ia(i),ia(i+1)-1
             jj = ja(j)
             if (mask(jj) == 0) then
                mask(jj) = nz
                ja(nz) = jj
                nz = nz + 1
             end if
          end do
          ib(i+1) = nz
          do j = ib(i),nz-1
             mask(ja(j)) = 0
          end do
       end do
    end if
    nz = nz - 1
    ia(1:m+1) = ib(1:m+1)
  end subroutine csr_matrix_clean_private


  subroutine csr_matrix_sort(matrix,info,stat)
    ! subroutine csr_matrix_sort(matrix,info)
    !
    ! sort all rows of a matrix into increasing column indices
    ! This is achieved by doing transpose twice.
    ! Repeated entries will be merged.

    ! matrix: is of the derived type ZD11_type, INTENT (inOUT),
    !         the sparse matrix in compressed sparse row format
    type (zd11_type), intent (inout) :: matrix

    ! info: integer scaler of INTENT (OUT).
    !      INFO = 0 if the subroutine completes successfully;
    !      INFO = MC65_ERR_MEMORY_ALLOC if memory
    !      allocation failed; and
    !      INFO = MC65_ERR_MEMORY_DEALLOC if memory
    !      deallocation failed.
    integer (kind = myint), intent (out) :: info

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! temporary copy of transpose of the matrix
    type (zd11_type) :: matrix_trans

    if (present(stat)) stat = 0

    if ( present(stat) ) then
       call csr_matrix_transpose(matrix,matrix_trans,info,merge=.true.,&
            stat=stat)
    else
       call csr_matrix_transpose(matrix,matrix_trans,info,merge=.true.)
    end if
    if (info < 0) return

    if ( present(stat) )then
       call csr_matrix_destruct(matrix,info,stat=stat)
    else
       call csr_matrix_destruct(matrix,info)
    end if
    if (info < 0) return

    if ( present(stat) ) then
       call csr_matrix_transpose(matrix_trans,matrix,info,stat=stat)
    else
       call csr_matrix_transpose(matrix_trans,matrix,info)
    end if
    if (info < 0) return

    if ( present(stat) ) then
       call csr_matrix_destruct(matrix_trans,info,stat=stat)
    else
       call csr_matrix_destruct(matrix_trans,info)
    end if

  end subroutine csr_matrix_sort

  subroutine csr_matrix_sum(matrix1,matrix2,result_matrix,&
       info,scaling,graph,stat)
    ! subroutine csr_matrix_sum(matrix1,matrix2,result_matrix,
    !       info[,scaling,graph])
    ! Subroutine {\tt MC65\_MATRIX\_SUM} adds two matrix
    !   {\tt MATRIX1} and {\tt MATRIX2} together,
    !    and return the resulting matrix in {\tt RESULT\_MATRIX}.
    !   If the {\tt OPTIONAL} argument {\tt GRAPH}
    !   is present and is {\tt .TRUE.},
    !   the {\tt RESULT\_MATRIX} is the sum of {\tt MATRIX1}
    !   and {\tt MATRIX2}, but with diagonal removed and
    !   all other entries set to 1.0
    !   (1.0D0 for {\tt HSL\_MC65\_DOUBLE}).
    !    Otherwise if the {\tt OPTIONAL} argument
    !   {\tt SCALING} is present,
    !    the scaled summation {\tt Matrix1+SCALING*MATRISX2}
    !    is calculated.

    ! matrix1: of the derived type ZD11_type, INTENT (IN),
    !         a sparse matrix in compressed sparse row format
    ! matrix2: of the derived type ZD11_type, INTENT (IN),
    !         a sparse matrix in compressed sparse row format
    type (zd11_type), intent (in) :: matrix1,matrix2

    ! \itt{result_matrix} is of the derived type
    !   {\tt ZD11\_type} with INTENT (inOUT),
    !   result_matrix is one of the following
    ! \begin{itemize}
    !     \item{} If {\tt GRAPH} is present and is
    !       {\tt .TRUE.}, {\tt RESULT\_MATRIX}
    !       is a matrix of type ``general'', and is set to the sum
    !       of {\tt MATRIX1} and {\tt MATRIX2},
    !       with the diagonal removed
    !       and with all entry values set to 1.0
    !     \item{} Otherwise {\tt RESULT\_MATRIX =
    !       MATRIX1 + MATRIX2}. The type of {\tt RESULT\_MATRIX}
    !       is set to ``pattern'' either {\tt MATRIX1}
    !       or {\tt MATRIX2} is of this type, or else it is
    !       set to the type of {\tt MATRIX1} if
    !       both {\tt MATRIX1} and {\tt MATRIX2}
    !       have the same type, otherwise
    !       it is set as ``general''.
    !     \item{} If the {\tt OPTIONAL} argument
    !       {\tt SCALING} is present,
    !       {\tt RESULT\_MATRIX = MATRIX1 + SCALING*MATRIX2}
    ! \end{itemize}
    type (zd11_type), intent (inout) :: result_matrix

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
    !       = MC65_ERR_SUM_DIM_MISMATCH if row (column)
    !         dimension of matrix1 and matrix2 does not match.
    integer (kind = myint), intent (out) :: info

    ! scaling:  matrix2 is to be scaled by *scaling is
    real (kind = myreal), intent (in), optional :: scaling

    ! graph: optional logical scaler of intent (in).
    !        When present and is .true.
    !        the diagonals of the summation is removed, and
    !        entry values set to one if the result_matrix is
    !        not pattern only.
    logical, intent(in), optional :: graph

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! m: row dimension of matrix1
    ! n: column dimension of matrix1
    ! nz: number of entries in the summation matrix1+matrix2
    integer (kind = myint) :: n,m,nz

    integer (kind = myint) :: ierr

    ! mask: mask array to indicate if a particular column index
    !       has been seen
    integer (kind = myint), allocatable, dimension (:) :: mask
    logical :: is_graph

    info = 0
    if (present(stat)) stat = 0

    is_graph = .false.
    if (present(graph)) then
       is_graph = graph
    end if

    m = matrix1%m
    n = matrix1%n
    if (m /= matrix2%m .or. n /= matrix2%n) then
       info = MC65_ERR_SUM_DIM_MISMATCH
       return
    end if

    allocate (mask(n), stat = ierr)
    if ( present(stat) ) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

    mask = 0

    ! work out the number of entries of the summation.
    if (is_graph) then
       call csr_matrix_sum_getnz(m,matrix1%ptr,matrix1%col,matrix2%ptr, &
                       matrix2%col,mask,nz,graph = is_graph)
    else
       call csr_matrix_sum_getnz(m,matrix1%ptr,matrix1%col,matrix2%ptr, &
                       matrix2%col,mask,nz)
    end if

    ! construct the new matrix
    ! if one of the matrix is pattern only, then
    ! the result matrix has to be pattern only. If
    ! a graph is to be formed, the matrix is always of type "general"
    ! and filled with values 1.0
    if (is_graph) then
       call csr_matrix_construct(result_matrix,m,nz,info,&
            n = n,type = "general",stat=ierr)
    else if (ZD11_get(matrix1%type) == "pattern".or.&
         &ZD11_get(matrix2%type) == "pattern") then
       call csr_matrix_construct(result_matrix,m,nz,info,&
            n = n,type = "pattern",stat=ierr)
    else if (ZD11_get(matrix1%type) == ZD11_get(matrix2%type)) then
       call csr_matrix_construct(result_matrix,m,nz,info,n = n,&
            type = ZD11_get(matrix1%type),stat=ierr)
    else
       call csr_matrix_construct(result_matrix,m,nz,info,&
            n = n,type = "general",stat=ierr)
    end if
    if (present(stat)) stat=ierr
    if (info < 0) then
       return
    end if

    if (ZD11_get(result_matrix%type) /= "pattern" ) then
       IF (PRESENT(SCALING).and.(.not.is_graph)) THEN
          call csr_matrix_sum_values(m,nz,matrix1%ptr,matrix1%col,matrix1%val, &
                   matrix2%ptr,matrix2%col,matrix2%val,result_matrix%ptr, &
                   result_matrix%col,result_matrix%val,mask,scaling = scaling)
       ELSE if (is_graph) then
          call csr_matrix_sum_graph(m,nz,matrix1%ptr,matrix1%col, &
                   matrix2%ptr,matrix2%col,result_matrix%ptr, &
                   result_matrix%col,result_matrix%val,mask)
       else
          call csr_matrix_sum_values(m,nz,matrix1%ptr,matrix1%col,matrix1%val, &
                   matrix2%ptr,matrix2%col,matrix2%val,result_matrix%ptr, &
                   result_matrix%col,result_matrix%val,mask)
       END IF
    else
       call csr_matrix_sum_pattern(m,nz,matrix1%ptr,matrix1%col, &
                   matrix2%ptr,matrix2%col,result_matrix%ptr, &
                   result_matrix%col,mask)
    end if

    deallocate(mask,stat = ierr)
    if ( present(stat) ) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
    end if

  end subroutine csr_matrix_sum


  subroutine csr_matrix_sum_getnz(m,ia,ja,ib,jb,mask,nz,graph)
    ! get the number of nonzeros in the sum of two matrix
    ! (ia,ja) and (ib,jb), and return that number in nz
    ! this is a private routine called by csr_matrix_sum
    !
    ! m: row dimension of A
    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    ! mask: masking array to show whether a particular entry
    !       in the summed matrix has already been see. Initialised
    !       to zero on entry, on exit <= 0.
    ! nz: number of nonzeros in the summed matrix
    ! GRAPH: optional. If present and is true, diagonals
    !        in the summation matrix is not counted
    integer (kind = myint), intent (in) :: m,ia(*),ib(*),ja(*),jb(*)
    integer (kind = myint), intent (inout) :: mask(*)
    integer (kind = myint), intent (out) :: nz
    logical, intent (in),optional :: graph

    logical :: is_graph

    integer (kind = myint) :: i,j,jj

    is_graph = .false.
    if (present(graph)) then
       is_graph = graph
    end if

    if (is_graph) then
       nz = 0
       do i = 1, m
          do j = ia(i), ia(i+1)-1
             jj = ja(j)
             if (jj == i) cycle
             if (mask(jj) /= -i) then
                mask(jj) = -i
                nz = nz + 1
             end if
          end do
          do j = ib(i), ib(i+1)-1
             jj = jb(j)
             if (jj == i) cycle
             if (mask(jj) /= -i) then
                mask(jj) = -i
                nz = nz + 1
             end if
          end do
       end do
    else

       nz = 0
       do i = 1, m
          do j = ia(i), ia(i+1)-1
             jj = ja(j)
             if (mask(jj) /= -i) then
                mask(jj) = -i
                nz = nz + 1
             end if
          end do
          do j = ib(i), ib(i+1)-1
             jj = jb(j)
             if (mask(jj) /= -i) then
                mask(jj) = -i
                nz = nz + 1
             end if
          end do
       end do
    end if
  end subroutine csr_matrix_sum_getnz


  subroutine csr_matrix_sum_values(m,nz,ia,ja,a,ib,jb,b,ic,jc,c,&
       mask,scaling)
    ! sum the two matrix A=(ia,ja,a) and B=(ib,jb,b) into C=(ic,jc,c)
    ! if scaling is present, B will be scaled by it before summing
    !
    !
    ! m: row dimension of A
    ! nz: number of nonzeros in the summation C
    ! ia: row pointer of A
    ! ja: column indices of A
    ! a: values of A
    ! ib: row pointer of B
    ! jb: column indices of B
    ! b: values of B
    ! ic: row pointer of C
    ! jc: column indices of C
    ! c: values of C
    ! mask: masking array to show whether a particular entry
    !       in the summed matrix has already been see. Initialised
    !       to <= 0 on entry
    integer (kind = myint), intent (in) :: m,ia(*),ib(*),ja(*),jb(*)
    real (kind = myreal), intent (in) :: a(*),b(*)
    integer (kind = myint), intent (out) :: ic(*),jc(*)
    real (kind = myreal), intent (out) :: c(*)
    integer (kind = myint), intent (inout) :: mask(*)
    real (kind = myreal), intent (in), optional :: scaling
    integer (kind = myint), intent (out) :: nz


    integer (kind = myint) :: i,j,jj


    ic(1) = 1
    nz = 1
    if (.not.present(scaling)) then
       do i = 1, m
          do j = ia(i), ia(i+1)-1
             jj = ja(j)
             if (mask(jj) <= 0) then
                mask(jj) = nz
                jc(nz) = jj
                c(nz) = a(j)
                nz = nz + 1
             else
                c(mask(jj)) = c(mask(jj)) +  a(j)
             end if
          end do
          do j = ib(i),ib(i+1)-1
             jj = jb(j)
             if (mask(jj) <= 0) then
                mask(jj) = nz
                jc(nz) = jj
                c(nz) = b(j)
                nz = nz + 1
             else
                c(mask(jj)) = c(mask(jj)) +  b(j)
             end if
          end do
          ic(i+1) = nz
          do j = ic(i),nz-1
             mask(jc(j)) = 0
          end do
       end do
    else
       do i = 1, m
          do j = ia(i), ia(i+1)-1
             jj = ja(j)
             if (mask(jj) <= 0) then
                mask(jj) = nz
                jc(nz) = jj
                c(nz) = a(j)
                nz = nz + 1
             else
                c(mask(jj)) = c(mask(jj)) +  a(j)
             end if
          end do
          do j = ib(i),ib(i+1)-1
             jj = jb(j)
             if (mask(jj) <= 0) then
                mask(jj) = nz
                jc(nz) = jj
                c(nz) = scaling*b(j)
                nz = nz + 1
             else
                c(mask(jj)) = c(mask(jj)) +  scaling*b(j)
             end if
          end do
          ic(i+1) = nz
          do j = ic(i),nz-1
             mask(jc(j)) = 0
          end do
       end do

    end if
    nz = nz-1
  end subroutine csr_matrix_sum_values



  subroutine csr_matrix_sum_graph(m,nz,ia,ja,ib,jb,ic,jc,c,mask)
    ! sum the two matrix A=(ia,ja,a) and B=(ib,jb,b) into C=(ic,jc,c)
    ! diagonals of C is ignored and
    ! all entries of C set to one
    !
    ! m: row dimension of A
    ! nz: number of nonzeros in the summation C
    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    ! ic: row pointer of C
    ! jc: column indices of C
    ! c: values of C
    ! mask: masking array to show whether a particular entry
    !       in the summed matrix has already been see. Initiualised
    !       to <= 0 on entry
    integer (kind = myint), intent (in) :: m,ia(*),ib(*),ja(*),jb(*)
    integer (kind = myint), intent (out) :: ic(*),jc(*)
    real (kind = myreal), intent (out) :: c(*)
    integer (kind = myint), intent (inout) :: mask(*)
    integer (kind = myint), intent (out) :: nz


    integer (kind = myint) :: i,j,jj


    ic(1) = 1
    nz = 1
    do i = 1, m
       do j = ia(i), ia(i+1)-1
          jj = ja(j)
          if (jj == i) cycle
          if (mask(jj) <= 0) then
             mask(jj) = nz
             jc(nz) = jj
             c(nz) =  1.0_myreal
             nz = nz + 1
          else
             c(mask(jj)) = 1.0_myreal
          end if
       end do
       do j = ib(i),ib(i+1)-1
          jj = jb(j)
          if (jj == i) cycle
          if (mask(jj) <= 0) then
             mask(jj) = nz
             jc(nz) = jj
             c(nz) = 1.0_myreal
             nz = nz + 1
          else
             c(mask(jj)) = 1.0_myreal
          end if
       end do
       ic(i+1) = nz
       do j = ic(i),nz-1
          mask(jc(j)) = 0
       end do
    end do
    nz = nz-1
  end subroutine csr_matrix_sum_graph


  subroutine csr_matrix_sum_pattern(m,nz,ia,ja,ib,jb,ic,jc,mask)
    ! sum the two matrix A=(ia,ja) and B=(ib,jb, into C=(ic,jc),
    !    patterns only.
    ! m: row dimension of A
    ! nz: number of nonzeros in the summation C
    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    ! ic: row pointer of C
    ! jc: column indices of C
    ! mask: masking array to show whether a particular entry
    !       in the summed matrix has already been see. Initiualised
    !       to <= 0 on entry
    ! graph: optional logical, when present and is .true.,
    !        the diagonals of the C will not be generated.
    integer (kind = myint), intent (in) :: m,ia(*),ib(*),ja(*),jb(*)
    integer (kind = myint), intent (out) :: ic(*),jc(*)
    integer (kind = myint), intent (inout) :: mask(*)
    integer (kind = myint), intent (out) :: nz

    integer (kind = myint) :: i,j,jj

!!$    is_graph = .false.
!!$    if (present(graph)) then
!!$       is_graph = graph
!!$    end if


    ic(1) = 1
    nz = 1

!!$    if (is_graph) then
!!$       do i = 1, m
!!$          do j = ia(i), ia(i+1)-1
!!$             jj = ja(j)
!!$             if (jj == i) cycle
!!$             if (mask(jj) <= 0) then
!!$                mask(jj) = nz
!!$                jc(nz) = jj
!!$                nz = nz + 1
!!$             end if
!!$          end do
!!$          do j = ib(i),ib(i+1)-1
!!$             jj = jb(j)
!!$             if (jj == i) cycle
!!$             if (mask(jj) <= 0) then
!!$                mask(jj) = nz
!!$                jc(nz) = jj
!!$                nz = nz + 1
!!$             end if
!!$          end do
!!$          ic(i+1) = nz
!!$          do j = ic(i),nz-1
!!$             mask(jc(j)) = 0
!!$          end do
!!$       end do
!!$    else
       do i = 1, m
          do j = ia(i), ia(i+1)-1
             jj = ja(j)
             if (mask(jj) <= 0) then
                mask(jj) = nz
                jc(nz) = jj
                nz = nz + 1
             end if
          end do
          do j = ib(i),ib(i+1)-1
             jj = jb(j)
             if (mask(jj) <= 0) then
                mask(jj) = nz
                jc(nz) = jj
                nz = nz + 1
             end if
          end do
          ic(i+1) = nz
          do j = ic(i),nz-1
             mask(jc(j)) = 0
          end do
       end do
!    end if
    nz = nz-1
  end subroutine csr_matrix_sum_pattern

  subroutine csr_matrix_symmetrize(matrix,info,graph,stat)
    ! The subroutine {\tt MC65\_MATRIX\_SYMMETRIZE} make the
    !    matrix symmetric by summing it with its transpose.
    !    If the optional argument {\tt GRAPH}
    !    is present and is set to {\tt .TRUE.}, the diagonal of the
    !    summation is not included, and whenever entry values are
    !    available, they are set to {\tt 1.0}(or {\tt 1.0D0}for
    !    {\tt HSL\_MC65\_double}.

    ! matrix: of the derived type ZD11_type, INTENT (INOUT),
    !         the sparse matrix in compressed sparse row format
    !         to be made symmetric by matrix+matrix^T
    type (zd11_type), intent (inout) :: matrix

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
    !       = MC65_ERR_SUM_DIM_MISMATCH if trying to symmetrize
    !          a nonsquare matrix.
    integer (kind = myint), intent (out) :: info

    ! graph: optional logical scaler. when present and is .true.,
    !        the diagonals of matrix+matrix^T will be removed,
    !        and entry values set to 1.0_myreal matrix is
    !        not pattern only.
    logical, intent (in), optional :: graph

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! matrix1: transposed of matrix
    ! matrix2: A+A^T
    type (zd11_type) :: matrix1, matrix2

    ! is_graph: whether we are generating a graph by the summation
    logical :: is_graph

    if ( present(stat) ) stat = 0

    is_graph = .false.
    if (present(graph)) is_graph = graph
    if (present(stat)) then
       call csr_matrix_transpose(matrix,matrix1,info,stat=stat)
    else
       call csr_matrix_transpose(matrix,matrix1,info)
    end if
    if (info < 0) return

    if (is_graph) then
       if (present(stat)) then
          call csr_matrix_sum(matrix,matrix1,matrix2,info,graph=.true.,&
               &stat=stat)
       else
          call csr_matrix_sum(matrix,matrix1,matrix2,info,graph=.true.)
       end if
    else
       if (present(stat)) then
          call csr_matrix_sum(matrix,matrix1,matrix2,info,stat=stat)
       else
          call csr_matrix_sum(matrix,matrix1,matrix2,info)
       end if
    end if

    if (info < 0) return

    if (present(stat)) then
       call csr_matrix_copy(matrix2,matrix,info,stat=stat)
    else
       call csr_matrix_copy(matrix2,matrix,info)
    end if

    if (info < 0) return

    if (present(stat)) then
       call csr_matrix_destruct(matrix2,info,stat=stat)
    else
       call csr_matrix_destruct(matrix2,info)
    end if

    if (info < 0) return

    if (present(stat)) then
       call csr_matrix_destruct(matrix1,info,stat=stat)
    else
       call csr_matrix_destruct(matrix1,info)
    end if

    if (info < 0) return

  end subroutine csr_matrix_symmetrize


  subroutine csr_matrix_getrow(a,i,col)

    ! get the column indices of the i-th row of matrix a

    ! a: of the derived type ZD11_type, INTENT (IN), TARGET
    !    the sparse matrix in compressed sparse row format
    type (zd11_type), intent (in), target :: a

    ! i: is an integer of intent (in), it is the
    ! the index of a row of A to be extracted.
    integer (kind = myint), intent (in) :: i

    ! col: is the array of column indices for the row
    integer (kind = myint), pointer, dimension (:) :: col

    col => a%col(a%ptr(i):a%ptr(i+1)-1)

  end subroutine csr_matrix_getrow


  subroutine csr_matrix_getrowval(a,i,val)

    ! get the values of the i-th row of matrix a

    ! a: of the derived type ZD11_type, INTENT (IN), TARGET
    !    the sparse matrix in compressed sparse row format
    type (zd11_type), intent (in), target :: a

    ! i: the row index of A which the column values are
    !     to be extracted.
    integer (kind = myint), intent (in) :: i

    ! val: the real (or double) array of
    !     values of the i-th row of the matrix
    real (kind = myreal), pointer, dimension (:) :: val

    val => a%val(a%ptr(i):a%ptr(i+1)-1)

  end subroutine csr_matrix_getrowval


  subroutine csr_matrix_is_symmetric(matrix,symmetric,info,&
       pattern,tol,stat)
    ! function csr_matrix_is_symmetric(matrix,symmetric,&
    !    info[,pattern,tol,stat])
    !
    ! The subroutine {\tt MC65\_MATRIX\_IS\_SYMMETRIC}
    !   check whether a matrix is symmetric or not.
    !   {\tt symmetric = .TRUE.} if the matrix is symmetric,
    !   and {\tt symmetric= .FALSE.} if not. When the matrix
    !   is of type ``pattern'', only structural symmetry is
    !   checked, otherwise the matrix is said to be symmetric
    !   if it is both structural symmetric
    !   and that for each entry {\tt MATRIX(I,J)}, the inequality
    !  {\tt MATRIX(I,J) - MATRIX(J,I) <= E} holds.
    !  Here {\tt E} is a tolerance set by default to
    !  {\tt E = 0.0} (or {\tt E = 0.0D0) for HSL\_MC65\_DOUBLE}.

    ! Two optional arguments are supplied. If the optional argument
    !   {\tt PATTERN} is present and is set to {\tt .TRUE.},
    !   only structural symmetry is checked. If the optional
    !   argument {\tt TOL} is present, the internal tolerance
    !   {\tt E} is set to {\tt TOL}.



    ! matrix:  is of the derived type ZD11_type with INTENT (IN).
    !     It is the sparse matrix whose symmetry is to be checked.
    type (zd11_type), intent (in) :: matrix

    ! symmetric: is logical of intent (out).
    !             It is set to {\tt .TRUE.} only in one of
    !            the following cases
    !            1) When {\tt PATTERN} is not present, or
    !               is present but {\tt PATTERN /= .TRUE.}
    !               1a) the matrix is not of type "pattern" and
    !                   is symmetric both structurally
    !                   and numerically.
    !               1b) the matrix is of type "pattern" and
    !                   is structurally symmetric.
    !            2) when {\tt PATTERN} is present and
    !               {\tt PATTERN = .TRUE.}, and
    !               the matrix is structurally symmetric
    logical, intent (out) :: symmetric


    ! info: integer scaler of INTENT (OUT).
    !      {\tt INFO = 0} if the subroutine returns successfully;
    !      {\tt INFO = MC65_ERR_MEMORY_ALLOC} if memory
    !         allocation failed; and
    !      {\tt INFO = MC65_ERR_MEMORY_DEALLOC} if memory
    !        deallocation failed
    integer (kind = myint), intent (out) :: info

    ! pattern:  is an optional logical scaler of INTENT (IN).
    !         when it is present and is {\tt .TRUE.},
    !         only structural symmetry is checked.
    logical, optional, intent (in) :: pattern

    ! tol: an optional real scaler of INTENT (IN).
    !      tol is the tolerance used to decide if two
    !      entries are the same
    real (kind = myreal), optional, intent (in) :: tol

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! matrix_trans: the transpose of matrix
    type (zd11_type) :: matrix_trans

    ! m: row dimension of matrix
    ! n: column dimension of matrix
    ! nz: number of nonzero entries in the matrix
    integer (kind = myint) :: n,m


    integer (kind = myint), dimension (:), allocatable :: mask
    integer (kind = myint) :: ierr
    real (kind = myreal), dimension (:), allocatable :: val

    ! structural: = .true. if only structural symmetry is checked
    logical :: structural

    ! e: error allowed for judging symmetry on matrix that
    !    is not pattern only.
    real (kind = myreal) :: e

    info = 0
    if (present(stat)) stat = 0

    e = real(0,myreal)
    if (present(tol)) then
       e = tol
    end if

    structural = .false.
    if (present(pattern)) then
       structural = pattern
    end if
    if (zd11_get(matrix%type) == "pattern") structural = .true.

    symmetric = .true.
    m = matrix%m
    n = matrix%n
    if (n /= m) then
       symmetric = .false.
       return
    end if

    ! transpose matrix
    if (structural) then
       if (present(stat)) then
          call csr_matrix_transpose(matrix,matrix_trans,info,&
               merge = .true.,pattern = .true., stat = stat)
       else
          call csr_matrix_transpose(matrix,matrix_trans,info,&
               merge = .true.,pattern = .true.)
       end if
    else
       if (present(stat)) then
          call csr_matrix_transpose(matrix,matrix_trans,info,&
               merge = .true., stat=stat)
       else
          call csr_matrix_transpose(matrix,matrix_trans,info,&
               merge = .true.)
       end if
    end if
    if (info < 0) return

    allocate(mask(m),stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
    end if

    if (structural) then
       symmetric = csr_matrix_is_same_pattern(m,n,&
            ia = matrix%ptr,ja = matrix%col,&
            ib = matrix_trans%ptr,jb = matrix_trans%col,mask = mask)
    else
       allocate(val(m),stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          info = MC65_ERR_MEMORY_ALLOC
       end if
       symmetric = csr_matrix_is_same_values(m,n,&
            ia = matrix%ptr,ja = matrix%col,a = matrix%val, &
            ib = matrix_trans%ptr,jb = matrix_trans%col,&
            b = matrix_trans%val, &
            mask = mask,val = val,tol = e)
       ! check only structural symmetry as well as values
       deallocate(val,stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          info = MC65_ERR_MEMORY_DEALLOC
          return
       end if
    end if

    deallocate(mask,stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if
    call csr_matrix_destruct(matrix_trans,info)

  end subroutine csr_matrix_is_symmetric


  function csr_matrix_is_same_pattern(m,n,ia,ja,ib,jb,mask) &
       result(thesame)
    ! check if two matrix, A=(ia,ja) and B = (ib,jb),
    ! is structurally the same. Here B is assumed to have
    ! no repeated entries. A and B is assumed to have the same row
    ! size.
    ! This is a private function called by csr_matrix_is_symmetric
    !
    ! m: row size of A
    ! n: column size of A
    ! ia: row pointer array of A
    ! ja: column indices of A
    ! ib: row pointer array of B
    ! jb: column indices of B
    ! mask: a mask which indicate if a certain column entry k
    !       has been seen in this row i (if mask(k) = i). size >= n
    ! thesame: A and B is the same pattern
    integer (kind = myint), intent (in) :: m,n,ia(*),ja(*),ib(*),jb(*)
    integer (kind = myint), intent (inout) :: mask(*)
    logical :: thesame
    !
    ! na: number of distinctive entries in a row of A
    integer (kind = myint) :: na
    integer (kind = myint) :: i,j,k,l1,l2

    thesame = .true.

    mask(1:n) = 0
    do i = 1, m
       na = 0
       do j = ia(i),ia(i+1)-1
          k = ja(j)
          if (mask(k) /= i) then
             mask(k) = i
             na = na + 1
          end if
       end do
       l1 = ib(i); l2 = ib(i+1) - 1
       ! make sure the number of distinctive entries
       ! of the A and B is the same. Since B is assumed to have
       ! no repeated entries, its number of distinctive entries
       ! is simply l2-l1+1
       if (na /= l2 - l1 + 1) then
          thesame = .false.
          return
       end if
       do j = l1,l2
          if (mask(jb(j)) /= i) then
             thesame = .false.
             return
          end if
       end do
    end do
  end function csr_matrix_is_same_pattern





  function csr_matrix_is_same_values(m,n,ia,ja,a,ib,jb,b,mask,&
       val,tol) result(thesame)
    ! check if two matrix, A=(ia,ja,a) and B = (ib,jb,b),
    ! is the same structurally as well as in values.
    ! Here B is assumed to have
    ! no repeated entries. A and B is assumed to have the same row
    ! size as well as column size.
    ! This is a private function called by csr_matrix_is_symmetric
    !
    ! m: row size of A
    ! n: column size of A.
    ! ia: row pointer array of A
    ! ja: column indices of A
    ! a: entry values of A
    ! ib: row pointer array of B
    ! jb: column indices of B
    ! b: entry values of B
    ! mask: a mask which indicate if a certain column entry k
    !       has been seen in this row i (if mask(k) = i). Initialised
    !       to zero. size >= n
    ! val: the entry values in a row of A. If A has repeated entries,
    !      val will contain the sum of the repeated entry values.
    ! thesame: A and B is the same pattern
    ! tol: error allowed for judging symmetry on matrix that
    !     is not pettern only.
    integer (kind = myint), intent (in) :: m,n,ia(*),ja(*),ib(*),jb(*)
    integer (kind = myint), intent (inout) :: mask(*)
    real (kind = myreal), intent (in) :: a(*),b(*)
    real (kind = myreal), intent (inout) :: val(*)
    logical :: thesame
    real (kind = myreal), intent (in) :: tol
    !
    ! na: number of distinctive entries in a row of A
    integer (kind = myint) :: na
    integer (kind = myint) :: i,j,k,l1,l2

    thesame = .true.
    mask(1:n) = 0
    do i = 1, m
       na = 0
       do j = ia(i),ia(i+1)-1
          k = ja(j)
          if (mask(k) /= i) then
             mask(k) = i
             na = na + 1
             val(k) = a(j)
          else
             val(k) = val(k) + a(j)
          end if
       end do
       l1 = ib(i); l2 = ib(i+1) - 1
       ! make sure the number of distinctive entries
       ! of the A and B is the same. Since B
       ! is assumed to have no repeated entries,
       ! its number of distinctive entries
       ! is simply l2-l1+1
       if (na /= l2 - l1 + 1) then
          thesame = .false.
          return
       end if
       do j = l1,l2
          k = jb(j)
          if (mask(k) /= i) then
             thesame = .false.
             return
          else
             if (abs(val(k)-b(j)) > tol) then
                thesame = .false.
                return
             end if
          end if
       end do
    end do
  end function csr_matrix_is_same_values


  subroutine csr_matrix_diff(MATRIX1,MATRIX2,diff,info,pattern,stat)
    !   subroutine csr_matrix_diff(MATRIX1,MATRIX2,diff,&
    !      info[,pattern,stat])

    ! return the difference between two matrix MATRIX1 and MATRIX2.
    !    Here MATRIX1 and MATRIX2 are assumed to have not
    !    repeated entries.
    !
    ! 1) if the dimension of the matrix or the number of entries
    !    is different, huge(0.0_myreal) is returned.
    ! 2) else if the two matrices have different patterns,
    !    huge(0.0_myreal) is returned
    ! 3) else if both matrices are not pattern only, f1 is returned
    ! here: f1=sum(abs((a-b)%val))
    ! 4) else if either are pattern only, 0.0_myreal is returned
    !
    ! If the optional argument pattern is present and is .true.,
    !   both matrices are taken as pattern only.

    ! MATRIX1,MATRIX2: of the derived type ZD11_type, INTENT (IN),
    !      the sparse matrix in compressed sparse row format
    !      Two matrices to be compared.
    type (zd11_type), intent (in) :: MATRIX1,MATRIX2

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
    integer (kind = myint), intent(out) :: info

    ! pattern: optional logical scaler of INTENT (IN). If present and
    !          is .true., only the difference in pattern is returned.
    logical, optional, INTENT (IN) :: pattern


    ! diff:  is a real (double precision real in HSL_MC65_double)
    !         scaler returned on successful completion of the
    !         function. It is difference in the matrices
    !         {\tt MATRIX1} and {\tt MATRIX2} defined as follows
    !   \begin{itemize}
    !   \item{} If the row or column dimension of {\tt MATRIX1}
    !           and {\tt MATRIX2} does not match,
    !           {\tt HUGE\{0.0\}} ({\tt HUGE\{0.0D0\}} for
    !           {\tt HSL_MC65_double}) is returned;
    !   \item{} Else if {\tt MATRIX1} and {\tt MATRIX2} has different
    !            pattern, {\tt HUGE\{0.0\}} ({\tt HUGE\{0.0D0\}} for
    !           {\tt HSL_MC65_double}) is returned;
    !   \item{} Else if either {\tt MATRIX1} or {\tt MATRIX2}
    !            is of type {\tt ``pattern''},
    !            0.0 (0.0D0 for {\tt HSL_MC65_double}) is returned.
    !   \item{} Else if neither {\tt MATRIX1} nor {\tt MATRIX2}
    !            is of the type
    !            {\tt ``pattern''}, the sum of the absolute value
    !             of the difference in entry value,
    !            {\tt SUM(ABS((MATRIX1\%VAL - MATRIX2\%VAL))},
    !            is returned.
    !   \end{itemize}

    real (kind = myreal), intent (out) :: diff

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat


    ! ma,ma: row dimension of a and b
    ! na,nb: column dimension of a and b
    ! nza,nzb: number of entries of a and b and c
    integer (kind = myint) :: ma,mb,na,nb,nza,nzb


    ! pattern_only: only find the difference in patterns
    logical :: pattern_only,pattern_the_same

    integer (kind = myint), allocatable :: mask(:)

    integer (kind = myint) :: ierr


    info = 0
    if (present(stat)) stat = 0

    pattern_only = .false.
    if (present(pattern)) then
       pattern_only = pattern
    end if
    if (ZD11_get(MATRIX1%type) == "pattern".or.&
         ZD11_get(MATRIX2%type) == "pattern") &
         pattern_only = .true.


    ma = MATRIX1%m
    mb = MATRIX2%m
    na = MATRIX1%n
    nb = MATRIX2%n
    nza = MATRIX1%ptr(ma+1)-1
    nzb = MATRIX2%ptr(mb+1)-1
    if (ma /= mb .or. na /= nb .or. nza /= nzb) then
       diff = huge(0.0_myreal)
       return
    end if

    allocate(mask(na),stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if
    if (pattern_only) then
       pattern_the_same = csr_matrix_is_same_pattern(m = ma,n = na,&
            ia = MATRIX1%ptr,ja = MATRIX1%col,&
            ib = MATRIX2%ptr,jb = MATRIX2%col,mask = mask)
       if (pattern_the_same) then
          diff = 0.0_myreal
       else
          diff = huge(0.0_myreal)
       end if
    else
       diff = csr_matrix_diff_values(m = ma,n = na,&
            ia = MATRIX1%ptr,ja = MATRIX1%col,a=MATRIX1%val,&
            ib = MATRIX2%ptr,jb = MATRIX2%col,b=MATRIX2%val,&
            mask = mask)
    end if


    deallocate(mask, stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if

  end subroutine csr_matrix_diff


  function csr_matrix_diff_values(m,n,ia,ja,a,ib,jb,b,mask) &
       result(diff)
    ! find the difference in values of two matrices,
    !   A=(ia,ja,a) and B = (ib,jb,b),
    ! If the two are of different pattern,
    !   the difference is huge(0.0_myreal).
    ! Here A and B are assumed to have
    ! no repeated entries, they are assumed to have the same row
    ! size as well as column size.
    ! This is a private function called by csr_matrix_diff
    !
    ! m: row size of A
    ! n: column size of A.
    ! ia: row pointer array of A
    ! ja: column indices of A
    ! a: entry values of A
    ! ib: row pointer array of B
    ! jb: column indices of B
    ! b: entry values of B
    ! mask: a mask which indicate if a certain column entry k
    !       has been seen in this row i (if mask(k) = i). Initialised
    !       to zero. size >= n
    ! diff: difference in values between A and B.
    !       If A and B is of different pattern,
    !       diff - huge(0.0_myreal),
    !       else diff = sum(abs(A-B))
    real (kind = myreal) :: diff


    integer (kind = myint), intent (in) :: m,n,ia(*),ja(*),ib(*),jb(*)
    integer (kind = myint), intent (inout) :: mask(*)
    real (kind = myreal), intent (in) :: a(*),b(*)

    !
    integer (kind = myint) :: i,j,k,la1,la2,lb1,lb2,jj

    mask(1:n) = 0
    diff = 0
    do i = 1, m
       la1 = ia(i);la2 = ia(i+1)-1
       do j = la1,la2
          k = ja(j)
          if (mask(k) == 0) then
             mask(k) = j
          end if
       end do

       ! of the A and B is the same. Since neither
       ! is assumed to have repeated entries,
       ! its number of distinctive entries
       ! is simply ia(i+1)-ia(i)
       lb1 = ib(i);lb2 = ib(i+1)-1
       if (la2-la1 /= lb2-lb1) then
          diff = huge(0.0_myreal)
          return
       end if
       do j = lb1,lb2
          k = jb(j)
          jj = mask(k)
          if (jj == 0) then
             diff = huge(0.0_myreal)
             return
          else
             diff = diff + abs(b(j) - a(jj))
          end if
       end do
       do j = la1,la2
          mask(ja(j)) = 0
       end do
    end do
  end function csr_matrix_diff_values





  subroutine csr_matrix_multiply(matrix1,matrix2,result_matrix,info,stat)
    ! subroutine csr_matrix_multiply(matrix1,matrix2,&
    !     result_matrix,info)
    !
    ! multiply two matrices matrix1 and matrix2,
    ! and put the resulting matrix in result_matrix
    ! The storage is allocated inside this routine for
    ! the sparse matrix object {tt RESULT\_MATRIX},
    ! and when result_matrix is no longer needed,
    ! this memory should be deallocated using
    ! {\tt MATRIX\_DESTRUCT}


    ! matrix1: is of the derived type ZD11_type, INTENT (IN),
    !         the sparse matrix in compressed sparse row format
    ! matrix2: is of the derived type ZD11_type, INTENT (IN),
    !         the sparse matrix in compressed sparse row format
    type (zd11_type), intent (in) :: matrix1,matrix2

    ! result_matrix: of the derived type ZD11_type, INTENT (INOUT)
    !                result_matrix=matrix1*matrix2
    type (zd11_type), intent (inout) :: result_matrix

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed.
    !       = MC65_ERR_MATMUL_DIM_MISMATCH dimension mismatch
    integer (kind = myint), intent (out) :: info


    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! ================ local variables =================
    ! m: row size of the matrix1
    ! nz: number of nonzeros in the product matrix
    ! n: column side of matrix2
    integer (kind = myint) :: n,nz,m

    ! mask: an array with a size the same as the
    !       column size of matrix2,
    !       used to indicating whether an entry in the product
    !       has been seen or not.
    integer (kind = myint), allocatable, dimension(:) :: mask

    integer (kind = myint) :: ierr
    if (present(stat)) stat = 0

    info = 0

    m = matrix1%m
    n = matrix2%n

    ! check if the dimension match
    if (matrix1%n /= matrix2%m) then
       info = MC65_ERR_MATMUL_DIM_MISMATCH
       return
    end if

    allocate(mask(n), stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

    ! now find the size
    call matmul_size(m,n,matrix1%ptr,matrix1%col,matrix2%ptr,matrix2%col,&
         nz,mask)

    ! construct the new matrix
    ! if one of the matrix is pattern only, then
    ! the result matrix has to be pettern only
    if (ZD11_get(matrix1%type) == "pattern".or.&
         ZD11_get(matrix2%type) == "pattern") then
       call csr_matrix_construct(result_matrix,m,nz,info,n=n, &
            type = "pattern",stat=ierr)
       if (present(stat)) stat = ierr
       if (info < 0) return
    else if (ZD11_get(matrix1%type) == ZD11_get(matrix2%type)) then
       call csr_matrix_construct(result_matrix,m,nz,info,n=n, &
            type = ZD11_get(matrix1%type),stat=ierr)
       if (present(stat)) stat = ierr
       if (info < 0) return
    else
       call csr_matrix_construct(result_matrix,m,nz,info,n=n, &
            type = "general",stat=ierr)
       if (present(stat)) stat = ierr
       if (info < 0) return
    end if

    if (ZD11_get(result_matrix%type) == "pattern")  then
       call matmul_noval(m,n,matrix1%ptr,matrix1%col,matrix2%ptr,matrix2%col, &
                       result_matrix%ptr,result_matrix%col,mask)
    else
       call matmul_normal(m,n,matrix1%ptr,matrix1%col,matrix1%val, &
                       matrix2%ptr,matrix2%col,matrix2%val,result_matrix%ptr, &
                       result_matrix%col,result_matrix%val,mask)
    end if

    deallocate(mask, stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if

  end subroutine csr_matrix_multiply


  subroutine matmul_size(m,n,ia1,ja1,ia2,ja2,nz,mask)
    ! work out the number of nonzeros in the
    ! product of two matrix A=(ia1,ja1) and B=(ia2,ja2),

    ! m: row size of A
    ! n: column size of B
    integer (kind = myint), intent (in) :: n,m
    ! ia1: row pointer of A
    ! ja1: column indices of A
    ! ia2: row pointer of B
    ! ja2: column indices of B
    integer (kind = myint), intent (in), dimension (*) :: &
         ia1,ja1,ia2,ja2

    ! nz: number of nonzeros in the product matrix
    integer (kind = myint), intent (inout) :: nz

    ! mask: working array
    integer (kind = myint), intent (inout) :: mask(*)
    ! i,j,k,neigh,icol_add: working integers
    integer (kind = myint) :: i,j,k,neigh,icol_add

    ! initialise the mask array which is an array that has
    ! value i if column index already exist in row i of
    ! the new matrix
    mask(1:n) = 0
    nz = 0
    do i=1,m
       do j = ia1(i),ia1(i+1)-1
          neigh = ja1(j)
          do k = ia2(neigh),ia2(neigh+1)-1
             icol_add = ja2(k)
             if (mask(icol_add) /= i) then
                nz = nz + 1
                mask(icol_add) = i     ! add mask
             end if
          end do
       end do
    end do
  end subroutine matmul_size


  subroutine matmul_normal(m,n,ia,ja,a,ib,jb,b,ic,jc,c,mask)
    ! calculate the product of two matrix A=(ia,ja,a) and
    ! B = (ib,jb,b).

    ! m: row size of A
    ! n: column size of B
    integer (kind = myint), intent (in) :: n,m

    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    integer (kind = myint), intent (in) :: ia(*),ib(*),ja(*),jb(*)
    ! ic: row pointer of C
    ! jc: column indices of C
    integer (kind = myint), intent (out) :: ic(*),jc(*)
    ! a: entry values of A
    ! b: entry values of B
    ! c: entry values of C
    real (kind = myreal), intent (in) :: a(*),b(*)
    real (kind = myreal), intent (out) :: c(*)
    ! working array
    integer (kind = myint), intent (inout) :: mask(*)
    ! working variables
    integer (kind = myint) :: nz,i,j,k,icol,icol_add,neigh
    real (kind = myreal) :: aij

    ! initialise the mask array which is an array that
    ! has non-zero value if the
    ! column index already exist, in which case
    ! the value is the index of that column
    mask(1:n) = 0
    nz = 0
    ic(1) = 1
    do i=1,m
       do j = ia(i),ia(i+1)-1
          aij = a(j)
          neigh = ja(j)
          do k = ib(neigh),ib(neigh+1)-1
             icol_add = jb(k)
             icol = mask(icol_add)
             if (icol == 0) then
                nz = nz + 1
                jc(nz) = icol_add
                c(nz) =  aij*b(k)
                mask(icol_add) = nz     ! add mask
             else
                c(icol) = c(icol) + aij*b(k)
             end if
          end do
       end do
       do j = ic(i),nz
          ! done this row i, so set mask to zero again
          mask(jc(j)) = 0
       end do
       ic(i+1) = nz+1
    end do
  end subroutine matmul_normal




  subroutine matmul_noval(m,n,ia,ja,ib,jb,ic,jc,mask)
    ! calculate the product of two matrix A=(ia,ja) and
    ! B = (ib,jb). Sparse patterns only is calculated.

    ! m: row size of A
    ! n: column size of B
    integer (kind = myint), intent (in) :: n,m
    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    integer (kind = myint), intent (in) :: ia(*),ib(*),ja(*),jb(*)
    ! ic: row pointer of C
    ! jc: column indices of C
    integer (kind = myint), intent (out) :: ic(*),jc(*)

    ! working array
    integer (kind = myint), intent (inout) :: mask(*)
    ! working variables
    integer (kind = myint) :: nz,i,j,k,icol,icol_add,neigh

    ! initialise the mask array which is an array
    ! that has non-zero value if the
    ! column index already exist, in which case
    ! the value is the index of that column
    mask(1:n) = 0
    nz = 0
    ic(1) = 1
    do i=1,m
       do j = ia(i),ia(i+1)-1
          neigh = ja(j)
          do k = ib(neigh),ib(neigh+1)-1
             icol_add = jb(k)
             icol = mask(icol_add)
             if (icol == 0) then
                nz = nz + 1
                jc(nz) = icol_add
                mask(icol_add) = nz     ! add mask
             else
                cycle
             end if
          end do
       end do
       do j = ic(i),nz
          ! done this row i, so set mask to zero again
          mask(jc(j)) = 0
       end do
       ic(i+1) = nz+1
    end do
  end subroutine matmul_noval

  subroutine csr_matrix_multiply_graph(matrix1,matrix2,&
       result_matrix,info,col_wgt,stat)
    ! subroutine csr_matrix_multiply_graph(matrix1,matrix2,&
    !    result_matrix,info[,col_wgt])
    !
    ! Given two matrix, find
    ! result_matrix = (\bar matrix1)*(\bar matrix2) where
    ! (\bar A) is the matrix A with all non-zero entries set to 1.
    ! The diagonal entries of the result_matrix will be discarded.
    ! the routine is used for finding the row connectivity graph
    ! A, which is (\bar A)*(\bar A)^T, with diagonal elements removed

    ! matrix1: of the derived type ZD11_type, INTENT (IN),
    !         the sparse matrix in compressed sparse row format
    ! matrix2: of the derived type ZD11_type, INTENT (IN),
    !         the sparse matrix in compressed sparse row format
    type (zd11_type), intent (in) :: matrix1,matrix2

    ! result_matrix: of the derived type ZD11_type, INTENT (INOUT)
    !                result_matrix=matrix1*matrix2
    type (zd11_type), intent (inout) :: result_matrix

    ! col_wgt: optional real array of INTENT (IN).
    !          the super-column weight for the columns of
    !          matrix1. this is used when calculating A*A^T
    !          with A having super-columns The effect is the same as
    !          result_matrix =
    !          (\bar matrix1)*DIAG(col_wgt)*(\bar matrix2)
    real (kind = myreal), dimension (matrix1%n), intent (in), &
         optional :: col_wgt

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed.
    !       = MC65_ERR_MATMULG_DIM_MISMATCH dimension mismatch
    integer (kind = myint), intent (out), optional :: info

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! ================= local variables =================
    ! m: row size of matrix1
    ! n: column size of matrix2
    ! nz: number of nonzeros in the product matrix
    integer (kind = myint) :: n,nz,m

    ! mask: an array with a size the same as the column
    !       size of matrix2, used to indicating whether an
    !       entry in the product has been seen or not.
    integer (kind = myint), allocatable, dimension(:) :: mask

    integer (kind = myint) :: ierr

    info = 0
    if (present(stat)) stat = 0

    m = matrix1%m

    n = matrix2%n

    ! Check if dimensions match
    if ( matrix1%n /= matrix2%m) then
       info = MC65_ERR_MATMULG_DIM_MISMATCH
       return
    end if

    allocate(mask(n), stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
    end if

    ! now find the size
    call matmul_size_graph(m,n,matrix1%ptr,matrix1%col,matrix2%ptr,&
         matrix2%col,nz,mask)

    ! initialise the new matrix
    call csr_matrix_construct(result_matrix,m,nz,info,n = n,&
         type = "general",stat=ierr)
    if (present(stat)) stat = ierr
    if (info < 0) return

    if (present(col_wgt)) then
       call matmul_wgraph(m,n,matrix1%ptr,matrix1%col,matrix2%ptr, &
                  matrix2%col,result_matrix%ptr,result_matrix%col, &
                  result_matrix%val,mask,col_wgt)
    else
       call matmul_graph(m,n,matrix1%ptr,matrix1%col,matrix2%ptr, &
                  matrix2%col,result_matrix%ptr,result_matrix%col, &
                  result_matrix%val,mask)
    end if

    deallocate(mask, stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if

  end subroutine csr_matrix_multiply_graph


  subroutine matmul_size_graph(m,n,ia1,ja1,ia2,ja2,nz,mask)
    ! work out the number of nonzeros in the
    ! product of two matrix A=(ia1,ja1) and B=(ia2,ja2),
    ! diagonal element of A*B is not counted
    ! because this routine is used in finding row graphs.
    ! (a graph is represented as a symm. matrix with no diag. entries)
    ! m: row size of A
    ! n: column size of B
    integer (kind = myint), intent (in) :: n,m
    ! ia1: row pointer of A
    ! ja1: column indices of A
    ! ia2: row pointer of B
    ! ja2: column indices of B
    integer (kind = myint), intent (in), dimension (*) :: &
         ia1,ja1,ia2,ja2

    ! nz: number of nonzeros in the product matrix
    integer (kind = myint), intent (inout) :: nz

    ! mask: working array
    integer (kind = myint), intent (inout) :: mask(*)

    ! i,j,k,neigh,icol_add: working integers
    integer (kind = myint) :: i,j,k,neigh,icol_add

    ! initialise the mask array which is an array that has
    !  value i if column index already exist in row i
    !  of the new matrix
    mask(1:n) = 0
    nz = 0
    do i=1,m
       do j = ia1(i),ia1(i+1)-1
          neigh = ja1(j)
          do k = ia2(neigh),ia2(neigh+1)-1
             icol_add = ja2(k)
             if (icol_add /= i .and. mask(icol_add) /= i) then
                nz = nz + 1
                mask(icol_add) = i     ! add mask
             end if
          end do
       end do
    end do
  end subroutine matmul_size_graph

  subroutine matmul_graph(m,n,ia,ja,ib,jb,ic,jc,c,mask)
    ! given two matrix, A=(ia,ja) and B = (ib,jb), find
    ! C = (\bar A)*(\bar B) where
    ! (\bar A) is the matrix A with all non-zero entries set to 1.
    ! The diagonal entries of the c matrix will be discarded.
    ! the routine is mainly used for finding the row
    ! connectivity graph
    ! of a matrix A, which is (\bar A)*(\bar A)^T
    ! it is assumed that the size of arrays (ic,jc,c)
    ! used to hold C is large enough.

    ! m: row size of A
    ! n: column size of B
    integer (kind = myint), intent (in) :: n,m
    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    integer (kind = myint), intent (in), dimension (*) :: &
         ia,ib,ja,jb
    ! ic: row pointer of C
    ! jc: column indices of C
    integer (kind = myint), intent (out), dimension (*)  :: ic,jc
    ! c: entry values of C
    real (kind = myreal), intent (out), dimension (*)  :: c
    ! working array
    integer (kind = myint), intent (inout) :: mask(*)
    ! working variables
    integer (kind = myint) :: nz,i,j,k,icol,icol_add,neigh

    ! initialise the mask array which is an array
    ! that has non-zero value if the
    ! column index already exist, in which case
    ! the value is the index of that column
    mask(1:n) = 0
    nz = 0
    ic(1) = 1
    do i=1,m
       do j = ia(i),ia(i+1)-1
          neigh = ja(j)
          do k = ib(neigh),ib(neigh+1)-1
             icol_add = jb(k)
             if (icol_add /= i) then ! do not consider diagonal
                icol = mask(icol_add)
                if (icol == 0) then
                   nz = nz + 1
                   jc(nz) = icol_add
                   c(nz) =  1.0_myreal
                   mask(icol_add) = nz     ! add mask
                else
                   c(icol) = c(icol) + 1.0_myreal
                end if
             end if
          end do
       end do
       do j = ic(i),nz
          ! done this row i, so set mask to zero again
          mask(jc(j)) = 0
       end do
       ic(i+1) = nz+1
    end do
  end subroutine matmul_graph


  subroutine matmul_wgraph(m,n,ia,ja,ib,jb,ic,jc,c,mask,col_wgt)
    ! given two matrix, find
    ! C = (\bar A)*W*(\bar B) where
    ! (\bar A) is the matrix A, with all non-zero entries set to 1.
    ! W is a diagonal matrix consists of the column weights of A
    ! The diagonal entries of the c matrix will be discarded.
    ! the routine is mainly used for finding the
    ! row connectivity graph of
    ! a matrix A with column weights, which is (\bar A)*W*(\bar A)^T
    ! it is assumed that the size of arrays (ic,jc,c)
    ! used to hold C is large enough.

    ! m: row size of A
    ! n: column size of B
    integer (kind = myint), intent (in) :: n,m

    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    integer (kind = myint), intent (in), dimension (*) :: &
         ia,ib,ja,jb
    ! col_wgt: the column weight, or the diagonal matrix W
    real(kind = myreal), intent (in), dimension (*) :: col_wgt
    ! ic: row pointer of C
    ! jc: column indices of C
    integer (kind = myint), intent (out), dimension (*)  :: ic,jc
    ! c: entry values of C
    real (kind = myreal), intent (out), dimension (*)  :: c
    ! working array
    integer (kind = myint), intent (inout) :: mask(*)
    ! working variables
    integer (kind = myint) :: nz,i,j,k,icol,icol_add,neigh

    ! initialise the mask array which is an array
    ! that has non-zero value if the
    ! column index already exist, in which case
    ! the value is the index of that column
    mask(1:n) = 0
    nz = 0
    ic(1) = 1
    do i=1,m
       do j = ia(i),ia(i+1)-1
          neigh = ja(j)
          do k = ib(neigh),ib(neigh+1)-1
             icol_add = jb(k)
             if (icol_add /= i) then ! do not consider diagonal
                icol = mask(icol_add)
                if (icol == 0) then
                   nz = nz + 1
                   jc(nz) = icol_add
                   c(nz) =  col_wgt(neigh)
                   mask(icol_add) = nz     ! add mask
                else
                   c(icol) = c(icol) + col_wgt(neigh)
                end if
             end if
          end do
       end do
       do j = ic(i),nz
          ! done this row i, so set mask to zero again
          mask(jc(j)) = 0
       end do
       ic(i+1) = nz+1
    end do
  end subroutine matmul_wgraph



  subroutine csr_matrix_multiply_rvector(matrix,x,y,info,trans)
    ! subroutine csr_matrix_multiply_rvector(matrix,x,y,info)
    !
    ! y = matrix*x where x and y are real vectors. Dimension of y
    ! is checked and returned if it is smaller than the row dimension
    ! of x
    !
    ! if trans is present and == 1;
    ! y = matrix^T*x where x and y are real vectors.


    ! matrix: of the derived type ZD11_type, INTENT (IN),
    !         the sparse matrix in compressed sparse row format
    type (zd11_type), intent (in) :: matrix

    ! x: real array of intent (IN), a vector to be multiplied
    !    with the matrix
    real (kind = myreal), intent (in), dimension (*) :: x

    ! y: real array of intent (OUT), the result of matrix*x
    !    or matrix^T*x
    real (kind = myreal), intent (out), dimension (*) :: y

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MATVEC_NOVALUE, &  if A is of type
    !         pattern only and x real.
    integer (kind = myint), intent (out) :: info

    ! trans: optional logical scalar. If present and trans = true,
    ! the transpose of the matrix is multiplied with the vector
    logical, optional, intent (in) :: trans
    ! ========================local variables=======================
    integer (kind = myint) :: m,i,l1,l2,j,n,jc
    real (kind = myreal) :: xx
    ! transpose: whether it is a matrix transpose multiplying a vector
    logical :: transpose

    transpose = .false.

    if (present(trans)) then
       if (trans) transpose = .true.
    end if

    info = 0

    m = matrix%m; n = matrix%n

    if (ZD11_get(matrix%type) == "pattern") then
       info = MC65_ERR_MATVEC_NOVALUE
       return
    end if
    if (transpose) then
       y(1:n) = 0
       do i = 1, m
          xx = x(i)
          do j = matrix%ptr(i),matrix%ptr(i+1)-1
             jc = matrix%col(j)
             y(jc) = y(jc) + matrix%val(j)*xx
          end do
       end do
    else
       do i = 1, m
          l1 = matrix%ptr(i); l2 = matrix%ptr(i+1)-1
          y(i) = dot_product(matrix%val(l1:l2),x(matrix%col(l1:l2)))
       end do
    end if

  end subroutine csr_matrix_multiply_rvector


  subroutine csr_matrix_multiply_ivector(matrix,x,y,info,trans)
    ! subroutine csr_matrix_multiply_ivector(matrix,x,y,info)
    !
    ! y = matrix*x where x and y are integer vectors. Entries of
    ! matrix is assumed to be one. Dimension of y
    ! is checked and returned if it is smaller than the row dimension
    ! of x
    !
    ! if trans is present and == 1;
    ! y = matrix^T*x where x and y are integers and entries of
    ! the matrix are assumed to have the values one.

    ! matrix: of the derived type ZD11_type, INTENT (IN),
    !         the sparse matrix in compressed sparse row format
    type (zd11_type), intent (in) :: matrix

    ! x: integer array of intent (IN), a vector to be
    !    multiplied with the matrix
    integer (kind = myint), intent (in), dimension (*) :: x

    ! y: integer array of intent (OUT), the result of
    !    matrix*x or matrix^T*x
    integer (kind = myint), intent (out), dimension (*) :: y

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    integer (kind = myint), intent (out) :: info

    ! trans: optional integer scalar. If present and trans = true,
    ! transpose of the matrix is multiplied with the vector
    logical, optional, intent (in) :: trans

    ! local ==========
    integer (kind = myint) :: m,n,i,l1,l2,j,xx,jc
    ! transpose: whether it is a matrix transpose
    !   multiplying a vector
    logical :: transpose


    info = 0

    transpose = .false.

    if (present(trans)) then
       if (trans) transpose = .true.
    end if

    m = matrix%m; n = matrix%n

    if (transpose) then
       y(1:n) = 0
       do i = 1, m
          xx = x(i)
          do j = matrix%ptr(i),matrix%ptr(i+1)-1
             jc = matrix%col(j)
             y(jc) = y(jc) + xx
          end do
       end do
    else
       do i = 1, m
          l1 = matrix%ptr(i); l2 = matrix%ptr(i+1)-1
          y(i) = sum(x(matrix%col(l1:l2)))
       end do
    end if
  end subroutine csr_matrix_multiply_ivector

  subroutine csr_to_csr_matrix(matrix,m,n,ptr,col,&
       info,val,type,checking,stat)
    !   subroutine csr_to_csr_matrix(matrix,m,n,ptr,col,
    !      info[,val,type,copying,checking])
    ! convert ptr, col and a arrays into the compressed sparse row
    ! matrix format. New storage is allocated for MATRIX and
    ! the content of ptr,col (and val) is copied into MATRIX.

    ! m: is an integer of intent (in).
    !    It is the row dimension of the sparse matrix
    ! n: is an integer of intent (in).
    !    It is the column dimension of the sparse matrix
    integer (kind = myint), intent (in) :: m,n

    ! MATRIX: is of the derived type {\tt ZD11\_type}
    !     with {\tt INTENT (INOUT)},
    type (zd11_type), intent (inout) :: matrix

    ! ptr: is an integer array of intent (in) of size M+1.
    !     It is the row pointer of the sparse matrix.
    ! col: is an integer array of intent (in) of size ptr(m+1)-1.
    !     It is the array of column indices.
    integer (kind = myint), INTENT (IN)  :: ptr(m+1), col(ptr(m+1)-1)


    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation fails.
    !       = MC65_ERR_MEMORY_deALLOC if memory deallocation fails.
    !       = MC65_ERR_RANGE_PTR if PTR(1:M+1) is not
    !         monotonically increasing
    !       = MC65_WARN_RANGE_COL: COL is not within the
    !         range of [1,n], such entries are removed
    !       = MC65_WARN_DUP_ENTRY: duplicate entries were
    !         found when cleaning the matrix and were summed
    integer (kind = myint), intent (out) :: info

    ! val: is an optional REAL array of intent (in) of
    !    size ptr(m+1)-1.
    !    It holds the entry values of the sparse matrix.
    !    If this argument is not present, the matrix is
    !    taken as pattern only.
    real (kind = myreal),  INTENT (IN), optional :: val(ptr(m+1)-1)


    ! type: an optional character array of intent (IN).
    !       If both type and A is present,
    !       the type of the MATRIX will be set as TYPE.
    !
    character (len = *), INTENT (IN), OPTIONAL :: type


    ! checking: is an optional integer scaler of INTENT (IN).
    !           If present and is
    !           1: the input data will be checked for consistency,
    !              the matrix will be cleaned by merging duplicate
    !              entries. If a column index is out of range,
    !              the entry will be excluded in the construction
    !              of the matrix.
    !              A warning will be recorded.
    !           < 0: no checking or cleaning is done
    !           Other values: current not used and has the
    !           same effect 1
    integer (kind = myint), parameter :: CHECK_RMV = 1
    integer (kind = myint), optional, INTENT (IN) :: checking

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat


    integer (kind = myint) :: nz,nzz,i,j,icol,ierr


    ! whether the matrix is pattern only
    logical :: pattern_only

    ! whether to check input data
    integer (kind = myint) :: to_check

    nzz = ptr(m+1)-1

    info = 0
    if (present(stat)) stat = 0
    to_check = 0
    if (present(checking)) then
       if (checking == CHECK_RMV) then
          to_check = CHECK_RMV
       else if (checking >= 0) then
          to_check = CHECK_RMV
       end if
    end if


    pattern_only = .true.
    if (present(val)) then
       pattern_only = .false.
       if (present(type)) then
          if (len_trim(type)>=7) then
             if (type(1:7) == "pattern".or.type(1:7) == "PATTERN") &
                  pattern_only = .true.
          end if
       end if
    end if

    if (to_check == CHECK_RMV) then
       do i = 1,m
          if (ptr(i) > ptr(i+1)) then
             info = MC65_ERR_RANGE_PTR
             return
          end if
       end do
       nz = 0
       do i = 1, nzz
          if (col(i) <= n.and.col(i) >= 1) then
             nz = nz + 1
          end if
       end do
    else
       nz = nzz
    end if


    if (pattern_only) then
       call csr_matrix_construct(matrix,m,nz,info,n = n,&
             type = "pattern",stat=ierr)
    else
       if (present(type)) then
          call csr_matrix_construct(matrix,m,nz,info,n = n,&
               type = type,stat=ierr)
       else
          call csr_matrix_construct(matrix,m,nz,info,n = n,&
               type = "general",stat=ierr)
       end if
    end if
    if (present(stat)) stat = ierr
    if (info < 0) return

    ! column index outside the range
    matrix%ptr(1) = 1
    if (nz /= nzz.and.to_check == CHECK_RMV) then
       nz = 1
       if (.not.pattern_only) then
          do i = 1, m
             do j = ptr(i), ptr(i+1)-1
                icol = col(j)
                if (icol >= 1.and.icol <= n) then
                   matrix%col(nz) = icol
                   matrix%val(nz) = val(j)
                   nz = nz + 1
                end if
             end do
             matrix%ptr(i+1) = nz
          end do
       else
          do i = 1, m
             do j = ptr(i), ptr(i+1)-1
                icol = col(j)
                if (icol >= 1.and.icol <= n) then
                   matrix%col(nz) = icol
                   nz = nz + 1
                end if
             end do
             matrix%ptr(i+1) = nz
          end do
       end if
       call csr_matrix_clean(matrix,info,stat=ierr)
       if (present(stat)) stat = ierr
       if (info < 0) return
       ! only set info for COL out of range, since info is reset when
       ! calling csr_matrix_construct or csr_matrix_clean
       if (info == MC65_WARN_DUP_ENTRY) then
          info = MC65_WARN_ENTRIES
       else
          info =  MC65_WARN_RANGE_COL
       end if
    else if (to_check == CHECK_RMV) then
       ! column index not outside the range but cleaning is needed
       matrix%ptr(1:m+1) = ptr(1:m+1)
       matrix%col(1:nz) = col(1:nz)
       if (.not.pattern_only) matrix%val(1:nz) = val(1:nz)
       call csr_matrix_clean(matrix,info,stat=ierr)
       if (present(stat)) stat = ierr
    else
       matrix%ptr(1:M+1) = ptr(1:m+1)
       matrix%col(1:nz) = col(1:nz)
       if (.not.pattern_only) matrix%val(1:nz) = val(1:nz)
    end if


  end subroutine csr_to_csr_matrix


  subroutine coo_to_csr_format(matrix,m,n,nz,irn,jcn,&
       info,val,type,checking,stat)
    !
    ! convert coordinate format (irn,jcn,val) to CSR matrix format

    ! matrix:  is of the derived type {\tt ZD11\_type}
    !          with {\tt INTENT (INOUT)}. It is
    !          the sparse matrix to be generated from the
    !          coordinated formatted
    !          input data.
    ! m: is an INTEGER scaler of INTENT (IN). It holds the row size
    ! n: is an INTEGER scaler of INTENT (IN). It holds the
    !    column size
    ! nz:  is an INTEGER scaler. It holds the number of
    !      nonzeros in the sparse matrix
    ! irn: is an INTEGER ARRAY with INTENT (IN), of size = NZ.
    !      It holds the row indices of the matrix
    ! jcn:  is an INTEGER ARRAY with INTENT (IN), of size = NZ.
    !       It holds column indices of the matrix
    ! info:  is an {\tt INTEGER} scaler of {\tt INTENT (OUT)}.
    !       INFO = 0 if the subroutine completes successful
    !       INFO = MC65_ERR_MEMORY_ALLOC if memory allocation failed.
    !       INFO = MC65_ERR_MEMORY_DEALLOC if memory
    !              deallocation failed.
    !       {\tt INFO = MC65\_WARN\_RANGE\_IRN} if {\tt IRN}
    !           is not within the
    !           range of {\tt [1,M]}. Related entries excluded
    !       {\tt INFO = MC65\_WARN\_RANGE\_JCN} if {\tt JCN}
    !           is not within the
    !           range of {\tt [1,N]}. Related entries excluded
    !       INDO = MC65_WARN_DUP_ENTRY: duplicate entries were
    !              found when cleaning the matrix
    !              and were summed
    ! val:  is an optional REAL (double precision REAL for
    !       HSL_MC65_double)
    !       ARRAY  with INTENT (IN), of size = NZ.
    !       It holds entry values of the matrix
    integer (kind=myint), intent (in) :: m, n, nz
    type (zd11_type), intent (inout) :: matrix
    integer (kind=myint), intent (in) :: irn(nz)
    integer (kind=myint), intent (in) :: jcn(nz)
    integer (kind=myint), intent (out) :: info
    real (kind = myreal), intent (in), optional :: val(nz)

    ! type: is an OPTIONAL CHARACTER array of unspecified
    !       size with INTENT (IN).
    !       If both val and type are present,
    !       the type of the MATRIX will be set as type.
    character (len=*), intent (in), optional :: type

    ! checking: is an optional integer scaler of INTENT (IN).
    !           If present and is
    !           1: the input data will be checked for consistency,
    !              the matrix will be cleaned by merging duplicate
    !              entries.
    !              If a column index is out of range, the entry will
    !              be excluded in the construction of the matrix.
    !              A warning will be recorded.
    !           < 0: no checking or cleaning is done
    !              Other values: current not used and has the same
    !              effect as 1
    integer (kind = myint), parameter :: CHECK_RMV = 1
    integer (kind = myint), optional, INTENT (IN) :: checking

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! nzz: actual number of nonzero

    integer (kind=myint) :: ierr,nzz,i
    integer (kind = myint), dimension (:), allocatable :: nrow,iia

    ! whether the matrix is pattern only
    logical :: pattern_only

    ! whether to check input data. 0 not, /= 0 yes
    integer (kind = myint) :: to_check

    logical irn_out_range,jcn_out_range

    irn_out_range = .false.
    jcn_out_range = .false.

    info = 0
    if (present(stat)) stat = 0
    to_check = 0
    if (present(checking)) then
       if (checking == CHECK_RMV) then
          to_check = CHECK_RMV
       else if (checking >= 0) then
          to_check = CHECK_RMV
       end if
    end if

    pattern_only = .true.
    if (present(val)) then
       pattern_only = .false.
       if (present(type)) then
          if (len_trim(type)>=7) then
             if (type(1:7) == "pattern".or.type(1:7) == "PATTERN") &
                  pattern_only = .true.
          end if
       end if
    end if

    if (to_check == CHECK_RMV) then
       nzz = 0
       do i = 1, nz
          if (jcn(i) <= n.and.jcn(i) >= 1) then
             if (irn(i) <= m.and.irn(i) >= 1) then
                nzz = nzz + 1
             else
                irn_out_range = .true.
             end if
          else
             jcn_out_range = .true.
          end if
       end do
    else
       nzz = nz
    end if

    if (pattern_only) then
       call csr_matrix_construct(matrix,m,nzz,info,n = n,&
            type = "pattern",stat=ierr)
    else
       if (present(type)) then
          call csr_matrix_construct(matrix,m,nzz,info,n = n,&
               type = type,stat=ierr)
       else
          call csr_matrix_construct(matrix,m,nzz,info,n = n,&
               type = "general",stat=ierr)
       end if
    end if
    if (present(stat)) stat = ierr
    if (info < 0) return

    allocate(nrow(m),iia(m+1), stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

    if (pattern_only) then
       call coo_to_csr_private(to_check,m,n,nz,irn,jcn,matrix%ptr,matrix%col,&
            nrow,iia,pattern_only)
    else
       call coo_to_csr_private(to_check,m,n,nz,irn,jcn,matrix%ptr,matrix%col,&
            nrow,iia,pattern_only,val = val,a = matrix%val)
    end if

    deallocate(nrow, iia, stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if

    if (to_check == CHECK_RMV) then
       call csr_matrix_clean(matrix,info,stat=ierr)
       if (present(stat)) stat = ierr
       if (info < 0) return
       if (irn_out_range .and. info==MC65_WARN_DUP_ENTRY) then
          info = MC65_WARN_ENTRIES
          return
       end if
       if (jcn_out_range .and. info==MC65_WARN_DUP_ENTRY) then
          info = MC65_WARN_ENTRIES
          return
       end if
       if (irn_out_range) info = MC65_WARN_RANGE_IRN
       if (jcn_out_range) info = MC65_WARN_RANGE_JCN
       if (irn_out_range .and. jcn_out_range) info = MC65_WARN_RANGE_BOTH
    end if

  end subroutine coo_to_csr_format


  subroutine coo_to_csr_private(to_check,m,n,nz,irn,jcn,ia,ja,&
       nrow,iia,pattern_only,val,a)
    ! wrapper routine for coo_to_csr_format for performance
    ! to_check: whether to check input data. 0 not, /= 0 yes
    ! m,n: row/column dimension
    ! ia: row pointer
    ! ja: column pointer
    ! nrow: number of entries in a row
    ! iia: row pointer
    ! a: entry values
    integer (kind = myint) :: to_check
    integer (kind = myint), intent (in) :: n,m,nz,irn(*),jcn(*)
    integer (kind = myint), intent (out) :: ia(*),ja(*),nrow(*),iia(*)
    real (kind = myreal), intent (OUT), optional :: a(*)
    real (kind = myreal), intent (IN), optional :: val(*)
    logical :: pattern_only

    integer (kind=myint) :: i
    integer (kind = myint) :: ii,jj

    if (to_check == 0) then
       nrow(1:m) = 0
       do i = 1, nz
          nrow(irn(i)) = nrow(irn(i)) + 1
       end do
       ia (1) = 1
       do i = 1, m
          ia(i+1) = ia(i) + nrow(i)
       end do

       iia(1:m+1) = ia(1:m+1)
       if (pattern_only) then
          do i = 1, nz
             ii = irn(i)
             jj = jcn (i)
             ja(iia(ii)) = jj
             iia(ii) = iia(ii) + 1
          end do
       else
          do i = 1, nz
             ii = irn(i)
             jj = jcn (i)
             ja(iia(ii)) = jj
             a(iia(ii)) = val(i)
             iia(ii) = iia(ii) + 1
          end do
       end if
    else
       ! checking range, exclude out of range entries
       nrow(1:m) = 0
       do i = 1, nz
          ii = irn(i)
          jj = jcn(i)
          if (ii >= 1.and.ii <= m.and.jj >= 1.and.jj <= n) &
               nrow(ii) = nrow(ii) + 1
       end do
       ia (1) = 1
       do i = 1, m
          ia(i+1) = ia(i) + nrow(i)
       end do

       iia(1:m+1) = ia(1:m+1)
       if (pattern_only) then
          do i = 1, nz
             ii = irn(i)
             jj = jcn (i)
             if (ii >= 1.and.ii <= m.and.jj >= 1.and.jj <= n) then
                ja(iia(ii)) = jj
                iia(ii) = iia(ii) + 1
             end if
          end do
       else
          do i = 1, nz
             ii = irn(i)
             jj = jcn (i)
             if (ii >= 1.and.ii <= m.and.jj >= 1.and.jj <= n) then
                ja(iia(ii)) = jj
                a(iia(ii)) = val(i)
                iia(ii) = iia(ii) + 1
             end if
          end do
       end if
    end if
  end subroutine coo_to_csr_private


  subroutine csr_matrix_to_coo(matrix,m,n,nz,irn,jcn,info,val,stat)
    ! convert the ZD11_type matrix format to coordinate form.
    ! Storage space for the array pointers  {\tt IRN},
    ! {\tt JCN} and {\tt VAL}
    ! are allocated inside the subroutine, and should be deallocated
    ! by the user once these array pointers are no longer used.
    ! being used. When the matrix is of type "pattern"
    ! and VAL is present, VAL will be allocated a size of 0.

    ! matrix:  is of the derived type {\tt ZD11\_TYPE} with
    !   {\tt INTENT (IN)}.
    type (zd11_type), intent (in) :: MATRIX

    ! m: is an INTEGER scaler of INTENT (OUT). It is the row
    !    dimension
    ! n: is an INTEGER scaler of INTENT (OUT). It is the column
    !    dimension
    ! nz: is an INTEGER scaler of INTENT (OUT). It is the number
    !    of nonzeros
    integer (kind = myint), intent (out) :: m,n,nz

    ! irn: is an INTEGER allocatable array. It holds the row indices
    ! jcn:  is an INTEGER allocatable array. It holds the column indices
    ! val:  is an OPTIONAL INTEGER allocatable array. When present,
    !       and the matrix is not of type "pattern", VAL holds
    !       the entry values
    !       when present and the matrix is of type "pattern"
    !       VAL will be allocated a size of 0.
    integer (kind = myint), allocatable :: irn(:)
    integer (kind = myint), allocatable :: jcn(:)
    real (kind = myreal), optional, allocatable :: val(:)

    ! info:  is an {\tt INTEGER} scaler of {\tt INTENT (OUT)}.
    !       INFO = 0 if the subroutine completes successful
    !       INFO = MC65_ERR_MEMORY_ALLOC if memory allocation failed.
    integer (kind = myint), intent (out) :: info

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! working integers
    integer (kind = myint) :: ierr

    ! pattern_only: whether the matrix is pattern only
    ! entry_values: whether val array is to be allocated and filled
    logical :: pattern_only,entry_values

    info = 0
    if (present(stat)) stat = 0
    pattern_only = .false.
    entry_values = .false.
    if (ZD11_get(matrix%type) == "pattern") pattern_only = .true.
    if ((.not.pattern_only).and.present(val)) entry_values = .true.

    m = matrix%m; n = matrix%n

    nz = matrix%ptr(m+1)-1

    if (nz < 0) return
    if (entry_values) then
       allocate (irn(1:nz),jcn(1:nz),val(1:nz),stat = ierr)
    else
       allocate (irn(1:nz),jcn(1:nz),stat = ierr)
       IF (PRESENT(VAL)) allocate(val(0),stat = ierr)
    end if
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       INFO = MC65_ERR_MEMORY_ALLOC
       return
    end if

    if (entry_values) then
       call csr_matrix_to_coo_private(m,nz,matrix%ptr,matrix%col,irn,jcn, &
                                      matrix%val,val)
    else
       call csr_matrix_to_coo_private(m,nz,matrix%ptr,matrix%col,irn,jcn)
    end if

  end subroutine csr_matrix_to_coo

  subroutine csr_matrix_to_coo_private(m,nz,ia,ja,irn,jcn,a,val)
    ! private wrapper routine called by csr_matrix_to_coo
    ! for efficient
    !
    ! m: row dimension
    ! nz: number of nonzeros
    ! ia: row pointer
    ! ja: column indices
    ! a: entry values, optional
    ! irn: row indices
    ! jcn: column indices
    ! val: entry values, optional
    integer (kind = myint), intent (in) :: m
    integer (kind = myint), intent (out) :: nz
    integer (kind = myint), dimension (*), intent (in) :: ia,ja
    real (kind = myreal), dimension (*), intent (in), optional :: a
    integer (kind = myint), dimension (*), intent (out) :: irn,jcn
    real (kind = myreal), dimension (*), intent (out), optional :: val

    integer (kind = myint) :: i,n1,i1,i2


    nz = 1
    if (present(a)) then
       do i=1,m
          i1 = ia(i); i2 = ia(i+1)-1
          n1 = i2-i1+1
          irn(nz:nz+n1-1) = i
          jcn(nz:nz+n1-1) = ja(i1:i2)
          val(nz:nz+n1-1) = a(i1:i2)
          nz = nz + n1
       end do
    else
       do i=1,m
          i1 = ia(i); i2 = ia(i+1)-1
          n1 = i2-i1+1
          irn(nz:nz+n1-1) = i
          jcn(nz:nz+n1-1) = ja(i1:i2)
          nz = nz + n1
       end do
    end if
    nz = nz -1
  end subroutine csr_matrix_to_coo_private

  subroutine csr_matrix_remove_diagonal(matrix)
    ! remove the diagonal entry of a matrix. This is used when
    ! the matrix represents a graph and there should not be
    ! diagonals in such a matrix

    ! matrix: is of the derived type {\tt ZD11\_TYPE} with
    !         {\tt INTENT (INOUT)}.
    !         It is the matrix whose diagonal is to be removed
    type (zd11_type), intent (inout) :: matrix

    ! m: row dimension
    integer (kind = myint) :: m

    m = matrix%m
    if (ZD11_get(matrix%type) /= "pattern" ) then
       call csr_matrix_remove_diag_private(m,matrix%ptr,matrix%col,matrix%val)
    else
       call csr_matrix_remove_diag_private(m,matrix%ptr,matrix%col)
    end if
  end subroutine csr_matrix_remove_diagonal

  subroutine csr_matrix_remove_diag_private(m,ia,ja,a)
    ! private  subroutine called by csr_matrix_remove_diag
    ! for performance reason
    !
    ! m: row dimension
    ! ia: row pointer
    ! ja: column indices
    ! a: entry values
    integer (kind = myint), intent (in) :: m
    integer (kind = myint), intent (inout) :: ia(*),ja(*)
    real (kind = myreal), optional, intent (inout) :: a(*)

    integer (kind = myint) :: i,j,l1,nz


    nz = 1
    l1 = 1
    if (present(a)) then
       do i = 1, m
          do j = l1, ia(i+1) - 1
             if (ja(j) /= i) then
                ja(nz) = ja(j)
                a(nz) = a(j)
                nz = nz + 1
             end if
          end do
          l1 = ia(i+1)
          ia(i+1) = nz
       end do
    else
       do i = 1, m
          do j = l1, ia(i+1) - 1
             if (ja(j) /= i) then
                ja(nz) = ja(j)
                nz = nz + 1
             end if
          end do
          l1 = ia(i+1)
          ia(i+1) = nz
       end do
    end if
  end subroutine csr_matrix_remove_diag_private


  subroutine csr_matrix_diagonal_first(matrix,info)
    ! subroutine csr_matrix_diagonal_first(matrix,info)
    !
    ! changes the storage arrange so that the diagonal of
    ! row i of the matrix is always
    ! at position matrix%ptr(i) of the matrix%col and matrix%val
    !    arrays

    ! matrix:  is of the derived type {\tt ZD11\_TYPE} with
    !    {\tt INTENT (INOUT)}.
    type (zd11_type), intent (inout) :: matrix

    ! INFO: an INTEGER scaler of INTENT (OUT).
    !       INFO = 0 if the subroutine completes successfully; and
    !       INFO = MC65_WARN_MV_DIAG_MISSING if some diagonal
    !       elements are missing.
    integer (kind = myint), intent (OUT) :: info

    integer (kind = myint) :: m

    info = 0

    m = matrix%m
    if (zd11_get(matrix%type) == "pattern") then
       call csr_matrix_diagonal_first_priv(m,matrix%ptr,matrix%col,info)
    else
       call csr_matrix_diagonal_first_priv(m,matrix%ptr,matrix%col,info, &
                                           matrix%val)
    end if
  end subroutine csr_matrix_diagonal_first

  subroutine csr_matrix_diagonal_first_priv(m,ia,ja,info,a)
    ! private routine for csr_matrix_diagonal_first for efficiency

    ! ia: row pointer to matrix
    ! ja: column indices
    integer (kind = myint), dimension (*), intent (inout) :: ia,ja
    ! a: entry values
    real (kind = myreal), dimension (*), intent (inout), optional :: a

    integer (kind = myint) :: m,i,j,jdiag,jf,info
    real (kind = myreal) :: tmp

    if (present(a)) then
       do i=1,m
          jdiag = -1
          jf = ia(i)
          do j = jf,ia(i+1)-1
             if (ja(j) == i) then
                jdiag = j
                exit
             end if
          end do
          if (jdiag < 0) then
             INFO = MC65_WARN_MV_DIAG_MISSING
             cycle
          end if
          ! swap the j-th and the jdiag-th entry in the row
          tmp = a(jdiag)
          a(jdiag) = a(jf)
          a(jf) = tmp
          ja(jdiag) = ja(jf)
          ja(jf) = i
       end do
    else
       do i=1,m
          jdiag = -1
          jf = ia(i)
          do j = jf,ia(i+1)-1
             if (ja(j) == i) then
                jdiag = j
                exit
             end if
          end do
          if (jdiag < 0) then
             INFO = MC65_WARN_MV_DIAG_MISSING
             cycle
          end if
          ja(jdiag) = ja(jf)
          ja(jf) = i
       end do
    end if
  end subroutine csr_matrix_diagonal_first_priv


  subroutine csr_matrix_write(matrix,file_name,info,form,row_ord,&
       col_ord,stat)
    ! write the matrix in the file <file_name> with
    ! various <form>. If row_ord or/and col_ord is present,
    ! the reordered matrix will be written
    !  (unless form = "unformatted"
    ! or form = "hypergraph", in which case no reordering will
    !  take place).

    ! matrix: is of the derived type {\tt ZD11\_TYPE} with
    !         {\tt INTENT (IN)}.
    !         It is the matrix to be written.
    type (zd11_type), intent (in) :: matrix

    ! file_name: is a CHARACTER array of assumed size with
    !   {\tt INTENT (IN)}.
    ! It is the file to write the matrix in.
    character (len = *), intent (in) :: file_name

    ! form: an OPTIONAL CHARACTER array of assumed size with
    !         {\tt INTENT (IN)}. It
    !         is the format in which the matrix is to be written.
    !         Not supplying <form>
    !         or a <form> that is not one of the following strings
    !         has the effect of
    !         form = "ija".
    !         In the follow let M, N and NZ be the row, column
    !         indices and the number of entries in the matrix.
    !         TYPE be the type of the matrix as a character array
    !         with trailing blanks removed.
    !         LEN_TYPE be the string length of TYPE.
    ! form =
    !       "gnuplot": the matrix is written in <file_name> and
    !                the gnuplot command in <file_name.com>,
    !                then in gnuplot, you can type
    !                "load 'file_name'" to view the matrix
    !       "hypergraph": the matrix is written in <file_name>
    !                 in hypergraph format.
    !                       N,M
    !                       row_indices_of_column_1
    !                       row_indices_of_column_2
    !                       row_indices_of_column_3
    !                       ...
    !
    !
    !       "unformatted": the matrix is written in <file_name>
    !            as unformatted data:
    !                         LEN_TYPE
    !                         TYPE
    !                         M N NZ
    !                         row_pointers
    !                         column_indices
    !                         entry_values (Only if the matrix is
    !                         not of type "pattern")
    !
    !       "ij": the matrix is written in <file_name> in
    !             coordinate format,
    !             with only row and column indices. That is
    !                         M N NZ
    !                         ....
    !                         a_row_index a_column_index
    !                         another_row_index another_column_index
    !                         ...
    !
    !       "ija" or "coordinate": the matrix is written in
    !              <file_name> in coordinate format.
    !              Note that entry values will be written only
    !              if the matrix is
    !              not of type "pattern"
    !               M N NZ
    !               ....
    !               row_index column_index a(row_index,column_index)
    !               ....
    !
    character (len = *), intent (in), optional :: form

    ! INFO: an INTEGER of INTENT (OUT).
    !       INFO = 0 if the subroutine completes successfully
    !       INFO = MC65_ERR_NO_VACANT_UNIT if no vacant unit has
    !              been found
    !              in MC65_MATRIX_WRITE
    !       INFO = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       INFO = MC65_ERR_MEMORY_DEALLOC if memory
    !              deallocation failed
    INTEGER (KIND = MYINT), INTENT (OUT) :: INFO

    ! row_ord: is an OPTIONAL INTEGER array of size M with
    !          INTENT (IN). When present, a row index
    !          I will become ROW_ORD(I) in the reordered
    !          matrix.
    ! col_ord:  is an OPTIONAL INTEGER array of size N with
    !          INTENT (IN). When present, a column index
    !          I will become COL_ORD(I) in the reordered
    !          matrix.
    integer (kind = myint), dimension (matrix%m), intent (in), &
         optional :: row_ord
    integer (kind = myint), dimension (matrix%n), intent (in), &
         optional :: col_ord

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat


    ! =============== local variables ==============

    ! pattern_only: whether the matrix is of type "pattern"
    logical :: pattern_only

    integer (kind = myint) :: m,n,nz,ierr

    if (zd11_get(matrix%type) == "pattern") then
       pattern_only = .true.
    else
       pattern_only = .false.
    end if

    info = 0;
    if (present(stat)) stat = 0

    m = matrix%m
    n = matrix%n
    nz = matrix%ptr(m+1)-1

    if (.not.present(form)) then
       goto 10
    else if (form == "gnuplot".or.form == "GNUPLOT") then
       if (present(row_ord).and.present(col_ord)) then
          call csr_matrix_write_gnuplot(file_name,m,n,matrix%ptr,matrix%col,&
               info,row_ord,col_ord)
       else if (present(row_ord)) then
          call csr_matrix_write_gnuplot(file_name,m,n,matrix%ptr,matrix%col,&
               info,row_ord = row_ord)
       else if (present(col_ord)) then
          call csr_matrix_write_gnuplot(file_name,m,n,matrix%ptr,matrix%col,&
               info,col_ord = col_ord)
       else
          call csr_matrix_write_gnuplot(file_name,m,n,matrix%ptr,matrix%col,&
               info)
       end if
       return
    else if (form == "hypergraph" .or. form == "HYPERGRAPH") then
       call csr_matrix_write_hypergraph(matrix,file_name,info,stat=ierr)
       if (present(stat)) stat = ierr
       return
    else if (form == "unformatted".or.form == "UNFORMATTED") then
       if (pattern_only) then
          call csr_matrix_write_unformatted(file_name,m,n,nz,&
               matrix%type,matrix%ptr,matrix%col,info)
       else
          call csr_matrix_write_unformatted(file_name,m,n,nz,&
               matrix%type,matrix%ptr,matrix%col,info,matrix%val)
       end if
       return
    else if (form == "ij".or.form == "IJ") then
       pattern_only = .true.
       goto 10
    else if (form == "ija".or.form == "coordinate".or.&
         form == "IJA".or.form == "COORDINATE") then
    end if

10  if (pattern_only) then
       if (present(row_ord).and.present(col_ord)) then
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,row_ord = row_ord,col_ord = col_ord)
       else if (present(row_ord)) then
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,row_ord = row_ord)
       else if (present(col_ord)) then
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,col_ord = col_ord)
       else
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info)
       end if
    else
       if (present(row_ord).and.present(col_ord)) then
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,matrix%val,row_ord = row_ord,col_ord = col_ord)
       else if (present(row_ord)) then
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,matrix%val,row_ord = row_ord)
       else if (present(col_ord)) then
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,matrix%val,col_ord = col_ord)
       else
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,matrix%val)
       end if
    end if

  end subroutine csr_matrix_write


  subroutine CSR_MATRIX_WRITE_ija(file_name,m,n,nz,ia,ja,info,&
       a,row_ord,col_ord)
    ! file_name: is a CHARACTER array of assumed size with
    !            {\tt INTENT (IN)}.
    ! It is the file to write the matrix in.
    ! m: row dimension
    ! n: column dimension
    ! nz: number of nonzeros
    ! ia: row pointer
    ! ja: column indices
    ! info: flag to show whether successful or not
    ! a: entry values
    ! row_ord: row ordering. NewIndex = row_ord(OldIndex)
    ! col_ord: col ordering. NewIndex = col_ord(OldIndex)
    character (len = *), intent (in) :: file_name
    integer (kind = myint), intent (in) :: m,n,nz,ia(*),ja(*)
    integer (kind = myint), intent (inout) :: info
    real (kind = myreal), intent (in), optional :: a(*)
    integer (kind = myint), intent (in), optional :: row_ord(*),&
         col_ord(*)
    integer (kind = myint) :: file_unit
    integer (kind = myint) :: i,k

    file_unit = vacant_unit()
    if (file_unit < 0) then
       info = MC65_ERR_NO_VACANT_UNIT
       return
    end if
    open(unit = file_unit, file = file_name, form = "formatted")
    write(file_unit,"(i10,i10,i10)") m, n, nz

    if (present(row_ord).and.present(col_ord)) then
       if (present(a)) then
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10,g12.4)") &
                  (row_ord(i),col_ord(ja(k)),a(k),k = ia(i),ia(i+1)-1)
          end do
       else
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10)") &
                  (row_ord(i),col_ord(ja(k)),k = ia(i),ia(i+1)-1)
          end do
       end if
    else if (present(row_ord)) then
       if (present(a)) then
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10,g12.4)") &
                  (row_ord(i),ja(k),a(k),k = ia(i),ia(i+1)-1)
          end do
       else
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10)") &
                  (row_ord(i),ja(k),k = ia(i),ia(i+1)-1)
          end do
       end if
    else if (present(col_ord)) then
       if (present(a)) then
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10,g12.4)") &
                  (i,col_ord(ja(k)),a(k),k = ia(i),ia(i+1)-1)
          end do
       else
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10)") &
                  (i,col_ord(ja(k)),k = ia(i),ia(i+1)-1)
          end do
       end if
    else
       if (present(a)) then
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10,g12.4)") &
                  (i,ja(k),a(k),k = ia(i),ia(i+1)-1)
          end do
       else
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10)") (i,ja(k),k = ia(i),ia(i+1)-1)
          end do
       end if
    end if
    close(file_unit)

  end subroutine CSR_MATRIX_WRITE_ija

  subroutine csr_matrix_write_unformatted(file_name,m,n,nz,&
       type,ia,ja,info,a)
    ! file_name: is a CHARACTER array of assumed size with
    !  {\tt INTENT (IN)}.
    ! It is the file to write the matrix in.
    ! m: row dimension
    ! n: column dimension
    ! nz: number of nonzeros
    ! type: character pointer of matrix type
    ! ia: row pointer
    ! ja: column indices
    ! info: flag to show whether successful or not
    ! a: entry values
    character (len = *), intent (in) :: file_name
    integer (kind = myint), intent (in) :: m,n,nz,ia(*),ja(*)
    character, allocatable, dimension (:) :: type
    integer (kind = myint), intent (inout) :: info
    real (kind = myreal), intent (in), optional :: a(*)
    integer (kind = myint) :: file_unit

    file_unit = vacant_unit()
    if (file_unit < 0) then
       info = MC65_ERR_NO_VACANT_UNIT
       return
    end if
    open(unit = file_unit, file = file_name, form = "unformatted")
    write(file_unit) len_trim(zd11_get(type))
    write(file_unit) zd11_get(type)
    write(file_unit) m,n,nz
    write(file_unit) ia(1:m+1)
    if (nz >= 0) then
       write(file_unit) ja(1:nz)
       if (present(a)) write(file_unit) a(1:nz)
    end if
    close(file_unit)
  end subroutine csr_matrix_write_unformatted

  subroutine csr_matrix_write_gnuplot(file_name,m,n,ia,ja,&
       info,row_ord,col_ord)
    ! write a matrix in a gnuplot data file
    ! "file_name" and the corresponding gnuplot command
    ! to view the matrix in "file_name.com"
    ! The matrix is reordered using col_ord and row_ord
    ! before writing




    ! file_name: is a CHARACTER array of assumed size with
    !   {\tt INTENT (IN)}.
    ! It is the file to write the matrix in.
    ! m: row dimension
    ! n: column dimension
    ! ia: row pointer
    ! ja: column indices
    ! info: info tag.
    ! row_ord: row ordering. NewIndex = row_ord(OldIndex)
    ! col_ord: col ordering. NewIndex = col_ord(OldIndex)
    character (len = *), intent (in) :: file_name
    integer (kind = myint), intent (in) :: m,n,ia(*),ja(*)

    ! col_ord, row_ord: column and row ordering. if a
    !  column col_ord(i) = 1,
    ! column i will be the first column in the ordered matrix
    integer (kind = myint), dimension (*), intent (in), optional :: &
         col_ord, row_ord
    integer (kind = myint), intent (OUT) :: info

    integer (kind = myint) :: i,j
    integer (kind = myint) :: unit1,unit2
    character (len = 120):: fmt_col

    unit1 = vacant_unit()
    if (unit1 < 0) then
       info = MC65_ERR_NO_VACANT_UNIT
       return
    end if
    open(unit = unit1,file = file_name)

    unit2 = vacant_unit()
    if (unit2 < 0) then
       info = MC65_ERR_NO_VACANT_UNIT
       return
    end if
    open(unit = unit2,file = file_name/&
         &/".com")

    write(unit2,"(a,i12,a)") "set xrange [0:",n+1,"]"
    write(unit2,"(a,i12,a)") "set yrange [0:",m+1,"]"
    write(unit2,"(a,i12,a,i12,a)") 'set ytics ("1" ',&
         m,', "',m,'"  1 )'
    write(unit2,"(a)") 'set noxtics'

    ! use in formating for printing number of columns
    fmt_col = "(a,a,i"/&
         &/int2str(digit(n),digit(digit(n)))/&
         &/",a,i12,a)"
    write(unit2,fmt_col) 'set x2tics ("1" 1', &
         ', "', &
         n, &
         '" ', &
         n, &
         ' )'
    write(unit2,"(a)") 'set nolabel '
    !    write(unit2,"(a)") 'set term post portrait '
    !    write(unit2,"(a,a)") 'set out "',file_name/&
    !         &/'.ps"'
    ! to set x and y equal scale
    write(unit2,"(a)") 'set size  square'

    write(unit2,"(a,a,a)") 'plot "',file_name,&
         '" using 1:2 notitle w p ps 0.3'

    if (present(row_ord).and.present(col_ord)) then
       write(unit1,'(2i10)') ((col_ord(ja(j)), m+1-row_ord(i), &
            j=ia(i), ia(i+1)-1), i = 1, m)
    else if (present(row_ord)) then
       write(unit1,'(2i10)') ((ja(j), m+1-row_ord(i), &
            j=ia(i), ia(i+1)-1), i = 1, m)
    else if (present(col_ord)) then
       write(unit1,'(2i10)') ((col_ord(ja(j)), m+1-i, &
            j=ia(i), ia(i+1)-1), i = 1, m)
    else
       write(unit1,'(2i10)') ((ja(j), m+1-i, j=ia(i), ia(i+1)-1), &
            i = 1, m)
    end if

    close(unit1)
    close(unit2)
  end subroutine csr_matrix_write_gnuplot



  subroutine csr_matrix_write_hypergraph(matrix,file_name,info,stat)
    ! write the matrix in hypergraph format.
    ! that is,
    ! row_dimension_of_A^T column_dimension_of_A^T
    ! column_indices_of_row1_of_A^T
    ! column_indices_of_row2_of_A^T
    ! column_indices_of_row3_of_A^T
    ! .....


    ! matrix: the matrix to be written
    type (zd11_type), intent (in) :: matrix

    ! file_name: the file name to write the data into
    character (len = *), intent (in) :: file_name

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    integer (kind = myint) :: i,n,m,l1,l2,ierr
    integer (kind = myint) :: unit1
    type (zd11_type) :: matrix2
    character*80 fmt
    integer (kind = myint), intent (inout) :: info

    if (present(stat)) stat = 0

    call csr_matrix_transpose(matrix,matrix2,info,stat=ierr)
    if (present(stat)) stat = ierr
    if (info < 0) return

    unit1 = vacant_unit()
    if (unit1 < 0) then
       info = MC65_ERR_NO_VACANT_UNIT
       return
    end if
    open(unit = unit1,file = file_name)

    m = matrix%m
    n = matrix%n
    write(unit1,"(2i12)") n,m
    do i = 1, n
       l1 = matrix2%ptr(i); l2 = matrix2%ptr(i+1)-1
       write(fmt,"('(',i12,'i8)')") max(1,l2-l1+1)
       write(unit1,fmt) matrix2%col(l1:l2)
    end do

    call csr_matrix_destruct(matrix2,info,stat=ierr)
    if (present(stat)) stat = ierr
    if (info < 0) return
    close(unit1)
  end subroutine csr_matrix_write_hypergraph








  subroutine csr_matrix_read(matrix,file_name,info,form,stat)
    ! read the matrix as dumped by csr_matrix_write,
    !  unformatted only

    ! matrix:  is of the derived type {\tt ZD11\_TYPE}
    !  with {\tt INTENT (INOUT)}.
    type (zd11_type), intent (inout) :: matrix

    ! file_name:  is a {\tt CHARACTER} array of assumed size
    !              with {\tt INTENT (IN)}.
    !              It is the file to write the matrix in.
    character (len = *), intent (in) :: file_name

    ! INFO: is a INTEGER of INTENT (OUT).
    !       INFO = 0 if the subroutine returns successfully
    ! info = MC65_mem_alloc if memory allocation failed
    ! info = MC65_mem_dealloc if memory deallocation failed
    !       INFO = MC65_ERR_READ_WRONGFORM if in
    !              MC65_MATRIX_READ,form is present
    !              but is of an unsupported form.
    !       INFO = MC65_ERR_READ_FILE_MISS if in MC65_MATRIX_READ,
    !              the file <file_name> does not exist.
    !       INFO = MC65_ERR_READ_OPEN if in MC65_MATRIX_READ,
    !              opening of the file <file_name> returns iostat /= 0
    !       INFO = MC65_ERR_NO_VACANT_UNIT if no vacant unit
    !              has been found in MC65_MATRIX_WRITE
    integer (kind = myint), intent (out) :: INFO


    ! form:  is an {\tt OPTIONAL CHARACTER} array of
    !        assumed size with
    !       {\tt INTENT (IN)}. It is the format of the input file
    !        Only form = "unformatted" is supported. If
    !        form is not present, by default
    !        the subroutine assume that the input file is unformatted.
    !        If form is present but form != "unformatted",
    !        the subroutine will return with info =
    !        MC65_ERR_READ_WRONGFORM
    character (len = *), intent (in), optional :: form

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat


    ! file_unit: open the file in this unit
    integer (kind = myint) :: file_unit
    ! m,n: row and column dimension
    ! nz: number of nonzeros
    integer (kind = myint) :: m,n,nz,ierr

    integer (kind = myint) :: ios
    logical :: existed, opened
    integer (kind = myint), parameter :: maxlen = 2000
    character (len=maxlen) :: type
    integer (kind = myint) :: len_type
    integer (kind = myint) :: iform, UNFORMAT = 1, UNSUPPORTED = 2

    info = 0
    if (present(stat)) stat = 0

    iform = UNFORMAT
    if (present(form)) then
       if (form == "unformatted") then
          iform = UNFORMAT
          !  else if ...
          ! .... other forms comes here
       else
          iform = UNSUPPORTED
          info = MC65_ERR_READ_WRONGFORM
          return
       end if
    end if


    file_unit = vacant_unit()
    if (file_unit < 0) then
       info = MC65_ERR_NO_VACANT_UNIT
       return
    end if

    inquire(file = file_name, exist = existed, opened = opened, &
         iostat = ios)
    if (.not.existed) then
       info = MC65_ERR_READ_FILE_MISS
       return
    end if
    if (opened) then
       rewind(file_unit)
    end if
    if (ios /= 0) then
       info = MC65_ERR_READ_OPEN
    end if

    if (iform == UNFORMAT) then
       open(unit = file_unit, file = file_name, form = "unformatted")
    else
       open(unit = file_unit, file = file_name, form = "formatted")
    end if


    if (iform == UNFORMAT) then
       type = ""
       read(file_unit) len_type
       if (len_type > maxlen) then
          INFO = MC65_ERR_READ_MAXLEN
          return
       end if
       read(file_unit) type(1:len_type)
       read(file_unit) m,n,nz
       call csr_matrix_construct(matrix,m,nz,info,n=n,type=type,stat=ierr)
       if (present(stat)) stat = ierr
       if (info < 0) then
          return
       end if
       read(file_unit) matrix%ptr(1:matrix%m+1)
       if (nz >= 1) then
          read(file_unit) matrix%col(1:nz)
          if (ZD11_get(matrix%type) /= "pattern" ) &
               read(file_unit) matrix%val(1:nz)
       end if

       ! else if (...)
       ! other format will be here ...
    end if

    close(file_unit)

  end subroutine csr_matrix_read



  subroutine csr_matrix_condense(matrix,info,col_wgt_real,&
       col_wgt_int,realloc,iremove,stat)
    !   subroutine csr_matrix_condense(matrix,info[,&
    !        col_wgt_real,col_wgt_int,realloc])
    ! given a matrix, finds all columns that has
    ! <= iremove (by default 0) element and delete them, also
    ! merge columns with exactly the same pattern.
    ! it only works on matrix with pattern only.
    ! The matrix should also be cleaned using
    ! matrix_clean to agglomerate repeated
    ! entries!

    ! MATRIX is of the derived type {\tt ZD11\_TYPE} with
    !        {\tt INTENT (INOUT)}. on exit all columns that
    !        has no more than one entry are deleted;
    !        columns that have the same pattern are merged into one.
    !        The column dimension of the matrix is changed
    !        accordingly.
    type (zd11_type), intent (inout) :: matrix


    ! INFO  is an {\tt INTEGER} of {\tt INTENT (OUT)}.
    !       INFO = 0 if the subroutine completes successfully
    !       INFO = MC65_ERR_CONDENSE_NOPAT if the matrix
    !              is not of type {\tt "pattern"}
    !       INFO = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       INFO = MC65_ERR_MEMORY_DEALLOC if memory
    !              deallocation failed.
    integer (kind = myint), intent (out) :: info

    ! col_wgt_real: is an OPTIONAL REAL array with INTENT (INOUT) of
    !          size N (the column dimension). When present,
    !          On entry, the first N elements contain the column
    !          weights of the matrix,
    !          on exit, the first NSUP elements of col_wgt
    !          holds the super-column weight of the condensed matrix.
    !          Here nsup is the new column dimension of the
    !          condensed matrix. Note that the two optional arguments
    !          col_wgt_real and col_wgt_int are supplied so that
    !          the user is free to treat column weight either as
    !          real arrays, or as integer arrays.
    !          The correct usage is therefore
    !          to supply, if appropriate, either col_wgt_real or
    !          col_wgt_int as a keyword argument
    !          ({\tt col_wgt_real = colwgt} or
    !          (\tt col_wgt_int = colwgt},
    !          depending on whether colwgt is real or integer)
    ! col_wgt_int: is an OPTIONAL INTEGER array with INTENT (INOUT)
    !          of size N  (the column dimension). When present,
    !          On entry, the first N elements contain the column
    !           weights of the matrix,
    !          on exit, the first NSUP elements of col_wgt
    !          holds the super-column weight of the condensed matrix.
    !          Here nsup is the new column dimension of the
    !          condensed matrix.
    !          Note that the two optional arguments
    !          col_wgt_real and col_wgt_int are supplied so that
    !          the user is free to treat column weight either as
    !          real arrays, or as integer arrays.
    !          The correct usage is therefore
    !          to supply, if appropriate, either col_wgt_real or
    !          col_wgt_int as a keyword argument
    !          ({\tt col_wgt_real = colwgt} or
    !          (\tt col_wgt_int = colwgt},
    !          depending on whether colwgt is real or integer).
    integer (kind = myint), dimension (matrix%n), optional, &
         intent (inout) :: col_wgt_int
    real (kind = myreal), dimension (matrix%n), optional, &
         intent (inout) :: col_wgt_real


    ! realloc:  is an optional real scaler of INTENT (IN).
    !   It is used to control the reallocation of storage.
    ! \begin{itemize}
    !  \item{} if {\tt REALLOC < 0}, no reallocation.
    !  \item{} if {\tt REALLOC == 0}, reallocation is carried out
    !  \item{} if {\tt REALLOC > 0}, reallocation iff memory
    !       saving is greater than
    !      {\tt REALLOC*100}\%. For example,
    !      if {\tt REALLOC = 0.5}, reallocation will be carried out
    !      if saving in storage is greater than 50\%
    !  \end{itemize}
    ! If {\tt REALLOC} is not present, no reallocation is carried out.
    real (kind = myreal), INTENT(IN), optional :: realloc

    ! iremove: an OPTIONAL INTEGER scalar. When present,
    !    any columns that have <= iremove
    !    entries are removed. If not present, only columns
    !    of zero entries are removed.
    integer (kind = myint), intent (in), optional :: iremove

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat
    ! ===================== local variables =======================

    ! col_wgt_new: working array. the new super-column weight.
    real (kind = myreal), dimension (:), allocatable :: col_wgt_new_r
    integer (kind = myint), dimension (:), allocatable :: &
         col_wgt_new_i

    ! f2c: the array which says which original column is
    !      merged into which super-column
    ! numsup: how many original columns make up the super-columns
    ! new: the new identity of a supercolumn
    ! flag: in which row the first appearance of a super-column?
    integer (kind = myint), dimension (:), allocatable :: &
         f2c,numsup,new,flag

    ! nsup: number of super-columns so far
    ! m,n: row and column size
    ! newid: new group of an old super-column
    ! i,j,jj: loop index
    ! scol: super-column id ofan column entry in the current row
    ! nz: number of nonzeros in the condensed matrix
    ! jjs: starting index of original row in the matrix
    integer (kind = myint) :: nsup,n,m,newid,i,j,jj,scol,nz,jjs,&
         nz_old,ierr

    ! number of super-columns after getting rid of columns
    ! with less than one entry
    integer (kind = myint) :: nsup_new

    ! col_count: a counter for (the number of entries - 1) in
    !            each column
    integer (kind = myint), dimension (:), allocatable :: col_count
    ! col_count_new: a counter for (the number of entries - 1)
    !                in each column in the condensed matrix
    !                (of super-columns)
    ! sup_to_sup: super-column index after condensing and that after
    ! getting rid of columns with only <= 1 entries
    integer (kind = myint), dimension (:), allocatable :: &
         col_count_new, sup_to_sup

    ! num_remove: any column of <= num_remove entries are
    !             removed. by default num_remove = 0
    integer (kind = myint) :: num_remove

    num_remove = 0
    if (present(iremove)) num_remove = iremove

    info = 0 ! by default the subroutine is successful.
    if (present(stat)) stat = 0

    n = matrix%n; m = matrix%m

    ! check that the matrix is pattern only!
    if (ZD11_get(matrix%type) /= "pattern" ) then
       info = MC65_ERR_CONDENSE_NOPAT
       return
    end if

    ! if no column no need to work further.
    if (n <= 0) return

    ! all columns assigned as super-column 1 first
    allocate(f2c(n),numsup(n),new(n),flag(n),col_count(n),stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       INFO = MC65_ERR_MEMORY_ALLOC
       return
    end if

    f2c = 1
    numsup(1) = n
    flag(1) = 0
    nsup = 1
    col_count = 0

    ! loop over each row, assign each column entry as a new
    !    super-column
    do i = 1, m
       ! column entries in this row will change their super-column
       ! number, so reduce the number of super-columns of their
       ! old group
       do jj = matrix%ptr(i),matrix%ptr(i+1)-1
          j = matrix%col(jj)
          col_count(j) = col_count(j) + 1
          scol = f2c(j)
          numsup(scol) = numsup(scol) - 1
       end do
       ! assign new groups
       do jj = matrix%ptr(i),matrix%ptr(i+1)-1
          j = matrix%col(jj)
          scol = f2c(j)
         ! has this super-column been seen before in row i?
          if (flag(scol) < i) then
             ! the first appearance of scol in row i
             flag(scol) = i
             if (numsup(scol) > 0) then
                nsup = nsup + 1
                flag(nsup) = i
                f2c(j) = nsup
                numsup(nsup) = 1
                new(scol) = nsup
             else
                ! if zero, the old group disappear anyway, there is
                ! no need to add a new group name, just use
                ! the old one
                numsup(scol) = 1
                new(scol) = scol
             end if
          else
             newid = new(scol)
             f2c(j) = newid
             numsup(newid) = numsup(newid) + 1
          end if
       end do
    end do

    ! now get rid of columns of <= num_remove entries by modifying f2c
    allocate(col_count_new(nsup),sup_to_sup(nsup),stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       INFO = MC65_ERR_MEMORY_ALLOC
       return
    end if

    ! The following array statement can't be used, as multiple columns
    ! (with the _same_ non-zero pattern) can be assigned to each supercolumn
    ! so we need to use a loop instead
    !col_count_new(f2c(1:n)) = col_count(1:n)
    do i = 1, n
       col_count_new(f2c(i)) = col_count(i)
    end do
    nsup_new = 0
    do i = 1, nsup
       if (col_count_new(i) > num_remove) then
          nsup_new = nsup_new + 1
          sup_to_sup(i) = nsup_new
       else
          sup_to_sup(i) = -1
       end if
    end do
    f2c(1:n) = sup_to_sup(f2c(1:n))
    nsup = nsup_new

    deallocate(col_count,col_count_new,sup_to_sup,stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       INFO = MC65_ERR_MEMORY_DEALLOC
       return
    end if

    ! condense the matrix
    flag = 0
    nz = 1
    do i = 1, m
       ! loop over each row
       jjs = matrix%ptr(i)
       matrix%ptr(i) = nz
       do jj = jjs,matrix%ptr(i+1)-1
          j = matrix%col(jj)
          scol = f2c(j)
          ! forget about columns with less than num_remove entries
          if (scol < 0) cycle
          ! if this the first time this group appear?
          ! if so this column will remain
          if (flag(scol) < i) then
             flag(scol) = i
             matrix%col(nz) = scol
             nz = nz + 1
          end if
       end do
    end do
    nz_old = size(matrix%col)
    matrix%ptr(m+1) = nz
    matrix%n = nsup

    nz = nz -1
    if (present(realloc)) then
       if (realloc == 0.0_myreal.or.(realloc > 0.and.&
            nz_old > (1.0_myreal + realloc)*nz)) then
          call csr_matrix_reallocate(matrix,nz,info)
          if (info < 0) return
       end if
    end if


    ! condense the column weight
    deallocate(numsup,new,flag,stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       INFO = MC65_ERR_MEMORY_DEALLOC
       return
    end if

    if (present(col_wgt_real)) then
       allocate(col_wgt_new_r(nsup),stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          INFO = MC65_ERR_MEMORY_ALLOC
          return
       end if

       col_wgt_new_r = 0
       do i = 1, n
          scol = f2c(i)
          ! forget about columns with less than num_remove entries
          if (scol < 0) cycle
          col_wgt_new_r(scol) = col_wgt_new_r(scol) + col_wgt_real(i)
       end do
       col_wgt_real(1:nsup) = col_wgt_new_r
       deallocate(col_wgt_new_r,stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          INFO = MC65_ERR_MEMORY_DEALLOC
       end if
    end if


    if (present(col_wgt_int)) then
       allocate(col_wgt_new_i(nsup),stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          INFO = MC65_ERR_MEMORY_ALLOC
          return
       end if

       col_wgt_new_i = 0
       do i = 1, n
          scol = f2c(i)
          ! forget about columns with less than num_remove entries
          if (scol < 0) cycle
          col_wgt_new_i(scol) = col_wgt_new_i(scol) + col_wgt_int(i)
       end do
       col_wgt_int(1:nsup) = col_wgt_new_i
       deallocate(col_wgt_new_i,stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          INFO = MC65_ERR_MEMORY_DEALLOC
       end if
    end if


    deallocate(f2c,stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       INFO = MC65_ERR_MEMORY_DEALLOC
    end if


  end subroutine csr_matrix_condense

  function csr_matrix_is_pattern(matrix) result (pattern)
    ! check is a matrix is of type "pattern"
    ! MATRIX is of the derived type {\tt ZD11\_TYPE}
    !        with {\tt INTENT (IN)}.
    !        On entry it is the matrix whose type is to
    !        be determined
    type (zd11_type), intent (in) :: matrix
    ! pattern: logical scalar. On exit
    ! is set to true if the matrix is of type "pattern"
    !  and false if not
    logical :: pattern
    if (ZD11_get(matrix%type) == "pattern")  then
       pattern = .true.
    else
       pattern = .false.
    end if
  end function csr_matrix_is_pattern

! ======================== private routines from here =============


  function vacant_unit()  result (file_unit)
    ! utility routine to find a vacant unit and
    ! returns a unit number of a unit that exists and is
    ! not connected

    ! max_unit: maximum number of units to check
    integer (kind = myint) :: max_unit
    parameter (max_unit = 500)
    ! file_unit: the vacant unit
    integer (kind = myint) :: file_unit

    logical :: existed, opened
    integer (kind = myint) :: ios

    do file_unit = 10, max_unit
       inquire (unit = file_unit, exist = existed, &
            opened = opened, iostat = ios)
       if (existed .and. .not. opened .and. ios == 0) return
    end do

    file_unit = -1

  end function vacant_unit




  subroutine expand1(p1,upbound_new1,info,ierr)
    ! reallocate memory for an integer or real allocatable array
    !     of up to 4 dimension.
    ! usage: expand(array,new_bound_1, ...,new_bound_n)
    ! a space of size "new_bound_1 X ... X new_bound_n"
    ! will be allocated and the content of array
    ! will be copied to the beginning of this
    ! memory.
    real (kind = myreal), allocatable:: p2(:),p1(:)
    integer (kind = myint) upbound_new1
    integer (kind = myint) ub(1),lb(1),upb(1)
    integer (kind = myint) info,i,ierr
    ierr=0

    upb(1)=upbound_new1

    lb=lbound(p1)
    ub=ubound(p1)


    allocate(p2(lb(1):upb(1)),stat=ierr)
    if (ierr == 0) then
       p2=0
       do i=lb(1),min(ub(1),upb(1))
          p2(i)=p1(i)
       end do
       deallocate(p1,stat = ierr)
       allocate(p1(lb(1):upb(1)),stat=ierr)
       if (ierr == 0) then
          do i=lb(1),min(ub(1),upb(1))
             p1(i)=p2(i)
          end do
       else
          info = MC65_ERR_MEMORY_DEALLOC
          return
       end if
    else
       info = MC65_ERR_MEMORY_ALLOC
    end if


  end subroutine expand1

  ! ============ for integer (kind = myint) arrays ================

  subroutine iexpand1(p1,upbound_new1,info,ierr)
    integer (kind = myint), allocatable:: p2(:),p1(:)
    integer (kind = myint) upbound_new1
    integer (kind = myint) ub(1),lb(1),upb(1)
    integer (kind = myint) info,i,ierr

    upb(1)=upbound_new1

    lb=lbound(p1)
    ub=ubound(p1)


    allocate(p2(lb(1):upb(1)),stat=ierr)
    if (ierr == 0) then
       p2=0
       do i=lb(1),min(ub(1),upb(1))
          p2(i)=p1(i)
       end do
       deallocate(p1,stat = ierr)
       allocate(p1(lb(1):upb(1)),stat=ierr)
       if (ierr == 0) then
          do i=lb(1),min(ub(1),upb(1))
             p1(i)=p2(i)
          end do
       else
          info = MC65_ERR_MEMORY_DEALLOC
          return
       end if
    else
       info = MC65_ERR_MEMORY_ALLOC
    end if

  end subroutine iexpand1



  function int2str(i,idigit) result (the_str)
    ! converting an integer to a string.
    integer (kind = myint):: i,idigit

    character (len = idigit) :: the_str

    character (len = 20) :: fmt_str

    fmt_str = "(I                 )"

    write(fmt_str(3:19),"(I6)") idigit

    write(the_str,fmt_str) i

  end function int2str

  function digit(i)
    ! return the number of digit of an integer
    integer (kind = myint) :: i,digit

    integer (kind = myint) :: absi


    absi = abs(i)

    digit = 1

    do while (absi >= 10)

       absi = absi / 10

       digit = digit + 1

    end do

    if (i < 0) digit = digit + 1

  end function digit





!!$  subroutine csr_matrix_reset_type(matrix,type,info)
!!$    ! subroutine csr_matrix_reset_type(matrix,type,info):
!!$    ! reset the type of the matrix
!!$
!!$    ! matrix: is of the derived type ZD11_type, INTENT (INOUT).
!!$    !         This is the matrix whose type is to be changed.
!!$    type (zd11_type), intent (inout) :: matrix
!!$
!!$    ! type: is a character array of unspecified length with
!!$    !       INTENT (IN).
!!$    !       It hold the new matrix type to be set.
!!$    character (len=*), intent (in) :: type
!!$
!!$    ! info: integer scaler of INTENT (OUT).
!!$    !       = 0 if the subroutine returned successfully
!!$    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
!!$    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
!!$    integer (kind = myint), intent (out) :: info
!!$
!!$    ! =============== local variables ===============
!!$
!!$    integer (kind =  myint) :: ierr
!!$
!!$    info = 0
!!$
!!$    if (allocated(matrix%type)) deallocate(matrix%type,stat = ierr)
!!$    if (ierr /= 0) then
!!$       info = MC65_ERR_MEMORY_DEALLOC
!!$       return
!!$    end if
!!$    call zd11_put(matrix%type,type,stat = ierr)
!!$    if (ierr /= 0) then
!!$       info = MC65_ERR_MEMORY_ALLOC
!!$    end if
!!$
!!$  end subroutine csr_matrix_reset_type
!!$
!!$
!!$  subroutine csr_matrix_fill(matrix,i,row_counter,fill_entry,&
!!$       fill_value)
!!$    ! fill an entry -- fill the current
!!$    ! entry with column index fill_entry and column value
!!$    ! fill_value, in row i of the matrix.
!!$    ! row_counter is a counter given
!!$    ! for the vacant position in the rows. This subroutine
!!$    ! is used when the matrix has been allocated storage, and also
!!$    ! the row pointer has already been established, and
!!$    ! one is at the stage of filling in the entry values
!!$    ! and column indices.
!!$
!!$    ! matrix: the matrix whose entries are to be filled.
!!$    type (zd11_type), intent (inout) :: matrix
!!$
!!$    ! i: the entry is to be filled at this i-th row.
!!$    integer (kind = myint), intent (in) :: i
!!$
!!$    ! row_counter: integer array, row_counter(i)
!!$    !              is number of entries filled in row i so far.
!!$    integer (kind = myint), dimension (*), intent (inout) :: &
!!$         row_counter
!!$
!!$    ! fill_entry: column index of the current entry
!!$    integer (kind = myint), intent (in) :: fill_entry
!!$    real (kind = myreal), intent (in), optional :: fill_value
!!$
!!$    integer (kind = myint), pointer, dimension (:) :: ja
!!$    real (kind = myreal), pointer, dimension (:) :: aa
!!$
!!$
!!$    ja => csr_matrix_getrow(matrix,i)
!!$    row_counter(i) = row_counter(i) + 1
!!$    ja(row_counter(i)) = fill_entry
!!$    if (present(fill_value)) then
!!$       aa => csr_matrix_getrowval(matrix,i)
!!$       aa(row_counter(i)) = fill_value
!!$    end if
!!$  end subroutine csr_matrix_fill
!!$
!!$
!!$  subroutine csr_matrix_component(matrix,root,mask,comp_nvtx,&
!!$       comp_nz,comp_vtx_list)
!!$    !
!!$    ! matrix_component: this subroutine finds the component
!!$    ! of the graph starting from the root. This component
!!$    ! is stored in the list of vertices "comp_vtx_list",
!!$    ! with "comp_nvtx" vertices and "comp_nz"/2 edges
!!$    ! (thus comp_nz nonzeros in the matrix describing
!!$    ! the component).
!!$
!!$    ! =================== arguments =========================
!!$
!!$    ! matrix: the graph to be inspected. This must be a symmetric
!!$    ! matrix with no diagonal
!!$    type (zd11_type), intent (in) :: matrix
!!$
!!$    ! root: the root of the current component
!!$    integer (kind = myint), intent (in) :: root
!!$
!!$    ! comp_nvtx: number of vertices in this component
!!$    ! comp_nz: number of edges*2 in the component
!!$    integer (kind = myint), intent (out) :: comp_nvtx,comp_nz
!!$
!!$    ! mask: has the dimension as the number of vertices.
!!$    ! It is zero if the vertex is not yet visited,
!!$    ! nonzero if visited. When a vertex is visited the first time,
!!$    ! its mask is set to the new index of this vertex in the component
!!$    ! comp_vtx_list: a list of the vertices in the current component.
!!$    ! mask and comp_vtx_list satisfies mask(comp_vtx(i)) = i,
!!$    ! i = 1,...,comp_nvtx
!!$    integer (kind = myint), dimension (:), intent (inout) :: mask, &
!!$         comp_vtx_list
!!$
!!$    ! ===================== local variables =====================
!!$    ! front_list: this array hold the list of vertices that
!!$    !             form the front in the breadth first search
!!$    !             starting from the root
!!$    integer (kind = myint), dimension (:), allocatable :: front_list
!!$
!!$    ! front_sta: the starting position in the front_list of the
!!$    !            current front.
!!$    ! front_sto: the ending position in the front_list of the current
!!$    ! front.
!!$    integer (kind = myint) :: front_sta, front_sto
!!$
!!$    ! n: number of vertices in the graph
!!$    ! ierr: error tag
!!$    integer (kind = myint) :: n,ierr
!!$    ! v: a vertex
!!$    ! u: neighbor of v
!!$    ! i: loop index
!!$    ! j: loop index
!!$    integer (kind = myint) :: i,j,u,v
!!$
!!$    n = matrix%m
!!$
!!$    ! the first front
!!$    allocate(front_list(n),stat = ierr)
!!$    if (ierr /= 0) stop "error allocating in matrix_component"
!!$    front_sta = 1
!!$    front_sto = 1
!!$    front_list(1) = root
!!$
!!$    comp_nvtx = 1
!!$    comp_nz = 0
!!$    comp_vtx_list(1) = root ! one vertex in the component so far
!!$    mask(root) = comp_nvtx ! mask the root
!!$
!!$    do while (front_sto-front_sta >= 0)
!!$       do i = front_sta, front_sto
!!$          v = front_list(i) ! pick a vertex from the front
!!$          ! count link to all neighbors as edges
!!$          comp_nz = comp_nz + matrix%ptr(v+1)-matrix%ptr(v)
!!$          do j = matrix%ptr(v),matrix%ptr(v+1)-1
!!$             u = matrix%col(j) ! pick its neighbor
!!$             if (mask(u) /= 0) cycle
!!$             comp_nvtx = comp_nvtx + 1 ! found a unmasked vertex
!!$             mask(u) = comp_nvtx ! mask this vertex
!!$             ! add this vertex to the component
!!$             comp_vtx_list(comp_nvtx) = u
!!$             front_list(comp_nvtx) = u ! also add it to the front
!!$          end do
!!$       end do
!!$       front_sta = front_sto + 1
!!$       front_sto = comp_nvtx
!!$    end do
!!$    deallocate(front_list,stat = ierr)
!!$    if (ierr /= 0) stop "error deallocating in matrix_component"
!!$
!!$
!!$  end subroutine csr_matrix_component
!!$
!!$  subroutine csr_matrix_crop(matrix,comp_nvtx,comp_vtx_list,&
!!$       mask,submatrix,comp_nz,info)
!!$    ! matrix_crop: generate a submatrix describing a component of
!!$    ! the graph (matrix).
!!$
!!$    ! =================== arguments =========================
!!$    ! matrix: the graph to be inspected. This must be a symmetric
!!$    ! matrix with no diagonal
!!$    type (zd11_type), intent (in) :: matrix
!!$
!!$    ! comp_nvtx: number of vertices in this component
!!$    ! comp_nz: number of edges*2 in the component (optional)
!!$    integer (kind = myint), intent (in) :: comp_nvtx
!!$    integer (kind = myint), intent (in), optional :: comp_nz
!!$    integer (kind = myint), intent (out) :: info
!!$
!!$    ! mask: has the dimension as the number of vertices.
!!$    ! The mask of a vertex is set to the new index
!!$    ! of this vertex in the component it belongs
!!$    ! comp_vtx_list: a list of the vertices in the current component
!!$    integer (kind = myint), dimension (:), intent (in) :: mask, &
!!$         comp_vtx_list
!!$
!!$    ! submatrix: the submatrix describing the component
!!$    type (zd11_type), intent (inout) :: submatrix
!!$
!!$    ! ===================== local variables =====================
!!$    ! i: loop index
!!$    ! j: loop index
!!$    ! nz: number of nonzeros in the submatrix
!!$    ! n: number of rows in the submatrix
!!$    ! m: number of columns in the submatrix
!!$    ! v: a vertex in its original index
!!$    ! l1: starting index
!!$    ! l2: stopping index
!!$    ! sl1: starting index
!!$    ! sl2: stopping index
!!$    integer (kind = myint) :: i,j,nz,n,m,v,l1,l2,sl1,sl2
!!$
!!$    ! ia: row pointer of the original graph
!!$    ! ja: column indices of the original graph
!!$    ! sia: row pointer of the subgraph
!!$    ! sja: column indices of the subgraph
!!$    integer (kind = myint), dimension (:), pointer :: ia,ja,sia,sja
!!$
!!$    ! a: entry values of the original graph
!!$    ! sa: entry values of the subgraph
!!$    real (kind = myreal), dimension (:), pointer :: a,sa
!!$
!!$    info = 0
!!$    ia => matrix%ptr
!!$    ja => matrix%col
!!$    if (.not.present(comp_nz)) then
!!$       nz = 0
!!$       do i = 1, comp_nvtx
!!$          v = comp_vtx_list(i)
!!$          nz = nz + ia(v+1) - ia(v)
!!$       end do
!!$    else
!!$       nz = comp_nz
!!$    end if
!!$
!!$    n = comp_nvtx; m = n;
!!$    call csr_matrix_construct(submatrix,n,nz,info,n = m,&
!!$         type = zd11_get(matrix%type))
!!$    sia => submatrix%ptr
!!$    sja => submatrix%col
!!$    if (ZD11_get(matrix%type) /= "pattern" ) then
!!$       sa => submatrix%val
!!$       a => matrix%val
!!$    end if
!!$    sia(1) = 1
!!$    do i = 1, comp_nvtx
!!$       v = comp_vtx_list(i)
!!$       l1 = ia(v); l2 = ia(v+1)-1
!!$       sia(i+1) = sia(i) + l2 - l1 + 1
!!$       sl1 = sia(i); sl2 = sia(i+1)-1
!!$       ! convert original index to new index in the component
!!$       sja(sl1:sl2) = mask(ja(l1:l2))
!!$       if (ZD11_get(matrix%type) /= "pattern" ) then
!!$          sa(sl1:sl2) = a(l1:l2)
!!$       end if
!!$    end do
!!$  end subroutine csr_matrix_crop
!!$
!!$
!!$
!!$
!!$  subroutine csr_matrix_crop_unsym(matrix,nrow,row_list,submatrix,info)
!!$    ! matrix_crop: generate a submatrix of the whole matrix
!!$    ! by taking all the nrow rows in the row_list.
!!$    ! column size of submatrix is assigned as the same as that
!!$    ! for matrix.
!!$
!!$    ! =================== arguments =========================
!!$    ! matrix: the matrix to be inspected.
!!$    type (zd11_type), intent (in) :: matrix
!!$
!!$    ! nrow: number of vertices in this component
!!$    integer (kind = myint), intent (in) :: nrow
!!$
!!$    ! row_list: a list of the vertices in the current component
!!$    integer (kind = myint), dimension (:), intent (in) :: row_list
!!$
!!$    ! submatrix: the submatrix describing the component
!!$    type (zd11_type), intent (inout) :: submatrix
!!$    integer (kind = myint), intent (out) :: info
!!$    ! ===================== local variables =====================
!!$    ! i: loop index
!!$    ! j: loop index
!!$    ! nz: number of nonzeros in the submatrix
!!$    ! n: number of rows in the submatrix
!!$    ! m: number of columns in the submatrix
!!$    ! v: a vertex in its original index
!!$    ! l1: starting index
!!$    ! l2: stopping index
!!$    ! sl1: starting index
!!$    ! sl2: stopping index
!!$    integer (kind = myint) :: i,j,nz,n,m,v,l1,l2,sl1,sl2
!!$
!!$    ! ia: row pointer of the original matrix
!!$    ! ja: column indices of the original matrix
!!$    ! sia: row pointer of the submatrix
!!$    ! sja: column indices of the submatrix
!!$    integer (kind = myint), dimension (:), pointer :: ia,ja,sia,sja
!!$
!!$    ! a: entry values of the original matrix
!!$    ! sa: entry values of the submatrix
!!$    real (kind = myreal), dimension (:), pointer :: a,sa
!!$
!!$    info = 0
!!$    ia => matrix%ptr
!!$    ja => matrix%col
!!$    nz = 0
!!$    do i = 1, nrow
!!$       v = row_list(i)
!!$       nz = nz + ia(v+1) - ia(v)
!!$    end do
!!$
!!$    n = nrow; m = matrix%n;
!!$    call csr_matrix_construct(submatrix,n,nz,info,n = m,&
!!$         type = zd11_get(matrix%type))
!!$
!!$    sia => submatrix%ptr
!!$    sja => submatrix%col
!!$    if (ZD11_get(matrix%type) /= "pattern" ) then
!!$       sa => submatrix%val
!!$       a => matrix%val
!!$    end if
!!$    sia(1) = 1
!!$    do i = 1, nrow
!!$       v = row_list(i)
!!$       l1 = ia(v); l2 = ia(v+1)-1
!!$       sia(i+1) = sia(i) + l2 - l1 + 1
!!$       sl1 = sia(i); sl2 = sia(i+1)-1
!!$       ! convert original index to new index in the component
!!$       sja(sl1:sl2) = ja(l1:l2)
!!$       if (ZD11_get(matrix%type) /= "pattern" ) then
!!$          sa(sl1:sl2) = a(l1:l2)
!!$       end if
!!$    end do
!!$  end subroutine csr_matrix_crop_unsym
!!$
!!$

end module HSL_MC65_single




