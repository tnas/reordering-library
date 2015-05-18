* COPYRIGHT (c) 1987 AEA Technology
* Original date 10 Feb 1993
C       Toolpack tool decs employed.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC34AD(N,IRN,JCOLST,YESA,A,IW)
C THIS SUBROUTINE ACCEPTS AS INPUT THE STANDARD DATA STRUCTURE FOR
C     A SYMMETRIC MATRIX STORED AS A LOWER TRIANGLE AND PRODUCES
C     AS OUTPUT THE SYMMETRIC MATRIX HELD IN THE SAME DATA
C     STRUCTURE AS A GENERAL MATRIX.
C N IS AN INTEGER VARIABLE THAT MUST BE SET BY THE USER TO THE
C     ORDER OF THE MATRIX. NOT ALTERED BY THE ROUTINE
C     RESTRICTION (IBM VERSION ONLY): N LE 32767.
C IRN IS AN INTEGER (INTEGER*2 IN IBM VERSION) ARRAY THAT
C     MUST BE SET BY THE USER TO HOLD THE ROW INDICES OF THE LOWER
C     TRIANGULAR PART OF THE SYMMETRIC MATRIX.  THE ENTRIES OF A
C     SINGLE COLUMN ARE CONTIGUOUS. THE ENTRIES OF COLUMN J
C     PRECEDE THOSE OF COLUMN J+1 (J_=_1, ..., N-1), AND THERE IS
C     NO WASTED SPACE BETWEEN COLUMNS. ROW INDICES WITHIN A COLUMN
C     MAY BE IN ANY ORDER.  ON EXIT IT WILL HAVE THE SAME MEANING
C     BUT WILL BE CHANGED TO HOLD THE ROW INDICES OF ENTRIES IN
C     THE EXPANDED STRUCTURE.  DIAGONAL ENTRIES NEED NOT BE
C     PRESENT. THE NEW ROW INDICES ADDED IN THE UPPER TRIANGULAR
C     PART WILL BE IN ORDER FOR EACH COLUMN AND WILL PRECEDE THE
C     ROW INDICES FOR THE LOWER TRIANGULAR PART WHICH WILL REMAIN
C     IN THE INPUT ORDER.
C JCOLST IS AN INTEGER ARRAY OF LENGTH N+1 THAT MUST BE SET BY
C     THE USER SO THAT JCOLST(J) IS THE POSITION IN ARRAYS IRN AND
C     A OF THE FIRST ENTRY IN COLUMN J (J_=_1, ..., N).
C     JCOLST(N+1) MUST BE SET TO ONE MORE THAN THE TOTAL NUMBER OF
C     ENTRIES.  ON EXIT, JCOLST(J) WILL HAVE THE SAME MEANING BUT
C     WILL BE CHANGED TO POINT TO THE POSITION OF THE FIRST ENTRY
C     OF COLUMN J IN THE EXPANDED STRUCTURE. THE NEW VALUE OF
C     JCOLST(N+1) WILL BE ONE GREATER THAN THE NUMBER OF ENTRIES
C     IN THE EXPANDED STRUCTURE.
C YESA IS A LOGICAL VARIABLE THAT MUST BE SET TO .TRUE. IF THE
C     USER DESIRES TO GENERATE THE EXPANDED FORM FOR THE VALUES ALSO.
C     IF YESA IS .FALSE., THE ARRAY A WILL NOT BE REFERENCED.  IT IS
C     NOT ALTERED BY THE ROUTINE.
C A IS A REAL (DOUBLE PRECISION IN THE D VERSION) ARRAY THAT
C     CAN BE SET BY THE USER SO THAT A(K) HOLDS THE VALUE OF THE
C     ENTRY IN POSITION K OF IRN, {K = 1, _..._ JCOLST(N+1)-1}.
C     ON EXIT, IF YESA IS .TRUE., THE ARRAY WILL HOLD THE VALUES
C     OF THE ENTRIES IN THE EXPANDED STRUCTURE CORRESPONDING TO
C     THE OUTPUT VALUES OF IRN.   IF YESA IS .FALSE., THE ARRAY IS
C     NOT ACCESSED BY THE SUBROUTINE.
C IW IS AN INTEGER (INTEGER*2 IN IBM VERSION) ARRAY OF LENGTH
C     N THAT WILL BE USED AS WORKSPACE.
C
C CKP1 IS A LOCAL VARIABLE USED AS A RUNNING POINTER.
C OLDTAU IS NUMBER OF ENTRIES IN SYMMETRIC STORAGE.
C     .. Scalar Arguments ..
      INTEGER N
      LOGICAL YESA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*)
      INTEGER IRN(*),IW(*),JCOLST(*)
C     ..
C     .. Local Scalars ..
      INTEGER CKP1,I,I1,I2,II,IPKP1,IPOS,J,JSTART,LENK,NDIAG,NEWTAU,
     +        OLDTAU
C     ..
C     .. Executable Statements ..
C
      OLDTAU = JCOLST(N+1) - 1
C INITIALIZE WORK ARRAY
      DO 5 I = 1,N
        IW(I) = 0
    5 CONTINUE
C
C IW(J) IS SET EQUAL TO THE TOTAL NUMBER OF ENTRIES IN COLUMN J
C     OF THE EXPANDED SYMMETRIC MATRIX.
C NDIAG COUNTS NUMBER OF DIAGONAL ENTRIES PRESENT
      NDIAG = 0
      DO 20 J = 1,N
        I1 = JCOLST(J)
        I2 = JCOLST(J+1) - 1
        IW(J) = IW(J) + I2 - I1 + 1
        DO 10 II = I1,I2
          I = IRN(II)
          IF (I.NE.J) THEN
            IW(I) = IW(I) + 1

          ELSE
            NDIAG = NDIAG + 1
          END IF

   10   CONTINUE
   20 CONTINUE
C
C NEWTAU IS NUMBER OF ENTRIES IN EXPANDED STORAGE.
      NEWTAU = 2*OLDTAU - NDIAG
C IPKP1 POINTS TO POSITION AFTER END OF COLUMN BEING CURRENTLY
C     PROCESSED
      IPKP1 = OLDTAU + 1
C CKP1 POINTS TO POSITION AFTER END OF SAME COLUMN IN EXPANDED
C     STRUCTURE
      CKP1 = NEWTAU + 1
C GO THROUGH THE ARRAY IN THE REVERSE ORDER PLACING LOWER TRIANGULAR
C     ELEMENTS IN THE APPROPRIATE SLOTS.
      DO 40 J = N,1,-1
        I1 = JCOLST(J)
        I2 = IPKP1
C LENK IS NUMBER OF ENTRIES IN COLUMN J OF ORIGINAL STRUCTURE
        LENK = I2 - I1
C JSTART IS RUNNING POINTER TO POSITION IN NEW STRUCTURE
        JSTART = CKP1
C SET IKP1 FOR NEXT COLUMN
        IPKP1 = I1
        I2 = I2 - 1
C RUN THROUGH COLUMNS IN REVERSE ORDER
C LOWER TRIANGULAR PART OF COLUMN MOVED TO END OF SAME COLUMN IN
C     EXPANDED FORM
        DO 30 II = I2,I1,-1
          JSTART = JSTART - 1
          IF (YESA) A(JSTART) = A(II)
          IRN(JSTART) = IRN(II)
   30   CONTINUE
C JCOLST IS SET TO POSITION OF FIRST ENTRY IN LOWER TRIANGULAR PART OF
C     COLUMN J IN EXPANDED FORM
        JCOLST(J) = JSTART
C SET CKP1 FOR NEXT COLUMN
        CKP1 = CKP1 - IW(J)
C RESET IW(J) TO NUMBER OF ENTRIES IN LOWER TRIANGLE OF COLUMN.
        IW(J) = LENK
   40 CONTINUE
C
C AGAIN SWEEP THROUGH THE COLUMNS IN THE REVERSE ORDER, THIS
C     TIME WHEN ONE IS HANDLING COLUMN J THE UPPER TRIANGULAR
C     ELEMENTS A(J,I) ARE PUT IN POSITION.
      DO 80 J = N,1,-1
        I1 = JCOLST(J)
        I2 = JCOLST(J) + IW(J) - 1
C RUN DOWN COLUMN IN ORDER
C NOTE THAT I IS ALWAYS GREATER THAN OR EQUAL TO J
        DO 60 II = I1,I2
          I = IRN(II)
          IF (I.EQ.J) GO TO 60
          JCOLST(I) = JCOLST(I) - 1
          IPOS = JCOLST(I)
          IF (YESA) A(IPOS) = A(II)
          IRN(IPOS) = J
   60   CONTINUE
   80 CONTINUE
      JCOLST(N+1) = NEWTAU + 1
      RETURN

      END
C COPYRIGHT (c) 1998 Council for the Central Laboratory
C                    of the Research Councils
C
C Version 1.1.0. For history see ChangeLog

      SUBROUTINE MC60AD(N,LIRN,IRN,ICPTR,ICNTL,IW,INFO)

C  Accepts the pattern of the lower-triangular part of the matrix
C  and constructs the pattern of the whole matrix.

      INTEGER N
C N has intent IN and holds the number of variables. N>=1.
      INTEGER LIRN
C LIRN has intent IN and holds the length of the array IRN.
      INTEGER IRN(LIRN)
C IRN has intent INOUT. On entry, the leading part must hold the
C     row indices of the entries in the lower triangle (including
C     those on the diagonal) of the matrix.
C     On return, holds the row indices of the whole matrix in the
C     same format.
      INTEGER ICPTR(N+1)
C ICPTR has intent INOUT. ICPTR(J) holds the position in
C     the array IRN of the first entry in column J (J=1,...,N) and
C     the last entry is in position ICPTR(N+1)-1.
      INTEGER ICNTL(2)
C ICNTL has intent IN and controls the actions.
C   ICNTL(1)=0: computation terminates if duplicated or out-of-range
C               indices are detected.
C   ICNTL(1)=1: any duplicated or out-of-range indices are ignored.
C   ICNTL(2) holds the unit number for printing or 0 if
C     no printing is required.
      INTEGER IW(N)
C IW is used for workspace.
      INTEGER INFO(4)
C INFO has intent OUT.
C  INFO(1) is used as an error flag.
C     Normal:
C     0  No out-of-range or duplicated indices.
C     1  Some out-of-range or duplicated indices.
C     Errors:
C     -1 N<1 or LIRN < ICPTR(N+1)-1.  Immediate return with
C        IRN and ICPTR unchanged.
C     -2 LIRN is too small. INFO(4) is set to the minimum value
C        that will suffice. If ICNTL(0)=1, any out-of-range or
C        duplicated variable indices will have been excluded from
C        IRN and ICPTR.  Whatever the value of ICNTL(0),
C        any diagonal entries will have been moved to the front of
C        their columns. Otherwise, IRN and ICPTR are unchanged.
C     -3 ICNTL = 0 and one or more variable indices either lies
C        outside the lower triangle  of the matrix or is duplicated
C        (see INFO(2) and INFO(3)). IRN and ICPTR are unchanged.
C   INFO(2) holds the number of variable indices in IRN found to be
C        out-of-range.
C   INFO(3) holds the number of duplicate variable indices found
C        in IRN.
C   INFO(4) holds the minimum value that will suffice for LIRN,
C        unless INFO(1)=-1.

C     .. Local Scalars ..
      INTEGER CKP1,I,I1,I2,II,IOUT,IPOS,IREP,J,KZ,LP,NDIAG,NEWTAU
C CKP1   Position after end of current column in new structure.
C I      Row index
C I1     Starting do index
C I2     Ending do index
C II     Do index
C IOUT   Number of entries outside the lower triangle.
C IPOS   Position of entry in new structure.
C IREP   Number of repeated entries.
C J      Column index
C KZ     Number of genuine entries.
C LP     Unit for printing
C NDIAG  Number of diagonal entries.
C NEWTAU Number of entries in expanded storage.

      LP = ICNTL(2)

C Initialise INFO
      DO 5 J = 1,4
         INFO(J) = 0
    5 CONTINUE

C  Check N.
      IF (N.LT.1) THEN
          INFO(1) = -1
          IF (LP.GT.0) WRITE (LP,'(/,A,I3/,A,I6)')
     *         ' MC60A/AD error: INFO(1) =', INFO(1), ' N =',N
          RETURN
      END IF

C  Check LIRN.
      IF (LIRN.LT.ICPTR(N+1)-1) THEN
          INFO(1) = -1
          IF (LP.GT.0) WRITE (LP,'(/,A,I3/,A)')
     *         ' MC60A/AD error:  INFO(1) =', INFO(1),
     *         ' LIRN is less than ICPTR(N+1)-1'
          RETURN
      END IF

C    Look for any repeated entries and entries with out-of-range
C    indices.
C    Remove such entries and issue a warning if ICNTL(1) = 1.
C IW(I) will be set to J if row I is encountered in column J.
      DO 10 I = 1,N
        IW(I) = 0
   10 CONTINUE
      IOUT = 0
      IREP = 0

      KZ = 0
      IF (ICNTL(1).EQ.1) THEN
         I1 = ICPTR(1)
         ICPTR(1) = 1
         DO 16 J = 1,N
            DO 15 II = I1,ICPTR(J+1) - 1
               I = IRN(II)
               IF (I.GT.N .OR. I.LT.J) THEN
                  IOUT = IOUT + 1
               ELSE IF (IW(I).EQ.J) THEN
                  IREP = IREP + 1
               ELSE
                  KZ = KZ + 1
                  IRN(KZ) = I
                  IW(I) = J
               END IF
   15       CONTINUE
            I1 = ICPTR(J+1)
            ICPTR(J+1) = KZ + 1
   16    CONTINUE

         IF (IOUT.GT.0) THEN
             INFO(1) = 1
             IF (LP.GT.0) WRITE (LP,'(/,A,I6,A)')
     *         ' MC60A/AD warning:',IOUT,' out-of-range entries ignored'
         END IF
         IF (IREP.GT.0) THEN
             INFO(1) = 1
             IF (LP.GT.0) WRITE (LP,'(/,A,I6,A)')
     *         ' MC60A/AD warning:',IREP,' duplicated entries ignored'
         END IF

         INFO(2) = IOUT
         INFO(3) = IREP

      ELSE

         I1 = ICPTR(1)
         DO 26 J = 1,N
            DO 25 II = I1,ICPTR(J+1) - 1
               I = IRN(II)
               IF (I.GT.N .OR. I.LT.J) THEN
                  IOUT = IOUT + 1
               ELSE IF (IW(I).EQ.J) THEN
                  IREP = IREP + 1
               ELSE
                  KZ = KZ + 1
                  IW(I) = J
               END IF
   25       CONTINUE
            I1 = ICPTR(J+1)
   26    CONTINUE

         IF (IOUT.GT.0 .OR. IREP.GT.0) THEN
            INFO(1) = -3
            IF (LP.GT.0) THEN
              WRITE (LP,'(/,A,I3)')
     *            ' MC60A/AD error:  INFO(1) =', INFO(1)
            IF (IOUT.GT.0) WRITE (LP,'(A,I6,A)')
     *         ' There are ',IOUT,' out-of-range entries'
            IF (IREP.GT.0) WRITE (LP,'(A,I6,A)')
     *         ' There are ',IREP,' duplicated entries'
            END IF
            INFO(2) = IOUT
            INFO(3) = IREP
            RETURN
         END IF

      END IF

C IW(J) is set equal to the total number of entries in column J
C     of the expanded symmetric matrix.
      DO 30 J = 1,N
        IW(J) = 0
   30 CONTINUE
      NDIAG = 0
      DO 40 J = 1,N
        I1 = ICPTR(J)
        I2 = ICPTR(J+1) - 1
        IW(J) = IW(J) + I2 - I1 + 1
        DO 35 II = I1,I2
          I = IRN(II)
          IF (I.NE.J) THEN
            IW(I) = IW(I) + 1
          ELSE
            NDIAG = NDIAG + 1
          END IF
   35   CONTINUE
   40 CONTINUE

      NEWTAU = 2*KZ - NDIAG
      INFO(4) = NEWTAU
      IF (NEWTAU.GT.LIRN) THEN
          INFO(1) = -2
          IF (LP.GT.0) WRITE (LP,'(/,A)')
     *         ' MC60A/AD error: LIRN is too small'
          RETURN
      END IF

      I1 = KZ + 1
C CKP1 points to position after end of same column in expanded
C     structure.
      CKP1 = NEWTAU + 1
C Go through the array in the reverse order placing lower triangular
C     elements in the appropriate slots.
      DO 60 J = N,1,-1
        I2 = I1 - 1
        I1 = ICPTR(J)
C IPOS is running pointer to position in new structure.
        IPOS = CKP1
C Run through columns in reverse order.
C Lower triangular part of column moved to end of same column in
C expanded form.
          DO 50 II = I2,I1,-1
            IPOS = IPOS - 1
            IRN(IPOS) = IRN(II)
   50     CONTINUE
C ICPTR is set to position of first entry in lower triangular part of
C     column J in expanded form.
        ICPTR(J) = IPOS
C Set CKP1 for next column
        CKP1 = CKP1 - IW(J)
C Reset IW(J) to number of entries in lower triangle of column.
        IW(J) = I2 - I1 + 1
   60 CONTINUE
C
C Again sweep through the columns in the reverse order.
C     This time when one is handling column J, the upper triangular
C     elements A(J,I) are put in position.
      DO 80 J = N,1,-1
        I1 = ICPTR(J)
        I2 = ICPTR(J) + IW(J) - 1
        IF(I1.LE.I2) THEN
C Run down column in order.
          DO 70 II = I1,I2
            I = IRN(II)
C Skip the diagonal entry, if present
            IF(I.EQ.J) GO TO 70
C Note that I is always greater than J, so ICPTR(I) will not be
C    needed in later cycles of the loop.
            ICPTR(I) = ICPTR(I) - 1
            IRN(ICPTR(I)) = J
   70     CONTINUE
        END IF
   80 CONTINUE
      ICPTR(N+1) = NEWTAU + 1

      END


      SUBROUTINE MC60BD(N,LIRN,IRN,ICPTR,NSUP,SVAR,VARS,IW)

C Given a symmetric sparsity pattern, this procedure finds
C supervariables and replaces the pattern by its compressed
C equivalent. The user must supply the pattern of the entries of
C the whole matrix, including the diagonal.

C No check is made on the validity of the data. There must be no
C duplicate entries or out-of-range indices. Columns with no
C entries are permitted.

      INTEGER N
C N has intent IN and holds the matrix order.
      INTEGER LIRN
C LIRN has intent IN and holds the length of the array IRN.
      INTEGER IRN(LIRN)
C IRN has intent INOUT. On entry it holds the row indices of
C     the entries in the matrix. On exit it holds the row indices of
C     the entries in the condensed matrix.
      INTEGER ICPTR(N+1)
C ICPTR has intent INOUT. On entry, ICPTR(J) holds the position in
C     the array IRN of the first entry in column J (J=1,...,N) and
C     the last entry is in position ICPTR(N+1)-1. On return, it holds
C     the corresponding data for the compressed matrix.
      INTEGER NSUP
C NSUP has intent OUT. On return, it holds the number of
C     supervariables.
      INTEGER SVAR(N)
C SVAR has intent OUT. On return, SVAR(I) is the supervariable
C     to which variable I belongs, I = 1, ..., N.
      INTEGER VARS(N)
C VARS has intent OUT. On return, VARS(IS) holds the number of
C     variables in supervariable IS, IS = 1, ..., NSUP.
      INTEGER IW(2*N)
C IW is used for workspace.

C     .. Local Scalars ..
      INTEGER FLAG

C     .. External Subroutines ..
      EXTERNAL MC60OD,MC60PD

C Look for supervariables
      FLAG = N+1
      CALL MC60OD(N,N,LIRN,IRN,ICPTR,SVAR,NSUP,IW,VARS,IW(FLAG))

C Adjust the permutation and compress the pattern
      CALL MC60PD(N,NSUP,LIRN,IRN,ICPTR,SVAR,VARS,IW,IW(FLAG))

      END

      SUBROUTINE MC60CD(N,NSUP,LIRN,IRN,ICPTR,VARS,JCNTL,
     +                  PERMSV,WEIGHT,PAIR,INFO,IW,W)
C
C    Permute a condensed matrix for small wavefront and profile
C
      INTEGER N,NSUP,LIRN,IRN(LIRN),ICPTR(NSUP+1),VARS(NSUP),JCNTL(2)
      INTEGER PERMSV(NSUP),PAIR(2,*),INFO(4),IW(3*NSUP+1)
      DOUBLE PRECISION WEIGHT(2),W(NSUP)
C N has intent IN and holds the order of the matrix.
C NSUP has intent IN. It holds the number of supervariables.
C LIRN has intent IN and holds the length of the array IRN.
C IRN has intent IN. It must hold the row indices of the
C     supervariable representation of the matrix, by columns.
C ICPTR has intent IN. It must be set so that ICPTR(J) holds the
C     position in the array IRN of the first entry in column J
C     (J=1,...,NSUP) and the last entry is in position ICPTR(NSUP+1)-1.
C VARS has intent IN. It must hold the number of variables in
C     supervariable IS, IS = 1, ..., NSUP.
C JCNTL has intent IN.
C   JCNTL(1) controls the choice of algorithm:
C      0 - Sloan's global priority.
C      1 - Reverse Cuthill-McKee ordering.
C   JCNTL(2) controls the algorithmic details:
C      0 - Automatic choices
C      1 - Pseudoperipheral pairs specified in PAIR.
C      2 - Global priority supplied in PERMSV.
C PERMSV has intent INOUT. If JCNTL(2)=2, it must be set to the global
C      priorities. During execution, negative values indicate ordered
C      variables and hold negation of new position.
C      On exit, the new index for supervariable I
C      is given by PERMSV(I), I = 1,...,NSUP.
C WEIGHT has intent IN. It holds the weights for the priorities.
C PAIR has intent INOUT.
C     If JCNTL(2)=0, PAIR need not set on entry and on return
C        PAIR(1,IC), PAIR(2,IC) hold the pseudoperipheral pair for
C        nontrivial component IC, IC = 1,2,...,INFO(1).
C        The first component is the largest.
C     If JCNTL(2)=1, PAIR must be set on entry to the
C        pseudoperipheral pairs of the components and is not
C        altered; the first component need not be the largest.
C     If JCNTL(2)=2, PAIR is not used.
C INFO has intent OUT.
C     INFO(1) - number of nontrivial components (2 or more nodes)
C        in the graph of the matrix
C     INFO(2) - number of variables in the largest component.
C     INFO(3) - number of level sets in the level-set structure of the
C        largest component
C     INFO(4) - width of the level-set structure of the largest
C        component
C IW is used for workspace.
C   IW(1),... holds lists on nodes at level 1, level 2, ...
C   IW(XLS+1),... holds the starting positions of these lists
C   IW(LIST+1),... holds a list of nodes.
C W is used for workspace.

      INTEGER DEGREE,I,IL,HINFO(6),J,K,LIST,LSTNUM,LWIDTH,LWDTH1,
     *        LZNUM,MAXPSV,MINPSV,NLVL,NLVL1,NODES,NSTOP,NSTRT,NVARS,XLS
      DOUBLE PRECISION NWGHT(2)
C DEGREE Degree, counting each supervariable as 1.
C I      Do index.
C IL     Do index.
C HINFO   Infomation from MC60H/HD
C J      Do index.
C K      Column start in IRN.
C LIST   See comment for IW.
C LSTNUM Last nonzero column ordered
C LWIDTH Width of level-set structure
C LWDTH1 Width of level-set structure from first end
C LZNUM  Last zero column ordered
C MAXPSV Max. value in PERMSV
C MINPSV Min. positive value in PERMSV
C NLVL   Number of levels in level-set structure
C NLVL1  Number of levels in level-set structure from first end
C NODES  Number of nodes in current component
C NSTOP  Index of the ending node of current component
C NSTRT  Index of the starting node of current component
C NVARS  Number of variables in current component
C NWGHT  Normalized weights
C XLS    See comment for IW.

C     .. External Subroutines ..
      EXTERNAL MC60HD,MC60JD,MC60LD

C Partition IW
      XLS = NSUP
      LIST = 2*NSUP + 1

C Initializations
      NWGHT(1) = WEIGHT(1)
      NWGHT(2) = WEIGHT(2)
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      INFO(4) = 0
      NVARS = 0
      LSTNUM = 0
      LZNUM = NSUP+1

      IF(JCNTL(2).NE.2) THEN
C  Set PERMSV = 1 to denote all nodes being visible.
         DO 5 I = 1,NSUP
            PERMSV(I) = 1
    5    CONTINUE
      END IF

C     Order all nodes of degree zero
      DO 6 I = 1,NSUP
         K = ICPTR(I)
         DEGREE = ICPTR(I+1) - K
         IF (DEGREE.LE.1) THEN
           IF (DEGREE.EQ.0) THEN
             LZNUM = LZNUM - 1
             PERMSV(I) = -LZNUM
           ELSE IF (IRN(K).EQ.I) THEN
              LSTNUM = LSTNUM + 1
              PERMSV(I) = -LSTNUM
            END IF
         END IF
    6 CONTINUE

C     Loop while some nodes remain unnumbered
      DO 30 I = 1, NSUP
        IF (LSTNUM.GE.LZNUM-1) GO TO 35
        IF(JCNTL(2).EQ.0) THEN
C       Find pseudoperipheral pair of nodes for this component
          CALL MC60HD(N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +               IW(XLS+1),IW(LIST+1),HINFO)
          NSTRT = HINFO(1)
          NSTOP = HINFO(2)
          PAIR(1,I) = HINFO(1)
          PAIR(2,I) = HINFO(2)
          NLVL = HINFO(3)
          LWIDTH = HINFO(4)
          NVARS = HINFO(5)
          NODES = HINFO(6)
        ELSE IF(JCNTL(2).EQ.1) THEN
C Pseudoperipheral pair specified. Try first end.
          NSTOP = PAIR(1,I)
          CALL MC60LD(NSTOP,N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +            IW(XLS+1),NLVL,LWIDTH,NVARS)
          LWDTH1 = LWIDTH
          NLVL1 = NLVL
          NODES = IW(XLS+NLVL+1)-1
C    Try other end.
          NSTRT = NSTOP
          NSTOP = PAIR(2,I)
          CALL MC60LD(NSTOP,N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +            IW(XLS+1),NLVL,LWIDTH,NVARS)
          IF (NLVL1.GT.NLVL .OR.
     +             (NLVL1.EQ.NLVL .AND. LWDTH1.LT.LWIDTH)) THEN
C     Go back to other way
             NSTRT = NSTOP
             NSTOP = PAIR(1,I)
             CALL MC60LD(NSTOP,N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +            IW(XLS+1),NLVL,LWIDTH,NVARS)
          END IF
        ELSE
C Find max. and min. given priority value
          MAXPSV = 0
          DO 10 J = 1,NSUP
            IF (PERMSV(J).GT.0) THEN
               IF (PERMSV(J).GT.MAXPSV) THEN
                  IF (MAXPSV.EQ.0) THEN
                    MINPSV = PERMSV(J)
                    NSTRT = J
                  END IF
                  MAXPSV = PERMSV(J)
               END IF
               IF (PERMSV(J).LT.MINPSV) THEN
                  MINPSV = PERMSV(J)
                  NSTRT = J
               END IF
            END IF
   10     CONTINUE
          CALL MC60LD(NSTRT,N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +            IW(XLS+1),NLVL,LWIDTH,NVARS)
          NODES = IW(XLS+NLVL+1)-1
          IF(MAXPSV.NE.MINPSV) THEN
            NWGHT(2) = (WEIGHT(2)*(NLVL-1))/(MAXPSV-MINPSV)
          ELSE
            NWGHT(2) = WEIGHT(2)
          END IF
        END IF

        INFO(1) = INFO(1) + 1
        IF (NVARS.GT.INFO(2)) THEN
           INFO(2) = NVARS
           INFO(3) = NLVL
           INFO(4) = LWIDTH
        END IF

        IF(JCNTL(1).EQ.1) THEN
C Renumber nodes in this component by RCM
C First regenerate level structure, ordering neighbours by increasing
C degree (last time level set constructed it was with NSTOP as root)
           CALL MC60QD(NSTOP,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +            IW(XLS+1),NLVL,NVARS,IW(LIST+1))
           DO 11 J = NODES,1,-1
             LSTNUM = LSTNUM + 1
             PERMSV(IW(J)) = -LSTNUM
   11      CONTINUE
        ELSE
          IF(JCNTL(2).NE.2) THEN
C     Store level-set indices as global priorities
            DO 15 IL = 1,NLVL
              DO 12 J = IW(XLS+IL), IW(XLS+IL+1) - 1
                PERMSV(IW(J)) = NLVL - IL
   12         CONTINUE
   15       CONTINUE
          END IF

C         Order nodes in this component
          CALL MC60JD(NSUP,LIRN,NODES,NSTRT,LSTNUM,IRN,ICPTR,
     +         VARS,PERMSV,NWGHT,IW,IW(LIST+1),IW(XLS+1),W)
        END IF

   30 CONTINUE

C     Set new node numbers to +ve values
   35 DO 40 I = 1,NSUP
          PERMSV(I) = -PERMSV(I)
   40 CONTINUE

      END


      SUBROUTINE MC60DD(N,NSUP,SVAR,VARS,PERMSV,PERM,POSSV)
C Given supervariable permutation, find corresponding variable
C     permutation.

      INTEGER N
C N has intent IN and must be set to the matrix order.
      INTEGER NSUP
C NSUP has intent IN. It holds the number of supervariables.
      INTEGER SVAR(N)
C SVAR has intent IN. SVAR(I) is the supervariable to which variable
C     I belongs, I = 1, ..., N
      INTEGER VARS(NSUP)
C VARS has intent IN and holds the numbers of variables in the
C     supervariables.
      INTEGER PERMSV(NSUP)
C PERMSV has intent IN. Positions of supervariables in
C     supervariable order.
      INTEGER PERM(N)
C PERM has intent OUT. Returns positions of variables in variable
C     order.
      INTEGER POSSV(NSUP)
C POSSV has intent OUT. Set to positions of supervariables
C     in variable order.
      INTEGER I,IS,L

C Find inverse permutation
      DO 10 IS = 1,NSUP
         PERM(PERMSV(IS)) = IS
   10 CONTINUE

C Accumulate lengths
      L = 1
      DO 20 I = 1,NSUP
         L = L + VARS(PERM(I))
         POSSV(PERM(I)) = L
   20 CONTINUE

C Construct full permutation
      DO 30 I = 1,N
         IS = SVAR(I)
         L = POSSV(IS)-1
         POSSV(IS) = L
         PERM(I) = L
   30 CONTINUE

      END


      SUBROUTINE MC60ED(N,NSUP,LIRN,IRN,ICPTR,SVAR,VARS,PERMSV,PERM,IW)

C Given a supervariable elimination order for a symmetric
C compressed matrix, find
C the corresponding row presentation order, as required by MA42.

C   N      - (IN) Integer variable. Must be set to the order of A.
C   NSUP   - (IN) Integer variable. Must be set to the order of
C            the condensed matrix (that is to the number of
C            supervariables). Set NSUP to N if A not condensed.
C   LIRN   - (IN) Integer variable. Length of array IRN.
C   IRN    - (IN) Integer array of length LIRN. On entry it
C            must hold the row indices of the entries
C            in the condensed matrix, ordered by columns.
C   ICPTR  - (IN) Integer array of length NSUP+1.
C            On entry, ICPTR(J) holds the position in the
C             array IRN of the first entry in column J (J=1,...,NSUP)
C            and the last entry is in position ICPTR(NSUP+1)-1.
C   SVAR   - (IN) Integer array SVAR(N).
C            If N.NE.NSUP, on entry SVAR(I) must hold
C            the supervariable to which variable
C            I belongs, I = 1, ..., N
C            If NSUP = N,  SVAR(I) is not used.
C   VARS   - (IN) Integer array of length NSUP.
C            If N.NE.NSUP, on entry, VARS(IS) holds the number of
C            variables in supervariable IS, IS = 1, ..., NSUP.
C            If NSUP = N,  VARS(IS) is not used.
C  PERMSV  - (INOUT) Integer array of length NSUP. On entry, must hold
C            the supervariable permutation (the new index for
C            supervariable I must be given by PERMSV(I),
C            I = 1,...,NSUP). On exit, the order in which the rows
C            of the condensed matrix should be
C            presented is PERMSV(1), PERMSV(2), ..., PERMSV(NSUP).
C  PERM    - (OUT) Integer array of length N.
C            On exit, the order in which the rows should be
C            presented is PERM(1), PERM(2), ..., PERM(N).
C   IW    -  (OUT) Integer array of length NSUP. Used as workspace.


C     .. Scalar Arguments ..
      INTEGER N,NSUP,LIRN

C     .. Array Arguments ..
      INTEGER IRN(LIRN),ICPTR(NSUP+1),PERM(N),PERMSV(NSUP),IW(NSUP),
     *        SVAR(N),VARS(NSUP)

C     .. Local Scalars ..
      INTEGER I,IS,JS,K,L,M
C I      Row index
C IS     Row index of compressed matrix
C JS     Column index of compressed matrix
C K      Assembly index
C L      Elimination index
C M      DO index


C First, convert the supervariable elimination order to a presentation
C order for the rows of the compressed matrix.

C Set IW to hold the supervariable elimination order
C (ie supervariables are eliminated in the order
C IW(1), IW(2), ...,IW(NSUP).)
C Also, set PERM(I) (I=1,...,NSUP) to 0.
      DO 10 IS = 1,NSUP
        JS = PERMSV(IS)
        IW(JS) = IS
        PERM(IS) = 0
   10 CONTINUE

C In the next loop PERM(JS) holds the presentation position of row JS
C of the compressed matrix.
      K = 0
C Loop over rows of the compressed matrix
      DO 30 L = 1,NSUP
        IS = IW(L)
C Loop over entries in row IS of the compressed matrix
        DO 20 M = ICPTR(IS), ICPTR(IS+1)-1
          JS = IRN(M)
          IF(PERM(JS).EQ.0) THEN
C Row JS of the compressed matrix has not yet been presented ...
C put it next in the presentation order
            K = K + 1
            PERM(JS) = K
            PERMSV(K) = JS
          END IF
   20   CONTINUE
        IF(PERM(IS).EQ.0) THEN
C Row IS of the compressed matrix has not yet been presented ...
C put it next in the presentation order
            K = K + 1
            PERM(IS) = K
            PERMSV(K) = IS
          END IF
   30 CONTINUE

C Deal with the special case N=NSUP.
C In this case we just need to set PERM
      IF (N.EQ.NSUP) THEN
         DO 40 I = 1,N
            PERM(I) = PERMSV(I)
   40    CONTINUE
         RETURN
      END IF

C Accumulate lengths  (so that IW(JS) points to the
C position after the presentation position of the last variable
C belonging to supervariable JS)
      L = 1
      DO 45 IS = 1,NSUP
         JS = PERMSV(IS)
C The number of variables belonging to supervariable JS is VARS(JS)
         L = L + VARS(JS)
         IW(JS) = L
   45 CONTINUE

C Construct full presentation order
      DO 50 I = 1,N
         IS = SVAR(I)
         L = IW(IS) - 1
         IW(IS) = L
         PERM(L) = I
   50 CONTINUE


      END


      SUBROUTINE MC60FD(N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,RINFO)

C Code to compute the profile, maximum and r.m.s wavefronts and
C bandwidth of a permutation of A.
C A must have a SYMMETRIC pattern (BUT no checks are made).
C A held in condensed form if NSUP.ne.N.

C   N      - (IN) Integer variable. Must be set to the order of A.
C   NSUP   - (IN) Integer variable. Must be set to the order of
C            the condensed matrix (that is to the number of
C            supervariables). Set NSUP to N if A not condensed.
C   LIRN   - (IN) Integer variable. Length of array IRN.
C   IRN    - (IN) Integer array of length LIRN. On entry it
C            must hold the row indices of the entries
C            in the condensed matrix, ordered by columns.
C   ICPTR  - (IN) Integer array of length NSUP+1.
C            ICPTR(J) holds the position in the array IRN
C            of the first entry in column J (J=1,...,NSUP)
C            and the last entry is in position ICPTR(NSUP+1)-1.
C   VARS   - (IN) Integer array of length NSUP.
C            VARS(IS) must hold the number of
C            variables in supervariable IS, IS = 1, ..., NSUP.
C   PERMSV - (IN) Integer array of length NSUP.  On entry, must hold the
C            positions of the supervariables. If the original order
C            is required, the user must set PERMSV(I) = I,
C            I = 1,...,NSUP.
C   IW    -  (OUT) Integer array of length 2*NSUP+1. Used as workspace.
C            ABS(IW(I)), I=1,NSUP holds the supervariable in position I.
C               IW(II) is negated if supervariable IS is in the front.
C            IW(NSUP+I) holds the number of variables ahead of
C               supervariable I in the order.
C  RINFO  -  (OUT) Real (DP) array length 4. On exit, RINFO(1) holds the
C            profile, RINFO(2) holds the maximum wavefront, RINFO(3)
C            holds the bandwidth, and RINFO(4) holds the r.m.s.
C            wavefront for permuted ordering.

C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)

C     .. Scalar Arguments ..
      INTEGER LIRN,N,NSUP

C     .. Array Arguments ..
      INTEGER IRN(LIRN),IW(2*NSUP+1),PERMSV(NSUP),ICPTR(NSUP+1),
     *        VARS(NSUP)
      DOUBLE PRECISION RINFO(4)

C     .. Local Scalars ..
      INTEGER I,IMIN,J,JSTOP,JSTRT,K,NACTIV,NBR,NV
C I      Do index
C IMIN   Lowest numbered neighbour
C J      Current supervariable in pivot order
C JSTOP  Last position in IRN of current column
C JSTRT  Start position in IRN of current column
C K      Index in IRN
C NACTIV Number of active (frontal) variables
C NBR    Neighbour of current supervariable J
C NV     Nuber of variables in current supervariable J


C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,MAX,MIN,SQRT

C     ..
C
C Initialise the frontsize and profile
      DO 10 I = 1,4
         RINFO(I) = ZERO
  10  CONTINUE

      NACTIV = 0
      IW(NSUP+1) = 0

C Set IW to be the inverse permutation i.e. IW(PERMSV(I)) = I
C (if IW(J) = I, original number for node J is I).

      DO 30 I = 1,NSUP
         J = PERMSV(I)
         IW(J) = I
   30 CONTINUE

C Main loop. Loop over all nodes
      DO 80 I = 1,NSUP
         J = ABS(IW(I))
C NV is number of variables in supervariable J
         NV = VARS(J)
         IW(NSUP+I+1) = IW(NSUP+I) + NV
         JSTRT = ICPTR(J)
         JSTOP = ICPTR(J+1) - 1

         IMIN = I + 1
C Loop over all neighbours of J
         DO 50 K = JSTRT,JSTOP
            NBR = IRN(K)
C Find lowest numbered neighbour
            IMIN = MIN(IMIN,PERMSV(NBR))
            IF (IW(NBR).GT.0) THEN
C NBR not yet in front. Put it in the front.
               NACTIV = NACTIV +  VARS(NBR)
               IW(NBR) = -IW(NBR)
            END IF
   50    CONTINUE

C Update statistics
         RINFO(3) = MAX(RINFO(3),DBLE(IW(NSUP+I+1)-IW(NSUP+IMIN)))
         IF (IW(J).GT.0) THEN
C J not yet in front. Treat is diagonal block as a diagonal matrix.
           IW(J) = -IW(J)
           RINFO(2) = MAX(RINFO(2),DBLE(NACTIV+1))
           RINFO(1) = RINFO(1) + NV*DBLE(NACTIV+1)
           RINFO(4) = RINFO(4) + NV*DBLE(NACTIV+1)**2
         ELSE
           RINFO(2) = MAX(RINFO(2),DBLE(NACTIV))
           DO 70 J = 1, NV
             RINFO(1) = RINFO(1) + DBLE(NACTIV)
             RINFO(4) = RINFO(4) + DBLE(NACTIV)**2
C Remove node from the front
             NACTIV = NACTIV - 1
   70      CONTINUE
         END IF
   80 CONTINUE

      RINFO(4) = SQRT(RINFO(4)/DBLE(N))

      END



      SUBROUTINE MC60GD(N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,RINFO)

C Given a row-by-row order, this subroutine computes the maximum
C and r.m.s frontsizes which will arise if MA42 is employed.
C The order in which the equations are input to MA42 is
C PERMSV(1),PERMSV(2),...,PERMSV(NSUP).
C A must have a SYMMETRIC pattern (BUT no checks are made).
C A held in condensed form if NSUP.ne.N.
C
C   N      - (IN) Integer variable. Must be set to the order of A.
C   NSUP   - (IN) Integer variable. Must be set to the order of
C            the condensed matrix (that is to the number of
C            supervariables). Set NSUP to N if A not condensed.
C   LIRN   - (IN) Integer variable. Length of array IRN.
C   IRN    - (IN) Integer array of length LIRN. On entry it
C            must hold the row indices of the entries
C            in the condensed matrix, ordered by columns.
C   ICPTR  - (IN) Integer array of length NSUP+1.
C            On entry, ICPTR(J) holds the position in
C            the array IRN of the first entry in column J (J=1,...,NSUP)
C            and the last entry is in position ICPTR(NSUP+1)-1.
C            Unchanged on exit.
C   VARS   - (INOUT) Integer array of length NSUP.
C            VARS(IS) must hold the number of
C            variables in supervariable IS, IS = 1, ..., NSUP.
C   PERMSV - (IN) Integer array of length NSUP.  On entry, must hold the
C            presentation order for the condensed matrix
C            (that is, the rows of the condensed matrix are to
C            be presented in the order PERMSV(1), PERMSV(2),...,
C            PERMSV(NSUP)). If the natural order is required,
C             the user must set PERMSV(I) = I, I = 1,...,NSUP.
C   IW    -  (OUT) Integer array of length NSUP. Used as workspace.
C  RINFO  -  (OUT) Real (DP) array length 4. On exit, RINFO(1)
C            and RINFO(2) hold the maximum row and column frontsizes,
C            and RINFO(3) holds root mean square row front size
C            and RINFO(4) holds the mean frontal matrix size.

C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)

C     .. Scalar Arguments ..
      INTEGER LIRN,N,NSUP

C     .. Array Arguments ..
      DOUBLE PRECISION RINFO(4)
      INTEGER ICPTR(NSUP+1),IRN(LIRN),IW(NSUP),PERMSV(NSUP),VARS(NSUP)

C     .. Local Scalars ..
      INTEGER I,IEQ,J,JSTOP,JSTRT,K,KFRNT,LFRNT,MFR,NV
C I
C IEQ
C J
C JSTOP  Last position in IRN of current column
C JSTRT  First position in IRN of current column
C K
C KFRNT
C LFRNT
C MFR
C NV

C     .. Intrinsic Functions ..
      INTRINSIC DBLE,MAX,SQRT

C Initialisation
      DO 10 I = 1,4
         RINFO(I) = ZERO
   10 CONTINUE

C Initialise IW(I) (I=1,...,NSUP) to 0.
      DO 20 I = 1,NSUP
         IW(I) = 0
   20 CONTINUE
C
C Initialise row and column frontsizes
      KFRNT = 0
      LFRNT = 0

C Find out when each supervariable appears for the last time.
C For supervariable MFR this is stored in IW(MFR)
      DO 40 IEQ = 1,NSUP
         I = PERMSV(IEQ)
         IW(I) = IEQ
         JSTRT = ICPTR(I)
         JSTOP = ICPTR(I+1) - 1
         DO 30 J = JSTRT,JSTOP
            MFR = IRN(J)
            IW(MFR) = IEQ
   30    CONTINUE
   40 CONTINUE

C Symbolic factorization
      DO 90 IEQ = 1,NSUP
C
C Assemble incoming equation
         I = PERMSV(IEQ)
         JSTRT = ICPTR(I)
         JSTOP = ICPTR(I+1) - 1
C Take no action for null row
         IF (JSTRT.GT.JSTOP) GO TO 90
C NV is number of variables in supervariable I
         NV = VARS(I)
C Update row front size.
         KFRNT = KFRNT + NV

         IF (IW(I).GE.0) THEN
C Supervariable I not yet in front. Add it.
            LFRNT = LFRNT + NV
            IW(I) = -IW(I)
         END IF
C Loop over supervariables in incoming eqn.
         DO 60 J = JSTRT,JSTOP
            MFR = IRN(J)
            IF (IW(MFR).GE.0) THEN
C Supervariable not yet in front. Add it.
               NV = VARS(MFR)
               LFRNT = LFRNT + NV
               IW(MFR) = -IW(MFR)
            END IF
   60    CONTINUE

C Update row/column front sizes.
         RINFO(1) = MAX(RINFO(1),DBLE(KFRNT))
         RINFO(2) = MAX(RINFO(2),DBLE(LFRNT))
         IF (-IW(I).EQ.IEQ) THEN
C Supervariable I can be eliminated.
            IW(I) = 0
            NV = VARS(I)
            DO 65 K = 1,NV
               RINFO(3) = RINFO(3) + DBLE(KFRNT)**2
               RINFO(4) = RINFO(4) + DBLE(KFRNT)*DBLE(LFRNT)
               LFRNT = LFRNT - 1
               KFRNT = KFRNT - 1
   65       CONTINUE
         END IF
         DO 80 J = JSTRT,JSTOP
            MFR = IRN(J)
            IF (-IW(MFR).EQ.IEQ) THEN
C Supervariable MFR can be eliminated.
               NV = VARS(MFR)
               DO 70 K = 1,NV
                  RINFO(3) = RINFO(3) + DBLE(KFRNT)**2
                  RINFO(4) = RINFO(4) + DBLE(KFRNT)*DBLE(LFRNT)
                  LFRNT = LFRNT - 1
                  KFRNT = KFRNT - 1
   70          CONTINUE
            END IF
   80    CONTINUE

   90 CONTINUE

C Compute r.m.s row front size.
      RINFO(3) = SQRT(RINFO(3)/DBLE(N))
C Compute mean frontal matrix size.
      RINFO(4) = RINFO(4)/DBLE(N)

      END




      SUBROUTINE MC60HD(N,NSUP,LIRN,IRN,ICPTR,VARS,MASK,LS,XLS,LIST,
     +                 INFO)
C Find pseudoperipheral pair of nodes for a component of the graph.

C     .. Scalar Arguments ..
      INTEGER N,NSUP,LIRN

C     .. Array Arguments ..
      INTEGER ICPTR(NSUP+1),IRN(LIRN),LIST(NSUP),LS(NSUP),
     +        MASK(NSUP),VARS(NSUP),XLS(NSUP+1),INFO(6)

C N has intent IN and holds the matrix order.
C NSUP has intent IN and holds the number of nodes.
C LIRN has intent IN and holds the length of the array IRN.
C IRN has intent IN and holds the adjacency lists for all the nodes.
C ICPTR has intent IN. ICPTR(J) holds the position in
C     IRN of the first entry in the list for node J (J=1,...,NSUP)
C     and the last entry is in position ICPTR(NSUP+1)-1.
C VARS has intent IN and holds the number of variables at each node.
C MASK has intent INOUT, but is always restored to its input value.
C     MASK(NODE) = 1 for each visible node, including the root.
C LS has intent OUT. Each time a level set is generated, it
C     is set to hold the indices on nodes at level 1,
C     followed by those at level 2, ...
C XLS has intent OUT. Each time a level set is generated,
C     XLS(I) is set to hold the position in
C     LS of the first entry for level I (J=1,...,NLVL)
C     and the last entry is in position ICPTR(NLVL+1)-1.
C     It also holds the degrees of the nodes on the bottom
C     level of the  rooted structure.
C LIST is workspace used to hold the set of possible start
C     nodes.
C INFO has intent OUT.
C   INFO(1) set to starting pseudoperipheral node.
C   INFO(2) set to ending pseudoperipheral node.
C   INFO(3) set to number of levels in final level set.
C   INFO(4) set to the width of the final level set.
C   INFO(5) set to number of variables in this component of graph.
C   INFO(6) set to number of nodes in this component of graph.

C     .. Local Scalars ..
      INTEGER DEGREE,I,J,LSIZE,LWIDTH,MAIN,MAXDEP,
     +        MINDEG,MINWID,NLSIZE,NLVL,NODE,NODES,NSTOP,NSTRT,NVARS
C DEGREE  Degree of node
C I       DO index
C J       DO index
C LSIZE   Size of the final level set
C LWIDTH  Width of level set
C MAIN    DO index of main loop.
C MAXDEP  Maximum level-set depth found so far.
C MINDEG  Minimum degree of node in final level set.
C MINWID  Minimum level-set width.
C NLSIZE  No. of nodes of differing degrees in final level set
C NLVL    No. levels in current level set
C NODE    Index of graph node.
C NODES   Number of nodes in this component of graph.
C NSTOP   Ending pseudoperipheral node.
C NSTRT   Starting pseudoperipheral node.
C NVARS   Number of variables in current component

C     .. External Subroutines ..
      EXTERNAL MC60LD

C     Choose first guess for starting node by min supervariable degree,
C     ignoring nodes that are invisible (MASK NE 1)
      MINDEG = N+1
      INFO(5) = 0
      DO 10 I = 1,NSUP
          IF (MASK(I).EQ.1) THEN
            INFO(5) = INFO(5) + VARS(I)
            DEGREE = ICPTR(I+1) - ICPTR(I)
            IF (DEGREE.LE.MINDEG) THEN
              IF (DEGREE.LT.MINDEG) THEN
                  NSTRT = I
                  MINDEG = DEGREE
              END IF
            END IF
          END IF
   10 CONTINUE

C     Generate level structure for node with min supervariable degree
      CALL MC60LD(NSTRT,N,NSUP,LIRN,IRN,ICPTR,VARS,MASK,LS,XLS,
     +            MAXDEP,LWIDTH,NVARS)

C     Store number of nodes in this component
      NODES = XLS(MAXDEP+1) - 1
      NSTOP = 0

C     Iterate with different start nodes
      DO 70 MAIN = 1, NODES

C Record the width  from NSTRT
        INFO(4) = LWIDTH

C Store  nodes in the final level set and their degrees
        LSIZE = 0
        DO 30 I = XLS(MAXDEP),XLS(MAXDEP+1) - 1
          NODE = LS(I)
          LSIZE = LSIZE + 1
          LIST(LSIZE) = NODE
          XLS(NODE) = ICPTR(NODE+1) - ICPTR(NODE)
   30   CONTINUE

C Choose at most 5 nodes

        DO 50 NLSIZE = 1,5
C Look for candidate with least degree
           MINDEG = N+1
           DO 41 I = NLSIZE,LSIZE
             IF(XLS(LIST(I)).LT.MINDEG)THEN
                J = I
                MINDEG = XLS(LIST(I))
             END IF
   41     CONTINUE
C Jump out of loop if no cnadidates left
          IF(MINDEG.EQ.N+1) GO TO 55
C Swap chosen candidate to next position
          NODE = LIST(J)
          LIST(J) = LIST(NLSIZE)
          LIST(NLSIZE) = NODE
C Rule out the neighbours of the chosen node
          DO 42 I = ICPTR(NODE), ICPTR(NODE+1)-1
             XLS(IRN(I)) = N+1
   42     CONTINUE
   50   CONTINUE
   55   NLSIZE = NLSIZE-1

C  Loop over nodes in list
        MINWID = N

        DO 60 I = 1,NLSIZE
          NODE = LIST(I)

C  Form rooted level structures NODE
          CALL MC60LD(NODE,MINWID,NSUP,LIRN,IRN,ICPTR,VARS,MASK,LS,XLS,
     +               NLVL,LWIDTH,NVARS)
          IF (LWIDTH.LT.MINWID) THEN
            IF (NLVL.GT.MAXDEP) THEN
C             Level structure of greater depth. Begin a new iteration.
              NSTRT = NODE
              MAXDEP = NLVL
              GO TO 70
            ELSE
C             Level structure of lesser width. Record it.
              NSTOP = NODE
              MINWID = LWIDTH
            END IF
          END IF
   60   CONTINUE
        GO TO 80
   70 CONTINUE

C Swap the ends if this reduces the width
   80 IF (INFO(4) .LT. MINWID) THEN
         INFO(1) = NSTRT
         NSTRT = NSTOP
         NSTOP = INFO(1)
      END IF

      IF(NSTOP.NE.NODE) CALL MC60LD(NSTOP,N,NSUP,LIRN,IRN,ICPTR,VARS,
     +               MASK,LS,XLS,NLVL,LWIDTH,NVARS)
      INFO(1) = NSTRT
      INFO(2) = NSTOP
      INFO(3) = MAXDEP
      INFO(4) = LWIDTH
      INFO(5) = NVARS
      INFO(6) = NODES

      END



      SUBROUTINE MC60JD(NSUP,LIRN,NODES,NSTRT,LSTNUM,IRN,ICPTR,VARS,
     +                  STATUS,WEIGHT,NLIST,QUEUE,DEG,PRIOR)

C     Renumber nodes in component of graph for small wavefront and
C     profile. Supervariable version of Sloan's code.

C     .. Scalar Arguments ..
      INTEGER LIRN,LSTNUM,NSUP,NODES,NSTRT

C     .. Array Arguments ..
      INTEGER DEG(NSUP),ICPTR(NSUP+1),IRN(LIRN),NLIST(NSUP),
     +        QUEUE(0:NODES-1),STATUS(NSUP),VARS(NSUP)
      DOUBLE PRECISION WEIGHT(2),PRIOR(NSUP)

C NSUP has intent IN. It holds the number of supervariables.
C LIRN  has intent IN. It holds the length of IRN.
C NODES has intent IN and holds the number of nodes in this component.
C NSTRT has intent IN and holds the index of the starting node.
C LSTNUM has intent INOUT and holds the number of reindexed nodes.
C IRN has intent IN. It must hold the row indices of the
C     supervariable representation of the matrix, by columns.
C ICPTR has intent IN. It must be set so that ICPTR(J) holds the
C     position in the array IRN of the first entry in column J
C     (J=1,...,NSUP) and the last entry is in position ICPTR(NSUP+1)-1.
C VARS has intent INOUT. It must hold the number of variables in
C     supervariable IS, IS = 1, ..., NSUP. It is altered, but restored
C     to its given value.
C STATUS has intent INOUT. If node I is in this component,
C     STATUS(I) holds:
C        global priority value on entry
C        2 if inactive
C        1 if preactive
C        0 if active (in front)
C        -(its new index) if postactive (renumbered).
C     Other components are not altered.
C WEIGHT has intent IN. It holds the weights for the priorities.
C NLIST has intent INOUT. On entry, it holds the list of nodes in
C     this component. If heap is used, overwritten by positions of
C     the nodes on the heap.
C QUEUE is workspace. Holds the queue of nodes that are currently
C     active or preactive.
C DEG is workspace. It holds the current degrees of the nodes.
C PRIOR is workspace. It holds the priorities of nodes in the queue.

C     .. Local Scalars ..
      INTEGER ADDRES,DEGREE,FATHER,I,ISTOP,ISTRT,J,JSTOP,JSTRT,J1,J2,
     +        K,L,NABOR,NBR,NEXT,NODE,NQ,QNODE,SON,THRESH
      PARAMETER (THRESH=100)
      DOUBLE PRECISION MAXPRT,PNODE,PRTY
C ADDRES  Position in queue of chosen node
C DEGREE  Degree of node
C FATHER  Father node in heap
C I       Do index
C ISTOP   End of index list in IRN
C ISTRT   Start of index list in IRN
C J       Do index
C JSTOP   End of index list in IRN
C JSTRT   Start of index list in IRN
C J1
C J2
C K
C L       Index for main loop
C MAXPRT  Max priority value
C NABOR   Neighbouring node to NBR
C NBR     Neighbouring node to NEXT
C NEXT    Next node to be renumbered
C NODE    Node in queue
C NQ      Size of queue
C PNODE   Priority value of QNODE
C PRTY    Priority value
C QNODE   Node to be inserted in heap
C SON     Son node in heap
C THRESH  Threshold for switching to the heap.

C     Initialise priorities and status for each node in this component
C     Priority is
C     -WEIGHT(1)*DEGREE - WEIGHT(2)*DIST
C     DIST   = distance of node from end node
C     DEGREE = initial current degree for node

      DO 10 I = 1,NODES
         NODE = NLIST(I)
         DEGREE = VARS(NODE)
         K = DEGREE
         VARS(NODE) = 0
         DO 7 J = ICPTR(NODE),ICPTR(NODE+1) - 1
            DEGREE = DEGREE + VARS(IRN(J))
    7    CONTINUE
         VARS(NODE) = K
         PRIOR(NODE) = -WEIGHT(1)*DEGREE-WEIGHT(2)*STATUS(NODE)
         STATUS(NODE) = 2
         DEG(NODE) = DEGREE
   10 CONTINUE

C   Insert starting node in queue, assigning it a preactive status
      NQ = 1
      QUEUE(NQ) = NSTRT
      QUEUE(0) = NSTRT
      STATUS(NSTRT) = 1
      PRIOR(NSTRT) = 1.0E30

C     Loop while queue is not empty
      DO 70 L = 1, NODES
C If queue is too long, switch to heap
         IF(NQ.GT.THRESH)GO TO 100

C       Pick node with in queue with max priority
         ADDRES = 1
         MAXPRT = PRIOR(QUEUE(1))
         DO 30 I = 2,NQ
            PRTY = PRIOR(QUEUE(I))
            IF (PRTY.GT.MAXPRT) THEN
               ADDRES = I
               MAXPRT = PRTY
            END IF
   30    CONTINUE
         NEXT = QUEUE(ADDRES)

C       Delete node NEXT from queue
         QUEUE(ADDRES) = QUEUE(NQ)
         NQ = NQ - 1
         ISTRT = ICPTR(NEXT)
         ISTOP = ICPTR(NEXT+1) - 1
C        If NEXT is preactive, make adjustments corresponding to
C        making it active.
         IF (STATUS(NEXT).EQ.1) THEN
            DO 40 I = ISTRT,ISTOP
C           Decrease current degree of neighbour
               NBR = IRN(I)
               PRIOR(NBR) = PRIOR(NBR) + WEIGHT(1)*VARS(NEXT)
               DEG(NBR) = DEG(NBR) - VARS(NEXT)
C  If NBR has no neighbours outside the front, eliminate it next
               IF(DEG(NBR).EQ.0) PRIOR(NBR) = 1.0E30
C           If NBR is inactive, add it to queue as preactive
               IF (STATUS(NBR).EQ.2) THEN
                  NQ = NQ + 1
                  QUEUE(NQ) = NBR
                  STATUS(NBR) = 1
                  PRIOR(NBR) = PRIOR(NBR)
               END IF
   40       CONTINUE
         END IF

C       Assign new index for node NEXT
         LSTNUM = LSTNUM + 1
         STATUS(NEXT) = -LSTNUM

C       Make each preactive neighbour of NEXT active and make
C       corresponding adjustments to degrees of neighbours.
         DO 60 I = ISTRT,ISTOP
            NBR = IRN(I)
            IF (STATUS(NBR).NE.1) GO TO 60
C       Decrease the current degree of NBR and give it active status.
            PRIOR(NBR) = PRIOR(NBR) + WEIGHT(1)*VARS(NBR)
            STATUS(NBR) = -1
            DEG(NBR) = DEG(NBR) - VARS(NBR)
C  If NBR has no neighbours outside the front, eliminate it next
            IF(DEG(NBR).EQ.0) PRIOR(NBR) = 1.0E30
C        Loop over nodes adjacent to NBR
            JSTRT = ICPTR(NBR)
            JSTOP = ICPTR(NBR+1) - 1
            DO 50 J = JSTRT,JSTOP
               NABOR = IRN(J)
               IF (STATUS(NABOR).LT.0) GO TO 50
C            Decrease current degree of NABOR
               PRIOR(NABOR) = PRIOR(NABOR) + WEIGHT(1)*VARS(NBR)
               DEG(NABOR) = DEG(NABOR) - VARS(NBR)
C  If NABOR has no neighbours outside the front, eliminate it next
               IF(DEG(NABOR).EQ.0) PRIOR(NABOR) = 1.0E30
               IF (STATUS(NABOR).EQ.2) THEN
C           NABOR is currently inactive, but this node is now
C           adjacent to a newly activated node. Insert NABOR in
C           queue and assign it a preactive status
                  NQ = NQ + 1
                  QUEUE(NQ) = NABOR
                  STATUS(NABOR) = 1
               END IF
   50       CONTINUE
            STATUS(NBR) = 0
   60    CONTINUE
   70 CONTINUE
      RETURN

C   Work in heap mode
  100  DO 120 I = 1, NQ
         NBR = QUEUE(I)
         PNODE = PRIOR(NBR)
         J1 = I
         DO 116 K = 1, J1
            J2 = J1/2
            FATHER = QUEUE(J2)
            IF(PRIOR(FATHER).GE.PNODE) GO TO 118
            QUEUE(J1) = FATHER
            NLIST(FATHER) = J1
            J1 = J2
  116   CONTINUE
  118   QUEUE(J1) = NBR
        NLIST(NBR) = J1
  120 CONTINUE

C     Loop while queue is not empty
      I = L
      DO 170 L =I, NODES
         NEXT = QUEUE(1)

C       Delete node 1 from queue
         QNODE = QUEUE(NQ)
         PNODE = PRIOR(QNODE)
         NQ = NQ - 1
         J = 2
         J2 = 1
         IF(NQ.GT.1) QUEUE(NQ+1) = QUEUE(NQ)
         DO 125 I = 2, NQ
            IF(J.GT.NQ) GO TO 130
            IF( PRIOR(QUEUE(J)).LT.PRIOR(QUEUE(J+1)) ) J=J+1
            SON = QUEUE(J)
            IF(PNODE.GE.PRIOR(SON)) GO TO 130
            QUEUE(J2) = SON
            NLIST(SON) = J2
            J2 = J
            J = J*2
  125    CONTINUE
  130    QUEUE(J2) = QNODE
         NLIST(QNODE) = J2

         ISTRT = ICPTR(NEXT)
         ISTOP = ICPTR(NEXT+1) - 1
C        If NEXT is preactive, making adjustments corresponding to
C        making it active.
         IF (STATUS(NEXT).EQ.1) THEN
            DO 140 I = ISTRT,ISTOP
C           Decrease current degree of neighbour
               NBR = IRN(I)
               IF (NBR.EQ.NEXT) GO TO 140
               PRIOR(NBR) = PRIOR(NBR) + WEIGHT(1)*VARS(NEXT)
               DEG(NBR) = DEG(NBR) - VARS(NEXT)
C  If NBR has no neighbours outside the front, eliminate it next
               IF(DEG(NBR).EQ.0) PRIOR(NBR) = 1.0E30
C           If NBR is inactive, add it to queue as preactive
               IF (STATUS(NBR).EQ.2) THEN
                  NQ = NQ + 1
                  QUEUE(NQ) = NBR
                  STATUS(NBR) = 1
                  NLIST(NBR) = NQ
               END IF

C Adjust heap
               PNODE = PRIOR(NBR)
               J = NLIST(NBR)
               DO 133 K = 1, NQ
                  J2 = J/2
                  FATHER = QUEUE(J2)
                  IF(PRIOR(FATHER).GE.PNODE) GO TO 137
                  QUEUE(J) = FATHER
                  NLIST(FATHER) = J
                  J = J2
  133          CONTINUE
  137          QUEUE(J) = NBR
               NLIST(NBR) = J
  140       CONTINUE
         END IF

C       Assign new index for node NEXT
         LSTNUM = LSTNUM + 1
         STATUS(NEXT) = -LSTNUM

C       Make each preactive neighbour of NEXT active and make
C       corresponding adjustments to degrees of neighbours.
         DO 160 I = ISTRT,ISTOP
            NBR = IRN(I)
            IF (STATUS(NBR).NE.1) GO TO 160
C       Decrease the current degree of NBR and give it active status.
            PRIOR(NBR) = PRIOR(NBR) + WEIGHT(1)*VARS(NBR)
            STATUS(NBR) = -1
            DEG(NBR) = DEG(NBR) - VARS(NBR)
C  If NBR has no neighbours outside the front, eliminate it next
            IF(DEG(NBR).EQ.0) PRIOR(NBR) = 1.0E30
C Adjust heap
            PNODE = PRIOR(NBR)
            J = NLIST(NBR)
            DO 142 K = 1, NQ
               J2 = J/2
               FATHER = QUEUE(J2)
               IF(PRIOR(FATHER).GE.PNODE) GO TO 144
               QUEUE(J) = FATHER
               NLIST(FATHER) = J
               J = J2
  142       CONTINUE
  144       QUEUE(J) = NBR
            NLIST(NBR) = J

C        Loop over nodes adjacent to NBR
            JSTRT = ICPTR(NBR)
            JSTOP = ICPTR(NBR+1) - 1
            DO 150 J = JSTRT,JSTOP
               NABOR = IRN(J)
               IF (STATUS(NABOR).LT.0) GO TO 150
C            Decrease current degree of NABOR
               PRIOR(NABOR) = PRIOR(NABOR) + WEIGHT(1)*VARS(NBR)
               DEG(NABOR) = DEG(NABOR) - VARS(NBR)
C  If NABOR has no neighbours outside the front, eliminate it next
               IF(DEG(NABOR).EQ.0) PRIOR(NABOR) = 1.0E30
               IF (STATUS(NABOR).EQ.2) THEN
C           NABOR is currently inactive, but this node is now
C           adjacent to a newly activated node. Insert NABOR in
C           queue and assign it a preactive status
                  NQ = NQ + 1
                  QUEUE(NQ) = NABOR
                  STATUS(NABOR) = 1
                  NLIST(NABOR) = NQ
               END IF
C Adjust heap
               PNODE = PRIOR(NABOR)
               J1 = NLIST(NABOR)
               J2 = J1/2
               FATHER = QUEUE(J2)
               IF(PRIOR(FATHER).GE.PNODE) GO TO 148
               QUEUE(J1) = FATHER
               NLIST(FATHER) = J1
               J1 = J2
               DO 146 K = 2, NQ
                  J2 = J1/2
                  FATHER = QUEUE(J2)
                  IF(PRIOR(FATHER).GE.PNODE) GO TO 148
                  QUEUE(J1) = FATHER
                  NLIST(FATHER) = J1
                  J1 = J2
  146          CONTINUE
  148          QUEUE(J1) = NABOR
               NLIST(NABOR) = J1
  150       CONTINUE
            STATUS(NBR) = 0
  160    CONTINUE
  170 CONTINUE

      END



      SUBROUTINE MC60LD(ROOT,MAXWID,NSUP,LIRN,IRN,ICPTR,VARS,MASK,LS,
     +                 XLS,NLVL,LWIDTH,NVARS)
C Generate rooted level structure of width less than MAXWID.

C     .. Scalar Arguments ..
      INTEGER LWIDTH,MAXWID,NSUP,NLVL,LIRN,ROOT,NVARS
C     ..
C     .. Array Arguments ..
      INTEGER ICPTR(NSUP+1),IRN(LIRN),LS(NSUP),MASK(NSUP),
     *        VARS(NSUP),XLS(NSUP+1)
C
C ROOT has intent IN and holds the root node for level structure.
C MAXWID has intent IN and holds the maximum permissible width of
C     rooted level structure.
C NSUP has intent IN and holds the number of nodes.
C LIRN has intent IN and holds the length of the array IRN.
C IRN has intent IN and holds the adjacency lists for all the nodes.
C ICPTR has intent IN. ICPTR(J) holds the position in
C     IRN of the first entry in the list for node J (J=1,...,NSUP)
C     and the last entry is in position ICPTR(NSUP+1)-1.
C VARS has intent IN and holds the number of variables at each node.
C MASK has intent INOUT, but is always restored to its input value.
C     MASK(NODE) > 0 for each visible node, including the root.
C LS has intent OUT. It is set to hold the indices on nodes at
C     level 1, followed by those at level 2, ...
C XLS has intent OUT. XLS(I) is set to hold the position in
C     LS of the first entry for level I (J=1,...,NLVL)
C     and the last entry is in position ICPTR(NLVL+1)-1.
C NLVL has intent OUT. It is set to the number of levels.
C LWIDTH has intent OUT. It is set to the width of the structure.
C NVARS has intent OUT. It is set to the number of variables in
C     current component

C     .. Local Scalars ..
      INTEGER I,J,LBEGIN,LNBR,LVLEND,LW,NBR,NODE
C  I      index in previous level
C  J      index in IRN
C  LBEGIN beginning of present level in LS
C  LNBR   next position in LS
C  LVLEND end of previous level in LS
C  LW     level width
C  NBR    neighbour node
C  NODE   current node

C     .. Intrinsic Functions ..
      INTRINSIC MAX

C Establish level 1
      MASK(ROOT) = -MASK(ROOT)
      LS(1) = ROOT
      LVLEND = 0
      NVARS = 0
      LNBR = 1
      LWIDTH = VARS(ROOT)
      DO 35 NLVL = 1,NSUP
C     Generate next level by finding all unmasked neighbours
C     of nodes in present level
          LBEGIN = LVLEND + 1
          LVLEND = LNBR
          XLS(NLVL) = LBEGIN
          LW = 0
          DO 30 I = LBEGIN,LVLEND
              NODE = LS(I)
              DO 20 J = ICPTR(NODE), ICPTR(NODE+1) - 1
                  NBR = IRN(J)
                  IF (MASK(NBR).GT.0) THEN
                      LNBR = LNBR + 1
                      LS(LNBR) = NBR
                      MASK(NBR) = -MASK(NBR)
                      LW = LW + VARS(NBR)
                  END IF
   20         CONTINUE
   30     CONTINUE
          LWIDTH = MAX(LW, LWIDTH)
          NVARS = NVARS + LW
C If no neighbours found, we are done.
          IF (LNBR.EQ.LVLEND) GO TO 40
C  Abort construction if level structure too wide.
          IF (LWIDTH.GE.MAXWID) GO TO 40
   35 CONTINUE

   40 XLS(NLVL+1) = LVLEND + 1
C     Reset MASK for nodes in the level structure
      DO 50 I = 1,LNBR
          MASK(LS(I)) = ABS(MASK(LS(I)))
   50 CONTINUE

      END


      SUBROUTINE MC60OD(N,NC,LIRN,IRN,ICPTR,SVAR,NSUP,NEW,VARS,FLAG)

C Find supervariables consisting of variables whose columns have
C    the same pattern. This is done by updating the supervariable
C    structure for the first J-1 columns to that for the first J
C    columns, J=1,2,...,N. If a supervariable IS is partially
C    involved in column J, the variables that are involved
C    are removed to make a new supervariable, leaving the rest
C    in IS.

C No check is made on the validity of the data.
C    There must be no duplicate entries or out-of-range indices.
C    Columns with no entries are permitted.

      INTEGER N
C N has intent IN and must be set to the number of rows.
      INTEGER NC
C NC has intent IN and must be set to the number of columns.
      INTEGER LIRN
C LIRN has intent IN and must be set to the length of IRN.
      INTEGER IRN(LIRN)
C IRN has intent IN and must be set to hold the row indices of the
C     entries by columns. They may be in any order within the column.
      INTEGER ICPTR(NC+1)
C ICPTR has intent IN. ICPTR(J) must be set to the position in
C     IRN of the first entry of column J, J = 1,..., NC and
C     the last entry is in position ICPTR(NC+1)-1.
      INTEGER SVAR(N)
C SVAR has intent OUT. On return, SVAR(I) is the supervariable
C     to which variable I belongs, I = 1, ..., N
      INTEGER NSUP
C NSUP has intent OUT. On return, it holds the number of
C     supervariables.
      INTEGER NEW(N)
C NEW is used as workspace. NEW(IS) is the new supervariable
C     constructed from variables of supervariable IS.
      INTEGER VARS(N)
C VARS is used as workspace. VARS(IS) is the number of
C     variables in supervariable IS.
      INTEGER FLAG(N)
C FLAG is used as workspace. FLAG(IS) = J if supervariable IS has
C     been encountered for column J.

C Local Scalars.
      INTEGER I,IS,J,JS,K,K1,K2
C I    Variable index.
C IS   Supervariable.
C JS   Supervariable.
C J    Column index.
C K    DO index.
C K1,K2   Temporary variables.

C Begin by setting SVAR and VARS to represent all variables
C     belonging to supervariable 1. Initialize NEW and FLAG.
      DO 10 I = 1,N
         SVAR(I) = 1
   10 CONTINUE
      VARS(1) = N
      FLAG(1) = 0
      NSUP = 1


C Scan the columns in turn, splitting the supervariables that are
C     involved.
      DO 40 J = 1,NC
         K1 = ICPTR(J)
         K2 = ICPTR(J+1) - 1
C Loop over variables in column J, decrementing the counts of
C      variables in supervariables.
         DO 20 K = K1,K2
            IS = SVAR(IRN(K))
            VARS(IS) = VARS(IS) - 1
   20    CONTINUE

C Loop over variables in column J, incrementing the count
C   and resetting VARS.
         DO 30 K = K1,K2
            I = IRN(K)
            IS = SVAR(I)
            IF (FLAG(IS).LT.J) THEN
C First occurrence of supervariable IS for column J
               FLAG(IS) = J
               IF (VARS(IS).GT.0) THEN
C Establish new supervariable
                  NSUP = NSUP + 1
                  VARS(NSUP) = 1
                  FLAG(NSUP) = J
                  NEW(IS) = NSUP
                  SVAR(I) = NSUP
               ELSE
C No new supervariable needed
                  VARS(IS) = 1
                  NEW(IS) = IS
                  SVAR(I) = IS
               END IF
            ELSE
C Subsequent occurrence of IS for column J
               JS = NEW(IS)
               VARS(JS) = VARS(JS) + 1
               SVAR(I) = JS
            END IF
   30    CONTINUE

   40 CONTINUE

      END



      SUBROUTINE MC60PD(N,NSUP,LIRN,IRN,ICPTR,SVAR,VARS,VAR,FLAG)
C Change the variable to supervariable mapping so that FIRST(IS),
C     IS=1,...,NSUP is monotonic where FIRST(IS) is the least index
C     of a variable in supervariable IS. Compress the pattern to
C     refer to supervariables.

      INTEGER N
C N has intent IN and must be set to the matrix order.
      INTEGER LIRN
C NSUP has intent IN. It holds the number of supervariables.
      INTEGER IRN(LIRN)
C LIRN has intent IN and must be set to the length of IRN.
      INTEGER NSUP
C IRN has intent INOUT. On entry, it holds the row indices of the
C     entries by columns. It is overwritten by the
C     corresponding data for the compressed matrix.
      INTEGER ICPTR(N+1)
C ICPTR has intent INOUT. ICPTR(I) must be set to the position in
C     IRN of the first entry of column J, J = 1,..., N and the last
C     entry is in position ICPTR(N+1)-1. It is overwritten by the
C     corresponding data for the compressed matrix.
      INTEGER SVAR(N)
C SVAR has intent INOUT. On entry, SVAR(I) is the supervariable to
C     which variable I belongs, I = 1, ..., N. Changed so that
C     mapping to earliest variable is monotonic.
      INTEGER VARS(NSUP)
C VARS has intent OUT. Returns numbers of variables in
C     supervariables.
      INTEGER VAR(NSUP)
C VAR has intent OUT. Returns supervariable to earliest variable
C     mapping.
      INTEGER FLAG(NSUP)
C FLAG is a workarray

      INTEGER I,IS,J,JS,K,K1,L

C Initialize FLAG and VARS
      DO 10 IS = 1,NSUP
         FLAG(IS) = -1
         VARS(IS) = 1
   10 CONTINUE

C Find supervariable mapping. FLAG holds old to new supervariable
C     map.
      L = 1
      DO 20 I = 1,N
         IS = SVAR(I)
         JS = FLAG(IS)
         IF(JS.GT.0)THEN
C Second or later occurrence of supervariable
            SVAR(I) = JS
            VARS(JS) = VARS(JS) + 1
         ELSE IF(JS.LT.0)THEN
C First occurrence of supervariable
            FLAG(IS) = L
            VAR(L) = I
            SVAR(I) = L
            L = L + 1
         END IF
   20 CONTINUE

C Re-initialize FLAG
      DO 30 IS = 1,NSUP
         FLAG(IS) = 0
   30 CONTINUE

C Do actual compression
      L = 1
      K1 = 1
      DO 60 JS = 1,NSUP
         J = VAR(JS)
         K1 = ICPTR(J)
         ICPTR(JS) = L
         DO 50 K = K1, ICPTR(J+1)-1
            IS = SVAR(IRN(K))
            IF(FLAG(IS).NE.JS)THEN
               FLAG(IS) = JS
               IRN(L) = IS
               L = L+1
            END IF
   50    CONTINUE
   60 CONTINUE
      ICPTR(JS) = L
      END

      SUBROUTINE MC60QD(ROOT,NSUP,LIRN,IRN,ICPTR,VARS,MASK,LS,
     +                 XLS,NLVL,NVARS,DEGREE)
C Generate rooted level structure.
C Additionally, order the entries within each level set by increasing 
C current degree

C     .. Scalar Arguments ..
      INTEGER NSUP,NLVL,LIRN,ROOT,NVARS
C     ..
C     .. Array Arguments ..
      INTEGER ICPTR(NSUP+1),IRN(LIRN),LS(NSUP),MASK(NSUP),
     *        VARS(NSUP),XLS(NSUP+1),DEGREE(NSUP)
C
C ROOT has intent IN and holds the root node for level structure.
C NSUP has intent IN and holds the number of nodes.
C LIRN has intent IN and holds the length of the array IRN.
C IRN has intent IN and holds the adjacency lists for all the nodes.
C ICPTR has intent IN. ICPTR(J) holds the position in
C     IRN of the first entry in the list for node J (J=1,...,NSUP)
C     and the last entry is in position ICPTR(NSUP+1)-1.
C VARS has intent IN and holds the number of variables at each node.
C MASK has intent INOUT, but is always restored to its input value.
C     MASK(NODE) > 0 for each visible node, including the root.
C LS has intent OUT. It is set to hold the indices on nodes at
C     level 1, followed by those at level 2, ...
C XLS has intent OUT. XLS(I) is set to hold the position in
C     LS of the first entry for level I (J=1,...,NLVL)
C     and the last entry is in position ICPTR(NLVL+1)-1.
C NLVL has intent OUT. It is set to the number of levels.
C NVARS has intent OUT. It is set to the number of variables in
C     current component

C     .. Local Scalars ..
      INTEGER FNBR,I,II,J,K,L,LBEGIN,LNBR,LPERM,LVLEND,LW,NBR,NODE
C  FNBR   position in LS of first neighbour of current node.
C  I      index in previous level
C  J      index in IRN
C  LBEGIN beginning of present level in LS
C  LNBR   position in LS of last neighbour of current node
C  LVLEND end of previous level in LS
C  LW     level width
C  NBR    neighbour node
C  NODE   current node

C Establish level 1
      MASK(ROOT) = -MASK(ROOT)
      LS(1) = ROOT
      LVLEND = 0
      NVARS = 0
      LNBR = 1
      DO 35 NLVL = 1,NSUP
C     Generate next level by finding all unmasked neighbours
C     of nodes in present level
          LBEGIN = LVLEND + 1
          LVLEND = LNBR
          XLS(NLVL) = LBEGIN
          LW = 0
          DO 30 I = LBEGIN,LVLEND
              NODE = LS(I)
              FNBR = LNBR + 1
              DO 20 J = ICPTR(NODE), ICPTR(NODE+1) - 1
                  NBR = IRN(J)
                  IF (MASK(NBR).GT.0) THEN
                      LNBR = LNBR + 1
                      LS(LNBR) = NBR
                      MASK(NBR) = -MASK(NBR)
                      LW = LW + VARS(NBR)
                  END IF
   20         CONTINUE
              IF (FNBR.EQ.LNBR+1) THEN
C No unmasked neighbours so skip to next node in this level
                 GO TO 30
              END IF
          
C Sort the neighbours of NODE in increasing order by 
C current degree, where current degree of a node is the number
C of unmasked neighbours.
C Loop over the nodes in the level set just constructed
C and compute their current degrees
              DO 16 II = FNBR,LNBR
                 NODE = LS(II)
                 DEGREE(NODE) = 0
                 DO 15 J = ICPTR(NODE),ICPTR(NODE+1)-1
                    NBR = IRN(J)
                    IF (MASK(NBR).GT.0) DEGREE(NODE) = DEGREE(NODE) + 1
   15            CONTINUE
   16         CONTINUE 

C Linear insertion is used (because we want to sort in place as no
C more storage available).
              DO 27 K = FNBR+1,LNBR
                 NBR = LS(K)
                 DO 24 L =  K-1,FNBR+1,-1
                    LPERM = LS(L)
                    IF (DEGREE(LPERM).LE.DEGREE(NBR)) GO TO 25
                    LS(L+1) = LPERM
   24            CONTINUE
   25            LS(L+1) = NBR
   27         CONTINUE

C finished with NODE
   30     CONTINUE

          NVARS = NVARS + LW
C If no neighbours found, we are done.
          IF (LNBR.EQ.LVLEND) GO TO 40
   35 CONTINUE

   40 XLS(NLVL+1) = LVLEND + 1
C     Reset MASK for nodes in the level structure
      DO 50 I = 1,LNBR
          MASK(LS(I)) = ABS(MASK(LS(I)))
   50 CONTINUE

      END
! COPYRIGHT (c) 1979 AEA Technology and
! Council for the Central Laboratory of the Research Councils
!
! Version 1.0.1
! See ChangeLog for version history
!
      DOUBLE PRECISION FUNCTION FA14AD(IX,I)
C         NEARLY PORTABLE RANDOM NUMBER GENERATOR USING THE RECURSION
C                       IX=IX*A MOD P
C
C    WHERE A=7**5
C    AND P=2**31-1.
C
C         THIS FUNCTION DOES NOT ADHERE TO THE ANSI STANDARD 1966 IN
C    TWO RESPECTS:
C      1) IT ASSUMES AN INTEGER WORD LENGTH OF AT LEAST 32 BITS (I.E.
C    INTEGERS WHICH LIE IN THE RANGE 1-2**31 TO 2**31-1 INCLUSIVE MUST
C    BE REPRESENTABLE);
C      2) IT ASSUMES THAT A POSITIVE INTEGER LESS THAN 2**16 MAY BE
C    FLOATED WITHOUT LOSS OF DIGITS.
C
C         THIS CODE IS BASED ON CODE PUBLISHED BY LINUS SCHRAGE IN
C    T.O.M.S. VOL.5 NO.2 JUNE 1979 (PP 132-138)
C
C
C
C       THE FUNCTION IS USED AS FOLLOWS:
C
C                      R=FA14AD(IX,I)
C
C       WHERE IX IS THE GENERATOR WORD
C             I IS AN INTEGER SET BY THE USER.
C
C
C       THE VALUE RETURNED BY FA14A/AD WILL LIE IN THE RANGE
C                (0.,1.)  IF I IS NON-NEGATIVE
C                (-1.,1.) IF I IS NEGATIVE.
C
C       THE METHOD EMPLOYED IS A MULTIPLICATIVE CONGRUENTIAL
C   ONE USING A MULTIPLIER OF 7**5 AND TAKING THE MODULO TO
C   2**31-1, I.E. THE GENERATOR NUMBER , G = IX, IS UPDATED ON
C   EACH CALL TO THE VALUE
C
C                  5          31
C               G*7  MODULO (2  -1)
C
C       THE RESULT RETURNED IS CALCULATED AS A DOUBLE
C  PRECISION NUMBER HAVING THE VALUE
C
C                      31
C                  G/(2   -1)    IF THE ARGUMENT IS
C                                NON-NEGATIVE
C           OR
C                      31
C                2*G/(2   -1)-1  IF THE ARGUMENT IS NEGATIVE
C
C
C 7**5, 2**15, 2**16, 2**31-1
C     .. Parameters ..
      INTEGER A,B15,B16,P
      PARAMETER (A=16807,B15=32768,B16=65536,P=2147483647)
C     ..
C     .. Scalar Arguments ..
      INTEGER IX,I
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION X
      INTEGER FHI,K,LEFTLO,XALO,XHI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC FLOAT
C     ..
C     .. Executable Statements ..
C
C GET 15 HI ORDER BITS OF IX
      XHI = IX/B16
C GET 16 LO BITS OF IX AND FORM LO PRODUCT
      XALO = (IX-XHI*B16)*A
C GET 15 HI ORDER BITS OF LO PRODUCT
      LEFTLO = XALO/B16
C     FORM THE 31 HIGHEST BITS OF FULL PRODUCT
      FHI = XHI*A + LEFTLO
C GET OVERFLOPAST 31ST BIT OF FULL PRODUCT
      K = FHI/B15
C ASSEMBLE ALL THE PARTS AND PRESUBTRACT P
C THE PARENTHESES ARE ESSENTIAL
      IX = (((XALO-LEFTLO*B16)-P)+ (FHI-K*B15)*B16) + K
C ADD P BACK IN IF NECCESSARY
      IF (IX.LT.0) IX = IX + P
C MULTIPLY BY 1/(2**31-1)
      XHI = IX/B16
      X = (FLOAT(XHI)*65536.0D0) + FLOAT(IX-XHI*B16)
      IF (I.GE.0) FA14AD = X*4.6566128752457969241D-10
      IF (I.LT.0) FA14AD = X*9.3132257504915938482D-10 - 1.0D0
      RETURN

      END
      SUBROUTINE FA14BD(IX,MAX,NRAND)
C         NEARLY PORTABLE RANDOM NUMBER GENERATOR USING THE RECURSION
C                       IX=IX*A MOD P
C
C    WHERE A=7**5
C    AND P=2**31-1.
C
C         THIS SUBROUTINE DOES NOT ADHERE TO THE ANSI STANDARD 1966
C    IN ONE RESPECT:
C         IT ASSUMES AN INTEGER WORD LENGTH OF AT LEAST 32 BITS (I.E.
C    INTEGERS WHICH LIE IN THE RANGE 1-2**31 TO 2**31-1 INCLUSIVE MUST
C    BE REPRESENTABLE).
C
C         THIS CODE IS BASED ON CODE PUBLISHED BY LINUS SCHRAGE IN
C    T.O.M.S. VOL.5 NO.2 JUNE 1979 (PP 132-138)
C
C
C       THE FUNCTION IS USED AS FOLLOWS:
C
C                  CALL FA14BD(IX,MAX,NRAND)
C
C       WHERE IX    IS THE GENERATOR WORD
C             MAX   IS AN INTEGER SET BY THE USER AND
C             NRAND IS AN INTEGER SET BY FA14B/BD.
C
C
C       THE VALUE OF NRAND RETURNED BY FA14B/BD WILL LIE IN THE
C   RANGE
C                        (1,MAX)
C
C       THE METHOD EMPLOYED IS A MULTIPLICATIVE CONGRUENTIAL
C   ONE USING A MULTIPLIER OF 7**5 AND TAKING THE MODULO TO
C   2**31-1, I.E. THE GENERATOR NUMBER , G = IX, IS UPDATED ON
C   EACH CALL TO THE VALUE
C
C                  5          31
C               G*7  MODULO (2  -1)
C
C       THE RESULT RETURNED IS AN INTEGER NUMBER
C   HAVING THE VALUE
C
C                        31
C   INT. PART( (MAX*G)/(2   -1) ) + 1
C
C
C 7**5, 2**15, 2**16, 2**31-1
C 2**30,  2**30-1
C     .. Parameters ..
      INTEGER A,B15,B16,P
      PARAMETER (A=16807,B15=32768,B16=65536,P=2147483647)
      INTEGER B30,Q
      PARAMETER (B30=1073741824,Q=1073741823)
C     ..
C     .. Scalar Arguments ..
      INTEGER IX,MAX,NRAND
C     ..
C     .. Local Scalars ..
      INTEGER BE1,BE2,C,D,F,FHI,G,K,LEFTLO,MHI,MLO,MU,NU,XALO,XHI,XLO
C     ..
C     .. Executable Statements ..
C
C GET 15 HI ORDER BITS OF IX
      XHI = IX/B16
C GET 16 LO BITS OF IX AND FORM LO PRODUCT
      XALO = (IX-XHI*B16)*A
C GET 15 HI ORDER BITS OF LO PRODUCT
      LEFTLO = XALO/B16
C     FORM THE 31 HIGHEST BITS OF FULL PRODUCT
      FHI = XHI*A + LEFTLO
C GET OVERFLOPAST 31ST BIT OF FULL PRODUCT
      K = FHI/B15
C ASSEMBLE ALL THE PARTS AND PRESUBTRACT P
C THE PARENTHESES ARE ESSENTIAL
      IX = (((XALO-LEFTLO*B16)-P)+ (FHI-K*B15)*B16) + K
C ADD P BACK IN IF NECCESSARY
      IF (IX.LT.0) IX = IX + P
C MULTIPLY BY MAX AND DIVIDE BY 2**31-1 IN INTEGER ARITHMETIC
C SPLIT IX AND MAX INTO HI AND LO PARTS
      XHI = IX/B15
      XLO = IX - B15*XHI
      MHI = MAX/B15
      MLO = MAX - B15*MHI
C CALCULATE INTERMEDIATE PRODUCT AND SPLIT INTO HI AND LO PARTS
C PRESUBTRACT P
      F = (XHI*MLO-P) + XLO*MHI
C F IS > 0 IF INTERMEDIATE PRODUCT WOULD HAVE OVERFLOWED
      IF (F.GT.0) GO TO 1
      F = F + P
      BE1 = F/B15
      BE2 = F - BE1*B15
      GO TO 2

    1 F = F - 1
      BE1 = F/B15
      BE2 = F - BE1*B15
      BE1 = BE1 + B16
C FORM PRODUCT OF LO PARTS AND ADD IN LO PART OF INTERMEDIATE PRODUCT
C TO GET LO PART OF COMPLETE PRODUCT
    2 G = B15*BE2 + XLO*MLO
C REPRESENT LO PART OF FULL PRODUCT IN BASE 2**30
      D = G/B30
      C = XHI/2
C CALCULATE FULL PRODUCT DIVIDED BY 2**30
      F = ((2* (C*MHI-Q)-1)+MHI* (XHI-2*C)) + D + BE1
C GET FULL PRODUCT DIVIDED IN BASE 2**31
      IF (F.GT.0) GO TO 3
      F = F + P
      NU = F/2
      MU = F - NU*2
      GO TO 4

    3 F = F - 1
      NU = F/2
      MU = F - 2*NU
      NU = NU + B30
C CALCULATE REMAINDER OF PRODUCT DIVIDED BY 2**31
    4 F = (B30*MU-P) + NU + (G-B30*D)
      NRAND = NU + 1
C  ADD ONE IF REMAINDER IS NOT < 2**31-1
      IF (F.GE.0) NRAND = NRAND + 1
      RETURN

      END
      SUBROUTINE FA14CD(IX,IGEN)
C        FA14CD IS A SUBROUTINE USED IN CONJUNCTION WITH FA14AD OR
C   FA14BD. IT PROVIDES THE USER WITH THE FACILITY OF SAVING THE
C   CURRENT VALUE OF THE GENERATOR NUMBER USED BY FA14AD AND FA14BD.
C
C        USE OF THE ROUTINE IS AS FOLLOWS:
C
C                       CALL FA14CD(IGEN)
C
C     WHERE IX   IS THE GENERATOR WORD
C           IGEN IS AN INTEGER WHICH IS SET BY FA14C/CD TO THE CURRENT
C                VALUE OF THE GENERATOR.
C
C
C     .. Scalar Arguments ..
      INTEGER IX,IGEN
C     ..
C     .. Executable Statements ..
      IGEN = IX
      RETURN

      END
      SUBROUTINE FA14DD(IX,IGEN)
C        FA14DD IS A SUBROUTINE USED IN CONJUNCTION WITH FA14AD OR
C   FA14BD. IT PROVIDES THE USER WITH THE FACILITY OF SETTING THE
C   CURRENT VALUE OF THE GENERATOR NUMBER USED BY FA14AD AND FA14BD.
C
C        USE OF THE ROUTINE IS AS FOLLOWS:
C
C                       CALL FA14DD(IGEN)
C
C    WHERE IX   IS THE GENERATOR WORD
C          IGEN IS AN INTEGER, SET BY THE USER TO THE VALUE TO WHICH
C               THE GENERATOR IS TO BE SET. IT IS RECOMMENDED THAT THIS
C               VALUE BE OBTAINED BY A PREVIOUS CALL TO FA14C/CD.
C
C     .. Scalar Arguments ..
      INTEGER IX,IGEN
C     ..
C     .. Executable Statements ..
      IX = IGEN
      RETURN

      END
      SUBROUTINE FA14ID(IX)
C        FA14ID IS A SUBROUTINE USED TO INITIALIZE THE GENERATOR WORD US
C   BY FA14AD AND FA14BD. IT MUST BE CALLED FIRST BEFORE ANY OF THE OTHE
C   ENTRIES ARE CALLED.
C
C        USE OF THE ROUTINE IS AS FOLLOWS:
C
C                       CALL FA14ID(IX)
C
C    WHERE IX   IS THE GENERATOR WORD
C
C     .. Scalar Arguments ..
      INTEGER IX
C     ..
C     .. Executable Statements ..
      IX = 1
      RETURN

      END
* COPYRIGHT (c) 1980 AEA Technology
* Original date 1 December 1993
C       Toolpack tool decs employed.
C       Arg dimensions changed to *.
C 1/4/99 Size of MARK increased to 100.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE KB07AD(COUNT,N,INDEX)
C
C             KB07AD      HANDLES DOUBLE PRECISION VARIABLES
C  THE WORK-SPACE 'MARK' OF LENGTH 100 PERMITS UP TO 2**50 NUMBERS
C  TO BE SORTED.

C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION COUNT(*)
      INTEGER INDEX(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AV,X
      INTEGER I,IF,IFK,IFKA,INT,INTEST,IP,IS,IS1,IY,J,K,K1,LA,LNGTH,M,
     +        MLOOP
C     ..
C     .. Local Arrays ..
      INTEGER MARK(100)
C     ..
C     .. Executable Statements ..
C  SET INDEX ARRAY TO ORIGINAL ORDER .
      DO 10 I = 1,N
        INDEX(I) = I
   10 CONTINUE
C  CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED .
      IF (N.EQ.1) GO TO 200
      IF (N.GE.1) GO TO 30
      WRITE (6,FMT=20)

   20 FORMAT (/,/,/,20X,' ***KB07AD***NO NUMBERS TO BE SORTED ** ',
     + 'RETURN TO CALLING PROGRAM')

      GO TO 200
C  'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER
C  THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.
   30 M = 12
C  SET UP INITIAL VALUES.
      LA = 2
      IS = 1
      IF = N
      DO 190 MLOOP = 1,N
C  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE .
        IFKA = IF - IS
        IF ((IFKA+1).GT.M) GO TO 70
C********* FINAL SORTING ***
C  ( A SIMPLE BUBBLE SORT )
        IS1 = IS + 1
        DO 60 J = IS1,IF
          I = J
   40     IF (COUNT(I-1).LT.COUNT(I)) GO TO 60
          IF (COUNT(I-1).GT.COUNT(I)) GO TO 50
          IF (INDEX(I-1).LT.INDEX(I)) GO TO 60
   50     AV = COUNT(I-1)
          COUNT(I-1) = COUNT(I)
          COUNT(I) = AV
          INT = INDEX(I-1)
          INDEX(I-1) = INDEX(I)
          INDEX(I) = INT
          I = I - 1
          IF (I.GT.IS) GO TO 40
   60   CONTINUE
        LA = LA - 2
        GO TO 170
C             *******  QUICKSORT  ********
C  SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS
C  THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENT'S
C  HIGHEST ADDRESS.
   70   IY = (IS+IF)/2
        X = COUNT(IY)
        INTEST = INDEX(IY)
        COUNT(IY) = COUNT(IF)
        INDEX(IY) = INDEX(IF)
C  THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END
C  OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE
C  OF X .
        K = 1
        IFK = IF
C  WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE
C  INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AND INDICES AS
C  NECESSARY, UNTIL THEY MEET .
        DO 110 I = IS,IF
          IF (X.GT.COUNT(I)) GO TO 110
          IF (X.LT.COUNT(I)) GO TO 80
          IF (INTEST.GT.INDEX(I)) GO TO 110
   80     IF (I.GE.IFK) GO TO 120
          COUNT(IFK) = COUNT(I)
          INDEX(IFK) = INDEX(I)
          K1 = K
          DO 100 K = K1,IFKA
            IFK = IF - K
            IF (COUNT(IFK).GT.X) GO TO 100
            IF (COUNT(IFK).LT.X) GO TO 90
            IF (INTEST.LE.INDEX(IFK)) GO TO 100
   90       IF (I.GE.IFK) GO TO 130
            COUNT(I) = COUNT(IFK)
            INDEX(I) = INDEX(IFK)
            GO TO 110

  100     CONTINUE
          GO TO 120

  110   CONTINUE
C  RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER
C  WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO
C  2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN OR EQUAL
C  TO ANY ELEMENT IN THE SECOND PART, AND THEY MAY NOW BE SORTED
C  INDEPENDENTLY .
  120   COUNT(IFK) = X
        INDEX(IFK) = INTEST
        IP = IFK
        GO TO 140

  130   COUNT(I) = X
        INDEX(I) = INTEST
        IP = I
C  STORE THE LONGER SUBDIVISION IN WORKSPACE.
  140   IF ((IP-IS).GT. (IF-IP)) GO TO 150
        MARK(LA) = IF
        MARK(LA-1) = IP + 1
        IF = IP - 1
        GO TO 160

  150   MARK(LA) = IP - 1
        MARK(LA-1) = IS
        IS = IP + 1
C  FIND THE LENGTH OF THE SHORTER SUBDIVISION.
  160   LNGTH = IF - IS
        IF (LNGTH.LE.0) GO TO 180
C  IF IT CONTAINS MORE THAN ONE ELEMENT SUPPLY IT WITH WORKSPACE .
        LA = LA + 2
        GO TO 190

  170   IF (LA.LE.0) GO TO 200
C  OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT
  180   IF = MARK(LA)
        IS = MARK(LA-1)
  190 CONTINUE
  200 RETURN

      END
C COPYRIGHT (c) 1998 Council for the Central Laboratory
*                    of the Research Councils

C Version 1.1.0. See ChangeLog for history

      SUBROUTINE MC61ID(ICNTL,CNTL)

C ICNTL - (OUT) INTEGER array of length 10. On  exit,
C         ICNTL contains default values.
C         If the user wishes to use values other
C         than the defaults, the corresponding entries in
C         ICNTL should be reset after the call to
C         MC61I/ID.

C ICNTL(1) is the stream number for error messages and has the
C          default value 6. Printing of error messages is
C          suppressed if ICNTL(1) < 0.

C ICNTL(2) is the stream number for warning messages.
C          It has the default value 6. Printing of warning
C          messages is suppressed if ICNTL(2) < 0.

C ICNTL(3) controls the action taken if duplicate or out-of-range
C          entries are detected . If ICNTL(3) = 0 and such entries are
C          detected, the computation terminates
C          with IRN and ICPTR unchanged. If ICNTL(3) = 1, a warning
C          is issued and the computation continues.
C          The default value is 0.

C ICNTL(4) controls whether supervariables are to be used
C          (a supervariable is a set of variables that
C          correspond to a set of identical columns).  If
C          ICNTL(4) = 0, supervariables are used. If
C          ICNTL(4) = 1, variables are used.  If the problem has
C          significantly fewer supervariables than variables,
C          using supervariables will reduce the execution time
C          significantly and will, in general, produce a permutation
C          or assembly order of comparable quality.
C          The default value is 0.

C ICNTL(5) indicates whether the user wishes to supply
C          a global priority function in PERM.
C          If ICNTL(5) = 0, no priority function is supplied;
C          if ICNTL(5) = 1, a priority function is supplied.
C          The default value is 0.

C ICNTL(6) indicates whether the user wishes to supply
C          the weights for the priority function.
C          If ICNTL(6) = 0, no weights are supplied
C                           and (2,1) and (16,1) tried;
C          If ICNTL(6) = 1, no weights are supplied
C                           and (1,2) and (16,1) tried;
C          if ICNTL(6) = 2, weights are supplied in CNTL(1), CNTL(2).
C          The default value is 0.

C ICNTL(7) to ICNTL(10)  are given default values of zero.
C          They are currently not used but may be used in a
C          later release of the code.

C CNTL  - (OUT) REAL (DP) array of length 5. On  exit,
C         CNTL contains default values.
C         If the user wishes to use values other
C         than the defaults, the corresponding entries in
C         CNTL should be reset after the call to
C         MC61I/ID.
C         CNTL holds the weights in the priority function.
C         They are given the default values 2.0 and 1.0.
C         They are only used if ICNTL(6) = 2.
C CNTL(3) to CNTL(5)  are given default values of zero.
C          They are currently not used but may be used in a
C          later release of the code.


      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0D0)
      DOUBLE PRECISION CNTL(5)
      INTEGER ICNTL(10)
      INTEGER I

      DO 10 I = 1,10
         ICNTL(I) = 0
   10 CONTINUE
      ICNTL(1) = 6
      ICNTL(2) = 6

      CNTL(1) = 2.0D0
      CNTL(2) = 1.0D0
      CNTL(3) = ZERO
      CNTL(4) = ZERO
      CNTL(5) = ZERO

      RETURN
      END



      SUBROUTINE MC61AD(JOB,N,LIRN,IRN,ICPTR,PERM,LIW,IW,W,ICNTL,
     +                  CNTL,INFO,RINFO)

C JOB   - (IN) INTEGER variable that must be set by the user to 1 if
C         a variable permutation to reduce the profile and wavefront of
C         the matrix is required, to 2 if a variable permutation to
C         reduce the bandwidth is required, and to 3 if an assembly
C         order is required.

C N     - (IN) INTEGER variable that must be set by the user to the
C         order of the matrix.

C LIRN  - (IN) INTEGER variable that must be set by the user to the
C         length of the array IRN, which must be large enough to
C         hold the pattern of the whole matrix.

C IRN   - (INOUT) INTEGER  array of length LIRN whose leading part must
C         be set by the user to hold the row indices of the
C         entries in the lower triangle of the matrix,
C         including the diagonal. The entries of each column must
C         be contiguous. The entries of column J must precede those
C         of column J+1 (J=1,...,N-1), and there must be no wasted
C         space between the columns.  Row indices within a
C         column may be in any order. On exit, holds the row
C         entries of the condensed matrix, using the same format.

C ICPTR - (INOUT) INTEGER array of length N+1 that must be set by the
C         user so that ICPTR(J) points to the position in the array
C         IRN of the first entry in column J (J=1,...,N), and
C         ICPTR(N+1)-1  must be the position of the last entry.
C         On exit, ICPTR holds corresponding data for the condensed
C         matrix.

C PERM  - (INOUT) INTEGER array of length N.  This array must only
C         be set on entry if ICNTL(5) = 1 (the default is ICNTL(5) = 0.
C         In this case, this array must be set by the user to
C         hold the global priority function.  On exit, the new ordering
C         is contained in PERM.  If a variable permutation is
C         requested (JOB = 1 or 2), the new index for variable
C         I is given by PERM(I) (I = 1,...,N).  If an
C         assembly order is requested (JOB = 3), the order
C         in which the rows should be assembled is
C         PERM(1), PERM(2),...,PERM(N).

C LIW   - (IN) INTEGER variable that must be set by the user to the
C         length of the array IW. Sufficient value is LIW = 8*N+2
C         (7N+2 if ICNTL(5)=1 and 6N+2 if ICNTL(5)=1, ICNTL(6)=2).

C IW    - (OUT) INTEGER array of length LIW that is used by
C         the routine as workspace.

C W     - (OUT) REAL (DP) array of length N that is used by
C         the routine as workspace.

C ICNTL - (IN) INTEGER array of length 10. Holds values of control
C         parameters. Default values may be set by calling
C         MC61I/ID.

C CNTL  - (IN) REAl (DP) array of length 5. Holds values of control
C         parameters. Default values may be set by calling
C         MC61I/ID.

C INFO  - (OUT) INTEGER  array of length 10 that need not be set by the
C         user. On each successful exit, INFO(1) is set to 0.
C         Negative values of INFO(1) indicate a fatal error has
C         been detected and positive values indicate a warning
C         has been issued.
C INFO(1) = -1 JOB is not equal to 1, 2 or 3.
C            Immediate return with input parameters unchanged.
C INFO(1) = -2  N < 1.  Immediate return with input parameters
C            unchanged.
C INFO(1) = -3  LIRN is too small.  INFO(6)  is set to the minimum value
C           which will suffice for LIRN.
C            If LIRN is at least as large as the input value
C            of ICPTR(N+1)-1 and ICNTL(3) = 1,
C            any out-of-range or duplicated variable indices
C            will have been excluded from IRN and ICPTR.
C            Otherwise, the input parameters  are unchanged.
C INFO(1) = -4 LIW is too small.  INFO(3)  is set to a value
C            which may suffice for LIW.
C INFO(1) = -5  One or more variable indices either lies outside the
C            lower triangle  of the matrix or is duplicated.
C            Further information is contained in
C            INFO(4)> and  INFO(5).
C
C INFO(2) holds, on successful exit, the total number of
C         supervariables in the problem.  If variables are not used
C         (ICNTL(4) = 1), INFO(3) is set to N.
C INFO(3) holds the amount of workspace used by the routine. If the
C         user has provided insufficient workspace (INFO(1) = -4),
C        INFO(3) is set to a value which may suffice for LIW.
C INFO(4) holds the number of variable indices
C         in IRN found to be out-of-range.
C INFO(5) holds the number of duplicate variable indices in IRN.
C INFO(6) holds the minimum value which will suffice for LIRN.
C INFO(7) holds the number of non-trivial components
C         in the graph of the matrix
C INFO(8) to INFO(10)  are currently not used but may be used in a
C later release of the code.


C RINFO - (OUT) REAL (DP) array of length 15 that need not be set by the
C          user.

C If JOB= 1 or 2,  on successful exit RINFO
C returns the following information:
C RINFO(1) holds the profile of the matrix A.
C RINFO(2) holds the maximum wavefront of the matrix A.
C RINFO(3) holds the bandwidth of the matrix A.
C RINFO(4) holds the root mean squared wavefront of the matrix A.
C RINFO(5) holds the profile of the permuted  matrix.
C RINFO(6) holds the maximum wavefront of the permuted matrix.
C RINFO(7) holds the bandwidth  of the permuted matrix.
C RINFO(8) holds the root mean squared wavefront of the permuted matrix.

C If JOB = 3,  on successful exit RINFO
C returns the following information:
C RINFO(1) holds the maximum row frontsize for
C the assembly order 1, 2, ..., N.
C RINFO(2) holds the maximum column frontsize for
C the assembly order 1, 2, ..., N.
C RINFO(3) holds the root mean square row front size
C for the assembly order 1, 2, ..., N.
C RINFO(4) holds the mean frontal matrix size
C for the assembly order 1, 2, ..., N.
C RINFO(5) holds the maximum row frontsize for
C the new assembly order PERM(1), PERM(2), ..., PERM(N).
C RINFO(6) holds the maximum column frontsize for
C the new assembly order PERM(1), PERM(2), ..., PERM(N).
C RINFO(7) holds the root mean square row front size
C for the new assembly order PERM(1), PERM(2), ..., PERM(N).
C RINFO(8) holds the mean frontal matrix size
C for the new assembly order PERM(1), PERM(2), ..., PERM(N).

C If JOB = 1 or 3,
C RINFO(9) and RINFO(10) hold weights used.

C RINFO(11) to RINFO(15)  are currently not used but may be used in a
C later release of the code.

      DOUBLE PRECISION ZERO,ONE,TWO,SIXTN
      PARAMETER (ZERO = 0.0D0, ONE=1.0D0, TWO=2.0D0, SIXTN=16.0D0)

      INTEGER JOB,N,LIW,LIRN

      DOUBLE PRECISION RINFO(15)
      DOUBLE PRECISION CNTL(5),W(N)
      INTEGER IRN(LIRN),ICPTR(N+1),INFO(10),ICNTL(10),IW(LIW),PERM(N)

      DOUBLE PRECISION RNFO5,RNFO6,RNFO7,RNFO8
      INTEGER I,PERMSV,IWORK,LP,MP,NSUP,SVAR,VARS,IRUN,NRUN,COPY,
     +        IPERM,J,ISUP,JSUP,PAIR
      LOGICAL LSWAP

      INTEGER ICON60(2),INFO60(4),JCNTL(2)

      EXTERNAL MC60AD,MC60BD,MC60CD,MC60DD,MC60ED,MC60FD,MC60GD

C Streams for messages
      LP = ICNTL(1)
      MP = ICNTL(2)
C Initialise
      DO 10 I = 1,10
         INFO(I) = 0
   10 CONTINUE
      DO 20 I = 1,15
         RINFO(I) = ZERO
   20 CONTINUE

C Check input data for errors

      IF (JOB.NE.1 .AND. JOB.NE.2 .AND. JOB.NE.3) GO TO 100
      IF (N.LT.1) GO TO 110
      IF (LIRN.LT.ICPTR(N+1)-1) GO TO 115

C Call MC60A/AD to construct pattern of whole matrix.
C MC60A/AD also checks data for errors.
C We do not want messages to be printed from MC60A/AD.
      ICON60(1) = ICNTL(3)
      ICON60(2) = 0

C Check workspace. The workspace has to be at least 5*N+2.
      INFO(3) = 5*N + 2
      IF (5*N+2.GT.LIW) THEN
C Not even minimum required ... return a sufficient value.
         INFO(3) = 8*N + 2
         IF (JOB.NE.2) THEN
           IF (ICNTL(5).EQ.1) INFO(3) = 7*N + 2
           IF (ICNTL(5).EQ.1 .AND. ICNTL(6).EQ.2) INFO(3) = 6*N + 2
         END IF
         GO TO 130
      END IF

      CALL MC60AD(N,LIRN,IRN,ICPTR,ICON60,IW,INFO60)

C Copy information returned from MC60A/AD into INFO
      INFO(4) = INFO60(2)
      INFO(5) = INFO60(3)
      INFO(6) = INFO60(4)

C Check for errors
C Errors  -2 and -3 are possible.
      IF (INFO60(1).EQ.-2) GO TO 120
      IF (INFO60(1).EQ.-3) GO TO 140

C Check for warnings (only possible if ICNTL(3)=1).
      IF (INFO60(1).GT.0) THEN
C Duplicated or out-of-range entries have been found.
         INFO(1) = INFO60(1)
         IF (MP.GT.0) THEN
            IF (INFO(4).GT.0) WRITE (MP,'(/,A,I8,A)')
     * ' MC61A/AD warning:',INFO(4),' out-of-range entries are ignored'
            IF (INFO(5).GT.0) WRITE (MP,'(/,A,I8,A)')
     * ' MC61A/AD warning:',INFO(5),' duplicated entries are ignored'
         END IF
      END IF

C Put initial permutation/assembly order into IW
C (we can't put it into PERM as we must not overwrite any global
C priority function supplied by the user).
C We do not use supervariables initially (so need VARS(I) = 1
C in MC60F/FD).
      SVAR = 1
      VARS = SVAR + N
      IPERM = VARS + N
      IWORK = IPERM + N
C IWORK is length 2*N+1 ... already checked sufficient space.
      DO 30 I = 1, N
         IW(VARS+I-1) = 1
         IW(IPERM+I-1) = I
   30 CONTINUE

      IF (JOB.EQ.1 .OR. JOB.EQ.2) THEN
C Use MC60F/FD to compute the profile etc for original matrix
         CALL MC60FD(N,N,LIRN,IRN,ICPTR,IW(VARS),IW(IPERM),IW(IWORK),
     *               RINFO)
      ELSE
C Use MC60G/GD to compute statistics for original assembly order
         CALL MC60GD(N,N,LIRN,IRN,ICPTR,IW(VARS),IW(IPERM),IW(IWORK),
     *               RINFO)
      END IF


      IF (ICNTL(4).EQ.0) THEN
C Supervariables are to be used.
C Call MC60B/BD to find supervariables and compress pattern.
         IWORK = VARS + N
C Workspace length 2*N+2 required ...already checked we have
C this
         CALL MC60BD(N,LIRN,IRN,ICPTR,NSUP,IW(SVAR),IW(VARS),IW(IWORK))
      ELSE
C Supervariables are NOT used.
         NSUP = N
      END IF
      INFO(2) = NSUP

C Find permutation for variables or supervariables.
      JCNTL(1) = 0
      IF (JOB.EQ.2) JCNTL(1) = 1

      JCNTL(2) = 0
      IF (ICNTL(5).EQ.1 .AND. JOB.NE.2) JCNTL(2) = 2

C Does the user wish to set the weights?
      IF (JOB.EQ.1 .OR. JOB.EQ.3) THEN
         IF (ICNTL(6).EQ.0)  THEN
C We will try the pairs (2,1) and (16,1) and choose the best
            NRUN = 2
            RINFO(9) = TWO
            RINFO(10) = ONE
         ELSE IF (ICNTL(6).EQ.1)  THEN
C We will try the pairs (1,2) and (16,1) and choose the best
            NRUN = 2
            RINFO(9) = ONE
            RINFO(10) = TWO
         ELSE IF (ICNTL(6).EQ.2) THEN
C Use the weights in CNTL
            NRUN = 1
            RINFO(9) = CNTL(1)
            RINFO(10) = CNTL(2)
         END IF
      ELSE
C No weights needed if JOB = 2 (RCM)
         NRUN = 1
      END IF

C Check workspace
C Workspace must be length 3*NSUP+1.
C We only need an array to copy into if we are going to try
C two sets of weights.
      PERMSV = VARS + NSUP
      IWORK = PERMSV + NSUP
      PAIR = IWORK + 3*NSUP + 1
      COPY = PAIR + NSUP
      IF (ICNTL(5).EQ.1 .AND. JOB.NE.2) COPY = PAIR
      IF (NRUN.EQ.1) INFO(3) = MAX(INFO(3),COPY)
      IF (NRUN.EQ.2) INFO(3) = MAX(INFO(3),COPY+N-1)
      IF (INFO(3).GT.LIW) GO TO 130

      DO 60 IRUN = 1,NRUN

         IF (IRUN.EQ.1) THEN
C If user has supplied global priority function, copy PERM(I) into
C IW(PERMSV+I-1).
C If we are working with supervariables, we have to convert given
C priority function into a supervariable priority function.
            IF (ICNTL(5).EQ.1) THEN
               IF (ICNTL(4).EQ.1) THEN
C Working with variables
                  DO 34 I = 1,N
                     J = PERM(I)
                     IW(PERMSV+I-1) = J
   34             CONTINUE
               ELSE
C Working with supervariables
C SVAR(I) is the supervariable to which variable I belongs.
                  DO 35 I = 1,N
                     J = PERM(I)
                     ISUP = IW(SVAR+I-1)
                     JSUP = IW(SVAR+J-1)
                     IW(PERMSV+ISUP-1) = JSUP
   35             CONTINUE
               END IF
            END IF
            IF (NRUN.EQ.2) THEN
C Take a copy of IW(PERMSV+I-1) for use with the second pair of weights.
               DO 36 I = 1,NSUP
                  IW(COPY+I-1) = IW(PERMSV+I-1)
   36          CONTINUE
            END IF
         END IF

C If we are trying a second pair of weights, we must take a copy
C of the permutation we obtained with the first and reset the weights
         IF (IRUN.EQ.2) THEN
C
            IF (ICNTL(5).EQ.1) THEN
C copy what is in IW(COPY+I-1) into IW(PERMSV+I-1)
               DO 37 I = 1,NSUP
                  IW(PERMSV+I-1) = IW(COPY+I-1)
   37          CONTINUE
            END IF
C Take a copy of PERM
            DO 40 I = 1,N
               IW(COPY+I-1) = PERM(I)
   40       CONTINUE
            RNFO5 = RINFO(5)
            RNFO6 = RINFO(6)
            RNFO7 = RINFO(7)
            RNFO8 = RINFO(8)

            RINFO(9) = SIXTN
            RINFO(10) = ONE

C On second run, the start and end nodes are already known
C (since returned on the first run) ... make use of this
C in the case when global priority function NOT supplied.
            IF (ICNTL(5).NE.1) JCNTL(2) = 1

         END IF

         CALL MC60CD(N,NSUP,LIRN,IRN,ICPTR,IW(VARS),JCNTL,IW(PERMSV),
     *              RINFO(9),IW(PAIR),INFO60,IW(IWORK),W)
         INFO(7) = INFO60(1)

         IF (JOB.EQ.1 .OR. JOB.EQ.2) THEN
C Variable permutation required.

C If supervariables have been used, find permutation for
C variables from supervariable permutation
            IF (ICNTL(4).EQ.0) THEN
               CALL MC60DD(N,NSUP,IW(SVAR),IW(VARS),IW(PERMSV),
     +                    PERM,IW(IWORK))
            ELSE
               DO 50 I = 1,N
                  PERM(I) = IW(PERMSV+I-1)
   50          CONTINUE
            END IF

C Use MC60F/FD to compute the profile etc for the permuted matrix
            CALL MC60FD(N,NSUP,LIRN,IRN,ICPTR,IW(VARS),IW(PERMSV),
     *                 IW(IWORK),RINFO(5))
         ELSE
C Assembly order required.
C Find row assembly order from (super)variable permutation
            CALL MC60ED(N,NSUP,LIRN,IRN,ICPTR,IW(SVAR),IW(VARS),
     *                 IW(PERMSV),PERM,IW(IWORK))
C Use MC60G/GD to compute statistics for new assembly order
C (compressed graph used for this).
            CALL MC60GD(N,NSUP,LIRN,IRN,ICPTR,IW(VARS),IW(PERMSV),
     *                 IW(IWORK),RINFO(5))
         END IF

   60 CONTINUE

C If we have done 2 runs, we now have to choose the better
C set of weights ... store in RINFO(9) and RINFO(10)
       LSWAP = .FALSE.
       IF (NRUN.EQ.2) THEN
          IF (JOB.EQ.1) THEN
C choose on basis of profile
             IF (RNFO5.LT.RINFO(5)) LSWAP = .TRUE.
          ELSE IF (JOB.EQ.3) THEN
C choose on basis of mean frontal matrix size
             IF (RNFO8.LT.RINFO(8)) LSWAP = .TRUE.
          END IF
          IF (LSWAP) THEN
C First set of weights gave best result
             RINFO(9) = TWO
             RINFO(10) = ONE
             IF (ICNTL(6).EQ.1)  THEN
                RINFO(9) = ONE
                RINFO(10) = TWO
             END IF
             DO 70 I = 1,N
               PERM(I) = IW(COPY+I-1)
   70       CONTINUE
            RINFO(5) = RNFO5
            RINFO(6) = RNFO6
            RINFO(7) = RNFO7
            RINFO(8) = RNFO8
         END IF
      END IF

C Issue a warning if new ordering is not better than old.
          IF (JOB.EQ.1 .AND. RINFO(1).LE.RINFO(5)) THEN
C choose on basis of profile
             IF (MP.GT.0) WRITE (MP,'(/,A)')
     *       ' MC61A/AD warning: Profile not reduced'
             INFO(1) = 1
          ELSE IF (JOB.EQ.2 .AND. RINFO(3).LE.RINFO(7)) THEN
C choose on basis of bandwidth
             IF (MP.GT.0) WRITE (MP,'(/,A)')
     *       ' MC61A/AD warning: Semibandwidth not reduced'
             INFO(1) = 1
          ELSE IF (JOB.EQ.3 .AND. RINFO(4).LE.RINFO(8)) THEN
C choose on basis of mean frontal size
             IF (MP.GT.0) WRITE (MP,'(/,A)')
     * ' MC61A/AD warning: mean frontal matrix size not reduced'
             INFO(1) = 1
          END IF

      GO TO 200
C Error returns

  100 INFO(1) = -1
      IF (LP.GT.0) THEN
         WRITE (LP,9000) INFO(1)
         WRITE (LP,9010) JOB
      END IF
      GO TO 200

  110 INFO(1) = -1
      IF (LP.GT.0) THEN
         WRITE (LP,9000) INFO(1)
         WRITE (LP,9020) N
      END IF
      GO TO 200

  115 INFO(1) = -1
      IF (LP.GT.0) THEN
         WRITE (LP,9000) INFO(1)
         WRITE (LP,9025)
      END IF
      GO TO 200

  120 INFO(1) = -2
      IF (LP.GT.0) THEN
         WRITE (LP,9000) INFO(1)
         WRITE (LP,9030) INFO(6)
      END IF
      GO TO 200

  130 INFO(1) = -3
      IF (LP.GT.0) THEN
         WRITE (LP,9000) INFO(1)
         WRITE (LP,9040) INFO(3)
      END IF
      GO TO 200

  140 INFO(1) = -4
      IF (LP.GT.0) THEN
         WRITE (LP,9000) INFO(1)
         IF (INFO(4).GT.0) WRITE (LP,'(A,I6,A)')
     *         ' There are ',INFO(4),' out-of-range entries'
         IF (INFO(5).GT.0) WRITE (LP,'(A,I6,A)')
     *         ' There are ',INFO(5),' duplicated entries'
      END IF
      GO TO 200

  200 RETURN

9000  FORMAT (/' MC61A/AD error:  INFO(1) = ',I3)
9010  FORMAT (' Value of JOB out-of-range.  JOB = ',I8)
9020  FORMAT (' Value of N out-of-range.  N = ',I8)
9025  FORMAT (' Value of LIRN is less than ICPTR(N+1)-1')
9030  FORMAT (' Value of LIRN too small.  Increase to at least = ',I8)
9040  FORMAT (' Value of LIW too small.  Sufficient value = ',I8)
      END


C COPYRIGHT (c) 2002 CCLRC Council for the Central Laboratory
*                    of the Research Councils
C Original date 26 February 2002
C 15 October 2002. ICNTL(7) termination added.

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC67ID(ICNTL)

      INTEGER ICNTL(10)

C ICNTL(1)  : stream number of error messages (suppressed if <0)
C ICNTL(2)  : stream number of warning messages (suppressed if <0)
C ICNTL(3)  : controls whether checks are carried out on user data
C             0  : no checks on data.
C             >0 : checks on data. If out-of-range or duplicates
C             found, they are removed.
C             <0 : checks on data. If out-of-range or duplicates
C             found, computation terminates.
C If ICNTL(3) is none zero, a check is also made that the
C pattern of the matrix is symmetric. The computation terminates if it
C is not SYMMETRIC.
C ICNTL(4)  : maximum number of times exchange algorithm applied
C             (1 is used if ICNTL(4) is negative or zero).
C             Default 5.
C ICNTL(5)  : controls whether down/up exchanges, up/down, only down
C             or only up.
C             0 = down/up (default)
C             1 = up/down
C             2 = down
C             3 = up
C (all values except 1,2,3 perform down/up exchanges)
C ICNTL(6)  : controls whether profile for given PERM computed
C             0 = not computed (default); otherwise computed
C ICNTL(7)  : controls termination. Terminates after a profile
C reduction that is not greater than ICNTL(7)% of the first reduction

      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 0
      ICNTL(4) = 5
      ICNTL(5) = 0
      ICNTL(6) = 0
      ICNTL(7) = 0
      ICNTL(8) = 0
      ICNTL(9) = 0
      ICNTL(10) = 0

      END


      SUBROUTINE MC67AD(N,LJCN,JCN,ROWPTR,PERM,LIW,IW,ICNTL,INFO,RINFO)

C Calling routine for Hager down/up exchanges.
C The matrix must have a symmetric pattern. The diagonal need NOT
C be present.
C The user's data is optionally checked for errors.

      INTEGER N
C Order of matrix is N. N.GE.1
      INTEGER LJCN
C Length of JCN is LJCN. LJCN.GE.ROWPTR(N+1)-1
      INTEGER JCN(LJCN)
C JCN (IN) holds column entries of WHOLE matrix, ordered by rows.
C The diagonal need not be present.
C JCN only changed if ICNTL(3) > 0 and duplicates/out-of-range entries
C found.
      INTEGER ROWPTR(N+1)
C ROWPTR (IN) holds row pointers for JCN
C ROWPTR only changed if ICNTL(3) > 0 and duplicates/out-of-range entrie
C found.
      INTEGER PERM(N)
C PERM (INOUT) used to hold permutation array. Must be set
C     on entry to initial permutation. PERM(I) is the
C     row in position I of the reordered matrix.
      INTEGER LIW
C LIW is length of work array IW
      INTEGER IW(LIW)
      INTEGER ICNTL(10)
C ICNTL holds control parameters
      INTEGER INFO(10)
C INFO(1) is used as an error flag
C INFO(2) holds number of times exchange algorithm applied
C         (at most ICNTL(4))
C INFO(3) holds number of out of range entries found
C         (ICNTL(3) nonzero only)
C INFO(4) holds number of duplicate entries found
C         (ICNTL(3) nonzero only)
      DOUBLE PRECISION RINFO(10)
C RINFO(1) holds the total profile reduction achieved.
C RINFO(2) holds profile reduction achieved through Hager down.
C RINFO(3) holds profile reduction achieved through Hager up.
C RINFO(4) holds the profile for initial PERM (ICNTL(6) nonzero only).
C RINFO(5) holds the profile for final PERM (ICNTL(6) nonzero only).

      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)

C Variables for splitting up workspace
      INTEGER DIAG,FIRST,GAINS,INVPM,NEXT,NUMFT,NUMSEC,STEPS,SECOND,
     +        START
C Print parameters
      INTEGER LP,MP
C Max. number of times exchange algorithm applied
      INTEGER MAXSWP,NPASS
C Accumulated gain, threshold for terminating
      DOUBLE PRECISION ACGAIN,THRESH
      INTEGER I,IOUT,IREP,II,I1,I2,J,KZ,SWEEP

      EXTERNAL MC67BD,MC67CD,MC67DD

C Initialise
      DO 10 I = 1,10
         INFO(I) = 0
         RINFO(I) = ZERO
  10  CONTINUE
      LP = ICNTL(1)
      MP = ICNTL(2)

C Basic error checks
      IF (N.LT.1) THEN
         INFO(1) = -1
         IF (LP.GE.0) WRITE (LP,'(A)')
     +     ' MC67A/AD error INFO(1) = -1 : N is out of range'
         RETURN
      ELSE IF (LJCN.LT.ROWPTR(N+1)-1) THEN
         INFO(1) = -2
         IF (LP.GE.0) WRITE (LP,'(A)')
     + ' MC67A/AD error INFO(1) = -2 : LJCN is too small'
         RETURN
      END IF

      MAXSWP = ICNTL(4)
      IF (MAXSWP.LE.0) MAXSWP = 1

C Requires workspace length 7N+2
C Check IW
      IF (LIW.LT.7*N+2) THEN
C Error
         INFO(1) = -3
         IF (LP.GE.0) WRITE (LP,'(A,I10)')
     +     ' MC67A/AD error INFO(1) = -3 : LIW must be at least ',7*N+2
         RETURN
      END IF

C Divide up workspace
      DIAG = 1
      NUMFT = DIAG + N
      INVPM = NUMFT + N
      SECOND = INVPM + N
      NUMSEC = SECOND + N
      START = NUMSEC + N+1
      NEXT = START + N
C NEXT is length N

      STEPS = INVPM + N
      FIRST = STEPS + N+1
      GAINS = FIRST + N
C GAINS is length N+1

C Initialise inverse permutation
      DO 20 I = 1,N
        IW(INVPM+I-1) = 0
   20 CONTINUE
      DO 25 I = 1,N
        J = PERM(I)
        IW(INVPM+J-1) = I
   25 CONTINUE
C Check we do have a permutation
      DO 30 I = 1,N
        IF (IW(INVPM+I-1).EQ.0) THEN
           INFO(1) = -4
           IF (LP.GE.0) WRITE (LP,'(A)')
     + ' MC67A/AD error INFO(1) = -4 : PERM does not hold a permutation'
           RETURN
        END IF
   30 CONTINUE
C Initialise DIAG (no need if Hager up algorithm not used)
      IF (ICNTL(5).NE.2) THEN
         DO 40 I = 1,N
            IW(DIAG+I-1) = 0
   40    CONTINUE
      END IF

C The user's data is optionally checked for errors
C Look for any repeated entries and entries with out-of-range
C indices. Remove such entries and issue a warning if ICNTL(3) > 0.
C IW(I) will be set to J if column I is encountered in row J.

      IF (ICNTL(3).NE.0) THEN
        DO 50 I = 1,N
           IW(NUMFT+I-1) = 0
   50    CONTINUE
         IOUT = 0
         IREP = 0
      END IF

      IF (ICNTL(3).GT.0) THEN
         KZ = 0
         I1 = ROWPTR(1)
         ROWPTR(1) = 1
         DO 70 J = 1,N
            DO 60 II = I1,ROWPTR(J+1) - 1
               I = JCN(II)
               IF (I.GT.N .OR. I.LT.1) THEN
C Out of range
                  IOUT = IOUT + 1
               ELSE IF (IW(NUMFT+I-1).EQ.J) THEN
C Duplicate
                  IREP = IREP + 1
               ELSE
                  IF (I.EQ.J) IW(DIAG+J-1) = 1
                  KZ = KZ + 1
                  JCN(KZ) = I
                  IW(NUMFT+I-1) = J
               END IF
   60       CONTINUE
            I1 = ROWPTR(J+1)
            ROWPTR(J+1) = KZ + 1
   70    CONTINUE

         IF (IOUT.GT.0) THEN
             INFO(1) = INFO(1) + 1
             IF (MP.GE.0) WRITE (MP,'(A,I6,A)')
     *         ' MC67A/AD warning:',IOUT,' out-of-range entries removed'
         END IF
         IF (IREP.GT.0) THEN
             IF (IOUT.EQ.0) INFO(1) = INFO(1) + 1
             IF (MP.GE.0) WRITE (MP,'(A,I6,A)')
     *         ' MC67A/AD warning:',IREP,' duplicated entries removed'
         END IF

         INFO(3) = IOUT
         INFO(4) = IREP

      ELSE IF (ICNTL(3).LT.0) THEN
C Any duplicates or out of range entries will terminate the computation
         DO 100 J = 1,N
            I1 = ROWPTR(J)
            I2 = ROWPTR(J+1) - 1
            DO 90 II = I1,I2
               I = JCN(II)
               IF (I.GT.N .OR. I.LT.1) THEN
                  IOUT = IOUT + 1
               ELSE IF (IW(NUMFT+I-1).EQ.J) THEN
                  IREP = IREP + 1
               ELSE
                  IF (I.EQ.J) IW(DIAG+J-1) = 1
                  IW(NUMFT+I-1) = J
               END IF
   90       CONTINUE
  100    CONTINUE

         IF (IOUT.GT.0 .OR. IREP.GT.0) THEN
            INFO(1) = -6
            IF (LP.GE.0) THEN
              WRITE (LP,'(/,A,I3)')
     *            ' MC67A/AD error:  INFO(1) =', INFO(1)
            IF (IOUT.GT.0) WRITE (LP,'(A,I6,A)')
     *         ' There are ',IOUT,' out-of-range entries'
            IF (IREP.GT.0) WRITE (LP,'(A,I6,A)')
     *         ' There are ',IREP,' duplicated entries'
            END IF
            INFO(3) = IOUT
            INFO(4) = IREP
            RETURN
         END IF

      ELSE IF (ICNTL(5).NE.2) THEN
C No checks on arrays; set DIAG as Hager up is to be called
         DO 115 J = 1,N
            I1 = ROWPTR(J)
            I2 = ROWPTR(J+1) - 1
            DO 110 II = I1,I2
               I = JCN(II)
               IF (I.EQ.J) THEN
                  IW(DIAG+J-1) = 1
                  GO TO 115
               END IF
  110       CONTINUE
  115    CONTINUE

      END IF

C Check symmetric pattern (we just check that row counts = column counts
C which is cheap but not a comprehensive check that matrix is symmetric)
      DO 120 I = 1,N
        IW(NUMFT+I-1) = 0
  120 CONTINUE
      DO 130 J = 1,N
         I1 = ROWPTR(J)
         I2 = ROWPTR(J+1) - 1
         DO 125 II = I1,I2
           I = JCN(II)
           IW(NUMFT+I-1) = IW(NUMFT+I-1) + 1
  125    CONTINUE
  130 CONTINUE
      DO 135 I = 1,N
        IF (IW(NUMFT+I-1).NE.ROWPTR(I+1) - ROWPTR(I)) THEN
C Col. count not equal to row count ... error return
           INFO(1) = -5
           IF (LP.GE.0) THEN
              WRITE (LP,'(/,A,I3)')
     +      ' MC67A/AD error:  INFO(1) =', INFO(1)
              WRITE (LP,'(A)')
     +      ' Matrix pattern not symmetric'
           END IF
           RETURN
        END IF
  135 CONTINUE

C Optionally compute profile for given permutation
      IF (ICNTL(6).NE.0)
     +   CALL MC67DD(N,LJCN,JCN,ROWPTR,PERM,IW(INVPM),RINFO(4))

      NPASS = 1
      IF (ICNTL(5).EQ.1) THEN
C Up/down
         DO 150 SWEEP = 1,MAXSWP
            ACGAIN = RINFO(1)
            CALL MC67CD(N,LJCN,JCN,ROWPTR,IW(DIAG),IW(GAINS),IW(NUMFT),
     &                  PERM,IW(INVPM),IW(STEPS),IW(FIRST),RINFO(3),
     &                  ICNTL,NPASS)
            CALL MC67BD(N,LJCN,JCN,ROWPTR,IW(SECOND),IW(NUMSEC),
     &                  IW(NUMFT),PERM,IW(INVPM),IW(START),IW(NEXT),
     &                  RINFO(2),ICNTL,NPASS)
            RINFO(1) = RINFO(2) + RINFO(3)
            IF (SWEEP.EQ.1) THRESH = MAX(ZERO,ICNTL(7)*RINFO(1)/100)
C Jump to return if improvement is too small
            IF (RINFO(1)-ACGAIN.LE.THRESH) THEN
C Set INFO(2) to hold number of sweeps performed
               INFO(2) = SWEEP
               GO TO 200
            END IF
  150    CONTINUE
      ELSE IF (ICNTL(5).EQ.2) THEN
C Down only
         CALL MC67BD(N,LJCN,JCN,ROWPTR,IW(SECOND),IW(NUMSEC),IW(NUMFT),
     &               PERM,IW(INVPM),IW(START),IW(NEXT),RINFO(2),
     &               ICNTL,MAXSWP)
      ELSE IF (ICNTL(5).EQ.3) THEN
C Up only
         CALL MC67CD(N,LJCN,JCN,ROWPTR,IW(DIAG),IW(GAINS),IW(NUMFT),
     &               PERM,IW(INVPM),IW(STEPS),IW(FIRST),
     &               RINFO(3),ICNTL,MAXSWP)
      ELSE
C Down/up (default)
         DO 180 SWEEP = 1,MAXSWP
            ACGAIN = RINFO(1)
            CALL MC67BD(N,LJCN,JCN,ROWPTR,IW(SECOND),IW(NUMSEC),
     &                  IW(NUMFT),PERM,IW(INVPM),IW(START),IW(NEXT),
     &                  RINFO(2),ICNTL,NPASS)
            CALL MC67CD(N,LJCN,JCN,ROWPTR,IW(DIAG),IW(GAINS),IW(NUMFT),
     &                  PERM,IW(INVPM),IW(STEPS),IW(FIRST),RINFO(3),
     &                  ICNTL,NPASS)
            RINFO(1) = RINFO(2) + RINFO(3)
            IF (SWEEP.EQ.1) THRESH = MAX(ZERO,ICNTL(7)*RINFO(1)/100)
C Jump to return if improvement is too small
            IF (RINFO(1)-ACGAIN.LE.THRESH) THEN
C Set INFO(2) to hold number of sweeps performed
               INFO(2) = SWEEP
               GO TO 200
            END IF
  180    CONTINUE
      END IF
C Set INFO(2) to hold number of sweeps performed
      INFO(2) = MAXSWP

  200 CONTINUE

      RINFO(1) = RINFO(2) + RINFO(3)

      IF (ICNTL(6).NE.0) RINFO(5) = RINFO(4) - RINFO(1)

      END

      SUBROUTINE MC67BD(N,LJCN,JCN,ROWPTR,SECOND,NUMSEC,NUMFT,
     &                  PERM,INVPM,START,NEXT,TGAIN,ICNTL,NPASS)

C This subroutine performs Hager downward exchanges.
C Row and column indices are for the original matrix, unless stated
C otherwise.
      INTEGER N
C N (IN) holds the order of the matrix.
      INTEGER LJCN
C LJCN (IN) holds the length of JCN.
      INTEGER JCN(LJCN)
C JCN (IN) holds column entries of WHOLE matrix, ordered by rows,
C     with no duplicates in a row. The diagonal is treated as if
c     present even if it is not.
      INTEGER ROWPTR(N+1)
C ROWPTR (IN) holds the positions of the row starts in JCN.
      INTEGER SECOND(N)
C SECOND is a work array. SECOND(J) is used to hold the permuted index
C     of the second entry in original col. J.
      INTEGER NUMSEC(N+1)
C NUMSEC is a work array. NUMSEC(I) is used to hold the number of
C     columns with first entry in permuted row K and second entry in
C     permuted row I.
      INTEGER NUMFT(N)
C NUMFT is a work array. NUMFT(I) is used to hold the number of first
C     entries in row I.
      INTEGER PERM(N)
C PERM (INOUT) used to hold permutation array. Must be set on entry
C     to initial permutation. PERM(I) is the row in position I of
C     the reordered matrix.
      INTEGER INVPM(N)
C INVPM (INOUT) holds the inverse permutation. Must be set on entry.
      INTEGER START(N)
C START is a work array. START(I) holds the original column index of
C     a first entry in original row I. All first entries in each row
C     are linked through NEXT.
      INTEGER NEXT(N)
C NEXT is a work array. NEXT(J) holds the original column index of
C     the next first entry in the same row.
      DOUBLE PRECISION TGAIN
C TGAIN (INOUT) must be set on entry to hold the total profile gain
C     so far achieved. The extra gain found is added.
      INTEGER ICNTL(10)
C ICNTL holds control parameters
      INTEGER NPASS
C NPASS (INOUT). On entry, must be set to a limit on the number of
C      sweeps. On exit, holds the number of times the matrix was swept.

      DOUBLE PRECISION ACGAIN,THRESH,ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER GAIN,I,II,IORG,J,K,KK,KORG,L,M,MORG
      INTEGER MXGAIN,NTGAIN,NXT,PREV,SEC,SWAPS,SWEEP
C ACGAIN Accumulated gain.
C GAIN  Gain in a single step.
C I     Row index.
C II    Temporary do index.
C IORG  Original index of row I.
C J     Column index.
C K     Row at start of cyclic permutation.
C KK    Temporary do index.
C KORG  Original index of row K.
C L     Row at end of cyclic permutation.
C M     Permuted row index
C MORG  Original index of row M.
C MXGAIN Maximum gain so far in this cyclic permutation.
C NTGAIN Net gain so far in this cyclic permutation
C NXT   Next entry in list for row K.
C PREV  Previous entry in list for row K.
C SEC   Permuted index of second entry of column.
C SWAPS Number of swaps in the sweep.
C SWEEP Do index for the complete sweeps.
C THRESH Threshold for terminating.

C Initialise
      DO 10 I = 1,N
        NUMFT(I) = 0
        NUMSEC(I) = 0
        NEXT(I) = 0
   10 CONTINUE
      NUMSEC(N+1) = 0

C Set the array NEXT to hold original indices of the
C first entries in the columns.
      DO 30 I = N,1,-1
        IORG = INVPM(I)
        NEXT(IORG) = IORG
        DO 20 II = ROWPTR(IORG),ROWPTR(IORG+1)-1
          J = JCN(II)
          NEXT(J) = IORG
   20   CONTINUE
   30 CONTINUE

C Link and count the first nodes in each row
      DO 40 J = 1,N
        I = NEXT(J)
        IF (I.GT.0) THEN
          NUMFT(I) = NUMFT(I) + 1
C Add J to its linked list
          IF (NUMFT(I).GT.1) NEXT(J) = START(I)
          START(I) = J
        END IF
   40 CONTINUE

C Continue sweeps until there is no gain
      DO 210 SWEEP = 1, NPASS
      SWAPS = 0
      ACGAIN = TGAIN
C Sweep from the bottom to top of the permuted matrix
C looking for a sequence of rows to permute cyclically.
      DO 200 K = N-1,1,-1
C Row K was row INVPM(K) in the original matrix.
        KORG = INVPM(K)
        GAIN =  NUMFT(KORG)
C Check if any gain is possible
        IF (GAIN.EQ.0) GO TO 200

C Loop over first entries in row K to find second entries in the cols
        J = START(KORG)
        DO 60 KK = 1, GAIN
          SEC = PERM(J)
          IF (SEC.EQ.K) SEC = N+1
          DO 50 II = ROWPTR(J),ROWPTR(J+1)-1
            I = PERM(JCN(II))
            IF (I.NE.K) SEC = MIN(SEC,I)
   50     CONTINUE
          SECOND(J) = SEC
          NUMSEC(SEC) = NUMSEC(SEC) + 1
          J = NEXT(J)
   60   CONTINUE

        NTGAIN = 0
        MXGAIN = 0
C Sweep through rows below K, looking for best swap
        DO 80 M = K+1,N
C Loop over first entries in row M (originally row INVPM(M))
          MORG = INVPM(M)
          GAIN = GAIN - NUMSEC(M)
C If GAIN = 0 then no need to examine any more rows as
C we cannot improve on MXGAIN
          IF (GAIN.EQ.0) GO TO 90
C All first entries in this row will lead to a loss of one
          NTGAIN = NTGAIN + GAIN - NUMFT(MORG)
C It is just possible that the above statement causes integer
C overflow, but this is very unlikely if the matrix has already
C been ordered sensibly. If integer overflow happens, we will
C still choose a very good permutation and later sweeps will
C pick up the gains that are missed this time.
          IF (NTGAIN.GT.MXGAIN) THEN
C M is best row found so far to swap with K
            MXGAIN = NTGAIN
            L = M
          END IF
   80   CONTINUE
C Restore NUMSEC to zero
   90   J = START(KORG)
        DO 100 KK = 1, NUMFT(KORG)
          SEC = SECOND(J)
          NUMSEC(SEC) = 0
          J = NEXT(J)
  100   CONTINUE

        IF (MXGAIN.EQ.0) GO TO 200
        SWAPS = SWAPS + L - K
C We have to update the linked list of first entries
        J = START(KORG)
        PREV = 0
        GAIN =  NUMFT(KORG)
        DO 160 KK = 1, GAIN
          NXT = NEXT(J)
          M = SECOND(J)
          IF (M.LE.L) THEN
C M has becomes first entry in column J instead of K.
             MORG = INVPM(M)
             NUMFT(KORG) = NUMFT(KORG) - 1
             NUMFT(MORG) = NUMFT(MORG) + 1
             IF (NUMFT(KORG).GT.0) THEN
C Remove J from its linked list
                IF (PREV.EQ.0) THEN
                   START(KORG) = NXT
                ELSE
                   NEXT(PREV) = NXT
                END IF
             END IF
C Add J to its new linked list
             IF (NUMFT(MORG).GT.1) NEXT(J) = START(MORG)
             START(MORG) = J
          ELSE
             PREV = J
          END IF
          J = NXT
 160    CONTINUE

C Update PERM and INVPM. We are swapping rows between K and L. Row K
C will move down to L, L moves up to L-1, L-1 to L-2 ,...
        DO 170 I = K,L-1
          IORG = INVPM(I+1)
          INVPM(I) = IORG
          PERM(IORG) = I
  170   CONTINUE
        INVPM(L) = KORG
        PERM(KORG) = L

        TGAIN = TGAIN + MXGAIN

  200 CONTINUE
      IF (SWEEP.EQ.1) THRESH = MAX(ZERO,ICNTL(7)*(TGAIN-ACGAIN)/100)
C Jump to return if improvement is too small
      IF (TGAIN-ACGAIN.LE.THRESH) THEN
         NPASS = SWEEP
         GO TO 220
      END IF
  210 CONTINUE

  220 RETURN

      END

      SUBROUTINE MC67CD(N,LJCN,JCN,ROWPTR,DIAG,GAINS,NUMFT,
     &                  PERM,INVPM,STEPS,FIRST,TGAIN,ICNTL,NPASS)

C This subroutine performs Hager upward exchanges.
C Row and column indices are for the original matrix, unless stated
C otherwise.
      INTEGER N
C N (IN) holds the order of the matrix.
      INTEGER LJCN
C LJCN (IN) holds the length of JCN.
      INTEGER JCN(LJCN)
C JCN (IN) holds column entries of WHOLE matrix, ordered by rows,
C     with no duplicates in a row. The diagonal is treated as if
C     present even if it is not.
      INTEGER ROWPTR(N+1)
C ROWPTR (IN) holds row pointers for JCN
      INTEGER DIAG(N)
C DIAG (IN). DIAG(I) has value 1 if row I includes a diagonal entry
C     and value 0 otherwise.
      INTEGER FIRST(N)
C FIRST is a work array. FIRST(J) is used to hold the row in which the
C     first entry in column J lies.
      INTEGER NUMFT(N)
C NUMFT is a work array. NUMFT(I) is used to hold number of first
C     entries in row I.
      INTEGER PERM(N)
C PERM (INOUT) used to hold permutation array. Must be set on entry
C     to initial permutation. PERM(I) is the row in position I of
C     the reordered matrix.
      INTEGER INVPM(N)
C INVPM (INOUT) holds the inverse permutation. Must be set on entry.
      INTEGER STEPS(0:N)
C STEPS is a work array. A row with a first entry in a column with an
C      entry in row K is called a 'step' row. STEPS holds the permuted
C      indices of the step rows, with duplicates if a row has more
C      than one such entry.
      INTEGER GAINS(N+1)
C GAINS is a work array. GAINS holds accumulated gains.
      DOUBLE PRECISION TGAIN
C TGAIN (INOUT) must be set on entry to hold the total profile gain
C     so far achieved. The extra gain found is added.
      INTEGER ICNTL(10)
C ICNTL holds control parameters
      INTEGER NPASS
C NPASS (INOUT). On entry, must be set to a limit on the number of
C      sweeps. On exit, holds the number of times the matrix was swept.

      DOUBLE PRECISION ACGAIN,THRESH,ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER I,II,IORG,J,K,KORG,L,LEN,M,NEXT,LOSS,LOSS1
      INTEGER MXGAIN,NTGAIN,NTGNQ,P,Q,QQ,SWAPS,SWEEP
C ACGAIN Accumulated gain.
C I     Row index.
C II    Temporary do index.
C IORG  Original index of row I.
C J     Column index.
C K     Row at start of cyclic permutation.
C KORG  Original index of row K.
C L     Row at end of cyclic permutation.
C LEN   Number of entres in row K (also number of step rows since
C       duplicates are included).
C M     Permuted row index
C NEXT  Next row to be considered when searching an interval.
C LOSS  Value of loss per row.
C LOSS1 Loss for first row.
C MXGAIN Maximum gain so far in this cyclic permutation.
C NTGAIN Net gain so far in this cyclic permutation
C NTGNQ Net gain at row Q
C NXT   Next entry in list for row K.
C P,Q   Search is between rows P-1 and Q.
C QQ    Upper limit for simple search.
C SWAPS Number of swaps in the sweep.
C SWEEP Do index for the complete sweeps.
C THRESH Threshold for terminating.

      EXTERNAL KB06AI

C Initialise
      DO 10 I = 1,N
        NUMFT(I) = 0
   10 CONTINUE

C Set the array FIRST to hold original indices of the
C first entries in the columns.
      DO 30 I = N,1,-1
        IORG = INVPM(I)
        FIRST(IORG) = IORG
        DO 20 II = ROWPTR(IORG),ROWPTR(IORG+1)-1
          J = JCN(II)
          FIRST(J) = IORG
   20   CONTINUE
   30 CONTINUE

C Count the first nodes in each row
      DO 40 J = 1,N
        I = FIRST(J)
        NUMFT(I) = NUMFT(I) + 1
   40 CONTINUE

C Accumulate potential gains
      GAINS(N+1) = 0
      DO 50 I = N,1,-1
        GAINS(I) = GAINS(I+1) + NUMFT(INVPM(I))
   50 CONTINUE

C Continue sweeps until there is no gain
      DO 210 SWEEP = 1, NPASS
      SWAPS = 0
      ACGAIN = TGAIN
C Sweep from the top to bottom of the matrix looking for rows to swap.
      DO 200 K = 2,N

C Row K was row INVPM(K) in the original matrix.
        KORG = INVPM(K)

C Set STEPS for first rows of all columns J with an entry in row K
        STEPS(0) = 0
        LEN = 0
        IF (DIAG(KORG).EQ.0) THEN
          LEN = 1
          STEPS(1) = PERM(FIRST(KORG))
        END IF
        DO 60 II = ROWPTR(KORG), ROWPTR(KORG+1) - 1
           J = FIRST(JCN(II))
           LEN = LEN + 1
           STEPS(LEN) = PERM(J)
   60   CONTINUE
C Sort STEPS (keeping duplicates).
        CALL KB06AI(STEPS,LEN+1)

        LOSS1 = NUMFT(KORG)
        IF (LOSS1.EQ.0) THEN
          Q = STEPS(0) + 1
          LOSS1 = 1
        ELSE
          Q = K
        END IF
        L = Q
        NTGNQ = GAINS(Q) - GAINS(K)
        MXGAIN = NTGNQ
        DO 80 LOSS = LOSS1, LEN
          P = Q
          Q = STEPS(LOSS) + 1
          QQ = Q
          NTGAIN = NTGNQ
          NTGNQ = NTGAIN+GAINS(Q)-GAINS(P)-LOSS*(P-Q)
C Outer loop.  Search the interval P-1 to Q.
C NTGAIN is the net gain at row P. The loss per row is LOSS.
C The gain from a single cyclic permutation cannot exceed
C N-1, so integer overflows cannot occur.

C Search exploiting the linearity of loss
          DO 66 II = 1,N
C We are sure that rows QQ+1 to Q, inclusive, are not advantageous.
            NEXT = P - (NTGAIN+GAINS(QQ)-GAINS(P)-MXGAIN)/LOSS
            IF (NEXT.GE.P) GO TO 80
            IF (NEXT.LE.QQ) GO TO 67
            QQ = NEXT
   66     CONTINUE

C Simple search between P-1 and QQ
   67     DO 70 M = P-1,QQ,-1
            NTGAIN = NTGAIN + NUMFT(INVPM(M)) - LOSS
            IF (NTGAIN.GT.MXGAIN) THEN
C Hold this as best swap so far found for row K
              L = M
              MXGAIN = NTGAIN
            END IF
            IF (NTGAIN+GAINS(QQ)-GAINS(M)-LOSS.LE.MXGAIN) GO TO 80
   70     CONTINUE

   80   CONTINUE

C If MXGAIN = 0, no row found to swap with K
        IF (MXGAIN.EQ.0) GO TO 200

        SWAPS = SWAPS + K - L

C We have to update FIRST
        IORG = FIRST(KORG)
        I = PERM(IORG)
        IF (I.GE.L) THEN
C The first entry that was in row I moves up to row K
           NUMFT(IORG) = NUMFT(IORG) - 1
           FIRST(KORG) = KORG
           NUMFT(KORG) = NUMFT(KORG) + 1
        END IF
        DO 160 II = ROWPTR(KORG),ROWPTR(KORG+1)-1
          J = JCN(II)
          IORG = FIRST(J)
          I = PERM(IORG)
          IF (I.GE.L) THEN
C The first entry that was in row I moves up to row K
             NUMFT(IORG) = NUMFT(IORG) - 1
             FIRST(J) = KORG
             NUMFT(KORG) = NUMFT(KORG) + 1
          END IF
 160    CONTINUE

C We are swapping rows K and L ie K will move up to L,
C       L moves down to L+1, L+1 to L+2 ,...
        DO 150 I = K,L+1,-1
          IORG = INVPM(I-1)
          INVPM(I) = IORG
          PERM(IORG) = I
          GAINS(I) = GAINS(I+1) + NUMFT(IORG)
  150   CONTINUE
        INVPM(L) = KORG
        PERM(KORG) = L

        TGAIN = TGAIN + MXGAIN

  200 CONTINUE
      IF (SWEEP.EQ.1) THRESH = MAX(ZERO,ICNTL(7)*(TGAIN-ACGAIN)/100)
C Jump to return if improvement is too small
      IF (TGAIN-ACGAIN.LE.THRESH) THEN
        NPASS = SWEEP
        GO TO 220
      END IF
  210 CONTINUE

  220 RETURN

      END

      SUBROUTINE MC67DD(N,LJCN,JCN,ROWPTR,PERM,INVPM,FRONT)
C This subroutine computes the profile of a given permutation
         INTEGER N,LJCN
         INTEGER JCN(LJCN),ROWPTR(N+1),PERM(N),INVPM(N)
         DOUBLE PRECISION FRONT
         INTEGER I,K,KORG,KMIN
         FRONT = N
         DO 20 K = 1,N
           KORG = INVPM(K)
           KMIN = K
           DO 10 I = ROWPTR(KORG),ROWPTR(KORG+1)-1
              KMIN = MIN(KMIN,PERM(JCN(I)))
   10      CONTINUE
           FRONT = FRONT + (K - KMIN)
   20   CONTINUE

      END
