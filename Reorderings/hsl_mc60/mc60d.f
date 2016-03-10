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
