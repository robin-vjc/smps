!-----------------------------------------------------------------------
!                     Copyright notice
!
! SMPSReader is copyright 1986-2002 by Horand I. Gassmann.
! All rights reserved.
!
! SMPSReader is protected by copyright.  It is not public-domain
! software or shareware, and it is not protected by a ``copyleft''
! agreement like the one used by the Free Software Foundation.
!
! SMPSReader is available free for any use in any field of endeavor.
! You may redistribute SMPSReader in whole or in part provided you
! acknowledge its source and include this copyright notice.  You may
! modify SMPSReader and create derived works, provided you retain this
! copyright notice, but the result may not be called SMPSReader without
! my written consent.
!
! You may sell a derived work, provided that all source code for your
! derived work is available, at no additional charge, to anyone who buys
! your derived work in any form.  You must give permisson for said
! source code to be used and modified under the terms of this license.
! You must state clearly that your work uses or is based on SMPSReader
! and that SMPSReader is available free of change.
!
!-----------------------------------------------------------------------
!
!                            No warranty
!
! This software is distributed in the hopes that it will be useful, but
! WITHOUT ANY WARRANTY; not even the implied warranty of MERCHANTIBILITY
! or FITNESS FOR A PARTICULAR PURPOSE.
!
!-----------------------------------------------------------------------
!
!     Please report bugs to hgassman@mgmt.dal.ca.
!
!                 Development history
!
!     Last change:  HG   7 Dec 2002    4:38 pm
!
!-----------------------------------------------------------------------
!
!                 Introductory remarks
!
! SMPSReader processes the input data for a stochastic linear program.
! It uses the SMPS format developed by Birge et al. (1987) and extended
! by Gassmann and Schweitzer (2001). A description of the format is also
! available on the world-wibe web at
! http://www.mgmt.da.ca/sba/profs/hgassmann/SMPS2.htm.
!
! The most general program that can be expressed in the SMPS format
! has the following mathematical form:
!
!                               -                      -
!                              |  Q_11  Q_12 ...  Q_1T  |
! Opt! 1/2 (x_1, x_2,...,x_T)' |  Q_21  Q_22 ...  Q_2T  |(x_1, x_2,...,x_T)
!                              |                        |
!                              |  Q_T1  Q_T2 ...  Q_TT  |
!                               -                      -
!
!      +   (c_1, c_2,...,c_T)' (x_1, x_2,...,x_T)
!
! s.t.
!        A_11 x_1 + A_12 x_2 + ... + A_1T x_T <> b_1
!        A_21 x_1 + A_22 x_2 + ... + A_2T x_T <> b_2
!           ...        ...              ...      ...
!        A_T1 x_1 + A_T2 x_2 + ... + A_TT x_T <> b_T
!
!               l_i <= x_i <= u_i,  i=1,...,T
!
!               x_i in R^{n^1_i} x Z^{n^2_i}.
!
! All the data items except A_11, Q_11, c_1, b_1, l_1 and u_1 can be
! stochastic. The symbol `<>' stands for an arbitrary relation (<=, =, >=).
!
! The constraint may be required to hold for every possible realization,
! or with probability one, or subject to a probabilistic constraint.
! Constraint matrices A_ij with i>j define so-called global (or linking)
! constraints. If A_ij = 0 for j > i+1, the problem is said to possess
! `staircase structure'.
!
! The input is contained in three files, called the time file, core file
! and stoch file. Each file is aranged into several sections, each containing
! a header record and zero or more data records in MPS format. Briefly,
! the MPS format for data files provides a four-character code field, three
! name fields and two numerical fields. The order of the sections is fixed,
! although some sections are optional and may be omitted.
!
! The purpose of the time file is to allow breaking down the core file
! scenario into nodes corresponding to the individual stages. If the
! core file is given in temporal order, this can be accomplished simply
! by placing markers for the first row and column of each stage.
! In addition each stage is given a name in the time file by which it
! can subsequently be referred.
!
! The core file lists all the deterministic information of the problem,
! in standard MPS format. Core file information includes the name and
! type of each constraint and variable, a column-ordered list of
! coefficients in the constraint matrix, right-hand side coefficients and
! any bounds and ranges on the variables and slacks. The core file also
! provides placeholders for all the stochastic elements (which must be
! mentioned in the core file and given a preliminary value
! that may or may not be meaningful). The core file thus represents a
! deterministic problem, which may be thought of as a typical scenario,
! an average case, a baseline case or similar.
!
! The stoch file allows the solver to build a deterministic
! equivalent, which requires information about the random variables.
! If all the random variables are finitely distributed, this task amounts
! to the construction of an event tree. The event tree can be described
! in three different ways: scenario by scenario, node by node, or
! implicitly using marginal distributions. Implicit descriptions are
! also available if the random variables are continuously distributed,
! where explicit descriptions would obviously fail.
!
! The internal representation of the data is very much node-oriented.
! Each node in the event tree contains a linear (or quadratic) subproblem
! made up of rows, columns and variables (columns plus logicals for the rows). All data items are stored in allocatable arrays;
! associated with every array are three pointers, for the currently
! allocated size, for the maximal allowable size and for the number
! of elements currently in use. This allows dynamic resizing of all arrays.
!
!---------------------------------------------------------------------

      MODULE STRING_LENGTH
      IMPLICIT NONE
!
!  If free form input is used, variable names can be longer than 8 characters.
!  NAMLEN designates the maximal length of a name. JCCLEN designates
!  the maximal length of the name of a multidimensional chance constraint.
!
      INTEGER NAMLEN, JCCLEN
      PARAMETER ( NAMLEN = 255 )
      PARAMETER ( JCCLEN =  31 )
!
      END MODULE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      MODULE IO_HANDLING
      USE STRING_LENGTH
      IMPLICIT NONE
!
!     File names and I/O channel assignments
!
      CHARACTER*(NAMLEN) TIMFIL,CORFIL,STOFIL,INBFIL,PARFIL,LOGFIL,      &
     &                   BASFIL,SUMFIL,SOLFIL,EVPFIL,SCRFIL,QUAFIL,      &
     &                   MPSFIL
      INTEGER  IOTIM, IOCOR, IOSTO, IOINB, IOPAR, IOLOG, IOBAS,      &
     &         IOSUM, IOSOL, IOEVP, IOSCR, IOQUA, IOMPS
!
      CONTAINS
!
      SUBROUTINE IOPREP(MODE)
!
!     THIS SUBROUTINE CONTAINS THE DEFAULT I/O UNIT ASSIGNMENTS AND
!     OPEN/CLOSE STATEMENTS SHOULD SUCH BE NECESSARY.
!
!          ***   PARAMETER DESCRIPTION   ***
!
!         MODE = 0   to set default
!         MODE = 1   for setup/open I/O channels
!         MODE = 2   shutdown/close I/O channels
!
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
!     ----------------------------------------------------------------
!
!     THE FOLLOWING I/O CHANNELS ARE USED:
!
!        IOTIM - Time file
!        IOCOR - Core file
!        IOQUA - File for QSECTIONs
!        IOSTO - Stoch file
!        IOINB - Input starting basis
!        IOPAR - Parameters
!        IOLOG - Detailed iteration log
!        IOBAS - Output final basis
!        IOSUM - Summary results
!        IOSOL - Optimal solution
!        IOEVP - Output file for EVPI report
!        IOSCR - Scratch file for specs file processing
!
!     ----------------------------------------------------------------
!
      IF (MODE .EQ. 1) THEN
         IOTIM = 1
         IOCOR = 2
         IOQUA = 2
         IOSTO = 3
         IOINB = 4
         IOPAR = 5
         IOLOG = 6
         IOBAS = 7
         IOSUM = 8
         IOSOL = 9
         IOEVP = 13
         IOSCR = 14
         IOMPS = 23
!
         TIMFIL = 'FORT.1'
         CORFIL = 'FORT.2'
         QUAFIL = 'FORT.2'
         STOFIL = 'FORT.3'
         INBFIL = 'FORT.4'
         PARFIL = 'FORT.5'
         LOGFIL = 'MSLiPlog.txt'
         BASFIL = 'FORT.7'
         SUMFIL = 'FORT.8'
         SOLFIL = 'FORT.9'
         EVPFIL = 'FORT.13'
         SCRFIL = 'FORT.14'
         MPSFIL = 'OUTPUT.MPS'
      ENDIF
      RETURN
      END SUBROUTINE
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CHOPEN(FILNAM,LUNIT,LSTAT,LFORM,IOERR)
!
!     This subroutine is used to open a file.
!
!     File names could be read from another file, or they could be
!     placed in the spec file (but the current record structure does
!     not allow for that). Alternatively, a dialog could be started
!     to inquire after the name the user wants. The final option
!     is to create a default file name based on the unit number,
!     in the way unix and VMS do.
!
!     FILNAM    is the file name. This can be a maximum of 255 characters
!               long, including the path and the extension.
!
!     LUNIT     is the logical unit number
!
!     LSTAT     = 1  means old;
!               = 2  means new;
!               = 3  means unknown;
!               = 4  means scratch.
!
!     LFORM     = 1  means unformatted;
!               = 2  means formatted.
!
!     IOERR     is output, corresponding to the error status
!
      CHARACTER*255 FILNAM
      CHARACTER*7  DSTAT(4)
      CHARACTER*11 DFORM(2)
      INTEGER LFORM,NFORM,LSTAT,NSTAT,LUNIT,IOERR
!
      DATA DSTAT/'OLD    ','NEW    ','UNKNOWN','SCRATCH'/
      DATA DFORM/'UNFORMATTED','FORMATTED  '/
!
      IF (FILNAM .EQ. '?') THEN
         WRITE (FILNAM, 1000) LUNIT
      ENDIF
!
      NSTAT = LSTAT
      NFORM = LFORM
      IF (LSTAT .EQ. 0) NSTAT = 3
      IF (LFORM .EQ. 0) NFORM = 2
!
      OPEN( UNIT=LUNIT, FILE=FILNAM, STATUS=DSTAT(NSTAT),    &
     &      FORM=DFORM(NFORM),IOSTAT=IOERR)
      IF (IOERR .GT. 0) THEN
         WRITE (*, 1200) FILNAM, LUNIT, IOERR
      ENDIF
!
      RETURN
!
 1000 FORMAT('fort.',I3.3)
 1200 FORMAT(' Error opening file ',A30,' on unit',I3,'.',       &
     & ' Error status = ',I3,'.')
      END SUBROUTINE
      END MODULE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      MODULE TIMER_STAT
!
!     Variables for timer
!
      INTEGER  FIRSTT
      LOGICAL  TIMESET
!
      CONTAINS
!
      SUBROUTINE STIMER ( TIME )
!
! -------------------------------------------------------------------------
!
!     This timer routine provides an interface between the system clock
!     on an Intel PC running under Windows and the calling program.
!     Identical interfaces should be provided on other platforms. This way,
!     only the interface has to be exchanged when going from one operating
!     system to the next.
!
!     This clock allows for one possible rollover. For Intel systems this
!     means that it is useful to measure time intervals up to 24 hours.
!
! -------------------------------------------------------------------------
!
!     *****     Parameters     *****
!
!     TIME  -   CPU seconds used by the program since the start of execution.
!
      REAL*8    TIME
      INTEGER   ITIM, ITICKS, MTICKS
!
      CALL SYSTEM_CLOCK(ITIM, ITICKS, MTICKS)
      IF (TIMESET) THEN
         TIME = (ITIM-FIRSTT)*(1.D0/ITICKS)
         IF (TIME .LT. 0.D0) TIME = TIME + (1.D0*MTICKS)/ITICKS
      ELSE
         TIMESET = .TRUE.
         FIRSTT  = ITIM
      ENDIF
!
      RETURN
!
      END SUBROUTINE
      END MODULE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      MODULE SMPS_DATA
      USE STRING_LENGTH
!
!     This is the module that sets up all the storage arrays needed to parse
!  the SMPS file, of which there are a great many. The program assumes that
!  the problem comes in the form of an event tree with a linear (or quadratic)
!  subproblem at each node, consisting of rows, columns and variables
!  (columns plus logicals for the rows). In addition there are blocks
!  associated with the constraint matrices and quadratic objectives.
!
!  Associated with every array are three pointers, for the currently
!  allocated size, for the maximal allowable size and for the number
!  of elements currently in use.
!
      IMPLICIT NONE
!
!  Each NODE contains the following parameters and variables:
!
      TYPE NODE
           REAL (KIND=8) :: PROB   ! Probability
           INTEGER       :: IANCTR ! Ancestor node
           INTEGER       :: IDESC  ! First descendant
           INTEGER       :: IBROTH ! Next sibling
           INTEGER       :: NDESC  ! # of descendants
           INTEGER       :: KCOL   ! Column pointer
           INTEGER       :: KROW   ! Row pointer
           INTEGER       :: KCOST  ! Cost pointer
           INTEGER       :: KRHS   ! RHS pointer
           INTEGER       :: KBOUND ! Bound pointer
           INTEGER       :: KNAMES ! Name pointer
           INTEGER       :: VARPTR ! Var pointer
           INTEGER       :: NCOL   ! # of columns
           INTEGER       :: NROW   ! # of rows
      END  TYPE NODE
!
      TYPE (NODE),   ALLOCATABLE, DIMENSION (:) :: NODEINFO
!
!  Pointers:
      INTEGER  NODES, IZNODE, MXNODE
!
!
! -----------------------------------------------------------------------
!     Basic constraint matrix elements.
!
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: A      ! A-matrix elems
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: IA     ! A-matrix rows
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: LA     ! A-matrix cols
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: KDATA  ! A-elem pointer
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: KCOLA  ! A-mtx pointer
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: KELMA  ! A-col pointer
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: NELMA  ! # A-elements
      LOGICAL,       ALLOCATABLE, DIMENSION (:,:) :: STOCHA
!
!     Usage pointers:
      INTEGER LASTA,  IZALMN, MXALMN   ! for A and IA
      INTEGER LASTCA, IZACOL, MXACOL   ! for LA
      INTEGER LASTBA, IZABLK, MXABLK   ! for KCOLA, KELMA, NELMA
!  There are no special pointers for KDATA (which is indexed over the event
!  nodes). The shape of the matrix is further described by MARKOV and HAVE_GC).
      LOGICAL MARKOV, HAVE_GC
!
!
! -----------------------------------------------------------------------
!  Distinct coefficients for the right-hand side:
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: RHS
      INTEGER LASTR, IZDRHS, MXDRHS
!
!
! -----------------------------------------------------------------------
!  Distinct cost coefficients:
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: COST
      INTEGER LASTC, IZCOST, MXCOST
!
!
! -----------------------------------------------------------------------
!  Distinct upper and lower bounds:
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: XLB    ! Lower bounds
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: XUB    ! Upper bounds
      INTEGER LASTBD, IZBNDS, MXBNDS
!
!
! -----------------------------------------------------------------------
!     Arrays for variable names and types
!
      INTEGER,     ALLOCATABLE, DIMENSION (:) :: NAME1
      CHARACTER*1, ALLOCATABLE, DIMENSION (:) :: QNAMES
      INTEGER,     ALLOCATABLE, DIMENSION (:) :: SOSTYP
      INTEGER,     ALLOCATABLE, DIMENSION (:) :: VARTYP
!
!  Pointers:
      INTEGER LASTNM, IZVNAM, MXVNAM   ! for NAME1
      INTEGER LASTCH, IZCHAR, MXCHAR   ! for QNAMES
      INTEGER NOFSOS, IZSTYP, MXSTYP   ! for SOSTYP
      INTEGER LASTVP, IZVTYP, MXVTYP   ! for VARTYP

!
! -----------------------------------------------------------------------
!     Arrays for Q-matrix are similar to the A-matrix arrays:
!
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: AQMTX  ! Q-matrix elems
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: IQMTX  ! Q-matrix rows
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: LQMTX  ! Q-matrix cols
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: KDATQ  ! Q-block pointer
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: IQOFF  ! Q-row pointer
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: LQOFF  ! Q-col pointer
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: NELMQ  ! # Q-elements
!
!     Pointers:
      INTEGER LASTQ,  IZQLMN, MXQLMN   ! for AQMTX and IQMTX
      INTEGER LASTQC, IZQCOL, MXQCOL   ! for LQMTX
      INTEGER LASTQB, IZQBLK, MXQBLK   ! for IQOFF, LQOFF, NELMQ
!  KDATQ is indexed over the event nodes, as its counterpart in the A-matrix.
      INTEGER NQPROF        ! describes the shape or profile of the Q-matrix
!
! -----------------------------------------------------------------------
!     Arrays for stochastic elements (D-matrix)
!
      REAL (KIND=8),      ALLOCATABLE, DIMENSION (:) :: DMTX    ! D-elements
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: IDMTX   ! D-rows
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: LDMTX   ! D-cols
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: MDIST   ! Distr.params
      CHARACTER*(NAMLEN), ALLOCATABLE, DIMENSION (:) :: SDIST   ! Distr.names
      REAL (KIND=8),      ALLOCATABLE, DIMENSION (:) :: RANDPAR ! Random params (REAL*8)
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: IRNDPAR ! Integer params
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: LRNDPAR ! Location params
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: KSTOCH  ! Stoch. elements
      CHARACTER*(NAMLEN), ALLOCATABLE, DIMENSION (:) :: ALIAS   ! Alternative names
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: IALIAS  ! Assoc. blocks/components
!
!     Pointers:
      INTEGER LASTD,  IZDLMN,  MXDLMN  ! for DMTX, IDMTX
      INTEGER LASTDC, IZDCOL,  MXDCOL  ! for LDMTX
      INTEGER NRVAR,  IZRAND,  MXRAND  ! for MDIST, SDIST
      INTEGER LRPAR,  IZRPAR,  MXRPAR  ! for RANDPAR
      INTEGER LIPAR,  IZIPAR,  MXIPAR  ! for IRNDPAR
      INTEGER LLPAR,  IZLPAR,  MXLPAR  ! for LRNDPAR
      INTEGER NSTELM, IZSTOC,  MXSTOC  ! for KSTOCH
      INTEGER NALIAS, IZALIAS, MXALIAS ! for ALIAS, IALIAS
!
!
! -----------------------------------------------------------------------
!     Arrays for individual and joint chance constraints
!
      TYPE CHANCE
           INTEGER       :: CTYPE  ! constraint type (+1 for >, -1 for <, etc.)
           REAL (KIND=8) :: PROB   ! required probability level
           INTEGER       :: SIZE   ! number of rows in this constraint
           INTEGER       :: NODE   ! node associated with this constraint
           INTEGER       :: ROWPTR ! pointer into RCHANCE for participating rows
           INTEGER       :: NAMPTR ! pointer into QCHANCE for multidimensional constraints
      END  TYPE CHANCE
!
      TYPE CHANCE_ROW
           INTEGER       :: STAGE  ! stage to which the row belongs
           INTEGER       :: RELROW ! relative row number
      END  TYPE CHANCE_ROW
!
      TYPE (CHANCE),      ALLOCATABLE, DIMENSION (:) :: XCHANCE ! constraint ptr
      TYPE (CHANCE_ROW),  ALLOCATABLE, DIMENSION (:) :: RCHANCE ! constraint row
      CHARACTER*(JCCLEN), ALLOCATABLE, DIMENSION (:) :: QCHANCE ! name array
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: NCHANCE ! number in each node
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: KCHANCE ! pointer into XCHANCE
!
!     Pointers:
      INTEGER LCHANCE, IZCCON, MXCCON  ! for XCHANCE
      INTEGER          IZCCRO, MXCCRO  ! for RCHANCE
      INTEGER NJOINT,  IZCCNM, MXCCNM  ! for QCHANCE
!
!     NCHANCE and KCHANCE are indexed over the time stages, so there is
!     no need for separate pointers.
!
!
! -----------------------------------------------------------------------
!     Arrays for integrated chance constraints
!
      TYPE ICC_DATA
           INTEGER       :: CTYPE  ! constraint type (+1 for >, -1 for <, etc.)
           REAL (KIND=8) :: LIMIT  ! limit on CVaR (RHS of ICC)
           INTEGER       :: NODE   ! node to which the row belongs
           INTEGER       :: RELROW ! relative row number
      END  TYPE ICC_DATA
!
      TYPE (ICC_DATA), ALLOCATABLE, DIMENSION (:) :: XICC ! ICC data
      INTEGER,         ALLOCATABLE, DIMENSION (:) :: NICC ! number in each node
      INTEGER,         ALLOCATABLE, DIMENSION (:) :: KICC ! pointer into XICC
      INTEGER LICC, IZNICC, MXNICC
!
!
! -----------------------------------------------------------------------
!     Arrays for piecewise linear-quadratic penalties
!
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: CPLQ    ! quadratic penalty
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: KPLQ    ! pointer to CPLQ
      INTEGER LASTPLQ, IZCPLQ, MXCPLQ  ! for CPLQ; KPLQ indexed over nodes
!
! -----------------------------------------------------------------------
!     Arrays for scenario names and leaf nodes
!
      INTEGER,     ALLOCATABLE, DIMENSION (:) :: LEAFND
      INTEGER,     ALLOCATABLE, DIMENSION (:) :: KSCNAM
      CHARACTER*1, ALLOCATABLE, DIMENSION (:) :: QSCNAM
      INTEGER LASTLF, IZSCEN, MXSCEN  ! for LEAFND and KSCNAM
      INTEGER NSCHAR, IZCHSC, MXCHSC  ! for QSCNAM
!
!
! -----------------------------------------------------------------------
!     Arrays of strings
!
      CHARACTER*(NAMLEN), ALLOCATABLE, DIMENSION (:) :: DTIME
      CHARACTER*(NAMLEN), ALLOCATABLE, DIMENSION (:) :: DTIMEC
      CHARACTER*(NAMLEN), ALLOCATABLE, DIMENSION (:) :: DTIMER
!
!
! =======================================================================
!
!     Logical variables for input
!
      LOGICAL QMAT, FREEFORM, PURELP, RESIZE, ALLOW_GC,      &
     &        HAVE_PLQP, CHNG_BLKS
!
!     QMAT     - used in RDSPEC to allow quadratic objective (checked in input routines)
!     FREEFORM - used to denote free/fixed format in the SMPS files
!     PURELP   - used to distinguish between pure LP and network problems
!     RESIZE   - used to allow dynamic resizing of arrays
!     ALLOW_GC - used to allow global constraints (and suppress some diagnostics)
!     HAVE_PLQP- signals the presence of piecewise linear-quadratic penalties
!     CHNG_BLKS- determines whether A-matrix coefficients stand alone or are
!                to be added to the core file values
!
!
!     Character variables for several problem components can be set
!     at runtime to distinguish between different problems in the same file
!
      CHARACTER*(NAMLEN) DRHS,  DBOUND,DRANGE,PROBNM,DCHANCE,      &
     &                   OBJNAM,DSIM,  DROB,  DPLQ,  DICC
!
!
!     Other array sizes (for temporary arrays)
!
      INTEGER IZTPER,  MXTPER  ! Number of time periods
      INTEGER IZROWP,  MXROWP  ! Number of rows per node
      INTEGER IZCOLP,  MXCOLP  ! Number of columns per node
      INTEGER IZROWS,  MXROWS  ! Number of rows altogether
      INTEGER IZCOLS,  MXCOLS  ! Number of columns altogether
      INTEGER IZANZB,  MXANZB  ! Number of nonzeroes in an A-matrix block
      INTEGER IZQNZB,  MXQNZB  ! Number of nonzeroes in an Q-matrix block
      INTEGER IZRTYP,  MXRTYP
      INTEGER IZCOPY,  MXCOPY
!
!
!     Input parameters
!
      REAL*8   PLINF,  ZTOLIN, DEFUB, DEFLB
      INTEGER  NECHO1
!
!     Input statistics
!
      INTEGER  TREEFMT, STFILE,  TIMSTAT, CORSTAT, STOSTAT, SLVCND,      &
     &         MAXROW,  MAXCOL,  MAXRHS,  LENOBJ,  IOBJ1,   NPSEEN,      &
     &         NPER,    IOBJ,    INODE,   IPER,    NROWS,   NCOLS,      &
     &         NGLOBC,  NDELEM,  NDCOLM
      INTEGER, ALLOCATABLE, DIMENSION (:) :: IRNGE0
!
!     Character variables
!
      CHARACTER*1  QA,QAST,QB,QBL,QC,QD,QE,QF,QG,QH,QI,QJ,QK,      &
     &             QL,QM,QN,QO,QP,QQ,QR,QS,QT,QU,QV,QX,QSTAT
!
      DATA QBL /' '/, QA  /'A'/, QB  /'B'/, QC  /'C'/, QD  /'D'/,      &
     &     QE  /'E'/, QF  /'F'/, QG  /'G'/, QH  /'H'/, QI  /'I'/,      &
     &     QJ  /'J'/, QK  /'K'/, QL  /'L'/, QM  /'M'/, QN  /'N'/,      &
     &     QO  /'O'/, QP  /'P'/, QQ  /'Q'/, QR  /'R'/, QS  /'S'/,      &
     &     QT  /'T'/, QU  /'U'/, QV  /'V'/, QX  /'X'/, QAST/'*'/
!
      CHARACTER*8 DOTS
      DATA DOTS /'  ...   '/
!
      END MODULE
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      MODULE SAMPLING_DATA
      INTEGER ISAMPL,ISMPRB,ISEED,IRAND,ISEED1,ISEED2(2),ISEED5(3),      &
     &        ISEED6,I97,J97,IJKL,I24,J24
      LOGICAL INITRG
      REAL (KIND=8) :: SEED3(97),SEED4(24),C,CD,CM,CARRY
      END MODULE
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      MODULE TEMP_DATA
      USE STRING_LENGTH
      IMPLICIT NONE
!
!     This module uses several temporary arrays. Set them up first.
!
      TYPE ROWINFO
           CHARACTER*(NAMLEN) NAME
           INTEGER            LENGTH
           INTEGER            ROWTYP
           REAL (KIND=8)      UPPER
           REAL (KIND=8)      LOWER
      END TYPE ROWINFO
!
      TYPE (ROWINFO), ALLOCATABLE, DIMENSION (:  )   :: T_ROWS
      TYPE (ROWINFO), ALLOCATABLE, DIMENSION (:,:)   :: T_VAR
!
      TYPE TIMENAME
            INTEGER PERIOD
            CHARACTER*(NAMLEN) NAME
      END TYPE TIMENAME
!
      TYPE (TIMENAME), ALLOCATABLE, DIMENSION (:)     :: TROWS
      TYPE (TIMENAME), ALLOCATABLE, DIMENSION (:)     :: TCOLS
      INTEGER IZTROW, MXTROW
      INTEGER IZTCOL, MXTCOL
!
      INTEGER,         ALLOCATABLE, DIMENSION (:  )   :: LMNS
      INTEGER,         ALLOCATABLE, DIMENSION (:,:)   :: IAUX
      INTEGER,         ALLOCATABLE, DIMENSION (:,:)   :: LAUX
      REAL (KIND=8),   ALLOCATABLE, DIMENSION (:,:)   :: AUX
!
      INTEGER,         ALLOCATABLE, DIMENSION (:,:  ) :: LMNS2
      INTEGER,         ALLOCATABLE, DIMENSION (:,:,:) :: IAUX2
      INTEGER,         ALLOCATABLE, DIMENSION (:,:,:) :: LAUX2
      REAL (KIND=8),   ALLOCATABLE, DIMENSION (:,:,:) :: AUX2
!
      INTEGER,         ALLOCATABLE, DIMENSION (:)     :: TCNROW
      REAL (KIND=8),   ALLOCATABLE, DIMENSION (:,:)   :: TCOST
!
      INTEGER,         ALLOCATABLE, DIMENSION (:  )   :: QLINKS
      INTEGER,         ALLOCATABLE, DIMENSION (:  )   :: QRTMP
      INTEGER,         ALLOCATABLE, DIMENSION (:  )   :: QCTMP
      REAL (KIND=8),   ALLOCATABLE, DIMENSION (:  )   :: QMTMP
      INTEGER,         ALLOCATABLE, DIMENSION (:,:)   :: ISTART
      integer, allocatable, dimension (:) :: NDUP, NBLK
!
!     These temporary arrays are used in INTREE
!
      CHARACTER*(NAMLEN), ALLOCATABLE, DIMENSION (:) :: DNODE
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: INUSE, KPATH
!
!     These arrays are used in MKNODE
!
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: ALINKS
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: ARTMP
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: ACTMP
      REAL (KIND=8),      ALLOCATABLE, DIMENSION (:) :: AMTMP
!
!     These arrays are used in INCONT
!
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: KREF
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: MARKER
      CHARACTER*(NAMLEN), ALLOCATABLE, DIMENSION (:) :: COPIED
!
!     These arrays are used in INSCEN
!
      CHARACTER*(NAMLEN), ALLOCATABLE, DIMENSION (:) :: DSCNAM
!
      END MODULE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      MODULE MEMORY_EXTENDER
!
      CONTAINS
!
      LOGICAL FUNCTION EXTMEM_ABLK(MNABLK,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER, ALLOCATABLE, DIMENSION (:) :: IT1,IT2,IT3
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNABLK .GT. MXABLK)      &
     &                 .OR. (IZABLK .GE. MXABLK)) GOTO 100
         ALLOCATE(IT1(IZABLK), IT2(IZABLK), IT3(IZABLK), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = KCOLA
         IT2 = KELMA
         IT3 = NELMA
         DEALLOCATE (KCOLA, KELMA, NELMA, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2ABLK = MAX(MNABLK, MIN(MXABLK,2*IZABLK))
         ALLOCATE( KCOLA(I2ABLK), KELMA(I2ABLK),      &
     &                            NELMA(I2ABLK), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         KCOLA(1:IZABLK) = IT1
         KELMA(1:IZABLK) = IT2
         NELMA(1:IZABLK) = IT3
         DEALLOCATE (IT1, IT2, IT3, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZABLK = I2ABLK
         EXTMEM_ABLK = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_ABLK = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_ACOL(MNACOL,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER, ALLOCATABLE, DIMENSION (:) :: IT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNACOL .GT. MXACOL)      &
     &                 .OR. (IZACOL .GE. MXACOL)) GOTO 100
         ALLOCATE(IT1(IZACOL), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = LA
         DEALLOCATE (LA, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2ACOL = MAX(MNACOL, MIN(MXACOL,2*IZACOL))
         ALLOCATE( LA(I2ACOL), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         LA(1:IZACOL) = IT1
         DEALLOCATE (IT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZACOL = I2ACOL
         EXTMEM_ACOL = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_ACOL = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_ALIAS(MNALIAS,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: IT1
      CHARACTER*(NAMLEN), ALLOCATABLE, DIMENSION (:) :: DT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNALIAS .GT. MXALIAS)      &
     &                 .OR. (IZALIAS .GE. MXALIAS)) GOTO 100
         ALLOCATE(IT1(IZALIAS), DT1(IZALIAS), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = IALIAS
         DT1 =  ALIAS
         DEALLOCATE (ALIAS, IALIAS, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2ALIAS = MAX(MNALIAS, MIN(MXALIAS,2*IZALIAS))
         ALLOCATE(IALIAS(I2ALIAS), ALIAS(I2ALIAS),  STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IALIAS(1:IZALIAS) = IT1
          ALIAS(1:IZALIAS) = DT1
         DEALLOCATE (IT1, DT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZALIAS = I2ALIAS
         EXTMEM_ALIAS = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_ALIAS = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_ALMN(MNALMN,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: IT1
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: RT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNALMN .GT. MXALMN)      &
     &                 .OR. (IZALMN .GE. MXALMN)) GOTO 100
         ALLOCATE(IT1(IZALMN), RT1(IZALMN), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = IA
         RT1 =  A
         DEALLOCATE (A, IA, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2ALMN = MAX(MNALMN, MIN(MXALMN,2*IZALMN))
         ALLOCATE(IA(I2ALMN), A(I2ALMN),  STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IA(1:IZALMN) = IT1
          A(1:IZALMN) = RT1
         DEALLOCATE (IT1, RT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZALMN = I2ALMN
         EXTMEM_ALMN = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_ALMN = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_ANZB(MNANZB,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER,       ALLOCATABLE, DIMENSION (:)     :: IT1, IT2, IT3
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:)     :: RT1
      INTEGER,       ALLOCATABLE, DIMENSION (:,:)   :: JT1
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:,:)   :: ST1
      INTEGER,       ALLOCATABLE, DIMENSION (:,:,:) :: KT1
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:,:,:) :: TT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNANZB .GT. MXANZB)      &
     &                 .OR. (IZANZB .GE. MXANZB)) GOTO 100
         IF (ALLOCATED(AUX)) THEN
            ALLOCATE(JT1(IZANZB,NPER), ST1(IZANZB,NPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            DO IP=1,NPER
               JT1(1:IZANZB,IP) = IAUX(1:IZANZB,IP)
               ST1(1:IZANZB,IP) =  AUX(1:IZANZB,IP)
            END DO
            DEALLOCATE (AUX, IAUX, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            I2ANZB = MAX(MNANZB, MIN(MXANZB,2*IZANZB))
            ALLOCATE(AUX(I2ANZB,NPER), IAUX(I2ANZB,NPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            DO IP=1,NPER
               IAUX(1:IZANZB,IP) = JT1(1:IZANZB,IP)
                AUX(1:IZANZB,IP) = ST1(1:IZANZB,IP)
            END DO
            DEALLOCATE (JT1, ST1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
!
         IF (ALLOCATED(AUX2)) THEN
            ALLOCATE(KT1(IZANZB,NPER,NPER),      &
     &               TT1(IZANZB,NPER,NPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            DO IP=1,NPER
               DO JP=1,NPER
                  KT1(1:IZANZB,IP,JP) = IAUX2(1:IZANZB,IP,JP)
                  TT1(1:IZANZB,IP,JP) =  AUX2(1:IZANZB,IP,JP)
               END DO
            END DO
            DEALLOCATE (AUX2, IAUX2, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            I2ANZB = MAX(MNANZB, MIN(MXANZB,2*IZANZB))
            ALLOCATE(AUX2(I2ANZB,NPER,NPER),      &
     &              IAUX2(I2ANZB,NPER,NPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            DO IP=1,NPER
               DO JP=1,NPER
                  IAUX2(1:IZANZB,IP,JP) = KT1(1:IZANZB,IP,JP)
                   AUX2(1:IZANZB,IP,JP) = TT1(1:IZANZB,IP,JP)
               END DO
            END DO
            DEALLOCATE (KT1, TT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
         IF (ALLOCATED(ALINKS)) THEN
            ALLOCATE(IT1(IZANZB), IT2(IZANZB), IT3(IZANZB),      &
     &                                         RT1(IZANZB), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            DEALLOCATE (ALINKS, ARTMP, ACTMP, AMTMP, STAT=IERR)
            IT1 = ALINKS
            IT2 = ARTMP
            IT3 = ACTMP
            RT1 = AMTMP
            IF (IERR .NE. 0) GOTO 100
            I2ANZB = MAX(MNANZB, MIN(MXANZB,2*IZANZB))
            ALLOCATE(ALINKS(I2ANZB), ARTMP(I2ANZB), ACTMP(I2ANZB),      &
     &                               AMTMP(I2ANZB), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            ALINKS(1:IZANZB) = IT1
            ARTMP (1:IZANZB) = IT2
            ACTMP (1:IZANZB) = IT3
            AMTMP (1:IZANZB) = RT1
            DEALLOCATE (IT1, IT2, IT3, RT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            ENDIF
         IZANZB = I2ANZB
         EXTMEM_ANZB = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_ANZB = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_BNDS(MNBNDS,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: RT1,RT2
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNBNDS .GT. MXBNDS)      &
     &                 .OR. (IZBNDS .GE. MXBNDS)) GOTO 100
         ALLOCATE(RT1(IZBNDS), RT2(IZBNDS), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         RT1 = XLB
         RT2 = XUB
         DEALLOCATE (XLB, XUB, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2BNDS = MAX(MNBNDS, MIN(MXBNDS,2*IZBNDS))
         ALLOCATE(XLB(I2BNDS), XUB(I2BNDS), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         XLB(1:IZBNDS) = RT1
         XUB(1:IZBNDS) = RT2
         DEALLOCATE (RT1, RT2, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZBNDS = I2BNDS
         EXTMEM_BNDS = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_BNDS = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_CCNM(MNCCNM,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(JCCLEN), ALLOCATABLE, DIMENSION (:) :: QT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNCCNM .GT. MXCCNM)      &
     &                 .OR. (IZCCNM .GE. MXCCNM)) GOTO 100
         ALLOCATE(QT1(IZCCNM), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         QT1 = QCHANCE
         DEALLOCATE (QCHANCE, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2CCNM = MAX(MNCCNM, MIN(MXCCNM,2*IZCCNM))
         ALLOCATE( QCHANCE(I2CCNM), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         QCHANCE(1:IZCCNM) = QT1
         DEALLOCATE (QT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZCCNM = I2CCNM
         EXTMEM_CCNM = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_CCNM = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_CCON(MNCCON,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      TYPE (CHANCE), ALLOCATABLE, DIMENSION (:) :: LT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNCCON .GT. MXCCON)      &
     &                 .OR. (IZCCON .GE. MXCCON)) GOTO 100
         ALLOCATE(LT1(IZCCON), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         LT1 = XCHANCE
         DEALLOCATE (XCHANCE, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2CCON = MAX(MNCCON, MIN(MXCCON,2*IZCCON))
         ALLOCATE( XCHANCE(I2CCON), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         XCHANCE(1:IZCCON) = LT1
         DEALLOCATE (LT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZCCON = I2CCON
         EXTMEM_CCON = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_CCON = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_CCRO(MNCCRO,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      TYPE (CHANCE_ROW), ALLOCATABLE, DIMENSION (:) :: LT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNCCRO .GT. MXCCRO)      &
     &                 .OR. (IZCCRO .GE. MXCCRO)) GOTO 100
         ALLOCATE(LT1(IZCCRO), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         LT1 = RCHANCE
         DEALLOCATE (RCHANCE, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2CCRO = MAX(MNCCRO, MIN(MXCCRO,2*IZCCRO))
         ALLOCATE( RCHANCE(I2CCRO), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         RCHANCE(1:IZCCRO) = LT1
         DEALLOCATE (LT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZCCRO = I2CCRO
         EXTMEM_CCRO = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_CCRO = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_CHAR(MNCHAR,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*1, ALLOCATABLE, DIMENSION (:) :: QT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNCHAR .GT. MXCHAR)      &
     &                 .OR. (IZCHAR .GE. MXCHAR)) GOTO 100
         ALLOCATE(QT1(IZCHAR), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         QT1 = QNAMES
         DEALLOCATE (QNAMES, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2CHAR = MAX(MNCHAR, MIN(MXCHAR,2*IZCHAR))
         ALLOCATE(QNAMES(I2CHAR), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         QNAMES(1:IZCHAR) = QT1
         DEALLOCATE (QT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZCHAR = I2CHAR
         EXTMEM_CHAR = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_CHAR = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_CHSC(MNCHSC,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*1, ALLOCATABLE, DIMENSION (:) :: QT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNCHSC .GT. MXCHSC)      &
     &                 .OR. (IZCHSC .GE. MXCHSC)) GOTO 100
         ALLOCATE(QT1(IZCHSC), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         QT1 = QSCNAM
         DEALLOCATE (QSCNAM, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2CHSC = MAX(MNCHSC, MIN(MXCHSC,2*IZCHSC))
         ALLOCATE(QSCNAM(I2CHSC), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         QSCNAM(1:IZCHSC) = QT1
         DEALLOCATE (QT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZCHSC = I2CHSC
         EXTMEM_CHSC = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_CHSC = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_COLP(MNCOLP,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER,        ALLOCATABLE, DIMENSION (:,:)   :: JT1
      INTEGER,        ALLOCATABLE, DIMENSION (:,:,:) :: KT1
      REAL (KIND=8),  ALLOCATABLE, DIMENSION (:,:)   :: RT1
      TYPE (ROWINFO), ALLOCATABLE, DIMENSION (:,:)   :: TT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNCOLP .GT. MXCOLP)      &
     &                 .OR. (IZCOLP .GE. MXCOLP)) GOTO 100
         I2COLP = MAX(MNCOLP, MIN(MXCOLP,2*IZCOLP))
         IF (ALLOCATED (LAUX)) THEN
            ALLOCATE(JT1(IZCOLP,NPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            DO IP=1,NPER
               JT1(1:IZCOLP,IP) = LAUX(1:IZCOLP,IP)
            END DO
            DEALLOCATE (LAUX, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            ALLOCATE(LAUX(I2COLP,NPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            DO IP=1,NPER
               LAUX(1:IZCOLP,IP) = JT1(1:IZCOLP,IP)
            END DO
            DEALLOCATE (JT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
         IF (ALLOCATED(LAUX2)) THEN
            ALLOCATE(KT1(IZCOLP,NPER,NPER), RT1(IZCOLP,NPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            DO IP=1,NPER
               RT1(1:IZCOLP,IP) = TCOST(1:IZCOLP,IP)
               DO JP=1,NPER
                  KT1(1:IZCOLP,IP,JP) = LAUX2(1:IZCOLP,IP,JP)
               END DO
            END DO
            DEALLOCATE (LAUX2, TCOST, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            ALLOCATE(LAUX2(I2COLP,NPER,NPER),      &
     &               TCOST(I2COLP,NPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            DO IP=1,NPER
               TCOST(1:IZCOLP,IP) = RT1(1:IZCOLP,IP)
               DO JP=1,NPER
                  LAUX2(1:IZCOLP,IP,JP) = KT1(1:IZCOLP,IP,JP)
               END DO
            END DO
            DEALLOCATE (KT1, RT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
         IF (ALLOCATED(T_VAR)) THEN
            ALLOCATE(TT1(IZCOLP,NPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            DO I=1,NPER
               TT1(1:IZCOLP,I) = T_VAR(1:IZCOLP,I)
            END DO
            DEALLOCATE (T_VAR, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            ALLOCATE(T_VAR(I2COLP,NPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            DO I=1,NPER
               T_VAR(1:IZCOLP,I) = TT1(1:IZCOLP,I)
            END DO
            DEALLOCATE (TT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
         IZCOLP = I2COLP
         EXTMEM_COLP = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_COLP = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_COPY(MNCOPY,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: IT1
      CHARACTER*(NAMLEN), ALLOCATABLE, DIMENSION (:) :: DT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNCOPY .GT. MXCOPY)      &
     &                 .OR. (IZCOPY .GE. MXCOPY)) GOTO 100
         ALLOCATE(IT1(IZCOPY), DT1(IZCOPY), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = MARKER
         DT1 = COPIED
         DEALLOCATE (COPIED, MARKER, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2COPY = MAX(MNCOPY, MIN(MXCOPY,2*IZCOPY))
         ALLOCATE(MARKER(I2COPY), COPIED(I2COPY),  STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         MARKER(1:IZCOPY) = IT1
         COPIED(1:IZCOPY) = DT1
         DEALLOCATE (IT1, DT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZCOPY = I2COPY
         EXTMEM_COPY = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_COPY = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_COST(MNCOST,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: RT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNCOST .GT. MXCOST)      &
     &                 .OR. (IZCOST .GE. MXCOST)) GOTO 100
         ALLOCATE(RT1(IZCOST), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         RT1 = COST
         DEALLOCATE (COST, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2COST = MAX(MNCOST, MIN(MXCOST,2*IZCOST))
         ALLOCATE(COST(I2COST), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         COST(1:IZCOST) = RT1
         DEALLOCATE (RT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZCOST = I2COST
         EXTMEM_COST = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_COST = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_CPLQ(MNCPLQ,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: RT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNCPLQ .GT. MXCPLQ)      &
     &                 .OR. (IZCPLQ .GE. MXCPLQ)) GOTO 100
         ALLOCATE(RT1(IZCPLQ), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         RT1 = CPLQ
         DEALLOCATE (CPLQ, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2CPLQ = MAX(MNCPLQ, MIN(MXCPLQ,2*IZCPLQ))
         ALLOCATE(CPLQ(I2CPLQ), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         CPLQ(1:IZCPLQ) = RT1
         DEALLOCATE (RT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZCPLQ = I2CPLQ
         EXTMEM_CPLQ = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_CPLQ = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_DCOL(MNDCOL,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER, ALLOCATABLE, DIMENSION (:) :: IT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNDCOL .GT. MXDCOL)      &
     &                 .OR. (IZDCOL .GE. MXDCOL)) GOTO 100
         ALLOCATE(IT1(IZDCOL), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = LDMTX
         DEALLOCATE (LDMTX, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2DCOL = MAX(MNDCOL, MIN(MXDCOL,2*IZDCOL))
         ALLOCATE( LDMTX(I2DCOL), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         LDMTX(1:IZDCOL) = IT1
         DEALLOCATE (IT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZDCOL = I2DCOL
         EXTMEM_DCOL = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_DCOL = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_DLMN(MNDLMN,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: IT1
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: RT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNDLMN .GT. MXDLMN)      &
     &                 .OR. (IZDLMN .GE. MXDLMN)) GOTO 100
         ALLOCATE(IT1(IZDLMN), RT1(IZDLMN), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = IDMTX
         RT1 =  DMTX
         DEALLOCATE (DMTX, IDMTX, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2DLMN = MAX(MNDLMN, MIN(MXDLMN,2*IZDLMN))
         ALLOCATE(IDMTX(I2DLMN), DMTX(I2DLMN),  STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IDMTX(1:IZDLMN) = IT1
          DMTX(1:IZDLMN) = RT1
         DEALLOCATE (IT1, RT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZDLMN = I2DLMN
         EXTMEM_DLMN = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_DLMN = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_DRHS(MNDRHS,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: RT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNDRHS .GT. MXDRHS)      &
     &                 .OR. (IZDRHS .GE. MXDRHS)) GOTO 100
         ALLOCATE(RT1(IZDRHS), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         RT1 = RHS
         DEALLOCATE (RHS, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2DRHS = MAX(MNDRHS, MIN(MXDRHS,2*IZDRHS))
         ALLOCATE(RHS(I2DRHS), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         RHS(1:IZDRHS) = RT1
         DEALLOCATE (RT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZDRHS = I2DRHS
         EXTMEM_DRHS = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_DRHS = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_IPAR(MNIPAR,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER, ALLOCATABLE, DIMENSION (:) :: IT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNIPAR .GT. MXIPAR)      &
     &                 .OR. (IZIPAR .GE. MXIPAR)) GOTO 100
         ALLOCATE(IT1(IZIPAR), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = IRNDPAR
         DEALLOCATE (IRNDPAR, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2IPAR = MAX(MNIPAR, MIN(MXIPAR,2*IZIPAR))
         ALLOCATE( IRNDPAR(I2IPAR), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IRNDPAR(1:IZIPAR) = IT1
         DEALLOCATE (IT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZIPAR = I2IPAR
         EXTMEM_IPAR = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_IPAR = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_LPAR(MNLPAR,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER, ALLOCATABLE, DIMENSION (:) :: IT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNLPAR .GT. MXLPAR)      &
     &                 .OR. (IZLPAR .GE. MXLPAR)) GOTO 100
         ALLOCATE(IT1(3*IZLPAR), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = LRNDPAR
         DEALLOCATE (LRNDPAR, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2LPAR = MAX(MNLPAR, MIN(MXLPAR,2*IZLPAR))
         ALLOCATE( LRNDPAR(3*I2LPAR), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         LRNDPAR(1:3*IZLPAR) = IT1
         DEALLOCATE (IT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZLPAR = I2LPAR
         EXTMEM_LPAR = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_LPAR = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_NICC(MNNICC,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      TYPE (ICC_DATA), ALLOCATABLE, DIMENSION (:) :: LT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNNICC .GT. MXNICC)      &
     &                 .OR. (IZNICC .GE. MXNICC)) GOTO 100
         ALLOCATE(LT1(IZNICC), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         LT1 = XICC
         DEALLOCATE (XICC, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2NICC = MAX(MNNICC, MIN(MXNICC,2*IZNICC))
         ALLOCATE( XICC(I2NICC), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         XICC(1:IZNICC) = LT1
         DEALLOCATE (LT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZNICC = I2NICC
         EXTMEM_NICC = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_NICC = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_NODE(MNNODE,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: IT1, IT2
      TYPE (NODE),        ALLOCATABLE, DIMENSION (:) :: NT1
      CHARACTER*(NAMLEN), ALLOCATABLE, DIMENSION (:) :: DT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNNODE .GT. MXNODE)      &
     &                 .OR. (IZNODE .GE. MXNODE)) GOTO 100
         IF (ALLOCATED(DNODE)) THEN
            ALLOCATE(DT1(IZNODE), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            DT1 = DNODE
            DEALLOCATE (DNODE, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            I2NODE = MAX(MNNODE, MIN(MXNODE,2*IZNODE))
            ALLOCATE(DNODE(I2NODE), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            DNODE(1:IZNODE) = DT1
            DEALLOCATE (DT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
!
         IF (ALLOCATED(KDATA)) THEN
            ALLOCATE(IT1(IZNODE), NT1(IZNODE), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            IT1 = KDATA
            NT1 = NODEINFO
            DEALLOCATE (KDATA, NODEINFO, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            I2NODE = MAX(MNNODE, MIN(MXNODE,2*IZNODE))
            ALLOCATE(KDATA(I2NODE), NODEINFO(I2NODE), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
               KDATA(1:IZNODE) = IT1
            NODEINFO(1:IZNODE) = NT1
            DEALLOCATE (IT1, NT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
!
         IF (ALLOCATED(KDATQ)) THEN
            ALLOCATE(IT1(IZNODE), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            IT1 = KDATQ
            DEALLOCATE (KDATQ, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            ALLOCATE(KDATQ(I2NODE), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            KDATQ(1:IZNODE) = IT1
            DEALLOCATE (IT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
!
         IF (ALLOCATED(KCHANCE)) THEN
            ALLOCATE(IT1(IZNODE), IT2(IZNODE), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            IT1 = KCHANCE
            IT2 = NCHANCE
            DEALLOCATE (KCHANCE, NCHANCE, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            ALLOCATE(KCHANCE(I2NODE), NCHANCE(I2NODE), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            KCHANCE(1:IZNODE) = IT1
            NCHANCE(1:IZNODE) = IT2
            DEALLOCATE (IT1, IT2, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
!
         IF (ALLOCATED(KICC)) THEN
            ALLOCATE(IT1(IZNODE), IT2(IZNODE), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            IT1 = KICC
            IT2 = NICC
            DEALLOCATE (KICC, NICC, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            ALLOCATE(KICC(I2NODE), NICC(I2NODE), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            KICC(1:IZNODE) = IT1
            NICC(1:IZNODE) = IT2
            DEALLOCATE (IT1, IT2, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
!
         IF (ALLOCATED(KPLQ)) THEN
            ALLOCATE(IT1(IZNODE), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            IT1 = KPLQ
            DEALLOCATE (KPLQ, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            ALLOCATE(KPLQ(I2NODE), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            KPLQ(1:IZNODE) = IT1
            DEALLOCATE (IT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
!
         IZNODE = I2NODE
         EXTMEM_NODE = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_NODE = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_QBLK(MNQBLK,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: IT1, IT2, IT3
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNQBLK .GT. MXQBLK)      &
     &                 .OR. (IZQBLK .GE. MXQBLK)) GOTO 100
         ALLOCATE(IT1(IZQBLK), IT2(IZQBLK), IT3(IZQBLK), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = IQOFF
         IT2 = LQOFF
         IT3 = NELMQ
         DEALLOCATE (IQOFF, LQOFF, NELMQ, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2QBLK = MAX(MNQBLK, MIN(MXQBLK,2*IZQBLK))
         ALLOCATE(IQOFF(I2QBLK), LQOFF(I2QBLK),      &
     &                           NELMQ(I2QBLK), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IQOFF(1:IZQBLK) = IT1
         LQOFF(1:IZQBLK) = IT2
         NELMQ(1:IZQBLK) = IT3
         DEALLOCATE (IT1, IT2, IT3, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZQBLK = I2QBLK
         EXTMEM_QBLK = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_QBLK = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_QCOL(MNQCOL,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER, ALLOCATABLE, DIMENSION (:) :: IT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNQCOL .GT. MXQCOL)      &
     &                 .OR. (IZQCOL .GE. MXQCOL)) GOTO 100
         ALLOCATE(IT1(IZQCOL), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = LQMTX
         DEALLOCATE (LQMTX, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2QCOL = MAX(MNQCOL, MIN(MXQCOL,2*IZQCOL))
         ALLOCATE( LQMTX(I2QCOL), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         LQMTX(1:IZQCOL) = IT1
         DEALLOCATE (IT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZQCOL = I2QCOL
         EXTMEM_QCOL = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_QCOL = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_QLMN(MNQLMN,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: IT1, IT2, IT3
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: RT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNQLMN .GT. MXQLMN)      &
     &                 .OR. (IZQLMN .GE. MXQLMN)) GOTO 100
         IF (ALLOCATED(QLINKS)) THEN
            ALLOCATE(IT1(IZQLMN), IT2(IZQLMN), IT3(IZQLMN),      &
     &                            RT1(IZQLMN), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            IT1 = QLINKS
            IT2 = QRTMP
            IT3 = QCTMP
            RT1 = QMTMP
            DEALLOCATE (QLINKS, QRTMP, QCTMP, QMTMP, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            I2QLMN = MAX(MNQLMN, MIN(MXQLMN,2*IZQLMN))
            ALLOCATE(QLINKS(I2QLMN), QRTMP(I2QLMN), QCTMP(I2QLMN),      &
     &               QMTMP(I2QLMN), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            QLINKS(1:IZQLMN) = IT1
            QRTMP (1:IZQLMN) = IT2
            QCTMP (1:IZQLMN) = IT3
            QMTMP (1:IZQLMN) = RT1
            DEALLOCATE (IT1, IT2, IT3, RT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
         IF (ALLOCATED(AQMTX)) THEN
            ALLOCATE(IT1(IZQLMN), RT1(IZQLMN), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            IT1 = IQMTX
            RT1 = AQMTX
            DEALLOCATE (AQMTX, IQMTX, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            I2QLMN = MAX(MNQLMN, MIN(MXQLMN,2*IZQLMN))
            ALLOCATE(IQMTX(I2QLMN), AQMTX(I2QLMN),  STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            IQMTX(1:IZQLMN) = IT1
            AQMTX(1:IZQLMN) = RT1
            DEALLOCATE (IT1, RT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
         IZQLMN = I2QLMN
         EXTMEM_QLMN = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_QLMN = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_RAND(MNRAND,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: IT1
      CHARACTER*(NAMLEN), ALLOCATABLE, DIMENSION (:) :: DT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNRAND .GT. MXRAND)      &
     &                 .OR. (IZRAND .GE. MXRAND)) GOTO 100
         ALLOCATE(IT1(8*IZRAND), DT1(IZRAND), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = MDIST
         DT1 = SDIST
         DEALLOCATE (MDIST, SDIST, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2RAND = MAX(MNRAND, MIN(MXRAND,2*IZRAND))
         ALLOCATE(MDIST(8*I2RAND), SDIST(I2RAND), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         MDIST(1:8*IZRAND) = IT1
         SDIST(1:  IZRAND) = DT1
         DEALLOCATE (IT1, DT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZRAND = I2RAND
         EXTMEM_RAND = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_RAND = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_RPAR(MNRPAR,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: RT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNRPAR .GT. MXRPAR)      &
     &                 .OR. (IZRPAR .GE. MXRPAR)) GOTO 100
         ALLOCATE(RT1(IZRPAR), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         RT1 = RANDPAR
         DEALLOCATE (RANDPAR, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2RPAR = MAX(MNRPAR, MIN(MXRPAR,2*IZRPAR))
         ALLOCATE( RANDPAR(I2RPAR), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         RANDPAR(1:IZRPAR) = RT1
         DEALLOCATE (RT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZRPAR = I2RPAR
         EXTMEM_RPAR = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_RPAR = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_RTYP(MNRTYP,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      TYPE (ROWINFO), ALLOCATABLE, DIMENSION (:) :: TT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNRTYP .GT. MXRTYP)      &
     &                 .OR. (IZRTYP .GE. MXRTYP)) GOTO 100
         ALLOCATE(TT1(IZRTYP), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         TT1 = T_ROWS
         DEALLOCATE (T_ROWS, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2RTYP = MAX(MNRTYP, MIN(MXRTYP,2*IZRTYP))
         ALLOCATE(T_ROWS(I2RTYP), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         T_ROWS(1:IZRTYP) = TT1
         DEALLOCATE (TT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZRTYP = I2RTYP
         EXTMEM_RTYP = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_RTYP = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_SCEN(MNSCEN,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER, ALLOCATABLE, DIMENSION (:) :: IT1, IT2
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNSCEN .GT. MXSCEN)      &
     &                 .OR. (IZSCEN .GE. MXSCEN)) GOTO 100
         ALLOCATE(IT1(IZSCEN), IT2(IZSCEN+1), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = LEAFND
         IT2 = KSCNAM
         DEALLOCATE (LEAFND, KSCNAM, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2SCEN = MAX(MNSCEN, MIN(MXSCEN,2*IZSCEN))
         ALLOCATE( KSCNAM(I2SCEN+1), LEAFND(I2SCEN), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         LEAFND(1:IZSCEN  ) = IT1
         KSCNAM(1:IZSCEN+1) = IT2
         DEALLOCATE (IT1, IT2, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZSCEN = I2SCEN
         EXTMEM_SCEN = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_SCEN = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_STOC(MNSTOC,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER, ALLOCATABLE, DIMENSION (:) :: IT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNSTOC .GT. MXSTOC)      &
     &                 .OR. (IZSTOC .GE. MXSTOC)) GOTO 100
         ALLOCATE(IT1(4*IZSTOC), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = KSTOCH
         DEALLOCATE (KSTOCH, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2STOC = MAX(MNSTOC, MIN(MXSTOC,2*IZSTOC))
         ALLOCATE( KSTOCH(4*I2STOC), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         KSTOCH(1:4*IZSTOC) = IT1
         DEALLOCATE (IT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZSTOC = I2STOC
         EXTMEM_STOC = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_STOC = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_STYP(MNSTYP,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER, ALLOCATABLE, DIMENSION (:) :: IT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNSTYP .GT. MXSTYP)      &
     &                 .OR. (IZSTYP .GE. MXSTYP)) GOTO 100
         ALLOCATE(IT1(IZSTYP), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = SOSTYP
         DEALLOCATE (SOSTYP, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2STYP = MAX(MNSTYP, MIN(MXSTYP,2*IZSTYP))
         ALLOCATE(SOSTYP(I2STYP), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         SOSTYP(1:IZSTYP) = IT1
         DEALLOCATE (IT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZSTYP = I2STYP
         EXTMEM_STYP = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_STYP = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_TCOL(MNTCOL,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      TYPE (TIMENAME),    ALLOCATABLE, DIMENSION (:) :: TT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNTCOL .GT. MXTCOL)      &
     &                 .OR. (IZTCOL .GE. MXTCOL)) GOTO 100
         I2TCOL = MAX(MNTCOL, MIN(MXTCOL,2*IZTCOL))
         IF (ALLOCATED(TCOLS)) THEN
            ALLOCATE(TT1(IZTCOL), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            TT1 = TCOLS
            DEALLOCATE (TCOLS, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            ALLOCATE(TCOLS(I2TCOL), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            TCOLS(1:IZTCOL) = TT1
            DEALLOCATE (TT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
!
         IZTCOL = I2TCOL
         EXTMEM_TCOL = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_TCOL = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_TPER(MNTPER,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER,            ALLOCATABLE, DIMENSION (:) :: IT1, IT2, IT3
      CHARACTER*(NAMLEN), ALLOCATABLE, DIMENSION (:) :: DT1, DT2, DT3
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNTPER .GT. MXTPER)      &
     &                 .OR. (IZTPER .GE. MXTPER)) GOTO 100
         I2TPER = MAX(MNTPER, MIN(MXTPER,2*IZTPER))
         IF (ALLOCATED(IRNGE0)) THEN
            ALLOCATE(IT1(IZTPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            IT1 = IRNGE0
            DEALLOCATE (IRNGE0, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            ALLOCATE( IRNGE0(I2TPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            IRNGE0(1:IZTPER) = IT1
            DEALLOCATE (IT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
!
         IF (ALLOCATED(INUSE)) THEN
            ALLOCATE(IT1(IZTPER), IT2(IZTPER), IT3(IZTPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            IT1 = INUSE
            IT2 = KPATH
            IT3 = KREF
            DEALLOCATE (INUSE, KPATH, KREF, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            ALLOCATE( INUSE(I2TPER), KPATH(I2TPER),      &
     &                               KREF (I2TPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            INUSE(1:IZTPER) = IT1
            KPATH(1:IZTPER) = IT2
            KREF (1:IZTPER) = IT3
            DEALLOCATE (IT1, IT2, IT3, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
!
         IF (ALLOCATED(DTIME)) THEN
            ALLOCATE(DT1(IZTPER), DT2(IZTPER), DT3(IZTPER), STAT=IER)
            IF (IER .NE. 0) GOTO 100
            DT1 = DTIME
            DT2 = DTIMEC
            DT3 = DTIMER
            DEALLOCATE(DTIME, DTIMEC, DTIMER, STAT=IER)
            IF (IER .NE. 0) GOTO 100
            I2TPER = MIN(MXTPER,2*IZTPER)
            ALLOCATE(DTIME(I2TPER), DTIMEC(I2TPER),      &
     &                              DTIMER(I2TPER), STAT=IER)
            IF (IER .NE. 0) GOTO 100
            DTIME (1:IZTPER) = DT1
            DTIMEC(1:IZTPER) = DT2
            DTIMER(1:IZTPER) = DT3
            DEALLOCATE(DT1, DT2, DT3, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
         IZTPER = I2TPER
         EXTMEM_TPER = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_TPER = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_TROW(MNTROW,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      TYPE (TIMENAME),    ALLOCATABLE, DIMENSION (:) :: TT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNTROW .GT. MXTROW)      &
     &                 .OR. (IZTROW .GE. MXTROW)) GOTO 100
         I2TROW = MAX(MNTROW, MIN(MXTROW,2*IZTROW))
         IF (ALLOCATED(TROWS)) THEN
            ALLOCATE(TT1(IZTROW), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            TT1 = TROWS
            DEALLOCATE (TROWS, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            ALLOCATE(TROWS(I2TROW), STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
            TROWS(1:IZTROW) = TT1
            DEALLOCATE (TT1, STAT=IERR)
            IF (IERR .NE. 0) GOTO 100
         ENDIF
!
         IZTROW = I2TROW
         EXTMEM_TROW = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_TROW = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_VNAM(MNVNAM,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER, ALLOCATABLE, DIMENSION (:) :: IT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNVNAM .GT. MXVNAM)      &
     &                 .OR. (IZVNAM .GE. MXVNAM)) GOTO 100
         ALLOCATE(IT1(IZVNAM), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = NAME1
         DEALLOCATE (NAME1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2VNAM = MAX(MNVNAM, MIN(MXVNAM,2*IZVNAM))
         ALLOCATE(NAME1(I2VNAM), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         NAME1(1:IZVNAM) = IT1
         DEALLOCATE (IT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZVNAM = I2VNAM
         EXTMEM_VNAM = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_VNAM = .FALSE.
         RETURN
         END FUNCTION
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      LOGICAL FUNCTION EXTMEM_VTYP(MNVTYP,IERR)
!
!     Auxiliary routine to allocate more memory in core file. The function
!     returns .TRUE. if everything went well and .FALSE. otherwise. The
!     error from ALLOCATE/DEALLOCATE is returned in IERR.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      INTEGER, ALLOCATABLE, DIMENSION (:) :: IT1
!
      IERR = 0
      IF (.NOT. RESIZE .OR. (MNVTYP .GT. MXVTYP)      &
     &                 .OR. (IZVTYP .GE. MXVTYP)) GOTO 100
         ALLOCATE(IT1(IZVTYP), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IT1 = VARTYP
         DEALLOCATE (VARTYP, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         I2VTYP = MAX(MNVTYP, MIN(MXVTYP,2*IZVTYP))
         ALLOCATE(VARTYP(I2VTYP), STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         VARTYP(1:IZVTYP) = IT1
         DEALLOCATE (IT1, STAT=IERR)
         IF (IERR .NE. 0) GOTO 100
         IZVTYP = I2VTYP
         EXTMEM_VTYP = .TRUE.
         RETURN
!
  100 CONTINUE
         EXTMEM_VTYP = .FALSE.
         RETURN
         END FUNCTION
!
      END MODULE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      MODULE SMPS_READER
!
      CONTAINS
!
      SUBROUTINE INPUT ( IERR )
!
!     This subroutine is the top level input routine. It first reads a
!     time file in the format laid out in Gassmann and Schweitzer and
!     described in broad terms below. It then calls further subroutines
!     to read the core file and the stoch file in one of the formats
!     described in the paper.
!
!     Routines called: INTIME, INCORE, INSTOC.
!
!   -----------------------------------------------------
!
!     A brief description of the input format follows:
!
!     All the information is contained in three input files, which are in
!     the order they are accessed:
!
!     - the  TIME FILE which breaks the rows and columns up into periods
!     - the  CORE FILE which contains information for a 'base scenario'
!     - the STOCH FILE which describes the stochastics of the problem
!
!     TIME FILE:
!
!     The first column and row of each period appear in the first two name
!     fields in standard MPSX format. The period is given a name in the
!     third name field. There is also an EXPLICIT form, signalled by the
!     presence of the keyword EXPLICIT on the PERIODS header. In this form
!     each row and column is explicitly assigned to a specific period in
!     separate ROW and COLUMN sections.
!
!
!     CORE FILE:
!
!     Standard MPSX format: The ROWS section lists all the rows for the
!     entire problem period by period, starting with period 1 and ending
!     with period T. The objective row is considered to be part of period 1.
!
!     The COLUMNS section is dealt with in the same way, columns are listed
!     period by period. The RHS, BOUNDS and RANGES sections follow as in the
!     MPSX standard. Segments for integer variables and special ordered sets
!     can be interspersed within the COLUMNS section.
!
!     NETGEN-type data cards for network components can also be accommodated.
!     An ARCS header switches to NETGEN format, and a COLUMNS header switches
!     back to MPS. SUPPLY and DEMAND are two additional headers which can be
!     interspersed in the RHS section.
!
!     Following the linear components, a QSECTION can be read to define
!     quadratic cost functions.
!
!
!     STOCH FILE:
!
!     Three mutually exclusive ways to specify random elements have been
!     implemented to date. All have as their goal to build the scenario
!     or event tree. There are two explicit specifications, scenario-wise
!     and node-wise, and one implicit specification (with many options).
!
!     The SCENario option allows dependence across time periods, but assumes
!     that all nodes in the decision tree belonging to the same time period
!     have identical problem dimensions and sparsity pattern in the
!     constraint matrix.
!
!     To facilitate varying problem dimensions, trap states, and the like,
!     one may use the TREE option, which requires a separate section for
!     each node in the decision tree. To avoid duplication, it is possible
!     to copy information from one node to another.
!
!     Implicit event tree specification assumes that the random elements
!     are independent from one period to the next (although dependence
!     within each period is permissible).
!
!     Independent random elements are specified with the keyword INDEP,
!     one element per data record. The system supports any distribution
!     desired, but only discrete, uniform, normal, gamma, beta and lognormal
!     distributions have been implemented. The user must provide the code
!     for any other distribution, if desired. The BLOCKS option is used
!     for blocks of random data which vary jointly but still exhibit
!     period-to-period independence. These two options can be mixed freely.
!     The only multivariate distributions implemented to date are discrete
!     and multivariate normal distributions; others are supported but code
!     must be supplied by the user.
!
!
!
!     Not all three files have to be present. The following table gives
!     all possible combinations:
!
!      -----------------------------------------------------------------
!     |  Time  |  Core  |  Stoch |              Remarks                 |
!     |  file  |  file  |  file  |                                      |
!      -----------------------------------------------------------------
!     |  yes   |  yes   |  yes   |  This is the normal case             |
!      -----------------------------------------------------------------
!     |  yes   |  yes   |   no   |  Problem is deterministic            |
!      -----------------------------------------------------------------
!     |  yes   |   no   |  yes   |  Only legal if NODES option is used  |
!      ----------------------------------------------------------------
!     |  yes   |   no   |   no   |  Not enough information to solve     |
!      -----------------------------------------------------------------
!     |   no   |  yes   |  yes   |  Only legal if NODES option is used  |
!      -----------------------------------------------------------------
!     |   no   |  yes   |   no   |  One-period deterministic problem    |
!      -----------------------------------------------------------------
!     |   no   |   no   |  yes   |  Only legal if NODES option is used  |
!      -----------------------------------------------------------------
!     |   no   |   no   |   no   |  Not enough information to solve     |
!      -----------------------------------------------------------------
!
!
!     In all cases, the program attempts to minimize storage by reducing
!     redundancies as much as possible.
!
!     -----------------------------------------------------------------
!
!     The internal representation is as follows.
!
!     DISCRETE distributions:
!
!     For each node N in the decision tree, N = 1,...,NODES,
!
!     find                         in array     with offset address
!     A matrix coefficients        A            KELMA (KDATA(N)+LMTX) (+)
!     A matrix locations           IA           KELMA (KDATA(N)+LMTX)
!     A matrix column pointers     LA           KCOLA (KDATA(N)+LMTX)
!     cost coefficients            COST         KCOST (N)
!     variable names               QNAMES       NAME1 (KNAMES(N)+K)   (++)
!     upper bounds                 XUB          KBOUND(N)
!     lower bound                  XLB          KBOUND(N)
!     right hand sides             RHS          KRHS  (N)
!     variable type                VARTYP       VARPTR(N)             (+++)
!     Q matrix coefficients        AQMTX        IQOFF (KDATQ(N)+LMTX) (++++)
!     Q matrix locations           IQMTX        IQOFF (KDATQ(N)+LMTX)
!     Q matrix column pointers     LQMTX        LQOFF (KDATQ(N)+LMTX)
!
!     decision variables           X            KROW  (N)
!     dual variables               YPI          KCOL  (N)
!
!     (+) LMTX = 1 for blocks on the main diagonal
!              = 2 for blocks immediately to the LEFT of the main diagonal
!              > 2 for blocks further away.
!         Staircase problems are indicated by MARKOV = .TRUE.
!
!     (++) Character strings may have differing lengths. The first character
!          of each string is marked in array NAME1, and the location of each
!          string is calculated from the offset for its node, recorded in
!          KNAMES. The length of the k'th string is NAME1(k+1) - NAME1(k).
!
!     (+++) This array identifies integer variables and special ordered sets.
!           VARTYP = 0 for continuous variables
!                  = 1 for integer variables
!                  < 0 for variables participating in special ordered set
!                      |VARTYP|. The type of the SOS (type 1 or 2) is stored
!                      in location SOSTYP(|VARTYP|).
!
!     (++++) LMTX = 1 for blocks on the main diagonal
!                 = 2 for blocks immediately to the LEFT of the main diagonal
!                 > 2 for blocks further away.
!            The number of blocks is indicated by the value of the variable
!            NQPROF. If NQPROF = 1, there is only one diagonal; if NQPROF = 2,
!            there is one diagonal and one subdiagonal block (linking a
!            period with the previous period); if NQPROF = 3, there are two
!            subdiagonal blocks, etc.
!
!
!     Problem dimensions are recorded in arrays NROW, NCOL (number of
!     columns including slacks) and NELMA.
!
!
!     Note that the identity matrix for the slack variables is *NOT* stored
!     as part of the A matrix and that the cost coefficients are separated,
!     even if costs are deterministic.
!
!
!     MSLiP solves the deterministic equivalent linear program, so all random
!     variables must eventually be discretized into an event tree. The input
!     format allows both explicit and implicit tree representation (including
!     continuous distributions). Different internal representations are used
!     for implicit and explicit trees.
!
!
!     In explicit form, the tree is represented by three pointer arrays IANCTR,
!     IDESC, IBROTH, which for each node give, respectively, the ancestor node,
!     the immediate successor node, and the next node in the same period.
!     If IBROTH > 0, then both nodes have the same ancestor, but it has
!     proven advantageous to link nodes in the same time period which have
!     different ancestors. This is indicated by a negative value for IBROTH,
!     and the next node in this case is given by ABS(IBROTH).
!
!
!     Each scenario has a name recorded in array QSCNAM. Since the length of
!     a scenario name may vary when free-form input is used, the offset of each
!     name is recorded in array KSCNAM. The association between the scenario
!     name and the nodes in the event tree is done through the pointer array
!     LEAFND, which for each scenario points to the leaf node on the path. The
!     other nodes can then be determined by iterating through the IANCTR array.
!     The number of scenarios is recorded in variable LASTLF.
!
!
!     ---------------------------------------------------------------------
!
!     The tree can also be given implicitly, in which case the deterministic
!     equivalent must be built (by calling routine MAKEDE) or a sample must
!     be generated (by calling routine SAMPLE). The implicit tree representation
!     is described below.
!
!     It is important to distinguish between STOCHASTIC ELEMENTS
!     which describe locations in the extended constraint matrix
!     (including RHS, cost, bounds, etc.) and RANDOM VARIABLES
!     which give information about their distributions.
!
!     The correspondence between random variables and stochastic elements
!     is usually one-to-one, but not if the LINTRAN option is used. The
!     linkage between random variables and stochastic elements is done in
!     the DISTRIBUTION MATRIX or D-matrix for short. It has block-diagonal
!     structure and is stored in sparse matrix form in arrays DMTX,
!     IDMTX and LDMTX. The random variables are the COLUMNS of this
!     matrix and the stochastic elements are the rows. Note that some
!     distributions define multi-dimensional random vectors. Such
!     random vectors occupy more than one column in the D-matrix.
!     Moreover, linear (affine) transformations usually involve more
!     than one random variable, including a degenerate one to describe
!     the translation component.
!
!
!     Data structures:
!     ----------------
!
!     Five data items must be stored for each random variable:
!
!     1. The distribution
!        In order to allow extensions and user defined routines, the
!        distributions are stored as strings in array SDIST.
!
!     2. The dimension
!     3. The number and location of integer parameters
!     4. The number and location of real parameters
!     5. The number and identity of parameters referenced by location
!        (including previous decision variables and historical data)
!     6. The number of the stage in which the random variable becomes known
!
!        These items are stored in array MDIST:
!        MDIST(8*IRV-7) gives the dimension of the random variable
!        MDIST(8*IRV-6) gives the offset for the first integer parameter
!                       (parameter values are in array IRNDPAR)
!        MDIST(8*IRV-5) gives the number of integer parameters
!        MDIST(8*IRV-4) gives the offset for the first real parameter
!                       (parameter values are in array RANDPAR)
!        MDIST(8*IRV-3) gives the number of real parameters
!        MDIST(8*IRV-2) gives the offset for the first location parameter
!                       (parameter values are in array XRNDPAR)
!        MDIST(8*IRV-3) gives the number of parameters referenced by location
!        MDIST(8*IRV  ) gives the number of the stage
!
!
!     The location and process mode of the stochastic elements are stored
!     in array KSTOCH in groups of four. The first entry gives the type
!     of random element, the second its block (or node), the third entry
!     describes the location within the relevant array, the fourth specifies
!     whether the values are to be treated as replacing the core file
!     information or if they are to be treated as additive or multiplicative
!     perturbations.
!
!     Types of random elements are as follows (can be extended as needed):
!
!        KSTOCH(4*NSTELM-3) = 1 for stochastic entry in an A-matrix block
!        KSTOCH(4*NSTELM-3) = 2 for stochastic cost
!        KSTOCH(4*NSTELM-3) = 3 for stochastic RHS
!        KSTOCH(4*NSTELM-3) = 4 for stochastic range
!        KSTOCH(4*NSTELM-3) = 5 for stochastic lower bound
!        KSTOCH(4*NSTELM-3) = 6 for stochastic upper bound
!        KSTOCH(4*NSTELM-3) = 7 for stochastic entry in a Q-matrix block
!        KSTOCH(4*NSTELM-3) = 11 for stochastic upper and lower bounds ('FX' type)
!        KSTOCH(4*NSTELM-3) = 12 for stochastic demand (network option)
!        KSTOCH(4*NSTELM-3) = 13 for stochastic supply (network option)
!
!     For stochastic A-matrix elements, the second entry is as in the
!     offset array KDATA, that is, depending on the value of MARKOV:
!
!        period           staircase problems      triangular problems
!          1                 1                       1
!          2                 3   2                   3   2
!          3                     5   4               6   5   4
!          4                         7  6           10   9   8   7
!         ...                          ...                ...
!
!
!     The location of each data item is recorded as the relative address
!     within the relevant array. Two utility routines, UNPACK and UNPACK2,
!     can be used to extract a column from one of the blocks making up the
!     constraint matrix.
!
!     For stochastic Q-matrix elements, IBLOCK is as in the offset array
!     KDATQ, that is, depending on the value of NQPROF:
!
!        period    NQPROF = 1       NQPROF = 2        NQPROF = 3      ...
!          1       1                1   4             1   4   9
!          2           2            3   2   7         3   2   7  14
!          3               3            6   5  10     8   6   5  12
!          4                   4            9   8        13  11  10
!         ...                ...             ...            ...
!
!
!     The processing mode is in location KSTOCH(4*NSTELM):
!
!        KSTOCH(4*NSTELM) = 1 if the core file info is to be replaced
!        KSTOCH(4*NSTELM) = 2 if the random information is additive
!        KSTOCH(4*NSTELM) = 3 if the information is multiplicative
!
!
!     Example:
!     ========
!
!     For a particular random element we have
!
!        KSTOCH(4*N + 1) = 1
!        KSTOCH(4*N + 2) = 5
!        KSTOCH(4*N + 3) = 8
!        KSTOCH(4*N + 4) = 2
!
!     The first entry identifies this as a A-random matrix element, so the core
!     file values are recorded in the array A (along with row indices in IA
!     and column indices in LA). The second entry specifies the block within
!     the matrix: the subdiagonal block linking periods 2 and 3 (independent
!     of whether the problem is Markovian or not). Offsets to this block are
!     then found in KELMA(5), and the particular element is in KELMA(5) + 8.
!     The fourth entry indicates that the stochastic value must be ADDED to the
!     core file value.
!
!     The linking between stochastic elements and random variables is done
!     in the `D'-matrix, which has block-diagonal structure and is stored
!     in sparse form in arrays DMTX, IDMTX and LDMTX. The stochastic
!     elements are the ROWS of this matrix, and the random variables are the
!     columns.
!
!
!     After the distributions of the random variables have been specified,
!     probabilistic constraints can be defined in a separate section. Chance
!     constraints are stored in three arrays ICHANCE, LCHANCE and PCHANCE.
!     The number of chance constraints is recorded in the variable NCHANCE.
!     PCHANCE gives the probability and direction of the probabilistic
!     inequality (>= 0 if the inequality represents a reliability level,
!     <= 0 for fault tolerances --- in the latter case, the reliability level
!     is ABS(PCHANCE(i))). LCHANCE acts as a pointer into the array ICHANCE.
!     For each probabilistic constraint i the number of rows participating
!     in the constraints are found (in absolute address form) in locations
!     ICHANCE(LCHANCE(i))... ICHANCE(LCHANCE(i+1)-1).
!
!     ----------------------------------------------------------------------
!
!     Development history:
!
!    21 Aug 2002: Translation to Fortran 90.
!
!    29 Sep 1997: Major rewrite to permit far greater sophistication
!                 in dealing with discrete distributions, as documented
!                 in the paper by Gassmann and Schweitzer.
!
!    29 April 89: First coding of the continuous distributions.
!
!    14 April 89: Restructured the construction of the decision tree for
!                 INDEP and BLOCKS options to detect coefficients which
!                 may become known in period t but do in fact belong to a
!                 later period. This is allowed, since information known at
!                 the outset is not confined to period 1, either.
!
!     1 Feb 1989: First attempts to detect the absence of input files in
!                 certain situations: For deterministic problems there is
!                 no need to read a stoch file. One-period LP problems
!                 could be specified by a core file alone, and the TREE
!                 structure does not require anything but a stoch file.
!
!  ---------------------------------------------------------------------
!
!                    ***DESCRIPTION OF PARAMETERS***
!
!        IOBJ1  = ORIGINAL OBJECTIVE ROW (MAY BE ZERO)
!                    THE OBJECTIVE ROW IS INTERCEPTED AND SWAPPED TO
!                    POSITION 1 FOR EASIER IDENTIFICATION IN SUBPROBLEMS
!                    IN PERIODS 2, 3, ..., T.
!
!        IERR   = ERROR STATUS
!                    = 0 - Normal termination
!                    > 0 - Error condition
!
!        STFILE = Form of the stoch file
!                    = 0 - none found (deterministic problem)
!                    = 1 - explicit event tree (TREE option)
!                    = 2 - explicit event tree (SCENARIO option)
!                    = 3 - implicit event tree (only DISCRETE distributions)
!                    = 4 - implicit event tree (continuous distributions)
!  ---------------------------------------------------------------------
!
!                  This version dated 21 August, 2002.
!
!  ---------------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE TIMER_STAT
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
      INTEGER  EXPL
!
      TIMESET = .FALSE.
      CALL STIMER(TIME0)
!
!     *********************************************
!     *****       PROCESS THE TIME FILE       *****
!     *********************************************
!
      CALL INTIME ( EXPL, IRTEMP, ICTEMP, IERR )
      IF (TIMSTAT .GE. 2) GOTO 9999
!
!     *********************************************
!     *****       PROCESS THE CORE FILE       *****
!     *********************************************
!
      CALL INCORE ( EXPL+1, IRTEMP, ICTEMP, IERR )
      IF (CORSTAT .GE. 2) GOTO 9999
!
!     *********************************************
!     *****       PROCESS THE STOCH-FILE      *****
!     *********************************************
!
      CALL INSTOC ( IERR )
      IF (STOSTAT .GE. 2) GOTO 9999
!
!     Normal return
!
      CALL STIMER(TIME1)
      WRITE (IOLOG, 1000) TIME1
      RETURN
!
!     Error return
!
 9999 CONTINUE
         IERR = 2
         RETURN
 1000 FORMAT(' Time to read input:',F10.3,' sec.')
!
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INTIME ( EXPL, IRTEMP, ICTEMP, IERR )
!
!     This routine reads the PERIODS section of the time file and if the
!     core file is in temporal order (signalled by the keyword EXPLICIT
!     on the PERIODS record) stores the first row and column of each period
!     in arrays for later processing in routine INCORE.
!
!
!     Routines called: GETNAM, GETMPS, CHOPEN
!
!     -----------------------------------------------------------------
!               This version dated 22 October 2002.
!     -----------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN)  DNAME(3), DBLANK, DTEMP, QLINE
      CHARACTER*8 DPER1
      CHARACTER*1 QBLANK(NAMLEN), QTEMP(NAMLEN)
      REAL*8      ATEMP(4)
      INTEGER*4   LENGTH(3),EXPL
!
      EQUIVALENCE (QBLANK,DBLANK), (QTEMP,DTEMP)
      DATA DPER1 /'PERIOD1 '/,  QBLANK/NAMLEN*' '/
!
! ----------------------------------------------------------------------
!
!     Allocate the name arrays
!
      ALLOCATE(DTIME(IZTPER), DTIMEC(IZTPER), DTIMER(IZTPER), STAT=IERR)
      IF (IERR .NE. 0) GOTO 9060
!
! ----------------------------------------------------------------------
!
!     Initializations
!
      NREC    = 0
      IENDAT  = 0
      TIMSTAT = 0
      PROBNM  = DBLANK
      NPER    = 0
      EXPL    = 0
      L       = 0
      IF (NECHO1 .GE. 1) WRITE (IOLOG, 1100)
!
!     Start by opening the time file. The subroutine FILE_WRAPPER is a
!     user-defined routine whose arguments have the following meaning:
!     TIMFIL - name of the time file,        INTENT(IN),    declared as CHARACTER*255.
!     IOTIM  - Fortran I/O channel,          INTENT(IN),    declared as INTEGER.
!     3        File status; 3 = 'UNKNOWN',   INTENT(IN),    declared as INTEGER.
!     2        File format; 2 = 'FORMATTED', INTENT(IN),    declared as INTEGER.
!     IOLOG  - File for error messages,      INTENT(IN),    declared as INTEGER.
!     IERR   - Error status                  INTENT(INOUT), declared as INTEGER.
!
      CALL FILE_WRAPPER(TIMFIL,IOTIM,3,2,IOLOG,IERR)
      IF (IERR .GT. 0) GOTO 110
!
!     Get the name of the problem first
!
  100 CONTINUE
      CALL GETNAM(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOTIM,NREC,LENGTH)
      IF (Q1     .EQ.  QAST) GOTO 100
      IF (IERR   .GT.    0 ) GOTO 110
      IF (Q1 .EQ. QT .AND. Q2 .EQ. QI) GOTO 120
         WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &          DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 1400)
         GOTO 9999
!
!     Error while reading the TIME file. Treat as missing and proceed.
!
  110 CONTINUE
         IF (NREC .GT. 0) WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                    DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3970)
         DTIMEC(1) = DOTS
         DTIMER(1) = DOTS
         DTIME(1)  = DPER1
         TIMSTAT = 1
         PROBNM  = DOTS
         NPER  = 1
         GOTO 900
!
!     We have found a TIME file and a name for our problem.
!
  120 CONTINUE
         IENDAT = 1
         PROBNM = DNAME(1)
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 1200)      &
     &      NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2)
!
!     Now read the rest of the TIME file
!
Read_Line: &
     & DO
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                        IERR,IOTIM,NREC,LENGTH)
         IF (Q1 .EQ. QAST ) CYCLE
!
!     Error while reading the TIME file.
!
         IF (IERR   .GT. 0) THEN
            IF (IENDAT .EQ. 1) GOTO 9050
               WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
               WRITE (IOLOG, 1700)
               GOTO 900
!
!     We have found the PERIODS card. The second name field might contain
!     information (IMPLICIT/EXPLICIT) which needs to be extracted.
!
         ELSEIF (Q1 .EQ. QP .AND. Q2 .EQ. QE) THEN
            IENDAT = 2
            IF (NECHO1 .GE. 2) WRITE (IOLOG, 1200)      &
     &                         NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2)
            IF (DNAME(2) .EQ. 'EXPLICIT') EXPL = 1
            L = 1 + EXPL
!
!     ROWS section header
!
         ELSEIF (Q1 .EQ. QR .AND. Q2 .EQ. QO) THEN
            IF (EXPL .EQ. 0) GOTO 9040
            IF (NECHO1 .GE. 1) WRITE (IOLOG, 1200)      &
     &         NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2)
            IZTROW = IZROWP*NPER
            IZTCOL = IZCOLP*NPER
            MXTROW = MXROWP*NPER
            MXTCOL = MXCOLP*NPER
            ALLOCATE(TROWS(IZTROW), TCOLS(IZTCOL), STAT=IERR)
            IF (IERR .NE. 0) GOTO 9060
            IRTEMP = 0
            ICTEMP = 0
            L = 3
!
!     COLS section header
!
         ELSEIF (Q1 .EQ. QC .AND. Q2 .EQ. QO) THEN
            IF (L .NE. 3) GOTO 9070
            IF (NECHO1 .GE. 1) WRITE (IOLOG, 1200)      &
     &         NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2)
            L = 4
!
!     ENDATA record
!
         ELSEIF (Q1 .EQ. QE .AND. Q2 .EQ. QN) THEN
            GOTO 270
!
!     Data card
!
         ELSEIF (Q1 .EQ. QBL) THEN
            SELECT CASE (L)
!
!     Regular data card. Extend storage arrays if necessary.
!
               CASE (1,2)
                  NPER = NPER + 1
                  IF (NPER .GT. IZTPER) THEN
                     IF (.NOT. EXTMEM_TPER(NPER,IERR)) GOTO 9030
                  ENDIF
!
!     Store the information
!
                  IF (DNAME(3) .NE. DBLANK) THEN
                     IF (NECHO1 .GE. 2) WRITE (IOLOG, 2500)      &
     &                                  NREC,DNAME(3),DNAME(2),DNAME(1)
                     DTIMEC(NPER) = DNAME(1)
                     DTIMER(NPER) = DNAME(2)
                     DTIME (NPER) = DNAME(3)
                  ELSE
                     IF (EXPL .EQ. 1) THEN
                        IF (NECHO1 .GE. 2)      &
     &                      WRITE (IOLOG, 2501) NREC,DNAME(1)
                        DTIME(NPER) = DNAME(1)
                     ELSE
                        GOTO 9080
                     ENDIF
                  ENDIF
!
!     Data record in ROWS section
!
               CASE (3)
                  IF (NECHO1 .GE. 3) WRITE (IOLOG, 1200)      &
     &                NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2)
                  DO IP=1,NPER
                     IF (DNAME(2) .EQ. DTIME(IP)) GOTO 210
                  END DO
                  GOTO 9120
!
 210             CONTINUE
                 IRTEMP = IRTEMP + 1
                 IF (IRTEMP .GT. IZTROW) THEN
                    IF (.NOT. EXTMEM_TROW(IRTEMP,IERR)) GOTO 9130
                 ENDIF
                 TROWS(IRTEMP)%PERIOD = IP
                 TROWS(IRTEMP)%NAME   = DNAME(1)
!
!     Data record in COLUMNS section
!
               CASE (4)
                  IF (NECHO1 .GE. 3) WRITE (IOLOG, 1200)      &
     &                NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2)
                  DO IP=1,NPER
                     IF (DNAME(2) .EQ. DTIME(IP)) GOTO 250
                  END DO
                  GOTO 9120
!
  250             CONTINUE
                  ICTEMP = ICTEMP + 1
                  IF (ICTEMP .GT. IZTCOL) THEN
                     IF (.NOT. EXTMEM_TROW(ICTEMP,IERR)) GOTO 9130
                  ENDIF
                  TCOLS(ICTEMP)%PERIOD = IP
                  TCOLS(ICTEMP)%NAME   = DNAME(1)
            END SELECT
!
!     Illegal record
!
         ELSE
            WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &             DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
            WRITE (IOLOG, 1400)
            GOTO 9999
         ENDIF
      END DO Read_Line
!
!     End of TIME file
!
  270 CONTINUE
         IF (L .EQ. 0) GOTO 9090
         IF (L .EQ. 2) GOTO 9100
         IF (L .EQ. 3) GOTO 9110
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 1200)      &
     &      NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2)
         GOTO 900
!
!     Come here if anything went wrong
!
 9030 CONTINUE
         WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3030) IZTPER
         GOTO 9999
!
 9040 CONTINUE
         WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3040) IZTPER
         GOTO 900
!
 9050 CONTINUE
         WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3050)
         GOTO 9999
!
 9060 CONTINUE
         WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3060)
         GOTO 9999
!
 9070 CONTINUE
         WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3070)
         GOTO 9999
!
 9080 CONTINUE
         WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3080)
         GOTO 9999
!
 9090 CONTINUE
         WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3090)
         GOTO 900
!
 9100 CONTINUE
         WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3100)
         GOTO 9999
!
 9110 CONTINUE
         WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3110)
         GOTO 9999
!
 9120 CONTINUE
         WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3120)
         GOTO 9999
!
 9130 CONTINUE
         WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3130)
         GOTO 9999
!
 9999 CONTINUE
         TIMSTAT = 2
  900 CONTINUE
         CLOSE (IOTIM)
         IERR = TIMSTAT
         RETURN
!
 1100 FORMAT(/,' Process TIME file:')
 1200 FORMAT(I8,4X,4A1,2X,A8,2X,A8,2X,F12.4,3X,A8,2X,F12.4)
 1400 FORMAT(' XXX -  FATAL  - Illegal record in TIME file')
 1700 FORMAT(' XXX - WARNING - Missing ENDATA card in TIME file')
 2500 FORMAT(I8,4X,' Period ',A32,/,      &
     &       ' -- first row ',A32,', first column ',A32)
 2501 FORMAT(I8,4X,' Period ',A32)
 3030 FORMAT(' XXX -  FATAL  - Too many periods specified: Use at most',      &
     &       I4)
 3040 FORMAT(' XXX -  INFO   - ROW/COLUMN section redundant - ignore',      &
     &       I4)
 3050 FORMAT(' XXX -  FATAL  - Detected EOF while reading TIME file')
 3060 FORMAT(' XXX -  FATAL  - Unable to allocate memory in INTIME')
 3070 FORMAT(' XXX -  FATAL  - Missing ROWS section in explicit TIME',      &
     &       ' file')
 3080 FORMAT(' XXX -  FATAL  - Missing name fields in implicit TIME',      &
     &       ' file')
 3090 FORMAT(' XXX - WARNING - Empty TIME file detected')
 3100 FORMAT(' XXX -  FATAL  - ROWS section missing in explicit TIME',      &
     &       ' file')
 3110 FORMAT(' XXX -  FATAL  - COLUMNS section missing in explicit',      &
     &       ' TIME file')
 3120 FORMAT(' XXX -  FATAL  - Name of stage not recognized')
 3130 FORMAT(' XXX -  FATAL  - Problem dimensions exceed storage',      &
     &       ' capacity')
 3970 FORMAT(' XXX - WARNING - Error during READ or non-existent TIME',      &
     &       ' file')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INCORE ( IFORM, IRTEMP, ICTEMP, IERR )
!
!     This subroutine reads the core file in the modified MPS format
!     described in the standards paper.
!
!     The task is broken into several pieces: the preamble and cleanup
!     are in this routine, which calls separate routines for ROWS,
!     COLUMNS and remaining sections, as well as a QUAD section if there
!     is one. ROWS and COLUMNS sections can be in temporal order or not;
!     the format is detected in INTIME and communicated via the parameter
!     IFORM: 1 - ROWS and COLUMNS are sorted in temporal order;
!            2 - explicit stage information is provided in the time file.
!
!     In addition to straight LP problems, problems can be specified in
!     the NETGEN format that is also part of the standards paper. It is
!     further possible to mix the two formats. This version of the reader
!     supports both fixed-format input (with eight-character name fields)
!     and free-format input (with name fields of arbitrary length).
!
!     Routines called: CHOPEN, GETNAM, GETMPS, GETNET, INQUAD.
!
!     --------------------------------------------------------------
!     This version dated 21 September 2002. Written by Gus Gassmann.
!     --------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      REAL*8 ATEMP(4)
      CHARACTER*(NAMLEN) DNAME(3), DBLANK, DROW, DCOL, QLINE
      CHARACTER*1 DB255(NAMLEN), QROW(NAMLEN), QCOL(NAMLEN)
      INTEGER*4   LENGTH(3)
!
      EQUIVALENCE (DBLANK,DB255), (DROW,QROW), (DCOL,QCOL)
!
      DATA        DB255 /NAMLEN*' '/
! ----------------------------------------------------------------------
!
!     Allocate the arrays
!
      ALLOCATE(A(IZALMN), IA(IZALMN), LA(IZACOL), KCOLA(IZABLK),      &
     &     KELMA(IZABLK), NELMA(IZABLK), KDATA(IZNODE), STAT=IERR)
      IF (IERR .NE. 0) GOTO 9060
      ALLOCATE (NODEINFO(IZNODE), STAT=IERR)
      IF (IERR .NE. 0) GOTO 9060
      ALLOCATE (XLB(IZBNDS), XUB(IZBNDS), COST(IZCOST),      &
     &          STAT=IERR)
      IF (IERR .NE. 0) GOTO 9060
      ALLOCATE (VARTYP(IZVTYP), SOSTYP(IZSTYP), NAME1(IZVNAM),      &
     &          QNAMES(IZCHAR), STAT=IERR)
      IF (IERR .NE. 0) GOTO 9060
      ALLOCATE (IRNGE0(IZTPER), STAT=IERR)
      IF (IERR .NE. 0) GOTO 9060
!
!     Initializations
!
         NODEINFO(1)%KBOUND = 0
         NODEINFO(1)%KROW   = 0
         NODEINFO(1)%KCOL   = 0
         NODEINFO(1)%KRHS   = 0
         NODEINFO(1)%KCOST  = 0
         NODEINFO(1)%KNAMES = 0
         NODEINFO(1)%IANCTR = 0
         NODEINFO(1)%IBROTH = 0
         NODEINFO(1)%VARPTR = 0
         NODEINFO(1)%PROB   = 1.D0
!
         IOBJ   = 1
         IOBJ1  = 0
         NPSEEN = 0
         INODE  = 1
         NROWS  = 0
         MAXROW = 0
!
         LASTA  = 0
         LASTBD = 0
         LASTC  = 0
         LASTCA = 0
         LASTD  = 0
         LASTNM = 0
         LASTR  = 0
         LASTCH = 0
         LASTVP = 0
!
         MARKOV = .TRUE.
         PURELP = .TRUE.
!
         KDATA(1) = 0
         KCOLA(1) = 0
         KELMA(1) = 0
         LA(1)    = 1
         NAME1(1) = 1
!
         NREC    = 0
         CORSTAT = 0
!
         IF (TIMSTAT .EQ.  2) GOTO 9960
         IF (TIMSTAT .EQ. -1) THEN
            WRITE (IOLOG, 2100)
            NPER = 1
         ENDIF
!
!     Open the core file
!
   60    CONTINUE
         CALL FILE_WRAPPER(CORFIL,IOCOR,3,2,IOLOG,IERR)
         IF (IERR .GT. 0) GOTO 80
!
         IF (NECHO1  .GE.  1) WRITE (IOLOG, 1000)
!
!     Read until we find the problem name
!
   70 CONTINUE
         CALL GETNAM(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOCOR,NREC,LENGTH)
         IF (Q1   .EQ. QAST) GOTO 70
         IF (IERR .GT.   0 ) GOTO 80
         IF (Q1   .EQ.  QN .AND. Q2 .EQ. QA) GOTO 90
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 2400)
            GOTO 9999
!
!     Error during read or missing CORE file. Note and keep going.
!
   80    CONTINUE
            IF (NREC .GT. 0) WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 3980)
            CORSTAT = 1
            IF (PROBNM .EQ. DBLANK) PROBNM = DOTS
            GOTO 999
!
!     Verify the problem name against the name in the time file.
!
   90 CONTINUE
         IF (PROBNM .EQ. DOTS) PROBNM = DNAME(1)
         IF (DNAME(1) .NE. PROBNM) GOTO 9160
            IF (NECHO1 .GE. 1) WRITE (IOLOG, 1400) NREC,QLINE
            CORSTAT = 0
            NPSEEN  = 1
!
!     Now read the sections of the core file in several subroutines.
!
  100 CONTINUE
         IF (IFORM .EQ. 1) THEN
            IZRTYP = IZTPER*IZROWP
            MXRTYP = 2147483647
            ALLOCATE(T_ROWS(IZRTYP), STAT=IERR)
            IF (IERR .NE. 0) GOTO 9060
!
!     Read ROWS section in temporal order
!
            CALL INROWS1(NREC, QLINE, NEXT, IENDAT, IERR)
            IF (CORSTAT .EQ. 2) THEN
               DEALLOCATE(T_ROWS, STAT=IERR)
               IF (IERR .GT. 0) GOTO 9060
                  GOTO 999
            ENDIF
            IF (IERR .GT. 0) THEN
            ENDIF
!
!     Read COLUMNS section in temporal order
!
            CALL INCOLS1(NREC, QLINE, NEXT, IENDAT, IERR)
            IF (IERR .GT. 0) THEN
            ENDIF
            DEALLOCATE(T_ROWS, STAT=IERR)
            IF (IERR    .GT. 0) GOTO 9060
            IF (CORSTAT .EQ. 2) GOTO 999
         ELSE
            ALLOCATE(T_VAR(IZCOLP,NPER), STAT=IERR)
            IF (IERR .NE. 0) GOTO 9060
!
!     Read ROWS section in arbitrary order
!
            CALL INROWS2 (NREC, QLINE, IRTEMP, NEXT, IENDAT, IERR)
            IF (CORSTAT .EQ. 2) THEN
               DEALLOCATE(T_VAR, STAT=IERR)
               IF (IERR .GT. 0) GOTO 9060
                  GOTO 999
            ENDIF
            IF (IERR .GT. 0) THEN
            ENDIF
!
!     Read COLUMNS section in arbitrary order
!
            CALL INCOLS2 (NREC, QLINE, ICTEMP, NEXT, IENDAT, IERR)
            IF (IERR .GT. 0) THEN
            ENDIF
            DEALLOCATE(T_VAR, STAT=IERR)
            IF (IERR    .GT. 0) GOTO 9060
            IF (CORSTAT .EQ. 2) GOTO 999
         ENDIF
!
!     Read RHS, BOUNDS and RANGES sections
!
         IF (NEXT .LE. 5) THEN
            CALL INRHS(NREC, QLINE, NEXT, IENDAT, IERR)
            IF (CORSTAT .EQ. 2) GOTO 999
         ENDIF
!
!     Read a QSECTION if so directed by the user
!
         IF ((NEXT .EQ. 6) .OR. QMAT) THEN
            CALL INQUAD ( NEXT, NREC, IERR )
            CORSTAT = MAX(IERR, CORSTAT)
            IF (CORSTAT .EQ. 2) GOTO 999
         ENDIF
!
!     End of core file. Compute the density of each main diagonal block.
!
         IF (NECHO1 .GE. 1) THEN
            WRITE (IOLOG, 1400) NREC,QLINE
            DO IP=1,NPER
               NROWS = NODEINFO(IP)%NROW
               NSCOL = NODEINFO(IP)%NCOL - NROWS
               RELEM = NELMA(KDATA(IP)+1)
               IF (NROWS * NSCOL .GT. 0) THEN
                  RDENS = RELEM / (NROWS * NSCOL)
               ELSE
                  RDENS = 0.D0
               ENDIF
               WRITE (IOLOG, 2200) IP, NROWS, NSCOL, RDENS
            END DO
         ENDIF
         GOTO 999

 9060 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3060)
         GOTO 9999
!
 9160 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3160)
         GOTO 9999
!
 9960 CONTINUE
         WRITE (IOLOG, 3960)
         GOTO 9999
!
 9999 CONTINUE
         CORSTAT = 2
!
!     Return to calling program
!
  999 CONTINUE
         CLOSE (IOCOR)
         IERR  = CORSTAT
         IF (CORSTAT .LT. 2) THEN
            NODES = NPER
         ENDIF
         RETURN
!
 1000 FORMAT(/,' Process CORE file:')
 1400 FORMAT(I8,2X,A80)
 2100 FORMAT(' XXX - WARNING - The time file must be read first.',      &
     &     /,' If you continue, the problem is assumed single-staged.')
 2200 FORMAT(' Period',I3,' has',I4,' rows and',I4,' columns.',      &
     &       ' Density of constraint matrix:',F6.3)
 2300 FORMAT(' XXX -  FATAL  - Illegal header card in CORE file.')
 2400 FORMAT(' XXX -  FATAL  - Name does not match info in time file.')
 3060 FORMAT(' XXX -  FATAL  - Unable to allocate memory in INCORE.')
 3160 FORMAT(' XXX -  FATAL  - Name does not match info in TIME file.')
 3960 FORMAT(' XXX -  FATAL  - Information from time file is in error.')
 3980 FORMAT(' XXX - WARNING - Error during read or CORE file missing.')
!
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INROWS1 ( NREC, QLINE, NEXT, IENDAT, IERR )
!
!     This routine reads the ROWS section in temporal form.
!
!     Routines called: CHOPEN, GETNAM, GETMPS, GETNET.
!
!     --------------------------------------------------------------
!     This version dated 21 September 2002. Written by Gus Gassmann.
!     --------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      REAL*8 ATEMP(4)
      CHARACTER*(NAMLEN) DNAME(3), QLINE, DBLANK, DNXROW
      CHARACTER*8 DUNCON
      CHARACTER*1 DB255(NAMLEN)
      INTEGER*4   LENGTH(3)
!
      EQUIVALENCE (DBLANK,DB255)
!
      DATA DUNCON/'''UNCONS'''/, DB255 /NAMLEN*' '/
! ----------------------------------------------------------------------
!
!     Store the name of the next period so it can be identified later.
!
         IENDAT = 0
         IROW   = 0
         DNXROW = DBLANK
         DO NXTPER=2,NPER
            IF (DTIMER(NXTPER) .NE. DUNCON) THEN
               DNXROW = DTIMER(NXTPER)
               GOTO 100
            ENDIF
         END DO
         NXTPER = NPER + 1
!
!     Read the next record, then distribute for further processing.
!     =============================================================
!
  100 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOCOR,NREC,LENGTH)
         IF (Q1 .EQ. QAST) GOTO 100
         IF (IERR .GT. 0 ) GOTO 9000
         IF (Q1 .EQ. QBL ) GOTO 120
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QO) GOTO 105
         IF (Q1 .EQ. QN .AND. Q2 .EQ. QO) GOTO 105
         IF (Q1 .EQ. QA .AND. Q2 .EQ. QR) GOTO 200
         IF (Q1 .EQ. QC .AND. Q2 .EQ. QO) GOTO 200
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QH) GOTO 9300
         IF (Q1 .EQ. QS .AND. Q2 .EQ. QU) GOTO 9300
         IF (Q1 .EQ. QD .AND. Q2 .EQ. QE) GOTO 9300
         IF (Q1 .EQ. QB .AND. Q2 .EQ. QO) GOTO 9300
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QA) GOTO 9300
         IF (Q1 .EQ. QQ .AND. Q2 .EQ. QS) GOTO 9300
         IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 9300
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 2300)
            GOTO 9999
!
!     ROW or NODE card. Make room for the objective row, even if temporary.
!
  105    CONTINUE
            IF (IENDAT .EQ. 0) THEN
               NROWS  = 1
               MAXROW = 1
               NODEINFO(1)%NROW = 1
               T_ROWS(1)%NAME   = OBJNAM
               T_ROWS(1)%LENGTH = LENOBJ
               T_ROWS(1)%ROWTYP = 2
               T_ROWS(1)%LOWER  = -PLINF
               T_ROWS(1)%UPPER  =  PLINF
               IENDAT = 1
            ENDIF
            IF (NECHO1 .GE. 1) WRITE (IOLOG, 1400) NREC,QLINE
            GOTO 100
!
!     Here we have a regular data card. Check whether it starts a new period.
!
  120    CONTINUE
            IF (NECHO1 .GE. 3) WRITE (IOLOG, 1400)      &
     &         NREC,QLINE
            IF (DNAME(1) .NE. DNXROW) GOTO 135
!
!      First row of a new time period. Enlarge arrays if needed.
!
  125       CONTINUE
            IF (NXTPER .GT. IZTPER) THEN
               IF (.NOT. EXTMEM_TPER(NXTPER,IERR)) GOTO 9060
            ENDIF
!
            IF (NXTPER .GT. IZNODE) THEN
               IF (.NOT. EXTMEM_NODE(NXTPER,IERR)) GOTO 9060
            ENDIF
!
         DO IP=NPSEEN+1,NXTPER
            CALL NEWNODE(IERR)
            IF (IERR .GT. 0) GOTO 9999
         END DO
!
         IROW   = NODEINFO(INODE)%KROW
         NROWS  = 1
         NPSEEN = NXTPER
         DNXROW = DBLANK
         DO NXTPER=NPSEEN+1,NPER
            IF (DTIMER(NXTPER) .NE. DUNCON) THEN
               DNXROW = DTIMER(NXTPER)
               GOTO 130
            ENDIF
         END DO
         NXTPER = NPER + 1
!
  130    CONTINUE
         IF (NPSEEN .LT. IZTPER) THEN
            DNXROW = DTIMER(NPSEEN+1)
         ELSE
            DNXROW = DBLANK
         ENDIF
!
!     Test row type. Expand memory if necessary.
!
  135    CONTINUE
            IF (NROWS  .GT. IZROWP) THEN
               IF (.NOT. RESIZE .OR. (IZROWP .GE. MXROWP)      &
     &                          .OR. (NROWS  .GT. MXROWP)) GOTO 9060
               IZROWP = NROWS
            ENDIF
            IF (MAXROW .GT. IZROWS) THEN
               IF (.NOT. RESIZE .OR. (IZROWS .GE. MXROWS)      &
     &                          .OR. (MAXROW .GT. MXROWS)) GOTO 9060
               IZROWS = MAXROW
            ENDIF
!
            IF (MAXROW .GE. IZVNAM) THEN
               IF (.NOT. EXTMEM_VNAM(MAXROW+1,IERR)) GOTO 9060
            ENDIF
!
            MNRTYP = IROW + NROWS + 1
            IF (MNRTYP .GT. IZRTYP) THEN
               IF (.NOT. EXTMEM_RTYP(MNRTYP,IERR)) GOTO 9060
            ENDIF
!
            IF ((Q2 .EQ. QE) .OR. (Q3 .EQ. QE) .OR.      &
     &         ((Q2 .EQ. QBL).AND.(Q3 .EQ. QBL))) THEN
               MAXROW = MAXROW + 1
               NROWS  = NROWS  + 1
               NODEINFO(INODE)%NROW = NROWS
               T_ROWS(IROW+NROWS)%NAME   = DNAME(1)
               T_ROWS(IROW+NROWS)%LENGTH = LENGTH(1)
               T_ROWS(IROW+NROWS)%ROWTYP = 0
               T_ROWS(IROW+NROWS)%LOWER  = 0.D0
               T_ROWS(IROW+NROWS)%UPPER  = 0.D0
!
            ELSEIF  ((Q2 .EQ. QG) .OR. (Q3 .EQ. QG)) THEN
               MAXROW = MAXROW + 1
               NROWS  = NROWS  + 1
               NODEINFO(INODE)%NROW = NROWS
               T_ROWS(IROW+NROWS)%NAME   = DNAME(1)
               T_ROWS(IROW+NROWS)%LENGTH = LENGTH(1)
               T_ROWS(IROW+NROWS)%ROWTYP = -1
               T_ROWS(IROW+NROWS)%LOWER  = -PLINF
               T_ROWS(IROW+NROWS)%UPPER  = 0.D0
!
            ELSEIF  ((Q2 .EQ. QL) .OR. (Q3 .EQ. QL)) THEN
               MAXROW = MAXROW + 1
               NROWS  = NROWS  + 1
               NODEINFO(INODE)%NROW = NROWS
               T_ROWS(IROW+NROWS)%NAME   = DNAME(1)
               T_ROWS(IROW+NROWS)%LENGTH = LENGTH(1)
               T_ROWS(IROW+NROWS)%ROWTYP = 1
               T_ROWS(IROW+NROWS)%LOWER  = 0.D0
               T_ROWS(IROW+NROWS)%UPPER  = PLINF
!
            ELSEIF  ((Q2 .EQ. QN) .OR. (Q3 .EQ. QN)) THEN
!
!     Unbounded row could be the objective. If not, treat like any other row.
!
               IF (OBJNAM .EQ. DOTS) THEN
                  OBJNAM = DNAME(1)
                  LENOBJ = MAX(LENGTH(1),8)
                  DO I=1,NPSEEN
                     IROBJ = NODEINFO(I)%KROW + 1
                     T_ROWS(IROBJ)%NAME   = OBJNAM
                     T_ROWS(IROBJ)%LENGTH = LENOBJ
                  END DO
               ELSEIF (OBJNAM .NE. DNAME(1)) THEN
                  MAXROW = MAXROW + 1
                  NROWS  = NROWS  + 1
                  NODEINFO(INODE)%NROW = NROWS
                  T_ROWS(IROW+NROWS)%NAME   = DNAME(1)
                  T_ROWS(IROW+NROWS)%LENGTH = LENGTH(1)
                  T_ROWS(IROW+NROWS)%ROWTYP = 2
                  T_ROWS(IROW+NROWS)%LOWER  = -PLINF
                  T_ROWS(IROW+NROWS)%UPPER  =  PLINF
               ENDIF
!
!     Unrecognized code in code field. Ignore this row
!
            ELSE
               WRITE (IOLOG, 1400) NREC,QLINE
               WRITE (IOLOG, 1800) Q2,Q3
            ENDIF
            GOTO 100
!
!     COLUMN header has been read. Flush unconstrained nodes if there are any.
!
  200 CONTINUE
         IF (NECHO1 .GE.    1) WRITE (IOLOG, 1400) NREC,QLINE
         DO IP=NPSEEN+1,NXTPER-1
            CALL NEWNODE(IERR)
            IF (IERR .GT. 0) GOTO 9999
         END DO
         NPSEEN = NXTPER - 1
         IF (Q1 .EQ. QC .AND. Q2 .EQ. QO) THEN
            NEXT = 1
         ELSE
            NEXT = 2
         ENDIF
!
         IF (NPSEEN .LT. NPER) THEN
            WRITE (IOLOG,2500) NPER, NPSEEN
            NPER = NPSEEN
         ENDIF
         IF (IENDAT .EQ.    0) GOTO 9700
            NODEINFO(INODE)%IDESC   = 0
            NODEINFO(INODE)%IBROTH  = 0
            NODEINFO(INODE)%IANCTR  = INODE - 1
            IRNGE0(NPSEEN) = INODE
            GOTO 999
!
!     Come here if anything went wrong
!
 9000 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         IF (IENDAT .EQ. 0) GOTO 9020
         IF (IENDAT .EQ. 1 .OR. NPER .NE. NPSEEN) GOTO 9010
            WRITE (IOLOG, 3000)
            GOTO 999
!
 9010    CONTINUE
            WRITE (IOLOG, 3010)
            GOTO 9999
!
 9020    CONTINUE
            CORSTAT = 1
            WRITE (IOLOG, 3020)
            GOTO 999
!
 9060 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3060)
         GOTO 9999
!
 9300 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3300)
         GOTO 9999
!
 9700 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3700)
         GOTO 9999
!
 9999 CONTINUE
         CORSTAT = 2
!
!     Return to calling program
!
  999 CONTINUE
         RETURN
!
 1400 FORMAT(I8,2X,A80)
 1800 FORMAT(' XXX - WARNING - Unrecognized code =',2A1,'= in ROWS',      &
     &       ' section. Row ignored.')
 2300 FORMAT(' XXX -  FATAL  - Illegal header card in CORE file.')
 2500 FORMAT(' XXX - WARNING - Time file had',I3,' periods; ROWS',      &
     &       ' section only finds',I3,'.')
 3000 FORMAT(' XXX - WARNING - Missing ENDATA card in CORE file.')
 3010 FORMAT(' XXX -  FATAL  - Detected EOF in CORE file.')
 3020 FORMAT(' XXX - WARNING - No information in CORE file.')
 3060 FORMAT(' XXX -  FATAL  - Unable to allocate memory in INCORE.')
! 3200 FORMAT(' XXX -  FATAL  - Number of periods misspecified in ROWS',
!     *       ' or COLUMNS section.')
 3300 FORMAT(' XXX -  FATAL  - No COLUMNS section specified.')
 3700 FORMAT(' XXX -  FATAL  - ROWS section is non-existent.')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INCOLS1 ( NREC, QLINE, NEXT, IENDAT, IERR )
!
!     This routine reads the COLUMNS section in temporal form.
!
!     Routines called: CHOPEN, GETNAM, GETMPS, GETNET.
!
!     --------------------------------------------------------------
!     This version dated 21 September 2002. Written by Gus Gassmann.
!     --------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      REAL*8 ATEMP(4)
      CHARACTER*(NAMLEN) DNAME(3), DBLANK, DROW, DCOL, QLINE
      CHARACTER*8 DUNCAP, DUNDIR, DINT1, DINT2, DINT3, DINT4
      CHARACTER*1 DB255(NAMLEN), QROW(NAMLEN), QCOL(NAMLEN), T
      INTEGER*4   LENGTH(3),S3ROW,S3BLK
      LOGICAL     MPSFORM
!
      EQUIVALENCE (DBLANK,DB255), (DROW,QROW), (DCOL,QCOL)
!
      DATA DUNCAP  /'UNCAP   '/, DUNDIR   /'UNDIR   '/,      &
     &     DINT1 /'''INTORG'''/, DINT2  /'''SOSORG'''/,      &
     &     DINT3 /'''INTEND'''/, DINT4  /'''SOSEND'''/,      &
     &     DB255 /NAMLEN*' '/
! ----------------------------------------------------------------------
!
!     Allocate temporary space for constraint matrix blocks.
!
         ALLOCATE( AUX(IZANZB,NPER), IAUX(IZANZB,NPER),      &
     &            LAUX(IZCOLP,NPER), LMNS(NPER), STAT=IER)
         IF (IER .GT. 0) GOTO 9060
!
!     Initialize the columns section
!
            IENDAT = 2
            IPER   = 1
            IROW   = 0
            ICOL   = 0
            INODE  = 1
            ICOLA  = 0
            ICOST  = 0
            IDATA  = 0
            IELMA  = 0
            INAMES = 0
            IBOUND = 0
            IVTYP  = 0
            NROWS  = NODEINFO(1)%NROW
            NCOLS  = NODEINFO(1)%NROW
            IROW1  = NROWS
!
            NOFSOS = 0
!
            MNABLK = 2*NPER - 1
            IF (MNABLK .GT. IZABLK) THEN
               IF (.NOT. EXTMEM_ABLK(MNABLK,IERR)) GOTO 9060
            ENDIF
!
            DO IP=2,NPER
               KDATA(IP) = 2*IP - 3
            END DO
!
            NELEM = 0
!
!     Copy the row information for period 1 (from temporary array T_ROWS).
!
            CALL CPROWS(1,IERR)
            IF (IERR .GT. 0) GOTO 9999
            IF (NEXT .EQ. 1) GOTO 235
               GOTO 240
!
!     Come back here to process the next record
!     =========================================
!
  230    CONTINUE
            IF (MPSFORM) THEN
               CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOCOR,NREC,LENGTH)
            ELSE
               CALL GETNET(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOCOR,NREC,LENGTH)
            ENDIF
            IF (Q1 .EQ. QAST ) GOTO 230
            IF (IERR   .GT. 0) GOTO 9980
            IF (Q1 .EQ. QBL  ) GOTO 250
            IF (Q1 .EQ. QC .AND. Q2 .EQ. QO) GOTO 235
            IF (Q1 .EQ. QA .AND. Q2 .EQ. QR) GOTO 240
            IF (Q1 .EQ. QR .AND. Q2 .EQ. QH) GOTO 600
            IF (Q1 .EQ. QD .AND. Q2 .EQ. QE) GOTO 610
            IF (Q1 .EQ. QS .AND. Q2 .EQ. QU) GOTO 620
            IF (Q1 .EQ. QB .AND. Q2 .EQ. QO) GOTO 630
            IF (Q1 .EQ. QR .AND. Q2 .EQ. QA) GOTO 640
            IF (Q1 .EQ. QQ .AND. Q2 .EQ. QS) GOTO 642
            IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 645
            WRITE (IOLOG, 1400) NREC,QLINE
               WRITE (IOLOG, 2300)
               GOTO 9999
!
!     Header starting a new COLUMNS segment
!
  235    CONTINUE
            MPSFORM = .TRUE.
            INTTYP = 0
            DEFUPB = DEFUB
            DEFLOB = DEFLB
            GOTO 230
!
!     Header starting a new ARCS segment
!
  240    CONTINUE
            MPSFORM = .FALSE.
            INTTYP  = 0
            IF     (DNAME(1) .EQ. DUNCAP .OR. DNAME(2) .EQ. DUNCAP) THEN
               DEFUPB = PLINF
               DEFLOB = 0.D0
            ELSEIF (DNAME(1) .EQ. DUNDIR .OR. DNAME(2) .EQ. DUNDIR) THEN
               DEFUPB =  PLINF
               DEFLOB = -PLINF
            ELSE
               DEFUPB = 0.D0
               DEFLOB = 0.D0
            ENDIF
            GOTO 230
!
!    Here we have a regular data card.
!    =================================
!
!    We must determine whether it starts a new period, whether it continues
!    a column already started, whether it begins or ends a segment of
!    integer variables, and whether it even contains relevant data
!
  250    CONTINUE
            IF (NECHO1 .GE. 3) WRITE (IOLOG, 1400) NREC,QLINE
            IF (MPSFORM) THEN
!
!     COLUMNS section.
!     ----------------
!     First match the column name to the last record read. Then check for
!     special markers to indicate integer variables and special ordered sets.
!
               JNM = 2
               IF (DNAME(1) .EQ. DCOL) THEN
                  IF (DABS(ATEMP(1)) .GT. ZTOLIN) GOTO 430
                  IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 230
                      JNM = 3
                      ATEMP(1) = ATEMP(2)
                      GO TO 430
               ELSEIF (DNAME(3) .EQ. DINT1) THEN
                  INTTYP = 1
                  GOTO 230
               ELSEIF (DNAME(3) .EQ. DINT2) THEN
                  IF (NOFSOS .GE. IZSTYP) THEN
                     IF (.NOT. EXTMEM_STYP(NOFSOS+1,IERR)) GOTO 9060
                  ENDIF
!
                  NOFSOS = NOFSOS + 1
                  INTTYP = -NOFSOS
                  IF     (Q2 .EQ. QS .AND. Q3 .EQ. '1') THEN
                     SOSTYP(NOFSOS) = 1
                  ELSEIF (Q2 .EQ. QS .AND. Q3 .EQ. '2') THEN
                     SOSTYP(NOFSOS) = 2
                  ELSE
!
!     SOS type 3. The first name field must contain a valid row name.
!     Each column participating in this set must be 0/1, the row must be
!     an equality row, with rhs = 1, and coefficients in this row are 1.
!
                     DCOL = DNAME(1)
                     DO JPER=1,NPER
                        KN = NODEINFO(JPER)%KNAMES
                        NN = NODEINFO(JPER)%NROW
                        IRES = NMSRCH(DCOL,KN+1,KN+NN)
                        IF (IRES .GT. 0) THEN
                           IF (T_ROWS(IRES)%ROWTYP .NE. 0) GOTO 258
                              SOSTYP(NOFSOS) = -IRES
                              S3BLK = JPER
                              S3ROW = IRES - KN
                              DEFUPB = 1.D0
                              DEFLOB = 0.D0
                              GOTO 230
                        ENDIF
                     END DO
!
!     Could be a linking constraint...
!
                     IF (ALLOW_GC) THEN
                        DO JPER=1,IPER-1
                           KN = NODEINFO(JPER)%KNAMES
                           NN = NODEINFO(JPER)%NROW
                           IRES = NMSRCH(DCOL,KN+1,KN+NN)
                           IF (IRES .GT. 0) THEN
                              IF (T_ROWS(IRES)%ROWTYP .NE. 0) GOTO 258
                                 SOSTYP(NOFSOS) = -IRES
                                 S3BLK = JPER
                                 S3ROW = IRES - KN
                                 DEFUPB = 1.D0
                                 DEFLOB = 0.D0
                                 GOTO 230
                           ENDIF
                        END DO
                     ENDIF
!
  258             CONTINUE
                     WRITE (IOLOG, 1400) NREC,QLINE
                     WRITE (IOLOG, 1900)
                     GOTO 9999
                  ENDIF
                  GOTO 230
!
               ELSEIF (DNAME(3) .EQ. DINT3 .OR.      &
     &                 DNAME(3) .EQ. DINT4) THEN
                  INTTYP = 0
                  DEFUPB = DEFUB
                  DEFLOB = DEFLB
                  GOTO 230
!
!     Column name not previously encountered.
!     Columns are named explicitly, and they appear in temporal order.
!     To find out if a new column still belongs to the current period,
!     we first check if the column name matches the first column of any
!     future period. If it does not, it must belong to the current period.
!
               ELSE
                  DO IPNEXT=IPER+1,NPER
                     IF (DNAME(1) .EQ. DTIMEC(IPNEXT))      &
     &                  GOTO 300
                  END DO
                  IF (DABS(ATEMP(1)) .GT. ZTOLIN) GOTO 370
                  IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 230
                      JNM = 3
                      ATEMP(1) = ATEMP(2)
                      GO TO 370
               ENDIF
            ELSE
!
!     ARCS section
!     ------------
!     Arcs are not named, so we infer the time period as follows:
!     Check if the third name field contains an explicit period name.
!     If not, use the earlier of the periods associated with the FROM-node
!     (IPFROM) and the TO-node (IPTO).
!
               PURELP = .FALSE.
               INTTYP = 0
               JNM    = 2
               IF (Q2 .NE. QBL) THEN
                  T = Q2
               ELSE
                  T = Q3
               ENDIF
               IF ((DNAME(1) .EQ. DCOL) .AND.      &
     &             (DNAME(2) .EQ. DROW) .AND. (T .EQ. QCODE)) THEN
                  GOTO 500
               ELSE
!
!     Determine FROM-node or use previous info
!
                  IF (DNAME(1) .NE. DCOL) THEN
                     DCOL = DNAME(1)
                     DO IPFROM=IPER,NPER
                        KN = NODEINFO(IPFROM)%KNAMES
                        NN = NODEINFO(IPFROM)%NROW
                        IRES = NMSRCH(DCOL,KN+1,KN+NN)
                        IF (IRES .GT. 0) GOTO 270
                     END DO
!
                     WRITE (IOLOG, 1400) NREC,QLINE
                     WRITE (IOLOG, 3800)
                     CORSTAT = 2
                     GOTO 999
!
  270                CONTINUE
                        IFROM = IRES - KN
                  ENDIF
!
!     Determine TO-node or use previous info
!
                  IF (DNAME(2) .NE. DROW) THEN
                     DROW = DNAME(2)
                     DO IPTO=IPER,NPER
                        KN = NODEINFO(IPTO)%KNAMES
                        NN = NODEINFO(IPTO)%NROW
                        IRES = NMSRCH(DROW,KN+1,KN+NN)
                        IF (IRES .GT. 0) GOTO 280
                     END DO
!
                     WRITE (IOLOG, 1400) NREC,QLINE
                     WRITE (IOLOG, 3900)
                     CORSTAT = 2
                     GOTO 999
!
  280                CONTINUE
                        ITO = IRES - KN
                  ENDIF
!
!     Check the code field for a marker to indicate parallel arcs
!
                  IF (Q2 .NE. QBL) THEN
                     QCODE = Q2
                  ELSE
                     QCODE = Q3
                  ENDIF
!
!     Finally check for an explicit period name in the third name field.
!     Due to the time order of the core file, only two cases are legal:
!     the arc belongs to the current period, or it belongs to the next period.
!
                  IF (DNAME(3) .EQ. DTIME(IPER)  ) THEN
                     DNAME(3) = DBLANK
                     GOTO 370
                  ENDIF
                  IF (DNAME(3) .EQ. DTIME(IPER+1)) THEN
                     DNAME(3) = DBLANK
                     GOTO 300
                  ENDIF
                  IPNEXT = MIN0(IPFROM,IPTO)
                  IF (IPNEXT .EQ. IPER) GOTO 370
                  IF (IPNEXT .GT. IPER) GOTO 300
                     WRITE (IOLOG, 1400) NREC, QLINE
                     WRITE (IOLOG, 3950)
                     CORSTAT = 2
                     GOTO 999
               ENDIF
            ENDIF
!
!     We should never, ever get here! Sum'pin's wrong!
!
            WRITE (IOLOG, 1400) NREC, QLINE
            WRITE (IOLOG, 3750)
            GOTO 9999
!
!     A new period is coming up. (The new period is IPER+1; old period is IPER)
!
  300 CONTINUE
         IF (IPER .GE. NPER) GOTO 9200
            NODEINFO(IPER)%NCOL = NCOLS
            NELMA(IDATA+1) = NELEM
            LASTA  = LASTA + NELEM
            LASTCA = LASTCA + NCOLS + 1 - NODEINFO(IPER)%NROW
!
!     Place the off-diagonal blocks for the current (old) node)
!
             CALL CPBLKS ( IERR )
             IF (IERR .GT. 0) GOTO 9999
!
!     Initial the theta column (of the current period)
!
  340    CONTINUE
            XINFTY = 1.D31
            KBD = NODEINFO(IPER)%KBOUND + NCOLS + 1
            XLB(KBD) = -XINFTY
            XUB(KBD) =  XINFTY
!
!     Now set the pointer values for the new period
!
            IPREV = IPER
            IPER  = IPER  + 1
            INODE = INODE + 1
            IDATA = KDATA(INODE)
            KCOLA(IDATA+1) = LASTCA
            KELMA(IDATA+1) = LASTA
            NODEINFO(IPER)%KCOL   = NODEINFO(IPREV)%KCOL   + NCOLS + 1
            NODEINFO(IPER)%KCOST  = NODEINFO(IPREV)%KCOST  + NCOLS      &
     &                            - NODEINFO(IPREV)%NROW
            NODEINFO(IPER)%KBOUND = NODEINFO(IPREV)%KBOUND + NCOLS + 1
            NODEINFO(IPER)%KNAMES = NODEINFO(IPREV)%KNAMES + NCOLS
            NODEINFO(IPER)%VARPTR = NODEINFO(IPREV)%VARPTR + NCOLS
!
!     Initial the new node
!
            ICOL   = NODEINFO(IPER)%KCOL
            IROW   = NODEINFO(IPER)%KROW
            ICOST  = NODEINFO(IPER)%KCOST
            INAMES = NODEINFO(IPER)%KNAMES
            IBOUND = NODEINFO(IPER)%KBOUND
            IVTYP  = NODEINFO(IPER)%VARPTR
            NROWS  = NODEINFO(IPER)%NROW
            NCOLS  = NODEINFO(IPER)%NROW
            ICOLA  = LASTCA
            IELMA  = LASTA
            NELEM  = 0
            LA(ICOLA+1) = 1
            DO JJ=1,NPER
               LMNS(JJ) = 0
               LAUX(1,JJ) = 1
            END DO
!
!     Extend memory if necessary
!
            IF (IROW + NROWS .GT. IZROWS) THEN
               IF (.NOT. RESIZE .OR. (IROW + NROWS .GE. MXROWS)      &
     &                          .OR. (IZROWS .GE. MXROWS)) GOTO 9060
               IZROWS = MAX(IROW + NROWS, MIN(MXROWS, 2*IZROWS))
            ENDIF
!
            MNVNAM = INAMES + NROWS + 1
            IF (MNVNAM .GT. IZVNAM) THEN
               IF (.NOT. EXTMEM_VNAM(MNVNAM,IERR)) GOTO 9060
            ENDIF
!
            MNBNDS = IBOUND + NROWS
            IF (MNBNDS .GT. IZBNDS) THEN
               IF (.NOT. EXTMEM_BNDS(MNBNDS,IERR)) GOTO 9060
            ENDIF
!
            MNVTYP = IVTYP  + NROWS
            IF (MNVTYP .GT. IZVTYP) THEN
               IF (.NOT. EXTMEM_VTYP(MNVTYP,IERR)) GOTO 9060
            ENDIF
!
            NCHARS = 0
            DO J=1,NROWS
               NCHARS = NCHARS + T_ROWS(IROW+J)%LENGTH
            END DO
            MNCHAR = LASTCH + NCHARS
            IF (MNCHAR .GT. IZCHAR) THEN
               IF (.NOT. EXTMEM_CHAR(MNCHAR,IERR)) GOTO 9060
            ENDIF
!
!     Copy row information for the new node (from temporary array T_ROWS).
!
            DO J=1,NROWS
               LENR            = T_ROWS(IROW+J)%LENGTH
               DROW            = T_ROWS(IROW+J)%NAME
               VARTYP(IVTYP+J) = T_ROWS(IROW+J)%ROWTYP
               XLB(IBOUND+J)   = T_ROWS(IROW+J)%LOWER
               XUB(IBOUND+J)   = T_ROWS(IROW+J)%UPPER
               NAME1(INAMES+J) = LASTCH + 1
               DO I=1,LENR
                  QNAMES(I+LASTCH) = QROW(I)
               END DO
               LASTCH = LASTCH + LENR
            END DO
            NAME1(INAMES+NROWS+1) = LASTCH + 1
!
            IF (IPER   .LT. IPNEXT) GOTO 300
            IF (IENDAT .GE.      3) GOTO 650
!
!     Start a new variable. In a COLUMNS section, we have the name directly
!
  370    CONTINUE
            IF (MPSFORM) THEN
               MNCHAR = LASTCH + LENGTH(1)
               IF (MNCHAR .GT. IZCHAR) THEN
                  IF (.NOT. EXTMEM_CHAR(MNCHAR,IERR)) GOTO 9060
               ENDIF
!
               DCOL = DNAME(1)
               DO I=1,LENGTH(1)
                  QNAMES(LASTCH+I) = QCOL(I)
               END DO
               LASTCH = LASTCH + LENGTH(1)
!
!     In an ARCS section, manufacture a column name:
!                                "{From-node}->{To-node}[:{code}]"
!
            ELSE
               IF (QCODE .EQ. QBL) THEN
                  MNCHAR = LASTCH + LENGTH(1) + LENGTH(2) + 2
               ELSE
                  MNCHAR = LASTCH + LENGTH(1) + LENGTH(2) + 4
               ENDIF
               IF (MNCHAR .GT. IZCHAR) THEN
               IF (.NOT. EXTMEM_CHAR(MNCHAR,IERR)) GOTO 9060
               ENDIF
!
               DO I=1,LENGTH(1)
                  QNAMES(LASTCH+I) = QCOL(I)
               END DO
               QNAMES(LASTCH+LENGTH(1)+1) = '-'
               QNAMES(LASTCH+LENGTH(1)+2) = '>'
               LASTCH = LASTCH + LENGTH(1) + 2
               DO I=1,LENGTH(2)
                  QNAMES(LASTCH+I) = QROW(I)
               END DO
               LASTCH = LASTCH + LENGTH(2)
               IF (QCODE .NE. QBL) THEN
                  QNAMES(LASTCH+1) = ':'
                  QNAMES(LASTCH+2) = QCODE
                  LASTCH = LASTCH + 2
               ENDIF
            ENDIF
!
!     Check array sizes and expand if necessary
!
            NCOLS  = NCOLS  + 1
            ICCA   = ICOLA  + NCOLS - NROWS
            JCOST  = ICOST  + NCOLS - NROWS
            JCOL   = ICOL   + NCOLS
            JNAMES = INAMES + NCOLS + 1
            JBOUND = IBOUND + NCOLS
            JVTYP  = IVTYP  + NCOLS
            LCC    = NCOLS  - NROWS
!
            IF (NCOLS .GT. IZCOLP) THEN
               IF (.NOT. EXTMEM_COLP(NCOLS,IERR)) GOTO 9060
            ENDIF
!
            IF (ICCA .GE. IZACOL) THEN
               IF (.NOT. EXTMEM_ACOL(ICCA+1,IERR)) GOTO 9060
            ENDIF
!
            IF (JCOST .GT. IZCOST) THEN
               IF (.NOT. EXTMEM_COST(JCOST,IERR)) GOTO 9060
            ENDIF
!
            IF (JCOL .GT. IZCOLS) THEN
               IF (.NOT. RESIZE      &
     &              .OR. (IZCOLS .GE. MXCOLS)      &
     &              .OR. (JCOL   .GT. MXCOLS)) GOTO 9060
               IZCOLS = MAX(ICOL,MIN(MXCOLS, 2*IZCOLS))
            ENDIF
!
            IF (JNAMES .GT. IZVNAM) THEN
               IF (.NOT. EXTMEM_VNAM(JNAMES,IERR)) GOTO 9060
            ENDIF
!
            IF (JBOUND .GT. IZBNDS) THEN
               IF (.NOT. EXTMEM_BNDS(JBOUND,IERR)) GOTO 9060
            ENDIF
!
            IF (JVTYP  .GT. IZVTYP) THEN
               IF (.NOT. EXTMEM_VTYP(JVTYP,IERR)) GOTO 9060
            ENDIF
!
!     Set column pointers and defaults for costs and bounds
!
            NAME1(JNAMES)        = LASTCH + 1
            COST(JCOST)          = 0.D0
            XLB(JBOUND)          = DEFLOB
            XUB(JBOUND)          = DEFUPB
            LA(ICCA)             = NELEM + 1
            LA(ICCA+1)           = NELEM + 1
!
            DO JP=IPER+1,NPER
                  LAUX(LCC+1,JP) = LMNS(JP) + 1
            END DO
            IF (ALLOW_GC) THEN
               DO JP=1,IPER-1
                  LAUX(LCC+1,JP) = LMNS(JP) + 1
               END DO
            ENDIF
!
!     Store the column type and, if SOS type 3, the implied matrix coefficient
!
            IF (INTTYP .GE. 0 .OR. SOSTYP(NOFSOS) .GT. 0) THEN
               VARTYP(JVTYP) = INTTYP
            ELSE
               VARTYP(JVTYP) = 1
               IF (S3BLK .EQ. IPER) THEN
                  NELEM = NELEM + 1
                  JELMA = IELMA + NELEM
                  IF (JELMA .GT. IZALMN) THEN
                     IF (.NOT. EXTMEM_ALMN(JELMA,IERR)) GOTO 9060
                  ENDIF
!
                  IA(JELMA) = S3ROW
                   A(JELMA) = 1.D0
!
!     Here the implied coefficient is in an off-diagonal block
!
               ELSE
                  JELEM = LMNS(S3BLK)
                  IF (JELEM .GT. IZANZB) THEN
                     IF (.NOT. EXTMEM_ANZB(JELEM,IERR)) GOTO 9060
                  ENDIF
!
                  LMNS(S3BLK)               = JELEM + 1
                  LAUX(NCOLS-NROWS+1,S3BLK) = JELEM + 2
                  IAUX(JELEM,        S3BLK) = S3ROW
                   AUX(JELEM,        S3BLK) = 1.D0
               ENDIF
            ENDIF
!
!        Test for row match
!
  430    CONTINUE
         IF (.NOT. MPSFORM) GOTO 500
         DROW = DNAME(JNM)
         IRES = NMSRCH(DROW,INAMES+1,INAMES+NROWS)
         IF (IRES .GT. 0) THEN
!
!    Matched a coefficient in the A-matrix. Expand array if needed, then store.
!
            I = IRES - INAMES
            IF (I .NE. IOBJ) THEN
               JELMA = IELMA + NELEM + 1
               IF (JELMA .GT. IZALMN) THEN
                  IF (.NOT. EXTMEM_ALMN(JELMA,IERR)) GOTO 9060
               ENDIF
!
               NELEM     = NELEM + 1
               IA(JELMA) = I
                A(JELMA) = ATEMP(1)
               LA(ICOLA+NCOLS-NROWS+1) = NELEM + 1
!
!     Cost coefficients (even if fixed) are not stored in the A-matrix
!
            ELSE
               COST(JCOST) = ATEMP(1)
            ENDIF
            GOTO 490
!
!     Search other time periods if not found
!
         ELSE
            DO JPER=IPER+1,NPER
               JROWS = NODEINFO(JPER)%NROW
               IROW0 = NODEINFO(JPER)%KROW
               DO IRES=2,JROWS
                  IF (DROW .EQ. T_ROWS(IROW0+IRES)%NAME) GOTO 470
               END DO
            END DO
!
!     Maybe it is a linking constraint or an alternative objective row?
!
            DO JPER=1,IPER-1
               J0 = NODEINFO(JPER)%KNAMES
               JR = NODEINFO(JPER)%NROW
               IRES = NMSRCH(DROW,J0+2,J0+JR)
               IF (IRES .GT. 0) THEN
                  IRES = IRES - J0
                  IF (.NOT. ALLOW_GC) THEN
                     WRITE (IOLOG, 1400) NREC,QLINE
                     WRITE (IOLOG, 2600) DROW
                  ENDIF
                  GOTO 470
               ENDIF
            END DO
!
!     No match!
!
  460       CONTINUE
               WRITE (IOLOG, 1400) NREC,QLINE
               WRITE (IOLOG, 1100) DROW
               GOTO 9999
!
!     We have found an element in an off-diagonal block
!
  470 CONTINUE
               JELEM = LMNS(JPER)
               IF (JELEM .GE. IZANZB) THEN
                  IF (.NOT. EXTMEM_ANZB(JELEM,IERR)) GOTO 9060
               ENDIF
!
               LMNS(JPER)               = JELEM + 1
               LAUX(NCOLS+1-NROWS,JPER) = JELEM + 2
               IAUX(LMNS(JPER),   JPER) = IRES
                AUX(LMNS(JPER),   JPER) = ATEMP(1)
         ENDIF
!
!     There may be another coefficient on this card. If not, get the next one
!
  490    CONTINUE
            IF (JNM .EQ. 3) GOTO 230
            IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 230
               JNM = 3
               ATEMP(1) = ATEMP(2)
               GOTO 430
!
!        Data lines in ARCS segment define arcs or side constraints,
!        depending on the value of the third name field.
!
  500    CONTINUE
         IF (DNAME(3) .NE. DBLANK) GOTO 530
!
!     First store the cost coefficient and the bounds (if specified)
!
            COST(JCOST) = ATEMP(1)
            IF (DABS(ATEMP(2)) .GT. ZTOLIN)      &
     &                              XUB(IBOUND+NCOLS) = ATEMP(2)
            IF (DABS(ATEMP(3)) .GT. ZTOLIN)      &
     &                              XLB(IBOUND+NCOLS) = ATEMP(3)
!
!     The coefficient in the FROM-node is -1.D0
!
            IF (IPFROM .EQ. IPER) THEN
               JELMA = NELEM + IELMA
               IF (JELMA .GE. IZALMN) THEN
                  IF (.NOT. EXTMEM_ALMN(JELMA,IERR)) GOTO 9060
               ENDIF
!
               IA(IELMA+NELEM+1) = IFROM
                A(IELMA+NELEM+1) = -1.D0
               NELEM = NELEM + 1
               LA(ICOLA+NCOLS-NROWS+1) = NELEM + 1
!
            ELSE
               JELEM = LMNS(JPER)
               IF (JELEM .GE. IZANZB) THEN
                  IF (.NOT. EXTMEM_ANZB(JELEM,IERR)) GOTO 9060
               ENDIF
!
               LMNS(JPER)               = JELEM + 1
               LAUX(NCOLS+1-NROWS,JPER) = JELEM + 2
               IAUX(LMNS(JPER),   JPER) = IFROM
                AUX(LMNS(JPER),   JPER) = -1.D0
            ENDIF
!
!      The coefficient in the TO-node could be a multiplier
!
            IF (IPTO .EQ. IPER) THEN
               JELMA = NELEM + IELMA
               IF (JELMA .GE. IZALMN) THEN
                  IF (.NOT. EXTMEM_ALMN(JELMA,IERR)) GOTO 9060
               ENDIF
!
                  IA(IELMA+NELEM+1) = ITO
                  IF (DABS(ATEMP(4)) .GT. ZTOLIN) THEN
                     A(IELMA+NELEM+1) = ATEMP(4)
                  ELSE
                     A(IELMA+NELEM+1) = 1.0
                  ENDIF
                  NELEM = NELEM + 1
                  LA(ICOLA+NCOLS-NROWS+1) = NELEM + 1
            ELSE
               JELEM = LMNS(JPER)
               IF (JELEM .GE. IZANZB) THEN
                  IF (.NOT. EXTMEM_ANZB(JELEM,IERR)) GOTO 9060
               ENDIF
!
                  LMNS(JPER)               = LMNS(JPER) + 1
                  LAUX(NCOLS+1-NROWS,JPER) = LMNS(JPER) + 2
                  IAUX(LMNS(JPER),   JPER) = ITO
                  IF (DABS(ATEMP(4)) .GT. ZTOLIN) THEN
                     AUX(LMNS(JPER), JPER) = ATEMP(4)
                  ELSE
                     AUX(LMNS(JPER), JPER) = 1.0
                  ENDIF
              ENDIF
              GOTO 230
!
!        Here we have a side constraint. Identify the row.
!
  530    CONTINUE
         DROW = DNAME(3)
         IRES = NMSRCH(DROW,INAMES+1,INAMES+NROWS)
         IF (IRES .GT. 0) THEN
!
!    Matched a coefficient in the A-matrix
!
            I = IRES - INAMES
            IF (I .NE. IOBJ) THEN
               JELMA = NELEM + IELMA
               IF (JELMA .GE. IZALMN) THEN
                  IF (.NOT. EXTMEM_ALMN(JELMA,IERR)) GOTO 9060
               ENDIF
!
                  NELEM           = NELEM + 1
                  IA(IELMA+NELEM) = I
                   A(IELMA+NELEM) = ATEMP(1)
                  LA(ICOLA+NCOLS-NROWS+1) = NELEM + 1
!
!     Error. Cost coefficient should be specified on a pure arc card
!
            ELSE
               WRITE (IOLOG, 1400) NREC,QLINE
               WRITE (IOLOG, 1600)
            ENDIF
!
!     Search other time periods
!
         ELSE
            DO JPER=IPER+1,NPER
               JROWS = NODEINFO(JPER)%NROW
               IROW0 = NODEINFO(JPER)%KROW
               DO IRES=2,JROWS
                  IF (DROW .EQ. T_ROWS(IROW0+IRES)%NAME) GOTO 590
               END DO
            END DO
!
!     Perhaps we have an alternative objective row or a linking constraint?
!
            DO JPER=1,IPER-1
               JROWS = NODEINFO(JPER)%NROW
               IROW1 = NODEINFO(JPER)%KNAMES
               IRES = NMSRCH(DROW,IROW1+1,IROW1+JROWS)
               IF (IRES .GT. 0) THEN
                  IRES = IRES - JROWS
                  IF (.NOT. ALLOW_GC) THEN
                     WRITE (IOLOG, 1400) NREC,QLINE
                     WRITE (IOLOG, 2600) DROW
                  ENDIF
                  GOTO 590
               ENDIF
            END DO
!
!     No match!
!
  570       CONTINUE
               WRITE (IOLOG, 1400) NREC,QLINE
               WRITE (IOLOG, 1100) DROW
               GOTO 9999
         ENDIF
         GOTO 230
!
!     We have found an element in an off-diagonal block
!
  590          CONTINUE
                  I = IRES - IROW1
               JELEM = LMNS(JPER)
               IF (JELEM .GE. IZANZB) THEN
                  IF (.NOT. EXTMEM_ANZB(JELEM,IERR)) GOTO 9060
               ENDIF
                      LMNS(JPER)               = LMNS(JPER) + 1
                      LAUX(NCOLS+1-NROWS,JPER) = LMNS(JPER) + 2
                      IAUX(LMNS(JPER),   JPER) = IRES
                       AUX(LMNS(JPER),   JPER) = ATEMP(1)
                      GOTO 230
!
!     The COLUMNS section is done. What is next?
!
!     RHS section
!
  600    CONTINUE
            NEXT = 1
            GOTO 650
!
!     DEMAND section
!
  610    CONTINUE
            NEXT = 2
            GOTO 650
!
!     SUPPLY section
!
  620    CONTINUE
            NEXT = 3
            GOTO 650
!
!     BOUNDS section
!
  630    CONTINUE
            NEXT = 4
            GOTO 650
!
!     RANGES section
!
  640    CONTINUE
            NEXT = 5
            GOTO 650
!
!     In-line QSECT section
!
  642    CONTINUE
            NEXT = 6
            GOTO 650
!
!     ENDATA record
!
  645    CONTINUE
            NEXT = 7
!
!     Before processing the other sections, we initialize the bounds for
!     the slacks. (This cannot be done properly until we know how many
!     columns there are and where in the XLB array we have the slacks.)
!     Also check if the last stage is empty and copy any remaining
!     linking constraints.
!
  650 CONTINUE
            NODEINFO(IPER)%NCOL = NCOLS
            NELMA(IDATA+1) = NELEM
            LASTA  = LASTA + NELEM
            LASTCA = LASTCA +  1  + NCOLS  - NROWS
!
            CALL CPBLKS(IERR)
            IF (IERR .GT. 0) GOTO 9999
!
            NROWS = NODEINFO(IPER)%NROW
            LASTC  = NODEINFO(IPER)%KCOST  + NCOLS - NROWS
            LASTR  = NODEINFO(IPER)%KRHS   + NROWS
            LASTVP = NODEINFO(IPER)%VARPTR + NCOLS
            LASTBD = NODEINFO(IPER)%KBOUND + NCOLS + 1
            LASTNM = NODEINFO(IPER)%KNAMES + NCOLS + 1
!
            IF (MARKOV) THEN
               IF (.NOT. HAVE_GC) THEN
                  LASTBA = 2 * IPER - 1
               ELSE
                  LASTBA = IPER * (IPER+1)/2 + IPER - 1
               ENDIF
            ELSE
               IF (.NOT. HAVE_GC) THEN
                  LASTBA = IPER * (IPER+1)/2
               ELSE
                  LASTBA = IPER**2
               ENDIF
            ENDIF
!
!     Initial the theta column as well...
!
            XLB(LASTBD) = -XINFTY
            XUB(LASTBD) =  XINFTY
!
            IENDAT = 3
            IF (IPER .LT. NPER) THEN
               IPNEXT = NPER
               GOTO 300
            ENDIF
         GOTO 999
!
!     Come here if anything went wrong
!
 9060 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3060)
         GOTO 9999
!
 9200 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3200)
         GOTO 9999
!
 9980 CONTINUE
         WRITE (IOLOG, 3980)
         GOTO 9999
!
 9999 CONTINUE
         CORSTAT = 2
!
!     Return to calling program
!
  999 CONTINUE
!
         DEALLOCATE(AUX,IAUX,LAUX,LMNS,STAT=IER)
         IF (IER .NE. 0) THEN
            WRITE (IOLOG, 3070)
            CORSTAT = 2
         ENDIF
         RETURN
!
 1100 FORMAT(' XXX -  FATAL  - Row "',A20,'" violates time order in',      &
     &       ' core file.')
 1400 FORMAT(I8,2X,A80)
 1600 FORMAT(' XXX - WARNING - Cost coefficient should be specified',      &
     &       ' on a pure arc card. Data ignored.')
 1900 FORMAT(' XXX -  FATAL  - Inappropriate row name in SOS type 3.')
 2300 FORMAT(' XXX -  FATAL  - Illegal header card in CORE file.')
 2600 FORMAT(' XXX - WARNING - Row "',A20,'" violates time order.',      &
     &       ' Enable linking constraints.')
 3060 FORMAT(' XXX -  FATAL  - Unable to allocate memory in INCORE')
 3070 FORMAT(' XXX -  FATAL  - Memory handler returned error in INCORE')
 3200 FORMAT(' XXX -  FATAL  - Number of periods misspecified in ROWS',      &
     &       ' or COLUMNS section.')
 3750 FORMAT(' XXX -  FATAL  - Error in COLUMNS section.')
 3800 FORMAT(' XXX -  FATAL  - FROM-node not previously specified.')
 3900 FORMAT(' XXX -  FATAL  - TO-node not previously specified.')
 3950 FORMAT(' XXX -  FATAL  - ARCS section violates temporal order.')
 3980 FORMAT(' XXX -  FATAL  - Error while reading CORE file.')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INROWS2 ( NREC, QLINE, IRTEMP, NEXT, IENDAT, IERR )
!
!     This routine reads the ROWS section in temporal form.
!
!     Routines called: CHOPEN, GETNAM, GETMPS, GETNET.
!
!     --------------------------------------------------------------
!     This version dated 21 September 2002. Written by Gus Gassmann.
!     --------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      REAL*8 ATEMP(4)
      CHARACTER*(NAMLEN) DNAME(3), DBLANK, QLINE
      CHARACTER*1 DB255(NAMLEN)
      INTEGER*4   LENGTH(3)
!
      EQUIVALENCE (DBLANK,DB255)
!
      DATA DB255 /NAMLEN*' '/
!
      ALLOCATE (TCNROW(NPER), STAT=IERR)
      IF (IERR .NE. 0) GOTO 9060
! ----------------------------------------------------------------------
!
!     Initializations.
!
         IENDAT = 0
!
         DO I=1,NPER
            TCNROW(I) = 0
         END DO
!
!     Read the next record, then distribute for further processing.
!     =============================================================
!
  100 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOCOR,NREC,LENGTH)
         IF (Q1 .EQ. QAST) GOTO 100
         IF (IERR .GT. 0 ) GOTO 9000
         IF (Q1 .EQ. QBL ) GOTO 120
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QO) GOTO 105
         IF (Q1 .EQ. QN .AND. Q2 .EQ. QO) GOTO 105
         IF (Q1 .EQ. QA .AND. Q2 .EQ. QR) GOTO 200
         IF (Q1 .EQ. QC .AND. Q2 .EQ. QO) GOTO 200
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QH) GOTO 9300
         IF (Q1 .EQ. QS .AND. Q2 .EQ. QU) GOTO 9300
         IF (Q1 .EQ. QD .AND. Q2 .EQ. QE) GOTO 9300
         IF (Q1 .EQ. QB .AND. Q2 .EQ. QO) GOTO 9300
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QA) GOTO 9300
         IF (Q1 .EQ. QQ .AND. Q2 .EQ. QS) GOTO 9300
         IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 9300
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 2300)
            GOTO 9999
!
!     ROW or NODE card. The core file is non-trivial.
!
  105    CONTINUE
            IF (IENDAT .EQ. 0) THEN
               IENDAT = 1
               MAXROW = 1
               IRT    = 1
            ENDIF
            IF (NECHO1 .GE. 1) WRITE (IOLOG, 1400) NREC,QLINE
            GOTO 100
!
!     Here we have a regular data card. Identify its period.
!
  120    CONTINUE
            IF (NECHO1 .GE. 3) WRITE (IOLOG, 1400) NREC,QLINE
            DO J=IRT+1,IRTEMP
               IF (DNAME(1) .EQ. TROWS(J)%NAME) GOTO 135
            END DO
            DO J=1,IRT
               IF (DNAME(1) .EQ. TROWS(J)%NAME) GOTO 135
            END DO
            GOTO 9160
!
!     Test row type. Expand memory if necessary.
!
  135    CONTINUE
            IRT = J
            INODE = TROWS(IRT)%PERIOD
            NROWS = TCNROW(INODE)
            IF (NROWS  .GT. IZROWP) THEN
               IF (.NOT. RESIZE .OR. (IZROWP .GE. MXROWP)      &
     &                          .OR. (NROWS  .GT. MXROWP)) GOTO 9060
               IZROWP = NROWS
            ENDIF
            IF (MAXROW .GT. IZROWS) THEN
               IF (.NOT. RESIZE .OR. (IZROWS .GE. MXROWS)      &
     &                          .OR. (MAXROW .GT. MXROWS)) GOTO 9060
               IZROWS = MAXROW
            ENDIF
!
            IF (MAXROW .GE. IZVNAM) THEN
               IF (.NOT. EXTMEM_VNAM(MAXROW+1,IERR)) GOTO 9060
            ENDIF
!
            MNCOLP = MAXROW + 1
            IF (MNCOLP .GT. IZCOLP) THEN
               IF (.NOT. EXTMEM_COLP(MNCOLP,IERR)) GOTO 9060
            ENDIF
!
            IF ((Q2 .EQ. QE) .OR. (Q3 .EQ. QE) .OR.      &
     &         ((Q2 .EQ. QBL).AND.(Q3 .EQ. QBL))) THEN
               MAXROW = MAXROW + 1
               NROWS  = NROWS  + 1
               TCNROW(INODE) = NROWS
               T_VAR(NROWS,INODE)%NAME   = DNAME(1)
               T_VAR(NROWS,INODE)%LENGTH = LENGTH(1)
               T_VAR(NROWS,INODE)%ROWTYP = 0
               T_VAR(NROWS,INODE)%LOWER  = 0.D0
               T_VAR(NROWS,INODE)%UPPER  = 0.D0
!
            ELSEIF  ((Q2 .EQ. QG) .OR. (Q3 .EQ. QG)) THEN
               MAXROW = MAXROW + 1
               NROWS  = NROWS  + 1
               TCNROW(INODE) = NROWS
               T_VAR(NROWS,INODE)%NAME   = DNAME(1)
               T_VAR(NROWS,INODE)%LENGTH = LENGTH(1)
               T_VAR(NROWS,INODE)%ROWTYP = -1
               T_VAR(NROWS,INODE)%LOWER  = -PLINF
               T_VAR(NROWS,INODE)%UPPER  = 0.D0
!
            ELSEIF  ((Q2 .EQ. QL) .OR. (Q3 .EQ. QL)) THEN
               MAXROW = MAXROW + 1
               NROWS  = NROWS  + 1
               TCNROW(INODE) = NROWS
               T_VAR(NROWS,INODE)%NAME   = DNAME(1)
               T_VAR(NROWS,INODE)%LENGTH = LENGTH(1)
               T_VAR(NROWS,INODE)%ROWTYP = 1
               T_VAR(NROWS,INODE)%LOWER  = 0.D0
               T_VAR(NROWS,INODE)%UPPER  = PLINF
!
            ELSEIF  ((Q2 .EQ. QN) .OR. (Q3 .EQ. QN)) THEN
!
!     Unbounded row could be the objective. If not, treat like any other row.
!
               IF (OBJNAM .EQ. DOTS) THEN
                  OBJNAM = DNAME(1)
                  LENOBJ = MAX(LENGTH(1),8)
               ELSEIF (OBJNAM .NE. DNAME(1)) THEN
                  MAXROW = MAXROW + 1
                  NROWS  = NROWS  + 1
                  TCNROW(INODE) = NROWS
                  T_VAR(NROWS,INODE)%NAME   = DNAME(1)
                  T_VAR(NROWS,INODE)%LENGTH = LENGTH(1)
                  T_VAR(NROWS,INODE)%ROWTYP = 2
                  T_VAR(NROWS,INODE)%LOWER  = -PLINF
                  T_VAR(NROWS,INODE)%UPPER  =  PLINF
               ENDIF
!
!     Unrecognized code in code field. Ignore this row
!
            ELSE
               WRITE (IOLOG, 1400) NREC,QLINE
               WRITE (IOLOG, 1800) Q2,Q3
            ENDIF
            GOTO 100
!
!     COLUMN header has been read.
!
  200 CONTINUE
         IF (NECHO1 .GE.    1) WRITE (IOLOG, 1400) NREC,QLINE
         IF (Q1 .EQ. QC .AND. Q2 .EQ. QO) THEN
            NEXT = 1
         ELSE
            NEXT = 2
         ENDIF
!
         IF (IENDAT .EQ.    0) GOTO 9700
            GOTO 999
!
!     Come here if anything went wrong
!
 9000 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         IF (IENDAT .EQ. 0) GOTO 9020
         IF (IENDAT .EQ. 1) GOTO 9010
            WRITE (IOLOG, 3000)
            GOTO 999
!
 9010    CONTINUE
            WRITE (IOLOG, 3010)
            GOTO 9999
!
 9020    CONTINUE
            CORSTAT = 1
            WRITE (IOLOG, 3020)
            goto 999     !!RETURN
!
 9060 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3060)
         GOTO 9999
!
 9160 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3160)
         GOTO 9999
!
 9300 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3300)
         GOTO 9999
!
 9700 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3700)
         GOTO 9999
!
 9999 CONTINUE
         CORSTAT = 2
!
!     Return to calling program
!
  999 CONTINUE
         RETURN
!
 1400 FORMAT(I8,2X,A80)
 1800 FORMAT(' XXX - WARNING - Unrecognized code =',2A1,'= in ROWS',      &
     &       ' section. Row ignored.')
 2300 FORMAT(' XXX -  FATAL  - Illegal header card in CORE file.')
 2500 FORMAT(' XXX - WARNING - Time file had',I3,' periods; ROWS',      &
     &       ' section only finds',I3,'.')
 3000 FORMAT(' XXX - WARNING - Missing ENDATA card in CORE file.')
 3010 FORMAT(' XXX -  FATAL  - Detected EOF in CORE file.')
 3020 FORMAT(' XXX - WARNING - No information in CORE file.')
 3060 FORMAT(' XXX -  FATAL  - Unable to allocate memory in INCORE.')
 3160 FORMAT(' XXX -  FATAL  - Name does not match info in TIME file.')
 3300 FORMAT(' XXX -  FATAL  - No COLUMNS section specified.')
 3700 FORMAT(' XXX -  FATAL  - ROWS section is non-existent.')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INCOLS2 ( NREC, QLINE, ICTEMP, NEXT, IENDAT, IERR )
!
!     This routine reads the COLUMNS section in temporal form.
!
!     Routines called: CHOPEN, GETNAM, GETMPS, GETNET.
!
!     --------------------------------------------------------------
!     This version dated 21 September 2002. Written by Gus Gassmann.
!     --------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      REAL*8 ATEMP(4)
      CHARACTER*(NAMLEN) DNAME(3), DBLANK, DROW, DCOL, DVNAM, QLINE
      CHARACTER*8 DUNCAP, DUNDIR, DINT1, DINT2, DINT3, DINT4
      CHARACTER*1 DB255(NAMLEN), QROW(NAMLEN), QCOL(NAMLEN), T
      INTEGER*4   LENGTH(3),S3ROW,S3BLK,TLEN
      LOGICAL     MPSFORM
!
      INTEGER,    ALLOCATABLE, DIMENSION (:) :: TCNCOL
!
      EQUIVALENCE (DBLANK,DB255), (DROW,QROW), (DCOL,QCOL)
!
      DATA DUNCAP  /'UNCAP   '/, DUNDIR   /'UNDIR   '/,      &
     &     DINT1 /'''INTORG'''/, DINT2  /'''SOSORG'''/,      &
     &     DINT3 /'''INTEND'''/, DINT4  /'''SOSEND'''/,      &
     &     DB255 /NAMLEN*' '/
! ----------------------------------------------------------------------
!
!     Allocate temporary space for constraint matrix blocks.
!
         ALLOCATE( AUX2(IZANZB,NPER,NPER), IAUX2(IZANZB,NPER,NPER),      &
     &            LAUX2(IZCOLP,NPER,NPER), LMNS2(NPER,NPER), STAT=IER)
         IF (IER .GT. 0) GOTO 9060
         ALLOCATE ( TCOST(IZCOLP,NPER), STAT=IER)
         IF (IER .GT. 0) GOTO 9060
         ALLOCATE ( TCNCOL(NPER), STAT=IER)
         IF (IER .GT. 0) GOTO 9060
!
!     Initialize the columns section
!
            ICT    = 1
            IENDAT = 2
            NOFSOS = 0
!
            DO I=1,NPER
               TCNCOL(I) = TCNROW(I)
               DO J=1,NPER
                  LMNS2(I,J) = 0
                  LAUX2(1,I,J) = 1
               END DO
            END DO
!
            IF (NEXT .EQ. 1) GOTO 235
               GOTO 240
!
!     Come back here to process the next record
!     =========================================
!
  230    CONTINUE
            IF (MPSFORM) THEN
               CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOCOR,NREC,LENGTH)
            ELSE
               CALL GETNET(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOCOR,NREC,LENGTH)
            ENDIF
            IF (Q1 .EQ. QAST ) GOTO 230
            IF (IERR   .GT. 0) GOTO 9980
            IF (Q1 .EQ. QBL  ) GOTO 250
            IF (Q1 .EQ. QC .AND. Q2 .EQ. QO) GOTO 235
            IF (Q1 .EQ. QA .AND. Q2 .EQ. QR) GOTO 240
            IF (Q1 .EQ. QR .AND. Q2 .EQ. QH) GOTO 600
            IF (Q1 .EQ. QD .AND. Q2 .EQ. QE) GOTO 610
            IF (Q1 .EQ. QS .AND. Q2 .EQ. QU) GOTO 620
            IF (Q1 .EQ. QB .AND. Q2 .EQ. QO) GOTO 630
            IF (Q1 .EQ. QR .AND. Q2 .EQ. QA) GOTO 640
            IF (Q1 .EQ. QQ .AND. Q2 .EQ. QS) GOTO 642
            IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 645
            WRITE (IOLOG, 1400) NREC,QLINE
               WRITE (IOLOG, 2300)
               GOTO 9999
!
!     Header starting a new COLUMNS segment
!
  235    CONTINUE
            MPSFORM = .TRUE.
            INTTYP = 0
            DEFUPB = DEFUB
            DEFLOB = DEFLB
            GOTO 230
!
!     Header starting a new ARCS segment
!
  240    CONTINUE
            MPSFORM = .FALSE.
            INTTYP  = 0
            IF     (DNAME(1) .EQ. DUNCAP .OR. DNAME(2) .EQ. DUNCAP) THEN
               DEFUPB = PLINF
               DEFLOB = 0.D0
            ELSEIF (DNAME(1) .EQ. DUNDIR .OR. DNAME(2) .EQ. DUNDIR) THEN
               DEFUPB =  PLINF
               DEFLOB = -PLINF
            ELSE
               DEFUPB = 0.D0
               DEFLOB = 0.D0
            ENDIF
            GOTO 230
!
!    Here we have a regular data card.
!    =================================
!
!    We must determine whether it continues
!    a column already started, whether it begins or ends a segment of
!    integer variables, and whether it even contains relevant data
!
  250    CONTINUE
            IF (NECHO1 .GE. 3) WRITE (IOLOG, 1400) NREC,QLINE
            IF (MPSFORM) THEN
!
!     COLUMNS section.
!     ----------------
!     First match the column name to the last record read. Then check for
!     special markers to indicate integer variables and special ordered sets.
!
               JNM = 2
               IF (DNAME(1) .EQ. DVNAM) THEN
                  IF (DABS(ATEMP(1)) .GT. ZTOLIN) GOTO 430
                  IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 230
                      JNM = 3
                      ATEMP(1) = ATEMP(2)
                      GO TO 430
               ELSEIF (DNAME(3) .EQ. DINT1) THEN
                  INTTYP = 1
                  GOTO 230
               ELSEIF (DNAME(3) .EQ. DINT2) THEN
                  IF (NOFSOS .GE. IZSTYP) THEN
                     IF (.NOT. EXTMEM_STYP(NOFSOS+1,IERR)) GOTO 9060
                  ENDIF
!
                  NOFSOS = NOFSOS + 1
                  INTTYP = -NOFSOS
                  IF     (Q2 .EQ. QS .AND. Q3 .EQ. '1') THEN
                     SOSTYP(NOFSOS) = 1
                  ELSEIF (Q2 .EQ. QS .AND. Q3 .EQ. '2') THEN
                     SOSTYP(NOFSOS) = 2
                  ELSE
!
!     SOS type 3. The first name field must contain a valid row name.
!     Each column participating in this set must be 0/1, the row must be
!     an equality row, with rhs = 1, and coefficients in this row are 1.
!
                     DCOL = DNAME(1)
                     DO JPER=1,NPER
                        DO J=1,TCNROW(JPER)
                           IF (DCOL .EQ. T_VAR(J,JPER)%NAME) THEN
                              IF (T_VAR(J,JPER)%ROWTYP .NE. 0) GOTO 258
                                 SOSTYP(NOFSOS) = -J
                                 S3BLK = JPER
                                 S3ROW = J
                                 DEFUPB = 1.D0
                                 DEFLOB = 0.D0
                                 GOTO 230
                           ENDIF
                        END DO
                     END DO
!
  258             CONTINUE
                     WRITE (IOLOG, 1400) NREC,QLINE
                     WRITE (IOLOG, 1900)
                     GOTO 9999
                  ENDIF
                  GOTO 230
!
               ELSEIF (DNAME(3) .EQ. DINT3 .OR.      &
     &                 DNAME(3) .EQ. DINT4) THEN
                  INTTYP = 0
                  DEFUPB = DEFUB
                  DEFLOB = DEFLB
                  GOTO 230
!
!     Column name not previously encountered. Determine the period.
!
               ELSE
                  DCOL = DNAME(1)
                  DO J=ICT,ICTEMP
                     IF (DCOL .EQ. TCOLS(J)%NAME) THEN
                        IPNEXT = TCOLS(J)%PERIOD
                        DVNAM = DCOL
                        TLEN  = LENGTH(1)
                        ICT   = J
                        GOTO 370
                     ENDIF
                  END DO
!
                  DO J=1,ICT-1
                     IF (DCOL .EQ. TCOLS(J)%NAME) THEN
                        IPNEXT = TCOLS(J)%PERIOD
                        DVNAM = DCOL
                        TLEN  = LENGTH(1)
                        ICT   = J
                        GOTO 370
                     ENDIF
                  END DO
!
                  WRITE (IOLOG, 1400) NREC,QLINE
                  WRITE (IOLOG, 3900)
                  CORSTAT = 2
                  GOTO 999
               ENDIF
            ELSE
!
!     ARCS section
!     ------------
!     Arcs are not named, so we infer the time period as follows:
!     Check if the third name field contains an explicit period name.
!     If not, use the earlier of the periods associated with the FROM-node
!     (IPFROM) and the TO-node (IPTO).
!
               PURELP = .FALSE.
               INTTYP = 0
               JNM    = 2
               IF (Q2 .NE. QBL) THEN
                  T = Q2
               ELSE
                  T = Q3
               ENDIF
               IF ((DNAME(1) .EQ. DCOL) .AND.      &
     &             (DNAME(2) .EQ. DROW) .AND. (T .EQ. QCODE)) THEN
                  GOTO 500
               ELSE
                  IF (DNAME(3) .NE. DBLANK) THEN
                     DO IPNEXT=1,NPER
                        IF (DNAME(3) .EQ. DTIME(IPNEXT)) GOTO 285
                     END DO
                     GOTO 9160
                  ENDIF
!
!     Determine FROM-node or use previous info
!
                  IF (DNAME(1) .NE. DCOL) THEN
                     DCOL = DNAME(1)
                     DO IPFROM=1,NPER
                        DO IFROM=1,TCNROW(IPFROM)
                           IF (DCOL .EQ. T_VAR(IFROM,IPFROM)%NAME)      &
     &                        GOTO 270
                        END DO
                     END DO
!
                     WRITE (IOLOG, 1400) NREC,QLINE
                     WRITE (IOLOG, 3800)
                     CORSTAT = 2
                     GOTO 999
                  ENDIF
!
!     Determine TO-node or use previous info
!
  270          CONTINUE
                  IF (DNAME(2) .NE. DROW) THEN
                     DROW = DNAME(2)
                     DO IPTO=1,NPER
                        DO ITO=1,TCNROW(IPTO)
                           IF (DROW .EQ. T_VAR(ITO,IPTO)%NAME) GOTO 280
                        END DO
                     END DO
!
                     WRITE (IOLOG, 1400) NREC,QLINE
                     WRITE (IOLOG, 3900)
                     CORSTAT = 2
                     GOTO 999
                  ENDIF
  280          CONTINUE
                  IPNEXT = MIN0(IPFROM,IPTO)
!
!     Check the code field for a marker to indicate parallel arcs
!
  285          CONTINUE
                  IF (Q2 .NE. QBL) THEN
                     QCODE = Q2
                  ELSE
                     QCODE = Q3
                  ENDIF
               ENDIF
!
!     Manufacture a column name:   "{From-node}->{To-node}[:{code}]"
!
               DVNAM = DCOL(1:LENGTH(1)) // '->' // DROW(1:LENGTH(2))
               TLEN = LENGTH(1) + LENGTH(2) + 2
!
               IF (QCODE .NE. QBL) THEN
                  DVNAM = DVNAM(1:TLEN) // ':' // QCODE
                  TLEN = TLEN + 2
               ENDIF
            ENDIF
!
!     Start a new variable.
!
  370    CONTINUE
            NCOLS = TCNCOL(IPNEXT) + 1
            NROWS = TCNROW(IPNEXT)
            JCOST  = NCOLS - NROWS
            TCOST(JCOST, IPNEXT) = 0.D0
            TCNCOL(IPNEXT)       = NCOLS
            T_VAR(NCOLS,IPNEXT)%NAME   = DVNAM
            T_VAR(NCOLS,IPNEXT)%LENGTH = TLEN
            T_VAR(NCOLS,IPNEXT)%LOWER  = DEFLOB
            T_VAR(NCOLS,IPNEXT)%UPPER  = DEFUPB
!
!     Check array sizes and expand if necessary
!
            IF (NCOLS .GT. IZCOLP) THEN
               IF (.NOT. EXTMEM_COLP(NCOLS,IERR)) GOTO 9060
            ENDIF
!
!     Set column pointers for all temporary matrix blocks
!
            LCC = NCOLS - NROWS
            DO JP=1,NPER
               LAUX2(LCC+1,JP,IPNEXT) = LMNS2(JP,IPNEXT) + 1
            END DO
!
!     Store the column type and, if SOS type 3, the implied matrix coefficient
!
            IF (INTTYP .GE. 0 .OR. SOSTYP(NOFSOS) .GT. 0) THEN
               T_VAR(NCOLS,IPNEXT)%ROWTYP = INTTYP
            ELSE
               T_VAR(NCOLS,IPNEXT)%ROWTYP = 1
               JELEM = LMNS2(S3BLK,IPNEXT)
               IF (JELEM .GE. IZANZB) THEN
                  IF (.NOT. EXTMEM_ANZB(JELEM+1,IERR)) GOTO 9060
               ENDIF
!
               LMNS2(      S3BLK,IPNEXT) = JELEM + 1
               LAUX2(LCC,  S3BLK,IPNEXT) = JELEM + 2
               IAUX2(JELEM,S3BLK,IPNEXT) = S3ROW
                AUX2(JELEM,S3BLK,IPNEXT) = 1.D0
            ENDIF
!
!        Test for row match
!
  430    CONTINUE
         IF (.NOT. MPSFORM) GOTO 500
            DROW = DNAME(JNM)
!
!     Cost coefficients (even if fixed) are not stored in the A-matrix
!
            IF (DROW .EQ. OBJNAM) THEN
               TCOST(JCOST, IPNEXT) = ATEMP(1)
               GOTO 490
            ENDIF
!
            DO JP=IPNEXT,NPER
               DO JR=1,TCNROW(JP)
                  IF (DROW .EQ. T_VAR(JR,JP)%NAME) GOTO 470
               END DO
            END DO
!
!     Maybe it is a linking constraint or an alternative objective row?
!
            DO JP=1,IPNEXT-1
               DO JR=1,TCNROW(JP)
                  IF (DROW .EQ. T_VAR(JR,JP)%NAME) THEN
                     IF (.NOT. ALLOW_GC) THEN
                        WRITE (IOLOG, 1400) NREC,QLINE
                        WRITE (IOLOG, 2600) DROW
                        GOTO 470
                     ENDIF
                  ENDIF
               END DO
            END DO
!
!     No match!
!
  460       CONTINUE
               WRITE (IOLOG, 1400) NREC,QLINE
               WRITE (IOLOG, 1100) DROW
               GOTO 9999
!
!     Matched a coefficient in the A-matrix. Expand array if needed, then store.
!
  470 CONTINUE
               JELEM = LMNS2(JP,IPNEXT) + 1
               IF (JELEM .GT. IZANZB) THEN
                  IF (.NOT. EXTMEM_ANZB(JELEM,IERR)) GOTO 9060
               ENDIF
!
               LMNS2(              JP,IPNEXT) = JELEM
               LAUX2(NCOLS+1-NROWS,JP,IPNEXT) = JELEM + 1
               IAUX2(JELEM,        JP,IPNEXT) = JR + 1
                AUX2(JELEM,        JP,IPNEXT) = ATEMP(1)
!
!     There may be another coefficient on this card. If not, get the next one
!
  490    CONTINUE
            IF (JNM .EQ. 3) GOTO 230
            IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 230
               JNM = 3
               ATEMP(1) = ATEMP(2)
               GOTO 430
!
!        Data lines in ARCS segment define arcs or side constraints,
!        depending on the value of the third name field.
!
  500    CONTINUE
         IF (DNAME(3) .NE. DBLANK) GOTO 530
!
!     First store the cost coefficient and the bounds (if specified)
!
            TCOST(JCOST, IPNEXT) = ATEMP(1)
            IF (DABS(ATEMP(2)) .GT. ZTOLIN)      &
     &                              T_VAR(NCOLS,IPNEXT)%UPPER = ATEMP(2)
            IF (DABS(ATEMP(3)) .GT. ZTOLIN)      &
     &                              T_VAR(NCOLS,IPNEXT)%LOWER = ATEMP(3)
!
!     The coefficient in the FROM-node is -1.D0
!
               JELEM = LMNS2(JPER,IPNEXT) + 1
               IF (JELEM .GT. IZANZB) THEN
                  IF (.NOT. EXTMEM_ANZB(JELEM,IERR)) GOTO 9060
               ENDIF
!
               LMNS2(JPER, IPNEXT)               = JELEM
               LAUX2(NCOLS+1-NROWS,JPER, IPNEXT) = JELEM + 1
               IAUX2(JELEM,        JPER, IPNEXT) = IFROM + 1
                AUX2(JELEM,        JPER, IPNEXT) = -1.D0
!
!      The coefficient in the TO-node could be a multiplier
!
               JELEM = LMNS2(JPER,IPNEXT) + 1
               IF (JELEM .GT. IZANZB) THEN
                  IF (.NOT. EXTMEM_ANZB(JELEM,IERR)) GOTO 9060
               ENDIF
!
               LMNS2(              JPER, IPNEXT) = JELEM
               LAUX2(NCOLS+1-NROWS,JPER, IPNEXT) = JELEM + 1
               IAUX2(JELEM,        JPER, IPNEXT) = ITO + 1
               IF (DABS(ATEMP(4))  .GT.  ZTOLIN)  THEN
                  AUX2(JELEM,      JPER, IPNEXT) = ATEMP(4)
               ELSE
                  AUX2(JELEM,      JPER, IPNEXT) = 1.0
              ENDIF
              GOTO 230
!
!        Here we have a side constraint. Identify the row.
!
  530    CONTINUE
         DROW = DNAME(3)
!
!     Cost coefficient should be specified on a pure arc card
!
            IF (DROW .EQ. OBJNAM) THEN
               WRITE (IOLOG, 1400) NREC,QLINE
               WRITE (IOLOG, 1600)
!
            ELSE
            DO JP=IPNEXT,NPER
               DO JR=2,TCNROW(JP)
                  IF (DROW .EQ. T_VAR(JR,JP)%NAME) GOTO 570
               END DO
            END DO
!
!     Maybe it is a linking constraint or an alternative objective row?
!
            DO JP=1,IPNEXT-1
               DO JR=2,TCNROW(JP)
                  IF (DROW .EQ. T_VAR(JR,JP)%NAME) THEN
                     IF (.NOT. ALLOW_GC) THEN
                        WRITE (IOLOG, 1400) NREC,QLINE
                        WRITE (IOLOG, 2600) DROW
                        GOTO 570
                     ENDIF
                  ENDIF
               END DO
            END DO
!
!     No match!
!
  560       CONTINUE
               WRITE (IOLOG, 1400) NREC,QLINE
               WRITE (IOLOG, 1100) DROW
               GOTO 9999
!
!     Matched a coefficient in the A-matrix. Expand array if needed, then store.
!
  570       CONTINUE
               JELEM = LMNS2(JP,IPNEXT) + 1
               IF (JELEM .GT. IZANZB) THEN
                  IF (.NOT. EXTMEM_ANZB(JELEM,IERR)) GOTO 9060
               ENDIF
!
               LMNS2(              JP,IPNEXT) = JELEM
               LAUX2(NCOLS+1-NROWS,JP,IPNEXT) = JELEM + 1
               IAUX2(JELEM,        JP,IPNEXT) = JR
                AUX2(JELEM,        JP,IPNEXT) = ATEMP(1)
            ENDIF
            GOTO 230
!
!     The COLUMNS section is done. What is next?
!
!     RHS section
!
  600    CONTINUE
            NEXT = 1
            GOTO 650
!
!     DEMAND section
!
  610    CONTINUE
            NEXT = 2
            GOTO 650
!
!     SUPPLY section
!
  620    CONTINUE
            NEXT = 3
            GOTO 650
!
!     BOUNDS section
!
  630    CONTINUE
            NEXT = 4
            GOTO 650
!
!     RANGES section
!
  640    CONTINUE
            NEXT = 5
            GOTO 650
!
!     In-line QSECT section
!
  642    CONTINUE
            NEXT = 6
            GOTO 650
!
!     ENDATA record
!
  645    CONTINUE
            NEXT = 7
!
!     Now we place all the information read into the data structure
!     and initial the event tree. First we determine the presence of
!     linking constraints and/or full block lower diagonal structure.
!
  650 CONTINUE
         MARKOV  = .TRUE.
         HAVE_GC = .FALSE.
         DO IP=1,NPER
            DO JP=IP+2,NPER
               IF (LMNS2(JP,IP) .GT. 0) THEN
                  MARKOV = .FALSE.
                  GOTO 652
               ENDIF
            END DO
         END DO
!
  652 CONTINUE
         DO IP=1,NPER
            DO JP=1,IP-1
               IF (LMNS2(JP,IP) .GT. 0) THEN
                  HAVE_GC = .TRUE.
                  GOTO 654
               ENDIF
            END DO
         END DO
!
!     We construct the event tree (path) one node at a time, allowing
!     for the possibility of empty periods, which are suppressed.
!
  654 CONTINUE
         INODE  = 0
         LASTA  = 0
         LASTBA = 0
         LASTBD = 0
         LASTC  = 0
         LASTCA = 0
         LASTCH = 0
         LASTD  = 0
         LASTNM = 0
         LASTR  = 0
         LASTVP = 0
         MAXROW = 0
         MAXCOL = 0
         XINFTY = 1.D+31
!
         DO IP=1,NPER
            IF (TCNCOL(IP) .GT. 0) THEN
               INODE = INODE + 1
               IRNGE0(INODE) = INODE
               IF (INODE .GT. IZNODE) THEN
                  IF (.NOT. EXTMEM_NODE(INODE,IERR)) GOTO 9060
               ENDIF
               NODEINFO(INODE)%PROB   = 1.0
               NODEINFO(INODE)%IANCTR = INODE - 1
               NODEINFO(INODE)%IDESC  = 0
               NODEINFO(INODE)%NDESC  = 0
               NODEINFO(INODE)%IBROTH = 0
               IF (INODE .GT. 1) THEN
                  NODEINFO(INODE-1)%IDESC = INODE
                  NODEINFO(INODE-1)%NDESC = 1
               END IF
!
               NODEINFO(INODE)%KCOL   = MAXCOL
               NODEINFO(INODE)%KROW   = MAXROW
               NODEINFO(INODE)%KCOST  = LASTC
               NODEINFO(INODE)%KRHS   = LASTR
               NODEINFO(INODE)%KBOUND = LASTBD
               NODEINFO(INODE)%KNAMES = LASTNM
               NODEINFO(INODE)%VARPTR = LASTVP
               NODEINFO(INODE)%NCOL   = TCNCOL(IP) + 1
               NODEINFO(INODE)%NROW   = TCNROW(IP) + 1
               KDATA(INODE) = LASTBA
!
!     Extend memory if necessary
!
               MNBNDS = LASTBD + TCNCOL(IP) + 2
               IF (MNBNDS .GT. IZBNDS) THEN
                  IF (.NOT. EXTMEM_BNDS(MNBNDS,IERR)) GOTO 9060
               ENDIF
!
               MNCOST = LASTC  + TCNCOL(IP) - TCNROW(IP)
               IF (MNCOST .GT. IZCOST) THEN
                  IF (.NOT. EXTMEM_COST(MNCOST,IERR)) GOTO 9060
               ENDIF
!
               MNVNAM = LASTNM + TCNCOL(IP) + 1
               IF (MNVNAM .GT. IZVNAM) THEN
                  IF (.NOT. EXTMEM_VNAM(MNVNAM,IERR)) GOTO 9060
               ENDIF
!
               MNVTYP = LASTVP + TCNCOL(IP) + 2
               IF (MNVTYP .GT. IZVTYP) THEN
                  IF (.NOT. EXTMEM_VTYP(MNVTYP,IERR)) GOTO 9060
               ENDIF
!
!     Copy objective name and default bounds
!
               MNCHAR = LASTCH + LENOBJ
               IF (MNCHAR .GT. IZCHAR) THEN
                  IF (.NOT. EXTMEM_CHAR(MNCHAR,IERR)) GOTO 9060
               ENDIF
               DROW = OBJNAM
               DO I=1,LENOBJ
                  QNAMES(LASTCH+I) = QROW(I)
               END DO
               NAME1(LASTNM+1)  = LASTCH + 1
               XLB(LASTBD+1)    = -PLINF
               XUB(LASTBD+1)    =  PLINF
               VARTYP(LASTVP+1) =  2
               LASTCH = LASTCH + LENOBJ
!
!     Copy variable names, types, default bounds
!
               DO IR=1,TCNCOL(IP)
                  TLEN = T_VAR(IR,IP)%LENGTH
                  MNCHAR = LASTCH + TLEN
                  IF (MNCHAR .GT. IZCHAR) THEN
                     IF (.NOT. EXTMEM_CHAR(MNCHAR,IERR)) GOTO 9060
                  ENDIF
                  DROW = T_VAR(IR,IP)%NAME
                  DO I=1,TLEN
                     QNAMES(LASTCH+I) = QROW(I)
                  END DO
                  NAME1(LASTNM+IR+1) = LASTCH + 1
                  XLB(LASTBD+IR+1)    = T_VAR(IR,IP)%LOWER
                  XUB(LASTBD+IR+1)    = T_VAR(IR,IP)%UPPER
                  VARTYP(LASTVP+IR+1) = T_VAR(IR,IP)%ROWTYP
                  LASTCH = LASTCH + TLEN
               END DO
               NAME1(LASTNM+TCNCOL(IP)+2) = LASTCH + 1
!
!     Initial the theta column
!
               ITHETA = LASTBD + TCNCOL(IP) + 2
               XLB(ITHETA) = -XINFTY
               XUB(ITHETA) =  XINFTY
!
!     Copy cost coefficients
!
               DO IR=1,TCNCOL(IP)-TCNROW(IP)
                  COST(LASTC+IR) = TCOST(IR,IP)
               END DO
!
!     Copy the A_matrix block on the main diagonal
!
               MNALMN = LASTA  + LMNS2(IP,IP)
               IF (MNALMN .GT. IZALMN) THEN
                  IF (.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 9060
               ENDIF
               NCOLA = TCNCOL(IP) - TCNROW(IP)
               MNACOL = LASTCA + NCOLA + 1
               IF (MNACOL .GT. IZACOL) THEN
                  IF (.NOT. EXTMEM_ACOL(MNACOL,IERR)) GOTO 9060
               ENDIF
!
               JBLK = KDATA(INODE) + 1
               IF (JBLK .GT. IZABLK) THEN
                  IF (.NOT. EXTMEM_ABLK(JBLK,IERR)) GOTO 9060
               ENDIF
               KCOLA(JBLK) = LASTCA
               KELMA(JBLK) = LASTA
               NELMA(JBLK) = LMNS2(IP,IP)
               DO JC=1,NCOLA
                  LL=LAUX2(JC,IP,IP)
                  KK=LAUX2(JC+1,IP,IP) - 1
                  LA(LASTCA+JC) = LL
                  DO JR=LL,KK
                     IA(LASTA+JR) = IAUX2(JR,IP,IP)
                      A(LASTA+JR) =  AUX2(JR,IP,IP)
                  END DO
               END DO
               LASTA  = LASTA  + LMNS2(IP,IP)
               LASTCA = LASTCA + NCOLA + 1
               LA(LASTCA) = KK + 1
!
!     Copy the A-matrix block(s) to the left of (below) the main diagonal
!
               IF (IP .GT. 1) THEN
                  DO JP=IP-1,1,-1
                     IF (TCNROW(JP) .GT. 0 .AND. TCNCOL(JP) .GT. 0) THEN
                        MNALMN = LASTA  + LMNS2(IP,JP)
                        IF (MNALMN .GT. IZALMN) THEN
                           IF (.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 9060
                        ENDIF
                        NCOLB = TCNCOL(JP) - TCNROW(JP)
                        MNACOL = LASTCA + NCOLB + 1
                        IF (MNACOL .GT. IZACOL) THEN
                           IF (.NOT. EXTMEM_ACOL(MNACOL,IERR)) GOTO 9060
                        ENDIF
!
                        JBLK = JBLK + 1
                        IF (JBLK .GT. IZABLK) THEN
                           IF (.NOT. EXTMEM_ABLK(JBLK,IERR)) GOTO 9060
                        ENDIF
                        KCOLA(JBLK) = LASTCA
                        KELMA(JBLK) = LASTA
                        NELMA(JBLK) = LMNS2(IP,JP)
                        DO JC=1,NCOLA
                           LL=LAUX2(JC,IP,JP)
                           KK=LAUX2(JC+1,IP,JP) - 1
                           LA(LASTCA+JC) = LL
                           DO JR=LL,KK
                              IA(LASTA+JR) = IAUX2(JR,IP,JP)
                               A(LASTA+JR) =  AUX2(JR,IP,JP)
                           END DO
                        END DO
                        LASTA  = LASTA  + LMNS2(IP,JP)
                        LASTCA = LASTCA + NCOLB + 1
                        LA(LASTCA) = KK + 1
                        IF (MARKOV) GOTO 670
                     ENDIF
                  END DO
!
!     Copy the blocks above the main diagonal (linking constraints) --- if any
!
  670          CONTINUE
                  IF (HAVE_GC) THEN
                     DO JP=1,IP-1
                        IF (TCNROW(JP) .GT. 0 .AND.      &
     &                      TCNCOL(JP) .GT. 0) THEN
                           MNALMN = LASTA  + LMNS2(JP,IP)
                           IF (MNALMN .GT. IZALMN) THEN
                              IF (.NOT. EXTMEM_ALMN(MNALMN,IERR))      &
     &                           GOTO 9060
                           ENDIF
                           MNACOL = LASTCA + NCOLA + 1
                           IF (MNACOL .GT. IZACOL) THEN
                              IF (.NOT. EXTMEM_ACOL(MNACOL,IERR))      &
     &                           GOTO 9060
                           ENDIF
!
                           JBLK = JBLK + 1
                           IF (JBLK .GT. IZABLK) THEN
                              IF (.NOT. EXTMEM_ABLK(JBLK,IERR))      &
     &                           GOTO 9060
                           ENDIF
                           KCOLA(JBLK) = LASTCA
                           KELMA(JBLK) = LASTA
                           NELMA(JBLK) = LMNS2(JP,IP)
                           DO JC=1,NCOLA
                              LL=LAUX2(JC,JP,IP)
                              KK=LAUX2(JC+1,JP,IP) - 1
                              LA(LASTCA+JC) = LL
                              DO JR=LL,KK
                                 IA(LASTA+JR) = IAUX2(JR,JP,IP)
                                  A(LASTA+JR) =  AUX2(JR,JP,IP)
                              END DO
                           END DO
                           LASTA  = LASTA  + LMNS2(JP,IP)
                           LASTCA = LASTCA + NCOLA + 1
                           LA(LASTCA) = KK + 1
                        ENDIF
                     END DO
                  ENDIF
               ENDIF
!
               LASTBA = JBLK
               LASTBD = LASTBD + TCNCOL(IP) + 2
               LASTC  = LASTC  + TCNCOL(IP) - TCNROW(IP)
               LASTNM = LASTNM + TCNCOL(IP) + 1
               LASTR  = LASTR  + TCNROW(IP) + 1
               LASTVP = LASTVP + TCNCOL(IP) + 1
               MAXROW = MAXROW + TCNROW(IP) + 1
               MAXCOL = MAXCOL + TCNCOL(IP) + 2
            ENDIF
         END DO
         NPER = INODE
!
         IENDAT = 3
         GOTO 999
!
!     Come here if anything went wrong
!
 9060 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3060)
         GOTO 9999
!
 9160 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3160)
         GOTO 9999
!
 9980 CONTINUE
         WRITE (IOLOG, 3980)
         GOTO 9999
!
 9999 CONTINUE
         CORSTAT = 2
!
!     Return to calling program
!
  999 CONTINUE
         DEALLOCATE(AUX2,IAUX2,LAUX2,LMNS2,STAT=IER)
         IF (IER .NE. 0) THEN
            WRITE (IOLOG, 3070)
            CORSTAT = 2
         ENDIF
         RETURN
!
 1100 FORMAT(' XXX -  FATAL  - Row "',A20,'" violates time order in',      &
     &       ' core file.')
 1400 FORMAT(I8,2X,A80)
 1600 FORMAT(' XXX - WARNING - Cost coefficient should be specified',      &
     &       ' on a pure arc card. Data ignored.')
 1900 FORMAT(' XXX -  FATAL  - Inappropriate row name in SOS type 3.')
 2300 FORMAT(' XXX -  FATAL  - Illegal header card in CORE file.')
 2600 FORMAT(' XXX - WARNING - Row "',A20,'" violates time order',      &
     &       ' Enable linking constraints.')
 3060 FORMAT(' XXX -  FATAL  - Unable to allocate memory in INCORE')
 3070 FORMAT(' XXX -  FATAL  - Memory handler returned error in INCORE')
 3160 FORMAT(' XXX -  FATAL  - Name does not match info in TIME file.')
 3750 FORMAT(' XXX -  FATAL  - Error in COLUMNS section.')
 3800 FORMAT(' XXX -  FATAL  - FROM-node not previously specified.')
 3900 FORMAT(' XXX -  FATAL  - TO-node not previously specified.')
 3950 FORMAT(' XXX -  FATAL  - ARCS section violates temporal order.')
 3980 FORMAT(' XXX -  FATAL  - Error while reading CORE file.')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INRHS(NREC, QLINE, NEXT, IENDAT, IERR)
!
!     This routine reads the remaining sections of the core file.
!
!     Routines called: CHOPEN, GETNAM, GETMPS, GETNET.
!
!     --------------------------------------------------------------
!     This version dated 21 September 2002. Written by Gus Gassmann.
!     --------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      REAL*8 ATEMP(4)
      CHARACTER*(NAMLEN) DNAME(3), QLINE, DBLANK, DROW
      CHARACTER*1 DB255(NAMLEN)
      INTEGER*4   LENGTH(3)
      LOGICAL     MPSFORM,SUPPLY
!
      EQUIVALENCE (DBLANK,DB255)
!
      DATA DB255 /NAMLEN*' '/
! ----------------------------------------------------------------------
!
!     The COLUMNS section is done. What is next?
!
!
!     Allocate and initial the RHS and auxiliary array STOCHA
!
  660 CONTINUE
               IF (MAXROW .GT. IZDRHS) THEN
                  IF (.NOT. RESIZE .OR. (IZDRHS .GE. MXDRHS)      &
     &                             .OR. (MAXROW .GT. MXDRHS)) GOTO 9060
                  IZDRHS = MAXROW
               ENDIF
               ALLOCATE(RHS(IZDRHS), STAT=IERR)
               IF (IERR .NE. 0) GOTO 9060
               DO I=1,NPER
                  JS = NODEINFO(I)%KRHS
                  NR = NODEINFO(I)%NROW
                  DO J=1,NR
                     RHS(JS+J) = 0.D0
                  END DO
               END DO
               IP0 = 1
               I0 = 1
               ALLOCATE(STOCHA(NPER,NPER), STAT=IERR)
               IF (IERR .NE. 0) GOTO 9060
               DO I=1,NPER
                  DO J=1,NPER
                     STOCHA(I,J) = .FALSE.
                  END DO
               END DO
            GOTO (600, 610, 620, 630, 640), NEXT
!
!     RHS section
!
  600    CONTINUE
            L    = 1
            SUPPLY  = .FALSE.
            MPSFORM = .TRUE.
            NAMFLD = 2
            GOTO 700
!
!     DEMAND section
!
  610    CONTINUE
            L = 1
            SUPPLY  = .FALSE.
            MPSFORM = .FALSE.
            NAMFLD = 1
            GOTO 700
!
!     SUPPLY section
!
  620    CONTINUE
            L = 1
            SUPPLY  = .TRUE.
            MPSFORM = .FALSE.
            NAMFLD = 1
            GOTO 700
!
!     BOUNDS section
!
  630    CONTINUE
            L = 2
            MPSFORM = .TRUE.
            GOTO 700
!
!     RANGES section
!
  640    CONTINUE
            L = 3
            MPSFORM = .TRUE.
!
!     Process right-hand sides (including supply and demand), bounds and ranges
!
  700 CONTINUE
         IF (MPSFORM) THEN
            CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                        IERR,IOCOR,NREC,LENGTH)
         ELSE
            CALL GETNET(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                        IERR,IOCOR,NREC,LENGTH)
         ENDIF
         IF (Q1 .EQ. QAST ) GOTO 700
         IF (IERR   .GT. 0) GOTO 710
         IF (Q1 .EQ. QBL  ) GOTO 720
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QH) GOTO 600
         IF (Q1 .EQ. QD .AND. Q2 .EQ. QE) GOTO 610
         IF (Q1 .EQ. QS .AND. Q2 .EQ. QU) GOTO 620
         IF (Q1 .EQ. QB .AND. Q2 .EQ. QO) GOTO 630
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QA) GOTO 640
         IF (Q1 .EQ. QQ .AND. Q2 .EQ. QS) GOTO 980
         IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 999
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 2300)
            GOTO 9999
!
  710    CONTINUE
            IF (IERR  .EQ. 2) GOTO 9000
               GOTO 9980
!
  720    CONTINUE
            IF (NECHO1 .GE. 3) WRITE (IOLOG, 1400) NREC,QLINE
            GOTO (750,800,850,999), L
!
!     RHS section (includes SUPPLY and DEMAND)
!
  750    CONTINUE
         J = NAMFLD
         IF (DABS(ATEMP(1)) .GT. ZTOLIN) GOTO 760
         IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 700
            J = 3
            ATEMP(1) = ATEMP(2)
!
!        Test for row match
!
  760    CONTINUE
         IF (DRHS .EQ. DOTS) DRHS = DNAME(1)
         IF (MPSFORM .AND. (DRHS .NE. DNAME(1))) GOTO 700
         DROW = DNAME(J)
         IP = IP0
         KN = NODEINFO(IP)%KNAMES
         NR = NODEINFO(IP)%NROW
         IRES = NMSRCH(DROW,KN+I0,KN+NR)
         IF (IRES .GT. 0) GOTO 790
         IRES = NMSRCH(DROW,KN+1 ,KN+I0)
         IF (IRES .GT. 0) GOTO 790
         DO IP=IP0+1,NPER
            KN = NODEINFO(IP)%KNAMES
            NR = NODEINFO(IP)%NROW
            IRES = NMSRCH(DROW,KN+1,KN+NR)
            IF (IRES .GT. 0) GOTO 790
         END DO
         DO IP=1,IP0-1
            KN = NODEINFO(IP)%KNAMES
            NR = NODEINFO(IP)%NROW
            IRES = NMSRCH(DROW,KN+1,KN+NR)
            IF (IRES .GT. 0) GOTO 790
         END DO
!
!     No match!
!
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 1200) DROW
         GOTO 9999
!
!     Matched
!
  790    CONTINUE
         IP0 = IP
         I0  = IRES - KN
         IF (SUPPLY) THEN
            RHS(NODEINFO(IP)%KRHS+I0) = -ATEMP(1)
         ELSE
            RHS(NODEINFO(IP)%KRHS+I0) =  ATEMP(1)
         ENDIF
         IF (J .EQ. 1 .OR. J .EQ. 3) GOTO 700
            IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 700
               J = 3
               ATEMP(1) = ATEMP(2)
               GOTO 760
!
!     BOUNDS Section. Match the column name.
!
  800 CONTINUE
         IF (DBOUND .EQ. DOTS     ) DBOUND = DNAME(1)
         IF (DBOUND .NE. DNAME(1) ) GOTO 700
         DROW = DNAME(2)
         IP = IP0
         KN = NODEINFO(IP)%KNAMES
         NC = NODEINFO(IP)%NCOL
         IRES = NMSRCH(DROW,KN+I0,KN+NC)
         IF (IRES .GT. 0) GOTO 830
         IRES = NMSRCH(DROW,KN+1 ,KN+I0)
         IF (IRES .GT. 0) GOTO 830
         DO IP=IP0+1,NPER
            KN = NODEINFO(IP)%KNAMES
            NC = NODEINFO(IP)%NCOL
            IRES = NMSRCH(DROW,KN+1,KN+NC)
            IF (IRES .GT. 0) GOTO 830
         END DO
         DO IP=1,IP0-1
            KN = NODEINFO(IP)%KNAMES
            NC = NODEINFO(IP)%NCOL
            IRES = NMSRCH(DROW,KN+1,KN+NC)
            IF (IRES .GT. 0) GOTO 830
         END DO
!
!     No match!
!
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 1300) DROW
            GOTO 9999
!
!     Matched. Now determine the bound type
!
  830    CONTINUE
         IP0 = IP
         I0  = IRES - KN
         IC  = NODEINFO(IP)%KBOUND + I0
!
         IF     (Q2 .EQ. QL .AND. Q3 .EQ. QO) THEN
            XLB(IC) = ATEMP(1)
         ELSEIF (Q2 .EQ. QU .AND. Q3 .EQ. QP) THEN
            XUB(IC) = ATEMP(1)
         ELSEIF (Q2 .EQ. QF .AND. Q3 .EQ. QX) THEN
            XLB(IC) = ATEMP(1)
            XUB(IC) = ATEMP(1)
         ELSEIF (Q2 .EQ. QF .AND. Q3 .EQ. QR) THEN
            XLB(IC) = -PLINF
            XUB(IC) =  PLINF
         ELSEIF (Q2 .EQ. QM .AND. Q3 .EQ. QI) THEN
            XLB(IC) = -PLINF
         ELSEIF (Q2 .EQ. QP .AND. Q3 .EQ. QL) THEN
            XUB(IC) =  PLINF
         ELSEIF (Q2 .EQ. QU .AND. Q3 .EQ. QI) THEN
            XUB(IC) =  DINT(ATEMP(1))
         ELSEIF (Q2 .EQ. QB .AND. Q3 .EQ. QV) THEN
            XLB(IC) =  0.D0
            XUB(IC) =  1.D0
         ELSE
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 1700)
            GOTO 9999
         ENDIF
         GOTO 700
!
!     RANGES Section. Match the row name.
!
  850    CONTINUE
         IF (DRANGE .EQ. DOTS     ) DRANGE = DNAME(1)
         IF (DRANGE .NE. DNAME(1) ) GOTO 700
         J = 2
         IF (DABS(ATEMP(1)) .GT. ZTOLIN) GOTO 860
         IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 700
            J = 3
            ATEMP(1) = ATEMP(2)
!
!        Test for row match
!
  860    CONTINUE
         DROW = DNAME(J)
         IP = IP0
         KN = NODEINFO(IP)%KNAMES
         NR = NODEINFO(IP)%NROW
         IRES = NMSRCH(DROW,KN+I0,KN+NR)
         IF (IRES .GT. 0) GOTO 890
         IRES = NMSRCH(DROW,KN+1 ,KN+I0)
         IF (IRES .GT. 0) GOTO 890
         DO IP=IP0+1,NPER
            KN = NODEINFO(IP)%KNAMES
            NR = NODEINFO(IP)%NROW
            IRES = NMSRCH(DROW,KN+1,KN+NR)
            IF (IRES .GT. 0) GOTO 890
         END DO
         DO IP=1,IP0-1
            KN = NODEINFO(IP)%KNAMES
            NR = NODEINFO(IP)%NROW
            IRES = NMSRCH(DROW,KN+1,KN+NR)
            IF (IRES .GT. 0) GOTO 890
         END DO
!
!     No match!
!
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 1200) DROW
         GOTO 9999
!
!     Matched
!
  890 CONTINUE
         I   = IRES - KN
         KB  = NODEINFO(IP)%KBOUND
         IV  = NODEINFO(IP)%VARPTR + I
         IT  = VARTYP(IV)
         I0  = I
         IP0 = IP
         IF (IT .EQ.  1) THEN
            XUB(KB+I) =  DABS(ATEMP(1))
         ELSEIF (IT .EQ. -1) THEN
            XLB(KB+I) = -DABS(ATEMP(1))
         ELSEIF (IT .EQ.  0) THEN
            IF (ATEMP(1) .GT. 0.) XUB(KB+I) = ATEMP(1)
            IF (ATEMP(1) .LT. 0.) XLB(KB+I) =-ATEMP(1)
         ELSE
            GOTO 9600
         ENDIF
         IF (J .EQ. 3) GOTO 700
         IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 700
            J = 3
            ATEMP(1) = ATEMP(2)
            GOTO 860
!
!     Inline QSECT section
!
  980 CONTINUE
         NEXT = 6
         GOTO 999
!
!     Come here if anything went wrong
!
 9000 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         IF (IENDAT .EQ. 0) GOTO 9020
         IF (IENDAT .EQ. 1) GOTO 9010
            WRITE (IOLOG, 3000)
            GOTO 999
!
 9010    CONTINUE
            WRITE (IOLOG, 3010)
            GOTO 9999
!
 9020    CONTINUE
            CORSTAT = 1
            WRITE (IOLOG, 3020)
            goto 999     !!RETURN
!
 9060 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3060)
         GOTO 9999
!
 9600 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3600)
         GOTO 9999
!
 9980 CONTINUE
         WRITE (IOLOG, 3980)
         GOTO 9999
!
 9999 CONTINUE
         CORSTAT = 2
!
!     Return to calling program
!
  999 CONTINUE
         RETURN
!
 1200 FORMAT(' XXX -  FATAL  - Unmatched row name "',A20,      &
     &       '" in RHS or RANGES section.')
 1300 FORMAT(' XXX -  FATAL  - Unmatched column name "',A20,      &
     &       '" in BOUNDS section.')
 1400 FORMAT(I8,2X,A80)
 1700 FORMAT(' XXX -  FATAL  - Unrecognized or missing code in BOUNDS',      &
     &       ' section.')
 2300 FORMAT(' XXX -  FATAL  - Illegal header card in CORE file.')
 3000 FORMAT(' XXX - WARNING - Missing ENDATA card in CORE file.')
 3010 FORMAT(' XXX -  FATAL  - Detected EOF in CORE file.')
 3020 FORMAT(' XXX - WARNING - No information in CORE file.')
 3060 FORMAT(' XXX -  FATAL  - Unable to allocate memory in INCORE')
 3600 FORMAT(' XXX -  FATAL  - Illegal row type in RANGES section.')
 3980 FORMAT(' XXX -  FATAL  - Error while reading CORE file.')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CPROWS(ND, IERR)
!
!     This routine copies the row information from temporary storage
!     to permanent array locations.
!
!     Routines called: CHOPEN, GETNAM, GETMPS, GETNET.
!
!     --------------------------------------------------------------
!     This version dated 21 September 2002. Written by Gus Gassmann.
!     --------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DBLANK, DROW
      CHARACTER*1 DB255(NAMLEN), QROW(NAMLEN)
!
      EQUIVALENCE (DBLANK,DB255), (DROW, QROW)
!
      DATA DB255 /NAMLEN*' '/
! ----------------------------------------------------------------------
!
            NCHARS = 0
            DO J=1,NROWS
               NCHARS = NCHARS + T_ROWS(IROW+J)%LENGTH
            END DO
!
            MNCHAR = LASTCH + NCHARS
            IF (MNCHAR .GT. IZCHAR) THEN
               IF (.NOT. EXTMEM_CHAR(MNCHAR,IERR)) GOTO 9060
            ENDIF
!
            IF (MAXROW .GE. IZBNDS) THEN
               IF (.NOT. EXTMEM_BNDS(MAXROW+1,IERR)) GOTO 9060
            ENDIF
!
            IF (MAXROW .GT. IZVTYP) THEN
               IF (.NOT. EXTMEM_VTYP(MAXROW,IERR)) GOTO 9060
            ENDIF
!
            IROW   = NODEINFO(ND)%KROW
            IVTYP  = NODEINFO(ND)%VARPTR
            IBOUND = NODEINFO(ND)%KBOUND
            INAMES = NODEINFO(ND)%KNAMES
            NROWS  = NODEINFO(ND)%NROW
!
            DO J=1,NROWS
               LENR            = T_ROWS(IROW+J)%LENGTH
               DROW            = T_ROWS(IROW+J)%NAME
               VARTYP(IVTYP+J) = T_ROWS(IROW+J)%ROWTYP
               XLB(IBOUND+J)   = T_ROWS(IROW+J)%LOWER
               XUB(IBOUND+J)   = T_ROWS(IROW+J)%UPPER
               NAME1(INAMES+J) = LASTCH + 1
               DO I=1,LENR
                  QNAMES(I+LASTCH) = QROW(I)
               END DO
               LASTCH = LASTCH + LENR
            END DO
            NAME1(INAMES+NROWS+1) = LASTCH + 1
!
            DO JJ=1,NPER
               LMNS(JJ) = 0
               LAUX(1,JJ) = 1
            END DO
            RETURN
!
 9060 CONTINUE
         WRITE (IOLOG, 3060)
         IERR = 2
         RETURN
!
 3060 FORMAT(' XXX -  FATAL  - Unable to allocate memory in INCORE')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE NEWNODE(IERR)
!
!     This subroutine adds a new node into the event tree and initials
!     the pointer array.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
            MNRTYP = NODEINFO(INODE)%KROW + NODEINFO(INODE)%NROW
            IF (MNRTYP .GT. IZRTYP) THEN
               IF (.NOT. EXTMEM_RTYP(MNRTYP,IERR)) GOTO 9060
            ENDIF
!
            IF (MAXROW .GE. IZVNAM) THEN
               IF (.NOT. EXTMEM_VNAM(MAXROW+1,IERR)) GOTO 9060
            ENDIF
!
!     Add the new node(s) into the event tree and initial pointers.
!
            IPREV         = INODE
            INODE         = INODE + 1
            IRNGE0(IPREV) = IPREV
!
            NODEINFO(IPREV)%IDESC  = INODE
            NODEINFO(IPREV)%NDESC  = 1
!
            NODEINFO(INODE)%PROB   = 1.0D0
            NODEINFO(INODE)%IANCTR = IPREV
            NODEINFO(INODE)%IDESC  = 0
            NODEINFO(INODE)%IBROTH = 0
            NODEINFO(INODE)%NDESC  = 0
            NODEINFO(INODE)%KROW   = NODEINFO(IPREV)%KROW + NROWS
            NODEINFO(INODE)%KRHS   = NODEINFO(IPREV)%KRHS + NROWS
!
            IROW   = NODEINFO(INODE)%KROW
            MAXROW = MAXROW + 1
            T_ROWS(IROW+1)%NAME   = OBJNAM
            T_ROWS(IROW+1)%LENGTH = LENOBJ
            T_ROWS(IROW+1)%ROWTYP = 2
            T_ROWS(IROW+1)%LOWER  = -PLINF
            T_ROWS(IROW+1)%UPPER  =  PLINF
            NODEINFO(INODE)%NROW  = 1
            NODEINFO(INODE)%NCOL  = 1
            RETURN
!
 9060 CONTINUE
         WRITE (IOLOG, 3060)
         IERR = 2
         RETURN
!
 3060 FORMAT(' XXX -  FATAL  - Unable to allocate memory in INCORE')
!
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CPBLKS(IERR)
!
!     Copy subdiagonal matrices of current node
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DBLANK, DROW
      CHARACTER*1 DB255(NAMLEN), QROW(NAMLEN)
!
      EQUIVALENCE (DBLANK,DB255), (DROW, QROW)
!
      DATA DB255 /NAMLEN*' '/
!
!     Check if the problem structure has changed.
!
               MVBLKS = 0
               IF (MARKOV) THEN
                  DO J=IPER+2,NPER
                     IF (LMNS(J) .GT. 0) MVBLKS = 1
                  END DO
               ENDIF
!
               IF (ALLOW_GC .AND..NOT. HAVE_GC) THEN
                  DO J=1,IPER-1
                     IF (LMNS(J) .GT. 0 .AND. MVBLKS .LT. 4)      &
     &                  MVBLKS = MVBLKS + 2
                  END DO
               ENDIF
!
               IF (MVBLKS .GT. 0) CALL SHIFTA(MVBLKS, IERR)
               IF (IERR .GT. 0) GOTO 9999
!
!     Copy subdiagonal matrices of the current node
!
  320          CONTINUE
                  NCOLA = NCOLS - NODEINFO(IPER)%NROW
                  IF (MARKOV .AND. IPER .LT. NPER) THEN
                     MPER = IPER + 1
                  ELSE
                     MPER = NPER
                  ENDIF
                  DO JP=IPER+1,MPER
                     JBLK = KDATA(JP) + JP + 1 - IPER
                     KCOLA(JBLK) = LASTCA
                     KELMA(JBLK) = LASTA
                     NELMA(JBLK) = LMNS(JP)
!
                     MNALMN = LASTA + LMNS(JP)
                     IF (MNALMN .GT. IZALMN) THEN
                        IF (.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 9060
                     ENDIF
!
                     MNACOL = LASTCA + NCOLA + 1
                     IF (MNACOL .GT. IZACOL) THEN
                        IF (.NOT. EXTMEM_ACOL(MNACOL,IERR)) GOTO 9060
                     ENDIF
!
                     DO JC=1,NCOLA
                        LL = LAUX(JC,JP)
                        KK = LAUX(JC+1,JP) - 1
                        LA(LASTCA+JC) = LL
                        DO JR=LL,KK
                           IA(LASTA+JR) = IAUX(JR,JP)
                            A(LASTA+JR) =  AUX(JR,JP)
                        END DO
                     END DO
                     LASTA  = LASTA  + LMNS(JP)
                     LASTCA = LASTCA + NCOLA + 1
                     LA(LASTCA) = KK + 1
                  END DO
!
!     If there are global constraints, add them as well.
!
                  IF (HAVE_GC) THEN
                     DO JP=1,IPER-1
                        IF (MARKOV) THEN
                           JBLK = KDATA(IPER) + JP + 2
                        ELSE
                           JBLK = KDATA(IPER) + JP + IPER
                        ENDIF
                        KCOLA(JBLK) = LASTCA
                        KELMA(JBLK) = LASTA
                        NELMA(JBLK) = LMNS(JP)
!
                        MNALMN = LASTA + LMNS(JP)
                        IF (MNALMN .GT. IZALMN) THEN
                        IF (.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 9060
                        ENDIF
!
                        MNACOL = LASTCA + NCOLA + 1
                        IF (MNACOL .GT. IZACOL) THEN
                        IF (.NOT. EXTMEM_ACOL(MNACOL,IERR)) GOTO 9060
                        ENDIF
!
                        DO JC=1,NCOLA
                           LL = LAUX(JC,JP)
                           KK = LAUX(JC+1,JP) - 1
                           LA(LASTCA+JC) = LL
                           DO JR=LL,KK
                              IA(LASTA+JR) = IAUX(JR,JP)
                               A(LASTA+JR) =  AUX(JR,JP)
                           END DO
                        END DO
                        LASTA  = LASTA  + LMNS(JP)
                        LASTCA = LASTCA + NCOLA + 1
                        LA(LASTCA) = KK + 1
                     END DO
                  END IF
  330          CONTINUE
               RETURN
!
 9060 CONTINUE
         WRITE (IOLOG, 3060)
 9999 CONTINUE
         IERR = 2
         RETURN
!
 3060 FORMAT(' XXX -  FATAL  - Unable to allocate memory in INCORE')
!
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE SHIFTA(MVBLKS, IERR)
!
!     This subroutine shifts the pointers for the blocks of the A-matrix
!     when the structure changes (from staircase to full lower-block
!     triangular or to include linking constraints or both).
!
!     --------------------------------------------------------------
!     This version dated 21 September 2002. Written by Gus Gassmann.
!     --------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DBLANK, DROW, DCOL
      CHARACTER*1 DB255(NAMLEN),QROW(NAMLEN),QCOL(NAMLEN)
!
      EQUIVALENCE (DBLANK,DB255), (DROW,QROW), (DCOL,QCOL)
!
      DATA DB255 /NAMLEN*' '/
               IF (MVBLKS .GT. 0) THEN
!
! -------------------------------------------------------------------------
!
!     Here we have to make the switch. This involves shifting pointers KDATA,
!     KCOLA, KELMA and NELMA. There are four possible combinations:
!
!        Markov, no global constraints      Non-markov, no global constraints
!                 1                               1
!                 3   2                           3   2
!                     5   4                       6   5   4
!                         7   6                   10  9   8   7
!                             ...                  ...        ...
!
!        Markov, global constraints         Non-markov, global constraints
!
!                 1   4   7   11  ...             1   4   8   14  ...
!                 3   2   8   12                  3   2   9   15
!                     6   5   13                  7   6   5   16
!                        10   9                   13  12  11  10
!                             ...                  ...        ...
!
!     Altogether this defines five different update paths, based on MVBLKS.
!
!        MVBLKS = 1: from staircase to non-Markov, no global constraints
!        MVBLKS = 2: from staircase to global constraints
!        MVBLKS = 3: from staircase to non-Markov, global constraints
!        MVBLKS = 4: add global constraints to non-Markovian problem
!        MVBLKS = 5: add lower triangle to global constraints
!
! -------------------------------------------------------------------------
!
!     Check array sizes and enlarge if necessary
!
                  IF     (MVBLKS .EQ. 1) THEN
                     MARKOV = .FALSE.
                     MINPER = 4
                     IF (HAVE_GC) THEN
                        MNABLK = NPER * NPER
                        MVBLKS = 5
                     ELSE
                        MNABLK = NPER*(NPER+1)/2
                     ENDIF
                  ELSEIF (MVBLKS .EQ. 2) THEN
                     HAVE_GC = .TRUE.
                     MINPER  = 3
                     IF (MARKOV) THEN
                        MNABLK = NPER*(NPER+3)/2 - 1
                     ELSE
                        MVBLKS = 4
                        MNABLK = NPER * NPER
                     ENDIF
                  ELSE
                     MARKOV  = .FALSE.
                     HAVE_GC = .TRUE.
                     MINPER  = 3
                     MNABLK  = NPER * NPER
                  ENDIF
                  IF (MNABLK .GT. IZABLK) THEN
                     IF (.NOT. EXTMEM_ABLK(MNABLK,IERR)) GOTO 9060
                  ENDIF
!
!     First we shift the existing pointers...
!
                     DO JP=NPER,MINPER,-1
                        IF (MVBLKS .EQ. 5) THEN
                           KOLD = JP*(JP+1)/2
                           KNEW = (JP-1)**2+JP
                           DO K=JP-1,1,-1
                              KCOLA(KNEW+K) = KCOLA(KOLD+K)
                              KELMA(KNEW+K) = KELMA(KOLD+K)
                              NELMA(KNEW+K) = NELMA(KOLD+K)
                           END DO
                        ENDIF
!
                        IF (MVBLKS .EQ. 1) THEN
                           KDATA(JP) = JP*(JP-1)/2
                           KOLD = 2*JP - 3
                           KNEW = JP*(JP-1)/2
                           KBLK = 2
                        ELSEIF (MVBLKS .EQ. 2) THEN
                           KDATA(JP) = JP*(JP+1)/2 - 2
                           KOLD = 2*JP - 3
                           KNEW = JP*(JP+1)/2 - 2
                           KBLK = 2
                        ELSEIF (MVBLKS .EQ. 3) THEN
                           KDATA(JP) = (JP-1)**2
                           KOLD = 2*JP - 3
                           KNEW = (JP-1)**2
                           KBLK = 2
                        ELSEIF (MVBLKS .EQ. 4) THEN
                           KDATA(JP) = (JP-1)**2
                           KOLD = JP*(JP-1)/2
                           KNEW = (JP-1)**2
                           KBLK = JP
                        ELSEIF (MVBLKS .EQ. 5) THEN
                           KDATA(JP) = (JP-1)**2
                           KOLD = JP*(JP+1)/2 - 2
                           KNEW = (JP-1)**2
                           KBLK = 2
                        ENDIF
!
                        DO K=KBLK,1,-1
                           KCOLA(KNEW+K) = KCOLA(KOLD+K)
                           KELMA(KNEW+K) = KELMA(KOLD+K)
                           NELMA(KNEW+K) = NELMA(KOLD+K)
                        END DO
                     END DO
!
!     ... then we put dummy values for the newly created blocks.
!     The main reason to put dummies is to guard against empty blocks
!     (which wouldn't be initialized otherwise). When an element of a
!     particular block is found, we could check if it was previously
!     initialized and save some column pointers, but in this version
!     we don't take the trouble.          GG 21 September 1998
!
!     For MVBLKS = 1, 3, 5 we put subdiagonal blocks
!
                  IF (MVBLKS .NE. 2 .AND. MVBLKS .NE. 4) THEN
                     DO IPT=1,IPER-1
                        NCOLA = 1 + NODEINFO(IPT)%NCOL      &
     &                            - NODEINFO(IPT)%NROW
                        DO JN=IPT+2,NPER
                           IF (MVBLKS .EQ. 1) THEN
                              K0 = JN*(JN-1)/2
                           ELSE
                              K0 = (JN-1)**2
                           ENDIF
                           K = JN + 1 - IPT
                              KCOLA(K0+K) = LASTCA
                              KELMA(K0+K) = LASTA
                              NELMA(K0+K) = 0
                              MNACOL = LASTCA + NCOLA
                              IF (MNACOL .GT. IZACOL) THEN
                                 IF (.NOT. EXTMEM_ACOL(MNACOL,IERR))      &
     &                              GOTO 9060
                              ENDIF
                              DO JC=1,NCOLA
                                 LA(LASTCA+JC) = 1
                              END DO
                              LASTCA = LASTCA + NCOLA
                           END DO
                        END DO
                     END IF
!
!     For MVBLKS = 2, 3, 4 we put superdiagonal blocks for linking constraints
!
                     IF (MVBLKS .NE. 1 .AND. MVBLKS .NE. 5) THEN
                        DO JN=2,IPER-1
                           IF (MVBLKS .EQ. 2) THEN
                              K0 = JN*(JN+1)/2
                           ELSE
                              K0 = (JN-1)**2 + JN
                           ENDIF
                           NCOLA = 1 + NODEINFO(IPER)%NCOL      &
     &                               - NODEINFO(IPER)%NROW
                           DO K=1,JN-1
                              KCOLA(K0+K) = LASTCA
                              KELMA(K0+K) = LASTA
                              NELMA(K0+K) = 0
                              MNACOL = LASTCA + NCOLA
                              IF (MNACOL .GE. IZACOL) THEN
                                 IF (.NOT. EXTMEM_ACOL(MNACOL,IERR))      &
     &                              GOTO 9060
                              ENDIF
                              DO JC=1,NCOLA
                                 LA(LASTCA+JC) = 1
                              END DO
                              LASTCA = LASTCA + NCOLA
                           END DO
                        END DO
                     END IF
               END IF
               RETURN
!
 9060 CONTINUE
         WRITE (IOLOG, 3060)
         IERR = 2
         RETURN
!
 3060 FORMAT(' XXX -  FATAL  - Unable to allocate memory in INCORE')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INQUAD ( NEXT, NREC, IERR )
!
!     Routine: INQUAD
!
!     Purpose: To read a QSECTION giving a quadratic objective.
!
!     Calling sequence: CALL INQUAD ( NEXT, NREC, IERR )
!
!           NEXT - = 6 if QSECT is inline, = 7 otherwise
!           NREC - Number of the current record
!           IERR - Error status
!
!     Routines called: GETNAM, GETMPS
!
!
!     Method:
!     ===============
!
!     This subroutine reads a QSECTION in the format specified in the
!     on-line OSL documentation (www.research.ibm.com/osl/ekkgc10.htm).
!     Althought the documentation does not make it clear, it seems
!     that the Q-matrix is supposed to be symmetric and only the
!     main diagonal and the triangle below are defined and stored.
!     The QSECTION has a name record that is checked against the other
!     problem components and may be in the core file following the
!     ENDATA record for the linear part, or it may be in a separate
!     file.
!
!     The data elements are first entered into an array of linked lists
!     and transferred to standard sparse column form later. We also keep
!     track of the "profile" of the Q-matrix. In order to maintain
!     separability, the profile should be 1, that is, the only blocks
!     should be along the main diagonal, but it is possible to specify
!     nonseparable problems. If the profile equals 2, there is linkage
!     between a period and its immediate predecessor and successor, if
!     the profile equals 3, there is linkage between periods two steps
!     apart, and so on. The maximal linkage is tracked by the variable
!     NQPROF. (That is, NQPROF = 3 means linkage up to two steps apart.)
!
!     Nonzero elements are stored in the Q-matrix block by block. For
!     each block we record the offsets in the column pointers (array
!     LQMTX) in LQOFF, and the offsets in the row pointers (IQMTX) and
!     the value array (AQMTX) in IQOFF. The order of the blocks is
!     similar to the order of the A-matrix blocks: For each node in
!     the event tree the pointer array KDATQ points at the first block
!     for that node. In order to allow stochastic quadratic costs (why
!     would anybody want to use those?), blocks off the main diagonal
!     are associated with the later of the two periods involved. The
!     first block (in position KDATQ(node)+1) is the block on the main
!     diagonal, followed by additional blocks below the diagonal ---
!     at most min(IPER,NQPROF) --- in increasing distance from the main
!     diagonal.
!
!     Example:
!     ========
!
!     NODE = 8, IPER = 4 (node 8 is in period 4), NQPROF=3, and KDATQ(8) = 5.
!
!     There are three blocks associated with this node, and their numbers
!     are found in the following block schematic:
!
!     period     1   2   3   4   5   ...
!               --- --- --- --- ---
!        1     |   |   |   |   |   | ...
!               --- --- --- --- ---
!        2     |   |   |   |   |   | ...
!               --- --- --- --- ---
!        3     |   |   |   |   |   | ...
!               --- --- --- --- ---
!        4     |   | 8 | 7 | 6 |   | ...
!               --- --- --- --- ---
!        5     |   |   |   |   |   | ...
!               --- --- --- --- ---
!       ...     ... ... ... ... ...
!
!     ---------------------------------------------------------------
!     This version dated September 28, 2002. Written by Gus Gassmann.
!     ---------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      REAL*8 ATEMP(4)
      CHARACTER*(NAMLEN) DNAME(3), QLINE, DBLANK, DCOL1, DCOL2
      CHARACTER*1 DB255(NAMLEN)
      INTEGER*4   LENGTH(3)
!
      EQUIVALENCE (DBLANK,DB255)
!
      DATA DB255/NAMLEN*' '/
!
!     Allocate temporary arrays
!
      ALLOCATE (QLINKS(IZQLMN), QRTMP(IZQLMN), QCTMP(IZQLMN),      &
     &          QMTMP(IZQLMN), ISTART(NPER,NPER), STAT=IER)
      IF (IER .NE. 0) GOTO 9060
!
!     INITIALIZE POINTERS
!
      DO J=1,NPER
         DO K=1,NPER
            ISTART(J,K) = 0
         END DO
      END DO
!
      DCOL1  = DBLANK
      DCOL2  = DBLANK
      IP0    = 1
      IC0    = 1
      LMN    = 0
      IERR   = 0
      IENDAT = 0
      IF (NEXT .EQ. 6) THEN
         IF (IOQUA .NE. IOCOR) THEN
            WRITE (IOLOG, 1500)
            IOQUA = IOCOR
         ENDIF
         GOTO 150
      ELSE
         NREC = 0
      ENDIF
!
!     Open the QSECTION file.
!
      IF (IOQUA .NE. IOCOR) THEN
         CALL FILE_WRAPPER(QUAFIL,IOQUA,3,2,IOLOG,IERR)
         IF (IERR .GT. 0) GOTO 142
      ENDIF
!
!     Get the name of the QSECTION first
!
  140 CONTINUE
         CALL GETNAM(QLINE,Q1,Q2,Q3,Q4,      &
     &                     DNAME,ATEMP,IERR,IOQUA,NREC,LENGTH)
         IF (Q1 .EQ. QAST ) GOTO 140
         IF (IERR   .GT. 0) GOTO 142
         IF (Q1 .EQ. QN .AND. Q2 .EQ. QA) GOTO 145
         IF (Q1 .EQ. QQ .AND. Q2 .EQ. QS) GOTO 9150
         IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 400
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 2300)
            GOTO 9999
!
!     Error during read. Treat as missing and continue
!
  142 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3000)
         GOTO 400
!
!     We have found a name. Check that it matches the core file.
!
  145 CONTINUE
         IF (DNAME(1) .NE. PROBNM) GOTO 9150
!
!     READ statement for the rest of the QSECTION
!
  150 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,      &
     &                     DNAME,ATEMP,IERR,IOQUA,NREC,LENGTH)
         IF (IERR .GT.  0) GOTO 152
         IF (Q1 .EQ. QAST) GOTO 150
         IF (Q1 .EQ. QBL ) GOTO 160
         IF (Q1 .EQ. QQ .AND. Q2 .EQ. QS) GOTO 151
         IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 400
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 2300)
            GOTO 9999
!
!     QSECT record. There could be more than one, in which case the rest
!     of the file is ignored
!
  151 CONTINUE
         IF (IENDAT .EQ. 1) GOTO 400
            IENDAT = 1
            IF (NECHO1 .GE. 1) WRITE (IOLOG, 1400) NREC,QLINE
            GOTO 150
!
!     Here we have an error after the QSECT header has been read.
!     This is recoverable, but the user must be informed. Treat like EOF.
!
  152 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3000)
         GOTO 400
!
!    Here we have a regular data card.
!
  160 CONTINUE
         IF (NECHO1 .GE. 2) WRITE (IOLOG, 1400) NREC,QLINE
         IF (DABS(ATEMP(1)) .LT. ZTOLIN) GOTO 150
         IF (DNAME(1) .EQ. DCOL1) GOTO 340
!
!     Match the first column name
!
         DCOL1 = DNAME(1)
         IP = IP0
         KN = NODEINFO(IP)%KNAMES
         NR = NODEINFO(IP)%NROW
         NC = NODEINFO(IP)%NCOL
         IRES = NMSRCH(DCOL1,KN+NR+IC0,KN+NC)
         IF (IRES .GT. 0) GOTO 330
         IRES = NMSRCH(DCOL1,KN+NR+1,KN+NR+IC0)
         IF (IRES .GT. 0) GOTO 330
!
         DO IP=IP0+1,NPER
            KN = NODEINFO(IP)%KNAMES
            NR = NODEINFO(IP)%NROW
            NC = NODEINFO(IP)%NCOL
            IRES = NMSRCH(DCOL1,KN+NR+1,KN+NC)
            IF (IRES .GT. 0) GOTO 330
         END DO
         DO IP=1,IP0-1
            KN = NODEINFO(IP)%KNAMES
            NR = NODEINFO(IP)%NROW
            NC = NODEINFO(IP)%NCOL
            IRES = NMSRCH(DCOL1,KN+NR+1,KN+NC)
            IF (IRES .GT. 0) GOTO 330
         END DO
!
!     Not found!
!
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 1200) DCOL1
         GOTO 9999
!
!     Matched
!
  330 CONTINUE
         IP0 = IP
         IC0 = IRES - KN - NR
!
!     Match the second column name
!
  340 CONTINUE
         DCOL2 = DNAME(2)
         IP = IP0
         KN = NODEINFO(IP)%KNAMES
         NR = NODEINFO(IP)%NROW
         NC = NODEINFO(IP)%NCOL
         IRES = NMSRCH(DCOL2,KN+NR+1,KN+NC)
         IF (IRES .GT. 0) GOTO 350
         DO IP=1,NPER
            KN = NODEINFO(IP)%KNAMES
            NR = NODEINFO(IP)%NROW
            NC = NODEINFO(IP)%NCOL
            IRES = NMSRCH(DCOL2,KN+NR+1,KN+NC)
            IF (IRES .GT. 0) GOTO 350
         END DO
!
!     Not found!
!
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 1200) DCOL2
         GOTO 9999
!
!   Matched. Verify that the element is in the lower triangular part,
!   then store in temporary data structure (linked lists). If in the
!   upper half, warn the user and ignore the entry.
!
  350 CONTINUE
      IR = IRES - KN - NR
      IF (IP0 .LT. IP .OR. (IP .EQ. IP0 .AND. IC0 .LE. IR)) THEN
         LMN = LMN + 1
         IF (LMN .GT. IZQLMN) THEN
            IF (.NOT. EXTMEM_QLMN(LMN,IERR)) GOTO 9060
         ENDIF
!
         QRTMP(LMN) = IR
         QCTMP(LMN) = IC0
         QMTMP(LMN) = ATEMP(1)
         IF (IABS(IP0-IP) .GE. NQPROF) NQPROF = IABS(IP0-IP) + 1
         IF (ISTART(IP0,IP) .EQ. 0) THEN
            ISTART(IP0,IP) = LMN
            QLINKS(LMN) = 0
         ELSE
            IL = ISTART(IP0,IP)
            IF (QCTMP(IL) .GE. IC0) THEN
               ISTART(IP0,IP) = LMN
               QLINKS(LMN) = IL
            ELSE
  360       CONTINUE
               IF (QLINKS(IL) .EQ. 0) THEN
                   QLINKS(IL)  = LMN
                   QLINKS(LMN) = 0
               ELSEIF (QCTMP(QLINKS(IL)) .GE. IC0) THEN
                   QLINKS(LMN) = QLINKS(IL)
                   QLINKS(IL)  = LMN
               ELSE
                   IL = QLINKS(IL)
                   GOTO 360
               ENDIF
            ENDIF
         ENDIF
      ELSE
         WRITE (IOLOG, 1400) NREC, QLINE
         WRITE (IOLOG, 1300)
      ENDIF
      GOTO 150
!
!     End of file. Copy data to permanent storage positions.
!     Start by computing the exact storage requirements. If resizing is
!     permitted, this may recover some storage.
!
  400 CONTINUE
         IF (NQPROF .GT. 0) THEN
            MNQBLK = NPER*NQPROF - NQPROF*(NQPROF-1)/2 + 1
            MNQLMN = LMN
            MNQCOL = 0
            DO I=1,NPER
               MNQCOL = MNQCOL + MIN(NQPROF, 1+NPER-I)*      &
     &                     (NODEINFO(I)%NCOL - NODEINFO(I)%NROW + 1)
            END DO
            IF (RESIZE) THEN
               IF ((MNQBLK .GT. MXQBLK) .OR. (MNQLMN .GT. MXQLMN)      &
     &                                  .OR. (MNQCOL .GT. MXQCOL)) THEN
                  GOTO 9060
               ELSE
                  IZQBLK = MNQBLK
                  IZQCOL = MNQCOL
                  IZQLMN = MNQLMN
               ENDIF
            ELSE
               IF ((MNQBLK .GT. IZQBLK) .OR. (MNQLMN .GT. IZQLMN)      &
     &                                  .OR. (MNQCOL .GT. IZQCOL)) THEN
                  GOTO 9060
               ENDIF
            ENDIF
!
            ALLOCATE(AQMTX(IZQLMN), IQMTX(IZQLMN), LQMTX(IZQCOL),      &
     &               KDATQ(IZNODE), IQOFF(IZQBLK), LQOFF(IZQBLK),      &
     &               NELMQ(IZQBLK), STAT=IER)
            IF (IER .NE. 0) GOTO 9060
!
            LOFF  = 0
            NQBLK = 0
            NQLMN = 0
!
!     Main diagonal is first
!
            DO I=1,NPER
               KDATQ(I) = NQBLK
               NQBLK    = NQBLK + 1
               IQOFF(NQBLK) = NQLMN
               LQOFF(NQBLK) = LOFF
               NELMQ(NQBLK) = NQLMN
               LC = 0
               IF (ISTART(I,I) .GT. 0) THEN
                  IL = ISTART(I,I)
  410          CONTINUE
                  IC = QCTMP(IL)
                  DO K=LC+1,IC
                     LQMTX(LOFF+K) = NQLMN + 1 - IQOFF(NQBLK)
                  END DO
                  NQLMN = NQLMN + 1
                  AQMTX(NQLMN) = QMTMP(IL)
                  IQMTX(NQLMN) = QRTMP(IL)
                  IL = QLINKS(IL)
                  LC = IC
                  IF (IL .GT. 0) GOTO 410
               ENDIF
               NELMQ(NQBLK)  = NQLMN - NELMQ(NQBLK)
               NC = NODEINFO(I)%NCOL - NODEINFO(I)%NROW
               DO K=LC,NC
                  LQMTX(LOFF+K+1) = NQLMN + 1 - IQOFF(NQBLK)
               END DO
               LOFF = LOFF + NC + 1
               J = I
!
!     Now deal with subdiagonals, if there are any
!
               IF (MIN(I,NQPROF) .GT. 1) THEN
                  J = I - 1
  425          CONTINUE
                  NQBLK = NQBLK + 1
                  IQOFF(NQBLK)  = NQLMN
                  LQOFF(NQBLK)  = LOFF
                  NELMQ(NQBLK)  = NQLMN
                  LQMTX(LOFF+1) = NQLMN + 1 - IQOFF(NQBLK)
                  LC = 0
                  IF (ISTART(J,I) .GT. 0) THEN
                     IL = ISTART(J,I)
  430             CONTINUE
                     IC = QCTMP(IL)
                     DO K=LC+1,IC
                        LQMTX(LOFF+K) = NQLMN + 1 - IQOFF(NQBLK)
                     END DO
                     NQLMN = NQLMN + 1
                     AQMTX(NQLMN) = QMTMP(IL)
                     IQMTX(NQLMN) = QRTMP(IL)
                     IL = QLINKS(IL)
                     LC = IC
                     IF (IL .GT. 0) GOTO 430
                  ENDIF
                  NELMQ(NQBLK)  = NQLMN - NELMQ(NQBLK)
                  NC = NODEINFO(J)%NCOL - NODEINFO(J)%NROW
                  DO K=LC,NC
                     LQMTX(LOFF+K+1) = NQLMN + 1 - IQOFF(NQBLK)
                  END DO
                  LOFF = LOFF + NC + 1
!
                  J = J - 1
                  IF (J .GT. 0 .AND. I-J .LT. NQPROF) GOTO 425
               ENDIF
            END DO
!
            LASTQB = NQBLK
            LASTQC = LOFF
            LASTQ  = NQLMN
            IQOFF(NQBLK+1) = NQLMN
         ENDIF
!
!     END OF FILE
!
  900 CONTINUE
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 1400) NREC,QLINE
         GOTO 999
!
 9060 CONTINUE
         WRITE (IOLOG, 3060)
         IERR = 2
         GOTO 999
!
 9150 CONTINUE
         WRITE (IOLOG, 3150)
         IERR = 2
         GOTO 999
!
 9999 CONTINUE
         IERR = 2
!
  999 CONTINUE
      DEALLOCATE(QLINKS, QRTMP, QCTMP, ISTART, QMTMP, STAT=IER)
      IF (IER .NE. 0) WRITE (IOLOG, 3061)
      CLOSE (IOQUA)
      RETURN
!
 1200 FORMAT(' XXX -  FATAL  - Unmatched column name ',A8,      &
     &       ' in QSECTION')
 1300 FORMAT(' XXX - WARNING - Q-matrix element in upper triangle.'      &
     &       ' Ignore.')
 1400 FORMAT(I8,2X,A80)
 1500 FORMAT(' XXX - WARNING - Inline QSECT must be read from IOCOR.',      &
     &       ' Reset I/O channel.')
 2300 FORMAT(' XXX -  FATAL  - Illegal header card in QSECTION.')
 3000 FORMAT(' XXX - WARNING - Missing ENDATA card in QSECTION.')
 3060 FORMAT(' XXX -  FATAL  - Memory allocation error in INQUAD.')
 3061 FORMAT(' XXX - WARNING - Memory allocation error in INQUAD.')
 3150 FORMAT(' XXX -  FATAL  - Problem name not matched.')
      END SUBROUTINE
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INSTOC ( IERR )
!
!     Routine: INSTOC
!
!     Purpose: To verify the name in the stoch file against the name in the time
!              and core files and to determine the format of the event tree.
!              Based on this format, further routines are called.
!
!     Calling sequence: CALL INSTOC ( IERR )
!
!        IERR   Error status
!         = 0   Normal termination
!         > 0   Error condition
!
!     Calls:    CHOPEN, GETNAM, GETMPS, INSCEN, INTREE, INCONT.
!
!     Method:   The name on the NAME card is compared to the info read in the
!               time and core files. Subsequent headers define the distributions
!               of the random variables that define the event tree, as set out
!               in the paper by Gassmann and Schweitzer.
!
!               The format could be explicit, as flagged by the header NODES
!               or SCENARIOS, or implicit if any other header is present.
!               Depending on the format, different routines are called to
!               read the rest of the stoch file.
!
!               The order of sections in the stoch file is as follows:
!               (these sections can be given in any order)
!                  SIMPLE    - for simple recourse
!                  ROBUST    - for robust recourse (quadratic penalty)
!                  PLINQUAD  - for piecewise linear-quadratic penalty
!                  CHANCE    - defines chance constraints
!                  ICC       - defines integrated chance constraints
!
!               These sections are followed by the definition of the event
!               tree. There are three mutually exclusive forms:
!                  NODES     - explicit tree in node form
!                  SCENARIOS - explicit tree in scenario form
!                  other     - implicit tree, using these three headers:
!                     INDEP   - independent (univariate) random variables
!                     BLOCKS  - random vector
!                     DISTRIB - distribution not linked directly to matrix
!                               coefficients
!
!     Remarks:
!
!     Before specifying the distributions of the random variables,
!     additional sections for chance constraints and simple recourse may
!     be processed.
!
!     Each data line in the SIMPLE section contains an identifier in the first
!     name field, a valid row name in the second, and the cost coefficients in
!     the two numeric fields; $q^-$ is first, followed by $q^+$. The penalty
!     costs $q^-$ and $q^+$ must satisfy the condition $q^+ -q^- \ge 0$.
!     Each occurrence of a row name in this section triggers the generation of
!     additional decision variables, identified by appending the suffix
!    `.minus' and `.plus' after the name of the row. The recourse coefficients
!     ought to be deterministic, but in a pinch the user can reference the
!     generated column names in a distribution section later. For that reason
!     the SIMPLE section precedes the INDEP, BLOCKS and DISTRIB sections.

!     A ROBUST section can be specified to set up quadratic penalties. The
!     syntax is identical to the SIMPLE section. In keeping with the conventions
!     for quadratic objectives, the penalty term includes an implicit factor
!     of 2, that is, if the specified coefficient is 1.5, the penalty function
!     is 0.75*(x^2).

!     Finally the program recognizes a PLINQUAD section for piecewise linear-
!     quadratic penalties. These are penalties of the form
!
!        p(t) = a/2 t^2          if 0 <= t <= b/a
!               bt - (b^2)/2a    if b/a < t
!
!     Each data line in this section lists an identifying name in the first name
!     field, a valid row name in the second, the quadratic parameter $a$ in the
!     first numeric field, and the linear parameter $b$ in the second. (If $b$
!     is negative, the penalty is assumed to be purely quadratic, if $a$ is
!     negative, the penalty is assumed to be purely linear.)
!
!     If the mentioned row is of type `L', the penalty is applied to positive
!     deviations, if the row is of type `G', the penalty applies to negative
!     deviations. In the case of equality constraints, the range can be
!     clarified in the code field. If the code field is empty, the penalty is
!     applied to both positive and negative deviations; the code `PL' denotes
!     positive deviations, `MI' is used to specify negative deviations.
!
!     Probabilistic constraints are stored as follows: For each node in the
!     event tree NCHANCE(node) records the number of probabilistic
!     constraints, KCHANCE(node) acts as a pointer into array XCHANCE.
!     XCHANCE(i) contains the type, probability and dimension of the ith
!     chance constraint, along with pointers into RCHANCE and QCHANCE.
!     RCHANCE gives a list of rows associated with the ith constraint,
!     QCHANCE gives an identifying name for multidimensional constraints.
!     The total number of chance constraints is recorded in LCHANCE.
!
!     Integrated chance constraints are stored in a similar fashion.
!     For each node in the tree NICC(node) records the number of ICCs,
!     KICC(node) acts as a pointer into array XICC, which contains the
!     type and right-hand side of the ith ICC. ICCs are one-dimensional,
!     so there is no need for row and name arrays. The total number
!     of ICCs is recorded in LICC.
!
!  --------------------------------------------------------------------
!     This version dated September 28, 2002. Written by Gus Gassmann.
!  --------------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME(3), DBLANK, QLINE, DTEMP
      CHARACTER*30 F1500
      CHARACTER*1 QBLANK(NAMLEN), QTEMP(NAMLEN)
      REAL*8      ATEMP(4)
      INTEGER*4   LENGTH(3)
!
      EQUIVALENCE (QBLANK,DBLANK), (QTEMP,DTEMP)
      DATA        QBLANK/NAMLEN*' '/
!
! ---------------------------------------------------------------------
!
!     Some initializations at the start
!
!
      IERR    = 0
      NREC    = 0
      NRVAR   = 0
      NSTELM  = 0
      STOSTAT = 0
      NDCOLM  = 0
      NDELEM  = 0
      NSTELM  = 0
      LCHANCE = 0
      NJOINT  = 0
!
      IF (CORSTAT .EQ. 2 .OR. TIMSTAT .EQ. 2) GOTO 9980
!
      IF (NECHO1  .GE. 1) WRITE (IOLOG, 1900)
!
!     Open the stoch file
!
      CALL FILE_WRAPPER(STOFIL,IOSTO,3,2,IOLOG,IERR)
      IF (IERR .GT. 0) GOTO 340
!
!     Look for the STOCH card and the problem name
!
  310 CONTINUE
         CALL GETNAM(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOSTO,NREC,LENGTH)
         IF (Q1 .EQ. QAST) GOTO 310
         IF (IERR  .GT. 0) GOTO 340
         IF (Q1 .EQ. QS .AND. Q2 .EQ. QT) GOTO 320
            WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &             DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
            WRITE (IOLOG, 1600)
            GOTO 9999
!
!     Verify the problem name
!
  320    CONTINUE
            IF (PROBNM .EQ. DOTS) PROBNM = DNAME(1)
            IF (DNAME(1) .NE. PROBNM) GOTO 9150
            IF (NECHO1 .GE. 1) WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &            DNAME(1),DNAME(2)
!
!     Read one more record to determine the next section
!
  330 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOSTO,NREC,LENGTH)
         IF (Q1 .EQ. QAST ) GOTO 330
         IF (IERR   .GT. 0) GOTO 340
!
!     Distribute for further processing
!
  335 CONTINUE
         IF (Q1 .EQ. QS .AND. Q2 .EQ. QI) GOTO 350
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QO) GOTO 355
         IF (Q1 .EQ. QP .AND. Q2 .EQ. QL) GOTO 360
         IF (Q1 .EQ. QC .AND. Q2 .EQ. QH) GOTO 365
         IF (Q1 .EQ. QI .AND. Q2 .EQ. QC) GOTO 370
         IF (Q1 .EQ. QS .AND. Q2 .EQ. QC) GOTO 400
         IF (Q1 .EQ. QN .AND. Q2 .EQ. QO) GOTO 500
         IF (Q1 .EQ. QI .AND. Q2 .EQ. QN) GOTO 600
         IF (Q1 .EQ. QB .AND. Q2 .EQ. QL) GOTO 600
         IF (Q1 .EQ. QD .AND. Q2 .EQ. QI) GOTO 600
         IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 900
            WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &             DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
            WRITE (IOLOG, 1600)
            GOTO 9999
!
!     Error during read or missing STOCH file - Keep going
!     (This means the problem is assumed to be deterministic)
!
  340    CONTINUE
            IF (NREC .GT. 0) WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
            WRITE (IOLOG, 3990)
            IF (CORSTAT .EQ. 1 .OR. CORSTAT .EQ. -1      &
     &                         .OR. PROBNM .EQ. DOTS) GOTO 9995
               IF (NPER .EQ. 0) NPER = 1
               NODES = NPER
               GOTO 900
!
!     Simple recourse section
!
  350    CONTINUE
            CALL INSIMP(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
            IF (IERR .NE. 0) GOTO 9999
            GOTO 335
!
!     Section for quadratic penalties (ROBUST)
!
  355    CONTINUE
            CALL INRBST(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
            IF (IERR .NE. 0) GOTO 9999
            GOTO 335
!
!     Section for piecewise linear-quadratic penalties
!
  360    CONTINUE
            CALL INPLQP(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
            IF (IERR .NE. 0) GOTO 9999
            GOTO 335
!
!     Chance-constraint section
!
  365    CONTINUE
            CALL INCHAN(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
            IF (IERR .NE. 0) GOTO 9999
            GOTO 335
!
!     Integrated chance-constraint section
!
  370    CONTINUE
            CALL INICC(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
            IF (IERR .NE. 0) GOTO 9999
            GOTO 335
!
!     Explicit event tree representation: SCENARIOS
!    ===============================================
!
  400    CONTINUE
            IF (CORSTAT .EQ. 1 .OR. CORSTAT .EQ. -1 .OR.      &
     &          TIMSTAT .EQ. 1 .OR. TIMSTAT .EQ. -1) GOTO 9995
            IF (NPSEEN  .NE. NPER) WRITE (IOLOG, 2000)
            NPER = NPSEEN
            IF (NECHO1  .GE.    2) WRITE (IOLOG, 1200)      &
     &         NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2)
!
!     Allocate arrays for scenario names and leaf nodes
!
            ALLOCATE (LEAFND(IZSCEN), KSCNAM(IZSCEN+1),      &
     &                QSCNAM(IZCHSC), STAT=IERR)
            IF (IERR .NE. 0) GOTO 9200
!
            CALL INSCEN( NREC, IERR)
            IF (IERR .GT. 1) GOTO 9999
               STOSTAT = IERR
               STFILE  = 2
               CHNG_BLKS = .TRUE.
!
!     Link together the nodes in the event tree (i.e., identify cousins)
!
      DO 430 IP=1,NPER
         ISC1 = IRNGE0(IP)
  410    CONTINUE
            IF (NODEINFO(ISC1)%IBROTH .GT. 0) THEN
               ISC1 = NODEINFO(ISC1)%IBROTH
               GOTO 410
            ENDIF
!
            IAN = NODEINFO(ISC1)%IANCTR
            IF (IAN .EQ. 0) GOTO 430
  420       CONTINUE
               IBRO = ABS(NODEINFO(IAN)%IBROTH)
               IF (IBRO .EQ. 0) GOTO 430
                  ISC2 = NODEINFO(IBRO)%IDESC
                  IF (ISC2 .GT. 0) THEN
                     NODEINFO(ISC1)%IBROTH = -ISC2
                     ISC1 = ISC2
                     GOTO 410
                  ELSE
                     IAN = IBRO
                     GOTO 420
                  ENDIF
  430 CONTINUE
!
!     Count descendants for each problem
!
      DO IP=1,NPER
         I = IRNGE0(IP)
  450    CONTINUE
            N0 = 0
            I0 = NODEINFO(I)%IDESC
  460       CONTINUE
               IF (I0 .GT. 0) THEN
                  I0 = NODEINFO(I0)%IBROTH
                  N0 = N0 + 1
                  GOTO 460
               ENDIF
               NODEINFO(I)%NDESC = N0
               I = ABS(NODEINFO(I)%IBROTH)
               IF (I .GT. 0) GOTO 450
      END DO
!
!     Find conditional probabilities
!
               DO IP=2,NPER
                  ISC1 = IRNGE0(NPER+2-IP)
  480             CONTINUE
                     IN = NODEINFO(ISC1)%IANCTR
                     PT = NODEINFO(ISC1)%PROB/NODEINFO(IN)%PROB
                     NODEINFO(ISC1)%PROB = PT
                     ISC1 = ABS(NODEINFO(ISC1)%IBROTH)
                     IF (ISC1 .GT. 0) GOTO 480
               END DO
               TREEFMT = 2
               GOTO 900
!
!     Explicit event tree representation: NODES
!    ===========================================
!
  500    CONTINUE
!
!     NODE representation is incompatible with global constraints
!     (because there is no way to recognize them in MKNODE and because
!     they don't seem to make sense in this context).
!
            IF (HAVE_GC) GOTO 9300
!
            NPER = NPSEEN
            MULTI = 1
            IF (CORSTAT .EQ. 1 .OR. CORSTAT .EQ. -1 .OR. NPER .EQ. 0)      &
     &         IERR = 1
            IF (NECHO1   .GE. 2) WRITE (IOLOG, 1200)      &
     &         NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2)
!
!     Allocate arrays for scenario names and leaf nodes
!
            ALLOCATE (LEAFND(IZSCEN), KSCNAM(IZSCEN+1),      &
     &                QSCNAM(IZCHSC), STAT=IERR)
            IF (IERR .NE. 0) GOTO 9200
!
            CALL INTREE( IERR, NREC)
            IF (IERR .GT. 1) GOTO 9999
               STOSTAT = IERR
               STFILE  = 1
               CHNG_BLKS = .FALSE.

               IF ((.NOT. RESIZE .AND. (NPER .GT. IZTPER))      &
     &                           .OR.  (NPER .GT. MXTPER) ) GOTO 9200
               IF (ALLOCATED(STOCHA)) THEN
                  DEALLOCATE(STOCHA, STAT=IERR)
                  IF (IERR .NE. 0) GOTO 9200
               ENDIF
               ALLOCATE(STOCHA(NPER,NPER), STAT=IERR)
               IF (IERR .NE. 0) GOTO 9200
!
               DO I=1,NPER
                  DO J=1,I
                     STOCHA(I,J) = .TRUE.
                  END DO
               END DO
!
!     Link together the nodes in the event tree (i.e., identify cousins)
!
      DO 550 IP=1,NPER
         IF (IP .EQ. 1) THEN
            ISC1 = IRNGE0(IP)
         ELSE
            J = IRNGE0(IP-1)
  525       CONTINUE
               ISC1 = NODEINFO(J)%IDESC
               IF (ISC1 .EQ. 0) THEN
                  J = ABS(NODEINFO(J)%IBROTH)
                  IF (J .GT. 0) GOTO 525
                     GOTO 9400
               ELSE
                  IRNGE0(IP) = ISC1
               ENDIF
         ENDIF
!
  530    CONTINUE
            IF (NODEINFO(ISC1)%IBROTH .GT. 0) THEN
               ISC1 = NODEINFO(ISC1)%IBROTH
               GOTO 530
            ENDIF
!
            IAN = NODEINFO(ISC1)%IANCTR
            IF (IAN .EQ. 0) GOTO 550
  540       CONTINUE
               IBRO = ABS(NODEINFO(IAN)%IBROTH)
               IF (IBRO .EQ. 0) GOTO 550
                  ISC2 = NODEINFO(IBRO)%IDESC
                  IF (ISC2 .GT. 0) THEN
                     NODEINFO(ISC1)%IBROTH = -ISC2
                     ISC1 = ISC2
                     GOTO 530
                  ELSE
                     IAN = IBRO
                     GOTO 540
                  ENDIF
  550 CONTINUE
!
!     Count descendants for each problem
!
      NSCHAR = 0
      LASTLF = 0
      KSCNAM(1) = 0
      DO IP=1,NPER
         I = IRNGE0(IP)
  560    CONTINUE
            N0 = 0
            I0 = NODEINFO(I)%IDESC
  570       CONTINUE
               IF (I0 .GT. 0) THEN
                  I0 = NODEINFO(I0)%IBROTH
                  N0 = N0 + 1
                  GOTO 570
               ENDIF
               NODEINFO(I)%NDESC = N0
               I = ABS(NODEINFO(I)%IBROTH)
               IF (I .GT. 0) GOTO 560
      END DO
!
!     Identify leaf nodes and generate a scenario name.
!
      ND = IRNGE0(1)
  585 CONTINUE
         IF (NODEINFO(ND)%IDESC .GT. 0) THEN
            ND = NODEINFO(ND)%IDESC
            GOTO 585
!
!     Store the information.
!     (The hardest part is to produce the FORMAT statement.)
!
         ELSE
            IF (LASTLF .GE. IZSCEN) THEN
               IF (.NOT. EXTMEM_SCEN(LASTLF+1,IERR)) GOTO 9200
            ENDIF
            LASTLF = LASTLF + 1
            LEAFND(LASTLF) = ND
            IF (LASTLF .LT. 10) THEN
               LEN = 8
               F1500 = '(''Path000'',I1)'
            ELSEIF (LASTLF .LT. 100) THEN
               LEN = 8
               F1500 = '(''Path00'',I2)'
            ELSEIF (LASTLF .LT. 1000) THEN
               LEN = 8
               F1500 = '(''Path0'',I3)'
            ELSEIF (LASTLF .LT. 10000) THEN
               LEN = 8
               F1500 = '(''Path'',I4)'
            ELSE
               LEN = IFIX(ALOG10(FLOAT(LASTLF))) + 1
               IF (LEN .LT. 10) THEN
                  WRITE (F1500, 1510) LEN
               ELSE
                  WRITE (F1500, 1520) LEN
               ENDIF
               LEN = LEN + 4
            ENDIF
            WRITE (DTEMP, F1500) LASTLF

            IF (NSCHAR+LEN .GT. IZCHSC) THEN
               IF (.NOT. EXTMEM_CHSC(NSCHAR+LEN,IERR)) GOTO 9200
            ENDIF
!
            DO J=1,LEN
               QSCNAM(NSCHAR+J) = QTEMP(J)
            END DO
            NSCHAR = NSCHAR + LEN
            KSCNAM(LASTLF+1) = NSCHAR
         ENDIF
!
!     Start the search for the next leaf
!
  595    CONTINUE
            IF (NODEINFO(ND)%IBROTH .GT. 0) THEN
               ND = NODEINFO(ND)%IBROTH
               GOTO 585
            ELSE
               ND = NODEINFO(ND)%IANCTR
               IF (ND .GT. 0) GOTO 595
            ENDIF
!
      IF (NECHO1 .EQ. 0) WRITE (IOLOG, 2800) PROBNM, NPER
      WRITE (IOLOG, 2100)
      TREEFMT = 2
      GOTO 900
!
!     Implicit event tree representation
!    ====================================
!
  600    CONTINUE
            IF ((TIMSTAT .EQ. 1 .OR. TIMSTAT .EQ. -1) .AND.      &
     &           CORSTAT .NE. 0) GOTO 9995
            IF (NECHO1 .GE. 1)      &
     &         WRITE (IOLOG, 2700) NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2)
!
!     Allocate arrays for scenario names and leaf nodes
!
            ALLOCATE (LEAFND(IZSCEN), KSCNAM(IZSCEN+1),      &
     &                QSCNAM(IZCHSC), STAT=IERR)
            IF (IERR .NE. 0) GOTO 9200
!
            CALL INCONT(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
            IF (IERR .GT. 1) GOTO 9999
               STOSTAT = IERR
               CHNG_BLKS = .TRUE.
!
!     Linking is not necessary, and counting is much easier for this
!     type of stochastic information, since there can only be one path.
!
               DO I=1,NPER-1
                  NODEINFO(I)%NDESC = 1
               END DO
               NODEINFO(NPER)%NDESC = 0
!
!     End of input
!    =============
!
  900    CONTINUE
            IF (NECHO1 .EQ. 0) WRITE (IOLOG, 2800) PROBNM, NPER
            WRITE (IOLOG, 2100)
!
!     Determine overall problem size
!
            IF (MULTI .EQ. -1) MULTI = 0
            NEXTND = NODES + 1
            NROWS  = NODEINFO(NODES)%NROW
            NCOLS  = NODEINFO(NODES)%NCOL
            MAXCOL = NODEINFO(NODES)%KCOL + NCOLS + 1
            MAXROW = NODEINFO(NODES)%KROW + NROWS
            MAXRHS = LASTR
!
            IF (MAXROW .GT. IZROWS) THEN
               IF (.NOT. RESIZE .OR. (IZROWS .GE. MXROWS)      &
     &                          .OR. (MAXROW .GT. MXROWS)) GOTO 9200
               IZROWS = MAXROW
            ENDIF
!
            IF (MAXCOL .GT. IZCOLS) THEN
               IF (.NOT. RESIZE .OR. (IZCOLS .GE. MXCOLS)      &
     &                          .OR. (MAXCOL .GT. MXCOLS)) GOTO 9200
               IZCOLS = MAXCOL
            ENDIF
!
            IF (NEXTND .GT. IZNODE) THEN
               IF (.NOT. EXTMEM_NODE(NEXTND,IERR)) GOTO 9200
            ENDIF
!
            NODEINFO(NEXTND)%KCOL   = MAXCOL
            NODEINFO(NEXTND)%KROW   = MAXROW
            NODEINFO(NEXTND)%KRHS   = LASTR
            NODEINFO(NEXTND)%KCOST  = LASTC
            NODEINFO(NEXTND)%KBOUND = LASTBD
            NODEINFO(NEXTND)%KNAMES = LASTNM
!
            GOTO 999
!
!     Come here if anything went wrong
!
 9150 CONTINUE
         WRITE (IOLOG, 1200) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3150)
         GOTO 9999
!
 9200 CONTINUE
         WRITE (IOLOG, 3200)
         GOTO 9999
!
 9300 CONTINUE
         WRITE (IOLOG, 3300)
         GOTO 9999
!
 9400 CONTINUE
         WRITE (IOLOG, 3400)
         GOTO 9999
!
 9980 CONTINUE
         WRITE (IOLOG, 3980)
         GOTO 9999
!
 9995 CONTINUE
         WRITE (IOLOG, 3995)
 9999 CONTINUE
         STOSTAT = 2
!
  999 CONTINUE
         CLOSE(IOSTO)
         IERR = STOSTAT
         RETURN
!
 1200 FORMAT(I8,4X,4A1,2X,A8,2X,A8,2X,F12.4,3X,A8,2X,F12.4)
 1510 FORMAT('(''Path'',I',i1,')')
 1520 FORMAT('(''Path'',I',i2,')')
 1600 FORMAT(' XXX -  FATAL  - Illegal record in STOCH file')
 1900 FORMAT(/,' Process STOCH file:')
 2000 FORMAT(' XXX - WARNING - Number of periods in CORE file does not',      &
     &       ' match information in TIME file')
 2100 FORMAT(' ')
 2600 FORMAT(' ***   Number of periods has been adjusted to',I3,'  ***')
 2700 FORMAT(I8,4X,4A1,A8,2X,A8)
 2800 FORMAT(' Solving problem ',A8,' -- ',I2,' periods')
 3150 FORMAT(' XXX -  FATAL  - Name does not match info in TIME file')
 3200 FORMAT(' XXX -  FATAL  - Memory allocation failed in routine',      &
     &       ' INSTOC.')
 3300 FORMAT(' XXX -  FATAL  - Global constraints are not allowed in',      &
     &       ' NODE format.')
 3400 FORMAT(' XXX -  FATAL  - Event nodes do not seem to form a tree.')
 3980 FORMAT(' XXX - WARNING - Error during READ or non-existent CORE',      &
     &       ' file')
 3990 FORMAT(' XXX - WARNING - Error during READ or non-existent STOCH',      &
     &       ' file')
 3995 FORMAT(' XXX -  FATAL  - Not enough information to solve the',      &
     &       ' problem')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INSIMP(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
!
!     Routine:   INSIMP
!
!     Purpose:
!
!         This subroutine reads the penalties (cost coefficients) for
!         problems with simple recourse.
!
!     Calling sequence:
!
!         CALL INSIMP(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
!
!            Q1     First character of the code field on an MPS-style data card
!            Q2     Second character of the code field
!            Q3     Third character of the code field
!            Q4     Fourth character of the code field
!            DNAME  Name fields on the data card (3 values)
!            ATEMP  Numeric fields on the data card (2 values)
!            NREC   Number of records processed
!            IERR   Error status
!
!     Routines called: GETMPS, IDROW, GENCOL
!
!     ------------------------------------------------------------------
!        This version dated 25 September 2002. Written by Gus Gassmann.
!     ------------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME(*), QLINE
      CHARACTER*8 DINT1, DINT2, DINT3, DINT4
!
      REAL*8      ATEMP(*)
      INTEGER*4   LENGTH(3)
!
      DATA DINT1 /'''INTORG'''/, DINT2  /'''SOSORG'''/,      &
     &     DINT3 /'''INTEND'''/, DINT4  /'''SOSEND'''/
!
      ALLOCATE (NDUP(NPER), NBLK(NPER), STAT=IERR)
      IF (IERR .GT. 0) GOTO 950
!
!     Read the next record
!
  890 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,      &
     &                     DNAME,ATEMP,IERR,IOSTO,NREC,LENGTH)
         IF (Q1     .EQ. QAST) GOTO 890
         IF (IERR   .GT.    0) GOTO 980
         IF (NECHO1 .GE.    2) WRITE (IOLOG, 1900) NREC, QLINE
         IF (Q1     .NE.  QBL) GOTO 900
!
!     Data card for simple recourse (ITYPE = 7).
!     The first name field should match the name of the simple recourse
!     section (if given in the specs file), the second name field
!     should match an existing row name. The two numeric fields must
!     satisfy num1 + num2 >= 0. Each occurrence of a row name in this section
!     triggers the generation of additional decision variables, identified by
!     appending the suffix `.minus' and `.plus' after the name of the row. If
!     the row represents an equality constraint, both variables are generated,
!     along with their coefficients in the A-matrix; if the row is of type `L',
!     only the shortage variable <rowname>.minus is created; if the row is of
!     type `G', only the surplus variable <rowname>.plus is added.
!
!     -----------------------------------------------------------------
!
!     First check if there is an integer marker in the row
!
  810 CONTINUE
         IF (    DNAME(3) .EQ. DINT1) THEN
            INTTYP = 1
            GOTO 890
         ELSEIF (DNAME(3) .EQ. DINT2) THEN
            WRITE (IOLOG, 1900) NREC, QLINE
            WRITE (IOLOG, 1600)
            GOTO 890
         ELSEIF (DNAME(3) .EQ. DINT3 .OR.      &
     &           DNAME(3) .EQ. DINT4) THEN
            INTTYP = 0
            GOTO 890
         ENDIF
!
!     Verify the name of the simple recourse section
!
         IF (DSIM .EQ. DOTS    ) DSIM = DNAME(1)
         IF (DSIM .NE. DNAME(1)) GOTO 890
!
!     Is the implied cost convex?
!
         IF (ATEMP(1) + ATEMP(2) .LT. 0.D0) THEN
            WRITE (IOLOG, 1900) NREC, QLINE
            WRITE (IOLOG, 3800)
            GOTO 890
         ENDIF
!
!     Identify the row
!
         DO IP=1,NPER
            IR = NMSRCH( DNAME(2),   NODEINFO(IP)%KNAMES+1,      &
     &         NODEINFO(IP)%KNAMES + NODEINFO(IP)%NROW )
            IF (IR .GT. 0) GOTO 814
         END DO
         WRITE (IOLOG, 3900) DNAME(2)
         GOTO 890
!
  814    CONTINUE
         FREEFORM = .TRUE.
         IR = IR - NODEINFO(IP)%KNAMES
         IV = IR + NODEINFO(IP)%VARPTR
         IF (VARTYP(IV) .EQ. 0) THEN
            IF (INTTYP .EQ. 0) THEN
               CALL GENCOL(DNAME(2),IP,IR,ATEMP,3,1)
            ELSE
               CALL GENCOL(DNAME(2),IP,IR,ATEMP,4,1)
            ENDIF
         ELSEIF (VARTYP(IV) .EQ. 1) THEN
             CALL GENCOL(DNAME(2),IP,IR,ATEMP,2,1)
         ELSEIF (VARTYP(IV) .EQ. -1) THEN
             CALL GENCOL(DNAME(2),IP,IR,ATEMP,1,1)
         ELSE
            WRITE (IOLOG, 3910) DNAME(2)
         ENDIF
         GOTO 890
!
!     Found the next header card
!
  900 CONTINUE
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 1900) NREC, QLINE
         GOTO 999
!
!     Could not allocate sufficient memory
!
  950 CONTINUE
      WRITE (IOLOG, 3990)
      STOSTAT = 2
      GOTO 999
!
!     End of file during read
!
  980 CONTINUE
         WRITE (IOLOG, 1900) NREC, QLINE
         WRITE (IOLOG, 1700)
  999 CONTINUE
         DEALLOCATE (NDUP, NBLK, STAT=IER)
         IF (IER .NE. 0) WRITE (IOLOG, 3980)
         RETURN
!
 1600 FORMAT(' XXX - WARNING - Integer directive ignored.')
 1700 FORMAT(' XXX - WARNING - Missing ENDATA card.')
 1900 FORMAT(I8,4X,A80)
 3400 FORMAT(' XXX - WARNING - Simple recourse and other penalty',      &
     &       ' functions must precede distributions.',/,      &
     &       ' Flush to the next header.')
 3800 FORMAT(' XXX - WARNING - Penalty parameters in simple recourse',      &
     &       ' section not convex. Data ignored.')
 3900 FORMAT(' XXX - WARNING - Row not found in SIMPLE section: ',A30)
 3910 FORMAT(' XXX - WARNING - Bad row type in SIMPLE section: ',A30)
 3980 FORMAT(' XXX - WARNING - Recoverable memory error in INSIMP.')
 3990 FORMAT(' XXX -  FATAL  - Insufficient memory in INSIMP.')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INRBST(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
!
!     Routine:   INRBST
!
!     Purpose:
!
!         This subroutine reads the penalties (entries in the Q-matrix)
!         for problems with quadratic penalties for constraint violations.
!         This also goes under the name "robust optimization".
!
!     Calling sequence:
!
!         CALL INRBST(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
!
!            Q1     First character of the code field on an MPS-style data card
!            Q2     Second character of the code field
!            Q3     Third character of the code field
!            Q4     Fourth character of the code field
!            DNAME  Name fields on the data card (3 values)
!            ATEMP  Numeric fields on the data card (2 values)
!            NREC   Number of records processed
!            IERR   Error status
!
!     Routines called: GETMPS, IDROW, GENCOL
!
!     ------------------------------------------------------------------
!        This version dated 25 September 2002. Written by Gus Gassmann.
!     ------------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME(*), QLINE
      CHARACTER*8 DINT1, DINT2, DINT3, DINT4
!
      REAL*8      ATEMP(*)
      INTEGER*4   LENGTH(3)
!
      DATA DINT1 /'''INTORG'''/, DINT2  /'''SOSORG'''/,      &
     &     DINT3 /'''INTEND'''/, DINT4  /'''SOSEND'''/
!
      ALLOCATE (NDUP(NPER), NBLK(NPER), STAT=IERR)
      IF (IERR .GT. 0) GOTO 950
!
      DEFPEN = ATEMP(1)
!
!     Read the next record
!
  890 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,      &
     &                     DNAME,ATEMP,IERR,IOSTO,NREC,LENGTH)
         IF (Q1     .EQ. QAST) GOTO 890
         IF (IERR   .GT.    0) GOTO 980
         IF (NECHO1 .GE.    2) WRITE (IOLOG, 1900) NREC, QLINE
         IF (Q1     .NE.  QBL) GOTO 900
!
!     Data card for robust optimization segment (ITYPE = 8).
!     The first name field should match the name of the robust optimization
!     section (if given in the specs file), the second name field should
!     match an existing row name. Each occurrence of a row name in this section
!     triggers the generation of additional decision variables, identified by
!     appending the suffix `.minus' and `.plus' after the name of the row. If
!     the row represents an equality constraint, both variables are generated,
!     along with their coefficients in the A-matrix; if the row is of type `L',
!     only the shortage variable <rowname>.minus is created; if the row is of
!     type `G', only the surplus variable <rowname>.plus is added.
!
!     -----------------------------------------------------------------
!
!     First check if there is an integer marker in the row
!
  810 CONTINUE
         IF (    DNAME(3) .EQ. DINT1) THEN
            INTTYP = 1
            GOTO 890
         ELSEIF (DNAME(3) .EQ. DINT2) THEN
            WRITE (IOLOG, 1900) NREC, QLINE
            WRITE (IOLOG, 1600)
            GOTO 890
         ELSEIF (DNAME(3) .EQ. DINT3 .OR.      &
     &           DNAME(3) .EQ. DINT4) THEN
            INTTYP = 0
            GOTO 890
         ENDIF
!
!     Verify the name of the robust optimization section
!
         IF (DROB .EQ. DOTS    ) DROB = DNAME(1)
         IF (DROB .NE. DNAME(1)) GOTO 890
!
!     Identify the row
!
         DO IP=1,NPER
            IR = NMSRCH( DNAME(2),   NODEINFO(IP)%KNAMES+1,      &
     &         NODEINFO(IP)%KNAMES + NODEINFO(IP)%NROW )
            IF (IR .GT. 0) GOTO 814
         END DO
         WRITE (IOLOG, 3900) DNAME(2)
         GOTO 890
!
  814    CONTINUE
         FREEFORM = .TRUE.
         IR = IR - NODEINFO(IP)%KNAMES
         IV = IR + NODEINFO(IP)%VARPTR
         IF (VARTYP(IV) .EQ. 0) THEN
            IF (ATEMP(1) .EQ. 0.D0) ATEMP(1) = DEFPEN
            IF (ATEMP(2) .EQ. 0.D0) ATEMP(2) = DEFPEN
            IF (INTTYP .EQ. 0) THEN
               CALL GENCOL(DNAME(2),IP,IR,ATEMP,3,2)
            ELSE
               CALL GENCOL(DNAME(2),IP,IR,ATEMP,4,2)
            ENDIF
         ELSEIF (VARTYP(IV) .EQ. 1) THEN
            IF (ATEMP(1) .EQ. 0.D0) ATEMP(1) = DEFPEN
            CALL GENCOL(DNAME(2),IP,IR,ATEMP,2,2)
         ELSEIF (VARTYP(IV) .EQ. -1) THEN
            IF (ATEMP(1) .EQ. 0.D0) ATEMP(1) = DEFPEN
            IF (ATEMP(2) .EQ. 0.D0) ATEMP(2) = ATEMP(1)
            CALL GENCOL(DNAME(2),IP,IR,ATEMP,1,2)
         ELSE
            WRITE (IOLOG, 3910) DNAME(2)
         ENDIF
         GOTO 890
!
!     Found the next header card
!
  900 CONTINUE
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 1900) NREC, QLINE
         GOTO 999
!
!     Could not allocate sufficient memory
!
  950 CONTINUE
      WRITE (IOLOG, 3990)
      STOSTAT = 2
      GOTO 999
!
!     End of file during read
!
  980 CONTINUE
         WRITE (IOLOG, 1900) NREC, QLINE
         WRITE (IOLOG, 1700)
  999 CONTINUE
         DEALLOCATE (NDUP, NBLK, STAT=IER)
         IF (IER .NE. 0) WRITE (IOLOG, 3980)
         RETURN
!
 1600 FORMAT(' XXX - WARNING - Integer directive ignored.')
 1700 FORMAT(' XXX - WARNING - Missing ENDATA card.')
 1900 FORMAT(I8,4X,A80)
 3400 FORMAT(' XXX - WARNING - Simple recourse and other penalty',      &
     &       ' functions must precede distributions.',/,      &
     &       ' Flush to the next header.')
 3800 FORMAT(' XXX - WARNING - Penalty parameters in simple recourse',      &
     &       ' section not convex. Data ignored.')
 3900 FORMAT(' XXX - WARNING - Row not found in SIMPLE section: ',A30)
 3910 FORMAT(' XXX - WARNING - Bad row type in SIMPLE section: ',A30)
 3980 FORMAT(' XXX - WARNING - Recoverable memory error in INRBST.')
 3990 FORMAT(' XXX -  FATAL  - Insufficient memory in INRBST.')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INPLQP(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
!
!     Routine:   INPLQP
!
!     Purpose:
!
!         This subroutine reads the cost and quadratic coefficients for
!         problems with piecewise lineaer quadratic penalties.
!
!     Calling sequence:
!
!         CALL INPLQP(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
!
!            Q1     First character of the code field on an MPS-style data card
!            Q2     Second character of the code field
!            Q3     Third character of the code field
!            Q4     Fourth character of the code field
!            DNAME  Name fields on the data card (3 values)
!            ATEMP  Numeric fields on the data card (2 values)
!            NREC   Number of records processed
!            IERR   Error status
!
!     Routines called: GETMPS, IDROW, GENCOL
!
!     ------------------------------------------------------------------
!        This version dated 25 September 2002. Written by Gus Gassmann.
!     ------------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME(*), QLINE
      CHARACTER*8 DINT1, DINT2, DINT3, DINT4
!
      REAL*8      ATEMP(*)
      INTEGER*4   LENGTH(3)
!
      DATA DINT1 /'''INTORG'''/, DINT2  /'''SOSORG'''/,      &
     &     DINT3 /'''INTEND'''/, DINT4  /'''SOSEND'''/
!
      ALLOCATE (NDUP(NPER), NBLK(NPER), STAT=IERR)
      IF (IERR .GT. 0) GOTO 950
!
!     Read the next record
!
  890 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,      &
     &                     DNAME,ATEMP,IERR,IOSTO,NREC,LENGTH)
         IF (Q1     .EQ. QAST) GOTO 890
         IF (IERR   .GT.    0) GOTO 980
         IF (NECHO1 .GE.    2) WRITE (IOLOG, 1900) NREC, QLINE
         IF (Q1     .NE.  QBL) GOTO 900
!
!     Data card for piecewise linear-quadratic penalty (ITYPE = 9).
!     If the row in name field 2 is an equality row, the card could have
!     a 'PL' or 'MI' code.
!
!     -----------------------------------------------------------------
!
!     First check if there is an integer marker in the row
!
  810 CONTINUE
         IF (    DNAME(3) .EQ. DINT1) THEN
            INTTYP = 1
            GOTO 890
         ELSEIF (DNAME(3) .EQ. DINT2) THEN
            WRITE (IOLOG, 1900) NREC, QLINE
            WRITE (IOLOG, 1600)
            GOTO 890
         ELSEIF (DNAME(3) .EQ. DINT3 .OR.      &
     &           DNAME(3) .EQ. DINT4) THEN
            INTTYP = 0
            GOTO 890
         ENDIF
!
!     Verify the name of the simple recourse section
!
         IF (DPLQ .EQ. DOTS    ) DPLQ = DNAME(1)
         IF (DPLQ .NE. DNAME(1)) GOTO 890
!
!     Identify the row
!
         DO IP=1,NPER
            IR = NMSRCH( DNAME(2),   NODEINFO(IP)%KNAMES+1,      &
     &         NODEINFO(IP)%KNAMES + NODEINFO(IP)%NROW )
            IF (IR .GT. 0) GOTO 814
         END DO
         WRITE (IOLOG, 3900) DNAME(2)
         GOTO 890
!
  814    CONTINUE
         FREEFORM = .TRUE.
         IR = IR - NODEINFO(IP)%KNAMES
         IV = IR + NODEINFO(IP)%VARPTR
         IF (VARTYP(IV) .EQ. 0) THEN
            IF (Q2 .EQ. QP .AND. Q3 .EQ. QL) THEN
               CALL GENCOL(DNAME(2),IP,IR,ATEMP,2,3)
            ELSEIF (Q2 .EQ. QM .AND. Q3 .EQ. QI) THEN
               CALL GENCOL(DNAME(2),IP,IR,ATEMP,1,3)
            ELSE
               IF (INTTYP .EQ. 0) THEN
                  CALL GENCOL(DNAME(2),IP,IR,ATEMP,3,3)
               ELSE
                  CALL GENCOL(DNAME(2),IP,IR,ATEMP,4,3)
               ENDIF
            ENDIF
         ELSEIF (VARTYP(IV) .EQ. 1) THEN
             CALL GENCOL(DNAME(2),IP,IR,ATEMP,2,3)
         ELSEIF (VARTYP(IV) .EQ. -1) THEN
             CALL GENCOL(DNAME(2),IP,IR,ATEMP,1,3)
         ELSE
            WRITE (IOLOG, 3950) DNAME(2)
         ENDIF
         GOTO 890
!
!     Found the next header card
!
  900 CONTINUE
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 1900) NREC, QLINE
         GOTO 999
!
!     Could not allocate sufficient memory
!
  950 CONTINUE
      WRITE (IOLOG, 3990)
      STOSTAT = 2
      GOTO 999
!
!     End of file during read
!
  980 CONTINUE
         WRITE (IOLOG, 1900) NREC, QLINE
         WRITE (IOLOG, 1700)
  999 CONTINUE
         DEALLOCATE (NDUP, NBLK, STAT=IER)
         IF (IER .NE. 0) WRITE (IOLOG, 3980)
         RETURN
!
 1600 FORMAT(' XXX - WARNING - Integer directive ignored.')
 1700 FORMAT(' XXX - WARNING - Missing ENDATA card.')
 1900 FORMAT(I8,4X,A80)
 3400 FORMAT(' XXX - WARNING - Simple recourse and other penalty',      &
     &       ' functions must precede distributions.',/,      &
     &       ' Flush to the next header.')
 3800 FORMAT(' XXX - WARNING - Penalty parameters in simple recourse',      &
     &       ' section not convex. Data ignored.')
 3900 FORMAT(' XXX - WARNING - Row not found in SIMPLE section: ',A30)
 3910 FORMAT(' XXX - WARNING - Bad row type in SIMPLE section: ',A30)
 3950 FORMAT(' XXX - WARNING - Bad row type in PLINQUAD section: ',A30)
 3980 FORMAT(' XXX -  FATAL  - Insufficient memory in INPLQP.')
 3990 FORMAT(' XXX - WARNING - Recoverable memory error in INPLQP.')
      END SUBROUTINE
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE GENCOL(ROWNAM,IP,ROWNBR,VALUE,MODE0,ITYPE)
!
!     This routine adds logical columns for simple recourse and other
!     forms of penalties.
!
!     Parameters:
!
!        ROWNAM - Name of the row mentioned in the stoch file
!        IP     - Period to which the row belongs
!        ROWNBR - Number of the row
!        VALUE  - Values of two coefficients
!        MODE0  - Determines which columns are to be generated:
!          1      <rowname>.minus (value is in VALUE(1))
!          2      <rowname>.plus  (value is in VALUE(2))
!          3      both columns generated (continuous variables)
!          4      both columns generated (integer variables)
!        ITYPE  - Determines the form of the penalty:
!          1      simple recourse (linear penalty)
!          2      robust optimization (quadratic penalty)
!          3      piecewise linear-quadratic penalty
!
!     Method:
!
!        The following steps are executed:
!            If ITYPE > 1
!                if there is no Q-matrix (NQPROF = 0)
!                    make a Q-matrix
!                    mark for insert: Q-matrix
!            For each column:
!                if defined before
!                    update cost (if ITYPE = 1 or ITYPE = 3)
!                    if IYPE = 2 or ITYPE = 3
!                       search Q-matrix for entry
!                       if found
!                          update Q-matrix
!                       else
!                          mark for insert: Q-matrix
!                else
!                    mark for insert: name, cost, bound, A-matrix
!                    if Q-matrix exists
!                        mark for insertion: Q-matrix columns
!            next column
!            insert marked names
!            insert marked bounds
!            insert marked cost (at zero level if ITYPE = 2)
!            if equality row and integer recourse
!                insert one row
!            insert marked A-matrix
!            if equality row and integer recourse
!                duplicate one constrainst
!            if Q-matrix exists
!                insert marked Q-matrix columns
!
!        The following arrays must be moved to make space for the insert:
!
!            COST,XLB,XUB,VARTYP,QNAMES,NAME1,LA,IA,A,LQMTX(,IQMTX,AQMTX).
!
!        The following arrays must be updated:
!
!            NCOL, NELMA, (NROW)
!
!        In addition, we must adjust these pointer arrays:
!
!            KCOL,KCOST,KNAMES,KBOUND,VARPTR,KCOLA,KELMA,LQOFF,IQOFF,LQOFF.
!
! ------------------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) ROWNAM,DBLANK,DTEMP,DROW
      CHARACTER*1        QROW(NAMLEN),QBLANK(NAMLEN),QTEMP(NAMLEN)
      REAL*8             VALUE(2)
      INTEGER*4          PROCESS(3),ROWNBR,INSERTC(2),INSERTQ(2)
!
      EQUIVALENCE (DROW,QROW),(DBLANK,QBLANK),(DTEMP,QTEMP)
!
      DATA QBLANK /NAMLEN*' '/
!
      IF (IERR .GT. 0) GOTO 770
      MODE = MODE0
      DROW = ROWNAM
      IR   = ROWNBR
      INSERTC(1) = 0
      INSERTC(2) = 0
      INSERTQ(1) = 0
      INSERTQ(2) = 0
      PROCESS(1) = 0
      PROCESS(2) = 0
      PROCESS(3) = 0
!
!     NSHIFT gives the number of characters by which the names have to be moved
!     ISHIFT gives the number of inserted columns
!     JSHIFT gives the number of inserted Q-matrix elements
!
      NSHIFT = 0
      ISHIFT = 0
      JSHIFT = 0
!
      KC = NODEINFO(IP)%KCOST
      KN = NODEINFO(IP)%KNAMES
      NR = NODEINFO(IP)%NROW
      NC = NODEINFO(IP)%NCOL
!
!     If the Q-matrix does not exist, make sure to initialize it
!
      IF (ITYPE .EQ. 2) THEN
         IF (.NOT. ALLOCATED(KDATQ)) THEN
            MNQBLK = MAX(IZQBLK,NPER)
            IF ((.NOT. RESIZE .AND. (MNQBLK .GT. IZQBLK))      &
     &                         .OR. (MNQBLK .GT. MXQBLK)) GOTO 760
            IZQBLK = MNQBLK
            MNQCOL = (MODE-1)/2
            DO I=1,NPER
               MNQCOL = MNQCOL + NODEINFO(I)%NCOL - NODEINFO(I)%NROW + 1
            END DO
            MNQCOL = MAX(IZQCOL,MNQCOL)
            IF ((.NOT. RESIZE .AND. (MNQCOL .GT. IZQCOL))      &
     &                         .OR. (MNQCOL .GT. MXQCOL)) GOTO 760
            IZQCOL = MNQCOL
            ALLOCATE(AQMTX(IZQLMN), IQMTX(IZQLMN), LQMTX(IZQCOL),      &
     &               KDATQ(IZNODE), IQOFF(IZQBLK), LQOFF(IZQBLK),      &
     &               NELMQ(IZQBLK), STAT=IER)
            IF (IER .NE. 0) GOTO 760
         ENDIF
         KQ = KDATQ(IP)
      ENDIF
!
      IF (ITYPE .EQ. 3) THEN
         IF (.NOT. ALLOCATED(CPLQ)) THEN
            IF (.NOT. RESIZE .AND. (IZCPLQ .LT. LASTC)) GOTO 760
            IF (LASTC .GT. MXCPLQ) GOTO 760
            IZCPLQ = MAX(IZCPLQ,LASTC)
            ALLOCATE(CPLQ(IZCPLQ), KPLQ(IZNODE), STAT=IER)
            IF (IER .NE. 0) GOTO 760
            HAVE_PLQP = .TRUE.
            CPLQ = 0.D0
            KPLQ(1:NPER) = NODEINFO(1:NPER)%KCOST
            LASTPLQ = LASTC
         ENDIF
      ENDIF
!
      IF (NQPROF .EQ. 0 .AND. ITYPE .EQ. 2) THEN
         NQPROF = 1
         DO I=1,NPER
            IQOFF(I) = LASTQ
            LQOFF(I) = LASTQC
            NELMQ(I) = 0
            KDATQ(I) = I - 1
            NCC = NODEINFO(I)%NCOL - NODEINFO(I)%NROW + 1
            DO J=1,NCC
               LQMTX(LASTQC+J) = 1
            END DO
            LASTQC = LASTQC + NCC
         END DO
         LASTQ = 0
      ENDIF
!
!     Determine the length of the row name.
!
      LEN1 = NAME1(KN+IR+1) - NAME1(KN+IR)
!
!     If the length is 8, there may be blanks at the end. If so, truncate.
!
      IF (LEN1 .EQ. 8) THEN
         DO LEN=LEN1,1,-1
            IF (QROW(LEN) .NE. QBL) GOTO 110
         END DO
      ENDIF
      LEN = LEN1
!
!     Check if the column names have been generated already.
!
  110 CONTINUE
      IF (MODE .GE. 2) THEN
         DTEMP = DBLANK
         DO I=1,LEN
            QTEMP(I) = QROW(I)
         END DO
         QTEMP(LEN+1) = '.'
         QTEMP(LEN+2) = 'm'
         QTEMP(LEN+3) = 'i'
         QTEMP(LEN+4) = 'n'
         QTEMP(LEN+5) = 'u'
         QTEMP(LEN+6) = 's'
         IR = NMSRCH(DTEMP,KN+NR,KN+NC)
         IF (IR .GT. 0) THEN
            WRITE (IOLOG, 1000) DTEMP
            PROCESS(1) = 2
            IF (ITYPE .EQ. 1) THEN
               COST(KC+IR) = VALUE(1)
            ELSEIF (ITYPE .EQ. 3) THEN
               COST(KC+IR) = DMAX1(VALUE(1),0.D0)
               CPLQ(KC+IR) = DMAX1(VALUE(2),0.D0)
            ELSE
               IQR = LQOFF(KQ+1) + IR - KN - NR
               DO J=LQMTX(IQR),LQMTX(IQR+1)-1
                  IF (IQMTX(IQOFF(IBLK)+J) .EQ. IR) THEN
                     AQMTX(IQOFF(IBLK)+J) = VALUE(1)
                     GOTO 140
                  ENDIF
               END DO
!
               INSERTC(1) = IR - KN - NR
               INSERTQ(1) = LQMTX(IQR+1) + JSHIFT
               JSHIFT = JSHIFT + 1
            ENDIF
         ELSE
            NSHIFT = NSHIFT + LEN + 6
            ISHIFT = ISHIFT + 1
            PROCESS(1) = 1
            IF (ITYPE .EQ. 2) THEN
               JSHIFT = JSHIFT + 1
               IBLK   = KQ + 1
               INSERTC(1) = NC - NR + JSHIFT
               INSERTQ(1) = LQMTX(LQOFF(IBLK)+NC-NR+1)      &
     &                         + JSHIFT - 1
            ENDIF
         ENDIF
      ENDIF
!
  140 CONTINUE
      IF (MODE .NE. 2) THEN
         DTEMP = DBLANK
         DO I=1,LEN
            QTEMP(I) = QROW(I)
         END DO
         QTEMP(LEN+1) = '.'
         QTEMP(LEN+2) = 'p'
         QTEMP(LEN+3) = 'l'
         QTEMP(LEN+4) = 'u'
         QTEMP(LEN+5) = 's'
         QTEMP(LEN+6) = ' '
         IR = NMSRCH(DTEMP,KN+NR,KN+NC)
         IF (IR .GT. 0) THEN
            WRITE (IOLOG, 1000) DTEMP
            PROCESS(2) = 2
            IF (ITYPE .EQ. 1) THEN
               IF (MODE .GE. 2 .OR. VALUE(2) .NE. 0.D0) THEN
                  COST(KC+IR) = VALUE(2)
               ELSE
                  COST(KC+IR) = VALUE(1)
               ENDIF
            ELSEIF (ITYPE .EQ. 3) THEN
               COST(KC+IR) = DMAX1(VALUE(1),0.D0)
               CPLQ(KC+IR) = DMAX1(VALUE(2),0.D0)
            ELSE
               IQR = LQOFF(KQ+1) + IR - KN - NR
               DO J=LQMTX(IQR),LQMTX(IQR+1)-1
                  IF (IQMTX(IQOFF(IBLK)+J) .EQ. IR) THEN
                     IF (MODE .GE. 2 .OR. VALUE(2) .NE. 0.D0) THEN
                        AQMTX(IQOFF(IBLK)+J) = VALUE(2)
                     ELSE
                        AQMTX(IQOFF(IBLK)+J) = VALUE(1)
                     ENDIF
                     GOTO 180
                  ENDIF
               END DO
!
               INSERTC(2) = IR - KN - NR
               INSERTQ(2) = LQMTX(IQR+1) + JSHIFT
               JSHIFT = JSHIFT + 1
            ENDIF
         ELSE
            NSHIFT = NSHIFT + LEN + 5
            ISHIFT = ISHIFT + 1
            PROCESS(2) = 1
            IF (ITYPE .EQ. 2) THEN
               JSHIFT = JSHIFT + 1
               IBLK   = KQ + 1
               INSERTC(2) = NC - NR + JSHIFT
               INSERTQ(2) = LQMTX(LQOFF(IBLK)+NC-NR+1)      &
     &                         + JSHIFT - 1
            ENDIF
!
!     If MODE = 4, then we have to duplicate a row (and create a new name)
!
            IF (MODE .EQ. 4) THEN
               QTEMP(LEN+1) = '_'
               QTEMP(LEN+2) = '2'
               QTEMP(LEN+3) = ' '
               QTEMP(LEN+4) = ' '
               QTEMP(LEN+5) = ' '
               MX = 0
!
  170          CONTINUE
                  IR = NMSRCH(DTEMP,KN+NR,KN+NC)
                  IF (IR .GT. 0) THEN
                     MX = MX + 1
                     QTEMP(LEN+MX+2) = 'x'
                     GOTO 170
                  ENDIF
!
                  PROCESS(3) = 1
                  NSHIFT = NSHIFT + LEN + 2 + MX
                  ISHIFT = ISHIFT + 1
            ENDIF
         ENDIF
      ENDIF
!
! --------------------------------------------------------------------------
!     If ISHIFT > 0, then we insert new columns. This will involve names,
!     cost vectors, bounds, and types, along with entries in the A-matrix.
! --------------------------------------------------------------------------
!
  180 CONTINUE
      IF (ISHIFT .GT. 0) THEN
!
!     Shift the names (in array QNAMES) to make room for the new column(s)
!
         MNCHAR = LASTCH + NSHIFT
         IF (MNCHAR .GT. IZCHAR) THEN
            IF(.NOT. EXTMEM_CHAR(MNCHAR,IERR)) GOTO 760
         ENDIF
         ICHAR0=NAME1(KN+NC+1) - 1
         IF (PROCESS(3) .EQ. 1) ICHAR1 = ICHAR0
         DO I=LASTCH,ICHAR0+1,-1
            QNAMES(I+NSHIFT) = QNAMES(I)
         END DO
         LASTCH = LASTCH + NSHIFT
!
!     Insert the names
!
         IF (PROCESS(3) .EQ. 1) ICHAR0 = ICHAR0 + LEN + 2 + MX
!
         IF (PROCESS(1) .EQ. 1) THEN
            DO I=1,LEN
               QNAMES(ICHAR0+I) = QTEMP(I)
            END DO
            ICHAR0 = ICHAR0 + LEN
            QNAMES(ICHAR0+1) = '.'
            QNAMES(ICHAR0+2) = 'm'
            QNAMES(ICHAR0+3) = 'i'
            QNAMES(ICHAR0+4) = 'n'
            QNAMES(ICHAR0+5) = 'u'
            QNAMES(ICHAR0+6) = 's'
            ICHAR0 = ICHAR0 + 6
         ENDIF
!
         IF (PROCESS(2) .EQ. 1) THEN
            DO I=1,LEN
               QNAMES(ICHAR0+I) = QTEMP(I)
            END DO
            ICHAR0 = ICHAR0 + LEN
            QNAMES(ICHAR0+1) = '.'
            QNAMES(ICHAR0+2) = 'p'
            QNAMES(ICHAR0+3) = 'l'
            QNAMES(ICHAR0+4) = 'u'
            QNAMES(ICHAR0+5) = 's'
            ICHAR0 = ICHAR0 + 5
         ENDIF

         IF (PROCESS(3) .EQ. 1) THEN
            ICHAR2 = NAME1(KN+NR+1) - 1
            DO I=ICHAR1,ICHAR2+1,-1
               QNAMES(I+LEN+2+MX) = QNAMES(I)
            END DO
            DO I=1,LEN
               QNAMES(ICHAR2+I) = QTEMP(I)
            END DO
            ICHAR2 = ICHAR2 + LEN
            QNAMES(ICHAR2+1) = '_'
            QNAMES(ICHAR2+2) = '2'
            ICHAR2 = ICHAR2 + 2
            DO I=1,MX
               QNAMES(ICHAR2+I) = 'x'
            END DO
         ENDIF
!
!     Adjust the name pointers (in array NAME1)
!
         MNVNAM = LASTNM + ISHIFT
         IF (MNVNAM .GT. IZVNAM) THEN
            IF (.NOT. EXTMEM_VNAM(MNVNAM,IERR)) GOTO 760
         ENDIF
!
         INAM0 = KN + NC + 1
         DO I=LASTNM,INAM0,-1
            NAME1(I+ISHIFT) = NAME1(I) + NSHIFT
         END DO
         IF (ISHIFT .EQ. 2) NAME1(INAM0+1) = NAME1(INAM0) + LEN + 6
         IF (ISHIFT .EQ. 3) THEN
            NAME1(INAM0+2) = NAME1(INAM0+3) - LEN - 5
            DO I=NC+1,NR+1,-1
               NAME1(KN+I+1) = NAME1(KN+I) + LEN + 2 + MX
            END DO
         ENDIF
         LASTNM = LASTNM + ISHIFT
!
!     Move and insert cost coefficients and piecewise linear quadratic penalty
!
         MNCOST = LASTC + ISHIFT
         IF (MNCOST .GT. IZCOST) THEN
            IF (.NOT. EXTMEM_COST(MNCOST,IERR)) GOTO 760
         ENDIF
!
         ICOST0 = KC + NC - NR + 1
         DO I=LASTC,ICOST0,-1
            COST(I+ISHIFT) = COST(I)
         END DO
         LASTC = LASTC + ISHIFT
!
         IF (ALLOCATED(CPLQ)) THEN
            IF (MNCOST .GT. IZCPLQ) THEN
               IF (.NOT. EXTMEM_CPLQ(MNCOST,IERR)) GOTO 760
            ENDIF
!
            DO I=LASTPLQ,ICOST0,-1
               CPLQ(I+ISHIFT) = CPLQ(I)
            END DO
            LASTPLQ = LASTPLQ + ISHIFT
         ENDIF
!
!
         IF (PROCESS(1) .EQ. 1) THEN
            IF (ITYPE .EQ. 1) THEN
               COST(ICOST0) = VALUE(1)
            ELSEIF (ITYPE .EQ. 2) THEN
               COST(ICOST0) = 0.D0
            ELSE
               COST(ICOST0) = DMAX1(VALUE(1),0.D0)
               CPLQ(ICOST0) = DMAX1(VALUE(2),0.D0)
            ENDIF
            ICOST0 = ICOST0 + 1
         ENDIF
         IF (PROCESS(2) .EQ. 1) THEN
            IF (ITYPE .EQ. 1) THEN
               IF (MODE .GE. 2 .OR. VALUE(2) .NE. 0.D0) THEN
                  COST(ICOST0) = VALUE(2)
               ELSE
                  COST(ICOST0) = VALUE(1)
               ENDIF
            ELSEIF (ITYPE .EQ. 2) THEN
               COST(ICOST0) = 0.D0
            ELSE
               COST(ICOST0) = DMAX1(VALUE(1),0.D0)
               CPLQ(ICOST0) = DMAX1(VALUE(2),0.D0)
            ENDIF
         ENDIF
!
!     If MODE = 4, then we must make room for one duplicated rhs
!
         IF (LASTR .GE. IZDRHS) THEN
            IF (.NOT. EXTMEM_DRHS(LASTR+1,IERR)) GOTO 760
         ENDIF
!
         IF (MODE .EQ. 4) THEN
            KR = NODEINFO(IP)%KRHS
            IRHS0 = KR + NR + 1
            DO I=LASTR,IRHS0,-1
               RHS(I+1) = RHS(I)
            END DO
            RHS(IRHS0) = RHS(KR+ROWNBR)
            LASTR = LASTR + 1
         ENDIF
!
!     Move and insert bounds
!
         MNBNDS = LASTBD + ISHIFT
         IF (MNBNDS .GT. IZBNDS) THEN
            IF (.NOT. EXTMEM_BNDS(MNBNDS,IERR)) GOTO 760
         ENDIF
!
         LSHIFT = MIN(ISHIFT,2)
         KB = NODEINFO(IP)%KBOUND
         IBND0 = KB + NC + 1
         DO I=LASTBD,IBND0,-1
            XLB(I+ISHIFT) = XLB(I)
            XUB(I+ISHIFT) = XUB(I)
         END DO
         LASTBD = LASTBD + ISHIFT
!
         DO I=ISHIFT-LSHIFT,ISHIFT-1
            XLB(IBND0+I) = 0.D0
            XUB(IBND0+I) = PLINF
         END DO
!
         IF (ISHIFT .GT. LSHIFT) THEN
            DO I=KB+NC,KB+NR+1,-1
               XLB(I+1) = XLB(I)
               XUB(I+1) = XUB(I)
            END DO
            XLB(KB+ROWNBR) = -PLINF
            XUB(KB+NR+1)   =  PLINF
         ENDIF
!
!     Move and insert variable type
!
         MNVTYP = LASTVP + ISHIFT
         IF (MNVTYP .GT. IZVTYP) THEN
            IF (.NOT. EXTMEM_VTYP(MNVTYP,IERR)) GOTO 760
         ENDIF
!
         KV = NODEINFO(IP)%VARPTR
         IVAR0 = KV + NC + 1
         DO I=LASTVP,IVAR0,-1
            VARTYP(I+ISHIFT) = VARTYP(I)
         END DO
         LASTVP = LASTVP + ISHIFT
!
         IF (MODE .EQ. 4) THEN
            IVT = 1
         ELSE
            IVT = 0
         ENDIF
         DO I=ISHIFT-LSHIFT,ISHIFT-1
            VARTYP(IVAR0+I) = IVT
         END DO
!
         IF (ISHIFT .GT. LSHIFT) THEN
            DO I=KV+NC,KV+NR+1,-1
               VARTYP(I+1) = VARTYP(I)
            END DO
            VARTYP(KV+ROWNBR) = -1
            VARTYP(KV+NR+1)   =  1
         ENDIF
!
!     Adjust variable pointers
!
         DO I=IP+1,NPER
            NODEINFO(I)%KCOL   = NODEINFO(I)%KCOL   + ISHIFT
            NODEINFO(I)%KCOST  = NODEINFO(I)%KCOST  + LSHIFT
            NODEINFO(I)%KNAMES = NODEINFO(I)%KNAMES + ISHIFT
            NODEINFO(I)%KBOUND = NODEINFO(I)%KBOUND + ISHIFT
            NODEINFO(I)%VARPTR = NODEINFO(I)%VARPTR + ISHIFT
         END DO
!
         IF (ALLOCATED(CPLQ)) THEN
            DO I=IP+1,NPER
               KPLQ(I) = NODEINFO(I)%KCOST
            END DO
         ENDIF
!
         IF (PROCESS(3) .EQ. 1) THEN
            DO I=IP+1,NPER
               NODEINFO(I)%KROW = NODEINFO(I)%KROW + 1
            END DO
         ENDIF
!
!     Determine the number of column pointers to add and see if this fits
!
         IF (MARKOV) THEN
            NBLKS = MIN(2,NPER-IP+1)
         ELSE
            NBLKS = NPER -IP + 1
         ENDIF
         IF (HAVE_GC) NBLKS = NBLKS + IP - 1
!
         MNACOL = LASTCA + LSHIFT*NBLKS
         IF (MNACOL .GT. IZACOL) THEN
            IF(.NOT. EXTMEM_ACOL(MNACOL,IERR)) GOTO 760
         ENDIF
!
!     Move and insert column pointers for A-matrix.
!
         IBLK = KDATA(IP) + 1
         ICA0 = KCOLA(IBLK) + NC - NR + 1
         DO I=LASTCA,ICA0+1,-1
            LA(I+LSHIFT) = LA(I)
         END DO
!
!     Move and insert row pointers and values.
!
!     The code is simpler if MODE < 4, where the only coefficients are +/-1
!
         IF (PROCESS(3) .EQ. 0) THEN
!
            LA(ICA0+1) = LA(ICA0) + 1
            IF (LSHIFT .EQ. 2) LA(ICA0+2) = LA(ICA0) + 2
            LASTCA = LASTCA + LSHIFT
!
            MNALMN = LASTA + LSHIFT
            IF (MNALMN .GT. IZALMN) THEN
               IF(.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 760
            ENDIF
!
            LMN0 = KELMA(IBLK) + NELMA(IBLK)
            DO I=LASTA,LMN0+1,-1
               IA(I+LSHIFT) = IA(I)
                A(I+LSHIFT) =  A(I)
            END DO
            LASTA = LASTA + LSHIFT
!
            IF (PROCESS(1) .EQ. 1) THEN
               LMN0 = LMN0 + 1
               IA(LMN0) = ROWNBR
                A(LMN0) = +1.D0
            ENDIF
            IF (PROCESS(2) .EQ. 1) THEN
               LMN0 = LMN0 + 1
               IA(LMN0) = ROWNBR
                A(LMN0) = -1.D0
            ENDIF
!
            LMN0 = LMN0 - LSHIFT
            DO I=1,LASTBA
               IF  (KCOLA(I) .GE. ICA0) KCOLA(I) = KCOLA(I) + LSHIFT
               IF ((KELMA(I) .GE. LMN0) .AND. (I .NE. IBLK))      &
     &              KELMA(I) = KELMA(I) + LSHIFT
            END DO
            NELMA(IBLK) = NELMA(IBLK) + LSHIFT
!
!     If MODE = 4, we split and copy an existing row, so we must be careful.
!     Start by counting the number of elements that must be duplicated;
!     starting with the subdiagonal matrices.
!
         ELSE
            IF (MARKOV .AND. IP .GT. 1) THEN
               NMTX = 2
            ELSE
               NMTX = IP
            ENDIF
!
            MDUPL = 0
            DO I=1,NMTX
               J = KDATA(IP) + I
               MDUP2 = 0
               DO K=KELMA(J)+1,KELMA(J)+NELMA(J)
                  IF (IA(K) .EQ. ROWNBR) MDUP2 = MDUP2 + 1
               END DO
               MDUPL = MDUPL + MDUP2
!
!     Sort the blocks intersecting the row into descending order of offset
!
               DO K=1,I-1
                  IF (KELMA(J) .GT. KELMA(NBLK(K))) GOTO 450
               END DO
               K = I
!
  450          CONTINUE
               DO L=I-1,K,-1
                  NBLK(L+1) = NBLK(L)
                  NDUP(L+1) = NDUP(L)
               END DO
               NBLK(K) = J
               NDUP(K) = MDUP2
            END DO
!
!     There might also be superdiagonal matrices to worry about.
!
            IF (HAVE_GC) THEN
               DO I=IP+1,NPER
                  IF (MARKOV) THEN
                     NSUBD = 2
                  ELSE
                     NSUBD = I
                  ENDIF
                  J = KDATA(I) + NSUBD + IP
!
                  MDUP2 = 0
                  DO K=KELMA(J)+1,KELMA(J)+NELMA(J)
                     IF (IA(K) .EQ. ROWNBR) MDUP2 = MDUP2 + 1
                  END DO
                  MDUPL = MDUPL + MDUP2
!
!     Sort the blocks intersecting the row into descending order of offset
!
                  DO K=1,I-1
                     IF (KELMA(J) .GT. KELMA(NBLK(K))) GOTO 480
                  END DO
                  K = I
!
  480             CONTINUE
                  DO L=I-1,K,-1
                     NBLK(L+1) = NBLK(L)
                     NDUP(L+1) = NDUP(L)
                  END DO
                  NBLK(K) = J
                  NDUP(K) = MDUP2
               END DO
            END IF
!
!     Make sure there is enough room to fit everything
!
            MSHIFT = MDUPL + 2
            MNALMN = LASTA + MSHIFT
            IF (MNALMN .GT. IZALMN) THEN
               IF(.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 760
            ENDIF
!
!     Now shift the matrix elements around. Blocks that don't intersect the
!     row are simply moved; other blocks are carefully moved and extended.
!
            LMN1  = LASTA
            LASTA = LASTA + MSHIFT
            DO II=1,NMTX
               JBLK = NBLK(II)
               LMNX = LMN1
               LMN1 = KELMA(JBLK) + NELMA(JBLK)
               DO I=LMNX,LMN1+1,-1
                  IA(I+MSHIFT) = IA(I)
                   A(I+MSHIFT) =  A(I)
               END DO
               DO I=1,LASTBA
                  IF (KELMA(I) .GT. KELMA(JBLK) .AND.      &
     &                KELMA(I) .LE. LMNX)  KELMA(I) = KELMA(I) + MSHIFT
               END DO
!
               NINS = NDUP(II)
               NELMA(JBLK) = NELMA(JBLK) + NINS
!
               IF (JBLK .EQ. IBLK) THEN
                  LA(ICA0+1) = LA(ICA0) + 1
                  LA(ICA0+2) = LA(ICA0) + 2
                  LASTCA = LASTCA + 2
                  MSHIFT = MSHIFT - 2
                  DO I=1,LASTBA
                     IF (KCOLA(I) .GE. ICA0) KCOLA(I) = KCOLA(I) + 2
                  END DO
!
                  IA(LMNX+MSHIFT+1) = ROWNBR
                  IA(LMNX+MSHIFT+2) = NR + 1
                   A(LMNX+MSHIFT+1) = +1.D0
                   A(LMNX+MSHIFT+2) = -1.D0
                  NELMA(IBLK) = NELMA(IBLK) + 2
               ENDIF
!
               IF (MSHIFT .GT. 0) THEN
                  LCC = NODEINFO(IP+I-II)%NCOL - NODEINFO(IP+I-II)%NROW
                  DO I=LCC,1,-1
                     L1 = KELMA(JBLK) + LA(KCOLA(JBLK)+I)
                     LZ = KELMA(JBLK) + LA(KCOLA(JBLK)+I+1) - 1
                     LA(KCOLA(JBLK)+I+1) = LA(KCOLA(JBLK)+I+1) + NINS
                     DO J=LZ,L1,-1
                        IF (IA(J) .EQ. ROWNBR) THEN
                           IA(J+MSHIFT) = NR + 1
                            A(J+MSHIFT) = A(J)
                           MSHIFT = MSHIFT - 1
                           NINS   = NINS   - 1
                        ENDIF
                        IA(J+MSHIFT) = IA(J)
                         A(J+MSHIFT) =  A(J)
                     END DO
                  END DO
               ENDIF
            END DO
         ENDIF
!
!     Insert empty columns in subdiagonal blocks
!
         IF (MARKOV .AND. IP .LT. NPER) THEN
            MAXP = IP + 1
         ELSE
            MAXP = NPER
         ENDIF
         DO IPC=IP+1,MAXP
            IBLK = KDATA(IPC)  + IPC - IP + 1
            ICA0 = KCOLA(IBLK) + NC  - NR + 1
            DO I=LASTCA,ICA0,-1
               LA(I+LSHIFT) = LA(I)
            END DO
            IF (ISHIFT .EQ. 2) LA(ICA0+1) = LA(ICA0)
            LASTCA = LASTCA + LSHIFT
!
            DO I=1,LASTBA
               IF (KCOLA(I) .GE. ICA0) KCOLA(I) = KCOLA(I) + ISHIFT
            END DO
         END DO
!
!     Insert empty columns in superdiagonal blocks
!
         IF (HAVE_GC) THEN
            IF (MARKOV) THEN
               IBLK0 = KDATA(IP) + 2
            ELSE
               IBLK0 = KDATA(IP) + IP
            ENDIF
            DO IPC=1,IP-1
               IBLK = IBLK0 + IPC
               ICA0 = KCOLA(IBLK) + NC - NR + 1
               DO I=LASTCA,ICA0,-1
                  LA(I+LSHIFT) = LA(I)
               END DO
               IF (ISHIFT .EQ. 2) LA(ICA0+1) = LA(ICA0)
               LASTCA = LASTCA + LSHIFT
!
               DO I=1,LASTBA
                  IF (KCOLA(I) .GE. ICA0) KCOLA(I) = KCOLA(I) + ISHIFT
               END DO
            END DO
         ENDIF
!
!     If Q-matrix exists, move and insert column pointers for Q-matrix.
!     We don't have to worry about splitting rows here (all the coefficients
!     refer to pairs of columns, so the code is simpler than for the A-matrix.)
!
         IF (NQPROF .GT. 0) THEN
            NBLKS = MIN(IP,NQPROF)
            MNQCOL = LASTQC + LSHIFT*NBLKS
            IF (MNQCOL .GT. IZQCOL) THEN
               IF(.NOT. EXTMEM_QCOL(MNQCOL,IERR)) GOTO 760
            ENDIF
!
            IBLK = KDATQ(IP) + 1
            ICQ0 = LQOFF(IBLK) + NC - NR + 1
            DO I=LASTQC,ICQ0,-1
               LQMTX(I+LSHIFT) = LQMTX(I)
            END DO
            IF (ISHIFT .EQ. 2) LQMTX(ICQ0+1) = LQMTX(ICQ0)
            LASTQC = LASTQC + LSHIFT
!
            DO I=1,LASTQB
               IF (LQOFF(I) .GE. ICQ0) LQOFF(I) = LQOFF(I) + LSHIFT
            END DO
!
!     Insert empty columns in subdiagonal blocks
!
            MAXP = MIN0(IP+NQPROF-1,NPER)
            DO IPC=IP+1,MAXP
               IBLK = KDATQ(IPC) + IPC - IP
               ICQ0 = LQOFF(IBLK) + NC - NR + 1
               DO I=LASTQC,ICQ0,-1
                  LQMTX(I+LSHIFT) = LQMTX(I)
               END DO
               IF (LSHIFT .EQ. 2) LQMTX(ICQ0+1) = LQMTX(ICQ0)
               LASTQC = LASTQC + LSHIFT
!
               DO I=1,LASTQB
                  IF (LQOFF(I) .GE. ICQ0) LQOFF(I) = LQOFF(I) + LSHIFT
               END DO
            END DO
         ENDIF
!
!     Here we insert the Q-matrix elements. First item of business:
!     to extend the columns containing the entries.
!
      IF (JSHIFT .GT. 0) THEN
         MNQLMN = LASTQ + JSHIFT
         IF (MNQLMN .GT. IZQLMN) THEN
            IF(.NOT. EXTMEM_QLMN(MNQLMN,IERR)) GOTO 760
         ENDIF
         IBLK = KQ + 1
         IF (INSERTC(1) .NE. 0) THEN
            ICQ0 = LQOFF(IBLK) + INSERTC(1)
         ELSE
            ICQ0 = LQOFF(IBLK) + INSERTC(2)
         ENDIF
!
         DO I=1,JSHIFT
            LQMTX(ICQ0+I) = LQMTX(ICQ0+I) + I
         END DO
!
         DO I=ICQ0+JSHIFT+1,LQOFF(IBLK)+NC-NR+1
            LQMTX(I) = LQMTX(I) + JSHIFT
         END DO
!
!     Move and insert row pointers and values
!
         LMN0 = IQOFF(IBLK) + LQMTX(ICQ0+1) - 1
         DO I=LASTQ,LMN0,-1
            IQMTX(I+JSHIFT) = IQMTX(I)
            AQMTX(I+JSHIFT) = AQMTX(I)
         END DO
         LASTQ = LASTQ + JSHIFT
!
         DO I=1,2
            IF (INSERTQ(I) .NE. 0) THEN
               IQMTX(IQOFF(IBLK)+INSERTQ(I)) = INSERTC(I)
               AQMTX(IQOFF(IBLK)+INSERTQ(I)) = VALUE(I)
            ENDIF
         END DO
         NELMQ(IBLK) = NELMQ(IBLK) + JSHIFT
!
         DO I=1,LASTQB
            IF ((IQOFF(I) .GE. LMN0-1) .AND. (I .NE. IBLK))      &
     &           IQOFF(I) = IQOFF(I) + JSHIFT
         END DO
      ENDIF
!
!     Finally adjust the number of columns (array NCOL)
!
         NODEINFO(IP)%NCOL = NODEINFO(IP)%NCOL + ISHIFT
         IF (MODE .EQ. 4) NODEINFO(IP)%NROW = NODEINFO(IP)%NROW + 1
      ENDIF
      GOTO 770
!
!     Error condition
!
  760 CONTINUE
      WRITE (IOLOG, 1100)
      IERR = 2
!
  770 CONTINUE
      RETURN
!
 1000 FORMAT(' XXX - WARNING - Column previously generated: ',A30)
 1100 FORMAT(' XXX -  FATAL  - Unable to allocate memory in GENCOL.')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INCHAN(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
!
!     Routine:   INCHAN
!
!     Purpose:
!
!     This routine processes probabilistic constraints, both individual
!     and joint. These data items are used to describe chance constraints:
!
!     LCHANCE - the number of chance constraints
!     NCHANCE - the number of chance constraints in each node
!     KCHANCE - pointer into XCHANCE for each node in the tree
!     XCHANCE - a pointer for each chance constraint that records the
!               type, probability, dimension, stage and offset.
!               CType:  like row type, +1 for >, -1 for <, etc.
!               Prob:   required probability level
!               Size:   number of rows participating in this constraint
!               Node:   node of the LATEST row participating in this constraint
!               ROWPTR: pointer into array RCHANCE containing the rows
!               NAMPTR: pointer into QCHANCE for identifying name of
!                       multidimensional chance constraints
!
!     RCHANCE - structure that for each chance constraint records the stage
!               and relative row number of every row participating in it.
!               Stage:  stage to which the row belongs
!               RelRow: row number (in relative address form)
!
!     Calling sequence:
!
!         CALL INCHAN(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
!
!            Q1     First character of the code field on an MPS-style data card
!            Q2     Second character of the code field
!            Q3     Third character of the code field
!            Q4     Fourth character of the code field
!            DNAME  Name fields on the data card (3 values)
!            ATEMP  Numeric fields on the data card (2 values)
!            NREC   Number of records processed
!            IERR   Error status
!
!     Routines called: GETMPS, IDROW, GENCOL
!
!     ------------------------------------------------------------------
!        This version dated 25 October 2002. Written by Gus Gassmann.
!     ------------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME(*), QLINE, DBLANK
      CHARACTER*1 QBLANK(NAMLEN)
!
      REAL*8        ATEMP(*)
      INTEGER*4     LENGTH(3)
      LOGICAL       KEEPCC
      TYPE (CHANCE) TCHANCE
!
      EQUIVALENCE (DBLANK,QBLANK)
!
      DATA QBLANK/NAMLEN*' '/
!
      IF (.NOT. ALLOCATED(NCHANCE)) THEN
         ALLOCATE(NCHANCE(IZNODE), XCHANCE(IZCCON), RCHANCE(IZCCRO),      &
     &            KCHANCE(IZNODE), QCHANCE(IZCCNM), STAT=IERR)
         IF (IERR .NE. 0) GOTO 950
      ENDIF
!
!     Read the next record
!
  890 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,      &
     &                     DNAME,ATEMP,IERR,IOSTO,NREC,LENGTH)
         IF (Q1     .EQ. QAST) GOTO 890
         IF (IERR   .GT.    0) GOTO 980
         IF (NECHO1 .GE.    2) WRITE (IOLOG, 1900) NREC, QLINE
         IF (Q1     .NE.  QBL) GOTO 900
!
!     Data card for chance-constraint segment (ITYPE = 6).
!     The card could have a 'G' or 'L' code (univariate) or a
!     'JG' or 'JL' code (multivariate), or it could contain just
!     a row name (in name field 1 or 2).
!
  800 CONTINUE
         IF (Q2 .EQ. QJ .AND. Q3 .EQ. QG) THEN
            IF (DCHANCE .EQ. DOTS) DCHANCE = DNAME(1)
            IF (DCHANCE .EQ. DNAME(1)) THEN
               LCHANCE = LCHANCE + 1
               IF (LCHANCE .GT. IZCCON) THEN
                  IF (.NOT. EXTMEM_CCON(LCHANCE,IERR)) GOTO 950
               ENDIF
               NJOINT = NJOINT + 1
               IF (NJOINT .GT. IZCCNM) THEN
                  IF (.NOT. EXTMEM_CCNM(NJOINT,IERR)) GOTO 950
               ENDIF
               XCHANCE(LCHANCE)%ROWPTR = NCELEM
               XCHANCE(LCHANCE)%NAMPTR = NJOINT
               XCHANCE(LCHANCE)%CTYPE  = -1
               XCHANCE(LCHANCE)%PROB   = ATEMP(1)
               XCHANCE(LCHANCE)%SIZE   = 0
               XCHANCE(LCHANCE)%NODE   = 0
               QCHANCE(NJOINT ) = DNAME(2)
               KEEPCC = .TRUE.
            ELSE
               KEEPCC = .FALSE.
            ENDIF
            GOTO 890
!
         ELSEIF (Q2 .EQ. QJ .AND. Q3 .EQ. QL) THEN
            IF (DCHANCE .EQ. DOTS) DCHANCE = DNAME(1)
            IF (DCHANCE .EQ. DNAME(1)) THEN
               LCHANCE = LCHANCE + 1
               IF (LCHANCE .GT. IZCCON) THEN
                  IF (.NOT. EXTMEM_CCON(LCHANCE,IERR)) GOTO 950
               ENDIF
               NJOINT = NJOINT + 1
               IF (NJOINT .GT. IZCCNM) THEN
                  IF (.NOT. EXTMEM_CCNM(NJOINT,IERR)) GOTO 950
               ENDIF
               XCHANCE(LCHANCE)%ROWPTR = NCELEM
               XCHANCE(LCHANCE)%NAMPTR = NJOINT
               XCHANCE(LCHANCE)%CTYPE  = 1
               XCHANCE(LCHANCE)%PROB   = ATEMP(1)
               XCHANCE(LCHANCE)%SIZE   = 0
               XCHANCE(LCHANCE)%NODE   = 0
               QCHANCE(NJOINT ) = DNAME(2)
               KEEPCC = .TRUE.
            ELSE
               KEEPCC = .FALSE.
            ENDIF
            GOTO 890
!
         ELSEIF (Q2 .EQ. QJ .AND. Q3 .EQ. QE) THEN
            IF (DCHANCE .EQ. DOTS) DCHANCE = DNAME(1)
            IF (DCHANCE .EQ. DNAME(1)) THEN
               LCHANCE = LCHANCE + 1
               IF (LCHANCE .GT. IZCCON) THEN
                  IF (.NOT. EXTMEM_CCON(LCHANCE,IERR)) GOTO 950
               ENDIF
               XCHANCE(LCHANCE)%ROWPTR = NCELEM
               XCHANCE(LCHANCE)%NAMPTR = NJOINT
               XCHANCE(LCHANCE)%CTYPE  = 0
               XCHANCE(LCHANCE)%PROB   = ATEMP(1)
               XCHANCE(LCHANCE)%SIZE   = 0
               XCHANCE(LCHANCE)%NODE   = 0
               KEEPCC = .TRUE.
            ELSE
               KEEPCC = .FALSE.
            ENDIF
            GOTO 890
!
         ELSEIF (Q2 .EQ. QG .OR. Q3 .EQ. QG) THEN
            IF (DCHANCE .EQ. DOTS) DCHANCE = DNAME(1)
            IF (DCHANCE .EQ. DNAME(1)) THEN
               CALL IDROW(DNAME(2),NR,IRP)
               IF (NR .GT. 0) THEN
                  LCHANCE = LCHANCE + 1
                  IF (LCHANCE .GT. IZCCON) THEN
                     IF (.NOT. EXTMEM_CCON(LCHANCE,IERR)) GOTO 950
                  ENDIF
                  XCHANCE(LCHANCE)%ROWPTR = NCELEM
                  XCHANCE(LCHANCE)%NAMPTR = 0
                  XCHANCE(LCHANCE)%CTYPE  = -1
                  XCHANCE(LCHANCE)%PROB   = ATEMP(1)
                  XCHANCE(LCHANCE)%SIZE   = 1
                  XCHANCE(LCHANCE)%NODE   = IRP
                  NCELEM  = NCELEM + 1
                  IF (NCELEM .GT. IZCCRO) THEN
                     IF (.NOT. EXTMEM_CCRO(NCELEM,IERR)) GOTO 950
                  ENDIF
                  RCHANCE(NCELEM)%STAGE  = IRP
                  RCHANCE(NCELEM)%RELROW = NR - NODEINFO(IRP)%KNAMES
               ELSE
                  WRITE (IOLOG, 3200)
               ENDIF
            ENDIF
            GOTO 890
!
         ELSEIF (Q2 .EQ. QL .OR. Q3 .EQ. QL) THEN
            IF (DCHANCE .EQ. DOTS) DCHANCE = DNAME(1)
            IF (DCHANCE .EQ. DNAME(1)) THEN
               CALL IDROW(DNAME(2),NR,IRP)
               IF (NR .GT. 0) THEN
                  LCHANCE = LCHANCE + 1
                  IF (LCHANCE .GT. IZCCON) THEN
                     IF (.NOT. EXTMEM_CCON(LCHANCE,IERR)) GOTO 950
                  ENDIF
                  XCHANCE(LCHANCE)%ROWPTR = NCELEM
                  XCHANCE(LCHANCE)%NAMPTR = 0
                  XCHANCE(LCHANCE)%CTYPE  = 1
                  XCHANCE(LCHANCE)%PROB   = ATEMP(1)
                  XCHANCE(LCHANCE)%SIZE   = 1
                  XCHANCE(LCHANCE)%NODE   = IRP
                  NCELEM  = NCELEM + 1
                  IF (NCELEM .GT. IZCCRO) THEN
                     IF (.NOT. EXTMEM_CCRO(NCELEM,IERR)) GOTO 950
                  ENDIF
                  RCHANCE(NCELEM)%STAGE  = IRP
                  RCHANCE(NCELEM)%RELROW = NR - NODEINFO(IRP)%KNAMES
               ELSE
                  WRITE (IOLOG, 3200)
               ENDIF
            ENDIF
            GOTO 890
!
         ELSEIF (Q2 .EQ. QE .OR. Q3 .EQ. QE) THEN
            IF (DCHANCE .EQ. DOTS) DCHANCE = DNAME(1)
            IF (DCHANCE .EQ. DNAME(1)) THEN
               CALL IDROW(DNAME(2),NR,IRP)
               IF (NR .GT. 0) THEN
                  LCHANCE = LCHANCE + 1
                  IF (LCHANCE .GT. IZCCON) THEN
                     IF (.NOT. EXTMEM_CCON(LCHANCE,IERR)) GOTO 950
                  ENDIF
                  XCHANCE(LCHANCE)%ROWPTR = NCELEM + 1
                  XCHANCE(LCHANCE)%NAMPTR = 0
                  XCHANCE(LCHANCE)%CTYPE  = 0
                  XCHANCE(LCHANCE)%PROB   = ATEMP(1)
                  XCHANCE(LCHANCE)%SIZE   = 1
                  XCHANCE(LCHANCE)%NODE   = IRP
                  NCELEM  = NCELEM + 1
                  IF (NCELEM .GT. IZCCRO) THEN
                     IF (.NOT. EXTMEM_CCRO(NCELEM,IERR)) GOTO 950
                  ENDIF
                  RCHANCE(NCELEM)%STAGE  = IRP
                  RCHANCE(NCELEM)%RELROW = NR - NODEINFO(IRP)%KNAMES
               ELSE
                  WRITE (IOLOG, 3200)
               ENDIF
            ENDIF
            GOTO 890
!
         ELSEIF (Q2 .EQ. QN .OR. Q3 .EQ. QN) THEN
            IF (DCHANCE .EQ. DOTS) DCHANCE = DNAME(1)
            IF (DCHANCE .EQ. DNAME(1)) THEN
               CALL IDROW(DNAME(2),NR,IRP)
               IF (NR .GT. 0) THEN
                  LCHANCE = LCHANCE + 1
                  IF (LCHANCE .GT. IZCCON) THEN
                     IF (.NOT. EXTMEM_CCON(LCHANCE,IERR)) GOTO 950
                  ENDIF
                  XCHANCE(LCHANCE)%ROWPTR = NCELEM + 1
                  XCHANCE(LCHANCE)%NAMPTR = 0
                  XCHANCE(LCHANCE)%CTYPE  = 2
                  XCHANCE(LCHANCE)%PROB   = ATEMP(1)
                  XCHANCE(LCHANCE)%SIZE   = 1
                  XCHANCE(LCHANCE)%NODE   = IRP
                  NCELEM  = NCELEM + 1
                  IF (NCELEM .GT. IZCCRO) THEN
                     IF (.NOT. EXTMEM_CCRO(NCELEM,IERR)) GOTO 950
                  ENDIF
                  RCHANCE(NCELEM)%STAGE  = IRP
                  RCHANCE(NCELEM)%RELROW = NR - NODEINFO(IRP)%KNAMES
               ELSE
                  WRITE (IOLOG, 3200)
               ENDIF
            ENDIF
            GOTO 890
!
         ELSEIF (Q2 .EQ. QBL .AND. Q3 .EQ. QBL) THEN
            IF (KEEPCC) THEN
               IF (DNAME(1) .EQ. DBLANK) DNAME(1) = DNAME(2)
               CALL IDROW(DNAME(1),NR,IRP)
               IF (NR .GT. 0) THEN
                  NCELEM = NCELEM + 1
                  IF (NCELEM .GT. IZCCRO) THEN
                     IF (.NOT. EXTMEM_CCRO(NCELEM,IERR)) GOTO 950
                  ENDIF
                  XCHANCE(LCHANCE)%SIZE   = XCHANCE(LCHANCE)%SIZE + 1
                  XCHANCE(LCHANCE)%NODE   =      &
     &                MAX(XCHANCE(LCHANCE)%NODE,IRP)
                  RCHANCE(NCELEM )%STAGE  = IRP
                  RCHANCE(NCELEM )%RELROW = NR - NODEINFO(IRP)%KNAMES
               ELSE
                  WRITE (IOLOG, 3200)
               ENDIF
            ENDIF
            GOTO 890
!
         ELSE
            WRITE (IOLOG, 3300)
            GOTO 850
         ENDIF
         GOTO 890
!
!     Flush to next header card
!
  850 CONTINUE
  860 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,      &
     &                     DNAME,ATEMP,IERR,IOSTO,NREC,LENGTH)
         IF (Q1 .EQ. QBL .OR. Q1 .EQ. QAST) GOTO 860
         IF (IERR .GT. 0) GOTO 980
!
!     Found the next header card. Sort the array XCHANCE into temporal order.
!     Also count the chance constraints for each node and initial pointers.
!
  900 CONTINUE
         LASTCC = 0
         DO I=1,NPER
            NCHANCE(I) = 0
            KCHANCE(I) = LASTCC
            DO J=LASTCC+1,LCHANCE
               IF (XCHANCE(J)%NODE .EQ. I) THEN
                  IF (J .GT. LASTCC+1) THEN
                     TCHANCE           = XCHANCE(LASTCC+1)
                     XCHANCE(LASTCC+1) = XCHANCE(J)
                     XCHANCE(J)        = TCHANCE
                  ENDIF
                  NCHANCE(I) = NCHANCE(I) + 1
                  LASTCC = LASTCC + 1
               ENDIF
            END DO
         END DO
!
         IF (LASTCC .LT. LCHANCE) THEN
            WRITE (IOLOG, 1300)
            LCHANCE = LASTCC
         ENDIF
!
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 1900) NREC, QLINE
         GOTO 999
!
!     Could not allocate sufficient memory
!
  950 CONTINUE
      WRITE (IOLOG, 3990)
      STOSTAT = 2
      GOTO 999
!
!     End of file during read
!
  980 CONTINUE
         WRITE (IOLOG, 1900) NREC, QLINE
         WRITE (IOLOG, 1700)
  999 CONTINUE
         RETURN
!
 1300 FORMAT(' XXX - WARNING - Found disassociated chance constraints',      &
     &       ' --- Removing.')
 1600 FORMAT(' XXX - WARNING - Integer directive ignored.')
 1700 FORMAT(' XXX - WARNING - Missing ENDATA card.')
 1900 FORMAT(I8,4X,A80)
 3200 FORMAT(' XXX - WARNING - Unidentified row name in',      &
     &       ' chance-constraint segment. Data ignored.')
 3300 FORMAT(' XXX - WARNING - Illegal code in chance-constraint',      &
     &       ' segment. Flush to the next header.')
 3400 FORMAT(' XXX - WARNING - Simple recourse and other penalty',      &
     &       ' functions must precede distributions.',/,      &
     &       ' Flush to the next header.')
 3800 FORMAT(' XXX - WARNING - Penalty parameters in simple recourse',      &
     &       ' section not convex. Data ignored.')
 3900 FORMAT(' XXX - WARNING - Row not found in SIMPLE section: ',A30)
 3910 FORMAT(' XXX - WARNING - Bad row type in SIMPLE section: ',A30)
 3950 FORMAT(' XXX - WARNING - Bad row type in PLINQUAD section: ',A30)
 3990 FORMAT(' XXX -  FATAL  - Insufficient memory in INCHAN.')
      END SUBROUTINE
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INICC(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
!
!     Routine:   INICC
!
!     Purpose:
!
!     This routine processes integrated chance constraints, sometimes
!     referred to as conditional value at risk (CVaR) constraints.
!
!     LICC - the number of integrated chance constraints
!     NICC - the number of integrated chance constraints in each node
!     KICC - pointer for each node to XICC
!     XICC - structure that for each integrated chance constraint records
!            the type and RHS, as well as the node and relative row number
!            of the participating row.
!            CType:  like row type, +1 for >, -1 for <, etc.
!            Limit:  RHS of the ICC
!            Node:   node to which the row belongs
!            RelRow: row number (in relative address form)
!
!     Calling sequence:
!
!         CALL INICC(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
!
!            Q1     First character of the code field on an MPS-style data card
!            Q2     Second character of the code field
!            Q3     Third character of the code field
!            Q4     Fourth character of the code field
!            DNAME  Name fields on the data card (3 values)
!            ATEMP  Numeric fields on the data card (2 values)
!            NREC   Number of records processed
!            IERR   Error status
!
!     Routines called: GETMPS, IDROW, GENCOL
!
!     ------------------------------------------------------------------
!        This version dated 25 October 2002. Written by Gus Gassmann.
!     ------------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME(*), QLINE, DBLANK
      CHARACTER*1 QBLANK(NAMLEN)
!
      REAL*8          ATEMP(*)
      INTEGER*4       LENGTH(3)
      TYPE (ICC_DATA) TICC
!
      EQUIVALENCE (DBLANK,QBLANK)
!
      DATA QBLANK/NAMLEN*' '/
!
      IF (.NOT. ALLOCATED(NICC)) THEN
         ALLOCATE(NICC(IZNODE), KICC(IZNODE), XICC(IZNICC), STAT=IERR)
         IF (IERR .GT. 0) GOTO 950
      ENDIF
!
!     Read the next record
!
  890 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,      &
     &                     DNAME,ATEMP,IERR,IOSTO,NREC,LENGTH)
         IF (Q1     .EQ. QAST) GOTO 890
         IF (IERR   .GT.    0) GOTO 980
         IF (NECHO1 .GE.    2) WRITE (IOLOG, 1900) NREC, QLINE
         IF (Q1     .NE.  QBL) GOTO 900
!
!     Data card for integrated chance-constraint segment.
!     The code field determines the type of the constraint.
!
  800 CONTINUE
         IF (Q2 .EQ. QG .OR. Q3 .EQ. QG) THEN
            IF (DICC .EQ. DOTS) DICC = DNAME(1)
            IF (DICC .EQ. DNAME(1)) THEN
               CALL IDROW(DNAME(2),NR,IRP)
               IF (NR .GT. 0) THEN
                  LICC = LICC + 1
                  IF (LICC .GT. IZNICC) THEN
                     IF (.NOT. EXTMEM_NICC(LICC,IERR)) GOTO 950
                  ENDIF
                  XICC(LICC)%CTYPE  = -1
                  XICC(LICC)%LIMIT  = ATEMP(1)
                  XICC(LICC)%NODE   = IRP
                  XICC(LICC)%RELROW = NR - NODEINFO(IRP)%KNAMES
               ELSE
                  WRITE (IOLOG, 3200)
               ENDIF
            ENDIF
            GOTO 890
!
         ELSEIF (Q2 .EQ. QL .OR. Q3 .EQ. QL) THEN
            IF (DICC .EQ. DOTS) DICC = DNAME(1)
            IF (DICC .EQ. DNAME(1)) THEN
               CALL IDROW(DNAME(2),NR,IRP)
               IF (NR .GT. 0) THEN
                  LICC = LICC + 1
                  IF (LICC .GT. IZNICC) THEN
                     IF (.NOT. EXTMEM_NICC(LICC,IERR)) GOTO 950
                  ENDIF
                  XICC(LICC)%CTYPE  = 1
                  XICC(LICC)%LIMIT  = ATEMP(1)
                  XICC(LICC)%NODE   = IRP
                  XICC(LICC)%RELROW = NR - NODEINFO(IRP)%KNAMES
               ELSE
                  WRITE (IOLOG, 3200)
               ENDIF
            ENDIF
            GOTO 890
!
         ELSEIF (Q2 .EQ. QE .OR. Q3 .EQ. QE) THEN
            IF (DICC .EQ. DOTS) DICC = DNAME(1)
            IF (DICC .EQ. DNAME(1)) THEN
               CALL IDROW(DNAME(2),NR,IRP)
               IF (NR .GT. 0) THEN
                  LICC = LICC + 1
                  IF (LICC .GT. IZNICC) THEN
                     IF (.NOT. EXTMEM_NICC(LICC,IERR)) GOTO 950
                  ENDIF
                  XICC(LICC)%CTYPE  = 0
                  XICC(LICC)%LIMIT  = ATEMP(1)
                  XICC(LICC)%NODE   = IRP
                  XICC(LICC)%RELROW = NR - NODEINFO(IRP)%KNAMES
               ELSE
                  WRITE (IOLOG, 3200)
               ENDIF
            ENDIF
            GOTO 890
!
         ELSEIF (Q2 .EQ. QN .OR. Q3 .EQ. QN) THEN
            IF (DICC .EQ. DOTS) DICC = DNAME(1)
            IF (DICC .EQ. DNAME(1)) THEN
               CALL IDROW(DNAME(2),NR,IRP)
               IF (NR .GT. 0) THEN
                  LICC = LICC + 1
                  IF (LICC .GT. IZNICC) THEN
                     IF (.NOT. EXTMEM_NICC(LICC,IERR)) GOTO 950
                  ENDIF
                  XICC(LICC)%CTYPE  = 2
                  XICC(LICC)%LIMIT  = ATEMP(1)
                  XICC(LICC)%NODE   = IRP
                  XICC(LICC)%RELROW = NR - NODEINFO(IRP)%KNAMES
               ELSE
                  WRITE (IOLOG, 3200)
               ENDIF
            ENDIF
            GOTO 890
!
         ELSE
            WRITE (IOLOG, 3300)
         ENDIF
         GOTO 890
!
!     Found the next header card. Sort the array XICC into temporal order.
!     Also count the ICCs for each node and initial pointers.
!
  900 CONTINUE
         LASTCC = 0
         DO I=1,NPER
            NICC(I) = 0
            KICC(I) = LASTCC
            DO J=LASTCC+1,LICC
               IF (XICC(J)%NODE .EQ. I) THEN
                  IF (J .GT. LASTCC+1) THEN
                     TICC           = XICC(LASTCC+1)
                     XICC(LASTCC+1) = XICC(J)
                     XICC(J)        = TICC
                  ENDIF
                  NICC(I) = NICC(I) + 1
                  LASTCC = LASTCC + 1
               ENDIF
            END DO
         END DO
!
         IF (LASTCC .LT. LICC) THEN
            WRITE (IOLOG, 1300)
            LICC = LASTCC
         ENDIF
!
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 1900) NREC, QLINE
         GOTO 999
!
!     Could not allocate sufficient memory
!
  950 CONTINUE
      WRITE (IOLOG, 3990)
      STOSTAT = 2
      GOTO 999
!
!     End of file during read
!
  980 CONTINUE
         WRITE (IOLOG, 1900) NREC, QLINE
         WRITE (IOLOG, 1700)
  999 CONTINUE
         RETURN
!
 1300 FORMAT(' XXX - WARNING - Found disassociated ICCs --- Removing.')
 1600 FORMAT(' XXX - WARNING - Integer directive ignored.')
 1700 FORMAT(' XXX - WARNING - Missing ENDATA card.')
 1900 FORMAT(I8,4X,A80)
 3200 FORMAT(' XXX - WARNING - Unidentified row name in',      &
     &       ' chance-constraint segment. Data ignored.')
 3300 FORMAT(' XXX - WARNING - Illegal code in ICC segment.',      &
     &       ' Record ignored.')
 3400 FORMAT(' XXX - WARNING - Simple recourse and other penalty',      &
     &       ' functions must precede distributions.',/,      &
     &       ' Flush to the next header.')
 3800 FORMAT(' XXX - WARNING - Penalty parameters in simple recourse',      &
     &       ' section not convex. Data ignored.')
 3900 FORMAT(' XXX - WARNING - Row not found in SIMPLE section: ',A30)
 3910 FORMAT(' XXX - WARNING - Bad row type in SIMPLE section: ',A30)
 3950 FORMAT(' XXX - WARNING - Bad row type in PLINQUAD section: ',A30)
 3990 FORMAT(' XXX -  FATAL  - Insufficient memory in INCHAN.')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INTREE ( IERR, NREC )
!
!--------------------------------------------------------------------------
!     Routine:   INTREE
!
!     Purpose:   This routine reads a stoch file in NODES format, where
!                every node in the decision tree is described explicitly.
!
!     Calling sequence:
!
!         CALL INTREE ( IERR, NREC )
!
!           IERR  Error status
!                   = 0  - no error
!                   > 0  - error condition
!
!           NREC  Number of records processed
!
!     Calls:  GETMPS, CPNODE, MKNODE
!
!     Method:
!
!         This program reads and processes stoch files in explicit NODES format.
!         The format is described in more detail in the paper by Gassmann and
!         Schweitzer. Each node in the event tree is explicitly mentioned and
!         the (conditional) branching probability is given along with the
!         ancestor node. Two options are available:
!
!             MK  to build an MPS file from scratch
!             CP  assemble an MPS file from prior info plus modifications
!
!     Remarks:  None
!
! ----------------------------------------------------------------------------
!               Version of December 17, 1997. Written by Gus Gassmann.
! ----------------------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME(3), DBLANK, DROW, QLINE
      CHARACTER*8 DROOT, DCORE
      LOGICAL     GOTONE, ERRCOR
      CHARACTER*1 QBLANK(NAMLEN)
      REAL*8      ATEMP(4)
      INTEGER*4   LENGTH(3)
!
      EQUIVALENCE (DBLANK,QBLANK)
!
      DATA QBLANK/NAMLEN*' '/, DROOT /'''ROOT'''/, DCORE /'''CORFIL'''/
!
!     Allocate temporary arrays
!
      ALLOCATE (DNODE(IZNODE), INUSE(IZTPER), KPATH(IZTPER),      &
     &                          KREF(IZTPER), STAT=IERR)
      IF (IERR .NE. 0) GOTO 9200
!
!     Initialization
!
      GOTONE   = .FALSE.
      ERRCOR   = (IERR .EQ. 1)
      IERR     = 0
      DROW     = DBLANK
      NODES    = NPER
      IPER0    = IZTPER
      NDINFO   = 0
      NPCORE   = NPER
      NAME1(1) = 1
      NODEINFO(1)%PROB = 1.0
!
!     Some of the nodes may have been read in the core file (if NPER > 0)
!
      IF (NPER .GT. 0) THEN
         DO I=1,NPER
            INUSE(I)  = 0
            IRNGE0(I) = 0
            DNODE(I)  = DBLANK
            NODEINFO(I)%IDESC = 0
         END DO
         LASTNM = LASTNM - 1
      ENDIF
!
!     Come here to process the next record
!
  800 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOSTO,NREC,LENGTH)
         IF (Q1   .EQ. QAST) GOTO 800
         IF (IERR .GT.    0) GOTO 910
         IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 900
         IF (Q2 .EQ. QC .AND. Q3 .EQ. QP) NDINFO = 1
         IF (Q2 .EQ. QM .AND. Q3 .EQ. QK) NDINFO = 2
         IF (NDINFO .EQ. 0) GOTO 9990
!
!     We have found the start of a new node
!
  802    CONTINUE
            IF (NECHO1 .GE. 1) WRITE (IOLOG, 1800) DNAME(1)
            DROW = DNAME(2)
            IF (DROW .EQ. DROOT .OR. DROW .EQ. DBLANK) THEN
!
!     This is the root node
!
            IF (GOTONE) GOTO 9805
               NPER = 1
               GOTONE = .TRUE.
               IP = 1
               J  = 0
               IF (DNAME(3) .EQ.  DCORE) GOTO 820
               IF (NODES    .GE. IZNODE) THEN
                  IF (.NOT. EXTMEM_NODE (NODES+1,IERR)) GOTO 9200
               ENDIF
!
                  IPER0 = 1
                  NODES = NODES + 1
                  DNODE(NODES) = DNAME(1)
                  IRNGE0(IP) = NODES
                  NODEINFO(NODES)%IANCTR = 0
                  IF (NODES .GT. 1) THEN
                     NODEINFO(NODES)%KCOL = NODEINFO(NODES-1)%NCOL      &
     &                                    + NODEINFO(NODES-1)%KCOL + 1
                     NODEINFO(NODES)%KROW = NODEINFO(NODES-1)%NROW      &
     &                                    + NODEINFO(NODES-1)%KROW
                  ENDIF
            ELSE
!
!     Not the root node; fix the period
!
               DO J=1,NODES
                  IF (DNODE(J) .EQ. DROW) GOTO 807
               END DO
               GOTO 9820

  807       CONTINUE
               IP = 1
               JAUX = J
               NCTR = J
  808       CONTINUE
               IP = IP + 1
               JAUX = NODEINFO(JAUX)%IANCTR
               IF (JAUX .GT. 0) GOTO 808
               IF (IP .GT. IZTPER) THEN
                  IF (.NOT. EXTMEM_TPER (IP,IERR)) GOTO 9200
               ENDIF
!
               IF (IP .GT. NPER ) NPER  = IP
               IF (IP .LT. IPER0) IPER0 = IP
               IF (DNAME(3) .EQ. DCORE .AND. INUSE(IP) .EQ. 0) GOTO 820
               IF (NODES    .GE. IZNODE) THEN
                  IF (.NOT. EXTMEM_NODE (NODES+1,IERR)) GOTO 9200
               ENDIF
!
               NODES = NODES + 1
               DNODE(NODES) = DNAME(1)
!
!     Place the node into the tree and adjust pointers
!
               IF (IRNGE0(IP) .EQ. 0) IRNGE0(IP) = NODES
               NODEINFO(NODES)%PROB   = ATEMP(1)
               NODEINFO(NODES)%IBROTH = 0
               NODEINFO(NODES)%IDESC  = 0
               NODEINFO(NODES)%NDESC  = 0
               NODEINFO(NODES)%IANCTR = NCTR
               IBRO1 = NODEINFO(NCTR)%IDESC
               NODEINFO(NCTR)%NDESC = NODEINFO(NCTR)%NDESC + 1
               IF (IBRO1 .EQ. 0) THEN
                  NODEINFO(NCTR)%IDESC = NODES
               ELSE
!
  811          CONTINUE
                  IF (NODEINFO(IBRO1)%IBROTH .NE. 0) THEN
                     IBRO1 = NODEINFO(IBRO1)%IBROTH
                     GOTO 811
                  ENDIF
!
                  NODEINFO(IBRO1)%IBROTH = NODES
               ENDIF
            ENDIF
!
!     Trace back ancestors to the root node
!
            KPATH(IP) = NODES
            NCTR = NODES
            JP = IP
  813       CONTINUE
               NCTR = NODEINFO(NCTR)%IANCTR
               JP = JP - 1
               IF (JP .GT. 0) THEN
                  KPATH(JP) = NCTR
                  GOTO 813
               ENDIF
!
!     Where does the information come from?
!
            IF (NDINFO .EQ. 1) THEN
!
!     Copy from the reference node (CP option)
!
               DROW = DNAME(3)
               DO J=1,NODES-1
                  IF (DROW .EQ. DNODE(J)) GOTO 817
               END DO
               IF (DROW .NE. DCORE) GOTO 9820
                  J = IP
!
  817       CONTINUE
               IREFND = J
               KREF(IP) = IREFND
               NCTR = IREFND
               JP = IP
  818       CONTINUE
               NCTR = NODEINFO(NCTR)%IANCTR
               JP = JP - 1
               IF (JP .GT. 0) THEN
                  KREF(JP) = NCTR
                  GOTO 818
               ENDIF
!
               ND = NODES
               IPER = IP
               CALL CPNODE ( ND, NDINFO, IERR, NREC, IREFND,      &
     &                       1, DNAME, ATEMP, Q1, Q2, Q3, Q4 )
               IF (IERR   .GT. 0) GOTO 9960
               IF (NDINFO .GT. 0) GOTO 802
            ELSE
!
!     Build the node (MK option)
!
               ND = NODES
               IPER = IP
               CALL MKNODE( ND, NDINFO, IERR, NREC, KPATH,      &
     &                      DNAME, ATEMP, Q1, Q2, Q3, Q4 )
               IF (IERR .GT. 0) GOTO 9960
               IF (NDINFO .GT. 0) GOTO 802
            ENDIF
            GOTO 900
!
!     Use information provided in the core file
!
  820    CONTINUE
            IF (ERRCOR) GOTO 9980
            IF (IP .GT. NPCORE) GOTO 9980
               DNODE(IP)  = DNAME(1)
               IRNGE0(IP) = IP
               INUSE(IP)  = 1
               NODEINFO(IP)%PROB = ATEMP(1)
               IF (J  .GT. 0) NODEINFO(IP)%IBROTH = NODEINFO(J)%IDESC
               IF (J  .GT. 0) NODEINFO(J )%IDESC  = IP
               IF (IP .GT. 1) NODEINFO(IP)%IANCTR = J
               IF (NODEINFO(IP)%IBROTH .EQ. IP) NODEINFO(IP)%IBROTH = 0
               IPER = IP
               IREFND = IP
               DO J=1,IP
                  KPATH(J) = J
                  KREF(J) = J
               END DO
               CALL CPNODE ( IP, NDINFO, IERR, NREC, IREFND,      &
     &                       3, DNAME, ATEMP, Q1, Q2, Q3, Q4 )
               IF (IERR   .GT. 0) GOTO 9960
               IF (NDINFO .GT. 0) GOTO 802
!
!     End of stoch file
!
  900 CONTINUE
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 1100)      &
     &      NREC,Q1,Q2,Q3,Q4,DNAME(1)
         GOTO 999
!
  910 CONTINUE
         WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 1200)
         GOTO 999
!
!     Come here if anything went wrong
!
 9200 CONTINUE
         WRITE (IOLOG, 3200)
         GOTO 9999
!
 9805 CONTINUE
         WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3805)
         GOTO 9999
!
 9820 CONTINUE
         WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3820)
         GOTO 9999
!
 9960 CONTINUE
         WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3960)
         GOTO 9999
!
 9980 CONTINUE
         WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3980)
         GOTO 9999
!
 9990 CONTINUE
         WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                 DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3990)
 9999 CONTINUE
         IERR = 2
  999 CONTINUE
         DEALLOCATE (DNODE, INUSE, KPATH, KREF, STAT=IER)
         IF (IER .NE. 0) WRITE (IOLOG, 3210)
         RETURN
!
 1100 FORMAT(I8,4X,4A1,A8,2X,A8,2X,F12.4,3X,A8,2X,F12.4)
 1200 FORMAT(' XXX - WARNING - Missing ENDATA card')
 1800 FORMAT(' Found node ',A8)
 3200 FORMAT(' XXX -  FATAL  - Unable to allocate memory in INTREE.')
 3210 FORMAT(' XXX - WARNING - Recoverable memory error in INTREE.')
 3805 FORMAT(' XXX -  FATAL  - There can be only ONE root node!')
 3820 FORMAT(' XXX -  FATAL  - Misspecified node name')
 3960 FORMAT(' XXX -  FATAL  - Error upon return from routine INNODE')
 3980 FORMAT(' XXX -  FATAL  - Error in CORE file detected while',      &
     &       ' processing STOCH FILE')
 3990 FORMAT(' XXX -  FATAL  - Error while reading STOCH FILE')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE MKNODE ( ND, NDINFO, IERR, NREC, KPREF,      &
     &                    DNAME, ATEMP, Q1, Q2, Q3, Q4 )
!
!     Routine:  MKNODE
!
!     Purpose:  This subroutine reads the data for one node in the modified
!               SMPS format described in the paper by Gassmann and Schweitzer.
!
!     Calling Sequence:
!
!           CALL MKNODE ( ND, NDINFO, IERR, NREC, KPREF,
!     *                   DNAME, ATEMP, Q1, Q2, Q3, Q4 )
!
!     Parameters:
!
!         ND      Number of the node to be processed
!         NDINFO  Format of the *next* node to be processed
!                 (has to be tracked because the end of a node may be
!                 defined implicitly by encountering the next node header)
!           = 0   END card has been detected; no format implied
!           = 1   CP header has been detected: copy the next node
!           = 2   MK header has been detected: build the next node
!         IERR    Error status
!           = 0   Normal termination
!           > 0   Error condition
!         NREC    Number of records read
!         KPREF   Array of reference nodes, contains all the ancestor nodes
!         DNAME   Name fields on the MPS-style data card
!         ATEMP   Numerical fields on the MPS-style data card
!         Q1      First character of the code field
!         Q2      Second character of the code field
!         Q3      Third character of the code field
!         Q4      Fourth character of the code field
!
!
!     Calls: GETMPS, GETNET
!
!     Method:
!
!     This subroutine reads the data for one node in the modified MPS
!     format described in the standards paper. There are a great many
!     possible sections, and the ENDATA card is ambiguous, as it can
!     denote the end of either the LP section or the QSECTION or the
!     end of the node information altogether.
!
!     The different header cards that may be encountered can be read off from
!     the table below:
!
!      ----------------------------------------------------------------------
!     | possible | before  | during  | during  | during  | during  | after   |
!     | headers  | ROWS    | ROWS    | COLS    | RHS etc.| QSECT   | QSECT   |
!     |          | section | section | section | section | section | section |
!      ----------------------------------------------------------------------
!     |          | ignore  |         | finish  | wait    |         |         |
!     | NAME     |(no name | illegal | COLS    | for     | illegal | illegal |
!     |          |   check)|         | start QS| QSECTION|         |         |
!      ----------------------------------------------------------------------
!     |          | start   |         |         |         |         |         |
!     | ROWS     | ROWS    | ignore  | illegal | illegal | illegal | illegal |
!     |          | section |         |         |         |         |         |
!      ----------------------------------------------------------------------
!     |          | start   |         |         |         |         |         |
!     | NODE     | ROWS    | ignore  | illegal | illegal | illegal | illegal |
!     |          | section |         |         |         |         |         |
!      ----------------------------------------------------------------------
!     |          |         | start   | switch  |         |         |         |
!     | COLUMNS  | illegal | COLS    | to COLS | illegal | illegal | illegal |
!     |          |         | section |         |         |         |         |
!      ----------------------------------------------------------------------
!     |          |         | ignore  | switch  |         |         |         |
!     | ARCS     | illegal | data    | to ARCS | illegal | illegal | illegal |
!     |          |         | cards   |         |         |         |         |
!      ----------------------------------------------------------------------
!     |          |         |         | finish  | process |         |         |
!     | RHS      | illegal | illegal | COLS    | and     | illegal | illegal |
!     |          |         |         |         | continue|         |         |
!      ----------------------------------------------------------------------
!     | SUPPLY   |         |         | finish  | process |         |         |
!     |   OR     | illegal | illegal | COLS    | and     | illegal | illegal |
!     | DEMAND   |         |         |         | continue|         |         |
!      ----------------------------------------------------------------------
!     | RANGES   |         |         | finish  |         |         |         |
!     |   OR     | illegal | illegal | COLS    | continue| illegal | illegal |
!     | BOUNDS   |         |         |         |         |         |         |
!      ----------------------------------------------------------------------
!     |          |         |         | finish  |         | finish  |         |
!     | ENDATA   | illegal | illegal | COLS    | continue| QSECT   | continue|
!     |          |         |         |         |         |         |         |
!      ----------------------------------------------------------------------
!     |          |         |         | finish  | start   |         |         |
!     | QSECTION | illegal | illegal | COLS    | QSECT   | illegal | illegal |
!     |          |         |         | start QS|         |         |         |
!      ----------------------------------------------------------------------
!     |          | print   |         | finish  |         | finish  |         |
!     | 'CP'/'MK'| warning;| illegal | COLS    | return  | QSECT   | return  |
!     |          | return  |         | return  |         | return  |         |
!      ----------------------------------------------------------------------
!     |          | print   |         | finish  |         | finish  |         |
!     | <EOF>    | warning;| illegal | COLS    | return  | QSECT   | return  |
!     |          | return  |         | return  |         | return  |         |
!      ----------------------------------------------------------------------
!
!     Remarks: None.
!
!     --------------------------------------------------------------
!     This version dated December 15, 1997. Written by Gus Gassmann.
!     --------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME(*), QLINE, DBLANK, DROW, DCOL, DCOL1,      &
     &            DCOL2
      CHARACTER*8 DINT1, DINT2,  DINT3, DINT4
      CHARACTER*1 DB255(NAMLEN), QROW(NAMLEN), QCOL(NAMLEN)
      INTEGER*4   KPREF(*)
      INTEGER*4   LENGTH(3),S3ROW,FIXCOL
      LOGICAL     MPSFORM,SUPPLY
      REAL*8      ATEMP(*)
      INTEGER, ALLOCATABLE, DIMENSION (:,:) :: MREF
!
      EQUIVALENCE (DBLANK,DB255), (DROW,  QROW), (DCOL,  QCOL)
!
      DATA DINT1 /'''INTORG'''/, DINT2  /'''SOSORG'''/,      &
     &     DINT3 /'''INTEND'''/, DINT4  /'''SOSEND'''/,      &
     &     DB255 /NAMLEN*' '/
!
!     Allocate temporary arrays
!
      ALLOCATE (ISTART(IZTPER,2),      &
     &   QLINKS(IZQLMN), QRTMP(IZQLMN), QCTMP(IZQLMN), QMTMP(IZQLMN),      &
     &   ALINKS(IZANZB), ARTMP(IZANZB), ACTMP(IZANZB), AMTMP(IZANZB),      &
     &   STAT=IERR)
      IF (IERR .GT. 0) GOTO 9200
!
!     Initialize pointers
!
      IF (ND .GT. 1) THEN
         NODEINFO(ND)%KROW   = NODEINFO(ND-1)%KROW + NODEINFO(ND-1)%NROW
         NODEINFO(ND)%KCOL   = NODEINFO(ND-1)%KCOL      &
     &                       + NODEINFO(ND-1)%NCOL + 1
         NODEINFO(ND)%KBOUND = LASTBD
         NODEINFO(ND)%NROW   = 0
         NODEINFO(ND)%NCOL   = 0
         NODEINFO(ND)%KRHS   = LASTR
         NODEINFO(ND)%KCOST  = LASTC
         NODEINFO(ND)%KNAMES = LASTNM
         NODEINFO(ND)%VARPTR = LASTVP
         KDATA(ND) = LASTBA
         IF (ALLOCATED(KDATQ)) KDATQ(ND) = LASTQB
      ENDIF
!
      DO I=1,IPER
         ISTART(I,1) = 0
      END DO
!
!     Other initializations
!
         NROWS  = 0
         IENDAT = 0
         IROW   = NODEINFO(ND)%KROW
         ICOL   = NODEINFO(ND)%KCOL
         ICOST  = NODEINFO(ND)%KCOST
         IBOUND = NODEINFO(ND)%KBOUND
         INAMES = NODEINFO(ND)%KNAMES
         IVTYP  = NODEINFO(ND)%VARPTR
         NAME1(INAMES+1) = LASTCH + 1
         LASTBA = LASTBA + 1
         IF (LASTBA .GT. IZABLK) THEN
            IF (.NOT. EXTMEM_ABLK(LASTBA,IERR)) GOTO 9200
         ENDIF
!
         KCOLA(LASTBA) = LASTCA
         KELMA(LASTBA) = LASTA
         IELMA = LASTA
         ICOLA = LASTCA
         DCOL  = DOTS
!
!     Now read the node information. Start with the ROWS section
!
  150 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOSTO,NREC,LENGTH)
         IF (Q1 .EQ. QAST) GOTO 150
         IF (IERR .GT.  0) GOTO 154
         IF (Q2 .EQ. QC .AND. Q3 .EQ. QP) GOTO 153
         IF (Q2 .EQ. QM .AND. Q3 .EQ. QK) GOTO 153
         IF (Q1 .EQ. QBL ) GOTO 160
         IF (Q1 .EQ. QN .AND. Q2 .EQ. QA) GOTO 155
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QO) GOTO 151
         IF (Q1 .EQ. QN .AND. Q2 .EQ. QO) GOTO 151
         IF (Q1 .EQ. QA .AND. Q2 .EQ. QR) GOTO 200
         IF (Q1 .EQ. QC .AND. Q2 .EQ. QO) GOTO 200
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QH) GOTO 9300
         IF (Q1 .EQ. QS .AND. Q2 .EQ. QU) GOTO 9300
         IF (Q1 .EQ. QD .AND. Q2 .EQ. QE) GOTO 9300
         IF (Q1 .EQ. QB .AND. Q2 .EQ. QO) GOTO 9300
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QA) GOTO 9300
         IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 156
         IF (Q1 .EQ. QQ .AND. Q2 .EQ. QS) GOTO 157
!
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 2300)
            GOTO 9999
!
!     ROW or NODE card. Store the objective as the first row
!
  151    CONTINUE
            IF (NECHO1 .GE. 1) WRITE (IOLOG, 1400) NREC,QLINE
            IF (IENDAT .EQ. 0) THEN
               IF (LASTCH + LENOBJ .GE. IZCHAR) THEN
                  IF (.NOT. EXTMEM_CHAR(LASTCH+LENOBJ,IERR)) GOTO 9200
               ENDIF
               NROWS  = 1
               MAXROW = MAXROW + 1
               NODEINFO(ND)%NROW = 1
               NODEINFO(ND)%NCOL = 1
               VARTYP(IVTYP+1)   = 2
               XUB(LASTBD+1)     = PLINF
               XLB(LASTBD+1)     =-PLINF
               DO J=1,LENOBJ
                  QNAMES(LASTCH+J) = QNAMES(J)
               END DO
               LASTCH = LASTCH + LENOBJ
               MNVNAM = INAMES + 2
               IF (MNVNAM .GT. IZVNAM) THEN
                  IF (.NOT. EXTMEM_VNAM(MNVNAM,IERR)) GOTO 9200
               ENDIF
               NAME1(INAMES+2) = LASTCH + 1
               IENDAT = 1
            ENDIF
            GOTO 150
!
!     Next node detected ('CP' or 'MK' code)
!
  153 CONTINUE
         IF (IENDAT .EQ. 1) GOTO 9300
            WRITE (IOLOG, 2600)
            GOTO 900
!
!     End of file or error during read
!
  154 CONTINUE
         IF (IENDAT .EQ. 1) GOTO 9300
         IF (IERR   .GT. 1) GOTO 9400
            WRITE (IOLOG, 2600)
            Q1 = QE
            Q2 = QN
            Q3 = QD
            GOTO 900
!
!     Found a NAME card
!
  155 CONTINUE
         IF (IENDAT .EQ. 1) GOTO 9300
            GOTO 150
!
!     Found an ENDATA card
!
  156 CONTINUE
         IF (IENDAT .EQ. 1) GOTO 9300
            WRITE (IOLOG, 2600)
            IENDAT = 4
            GOTO 450
!
!     Found a QSECTION card
!
  157 CONTINUE
         IF (IENDAT .EQ. 1) GOTO 9300
            IF (IENDAT .NE. 4) WRITE (IOLOG, 2600)
            IENDAT = 5
            GOTO 449
!
!     Here we have a regular data card. Store the row info.
!
  160    CONTINUE
            IF (NECHO1 .GE. 2) WRITE (IOLOG, 1400)      &
     &         NREC,QLINE
            IF (NROWS  .GE. IZROWP) THEN
               IF (.NOT. RESIZE .OR. (IZROWP .GE. MXROWP)      &
     &                          .OR. (NROWS  .GE. MXROWP)) GOTO 9200
               IZROWP = NROWS + 1
            ENDIF
            IF (MAXROW .GE. IZROWS) THEN
               IF (.NOT. RESIZE .OR. (IZROWS .GE. MXROWS)      &
     &                          .OR. (MAXROW .GE. MXROWS)) GOTO 9200
               IZROWS = MAXROW + 1
            ENDIF
            IF (IVTYP + NROWS .GE. IZVTYP) THEN
               IF (.NOT. EXTMEM_VTYP(IVTYP+NROWS+1,IERR)) GOTO 9200
            ENDIF
            IF (LASTNM + NROWS + 1.GE. IZVNAM) THEN
               IF (.NOT. EXTMEM_VNAM(LASTNM+NROWS+2,IERR)) GOTO 9200
            ENDIF
            IF (LASTBD + NROWS .GE. IZBNDS) THEN
               IF (.NOT. EXTMEM_BNDS(LASTBD+NROWS+1,IERR)) GOTO 9200
            ENDIF
            IF (LASTCH + LENGTH(1) .GT. IZCHAR) THEN
               IF (.NOT. EXTMEM_CHAR(LASTCH+LENGTH(1),IERR)) GOTO 9200
            ENDIF
!
            IF ((Q2 .EQ. QE) .OR. (Q3 .EQ. QE) .OR.      &
     &         ((Q2 .EQ. QBL).AND.(Q3 .EQ. QBL))) THEN
               MAXROW              = MAXROW + 1
               NROWS               = NROWS + 1
               NODEINFO(ND)%NROW   = NROWS
               NODEINFO(ND)%NCOL   = NROWS
               XUB(LASTBD+NROWS)   = 0.D0
               XLB(LASTBD+NROWS)   = 0.D0
               VARTYP(IVTYP+NROWS) = 0
               DROW = DNAME(1)
               DO I=1,LENGTH(1)
                  QNAMES(LASTCH+I) = QROW(I)
               END DO
               LASTCH = LASTCH + LENGTH(1)
               NAME1(INAMES+NROWS+1) = LASTCH + 1
!
           ELSEIF  ((Q2 .EQ. QG) .OR. (Q3 .EQ. QG)) THEN
               MAXROW              = MAXROW + 1
               NROWS               = NROWS + 1
               NODEINFO(ND)%NROW   = NROWS
               NODEINFO(ND)%NCOL   = NROWS
               XUB(LASTBD+NROWS)   = 0.D0
               XLB(LASTBD+NROWS)   = -PLINF
               VARTYP(IVTYP+NROWS) = -1
               DROW = DNAME(1)
               DO I=1,LENGTH(1)
                  QNAMES(LASTCH+I) = QROW(I)
               END DO
               LASTCH = LASTCH + LENGTH(1)
               NAME1(INAMES+NROWS+1) = LASTCH + 1
!
            ELSEIF  ((Q2 .EQ. QL) .OR. (Q3 .EQ. QL)) THEN
               MAXROW              = MAXROW + 1
               NROWS               = NROWS + 1
               NODEINFO(ND)%NROW   = NROWS
               NODEINFO(ND)%NCOL   = NROWS
               XUB(LASTBD+NROWS)   = PLINF
               XLB(LASTBD+NROWS)   = 0.D0
               VARTYP(IVTYP+NROWS) = 1
               DROW = DNAME(1)
               DO I=1,LENGTH(1)
                  QNAMES(LASTCH+I) = QROW(I)
               END DO
               LASTCH = LASTCH + LENGTH(1)
               NAME1(INAMES+NROWS+1) = LASTCH + 1
!
            ELSEIF  ((Q2 .EQ. QN) .OR. (Q3 .EQ. QN)) THEN
!
!     Unbounded row could be the objective. If not, treat like any other row.
!
               IF (OBJNAM .EQ. DOTS) THEN
                  DCOL = DNAME(1)
                  IF (LENGTH(1) .LE. 8) THEN
                     DO I=1,NODES
                        N1 = NAME1(NODEINFO(I)%KNAMES+1) - 1
                        DO J=1,8
                           QNAMES(N1+J) = QCOL(J)
                        END DO
                     END DO
                  ELSE
                     ISHIFT = LENGTH(1) - 8
                     MNCHAR = LASTCH + ISHIFT
                     IF (MNCHAR .GT. IZCHAR) THEN
                        IF (.NOT. EXTMEM_CHAR(MNCHAR,IERR)) GOTO 9200
                     ENDIF
                     DO I=1,NODES
                        INAMES = NODEINFO(I)%KNAMES
                        DO J=LASTCH,NAME1(INAMES+2),-1
                           QNAMES(J+ISHIFT) = QNAMES(J)
                        END DO
                        DO J=INAMES+2,INAMES+NROWS+1
                           NAME1(J) = NAME1(J) + ISHIFT
                        END DO
                        N1 = NAME1(INAMES+1) - 1
                        DO J=1,LENGTH(1)
                           QNAMES(N1+J) = QCOL(J)
                        END DO
                        LASTCH = LASTCH + ISHIFT
                        LENOBJ = LENGTH(1)
                     END DO
                  ENDIF
                  OBJNAM = DNAME(1)
               ELSEIF (OBJNAM .NE. DNAME(1)) THEN
                  MAXROW              = MAXROW + 1
                  NROWS               = NROWS + 1
                  NODEINFO(ND)%NROW   = NROWS
                  NODEINFO(ND)%NCOL   = NROWS
                  XUB(LASTBD+NROWS)   = PLINF
                  XLB(LASTBD+NROWS)   =-PLINF
                  VARTYP(IVTYP+NROWS) = 2
                  DROW = DNAME(1)
                  IF (LASTCH + LENGTH(1) .GT. IZCHAR) THEN
                     IF (.NOT. EXTMEM_CHAR(LASTCH+LENGTH(1),IERR))      &
     &                   GOTO 9200
                  ENDIF
                  DO I=1,LENGTH(1)
                     QNAMES(LASTCH+I) = QROW(I)
                  END DO
                  LASTCH = LASTCH + LENGTH(1)
                  NAME1(INAMES+NROWS+1) = LASTCH + 1
               ENDIF
!
!     Unrecognized code in code field. Ignore this row
!
            ELSE
               WRITE (IOLOG, 1400) NREC,QLINE
               WRITE (IOLOG, 1800) Q2,Q3
            ENDIF
            GOTO 150
!
!     COLUMNS or ARCS card. Initialize the columns section.
!
  200 CONTINUE
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 1400) NREC,QLINE
         IF (IENDAT .EQ. 0) GOTO 9700
            FIXCOL = 0
            IENDAT = 2
            NROWS  = NODEINFO(ND)%NROW
            NCOLS  = NODEINFO(ND)%NROW
            NELEM  = 0
            KMTX   = 1
            LDIAG  = 0
            LMN    = 0
            IF (Q1 .EQ. QC .AND. Q2 .EQ. QO) GOTO 203
               GOTO 204
!
!     Come back here to process the next record
!
  202    CONTINUE
            IF (MPSFORM) THEN
               CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOSTO,NREC,LENGTH)
            ELSE
               CALL GETNET(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOSTO,NREC,LENGTH)
            ENDIF
            IF (Q1 .EQ. QAST  ) GOTO 202
            IF (IERR   .GT.  0) GOTO 284
            IF (Q2 .EQ. QC .AND. Q3 .EQ. QP) GOTO 283
            IF (Q2 .EQ. QM .AND. Q3 .EQ. QK) GOTO 283
            IF (Q1 .EQ. QBL   ) GOTO 205
            IF (Q1 .EQ. QC .AND. Q2 .EQ. QO) GOTO 203
            IF (Q1 .EQ. QA .AND. Q2 .EQ. QR) GOTO 204
            IF (Q1 .EQ. QR .AND. Q2 .EQ. QH) GOTO 260
            IF (Q1 .EQ. QD .AND. Q2 .EQ. QE) GOTO 264
            IF (Q1 .EQ. QS .AND. Q2 .EQ. QU) GOTO 265
            IF (Q1 .EQ. QB .AND. Q2 .EQ. QO) GOTO 270
            IF (Q1 .EQ. QR .AND. Q2 .EQ. QA) GOTO 280
            IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 281
            IF (Q1 .EQ. QN .AND. Q2 .EQ. QA) GOTO 281
            IF (Q1 .EQ. QQ .AND. Q2 .EQ. QS) GOTO 282
               WRITE (IOLOG, 1400) NREC,QLINE
               WRITE (IOLOG, 2300)
               GOTO 9999
!
!     COLUMNS segment
!
  203    CONTINUE
            MPSFORM = .TRUE.
            INTTYP = 0
            DEFUPB = DEFUB
            DEFLOB = DEFLB
            GOTO 202
!
!     ARCS segment
!
  204    CONTINUE
            MPSFORM = .FALSE.
            WRITE (IOLOG, 2000)
            GOTO 202
!
!     Here we have a regular data card.
!     =================================
!
!     We must determine whether it continues a column already started,
!     which node it belongs to, whether it begins or ends a segment of
!     integer variables, and whether it even contains relevant data.
!
  205    CONTINUE
            IF (NECHO1 .GE. 2) WRITE (IOLOG, 1400) NREC,QLINE
            IF (MPSFORM) THEN
!
!     COLUMNS section.
!     ----------------
!
!     First check for integer and SOS markers.
!
               IF (DNAME(3) .EQ. DINT1) THEN
                  INTTYP = 1
                  GOTO 202
               ELSEIF (DNAME(3) .EQ. DINT2) THEN
                  IF (NOFSOS .GE. IZSTYP) THEN
                     IF (.NOT. EXTMEM_STYP(NOFSOS+1,IERR)) GOTO 9200
                  ENDIF
                  NOFSOS = NOFSOS + 1
                  INTTYP = -NOFSOS
                  IF     (Q2 .EQ. QS .AND. Q3 .EQ. '1') THEN
                     SOSTYP(NOFSOS) = 1
                  ELSEIF (Q2 .EQ. QS .AND. Q3 .EQ. '2') THEN
                     SOSTYP(NOFSOS) = 2
                  ELSE
!
!     SOS type 3. The first name field must contain a valid row name.
!     Each column participating in this set must be 0/1, the row must be
!     an equality row, with rhs = 1, coefficients in this row are all 1.
!
                     DCOL = DNAME(1)
                     KN = NODEINFO(ND)%KNAMES
                     NR = NODEINFO(ND)%NROW
                     IRES = NMSRCH(DCOL,KN+1,KN+NR)
                     IF (IRES .GT. 0) THEN
                        IRES = IRES - KN
                        IF (VARTYP(IRES+IVTYP) .NE. 0) GOTO 210
                           SOSTYP(NOFSOS) = - (IRES + IROW)
                           S3ROW  = IRES
                           DEFUPB = 1.D0
                           DEFLOB = 0.D0
                           GOTO 202
                     ENDIF
!
  210            CONTINUE
                     WRITE (IOLOG, 1400) NREC,QLINE
                     WRITE (IOLOG, 1900)
                     GOTO 9999
                  ENDIF
                  GOTO 202
!
               ELSEIF (DNAME(3) .EQ. DINT3 .OR.      &
     &                 DNAME(3) .EQ. DINT4) THEN
                  INTTYP = 0
                  DEFUPB = DEFUB
                  DEFLOB = DEFLB
                  GOTO 202
               ELSE
!
!     Here we have a regular MPS style column card. Test for row match.
!
                  JNM = 1
  220          CONTINUE
                  IF (DABS(ATEMP(JNM)) .LE. ZTOLIN) THEN
                      IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 202
                         JNM = 2
                  ENDIF
                  DROW = DNAME(JNM+1)
                  IRES = NMSRCH(DROW,INAMES+1,INAMES+NROWS)
                  IF (IRES .EQ. 0) THEN
!
!     No match!
!
                     WRITE (IOLOG, 1400) NREC,QLINE
                     WRITE (IOLOG, 1100) DROW
                     GOTO 9999
                  ELSE
!
!     Row was found in current period.
!
                     IR = IRES - INAMES
                     IF (DNAME(1) .EQ. DCOL) THEN
                        IF (IPC .EQ. IPER) GOTO 235
                           GO TO 241
                     ELSE
                        DCOL = DNAME(1)
                        DO IPC=1,IPER-1
                           IREF = KPREF(IPC)
                           KN = NODEINFO(IREF)%KNAMES
                           NR = NODEINFO(IREF)%NROW
                           NC = NODEINFO(IREF)%NCOL
                           IRES = NMSRCH(DCOL, KN+NR, KN+NC)
                           IF (IRES .GT. 0) GOTO 240
                        END DO
!
!      We have not seen the column before, so it must be new. Store the name.
!
                        IPC = IPER
                        IF (LASTCH + LENGTH(1) .GT. IZCHAR) THEN
                           IF (.NOT. EXTMEM_CHAR(LASTCH+LENGTH(1),IERR))      &
     &                         GOTO 9200
                        ENDIF
                        DO I=1,LENGTH(1)
                           QNAMES(LASTCH+I) = QCOL(I)
                        END DO
                        LASTCH = LASTCH + LENGTH(1)
!
!     Set column pointers and defaults for costs and bounds
!
                        NCOLS = NODEINFO(ND)%NCOL + 1
                        ICCA  = ICOLA + NCOLS - NROWS
                        JCOST = ICOST + NCOLS - NROWS
                        JVPTR = NODEINFO(ND)%VARPTR + NCOLS
                        NODEINFO(ND)%NCOL = NCOLS
!
!     Verify array sizes and expand memory if necessary/possible
!
                        IF (NCOLS .GT. IZCOLP) THEN
                           IF (.NOT. EXTMEM_COLP(NCOLS,IERR)) GOTO 9200
                        ENDIF
!
                        IF (ICCA .GE. IZACOL) THEN
                           IF (.NOT. EXTMEM_ACOL(ICCA+1,IERR)) GOTO 9200
                        ENDIF
!
                        IF (JCOST .GT. IZCOST) THEN
                           IF (.NOT. EXTMEM_COST(JCOST,IERR)) GOTO 9200
                        ENDIF
!
                        JCOL = ICOL + NCOLS
                        IF (JCOL .GT. IZCOLS) THEN
                           IF (.NOT. RESIZE      &
     &                          .OR. (IZCOLS .GE. MXCOLS)      &
     &                          .OR. (JCOL   .GT. MXCOLS)) GOTO 9200
                           IZCOLS = MAX(ICOL,MIN(MXCOLS, 2*IZCOLS))
                        ENDIF
!
                        IF (LASTNM + NCOLS .GE. IZVNAM) THEN
                           IF (.NOT. EXTMEM_VNAM(LASTNM+NCOLS+1,IERR))      &
     &                         GOTO 9200
                        ENDIF
                        NAME1(INAMES+NCOLS+1) = LASTCH + 1
!
                        JNAMES = INAMES + NCOLS
                        IF (JNAMES .GT. IZVNAM) THEN
                           IF (.NOT. EXTMEM_VNAM(JNAMES,IERR)) GOTO 9200
                        ENDIF
!
                        JBOUND = IBOUND + NCOLS
                        IF (JBOUND .GT. IZBNDS) THEN
                           IF (.NOT. EXTMEM_BNDS(JBOUND,IERR)) GOTO 9200
                        ENDIF
!
                        IF (JVPTR  .GT. IZVTYP) THEN
                           IF (.NOT. EXTMEM_VTYP(JVPTR,IERR)) GOTO 9200
                        ENDIF
!
                        COST(JCOST) = 0.D0
                        XLB(JBOUND) = DEFLOB
                        XUB(JBOUND) = DEFUPB
                        LA(ICCA)    = NELEM + 1
                        LA(ICCA+1)  = NELEM + 1
!
!     Store the column type and, if SOS type 3, the implied matrix coefficient
!
                        IF (INTTYP .GE. 0 .OR.      &
     &                             SOSTYP(NOFSOS) .GT. 0) THEN
                           VARTYP(JVPTR) = INTTYP
                        ELSE
                           VARTYP(JVPTR) = 1
                           NELEM = NELEM + 1
                           JELMA = IELMA + NELEM
                           IF (JELMA .GT. IZALMN) THEN
                              IF (.NOT. EXTMEM_ALMN(JELMA,IERR))      &
     &                           GOTO 9200
                           ENDIF
!
                           IA(JELMA)  = S3ROW
                            A(JELMA)  = 1.D0
                           LA(ICCA+1) = NELEM + 1
                        ENDIF
                     ENDIF
!
  235             CONTINUE
                     IF (IR .NE. IOBJ) THEN
!
!    Matched a coefficient in the A-matrix
!
                        NELEM = NELEM + 1
                        JELMA = IELMA + NELEM
                        IF (JELMA .GT. IZALMN) THEN
                           IF (.NOT. EXTMEM_ALMN(JELMA,IERR)) GOTO 9200
                        ENDIF
                        JCOLA = ICOLA+NCOLS-NROWS+1
                        IF (JCOLA .GT. IZACOL) THEN
                           IF (.NOT. EXTMEM_ACOL(JCOLA,IERR)) GOTO 9200
                        ENDIF
!
                        IA(JELMA) = IR
                         A(JELMA) = ATEMP(JNM)
                        LA(JCOLA) = NELEM + 1
                     ELSE
!
!     Cost coefficients (even if fixed) are not stored in the A-matrix
!
                        COST(JCOST) = ATEMP(JNM)
                     ENDIF
!
!     There might be another coefficient on this card.
!
                     IF (JNM .EQ. 2) GOTO 202
                     IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 202
                        JNM = 2
                        GO TO 220
!
!      Subdiagonal matrix entry. Identify the block this column belongs to.
!
  240             CONTINUE
                     IC0 = IRES - NODEINFO(IREF)%KNAMES      &
     &                          - NODEINFO(IREF)%NROW
  241             CONTINUE
                     IF (IR .NE. IOBJ) THEN
                        IF (LMN .GE. IZANZB) THEN
                           IF (.NOT. EXTMEM_ANZB(LMN+1,IERR)) GOTO 9200
                        ENDIF
                        LMN = LMN + 1
                        ARTMP(LMN) = IR
                        ACTMP(LMN) = IC0
                        AMTMP(LMN) = ATEMP(JNM)
                        IF (IPER-IPC .GT. LDIAG) LDIAG = IPER - IPC
                        IF (ISTART(IPER-IPC,1) .EQ. 0) THEN
                            ISTART(IPER-IPC,1) = LMN
                            ALINKS(LMN) = 0
                        ELSE
                            IL = ISTART(IPER-IPC,1)
                            IF (ACTMP(IL) .GE. IC0) THEN
                               ISTART(IPER-IPC,1) = LMN
                               ALINKS(LMN) = IL
                            ELSE
  250                       CONTINUE
                               IF (ALINKS(IL) .EQ. 0) THEN
                                   ALINKS(IL)  = LMN
                                   ALINKS(LMN) = 0
                               ELSEIF (ACTMP(ALINKS(IL)) .GE. IC0) THEN
                                   ALINKS(LMN) = ALINKS(IL)
                                   ALINKS(IL)  = LMN
                               ELSE
                                   IL = ALINKS(IL)
                                   GOTO 250
                               ENDIF
                           ENDIF
                        ENDIF
                     ELSE
!
!     Cost coefficient in a previous period --- should be ignored.
!
                        WRITE (IOLOG, 1400) NREC, QLINE
                        WRITE (IOLOG, 2700)
                     ENDIF
!
!     There may be another coefficient on this card. If not, get the next one
!
                     IF (JNM .EQ. 2) GOTO 202
                     IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 202
                        JNM = 2
                        GOTO 220
                  ENDIF
               ENDIF
            ELSE
!
!     ARCS section. Ignore the card.
!     ------------------------------
!
            ENDIF
            GOTO 202
!
!     The COLUMNS section is done. What is next?
!     ==========================================
!
!     RHS section
!
  260    CONTINUE
            L = 1
            SUPPLY  = .FALSE.
            MPSFORM = .TRUE.
            NAMFLD = 2
            GOTO 285
!
!     Demand segment
!
  264    CONTINUE
            L = 1
            SUPPLY  = .FALSE.
            MPSFORM = .FALSE.
            NAMFLD = 1
            GOTO 285
!
!     Supply segment
!
  265    CONTINUE
            L = 1
            SUPPLY  = .TRUE.
            MPSFORM = .FALSE.
            NAMFLD = 1
            GOTO 285
!
!     RANGES section
!
  270    CONTINUE
            L = 2
            MPSFORM = .TRUE.
            GOTO 285
!
!     BOUNDS section
!
  280    CONTINUE
            L = 3
            MPSFORM = .TRUE.
            GOTO 285
!
!     ENDATA card
!
  281    CONTINUE
            IENDAT = 4
            L = 4
            GOTO 285
!
!     QSECTION card
!
  282    CONTINUE
            IENDAT = 5
            L = 4
            GOTO 285
!
!     'CP' or 'MK' code (next node found)
!
  283    CONTINUE
            L = 5
            GOTO 285
!
!     End of file or error during read
!
  284    CONTINUE
            IF (IERR .GT. 1) GOTO 9400
               L = 5
               MPSFORM = .TRUE.
               Q1 = QE
               Q2 = QN
               Q3 = QD
!
!     Update the pointers
!
  285 CONTINUE
         IF (FIXCOL .EQ. 0) THEN
            FIXCOL = 1
            LASTA  = LASTA  + NELEM
            LASTNM = LASTNM + NODEINFO(ND)%NCOL + 1
            LASTC  = LASTC  - NODEINFO(ND)%NROW + NCOLS
            LASTCA = LASTCA - NODEINFO(ND)%NROW + NCOLS + 1
            LASTBD = LASTBD + NODEINFO(ND)%NCOL + 1
            LASTR  = LASTR  + NODEINFO(ND)%NROW
            LASTVP = LASTVP + NODEINFO(ND)%NCOL
            NODEINFO(ND)%NCOL = NCOLS
            NELMA(LASTBA) = NELEM
!
            IF (MARKOV .AND. IPER .GT. 1) THEN
               KMTX = 1
            ELSE
               KMTX = IPER - 1
            ENDIF
!
!     Switch to lower block-triangular form if necessary.
!
         IF (LDIAG .GT. KMTX) THEN
            KMTX = LDIAG
            MARKOV = .FALSE.
            DO J=NODES-1,1,-1
               NCTR = J
               JP   = 0
  286          CONTINUE
                  JP = JP + 1
                  NCTR = NODEINFO(NCTR)%IANCTR
                  IF (NCTR .GT. 0) GOTO 286
!
!     If node J is in stage 3 or later, we need to shift pointers around
!
               IF (JP .GT. 2) THEN
                  KSHIFT = JP - 2
                  MNABLK = LASTBA + KSHIFT
                  IF (MNABLK .GT. IZABLK) THEN
                     IF (.NOT. EXTMEM_ABLK(MNABLK,IERR)) GOTO 9200
                  ENDIF
!
                  DO K=LASTBA,KDATA(J)+2,-1
                     KCOLA(K+KSHIFT) = KCOLA(K)
                     KELMA(K+KSHIFT) = KELMA(K)
                     NELMA(K+KSHIFT) = NELMA(K)
                  END DO
                  LASTBA = LASTBA + KSHIFT
!
                  DO K=1,NODES
                     IF (KDATA(K) .GT. KDATA(J))      &
     &                   KDATA(K) = KDATA(K) + KSHIFT
                  END DO
!
                  NCOLA = NODEINFO(J)%NCOL + 1 - NODEINFO(J)%NROW
                  MNACOL = LASTCA + (JP-3)*NCOLA
                  IF (MNACOL .GT. IZACOL) THEN
                     IF (.NOT. EXTMEM_ACOL(MNACOL,IERR)) GOTO 9200
                  ENDIF
!
                  DO K=3,JP
                     KCOLA(KDATA(J)+K) = LASTCA
                     KELMA(KDATA(J)+K) = LASTA
                     NELMA(KDATA(J)+K) = 0
                     DO JC=1,NCOLA
                        LA(LASTCA+JC) = 1
                     END DO
                     LASTCA = LASTCA + NCOLA
                  END DO
!
               ENDIF
            END DO
            KMTX = IPER - 1
         ENDIF
!
!     Copy subdiagonal blocks to permanent positions.
!
            LOFF  = LASTCA
            NABLK = LASTBA
            NALMN = LASTA
            MNABLK = LASTBA + KMTX
            IF (MNABLK .GT. IZABLK) THEN
              IF (.NOT. EXTMEM_ABLK(MNABLK,IERR)) GOTO 9200
            ENDIF
!
            MNALMN = NELEM + LMN
            IF (MNALMN .GT. IZALMN) THEN
               IF (.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 9200
            ENDIF
!
            DO J=1,KMTX
               JJ = KPREF(IPER-J)
               NC = NODEINFO(JJ)%NCOL - NODEINFO(JJ)%NROW
               MNACOL = LASTCA + NC + 1
               IF (MNACOL .GT. IZACOL) THEN
                  IF (.NOT. EXTMEM_ACOL(MNACOL,IERR)) GOTO 9200
               ENDIF
!
               NABLK = NABLK + 1
               KELMA(NABLK)  = NALMN
               KCOLA(NABLK)  = LOFF
               LA(LOFF+1) = NALMN + 1 - KELMA(NABLK)
               LC = 0
               IF (ISTART(J,1) .GT. 0) THEN
                  IL = ISTART(J,1)
!
  296          CONTINUE
                  IC = ACTMP(IL)
                  DO K=LC+1,IC
                     LA(LOFF+K) = NALMN + 1 - KELMA(NABLK)
                  END DO
                  NALMN = NALMN + 1
!
                   A(NALMN) = AMTMP(IL)
                  IA(NALMN) = ARTMP(IL)
                  IL = ALINKS(IL)
                  LC = IC
                  IF (IL .GT. 0) GOTO 296
               ENDIF
!
               DO K=LC,NC
                  LA(LOFF+K+1) = NALMN + 1 - KELMA(NABLK)
               END DO
               NELMA(NABLK) = NALMN - KELMA(NABLK)
               LOFF   = LOFF + NC + 1
               LASTCA = LOFF
               LASTBA = NABLK
               LASTA  = NALMN
            END DO
!
!    Initial the RHS
!
            IF (.NOT. ALLOCATED(RHS)) THEN
               ALLOCATE(RHS(IZDRHS), STAT=IERR)
               IF (IERR .NE. 0) GOTO 9200
            ENDIF
            KR = NODEINFO(ND)%KRHS
            NR = NODEINFO(ND)%NROW
            MNDRHS = KR + NR
            IF (MNDRHS .GT. IZDRHS) THEN
               IF (.NOT. EXTMEM_DRHS(MNDRHS,IERR)) GOTO 9200
            ENDIF
            DO J=1,NR
               RHS(KR+J) = 0.D0
            END DO
!
            IF (L .EQ. 4 .AND. IENDAT .EQ. 4) GOTO 450
            IF (L .EQ. 4 .AND. IENDAT .EQ. 5) GOTO 449
            IF (L .EQ. 5) GOTO 900
         ENDIF
!
         I0  = 1
!
!     Process right-hand sides (including supply and demand), bounds and ranges
!
  301 CONTINUE
         IF (MPSFORM) THEN
            CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                        IERR,IOSTO,NREC,LENGTH)
         ELSE
            CALL GETNET(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                        IERR,IOSTO,NREC,LENGTH)
         ENDIF
         IF (Q1 .EQ. QAST ) GOTO 301
         IF (IERR   .GT. 0) GOTO 305
         IF (Q2 .EQ. QC .AND. Q3 .EQ. QP) GOTO 450   
         IF (Q2 .EQ. QM .AND. Q3 .EQ. QK) GOTO 450
         IF (Q1 .EQ. QBL  ) GOTO 309
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QH) GOTO 260
         IF (Q1 .EQ. QD .AND. Q2 .EQ. QE) GOTO 264
         IF (Q1 .EQ. QS .AND. Q2 .EQ. QU) GOTO 265
         IF (Q1 .EQ. QB .AND. Q2 .EQ. QO) GOTO 270
         IF (Q1 .EQ. QR .AND. Q2 .EQ. QA) GOTO 280
         IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 306
         IF (Q1 .EQ. QN .AND. Q2 .EQ. QA) GOTO 307
         IF (Q1 .EQ. QQ .AND. Q2 .EQ. QS) GOTO 308
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 2300)
            GOTO 9999
!
!    Error during read or end of file detected
!
  305    CONTINUE
            IF (IERR .GT. 1) GOTO 9400
               Q1 = QE
               Q2 = QN
               Q3 = QD
               GOTO 450
!
!     END record
!
  306    CONTINUE
            IENDAT = 4
            GOTO 301
!
!     NAME record starts the next node
!
  307    CONTINUE
            IENDAT = 4
            GOTO 450
!
!     QSECT record
!
  308    CONTINUE
            IENDAT = 5
            GOTO 449
!
  309    CONTINUE
            IF (NECHO1 .GE. 2) WRITE (IOLOG, 1400) NREC,QLINE
            GOTO (310,350,400,450), L
!
!     RHS section (includes SUPPLY and DEMAND)
!
  310    CONTINUE
         J = NAMFLD
         IF (DABS(ATEMP(1)) .GT. ZTOLIN) GOTO 312
         IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 301
            J = 3
            ATEMP(1) = ATEMP(2)
!
!        Test for row match
!
  312    CONTINUE
         IF (DRHS .EQ. DOTS) DRHS = DNAME(1)
         IF (MPSFORM .AND. (DRHS .NE. DNAME(1))) GOTO 301
         DROW = DNAME(J)
         IREF = KPREF(IPER)
         IRES = NMSRCH(DROW,INAMES+I0,INAMES+NODEINFO(IREF)%NROW)
         IF (IRES .GT. 0) GOTO 330
         IRES = NMSRCH(DROW,INAMES+1 ,INAMES+         I0        )
         IF (IRES .GT. 0) GOTO 330
!
!     No match!
!
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 1200) DROW
         GOTO 9999
!
!     Matched
!
  330    CONTINUE
         I0 = IRES - NODEINFO(IREF)%KNAMES
         IF (SUPPLY) THEN
            RHS(NODEINFO(IREF)%KRHS+I0) = -ATEMP(1)
         ELSE
            RHS(NODEINFO(IREF)%KRHS+I0) =  ATEMP(1)
         ENDIF
         IF (J .EQ. 1 .OR. J .EQ. 3) GOTO 301
            IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 301
               J = 3
               ATEMP(1) = ATEMP(2)
               GOTO 312
!
!     BOUNDS Section. Match the column name.
!
  350 CONTINUE
         IF (DBOUND .EQ. DOTS     ) DBOUND = DNAME(1)
         IF (DBOUND .NE. DNAME(1) ) GOTO 301
         DROW = DNAME(2)
         IREF = KPREF(IPER)
         IRES = NMSRCH(DROW,INAMES+I0,INAMES+NODEINFO(IREF)%NCOL)
         IF (IRES .GT. 0) GOTO 360
         IRES = NMSRCH(DROW,INAMES+1 ,INAMES+I0        )
         IF (IRES .GT. 0) GOTO 360
!
!     No match!
!
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 1300) DROW
            GOTO 9999
!
!     Matched. Now determine the bound type
!
  360    CONTINUE
         I0  = IRES - NODEINFO(IREF)%KNAMES
         IC  = NODEINFO(IREF)%KBOUND + I0
         IF (Q2 .EQ. QL .AND. Q3 .EQ. QO) GOTO 361
         IF (Q2 .EQ. QU .AND. Q3 .EQ. QP) GOTO 366
         IF (Q2 .EQ. QF .AND. Q3 .EQ. QX) GOTO 365
         IF (Q2 .EQ. QF .AND. Q3 .EQ. QR) GOTO 370
         IF (Q2 .EQ. QM .AND. Q3 .EQ. QI) GOTO 368
         IF (Q2 .EQ. QP .AND. Q3 .EQ. QL) GOTO 372
         IF (Q2 .EQ. QU .AND. Q3 .EQ. QI) GOTO 374
         IF (Q2 .EQ. QB .AND. Q3 .EQ. QV) GOTO 376
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 1700)
            GOTO 9999
!
  361    CONTINUE
            XLB(IC) = ATEMP(1)
            GOTO 301
  365    CONTINUE
            XLB(IC) = ATEMP(1)
  366    CONTINUE
            XUB(IC) = ATEMP(1)
            GOTO 301
  368    CONTINUE
            XLB(IC) = -PLINF
            GOTO 301
  370    CONTINUE
            XLB(IC) = -PLINF
  372    CONTINUE
            XUB(IC) =  PLINF
            GOTO 301
  374    CONTINUE
            XUB(IC) =  DINT(ATEMP(1))
            GOTO 301
  376    CONTINUE
            XLB(IC) =  0.D0
            XUB(IC) =  1.D0
            GOTO 301
!
!     RANGES Section. Match the row name.
!
  400    CONTINUE
         IF (DRANGE .EQ. DOTS     ) DRANGE = DNAME(1)
         IF (DRANGE .NE. DNAME(1) ) GOTO 301
         J = 2
         IF (DABS(ATEMP(1)) .GT. ZTOLIN) GOTO 412
         IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 301
            J = 3
            ATEMP(1) = ATEMP(2)
!
!        Test for row match
!
  412    CONTINUE
         DROW = DNAME(J)
         IREF = KPREF(IPER)
         IRES = NMSRCH(DROW,INAMES+I0,INAMES+NODEINFO(IREF)%NROW)
         IF (IRES .GT. 0) GOTO 430
         IRES = NMSRCH(DROW,INAMES+1 ,INAMES   +      I0        )
         IF (IRES .GT. 0) GOTO 430
!
!     No match!
!
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 1200) DROW
         GOTO 9999
!
!     Matched
!
  430    CONTINUE
         I0  = IRES - NODEINFO(IREF)%KNAMES
         IR  = NODEINFO(IREF)%KBOUND + I0
         IT  = VARTYP(IVTYP+I0)
         IF (IT .EQ.  1) GOTO 435
         IF (IT .EQ. -1) GOTO 440
         IF (IT .NE.  0) GOTO 9600
            IF (ATEMP(1) .GT. 0.) XUB(IR) = ATEMP(1)
            IF (ATEMP(1) .LT. 0.) XLB(IR) =-ATEMP(1)
            GOTO 442
  435    CONTINUE
            XUB(IR) =  DABS(ATEMP(1))
            GOTO 442
  440    CONTINUE
            XLB(IR) = -DABS(ATEMP(1))
  442    CONTINUE
         IF (J .EQ. 3) GOTO 301
            IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 301
               J = 3
               ATEMP(1) = ATEMP(2)
               GOTO 412
!
!     End of core file
!
  449 CONTINUE
         IF (.NOT. QMAT) THEN
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 2500)
            QMAT = .TRUE.
            ALLOCATE(AQMTX(IZQLMN), IQMTX(IZQLMN), LQMTX(IZQCOL),      &
     &               KDATQ(IZNODE), IQOFF(IZQBLK), LQOFF(IZQBLK),      &
     &               NELMQ(IZQBLK), STAT=IER)
            IF (IER .NE. 0) GOTO 9200
            DO I=1,NODES
               KDATQ(I) = 0
            END DO
            LASTQB = 0
         ENDIF
  450 CONTINUE
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 1400) NREC,QLINE
!
!     Update the pointers
!
         LASTR  = NODEINFO(ND)%KRHS   + NODEINFO(ND)%NROW
         LASTBD = NODEINFO(ND)%KBOUND + NCOLS + 1
         LASTVP = NODEINFO(ND)%VARPTR + NCOLS
!
!     Store bounds for the theta column
!
         XINFTY = 1.D31
         XLB(LASTBD) = -XINFTY
         XUB(LASTBD) =  XINFTY
!
!     Read a QSECTION if necessary
!
!     INITIALIZE POINTERS
!
      IF (Q1 .EQ. QQ .AND. Q2 .EQ. QS) THEN
      DO I=1,IPER
         ISTART(I,2) = 0
      END DO
!
      LDIAG  = 0
      DCOL1  = DBLANK
      DCOL2  = DBLANK
      IP0 = 1
      IC0 = 1
      LMN = 0
!
  550 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,      &
     &                     DNAME,ATEMP,IERR,IOSTO,NREC,LENGTH)
         IF (Q1 .EQ. QAST  ) GOTO 550
         IF (IERR   .GT.  0) GOTO 554
         IF (Q2 .EQ. QC .AND. Q3 .EQ. QP) GOTO 780
         IF (Q2 .EQ. QM .AND. Q3 .EQ. QK) GOTO 780
         IF (Q1 .EQ. QBL   ) GOTO 560
         IF (Q1 .EQ. QN .AND. Q2 .EQ. QA) GOTO 550
         IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 550
         IF (Q1 .EQ. QQ .AND. Q2 .EQ. QS) GOTO 550
            WRITE (IOLOG, 1400) NREC,QLINE
            WRITE (IOLOG, 2100)
            GOTO 9999
!
!     End of file or error during read.
!
  554 CONTINUE
         IF (IERR .GT. 1) GOTO 9400
            Q1 = QE
            Q2 = QN
            Q3 = QD
            GOTO 780
!
!    Here we have a regular data card.
!
  560 CONTINUE
         IF (NECHO1 .GE. 2) WRITE (IOLOG, 1400) NREC,QLINE
         IF (DABS(ATEMP(1)) .LT. ZTOLIN) GOTO 550
         IF (DNAME(1) .EQ. DCOL1) GOTO 740
!
!     Match the first column name
!
         DCOL1 = DNAME(1)
         IP = IP0
         IR = KPREF(IP)
         KN = NODEINFO(IR)%KNAMES
         NR = NODEINFO(IR)%NROW
         NC = NODEINFO(IR)%NCOL
         IRES = NMSRCH(DCOL1,KN+NR+IC0,KN+NC)
         IF (IRES .GT. 0) GOTO 730
         IRES = NMSRCH(DCOL1,KN+NR+1,  KN+NR+IC0)
         IF (IRES .GT. 0) GOTO 730
!
         DO IP=IP0+1,IPER
            IR = KPREF(IP)
            KN = NODEINFO(IR)%KNAMES
            NR = NODEINFO(IR)%NROW
            NC = NODEINFO(IR)%NCOL
            IRES = NMSRCH(DCOL1,KN+NR+1,KN+NR+NC)
            IF (IRES .GT. 0) GOTO 730
         END DO
         DO IP=1,IP0-1
            IR = KPREF(IP)
            KN = NODEINFO(IR)%KNAMES
            NR = NODEINFO(IR)%NROW
            NC = NODEINFO(IR)%NCOL
            IRES = NMSRCH(DCOL1,KN+NR+1,KN+NR+NC)
            IF (IRES .GT. 0) GOTO 730
         END DO
!
!     Not found!
!
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 2200) DCOL1
         GOTO 9999
!
!     Matched
!
  730 CONTINUE
         IP0 = IP
         IC0 = IRES - NODEINFO(IR)%KNAMES - NODEINFO(IR)%NROW
!
!     Match the second column name
!
  740 CONTINUE
         DCOL2 = DNAME(2)
         IP = IP0
         IR = KPREF(IP)
         KN = NODEINFO(IR)%KNAMES
         NR = NODEINFO(IR)%NROW
         NC = NODEINFO(IR)%NCOL
         IRES = NMSRCH(DCOL2,KN+NR+1,KN+NC)
         IF (IRES .GT. 0) GOTO 750
         DO IP=1,IPER
            IR = KPREF(IP)
            KN = NODEINFO(IR)%KNAMES
            NR = NODEINFO(IR)%NROW
            NC = NODEINFO(IR)%NCOL
            IRES = NMSRCH(DCOL2,KN+NR+1,KN+NR+NC)
            IF (IRES .GT. 0) GOTO 750
         END DO
!
!     Not found!
!
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 2200) DCOL2
         GOTO 9999
!
!   Matched. Verify that the element is in the lower triangular part,
!   then store in temporary data structure (linked lists). If in the
!   upper half, warn the user and ignore the entry.
!
  750 CONTINUE
      IF ((IP0 .LT. IP .OR. (IP .EQ. IP0 .AND. IC0 .LE. IR)) .AND.      &
     &                      (IP .EQ. IPER) ) THEN
         IR = IRES - KN - NR
         LMN = LMN + 1
         IF (LMN .GT. IZQLMN) THEN
            IF (.NOT. EXTMEM_QLMN(LMN,IERR)) GOTO 9200
         ENDIF
!
         QRTMP(LMN) = IR
         QCTMP(LMN) = IC0
         QMTMP(LMN) = ATEMP(1)
         IF (IABS(IP0-IP) .GE. LDIAG) LDIAG = IABS(IP0-IP) + 1
         IF (ISTART(IP0,2) .EQ. 0) THEN
            ISTART(IP0,2) = LMN
            QLINKS(LMN) = 0
         ELSE
            IL = ISTART(IP0,2)
            IF (QCTMP(IL) .GE. IC0) THEN
               ISTART(IP0,2) = LMN
               QLINKS(LMN) = IL
            ELSE
  760       CONTINUE
               IF (QLINKS(IL) .EQ. 0) THEN
                   QLINKS(IL)  = LMN
                   QLINKS(LMN) = 0
               ELSEIF (QCTMP(QLINKS(IL)) .GE. IC0) THEN
                   QLINKS(LMN) = QLINKS(IL)
                   QLINKS(IL)  = LMN
               ELSE
                   IL = QLINKS(IL)
                   GOTO 760
               ENDIF
            ENDIF
         ENDIF
      ELSE
         WRITE (IOLOG, 1400) NREC, QLINE
         WRITE (IOLOG, 2800)
      ENDIF
         GOTO 550
!
!     End of file. Copy data to permanent storage positions
!
  780 CONTINUE
         IF (LDIAG .GT. NQPROF) THEN
!
!     The Q-matrix must be expanded; this involves blocks as well as column
!     pointers. First count the number of array elements to add.
!
            ALLOCATE(MREF(NODES,2), STAT=IERR)
            IF (IERR .NE. 0) GOTO 9200
            KSHIFT = 0
            LSHIFT = 0
            DO J=NODES-1,1,-1
               NCTR = J
               JP   = 0
               MREF(J,2) = 0
  781          CONTINUE
                  JP = JP + 1
                  IF (JP .GT. NQPROF .AND. JP .LE. LDIAG) THEN
                     KSHIFT = KSHIFT + 1
                     LSHIFT = LSHIFT + NODEINFO(NCTR)%NCOL      &
     &                               - NODEINFO(NCTR)%NROW + 1
                     IF (JP .EQ. NQPROF+1) MREF(J,2) = NCTR
                  ENDIF
                  NCTR = NODEINFO(NCTR)%IANCTR
                  IF (NCTR .GT. 0) GOTO 781
                     MREF(J,1) = JP
            END DO
            KNEW = 0
            LNEW = 0
            NCTR = NODES
            JP   = 0
 782        CONTINUE
               JP = JP + 1
               IF (JP .GT. NQPROF .AND. JP .LE. LDIAG) THEN
                  KNEW = KNEW + 1
                  LNEW = LNEW + NODEINFO(NCTR)%NCOL      &
     &                        - NODEINFO(NCTR)%NROW + 1
               ENDIF
               NCTR = NODEINFO(NCTR)%IANCTR
               IF (NCTR .GT. 0) GOTO 782
!
!     Make sure there is enough room
!
            MNQBLK = LASTQB + KSHIFT + KNEW
            IF (MNQBLK .GT. IZQBLK) THEN
               IF (.NOT. EXTMEM_QBLK(MNQBLK,IERR)) GOTO 9200
            ENDIF
            MNQCOL = LASTQC + LSHIFT + LNEW
            IF (MNQCOL .GT. IZQCOL) THEN
               IF (.NOT. EXTMEM_QCOL(MNQCOL,IERR)) GOTO 9200
            ENDIF
!
!     Now shift the blocks
!
            JSHIFT = KSHIFT
            DO J=NODES-1,1,-1
               JP = MIN(LDIAG,MREF(J,1))
               JSHIFT = JSHIFT - MAX(0,JP-NQPROF)
               KDATQ(J) = KDATQ(J) + JSHIFT
               JP0 = MIN(MREF(J,1),NQPROF)
               DO KK=KDATQ(J)+JP0,KDATQ(J)+1,-1
                  IQOFF(KK) = IQOFF(KK-JSHIFT)
                  LQOFF(KK) = LQOFF(KK-JSHIFT)
                  NELMQ(KK) = NELMQ(KK-JSHIFT)
               END DO
!
!     Insert dummy blocks
!
               IF (MREF(J,1) .GT. NQPROF) THEN
                  NCTR = MREF(J,2)
                  DO KK=NQPROF+1,JP
                     IQOFF(KDATQ(J)+KK) = LASTQ
                     LQOFF(KDATQ(J)+KK) = LASTQC
                     NELMQ(KDATQ(J)+KK) = 0
                     NLC = NODEINFO(NCTR)%NCOL - NODEINFO(NCTR)%NROW + 1
                     DO LL=1,NLC
                        LQMTX(LASTQC+LL) = 1
                     END DO
                     LASTQC = LASTQC + NLC
                     NCTR = NODEINFO(NCTR)%IANCTR
                  END DO
               ENDIF
            END DO
            LASTQB = LASTQB + KSHIFT
            NQPROF = LDIAG
            DEALLOCATE(MREF, STAT=IERR)
         ENDIF
!
!     Place the info for the current node into the data structure
!
         IF (LDIAG .GT. 0) THEN
            LOFF  = LASTQC
            NQBLK = LASTQB
            NQLMN = LASTQ
!
!     Main diagonal is first
!
            KDATQ(ND)    = NQBLK
            NQBLK        = NQBLK + 1
!
            NCOLQ = NODEINFO(ND)%NCOL-NODEINFO(ND)%NROW + 1
            MNQCOL = LASTQC + NCOLQ
!
            IQOFF(NQBLK) = NQLMN
            LQOFF(NQBLK) = LOFF
            NELMQ(NQBLK) = NQLMN
            LC = 0
            IF (ISTART(IPER,2) .GT. 0) THEN
               IL = ISTART(IPER,2)
  810       CONTINUE
               IC = QCTMP(IL)
               DO K=LC+1,IC
                  LQMTX(LOFF+K) = NQLMN + 1 - IQOFF(NQBLK)
               END DO
               NQLMN = NQLMN + 1
               AQMTX(NQLMN) = QMTMP(IL)
               IQMTX(NQLMN) = QRTMP(IL)
               IL = QLINKS(IL)
               LC = IC
               IF (IL .GT. 0) GOTO 810
            ENDIF
            NELMQ(NQBLK) = NQLMN - NELMQ(NQBLK)
            DO K=LC,NCOLQ-1
               LQMTX(LOFF+K+1) = NQLMN + 1 - IQOFF(NQBLK)
            END DO
            LOFF = LOFF + NCOLQ
            J = IPER
!
!     Now deal with sub-diagonals, if there are any
!
  825       CONTINUE
               J = J - 1
               IF (J .LE. 0 .OR. IPER-J .GE. NQPROF) GOTO 890
!
!     Here we copy the subdiagonal matrix
!
                  NQBLK = NQBLK + 1
!
                  JJ = KPREF(J)
                  NCOLQ = NODEINFO(JJ)%NCOL - NODEINFO(JJ)%NROW + 1
                  MNQCOL = LASTQC + NCOLQ
!
                  IQOFF(NQBLK)  = NQLMN
                  LQOFF(NQBLK)  = LOFF
                  NELMQ(NQBLK)  = NQLMN
                  LQMTX(LOFF+1) = NQLMN + 1 - IQOFF(NQBLK)
                  LC = 0
                  IF (ISTART(J,2) .GT. 0) THEN
                     IL = ISTART(J,2)
  830                CONTINUE
                        IC = QCTMP(IL)
                        DO K=LC+1,IC
                           LQMTX(LOFF+K) = NQLMN + 1 - IQOFF(NQBLK)
                        END DO
                        NQLMN = NQLMN + 1
                        AQMTX(NQLMN) = QMTMP(IL)
                        IQMTX(NQLMN) = QRTMP(IL)
                        IL = QLINKS(IL)
                        LC = IC
                        IF (IL .GT. 0) GOTO 830
                  ENDIF
                  NELMQ(NQBLK)  = NQLMN - NELMQ(NQBLK)
                  JJ = KPREF(J)
                  DO K=LC,NCOLQ-1
                     LQMTX(LOFF+K+1) = NQLMN + 1 - IQOFF(NQBLK)
                  END DO
                  LOFF = LOFF + NCOLQ
                  GOTO 825
!
  890       CONTINUE
               LASTQC = LOFF
               LASTQ  = NQLMN
               LASTQB = NQBLK
            ENDIF
!
!     If there is no QSECT in this node, but one was defined earlier,
!     adjust the pointers
!
      ELSEIF (NQPROF .GT. 0) THEN
         NQBLK = MIN(IPER,NQPROF)
         NLCOL = 0
         NCTR  = ND
         DO JJ=1,NQBLK
            NLCOL = NLCOL + NODEINFO(NCTR)%NCOL      &
     &                    - NODEINFO(NCTR)%NROW + 1
            NCTR  = NODEINFO(NCTR)%IANCTR
         END DO
!
         MNQBLK = LASTQB + NQBLK
         IF (MNQBLK .GT. IZQBLK) THEN
            IF (.NOT. EXTMEM_QBLK(MNQBLK,IERR)) GOTO 9200
         ENDIF
         MNQCOL = LASTQC + NLCOL
         IF (MNQCOL .GT. IZQCOL) THEN
            IF (.NOT. EXTMEM_QCOL(MNQCOL,IERR)) GOTO 9200
         ENDIF
!
         NCTR = ND
         DO JJ=1,NQBLK
            IQOFF(KDATQ(ND)+JJ) = LASTQ
            LQOFF(KDATQ(ND)+JJ) = LASTQC
            NELMQ(KDATQ(ND)+JJ) = 0
            NLCOL = NODEINFO(NCTR)%NCOL - NODEINFO(NCTR)%NROW + 1
            NCTR  = NODEINFO(NCTR)%IANCTR
            DO JL=1,NLCOL
               LQMTX(LASTQC+JL) = 1
            END DO
            LASTQC = LASTQC + NLCOL
         END DO
         LASTQB = LASTQB + NQBLK
      ENDIF
!
!     Next node card or end of file
!
  900 CONTINUE
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 1400) NREC,QLINE
!
         IF (Q2 .EQ. QN .AND. Q3 .EQ. QD) NDINFO = 0
         IF (Q2 .EQ. QC .AND. Q3 .EQ. QP) NDINFO = 1
         IF (Q2 .EQ. QM .AND. Q3 .EQ. QK) NDINFO = 2
         IF (IERR .EQ. 1) IERR = 0
         GOTO 999
!
!     Come here if anything went wrong
!
 9200 CONTINUE
         WRITE (IOLOG, 3200)
         GOTO 9999
!
 9300 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3300)
         GOTO 9999
!
 9400 CONTINUE
         WRITE (IOLOG, 3400)
         GOTO 9999
!
 9600 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3600)
         GOTO 9999
!
 9700 CONTINUE
         WRITE (IOLOG, 1400) NREC,QLINE
         WRITE (IOLOG, 3700)
         GOTO 9999
!
 9999 CONTINUE
         IERR = 2
  999 CONTINUE
         DEALLOCATE(ISTART, QLINKS, QRTMP, QCTMP, QMTMP,      &
     &                      ALINKS, ARTMP, ACTMP, AMTMP, STAT=IER)
         IF (IER .NE. 0) IERR = 2
         IF (ALLOCATED(MREF)) THEN
            DEALLOCATE(MREF, STAT=IER)
            IF (IER .NE. 0) IERR = 2
         ENDIF
         RETURN
!
 1100 FORMAT(' XXX -  FATAL  - Row "',A20,'" was never defined in',      &
     &       ' ROWS section.')
 1200 FORMAT(' XXX -  FATAL  - Unmatched row name "',A20,      &
     &       '" in RHS or RANGES section.')
 1300 FORMAT(' XXX -  FATAL  - Unmatched column name "',A20,      &
     &       '" in BOUNDS section.')
 1400 FORMAT(I8,2X,A80)
 1700 FORMAT(' XXX -  FATAL  - Error in BOUNDS section.')
 1800 FORMAT(' XXX - WARNING - Unrecognized code =',2A1,'= in ROWS',      &
     &       ' section. Row ignored.')
 1900 FORMAT(' XXX -  FATAL  - Inappropriate row name in SOS type 3.')
 2000 FORMAT(' XXX - WARNING - Node format cannot process ARCS cards.',      &
     &       ' Flush to next header.')
 2100 FORMAT(' XXX -  FATAL  - Illegal header card in QSECTION.')
 2200 FORMAT(' XXX -  FATAL  - Unmatched column name ',A8,      &
     &       ' in QSECTION')
 2300 FORMAT(' XXX -  FATAL  - Illegal header card in CORE file.')
 2500 FORMAT(' XXX - WARNING - QSECTION detected but not present in',      &
     &       ' core file.')
 2600 FORMAT(' XXX - WARNING - Node contains no information.')
 2700 FORMAT(' XXX - WARNING - Cost coefficient in previous period.',      &
     &       ' Data ignored.')
 2800 FORMAT(' XXX - WARNING - Q-matrix element in upper triangle.',      &
     &       ' Ignore.')
 3200 FORMAT(' XXX -  FATAL  - Memory allocation failed in MKNODE.')
 3300 FORMAT(' XXX -  FATAL  - No COLUMNS section specified.')
 3400 FORMAT(' XXX -  FATAL  - Error while rading stoch file.')
 3600 FORMAT(' XXX -  FATAL  - Illegal row type in RANGES section.')
 3700 FORMAT(' XXX -  FATAL  - ROWS section is non-existent.')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CPNODE ( ND, NDINFO, IERR, NREC, IREFND,      &
     &                    MODE, DNAME, ATEMP, Q1, Q2, Q3, Q4 )
!
!     Routine:  CPNODE
!
!     Purpose:  This subroutine copies the data for one node in the explicitly
!               given event tree to another. The input format is described in
!               further detail in the paper by Gassmann and Schweitzer.
!
!     Calling Sequence:
!
!           CALL CPNODE ( ND, NDINFO, IERR, NREC,
!     *                   MODE, DNAME, ATEMP, Q1, Q2, Q3, Q4 )
!
!     Parameters:
!
!         ND      Number of the node to be processed
!         NDINFO  Format of the *next* node to be processed
!                 (has to be tracked because the end of a node may be
!                 defined implicitly by encountering the next node header)
!           = 0   END card has been detected; no format implied
!           = 1   CP header has been detected: copy the next node
!           = 2   MK header has been detected: build the next node
!         IERR    Error status
!           = 0   Normal termination
!           > 0   Error condition
!         NREC    Number of records read
!         IREFND  Reference node from which to copy the information
!         MODE    Descibes the location of the node information
!           = 1   Information is found in a previously recorded node
!           = 3   Information is found in the core file
!         DNAME   Name fields on the MPS-style data card
!         ATEMP   Numerical fields on the MPS-style data card
!         Q1      First character of the code field
!         Q2      Second character of the code field
!         Q3      Third character of the code field
!         Q4      Fourth character of the code field
!
!
!     Calls:    GETMPS, IDENTIFY
!
!     Method:   First copy all the information from the reference node,
!               then process the stoch file one record at a time and
!               replace the information.
!
!     Remarks:  None
!
!     -----------------------------------------------------
!     Version of 14 October 2002. Written by Gus Gassmann.
!     -----------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME(*), DBLANK, DROW, DCOL, QLINE
      CHARACTER*1 QBLANK(NAMLEN)
      INTEGER*4   LENGTH(3)
      REAL*8      ATEMP(*)
!
      INTEGER, ALLOCATABLE, DIMENSION (:) :: KPREF
      EQUIVALENCE (QBLANK,DBLANK)
!
      DATA QBLANK/NAMLEN*' '/
!
!
!     ALLOCATE TEMPORARY ARRAY
!
      ALLOCATE (KPREF(NPER), STAT=IERR)
      IF (IERR .NE. 0) GOTO 9200
!
!     Initialize pointers
!
      IF (ND .NE. IREFND) THEN
         KDATA(ND) = LASTBA
         IF (ALLOCATED(KDATQ)) KDATQ(ND) = LASTQB
         NODEINFO(ND)%KBOUND = NODEINFO(IREFND)%KBOUND
         NODEINFO(ND)%NROW   = NODEINFO(IREFND)%NROW
         NODEINFO(ND)%NCOL   = NODEINFO(IREFND)%NCOL
         NODEINFO(ND)%KRHS   = NODEINFO(IREFND)%KRHS
         NODEINFO(ND)%KCOST  = NODEINFO(IREFND)%KCOST
         NODEINFO(ND)%KNAMES = NODEINFO(IREFND)%KNAMES
         NODEINFO(ND)%VARPTR = NODEINFO(IREFND)%VARPTR
         NODEINFO(ND)%KROW   = NODEINFO(ND-1)%KROW      &
     &                       + NODEINFO(ND-1)%NROW
         NODEINFO(ND)%KCOL   = NODEINFO(ND-1)%KCOL      &
     &                       + NODEINFO(ND-1)%NCOL + 1
         IF (NODEINFO(ND)%KROW .GT. IZROWS) THEN
            IF (.NOT. RESIZE .OR. (NODEINFO(ND)%KROW .GT. MXROWS))      &
     &         GOTO 9200
            IZROWS = MIN(MXROWS,MAX(2*IZROWS,NODEINFO(ND)%KROW))
         ENDIF
         IF (NODEINFO(ND)%KCOL .GT. IZCOLS) THEN
            IF (.NOT. RESIZE .OR. (NODEINFO(ND)%KCOL .GT. MXCOLS))      &
     &         GOTO 9200
            IZCOLS = MIN(MXCOLS,MAX(2*IZCOLS,NODEINFO(ND)%KCOL))
         ENDIF
!
         NMTX = 1
         IF (IPER .GT. 1) THEN
            IF (MARKOV) THEN
               NMTX = NMTX + 1
            ELSE
               NMTX = NMTX + IPER - 1
            ENDIF
            IF (HAVE_GC) THEN
               NMTX = NMTX + IPER - 1
            ENDIF
         ENDIF
!
         MNABLK = LASTBA + NMTX
         IF (MNABLK .GT. IZABLK) THEN
            IF (.NOT. EXTMEM_ABLK(MNABLK,IERR)) GOTO 9200
         ENDIF
!
         DO I=1,NMTX
            KCOLA(KDATA(ND)+I) = KCOLA(KDATA(IREFND)+I)
            KELMA(KDATA(ND)+I) = KELMA(KDATA(IREFND)+I)
            NELMA(KDATA(ND)+I) = NELMA(KDATA(IREFND)+I)
         END DO
         LASTBA = LASTBA + NMTX
!
         IF (NQPROF .GT. 0) THEN
            NMTX = MIN(IPER,NQPROF)
            MNQBLK = LASTQB + NMTX
            IF (MNQBLK .GT. IZQBLK) THEN
               IF (.NOT. EXTMEM_QBLK(MNQBLK,IERR)) GOTO 9200
            ENDIF
!
            DO I=1,NMTX
               IQOFF(KDATQ(ND)+I) = IQOFF(KDATQ(IREFND)+I)
               LQOFF(KDATQ(ND)+I) = LQOFF(KDATQ(IREFND)+I)
               NELMQ(KDATQ(ND)+I) = NELMQ(KDATQ(IREFND)+I)
            END DO
            LASTQB = LASTQB + NMTX
         ENDIF
      ENDIF
!
!     Other initializations
!
      DROW = DBLANK
      DCOL = DBLANK
      JERR = 0
!
      J = IPER
      K = IREFND
   99 CONTINUE
         KPREF(J) = K
         J = J - 1
         K = NODEINFO(K)%IANCTR
         IF (J .GT. 0) GOTO 99
!
  100 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOSTO,NREC,LENGTH)
         IF (Q1   .EQ. QAST) GOTO 100
         IF (IERR .GT.    0) GOTO 120
         IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 1900
         IF (Q2 .EQ. QM .AND. Q3 .EQ. QK) GOTO 1900
         IF (Q2 .EQ. QC .AND. Q3 .EQ. QP) GOTO 1900
         IF (Q1 .EQ. QBL ) GOTO 300
            WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,      &
     &             DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
            WRITE (IOLOG, 2600)
            GOTO 9999
!
!     Error during read. Probably end of file
!
  120    CONTINUE
            Q1 = QE
            Q2 = QN
            Q3 = QD
            GOTO 1900
!
!     Here we have a regular data line. Identify the location.
!     --------------------------------------------------------
!
  300 CONTINUE
         IF (NECHO1 .GE. 2) WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,      &
     &                      DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
            DCOL = DNAME(1)
            DROW = DNAME(2)
!
  310 CONTINUE
         CALL IDENTIFY(NREC,Q1,Q2,Q3,Q4,DCOL,DROW,ATEMP(1),DNAME(3),      &
     &                 ATEMP(2),ISTAGE,ISTOCH,IBLOCK,IRPOS,MROW,MCOL)
         IF (ISTOCH .EQ.  1) GOTO 700
         IF (ISTOCH .EQ.  2) GOTO 400
         IF (ISTOCH .EQ.  3) GOTO 500
         IF (ISTOCH .EQ.  4) GOTO 610
         IF (ISTOCH .EQ.  5) GOTO 600
         IF (ISTOCH .EQ.  6) GOTO 602
         IF (ISTOCH .EQ.  7) GOTO 800
         IF (ISTOCH .EQ. 11) GOTO 604
         IF (ISTOCH .EQ. 12) GOTO 500
         IF (ISTOCH .EQ. 13) GOTO 500
         IF (ISTOCH .EQ. -1) GOTO 700
         IF (ISTOCH .EQ. -7) GOTO 800
            IF (ISTOCH .NE. 99) THEN
               WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,DNAME(1),      &
     &                   DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
               WRITE (IOLOG, 2900)
               JERR = 1
            ENDIF
            GOTO 890
!
!     Here we have a random cost coefficient. Copy info if necessary.
!
  400          CONTINUE
                  IF (ISTAGE .EQ. IPER) THEN
                     IF (MODE .EQ. 3) THEN
                        NCURR = ISTAGE
                        GOTO 440
                     ENDIF
                     NREF = KPREF(ISTAGE)
                     NCURR = ND
                     IF (NODEINFO(NCURR)%KCOST .NE.      &
     &                    NODEINFO(NREF)%KCOST) GOTO 440
                        NCOEFF = NODEINFO(ISTAGE)%NCOL      &
     &                         - NODEINFO(ISTAGE)%NROW
                        MNCOST = LASTC + NCOEFF
                        IF (MNCOST .GT. IZCOST) THEN
                           IF (.NOT. EXTMEM_COST(MNCOST,IERR)) GOTO 9200
                        ENDIF
!
                        DO J=1,NCOEFF
                           COST(LASTC+J) = COST(NODEINFO(NREF)%KCOST+J)
                        END DO
                        COST(LASTC+IRPOS) = ATEMP(1)
                        NODEINFO(NCURR)%KCOST = LASTC
                        LASTC = LASTC + NCOEFF
                        GOTO 890
!
  440                CONTINUE
                        COST(NODEINFO(NCURR)%KCOST+IRPOS) = ATEMP(1)
                  ELSE
                     WRITE (IOLOG, 2400) NREC, QLINE
                     WRITE (IOLOG, 2500)
                  ENDIF
                  GOTO 890
!
!     Here we have a random RHS
!
  500          CONTINUE
                  IF (ISTAGE .EQ. IPER) THEN
                     IF (MODE .EQ. 3) THEN
                        NCURR = ISTAGE
                        GOTO 520
                     ENDIF
                     NREF = KPREF(ISTAGE)
                     NCURR = ND
                     IF (NODEINFO(NCURR)%KRHS .NE.      &
     &                    NODEINFO(NREF)%KRHS) GOTO 520
                     MNDRHS = NODEINFO(ISTAGE)%NROW + LASTR
                     IF (MNDRHS .GT. IZDRHS) THEN
                        IF (.NOT. EXTMEM_DRHS(MNDRHS,IERR)) GOTO 9200
                     ENDIF
!
                     DO J=1,NODEINFO(ISTAGE)%NROW
                        RHS(LASTR+J) = RHS(NODEINFO(NREF)%KRHS+J)
                     END DO
                     RHS(LASTR+IRPOS) = ATEMP(1)
                     NODEINFO(NCURR)%KRHS = LASTR
                     LASTR = LASTR + NODEINFO(ISTAGE)%NROW
                     GOTO 890
!
  520             CONTINUE
                     RHS(NODEINFO(NCURR)%KRHS+IRPOS) = ATEMP(1)
                  ELSE
                     WRITE (IOLOG, 2400) NREC, QLINE
                     WRITE (IOLOG, 2700)
                  ENDIF
                  GOTO 890
!
!     RANDOM BOUND ON A DECISION VARIABLE.
!     Case 1: Lower bound
!
  600          CONTINUE
                  JL = 1
                  JU = 0
                  TMPL = ATEMP(1)
                  GOTO 650
!
!     Case 2: Upper bound
!
  602          CONTINUE
                  JL = 0
                  JU = 1
                  TMPU = ATEMP(1)
                  GOTO 650
!
!     Case 3: Both upper and lower bound (type 'FX')
!
  604         CONTINUE
                  JL = 1
                  JU = 1
                  TMPL = ATEMP(1)
                  TMPU = ATEMP(1)
                  GOTO 650
!
!     Stochastic range for one of the rows
!
  610          CONTINUE
                  JL = 0
                  JU = 0
                  IVAR = NODEINFO(NREF)%VARPTR
                  IF     (VARTYP(IVAR+IRPOS) .EQ.  1) THEN
                     JU = 1
                     TMPU = DABS(ATEMP(1))
                  ELSEIF (VARTYP(IVAR+IRPOS) .EQ. -1) THEN
                     JL = 1
                     TMPL = -DABS(ATEMP(1))
                  ELSEIF (VARTYP(IVAR+IRPOS) .EQ.  0) THEN
                     IF (ATEMP(1) .GT. 0.D0) THEN
                        JU = 1
                        TMPU = ATEMP(1)
                     ELSE
                        JL = 1
                        TMPL = ATEMP(1)
                     ENDIF
                  ELSE
                     WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,      &
     &                  DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
                     WRITE (IOLOG, 3070)
                  ENDIF
!
!     Store information -- same code for BOUNDS and RANGES
!
  650          CONTINUE
                  IF (ISTAGE .EQ. IPER) THEN
                     IF (MODE .EQ. 3) THEN
                        NCURR = ISTAGE
                        GOTO 670
                     ENDIF
                     NREF = KPREF(ISTAGE)
                     NCURR = ND
                     KBCUR = NODEINFO(NCURR)%KBOUND
                     KBREF = NODEINFO(NREF )%KBOUND
                     IF (KBCUR .NE. KBREF) GOTO 670
                        MNBNDS = LASTBD + NODEINFO(ISTAGE)%NCOL + 1
                        IF (MNBNDS .GT. IZBNDS) THEN
                           IF (.NOT. EXTMEM_BNDS(MNBNDS,IERR)) GOTO 9200
                        ENDIF
!
                        DO J=1,NODEINFO(ISTAGE)%NCOL+1
                           XLB(LASTBD+J) = XLB(KBREF+J)
                           XUB(LASTBD+J) = XUB(KBREF+J)
                        END DO
                        IF (JL .EQ. 1) XLB(LASTBD+IRPOS) = TMPL
                        IF (JU .EQ. 1) XUB(LASTBD+IRPOS) = TMPU
                        NODEINFO(NCURR)%KBOUND = LASTBD
                        LASTBD = LASTBD + NODEINFO(ISTAGE)%NCOL + 1
                        GOTO 680
!
  670                CONTINUE
                        IF (JL .EQ. 1) XLB(KBCUR+IRPOS) = TMPL
                        IF (JU .EQ. 1) XUB(KBCUR+IRPOS) = TMPU
                  ELSE
                     WRITE (IOLOG, 2400) NREC, QLINE
                     WRITE (IOLOG, 2800)
                  ENDIF
  680             CONTINUE
                     IF (IRPOS .GT. NODEINFO(ISTAGE)%NROW) GOTO 100
                     GOTO 890
!
!     Here we have a random coefficient in the A-matrix
!
  700          CONTINUE
                  IF (IGROUP .EQ.  1) IGROUP = 4
                  IF (MULTI  .EQ. -1) MULTI  = 1
                  DO IPP=NPER,1,-1
                     IF (IBLOCK .GT. KDATA(IPP)) GOTO 720
                  ENDDO
!
  720             CONTINUE
                  NREF = KREF(IPP)
                  IF (MODE .EQ. 3) THEN
                     IF (ISTOCH .EQ. 1) THEN
                        A(KELMA(IBLOCK)+IRPOS) = ATEMP(1)
                        GOTO 890
                     ELSE
                        WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,DNAME(1),      &
     &                        DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
                        WRITE (IOLOG, 2200)
                        GOTO 9999
                     ENDIF
                  ENDIF
!
                  JMTX = IBLOCK - KDATA(IPP)
                  IF (MARKOV) THEN
                     NMTX = MIN(IPP,2)
                  ELSE
                     NMTX = IPP
                  ENDIF
                  IF (JMTX .LE. NMTX) THEN
                     STOCHA(IPP,JMTX) = .TRUE.
                  ELSE
                     STOCHA(JMTX-NMTX,IPP) = .TRUE.
                  ENDIF
!
                  NCURR = ND
                  KDREF = KDATA(NREF)
                  KDCUR = KDATA(NCURR)
                  IF (KDCUR .EQ. KDREF) THEN
                     IF (HAVE_GC) NMTX = NMTX + IPP - 1
                     MNABLK = LASTBA + NMTX
                     IF (MNABLK .GT. IZABLK) THEN
                        IF (.NOT. EXTMEM_ABLK(MNABLK,IERR)) GOTO 9200
                     ENDIF
!
                     KDATA(NCURR) = LASTBA
                     DO I=1,NMTX
                        KCOLA(LASTBA+I) = KCOLA(KDREF+I)
                        KELMA(LASTBA+I) = KELMA(KDREF+I)
                        NELMA(LASTBA+I) = NELMA(KDREF+I)
                     END DO
                     KDCUR  = LASTBA
                     LASTBA = LASTBA + NMTX
                  ENDIF
!
                  IAREF = KDREF + JMTX
                  IACUR = KDCUR + JMTX
                  JCPER = IPP + 1 - JMTX
                  KEREF = KELMA(IAREF)
                  KCREF = KCOLA(IAREF)
                  NELMS = NELMA(IAREF)
                  NCOLS = NODEINFO(KREF(JCPER))%NCOL      &
     &                  - NODEINFO(KREF(JCPER))%NROW
!
                  IF (ISTOCH .EQ. 1) THEN
                     IF (KELMA(IACUR) .EQ. KELMA(IAREF)) THEN
!
!    Copy the matrix coefficients
!
                        MNALMN = LASTA + NELMS
                        IF (MNALMN .GT. IZALMN) THEN
                           IF (.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 9200
                        ENDIF
!
                        DO JCOEF=1,NELMS
                            A(LASTA+JCOEF) =  A(KEREF+JCOEF)
                           IA(LASTA+JCOEF) = IA(KEREF+JCOEF)
                        END DO
                        A(LASTA+IRPOS) = ATEMP(1)
                        KELMA(IACUR) = LASTA
                        NELMA(IACUR) = NELMS
                        LASTA = LASTA + NELMS
                        GOTO 890
!
                     ELSE
                        A(KELMA(IACUR)+IRPOS) = ATEMP(1)
                        GOTO 890
                     ENDIF
                  ELSE
!
!     Here the random coefficient was not mentioned in the reference node
!     We may have to copy both matrix coefficients and column pointers
!
                     WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,DNAME(1),      &
     &                     DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
                     WRITE (IOLOG, 2100)
                     IF (KELMA(IACUR) .EQ. KELMA(IAREF)) THEN
!
!     The matrix has not been copied yet; copy matrix and column pointers
!
                        MNALMN = LASTA  + NELMS + 1
                        IF (MNALMN .GT. IZALMN) THEN
                           IF (.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 9200
                        ENDIF
!
                        MNACOL = LASTCA + NCOLS + 1
                        IF (MNACOL .GT. IZACOL) THEN
                           IF (.NOT. EXTMEM_ACOL(MNACOL,IERR)) GOTO 9200
                        ENDIF
!
                        MOVE  = 0
                        DO JC=1,NCOLS
                           LL = LA(KCREF+JC)
                           KK = LA(KCREF+JC+1) - 1
                           LA(LASTCA+JC) = LA(KCREF+JC) + MOVE
                           DO JCOEF=LL,KK
                               A(LASTA+JCOEF+MOVE) =  A(KEREF+JCOEF)
                              IA(LASTA+JCOEF+MOVE) = IA(KEREF+JCOEF)
                           ENDDO
                           IF (JC .EQ. MCOL) MOVE = 1
                        ENDDO
!
                        IRPOS = LA(KCREF+MCOL+1)
                         A(LASTA +IRPOS)   = ATEMP(1)
                        IA(LASTA +IRPOS)   = MROW
                        LA(LASTCA+NCOLS+1) = NELMS + 2
                        NELMA(KDATA(NCURR)+JMTX) = NELMS + 1
                        KELMA(KDATA(NCURR)+JMTX) = LASTA
                        KCOLA(KDATA(NCURR)+JMTX) = LASTCA
                        LASTA  = LASTA  + NELMS + 1
                        LASTCA = LASTCA + NCOLS + 1
                        GOTO 890
                     ELSE
!
!     The matrix has been copied before, so we just need to make room
!     for the new element: we shift the tail by one position
!
                        MNALMN = LASTA + 1
                        IF (MNALMN .GT. IZALMN) THEN
                           IF (.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 9200
                        ENDIF
!
                        IRPOS = KELMA(IACUR) + LA(KCREF+MCOL+1)
                        DO I=LASTA,IRPOS,-1
                            A(I+1) =  A(I)
                           IA(I+1) = IA(I)
                        END DO
                         A(IRPOS) = ATEMP(1)
                        IA(IRPOS) = MROW
                        LASTA = LASTA  + 1
                        NELMA(IACUR) = NELMA(IACUR) + 1
!
!     This node has different zero pattern from the reference node, so the
!     column pointers can't be shared, either. Copy if necessary and adjust.
!
                        MOVE = 0
                        IF (KCOLA(IACUR) .EQ. KCOLA(IAREF)) THEN
                           MNACOL = LASTCA + NCOLS + 1
                           IF (MNACOL .GT. IZACOL) THEN
                              IF (.NOT. EXTMEM_ACOL(MNACOL,IERR))      &
     &                           GOTO 9200
                           ENDIF
!
                           DO JC=1,NCOLS+1
                              LA(LASTCA+JC) = LA(KCREF+JC) + MOVE
                              IF (JC .EQ. MCOL) MOVE = 1
                           ENDDO
                           KCOLA(KDATA(NCURR)+JMTX) = LASTCA
                           LASTCA = LASTCA + NCOLS + 1
                        ELSE
                           KCURR = KCOLA(KDATA(NCURR)+JMTX)
                           DO JC=1,NCOLS+1
                              LA(KCURR+JC) = LA(KCURR+JC) + MOVE
                              IF (JC .EQ. MCOL) MOVE = 1
                           END DO
                        ENDIF
                     ENDIF
                  ENDIF
                  GOTO 890
!
!     Here we have a random coefficient in the Q-matrix
!
  800       CONTINUE
               DO IPP=NPER,1,-1
                  IF (IBLOCK .GT. KDATQ(IPP)) GOTO 820
               ENDDO
!
  820          CONTINUE
               NREF = KREF(IPP)
               IF (MODE .EQ. 3) THEN
                  IF (ISTOCH .EQ. 7) THEN
                     AQMTX(IQOFF(IBLOCK)+IRPOS) = ATEMP(1)
                     GOTO 890
                  ELSE
                     WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,DNAME(1),      &
     &                         DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
                     WRITE (IOLOG, 2200)
                     GOTO 9999
                  ENDIF
               ENDIF
!
               JMTX = IBLOCK - KDATQ(IPP)
!
               NCURR = ND
               KDREF = KDATQ(NREF)
               KDCUR = KDATQ(NCURR)
               IF (KDCUR .EQ. KDREF) THEN
                  NMTX = MIN(IPP,NQPROF)
                  MNQBLK = LASTQB + NMTX
                  IF (MNQBLK .GT. IZQBLK) THEN
                     IF (.NOT. EXTMEM_QBLK(MNQBLK,IERR)) GOTO 9200
                  ENDIF
!
                  KDATQ(NCURR) = LASTQB
                  DO I=1,NMTX
                     IQOFF(LASTQB+I) = IQOFF(KDREF+I)
                     LQOFF(LASTQB+I) = LQOFF(KDREF+I)
                     NELMQ(LASTQB+I) = NELMQ(KDREF+I)
                  END DO
                  KDCUR  = LASTQB
                  LASTQB = LASTQB + NMTX
               ENDIF
!
               IAREF = KDREF + JMTX
               IACUR = KDCUR + JMTX
               JCPER = IPP + 1 - JMTX
               IQREF = IQOFF(IAREF)
               LQREF = LQOFF(IAREF)
               NELMS = NELMQ(IAREF)
               NCOLS = NODEINFO(KREF(JCPER))%NCOL      &
     &               - NODEINFO(KREF(JCPER))%NROW
!
               IF (ISTOCH .EQ. 7) THEN
                  IF (IQOFF(IACUR) .EQ. IQOFF(IAREF)) THEN
!
!     Copy the matrix coefficients
!
                     MNQLMN = LASTQ + NELMS
                     IF (MNQLMN .GT. IZQLMN) THEN
                        IF (.NOT. EXTMEM_QLMN(MNQLMN,IERR)) GOTO 9200
                     ENDIF
!
                     DO JCOEF=1,NELMS
                        AQMTX(LASTQ+JCOEF) = AQMTX(IQREF+JCOEF)
                        IQMTX(LASTQ+JCOEF) = IQMTX(IQREF+JCOEF)
                     END DO
                     AQMTX(LASTQ+IRPOS) = ATEMP(1)
                     IQOFF(IACUR)  = LASTQ
                     NELMQ(IACUR)  = NELMS
                     LASTQ = LASTQ + NELMS
                     GOTO 890
!
                  ELSE
                     AQMTX(IQOFF(IACUR)+IRPOS) = ATEMP(1)
                     GOTO 890
                  ENDIF
               ELSE
!
!     Here the random coefficient was not mentioned in the reference node
!     We may have to copy both matrix coefficients and column pointers
!
                     WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,DNAME(1),      &
     &                     DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
                     WRITE (IOLOG, 2100)
                     IF (IQOFF(IACUR) .EQ. IQOFF(IAREF)) THEN
!
!     The matrix has not been copied yet; copy matrix and column pointers
!
                        MNQLMN = LASTQ + NELMS + 1
                        IF (MNQLMN .GT. IZQLMN) THEN
                           IF (.NOT. EXTMEM_QLMN(MNQLMN,IERR)) GOTO 9200
                        ENDIF
!
                        MNQCOL = LASTQC + NCOLS + 1
                        IF (MNQCOL .GT. IZQCOL) THEN
                           IF (.NOT. EXTMEM_QCOL(MNQCOL,IERR)) GOTO 9200
                        ENDIF
!
                        MOVE  = 0
                        DO JC=1,NCOLS
                           LL = LQMTX(LQREF+JC)
                           KK = LQMTX(LQREF+JC+1) - 1
                           LQMTX(LASTQC+JC) = LQMTX(LQREF+JC) + MOVE
                           DO JCOEF=LL,KK
                              AQMTX(LASTQ+JCOEF+MOVE) =      &
     &                                    AQMTX(IQREF+JCOEF)
                              IQMTX(LASTQ+JCOEF+MOVE) =      &
     &                                    IQMTX(IQREF+JCOEF)
                           ENDDO
                           IF (JC .EQ. MCOL) MOVE = 1
                        ENDDO
!
                        IRPOS = LQMTX(LQREF+MCOL+1)
                        AQMTX(LASTQ +IRPOS)   = ATEMP(1)
                        IQMTX(LASTQ +IRPOS)   = MROW
                        LQMTX(LASTQC+NCOLS+1) = NELMS + 2
                        IQOFF(KDATQ(NCURR)+JMTX) = LASTQ
                        LQOFF(KDATQ(NCURR)+JMTX) = LASTQC
                        NELMQ(KDATQ(NCURR)+JMTX) = NELMS + 1
                        LASTQ  = LASTQ  + NELMS + 1
                        LASTQC = LASTQC + NCOLS + 1
                        GOTO 890
                     ELSE
!
!     The matrix has been copied before, so we just need to make room
!     for the new element: we shift the tail by one position
!
                        MNQLMN = LASTQ + 1
                        IF (MNQLMN .GT. IZQLMN) THEN
                           IF (.NOT. EXTMEM_QLMN(MNQLMN,IERR)) GOTO 9200
                        ENDIF
!
                        IRPOS = IQOFF(IACUR) + LQMTX(KCREF+MCOL+1)
                        DO I=LASTQ,IRPOS,-1
                           AQMTX(I+1) = AQMTX(I)
                           IQMTX(I+1) = IQMTX(I)
                        END DO
                        AQMTX(IRPOS) = ATEMP(1)
                        IQMTX(IRPOS) = MROW
                        LASTQ  = LASTQ  + 1
                        NELMQ(IACUR) = NELMQ(IACUR) + 1
!
!     This node has different zero pattern from the reference node, so the
!     column pointers can't be shared, either. Copy if necessary and adjust.
!
                        MOVE = 0
                        IF (LQOFF(IACUR) .EQ. LQOFF(IAREF)) THEN
                           MNQCOL = LASTQC + NCOLS + 1
                           IF (MNQCOL .GT. IZQCOL) THEN
                              IF (.NOT. EXTMEM_QCOL(MNACOL,IERR))      &
     &                           GOTO 9200
                           ENDIF
!
                           DO JC=1,NCOLS+1
                              LQMTX(LASTQC+JC) = LQMTX(LQREF+JC) + MOVE
                              IF (JC .EQ. MCOL) MOVE = 1
                           ENDDO
                           LQOFF(KDATQ(NCURR)+JMTX) = LASTQC
                           LASTQC = LASTQC + NCOLS + 1
                        ELSE
                           KCURR = LQOFF(KDATQ(NCURR)+JMTX)
                           DO JC=1,NCOLS+1
                              LQMTX(KCURR+JC) = LQMTX(KCURR+JC) + MOVE
                              IF (JC .EQ. MCOL) MOVE = 1
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDIF
                  GOTO 890
!
!     The third name field might contain more information
!
  890             CONTINUE
                     IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 100
                        DROW = DNAME(3)
                        DNAME(3) = DBLANK
                        ATEMP(1) = ATEMP(2)
                        ATEMP(2) = 0.D0
                        GOTO 310
!
!     Next node card or end of file
!
 1900 CONTINUE
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 2000)      &
     &      NREC,Q1,Q2,Q3,Q4,DNAME(1)
!
         IF (Q2 .EQ. QN .AND. Q3 .EQ. QD) NDINFO = 0
         IF (Q2 .EQ. QC .AND. Q3 .EQ. QP) NDINFO = 1
         IF (Q2 .EQ. QM .AND. Q3 .EQ. QK) NDINFO = 2
         IF (IERR .EQ. 1) IERR = 0
         GOTO 999
!
!     Come here if anything went wrong
!
 9200 CONTINUE
         WRITE (IOLOG, 3200)
         GOTO 9999
!
 9999 CONTINUE
         IERR = 2
  999 CONTINUE
         DEALLOCATE(KPREF, STAT=IER)
         IF (IER .NE. 0 .OR. JERR .NE. 0) IERR = 2
         RETURN
!
 2000 FORMAT(I8,4X,4A1,A8,2X,A8,2X,F12.4,3X,A8,2X,F12.4)
 2100 FORMAT(' XXX - WARNING - Nonzero pattern not as in reference',      &
     &       ' node')
 2200 FORMAT(' XXX -  FATAL  - Nonzero pattern not as in core file')
 2400 FORMAT(I8,2X,A80)
 2500 FORMAT(' XXX - WARNING - Cost coefficient not in current node.',      &
     &       ' Data ignored.')
 2600 FORMAT(' XXX -  FATAL  - Illegal header card in STOCH file')
 2700 FORMAT(' XXX - WARNING - RHS coefficient not in current node.',      &
     &       ' Data ignored.')
 2800 FORMAT(' XXX - WARNING - Range or bound not in current node.',      &
     &       ' Data ignored.')
 2900 FORMAT(' XXX -  FATAL  - Illegal type of random element')
 3070 FORMAT(' XXX - WARNING - Illegal row type in stochastic RANGES',      &
     &       ' section. Data ignored.')
 3200 FORMAT(' XXX -  FATAL  - Memory allocation failed in CPNODE.')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INSCEN ( NREC, IERR )
!
!     Routine: INSCEN
!
!     Purpose:
!
!         To input stoch file in SCENARIO format (See Gassmann and Schweitzer)
!
!     Calling sequence:
!
!         CALL INSCEN (NREC, IERR)
!
!         NREC - Number of records processed
!         IERR - Error status
!
!     Routines called: GETMPS, IDENTIFY
!
!     Remarks: The input format is described in more detail in the paper by
!              Gassmann and Schweitzer.
!
!     ---------------------------------------------------
!     Version of 5 October 2002. Written by Gus Gassmann.
!     ---------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME(3), DBLANK, DROW, DCOL, QLINE, DTEMP
      CHARACTER*1 QBLANK(NAMLEN),QTEMP(NAMLEN)
      INTEGER*4   LENGTH(3)
      REAL*8      ATEMP(4)
!
      EQUIVALENCE (QBLANK,DBLANK),(QTEMP,DTEMP)
!
      DATA QBLANK/NAMLEN*' '/
!
!     Allocate temporary arrays
!
      ALLOCATE (KREF(NPER), KPATH(NPER), STAT=IERR)
      IF (IERR .NE. 0) GOTO 9200
!
      IF (IGROUP .EQ. 1) IGROUP = 4
      DROW      = DBLANK
      DCOL      = DBLANK
      NODES     = NPER
      IIPER     = 1
      IERR      = 0
      JERR      = 0
      LASTLF    = 0
      NSCHAR    = 0
      KSCNAM(1) = 0
      NODEINFO(1)%PROB = 1.D0
!
  100 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                     IERR,IOSTO,NREC,LENGTH)
         IF (Q1   .EQ. QAST) GOTO 100
         IF (IERR .GT.    0) GOTO 910
         IF (Q1 .EQ. QE .AND. Q2 .EQ. QN) GOTO 1900
         IF (Q2 .EQ. QS .AND. Q3 .EQ. QC) GOTO 120
         IF (Q1 .EQ. QBL ) GOTO 300
            WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,      &
     &             DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
            WRITE (IOLOG, 2300)
            GOTO 9999
!
!     Here we have a new scenario (code SC)
!
  120    CONTINUE
            IF (NECHO1 .GE. 1)      &
     &         WRITE (IOLOG, 2400) NREC,DNAME(1),DNAME(2),DNAME(3)
            IF (LASTLF .GT. 0) GOTO 140
!
!     The current scenario is scenario 1; store the information
!
               DO I=1,NPER
                  NODEINFO(I)%PROB = ATEMP(1)
                  KREF(I)  = I
                  KPATH(I) = I
               END DO
!
               LASTLF = LASTLF + 1
               MNCHSC = NSCHAR+LENGTH(1)
               IF (MNCHSC .GT. IZCHSC) THEN
                  IF (.NOT. EXTMEM_CHSC(MNCHSC,IERR)) GOTO 9200
               ENDIF
               IF (LASTLF .GE. IZSCEN) THEN
                  IF (.NOT. EXTMEM_SCEN(LASTLF+1,IERR)) GOTO 9200
               ENDIF
!
               DTEMP = DNAME(1)
               DO I=1,LENGTH(1)
                  QSCNAM(NSCHAR+I) = QTEMP(I)
               END DO
               NSCHAR = NSCHAR + LENGTH(1)
               KSCNAM(LASTLF+1) = NSCHAR
               LEAFND(LASTLF) = NPER
               GOTO 100
!
!     This is not scenario 1; find the scenario it branches from
!
  140       CONTINUE
               DO I=1,LASTLF
                  LEN = KSCNAM(I+1) - KSCNAM(I)
                  KSC = KSCNAM(I)
                  DO J=1,LEN
                     QTEMP(J) = QSCNAM(KSC+J)
                  END DO
                  IF (DNAME(2) .EQ. DTEMP) GOTO 200
               END DO
!
               WRITE (IOLOG, 2000) NREC, Q1,Q2,Q3,Q4,DNAME(1),      &
     &            DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
               WRITE (IOLOG, 2600)
               GOTO 9999
!
!     Found the branching scenario. Store the name.
!
  200       CONTINUE
               LASTLF = LASTLF + 1
               IF (LASTLF .GE. IZSCEN) THEN
                  IF (.NOT. EXTMEM_SCEN(LASTLF+1,IERR)) GOTO 9200
               ENDIF
               MNCHSC = NSCHAR+LENGTH(1)
               IF (MNCHSC .GT. IZCHSC) THEN
                  IF (.NOT. EXTMEM_CHSC(MNCHSC,IERR)) GOTO 9200
               ENDIF
!
               DTEMP  = DNAME(1)
               DO J=1,LENGTH(1)
                  QSCNAM(NSCHAR+J) = QTEMP(J)
               END DO
               NSCHAR = NSCHAR + LENGTH(1)
               KSCNAM(LASTLF+1) = NSCHAR
               LASTN = LEAFND(I)
               IP = NPER
!
!     Identify the branching period (by counting back)
!
  220       CONTINUE
               KREF(IP) = LASTN
               IF (DTIME(IP) .EQ. DNAME(3)) GOTO 230
               IF (IP .EQ. 1) GOTO 9040
                  IP = IP - 1
                  LASTN = NODEINFO(LASTN)%IANCTR
                  GOTO 220
!
!     Fix the probabilities of all the ancestors
!
  230       CONTINUE
               IIPER = IP
               IBRO1 = LASTN
  240       CONTINUE
               LASTN = NODEINFO(LASTN)%IANCTR
               IP = IP - 1
               IF (LASTN .EQ. 0) GOTO 250
                  NODEINFO(LASTN)%PROB = NODEINFO(LASTN)%PROB + ATEMP(1)
                  KPATH(IP) = LASTN
                  KREF(IP)  = LASTN
                  GOTO 240
!
!     Set pointers for all nodes on this scenario
!
  250       CONTINUE
               MNNODE = NODES + NPER - IIPER + 1
               IF (MNNODE .GT. IZNODE) THEN
                  IF (.NOT. EXTMEM_NODE(MNNODE,IERR)) GOTO 9200
               ENDIF
!
               NODEINFO(NODES+1)%IANCTR = KREF(IIPER-1)
               DO WHILE ( NODEINFO(IBRO1)%IBROTH .NE. 0 )
                  IBRO1 = NODEINFO(IBRO1)%IBROTH
               END DO
!
               NODEINFO(IBRO1)%IBROTH = NODES + 1
               DO IP=IIPER,NPER
                  NODES = NODES + 1
                  LASTN = KREF(IP)
                  IF (NECHO1 .GE. 2) WRITE (IOLOG, 2200) NODES
                  KPATH(IP) = NODES
                  NODEINFO(NODES)%IBROTH = 0
                  IF (IP .EQ. NPER) THEN
                     NODEINFO(NODES)%IDESC = 0
                  ELSE
                     NODEINFO(NODES  )%IDESC  = NODES + 1
                     NODEINFO(NODES+1)%IANCTR = NODES
                  ENDIF
                  KDATA(NODES) = KDATA(LASTN)
                  NODEINFO(NODES)%PROB   = ATEMP(1)
                  NODEINFO(NODES)%KROW   = NODEINFO(NODES-1)%KROW +      &
     &                                     NODEINFO(NODES-1)%NROW
                  NODEINFO(NODES)%KCOL   = NODEINFO(NODES-1)%KCOL +      &
     &                                     NODEINFO(NODES-1)%NCOL + 1
                  NODEINFO(NODES)%KCOST  = NODEINFO(LASTN)%KCOST
                  NODEINFO(NODES)%KRHS   = NODEINFO(LASTN)%KRHS
                  NODEINFO(NODES)%KBOUND = NODEINFO(LASTN)%KBOUND
                  NODEINFO(NODES)%KNAMES = NODEINFO(LASTN)%KNAMES
                  NODEINFO(NODES)%VARPTR = NODEINFO(LASTN)%VARPTR
                  NODEINFO(NODES)%NROW   = NODEINFO(LASTN)%NROW
                  NODEINFO(NODES)%NCOL   = NODEINFO(LASTN)%NCOL
                  IF (ALLOCATED(KDATQ)) KDATQ(NODES) = KDATQ(LASTN)
                  IF (LCHANCE .GT. 0) THEN
                     KCHANCE(NODES) = KCHANCE(LASTN)
                     NCHANCE(NODES) = NCHANCE(LASTN)
                  ENDIF
                  IF (LICC    .GT. 0) THEN
                     KICC(NODES) = KICC(LASTN)
                     NICC(NODES) = NICC(LASTN)
                  ENDIF
                  IF (ALLOCATED(KPLQ)) KPLQ(NODES) = KPLQ(LASTN)
!
                  MNROWS = NODEINFO(NODES)%NROW + NODEINFO(NODES)%KROW
                  IF (MNROWS .GT. IZROWS) THEN
                     IF (.NOT. RESIZE .OR. (IZROWS .GE. MXROWS)      &
     &                                .OR. (MNROWS .GT. MXROWS))      &
     &                  GOTO 9200
                     IZROWS = MNROWS
                  ENDIF
!
                  MNCOLS = NODEINFO(NODES)%NCOL + NODEINFO(NODES)%KCOL
                  IF (MNCOLS .GT. IZCOLS) THEN
                     IF (.NOT. RESIZE .OR. (IZCOLS .GE. MXCOLS)      &
     &                                .OR. (MNCOLS .GT. MXCOLS))      &
     &                  GOTO 9200
                     IZCOLS = MNCOLS
                  ENDIF
!
               END DO
               LEAFND(LASTLF) = NODES
               GOTO 100
!
!     Here we have a regular data line. Identify the location.
!     --------------------------------------------------------
!
  300 CONTINUE
         IF (NECHO1 .GE. 3) WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,      &
     &                      DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
            DCOL = DNAME(1)
            DROW = DNAME(2)
!
  310 CONTINUE
         CALL IDENTIFY(NREC,Q1,Q2,Q3,Q4,DCOL,DROW,ATEMP(1),DNAME(3),      &
     &                 ATEMP(2),ISTAGE,ISTOCH,IBLOCK,IRPOS,MROW,MCOL)
         IF (ISTOCH .EQ.  1) GOTO 700
         IF (ISTOCH .EQ.  2) GOTO 400
         IF (ISTOCH .EQ.  3) GOTO 500
         IF (ISTOCH .EQ.  4) GOTO 610
         IF (ISTOCH .EQ.  5) GOTO 600
         IF (ISTOCH .EQ.  6) GOTO 602
         IF (ISTOCH .EQ.  7) GOTO 800
         IF (ISTOCH .EQ. 11) GOTO 604
         IF (ISTOCH .EQ. 12) GOTO 500
         IF (ISTOCH .EQ. 13) GOTO 500
         IF (ISTOCH .EQ. -1) GOTO 700
         IF (ISTOCH .EQ. -7) GOTO 800
         IF (ISTOCH .EQ.  8) GOTO 900
         IF (ISTOCH .EQ.  9) GOTO 1000
         IF (ISTOCH .EQ. 10) GOTO 1100
            IF (ISTOCH .NE. 99) THEN
               WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,DNAME(1),      &
     &                   DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
               WRITE (IOLOG, 2500)
               JERR = 1
            ENDIF
            GOTO 890
!
!     Here we have a random cost coefficient. Copy info if necessary.
!
  400    CONTINUE
            IF (LASTLF .GT. 1) THEN
               NREF = KREF(IBLOCK)
               NCURR = NODES + IBLOCK - NPER
!
!     Check whether NCURR is indeed on the current path.
!
               DO JP=1,NPER
                  IF (NCURR .EQ. KPATH(JP)) GOTO 401
               END DO
               GOTO 9100
!
  401       CONTINUE
               IF (NODEINFO(NCURR)%KCOST .EQ. NODEINFO(NREF)%KCOST)      &
     &               THEN
                  NCOEFF = NODEINFO(NREF)%NCOL - NODEINFO(NREF)%NROW
                  MNCOST = LASTC + NCOEFF
                  IF (MNCOST .GT. IZCOST) THEN
                     IF (.NOT. EXTMEM_COST(MNCOST,IERR)) GOTO 9200
                  ENDIF
!
                  KCREF = NODEINFO(NREF)%KCOST
                  DO J=1,NCOEFF
                     COST(LASTC+J) = COST(KCREF+J)
                  END DO
                  COST(LASTC+IRPOS) = ATEMP(1)
                  NODEINFO(NCURR)%KCOST = LASTC
                  LASTC = LASTC + NCOEFF
               ELSE
                  COST(NODEINFO(NCURR)%KCOST+IRPOS) = ATEMP(1)
               ENDIF
            ELSE
                  NCURR = IBLOCK
                  COST(NODEINFO(NCURR)%KCOST+IRPOS) = ATEMP(1)
            ENDIF
            GOTO 890
!
!     Here we have a random RHS
!
  500    CONTINUE
            IF (LASTLF .GT. 1) THEN
               NREF = KREF(IBLOCK)
               NCURR = NODES + IBLOCK - NPER
!
!     Check whether NCURR is indeed on the current path.
!
               DO JP=1,NPER
                  IF (NCURR .EQ. KPATH(JP)) GOTO 501
               END DO
               GOTO 9100
!
  501       CONTINUE
               IF (NODEINFO(NCURR)%KRHS .EQ. NODEINFO(NREF)%KRHS)      &
     &               THEN
                  NCOEFF = NODEINFO(NREF)%NROW
                  MNDRHS = LASTR + NCOEFF
                  IF (MNDRHS .GT. IZDRHS) THEN
                     IF (.NOT. EXTMEM_DRHS(MNDRHS,IERR)) GOTO 9200
                  ENDIF
!
                  KRREF = NODEINFO(NREF)%KRHS
                  DO J=1,NCOEFF
                     RHS(LASTR+J) = RHS(KRREF+J)
                  END DO
                  RHS(LASTR+IRPOS) = ATEMP(1)
                  NODEINFO(NCURR)%KRHS = LASTR
                  LASTR = LASTR + NODEINFO(IBLOCK)%NROW
               ELSE
                  RHS(NODEINFO(NCURR)%KRHS+IRPOS) = ATEMP(1)
               ENDIF
            ELSE
                  NCURR = IBLOCK
                  RHS(NODEINFO(NCURR)%KRHS+IRPOS) = ATEMP(1)
            ENDIF
            GOTO 890
!
!     RANDOM BOUND ON A DECISION VARIABLE.
!     Case 1: Lower bound
!
  600    CONTINUE
            JL = 1
            JU = 0
            TMPL = ATEMP(1)
            GOTO 650
!
!     Case 2: Upper bound
!
  602    CONTINUE
            JL = 0
            JU = 1
            TMPU = ATEMP(1)
            GOTO 650
!
!     Case 3: Both upper and lower bound (type 'FX')
!
  604    CONTINUE
            JL = 1
            JU = 1
            TMPL = ATEMP(1)
            TMPU = ATEMP(1)
            GOTO 650
!
!     Stochastic range for one of the rows. Treatment depends on row type.
!
  610    CONTINUE
            JL = 0
            JU = 0
            IT = VARTYP(NODEINFO(IIPER)%VARPTR + IRPOS)
            IF (IT .EQ. -1) THEN
               JL = 1
               TMPL = -DABS(ATEMP(1))
            ELSEIF (IT .EQ.  1) THEN
               JU = 1
               TMPU = DABS(ATEMP(1))
            ELSEIF (IT .EQ. 0) THEN
               IF (ATEMP(1) .GT. 0.D0) THEN
                  JU = 1
                  TMPU = ATEMP(1)
               ELSE
                  JL = 1
                  TMPL = ATEMP(1)
               ENDIF
            ELSE
               GOTO 9070
            ENDIF
!
!     Store information -- same code for BOUNDS and RANGES
!
  650       CONTINUE
               IF (LASTLF .GT. 1) THEN
                  NREF = KREF(IBLOCK)
                  NCURR = NODES + IBLOCK - NPER
!
!     Check whether NCURR is indeed on the current path.
!
               DO JP=1,NPER
                  IF (NCURR .EQ. KPATH(JP)) GOTO 651
               END DO
               GOTO 9100
!
  651       CONTINUE
                  IF (NODEINFO(NCURR)%KBOUND .EQ. NODEINFO(NREF)%KBOUND)      &
     &                  THEN
                     NCOEFF = NODEINFO(NREF)%NCOL + 1
                     MNBNDS = LASTBD + NCOEFF
                     IF (MNBNDS .GT. IZBNDS) THEN
                        IF (.NOT. EXTMEM_BNDS(MNBNDS,IERR)) GOTO 9200
                     ENDIF
!
                     KBREF = NODEINFO(NREF)%KBOUND
                     DO J=1,NCOEFF
                        XLB(LASTBD+J) = XLB(KBREF+J)
                        XUB(LASTBD+J) = XUB(KBREF+J)
                     END DO
                     IF (JL .EQ. 1) XLB(LASTBD+IRPOS) = TMPL
                     IF (JU .EQ. 1) XUB(LASTBD+IRPOS) = TMPU
                     NODEINFO(NCURR)%KBOUND = LASTBD
                     LASTBD = LASTBD + NCOEFF
!
                  ELSE
                     KBCUR = NODEINFO(NCURR)%KBOUND
                     IF (JL .EQ. 1) XLB(KBCUR+IRPOS) = TMPL
                     IF (JU .EQ. 1) XUB(KBCUR+IRPOS) = TMPU
                  ENDIF
               ELSE
                     NCURR = IBLOCK
                     KBCUR = NODEINFO(NCURR)%KBOUND
                     IF (JL .EQ. 1) XLB(KBCUR+IRPOS) = TMPL
                     IF (JU .EQ. 1) XUB(KBCUR+IRPOS) = TMPU
               ENDIF
               IF (IRPOS .GT. NODEINFO(IBLOCK)%NROW) GOTO 100
                  GOTO 890
!
!     Here we have a random coefficient in the A-matrix
!
  700          CONTINUE
                  IF (IGROUP .EQ.  1) IGROUP = 4
                  IF (MULTI  .EQ. -1) MULTI  = 1
                  DO IPP=NPER,1,-1
                     IF (IBLOCK .GT. KDATA(IPP)) GOTO 720
                  ENDDO
!
  720             CONTINUE
                  NREF = KREF(IPP)
                  IF (LASTLF .EQ. 1) THEN
                     IF (ISTOCH .EQ. 1) THEN
                        A(KELMA(IBLOCK)+IRPOS) = ATEMP(1)
                        GOTO 890
                     ELSE
                        WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,DNAME(1),      &
     &                        DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
                        WRITE (IOLOG, 2800)
                        GOTO 9999
                     ENDIF
                  ENDIF
!
                  JMTX = IBLOCK - KDATA(IPP)
                  IF (MARKOV) THEN
                     NMTX = MIN(IPP,2)
                  ELSE
                     NMTX = IPP
                  ENDIF
                  IF (JMTX .LE. NMTX) THEN
                     STOCHA(IPP,JMTX) = .TRUE.
                     JCPER = IPP + 1 - JMTX
                  ELSE
                     STOCHA(JMTX-NMTX,IPP) = .TRUE.
                     JCPER = 5
                  ENDIF
!
                  NCURR = NODES + IPP - NPER
!
!     Check whether NCURR is indeed on the current path.
!
                  DO JP=1,NPER
                     IF (NCURR .EQ. KPATH(JP)) GOTO 721
                  END DO
                  GOTO 9100
!
  721          CONTINUE
                  KDREF = KDATA(NREF)
                  KDCUR = KDATA(NCURR)
                  IF (KDCUR .EQ. KDREF) THEN
                     IF (HAVE_GC) NMTX = NMTX + IPP - 1
                     MNABLK = LASTBA + NMTX
                     IF (MNABLK .GT. IZABLK) THEN
                        IF (.NOT. EXTMEM_ABLK(MNABLK,IERR)) GOTO 9200
                     ENDIF
!
                     KDATA(NCURR) = LASTBA
                     DO I=1,NMTX
                        KCOLA(LASTBA+I) = KCOLA(KDREF+I)
                        KELMA(LASTBA+I) = KELMA(KDREF+I)
                        NELMA(LASTBA+I) = NELMA(KDREF+I)
                     END DO
                     KDCUR  = LASTBA
                     LASTBA = LASTBA + NMTX
                  ENDIF
!
                  IAREF = KDREF + JMTX
                  IACUR = KDCUR + JMTX
                  KEREF = KELMA(IAREF)
                  KCREF = KCOLA(IAREF)
                  NELMS = NELMA(IAREF)
                  NCOLS = NODEINFO(KREF(JCPER))%NCOL      &
     &                  - NODEINFO(KREF(JCPER))%NROW
!
                  IF (ISTOCH .EQ. 1) THEN
                     IF (KELMA(IACUR) .EQ. KELMA(IAREF)) THEN
!
!     Copy the matrix coefficients
!
                        MNALMN = LASTA + NELMS
                        IF (MNALMN .GT. IZALMN) THEN
                           IF (.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 9200
                        ENDIF
                        DO JCOEF=1,NELMS
                            A(LASTA+JCOEF) =  A(KEREF+JCOEF)
                           IA(LASTA+JCOEF) = IA(KEREF+JCOEF)
                        END DO
                        A(LASTA+IRPOS) = ATEMP(1)
                        KELMA(IACUR) = LASTA
                        NELMA(IACUR) = NELMS
                        LASTA = LASTA + NELMS
                        GOTO 890
!
                     ELSE
                        A(KELMA(IACUR)+IRPOS) = ATEMP(1)
                        GOTO 890
                     ENDIF
                  ELSE
!
!     Here the random coefficient was not mentioned in the reference node
!     We may have to copy both matrix coefficients and column pointers
!
                     WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,DNAME(1),      &
     &                     DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
                     WRITE (IOLOG, 2700)
                     IF (KELMA(IACUR) .EQ. KELMA(IAREF)) THEN
!
!     The matrix has not been copied yet; copy matrix and column pointers
!
                        MNALMN = LASTA + NELMS + 1
                        IF (MNALMN .GT. IZALMN) THEN
                           IF (.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 9200
                        ENDIF
!
                        MNACOL = LASTCA + NCOLS + 1
                        IF (MNACOL .GT. IZACOL) THEN
                           IF (.NOT. EXTMEM_ACOL(MNACOL,IERR)) GOTO 9200
                        ENDIF
!
                        MOVE  = 0
                        DO JC=1,NCOLS
                           LL = LA(KCREF+JC)
                           KK = LA(KCREF+JC+1) - 1
                           LA(LASTCA+JC) = LA(KCREF+JC) + MOVE
                           DO JCOEF=LL,KK
                               A(LASTA+JCOEF+MOVE) =  A(KEREF+JCOEF)
                              IA(LASTA+JCOEF+MOVE) = IA(KEREF+JCOEF)
                           ENDDO
                           IF (JC .EQ. MCOL) MOVE = 1
                        ENDDO
!
                        IRPOS = LA(KCREF+MCOL+1)
                         A(LASTA +IRPOS)   = ATEMP(1)
                        IA(LASTA +IRPOS)   = MROW
                        LA(LASTCA+NCOLS+1) = NELMS + 2
                        NELMA(KDATA(NCURR)+JMTX) = NELMS + 1
                        KELMA(KDATA(NCURR)+JMTX) = LASTA
                        KCOLA(KDATA(NCURR)+JMTX) = LASTCA
                        LASTA  = LASTA  + NELMS + 1
                        LASTCA = LASTCA + NCOLS + 1
                        GOTO 890
                     ELSE
!
!     The matrix has been copied before, so we just need to make room
!     for the new element: we shift the tail by one position
!
                        MNALMN = LASTA + 1
                        IF (MNALMN .GT. IZALMN) THEN
                           IF (.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 9200
                        ENDIF
!
                        IRPOS = KELMA(IACUR) + LA(KCREF+MCOL+1)
                        DO I=LASTA,IRPOS,-1
                            A(I+1) =  A(I)
                           IA(I+1) = IA(I)
                        ENDDO
                         A(IRPOS) = ATEMP(1)
                        IA(IRPOS) = MROW
                        LASTA = LASTA  + 1
                        NELMA(IACUR) = NELMA(IACUR) + 1
!
!     This node has different zero pattern from the reference node, so the
!     column pointers can't be shared, either. Copy if necessary and adjust.
!
                        MOVE = 0
                        IF (KCOLA(IACUR) .EQ. KCOLA(IAREF)) THEN
                           MNACOL = LASTCA + NCOLS + 1
                           IF (MNACOL .GT. IZACOL) THEN
                              IF (.NOT. EXTMEM_ACOL(MNACOL,IERR))      &
     &                           GOTO 9200
                           ENDIF
!
                           DO JC=1,NCOLS+1
                              LA(LASTCA+JC) = LA(KCREF+JC) + MOVE
                              IF (JC .EQ. MCOL) MOVE = 1
                           ENDDO
                           KCOLA(KDATA(NCURR)+JMTX) = LASTCA
                           LASTCA = LASTCA + NCOLS + 1
                        ELSE
                           KCURR = KCOLA(KDATA(NCURR)+JMTX)
                           DO JC=1,NCOLS+1
                              LA(KCURR+JC) = LA(KCURR+JC) + MOVE
                              IF (JC .EQ. MCOL) MOVE = 1
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDIF
                  GOTO 890
!
!     Here we have a random coefficient in the Q-matrix
!
  800       CONTINUE
               DO IPP=NPER,1,-1
                  IF (IBLOCK .GT. KDATQ(IPP)) GOTO 820
               ENDDO
!
  820          CONTINUE
               NREF = KREF(IPP)
               IF (LASTLF .EQ. 1) THEN
                  IF (ISTOCH .EQ. 7) THEN
                     AQMTX(IQOFF(IBLOCK)+IRPOS) = ATEMP(1)
                     GOTO 890
                  ELSE
                     WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,DNAME(1),      &
     &                         DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
                     WRITE (IOLOG, 2800)
                     GOTO 9999
                  ENDIF
               ENDIF
!
               JMTX  = IBLOCK - KDATQ(IPP)
               NCURR = NODES + IPP - NPER
!
!     Check whether NCURR is indeed on the current path.
!
               DO JP=1,NPER
                  IF (NCURR .EQ. KPATH(JP)) GOTO 821
               END DO
               GOTO 9100
!
  821       CONTINUE
               KDREF = KDATQ(NREF)
               KDCUR = KDATQ(NCURR)
               IF (KDCUR .EQ. KDREF) THEN
                  NMTX = MIN(IPP,NQPROF)
                  MNQBLK = LASTQB + NMTX
                  IF (MNQBLK .GT. IZQBLK) THEN
                     IF (.NOT. EXTMEM_QBLK(MNQBLK,IERR)) GOTO 9200
                  ENDIF
!
                  KDATQ(NCURR) = LASTQB
                  DO I=1,NMTX
                     IQOFF(LASTQB+I) = IQOFF(KDREF+I)
                     LQOFF(LASTQB+I) = LQOFF(KDREF+I)
                     NELMQ(LASTQB+I) = NELMQ(KDREF+I)
                  ENDDO
                  KDCUR  = LASTQB
                  LASTQB = LASTQB + NMTX
               ENDIF
!
               IAREF = KDREF + JMTX
               IACUR = KDCUR + JMTX
               JCPER = IPP + 1 - JMTX
               IQREF = IQOFF(IAREF)
               LQREF = LQOFF(IAREF)
               NELMS = NELMQ(IAREF)
               NCOLS = NODEINFO(KREF(JCPER))%NCOL      &
     &               - NODEINFO(KREF(JCPER))%NROW
!
               IF (ISTOCH .EQ. 7) THEN
                  IF (IQOFF(IACUR) .EQ. IQOFF(IAREF)) THEN
!
!     Copy the matrix coefficients
!
                     MNQLMN = LASTQ + NELMS
                     IF (MNQLMN .GT. IZQLMN) THEN
                        IF (.NOT. EXTMEM_QLMN(MNQLMN,IERR)) GOTO 9200
                     ENDIF
!
                     DO JCOEF=1,NELMS
                        AQMTX(LASTQ+JCOEF) = AQMTX(IQREF+JCOEF)
                        IQMTX(LASTQ+JCOEF) = IQMTX(IQREF+JCOEF)
                     END DO
                     AQMTX(LASTQ+IRPOS) = ATEMP(1)
                     IQOFF(IACUR) = LASTQ
                     NELMQ(IACUR) = NELMS
                     LASTQ = LASTQ + NELMS
                     GOTO 890
!
                  ELSE
                     AQMTX(IQOFF(KDATQ(NCURR)+JMTX)+IRPOS) = ATEMP(1)
                     GOTO 890
                  ENDIF
               ELSE
!
!     Here the random coefficient was not mentioned in the reference node
!     We may have to copy both matrix coefficients and column pointers
!
                     WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,DNAME(1),      &
     &                     DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
                     WRITE (IOLOG, 2700)
                     IF (IQOFF(IACUR) .EQ. IQOFF(IAREF)) THEN
!
!     The matrix has not been copied yet; copy matrix and column pointers
!
                        MNQLMN = LASTQ + NELMS + 1
                        IF (MNQLMN .GT. IZQLMN) THEN
                           IF (.NOT. EXTMEM_QLMN(MNQLMN,IERR)) GOTO 9200
                        ENDIF
!
                        MNQCOL = LASTQC + NCOLS + 1
                        IF (MNQCOL .GT. IZQCOL) THEN
                           IF (.NOT. EXTMEM_QCOL(MNQCOL,IERR)) GOTO 9200
                        ENDIF
!
                        MOVE  = 0
                        DO JC=1,NCOLS
                           LL = LQMTX(LQREF+JC)
                           KK = LQMTX(LQREF+JC+1) - 1
                           LQMTX(LASTQC+JC) = LQMTX(LQREF+JC) + MOVE
                           DO JCOEF=LL,KK
                              AQMTX(LASTQ+JCOEF+MOVE) =      &
     &                                       AQMTX(IQREF+JCOEF)
                              IQMTX(LASTQ+JCOEF+MOVE) =      &
     &                                       IQMTX(IQREF+JCOEF)
                           ENDDO
                           IF (JC .EQ. MCOL)  MOVE = 1
                        ENDDO
!
                        IRPOS = LQMTX(LQREF+MCOL+1)
                        AQMTX(LASTQ +IRPOS)   = ATEMP(1)
                        IQMTX(LASTQ +IRPOS)   = MROW
                        LQMTX(LASTQC+NCOLS+1) = NELMS + 2
                        IQOFF(KDATQ(NCURR)+JMTX) = LASTQ
                        LQOFF(KDATQ(NCURR)+JMTX) = LASTQC
                        NELMQ(KDATQ(NCURR)+JMTX) = NELMS + 1
                        LASTQ  = LASTQ  + NELMS + 1
                        LASTQC = LASTQC + NCOLS + 1
                        GOTO 890
                     ELSE
!
!     The matrix has been copied before, so we just need to make room
!     for the new element: we shift the tail by one position
!
                        MNQLMN = LASTQ + 1
                        IF (MNQLMN .GT. IZQLMN) THEN
                           IF (.NOT. EXTMEM_QLMN(MNQLMN,IERR)) GOTO 9200
                        ENDIF
!
                        IRPOS = IQOFF(IACUR) + LQMTX(KCREF+MCOL+1)
                        DO I=LASTQ,IRPOS,-1
                           AQMTX(I+1) = AQMTX(I)
                           IQMTX(I+1) = IQMTX(I)
                        ENDDO
                        AQMTX(IRPOS) = ATEMP(1)
                        IQMTX(IRPOS) = MROW
                        LASTQ  = LASTQ  + 1
                        NELMQ(IACUR) = NELMQ(IACUR) + 1
!
!     This node has different zero pattern from the reference node, so the
!     column pointers can't be shared, either. Copy if necessary and adjust.
!
                        MOVE = 0
                        IF (LQOFF(IACUR) .EQ. LQOFF(IAREF)) THEN
                           MNQCOL = LASTQC + NCOLS + 1
                           IF (MNQCOL .GT. IZQCOL) THEN
                              IF (.NOT. EXTMEM_QCOL(MNQCOL,IERR))      &
     &                           GOTO 9200
                           ENDIF
!
                           DO JC=1,NCOLS+1
                              LQMTX(LASTQC+JC) = LQMTX(LQREF+JC) + MOVE
                              IF (JC .EQ. MCOL) MOVE = 1
                           ENDDO
                           LQOFF(KDATQ(NCURR)+JMTX) = LASTQC
                           LASTQC = LASTQC + NCOLS + 1
                        ELSE
                           KCURR = LQOFF(KDATQ(NCURR)+JMTX)
                           DO JC=1,NCOLS+1
                              LQMTX(KCURR+JC) = LQMTX(KCURR+JC) + MOVE
                              IF (JC .EQ. MCOL) MOVE = 1
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDIF
                  GOTO 890
!
!     Random RHS in a probabilistic constraint. Copy info if necessary.
!
  900    CONTINUE
            IF (LASTLF .GT. 1) THEN
               NREF = KREF(IBLOCK)
               NCURR = NODES + IBLOCK - NPER
!
!     Check whether NCURR is indeed on the current path.
!
               DO JP=1,NPER
                  IF (NCURR .EQ. KPATH(JP)) GOTO 901
               END DO
               GOTO 9100
!
  901       CONTINUE
               IF (KCHANCE(NCURR) .EQ. KCHANCE(NREF)) THEN
                  NCC = NCHANCE(NREF)
                  MNCCON = LCHANCE + NCC
                  IF (MNCCON .GT. IZCCON) THEN
                     IF (.NOT. EXTMEM_CCON(MNCCON,IERR)) GOTO 9200
                  ENDIF
!
                  KCREF = KCHANCE(NREF)
                  DO J=1,NCC
                     XCHANCE(LCHANCE+J) = XCHANCE(KCREF+J)
                  END DO
                  XCHANCE(LCHANCE+IRPOS)%PROB = ATEMP(1)
                  KCHANCE(NCURR) = LCHANCE
                  LCHANCE = LCHANCE + NCC
               ELSE
                  XCHANCE(KCHANCE(NCURR)+IRPOS)%PROB = ATEMP(1)
               ENDIF
            ELSE
                  NCURR = IBLOCK
                  XCHANCE(KCHANCE(NCURR)+IRPOS)%PROB = ATEMP(1)
            ENDIF
            GOTO 890
!
!     Random RHS in an integrated chance constraint. Copy info if necessary.
!
 1000    CONTINUE
            IF (LASTLF .GT. 1) THEN
               NREF = KREF(IBLOCK)
               NCURR = NODES + IBLOCK - NPER
!
!     Check whether NCURR is indeed on the current path.
!
               DO JP=1,NPER
                  IF (NCURR .EQ. KPATH(JP)) GOTO 1001
               END DO
               GOTO 9100
!
 1001       CONTINUE
               IF (KICC(NCURR) .EQ. KICC(NREF)) THEN
                  NCC = NICC(NREF)
                  MNNICC = LICC + NCC
                  IF (MNNICC .GT. IZNICC) THEN
                     IF (.NOT. EXTMEM_NICC(MNNICC,IERR)) GOTO 9200
                  ENDIF
!
                  KCREF = KICC(NREF)
                  DO J=1,NCC
                     XICC(LICC+J) = XICC(KCREF+J)
                  END DO
                  XICC(LICC+IRPOS)%LIMIT = ATEMP(1)
                  KICC(NCURR) = LICC
                  LICC = LICC + NCC
               ELSE
                  XICC(KICC(NCURR)+IRPOS)%LIMIT = ATEMP(1)
               ENDIF
            ELSE
                  NCURR = IBLOCK
                  XICC(KICC(NCURR)+IRPOS)%LIMIT = ATEMP(1)
            ENDIF
            GOTO 890
!
!     Curvature coefficient for piecewise linear-quadratic penalty
!
 1100    CONTINUE
            IF (LASTLF .GT. 1) THEN
               NREF = KREF(IBLOCK)
               NCURR = NODES + IBLOCK - NPER
!
!     Check whether NCURR is indeed on the current path.
!
               DO JP=1,NPER
                  IF (NCURR .EQ. KPATH(JP)) GOTO 1101
               END DO
               GOTO 9100
!
 1101       CONTINUE
               IF (KPLQ(NCURR) .EQ. KPLQ(NREF)) THEN
                  NCOEFF = NODEINFO(NREF)%NCOL - NODEINFO(NREF)%NROW
                  MNCPLQ = LASTPLQ + NCOEFF
                  IF (MNCPLQ .GT. IZCPLQ) THEN
                     IF (.NOT. EXTMEM_CPLQ(MNCPLQ,IERR)) GOTO 9200
                  ENDIF
!
                  KCREF = KPLQ(NREF)
                  DO J=1,NCOEFF
                     CPLQ(LASTC+J) = CPLQ(KCREF+J)
                  END DO
                  CPLQ(LASTPLQ+IRPOS) = ATEMP(1)
                  KPLQ(NCURR) = LASTPLQ
                  LASTPLQ = LASTPLQ + NCOEFF
               ELSE
                  CPLQ(KPLQ(NCURR)+IRPOS) = ATEMP(1)
               ENDIF
!
            ELSE
                  NCURR = IBLOCK
                  CPLQ(KPLQ(NCURR)+IRPOS) = ATEMP(1)
            ENDIF
            GOTO 890
!
!     The third name field might contain more information
!
  890 CONTINUE
         IF (DABS(ATEMP(2)) .LE. ZTOLIN) GOTO 100
            DROW = DNAME(3)
            DNAME(3) = DBLANK
            ATEMP(1) = ATEMP(2)
            ATEMP(2) = 0.D0
            GOTO 310
!
!     End of stoch file
!
 1900 CONTINUE
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 2000)      &
     &      NREC,Q1,Q2,Q3,Q4,DNAME(1)
         GOTO 999
!
  910 CONTINUE
         WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,      &
     &          DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 2100)
         GOTO 999
!
!     Come here if anything went wrong
!
 9040 CONTINUE
         WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,      &
     &          DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3040)
         GOTO 9999
!
 9070 CONTINUE
         WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,      &
     &          DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3070)
         GOTO 9999
!
 9100 CONTINUE
         WRITE (IOLOG, 2000) NREC,Q1,Q2,Q3,Q4,      &
     &          DNAME(1),DNAME(2),ATEMP(1),DNAME(3),ATEMP(2)
         WRITE (IOLOG, 3100)
         GOTO 9999
!
 9200 CONTINUE
         WRITE (IOLOG, 3200)
         GOTO 9999
!
!
 9999 CONTINUE
         IERR = 2
  999 CONTINUE
       DEALLOCATE (KREF, KPATH, STAT=IER)
         IF (IER .NE. 0 .OR. JERR .NE. 0) IERR = 2
         IF (JERR .NE. 0) IERR = 2
         RETURN
!
 2000 FORMAT(I8,4X,4A1,A8,2X,A8,2X,F12.4,3X,A8,2X,F12.4)
 2100 FORMAT(' XXX - WARNING - Missing ENDATA card')
 2200 FORMAT(' Creating node',I8)
 2300 FORMAT(' XXX -  FATAL  - Illegal header card in STOCH file')
 2400 FORMAT(I8,4X,' Found scenario ',A8,' branching from ',A8,      &
     &       ' in stage ',A8)
 2500 FORMAT(' XXX -  FATAL  - Illegal type of random element')
 2600 FORMAT(' XXX -  FATAL  - Misspecified branch in decision tree')
 2700 FORMAT(' XXX - WARNING - Nonzero pattern not as in reference',      &
     &       ' node')
 2800 FORMAT(' XXX -  FATAL  - Nonzero pattern not as in core file')
 2900 FORMAT(I8,4X,A80)
 3040 FORMAT(' XXX -  FATAL  - Period could not be found')
 3070 FORMAT(' XXX -  FATAL  - Illegal row type in stochastic RANGES',      &
     &       ' section')
 3100 FORMAT(' XXX -  FATAL  - Stochastic element violates time',      &
     &       ' structure.')
 3200 FORMAT(' XXX -  FATAL  - Insufficient memory in routine INSCEN.')
 3800 FORMAT(' XXX - WARNING - Penalty parameters in simple recourse',      &
     &       ' section not convex. Data ignored.')
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE INCONT(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
!
!     Routine:   INCONT
!
!     Purpose:
!
!         This subroutine reads the distributions of continuously
!         distributed stochastic elements and stores them internally for
!         future reference. The user will have to find his/her own way of
!         dealing with them in a solver routine.
!
!     Calling sequence:
!
!         CALL INCONT(Q1,Q2,Q3,Q4, DNAME, ATEMP, NREC, IERR)
!
!            Q1     First character of the code field on an MPS-style data card
!            Q2     Second character of the code field
!            Q3     Third character of the code field
!            Q4     Fourth character of the code field
!            DNAME  Name fields on the data card
!            ATEMP  Numeric fields on the data card
!            NREC   Number of records processed
!            IERR   Error status
!
!     Routines called: GETMPS, GETPAR, IDENTIFY, IDROW, GENCOL
!
!     --------------------------------------------------------------
!        This version dated 25 August 1998. Written by Gus Gassmann.
!     --------------------------------------------------------------
!                                                 X
!     It is important to distinguish between STOCHASTIC ELEMENTS
!     which describe locations in the extended constraint matrix
!     (including RHS, cost, bounds, etc.) and RANDOM VARIABLES
!     which give information about their distributions.
!
!     In the LINTRAN option, the correspondence between random variables
!     and stochastic elements is not one-to-one; the linkage between
!     random variables and stochastic elements is done in the
!     DISTRIBUTION MATRIX or D-matrix for short. It has block-diagonal
!     structure and is stored in sparse matrix form in arrays DMTX,
!     IDMTX and LDMTX. The random variables are the COLUMNS of this
!     matrix and the stochastic elements are the rows. Note that some
!     distributions define multi-dimensional random vectors. Such
!     random vectors occupy more than one column in the D-matrix.
!     Moreover, linear (affine) transformations usually involve more
!     than one random variable, including a degenerate one to describe
!     the translation component.
!
!
!     Data structures:
!     ----------------
!
!     Several data items must be stored for each random variable:
!
!     1. The distribution
!        In order to allow extensions and user defined routines, the
!        distributions are stored as strings in array SDIST.
!
!     2. The dimension
!     3. The number and location of integer parameters
!     4. The number and location of real parameters
!     5. The number and identity of parameters referenced by location
!     6. The number of the stage in which the random variable becomes known
!
!        These items are stored in array MDIST:
!        MDIST(8*IRV-7) gives the dimension of the random variable IRV
!        MDIST(8*IRV-6) gives the offset for the first integer parameter
!                       (parameter values are in array IRNDPAR)
!        MDIST(8*IRV-5) gives the number of integer parameters
!        MDIST(8*IRV-4) gives the offset for the first real parameter
!                       (parameter values are in array RANDPAR)
!        MDIST(8*IRV-3) gives the number of real parameters
!        MDIST(8*IRV-2) gives the offset for the first parameter referenced
!                       by location (parameter values are in array LRNDPAR)
!        MDIST(8*IRV-1) gives the number of location parameters
!        MDIST(8*IRV  ) gives the number of the stage
!
!
!     The location and process mode of the stochastic elements are stored
!     in array KSTOCH in groups of four. The first entry gives the type
!     of random element, the second its block (or node), the third entry
!     describes the location within the relevant array, the fourth specifies
!     whether the values are to be treated as replacing the core file
!     information or if they are to be treated as additive or multiplicative
!     perturbations.
!
!     Types of random elements are as follows (can be extended as needed):
!
!        KSTOCH(4*NSTELM-3) = 1 for stochastic entry in an A-matrix block
!        KSTOCH(4*NSTELM-3) = 2 for stochastic cost
!        KSTOCH(4*NSTELM-3) = 3 for stochastic RHS
!        KSTOCH(4*NSTELM-3) = 4 for stochastic range
!        KSTOCH(4*NSTELM-3) = 5 for stochastic lower bound
!        KSTOCH(4*NSTELM-3) = 6 for stochastic upper bound
!        KSTOCH(4*NSTELM-3) = 7 for stochastic entry in a Q-matrix block
!        KSTOCH(4*NSTELM-3) = 11 for stochastic upper/lower bounds ('FX' type)
!        KSTOCH(4*NSTELM-3) = 12 for stochastic demand (network option)
!        KSTOCH(4*NSTELM-3) = 13 for stochastic supply (network option)
!
!     For stochastic A-matrix elements, the second entry is as in the
!     offset array KDATA, that is, depending on the value of MARKOV
!     and the presence or absence of global (linking) variables:
!
!      a) no global variables
!        period           staircase problems      triangular problems
!          1                 1                       1
!          2                 3   2                   3   2
!          3                     5   4               6   5   4
!          4                         7  6           10   9   8   7
!         ...                          ...                ...

!      b) with global variables
!        period           staircase problems      triangular problems
!          1                 1   4   7  11           1   4   8  14
!          2                 3   2   8  12           3   2   9  15
!          3                     6   5  13           7   6   5  16
!          4                        10  14          13  12  11  10
!         ...                          ...                ...
!
!
!     The location of each data item is recorded as the relative address
!     within the relevant array.
!
!     For stochastic Q-matrix elements, IBLOCK is as in the offset array
!     KDATQ, that is, depending on the value of NQPROF:
!
!        period    NQPROF = 1       NQPROF = 2        NQPROF = 3      ...
!          1       1                1                 1
!          2           2            3   2             3   2
!          3               3            5   4         6   5   4
!          4                   4            7   6         9   8   7
!         ...                ...             ...            ...
!
!
!     The processing mode is in location KSTOCH(4*NSTELM):
!
!        KSTOCH(4*NSTELM) = 1 if the core file info is to be replaced
!        KSTOCH(4*NSTELM) = 2 if the random information is additive
!        KSTOCH(4*NSTELM) = 3 if the information is multiplicative
!
!
!     Example:
!     ========
!
!     For a particular random element we have
!
!        KSTOCH(4*N + 1) = 1
!        KSTOCH(4*N + 2) = 5
!        KSTOCH(4*N + 3) = 8
!        KSTOCH(4*N + 4) = 2
!
!     The first entry identifies this as a A-random matrix element, so the core
!     file values are recorded in the array A (along with row indices in IA
!     and column indices in LA). The second entry specifies the block within
!     the matrix: the subdiagonal block linking periods 2 and 3 (independent
!     of whether the problem is markovian or not). Offsets to this block are
!     then found in KELMA(5), and the particular element is in KELMA(5) + 8.
!     The fourth entry indicates that the stochastic value must be ADDED to the
!     core file value.
!
! ============================================================================
!
!     Aliases
!     =======
!
!     Aliases can arise in two situations: In a DISTRIB segment, random
!     variables (and vectors) are named for future use, and in an MVNORM
!     segment aliases are assigned to stochastic elements for later
!     processing of the correlation/covariance matrix. Since there is
!     a one-to-one correspondence between stochastic elements and random
!     variables in the latter case, we can treat both the same way.
!
!     For univariate random variables, the alias itself is stored in a
!     string array ALIAS, indexed by NALIAS. The number of the random
!     variable is recorded in a corresponding array IALIAS. For random
!     vectors, the vector as well as each component get a name recorded
!     in the ALIAS array. The IALIAS array contains the number of the
!     random variable in the location corresponding to the vector and
!     an offset (which is a negative number) describing where the random
!     variable's number is to be found.
!
!
!     Example:
!     ========
!
!      ALIAS: ...  BLOCK1  VAR1    VAR2    VAR3     INDEP1  ...
!     IALIAS: ...  11      -1      -2      -3       12
!
!     If a correlation matrix makes reference to the component INDEP1,
!     the corresponding entry in the IALIAS array is 12, indicating that
!     we are dealing with a univariate random variable, whose number is 12.
!     Further information is then found in array MDIST (in locations 89-96).
!     If the reference is to the component VAR2, then the negative entry
!     in the IALIAS array indicates that we have to look two places to the
!     left for more. The entry corresponding to BLOCK1 then directs us to
!     locations 81 to 88 in the MDIST array.
!
! ============================================================================
!
!     In LINTR segments, reference can be made to random variables set up
!     earlier in a DISTRIB segment. Such variables can be reused in different
!     LINTR segments, and therefore we must make a copy for each occurrence.
!     MXCOPY gives a limit to the number of distributions that can be handled
!     in a single LINTR. (There should ordinarily be no more than one
!     distribution, but it's best to give the user a little leeway...)
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME(*),DN1,DN2,DTYPE,BLNAME,DBLANK,QLINE,AKA
      CHARACTER*8  DCONST, DISCR,  DLIN1,  DLIN2,  DMODE,  DMNORM,      &
     &             DADD,   DMULT1, DMULT2, DREP
!
      CHARACTER*2 DCODE
      CHARACTER*1 QBLANK(NAMLEN),QN2(NAMLEN)
      REAL*8      ATEMP(*)
      INTEGER     LENGTH(3)
!
      EQUIVALENCE (DBLANK,QBLANK),(DN2,QN2)
!
      DATA DCONST/'CONSTANT'/, DISCR /'DISCRETE'/, DLIN1 /'LINTR   '/,      &
     &     DLIN2 /'LINTRAN '/, DMNORM/'MVNORMAL'/, DADD  /'ADD     '/,      &
     &     DMULT1/'MULTIPLY'/, DREP  /'REPLACE '/, DMULT2/'MULT    '/
!
      DATA QBLANK/NAMLEN*' '/
!
!     Allocate storage arrays for random variables and stochastic elements
!
      ALLOCATE (DMTX(IZDLMN),   IDMTX(IZDLMN),  LDMTX(IZDCOL),      &
     &          RANDPAR(IZRPAR),IRNDPAR(IZIPAR),LRNDPAR(3*IZLPAR),      &
     &          MDIST(8*IZRAND),SDIST(IZRAND),  KSTOCH(4*IZSTOC),      &
     &          STAT=IERR)
      IF (IERR .NE. 0) GOTO 9200
!
!     Allocate temporary arrays
!
      ALLOCATE (ALIAS(IZALIAS), IALIAS(IZALIAS), MARKER(IZCOPY),      &
     &          COPIED(IZCOPY), KPATH(NPER), KREF(NPER), STAT=IERR)
      IF (IERR .NE. 0) GOTO 9200
!
!     Initialize various pointers
!
      DN1      = DBLANK
      DN2      = DBLANK
      BLNAME   = DBLANK
      NDCOLM   = 0
      NDELEM   = 0
      NSTELM   = 0
      NALIAS   = 0
      NPREV    = 0
      NRVAR    = 0
      LIPAR    = 0
      LRPAR    = 0
      LLPAR    = 0
      LDMTX(1) = 1
      DO I=1,NPER
         KREF(I) = I
         KPATH(I) = I
      END DO
!
!     The first data record is read in the top-level input routine INPUT
!     and is passed to this routine as part of the parameter list. For this
!     reason the READ statement is at the END of the routine.
!
  100 CONTINUE
      IF (Q1 .EQ. QI  .AND.  Q2   .EQ. QN) GOTO 110
      IF (Q1 .EQ. QB  .AND.  Q2   .EQ. QL) GOTO 120
      IF (Q1 .EQ. QD  .AND.  Q2   .EQ. QI) GOTO 130
      IF (Q1 .EQ. QE  .AND.  Q2   .EQ. QN) GOTO 900
!
      IF (Q1 .EQ. QBL .AND. ITYPE .EQ.  1) GOTO 200
      IF (Q1 .EQ. QBL .AND. ITYPE .EQ.  2) GOTO 300
      IF (Q1 .EQ. QBL .AND. ITYPE .EQ.  3) GOTO 600
      IF (Q1 .EQ. QBL .AND. ITYPE .EQ.  4) GOTO 650
      IF (Q1 .EQ. QBL .AND. ITYPE .EQ.  5) GOTO 599
!
      IF (Q1 .EQ. QAST) GOTO 890
!
         WRITE (IOLOG, 1900) NREC,QLINE
         WRITE (IOLOG, 1000)
         IERR = 2
         GOTO 999
!
!     ================================================================
!     Start of a new segment of independent distributions (INDEP)
!     ================================================================
!
  110 CONTINUE
         LISTPAR = 0
         ITYPE = 1
         DTYPE = DNAME(2)
         DMODE = DNAME(3)
         IF (DTYPE .EQ. DISCR) THEN
            IF (STFILE .LT. 3) STFILE = 3
         ELSE
             STFILE = 4
         ENDIF
         IF (DMODE .EQ. DADD) THEN
            METHOD = 2
         ELSEIF (DMODE .EQ. DMULT1 .OR. DMODE .EQ. DMULT2) THEN
            METHOD = 3
         ELSE
            METHOD = 1
            IF (DMODE .NE. DBLANK .AND. DMODE .NE. DREP) THEN
               WRITE (IOLOG, 1900) NREC,QLINE
               WRITE (IOLOG, 2100)
            ENDIF
         ENDIF
         GOTO 890
!
!     Data card for independent distribution (ITYPE = 1).
!     ----------------------------------------------------
!     The card could have a 'PR' code, or it could be a repeat realization of
!     a discrete distribution.
!
  200 CONTINUE
         IF (Q2       .EQ. QP  .AND. Q3       .EQ. QR)    GOTO 280
         IF (DNAME(1) .EQ. DN1 .AND. DNAME(2) .EQ. DN2      &
     &                         .AND. DTYPE    .EQ. DISCR) GOTO 270
         IF (LISTPAR  .EQ. 1)   GOTO 285
!
!     Identify the stochastic element.
!
!     ISTOCH gives the type of random element and is used for diagnostics:
!     If ISTOCH = 0, a row or column name was not found. Flush to next header.
!     If ISTOCH < 0, we make room for the coefficient within the A or Q matrix.
!     If ISTOCH = 99, a recoverable error occurred. Ignore this record
!
            CALL IDENTIFY(NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2),ATEMP(1),      &
     &                    DNAME(3),ATEMP(2),ISTAGE,ISTOCH,IBLOCK,IRPOS,      &
     &                    IROW,ICOL)
!
            IF (ISTOCH .EQ.  0) GOTO 850
            IF (ISTOCH .EQ. 99) GOTO 890
               IF (ISTOCH .LT. 0) THEN
!
!     Enlarge the A-matrix
!
                  IF (ISTOCH .EQ. -1) THEN
                     IF (LASTA .GE. IZALMN) THEN
                        IF (.NOT. EXTMEM_ALMN(LASTA+1,IERR)) GOTO 9200
                     ENDIF
!
                     DO I=LASTA,IRPOS,-1
                         A(I+1) =  A(I)
                        IA(I+1) = IA(I)
                     END DO
                      A(IRPOS) = 0.D0
                     IA(IRPOS) = IROW
                     LASTA = LASTA + 1
!
                     ISTOCH = 1
                     DO I=1,LASTBA
                        IF (I .NE. IBLOCK .AND. KELMA(I) .GE. IRPOS-1)      &
     &                               KELMA(I) = KELMA(I) + 1
                     END DO
!
                     DO I=1,NPER
                        IF (IBLOCK .GT. KDATA(I) .AND.      &
     &                      IBLOCK .LE. KDATA(I+1) ) GOTO 215
                     END DO
!
!     This should not happen!
!
                     WRITE (IOLOG, 3500)
                     IERR = 2
                     GOTO 999
!
  215             CONTINUE
                     IF (MARKOV) THEN
                        NSUB = MIN(I,2)
                     ELSE
                        NSUB = I
                     ENDIF
                     IBL = IBLOCK - KDATA(I)
                     IF (IBL .LE. NSUB) THEN
                        IP = I - IBL + 1
                     ELSE
                        IP = I
                     ENDIF
                     LCOL = NODEINFO(IP)%NCOL - NODEINFO(IP)%NROW + 1
                     LL = KCOLA(IBLOCK)
                     DO I=ICOL+1,LCOL
                        LA(LL+I) = LA(LL+I) + 1
                     END DO
                     IRPOS = IRPOS - KELMA(IBLOCK)
                     NELMA(IBLOCK) = NELMA(IBLOCK) + 1
!
!     Update the stochastic elements found so far
!
                     DO I=1,NSTELM
                        IF (KSTOCH(4*I-3) .EQ.   1    .AND.      &
     &                      KSTOCH(4*I-2) .EQ. IBLOCK .AND.      &
     &                      KSTOCH(4*I-1) .GE. IRPOS )      &
     &                      KSTOCH(4*I-1) = KSTOCH(4*I-1) + 1
                     END DO
!
!     Enlarge the Q-matrix
!
                  ELSEIF (ISTOCH .EQ. -7) THEN
                     IF (LASTQ .GE. IZQLMN) THEN
                        IF (.NOT. EXTMEM_QLMN(LASTQ+1,IERR)) GOTO 9200
                     ENDIF
!
                     DO I=LASTQ,IRPOS,-1
                        AQMTX(I+1) = AQMTX(I)
                        IQMTX(I+1) = IQMTX(I)
                     END DO
                     AQMTX(IRPOS) = 0.D0
                     IQMTX(IRPOS) = IROW
                     LASTQ = LASTQ + 1
!
                     ISTOCH = 7
                     DO I=1,LASTQB
                        IF (I .NE. IBLOCK .AND. IQOFF(I) .GE. IRPOS-1)      &
     &                               IQOFF(I) = IQOFF(I) + 1
                     END DO
!
                     DO I=1,NPER
                        IBL = IBLOCK - KDATQ(I)
                        IF (IBL .GT. 0 .AND. IBL .LE. NQPROF ) GOTO 216
                     END DO
                     WRITE (IOLOG, 3700)
                     IERR = 2
                     GOTO 999
!
  216             CONTINUE
                     IP = I - IBL + 1
                     LCOL = NODEINFO(IP)%NCOL - NODEINFO(IP)%NROW + 1
                     LQB = LQOFF(IBLOCK)
                     DO I=ICOL+1,LCOL
                        LQMTX(LQB+I) = LQMTX(LQB+I) + 1
                     END DO
                     IRPOS = IRPOS - IQOFF(IBLOCK)
                     NELMQ(IBLOCK) = NELMQ(IBLOCK) + 1
!
!     Update the stochastic elements found so far
!
                     DO I=1,NSTELM
                        IF (KSTOCH(4*I-3) .EQ.   7    .AND.      &
     &                      KSTOCH(4*I-2) .EQ. IBLOCK .AND.      &
     &                      KSTOCH(4*I-1) .GE. IRPOS )      &
     &                      KSTOCH(4*I-1) = KSTOCH(4*I-1) + 1
                     END DO
!
!     Something happened to ISTOCH. Stop the program.
!
                  ELSE
                     WRITE (IOLOG, 3700)
                     IERR = 2
                     GOTO 999
                  ENDIF
               ENDIF
!
!     Store the information.
!
               IF (NRVAR .GE. IZRAND) THEN
                  IF (.NOT. EXTMEM_RAND(NRVAR+1,IERR)) GOTO 9200
               ENDIF
!
               NRVAR = NRVAR  + 1
               SDIST(NRVAR) = DTYPE
               MDIST(8*NRVAR-7) = 1
               MDIST(8*NRVAR-6) = LIPAR
               MDIST(8*NRVAR-5) = 0
               MDIST(8*NRVAR-4) = LRPAR
               MDIST(8*NRVAR-3) = 2
               MDIST(8*NRVAR-2) = LLPAR
               MDIST(8*NRVAR-1) = 0
               MDIST(8*NRVAR  ) = ISTAGE
!
               IF (LRPAR+2 .GT. IZRPAR) THEN
                  IF (.NOT. EXTMEM_RPAR(LRPAR+2,IERR)) GOTO 9200
               ENDIF
               IF (DTYPE .EQ. DISCR) THEN
                  RANDPAR(LRPAR+1) = ATEMP(2)
                  RANDPAR(LRPAR+2) = ATEMP(1)
               ELSE
                  RANDPAR(LRPAR+1) = ATEMP(1)
                  RANDPAR(LRPAR+2) = ATEMP(2)
               ENDIF
               NINTPAR = 0
               NREAPAR = 2
               NLOCPAR = 0
               LRPAR = LRPAR + 2
!
               IF (NSTELM .GE. IZSTOC) THEN
                  IF (.NOT. EXTMEM_STOC(NSTELM+1,IERR)) GOTO 9200
               ENDIF
!
               NSTELM = NSTELM + 1
               KSTOCH(4*NSTELM-3) = ISTOCH
               KSTOCH(4*NSTELM-2) = IBLOCK
               KSTOCH(4*NSTELM-1) = IRPOS
               KSTOCH(4*NSTELM  ) = METHOD
!
               MNDCOL = NDCOLM + 2
               IF (MNDCOL .GT. IZDCOL) THEN
                  IF (.NOT. EXTMEM_DCOL(MNDCOL,IERR)) GOTO 9200
               ENDIF
!
               IF (NDELEM .GE. IZDLMN) THEN
                  IF (.NOT. EXTMEM_DLMN(NDELEM+1,IERR)) GOTO 9200
               ENDIF
!
               NDCOLM = NDCOLM + 1
               NDELEM = NDELEM + 1
                DMTX(NDELEM)   = 1.D0
               IDMTX(NDELEM)   = NSTELM
               LDMTX(NDCOLM+1) = NDELEM + 1
               GOTO 890
!
!     Repeat realization of a discrete random variable.
!
  270 CONTINUE
         IF (LRPAR+2 .GT. IZRPAR) THEN
            IF (.NOT. EXTMEM_RPAR(LRPAR+2,IERR)) GOTO 9200
         ENDIF
!
         RANDPAR(LRPAR + 1) = ATEMP(2)
         RANDPAR(LRPAR + 2) = ATEMP(1)
         NREAPAR = NREAPAR + 2
         MDIST(8*NRVAR-3) = MDIST(8*NRVAR-3) + 2
         LRPAR = LRPAR + 2
         GOTO 890
!
!     'PR' card
!
  280 CONTINUE
         LISTPAR = 1
         NREAPAR = 0
         NINTPAR = 0
         NLOCPAR = 0
         LRPAR = LRPAR - 2
         MDIST(8*NRVAR-3) = 0
!
!     Processing the parameter list
!
  285 CONTINUE
         CALL GETPAR(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,XPAR,IPAR,MTYPE,      &
     &               IBLOCK,IERR,IOSTO,NREC)
         IF (Q1 .EQ. QAST) GOTO 285
         IF (IERR .GT. 0) THEN
            IF (IERR .EQ. 3) THEN
               LISTPAR = 0
               GOTO 100
            ENDIF
         ENDIF
         IF (MTYPE .EQ. 1) THEN
            NINTPAR = NINTPAR + 1
            IF (LIPAR+2 .GT. IZIPAR) THEN
               IF (.NOT. EXTMEM_IPAR(LIPAR+2,IERR)) GOTO 9200
            ENDIF
!
            LIPAR = LIPAR + 1
            IRNDPAR(LIPAR) = IPAR
            MDIST(8*NRVAR-5) = MDIST(8*NRVAR-5) + 1
         ELSEIF (MTYPE .EQ. 2) THEN
            NREAPAR = NREAPAR + 1
            IF (LRPAR+2 .GT. IZRPAR) THEN
               IF (.NOT. EXTMEM_RPAR(LRPAR+2,IERR)) GOTO 9200
            ENDIF
!
            LRPAR = LRPAR + 1
            RANDPAR(LRPAR) = XPAR
            MDIST(8*NRVAR-3) = MDIST(8*NRVAR-3) + 1
         ELSEIF (MTYPE .LT. 0) THEN
            NLOCPAR = NLOCPAR + 1
            IF (LLPAR+2 .GT. IZLPAR) THEN
               IF (.NOT. EXTMEM_LPAR(LLPAR+2,IERR)) GOTO 9200
            ENDIF
!
            LLPAR = LLPAR + 1
            LRNDPAR(3*LLPAR-2) = -MTYPE
            LRNDPAR(3*LLPAR-1) =  IBLOCK
            LRNDPAR(3*LLPAR  ) =  IPAR
            MDIST(8*NRVAR-1) = MDIST(8*NRVAR-1) + 1
         ENDIF
         GOTO 285
!
!     ======================================================================
!     Start of a new segment of joint distributions (BLOCK)
!
!     There are four different BLOCK types (DISCRETE, MVNORMAL, LINTR, SUB)
!     Different types allow different codes to appear in the code field,
!     but all types must start with a 'BL' type card to fix the locations.
!     We keep track of the last code card read in variable DCODE.
!
!     For the record, here are the legal values for the code:
!
!        DISCRETE:           BL
!        MVNORMAL:           BL, CV, CR
!        LINTR:              BL, RV
!        Other (subroutine): BL, PR
!     ======================================================================
!
  120 CONTINUE
         LISTPAR = 0
         DCODE = 'XX'
         ITYPE = 2
         DTYPE = DNAME(2)
         IF (DTYPE .EQ. DLIN2) DTYPE = DLIN1
         IF (DTYPE .EQ. DISCR) THEN
            IF (STFILE .LT. 3) STFILE = 3
         ELSEIF (DTYPE .NE. DLIN1) THEN
            STFILE = 4
         ENDIF
         DMODE = DNAME(3)
         IF (DMODE .EQ. DADD) THEN
            METHOD = 2
         ELSEIF (DMODE .EQ. DMULT1 .OR. DMODE .EQ. DMULT2) THEN
            METHOD = 3
         ELSE
            METHOD = 1
            IF (DMODE .NE. DBLANK .AND. DMODE .NE. DREP) THEN
               WRITE (IOLOG, 1900) NREC,QLINE
               WRITE (IOLOG, 2100)
            ENDIF
         ENDIF
         GOTO 890
!
!     Data card for BLOCK segment. First check the code field:
!     The first non-comment card must be a 'BL' card.
!
  300 CONTINUE
         IF (Q2 .EQ. QB .AND. Q3 .EQ. QL) GOTO 400
         IF (DCODE      .EQ.        'XX') GOTO 970
         IF (Q2 .EQ. QR .AND. Q3 .EQ. QV) GOTO 450
         IF (Q2 .EQ. QC .AND. Q3 .EQ. QR) GOTO 500
         IF (Q2 .EQ. QC .AND. Q3 .EQ. QV) GOTO 550
         IF (Q2 .EQ. QP .AND. Q3 .EQ. QR) GOTO 570
            GOTO 350
!
!     A 'BL' card starts the definition of locations. Keep the period
!     ---------------------------------------------------------------
!
  400    CONTINUE
            DCODE = 'BL'
            IF (BLNAME .EQ. DNAME(1)) GOTO 440
               BLNAME = DNAME(1)
               DO ISTAGE=1,NPER
                  IF (DNAME(2) .EQ. DTIME(ISTAGE)) GOTO 430
               END DO
!
!     The stage information is missing. If we have a two-stage problem
!     or a one-period chance-constrained problem, we can recover.
!
               IF (NPER .GT. 2) THEN
                  WRITE (IOLOG, 1900) NREC,QLINE
                  WRITE (IOLOG, 1200) DNAME(2)
                  IERR = 2
                  GOTO 999
               ELSE
                  WRITE (IOLOG, 1900) NREC, QLINE
                  WRITE (IOLOG, 1250) NPER
                  ISTAGE = NPER
               ENDIF
!
!     Keep the information
!
  430       CONTINUE
               IF (NRVAR .GE. IZRAND) THEN
                  IF (.NOT. EXTMEM_RAND(NRVAR+1,IERR)) GOTO 9200
               ENDIF
!
               NRVAR  = NRVAR  + 1
               MDIST(8*NRVAR-6) = LIPAR
               MDIST(8*NRVAR-5) = 0
               MDIST(8*NRVAR-4) = LRPAR
               MDIST(8*NRVAR-2) = LLPAR
               MDIST(8*NRVAR-1) = 0
               MDIST(8*NRVAR  ) = ISTAGE
               NINTPAR = 0
               NLOCPAR = 0
               IF (DTYPE .EQ. DISCR) THEN
                  SDIST(NRVAR) = DTYPE
                  MDIST(8*NRVAR-7) = 0
                  NREAPAR = 1
                  IF (LRPAR .GE. IZRPAR) THEN
                     IF (.NOT. EXTMEM_RPAR(LRPAR+1,IERR)) GOTO 9200
                  ENDIF
!
                  LRPAR = LRPAR + 1
                  RANDPAR(LRPAR) = ATEMP(1)
                  MDIST(8*NRVAR-3) = 1
                  NREAL = 1
               ELSEIF (DTYPE .EQ. DMNORM) THEN
                  SDIST(NRVAR) = DTYPE
                  MDIST(8*NRVAR-7) = 0
                  MDIST(8*NRVAR-3) = 0
                  NREAPAR = 0
                  NINTPAR = 1
                  MDIST(8*NRVAR-5) = 1
                  IF (LIPAR .GE. IZIPAR) THEN
                     IF (.NOT. EXTMEM_IPAR(LIPAR+1,IERR)) GOTO 9200
                  ENDIF
!
                  LIPAR = LIPAR + 1
                  IRNDPAR(LIPAR) = 0
                  IF (NALIAS .GE. IZALIAS) THEN
                     IF (.NOT. EXTMEM_ALIAS(NALIAS+1,IERR)) GOTO 9200
                  ENDIF
!
                  NALIAS = NALIAS + 1
                   ALIAS(NALIAS) = BLNAME
                  IALIAS(NALIAS) = NRVAR + 1
               ELSEIF (DTYPE .EQ. DLIN1) THEN
!
!     The translation component c of the affine transformation v = HU + c
!     is stored as a separate column in the D-matrix, which means that
!     it must be treated like a (degenerate) random variable. We use
!     the reserved name 'CONSTANT' for this purpose.
!
                  SDIST(NRVAR) = DCONST
                  MDIST(8*NRVAR-7) = 1
                  NREAPAR = 1
                  IF (LRPAR+1 .GT. IZRPAR) THEN
                     IF (.NOT. EXTMEM_RPAR(LRPAR+1,IERR)) GOTO 9200
                  ENDIF
!
                  LRPAR = LRPAR + 1
                  RANDPAR(LRPAR) = 1.D0
                  MDIST(8*NRVAR-3) = 1
                  KRV     = 0
                  NPREV   = NSTELM
                  NDCOLM  = NDCOLM + 1
                  NCOPIED = 0
                  LDMTX(NDCOLM+1) = NDELEM + 1
               ELSE
                  SDIST(NRVAR) = DTYPE
                  MDIST(8*NRVAR-7) = 0
                  MDIST(8*NRVAR-3) = 0
                  NREAPAR = 0
               ENDIF
               GOTO 890
!
!     This block has the same name as the last. Only legal for DISCRETE
!     distributions. Record the probability and copy base values.
!
  440       CONTINUE
               IF (DTYPE .NE. DISCR) GOTO 950
                  NREAL = NREAL + 1
                  MNRPAR = LRPAR + NREAPAR
                  IF (MNRPAR .GT. IZRPAR) THEN
                     IF (.NOT. EXTMEM_RPAR(MNRPAR,IERR)) GOTO 9200
                  ENDIF
!
                  RANDPAR(LRPAR+1) = ATEMP(1)
                  DO I=2,NREAPAR
                     RANDPAR(LRPAR + I) = RANDPAR(MDIST(8*NRVAR-4) + I)
                  END DO
                  MDIST(8*NRVAR-3) = MDIST(8*NRVAR-3) + NREAPAR
                  LRPAR = LRPAR + NREAPAR
                  GOTO 890
!
!     An 'RV' card sets up a new auxiliary random variable. LINTR ONLY.
!     -----------------------------------------------------------------
!
  450       CONTINUE
            IF (DTYPE .NE. DLIN1) GOTO 960
               DCODE = 'RV'
               IF (DNAME(2) .NE. DBLANK) THEN
!
!     If the second name field contains a name, then this name is taken to
!     describe the (univariate) distribution of this random variable.
!
                  IF (DNAME(2) .EQ. DISCR) THEN
                     IF (STFILE .LT. 3) STFILE = 3
                  ELSE
                     STFILE = 4
                  ENDIF
                  IF (NRVAR .GE. IZRAND) THEN
                     IF (.NOT. EXTMEM_RAND(NRVAR+1,IERR)) GOTO 9200
                  ENDIF
!
                  NRVAR = NRVAR  + 1
                  SDIST(NRVAR) = DNAME(2)
                  MDIST(8*NRVAR-7) = 1
                  MDIST(8*NRVAR-6) = LIPAR
                  MDIST(8*NRVAR-5) = 0
                  MDIST(8*NRVAR-4) = LRPAR
                  MDIST(8*NRVAR-3) = 2
                  MDIST(8*NRVAR-2) = LLPAR
                  MDIST(8*NRVAR-1) = 0
                  MDIST(8*NRVAR  ) = ISTAGE
                  NINTPAR = 0
                  NREAPAR = 2
                  NLOCPAR = 0
                  IF (LRPAR+2 .GT. IZRPAR) THEN
                     IF (.NOT. EXTMEM_RPAR(LRPAR+2,IERR)) GOTO 9200
                  ENDIF
!
                  RANDPAR(LRPAR+1) = ATEMP(1)
                  RANDPAR(LRPAR+2) = ATEMP(2)
                  LRPAR = LRPAR + 2
                  NDCOLM = NDCOLM + 1
                  LDMTX(NDCOLM+1) = NDELEM + 1
                  INSERT = NDCOLM + 1
               ELSE
!
!     If the second name field is empty, the random variable is assumed
!     to be the name of a (component of a ) random variable described
!     in a DISTRIB section earlier. Such distributions can be reused
!     several times, so we must make a copy first.
!
                  DO I=1,NALIAS
                     IF (ALIAS(I) .EQ. DNAME(1)) GOTO 461
                  END DO
!
                  WRITE (IOLOG, 1900) NREC,QLINE
                  WRITE (IOLOG, 2500)
                  GOTO 890
!
  461          CONTINUE
                  IF (IALIAS(I) .GE. 0) THEN
!
!     Found a univariate distribution.
!
                     INC = 1
                     KRV = IALIAS(I)
                     AKA =  ALIAS(I)
                  ELSE
!
!     Found a multivariate distribution.
!
                     INC = -IALIAS(I)
                     KRV =  IALIAS(I + IALIAS(I))
                     AKA =   ALIAS(I + IALIAS(I))
                  ENDIF
!
!     Check if the random variable was accessed before in this block.
!
                  DO J=1,NCOPIED
                     IF (COPIED(J) .EQ. AKA) GOTO 465
                  END DO
!
!     Random variable was not copied before. Copy and allocate columns,
!     then mark insertion point before the start of the next column.
!
                  IF (NRVAR .GE. IZRAND) THEN
                     IF (.NOT. EXTMEM_RAND(NRVAR+1,IERR)) GOTO 9200
                  ENDIF
!
                  NRVAR = NRVAR + 1
                  SDIST(NRVAR) = SDIST(KRV)
                  MDIST(8*NRVAR-7) = MDIST(8*KRV-7)
                  MDIST(8*NRVAR-6) = MDIST(8*KRV-6)
                  MDIST(8*NRVAR-5) = MDIST(8*KRV-5)
                  MDIST(8*NRVAR-4) = MDIST(8*KRV-4)
                  MDIST(8*NRVAR-3) = MDIST(8*KRV-3)
                  MDIST(8*NRVAR-2) = MDIST(8*KRV-2)
                  MDIST(8*NRVAR-1) = MDIST(8*KRV-1)
                  MDIST(8*NRVAR  ) = ISTAGE
!
                  IF (NCOPIED .GE. IZCOPY) THEN
                     IF (.NOT. EXTMEM_COPY(NCOPIED+1,IERR)) GOTO 9200
                  ENDIF
!
                  NCOPIED = NCOPIED + 1
                  COPIED(NCOPIED) = AKA
                  MARKER(NCOPIED) = NDCOLM + 1
                  INSERT = NDCOLM + INC + 1
                  MNDCOL = NDCOLM + MDIST(8*KRV-7)
                  IF (MNDCOL .GT. IZDCOL) THEN
                     IF (.NOT. EXTMEM_DCOL(MNDCOL,IERR)) GOTO 9200
                  ENDIF
!
                  DO I=1,MDIST(8*KRV-7)
                     LDMTX(NDCOLM+I+1) = NDELEM + 1
                  END DO
                  NDCOLM = NDCOLM + MDIST(8*KRV-7)
                  GOTO 890
!
!     Random variable was copied before. Find the insertion point.
!
  465             CONTINUE
                     INSERT = MARKER(J) + INC
               ENDIF
               GOTO 890
!
!     'CR' card - MVNORMAL only
!     -------------------------
!
  500 CONTINUE
         IF (DTYPE .NE. DMNORM) GOTO 960
            DCODE = 'CR'
            GOTO 890
!
!     'CV' card - MVNORMAL only
!     -------------------------
!
  550 CONTINUE
         IF (DTYPE .NE. DMNORM) GOTO 960
            DCODE = 'CV'
            GOTO 890
!
!     'PR' card - SUB option only
!     ---------------------------
!
  570 CONTINUE
         DCODE = 'PR'
         LISTPAR = 1
         NREAPAR = 0
         NINTPAR = 0
         NLOCPAR = 0
!
!     Process the parameter list
!
  585 CONTINUE
         CALL GETPAR(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,XPAR,IPAR,MTYPE,      &
     &               IBLOCK,IERR,IOSTO,NREC)
         IF (Q1 .EQ. QAST) GOTO 585
         IF (IERR .GT. 0) THEN
            IF (IERR .EQ. 3) THEN
               LISTPAR = 0
               GOTO 100
            ENDIF
         ENDIF
         IF (MTYPE .EQ. 1) THEN
            IF (LIPAR+2 .GT. IZIPAR) THEN
               IF (.NOT. EXTMEM_IPAR(LIPAR+2,IERR)) GOTO 9200
            ENDIF
!
            NINTPAR = NINTPAR + 1
            LIPAR = LIPAR + 1
            IRNDPAR(LIPAR) = IPAR
            MDIST(8*NRVAR-5) = MDIST(8*NRVAR-5) + 1

         ELSEIF (MTYPE .EQ. 2) THEN
            IF (LRPAR+2 .GT. IZRPAR) THEN
               IF (.NOT. EXTMEM_RPAR(LRPAR+2,IERR)) GOTO 9200
            ENDIF
!
            NREAPAR = NREAPAR + 1
            LRPAR = LRPAR + 1
            RANDPAR(LRPAR) = XPAR
            MDIST(8*NRVAR-3) = MDIST(8*NRVAR-3) + 1

         ELSEIF (MTYPE .LT. 0) THEN
            IF (LLPAR+2 .GT. IZLPAR) THEN
               IF (.NOT. EXTMEM_LPAR(LLPAR+2,IERR)) GOTO 9200
            ENDIF
!
            NLOCPAR = NLOCPAR + 1
            LLPAR = LLPAR + 1
            LRNDPAR(3*LLPAR-2) = -MTYPE
            LRNDPAR(3*LLPAR-1) =  IBLOCK
            LRNDPAR(3*LLPAR  ) =  IPAR

            MDIST(8*NRVAR-1) = MDIST(8*NRVAR-1) + 1
         ENDIF
         GOTO 585
!
!     --------------------------------------------------------------------
!     Here we have a regular data card.
!     Process the card according to the segment it's in (BL/CR/CV/PR/RV).
!     --------------------------------------------------------------------
!
  350 CONTINUE
         IF (DCODE .EQ. 'BL') THEN
!
!     Setting up a new block (or perhaps repeat if DISCRETE)
!     ------------------------------------------------------
!
!     Identify the stochastic element and store the information, by type:
!        DISCRETE: value
!        MVNORMAL: mean, variance, alias; allot space for correlation matrix
!        LINTR:    component of c in transformation v = Hu + c
!        other:    nothing
!
         CALL IDENTIFY(NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2),ATEMP(1),      &
     &                    DNAME(3),ATEMP(2),ISTAGE,ISTOCH,IBLOCK,IRPOS,      &
     &                    IROW,ICOL)
         IF (DTYPE .EQ. DISCR .AND. NREAL .GT. 1) GOTO 370
!
!     ISTOCH gives the type of random element and is used for diagnostics:
!     If ISTOCH = 0, a row or column name was not found. Flush to next header.
!     If ISTOCH < 0, we make room for the coefficient within the A or Q matrix.
!     In addition we need to verify the time order: The stochastic element
!     must be in a stage no earlier than the current random vector.
!
            IF (ISTAGE .LT. MDIST(8*NRVAR)) THEN
               WRITE (IOLOG, 1400)
               GOTO 890
            ENDIF
            IF (ISTOCH .EQ.  0) GOTO 850
            IF (ISTOCH .EQ. 99) GOTO 890
               IF (ISTOCH .LT. 0) THEN
!
!     Enlarge the A-matrix
!
                  IF (ISTOCH .EQ. -1) THEN
                     IF (LASTA .GE. IZALMN) THEN
                        IF (.NOT. EXTMEM_ALMN(LASTA+1,IERR)) GOTO 9200
                     ENDIF
!
                     DO I=LASTA,IRPOS,-1
                         A(I+1) =  A(I)
                        IA(I+1) = IA(I)
                     END DO
                      A(IRPOS) = 0.D0
                     IA(IRPOS) = IROW
                     LASTA = LASTA + 1
!
                     ISTOCH = 1
                     DO I=1,LASTBA
                        IF (I .NE. IBLOCK .AND. KELMA(I) .GE. IRPOS-1)      &
     &                               KELMA(I) = KELMA(I) + 1
                     END DO
!
                     DO I=1,NPER
                        IF (IBLOCK .GT. KDATA(I) .AND.      &
     &                      IBLOCK .LE. KDATA(I+1) ) GOTO 353
                     END DO
!
!     This should not happen!
!
                     WRITE (IOLOG, 3500)
                     IERR = 2
                     GOTO 999
!
  353             CONTINUE
                     IF (MARKOV) THEN
                        NSUB = MIN(I,2)
                     ELSE
                        NSUB = I
                     ENDIF
                     IBL = IBLOCK - KDATA(I)
                     IF (IBL .LE. NSUB) THEN
                        IP = I - IBL + 1
                     ELSE
                        IP = I
                     ENDIF
                     LCOL = NODEINFO(IP)%NCOL - NODEINFO(IP)%NROW + 1
                     LL = KCOLA(IBLOCK)
                     DO I=ICOL+1,LCOL
                        LA(LL+I) = LA(LL+I) + 1
                     END DO
                     IRPOS = IRPOS - KELMA(IBLOCK)
                     NELMA(IBLOCK) = NELMA(IBLOCK) + 1
!
!     Update the stochastic elements found so far
!
                     DO I=1,NSTELM
                        IF (KSTOCH(4*I-3) .EQ.   1    .AND.      &
     &                      KSTOCH(4*I-2) .EQ. IBLOCK .AND.      &
     &                      KSTOCH(4*I-1) .GE. IRPOS )      &
     &                      KSTOCH(4*I-1) = KSTOCH(4*I-1) + 1
                     END DO
!
!     Enlarge the Q-matrix
!
                  ELSEIF (ISTOCH .EQ. -7) THEN
                     IF (LASTQ .GE. IZQLMN) THEN
                        IF (.NOT. EXTMEM_QLMN(LASTQ+1,IERR)) GOTO 9200
                     ENDIF
!
                     DO I=LASTQ,IRPOS,-1
                        AQMTX(I+1) = AQMTX(I)
                        IQMTX(I+1) = IQMTX(I)
                     END DO
                     AQMTX(IRPOS) = 0.D0
                     IQMTX(IRPOS) = IROW
                     LASTQ = LASTQ + 1
!
                     ISTOCH = 7
                     DO I=1,LASTQB
                        IF (I .NE. IBLOCK .AND. IQOFF(I) .GE. IRPOS-1)      &
     &                               IQOFF(I) = IQOFF(I) + 1
                     END DO
!
                     DO I=1,NPER
                        IBL = IBLOCK - KDATQ(I)
                        IF (IBL .GT. 0 .AND. IBL .LE. NQPROF ) GOTO 358
                     END DO
                     WRITE (IOLOG, 3700)
                     IERR = 2
                     GOTO 999
!
  358             CONTINUE
                     IP = I - IBL + 1
                     LCOL = NODEINFO(IP)%NCOL - NODEINFO(IP)%NROW + 1
                     LQB = LQOFF(IBLOCK)
                     DO I=ICOL+1,LCOL
                        LQMTX(LQB+I) = LQMTX(LQB+I) + 1
                     END DO
                     IRPOS = IRPOS - IQOFF(IBLOCK)
                     NELMQ(IBLOCK) = NELMQ(IBLOCK) + 1
!
!     Update the stochastic elements found so far
!
                     DO I=1,NSTELM
                        IF (KSTOCH(4*I-3) .EQ.   7    .AND.      &
     &                      KSTOCH(4*I-2) .EQ. IBLOCK .AND.      &
     &                      KSTOCH(4*I-1) .GE. IRPOS )      &
     &                      KSTOCH(4*I-1) = KSTOCH(4*I-1) + 1
                     END DO
!
!     Something happened to ISTOCH. Stop the program.
!
                  ELSE
                     WRITE (IOLOG, 3700)
                     IERR = 2
                     GOTO 999
                  ENDIF
               ENDIF
!
!     Store the information.
!
            IF (NSTELM .GE. IZSTOC) THEN
               IF (.NOT. EXTMEM_STOC(NSTELM+1,IERR)) GOTO 9200
            ENDIF
!
            NSTELM = NSTELM + 1
            KSTOCH(4*NSTELM-3) = ISTOCH
            KSTOCH(4*NSTELM-2) = IBLOCK
            KSTOCH(4*NSTELM-1) = IRPOS
            KSTOCH(4*NSTELM  ) = METHOD
            IF (DTYPE .EQ. DISCR) THEN
               IF (LRPAR .GE. IZRPAR) THEN
                  IF (.NOT. EXTMEM_RPAR(LRPAR+1,IERR)) GOTO 9200
               ENDIF
!
               IF (NDELEM .GE. IZDLMN) THEN
                  IF (.NOT. EXTMEM_DLMN(NDELEM+1,IERR)) GOTO 9200
               ENDIF
!
               IF (NDCOLM+1 .GE. IZDCOL) THEN
                  IF (.NOT. EXTMEM_DCOL(NDCOLM+2,IERR)) GOTO 9200
               ENDIF
!
               MDIST(8*NRVAR-7) = MDIST(8*NRVAR-7) + 1
               NREAPAR = NREAPAR + 1
               MDIST(8*NRVAR-3) = MDIST(8*NRVAR-3) + 1
!
               LRPAR = LRPAR + 1
               RANDPAR(LRPAR) = ATEMP(1)
!
               NDELEM = NDELEM + 1
               NDCOLM = NDCOLM + 1
                DMTX(NDELEM) = 1.D0
               IDMTX(NDELEM) = NSTELM
               LDMTX(NDCOLM+1) = NDELEM + 1
            ELSEIF (DTYPE .EQ. DLIN1) THEN
               IF (DABS(ATEMP(1)) .GT. ZTOLIN) THEN
                  IF (NDELEM .GE. IZDLMN) THEN
                     IF (.NOT. EXTMEM_DLMN(NDELEM+1,IERR)) GOTO 9200
                  ENDIF
                  IF (NDCOLM+1 .GE. IZDCOL) THEN
                     IF (.NOT. EXTMEM_DCOL(NDCOLM+2,IERR)) GOTO 9200
                  ENDIF
!
                  NDELEM = NDELEM + 1
                   DMTX(NDELEM) = ATEMP(1)
                  IDMTX(NDELEM) = NSTELM
                  LDMTX(NDCOLM+1) = NDELEM + 1
               ENDIF
            ELSEIF (DTYPE .EQ. DMNORM) THEN
               MDIST(8*NRVAR-7) = MDIST(8*NRVAR-7) + 1
               MNRPAR = LRPAR + MDIST(8*NRVAR-7) + 1
               IF (MNRPAR .GT. IZRPAR) THEN
                  IF (.NOT. EXTMEM_RPAR(MNRPAR,IERR)) GOTO 9200
               ENDIF
!
               RANDPAR(LRPAR+1) = ATEMP(1)
               RANDPAR(LRPAR+2) = ATEMP(2)
               DO I=2, MDIST(8*NRVAR-7)
                  RANDPAR(LRPAR+I+1) = 0.D0
               END DO
               MDIST(8*NRVAR-3) = MDIST(8*NRVAR-3)      &
     &                               + MDIST(8*NRVAR-7) + 1
               NREAPAR = NREAPAR + MDIST(8*NRVAR-7) + 1
               LRPAR = LRPAR + MDIST(8*NRVAR-7) + 1
               IF (NALIAS .GE. IZALIAS) THEN
                  IF (.NOT. EXTMEM_ALIAS(NALIAS+1,IERR)) GOTO 9200
               ENDIF
!
               IF (NDELEM .GE. IZDLMN) THEN
                  IF (.NOT. EXTMEM_DLMN(NDELEM+1,IERR)) GOTO 9200
               ENDIF
!
               IF (NDCOLM+1 .GE. IZDCOL) THEN
                  IF (.NOT. EXTMEM_DCOL(NDCOLM+2,IERR)) GOTO 9200
               ENDIF
!
               NALIAS = NALIAS + 1
                ALIAS(NALIAS) = DNAME(3)
               IALIAS(NALIAS) = - MDIST(8*NRVAR-7)
               NDELEM = NDELEM + 1
               NDCOLM = NDCOLM + 1
                DMTX(NDELEM) = 1.D0
               IDMTX(NDELEM) = NSTELM
               LDMTX(NDCOLM+1) = NDELEM + 1
            ELSE
               MDIST(8*NRVAR-7) = MDIST(8*NRVAR-7) + 1
               IF (NDELEM .GE. IZDLMN) THEN
                  IF (.NOT. EXTMEM_DLMN(NDELEM+1,IERR)) GOTO 9200
               ENDIF
!
               IF (NDCOLM+1 .GE. IZDCOL) THEN
                  IF (.NOT. EXTMEM_DCOL(NDCOLM+2,IERR)) GOTO 9200
               ENDIF
!
               NDELEM = NDELEM + 1
               NDCOLM = NDCOLM + 1
                DMTX(NDELEM) = 1.D0
               IDMTX(NDELEM) = NSTELM
               LDMTX(NDCOLM+1) = NDELEM + 1
            ENDIF
            GOTO 890
!
!     Repeat realization of a variable in a discrete block. Verify the
!     location by comparing ISTOCH, IBLOCK, IRPOS to the first realization.
!
  370    CONTINUE
            IDIM = MDIST(8*NRVAR-7)
            DO I=1,IDIM
               J = NSTELM - IDIM + I - 1
               IF (KSTOCH(4*J+1) .EQ. ISTOCH .AND.      &
     &             KSTOCH(4*J+2) .EQ. IBLOCK .AND.      &
     &             KSTOCH(4*J+3) .EQ. IRPOS         ) GOTO 377
            END DO
            WRITE (IOLOG, 1900) NREC, QLINE
            WRITE (IOLOG, 2400)
            IERR = 2
            GOTO 999
!
  377    CONTINUE
            RANDPAR(MDIST(8*NRVAR-4) + MDIST(8*NRVAR-3) - IDIM + I)      &
     &                                                        = ATEMP(1)
            GOTO 890
!
         ELSEIF (DCODE .EQ. 'RV') THEN
!
!     Accessing a new random variable (in LINTR segment)
!     --------------------------------------------------
!
            CALL IDENTIFY(NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2),ATEMP(1),      &
     &                    DNAME(3),ATEMP(2),ISTAGE,ISTOCH,IBLOCK,IRPOS,      &
     &                    IROW,ICOL)
            DO J=NPREV+1,NSTELM
               IF (ISTOCH .EQ. KSTOCH(4*J-3) .AND.      &
     &             IBLOCK .EQ. KSTOCH(4*J-2) .AND.      &
     &             IRPOS  .EQ. KSTOCH(4*J-1)      ) GOTO 379
            END DO
            WRITE (IOLOG, 1900) NREC, QLINE
            WRITE (IOLOG, 1300)
            IERR = 2
            GOTO 999
!
!     Write the information into the D matrix
!
  379       CONTINUE
               IF (DABS(ATEMP(1)) .GT. ZTOLIN) THEN
                  IF (NDELEM .GE. IZDLMN) THEN
                     IF (.NOT. EXTMEM_DLMN(NDELEM+1,IERR)) GOTO 9200
                  ENDIF
!
                  DO I=NDELEM,LDMTX(INSERT),-1
                      DMTX(I+1) =  DMTX(I)
                     IDMTX(I+1) = IDMTX(I)
                  END DO
                  IDMTX(LDMTX(INSERT)) = J
                   DMTX(LDMTX(INSERT)) = ATEMP(1)
                  DO I=INSERT,NDCOLM+1
                     LDMTX(I) = LDMTX(I) + 1
                  END DO
                  NDELEM = NDELEM + 1
               ENDIF
               GOTO 890
!
         ELSEIF (DCODE .EQ. 'CV' .OR. DCODE .EQ. 'CR') THEN
!
!     Correlation and covariance matrices are handled in similar fashion.
!     Start by verifying the two aliases.
!
            DO I1=1,NALIAS
               IF (ALIAS(I1) .EQ. DNAME(1)) GOTO 391
            END DO
!
            WRITE (IOLOG, 1900) NREC, QLINE
            WRITE (IOLOG, 2500)
            GOTO 890
!
  391    CONTINUE
            IF (IALIAS(I1) .GE. 0) THEN
               WRITE (IOLOG, 1900) NREC, QLINE
               WRITE (IOLOG, 2600)
               GOTO 890
            ENDIF
!
            NRV1 = IALIAS(I1 + IALIAS(I1))
            DO I2=1,NALIAS
               IF (ALIAS(I2) .EQ. DNAME(2)) GOTO 396
            END DO
!
            WRITE (IOLOG, 1900) NREC, QLINE
            WRITE (IOLOG, 2500)
            GOTO 890
!
  396    CONTINUE
            IF (IALIAS(I2) .GE. 0) THEN
               WRITE (IOLOG, 1900) NREC, QLINE
               WRITE (IOLOG, 2600)
               GOTO 890
            ENDIF
!
            NRV2 = IALIAS(I2 + IALIAS(I2))
            IF (NRV1 .NE. NRV2) THEN
               WRITE (IOLOG, 1900) NREC, QLINE
               WRITE (IOLOG, 2700)
               GOTO 890
            ENDIF
!
!     All is well: Random variables are two components of the same vector
!
            N1 = - IALIAS(I1)
            N2 = - IALIAS(I2)
            IF (N1 .LT. N2) THEN
               LCOV = MDIST(8*NRVAR-4) + N2*(N2+1)/2 + N1 + 1
               IF (DCODE .EQ. 'CV') THEN
                  RANDPAR(LCOV) = ATEMP(1)
               ELSE
                  LVAR1 = MDIST(8*NRVAR-4) + N1*(N1+1)/2 + 1
                  LVAR2 = MDIST(8*NRVAR-4) + N2*(N2+1)/2 + 1
                  RANDPAR(LCOV)      &
     &                 = ATEMP(1) * DSQRT(RANDPAR(LVAR1)*RANDPAR(LVAR2))
                  GOTO 890
               ENDIF
            ELSEIF (N1 .GT. N2) THEN
               LCOV = MDIST(8*NRVAR-4) + N1*(N1+1)/2 + N2 + 1
               IF (DCODE .EQ. 'CV') THEN
                  RANDPAR(LCOV) = ATEMP(1)
               ELSE
                  LVAR1 = MDIST(8*NRVAR-4) + N1*(N1+1)/2 + 1
                  LVAR2 = MDIST(8*NRVAR-4) + N2*(N2+1)/2 + 1
                  RANDPAR(LCOV)      &
     &                 = ATEMP(1) * DSQRT(RANDPAR(LVAR1)*RANDPAR(LVAR2))
                  GOTO 890
               ENDIF
            ELSE
               WRITE (IOLOG, 1900) NREC, QLINE
               WRITE (IOLOG, 2800) I1,I2,N1,N2
               GOTO 890
            ENDIF
!
         ENDIF
!
!     ======================================================================
!     Start of a new DISTRIB segment.
!
!     Univariate and multivariate distributions are distinguished by the
!     absence or presence of a 'BL' type card. Different multivariate types
!     allow different codes to appear in the code field, but the first code
!     must be of 'BL' type. We track the last code read in variable DCODE.
!
!     For the record, here are the legal values for the code field:
!
!        MVNORMAL:           BL, CV, CR
!        Other (subroutine): BL, PR
!     ======================================================================
!
  130 CONTINUE
         LISTPAR = 0
         DCODE = 'XX'
         ITYPE = 5
         DTYPE = DNAME(2)
         IF (DTYPE .EQ. DISCR) THEN
            IDISCR = 1
            IF (STFILE .LT. 3) STFILE = 3
         ELSE
            IDISCR = 0
            STFILE = 4
         ENDIF
         DMODE = DNAME(3)
         METHOD = 1
         IF (DMODE .NE. DBLANK .AND. DMODE .NE. DREP) THEN
            WRITE (IOLOG, 1900) NREC, QLINE
            WRITE (IOLOG, 2100)

         ENDIF
         GOTO 890
!
!     -------------------------------------------------------------------------
!     Data card for DISTRIB section. If the next card has code 'BL', the random
!     variable is multivariate (NOT DISCRETE!); else it must be univariate.
!     -------------------------------------------------------------------------
!
  599 CONTINUE
         IF (Q2 .EQ. QB .AND. Q3 .EQ. QL) THEN
            IF (IDISCR .EQ. 0) THEN
                ITYPE = 4
                GOTO 700
            ELSE
               WRITE (IOLOG, 1900) NREC, QLINE
               WRITE (IOLOG, 3000)
               GOTO 850
            ENDIF
         ELSE
            ITYPE = 3
         ENDIF
!
!     Here the distribution is univariate.
!     Check if repeat realization of a discrete distribution.
!
  600 CONTINUE
         IF (DNAME(1) .EQ. DN1 .AND. DTYPE .EQ. DISCR) GOTO 617
!
!     The random variable is not associated with any stochastic element (yet),
!     so there is much less information to process and store.
!
            IF (NRVAR .GE. IZRAND) THEN
               IF (.NOT. EXTMEM_RAND(NRVAR+1,IERR)) GOTO 9200
            ENDIF
!
            NRVAR = NRVAR  + 1
            SDIST(NRVAR) = DTYPE
            MDIST(8*NRVAR-7) = 1
            MDIST(8*NRVAR-6) = LIPAR
            MDIST(8*NRVAR-5) = 0
            MDIST(8*NRVAR-4) = LRPAR
            MDIST(8*NRVAR-3) = 2
            MDIST(8*NRVAR-2) = LLPAR
            MDIST(8*NRVAR-1) = 0
            MDIST(8*NRVAR  ) = 0
            NINTPAR = 0
            NREAPAR = 2
            NLOCPAR = 0
            IF (LRPAR+2 .GT. IZRPAR) THEN
               IF (.NOT. EXTMEM_RPAR(LRPAR+2,IERR)) GOTO 9200
            ENDIF
!
            IF (DTYPE .EQ. DISCR) THEN
               RANDPAR(LRPAR+1) = ATEMP(2)
               RANDPAR(LRPAR+2) = ATEMP(1)
            ELSE
               RANDPAR(LRPAR+1) = ATEMP(1)
               RANDPAR(LRPAR+2) = ATEMP(2)
            ENDIF
            LRPAR = LRPAR + 2
            IF (NALIAS .GE. IZALIAS) THEN
               IF (.NOT. EXTMEM_ALIAS(NALIAS+1,IERR)) GOTO 9200
            ENDIF
!
            NALIAS = NALIAS + 1
             ALIAS(NALIAS) = DNAME(1)
            IALIAS(NALIAS) = NRVAR
            IF (NDCOLM+1 .GE. IZDCOL) THEN
               IF (.NOT. EXTMEM_DCOL(NDCOLM+2,IERR)) GOTO 9200
            ENDIF
!
            NDCOLM = NDCOLM + 1
            LDMTX(NDCOLM+1) = NDELEM + 1
            GOTO 890
!
!     Repeat realization of a discrete random variable.
!
  617 CONTINUE
         IF (LRPAR+2 .GT. IZRPAR) THEN
            IF (.NOT. EXTMEM_RPAR(LRPAR+2,IERR)) GOTO 9200
         ENDIF
!
         RANDPAR(LRPAR + 1) = ATEMP(2)
         RANDPAR(LRPAR + 2) = ATEMP(1)
         NREAPAR = NREAPAR + 2
         LRPAR = LRPAR + 2
         MDIST(8*NRVAR-3) = NREAPAR
         GOTO 890
!
!     ---------------------------------------------------------------
!     Data card for multivariate distribution. Check the code field:
!     The first non-comment card must be a 'BL' card.
!     ---------------------------------------------------------------
!
  650 CONTINUE
         IF (Q2 .EQ. QB .AND. Q3 .EQ. QL) GOTO 700
         IF (Q2 .EQ. QC .AND. Q3 .EQ. QR) GOTO 750
         IF (Q2 .EQ. QC .AND. Q3 .EQ. QV) GOTO 760
         IF (Q2 .EQ. QP .AND. Q3 .EQ. QR) GOTO 770
            GOTO 660
!
!     A 'BL' card starts the definition of locations. Keep the info.
!     --------------------------------------------------------------
!
  700    CONTINUE
            DCODE = 'BL'
            NPREV  = NSTELM
            IF (BLNAME .EQ. DNAME(1)) GOTO 950
               IF (NRVAR .GE. IZRAND) THEN
                  IF (.NOT. EXTMEM_RAND(NRVAR+1,IERR)) GOTO 9200
               ENDIF
!
               BLNAME = DNAME(1)
               NRVAR  = NRVAR  + 1
               SDIST(NRVAR) = DTYPE
               MDIST(8*NRVAR-7) = 0
               MDIST(8*NRVAR-6) = LIPAR
               MDIST(8*NRVAR-5) = 0
               MDIST(8*NRVAR-4) = LRPAR
               MDIST(8*NRVAR-3) = 0
               MDIST(8*NRVAR-2) = LLPAR
               MDIST(8*NRVAR-1) = 0
               MDIST(8*NRVAR  ) = 0
               IF (DTYPE .EQ. DMNORM) THEN
                  IF (LIPAR .GE. IZIPAR) THEN
                     IF (.NOT. EXTMEM_IPAR(LIPAR+1,IERR)) GOTO 9200
                  ENDIF
!
                  NINTPAR = 1
                  MDIST(8*NRVAR-5) = 1
                  LIPAR = LIPAR + 1
                  IRNDPAR(LIPAR) = 0
               ELSE
                  NINTPAR = 0
               ENDIF
               NREAPAR = 0
               NLOCPAR = 0
               IF (NALIAS .GE. IZALIAS) THEN
                  IF (.NOT. EXTMEM_ALIAS(NALIAS+1,IERR)) GOTO 9200
               ENDIF
!
               NALIAS = NALIAS + 1
                ALIAS(NALIAS) = BLNAME
               IALIAS(NALIAS) = NRVAR
               GOTO 890
!
!     'CR' card - MVNORMAL only
!     -------------------------
!
  750 CONTINUE
         IF (DTYPE .NE. DMNORM) GOTO 960
            DCODE = 'CR'
            GOTO 890
!
!     'CV' card - MVNORMAL only
!     -------------------------
!
  760 CONTINUE
         IF (DTYPE .NE. DMNORM) GOTO 960
            DCODE = 'CV'
            GOTO 890
!
!     'PR' card - SUB option only
!     ---------------------------
!
  770 CONTINUE
         DCODE = 'PR'
         LISTPAR = 1
         NREAPAR = 0
         NINTPAR = 0
         NLOCPAR = 0
!
!     Process the parameter list
!
  785 CONTINUE
         CALL GETPAR(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,XPAR,IPAR,MTYPE,      &
     &               IBLOCK,IERR,IOSTO,NREC)
         IF (Q1 .EQ. QAST) GOTO 785
         IF (IERR .GT. 0) THEN
            IF (IERR .EQ. 3) THEN
               LISTPAR = 0
               GOTO 100
            ENDIF
         ENDIF
         IF (MTYPE .EQ. 1) THEN
            IF (LIPAR+2 .GT. IZIPAR) THEN
               IF (.NOT. EXTMEM_IPAR(LIPAR+2,IERR)) GOTO 9200
            ENDIF
!
            NINTPAR = NINTPAR + 1
            LIPAR = LIPAR + 1
            IRNDPAR(LIPAR) = IPAR
            MDIST(8*NRVAR-5) = MDIST(8*NRVAR-5) + 1

         ELSEIF (MTYPE .EQ. 2) THEN
            IF (LRPAR+2 .GT. IZRPAR) THEN
               IF (.NOT. EXTMEM_RPAR(LRPAR+2,IERR)) GOTO 9200
            ENDIF
!
            NREAPAR = NREAPAR + 1
            LRPAR = LRPAR + 1
            RANDPAR(LRPAR) = XPAR
            MDIST(8*NRVAR-3) = MDIST(8*NRVAR-3) + 1

         ELSEIF (MTYPE .LT. 0) THEN
            IF (LLPAR+2 .GT. IZLPAR) THEN
               IF (.NOT. EXTMEM_LPAR(LLPAR+2,IERR)) GOTO 9200
            ENDIF
!
            NLOCPAR = NLOCPAR + 1
            LLPAR = LLPAR + 1
            LRNDPAR(3*LLPAR-2) = -MTYPE
            LRNDPAR(3*LLPAR-1) =  IBLOCK
            LRNDPAR(3*LLPAR  ) =  IPAR
            MDIST(8*NRVAR-1) = MDIST(8*NRVAR-1) + 1

         ENDIF
         GOTO 785
!
!     ---------------------------------------------------
!     Here we have a regular data card.
!     Process according to the section it is in.
!     ---------------------------------------------------
!
  660 CONTINUE
         IF (DCODE .EQ. 'BL') THEN
!
!        Store the information, depending on type:
!        MVNORMAL: mean, variance, alias; also allot space for correlation matrix
!        other:    alias only
!
            MDIST(8*NRVAR-7) = MDIST(8*NRVAR-7) + 1
            IF (DTYPE .EQ. DMNORM) THEN
               MNRPAR = LRPAR + MDIST(8*NRVAR-7) + 1
               IF (MNRPAR .GT. IZRPAR) THEN
                  IF (.NOT. EXTMEM_RPAR(IZRPAR,IERR)) GOTO 9200
               ENDIF
!
               RANDPAR(LRPAR+1) = ATEMP(1)
               RANDPAR(LRPAR+2) = ATEMP(2)
               DO I=2, MDIST(8*NRVAR-7)
                  RANDPAR(LRPAR+I+1) = 0.D0
               END DO
               MDIST(8*NRVAR-3) = MDIST(8*NRVAR-3)      &
     &                               + MDIST(8*NRVAR-7) + 1
               NREAPAR = NREAPAR + MDIST(8*NRVAR-7) + 1
               LRPAR = LRPAR + MDIST(8*NRVAR-7) + 1
            ENDIF
!
            IF (NALIAS .GE. IZALIAS) THEN
               IF (.NOT. EXTMEM_ALIAS(NALIAS+1,IERR)) GOTO 9200
            ENDIF
!
            IF (NDCOLM+1 .GE. IZDCOL) THEN
               IF (.NOT. EXTMEM_DCOL(NDCOLM+2,IERR)) GOTO 9200
            ENDIF
!
            NALIAS = NALIAS + 1
             ALIAS(NALIAS) = DNAME(1)
            IALIAS(NALIAS) = - MDIST(8*NRVAR-7)
            NDCOLM = NDCOLM + 1
            LDMTX(NDCOLM+1) = NDELEM + 1
            GOTO 890
         ELSEIF (DCODE .EQ. 'CV' .OR. DCODE .EQ. 'CR') THEN
!
!     Correlation and covariance matrices are handled in similar fashion.
!     Start by verifying the two aliases.
!
            DO I1=1,NALIAS
               IF (ALIAS(I1) .EQ. DNAME(1)) GOTO 681
            END DO
!
            WRITE (IOLOG, 1900) NREC, QLINE
            WRITE (IOLOG, 2500)
            GOTO 890
!
  681    CONTINUE
            IF (IALIAS(I1) .GE. 0) THEN
               WRITE (IOLOG, 1900) NREC, QLINE
               WRITE (IOLOG, 2600)
               GOTO 890
            ENDIF
!
            NRV1 = IALIAS(I1 + IALIAS(I1))
            DO I2=1,NALIAS
               IF (ALIAS(I2) .EQ. DNAME(2)) GOTO 691
            END DO
!
            WRITE (IOLOG, 1900) NREC, QLINE
            WRITE (IOLOG, 2500)
            GOTO 890
!
  691    CONTINUE
            IF (IALIAS(I2) .GE. 0) THEN
               WRITE (IOLOG, 1900) NREC, QLINE
               WRITE (IOLOG, 2600)
               GOTO 890
            ENDIF
!
            NRV2 = IALIAS(I2 + IALIAS(I2))
            IF (NRV1 .NE. NRV2) THEN
               WRITE (IOLOG, 1900) NREC, QLINE
               WRITE (IOLOG, 2700)
               GOTO 890
            ENDIF
!
!     All is well: Random variables are two components of the same vector
!
            N1 = - IALIAS(I1)
            N2 = - IALIAS(I2)
            IF (N1 .LT. N2) THEN
               LCOV = MDIST(8*NRVAR-4) + N2*(N2+1)/2 + N1 + 1
               IF (DCODE .EQ. 'CV') THEN
                  RANDPAR(LCOV) = ATEMP(1)
               ELSE
                  LVAR1 = MDIST(8*NRVAR-4) + N1*(N1+1)/2 + 1
                  LVAR2 = MDIST(8*NRVAR-4) + N2*(N2+1)/2 + 1
                  RANDPAR(LCOV)      &
     &                 = ATEMP(1) * DSQRT(RANDPAR(LVAR1)*RANDPAR(LVAR2))
               ENDIF
            ELSEIF (N1 .GT. N2) THEN
               LCOV = MDIST(8*NRVAR-4) + N1*(N1+1)/2 + N2 + 1
               IF (DCODE .EQ. 'CV') THEN
                  RANDPAR(LCOV) = ATEMP(1)
               ELSE
                  LVAR1 = MDIST(8*NRVAR-4) + N1*(N1+1)/2 + 1
                  LVAR2 = MDIST(8*NRVAR-4) + N2*(N2+1)/2 + 1
                  RANDPAR(LCOV)      &
     &                 = ATEMP(1) * DSQRT(RANDPAR(LVAR1)*RANDPAR(LVAR2))
               ENDIF
            ELSE
               WRITE (IOLOG, 1900) NREC, QLINE
               WRITE (IOLOG, 2800) I1,I2,N1,N2
            ENDIF
            GOTO 890
         ELSE
            WRITE (IOLOG, 1900) NREC, QLINE
            WRITE (IOLOG, 2900)
            GOTO 850
         ENDIF
!
!     Flush to next header card
!
  850 CONTINUE
  860 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,      &
     &                     DNAME,ATEMP,IERR,IOSTO,NREC,LENGTH)
         IF (Q1 .EQ. QBL .OR. Q1 .EQ. QAST) GOTO 860
         IF (IERR .GT. 0) GOTO 900
            GOTO 100
!
!     Read the next record
!
  890 CONTINUE
         DN1 = DNAME(1)
         DN2 = DNAME(2)
  891 CONTINUE
         CALL GETMPS(QLINE,Q1,Q2,Q3,Q4,      &
     &                     DNAME,ATEMP,IERR,IOSTO,NREC,LENGTH)
         IF (Q1 .EQ. QAST) GOTO 891
         IF (NECHO1 .GE. 2) WRITE (IOLOG, 1900) NREC, QLINE
         IF (IERR .GT. 0) GOTO 980
            GOTO 100
!
!     Found an ENDATA card
!
  900 CONTINUE
         IF (NECHO1 .GE. 1) WRITE (IOLOG, 1900) NREC, QLINE
         GOTO 999
!
!     Insufficient memory
!
 9200 CONTINUE
      WRITE (IOLOG, 3200)
      IERR = 2
      GOTO 999
!
!     Other error
!
  950 CONTINUE
  960 CONTINUE
      WRITE (IOLOG, 1900) NREC, QLINE
      WRITE (IOLOG, 3960)
      IERR = 2
      GOTO 999
!
!     BLOCK segment does not contain 'BL' card
!
  970 CONTINUE
         WRITE (IOLOG, 1900) NREC, QLINE
         WRITE (IOLOG, 2300)
         IERR = 2
         GOTO 999
!
!     End of file during read
!
  980 CONTINUE
         WRITE (IOLOG, 1900) NREC, QLINE
         WRITE (IOLOG, 1700)
  999 CONTINUE
         DEALLOCATE(MARKER,COPIED,KPATH,KREF,STAT=IER)
!CC         DEALLOCATE(ALIAS,IALIAS,MARKER,COPIED,KPATH,KREF,STAT=IER)
         IF (IER .NE. 0) IERR = 2
         RETURN
!
 1000 FORMAT(' XXX -  FATAL  - Unidentifiable header card. ')
 1200 FORMAT(' XXX -  FATAL  - Variable name',A8,' not matched.')
 1250 FORMAT(' XXX - WARNING - Stage not given; assuming stage',I2)
 1300 FORMAT(' XXX -  FATAL  - Misspecified matrix element')
 1400 FORMAT(' XXX - WARNING - Inconsistent stage information.',      &
     &       ' Record ignored.')
 1700 FORMAT(' XXX - WARNING - Missing ENDATA card.')
 1900 FORMAT(I8,4X,A80)
 2100 FORMAT(' XXX - WARNING - Instruction not recognized.',      &
     &       ' Assume ''Replace.''')
 2300 FORMAT(' XXX -  FATAL  - BLOCK segment without ''BL'' card.')
 2400 FORMAT(' XXX -  FATAL  - Location not defined in first',      &
     &       ' appearance of this block.')
 2500 FORMAT(' XXX - WARNING - Referencing non-existent alias.',      &
     &       ' Data ignored.')
 2600 FORMAT(' XXX - WARNING - Referencing random vector instead',      &
     &       ' of its component. Data ignored.')
 2700 FORMAT(' XXX - WARNING - Two aliases reference different random',      &
     &       ' variables. Data ignored.')
 2800 FORMAT(' XXX - WARNING - Two aliases reference the same',      &
     &       ' component. Data ignored.',/,4I8)
 2900 FORMAT(' XXX - WARNING - Unrecognized code in code field.')
 3000 FORMAT(' XXX - WARNING - Unsupported distribution.',      &
     &       ' Flush to the next header.')
 3200 FORMAT(' XXX -  FATAL  - Memory allocation failed in INCONT.')
 3500 FORMAT(' XXX -  FATAL  - Cannot extend A-matrix.')
 3600 FORMAT(' XXX -  FATAL  - Cannot extend Q-matrix.')
 3700 FORMAT(' XXX -  FATAL  - Illegal value for ISTOCH.')
 3960 FORMAT(' XXX -  FATAL  - Header value or block name not valid in',      &
     &       ' this context.')
!
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE IDENTIFY(NREC,   Q1,Q2,Q3,Q4, DNAME1, DNAME2, ATEMP1,      &
     &                    DNAME3, ATEMP2, ISTAGE, ISTOCH,      &
     &                    IBLOCK, IRPOS, IROW,   ICOL)
!
!     This routine takes as input a recently read record from the stoch file
!     and checks whether the names in fields 1 and 2 describe a valid
!     location. If yes, it returns the stage, type and location in variables
!     ISTAGE, ISTOCH, IBLOCK, IRPOS, IROW and ICOL according to these rules:
!
!        ISTAGE is the stage in which the stochastic element becomes known.
!        This information may be given explicitly in the stoch file, or it
!        may be derived from the location of the stochastic element. The
!        convention used is that a stage consists of (event - decision),
!        that is, the stochastic elements associated with a stage are known
!        before the decisions are taken. This means that the timing of
!        probabilistic constraints must be adjusted.
!
!        ISTOCH =  1 for stochastic entry in one of the A-matrix blocks
!        ISTOCH = -1 for stochastic A-matrix entry not listed in core file
!        ISTOCH =  2 for stochastic cost
!        ISTOCH =  3 for stochastic RHS
!        ISTOCH =  4 for stochastic range (later converted to upper or lower bound)
!        ISTOCH =  5 for stochastic lower bound
!        ISTOCH =  6 for stochastic upper bound
!        ISTOCH =  7 for stochastic entry in one of the Q-matrix blocks
!        ISTOCH = -7 for stochastic Q-matrix entry not listed in core file
!        ISTOCH =  8 for stochastic RHS in an ordinary   chance-constraint
!        ISTOCH =  9 for stochastic RHS in an integrated chance-constraint
!        ISTOCH = 10 for stochastic curvature in piecewise linear-quadratic penalty
!        ISTOCH = 11 for fixed stochastic bound (affects upper and lower)
!        ISTOCH = 12 for stochastic demand (network option)
!        ISTOCH = 13 for stochastic supply (network option)
!        ISTOCH = 18 for reference to a location in the X array (primal value)
!        ISTOCH = 19 for reference to a location in the YPI array (dual value)
!        ISTOCH =  0 for unidentifiable stochastic elements
!        ISTOCH = 99 for stochastic elements that should be ignored
!
!     ISTOCH = 18 and ISTOCH = 19 are used in the parameter section of a
!     user-defined subroutine only.
!
!     Most quantities that can be specified in a network description have
!     a straightforward analog in the LP description (i.e., costs, bounds
!     and multipliers), but supply and demand are special forms of right-
!     hand sides (with an implied sign) that need to be tracked separately.
!
!     The nature of the stochastic element is determined from information found
!     in the first two name fields, and possibly in the code field. The first
!     name field can contain a tag (i.e., the identifier for a right-hand side,
!     ranges or bounds section), or a valid row name, or a valid column name.
!     The second name field can contain a valid row name or a valid column name.
!     If the card defines random supply or demand (as indicated in the code
!     field, the second name field can be empty. Finally, the second name
!     field may contain the tag 'PLQUAD' (including quotes) to define the
!     stochastic curvature in a piecewise linear-quadratic penalty. This
!     gives the following possibilities:
!
!       --------------------------------------------------------------------
!      | First ||                     Second name field                     |
!      | name  ||-----------------------------------------------------------
!      | field ||       row name    |    column name    |      (empty)      |
!       ====================================================================
!      |       ||         rhs       |      bound        |  (can't happen)   |
!      |  tag  ||        range      | (need code field) |                   |
!      |       || chance constraint |                   |                   |
!      |       ||         ICC       |                   |                   |
!       --------------------------------------------------------------------
!      |  row  ||   network card    |  (can't happen)   |   network card    |
!      |       || (need code field) |                   | (need code field) |
!       --------------------------------------------------------------------
!      |       ||         cost      |     Q-coeff       |  (can't happen)   |
!      |  col  ||       A-coeff     |                   |                   |
!      |       ||                   |                   |     'PLQUAD'      |
!      |       ||                   |                   |     curvature     |
!       --------------------------------------------------------------------
!
!     For stochastic elements other than matrix entries, IBLOCK gives
!     the number of the period containing the row or column in question.
!
!     For stochastic A-matrix elements, IBLOCK is as in the offset array
!     KDATA, that is, depending on the value of MARKOV and HAVE_GC:
!
!        period         staircase problems      triangular problems
!          1                 1                       1
!          2                 3   2                   3   2
!          3                     5   4               6   5   4
!          4                         7  6           10   9   8   7
!         ...                          ...                ...
!
!                       staircase with          triangular problems
!        period         linking constraints     with linking constraints
!          1                 1   4   7               1   4   8  14
!          2                 3   2   8               3   2   9  15
!          3                     6   5               7   6   5  16
!          4                        10   9          13  12  11  10
!         ...                          ...                ...
!
!     For stochastic Q-matrix elements, IBLOCK is as in the offset array
!     KDATQ, that is, depending on the value of NQPROF:
!
!        period    NQPROF = 1       NQPROF = 2        NQPROF = 3      ...
!          1       1                1                 1
!          2           2            3   2             3   2
!          3               3            5   4         6   5   4
!          4                   4            7   6         9   8   7
!         ...                ...             ...            ...
!
!
!     The location of the data item is recorded as the relative address
!     within the block of the relevant array in IRPOS.
!
!     Unlike other data items, matrix elements are stored in sparse form.
!     The user could therefore specify a legal location (row and column names
!     check out), that nevertheless does not have a storage location
!     associated with it (most likely because the user forgot to record
!     a zero value in the core file). This should generate an error, but it
!     ought to be recoverable --- the matrix element occurs in row IROW and
!     can be added at the end of the current column ICOL (i.e., in position
!     LA(ICOL+1)). We record the insertion point in this routine and flag the
!     situation by returning a NEGATIVE value in ISTOCH, ISTOCH = -1 for
!     A-matrix elements and ISTOCH = -7 for Q-matrix elements. IRPOS is
!     meaningless in this case.
!
!     If the name fields indicate a nonexisting row or column, the program
!     generates a warning error and returns the value ISTOCH = 0.
!
!     -------------------------------------------------------------------
!     This version dated 25 November 2002.
!     -------------------------------------------------------------------
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME1,DNAME2,DNAME3,DTPER,DCOL,DROW,DBLANK,DT
      CHARACTER*8   DPRIMAL,DDUAL,DPLQUAD
      CHARACTER*1   QBLANK(NAMLEN),QTEMP(NAMLEN)
      TYPE (CHANCE_ROW) TR
!
      EQUIVALENCE (QBLANK,DBLANK), (QTEMP,DT)
!
      DATA QBLANK /NAMLEN*' '/, DPRIMAL/'''PRIMAL'''/,DDUAL/'''DUAL'''/,      &
     &                          DPLQUAD/'''PLQUAD'''/
!
      ISTAGE = 0
      ISTOCH = 0
      IBLOCK = 0
      IRPOS  = 0
      IROW   = 0
!
!     First check if the stage is given explicitly.
!
         DTPER = DNAME3
         IF (DTPER .EQ. DBLANK) GOTO 110
            DO ISTAGE=1,NPER
               IF (DTPER .EQ. DTIME(ISTAGE)) GOTO 110
            END DO
            ISTAGE = 0
!
!     Next find location of the stochastic element (rhs, cost vector, etc.)
!
  110    CONTINUE
         DCOL = DNAME1
         DROW = DNAME2
         ISTOCH = 1
         IF (DCOL .EQ. DRHS .OR. DCOL .EQ. DBOUND .OR. DCOL .EQ. DRANGE      &
     &            .OR. DCOL .EQ. DICC .OR. DCOL   .EQ. DCHANCE      &
     &            .OR. DCOL .EQ. DPRIMAL .OR. DCOL .EQ. DDUAL) THEN
!
!     First name field is a tag. (Right-hand side, bound or variable reference.)
!     ==========================================================================
!
            IF (DCOL .EQ. DRHS   ) ISTOCH = 3
            IF (DCOL .EQ. DRANGE ) ISTOCH = 4
            IF (DCOL .EQ. DPRIMAL) ISTOCH = 18
            IF (DCOL .EQ. DDUAL  ) ISTOCH = 19
            IF (DCOL .EQ. DCHANCE) ISTOCH = 8
            IF (DCOL .EQ. DICC   ) ISTOCH = 9
!
!     For bounds, use the code field to distinguish the bound type
!
            IF (DCOL .EQ. DBOUND ) THEN
               IF     (Q2 .EQ. QL .AND. Q3 .EQ. QO) THEN
                  ISTOCH = 5
               ELSEIF (Q2 .EQ. QU .AND. Q3 .EQ. QP) THEN
                  ISTOCH = 6
               ELSEIF (Q2 .EQ. QU .AND. Q3 .EQ. QI) THEN
                  ISTOCH = 6
                  ATEMP1 = DINT(ATEMP1)
               ELSEIF (Q2 .EQ. QF .AND. Q3 .EQ. QX) THEN
                  ISTOCH = 11
               ELSE
!
!     Inappropriate code in code field.
!
                  WRITE (IOLOG, 1100)      &
     &               NREC,Q1,Q2,Q3,Q4,DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                  WRITE (IOLOG, 1200) Q2,Q3
                  ISTOCH = 0
                  RETURN
               ENDIF
            ENDIF
!
!     Identify the second name field.
!
            DO IRP=1,NPER
               JND = KPATH(IRP)
               KN  = NODEINFO(JND)%KNAMES
               NC  = NODEINFO(JND)%NCOL
               IRPOS = NMSRCH(DROW,KN+1,KN+NC)
               IF (IRPOS .GT. 0) GOTO 116
            ENDDO
!
!     Not matched. Could still be the name of a joint chance constraint...
!
            IF (ISTOCH .EQ. 8) THEN
               DO I=1,LCHANCE
                  IF (XCHANCE(I)%NAMPTR .GT. 0) THEN
                     IF (QCHANCE(XCHANCE(I)%NAMPTR) .EQ. DROW) THEN
                        IBLOCK = XCHANCE(I)%NODE
                        IRPOS  = I - KCHANCE(IBLOCK)
                        IF (ISTAGE .EQ. 0) THEN
                           ISTAGE = XCHANCE(I)%NODE
                        ELSEIF (ISTAGE .LT. XCHANCE(I)%NODE) THEN
                           WRITE (IOLOG, 2000)
                           ISTOCH = 0
                        ENDIF
                        RETURN
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
!
            WRITE (IOLOG, 1100)      &
     &         NREC,Q1,Q2,Q3,Q4,DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
            WRITE (IOLOG, 1000) DROW
            ISTOCH = 0
            RETURN
!
!     Matched. Perform a sanity check: bounds and primal variables
!     need a column name, all other tags need a row name.
!
  116       CONTINUE
               IRPOS = IRPOS - KN
               IF (ISTOCH .EQ. 5 .OR. ISTOCH .EQ. 6 .OR. ISTOCH .EQ. 11      &
     &                           .OR. ISTOCH .EQ. 18) THEN
                  IF (IRPOS .LE. NODEINFO(JND)%NROW) THEN
                     WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                        DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                     WRITE (IOLOG, 1400) DROW
                     ISTOCH = 0
                     RETURN
                  ENDIF
               ELSE
                  IF (IRPOS .GT. NODEINFO(JND)%NROW) THEN
                     WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                        DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                     WRITE (IOLOG, 1300) DROW
                     ISTOCH = 0
                     RETURN
                  ENDIF
               ENDIF
               IBLOCK  = IRP
!
!     For stochastic ranges and RHS pertaining to rows in chance constraints
!     we must adjust the stage information if it is not given explicitly.
!
               IF (ISTAGE .EQ. 0) THEN
                  ISTAGE = IRP
                  IF (ISTOCH .EQ. 3 .OR. ISTOCH .EQ. 4) THEN
                     IF (LCHANCE .GT. 0) THEN
                        DO JRP=IRP,NPER
                           DO IC=1,NCHANCE(JRP)
                              JCR = XCHANCE(IC)%ROWPTR
                              DO ICR=1,XCHANCE(IC)%SIZE
                                 TR = RCHANCE(ICR+JCR)
                                 IF (TR%RELROW .EQ. IRPOS .AND.      &
     &                               TR%STAGE  .EQ. IRP  ) THEN
                                    ISTAGE = JRP + 1
                                    GOTO 120
                                 END IF
                              END DO
                           END DO
                        END DO
                     ENDIF
                     IF (LICC .GT. 0) THEN
                        DO IC=1,NICC(IRP)
                           IF (XICC(IC)%RELROW .EQ. IRPOS) THEN
                              ISTAGE = IRP + 1
                              GOTO 120
                           END IF
                        END DO
                     ENDIF
                  ENDIF
!
!     Check nonanticipativity
!
               ELSE
                  IF (ISTAGE .LT. IRP) THEN
                     WRITE (IOLOG, 2000)
                     ISTOCH = 0
                     RETURN
                  ENDIF
               ENDIF
!
!     If we are dealing with a stochastic range, determine whether
!     the upper or lower bound is affected and reset ISTOCH.
!
  120       CONTINUE
               IF (ISTOCH .EQ. 4) THEN
                  IT = NODEINFO(IRP)%VARPTR + IRPOS
                  IF (VARTYP(IT) .EQ. 1) THEN
                     ISTOCH = 5
                  ELSEIF (VARTYP(IT) .EQ. -1) THEN
                     ISTOCH = 6
                  ELSEIF (VARTYP(IT) .EQ.  0) THEN
                     IF (ATEMP1 .GT. 0.D0) THEN
                        ISTOCH = 6
                     ELSE
                        ISTOCH = 5
                     ENDIF
                  ENDIF
!
!     If we have a univariate chance constraint, locate it.
!
               ELSEIF (ISTOCH .EQ. 8) THEN
                  DO I=1,LCHANCE
                     IF (XCHANCE(I)%SIZE .EQ. 1) THEN
                        J = XCHANCE(I)%SIZE
                        IF (RCHANCE(J)%RELROW .EQ. IRPOS .AND.      &
     &                      RCHANCE(J)%STAGE  .EQ. IRP  ) THEN
                           IRPOS = I - KCHANCE(IRP)
                           RETURN
                        ENDIF
                     ENDIF
                  ENDDO
!
!     Not found!
!
                  WRITE (IOLOG, 1100)      &
     &               NREC,Q1,Q2,Q3,Q4,DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                  WRITE (IOLOG, 1000) DROW
                  ISTOCH = 0
                  RETURN
!
!     If we have an integrated chance constraint, locate it.
!
               ELSEIF (ISTOCH .EQ. 9) THEN
                  DO I=1,LICC
                     IF (XICC(I)%RELROW .EQ. IRPOS .AND.      &
     &                   XICC(I)%NODE   .EQ. IRP  ) THEN
                        IRPOS = I - KICC(IRP)
                        RETURN
                     ENDIF
                  ENDDO
!
!     Not found!
!
                  WRITE (IOLOG, 1100)      &
     &               NREC,Q1,Q2,Q3,Q4,DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                  WRITE (IOLOG, 1000) DROW
                  ISTOCH = 0
                  RETURN
               ENDIF
!
!     First name field is a variable name. Identify.
!     ==============================================
!
         ELSE
            DO ICP=1,NPER
               JND = KPATH(ICP)
               KN  = NODEINFO(JND)%KNAMES
               NC  = NODEINFO(JND)%NCOL
               ICLOC = NMSRCH(DCOL,KN+1,KN+NC)
               IF (ICLOC .GT. 0) GOTO 130
            END DO
!
!     Not matched!
!
            WRITE (IOLOG, 1100)      &
     &         NREC,Q1,Q2,Q3,Q4,DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
            WRITE (IOLOG, 1000) DCOL
            ISTOCH = 0
            RETURN
!
!     Matched. Identify the name in the second name field.
!     (The second name field can be blank for supply/demand cards
!     or contain the tag 'PLQUAD' for the stochastic curvature of
!     a piecewise linear-quadratic penalty.)
!
  130 CONTINUE
         IF (DROW .EQ. DBLANK) THEN
            IF     (Q2 .EQ. QD .AND. Q3 .EQ. QE) THEN
               ISTOCH  = 12
               IBLOCK  = ICP
               IRPOS   = ICLOC - NODEINFO(JND)%KNAMES
               IF (ISTAGE .EQ. 0) THEN
                  ISTAGE = ICP
               ELSEIF (ISTAGE .LT. ICP) THEN
                  WRITE (IOLOG, 2000)
                  ISTOCH = 0
                  RETURN
               ENDIF
               IF (IRPOS .GT. NODEINFO(JND)%NROW) THEN
                  WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                                DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                  WRITE (IOLOG, 1300) DCOL
                  ISTOCH = 0
                  RETURN
               ENDIF
            ELSEIF (Q2 .EQ. QS .AND. Q3 .EQ. QU) THEN
               ISTOCH  = 13
               IBLOCK  = ICP
               IRPOS   = ICLOC - NODEINFO(JND)%KNAMES
               IF (ISTAGE .EQ. 0) THEN
                  ISTAGE = ICP
               ELSEIF (ISTAGE .LT. ICP) THEN
                  WRITE (IOLOG, 2000)
                  ISTOCH = 0
                  RETURN
               ENDIF
               IF (IRPOS .GT. NODEINFO(JND)%NROW) THEN
                  WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                                DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                  WRITE (IOLOG, 1300) DCOL
                  ISTOCH = 0
                  RETURN
               ENDIF
            ELSE
               WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                             DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
               WRITE (IOLOG, 1200) Q2,Q3
               ISTOCH = 0
               RETURN
            ENDIF
!
!     Second name field contains a tag for piecewise linear-quadratic penalty.
!
         ELSEIF (DROW .EQ. DPLQUAD) THEN
            IRPOS = ICLOC - NODEINFO(JND)%KNAMES - NODEINFO(JND)%NROW
            IF (IRPOS .LE. 0) THEN
               WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                  DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
               WRITE (IOLOG, 1400) DROW
               ISTOCH = 0
            ELSE
               ISTOCH = 10
               IBLOCK = ICP
               IF (ISTAGE .EQ. 0) THEN
                  ISTAGE = IBLOCK
               ELSEIF (ISTAGE .LT. IBLOCK) THEN
                  WRITE (IOLOG, 2000)
                  ISTOCH = 0
               ENDIF
            ENDIF
            RETURN
!
!     Second name field contains a name. Identify.
!
         ELSE
            DO IRP=1,NPER
               JND = KPATH(IRP)
               KN  = NODEINFO(JND)%KNAMES
               NC  = NODEINFO(JND)%NCOL
               IRLOC = NMSRCH(DROW,KN+1,KN+NC)
               IF (IRLOC .GT. 0) GOTO 150
            END DO
!
!     No match for second name field! If the card does not contain a
!     numeric field and has 'DE' or 'SU' code, then it is possible that
!     the second name field contains the name of a stage. Make sure.
!
            IF ( (Q2 .EQ. QS .AND. Q3 .EQ. QU) .OR.      &
     &           (Q2 .EQ. QD .AND. Q3 .EQ. QE) ) THEN
               DTPER = DROW
               DROW  = DBLANK
               DO ISTAGE=1,NPER
                  IF (DTPER .EQ. DTIME(ISTAGE)) GOTO 140
               END DO
               ISTAGE = 0
            ENDIF
            WRITE (IOLOG, 1100)      &
     &         NREC,Q1,Q2,Q3,Q4,DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
            WRITE (IOLOG, 1000) DROW
            ISTOCH = 0
            RETURN
!
!     Store the information
!
  140       CONTINUE
               IBLOCK = ICP
               IF (ISTAGE .EQ. 0) THEN
                  ISTAGE = ICP
               ELSEIF (ISTAGE .LT. ICP) THEN
                  WRITE (IOLOG, 2000)
                  ISTOCH = 0
                  RETURN
               ENDIF
               IRPOS  = ICLOC - NODEINFO(JND)%KNAMES
               IF (IRPOS .GT. NODEINFO(JND)%NROW) THEN
                  WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                                DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                  WRITE (IOLOG, 1300) DCOL
                  ISTOCH = 0
                  RETURN
               ENDIF
               IF (Q2 .EQ. QD .AND. Q3 .EQ. QE) THEN
                  ISTOCH = 12
               ELSE
                  ISTOCH = 13
               ENDIF
               RETURN
!
!     Got it!
!
  150    CONTINUE
            NCTEMP = NODEINFO(KPATH(ICP))%KNAMES      &
     &             + NODEINFO(KPATH(ICP))%NROW
            NRTEMP = NODEINFO(KPATH(IRP))%KNAMES      &
     &             + NODEINFO(KPATH(IRP))%NROW
            IF (ICLOC .LE. NCTEMP) THEN
               IF (IRLOC .LE. NRTEMP) THEN
!
!     Two row names: Must be a network card. Build the column name:
!     "{FROM-node}->{TO-node}[:<code>]" and identify.
!
                  DT = DBLANK
                  LEN1 = NAME1(ICLOC+1) - NAME1(ICLOC)
                  DO I=1,LEN1
                     QTEMP(I) = QNAMES(NAME1(ICLOC)+I-1)
                  END DO
                  QTEMP(LEN1+1) = '-'
                  QTEMP(LEN1+2) = '>'
                  LEN2 = NAME1(IRLOC+1) - NAME1(IRLOC)
                  DO I=1,LEN2
                     QTEMP(LEN1+2+I) = QNAMES(NAME1(IRLOC)+I-1)
                  END DO
                  IF (Q3 .NE. QBL) THEN
                     QTEMP(LEN1+LEN2+3) = ':'
                     QTEMP(LEN1+LEN2+4) = Q3
                  ENDIF
!
                  DO ICP=1,NPER
                     JND = KPATH(ICP)
                     KN  = NODEINFO(JND)%KNAMES
                     NC  = NODEINFO(JND)%NCOL
                     NR  = NODEINFO(JND)%NROW
                     IRES = NMSRCH(DT,KN+NR+1,KN+NC)
                     IF (IRES .GT. 0) GOTO 155
                  END DO
!
!     No match!
!
                  WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                   DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                  WRITE (IOLOG, 1000) DROW
                  ISTOCH = 0
                  RETURN
!
!     Found it! Check the code for more.
!
  155          CONTINUE
                  DT = DNAME3
                  IF (Q2 .EQ. QM) THEN
!
!     Multiplier, defines a location in the A-matrix
!
                     ISTOCH = 1
                     IROW = IRLOC - NODEINFO(KPATH(IRP))%KNAMES
                     IF (IRP .GE. ICP) THEN
                        IBLOCK = KDATA(IRP) + IRP - ICP + 1
                        JBLOCK = IBLOCK - KDATA(IRP) + KDATA(KREF(IRP))
                        IF (IRP-ICP .GT. 1 .AND. MARKOV) GOTO 284
                     ELSEIF (MARKOV) THEN
                        IBLOCK = KDATA(ICP) + IRP + MIN(2,ICP)
                        JBLOCK = IBLOCK - KDATA(ICP) + KDATA(KREF(ICP))
                        IF (.NOT. HAVE_GC) GOTO 284
                     ELSE
                        IBLOCK = KDATA(ICP) + IRP + ICP
                        JBLOCK = IBLOCK - KDATA(ICP) + KDATA(KREF(ICP))
                        IF (.NOT. HAVE_GC) GOTO 284
                     ENDIF
!
!     For rows in chance constraints adjust the stage information if not given explicitly.
!
                     IF (ISTAGE .EQ. 0) THEN
                        IF (LCHANCE .GT. 0) THEN
                           DO JRP=IRP,NPER
                              DO IC=1,NCHANCE(JRP)
                                 JCR = XCHANCE(IC)%ROWPTR
                                 DO ICR=1,XCHANCE(IC)%SIZE
                                    TR = RCHANCE(ICR+JCR)
                                    IF (TR%STAGE  .EQ. IRP   .AND.      &
     &                                  TR%RELROW .EQ. IRPOS) THEN
                                       ISTAGE = JRP + 1
                                       GOTO 156
                                    END IF
                                 END DO
                              END DO
                           END DO
                        ENDIF
!
                        IF (LICC .GT. 0) THEN
                           DO IC=1,NICC(IRP)
                              IF (XICC(IC)%RELROW .EQ. IRPOS) THEN
                                 ISTAGE = IRP + 1
                                 GOTO 156
                              END IF
                           END DO
                        ENDIF
                     ENDIF
!
  156             CONTINUE
                     II = IRES -  NODEINFO(KPATH(ICP))%KNAMES      &
     &                         -  NODEINFO(KPATH(ICP))%NROW
                     LL = KELMA(JBLOCK) + LA(KCOLA(JBLOCK)+II)
                     KK = KELMA(JBLOCK) + LA(KCOLA(JBLOCK)+II+1) - 1
                     DO IR=LL,KK
                        IF (IA(IR) .EQ. IROW) THEN
                           IRPOS = IR - KELMA(JBLOCK)
                           RETURN
                        ENDIF
                     END DO
!
!     Coefficient not defined in core file! Check that the structure is
!     right, at least. If not, issue a warning: this looks like an error.
!
                     WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                      DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                     IF (((IRP .LT. ICP)     .AND. .NOT. HAVE_GC) .OR.      &
     &                   ((IRP .GT. ICP + 1) .AND.       MARKOV)) THEN
                        WRITE (IOLOG, 1700) DROW,DCOL
                        ISTOCH = 0
                        RETURN
                     ELSE
                        WRITE (IOLOG, 1600) DROW,DCOL
                        ISTOCH = -1
                        IRPOS  = KK + 1
                        ICOL   = II
                        RETURN
                     ENDIF
!
!     Cost coefficient
!
                  ELSEIF (Q2 .EQ. QC) THEN
                     JND = KPATH(ICP)
                     IRPOS  = IRES - NODEINFO(JND)%KNAMES      &
     &                             - NODEINFO(JND)%NROW
                     IBLOCK = ICP
                     ISTOCH = 2
                     IF (ISTAGE .EQ. 0) THEN
                        ISTAGE = ICP
                     ELSEIF (ISTAGE .LT. ICP) THEN
                        WRITE (IOLOG, 2000)
                        ISTOCH = 0
                     ENDIF
                     RETURN
!
!     Random bound
!
                  ELSEIF (Q2 .EQ. QL .OR. Q2 .EQ. QU      &
     &                               .OR. Q2 .EQ. QX) THEN
                     IRPOS  = IRES - NODEINFO(KPATH(ICP))%KNAMES
                     IBLOCK = ICP
                     IF (Q2     .EQ. QL) ISTOCH = 5
                     IF (Q2     .EQ. QU) ISTOCH = 6
                     IF (Q2     .EQ. QX) ISTOCH = 11
                     IF (ISTAGE .EQ.  0) THEN
                        ISTAGE = ICP
                     ELSEIF (ISTAGE .LT. ICP) THEN
                        WRITE (IOLOG, 2000)
                        ISTOCH = 0
                     ENDIF
                     RETURN
                  ELSEIF (DNAME3 .NE. DBLANK) THEN
!
!     Could be a side constraint if third name field is a valid row
!
                     DO IRP=1,NPER
                        JND = KPATH(IRP)
                        KN  = NODEINFO(JND)%KNAMES
                        NR  = NODEINFO(JND)%NROW
                        IROW = NMSRCH(DT,KN+1,KN+NR)
                        IF (IROW .GT. 0) GOTO 164
                     END DO
!
!     No match!
!
                     WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                      DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                     WRITE (IOLOG, 1000) DROW
                     ISTOCH = 0
                     RETURN
!
!     A-matrix coefficient
!
  164             CONTINUE
                     IROW =  IROW - NODEINFO(KPATH(IRP))%KNAMES
                     ISTOCH = 1
                     IF (IRP .GE. ICP) THEN
                        IBLOCK = KDATA(IRP) + IRP - ICP + 1
                        JBLOCK = IBLOCK - KDATA(IRP) + KDATA(KREF(IRP))
                        IF (IRP-ICP .GT. 1 .AND. MARKOV) GOTO 284
                     ELSEIF (MARKOV) THEN
                        IBLOCK = KDATA(ICP) + IRP + MIN(2,ICP)
                        JBLOCK = IBLOCK - KDATA(ICP) + KDATA(KREF(ICP))
                        IF (.NOT. HAVE_GC) GOTO 284
                     ELSE
                        IBLOCK = KDATA(ICP) + IRP + ICP
                        JBLOCK = IBLOCK - KDATA(ICP) + KDATA(KREF(ICP))
                        IF (.NOT. HAVE_GC) GOTO 284
                     ENDIF
!
!     For rows in chance constraints adjust the stage information if not given explicitly.
!
                     IF (ISTAGE .EQ. 0) THEN
                        IF (LCHANCE .GT. 0) THEN
                           DO JRP=IRP,NPER
                              DO IC=1,NCHANCE(JRP)
                                 JCR = XCHANCE(IC)%ROWPTR
                                 DO ICR=1,XCHANCE(IC)%SIZE
                                    TR = RCHANCE(ICR+JCR)
                                    IF (TR%STAGE  .EQ. ICP   .AND.      &
     &                                  TR%RELROW .EQ. IRPOS) THEN
                                       ISTAGE = JRP + 1
                                       GOTO 165
                                    END IF
                                 END DO
                              END DO
                           END DO
                        ENDIF
!
                        IF (LICC .GT. 0) THEN
                           DO IC=1,NICC(IRP)
                              IF (XICC(IC)%RELROW .EQ. IRPOS) THEN
                                 ISTAGE = IRP + 1
                                 GOTO 165
                              END IF
                           END DO
                        ENDIF
                     ENDIF
!
  165             CONTINUE
                     IF (ISTAGE .EQ. 0) THEN
                        ISTAGE = ICP
                     ELSEIF (ISTAGE .LT. ICP) THEN
                        WRITE (IOLOG, 2000)
                        ISTOCH = 0
                        RETURN
                     ENDIF
                     II = IRES - NODEINFO(KPATH(ICP))%KNAMES      &
     &                         - NODEINFO(KPATH(ICP))%NROW
                     LL = KELMA(JBLOCK) + LA(KCOLA(JBLOCK)+II)
                     KK = KELMA(JBLOCK) + LA(KCOLA(JBLOCK)+II+1) - 1
                     DO IR=LL,KK
                        IF (IA(IR) .EQ. IROW) THEN
                           IRPOS = IR - KELMA(JBLOCK)
                           RETURN
                        ENDIF
                     END DO
!
!     Coefficient not defined in core file!
!
                     WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                      DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                     IF (((IRP .LT. ICP)     .AND. .NOT. HAVE_GC) .OR.      &
     &                   ((IRP .GT. ICP + 1) .AND.       MARKOV)) THEN
                        WRITE (IOLOG, 1700) DROW,DCOL
                        ISTOCH = 0
                        RETURN
                     ELSE
                        WRITE (IOLOG, 1600) DROW,DCOL
                        ISTOCH = -1
                        IRPOS  = KK + 1
                        ICOL   = II
                        RETURN
                     ENDIF
                  ELSE
!
!     Code not appropriate for this situation
!
                     WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                             DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                     WRITE (IOLOG, 1200) Q2,Q3
                     ISTOCH = 0
                     RETURN
                  ENDIF
               ELSE
!
!     First field is a row, second is a column; this can't happen!
!
                  WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                                DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                  WRITE (IOLOG, 1500) DROW,DCOL
                  ISTOCH = 0
                  RETURN
               ENDIF
            ELSE
!
!     First field is a column name. Here it could be an entry in the
!     A or Q matrix, or a cost coefficient.
!
               IF (IRLOC .EQ. 1) THEN
!
!     Stochastic cost coefficient
!
                  ISTOCH = 2
                  IBLOCK = ICP
                  IRPOS  = ICLOC - NODEINFO((KPATH(ICP)))%KNAMES      &
     &                           - NODEINFO((KPATH(ICP)))%NROW
                  IF (ISTAGE .EQ. 0) THEN
                     ISTAGE = ICP
                  ELSEIF (ISTAGE .LT. ICP) THEN
                     WRITE (IOLOG, 2000)
                     ISTOCH = 0
                     RETURN
                  ENDIF
               ELSEIF (IRLOC .LE.  NODEINFO((KPATH(IRP)))%KNAMES      &
     &                           + NODEINFO((KPATH(IRP)))%NROW) THEN
!
!     A-matrix coefficient
!
                  IROW =  IRLOC - NODEINFO(JND)%KNAMES
                  ISTOCH = 1
                  IF (IRP .GE. ICP) THEN
                     IBLOCK = KDATA(IRP) + IRP - ICP + 1
                     JBLOCK = IBLOCK - KDATA(IRP) + KDATA(KREF(IRP))
                     IF (IRP-ICP .GT. 1 .AND. MARKOV) GOTO 284
                        JSTAGE = IRP
                  ELSEIF (MARKOV) THEN
                     IBLOCK = KDATA(ICP) + IRP + MIN(2,ICP)
                     JBLOCK = IBLOCK - KDATA(ICP) + KDATA(KREF(ICP))
                     IF (.NOT. HAVE_GC) GOTO 284
                        JSTAGE = ICP
                  ELSE
                     IBLOCK = KDATA(ICP) + IRP + ICP
                     JBLOCK = IBLOCK - KDATA(ICP) + KDATA(KREF(ICP))
                     IF (.NOT. HAVE_GC) GOTO 284
                        JSTAGE = ICP
                  ENDIF
!
!     For rows in chance constraints adjust the stage information if not given explicitly.
!
                  IF (ISTAGE .EQ. 0) THEN
                     IF (LCHANCE .GT. 0) THEN
                        DO JRP=IRP,NPER
                           DO IC=1,NCHANCE(JRP)
                              JCR = XCHANCE(IC)%ROWPTR
                              DO ICR=1,XCHANCE(IC)%SIZE
                                 TR = RCHANCE(ICR+JCR)
                                 IF (TR%STAGE  .EQ. JSTAGE .AND.      &
     &                               TR%RELROW .EQ. IRPOS) THEN
                                    ISTAGE = MAX(JRP,JSTAGE) + 1
                                    GOTO 168
                                 END IF
                              END DO
                           END DO
                        END DO
                     ENDIF
!
                     IF (LICC .GT. 0) THEN
                        DO IC=1,NICC(IRP)
                           IF (XICC(IC)%RELROW .EQ. IRPOS) THEN
                              ISTAGE = JSTAGE + 1
                              GOTO 168
                           END IF
                        END DO
                     ENDIF
                  ENDIF
!
  168          CONTINUE
                  IF (ISTAGE .EQ. 0) THEN
                     ISTAGE = JSTAGE
                  ELSEIF (ISTAGE .LT. JSTAGE) THEN
                     WRITE (IOLOG, 2000)
                     ISTOCH = 0
                     RETURN
                  ENDIF
                  II = ICLOC - NODEINFO(KPATH(ICP))%KNAMES      &
     &                       - NODEINFO(KPATH(ICP))%NROW
                  LL = KELMA(JBLOCK) + LA(KCOLA(JBLOCK)+II)
                  KK = KELMA(JBLOCK) + LA(KCOLA(JBLOCK)+II+1) - 1
                  DO IR=LL,KK
                     IF (IA(IR) .EQ. IROW) THEN
                        IRPOS = IR - KELMA(JBLOCK)
                        RETURN
                     ENDIF
                  END DO
!
!     Coefficient not defined in core file!
!
                  WRITE (IOLOG, 1100)      &
     &               NREC,Q1,Q2,Q3,Q4,DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                     IF (((IRP .LT. ICP)     .AND. .NOT. HAVE_GC) .OR.      &
     &                   ((IRP .GT. ICP + 1) .AND.       MARKOV)) THEN
                     WRITE (IOLOG, 1700) DROW,DCOL
                     ISTOCH = 0
                     RETURN
                  ELSE
                     WRITE (IOLOG, 1600) DROW,DCOL
                     ISTOCH = -1
                     IRPOS  = KK + 1
                     ICOL   = II
                     RETURN
                  ENDIF
!
!     Two column names: Q-matrix coefficient
!
               ELSE
                  IROW = IRLOC - NODEINFO(KPATH(IRP))%KNAMES      &
     &                         - NODEINFO(KPATH(IRP))%NROW
                  ISTOCH = 7
                  IF (IRP .GE. ICP) THEN
                     IF (IRP-ICP .GE. NQPROF) GOTO 285
                        IBLOCK = KDATQ(IRP) + IRP - ICP + 1
                        JBLOCK = IBLOCK - KDATQ(IRP) + KDATQ(KREF(IRP))
                        IF (ISTAGE .EQ. 0) THEN
                           ISTAGE = IRP
                        ELSEIF (ISTAGE .LT. IRP) THEN
                           WRITE (IOLOG, 2000)
                           ISTOCH = 0
                           RETURN
                        ENDIF
                  ELSE
                     WRITE (IOLOG, 1100) NREC,Q1,Q2,Q3,Q4,      &
     &                      DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                     WRITE (IOLOG, 1900)
                     ISTOCH = 99
                     RETURN
                  ENDIF
                  II = ICLOC - NODEINFO(KPATH(ICP))%KNAMES      &
     &                       - NODEINFO(KPATH(ICP))%NROW
                  LL = IQOFF(JBLOCK) + LQMTX(LQOFF(JBLOCK)+II)
                  KK = IQOFF(JBLOCK) + LQMTX(LQOFF(JBLOCK)+II+1) - 1
                  DO IR=LL,KK
                     IF (IQMTX(IR) .EQ. IROW) THEN
                        IRPOS = IR - IQOFF(JBLOCK)
                        RETURN
                     ENDIF
                  END DO
!
!     Coefficient not defined in core file!
!
                  WRITE (IOLOG, 1100)      &
     &               NREC,Q1,Q2,Q3,Q4,DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                  WRITE (IOLOG, 1600) DROW,DCOL
                  ISTOCH = -7
                  IRPOS  = KK + 1
                  ICOL   = II
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      RETURN
!
!     A-Coefficient is in a block not previously defined. This is illegal.
!
  284          CONTINUE
                  WRITE (IOLOG, 1100)      &
     &               NREC,Q1,Q2,Q3,Q4,DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                  WRITE (IOLOG, 1700) DROW,DCOL
                  ISTOCH = 0
                  RETURN
!
!     Q-Coefficient is in a block not previously defined. This is illegal.
!
  285          CONTINUE
                  WRITE (IOLOG, 1100)      &
     &               NREC,Q1,Q2,Q3,Q4,DNAME1,DNAME2,ATEMP1,DNAME3,ATEMP2
                  IF (NQPROF .EQ. 0) THEN
                     WRITE (IOLOG, 1800)
                  ELSE
                     WRITE (IOLOG, 1700) DROW,DCOL
                  ENDIF
                  ISTOCH = 0
                  RETURN
!
 1000 FORMAT(' XXX - WARNING - Variable name not matched: ',A30)
 1100 FORMAT(I8,4X,4A1,A8,2X,A8,2X,F12.4,3X,A8,F12.4)
 1200 FORMAT(' XXX - WARNING - Inappropriate code detected: ',2A1)
 1300 FORMAT(' XXX - WARNING - Column name where row expected: ',A30)
 1400 FORMAT(' XXX - WARNING - Row name where column expected: ',A30)
 1500 FORMAT(' XXX - WARNING - Row and column names in wrong order:',      &
     &       /,2X,A20,2X,A20)
 1600 FORMAT(' XXX - WARNING - Location not defined in core file:',      &
     &       /,2X,A20,2X,A20)
 1700 FORMAT(' XXX -  FATAL  - Block not defined in core file for:',      &
     &       /,2X,A20,2X,A20)
 1800 FORMAT(' XXX -  FATAL  - Q-matrix was never read from core file!')
 1900 FORMAT(' XXX - WARNING - Q-matrix element in upper triangle.'      &
     &       ' Ignore.')
 2000 FORMAT(' XXX - WARNING - Inconsistent timing information.',      &
     &       ' Record ignored.')
      END SUBROUTINE
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE IDROW(DNAME,NR,IRP)
!
!     This routine takes as input a name field from a record in a
!     chance-constraint segment and checks if the name constitutes
!     a valid row. If yes, it returns the row number in variable NR
!     (in absolute row address form) and the period (for consistency
!     checking) in variable IRP.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME
!
      NR = 0
      DO I=1,NPER
         KN = NODEINFO(I)%KNAMES
         NN = NODEINFO(I)%NROW
         NR = NMSRCH(DNAME,KN+1,KN+NN)
         IF (NR .GT. 0) THEN
            IRP = I
            RETURN
         ENDIF
      END DO
      WRITE (IOLOG, 1000) DNAME
      NR  = 0
      IRP = 0
      RETURN
!
 1000 FORMAT(' XXX - WARNING - Variable name ',A8,' not matched')
      END SUBROUTINE
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      FUNCTION NMSRCH(S1,N1,N2)
!
!     This function searches the name strings in positions N1 through N2
!     for the occurrence of target string S1. If a match is encountered,
!     NMSRCH is set to the position where the match occurred, otherwise
!     NMSRCH returns 0. The logical variable FIX switches between two
!     search modes: If FIX is true, all name fields are eight characters
!     long, and the search is considerably simplified. Otherwise the
!     name fields may have variable lengths and must be copied first.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) S1,S2,EMPTY
      CHARACTER*1 ST(NAMLEN),SBLANK(NAMLEN)
!
      EQUIVALENCE (S2,ST), (SBLANK,EMPTY)
!
      DATA SBLANK /255*' '/
!
      DO I=N1,N2
         S2 = EMPTY
         ST(1:NAME1(I+1)-NAME1(I)) = QNAMES(NAME1(I):NAME1(I+1)-1)
         IF (S1 .EQ. S2) GOTO 130
      END DO
!
!     No match
!
      NMSRCH = 0
      RETURN
!
!     Matched
!
  130 CONTINUE
      NMSRCH = I
      RETURN
!
      END FUNCTION
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE LINKTREE
!
!     This routine links together the nodes in the event tree and
!     counts the number of descendants.
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      DO 120 IP=1,NPER
         ISC1 = IRNGE0(IP)
  100    CONTINUE
            IF (NODEINFO(ISC1)%IBROTH .GT. 0) THEN
               ISC1 = NODEINFO(ISC1)%IBROTH
               GOTO 100
            ENDIF
!
            IAN = NODEINFO(ISC1)%IANCTR
            IF (IAN .EQ. 0) GOTO 120
  110       CONTINUE
               IBRO = ABS(NODEINFO(IAN)%IBROTH)
               IF (IBRO .EQ. 0) GOTO 120
                  ISC2 = NODEINFO(IBRO)%IDESC
                  IF (ISC2 .GT. 0) THEN
                     NODEINFO(ISC1)%IBROTH = -ISC2
                     ISC1 = ISC2
                     GOTO 100
                  ENDIF
!
                  IAN = IBRO
                  GOTO 110
  120 CONTINUE
!
!     Count descendants for each problem
!
      DO I=1,NODES
         N0 = 0
         I0 = NODEINFO(I)%IDESC
  130    CONTINUE
            IF (I0 .GT. 0) THEN
               I0 = NODEINFO(I0)%IBROTH
               N0 = N0 + 1
               GOTO 130
            ENDIF
            NODEINFO(I)%NDESC = N0
         END DO
!
      RETURN
      END SUBROUTINE
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE GETMPS(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                                    IERR,IONUM,NREC,LENGTH)
!
!     This routine reads one line of an MPS file and extracts from
!     each record a 4-character key, up to three non-numeric modifiers,
!     and one or two numeric fields (real numbers). If the line contains
!     no information (flagged by a comment marker in the code field),
!     another line is read.
!
!     Adapted from MINOS routine M3SPC1 and ob1 routine hspec1.
!
!     Assumptions about the MPS file:
!
!     1. Only code information in the first four character positions.
!     2. The sequence of elements is fixed: Each line contains four
!        character keys, two name fields, one numeric field, possibly
!        followed by another name and number and a comment, in that order.
!     3. Cards whose first non-blank character is an asterisk or an
!        exclamation mark are considered comment cards.
!     4. Data items are separated by at least one space.
!     5. Real numbers in scientific notation have no imbedded blanks
!     6. Non-numeric strings have no imbedded blanks. The first
!        character is not a digit, nor does the string start with
!        a period followed by a digit. (BP47, F.04, and .COSTA will
!        all be parsed correctly, but 3BL and .5COST would be illegal.)
!     7. The line length is limited to 255 characters.
!
!     ==================================================================
!
!     This version dated 5 November 2002
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*1 QLOWER(26),QUPPER(26),QDIGIT(12),      &
     &            QPUNCT(4),QNUM(NAMLEN),QNAME(NAMLEN),QB1(NAMLEN),      &
     &            QCOPY(NAMLEN)
      CHARACTER*4 QCODE
      CHARACTER*(NAMLEN) QNUM1,DNAME(*),QNAME1,QBLANK,CCOPY,QLINE
      LOGICAL*4   REALNO,SIGNED,EXPONT,SIGNEX
      REAL*8      ATEMP(*)
      INTEGER*4   LENGTH(*)
!
      EQUIVALENCE (QNUM1,QNUM),(QNAME,QNAME1),(QBLANK,QB1),(QCOPY,CCOPY)
!
      DATA QLOWER/'a','b','c','d','e','f','g','h','i','j','k','l','m',      &
     &            'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      DATA QUPPER/'A','B','C','D','E','F','G','H','I','J','K','L','M',      &
     &            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      DATA QDIGIT/'1','2','3','4','5','6','7','8','9','0','+','-'/
      DATA QPUNCT/'!','*','.','='/
      DATA QB1 /NAMLEN*' '/
!
  100 CONTINUE
      READ (IONUM, 1000, END=910, ERR=920) QCOPY
      NREC = NREC + 1
      IERR = 0
      QLINE = CCOPY
!
!     If the file is not free format, try to read the line in fixed format.
!     Switch to free format if there is an error.
!
      IF (.NOT. FREEFORM) THEN
         READ (CCOPY, 1300, ERR=110) Q1,Q2,Q3,Q4,DNAME(1),DNAME(2),      &
     &                                  ATEMP(1),DNAME(3),ATEMP(2)
!
!     Check for trivial comment line
!
      IF (Q1 .EQ. QPUNCT(1) .OR. Q1 .EQ. QPUNCT(2) ) GOTO 100
         QCODE = Q1 // Q2 // Q3 //Q4
         IF (QCODE .EQ. ' !'     .OR.      &
     &       QCODE .EQ. '  !'    .OR.      &
     &       QCODE .EQ. '   !' ) GOTO 100
!
         LENGTH(1) = 8
         LENGTH(2) = 8
         LENGTH(3) = 8
         RETURN
      ENDIF
!
!     Read the line in free format. Access the fields one at a time.
!
  110 CONTINUE
      FREEFORM = .TRUE.
      NNAME = 0
      NNUM  = 0
      ATEMP(1) = 0.D0
      ATEMP(2) = 0.D0
      DNAME(1) = QBLANK
      DNAME(2) = QBLANK
      DNAME(3) = QBLANK
      LENGTH(1) = 0
      LENGTH(2) = 0
      LENGTH(3) = 0
!
      Q1 = QCOPY(1)
      Q2 = QCOPY(2)
      Q3 = QCOPY(3)
      Q4 = QCOPY(4)
!
!     Check for trivial comment line
!
      IF (Q1 .EQ. QPUNCT(1) .OR. Q1 .EQ. QPUNCT(2) ) GOTO 100
         QCODE = Q1 // Q2 // Q3 //Q4
         IF (QCODE .EQ. ' !'     .OR.      &
     &       QCODE .EQ. '  !'    .OR.      &
     &       QCODE .EQ. '   !' ) GOTO 100
!
         I = 4
!
!     Skip all blanks to the next item -- if there is one
!
  190    CONTINUE
            J = I
            DO I=J+1,NAMLEN
               IF (QCOPY(I) .NE. QBL) GOTO 210
            END DO
            GOTO 400
!
!     Ignore comments
!
  210 CONTINUE
         IF ( QCOPY(I) .EQ. QPUNCT(1)      &
     &                 .OR. QCOPY(I) .EQ. QPUNCT(2)) GOTO 400
!
!     Check for digit
!
         DO JDIGIT=1,10
            IF (QCOPY(I) .EQ. QDIGIT(JDIGIT)) GOTO 300
         END DO
         IF (QCOPY(I) .EQ. QDIGIT(11)      &
     &                .OR. QCOPY(I) .EQ. QDIGIT(12)) GOTO 300
         IF (QCOPY(I) .NE. QPUNCT(3)) GOTO 240
!
!     First character is a period. Must be a number if the next is a digit
!
         DO JDIGIT=1,10
            IF (QCOPY(I+1) .EQ. QDIGIT(JDIGIT)) GOTO 300
         END DO
!
!     Not a digit - treat as a name
!
  240    CONTINUE
            NNAME = NNAME + 1
            IF (NNAME .EQ. 3) NNUM = 1
            IF (NNAME .EQ. 4) GOTO 400
            QNAME1 = QBLANK
            DO J=I,NAMLEN
               QCHAR = QCOPY(J)
               IF (QCHAR .EQ. QBL) GOTO 260
                  QNAME(J+1-I) = QCHAR
            END DO
            J = NAMLEN + 1
!
  260    CONTINUE
            IF (FREEFORM) THEN
               LENGTH(NNAME) = J - I
            ELSEIF ((J-I) .LE. 8) THEN
               LENGTH(NNAME) = 8
            ELSE
               FREEFORM = .TRUE.
               WRITE (IOLOG, 1600)
               LENGTH(NNAME) = J - I
            ENDIF
            I = J
            READ (QNAME1, 1200) DNAME(NNAME)
            IF (NNUM .LT. 2) GOTO 190
               GOTO 400
!
!     Looks like a number...
!
  300    CONTINUE
            NNUM = NNUM + 1
            IF (NNUM .EQ. 2) NNAME = 3
            QNUM1 = QBLANK
!
            REALNO = .FALSE.
            SIGNED = .FALSE.
            EXPONT = .FALSE.
            SIGNEX = .FALSE.
            DO 360 J=I,NAMLEN
               QCHAR = QCOPY(J)
               IF (QCHAR .EQ. QBL) GOTO 370
               DO JCHAR=1,10
                  IF (QCHAR .EQ. QDIGIT(JCHAR)) GOTO 360
               END DO
               IF (QCHAR  .EQ. QDIGIT(11)      &
     &                    .OR. QCHAR .EQ. QDIGIT(12)) GOTO 340
               IF (QCHAR  .EQ. QUPPER( 5)      &
     &                    .OR. QCHAR .EQ. QLOWER( 5)) GOTO 350
               IF (QCHAR  .EQ. QUPPER( 4)      &
     &                    .OR. QCHAR .EQ. QLOWER( 4)) GOTO 350
               IF (QCHAR  .NE. QPUNCT( 3)) GOTO 930
               IF (REALNO .OR. EXPONT) GOTO 930
                  REALNO = .TRUE.
                  GOTO 360
!
  340          CONTINUE
                  IF (.NOT. EXPONT) THEN
                     IF (SIGNED) GOTO 930
                         SIGNED = .TRUE.
                  ELSE
                     IF (SIGNEX) GOTO 930
                         SIGNEX = .TRUE.
                  ENDIF
                  GOTO 360
!
  350          CONTINUE
                  IF (EXPONT) GOTO 930
                      EXPONT = .TRUE.
                      SIGNED = .FALSE.
  360       CONTINUE
            J = NAMLEN + 1
!
!     Number has been read. Copy info
!
  370    CONTINUE
            DO K=I,J-1
               QNUM(K+1-I) = QCOPY(K)
            END DO
            IF (.NOT. REALNO .AND. .NOT. EXPONT)      &
     &         QNUM(J+1-I) = '.'
            I = J
            READ (QNUM1, 1100) ATEMP(NNUM)
            IF (NNUM .LT. 2) GOTO 190
!
  400 CONTINUE
         RETURN
!
!     END card is missing
!
  910 CONTINUE
         IERR = 1
         WRITE (IOLOG, 1700) IONUM
         RETURN
!
!     Error after start. This should not happen!
!
  920 CONTINUE
         IERR = 2
         WRITE (IOLOG, 1000) QCOPY
         WRITE (IOLOG, 1800) IONUM
         RETURN
!
!     Illegal format in a numeric field
!
  930 CONTINUE
         WRITE (IOLOG, 1900) QCOPY
         GOTO 190
!
 1000 FORMAT(255(A1))
 1100 FORMAT(G255.20)
 1200 FORMAT(A255)
 1300 FORMAT(4A1,A8,2X,A8,2X,F12.6,3X,A8,2X,F12.6)
 1600 FORMAT(' XXX - WARNING - Switching to free-format input.')
 1700 FORMAT(' XXX - WARNING - End of file detected on I/O',      &
     &       ' channel',i3)
 1800 FORMAT(' XXX -  FATAL  - System error in routine RDMPS,',      &
     &       ' attached to channel',i4)
 1900 FORMAT(' XXX - WARNING - Illegal format in numeric field.',      &
     &       ' Ignore this record:',/,80A1)
      END SUBROUTINE
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE GETNET(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                                    IERR,IONUM,NREC,LENGTH)
!
!     This routine reads one line of a NETGEN file and extracts from
!     each record a 4-character key, up to three non-numeric modifiers,
!     and up to four numeric fields (real numbers).
!
!     Adapted from MINOS routine M3SPC1 and ob1 routine hspec1.
!
!     Assumptions about the NETGEN format:
!
!     1. Only code information in the first four character positions.
!     2. The sequence of elements is fixed: Each line contains four
!        character keys, two name fields, one to four numeric fields,
!        possibly followed by another name and a comment, in that order.
!     3. Cards whose first non-blank character is an asterisk or an
!        exclamation mark are considered comment cards.
!     4. Data items are separated by at least one space.
!     5. Real numbers in scientific notation have no imbedded blanks
!     6. Non-numeric strings have no imbedded blanks. The first
!        character is not a digit, nor does the string start with
!        a period followed by a digit. (BP47, F.04, and .COSTA will
!        all be parsed correctly, but 3BL and .5COST would be illegal.)
!     7. The line length is limited to 255 characters.
!
!     ==================================================================
!
!     This version dated 3 October 1997
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) QNUM1,DNAME(*),QNAME1,QBLANK,CCOPY,QLINE
      CHARACTER*4 QCODE
      CHARACTER*1 QLOWER(26),QUPPER(26),QDIGIT(12),QCOPY(NAMLEN),      &
     &            QPUNCT(4),QNUM(NAMLEN),QNAME(NAMLEN),QB1(NAMLEN)
      LOGICAL*4   REALNO,SIGNED,EXPONT,SIGNEX
      REAL*8      ATEMP(*)
      INTEGER*4   LENGTH(*)
!
      EQUIVALENCE (QNUM1,QNUM),(QNAME,QNAME1),(QBLANK,QB1),(QCOPY,CCOPY)
!
      DATA QLOWER/'a','b','c','d','e','f','g','h','i','j','k','l','m',      &
     &            'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      DATA QUPPER/'A','B','C','D','E','F','G','H','I','J','K','L','M',      &
     &            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      DATA QDIGIT/'1','2','3','4','5','6','7','8','9','0','+','-'/
      DATA QPUNCT/'!','*','.','='/
      DATA QB1 /NAMLEN*' '/
!
  100 CONTINUE
      READ (IONUM, 1000, END=910, ERR=920) QCOPY
      NREC = NREC + 1
      IERR = 0
      QLINE = CCOPY
!
!     If the file is not free format, try to read the line in fixed format.
!     Switch to free format if there is an error.
!
      IF (.NOT. FREEFORM) THEN
         READ (CCOPY, 1300, ERR=110) Q1,Q2,Q3,Q4,DNAME(1),DNAME(2),      &
     &                ATEMP(1),ATEMP(2),ATEMP(3),ATEMP(4),DNAME(3)
!
!     Check for trivial comment line
!
      IF (Q1 .EQ. QPUNCT(1) .OR. Q1 .EQ. QPUNCT(2) ) GOTO 100
         QCODE = Q1 // Q2 // Q3 //Q4
         IF (QCODE .EQ. ' !'     .OR.      &
     &       QCODE .EQ. '  !'    .OR.      &
     &       QCODE .EQ. '   !' ) GOTO 100
!
         LENGTH(1) = 8
         LENGTH(2) = 8
         LENGTH(3) = 8
         RETURN
      ENDIF
!
!     Read the line in free format. Access the fields one at a time.
!
  110 CONTINUE
      FREEFORM = .TRUE.
      NNAME = 0
      NNUM  = 0
      ATEMP(1) = 0
      ATEMP(2) = 0
      ATEMP(3) = 0
      ATEMP(4) = 0
      DNAME(1) = QBLANK
      DNAME(2) = QBLANK
      DNAME(3) = QBLANK
      LENGTH(1) = 0
      LENGTH(2) = 0
      LENGTH(3) = 0
!
      Q1 = QCOPY(1)
      Q2 = QCOPY(2)
      Q3 = QCOPY(3)
      Q4 = QCOPY(4)
!
!     Check for trivial comment line
!
      IF (Q1 .EQ. QPUNCT(1) .OR. Q1 .EQ. QPUNCT(2) ) GOTO 100
         QCODE = Q1 // Q2 // Q3 //Q4
         IF (QCODE .EQ. ' !'     .OR.      &
     &       QCODE .EQ. '  !'    .OR.      &
     &       QCODE .EQ. '   !' ) GOTO 100
!
         I = 4
!
!     Skip all blanks to the next item -- if there is one
!
  190    CONTINUE
            J = I
            DO I=J+1,NAMLEN
               IF (QCOPY(I) .NE. QBL) GOTO 210
            END DO
            GOTO 400
!
!     Ignore comments
!
  210 CONTINUE
         IF ( QCOPY(I) .EQ. QPUNCT(1) .OR.      &
     &        QCOPY(I) .EQ. QPUNCT(2)) GOTO 400
!
!     Check for digit
!
         REALNO = .FALSE.
         SIGNED = .FALSE.
         EXPONT = .FALSE.
         SIGNEX = .FALSE.
!
         DO JDIGIT=1,10
            IF (QCOPY(I) .EQ. QDIGIT(JDIGIT)) GOTO 300
         END DO
         IF (QCOPY(I) .EQ. QDIGIT(11) .OR.      &
     &       QCOPY(I) .EQ. QDIGIT(12)) GOTO 300
         IF (QCOPY(I) .NE. QPUNCT(3))  GOTO 240
!
!     First character is a period. Must be a number if the next is a digit
!
         DO JDIGIT=1,10
            IF (QCOPY(I+1) .EQ. QDIGIT(JDIGIT)) GOTO 300
         END DO
!
!     Not a digit - treat as a name
!
  240    CONTINUE
            NNAME = NNAME + 1
            IF (NNAME .EQ. 3) NNUM = 4
            QNAME1 = QBLANK
            DO J=I,NAMLEN
               QCHAR = QCOPY(J)
               IF (QCHAR .EQ. QBL) GOTO 260
                  QNAME(J+1-I) = QCHAR
            END DO
            J = NAMLEN + 1
!
  260    CONTINUE
            IF (FREEFORM) THEN
               LENGTH(NNAME) = J - I
            ELSEIF ((J-I) .LE. 8) THEN
               LENGTH(NNAME) = 8
            ELSE
               FREEFORM = .TRUE.
               WRITE (IOLOG, 1600)
               LENGTH(NNAME) = J - I
            ENDIF
            I = J
            READ (QNAME1, 1200) DNAME(NNAME)
            IF (NNAME .LT. 3) GOTO 190
               GOTO 400
!
!     Looks like a number...
!
  300    CONTINUE
            NNUM = NNUM + 1
            IF (NNUM .EQ. 4) NNAME = 3
            QNUM1 = QBLANK
            DO 360 J=I,NAMLEN
               QCHAR = QCOPY(J)
               IF (QCHAR .EQ. QBL) GOTO 370
               DO JCHAR=1,10
                  IF (QCHAR .EQ. QDIGIT(JCHAR)) GOTO 360
               END DO
               IF (QCHAR  .EQ. QDIGIT(11)      &
     &                    .OR. QCHAR .EQ. QDIGIT(12)) GOTO 340
               IF (QCHAR  .EQ. QUPPER( 5)      &
     &                    .OR. QCHAR .EQ. QLOWER( 5)) GOTO 350
               IF (QCHAR  .EQ. QUPPER( 4)      &
     &                    .OR. QCHAR .EQ. QLOWER( 4)) GOTO 350
               IF (QCHAR  .NE. QPUNCT( 3)) GOTO 930
               IF (REALNO .OR. EXPONT) GOTO 930
                  REALNO = .TRUE.
                  GOTO 360
!
  340          CONTINUE
                  IF (.NOT. EXPONT) THEN
                     IF (SIGNED) GOTO 930
                         SIGNED = .TRUE.
                  ELSE
                     IF (SIGNEX) GOTO 930
                         SIGNEX = .TRUE.
                  ENDIF
                  GOTO 360
!
  350          CONTINUE
                  IF (EXPONT) GOTO 930
                      EXPONT = .TRUE.
                      SIGNED = .FALSE.
  360       CONTINUE
            J = NAMLEN + 1
!
!     Number has been read. Copy info
!
  370    CONTINUE
            DO K=I,J-1
               QNUM(K+1-I) = QCOPY(K)
            END DO
            IF (.NOT. REALNO .AND. .NOT. EXPONT)      &
     &         QNUM(J+1-I) = '.'
            I = J
            READ (QNUM1, 1100) ATEMP(NNUM)
            IF (NNAME .LT. 3) GOTO 190
!
  400 CONTINUE
         RETURN
!
!     END card is missing
!                  
  910 CONTINUE
         IERR = 1
         WRITE (IOLOG, 1700)
         RETURN
!
!     Error after start. This should not happen!
!
  920 CONTINUE
         IERR = 2
         WRITE (IOLOG, 1000) QCOPY
         WRITE (IOLOG, 1800)
         RETURN
!
!     Illegal format in a numeric field
!
  930 CONTINUE
         WRITE (IOLOG, 1900) QCOPY
         GOTO 190
!
 1000 FORMAT(255(A1))
 1100 FORMAT(G255.20)
 1200 FORMAT(A255)
 1300 FORMAT(4A1,2X,2A6,2X,4F10.4,2X,A8)
 1600 FORMAT(' XXX - WARNING - Switching to free-format input.')
 1700 FORMAT(' XXX - WARNING - ENDATA card missing on I/O',      &
     &       ' channel',i3)
 1800 FORMAT(' XXX -  FATAL  - System error in routine RDARC.',      &
     &       ' Try again.')
 1900 FORMAT(' XXX - WARNING - Illegal format in numeric field.',      &
     &       ' Ignore this record:',/,80A1)
      END SUBROUTINE
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE GETNAM(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,      &
     &                                    IERR,IONUM,NREC,LENGTH)
!
!     This routine reads one record of an MPS file and extracts from it
!     a 4-character key and a single name.
!
!     Adapted from MINOS routine M3SPC1 and ob1 routine hspec1.
!
!     Assumptions about the MPS file:
!
!     1. Only code information in the first four character positions.
!     2. Following the codes there is one name field, possibly
!        followed by a comment. In free input format, the name field
!        begins in column 6 or later.
!     3. Cards whose first non-blank character is an asterisk or an
!        exclamation mark are considered comment cards.
!     4. Data items are separated by at least one space.
!     5. Strings can have embedded blanks and even the comment
!        characters, but a blank followed by a comment character marks
!        the start of the comment. ("Problem 1" and "STORM; version 1.1"
!        are valid names, but "STORM ; version 1.1" is not. The latter
!        will return just "STORM" in the name field.
!     6. The line length is limited to 255 characters.
!
!     ==================================================================
!
!     This version dated 7 November 1997
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) DNAME(*),QNAME1,QBLANK,QLINE,CCOPY
      CHARACTER*4 QCODE
      CHARACTER*1 QCOPY(NAMLEN),QPUNCT(4),QNAME(NAMLEN),QB1(NAMLEN)
      REAL*8      ATEMP(*)
      INTEGER*4   LENGTH(*)
      EQUIVALENCE (QNAME,QNAME1),(QBLANK,QB1), (QCOPY, CCOPY)
!
      DATA QPUNCT/'!','*','.','='/
      DATA QB1 /NAMLEN*' '/
!
  100 CONTINUE
      READ (IONUM, 1000, END=910, ERR=920) QCOPY
      NREC = NREC + 1
      IERR = 0
      DNAME(1) = QBLANK
      DNAME(2) = QBLANK
      QLINE = CCOPY
!
      Q1 = QCOPY(1)
      Q2 = QCOPY(2)
      Q3 = QCOPY(3)
      Q4 = QCOPY(4)
!
!     Check for trivial comment line
!
      IF (Q1 .EQ. QPUNCT(1) .OR. Q1 .EQ. QPUNCT(2) ) GOTO 100
         QCODE = Q1 // Q2 // Q3 //Q4
         IF (QCODE .EQ. ' !'     .OR.      &
     &       QCODE .EQ. '  !'    .OR.      &
     &       QCODE .EQ. '   !' ) GOTO 100
!
!     Skip all blanks to the next item -- if there is one
!
            DO I=6,NAMLEN
               IF (QCOPY(I) .NE. QBL) GOTO 210
            END DO
            GOTO 400
!
!     Ignore comment cards
!
  210 CONTINUE
         IF ( QCOPY(I) .EQ. QPUNCT(1) .OR.      &
     &        QCOPY(I) .EQ. QPUNCT(2)) GOTO 400
            QNAME1 = QBLANK
            DO J=I,NAMLEN
               QCHAR = QCOPY(J)
               IF ( (QCHAR .EQ. QPUNCT(1) .OR. QCHAR .EQ. QPUNCT(2) )      &
     &                    .AND. QCOPY(J-1) .EQ. QBL) GOTO 260
                  QNAME(J+1-I) = QCHAR
            END DO
!
  260    CONTINUE
            READ (QNAME1, 1200) DNAME(1)
!
  400 CONTINUE
         RETURN
!
!     END card is missing
!                  
  910 CONTINUE
         IERR = 1
         WRITE (IOLOG, 1100)
         RETURN
!
!     Error after start. This should not happen!
!
  920 CONTINUE
         IERR = 2
         WRITE (IOLOG, 1000) QCOPY
         WRITE (IOLOG, 1100)
         RETURN
!
 1000 FORMAT(255(A1))
 1100 FORMAT(' XXX - WARNING - End of file or file not found.')
 1200 FORMAT(A255)
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE GETPAR(QLINE,Q1,Q2,Q3,Q4,DNAME,ATEMP,XPAR,IPAR,MTYPE,      &
     &                  IBLOCK,IERR,IONUM,NREC)
!
!     This routine reads one line of the stoch file and extracts from
!     the record one numeric field, in standard MPS position. It
!     further determines the type of the number represented (real or
!     integer) and recovers gracefully if the record read is a header card.
!     Alternatively, the name fields are searched to determine the location
!     of a problem coefficient or decision variable.
!
!     Information is returned to the calling routine in variables XPAR,
!     IPAR,MTYPE,and IBLOCK. If the parameter is integer or real, XPAR
!     and IPAR will, respectively, contain the parameter value; IBLOCK
!     is not used. If the parameters references a location, IBLOCK gives
!     the block containing the reference, IPAR gives the location
!     within the block.
!
!     ==================================================================
!     This version dated 3 December 1997.
!     Written by Gus Gassmann.
!     ==================================================================
!
      USE SMPS_DATA
      USE TEMP_DATA
      USE IO_HANDLING
      IMPLICIT REAL*8 (A-H,O,P,R-Z), INTEGER*4 (I-N), CHARACTER*1 (Q)
!
      CHARACTER*(NAMLEN) QNUM1, DNAME(*), QNAME1, QBLANK, DTEMP(2),      &
     &            CCOPY, QLINE
      CHARACTER*4 QCODE
      CHARACTER*1 QCOPY(NAMLEN),QLOWER(26),QUPPER(26),QDIGIT(12),      &
     &            QPUNCT(4),QNUM(NAMLEN),QNAME(NAMLEN),QB1(NAMLEN)
      CHARACTER*8 DBLANK
      LOGICAL*4   REALNO,SIGNED,EXPONT,SIGNEX,PUTDOT(2)
      REAL*8      ATEMP(*)
!
      EQUIVALENCE (QNUM1,QNUM),(QNAME,QNAME1),(QBLANK,QB1),(QCOPY,CCOPY)
!
      DATA QLOWER/'a','b','c','d','e','f','g','h','i','j','k','l','m',      &
     &            'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      DATA QUPPER/'A','B','C','D','E','F','G','H','I','J','K','L','M',      &
     &            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      DATA QDIGIT/'1','2','3','4','5','6','7','8','9','0','+','-'/
      DATA QPUNCT/'!','*','.','='/
      DATA QB1 /NAMLEN*' '/
      DATA DBLANK  /'        '/
!
      XPAR = 0.D0
      IPAR = 0
!
!     Read the next line; ignore if comment card
!
  100 CONTINUE
         READ (IONUM, 1000, END=910, ERR=920) QCOPY
         NREC = NREC + 1
         IF (Q1 .EQ. QPUNCT(1) .OR. Q1 .EQ. QPUNCT(2) ) GOTO 100
            QCODE = Q1 // Q2 // Q3 //Q4
            IF (QCODE .EQ. ' !'     .OR.      &
     &          QCODE .EQ. '  !'    .OR.      &
     &          QCODE .EQ. '   !' ) GOTO 100
!
         IERR = 0
         QLINE = CCOPY
!
!     Always assume free-format input. Process the fields one at a time.
!
  105 CONTINUE
      NNAME = 0
      NNUM  = 0
      ATEMP(1) = 0.D0
      ATEMP(2) = 0.D0
      DNAME(1) = QBLANK
      DNAME(2) = QBLANK
      DNAME(3) = QBLANK
      DTEMP(1) = DBLANK
      DTEMP(2) = DBLANK
!
      Q1 = QCOPY(1)
      Q2 = QCOPY(2)
      Q3 = QCOPY(3)
      Q4 = QCOPY(4)
      IF (Q1 .EQ. QPUNCT(1)) Q1 = QPUNCT(2)
      IF (Q1 .EQ. QPUNCT(2) .OR. Q2 .EQ. QPUNCT(1) .OR.      &
     &    Q3 .EQ. QPUNCT(1) .OR. Q4 .EQ. QPUNCT(1) ) GOTO 100
         I = 4
!
!     Skip all blanks to the next item -- if there is one
!
  190    CONTINUE
            J = I
            DO I=J+1,NAMLEN
               IF (QCOPY(I) .NE. QBL) GOTO 210
            END DO
            GOTO 400
!
!     Ignore comments
!
  210 CONTINUE
         IF ( QCOPY(I) .EQ. QPUNCT(1) .OR.      &
     &        QCOPY(I) .EQ. QPUNCT(2)) GOTO 400
!
!     Check for digit
!
         DO JDIGIT=1,10
            IF (QCOPY(I) .EQ. QDIGIT(JDIGIT)) GOTO 300
         END DO
         IF (QCOPY(I) .EQ. QDIGIT(11) .OR.      &
     &       QCOPY(I) .EQ. QDIGIT(12)) GOTO 300
         IF (QCOPY(I) .NE. QPUNCT(3))  GOTO 240
!
!     First character is a period. Must be a number if the next is a digit
!
         DO JDIGIT=1,10
            IF (QCOPY(I+1) .EQ. QDIGIT(JDIGIT)) GOTO 300
         END DO
!
!     Not a digit -- must be a name
!
  240    CONTINUE
            NNAME = NNAME + 1
            IF (NNAME .EQ. 3) NNUM = 1
            QNAME1 = QBLANK
            DO J=I,NAMLEN
               QCHAR = QCOPY(J)
               IF (QCHAR .EQ. QBL) GOTO 260
                  QNAME(J+1-I) = QCHAR
            END DO
            J = NAMLEN + 1
!
  260    CONTINUE
            I = J
            READ (QNAME1, 1200) DNAME(NNAME)
            IF (NNUM .LT. 2) GOTO 190
               GOTO 400
!
!     Looks like a number...
!             
  300    CONTINUE
            NNUM = NNUM + 1
            IF (NNUM .EQ. 2) NNAME = 3
            QNUM1 = QBLANK
!
            REALNO = .FALSE.
            SIGNED = .FALSE.
            EXPONT = .FALSE.
            SIGNEX = .FALSE.
            DO 360 J=I,NAMLEN
               QCHAR = QCOPY(J)
               IF (QCHAR .EQ. QBL) GOTO 370
               DO JCHAR=1,10
                  IF (QCHAR .EQ. QDIGIT(JCHAR)) GOTO 360
               END DO
               IF (QCHAR  .EQ. QDIGIT(11)      &
     &                    .OR. QCHAR .EQ. QDIGIT(12)) GOTO 340
               IF (QCHAR  .EQ. QUPPER( 5)      &
     &                    .OR. QCHAR .EQ. QLOWER( 5)) GOTO 350
               IF (QCHAR  .EQ. QUPPER( 4)      &
     &                    .OR. QCHAR .EQ. QLOWER( 4)) GOTO 350
               IF (QCHAR  .NE. QPUNCT( 3)) GOTO 930
               IF (REALNO .OR. EXPONT) GOTO 930
                  REALNO = .TRUE.
                  GOTO 360
!
  340          CONTINUE
                  IF (.NOT. EXPONT) THEN
                     IF (SIGNED) GOTO 930
                         SIGNED = .TRUE.
                  ELSE
                     IF (SIGNEX) GOTO 930
                         SIGNEX = .TRUE.
                  ENDIF
                  GOTO 360
!
  350          CONTINUE
                  IF (EXPONT) GOTO 930
                      EXPONT = .TRUE.
                      SIGNED = .FALSE.
  360       CONTINUE
            J = NAMLEN + 1
!
!     Number has been read. Copy info
!
  370    CONTINUE
            DO K=I,J-1
               QNUM(K+NAMLEN-J) = QCOPY(K)
            END DO
            I = J
            DTEMP(NNUM) = QNUM1
            PUTDOT(NNUM) = (.NOT. REALNO .AND. .NOT. EXPONT)
            IF (NNUM .EQ. 1) NNAME = 2
            IF (NNUM .LT. 2) GOTO 190
!
!     The fields have been processed; interpret the information.
!
  400 CONTINUE
         IF (Q1 .NE. QBL) GOTO 940
         IF (Q2 .EQ. QB  .AND. Q3 .EQ. QL) GOTO 940
         IF (DNAME(1) .NE. DBLANK .OR. DNAME(2) .NE. DBLANK) THEN
            DO I=1,NPER
               KREF(I) = I
               KPATH(I) = I
            END DO
            CALL IDENTIFY(NREC,Q1,Q2,Q3,Q4,DNAME(1),DNAME(2),ATEMP(1),      &
     &                    DNAME(3),ATEMP(2),ISTAGE,ISTOCH,IBLOCK,IRPOS,      &
     &                    IROW,ICOL)
            IF (ISTOCH .LT. 0) THEN
               MTYPE = 0
            ELSE
               MTYPE = -ISTOCH
               IPAR  =  IRPOS
            ENDIF
         ELSEIF (NNUM .EQ. 0) THEN
            MTYPE = 0
         ELSEIF (.NOT. PUTDOT(1)) THEN
            READ (DTEMP(1), 1100) XPAR
            MTYPE = 2
         ELSEIF (DNAME(3) .EQ. '''DOUBLE''') THEN
            QNUM1 = DTEMP(1)
            QNUM(NAMLEN) = '.'
            READ (QNUM1, 1100) XPAR
            MTYPE = 2
         ELSE
            READ (DTEMP(1), 1400) IPAR
            MTYPE = 1
         ENDIF
         GOTO 999
!
!     END card is missing
!
  910 CONTINUE
         IERR = 1
         WRITE (IOLOG, 1000) QCOPY
         WRITE (IOLOG, 1700)
         GOTO 999
!
!     Error after start. This should not happen!
!
  920 CONTINUE
         IERR = 2
         WRITE (IOLOG, 1000) QCOPY
         WRITE (IOLOG, 1800)
         GOTO 999
!
!     Illegal format in a numeric field
!
  930 CONTINUE
         WRITE (IOLOG, 1900) QCOPY
         GOTO 190
!
!     Found another header or code card
!
  940 CONTINUE
         IERR = 3
         DO I=1,NNUM
            QNUM1 = DTEMP(I)
            IF (PUTDOT(I)) QNUM(NAMLEN) = '.'
            READ (QNUM1, 1100, ERR=920) ATEMP(I)
         END DO
!
  999 CONTINUE
         RETURN
!
 1000 FORMAT(255(A1))
 1100 FORMAT(G255.20)
 1200 FORMAT(A255)
 1300 FORMAT(4A1,A8,2X,A8,2X,F12.6,3X,A8,2X,F12.6)
 1400 FORMAT(I254)
 1500 FORMAT(' XXX -  FATAL  - Memory handling error in GETPAR.')
 1700 FORMAT(' XXX - WARNING - ENDATA card missing on I/O',      &
     &       ' channel',i3)
 1800 FORMAT(' XXX -  FATAL  - System error in routine RDMPS.',      &
     &       ' Try again.')
 1900 FORMAT(' XXX - WARNING - Illegal format in numeric field.',      &
     &       ' Ignore this record:',/,80A1)
      END SUBROUTINE
!
      END MODULE
!
! ======================================================================
!
      MODULE SAMPUTIL
!
!     This collection of routines aids the work with random variables and
!     stochastic elements. It contains the universal interface to all
!     distributions described in the SMPS standard. This routine can be
!     extended by the user to include calls to user-defined subroutines.
!
!     Also part of this collection are the distributions themselves. This
!     part can be expanded as the need arises. For now the system provides
!     facilities for the discrete, uniform, normal, beta, gamma, and
!     lognormal distribution. For each distribution there is a routine
!     to generate a random sample and a routine to compute conditional
!     probabilities and conditional means over regions of the form
!     {x | a < x < b} (where a or b may be infinity).
!
!     The uniform random number generator RNDGET underlying all this can
!     be selected in routine RDSPEC. (It may be useful to have several
!     generators, just in case one of them gets one into trouble.)
!
! ======================================================================
!
      CONTAINS
!
! ======================================================================
!
      SUBROUTINE CALLDIST(DNAME, INTPAR,  NINTPAR, DBLPAR, NDBLPAR,      &
     &                           XLOCPAR, NLOCPAR, X,      NDIMX,      &
     &                           MODE,    IERR,      &
     &                           PROB,    A,       B,      INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   CALLDIST
!
!       Purpose:   This routine provides the interface to the random number
!                  generators and integration routines.
!
!       Calling sequence:
!
!           CALL CALLDIST(DNAME, INTPAR,  NINTPAR, DBLPAR, NDBLPAR,
!    *                           XLOCPAR, NLOCPAR, X,      NDIMX,
!    *                           MODE,    IERR,
!    *                           PROB,    A,       B,      INF)
!
!      DNAME    The name of the routine to be called (character string)
!      INTPAR   Array containing integer parameters (in the order specified
!               in the stoch file)
!      NINTPAR  Number of integer parameters specified in the stoch file
!      DBLPAR   Array for real-valued parameters (in the order specified
!               in the stoch file)
!      NDBLPAR  Number of real-valued parameters specified in the stoch file
!      XLOCPAR  Array containing the VALUES of parameters referenced
!               by location (in the order specified in the stoch file)
!      NLOCPAR  Number of parameters referenced by location
!
!      X        Array which on return contains the realization or expectation
!      NDIMX    Number of elements in the random vector
!
!      MODE     Information about the value that is to be returned
!               0 - generate a random realization
!               1 - compute a conditional probability
!               2 - compute a conditional probability and conditional mean
!
!      IERR     Returns the error status (0 - no error; > 0 - error condition
!
!      PROB     Returns the conditional probability (if MODE > 0)
!
!      A        Lower endpoints of the integration region (if MODE > 0)
!      B        Upper endpoints of the integration region (if MODE > 0)
!      INF      Shape of the integration region (if MODE > 0):
!               If INF(k) =  2, then kth dimension is   A(k) < X <   B(k)
!               If INF(k) =  1, then kth dimension is   A(k) < X < +infty
!               If INF(k) =  2, then kth dimension is -infty < X <   B(k)
!               If INF(k) = -1, then kth dimension is -infty < X < +infty
!
!      The last four arguments are optional.
!
!
!      Calls:  DISCRETE, UNIFORM, BETA, GAMA, NORMAL, LOGNORM, MVNORMAL,
!              CONSTANT, plus any routines defined by the user.
!
!      Method:
!
!          We use a simple jump table based on the value of the parameter DNAME.
!
!      Remarks:
!
!          Since the standard is self-expanding, this routine must be maintained
!          by the user if he/she defines special subroutines in the stoch file.
!
!          To use this routine, the calling program must contain the line
!
!          USE SAMPUTIL
!
!          immediately following the PROGRAM or SUBROUTINE statement.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      INTEGER     INTPAR(*), MODE,NINTPAR,NDBLPAR,NLOCPAR,NDIMX,IERR
      REAL*8      DBLPAR(*),X(*),XLOCPAR(*)
      CHARACTER*8 DNAME
      REAL (KIND=8), OPTIONAL :: A(*), B(*), PROB
      INTEGER,OPTIONAL :: INF(*)
!
      IERR = 0
!
! -------------------------------------------------------------------
      IF     (DNAME .EQ. 'DISCRETE') THEN
         IF (MODE .EQ. 0) THEN
               CALL DISCRETE(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR)
         ELSE
            IF (PRESENT(PROB) .AND. PRESENT(A) .AND.      &
     &          PRESENT(INF ) .AND. PRESENT(B) ) THEN
               CALL DISCRETE(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR,      &
     &         PROB,A,B,INF)
            ELSE
               IERR = 1
               RETURN
            ENDIF
         ENDIF
!
! -------------------------------------------------------------------
      ELSEIF (DNAME .EQ. 'UNIFORM ') THEN
         IF (MODE .EQ. 0) THEN
               CALL UNIFORM(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR)
         ELSE
            IF (PRESENT(PROB) .AND. PRESENT(A) .AND.      &
     &          PRESENT(INF ) .AND. PRESENT(B) ) THEN
               CALL UNIFORM(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR,      &
     &         PROB,A,B,INF)
            ELSE
               IERR = 1
               RETURN
            ENDIF
         ENDIF
!
! -------------------------------------------------------------------
      ELSEIF (DNAME .EQ. 'NORMAL  ') THEN

         IF (MODE .EQ. 0) THEN
               CALL NORMAL(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR)
         ELSE
            IF (PRESENT(PROB) .AND. PRESENT(A) .AND.      &
     &          PRESENT(INF ) .AND. PRESENT(B) ) THEN
               CALL NORMAL(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR,      &
     &         PROB,A,B,INF)
            ELSE
               IERR = 1
               RETURN
            ENDIF
         ENDIF
!
! -------------------------------------------------------------------
      ELSEIF (DNAME .EQ. 'GAMMA   ') THEN
         IF (MODE .EQ. 0) THEN
               CALL GAMA(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR)
         ELSE
            IF (PRESENT(PROB) .AND. PRESENT(A) .AND.      &
     &          PRESENT(INF ) .AND. PRESENT(B) ) THEN
               CALL GAMA(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR,      &
     &         PROB,A,B,INF)
            ELSE
               IERR = 1
               RETURN
            ENDIF
         ENDIF
!
! -------------------------------------------------------------------
      ELSEIF (DNAME .EQ. 'LOGNORM ') THEN
         IF (MODE .EQ. 0) THEN
               CALL LOGNORM(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR)
         ELSE
            IF (PRESENT(PROB) .AND. PRESENT(A) .AND.      &
     &          PRESENT(INF ) .AND. PRESENT(B) ) THEN
               CALL LOGNORM(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR,      &
     &         PROB,A,B,INF)
            ELSE
               IERR = 1
               RETURN
            ENDIF
         ENDIF
!
! -------------------------------------------------------------------
      ELSEIF (DNAME .EQ. 'BETA    ') THEN
         IF (MODE .EQ. 0) THEN
               CALL BETA(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR)
         ELSE
            IF (PRESENT(PROB) .AND. PRESENT(A) .AND.      &
     &          PRESENT(INF ) .AND. PRESENT(B) ) THEN
               CALL BETA(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR,      &
     &         PROB,A,B,INF)
            ELSE
               IERR = 1
               RETURN
            ENDIF
         ENDIF
!
! -------------------------------------------------------------------
      ELSEIF (DNAME .EQ. 'MVNORMAL') THEN
         IF (MODE .EQ. 0) THEN
               CALL MVNORMAL(DBLPAR(1),NDBLPAR,INTPAR(1),NINTPAR,      &
     &                       X,NDIMX,MODE,IERR)
         ELSE
            IF (PRESENT(PROB) .AND. PRESENT(A) .AND.      &
     &          PRESENT(INF ) .AND. PRESENT(B) ) THEN
               CALL MVNORMAL(DBLPAR(1),NDBLPAR,INTPAR(1),NINTPAR,      &
     &                       X,NDIMX,MODE,IERR,PROB,A,B,INF)
            ELSE
               IERR = 1
               RETURN
            ENDIF
         ENDIF
!
! -------------------------------------------------------------------
      ELSEIF (DNAME .EQ. 'CONSTANT') THEN
         IF (MODE .EQ. 0) THEN
               CALL CONSTANT(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR)
         ELSE
            IF (PRESENT(PROB) .AND. PRESENT(A) .AND.      &
     &          PRESENT(INF ) .AND. PRESENT(B) ) THEN
               CALL CONSTANT(DBLPAR(1),NDBLPAR,X,NDIMX,MODE,IERR,      &
     &         PROB,A,B,INF)
            ELSE
               IERR = 1
               RETURN
            ENDIF
         ENDIF
!
! -------------------------------------------------------------------
      ELSEIF (DNAME .EQ. 'user_def') THEN
         IF (MODE .EQ. 0) THEN
               CALL user_def(DBLPAR(1),NDBLPAR,INTPAR(1),NINTPAR,      &
     &            XLOCPAR(1),NLOCPAR,X,NDIMX,MODE,IERR)
         ELSE
            IF (PRESENT(PROB) .AND. PRESENT(A) .AND.      &
     &          PRESENT(INF ) .AND. PRESENT(B) ) THEN
               CALL user_def(DBLPAR(1),NDBLPAR,INTPAR(1),NINTPAR,      &
     &            XLOCPAR(1),NLOCPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
            ELSE
               IERR = 1
               RETURN
            ENDIF
         ENDIF
!
! -------------------------------------------------------------------
!      ELSEIF (DNAME .EQ. '...') THEN
!         CALL ...
!
! -------------------------------------------------------------------
      ELSE
         WRITE(IOLOG, 1000) DNAME
         IERR = 2
         X(1:NDIMX) = 0.D0
      ENDIF
      RETURN
!
 1000 FORMAT('XXX -  FATAL  - Distribution routine ',A8,' not found.')
 1100 FORMAT('XXX -  FATAL  - Too many parameters called by location.',      &
     & /,' Increase parameter MXLOC (routine CALLDIST) to at least',I7)
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     This routine is added as a simple example of a user-defined
!     distribution.
!
      SUBROUTINE user_def(PAR,NPAR,IPAR,NIPAR,XPAR,LPAR,      &
     &            X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
! --------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8  PAR(*), XPAR(*), X(*)
      INTEGER NPAR, NIPAR, LPAR, NDIMX, MODE, IERR
      REAL (KIND=8), OPTIONAL :: PROB, A(*), B(*)
      INTEGER,OPTIONAL :: IPAR(*), INF(*)
!
      IERR = 0
      IF (NPAR .NE. 2 .OR. NIPAR .NE. 2 .OR. LPAR .NE. 1) THEN
         WRITE (IOLOG, 1000) NPAR
         IERR = 2
         RETURN
      ENDIF
!
      XVAL = PAR(1)
      IF (MODE .EQ. 0) THEN
         DO 100 I=1,NDIMX
            X(I) = XVAL
  100    CONTINUE
!
      ELSEIF (MODE .EQ. 1 .OR. MODE .EQ. 2) THEN
         IF (INF(1) .EQ. -1 .OR. INF(1) .EQ. 1) THEN
            D = XVAL
         ELSE
            D = MIN(XVAL,B(1))
         ENDIF
!
         IF (INF(1) .EQ. -1 .OR. INF(1) .EQ. 0) THEN
            C = XVAL
         ELSE
            C = MAX(XVAL,A(1))
         ENDIF
!
         IF (C .LE. XVAL .AND. D .GE. XVAL) THEN
            PROB = 1.D0
         ELSE
            PROB = 0.D0
         ENDIF
         IF (MODE .EQ. 2) X(1) = XVAL
      ELSE
         WRITE (IOLOG, 1200) MODE
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX -  FATAL  - Incorrect number of parameters for user',     &
     &       ' defined distribution.')
 1200 FORMAT(' XXX - WARNING - Unrecognized value for MODE:',I12)
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE UNIFORM(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   UNIFORM
!
!       Purpose:   Interface for the continuous uniform distribution on [a,b]
!
!       Calling sequence:
!
!          CALL UNIFORM(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!             PAR   --- Array containing the lower and upper endpoints
!                       PAR(1) = a
!                       PAR(2) = b
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns one realization of the random vector or the
!                       conditional mean (if MODE > 0)
!             NDIMX --- Dimension of the random vector X; number of normal
!                       variates to generate
!
!             MODE  --- Information about the value that is to be returned
!               0 - generate a random value/vector
!               1 - compute conditional probability
!               2 - compute conditional probability and mean
!
!             IERR     Error status
!                  = 0 No error
!                  > 0 Error condition
!
!             PROB     Returns the conditional probability (if MODE > 0)
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(k) =  2, then kth dimension is   A(k) < X <   B(k)
!               If INF(k) =  1, then kth dimension is   A(k) < X < +infty
!               If INF(k) =  2, then kth dimension is -infty < X <   B(k)
!               If INF(k) = -1, then kth dimension is -infty < X < +infty
!
!             The last four arguments are optional.
!
!       Calls:  RNDUNIF, CNDUNIF
!
!       Method:
!
!          This routine calls RNDUNIF or CNDUNIF, based on the value of MODE.
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8  PAR(*), X(*)
      INTEGER NPAR, NDIMX, MODE, IERR
      REAL (KIND=8), OPTIONAL :: PROB, A(*), B(*)
      INTEGER,OPTIONAL :: INF(*)
!
      IERR = 0
      IF (NPAR .LT. 2) THEN
         WRITE (IOLOG, 1000) NPAR
         IERR = 2
         RETURN
      ELSEIF (NPAR .GT. 2) THEN
         WRITE (IOLOG, 1100) NPAR
      ENDIF
!
      IF (MODE .EQ. 0) THEN
         CALL RNDUNIF(PAR,X,NDIMX,IERR)
      ELSE
         CALL CNDUNIF(PAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX -  FATAL  - Missing real parameters for uniform',      &
     &       ' distribution. Need 2, have',I2,'.')
 1100 FORMAT(' XXX - WARNING - Too many real parameters for uniform',      &
     &       ' distribution.',/,'    Need 2, have',I3,'.',      &
     &       ' Extra parameters ignored.')
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE RNDUNIF(PAR,X,NDIMX,IERR)
!
!--------------------------------------------------------------------------
!
!       Routine:   RNDUNIF
!
!       Purpose:   Generator for the continuous uniform distribution on [a,b]
!
!       Calling sequence:
!
!          CALL RNDUNIF(PAR,X,NDIMX,IERR)
!
!             PAR   --- Array containing the lower and upper endpoints
!                       PAR(1) = a
!                       PAR(2) = b
!             X     --- Returns one realization of the random vector
!             NDIMX --- Dimension of the random vector X
!             IERR  --- Error status
!
!       Calls:  RNDGET
!
!       Method:
!
!          Each component of the random vector X is generated as a pseudo-random
!          number Y uniformly distributed on [0,1] and translated as follows:
!
!              X(i) = a + Y * ( b - a )
!
!
!       Remarks:
!
!          This routine generates a univariate random variable. If NDIMX > 1,
!          we generate NDIMX independent uniform variates.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8 PAR(*), X(*), XHI, XLO, DIF
!
      IERR = 0
!
      XLO = PAR(1)
      XHI = PAR(2)
      DIF = XHI - XLO
!
      CALL RNDGET(X,NDIMX,IERR)
      DO I=1,NDIMX
         X(I) = XLO + X(I) * DIF
      END DO
!
      RETURN
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CNDUNIF(PAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   CNDUNIF
!
!       Purpose:   Computes conditional probabilities and conditional means
!                  for the continuous uniform distribution on [a,b]
!
!       Calling sequence:
!
!          CALL CNDUNIF(PAR,NPAR,X,NDIMX MODE,PROB,A,B,INF)
!
!             PAR   --- Array containing the lower and upper endpoints
!                       PAR(1) = a
!                       PAR(2) = b
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns the conditional expectation of the random vector
!             NDIMX --- Dimension of the random vector X
!             IERR  --- Error status
!             MODE  --- Information about the value that is to be returned
!               1 - compute conditional probability
!               2 - compute conditional probability and mean
!
!             PROB     Returns the conditional probability
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(1) =  2, then region is   A(1) < X <   B(1)
!               If INF(1) =  1, then region is   A(1) < X < +infty
!               If INF(1) =  2, then region is -infty < X <   B(1)
!               If INF(1) = -1, then region is -infty < X < +infty
!
!       Calls:  None.
!
!       Method:
!
!          Prob [c < X] = (c-a)/(b-a) if c in [a,b]
!          E[X | c < X] = (a+c)/2     if c in [a,b]
!          etc.
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      INTEGER MODE, INF(*)
      REAL*8 PAR(*), X(*), XHI, XLO, C, D, PROB, A(*), B(*)
!
      XLO = PAR(1)
      XHI = PAR(2)
!
      IF (INF(1) .EQ. -1 .OR. INF(1) .EQ. 1) THEN
         D = XHI
      ELSE
         D = MIN(XHI,B(1))
      ENDIF
!
      IF (INF(1) .EQ. -1 .OR. INF(1) .EQ. 0) THEN
         C = XLO
      ELSE
         C = MAX(XLO,A(1))
      ENDIF
!
      IF (C .GE. D) THEN
         PROB = 0.D0
         IF (MODE .EQ. 2) X(1) = XLO
      ELSEIF (D .LT. XLO) THEN
         PROB = 0.D0
         IF (MODE .EQ. 2) X(1) = XLO
      ELSEIF (C .GT. XHI) THEN
         PROB = 0.D0
         IF (MODE .EQ. 2) X(1) = XHI
      ELSE
         IF (XHI .GT. XLO) THEN
            PROB = (D - C) / (XHI - XLO)
            IF (MODE .EQ. 2) X(1) = (C + D) * 0.5D0
         ELSEIF (XHI .LT. XLO) THEN
            WRITE (IOLOG, 1000)
            IERR = 1
         ELSE
            IF (C .LE. XLO .AND. D .GE. XHI) THEN
               PROB = 1.D0
               IF (MODE .EQ. 2) X(1) = XLO
            ELSE
               PROB = 0.D0
               IF (MODE .EQ. 2) X(1) = XLO
            ENDIF
         ENDIF
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX - WARNING - Parameter error in routine CNDUNIF')
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE NORMAL(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   NORMAL
!
!       Purpose:   Interface for the normal distribution with mean MU and
!                  variance SIGMA^2.
!
!       Calling sequence:
!
!          CALL NORMAL(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!             PAR   --- Array containing the mean and variance
!                       PAR(1) = MU
!                       PAR(2) = SIGMA^2
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns one realization of the random vector or the
!                       conditional mean (if MODE > 0)
!             NDIMX --- Dimension of the random vector X; number of normal
!                       variates to generate
!
!             MODE  --- Information about the value that is to be returned
!               0 - generate a random value/vector
!               1 - compute Prob[X < c    ]
!               2 - compute Prob[    X > d]
!               3 - compute Prob[c < X < d]
!               4 - compute Prob[X < c    ] and E[X| X < c    ]
!               5 - compute Prob[    X > d] and E[X|     X > d]
!               6 - compute Prob[c < X < d] and E[X| c < X < d]
!
!             IERR     Error status: 0 - no error; > 0 - error condition
!
!             PROB     Returns the conditional probability (if MODE > 0)
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(k) =  2, then kth dimension is   A(k) < X <   B(k)
!               If INF(k) =  1, then kth dimension is   A(k) < X < +infty
!               If INF(k) =  2, then kth dimension is -infty < X <   B(k)
!               If INF(k) = -1, then kth dimension is -infty < X < +infty
!
!             The last four arguments are optional.
!
!       Calls:  RNDNORM, CNDNORM
!
!       Method:
!
!          This routine calls RNDNORM or CNDNORM, based on the value of MODE.
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8  PAR(*), X(*)
      INTEGER NPAR, NDIMX, MODE, IERR
      REAL (KIND=8), OPTIONAL :: PROB, A(*), B(*)
      INTEGER,OPTIONAL :: INF(*)
!
      IERR = 0
      IF (NPAR .LT. 2) THEN
         WRITE (IOLOG, 1000) NPAR
         IERR = 2
         RETURN
      ELSEIF (NPAR .GT. 2) THEN
         WRITE (IOLOG, 1100) NPAR
         IERR = 1
      ENDIF
!
      IF (MODE .EQ. 0) THEN
         CALL RNDNORM(PAR,NPAR,X,NDIMX,IERR)
      ELSE
         CALL CNDNORM(PAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX -  FATAL  - Missing real parameters for normal',      &
     &       ' distribution. Need 2, have',I2,'.')
 1100 FORMAT(' XXX - WARNING - Too many real parameters for normal',      &
     &       ' distribution.',/,'    Need 2, have',I3,'.',      &
     &       ' Extra parameters ignored.')
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE RNDNORM(PAR,NPAR,X,NDIMX,IERR)
!
!--------------------------------------------------------------------------
!
!       Routine:   RNDNORM
!
!       Purpose:   Generator for the normal distribution with mean MU and
!                  variance SIGMA^2.
!
!       Calling sequence:
!
!          CALL RNDNORM(PAR,NPAR,X,NDIMX)
!
!             PAR   --- Array containing the mean and variance
!                       PAR(1) = MU
!                       PAR(2) = SIGMA^2
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns one realization of the random vector
!             NDIMX --- Dimension of the random vector X; number of normal
!                       variates to generate
!
!       Calls:  RANDNORM
!
!       Method:
!
!          Each component of the random vector X is generated as a standard
!          normally distributed pseudo-random number Y, translated as follows:
!
!              X(i) = MU + SIGMA * NINV2(Y)
!
!
!       Remarks:
!
!          This routine generates a univariate random variable. If NDIMX > 1,
!          we generate NDIMX independent normal variates.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8 MU,SIGMA,X(*),PAR(*)
!
      MU    = PAR(1)
      IF (PAR(2) .LT. 0.D0) THEN
         WRITE (IOLOG, 1000)
         IERR = 1
         RETURN
      ENDIF
!
      IERR = 0
      SIGMA = DSQRT(PAR(2))
!
      DO I=1,NDIMX
         X(I) = MU + SIGMA*RANDNORM()
      END DO
!
      RETURN
!
 1000 FORMAT(' XXX - WARNING - Illegal parameter value: negative', &
     &       ' variance specified.')
!
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      REAL*8 FUNCTION RANDNORM()
!
!************************************************************************
!     Normal-Distribution Random Number Generator.
!
!     Computes a sample value from a normal distribution characterized by
!     zero mean and unit standard deviation.  Uses "Box-Muller polar method"
!     for calculation of normal deviates,
!
!     Box, G. E. P., and Muller, M. E., "A Note on the Generation of Normal
!     Deviates," ANNALS OF MATH. STAT., Vol. 29, pp. 610-611, 1958.
!
!     Refined in
!
!     Knuth, D. E., FUNDAMENTAL ALGORITHMS, Vol II, "Seminumerical Algorithms,"
!     Addison-Wesley Pub. Co., page 117.
!+---------------------------------------------------------------------------
!
      REAL*8 D,X(2)
      LOGICAL*4 EMPTY
      DATA EMPTY/.TRUE./
!
!     SAVE X, EMPTY
      SAVE X
!
      IF ( EMPTY ) THEN
  100 CONTINUE
         CALL RNDGET(X,2,IERR)
         X(1) = 2.D0*X(1) - 1.D0
         X(2) = 2.D0*X(2) - 1.D0
         D  = X(1)**2 + X(2)**2
         IF (D .GT. 1.D0 .OR. D .LE. 0) GOTO 100
            D  = DSQRT(-2.D0 * DLOG(D) / D)
            X(1) = X(1) * D
            X(2) = X(2) * D
            EMPTY = .FALSE.
            RANDNORM = X(1)
      ELSE
            EMPTY = .TRUE.
            RANDNORM = X(2)
      ENDIF
      RETURN
      END function
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CNDNORM(PAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   CNDNORM
!
!       Purpose:   Computes conditional probabilities and conditional means
!                  for the normal distribution with mean MU and variance SIGMA^2.
!
!       Calling sequence:
!
!          CALL CNDNORM(PAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
!
!             PAR   --- Array containing the mean and variance
!                       PAR(1) = MU
!                       PAR(2) = SIGMA^2
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns the conditional expectation of the random vector
!             NDIMX --- Dimension of the random vector X
!             IERR  --- Error status
!             MODE  --- Information about the value that is to be returned
!               1 - compute conditional probability
!               2 - compute conditional probability and mean
!
!             PROB     Returns the conditional probability
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(1) =  2, then region is   A(1) < X <   B(1)
!               If INF(1) =  1, then region is   A(1) < X < +infty
!               If INF(1) =  2, then region is -infty < X <   B(1)
!               If INF(1) = -1, then region is -infty < X < +infty
!
!       Calls:  SIGNORM
!
!       Method:
!
!          Conditional probabilities are computed by calling routine SIGNORM,
!          Conditional means are computed by analytic integration:
!             Integral(x*exp[(-x^2)/2] dx) = -exp[(-x^2)/2] + C
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      INTEGER MODE, INF(*)
      REAL*8 MU,SIGMA,X(*),PAR(*),C,D,PROB,A(*),B(*),      &
     &       ONEBY2PI
      DATA ONEBY2PI/0.3989422804014D0/
!
      MU    = PAR(1)
      IF (PAR(2) .LT. 0.D0) THEN
         WRITE (IOLOG, 1000)
         IERR = 1
         RETURN
      ENDIF
!
      IERR = 0
      SIGMA = DSQRT(PAR(2))
!
      IF     (INF(1) .EQ. 0) THEN
          C = (B(1) - MU) / SIGMA
          IF (C .LT. 0.D0) THEN
             PROB = SIGNORM( C )
          ELSE
             PROB = 1.D0 - SIGNORM( C )
          ENDIF
          IF (MODE .EQ. 2)      &
     &        X(1) = MU - SIGMA * ONEBY2PI * EXP(-0.5*(C**2)) / PROB
!
      ELSEIF (INF(1) .EQ. 1) THEN
          D = (A(1) - MU) / SIGMA
          IF (D .GT. 0.D0) THEN
             PROB = SIGNORM( D )
          ELSE
             PROB = 1.D0 - SIGNORM( D )
          ENDIF
          IF (MODE .EQ. 2)      &
     &        X(1) = MU + SIGMA * ONEBY2PI * EXP(-0.5*(D**2)) / PROB
!
      ELSEIF (INF(1) .EQ. 2) THEN
          C = (A(1) - MU) / SIGMA
          D = (B(1) - MU) / SIGMA
          IF (C .GT. D) THEN
              PROB = 0.D0
              X(1) = 0.D0
          ELSE
             IF (D .LT. 0.D0) THEN
                PROB = SIGNORM( D ) - SIGNORM( C )
             ELSEIF (C .GT. 0.D0) THEN
                PROB = SIGNORM( C ) - SIGNORM( D )
             ELSE
                PROB = 1.D0 - SIGNORM( C ) - SIGNORM( D )
             ENDIF
             IF (MODE .EQ. 2) X(1) = MU + SIGMA * ONEBY2PI *      &
     &                 (EXP(-0.5*(C**2)) - EXP(-0.5*(D**2))) / PROB
          ENDIF
!
      ELSEIF (INF(1) .EQ. -1) THEN
         PROB = 1.D0
         IF (MODE .EQ. 2) X(1) = MU
!
      ELSE
         WRITE (IOLOG, 1100) MODE
         IERR = 1
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX - WARNING - Illegal parameter value: negative', &
     &       ' variance specified.')
 1100 FORMAT(' XXX - WARNING - Unrecognized value for MODE:',I12)
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE MVNORMAL(RPAR,NRPAR,IPAR,NIPAR,      &
     &                    X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   MVNORMAL
!
!       Purpose:   Interface for the multivariate normal distribution
!
!       Calling sequence:
!
!          CALL MVNORMAL(RPAR,NRPAR,IPAR,NIPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!             RPAR  --- Array containing the mean and covariance matrix or
!                       its Cholesky factorization
!
!       Data structure for covariance matrix:
!       Mu(1),Var(1),Mu(2),Var(2),Cov(1,2),Mu(3),Var(3),Cov(1,3),Cov(2,3),...
!
!       The Cholesky factor L is stored in the same positions
!       (i.e. L(i,j) is stored in the same position as Cov(i,j))
!
!             NRPAR --- Number of parameters (should be NDIMX*(NDIMX+1))
!             IPAR  --- Array containing the integer parameter
!                       IPAR(1) = 0: RPAR contains the covariance matrix
!                       IPAR(1) = 1: RPAR contains the Cholesky factorization
!             X     --- Returns one realization of the random vector or the
!                       conditional mean (if MODE > 0)
!             NDIMX --- Dimension of the random vector X; number of normal
!                       variates to generate
!
!             MODE     Information about the value that is to be returned
!                      0 - generate a random realization
!                      1 - compute conditional probability
!                      2 - compute conditional probability and mean
!
!             IERR     Returns THE error status (0 - no error; > 0 - error condition
!
!             PROB     Returns the conditional probability (if MODE > 0)
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(k) =  2, then kth dimension is   A(k) < X <   B(k)
!               If INF(k) =  1, then kth dimension is   A(k) < X < +infty
!               If INF(k) =  2, then kth dimension is -infty < X <   B(k)
!               If INF(k) = -1, then kth dimension is -infty < X < +infty
!
!             The last four arguments are optional.
!
!       Calls:  RNDMNML, CNDMNML
!
!       Method:
!
!          This routine calls RNDMNML or CNDMNML, based on the value of MODE.
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8  RPAR(*), X(*)
      INTEGER NRPAR, NPAR2, NDIMX, MODE, IERR, IPAR(*)
      REAL (KIND=8), OPTIONAL :: PROB, A(*), B(*)
      INTEGER,OPTIONAL :: INF(*)
!
!     Verify the number of parameters
!
      IERR = 0
      NPAR2 = NDIMX*(NDIMX+3)/2
      IF (NRPAR .LT. NPAR2) THEN
         WRITE (IOLOG, 1000) NRPAR,NPAR2
         IERR = 2
         RETURN
      ELSEIF (NRPAR .GT. NPAR2) THEN
         WRITE (IOLOG, 1100) NRPAR,NPAR2
         IERR = 1
      ENDIF
!
      IF (NIPAR .LT. 1) THEN
         WRITE (IOLOG, 1200)
         IERR = 2
         RETURN
      ELSEIF (NIPAR .GT. 1) THEN
         WRITE (IOLOG, 1300) NIPAR
         IERR = 1
      ENDIF
!
      IF (MODE .EQ. 0) THEN
         CALL RNDMNML(RPAR,NRPAR,IPAR,NIPAR,X,NDIMX,IERR)
      ELSE
         CALL CNDMNML(RPAR,NRPAR,IPAR,NIPAR,X,NDIMX,IERR,MODE,      &
     &                PROB,A,B,INF)
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX -  FATAL  - Missing real parameters for',      &
     &       ' multinormal distribution. Have',I4,' need',I4,'.')
 1100 FORMAT(' XXX - WARNING - Too many real parameters for',      &
     &       ' multinormal distribution.',/,'   Have',I4,' need',I4,      &
     &      '. Extra parameters ignored.')
 1200 FORMAT(' XXX -  FATAL  - Missing integer parameter for',      &
     &       ' multinormal distribution.')
 1300 FORMAT(' XXX - WARNING - Too many integer parameters for',      &
     &       ' multinormal distribution.',/,'   Have',I4,' need  1.',      &
     &       ' Extra parameters ignored.')
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE RNDMNML(RPAR,NRPAR,IPAR,NIPAR,X,NDIMX,IERR)
!
!--------------------------------------------------------------------------
!
!       Routine:   RNDMNML
!
!       Purpose:   Generator for the multivariate normal distribution
!
!       Calling sequence:
!
!          CALL RNDMNML(RPAR,NRPAR,IPAR,NIPAR,X,NDIMX)
!
!             RPAR  --- Array containing the mean and covariance matrix or
!                       its Cholesky factorization
!
!       Data structure for covariance matrix:
!       Mu(1),Var(1),Mu(2),Var(2),Cov(1,2),Mu(3),Var(3),Cov(1,3),Cov(2,3),...
!
!       The Cholesky factor L is stored in the same positions
!       (i.e. L(i,j) is stored in the same position as Cov(i,j))
!
!             NRPAR --- Number of parameters (should be NDIMX*(NDIMX+1))
!             IPAR  --- Array containing the integer parameter
!                       IPAR(1) = 0: RPAR contains the covariance matrix
!                       IPAR(1) = 1: RPAR contains the Cholesky factorization
!             X     --- Returns one realization of the random vector
!             NDIMX --- Dimension of the random vector X
!             IERR  --- Error status
!
!       Calls:  RNDGET, NINV2, CHOLES
!
!       Method:
!
!          If necessary, the covariance matrix is factored first, using the
!          Cholesky factorization. At that time we also check whether the matrix
!          is positive definite (or positive semi-definite).
!          We then generate a vector Y of (0,1) normal variates and compute X as
!
!             X = L * Y + MU
!
!       Remarks:
!
!          This works because of the well-known linear transformation property
!          of the normal random variable: If U is a normal random variable with
!          mean M and covariance C, then S*U + T is normal with mean S * M + T
!          and covariance matrix S*C*S(trans).
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8    X(*),RPAR(*),TEMP
      INTEGER*4 IPAR(*)
!
!     Find the Cholesky factorization if necessary
!
      IF (IPAR(1) .EQ. 0) THEN
         CALL CHOLES(RPAR,NDIMX,IERR)
         IF (IERR .GT. 0) RETURN
         IPAR(1) = 1
      ENDIF
!
!     Get NDIMX independent standard normal variates
!
      CALL RNDGET(X,NDIMX,IERR)
      DO I=1,NDIMX
         TEMP = NINV2(X(I))
         X(I) = TEMP
      END DO
!
!     Compute the affine transformation X = L*Y + MU
!
      CALL CHMULT(RPAR,X,NDIMX)
      RETURN
!
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CHOLES(RPAR,NDIM,IERR)
!
!--------------------------------------------------------------------------
!
!       Routine:   CHOLES
!
!       Purpose:   Performs a Cholesky factorization of the covariance matrix
!
!       Calling sequence:
!
!          CALL CHOLES(RPAR,NDIM,IERR)
!
!             RPAR --- Array containing the mean and covariance matrix on call
!                       and the Cholesky factorization on return
!
!       Data structure for covariance matrix:
!       Mu(1),Var(1),Mu(2),Var(2),Cov(1,2),Mu(3),Var(3),Cov(1,3),Cov(2,3),...
!
!       The Cholesky factor L is stored in the same positions
!       (i.e. L(i,j) is stored in the same position as Cov(i,j))
!
!             NDIM --- Dimension of the covariance matrix
!             IERR --- Error status
!
!       Calls:  None
!
!       Method:
!
!          The Cholesky factor is found in a triple loop:
!
!          For I = 1 to NDIM
!              Xii = SQRT(Xii)
!              For J = I+1 to NDIM
!                  Xij = Xij / Xii
!                  For K = I+1 to J
!                      Xjk = Xjk - Xij * Xik
!
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8    RPAR(*),ZTOLZE
      INTEGER*4 ISTEP,JSTEP,KSTEP
!
      DATA ZTOLZE/1.D-6/
!
      IERR  = 0
      II    = 2
      ISTEP = 2
      DO I=1,NDIM
         IF (RPAR(II) .GT. ZTOLZE) THEN
            RPAR(II) = DSQRT(RPAR(II))
            JSTEP = ISTEP + 1
            IJ    = II    + JSTEP
            JJ    = II    + JSTEP - 1
            DO J=I+1,NDIM
               RPAR(IJ) = RPAR(IJ) / RPAR(II)
               RPAR(JJ) = RPAR(JJ) - RPAR(IJ)**2
               IK = II + ISTEP  + 1
               JK = JJ + I + 1
               KSTEP  = ISTEP  + 1
               DO K=I+1,J-1
                  RPAR(JK) = RPAR(JK) - RPAR(IJ) * RPAR(IK)
                  JK = JK + 1
                  IK = IK + KSTEP
                  KSTEP  = KSTEP + 1
               END DO
               IJ = IJ + JSTEP
               JJ = JJ + JSTEP
               JSTEP = JSTEP + 1
            END DO
         ELSEIF (RPAR(II) .GT. -ZTOLZE) THEN
            WRITE (IOLOG, 1000)
            RPAR(II) = 0.D0
            JSTEP = ISTEP + 1
            DO J=I+1,NDIM
               RPAR(IJ) = 0.D0
               IJ = IJ + JSTEP
               JSTEP = JSTEP + 1
            END DO
         ELSE
            WRITE (IOLOG, 1100)
            IERR = 2
            RETURN
         ENDIF
         II = II + ISTEP
         ISTEP = ISTEP + 1
      END DO
!
      RETURN
!
 1000 FORMAT(' Warning: Cholesky factor is singular.')
 1100 FORMAT(' Fatal:   Cholesky factor cannot be computed.')
!
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CHMULT(RPAR,X,NDIM)
!
!--------------------------------------------------------------------------
!
!       Routine:   CHMULT
!
!       Purpose:   Multiply the Cholesky factorization of the covariance matrix
!                  by a randomly generated vecor of standard normal variates.
!
!       Calling sequence:
!
!          CALL CHMULT(RPAR,X,NDIM)
!
!             RPAR --- Array containing the means and the Cholesky factorization
!
!       Data structure for covariance matrix:
!       Mu(1),L(1,1),Mu(2),L(2,2),L(1,2),Mu(3),L(3,3),L(1,3),L(2,3),...
!
!             X    --- On entry contains the vector of standard normal variates,
!                      on return contains the transformed values L*X + MU
!
!             NDIM --- Dimension of the Cholesky factor and the random vector X
!
!       Calls:  None
!
!       Method:
!
!          The product is found in a double loop:
!
!          For I = NDIM to 1
!              X(i) = L(i,i)*X(i) + MU(i)
!              For J = 1 to I-1
!                  X(i) = X(i) + L(i,j)*X(j)
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      REAL*8    X(*),RPAR(*)
      INTEGER*4 ISTEP,II
!
      II    = NDIM*(NDIM+1)/2 + 1
      ISTEP = NDIM
      DO I=NDIM,1,-1
         X(I) = RPAR(II)*X(I) + RPAR(II-1)
         DO J=1,I-1
            X(I) = X(I) + RPAR(II+J)*X(J)
         END DO
         II = II - ISTEP
         ISTEP = ISTEP - 1
      END DO
!
      RETURN
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CNDMNML(RPAR, NRPAR, IPAR, NIPAR, X, NDIMX, IERR,      &
     &                   MODE, PROB,  A, B, INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   CNDMNML
!
!       Purpose:   Computes conditional probabilities and conditional means
!                  for the multivariate normal distribution
!
!       Calling sequence:
!
!          CALL CNDMNML(RPAR,NRPAR,IPAR,NIPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
!
!             RPAR  --- Array containing the mean and covariance matrix or
!                       its Cholesky factorization
!
!       Data structure for covariance matrix:
!       Mu(1),Var(1),Mu(2),Var(2),Cov(1,2),Mu(3),Var(3),Cov(1,3),Cov(2,3),...
!
!       The Cholesky factor L is stored in the same positions
!       (i.e. L(i,j) is stored in the same position as Cov(i,j))
!
!             NRPAR --- Number of parameters (should be NDIMX*(NDIMX+1))
!             IPAR  --- Array containing the integer parameter
!                       IPAR(1) = 0: RPAR contains the covariance matrix
!                       IPAR(1) = 1: RPAR contains the Cholesky factorization
!             X     --- Returns the conditional expectation of the random vector
!             NDIMX --- Dimension of the random vector X
!             IERR  --- Error status
!             MODE  --- Information about the value that is to be returned
!               1 - compute conditional probability
!               2 - compute conditional probability and mean
!
!             IERR     Error status
!                  = 0 No error
!                  > 0 Error condition
!
!             PROB     Returns the conditional probability (if MODE > 0)
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(k) =  2, then kth dimension is   A(k) < X <   B(k)
!               If INF(k) =  1, then kth dimension is   A(k) < X < +infty
!               If INF(k) =  2, then kth dimension is -infty < X <   B(k)
!               If INF(k) = -1, then kth dimension is -infty < X < +infty
!
!       Calls:  RNDGET, CHOLES
!
!       Method:
!
!          If necessary, the Covariance matrix is factored first, using the
!          Cholesky factorization. At that time we also check whether the matrix
!          is positive definite (or positive semi-definite).
!          The computation of the conditional probability and the conditional mean
!          depends on the dimension NDIMX of the random vector: if NDIMX = 2, the
!          probability is computed by series expansion, and the mean is found by
!          analytic integration, if NDIMX > 2, we use a method described in
!            H.I. Gassmann, "Conditional probability and conditional expectation
!            of a random vector", in: Yu. Ermoliev and R.J-B Wets (eds.),
!            Numerical Techniques for Stochastic Optimization, Springer Verlag,
!            1988.
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8    X(*),RPAR(*), PROB, A(*), B(*)
      INTEGER*4 IPAR(*),MODE, INF(*)
!
!     Find the Cholesky factorization if necessary
!
      IF (IPAR(1) .EQ. 0) THEN
         CALL CHOLES(RPAR,NDIMX,IERR)
         IF (IERR .GT. 0) RETURN
         IPAR(1) = 1
      ENDIF
!
      IF (MODE .GT. 0 .AND. MODE .LE. 6) THEN
         IF (NDIMX .EQ. 2) THEN
            WRITE (IOLOG, 1000)
         ELSE
            WRITE (IOLOG, 1100)
         ENDIF
      ELSE
         WRITE (IOLOG, 1000) MODE
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX - WARNING - Unrecognized value for MODE:',I12)
 1100 FORMAT(' XXX - WARNING - Method not implemented yet.')
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE BETA(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   BETA
!
!       Purpose:   Interface for the beta distribution with shape
!                  parameters b and c.
!
!       Calling sequence:
!
!          CALL BETA(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!             PAR   --- Array containing the scale and shape
!                       PAR(1) = scale
!                       PAR(2) = shape
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns one realization of the random vector or the
!                       conditional mean (if MODE > 0)
!
!             NDIMX --- Dimension of the random vector X (should be 1)
!             MODE  --- Information about the value that is to be returned
!               0 - generate a random value/vector
!               1 - compute conditional probability
!               2 - compute conditional probability and mean
!
!             IERR     Error status
!                  = 0 No error
!                  > 0 Error condition
!
!             PROB     Returns the conditional probability (if MODE > 0)
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(k) =  2, then kth dimension is   A(k) < X <   B(k)
!               If INF(k) =  1, then kth dimension is   A(k) < X < +infty
!               If INF(k) =  2, then kth dimension is -infty < X <   B(k)
!               If INF(k) = -1, then kth dimension is -infty < X < +infty
!
!             The last four arguments are optional.
!
!       Calls:  RNDBETA, CNDBETA
!
!       Method:
!
!          This routine calls RNDBETA or CNDBETA, based on the value of MODE.
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8  PAR(*), X(*)
      INTEGER NPAR, NDIMX, MODE, IERR
      REAL (KIND=8), OPTIONAL :: PROB, A(*), B(*)
      INTEGER,OPTIONAL :: INF(*)
!
      IERR = 0
      IF (NPAR .LT. 2) THEN
         WRITE (IOLOG, 1000) NPAR
         IERR = 2
         RETURN
      ELSEIF (NPAR .GT. 2) THEN
         WRITE (IOLOG, 1100) NPAR
         IERR = 1
      ENDIF
!
      IF (MODE .EQ. 0) THEN
         CALL RNDBETA(PAR,NPAR,X,NDIMX,IERR)
      ELSE
         CALL CNDBETA(PAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX -  FATAL  - Missing real parameters for beta',      &
     &       ' distribution. Need 2, have',I2,'.')
 1100 FORMAT(' XXX - WARNING - Too many real parameters for beta',      &
     &       ' distribution.',/,'    Need 2, have',I3,'.',      &
     &       ' Extra parameters ignored.')
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE RNDBETA(PAR,NPAR,X,NDIMX,IERR)
!
!--------------------------------------------------------------------------
!
!       Routine:   RNDBETA
!
!       Purpose:   Generator for the beta distribution with shape
!                  parameters b and c.
!
!       Calling sequence:
!
!          CALL RNDBETA(PAR,NPAR,X,NDIMX)
!
!             PAR   --- Array containing the scale and shape
!                       PAR(1) = scale
!                       PAR(2) = shape
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns one realization of the random vector
!             NDIMX --- Dimension of the random vector X (should be 1)
!
!       Calls: RNDGAMA.
!
!       Method:
!
!          The algorithm uses an approach described in I. Deak, Random
!          Number Generators and Simulation, Akademiai Kiado, Budapest, 1990.
!
!          If X and Y are two independent standard Gamma variates with
!          shape parameter b and c, respectively, then X/(X+Y) is beta
!          distributed with shape parameters b and c.
!
!       Remarks:
!
!          This routine generates a univariate random variable. If NDIMX > 1,
!          we generate NDIMX independent beta variates.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8 X(*),PAR(*),WORK1(2),WORK2(2),WORK3(10)
      INTEGER*4 N2,K0
!
      WORK1(1) = 1.D0
      WORK1(2) = PAR(1)
      WORK2(1) = 1.D0
      WORK2(2) = PAR(2)
!
      CALL RNDGAMA(WORK1,2,X,NDIMX,IERR)
      DO I=1,NDIMX,10
         K0 = 10*(I-1)
         IF (K0+10 .GT. NDIMX) THEN
            N2 = NDIMX - K0
         ELSE
            N2 = 10
         ENDIF
         CALL RNDGAMA(WORK2,2,WORK3,N2,IERR)
         DO K=1,N2
            X(K0+K) = X(K0+K)/(X(K0+K)+WORK3(K))
         END DO
      END DO
!
      RETURN
!
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CNDBETA(PAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   CNDBETA
!
!       Purpose:   Computes conditional probabilities and conditional means
!                  for the beta distribution with shape parameters b and c.
!
!       Calling sequence:
!
!          CALL CNDBETA(PAR,NPAR,X,NDIMX,MODE,PROB,A,B,INF)
!
!             PAR   --- Array containing the scale and shape
!                       PAR(1) = scale
!                       PAR(2) = shape
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns the conditional expectation of the random vector
!             NDIMX --- Dimension of the random vector X (should be 1)
!             MODE  --- Information about the value that is to be returned
!               1 - compute conditional probability
!               2 - compute conditional probability and mean
!
!             PROB     Returns the conditional probability
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(1) =  2, then region is   A(1) < X <   B(1)
!               If INF(1) =  1, then region is   A(1) < X < +infty
!               If INF(1) =  2, then region is -infty < X <   B(1)
!               If INF(1) = -1, then region is -infty < X < +infty
!
!       Calls:   CDBETA.
!
!       Method:  Uses the fact that int( t * t**(a-1) * (1-t)**(b-1)) dt
!                is a scalar multiple of the incomplete beta function with
!                parameters (a+1) and b.
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8 X(*),PAR(*),T, C, D,PROB,A(*),B(*)
      INTEGER*4 MODE,IDIG,ITER,IFAULT, INF(*)
!
      MDIG  = 12
      MITER = 100
!
      IF     (INF(1) .EQ. 0) THEN
         IF     (B(1) .LE. 0.D0) THEN
            PROB = 0.D0
            IF (MODE .EQ. 2) X(1) = 0.D0
         ELSEIF (B(1) .GE. 1.D0) THEN
            PROB = 1.D0
            IF (MODE .EQ. 2) X(1) = PAR(1)/(PAR(1)+PAR(2))
         ELSE
            PROB = CDBETA(B(1),PAR(1),PAR(2),MDIG,MITER,        &
     &                    IDIG,ITER,IFAULT)
            IF (IFAULT .GT. 0) THEN
               WRITE (IOLOG, 1100) MDIG,MITER
               IERR = 1
            ENDIF
            IF (MODE .EQ. 2) THEN
               T = CDBETA(B(1),PAR(1)+1,PAR(2),MDIG,MITER,      &
     &                    IDIG,ITER,IFAULT)
               IF (IFAULT .GT. 0) THEN
                  WRITE (IOLOG, 1100) MDIG,MITER
                  IERR = 1
               ENDIF
               X(1) = PAR(1)/(PAR(1)+PAR(2)) * T / PROB
            ENDIF
         ENDIF
      ELSEIF (INF(1) .EQ. 1) THEN
         IF     (A(1) .GE. 1.D0) THEN
            PROB = 0.D0
            IF (MODE .EQ. 2) X(1) = 1.D0
         ELSEIF (A(1) .LE. 0.D0) THEN
            PROB = 1.D0
            IF (MODE .EQ. 2) X(1) = PAR(1)/(PAR(1)+PAR(2))
         ELSE
            PROB = 1.D0 - CDBETA(A(1),PAR(1),PAR(2),MDIG,MITER,      &
     &                    IDIG,ITER,IFAULT)
            IF (IFAULT .GT. 0) THEN
               WRITE (IOLOG, 1100) MDIG,MITER
               IERR = 1
            ENDIF
            IF (MODE .EQ. 2) THEN
               T = 1.D0 - CDBETA(A(1),PAR(1)+1,PAR(2),MDIG,MITER,      &
     &                    IDIG,ITER,IFAULT)
               IF (IFAULT .GT. 0) THEN
                  WRITE (IOLOG, 1100) MDIG,MITER
                  IERR = 1
               ENDIF
               X(1) = PAR(1)/(PAR(1)+PAR(2)) * T / PROB
            ENDIF
         ENDIF
      ELSEIF (INF(1) .EQ. 2) THEN
         IF (A(1) .GE. B(1)) THEN
            PROB = 0.D0
            IF (MODE .EQ. 2) X(1) = 0.D0
         ELSEIF (A(1) .GE. 1.D0) THEN
            PROB = 0.D0
            IF (MODE .EQ. 2) X(1) = 1.D0
         ELSEIF (B(1) .LE. 0.D0) THEN
            PROB = 0.D0
            IF (MODE .EQ. 2) X(1) = 0.D0
         ELSE
            C = A(1)
            D = B(1)
            IF (C .LT. 0.D0) C = 0.D0
            IF (D .GT. 1.D0) D = 1.D0
            PROB = CDBETA(D,PAR(1),PAR(2),MDIG,MITER,IDIG,ITER,IFAULT)      &
     &           - CDBETA(C,PAR(1),PAR(2),MDIG,MITER,IDIG,ITER,IFAULT)
            IF (IFAULT .GT. 0) THEN
               WRITE (IOLOG, 1100) MDIG,MITER
               IERR = 1
            ENDIF
            IF (MODE .EQ. 2) THEN
               T = CDBETA(D,PAR(1)+1,PAR(2),MDIG,MITER,IDIG,ITER,IFAULT)      &
     &           - CDBETA(C,PAR(1)+1,PAR(2),MDIG,MITER,IDIG,ITER,IFAULT)
               IF (IFAULT .GT. 0) THEN
                  WRITE (IOLOG, 1100) MDIG,MITER
                  IERR = 1
               ENDIF
               X(1) = PAR(1)/(PAR(1)+PAR(2)) * T / PROB
            ENDIF
         ENDIF
!
      ELSEIF (INF(1) .EQ. -1) THEN
         PROB = 1.D0
         IF (MODE .EQ. 2) X(1) = PAR(1)/(PAR(1)+PAR(2))
      ELSE
         WRITE (IOLOG, 1000) MODE
      ENDIF
      RETURN
!
 1000 FORMAT(' XXX - WARNING - Unrecognized value for MODE:',I12)
 1100 FORMAT(' XXX - WARNING - Error in routine DNDBETA.',/,' Increase',      &
     &       ' maximum iteration count (currently',i6,') or decrease',      &
     &       ' accuracy (currently',i4,')')
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE GAMA(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   GAMMA
!
!       Purpose:   Interface for the Gamma distribution with scale
!                  parameter b and shape parameter c.
!
!       Calling sequence:
!
!          CALL GAMA(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!             PAR   --- Array containing the scale and shape
!                       PAR(1) = scale
!                       PAR(2) = shape
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns one realization of the random vector or the
!                       conditional mean (if MODE > 0)
!             NDIMX --- Dimension of the random vector X
!             MODE  --- Information about the value that is to be returned
!               0 - generate a random value/vector
!               1 - compute conditional probability
!               2 - compute conditional probability and mean
!
!             IERR     Error status
!                  = 0 No error
!                  > 0 Error condition
!
!             PROB     Returns the conditional probability (if MODE > 0)
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(k) =  2, then kth dimension is   A(k) < X <   B(k)
!               If INF(k) =  1, then kth dimension is   A(k) < X < +infty
!               If INF(k) =  2, then kth dimension is -infty < X <   B(k)
!               If INF(k) = -1, then kth dimension is -infty < X < +infty
!
!             The last four arguments are optional.
!
!       Calls:  RNDGAMA, CNDGAMA
!
!       Method:
!
!          This routine calls RNDGAMA or CNDGAMA, based on the value of MODE.
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8  PAR(*), X(*)
      INTEGER NPAR, NDIMX, MODE, IERR
      REAL (KIND=8), OPTIONAL :: PROB, A(*), B(*)
      INTEGER,OPTIONAL :: INF(*)
!
      IERR = 0
      IF (NPAR .LT. 2) THEN
         WRITE (IOLOG, 1000) NPAR
         IERR = 2
         RETURN
      ELSEIF (NPAR .GT. 2) THEN
         WRITE (IOLOG, 1100) NPAR
         IERR = 1
      ENDIF
!
      IF (MODE .EQ. 0) THEN
         CALL RNDGAMA(PAR,NPAR,X,NDIMX,IERR)
      ELSE
         CALL CNDGAMA(PAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX -  FATAL  - Missing real parameters for Gamma',      &
     &       ' distribution. Need 2, have',I2,'.')
 1100 FORMAT(' XXX - WARNING - Too many real parameters for Gamma',      &
     &       ' distribution.',/,'    Need 2, have',I3,'.',      &
     &       ' Extra parameters ignored.')
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE RNDGAMA(PAR,NPAR,X,NDIMX,IERR)
!
!--------------------------------------------------------------------------
!
!       Routine:   RNDGAMA
!
!       Purpose:   Generator for the Gamma distribution with scale
!                  parameter b and shape parameter c.
!
!       Calling sequence:
!
!          CALL RNDGAMA(PAR,NPAR,X,NDIMX,IERR)
!
!             PAR   --- Array containing the scale and shape
!                       PAR(1) = scale
!                       PAR(2) = shape
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns one realization of the random vector
!             NDIMX --- Dimension of the random vector X (should be 1)
!
!       Calls:  RNDGET
!
!       Method:
!
!          The algorithm uses an approach described in I. Deak, Random
!          Number Generators and Simulation, Akademiai Kiado, Budapest, 1990.
!
!          If the shape parameter is greater than 1, we use acceptance-rejection
!          based on the Cauchy density; if the shape parameter is less than 1,
!          acceptance-rejection is based on a different function; if c = 1, the
!          Gamma is an exponential distribution which can be generated directly
!          by inverting the density.
!
!          A standard Gamma variate Y (scale = 1) is generated, then we return
!          X = b * Y as Gamma(b,c).
!
!       Remarks:
!
!          This routine generates a univariate random variable. If NDIMX > 1,
!          we generate NDIMX independent Gamma variates.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8 SCALE,SHAPE,X(*),PAR(*),E,WORK(2),T,PI,P1,P2,P3
      DATA PI/ 3.1415926535897932385D0 /
!
      E = DEXP(1.D0)
!
      SCALE = PAR(1)
      SHAPE = PAR(2)
!
!     Shape parameter = 1; this is just an exponential distribution
!
      IF (SHAPE .EQ. 1.D0) THEN
         CALL RNDGET(X,NDIMX,IERR)
         DO I=1,NDIMX
            X(I) = -DLOG(X(I))
         END DO
!
!     Shape parameter < 1; this is algorithm G2 from Deak's book
!
      ELSEIF (SHAPE .LT. 1.D0) THEN
         DO I=1,NDIMX
  110       CONTINUE
            CALL RNDGET(WORK,2,IERR)
            P = WORK(1)*(E+SHAPE)/E
            IF (P .GT. 1.D0) THEN
               T = P**(1/SHAPE)
               IF (WORK(2) .GT. DEXP(-T)) GOTO 110
            ELSE
               T = -DLOG(((SHAPE+E)/E-P)/SHAPE)
               IF (WORK(2) .GT. T**(SHAPE-1)) GOTO 110
            ENDIF
            X(I) = T
         END DO
!
!     Shape parameter > 1; this is algorithm G4 from Deak's book
!
      ELSE
         P1 = SHAPE - 1.D0
         P2 = 2.D0*SHAPE - 1.D0
         P3 = DSQRT(P2)
         DO I=1,NDIMX
  130       CONTINUE
            CALL RNDGET(WORK,2,IERR)
            T = P3*DTAN(PI*(WORK(1)-0.5D0))
            Y = T + P1
            IF (Y .LT. 0.D0) GOTO 130
               S = (1+(T**2)/P2)*DEXP(P1*DLOG(Y/P1)-T)
               IF (S .LT. WORK(2)) GOTO 130
                  X(I) = Y
         END DO
      ENDIF
!
      DO I=1,NDIMX
         X(I) = SCALE*X(I)
      END DO
!
      RETURN
!
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CNDGAMA(PAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   CNDGAMA
!
!       Purpose:   Computes conditional probabilities and conditional means
!                  for the continuous uniform distribution on [a,b]
!
!       Calling sequence:
!
!          CALL CNDGAMA(PAR,NPAR,X,NDIMX,MODE,PROB,A,B,INF)
!
!             PAR   --- Array containing the scale and shape
!                       PAR(1) = scale
!                       PAR(2) = shape
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns the coditional expectation of the random vector
!             NDIMX --- Dimension of the random vector X (should be 1)
!             MODE  --- Information about the value that is to be returned
!               1 - compute conditional probability
!               2 - compute conditional probability and mean
!
!             PROB     Returns the conditional probability
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(1) =  2, then region is   A(1) < X <   B(1)
!               If INF(1) =  1, then region is   A(1) < X < +infty
!               If INF(1) =  2, then region is -infty < X <   B(1)
!               If INF(1) = -1, then region is -infty < X < +infty
!
!       Calls:   GAMAIN.
!
!       Method:  Uses the fact that int( t * t**(c-1) * e**(-c)) dt
!                is a scalar multiple of the incomplete gamma function
!                with shape parameter (c+1).
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8  C,D,X(*),PAR(*),PROB,A(*),B(*)
      INTEGER INF(*)
!
      MDIG  = 12
      MITER = 100
!
      IF     (INF(1) .EQ. 0) THEN
         IF     (B(1) .LE. 0.D0) THEN
            PROB = 0.D0
            IF (MODE .EQ. 2) X(1) = 0.D0
         ELSE
            C = B(1)/PAR(2)
            PROB = GAMAIN(C,PAR(1),MDIG,MITER,IDIG,ITER,IER)
            IF (IER  .GT. 0) THEN
               WRITE (IOLOG, 1100) MDIG,MITER
               IERR = 1
            ENDIF
            IF (MODE .EQ. 2) THEN
               T = GAMAIN(C,PAR(1)+1,MDIG,MITER,IDIG,ITER,IER)
               IF (IER .GT. 0) WRITE (IOLOG, 1100) MDIG,MITER
               X(1) = PAR(1) * PAR(2) * T / PROB
            ENDIF
         ENDIF
      ELSEIF (INF(1) .EQ. 1) THEN
         IF     (A(1) .LE. 0.D0) THEN
            PROB = 1.D0
            IF (MODE .EQ. 2) X(1) = PAR(1) * PAR(2)
         ELSE
            C = A(1)/PAR(2)
            PROB = 1.D0 - GAMAIN(C,PAR(1),MDIG,MITER,IDIG,ITER,IER)
            IF (IER  .GT. 0) THEN
               WRITE (IOLOG, 1100) MDIG,MITER
               IERR = 1
            ENDIF
            IF (MODE .EQ. 2) THEN
               T = 1.0 - GAMAIN(C,PAR(1)+1,MDIG,MITER,IDIG,ITER,IER)
               IF (IER  .GT. 0) THEN
                  WRITE (IOLOG, 1100) MDIG,MITER
                  IERR = 1
               ENDIF
               X(1) = PAR(1) * PAR(2) * T / PROB
            ENDIF
         ENDIF
      ELSEIF (INF(1) .EQ. 2) THEN
         IF (A(1) .GE. B(1)) THEN
            PROB = 0.D0
            IF (MODE .EQ. 2) X(1) = 0.D0
         ELSEIF (B(1) .LE. 0.D0) THEN
            PROB = 0.D0
            IF (MODE .EQ. 2) X(1) = 0.D0
         ELSE
            C = A(1) / PAR(2)
            D = B(1) / PAR(2)
            IF (C .LT. 0.D0) C = 0.D0
            PROB = GAMAIN(D,PAR(1),MDIG,MITER,IDIG,ITER,IER)      &
     &           - GAMAIN(C,PAR(1),MDIG,MITER,IDIG,ITER,IER)
            IF (IER  .GT. 0) THEN
               WRITE (IOLOG, 1100) MDIG,MITER
               IERR = 1
            ENDIF
            IF (MODE .EQ. 2) THEN
               T = GAMAIN(D,PAR(1)+1,MDIG,MITER,IDIG,ITER,IER)      &
     &           - GAMAIN(C,PAR(1)+1,MDIG,MITER,IDIG,ITER,IER)
               IF (IER  .GT. 0) THEN
                  WRITE (IOLOG, 1100) MDIG,MITER
                  IERR = 1
               ENDIF
               X(1) = PAR(1) * PAR(2) * T / PROB
            ENDIF
         ENDIF
!
      ELSEIF (INF(1) .EQ. -1) THEN
         PROB = 1.D0
         IF (MODE .EQ. 2) X(1) = PAR(1) * PAR(2)
      ELSE
         WRITE (IOLOG, 1000) MODE
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX - WARNING - Unrecognized value for MODE:',I12)
 1100 FORMAT(' XXX - WARNING - Requested accuracy cannot be achieved',      &
     &       '; Digits: ',I5,', Iterations:',I6)
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE LOGNORM(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   LOGNORM
!
!       Purpose:   Interface for the lognormal distribution with parameters MU
!                  and SIGMA^2.
!
!       Calling sequence:
!
!          CALL LOGNORM(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!             PAR   --- Array containing the mean and variance
!                       PAR(1) = MU
!                       PAR(2) = SIGMA^2
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns one realization of the random vector or the
!                       conditional mean (if MODE > 0)
!             NDIMX --- Dimension of the random vector X
!             MODE  --- Information about the value that is to be returned
!               0 - generate a random value/vector
!               1 - compute Prob[X < c    ]
!               2 - compute Prob[    X > d]
!               3 - compute Prob[c < X < d]
!               4 - compute Prob[X < c    ] and E[X| X < c    ]
!               5 - compute Prob[    X > d] and E[X|     X > d]
!               6 - compute Prob[c < X < d] and E[X| c < X < d]
!
!             IERR     Error status: 0 - no error; > 0 - error condition
!
!             PROB     Returns the conditional probability (if MODE > 0)
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(k) =  2, then kth dimension is   A(k) < X <   B(k)
!               If INF(k) =  1, then kth dimension is   A(k) < X < +infty
!               If INF(k) =  2, then kth dimension is -infty < X <   B(k)
!               If INF(k) = -1, then kth dimension is -infty < X < +infty
!
!             The last four arguments are optional.
!
!       Calls:  RNDLOGN, CNDLOGN
!
!       Method:
!
!          This routine calls RNDLOGN or CNDLOGN, based on the value of MODE.
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8  PAR(*), X(*)
      INTEGER NPAR, NDIMX, MODE, IERR
      REAL (KIND=8), OPTIONAL :: PROB, A(*), B(*)
      INTEGER,OPTIONAL :: INF(*)
!
      IERR = 0
      IF (NPAR .LT. 2) THEN
         WRITE (IOLOG, 1000) NPAR
         IERR = 2
         RETURN
      ELSEIF (NPAR .GT. 2) THEN
         WRITE (IOLOG, 1100) NPAR
         IERR = 1
      ENDIF
!
      IF (MODE .EQ. 0) THEN
         CALL RNDLOGN(PAR,NPAR,X,NDIMX,IERR)
      ELSE
         CALL CNDLOGN(PAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX -  FATAL  - Missing real parameters for lognormal',      &
     &       ' distribution. Need 2, have',I2,'.')
 1100 FORMAT(' XXX - WARNING - Too many real parameters for lognormal',      &
     &       ' distribution.',/,'    Need 2, have',I3,'.',      &
     &       ' Extra parameters ignored.')
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE RNDLOGN(PAR,NPAR,X,NDIMX,IERR)
!
!--------------------------------------------------------------------------
!
!       Routine:   RNDLOGN
!
!       Purpose:   Generator for the lognormal distribution with parameters MU
!                  and SIGMA^2. (The natural logarithm of this distribution
!                  is normally distributed with mean MU and variance SIGMA^2.)
!
!       Calling sequence:
!
!          CALL RNDLOGN(PAR,X,NDIMX,IERR)
!
!             PAR   --- Array containing the mean and variance
!                       PAR(1) = MU
!                       PAR(2) = SIGMA^2
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns one realization of the random vector
!             NDIMX --- Dimension of the random vector X; number of normal
!                       variates to generate
!
!       Calls:  RANDNORM
!
!       Method:
!
!          Each component of the random vector X is generated as a standard
!          normally distributed pseudo-random number Y, translated as follows:
!
!              X(i) = EXP ( MU + SIGMA * NINV2(Y) )
!
!
!       Remarks:
!
!          This routine generates a univariate random variable. If NDIMX > 1,
!          we generate NDIMX independent normal variates.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8 MU,SIGMA,X(*),PAR(*)
!
      MU    = PAR(1)
!
      IF (PAR(2) .GE. 0.D0) THEN
         SIGMA = DSQRT(PAR(2))
!
         DO I=1,NDIMX
            X(I) = DEXP(MU + SIGMA*RANDNORM())
         END DO
      ELSE
         WRITE (IOLOG, 1000) PAR(2)
         IERR = 1
      ENDIF
!
 1000 FORMAT(' Value for PAR(2) must be nonnegative instead of',F20.10)
      RETURN
!
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CNDLOGN(PAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   CNDLOGN
!
!       Purpose:   Computes conditional probabilities and conditional means
!                  for the lognormal distribution with parameters MU and SIGMA^2.
!
!       Calling sequence:
!
!          CALL CNDLOGN(PAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
!
!             PAR   --- Array containing the parameter values
!                       PAR(1) = MU
!                       PAR(2) = SIGMA^2
!             NPAR  --- Number of parameters (should be 2)
!             X     --- Returns the conditional expectation of the random vector
!             NDIMX --- Dimension of the random vector X (should be 1)
!             MODE  --- Information about the value that is to be returned
!               1 - compute conditional probability
!               2 - compute conditional probability and mean
!
!             PROB     Returns the conditional probability
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(1) =  2, then region is   A(1) < X <   B(1)
!               If INF(1) =  1, then region is   A(1) < X < +infty
!               If INF(1) =  2, then region is -infty < X <   B(1)
!               If INF(1) = -1, then region is -infty < X < +infty
!
!       Calls:  SIGNORM
!
!       Method:
!
!          Conditional probabilities are computed by calling routine SIGNORM,
!          Conditional means are derived by analytic integration:
!             Integral(exp[(-LN x^2)/2] dx) = Integral(exp[u-(u^2)/2] du) =
!                                = exp(0.5) * Integral(exp[-((u-0.5)^2)/2] du)
!
!       Remarks: None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      INTEGER MODE, INF(*)
      REAL*8  MU,SIGMA,X(*),PAR(*),C,D,PROB,A(*),B(*)
!
      MU    = PAR(1)
      IF (PAR(2) .LT. 0.D0) THEN
         WRITE (IOLOG, 1000) PAR(2)
         IERR = 1
         RETURN
      ENDIF
!
      SIGMA = DSQRT(PAR(2))
!
      IF     (INF(1) .EQ. 0) THEN
          IF (B(1) .LE. 0.D0) THEN
             PROB = 0.D0
             IF (MODE .EQ. 2) X(1) = 0.D0
          ELSE
             C = (DLOG(B(1)) - MU) / SIGMA
             IF (C .LT. 0.D0) THEN
                PROB = SIGNORM( C )
             ELSE
                PROB = 1.D0 - SIGNORM( C )
             ENDIF
             IF (MODE .EQ. 2) THEN
                C = C - SIGMA
                IF (C .LT. 0.D0) THEN
                   T = SIGNORM( C )
                ELSE
                   T = 1.D0 - SIGNORM( C )
                ENDIF
                X(1) = EXP(MU + 0.5*SIGMA**2) * T / PROB
             ENDIF
          ENDIF
!
      ELSEIF (INF(1) .EQ. 1) THEN
         IF (A(1) .LE. 0.D0) THEN
             PROB = 1.D0
             IF (MODE .EQ. 2) X(1) = DEXP(MU+0.5*SIGMA**2)
          ELSE
             C = (DLOG(A(1)) - MU) / SIGMA
             IF (C .GT. 0.D0) THEN
                PROB = SIGNORM( C )
             ELSE
                PROB = 1.D0 - SIGNORM( C )
             ENDIF
             IF (MODE .EQ. 2) THEN
                C = C - SIGMA
                IF (C .GT. 0.D0) THEN
                   T = SIGNORM( C )
                ELSE
                   T = 1.D0 - SIGNORM( C )
                ENDIF
                X(1) = EXP(MU + 0.5*SIGMA**2) * T / PROB
             ENDIF
          ENDIF
!
      ELSEIF (INF(1) .EQ. 2) THEN
          C = (A(1) - MU) / SIGMA
          D = (B(1) - MU) / SIGMA
          IF (A(1) .GT. B(1)) THEN
             PROB = 0.D0
             IF (MODE .EQ. 2) X(1) = 0.D0
          ELSEIF (B(1) .LE. 0.D0) THEN
             PROB = 0.D0
             IF (MODE .EQ. 2) X(1) = 0.D0
          ELSEIF (A(1) .LT. 0.D0) THEN
             C = (DLOG(B(1)) - MU) / SIGMA
             IF (C .GT. 0.D0) THEN
                PROB = SIGNORM( C )
             ELSE
                PROB = 1.D0 - SIGNORM( C )
             ENDIF
             IF (MODE .EQ. 2) THEN
                C = C - SIGMA
                IF (C .GT. 0.D0) THEN
                   T = SIGNORM( C )
                ELSE
                   T = 1.D0 - SIGNORM( C )
                ENDIF
                X(1) = EXP(MU + 0.5*SIGMA**2) * T / PROB
             ENDIF
          ELSE
             C = (DLOG(A(1)) - MU) / SIGMA
             D = (DLOG(B(1)) - MU) / SIGMA
             IF (D .LT. 0.D0) THEN
                PROB = SIGNORM( D ) - SIGNORM( C)
             ELSEIF (C .GT. 0.D0) THEN
                PROB = SIGNORM( C ) - SIGNORM( D )
             ELSE
                PROB = 1.D0 - SIGNORM( C ) - SIGNORM( D )
             ENDIF
             IF (MODE .EQ. 2) THEN
                C = C - SIGMA
                D = D - SIGMA
                IF (D .LT. 0.D0) THEN
                   T = SIGNORM( D ) - SIGNORM( C)
                ELSEIF (C .GT. 0.D0) THEN
                   T = SIGNORM( C ) - SIGNORM( D )
                ELSE
                   T = 1.D0 - SIGNORM( C ) - SIGNORM( D )
                ENDIF
                X(1) = EXP(MU + 0.5*SIGMA**2) * T / PROB
             ENDIF
          ENDIF
!
      ELSEIF (INF(1) .EQ. -1) THEN
         PROB = 1.D0
         IF (MODE .EQ. 2) X(1) = DEXP(MU+0.5*SIGMA**2)
      ELSE
         WRITE (IOLOG, 1100) INF(1)
      ENDIF
!
      RETURN
!
 1000 FORMAT(' Value for PAR(2) must be nonnegative instead of',F20.10)
 1100 FORMAT(' XXX - WARNING - Unrecognized value for range:',I12)
!
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CONSTANT(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   CONSTANT
!
!       Purpose:   Auxiliary routine making available a degenerate distribution,
!                  which may be needed, for example, in a linear transformation.
!
!       Calling sequence:
!
!          CALL CONSTANT(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!             PAR   --- Array containing the value of the constant
!                       PAR(1) = XVAL
!             NPAR  --- Number of parameters (should be 1)
!             X     --- Returns one realization of the random vector
!             NDIMX --- Dimension of the random vector X
!             MODE  --- Information about the value that is to be returned
!               0 - generate a random value/vector
!               1 - compute conditional probability
!               2 - compute conditional probability and mean
!
!             IERR     Error status
!                  = 0 No error
!                  > 0 Error condition
!
!             PROB     Returns the conditional probability (if MODE > 0)
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(k) =  2, then kth dimension is   A(k) < X <   B(k)
!               If INF(k) =  1, then kth dimension is   A(k) < X < +infty
!               If INF(k) =  2, then kth dimension is -infty < X <   B(k)
!               If INF(k) = -1, then kth dimension is -infty < X < +infty
!
!             The last four arguments are optional.
!
!       Calls:     None.
!
!       Method:
!
!          The constant XVAL is stored in every component of the vector X.
!
!       Remarks:   None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8  PAR(*), X(*)
      INTEGER NPAR, NDIMX, MODE, IERR
      REAL (KIND=8), OPTIONAL :: PROB, A(*), B(*)
      INTEGER,OPTIONAL :: INF(*)
!
      IERR = 0
      IF (NPAR .LT. 1) THEN
         WRITE (IOLOG, 1000) NPAR
         IERR = 2
         RETURN
      ELSEIF (NPAR .GT. 1) THEN
         WRITE (IOLOG, 1100) NPAR
         IERR = 1
      ENDIF
!
      XVAL = PAR(1)
      IF (MODE .EQ. 0) THEN
         X(1:NDIMX) = XVAL
!
      ELSEIF (MODE .EQ. 1 .OR. MODE .EQ. 2) THEN
         IF (INF(1) .EQ. -1 .OR. INF(1) .EQ. 1) THEN
            D = XVAL
         ELSE
            D = MIN(XVAL,B(1))
         ENDIF
!
         IF (INF(1) .EQ. -1 .OR. INF(1) .EQ. 0) THEN
            C = XVAL
         ELSE
            C = MAX(XVAL,A(1))
         ENDIF
!
         IF (C .LE. XVAL .AND. D .GE. XVAL) THEN
            PROB = 1.D0
         ELSE
            PROB = 0.D0
         ENDIF
         IF (MODE .EQ. 2) X(1) = XVAL
      ELSE
         WRITE (IOLOG, 1200) MODE
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX -  FATAL  - Missing real parameters for degenerate',      &
     &       ' distribution. Need 1, have',I2,'.')
 1100 FORMAT(' XXX - WARNING - Too many real parameters for degenerate',      &
     &       ' distribution.',/,'    Need 1, have',I3,'.',      &
     &       ' Extra parameters ignored.')
 1200 FORMAT(' XXX - WARNING - Unrecognized value for MODE:',I12)
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE DISCRETE(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   DISCRETE
!
!       Purpose:   Interface for a discrete random variable
!
!       Calling sequence:
!
!          CALL DISCRETE(PAR,NPAR,X,NDIMX,MODE,IERR,PROB,A,B,INF)
!
!             PAR   --- Array containing the probabilities and values
!                       PAR(1) = first probability
!                       PAR(2)..PAR(NDIMX+1) = values for the first realization
!                       PAR(NDIMX+2) = second probability
!                       ...
!             NPAR  --- Number of parameters ( = NREAL * (NDIMX+1) )
!             X     --- Returns one realization of the random vector or the
!                       conditional mean (if MODE > 0)
!             NDIMX --- Dimension of the random vector X
!
!             MODE  --- Information about the value that is to be returned
!               0 - generate a random value/vector
!               1 - compute conditional probability
!               2 - compute conditional probability and mean
!
!             IERR     Error status
!                  = 0 No error
!                  > 0 Error condition
!
!             PROB     Returns the conditional probability (if MODE > 0)
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(k) =  2, then kth dimension is   A(k) < X <   B(k)
!               If INF(k) =  1, then kth dimension is   A(k) < X < +infty
!               If INF(k) =  2, then kth dimension is -infty < X <   B(k)
!               If INF(k) = -1, then kth dimension is -infty < X < +infty
!
!             The last four arguments are optional.
!
!       Calls:     RNDDISC, CNDDISC.
!
!       Method:
!
!          This routine calls RNDDISC or CNDDISC, based on the value of MODE.
!
!       Remarks:   None.
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8  PAR(*), X(*)
      INTEGER NPAR, NDIMX, MODE, IERR
      REAL (KIND=8), OPTIONAL :: PROB, A(*), B(*)
      INTEGER,OPTIONAL :: INF(*)
!
      IERR = 0
      NREAL = NPAR/(NDIMX+1)
      IF (NPAR .NE. NREAL*(NDIMX+1)) THEN
         WRITE (IOLOG, 1000) NDIMX+1
         IERR = 2
         RETURN
      ENDIF
!
      IF (MODE .EQ. 0) THEN
         CALL RNDDISC(PAR,NPAR,X,NDIMX,IERR)
      ELSE
         CALL CNDDISC(PAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX -  FATAL  - Missing or extra real parameters for',      &
     &       ' discrete distribution.',/,      &
     &       ' Should have ',i4,' for each realization.')
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE RNDDISC(RPAR,NPAR,X,NDIMX,IERR)
!
!--------------------------------------------------------------------------
!
!       Routine:   RNDDISC
!
!       Purpose:   This routine returns the value of a discrete random variable
!
!       Calling sequence:
!
!          CALL RNDDISC(RPAR,NPAR,X,NDIMX)
!
!             RPAR  --- Array containing the probabilities and values
!                       RPAR(1) = first probability
!                       RPAR(2)..RPAR(NDIMX+1) = values for first realization
!                       RPAR(NDIMX+2) = second probability
!                       ...
!             NPAR  --- Number of parameters (=NREAL * NDIMX)
!             X     --- Returns one realization of the random vector
!             NDIMX --- Dimension of the random vector X
!             IERR  --- Error status
!
!       Calls:     RNDGET.
!
!       Method:
!
!          A uniform random variable U is generated and compared to the
!          cumulative probability for each realization. If
!
!                   qprob(i) < U < qprob(i+1),
!
!          return the values of realization i.
!
!
!       Remarks:   None.
!
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      REAL*8 X(*),TEMP,RPAR(*),QPROB,ZTOLZE,RNDT(1)
      DATA ZTOLZE/1.D-6/
!
      IERR  = 0
      QPROB = 0.D0
      CALL RNDGET(RNDT,1,IERR)
      TEMP = RNDT(1)
!
      DO I0=1,NPAR,NDIMX+1
         QPROB = QPROB + RPAR(I0)
         IF (QPROB .GE. TEMP) THEN
            DO I=1,NDIMX
               X(I) = RPAR(I0+I)
            END DO
            RETURN
         ENDIF
      END DO
!
      IF (DABS(QPROB-1.D0) .GT. ZTOLZE) THEN
         IF (QPROB .GT. 0.D0) THEN
            WRITE (IOLOG, 1100) QPROB
            DO I0=1,NPAR,NDIMX+1
               RPAR(I0) = RPAR(I0) / QPROB
            END DO
         ELSE
            WRITE (IOLOG, 1200) QPROB
            IERR = 2
         ENDIF
      ENDIF
!
      RETURN
!
 1100 FORMAT(' XXX - WARNING - Discrete probabilities sum to ',F12.8,      &
     &       ' instead of 1. Normalizing.')
 1200 FORMAT(' XXX -  FATAL  - Cumulative probability equals ',F12.8,      &
     &       ' instead of 1. This is illegal.')
!
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE CNDDISC(RPAR,NPAR,X,NDIMX,IERR,MODE,PROB,A,B,INF)
!
!--------------------------------------------------------------------------
!
!       Routine:   CNDDISC
!
!       Purpose:   This routine computes conditional probabilities and
!                  conditional means for a discrete random variable
!
!       Calling sequence:
!
!          CALL CNDDISC(PAR,NPAR,X,NDIMX,MODE,PROB,A,B,INF)
!
!             RPAR  --- Array containing the probabilities and values
!                       RPAR(1) = first probability
!                       RPAR(2)..PAR(NDIMX+1) = values for the first realization
!                       RPAR(NDIMX+2) = second probability
!                       ...
!             NPAR  --- Number of parameters (=NREAL * NDIMX)
!             X     --- Returns the conditional mean (if MODE = 2)
!             NDIMX --- Dimension of the random vector X
!             MODE  --- Information about the value that is to be returned
!               1 - compute conditional probability
!               2 - compute conditional probability and mean
!
!             PROB     Returns the conditional probability
!
!             A        Lower endpoints of the region (if MODE > 0)
!             B        Upper endpoints of the region (if MODE > 0)
!             INF      Shape of the integration region (if MODE > 0):
!               If INF(k) =  2, then region is   A(k) < X <   B(k)
!               If INF(k) =  1, then region is   A(k) < X < +infty
!               If INF(k) =  2, then region is -infty < X <   B(k)
!               If INF(k) = -1, then region is -infty < X < +infty
!
!       Calls:     None.
!
!       Method:
!
!          Prob [c < X] = sum[ p(i) | c < X]
!          E[X | c < X] = sum[ X(i) | c < X] / Prob [c < X]
!          etc.
!
!       Remarks:   None.
!
!
!--------------------------------------------------------------------------
!
      USE IO_HANDLING
!
      INTEGER MODE, INF(*)
      REAL*8 X(*),RPAR(*),PROB,A(*),B(*)
!
      PROB  = 0.D0
!
      IF (MODE .EQ. 2) THEN
         X(1:NDIMX) = 0.D0
      ENDIF
!
      DO 160 I0=1,NPAR,NDIMX+1
         DO I=1,NDIMX
            IF     (INF(I) .EQ. 0) THEN
               IF (RPAR(I0+I) .GT. B(I)) GOTO 160
            ELSEIF (INF(I) .EQ. 1) THEN
               IF (RPAR(I0+I) .LT. A(I)) GOTO 160
            ELSEIF (INF(I) .EQ. 2) THEN
               IF (RPAR(I0+I) .LT. A(I) .OR.      &
     &             RPAR(I0+I) .GT. B(I)) GOTO 160
            ELSEIF (INF(I) .NE. -1) THEN
               WRITE (IOLOG, 1000) I,INF(I)
               IERR = 1
               RETURN
            ENDIF
         END DO
         PROB = PROB + RPAR(I0)
         IF (MODE .EQ. 2) THEN
            DO I=1,NDIMX
               X(I) = X(I) + RPAR(I0)*RPAR(I0+I)
            END DO
         ENDIF
  160 CONTINUE
!
      IF (MODE .EQ. 2) THEN
         DO I=1,NDIMX
            X(I) = X(I) / PROB
         END DO
      ENDIF
!
      RETURN
!
 1000 FORMAT(' XXX - WARNING - Unrecognized range for dimension',I5,   &
     & ':',I12)
      END subroutine
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE RNDGET(RVEC,NRAND,IERR)
!
!     This subroutine creates a vector of uniform deviates on [0,1].
!     The random deviates are generated as 32 bit numbers and stored in
!     REAL*8 array X. (This is not quite clean, but I have not found
!     a decent 64 bit random number generator anywhere.)
!
!     ITYPE describes the actual generator to use. So far I have
!     implemented six generators, there could be more if necessary.
!
!     All of the generators are described in a paper by F. James,
!     A review of pseudorandom number generators, Computer Physics
!     Communications, Vol. 60, 1990, p. 329--344.
!
!     ----------------------------------------------------------------
!
!     ITYPE = 1 uses a standard multiplicative linear congruential
!               generator, for which
!
!                         Z(i+1) = (a * Z(i) + c) mod m
!
!               with a = 16807, c = 0, m = 2147483647
!
!     ITYPE = 2 uses a combination of two MLCGs developed by L'Ecuyer.
!               This generator has a period of 2**60.
!
!     ITYPE = 3 uses a Fibonacci generator developed by Marsaglia,
!               Zaman and Tsang. It has a period of 2**144.
!
!     ITYPE = 4 uses another generator written by Marsaglia and Zaman
!               with an extremely long period of 2**568.
!
!     ITYPE = 5 uses a combination of three MLCGs taken from
!               BYTE magazine. Its period is given as 6.95*10**12.
!
!     ITYPE = 6 is my implementation of the RANDU generator.
!               This is a MLCG with a = 65539 and m = 536870912.
!     -----------------------------------------------------------------
!    |    WARNING:  THIS GENERATOR IS KNOWN TO BE SERIOUSLY FLAWED.    |
!    |                USE ONLY IF YOU KNOW WHAT YOU ARE DOING!!!       |
!     -----------------------------------------------------------------
!
      USE IO_HANDLING
      USE SAMPLING_DATA
!
      PARAMETER (TWOM24 = 2.0**(-24))
!
      REAL*8    RVEC(*),ZZ
!
      IF (.NOT. INITRG) THEN
         CALL RNDSET()
         INITRG = .TRUE.
      ENDIF
!
      SELECT CASE (IRAND)
!
!
!     --------------------------------------------------------------------
!     Multiplicative linear congruential Random Number Generator.
!
!     Reference:
!           Stephen K. Park and Keith W. Miller -
!           "Random Number Generators: Good Ones Are Hard To Find",
!           Communications of the ACM, October 1988, Vol 31, No.10.
!     -------------------------------------------------------------------
!
      CASE (1)
         DO I=1,NRAND
            K = ISEED1 / 127773
            ISEED1 = 16807 * (ISEED1 - K*127773) - K * 2836
            IF (ISEED1 .LE. 0) ISEED1 = ISEED1 + 2147483647
            RVEC(I) = REAL(ISEED1)/2147483647
         END DO
         RETURN
!
!     --------------------------------------------------------------------
!     Combined multiplicative linear congruential generator
!
!     Reference:
!            P. L'Ecuyer,
!            "Efficient and portable combined random number generators",
!            Communications of the ACM, Vol. 31, 1988, p.742
!     -------------------------------------------------------------------
!
      CASE (2)
         DO I=1,NRAND
            IS1 = ISEED2(1)
            K   = IS1/53668
            IS1 = 40014*(IS1 - K*53668) - K*12211
            IF (IS1 .LT. 0) IS1 = IS1 + 2147483563
            ISEED2(1) = IS1
!
            IS2 = ISEED2(2)
            K   = IS2/52774
            IS2 = 40692*(IS2 - K*52774) - K* 3791
            IF (IS2 .LT. 0) IS2 = IS2 + 2147483399
            ISEED2(2) = IS2
!
            IZ = IS1 - IS2
            IF (IZ .LT. 1)  IZ = IZ + 2147483562
!
            RVEC(I) = REAL(IZ) * 4.656613E-10
         END DO
!
!     --------------------------------------------------------------------
!     Fibonacci random number generator
!
!     Reference:
!            G. Marsaglia, A. Zaman and W.-W. Tsang,
!            "Toward a universal random number generator",
!            Statistics and Probability Letters, Vol. 9, 1990 p.35-39.
!     -------------------------------------------------------------------
!
      CASE (3)
         DO IVEC=1,NRAND
            UNI = SEED3(I97) - SEED3(J97)
            IF (UNI .LT. 0.)  UNI = UNI + 1.0
            SEED3(I97) = UNI
            I97 = I97 - 1
            IF (I97 .EQ. 0)  I97 = 97
            J97 = J97 - 1
            IF (J97 .EQ. 0)  J97 = 97
            C = C - CD
            IF (C .LT. 0.)  C = C + CM
            UNI = UNI - C
            IF (UNI .LT. 0.) UNI = UNI + 1.0
            RVEC(IVEC) = UNI
!
!     Replace exact zeros by uniform distribution *2**-24
!
            IF (UNI .EQ. 0.)  THEN
               ZUNI = TWOM24*SEED3(2)
!
!     An exact zero here is very unlikely, but let's be safe.
!
               IF (ZUNI .EQ. 0.) ZUNI= TWOM24*TWOM24
               RVEC(IVEC) = ZUNI
            ENDIF
         END DO
!
!     --------------------------------------------------------------------
!     Add-and-carry random number generator
!
!     Reference:
!            G. Marsaglia and A. Zaman,
!            "A new class of random number generators",
!            The Annals of Applied Probability, Vol. 1, p. 462-480.
!     -------------------------------------------------------------------
!
      CASE (4)
         DO IVEC=1,NRAND
            UNI = SEED4(J24) - SEED4(I24) - CARRY
            IF (UNI .LT. 0.)  THEN
               UNI = UNI + 1.0
               CARRY = TWOM24
            ELSE
               CARRY = 0.
            ENDIF
            SEED4(I24) = UNI
            I24 = I24 - 1
            IF (I24 .EQ. 0)  I24 = 24
            J24 = J24 - 1
            IF (J24 .EQ. 0)  J24 = 24
            RVEC(IVEC) = UNI
         END DO
!
!     --------------------------------------------------------------------
!     Triple multiplicative linear congruential generator
!
!     Reference:
!            B. Wichmann and D. Hill,
!            "Building a random number generator",
!            BYTE Magazine, March 1987, p. 127-128.
!     -------------------------------------------------------------------
!
      CASE (5)
         DO I=1,NRAND
            IS1 = ISEED5(1)
            K   = ISEED5(1)/177
            IS1 = 171*(IS1 - K*177) - K* 2
            IF (IS1 .LT. 0) IS1 = IS1 + 30269
            ISEED5(1) = IS1
!
            IS2 = ISEED5(2)
            K   = ISEED5(2)/176
            IS2 = 172*(IS2 - K*176) - K*35
            IF (IS2 .LT. 0) IS2 = IS2 + 30307
            ISEED5(2) = IS2
!
            IS3 = ISEED5(3)
            K   = ISEED5(3)/178
            IS3 = 170*(IS3 - K*178) - K*63
            IF (IS3 .LT. 0) IS3 = IS3 + 30323
            ISEED5(3) = IS3
!
            ZZ = IS1/30269.D0 + IS2/30307.D0 + IS3/30323.D0
!
            RVEC(I) = ZZ - DINT(ZZ)
         END DO
!
!     --------------------------------------------------------------------
!     My implementation of RANDU...
!
!     -----------------------------------------------------------------
!    |    WARNING:  THIS GENERATOR IS KNOWN TO BE SERIOUSLY FLAWED.    |
!    |                USE ONLY IF YOU KNOW WHAT YOU ARE DOING!!!       |
!     -----------------------------------------------------------------
!
      CASE (6)
         DO I=1,NRAND
            I1 = ISEED6 / 32768
            I2 = ISEED6 - 32768 * I1
            ISEED6 = MOD(3*I1 + MOD(2*I1+3*I2,16483), 536870912)
            RVEC(I) = REAL(ISEED6)/536870912.D0
         END DO
!
         CASE DEFAULT
            WRITE (IOLOG, 1000)
            RVEC(1:NRAND) = 0.D0
            IERR = 1
         END SELECT
!
         RETURN
!
 1000 FORMAT('RNDGET: Requested non-existent random number generator!')
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE RNDSET()
!
!     This subroutine initials the random number generator, either
!     with a default seed (if ISEED = 0) or the seed supplied by
!     the user. To be on the safe side, all random number generators
!     ought to be initialed, just in case the user changes his mind
!     midstream about which one he wanted!
!
!     ----------------------------------------------------------------
!
      USE SAMPLING_DATA
!
      PARAMETER (TWOM24=2.0**(-24))
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
!
!     First the linear congruential generator
!
      IF (ISEED .EQ. 0) THEN
          ISEED1 = 12345
      ELSE
          ISEED1 = ISEED
      ENDIF
!
!     Next the combined generator of L'Ecuyer.
!     In the original version this generator requires two seeds. To be
!     consistent with the other generators, the second seed is produced
!     from the first using a (hopefully!) independent generator.
!
      IF (ISEED .EQ. 0) THEN
          ISEED2(1) = 12345
          ISEED2(2) = 98765
      ELSE
          ISEED2(1) = ISEED
          KTEMP = ISEED / 127773
          ISEED2(2) = 16807 * (ISEED - KTEMP*127773) - KTEMP * 2836
      ENDIF
!
!     Marsaglia, Zaman and Tsang's Fibonacci generator
!     The default seed gives the values used in their paper
!
      IF (ISEED .EQ. 0) THEN
          IJKL = 54217137
      ELSE
          IJKL = ISEED
      ENDIF
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I  = MOD(IJ/177, 177) + 2
      J  = MOD(IJ, 177)     + 2
      K  = MOD(KL/169, 178) + 1
      L  = MOD(KL, 169)
      DO II=1,97
         S = 0.
         T = .5
         DO JJ=1,24
            M = MOD(MOD(I*J,179)*K, 179)
            I = J
            J = K
            K = M
            L = MOD(53*L+1, 169)
            IF (MOD(L*M,64) .GE. 32)  S = S + T
            T = 0.5*T
         END DO
         SEED3(II) = S
      END DO
      C  =   362436.*TWOM24
      CD =  7654321.*TWOM24
      CM = 16777213.*TWOM24
      I97 = 97
      J97 = 33
!
!     Marsaglia and Zaman's universal generator with the long period
!
      IF (ISEED .EQ. 0) THEN
          JSEED = 314159265
      ELSE
          JSEED = ISEED
      ENDIF
      DO I=1,24
         K = JSEED/53668
         JSEED = 40014 * (JSEED - K*53668) - K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED + ICONS
         ITEMP = MOD(JSEED,ITWO24)
         SEED4(I) = REAL(ITEMP)*TWOM24
      END DO
      I24 = 24
      J24 = 10
      CARRY = 0.
      IF (SEED4(24) .LT. SEED4(14)) CARRY = TWOM24
!
!     Next the triple generator of Wichmann and Hill.
!     In the original version this generator requires three seeds.
!     Only the first one is selected by the user, the others are
!     generated from the first. The defaults were given to me
!     by Adam Berger (Princeton University).
!
      IF (ISEED .EQ. 0) THEN
          ISEED5(1) = 2907
          ISEED5(2) = 2912
          ISEED5(3) =  706
      ELSE
          KTEMP1 = ISEED / 127773
          KTEMP1 = 16807 * (ISEED  - KTEMP1*127773) - KTEMP1 * 2836
          KTEMP2 = KTEMP1 / 127773
          KTEMP2 = 16807 * (KTEMP1 - KTEMP2*127773) - KTEMP2 * 2836
          ISEED5(1) = MOD( ISEED,30000)
          ISEED5(2) = MOD(KTEMP1,30000)
          ISEED5(3) = MOD(KTEMP2,30000)
      ENDIF
!
!     Then we get to RANDU...
!
      IF (ISEED .EQ. 0) THEN
          ISEED6 = 12345
      ELSE
          ISEED6 = ISEED
      ENDIF
!
      RETURN
!
      END subroutine
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! ******************************************************************************
! *                                                                            *
! *   The following routines are translated from a Pascal library I have used  *
! *   over the years. It came from the net, but I cannot recall details...     *
! *                                                                            *
! ******************************************************************************
!
      REAL*8 FUNCTION NINV( P )
!
!-------------------------------------------------------------------------
!
!       Function:  Ninv
!
!       Purpose:   Finds percentage point of normal distribution
!                  (low accuracy)
!
!       Calling Sequence:
!
!            Perc  := Ninv( P )
!
!                 P      --- input probability
!                 Perc   --- resultant percentage point
!
!       Calls:  None
!
!       Remark:  This method provides about 6.5 decimal digits of
!                accuracy.
!
!-------------------------------------------------------------------------
!
      REAL*8 LIM, Y, P, PR, NV, NUM, DENOM, AUX, PN(5), QN(5)

      DATA LIM/ 1.0E-20/
      DATA PN/-0.322232431088,   -1.0,  -0.342242088547,      &
     &        -0.0204231210245,         -0.453642210148E-4/
      DATA QN/ 0.0993484626060, 0.588581570495, 0.531103462366,      &
     &         0.103537752850 , 0.38560700634E-2/
!
      NINV   = 0.0D0
!
      IF( P .GT. 0.5D0 ) THEN
         PR = 1.0D0 - P
      ELSE
         PR = P
      ENDIF
!
      IF ( ( PR .GE. LIM ) .AND. ( PR .NE. 0.5D0 ) ) THEN
         Y = DSQRT ( -2.0D0 * DLOG( PR ) )
         AUX = Y * PN( 5 ) + PN( 4 )
         AUX = Y *  AUX    + PN( 3 )
         AUX = Y *  AUX    + PN( 2 )
         AUX = Y *  AUX    + PN( 1 )
         NUM = AUX
         AUX = Y * QN( 5 ) + QN( 4 )
         AUX = Y *  AUX    + QN( 3 )
         AUX = Y *  AUX    + QN( 2 )
         AUX = Y *  AUX    + QN( 1 )
         DENOM = AUX
!
         NV = Y + (NUM / DENOM)
!
         IF( P .LT. 0.5D0 ) THEN
            NINV = -NV
         ELSE
            NINV = NV
         ENDIF
      ENDIF
      RETURN
      END function
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      REAL*8 FUNCTION NINV2( P )
!
!-------------------------------------------------------------------------
!
!       Function:  Ninv2
!
!       Purpose:   Finds percentage point of normal distribution
!                  (high accuracy)
!
!       Calling Sequence:
!
!            Perc  := Ninv2( P )
!
!                 P      --- input probability
!                 Perc   --- resultant percentage point
!
!       Calls:   NINV, SIGNORM
!
!       Remark:  This method provides about 12-13 decimal digits of
!                accuracy.  The method used is to improve the
!                approximation produced by NINV using a Taylor series
!                on the approximation error.
!
!-------------------------------------------------------------------------
!
      REAL*8 P, XP,P1,Z,X1,X2,X3,PHI,PI
      DATA PI/ 3.1415926535897932385 /
!
      XP  = NINV( P )
      P1  = 1.D0 - SIGNORM( XP )
      PHI = DSQRT( 1.D0 / ( 2.D0 * PI ) ) * DEXP( -( XP * XP ) / 2.D0 )
      Z   = ( P - P1 ) / PHI
!
      X3  = ( 2.D0 * ( XP * XP ) + 1.D0 ) * Z / 3.D0
      X2  = ( X3 + XP ) * Z / 2.D0
      X1  = ( ( X2 + 1.D0 ) * Z )
!
      NINV2 = XP + X1
!
      RETURN
      END function
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      REAL*8 FUNCTION SIGNORM( X )
!
!-------------------------------------------------------------------------
!
!       Function:  SigNorm
!
!       Purpose:   Evaluates normal distribution probability
!
!       Calling Sequence:
!
!            P  = SigNorm( X )
!
!                 X      --- ordinate of normal distribution
!
!                 P      --- Resultant tail probability
!
!       Calls:
!
!            Erf
!
!       Method:
!
!            The simple relationship between the error function and the
!            normal distribution is used.
!
!-------------------------------------------------------------------------
!
      REAL*8 X,SQRT2
      DATA SQRT2 / 1.4142135623730950 /
!
      IF (X .GE. 0.D0) THEN
         SIGNORM = 1.D0 - ( 1.D0 + ERF(  X / SQRT2 ) ) / 2.D0
      ELSE
         SIGNORM = 1.D0 - ( 1.D0 - ERF( -X / SQRT2 ) ) / 2.D0
      ENDIF
      RETURN
      END function
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      REAL*8 FUNCTION ERF( Z )
!--------------------------------------------------------------------------
!
!       Function:  Erf
!
!       Purpose:   Calculates the error function.
!
!       Calling sequence:
!
!          Y = Erf( Z)
!
!             Z --- The argument to the error function
!             Y --- The resultant value of the error function
!
!       Calls:  None
!
!       Method:
!
!          A rational function approximation is used, adjusted for the
!          argument size.  The approximation gives 13-14 digits of
!          accuracy.
!
!       Remarks:
!
!          The error function can be used to find normal integral
!          probabilities because of the simple relationship between
!          the error function and the normal distribution:
!
!              Norm( X ) := ( 1 + Erf(  X / Sqrt( 2 ) ) ) / 2, X >= 0;
!              Norm( X ) := ( 1 - Erf( -X / Sqrt( 2 ) ) ) / 2, X < 0.
!
!--------------------------------------------------------------------------
!
      REAL*8 A(14), B(12),Z, U, X, S, NUM, DENOM, TEMP
!
      DATA A/ 1.1283791670955,      0.34197505591854,      &
     &        0.86290601455206E-1,  0.12382023274723E-1,      &
     &        0.11986242418302E-2,  0.76537302607825E-4,      &
     &        0.25365482058342E-5, -0.99999707603738,      &
     &       -1.4731794832805,     -1.0573449601594,      &
     &       -0.44078839213875,    -0.10684197950781,      &
     &       -0.12636031836273E-1, -0.1149393366616E-8/
!
      DATA B/-0.36359916427762,     0.52205830591727E-1,      &
     &       -0.30613035688519E-2, -0.46856639020338E-4,      &
     &        0.15601995561434E-4, -0.62143556409287E-6,      &
     &        2.6015349994799,      2.9929556755308,      &
     &        1.9684584582884,      0.79250795276064,      &
     &        0.18937020051337,     0.22396882835053E-1/
!
!     Get absolute value of the argument and remember the sign
!
      X = DABS( Z )
      IF (Z .GE. 0.0) THEN
         S = 1.0
      ELSE
         S = -1.0
      ENDIF
!
!     Check if the argument is in the approximation range.
!
      IF ( Z .EQ. 0.0 ) THEN
         ERF = 0.D0
      ELSEIF( X .GE. 5.5 ) THEN
         ERF = S
      ELSE
         U = X * X
!
!      Approx. for arg <= 1.5
!
         IF ( X .LE. 1.5 ) THEN
            NUM   = A(6) + U * A(7)
            NUM   = A(5) + U * NUM
            NUM   = A(4) + U * NUM
            NUM   = A(3) + U * NUM
            NUM   = A(2) + U * NUM
            NUM   = A(1) + U * NUM
            NUM   = X * DEXP( - U) * NUM
            DENOM = B(5) + U * B(6)
            DENOM = B(4) + U * DENOM
            DENOM = B(3) + U * DENOM
            DENOM = B(2) + U * DENOM
            DENOM = B(1) + U * DENOM
            DENOM = 1.D0 + U * DENOM
            ERF   = NUM / DENOM * S
!
!      Approx. for arg > 1.5
!
         ELSE
            NUM   = A(13) + X * A(14)
            NUM   = A(12) + X * NUM
            NUM   = A(11) + X * NUM
            NUM   = A(10) + X * NUM
            NUM   = A(9)  + X * NUM
            NUM   = A(8)  + X * NUM
            DENOM = B(11) + X * B(12)
            DENOM = B(10) + X * DENOM
            DENOM = B(9)  + X * DENOM
            DENOM = B(8)  + X * DENOM
            DENOM = B(7)  + X * DENOM
            DENOM = 1.0   + X * DENOM
            TEMP  = DEXP(-U)
            ERF   = ( TEMP * NUM / DENOM + 1.0 ) * S
         ENDIF
      ENDIF
      RETURN
      END function
!-------------------------------------------------------------------------
!               CDBeta  -- Cumulative Beta Distribution
!-------------------------------------------------------------------------
!
      REAL*8 FUNCTION CDBETA( X1, Alpha, Beta, Dprec, MaxIter,      &
     &                        Cprec, Iter, Ifault)
!
!-------------------------------------------------------------------------
!
!       Function:  CDBeta
!
!       Purpose:   Evaluates CPDF of Incomplete Beta Function
!
!       Calling Sequence:
!
!            P     = CDBeta( X1, Alpha, Beta, Dprec, MaxIter,
!     *                        Cprec, Iter, Ifault)
!
!                 X1     --- Upper percentage point of PDF
!                 Alpha  --- First shape parameter
!                 Beta   --- Second shape parameter
!                 Dprec  --- Number of digits of precision required
!                 Maxitr --- Maximum number of iterations
!                 Cprec  --- Actual resulting precision
!                 Iter   --- Iterations actually used
!                 Ifault --- error indicator
!                            = 0:  no error
!                            = 1:  argument error
!
!                 P      --- Resultant probability
!
!       Calls:
!
!            ALGama
!
!       Method:
!
!            The continued fraction expansion as given by
!            Abramowitz and Stegun (1964) is used.  This
!            method works well unless the minimum of (Alpha, Beta)
!            exceeds about 70000.
!
!            An error in the input arguments results in a returned
!            probability of -1.
!
!-------------------------------------------------------------------------
!
      PARAMETER (MAXPREC = 16)
      PARAMETER (RSMALL  = 4.19E-36)
!
      REAL*8    X,X1,ALPHA,BETA,EPSZ,A,B,C,F,FX,APB,ZM,ALO,AHI,      &
     &          BLO,BHI,BOD,BEV,ZM1,D1,AEV,AOD
      INTEGER*4 CPREC,DPREC,MAXITER,ITER,IFAULT
      LOGICAL*4 QSWAP,QCONV
!
!     Initialize
!
      IF (Dprec .GT. MaxPrec) THEN
         Dprec = MaxPrec
      ELSEIF (Dprec .LE. 0) THEN
         Dprec = 1
      ENDIF
!
      Cprec  = Dprec
!
      Epsz   = 10**( -Dprec )
!
      X      = X1
      A      = Alpha
      B      = Beta
      QSwap  = .FALSE.
      CDBeta = -1.0
!
!     Check arguments
!     Error if:
!        X <= 0
!        A <= 0
!        B <= 0
!
      Ifault = 1
!
      IF( X .LE. 0.0 ) GOTO 9000
!
      IF( ( A .LE. 0.0 ) .OR. ( B .LE. 0.0 ) ) GOTO 9000
!
      CDBeta = 1.0
      Ifault = 0
!
!     If X >= 1, return 1.0 as prob
!
      IF( X .GE. 1.0 ) GOTO 9000
!
!     If X > A / ( A + B ) then swap A, B for more efficient evaluation
!
      IF( X .GT. ( A / ( A + B ) ) ) THEN
         X      = 1.0 - X
         A      = Beta
         B      = Alpha
         QSwap  = .TRUE.
      ENDIF
!
!     Check for extreme values
!
      IF( ( X .EQ. A ) .OR. ( X .EQ. B ) ) GOTO 120
      IF( DABS( A - ( X * ( A + B ) ) ) .LE. Epsz ) GOTO 120
!
         C     = ALGama( A + B ) + A * DLOG( X ) +      &
     &            B * DLOG( 1.0 - X ) - ALGama( A ) - ALGama( B ) -      &
     &            DLOG( A - X * ( A + B ) )
!
         IF( ( C .LT. -36.0 ) .AND. QSwap ) GOTO 9000
!
            CDBeta = 0.0
            IF(  C .LT. -180.0 ) GOTO 9000
!
!     Set up continued fraction expansion evaluation.
!
  120 CONTINUE
         Apb    = A + B
         Zm     = 0.0
         Alo    = 0.0
         Bod    = 1.0
         Bev    = 1.0
         Bhi    = 1.0
         Blo    = 1.0
!
       Ahi    = EXP( ALGama( Apb ) + A * DLOG( X ) +      &
     &                  B * DLOG( 1.0 - X ) - ALGama( A + 1.0 ) -      &
     &                  ALGama( B ) )
!
         F      = Ahi
!
         Iter   = 0
!
!     Continued fraction loop begins here. Evaluation continues until
!     maximum iterations are exceeded, or convergence achieved.
!
         Qconv  = .FALSE.
!
  150    CONTINUE
            Fx  = F
            Zm1 = Zm
            Zm  = Zm + 1.0
            D1  = A + Zm + Zm1
            Aev = -( A + Zm1 ) * ( Apb + Zm1 ) * X / D1 / ( D1 - 1.0 )
            Aod = Zm * ( B - Zm ) * X / D1 / ( D1 + 1.0 )
            Alo = Bev * Ahi + Aev * Alo
            Blo = Bev * Bhi + Aev * Blo
            Ahi = Bod * Alo + Aod * Ahi
            Bhi = Bod * Blo + Aod * Bhi
!
            IF (ABS( Bhi ) .LT. Rsmall) Bhi = 0.0
!
            IF( Bhi .NE. 0.0 ) THEN
               F     = Ahi / Bhi
               Qconv = ( DABS( ( F - Fx ) / F ) .LT. Epsz )
            ENDIF
            Iter  = Iter + 1
            IF ( ( Iter .LE. MaxIter ) .AND. .NOT. Qconv ) GOTO 150
!
!     Arrive here when convergence achieved, or maximum iterations exceeded.
!
      IF ( Qswap ) THEN
         CDBeta = 1.0 - F
      ELSE
         CDBeta = F
      ENDIF
!
!     Calculate precision of result
!
      IF (DABS( F - Fx ) .NE. 0.0) THEN
         Cprec = -DLOG10( DABS( F - Fx ) )
      ELSE
         Cprec = MaxPrec
      ENDIF
!
 9000 CONTINUE
         RETURN
!
      END function
!
!-------------------------------------------------------------------------
!                GamaIn -- Incomplete Gamma integral
!-------------------------------------------------------------------------
!
      REAL*8 FUNCTION GamaIn(Y, P, Dprec, MaxIter, Cprec, Iter, Ifault)
!
!-------------------------------------------------------------------------
!
!       Function:  GamaIn
!
!       Purpose:   Evaluates Incomplete Gamma integral
!
!       Calling Sequence:
!
!            PR     = GamaIn(Y, P, Dprec, MaxIter, Cprec, Iter, Ifault)
!
!                 Y      --- Gamma distrib. value
!                 P      --- Degrees of freedom
!                 Ifault --- error indicator
!
!                 PR     --- Resultant probability
!
!       Calls:
!
!            ALGama
!
!       Remarks:
!
!            Either an infinite series summation or a continued fraction
!            approximation is used, depending upon the argument range.
!
!-------------------------------------------------------------------------
!
      PARAMETER (MAXPREC = 16)
      REAL*8    Y,P,OFLO,MINEXP,F,A,B,C,TERM,PN(6),GIN,AN,RN,      &
     &          DIF,EPS
      INTEGER*4 CPREC,DPREC,MAXITER,ITER,IFAULT
      LOGICAL*4 DONE
!
      DATA OFLO/1.0D+37/,MINEXP/-87.0/
!
!     Check arguments
!
      Ifault  = 1
      GamaIn = 1.0
!
      IF ( ( Y .LE. 0.0 ) .OR. ( P .LE. 0.0 ) ) GOTO 9000
!
!     Check value of F
!
      Ifault = 0
!
      F      = P * DLOG( Y ) - ALGama( P + 1.0 ) - Y
      IF ( F .LT. MinExp ) GOTO 9000
         F = EXP( F )
         IF( F .EQ. 0.0 ) GOTO 9000
!
!     Set precision
!
         IF (Dprec .GT. MaxPrec) THEN
            Dprec = MaxPrec
         ELSEIF (Dprec .LE. 0) THEN
            Dprec = 1
         ENDIF
!
         Cprec  = Dprec
!
         Eps    = 10**( -Dprec )
!
!     Choose infinite series or continued fraction.
!
         IF ( ( Y .GT. 1.0 ) .AND. ( Y .GE. P ) ) THEN
!
!     Continued Fraction
!
            A       = 1.0 - P
            B       = A + Y + 1.0
            Term    = 0.0
            Pn( 1 ) = 1.0
            Pn( 2 ) = Y
            Pn( 3 ) = Y + 1.0
            Pn( 4 ) = Y * B
            Gin     = Pn( 3 ) / Pn( 4 )
            Done    = .FALSE.
            Iter    = 0
!
  100    CONTINUE
            Iter   = Iter + 1
            A      = A + 1.0
            B      = B + 2.0
            Term   = Term + 1.0
            An     = A * Term

            Pn( 5 ) = B * Pn( 3 ) - An * Pn( 1 )
            Pn( 6 ) = B * Pn( 4 ) - An * Pn( 2 )

            IF ( Pn( 6 ) .NE. 0.0 ) THEN
                  Rn     = Pn( 5 ) / Pn( 6 )
                  Dif    = ABS( Gin - Rn )
                  IF ( Dif .LE. Eps ) THEN
                     IF ( Dif .LE. ( Eps * Rn ) ) Done = .TRUE.
                  ENDIF
                  Gin    = Rn
            ENDIF
!
            Pn( 1 ) = Pn( 3 )
            Pn( 2 ) = Pn( 4 )
            Pn( 3 ) = Pn( 5 )
            Pn( 4 ) = Pn( 6 )
!
            IF( ABS( Pn( 5 ) ) .GE. Oflo ) THEN
                Pn( 1 ) = Pn( 1 ) / Oflo
                Pn( 2 ) = Pn( 2 ) / Oflo
                Pn( 3 ) = Pn( 3 ) / Oflo
                Pn( 4 ) = Pn( 4 ) / Oflo
            ENDIF
            IF (ITER .LE. MAXITER .AND. .NOT. DONE) GOTO 100
!
         Gin    = 1.0 - ( F * Gin * P )
         GamaIn = Gin
!
!     Calculate precision of result
!
         IF (Dif .NE. 0.0) THEN
            Cprec = -DLOG10( Dif )
         ELSE
            Cprec = MaxPrec
         ENDIF
      ELSE
!
!     Infinite series
!
         Iter    = 0
         Term    = 1.0
         C       = 1.0
         A       = P
         Done    = .FALSE.
!
  150    CONTINUE
            A       = A + 1.0
            Term    = Term * Y / A
            C       = C + Term
            Iter    = Iter + 1
            IF (( Iter .LE. MaxIter ) .AND. ( ( Term / C ) .GT. Eps ))      &
     &         GOTO 150
!
         GamaIn = C * F
      ENDIF
!
!     Calculate precision of result
!
         Cprec = -DLOG10( Term / C )
!
 9000 CONTINUE
      RETURN
!
      END function
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!-------------------------------------------------------------------------
!               ALGama  -- Logarithm of Gamma Distribution
!-------------------------------------------------------------------------
!
      REAL*8 FUNCTION ALGama( Arg )
!
!-------------------------------------------------------------------------
!
!       Function:  ALGama
!
!       Purpose:  Calculates Log (base E) of Gamma function
!
!       Calling Sequence:
!
!            Val = ALGama( Arg )
!
!                 Arg   --- Gamma distribution parameter (Input)
!                 Val   --- output Log Gamma value
!
!       Calls:  None
!
!       Called By:  Many (CDBeta, etc.)
!
!       Remarks:
!
!            Minimax polynomial approximations are used over the
!            intervals (-inf,0), (0,.5), (.5,1.5), (1.5,4.0),
!            (4.0,12.0), (12.0,+inf).
!
!            See Hart et al, "Computer Approximations",
!            Wiley(1968), p. 130F, and also
!
!            Cody and Hillstrom, "Chebyshev approximations for
!            the natural logarithm of the Gamma function",
!            Mathematics of Computation, 21, April, 1967, P. 198F.
!
!
!            There are some machine-dependent constants used --
!
!               Rmax   --- Largest value for which ALGamma
!                          can be safely computed.
!               Rinf   --- Largest floating-point number.
!               Zeta   --- Smallest floating-point number
!                          such that (1 + Zeta) = 1.
!
!-------------------------------------------------------------------------
!
      REAL*8    ARG, RARG, ALINC, SCALE, TOP, BOT, FRAC, ALGVAL,      &
     &          PI,P(29),Q(24),RMAX,RINF,ZETA,Xln2sp
      INTEGER*4 I,IOF,ILO,IHI
      LOGICAL*4 QMINUS,QDOIT
!
      DATA PI/ 3.1415926535897932385 /
!
      DATA P/      4.12084318584770E+00 ,  8.56898206283132E+01 ,      &
     &             2.43175243524421E+02 , -2.61721858385614E+02 ,      &
     &            -9.22261372880152E+02 , -5.17638349802321E+02 ,      &
     &            -7.74106407133295E+01 , -2.20884399721618E+00 ,      &
     &             5.15505761764082E+00 ,  3.77510679797217E+02 ,      &
     &             5.26898325591498E+03 ,  1.95536055406304E+04 ,      &
     &             1.20431738098716E+04 , -2.06482942053253E+04 ,      &
     &            -1.50863022876672E+04 , -1.51383183411507E+03 ,      &
     &            -1.03770165173298E+04 , -9.82710228142049E+05 ,      &
     &            -1.97183011586092E+07 , -8.73167543823839E+07 ,      &
     &             1.11938535429986E+08 ,  4.81807710277363E+08 ,      &
     &            -2.44832176903288E+08 , -2.40798698017337E+08 ,      &
     &             8.06588089900001E-04 , -5.94997310888900E-04 ,      &
     &             7.93650067542790E-04 , -2.77777777688189E-03 ,      &
     &             8.33333333333330E-02   /
!
      DATA Q /     1.00000000000000E+00 ,  4.56467718758591E+01 ,      &
     &             3.77837248482394E+02 ,  9.51323597679706E+02 ,      &
     &             8.46075536202078E+02 ,  2.62308347026946E+02 ,      &
     &             2.44351966250631E+01 ,  4.09779292109262E-01 ,      &
     &             1.00000000000000E+00 ,  1.28909318901296E+02 ,      &
     &             3.03990304143943E+03 ,  2.20295621441566E+04 ,      &
     &             5.71202553960250E+04 ,  5.26228638384119E+04 ,      &
     &             1.44020903717009E+04 ,  6.98327414057351E+02 ,      &
     &             1.00000000000000E+00 , -2.01527519550048E+03 ,      &
     &            -3.11406284734067E+05 , -1.04857758304994E+07 ,      &
     &            -1.11925411626332E+08 , -4.04435928291436E+08 ,      &
     &            -4.35370714804374E+08 , -7.90261111418763E+07   /
!
      DATA RMAX/ 1.67E+38 /, RINF/ 1.67E+38 /, ZETA/ 1.0E-16 /,        &
     &     Xln2sp / 9.18938533204673E-01   /
!
!     Initialize
!
      Algval = Rinf
      Scale  = 1.0
      Alinc  = 0.0
      Frac   = 0.0
      Rarg   = Arg
      Iof    = 1
      Qminus = .FALSE.
      Qdoit  = .TRUE.
!
!     Adjust for negative argument
!
      IF ( Rarg .LT. 0.0 ) THEN
         Qminus = .TRUE.
         Rarg   = -Rarg
         Top    = Int( Rarg )
         Bot    = 1.0
!
         IF( ( INT( Top / 2.0 ) * 2.0 ) .EQ. 0.0 ) Bot = -1.0
!
         Top    = Rarg - Top
!
         IF( Top .EQ. 0.0 ) THEN
            Qdoit = .FALSE.
         ELSE
            Frac   = Bot * PI / SIN( Top * PI )
            Rarg   = Rarg + 1.0
            Frac   = DLOG( ABS( Frac ) )
         ENDIF
      ENDIF
!
!     Choose approximation interval based upon argument range
!
      IF ( Rarg .EQ. 0.0 ) THEN
         Qdoit = .FALSE.
!
      ELSEIF ( Rarg .LE. 0.5 ) THEN
         Alinc  = -DLOG( Rarg )
         Scale  = Rarg
         Rarg   = Rarg + 1.0
         IF( Scale .LT. Zeta ) THEN
            Algval = Alinc
            Qdoit  = .FALSE.
         ENDIF
      ELSEIF ( Rarg .LE. 1.5 ) THEN
         Scale = Rarg - 1.0
      ELSEIF ( Rarg .LE. 4.0 ) THEN
         Scale  = Rarg - 2.0
         Iof    = 9
      ELSEIF ( Rarg .LE. 12.0 ) THEN
         Iof = 17
      ELSEIF ( Rarg .LE. RMAX   ) THEN
         Alinc  = ( Rarg - 0.5 ) * DLOG( Rarg ) - Rarg + Xln2sp
         Scale  = 1.0 / Rarg
         Rarg   = Scale * Scale
         Top    = P( 25 )
         DO I=26,29
             Top    = Top * Rarg + P( I )
         END DO
         Algval = Scale * Top + Alinc
         Qdoit  = .FALSE.
      ENDIF
!
!     Common evaluation code for Arg <= 12.  Horner's method is used,
!     which seems to give better accuracy than continued fractions.
!
      IF (Qdoit) THEN
         Ilo = Iof + 1
         Ihi = Iof + 7
         Top = P( Iof )
         Bot = Q( Iof )
!
         DO I=Ilo,Ihi
            Top = Top * Rarg + P( I )
            Bot = Bot * Rarg + Q( I )
         END DO
!
         Algval = Scale * ( Top / Bot ) + Alinc
!
      ENDIF
!
      IF( Qminus ) Algval = Frac - Algval
      ALGama = Algval
      RETURN
!
      END function
!
! ==================================================================
!
      END MODULE SAMPUTIL
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      MODULE TREE_BUILDER
!
!     This module complements the SMPS reader. Its purpose is to build
!     a scenario tree from distribution information in one of three ways:
!       o by replication of data items for discrete random variables
!       o by sampling from continuous random variables
!       o by simplifying an existing scenario tree
!
! ---------------------------------------------------------------
!
!     Start by defining some work arrays that are shared between routines
!
      TYPE SMPS_SCALAR
          LOGICAL*1 RHS, BOUNDS, COST, PLQ, CHANCE, ICC, VAR
      END TYPE SMPS_SCALAR
!
      TYPE (SMPS_SCALAR), ALLOCATABLE, DIMENSION (:) :: COPIED_SCALARS
!
      TYPE SMPS_MATRIX
          LOGICAL*1 A_MATRIX, Q_MATRIX
      END TYPE SMPS_MATRIX
!
      TYPE (SMPS_MATRIX), ALLOCATABLE, DIMENSION(:,:) :: COPIED_MATRICES
!
      REAL (KIND=8),      ALLOCATABLE, DIMENSION (:) :: HIST
      REAL (KIND=8),      ALLOCATABLE, DIMENSION (:) :: WORK
      REAL (KIND=8),      ALLOCATABLE, DIMENSION (:) :: CPROB
      REAL (KIND=8),      ALLOCATABLE, DIMENSION (:) :: QPROB
      REAL (KIND=8),      ALLOCATABLE, DIMENSION (:) :: XRNDPAR
!
      INTEGER, ALLOCATABLE, DIMENSION (:,:) :: NBLK
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: INSTAGE, KREF, ITAKE
!
      CONTAINS
!
      SUBROUTINE MAKEDE(IERR)
!
!--------------------------------------------------------------------------
!
!       Routine:   MAKEDE
!
!       Purpose:   This routine builds a deterministic equivalent
!                  from a number of discrete distributions.
!
!       Calling sequence:
!
!           CALL MAKEDE(IERR)
!
!       IERR      Return status
!            = 0  Normal termination
!            = 1  Error condition
!
!       Calls:  SNGLBRCH.
!
!       Method:
!
!          Deterministic equivalents exist for a large class of problems.
!          In this first implementation only recourse problems with finite
!          (discrete) distributions are considered. We check for the
!          presence (or absence) of probabilistic constraints and count
!          the number of scenarios (to warn the user in case the situation
!          threatens to get out of hand). Values for the first realization
!          of each random variable are stored into the base scenario;
!          subsequent realizations are stored into replicated scenarios;
!          the number of replications depends on the number of realizations
!          of previously defined random variables in the same period.
!
!       Remarks: None.
!
! --------------------------------------------------------------------------
!          This version dated 4 December 2002. Written by Gus Gassmann.
! --------------------------------------------------------------------------
!
      USE SMPS_DATA
      USE IO_HANDLING
      USE TIMER_STAT
      USE MEMORY_EXTENDER
      implicit real*8 (a-h,o,p,r-z), integer*4(i-n), character*1(q)
!
      CHARACTER*(NAMLEN) DTEMP, DBLANK
      CHARACTER*8 DISCRETE, DCONST, DTEMP8, DBLANK8, F1100
      CHARACTER*1 QTEMP(NAMLEN), QBLANK(NAMLEN), QTEMP8(8)
      REAL*8      LN2
!
      EQUIVALENCE (DBLANK,QBLANK),(DTEMP,QTEMP),(QTEMP8,DTEMP8)
!
      DATA DISCRETE/'DISCRETE'/, DCONST/'CONSTANT'/, DBLANK8/'        '/
      DATA QBLANK/NAMLEN*' '/
!
!     Set up and allocate temporary arrays
!
      REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: REALS
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: NRCUR, NREAL, RVPOS
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: LNKFWD,  LNKBAK
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: LSTRV,   FSTRV
      INTEGER,       ALLOCATABLE, DIMENSION (:) :: NREP
!
      ALLOCATE (REALS(NPER), CPROB(NPER), INSTAGE(NPER),      &
     &          LSTRV(NPER), FSTRV(NPER), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
      ALLOCATE (NRCUR(NRVAR),  RVPOS(NRVAR),  NREAL(NRVAR),      &
     &          LNKFWD(NRVAR), LNKBAK(NRVAR), NREP(NRVAR), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
      ALLOCATE (COPIED_SCALARS(NPER), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
      ALLOCATE (COPIED_MATRICES(NPER,NPER), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
      ALLOCATE (WORK(NSTELM), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
      ALLOCATE (NBLK(0:NPER,2), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
      ALLOCATE (KREF(NPER), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
!
!     Initialization
!
      TIMESET = .FALSE.
      CALL STIMER(TIME0)
!
      DO I=1,NPER
         REALS(I) = 0.D0
         FSTRV(I) = 0
         LSTRV(I) = 0
         CPROB(I) = 1.D0
         NBLK(I-1,1) = KDATA(I)
         IF (ALLOCATED(KDATQ)) NBLK(I-1,2) = KDATQ(I)
      END DO
      NBLK(NPER,1) = LASTBA
      NBLK(NPER,2) = LASTQB
!
      NSCHAR = 0
      KSCNAM(1) = 0
!
!     Next we go through the random variables one at a time.
!     The types are checked: Only DISCRETE and CONSTANT entries are allowed;
!     as soon as anything else is encountered, pull the plug. Estimate
!     storage requirements along the way and set up the history vector.
!
      NRCOMP = 0
      DO I=1,NRVAR
         IF (SDIST(I) .EQ. DISCRETE .OR. SDIST(I) .EQ. DCONST) THEN
            NRCOMP = NRCOMP + MDIST(8*I-7)
         ELSE
            GOTO 930
         ENDIF
      END DO
      IF (NRCOMP .EQ. 0) GOTO 990
!
      ALLOCATE (HIST(NRCOMP), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
!
!     Process the random variables again. Discrete
!     random variables are put in a doubly linked list, sorted by period
!     (last period first), the first realization is stored, and the position
!     of the random variable is recorded. Constant values are simply stored.
!     We also compute the number of branches in each period. This can easily
!     be a very large number; thus we use a log transformation (base 2).
!
      LN2 = DLOG(2.D0)
      KDIM = 0
      INBAK = 0
      DO I=1,NRVAR
         IF (SDIST(I) .EQ. DISCRETE) THEN
            NR = MDIST(8*I-3) / (MDIST(8*I-7)+1)
            IP = MDIST(8*I)
            IF (IP .GT. 0) THEN
               REALS(IP) = REALS(IP) + LOG(REAL(NR))/LN2
               IF (REALS(IP) .GT. 32.D0) GOTO 920
                  NRCUR(I) = 1
                  DO K=1,MDIST(8*I-7)
                     HIST(KDIM+K) = RANDPAR(MDIST(8*I-4)+K+1)
                  END DO
                  CPROB(IP) = CPROB(IP) * RANDPAR(MDIST(8*I-4)+1)
!
                  IF (FSTRV(IP) .NE. 0) THEN
                     LNKFWD(I) = LNKFWD(LSTRV(IP))
                     LNKFWD(LSTRV(IP)) = I
                     LNKBAK(I) = LSTRV(IP)
                     LSTRV(IP) = I
                     IF (LNKFWD(I) .NE. 0) LNKBAK(LNKFWD(I)) = I
                  ELSE
                     FSTRV(IP) = I
                     LSTRV(IP) = I
                     IF (INBAK .LT. IP) THEN
                        IF (INBAK .EQ. 0) THEN
                           LNKFWD(I) = 0
                           LNKBAK(I) = 0
                           INFWD     = I
                        ELSE
                           LNKBAK(I) = LSTRV(INBAK)
                           LNKFWD(LNKBAK(I)) = I
                           LNKFWD(I) = 0
                        ENDIF
                        INBAK = IP
                     ELSE
                        DO JP=IP+1,NPER
                           IF (FSTRV(JP) .GT. 0) GOTO 130
                        END DO
!
  130                CONTINUE
                        LNKBAK(I) = LNKBAK(FSTRV(JP))
                        LNKBAK(FSTRV(JP)) = I
                        LNKFWD(I) = FSTRV(JP)
                        IF (LNKBAK(I) .NE. 0) THEN
                           LNKFWD(LNKBAK(I)) = I
                        ELSE
                           INFWD = I
                        ENDIF
                     ENDIF
                  ENDIF
                  NREAL(I) = NR
                  RVPOS(I) = KDIM
            ENDIF
!
         ELSE
            HIST(KDIM+1) = RANDPAR(MDIST(8*I-4)+1)
         ENDIF
         KDIM = KDIM + MDIST(8*I-7)
      END DO
!
!     Estimate the storage requirements; if too large, pull the plug.
!
      ND = NPER + 1
      NC = MAXCOL + NODEINFO(1)%NCOL
      NR = MAXROW + NODEINFO(1)%NROW
      IF (MARKOV) THEN
         NB = 2*NPER
      ELSE
         NB = NPER*(NPER+1)/2 + 1
      ENDIF
      IF (HAVE_GC) NB = NB + NPER*(NPER-1)/2
!
      DO I=2,NPER
         REALS(I) = REALS(I) + REALS(I-1)
         IF (REALS(I) .GT. 32.D0) GOTO 920
         NREALS = NINT(2**REALS(I))
         IF (NREALS .GT. MXNODE) GOTO 940
         ND = ND + NREALS
         RO = REALS(I) + LOG(REAL(NODEINFO(I)%NROW))/LN2
         CO = REALS(I) + LOG(REAL(NODEINFO(I)%NCOL))/LN2
         IF (RO .GT. 32.D0 .OR. CO .GT. 32.D0) GOTO 980
         NR = NR + NREALS*NODEINFO(I)%NROW
         NC = NC + NREALS*NODEINFO(I)%NCOL
         IF (MARKOV) THEN
            NB = NB + 2*NREALS
         ELSE
            NB = NB + I*NREALS
         ENDIF
         IF (HAVE_GC) NB = NB + NREALS*(I-1)
         IF (NECHO1 .GE. 2) THEN
            WRITE (IOLOG, 1470) I,NREALS,ND,NR,NC,NB
         ENDIF
      END DO
!
      IF (NR .GT. MXROWS) THEN
         WRITE (IOLOG, 1440) NR
         IERR = 1
         GOTO 999
      ENDIF
      IF (NC .GT. MXCOLS) THEN
         WRITE (IOLOG, 1450) NC
         IERR = 1
         GOTO 999
      ENDIF
      IF (ND .GT. MXNODE) THEN
         WRITE (IOLOG, 1460) ND
         IERR = 1
         GOTO 999
      ENDIF
      IF (NB .GT. MXABLK) THEN
         WRITE (IOLOG, 1490) NB
         IERR = 1
         GOTO 999
      ENDIF
!
!     Place the first scenario, which must duplicate the entire base path.
!     ====================================================================
!
      KSCEN = 1
      DO J=1,NPER-1
         INSTAGE(J) = NODES + J
      END DO
      DO J=1,NPER
         KREF(J) = NODES + J
      END DO
      CALL SNGLBRCH(1,0,IERR)
!
!     Generate the name for this scenario
!
      NOFF = 1
      DTEMP = DBLANK
      NV = INFWD
      QTEMP(1) = '1'
  220 CONTINUE
         NV = LNKFWD(NV)
         IF (NV .GT. 0) THEN
            QTEMP(NOFF+1) = '.'
            QTEMP(NOFF+2) = '1'
            NOFF = NOFF + 2
            GOTO 220
         ENDIF
!
         IF (NOFF .GT. IZCHSC) THEN
            IF (.NOT. EXTMEM_CHSC(NOFF,IERR)) GOTO 960
         ENDIF
         DO I=1,NOFF
            QSCNAM(I) = QTEMP(I)
         END DO
         NSCHAR = + NOFF
         LEAFND(1) = NODES
         KSCNAM(2) = NSCHAR
!
!     Now generate the remaining scenarios
!     ====================================
!
         IF (INBAK .EQ. 0) GOTO 900
!
!     Each realization of a random variable leads to multiple replications.
!     The number of replications is computed and stored in NREP.
!
         KRV   = LSTRV(INBAK)
         NREPS = 1
  205 CONTINUE
         NREP(KRV) = NREPS
         KRV = LNKBAK(KRV)
         IF (KRV .GT. 0) THEN
            NREPS = NREPS*NREAL(KRV)
            GOTO 205
         ENDIF
!
      NV = LSTRV(INBAK)
  300 CONTINUE
         IF (NRCUR(NV) .EQ. NREAL(NV)) THEN
            NV = LNKBAK(NV)
            IF (NV .EQ. 0) GOTO 900
               GOTO 300
         ELSE
            DO K=1,MDIST(8*NV-7)
               HIST(RVPOS(NV)+K) = RANDPAR(MDIST(8*NV-4) + K + 1 +      &
     &                                    (MDIST(8*NV-7)+1)*NRCUR(NV))
            END DO
            NRCUR(NV) = NRCUR(NV) + 1
            IPER = MDIST(8*NV) - 1
            IF (IPER .EQ. 0) GOTO 950
         ENDIF
         NX = LSTRV(INBAK)
!
  320 CONTINUE
      IF (NX .NE. NV) THEN
         NRCUR(NX) = 1
         IF (SDIST(NX) .EQ. DISCRETE) THEN
            DO K=1,MDIST(8*NX-7)
               HIST(RVPOS(NX)+K) = RANDPAR(MDIST(8*NX-4)+K+1)
            END DO
         ENDIF
         NX = LNKBAK(NX)
         GOTO 320
      ENDIF
!
!     The next scenario has been identified
!
            NV = LSTRV(INBAK)
            KSCEN = KSCEN + 1
!
!     Identify the branching node
!
            NODE0 = INSTAGE(IPER)
            DO J=IPER+1,NPER-1
               INSTAGE(J) = NODES + J - IPER
               CPROB(J) = 1.D0
            END DO
!
!     Recompute the conditional probabilities
!
            NX = LSTRV(INBAK)
            CPROB(IPER+1:NPER) = 1.D0
!
  360    CONTINUE
            IP = MDIST(8*NX)
            CPROB(IP) = CPROB(IP) * RANDPAR(MDIST(8*NX-4) + 1 +      &
     &                                (MDIST(8*NX-7)+1)*(NRCUR(NX)-1))
            IF (NX .NE. FSTRV(IPER+1)) THEN
               NX = LNKBAK(NX)
               GOTO 360
            ENDIF
!
!     Compute the reference nodes and store the scenario
!
            KRV = INFWD
  365       CONTINUE
               IF (NRCUR(KRV) .EQ. 1) THEN
                  IF (KRV .EQ. 0) THEN
                     IERR = 1
                     GOTO 999
                  ENDIF
                  KRV = LNKFWD(KRV)
                  GOTO 365
               ENDIF
               LEAF0 = 1
  366          CONTINUE
                  KRV = LNKFWD(KRV)
                  IF (KRV .GT. 0) THEN
                     LEAF0 = LEAF0 + (NRCUR(KRV)-1)*NREP(KRV)
                     GOTO 366
                  ENDIF
!
               IP = NPER
               NC = LEAFND(LEAF0)
  367       CONTINUE
               KREF(IP) = NC
               IF (IP .GT. IPER) THEN
                  IP = IP - 1
                  NC = NODEINFO(NC)%IANCTR
                  GOTO 367
               ENDIF
!
            CALL SNGLBRCH(IPER+1,NODE0,IERR)
!
!     Generate the name for this scenario
!
            NOFF = 0
            DTEMP = DBLANK
            DTEMP8 = DBLANK8
            NX = INFWD
  370    CONTINUE
            NDIGIT = IFIX(ALOG10(FLOAT(NRCUR(NX)))) + 1
            WRITE (F1100, 1200) NDIGIT
            WRITE (DTEMP8, F1100) NRCUR(NX)
!
            DO I=1,NDIGIT
               QTEMP(NOFF+I) = QTEMP8(I)
            END DO
            NOFF = NOFF + NDIGIT
            NX = LNKFWD(NX)
            IF (NX .GT. 0) THEN
               NOFF = NOFF + 1
               QTEMP(NOFF) = '.'
               GOTO 370
            ENDIF
!
            IF (KSCEN .GE. IZSCEN) THEN
               IF (.NOT. EXTMEM_SCEN(KSCEN+1,IERR)) GOTO 960
            ENDIF
!
            MNCHSC = NOFF + NSCHAR
            IF (MNCHSC .GT. IZCHSC) THEN
               IF (.NOT. EXTMEM_CHSC(MNCHSC,IERR)) GOTO 960
            ENDIF
            DO I=1,NOFF
               QSCNAM(NSCHAR+I) = QTEMP(I)
            END DO
            NSCHAR = NSCHAR + NOFF
            LEAFND(KSCEN) = NODES
            KSCNAM(KSCEN+1) = NSCHAR
!
!     Get the next scenario
!
            GOTO 300
!
!     We have processed all scenarios
!
  900 CONTINUE
      MAXRHS  = LASTR
      TREEFMT = 2
!
      CALL STIMER(TIME1)
      WRITE (IOLOG, 1000) TIME1
      GOTO 999
!
  920 CONTINUE
         WRITE (IOLOG, 1410)
         IERR = 1
         GOTO 999
!
  930 CONTINUE
         WRITE (IOLOG, 1420)
         IERR = 1
         GOTO 999
!
  940 CONTINUE
         WRITE (IOLOG, 1430) I, NREALS
         IERR = 1
         GOTO 999
!
  950 CONTINUE
         WRITE (IOLOG, 1600)
         IERR = 1
         GOTO 999
!
  960 CONTINUE
         WRITE (IOLOG, 1700)
         IERR = 1
         GOTO 999
!
  980 CONTINUE
         WRITE (IOLOG, 1480)
         IERR = 1
         GOTO 999
!
  990 CONTINUE
         WRITE (IOLOG, 1300)
  999 CONTINUE
         DEALLOCATE(REALS,CPROB,HIST,NRCUR,RVPOS,LNKFWD,LNKBAK,LSTRV,        &
     &              FSTRV,NREAL,INSTAGE,COPIED_SCALARS,COPIED_MATRICES,      &
     &              NBLK,WORK,KREF,STAT=IER)
         IF (IER .NE. 0) THEN
            WRITE (IOLOG, 1700)
            IERR = 1
         ENDIF
         RETURN
!
 1000 FORMAT(' Time to make deterministic equivalent:',F10.3,' sec.')
 1200 FORMAT('(I',I1,')')
 1300 FORMAT(' XXX - WARNING - Problem is deterministic!')
 1400 FORMAT(' XXX - WARNING - Can not find deterministic equivalent:',      &
     &     /,'                 Problem contains chance constraints.')
 1410 FORMAT(' XXX - WARNING - Can not find deterministic equivalent:',      &
     &     /,'                 More than 2**32 scenarios.')
 1420 FORMAT(' XXX - WARNING - Can not find deterministic equivalent:',      &
     &     /,'         Problem contains continuous random variables.')
 1430 FORMAT(' XXX - WARNING - Can not find deterministic equivalent:',      &
     &     /,'                 Period',i3,' contains',i8,' nodes.')
 1440 FORMAT(' XXX - WARNING - Can not find deterministic equivalent:',      &
     &       '                 Problem contains too many rows:',I8,'.')
 1450 FORMAT(' XXX - WARNING - Can not find deterministic equivalent:',      &
     &       '                 Problem contains too many cols:',I8,'.')
 1460 FORMAT(' XXX - WARNING - Can not find deterministic equivalent:',      &
     &       '                 Problem contains too many nodes:',I8,'.')
 1470 FORMAT(' Period',I4,' contains',i8,' nodes; total nodes =',i8,/,       &
     &       ' Number of rows so far:',i8,'; number of columns:',i8,/,       &
     &       ' Number of A-matrix blocks so far:',i8)
 1480 FORMAT(' XXX - WARNING - Can not find deterministic equivalent:',      &
     &     /,'                 More than 2**32 rows or columns.')
 1490 FORMAT(' XXX - WARNING - Can not find deterministic equivalent:',      &
     &     /,'                 Number of A-matrix blocks:',I8,'.')
 1500 FORMAT(' XXX - WARNING - Unrecognized type of stochastic element')
 1600 FORMAT(' XXX - WARNING - Can not find deterministic equivalent:',      &
     &     /,'                 Stochastic element in stage 1.')
 1700 FORMAT(' XXX -  FATAL  - Memory handling error in MAKEDE.')
      END SUBROUTINE
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE SNGLBRCH(IPER1,NODE0,IERR)
!
!--------------------------------------------------------------------------
!     Routine:   SNGLBRCH
!
!     Purpose:   To generate a single scenario and add it to the event tree.
!
!     Calling sequence:
!
!         CALL SNGLBRCH(IPER1,NODE0,IERR)
!
!              IPER1  - Starting period for this scenario
!              NODE0  - Ancestor node of this scenario
!              IERR   - Error status
!
!     Calls:   None.
!
!     Method:
!
!     This subroutine generates a single scenario, using information
!     read earlier from the input files and adds it to the event tree.
!     The new scenario branches from an existing scenario in NODE0 and starts
!     in IPER1. In other words, the first node on this scenario has NODE0
!     as its ancestor and is part of period IPER1.
!
!     Representation of random variables and stochastic elements is as in
!     routine INCONT, that is, using arrays MDIST and SDIST to describe the
!     distributions of the random variables, KSTOCH to describe the location
!     of the stochastic elements, and the D-matrix (DMTX,IDMTX and LDMTX)
!     to hold the two together.
!
!     Data structures:
!     ----------------
!
!     1. The distributions are stored as strings in array SDIST.
!
!     2. The dimension, number and locations of parameters and the stage
!        in which the random variable becomes known are stored in array MDIST:
!        MDIST(8*IRV-7) gives the dimension of the random variable
!        MDIST(8*IRV-6) gives the offset for the first integer parameter
!                       (parameter values are in array IRNDPAR)
!        MDIST(8*IRV-5) gives the number of integer parameters
!        MDIST(8*IRV-4) gives the offset for the first real parameter
!                       (parameter values are in array RANDPAR)
!        MDIST(8*IRV-3) gives the number of real parameters
!        MDIST(8*IRV-2) gives the offset for the first parameter referenced
!                       by location (parameter values are in array LRNDPAR)
!        MDIST(8*IRV-1) gives the number of parameters referenced by location
!        MDIST(8*IRV  ) gives the number of the stage
!
!
!     3. The location and process mode of the stochastic elements are stored
!        in array KSTOCH in groups of four. The first entry gives the type
!        of random element, the second its block (or node), the third entry
!        specifies describes the location within the relevant array, the
!        fourth whether the values are to be treated as replacing the core
!        file information or if they are to be treated as additive or
!        multiplicative perturbations.
!
!     Types of random elements are as follows (can be extended as needed):
!
!        KSTOCH(4*NSTELM-3) = 1 for stochastic entry in an A-matrix block
!        KSTOCH(4*NSTELM-3) = 2 for stochastic cost
!        KSTOCH(4*NSTELM-3) = 3 for stochastic RHS
!        KSTOCH(4*NSTELM-3) = 4 for stochastic range
!        KSTOCH(4*NSTELM-3) = 5 for stochastic lower bound
!        KSTOCH(4*NSTELM-3) = 6 for stochastic upper bound
!        KSTOCH(4*NSTELM-3) = 7 for stochastic entry in a Q-matrix block
!        KSTOCH(4*NSTELM-3) = 11 for stochastic upper and lower bounds ('FX' type)
!        KSTOCH(4*NSTELM-3) = 12 for stochastic demand (network option)
!        KSTOCH(4*NSTELM-3) = 13 for stochastic supply (network option)
!
!     ---------------------------------------------------------------
!     This version dated 7 December 2002. Written by Gus Gassmann.
!     ---------------------------------------------------------------
!
      USE SMPS_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      implicit real*8 (a-h,o,p,r-z), integer*4(i-n), character*1(q)
!
      LOGICAL   COPIED
!
! ----------------------------------------------------------------------
!
!     Start with some initializations
!
      DO I=1,NPER
         COPIED_SCALARS(I)%RHS    = .FALSE.
         COPIED_SCALARS(I)%BOUNDS = .FALSE.
         COPIED_SCALARS(I)%COST   = .FALSE.
         COPIED_SCALARS(I)%PLQ    = .FALSE.
         COPIED_SCALARS(I)%CHANCE = .FALSE.
         COPIED_SCALARS(I)%ICC    = .FALSE.
         COPIED_SCALARS(I)%VAR    = .FALSE.
         DO J=1,NPER
            COPIED_MATRICES(I,J)%A_MATRIX = .FALSE.
            COPIED_MATRICES(I,J)%Q_MATRIX = .FALSE.
         END DO
      END DO
!
!     Now compute the product V = D*Y
!
      DO J=1,NSTELM
         WORK(J) = 0.D0
      END DO
!
      KDIM = 0
      DO J=1,NRVAR
         DO K=1,MDIST(8*J-7)
            DO L=LDMTX(KDIM+K),LDMTX(KDIM+K+1)-1
               WORK(IDMTX(L)) = WORK(IDMTX(L)) + DMTX(L)*HIST(KDIM+K)
            END DO
         END DO
         KDIM = KDIM + MDIST(8*J-7)
      END DO
!
      MNNODE = NODES + NPER - IPER1 + 1
      IF (MNNODE .GT. IZNODE) THEN
         IF (.NOT. EXTMEM_NODE(MNNODE,IERR)) GOTO 9200
      ENDIF
!
!     Now store all the information into the data structure
!
      DO I=IPER1,NPER
!
!     Start by creating the node and placing it into the event tree
!
         NODES = NODES + 1
         I0 = KREF(I)
!
!     The first scenario must be copied in its entirety,
!     so as to allow node aggregation to work properly
!
         IF (NODE0 .EQ. 0) THEN
            IF (I .GT. IPER1) THEN
               NODEINFO(NODES  )%IANCTR = NODES - 1
               NODEINFO(NODES-1)%NDESC  = 1
               NODEINFO(NODES-1)%IDESC  = NODES
            ELSE
               NODEINFO(NODES  )%IANCTR = 0
            ENDIF
            NODEINFO(NODES)%IBROTH = 0
            NODEINFO(NODES)%NDESC  = 0
            NODEINFO(NODES)%IDESC  = 0
            NODEINFO(NODES)%PROB   = CPROB(I)
            IRNGE0(I)     = NODES
!
!     Initialize pointers
!
            KDATA(NODES) = LASTBA
            NODEINFO(NODES)%KROW   = NODEINFO(NODES-1)%KROW +      &
     &                               NODEINFO(NODES-1)%NROW
            NODEINFO(NODES)%KCOL   = NODEINFO(NODES-1)%KCOL +      &
     &                               NODEINFO(NODES-1)%NCOL + 1
            NODEINFO(NODES)%KCOST  = NODEINFO(I)%KCOST
            NODEINFO(NODES)%KRHS   = NODEINFO(I)%KRHS
            NODEINFO(NODES)%KBOUND = NODEINFO(I)%KBOUND
            NODEINFO(NODES)%KNAMES = NODEINFO(I)%KNAMES
            NODEINFO(NODES)%VARPTR = NODEINFO(I)%VARPTR
            NODEINFO(NODES)%NROW   = NODEINFO(I)%NROW
            NODEINFO(NODES)%NCOL   = NODEINFO(I)%NCOL
            IF (ALLOCATED(KDATQ)) KDATQ(NODES) = LASTQB
            IF (LCHANCE .GT. 0) THEN
               KCHANCE(NODES) = KCHANCE(I)
               NCHANCE(NODES) = NCHANCE(I)
            ENDIF
            IF (LICC    .GT. 0) THEN
               KICC(NODES) = KICC(I)
               NICC(NODES) = NICC(I)
            ENDIF
            IF (ALLOCATED(KPLQ)) KPLQ(NODES) = KPLQ(I)
!
            MNROWS = NODEINFO(NODES)%NROW + NODEINFO(NODES)%KROW
            IF (MNROWS .GT. IZROWS) THEN
               IF (.NOT. RESIZE .OR. (IZROWS .GE. MXROWS)      &
     &                          .OR. (MNROWS .GT. MXROWS)) GOTO 9200
               IZROWS = MNROWS
            ENDIF
!
            MNCOLS = NODEINFO(NODES)%NCOL + NODEINFO(NODES)%KCOL
            IF (MNCOLS .GT. IZCOLS) THEN
               IF (.NOT. RESIZE .OR. (IZCOLS .GE. MXCOLS)      &
     &                          .OR. (MNCOLS .GT. MXCOLS)) GOTO 9200
               IZCOLS = MNCOLS
            ENDIF
!
!     Copy the cost vector
!
            NCC = NODEINFO(NODES)%NCOL - NODEINFO(NODES)%NROW
            IF (LASTC + NCC .GT. IZCOST) THEN
               IF (.NOT. EXTMEM_COST(LASTC+NCC,IERR)) GOTO 9200
            ENDIF
!
            NODEINFO(NODES)%KCOST = LASTC
            DO K=1,NCC
               COST(NODEINFO(NODES)%KCOST+K) = COST(NODEINFO(I)%KCOST+K)
            END DO
            COPIED_SCALARS(I)%COST = .TRUE.
            LASTC = LASTC + NCC
!
!     Copy the RHS vector
!
            NR = NODEINFO(NODES)%NROW
            IF (LASTR + NR .GT. IZDRHS) THEN
               IF (.NOT. EXTMEM_DRHS(LASTR+NR,IERR)) GOTO 9200
            ENDIF
            NODEINFO(NODES)%KRHS = LASTR
            DO K=1,NR
               RHS(NODEINFO(NODES)%KRHS+K) = RHS(NODEINFO(I)%KRHS+K)
            END DO
            COPIED_SCALARS(I)%RHS = .TRUE.
            LASTR = LASTR + NR
!
!     Copy the upper and lower bounds
!
            NB = NODEINFO(NODES)%NCOL + 1
            IF (LASTBD + NB .GT. IZBNDS) THEN
               IF (.NOT. EXTMEM_BNDS(LASTBD+NB,IERR)) GOTO 9200
            ENDIF
            NODEINFO(NODES)%KBOUND = LASTBD
            DO K=1,NB
               XLB(NODEINFO(NODES)%KBOUND+K) = XLB(NODEINFO(I)%KBOUND+K)
               XUB(NODEINFO(NODES)%KBOUND+K) = XUB(NODEINFO(I)%KBOUND+K)
            END DO
            COPIED_SCALARS(I)%BOUNDS = .TRUE.
            LASTBD = LASTBD + NB
!
!     Copy variable types
!
            NV = NODEINFO(NODES)%NCOL
            NODEINFO(NODES)%VARPTR = LASTVP
            IF (LASTVP + NV .GT. IZVTYP) THEN
               IF (.NOT. EXTMEM_VTYP(LASTVP + NV,IERR)) GOTO 9200
            ENDIF
            DO K=1,NV
               VARTYP(LASTVP+K) = VARTYP(NODEINFO(I)%VARPTR+K)
            END DO
            COPIED_SCALARS(I)%VAR = .TRUE.
            LASTVP = LASTVP + NV
!
!     Copy the piecewise linear-quadratic penalty term
!
            IF (ALLOCATED(KPLQ)) THEN
               NCC = NODEINFO(NODES)%NCOL - NODEINFO(NODES)%NROW
               IF (LASTPLQ + NCC .GT. IZCPLQ) THEN
                  IF (.NOT. EXTMEM_CPLQ(LASTPLQ+NCC,IERR)) GOTO 9200
               ENDIF
!
               KPLQ(NODES) = LASTPLQ
               DO K=1,NCC
                  CPLQ(KPLQ(NODES)+K) = CPLQ(KPLQ(I)+K)
               END DO
               COPIED_SCALARS(I)%PLQ = .TRUE.
               LASTPLQ = LASTPLQ + NCC
            ENDIF
!
!     Copy the chance constraints
!
            IF (LCHANCE .GT. 0) THEN
               NCC = NCHANCE(NODES)
               IF (LCHANCE + NCC .GT. IZCCON) THEN
                  IF (.NOT. EXTMEM_CCON(LCHANCE+NCC,IERR)) GOTO 9200
               ENDIF
!
               KCHANCE(NODES) = LCHANCE
               DO K=1,NCC
                  XCHANCE(KCHANCE(NODES)+K) = XCHANCE(KCHANCE(I)+K)
               END DO
               COPIED_SCALARS(I)%CHANCE = .TRUE.
               LCHANCE = LCHANCE + NCC
            ENDIF
!
!     Copy the integrated chance constraints
!
            IF (LICC .GT. 0) THEN
               NCC = NICC(NODES)
               IF (LICC + NCC .GT. IZNICC) THEN
                  IF (.NOT. EXTMEM_NICC(LICC+NCC,IERR)) GOTO 9200
               ENDIF
!
               KICC(NODES) = LICC
               DO K=1,NCC
                  XICC(KICC(NODES)+K) = XICC(KICC(I)+K)
               END DO
               COPIED_SCALARS(I)%ICC = .TRUE.
               LICC = LICC + NCC
            ENDIF
!
!     Copy the A-matrix coefficients
!
            IF (MARKOV .AND. I .GT. 1) THEN
               NMTX = 2
            ELSE
               NMTX = I
            ENDIF
            NSUB = NMTX
            IF (HAVE_GC) NMTX = NMTX + I - 1
!
            MNABLK = LASTBA + NMTX
            IF (MNABLK .GT. IZABLK) THEN
               IF (.NOT. EXTMEM_ABLK(MNABLK,IERR)) GOTO 9200
            ENDIF
!
            DO J=1,NMTX
               KCOLA(KDATA(NODES)+J) = KCOLA(KDATA(I)+J)
               KELMA(KDATA(NODES)+J) = KELMA(KDATA(I)+J)
               NELMA(KDATA(NODES)+J) = NELMA(KDATA(I)+J)
            END DO
            LASTBA = LASTBA + NMTX
!
            DO JMTX=1,NMTX
               KAREF = KELMA(KDATA(I)+JMTX)
               LMN = NELMA(KDATA(NODES)+JMTX)
               IF (LASTA+LMN .GT. IZALMN) THEN
                  IF (.NOT. EXTMEM_ALMN(LASTA+LMN,IERR)) GOTO 9200
               ENDIF
               KELMA(KDATA(NODES)+JMTX) = LASTA
               DO K=1,LMN
                   A(LASTA+K) =  A(KAREF+K)
                  IA(LASTA+K) = IA(KAREF+K)
               END DO
               IF (JMTX .LE. NSUB) THEN
                  COPIED_MATRICES(I+1-JMTX, I)%A_MATRIX = .TRUE.
               ELSE
                  COPIED_MATRICES(I,JMTX-NSUB)%A_MATRIX = .TRUE.
               ENDIF
               LASTA = LASTA + LMN
!
               KCREF = KCOLA(KDATA(I)+JMTX)
               IF (JMTX .LE. NSUB) THEN
                  J = I + 1 - JMTX
                  NCA = NODEINFO(J)%NCOL - NODEINFO(J)%NROW + 1
               ELSE
                  NCA = NODEINFO(I)%NCOL - NODEINFO(I)%NROW + 1
               ENDIF
               IF (LASTCA+NCA .GT. IZACOL) THEN
                  IF (.NOT. EXTMEM_ACOL(LASTCA+NCA,IERR)) GOTO 9200
               ENDIF
               KCOLA(KDATA(NODES)+JMTX) = LASTCA
               DO K=1,NCA
                  LA(LASTCA+K) =  LA(KCREF+K)
               END DO
               LASTCA = LASTCA + NCA
            END DO
!
            IF (NQPROF .GT. 0) THEN
               NQMTX = MIN(NQPROF,I)
               MNQBLK = LASTQB + NQMTX
               IF (MNQBLK .GT. IZQBLK) THEN
                  IF (.NOT. EXTMEM_QBLK(MNQBLK,IERR)) GOTO 9200
               ENDIF
!
               DO J=1,NQMTX
                  IQOFF(KDATQ(NODES)+J) = IQOFF(KDATQ(I)+J)
                  LQOFF(KDATQ(NODES)+J) = LQOFF(KDATQ(I)+J)
                  NELMQ(KDATQ(NODES)+J) = NELMQ(KDATQ(I)+J)
               END DO
               LASTQB = LASTQB + NQMTX
!
               DO JMTX=1,NQMTX
                  KQREF = IQOFF(KDATQ(I)+JMTX)
                  LMN = NELMQ(KDATQ(NODES)+JMTX)
                  IF (LASTQ + LMN .GT. IZQLMN) THEN
                     IF (.NOT. EXTMEM_QLMN(LASTQ+LMN,IERR)) GOTO 9200
                  ENDIF
                  IQOFF(KDATQ(NODES)+JMTX) = LASTQ
                  DO K=1,LMN
                     AQMTX(LASTQ+K) = AQMTX(KQREF+K)
                     IQMTX(LASTQ+K) = IQMTX(KQREF+K)
                  END DO
                  COPIED_MATRICES(I+1-JMTX,I)%Q_MATRIX = .TRUE.
                  LASTQ = LASTQ + LMN
!
                  KCREF = LQOFF(KDATQ(I)+JMTX)
                  J = I + 1 - JMTX
                  NCA = NODEINFO(J)%NCOL - NODEINFO(J)%NROW + 1
                  IF (LASTQC+NCA .GT. IZQCOL) THEN
                     IF (.NOT. EXTMEM_QCOL(LASTQC+NCA,IERR)) GOTO 9200
                  ENDIF
                  LQOFF(KDATQ(NODES)+JMTX) = LASTQC
                  DO K=1,NCA
                     LQMTX(LASTQC+K) =  LQMTX(KCREF+K)
                  END DO
                  LASTQC = LASTQC + NCA
               END DO
            ENDIF
!
!     Here we have a scenario other than the root scenario
!
         ELSE
            IF (I .EQ. IPER1) THEN
               NODEINFO(NODE0)%NDESC  = NODEINFO(NODE0)%NDESC + 1
               NODEINFO(NODES)%IANCTR = NODE0
               IF (NODEINFO(NODE0)%IDESC .GT. 0) THEN
                  IBRO1 = NODEINFO(NODE0)%IDESC
  330             CONTINUE
                     IF (NODEINFO(IBRO1)%IBROTH .GT. 0) THEN
                        IBRO1 = NODEINFO(IBRO1)%IBROTH
                        GOTO 330
                     ENDIF
!
                  NODEINFO(NODES)%IBROTH = NODEINFO(IBRO1)%IBROTH
                  NODEINFO(IBRO1)%IBROTH = NODES
                  NODEINFO(NODES)%IDESC  = 0
                  NODEINFO(NODES)%NDESC  = 0
                  NODEINFO(NODES)%PROB   = CPROB(I)
               ELSE
                  NODEINFO(NODE0)%IDESC = NODES
                  NODEINFO(NODE0)%NDESC = 1
                  IBRO1 = IRNGE0(IPER1)
  340             CONTINUE
                     IF (NODEINFO(IBRO1)%IBROTH .NE. 0) THEN
                        IBRO1 = ABS(NODEINFO(IBRO1)%IBROTH)
                        GOTO 340
                     ENDIF
!
                  NODEINFO(NODES)%IBROTH = 0
                  NODEINFO(IBRO1)%IBROTH = -NODES
                  NODEINFO(NODES)%IDESC  = 0
                  NODEINFO(NODES)%NDESC  = 0
                  NODEINFO(NODES)%PROB   = CPROB(I)
               ENDIF
            ELSE
               NODEINFO(NODES)%IANCTR  = NODES - 1
               NODEINFO(NODES)%IDESC   = 0
               NODEINFO(NODES)%NDESC   = 0
               NODEINFO(NODES-1)%IDESC = NODES
               NODEINFO(NODES-1)%NDESC = 1
               NODEINFO(NODES)%PROB    = CPROB(I)
               IBRO1 = IRNGE0(I)
  350          CONTINUE
                  IF (NODEINFO(IBRO1)%IBROTH .NE. 0) THEN
                     IBRO1 = ABS(NODEINFO(IBRO1)%IBROTH)
                     GOTO 350
                  ENDIF
!
                  NODEINFO(NODES)%IBROTH = 0
                  NODEINFO(IBRO1)%IBROTH = -NODES
                  NODEINFO(NODES)%IDESC  = 0
                  NODEINFO(NODES)%NDESC  = 0
            ENDIF
!
!     Initialize pointers
!
            KDATA(NODES) = LASTBA
            NODEINFO(NODES)%KROW   = NODEINFO(NODES-1)%KROW +      &
     &                               NODEINFO(NODES-1)%NROW
            NODEINFO(NODES)%KCOL   = NODEINFO(NODES-1)%KCOL +      &
     &                               NODEINFO(NODES-1)%NCOL + 1
            NODEINFO(NODES)%KCOST  = NODEINFO(I0)%KCOST
            NODEINFO(NODES)%KRHS   = NODEINFO(I0)%KRHS
            NODEINFO(NODES)%KBOUND = NODEINFO(I0)%KBOUND
            NODEINFO(NODES)%KNAMES = NODEINFO(I0)%KNAMES
            NODEINFO(NODES)%VARPTR = NODEINFO(I0)%VARPTR
            NODEINFO(NODES)%NROW   = NODEINFO(I0)%NROW
            NODEINFO(NODES)%NCOL   = NODEINFO(I0)%NCOL
            IF (ALLOCATED(KDATQ)) KDATQ(NODES) = LASTQB
            IF (LCHANCE .GT. 0) THEN
               KCHANCE(NODES) = KCHANCE(I0)
               NCHANCE(NODES) = NCHANCE(I0)
            ENDIF
            IF (LICC    .GT. 0) THEN
               KICC(NODES) = KICC(I0)
               NICC(NODES) = NICC(I0)
            ENDIF
            IF (ALLOCATED(KPLQ)) KPLQ(NODES) = KPLQ(I0)
!
            MNROWS = NODEINFO(NODES)%NROW + NODEINFO(NODES)%KROW
            IF (MNROWS .GT. IZROWS) THEN
               IF (.NOT. RESIZE .OR. (IZROWS .GE. MXROWS)      &
     &                          .OR. (MNROWS .GT. MXROWS)) GOTO 9200
               IZROWS = MNROWS
            ENDIF
!
            MNCOLS = NODEINFO(NODES)%NCOL + NODEINFO(NODES)%KCOL
            IF (MNCOLS .GT. IZCOLS) THEN
               IF (.NOT. RESIZE .OR. (IZCOLS .GE. MXCOLS)      &
     &                          .OR. (MNCOLS .GT. MXCOLS)) GOTO 9200
               IZCOLS = MNCOLS
            ENDIF
!
!     Initial the matrix pointers (which will be needed in all cases)
!
            IF (MARKOV .AND. I .GT. 1) THEN
               NMTX = 2
            ELSE
               NMTX = I
            ENDIF
            NSUB = NMTX
            IF (HAVE_GC) NMTX = NMTX + I - 1
!
            MNABLK = LASTBA + NMTX
            IF (MNABLK .GT. IZABLK) THEN
               IF (.NOT. EXTMEM_ABLK(MNABLK,IERR)) GOTO 9200
            ENDIF
!
            DO J=1,NMTX
               KCOLA(KDATA(NODES)+J) = KCOLA(KDATA(I0)+J)
               KELMA(KDATA(NODES)+J) = KELMA(KDATA(I0)+J)
               NELMA(KDATA(NODES)+J) = NELMA(KDATA(I0)+J)
            END DO
            LASTBA = LASTBA + NMTX
!
            IF (NQPROF .GT. 0) THEN
               NQMTX = MIN(NQPROF,I)
               MNQBLK = LASTQB + NQMTX
               IF (MNQBLK .GT. IZQBLK) THEN
                  IF (.NOT. EXTMEM_QBLK(MNQBLK,IERR)) GOTO 9200
               ENDIF
!
               DO J=1,NQMTX
                  IQOFF(KDATQ(NODES)+J) = IQOFF(KDATQ(I0)+J)
                  LQOFF(KDATQ(NODES)+J) = LQOFF(KDATQ(I0)+J)
                  NELMQ(KDATQ(NODES)+J) = NELMQ(KDATQ(I0)+J)
               END DO
               LASTQB = LASTQB + NQMTX
            ENDIF
         ENDIF
!
!     Store the result for every stochastic element belonging to this period
!
         DO J=1,NSTELM
            JTYPE = KSTOCH(4*J-3)
            JBLOK = KSTOCH(4*J-2)
            JLOCN = KSTOCH(4*J-1)
            JMODE = KSTOCH(4*J)
!
!     Random cost coefficient
!
            IF (JTYPE .EQ. 2) THEN
               IF (JBLOK .EQ. I) THEN
                  REFVAL = COST(NODEINFO(I0)%KCOST+JLOCN)
                  IF (JMODE .EQ. 1) THEN
                     TARGET = WORK(J)
                  ELSEIF (JMODE .EQ. 2) THEN
                     TARGET = WORK(J) + COST(NODEINFO(I)%KCOST+JLOCN)
                  ELSE
                     TARGET = WORK(J) * COST(NODEINFO(I)%KCOST+JLOCN)
                  ENDIF
                  IF (ABS(REFVAL - TARGET) .GT. ZTOLIN) THEN
                     IF (.NOT. COPIED_SCALARS(JBLOK)%COST) THEN
                        NCC = NODEINFO(NODES)%NCOL                         &
     &                      - NODEINFO(NODES)%NROW
                        IF (LASTC + NCC .GT. IZCOST) THEN
                           IF (.NOT. EXTMEM_COST(LASTC+NCC,IERR))          &
     &                                                      GOTO 9200
                        ENDIF
!
                        NODEINFO(NODES)%KCOST = LASTC
                        DO K=1,NCC
                           COST(LASTC+K) = COST(NODEINFO(I0)%KCOST+K)
                        END DO
                        COPIED_SCALARS(I)%COST = .TRUE.
                        LASTC = LASTC + NCC
                     ENDIF
                     KC = NODEINFO(NODES)%KCOST
                     COST(KC+JLOCN) = TARGET
                  ENDIF
               ENDIF
!
!     Random right-hand side
!
            ELSEIF (JTYPE .EQ.  3) THEN
               IF (JBLOK .EQ. I) THEN
                  REFVAL = RHS(NODEINFO(I0)%KRHS+JLOCN)
                  IF (JMODE .EQ. 1) THEN
                     TARGET = WORK(J)
                  ELSEIF (JMODE .EQ. 2) THEN
                     TARGET = WORK(J) + RHS(NODEINFO(I)%KRHS+JLOCN)
                  ELSE
                     TARGET = WORK(J) * RHS(NODEINFO(I)%KRHS+JLOCN)
                  ENDIF
                  IF (ABS(REFVAL - TARGET) .GT. ZTOLIN) THEN
                     IF (.NOT. COPIED_SCALARS(JBLOK)%RHS) THEN
                        NR = NODEINFO(NODES)%NROW
                        IF (LASTR + NR .GT. IZDRHS) THEN
                           IF (.NOT. EXTMEM_DRHS(LASTR+NR,IERR))      &
     &                                                     GOTO 9200
                        ENDIF
                        NODEINFO(NODES)%KRHS = LASTR
                        DO K=1,NR
                           RHS(LASTR+K) = RHS(NODEINFO(I0)%KRHS+K)
                        END DO
                        COPIED_SCALARS(I)%RHS = .TRUE.
                        LASTR = LASTR + NR
                     ENDIF
                     KR = NODEINFO(NODES)%KRHS
                     RHS(KR+JLOCN) = TARGET
                  ENDIF
               ENDIF
!
!     Random bound (upper or lower or both)
!
            ELSEIF (JTYPE .EQ. 5 .OR. JTYPE .EQ. 6                       &
     &                           .OR. JTYPE .EQ. 11) THEN
               IF (JBLOK .EQ. I) THEN
                  IF (JTYPE .EQ. 5) THEN
                     REFVAL1 = XLB(NODEINFO(I0)%KBOUND+JLOCN)
                     IF (JMODE .EQ. 1) THEN
                        TARGET1 = WORK(J)
                     ELSEIF (JMODE .EQ. 2) THEN
                        TARGET1 = WORK(J)+XLB(NODEINFO(I)%KBOUND+JLOCN)
                     ELSE
                        TARGET1 = WORK(J)*XLB(NODEINFO(I)%KBOUND+JLOCN)
                     ENDIF
                     TARGET2 = 0.D0
                     REFVAL2 = 0.D0
                  ELSEIF (JTYPE .EQ. 6) THEN
                     REFVAL2 = XUB(NODEINFO(I0)%KBOUND+JLOCN)
                     IF (JMODE .EQ. 1) THEN
                        TARGET2 = WORK(J)
                     ELSEIF (JMODE .EQ. 2) THEN
                        TARGET2 = WORK(J)+XUB(NODEINFO(I)%KBOUND+JLOCN)
                     ELSE
                        TARGET2 = WORK(J)*XUB(NODEINFO(I)%KBOUND+JLOCN)
                     ENDIF
                     TARGET1 = 0.D0
                     REFVAL1 = 0.D0
                  ELSE
                     REFVAL1 = XLB(NODEINFO(I0)%KBOUND+JLOCN)
                     REFVAL2 = XUB(NODEINFO(I0)%KBOUND+JLOCN)
                     IF (JMODE .EQ. 1) THEN
                        TARGET1 = WORK(J)
                        TARGET2 = WORK(J)
                     ELSEIF (JMODE .EQ. 2) THEN
                        TARGET1 = WORK(J)+XLB(NODEINFO(I)%KBOUND+JLOCN)
                        TARGET2 = WORK(J)+XUB(NODEINFO(I)%KBOUND+JLOCN)
                     ELSE
                        TARGET1 = WORK(J)*XLB(NODEINFO(I)%KBOUND+JLOCN)
                        TARGET2 = WORK(J)*XUB(NODEINFO(I)%KBOUND+JLOCN)
                     ENDIF
                  ENDIF
                  IF (ABS(REFVAL1 - TARGET1) .GT. ZTOLIN .OR.          &
     &                ABS(REFVAL2 - TARGET2) .GT. ZTOLIN) THEN
                     IF (.NOT. COPIED_SCALARS(JBLOK)%BOUNDS) THEN
                        NB = NODEINFO(NODES)%NCOL + 1
                        IF (LASTBD + NB .GT. IZBNDS) THEN
                           IF (.NOT. EXTMEM_BNDS(LASTBD+NB,IERR))      &
     &                                                      GOTO 9200
                        ENDIF
                        NODEINFO(NODES)%KBOUND = LASTBD
                        DO K=1,NB
                           XLB(LASTBD+K) = XLB(NODEINFO(I0)%KBOUND+K)
                           XUB(LASTBD+K) = XUB(NODEINFO(I0)%KBOUND+K)
                        END DO
                        COPIED_SCALARS(I)%BOUNDS = .TRUE.
                        LASTBD = LASTBD + NB
                     ENDIF
!
                     KB = NODEINFO(NODES)%KBOUND
                     IF (JTYPE .EQ. 5) THEN
                        XLB(KB+JLOCN) = TARGET1
                     ELSEIF (JTYPE .EQ. 6) THEN
                        XUB(KB+JLOCN) = TARGET2
                     ELSE
                        XLB(KB+JLOCN) = TARGET1
                        XUB(KB+JLOCN) = TARGET2
                     ENDIF
                  ENDIF
               ENDIF
!
!     Random demand or supply
!
            ELSEIF (JTYPE .EQ. 12 .OR. JTYPE .EQ. 13) THEN
               IF (JBLOK .EQ. I) THEN
                  REFVAL = RHS(NODEINFO(I0)%KRHS+JLOCN)
                  IF (JMODE .EQ. 1) THEN
                     TARGET = WORK(J)
                  ELSEIF (JMODE .EQ. 2) THEN
                     TARGET = WORK(J) + RHS(NODEINFO(I)%KRHS+JLOCN)
                  ELSE
                     TARGET = WORK(J) * RHS(NODEINFO(I)%KRHS+JLOCN)
                  ENDIF
                  IF (ABS(REFVAL - TARGET) .GT. ZTOLIN) THEN
                     IF (.NOT. COPIED_SCALARS(JBLOK)%RHS) THEN
                        NR = NODEINFO(NODES)%NROW
                        IF (LASTR + NR .GT. IZDRHS) THEN
                           IF (.NOT. EXTMEM_DRHS(LASTR+NR,IERR))      &
     &                                                     GOTO 9200
                        ENDIF
                        NODEINFO(NODES)%KRHS = LASTR
                        DO K=1,NR
                           RHS(LASTR+K) = RHS(NODEINFO(I0)%KRHS+K)
                        END DO
                        COPIED_SCALARS(I)%RHS = .TRUE.
                        LASTR = LASTR + NR
                     ENDIF
!
                     KR = NODEINFO(NODES)%KRHS
                     IF (JTYPE .EQ. 12) THEN
                        RHS(KR+JLOCN) =  TARGET
                     ELSE
                        RHS(KR+JLOCN) = -TARGET
                     ENDIF
                  ENDIF
               ENDIF
!
!     Random probability limit in a chance constraint
!
            ELSEIF (JTYPE .EQ.  8) THEN
               IF (JBLOK .EQ. I) THEN
                  REFVAL = XCHANCE(KCHANCE(I0)+JLOCN)%PROB
                  IF (JMODE .EQ. 1) THEN
                     TARGET = WORK(J)
                  ELSEIF (JMODE .EQ. 2) THEN
                     TARGET = WORK(J) + XCHANCE(KCHANCE(I)+JLOCN)%PROB
                  ELSE
                     TARGET = WORK(J) * XCHANCE(KCHANCE(I)+JLOCN)%PROB
                  ENDIF
                  IF (ABS(REFVAL - TARGET) .GT. ZTOLIN) THEN
                     IF (.NOT. COPIED_SCALARS(JBLOK)%CHANCE) THEN
                        NC = NCHANCE(NODES)
                        MNCCON = LCHANCE + NC
                        IF (MNCCON .GT. IZCCON) THEN
                           IF (.NOT. EXTMEM_CCON(MNCCON,IERR)) GOTO 9200
                        ENDIF
                        KCHANCE(NODES) = LCHANCE
                        DO K=1,NC
                           XCHANCE(LCHANCE+K) = XCHANCE(KCHANCE(I0)+K)
                        END DO
                        COPIED_SCALARS(I)%CHANCE = .TRUE.
                        LCHANCE = LCHANCE + NC
                     ENDIF
                     KC = KCHANCE(NODES) + JLOCN
                     XCHANCE(KC)%PROB = TARGET
                  ENDIF
               ENDIF
!
!     Random CVaR limit in an integrated chance constraint
!
            ELSEIF (JTYPE .EQ.  9) THEN
               IF (JBLOK .EQ. I) THEN
                  REFVAL = XICC(KICC(I0)+JLOCN)%LIMIT
                  IF (JMODE .EQ. 1) THEN
                     TARGET = WORK(J)
                  ELSEIF (JMODE .EQ. 2) THEN
                     TARGET = WORK(J) + XICC(KICC(I)+JLOCN)%LIMIT
                  ELSE
                     TARGET = WORK(J) * XICC(KICC(I)+JLOCN)%LIMIT
                  ENDIF
                  IF (ABS(REFVAL - TARGET) .GT. ZTOLIN) THEN
                     IF (.NOT. COPIED_SCALARS(JBLOK)%ICC) THEN
                        NC = NICC(NODES)
                        IF (LICC + NC .GT. IZNICC) THEN
                           IF (.NOT. EXTMEM_NICC(LICC+NC,IERR))      &
     &                                                    GOTO 9200
                        ENDIF
                        KICC(NODES) = LICC
                        DO K=1,NC
                           XICC(LICC+K) = XICC(KICC(I0)+K)
                        END DO
                        COPIED_SCALARS(I)%ICC = .TRUE.
                        LICC = LICC + NC
                     ENDIF
                     KC = KICC(NODES) + JLOCN
                     XICC(KC)%LIMIT = TARGET
                  ENDIF
               ENDIF
!
!     Random curvature in a piecewise linear-quadratic penalty
!
            ELSEIF (JTYPE .EQ. 10) THEN
               IF (JBLOK .EQ. I) THEN
                  REFVAL = CPLQ(KPLQ(I0)+JLOCN)
                  IF (JMODE .EQ. 1) THEN
                     TARGET = WORK(J)
                  ELSEIF (JMODE .EQ. 2) THEN
                     TARGET = WORK(J) + CPLQ(KPLQ(I)+JLOCN)
                  ELSE
                     TARGET = WORK(J) * CPLQ(KPLQ(I)+JLOCN)
                  ENDIF
                  IF (ABS(REFVAL - TARGET) .GT. ZTOLIN) THEN
                     IF (.NOT. COPIED_SCALARS(JBLOK)%PLQ) THEN
                        NCC = NODEINFO(NODES)%NCOL      &
     &                      - NODEINFO(NODES)%NROW
                        MNPLQ = LASTPLQ + NCC
                        IF (MNPLQ .GT. IZCPLQ) THEN
                           IF (.NOT. EXTMEM_CPLQ(MNPLQ,IERR)) GOTO 9200
                        ENDIF
!
                        KPLQ(NODES) = LASTPLQ
                        DO K=1,NCC
                           COST(LASTPLQ+K) = CPLQ(KPLQ(I0)+K)
                        END DO
                        COPIED_SCALARS(I)%PLQ = .TRUE.
                        LASTPLQ = LASTPLQ + NCC
                     ENDIF
                     KC = KPLQ(NODES)
                     CPLQ(KC+JLOCN) = TARGET
                  ENDIF
               ENDIF
!
!     Random coefficient in the A-matrix
!
            ELSEIF (JTYPE .EQ.  1) THEN
               IF (JBLOK .GT. NBLK(I-1,1) .AND.       &
     &             JBLOK .LE. NBLK(I,  1)) THEN
                  JMTX = JBLOK - NBLK(I-1,1)
                  KAREF = KELMA(KDATA(I0)+JMTX)
                  REFVAL = A(KAREF+JLOCN)
                  IF (JMODE .EQ. 1) THEN
                     TARGET = WORK(J)
                  ELSEIF (JMODE .EQ. 2) THEN
                     TARGET = WORK(J) + A(KELMA(KDATA(I)+JMTX)+JLOCN)
                  ELSE
                     TARGET = WORK(J) * A(KELMA(KDATA(I)+JMTX)+JLOCN)
                  ENDIF
                  IF (ABS(REFVAL - TARGET) .GT. ZTOLIN) THEN
                     IF (MARKOV) THEN
                        NSUB = MIN(I,2)
                     ELSE
                        NSUB = I
                     ENDIF
                     IF (JMTX .LE. NSUB) THEN
                        COPIED = COPIED_MATRICES(I+1-JMTX, I)%A_MATRIX
                        COPIED_MATRICES(I+1-JMTX, I)%A_MATRIX = .TRUE.
                        STOCHA(I+1-JMTX, I) = .TRUE.
                     ELSE
                        COPIED = COPIED_MATRICES(I,JMTX-NSUB)%A_MATRIX
                        COPIED_MATRICES(I,JMTX-NSUB)%A_MATRIX = .TRUE.
                        STOCHA(I,JMTX-NSUB) = .TRUE.
                     ENDIF
                     IF (.NOT. COPIED) THEN
                        LMN = NELMA(KDATA(NODES)+JMTX)
                        MNALMN = LASTA + LMN
                        IF (MNALMN .GT. IZALMN) THEN
                           IF (.NOT. EXTMEM_ALMN(MNALMN,IERR)) GOTO 9200
                        ENDIF
                        KELMA(KDATA(NODES)+JMTX) = LASTA
                        DO K=1,LMN
                            A(LASTA+K) =  A(KAREF+K)
                           IA(LASTA+K) = IA(KAREF+K)
                        END DO
                        LASTA = LASTA + LMN
                     ENDIF
                     A(KELMA(KDATA(NODES)+JMTX)+JLOCN) = TARGET
                  ENDIF
               ENDIF
!
!     Random coefficient in the Q-matrix
!
            ELSEIF (JTYPE .EQ.  7) THEN
               IF (JBLOK .GT. NBLK(I-1,2) .AND.       &
     &             JBLOK .LE. NBLK(I,  2)) THEN
                  JMTX = JBLOK - NBLK(I-1,2)
                  KQREF = IQOFF(KDATQ(I0)+JMTX)
                  REFVAL = AQMTX(KQREF+JLOCN)
                  IF (JMODE .EQ. 1) THEN
                     TARGET = WORK(J)
                  ELSEIF (JMODE .EQ. 2) THEN
                     TARGET = WORK(J)+AQMTX(IQOFF(KDATQ(I)+JMTX)+JLOCN)
                  ELSE
                     TARGET = WORK(J)*AQMTX(IQOFF(KDATQ(I)+JMTX)+JLOCN)
                  ENDIF
                  IF (ABS(REFVAL - TARGET) .GT. ZTOLIN) THEN
                     COPIED = COPIED_MATRICES(I+1-JMTX, I)%Q_MATRIX
                     COPIED_MATRICES(I+1-JMTX, I)%Q_MATRIX = .TRUE.
                     IF (.NOT. COPIED) THEN
                        LMN = NELMQ(KDATQ(NODES)+JMTX)
                        MNQLMN = LASTQ + LMN
                        IF (MNQLMN .GT. IZQLMN) THEN
                           IF (.NOT. EXTMEM_QLMN(MNQLMN,IERR)) GOTO 9200
                        ENDIF
                        IQOFF(KDATQ(NODES)+JMTX) = LASTQ
                        DO K=1,LMN
                           AQMTX(LASTQ+K) = AQMTX(KQREF+K)
                           IQMTX(LASTQ+K) = IQMTX(KQREF+K)
                        END DO
                        LASTQ = LASTQ + LMN
                     ENDIF
                     AQMTX(IQOFF(KDATQ(NODES)+JMTX)+JLOCN) = TARGET
                  ENDIF
               ENDIF
            ENDIF
         END DO
      END DO
      RETURN
!
 9200 CONTINUE
      WRITE (IOLOG, 1000)
      IERR = 1
      RETURN
!
 1000 FORMAT(' XXX -  FATAL  - Memory error in routine SNGLBRCH.')
! 1100 FORMAT(' XXX -  FATAL  - Too many parameters called by location.',
!     * /,' Increase parameter MXLOC (routine CALLDIST) to at least',I7)
! 3100 FORMAT(' XXX -  FATAL  - More than',I6,' cost coefficients.',
!     *       ' Increase parameter MXCOST.')
! 3200 FORMAT(' XXX -  FATAL  - More than',I6,' RHS coefficients.',
!     *       ' Increase parameter MXDRHS.')
! 3300 FORMAT(' XXX -  FATAL  - More than',I6,' bound coefficients.',
!     *       ' Increase parameter MXBNDS.')
! 3400 FORMAT(' XXX -  FATAL  - More than',I6,' A-matrix coefficients.',
!     *       ' Increase parameter MXALMN.')
! 3410 FORMAT(' XXX -  FATAL  - More than',I6,' A-matrix blocks.',
!     *       ' Increase parameter MXABLK.')
! 3450 FORMAT(' XXX -  FATAL  - More than',I6,' A-matrix columns.',
!     *       ' Increase parameter MXACOL.')
! 3500 FORMAT(' XXX -  FATAL  - More than',I6,' Q-matrix coefficients.',
!     *       ' Increase parameter MXQLMN.')
! 3510 FORMAT(' XXX -  FATAL  - More than',I6,' Q-matrix blocks.',
!     *       ' Increase parameter MXQBLK.')
! 3550 FORMAT(' XXX -  FATAL  - More than',I6,' Q-matrix columns.',
!     *       ' Increase parameter MXQBLK.')
! 3600 FORMAT(' XXX -  FATAL  - More than',I6,' variables.',
!     *       ' Increase parameter MXVTYP.')
! 3700 FORMAT(' XXX -  FATAL  - More than',I6,' variable names.',
!     *       ' Increase parameter MXVNAM.')
      END SUBROUTINE
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE SAMPLE(NSAMPL,WEIGHT,IERR)
!
!     This subroutine samples NSAMPL scenarios, one at a time. The scenarios
!     are placed into an increasing event tree. For each scenario we have to
!     find a node in the existing tree from which the new scenario branches.
!     There are several ways to select the branching node. One could use any
!     existing non-terminal node with equal probability, or one could use a
!     more complicated scheme. Randomly selecting a node may introduce bias
!     because the longer the algorithm has been run, the more scenarios there
!     will be in later stages, and the larger the probability of introducing
!     a branch late in the tree.
!
!     We opt for a two-step process: First the stage is generated randomly,
!     and then a node is selected within the stage. If we select the node
!     randomly, there is the same problem with biasing, since nodes that were
!     generated early in the process are more likely to be selected as branching
!     nodes than others generated later. In order to obtain a balanced tree,
!     we only use ONE NODE per scenario to branch from. Once another node in
!     this period has been generated, it becomes the new branching node.
!
!     The problem with this approach is that the expected number of branches
!     per node depends not only on the number of samples generated, but also
!     on the probability with which each stage is selected. The larger the
!     sample, the smaller we should set the probability of selecting a stage
!     early in the process. If the stage probabilities are set carefully
!     (they will depend on the sample size), then a sufficiently large sample
!     can approximate the original data process arbitrarily closely.
!
!     We set the stage probabilities in a separate routine SETWT, since
!     a knowledgeable user may be able to calibrate the selection process
!     better than a general-purpose system.
!
!     ----------------------------------------------------------------------
!          This version dated 1 June 1998. Written by Gus Gassmann.
!     ----------------------------------------------------------------------
!
      USE SAMPUTIL, ONLY: RNDGET
      USE SMPS_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      implicit real*8 (a-h,o,p,r-z), integer*4(i-n), character*1(q)
!
      CHARACTER*(NAMLEN) DTEMP,DBLANK
      CHARACTER*40 F1500
      CHARACTER*1  QTEMP(NAMLEN),QBLANK(NAMLEN)
      REAL*8       WEIGHT(*),RVEC(2)
!
      EQUIVALENCE (QTEMP,DTEMP), (QBLANK,DBLANK)
!
      DATA QBLANK/NAMLEN*' '/
!
      IERR = 0
      NSCHAR = 0
!
      IF (NPER .EQ. 1) GOTO 200
!
!     Allocate temporary arrays
!
      ALLOCATE(INSTAGE(NPER), CPROB(NPER), QPROB(NPER), STAT=IER)
      IF (IER .GT. 0) GOTO 960
      ALLOCATE (COPIED_SCALARS(NPER), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
      ALLOCATE (COPIED_MATRICES(NPER,NPER), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
      ALLOCATE (WORK(NSTELM), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
      ALLOCATE (NBLK(0:NPER,2), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
      ALLOCATE (XRNDPAR(LLPAR), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
      ALLOCATE (KREF(NPER), STAT=IERR)
      IF (IERR .GT. 0) GOTO 960
!
      DO I=1,NPER
         NBLK(I-1,1) = KDATA(I)
         NBLK(I-1,2) = KDATQ(I)
      END DO
      NBLK(NPER,1) = LASTBA
      NBLK(NPER,2) = LASTQB
!
      NRCOMP = 0
      DO I=1,NRVAR
         NRCOMP = NRCOMP + MDIST(8*I-7)
      END DO
      IF (NRCOMP .EQ. 0) GOTO 999
!
      ALLOCATE (HIST(NRCOMP), STAT=IER)
      IF (IER .GT. 0) GOTO 960
      HIST(1:NRCOMP) = 0.D0
!
!     Set up breakpoints first (used to select period)
!
      QPROB(1) = WEIGHT(1)
      DO I=2,NPER-1
         QPROB(I) = WEIGHT(I) + QPROB(I-1)
      END DO
!
!     Now sample the scenarios one at a time.
!
      DO I=1,NSAMPL
!
!     The first sample must duplicate the entire base scenario
!
         IF (I .EQ. 1) THEN
            DO J=1,NPER-1
               INSTAGE(J) = NODES + J
            END DO
            CALL SAMPLE1(1,0,IERR)
            IF (IERR .GT. 0) RETURN
!
         ELSE
            CALL RNDGET(RVEC,1,IERR)
            DO J=1,NPER-1
               IF (RVEC(1) .LT. QPROB(J)) GOTO 150
            END DO
!
  150       CONTINUE
               IPER = J
!
!     Identify the branching node
!
               NODE0 = INSTAGE(IPER)
               DO J=IPER+1,NPER-1
                  INSTAGE(J) = NODES + J - IPER
               END DO
!
               CALL SAMPLE1(IPER+1,NODE0,IERR)
               IF (IERR .GT. 0) RETURN
         ENDIF
!
!     Give a name to the scenario
!
         IF (I .GE. IZSCEN) THEN
            IF (.NOT. EXTMEM_SCEN(I+1,IERR)) GOTO 960
         ENDIF
         LEAFND(I) = NODES
         IF (I .LT. 10000) THEN
            LEN = 8
            F1500 = '(''Smpl'',I4.4)'
         ELSE
            LEN = IFIX(ALOG10(REAL(I))) + 1
            IF (LEN .LT. 10) THEN
               WRITE (F1500, 1510) LEN
            ELSE
               WRITE (F1500, 1520) LEN
            ENDIF
            LEN = LEN + 4
         ENDIF
         WRITE (DTEMP, F1500) I
         IF (NSCHAR+LEN .GT. IZCHSC) THEN
            IF (.NOT. EXTMEM_CHSC(NSCHAR+LEN,IER)) GOTO 960
         ENDIF
         DO J=1,LEN
            QSCNAM(NSCHAR+J) = QTEMP(J)
         END DO
         NSCHAR = NSCHAR + LEN
         KSCNAM(I+1) = NSCHAR
      END DO
!
      MAXRHS  = LASTR
      TREEFMT = 1
!
!     Once the tree has been constructed, find the probabilities
!
      CALL SETPRB(IERR)
      RETURN
!
  200 CONTINUE
      WRITE (IOLOG, 1000)
      GOTO 999
!
  960 CONTINUE
      WRITE (IOLOG, 3200)
      IERR = 2
!
  999 CONTINUE
      DEALLOCATE (CPROB, QPROB, HIST, INSTAGE, KREF, COPIED_SCALARS,      &
     &            COPIED_MATRICES, XRNDPAR, NBLK, WORK, STAT=IER)
      IF (IER .NE. 0) THEN
         WRITE (IOLOG, 3200)
         IERR = 2
      ENDIF
      RETURN
!
 1000 FORMAT(' This feature still under construction...')
 1510 FORMAT('(''Smpl'',I',i1,')')
 1520 FORMAT('(''Smpl'',I',i2,')')
 3200 FORMAT(' XXX -  FATAL  - Memory handling error in SAMPLE.')
 3700 FORMAT(' XXX -  FATAL  - Too many scenario names.',      &
     &       ' Increase parameter MXCHSC.')
      END SUBROUTINE
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE SETWT(NPER,NSAMPL,WEIGHT)
!
!     This routine sets default weights for the sampling procedure. The idea
!     is to make a tree that is balanced, i.e., on average, the number of
!     branches should be the same for each nonterminal node. In order for this
!     to work, the selection probability for each stage is a multiple ALPHA
!     times that of the previous stage, where ALPHA = (NSAMPL)**(1/(NPER-1)).
!
      REAL*8 ALPHA, WEIGHT(*), TEMP, SUM
!
      ALPHA = REAL(NSAMPL)**(1.D0/(REAL(NPER-1)))
      SUM   = 0.D0
      TEMP  = 1.D0
      DO I=1,NPER-1
         WEIGHT(I) = TEMP
         SUM = SUM + TEMP
         TEMP = TEMP * ALPHA
      END DO
!
      DO I=1,NPER-1
         WEIGHT(I) = WEIGHT(I) / SUM
      END DO
      RETURN
      END SUBROUTINE
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE SAMPLE1(IPER1,NODE0,IERR)
!
!     This subroutine samples a single scenario, using information
!     read earlier from the input files and adds it to the event tree.
!     The new scenario branches from an existing scenario in NODE0 and starts
!     in IPER1. In other words, the first node on this scenario has NODE0
!     as its ancestor and is part of period IPER1.
!
!     Representation of random variables and stochastic elements is as in
!     routine INCONT, that is, using arrays MDIST and SDIST to describe the
!     distributions of the random variables, KSTOCH to describe the location
!     of the stochastic elements, and the D-matrix (DMTX,IDMTX and LDMTX)
!     to hold the two together.
!
!     Data structures:
!     ----------------
!
!     1. The distributions are stored as strings in array SDIST.
!
!     2. The dimension, number and locations of parameters and the stage
!        in which the random variable becomes known are stored in array MDIST:
!        MDIST(8*IRV-7) gives the dimension of the random variable
!        MDIST(8*IRV-6) gives the offset for the first integer parameter
!                       (parameter values are in array IRNDPAR)
!        MDIST(8*IRV-5) gives the number of integer parameters
!        MDIST(8*IRV-4) gives the offset for the first real parameter
!                       (parameter values are in array RANDPAR)
!        MDIST(8*IRV-3) gives the number of real parameters
!        MDIST(8*IRV-2) gives the offset for the first parameter referenced
!                       by location (parameter values are in array LRNDPAR)
!        MDIST(8*IRV-1) gives the number of parameters referenced by location
!        MDIST(8*IRV  ) gives the number of the stage
!
!
!     3. The location and process mode of the stochastic elements are stored
!        in array KSTOCH in groups of four. The first entry gives the type
!        of random element, the second its block (or node), the third entry
!        describes the location within the relevant array, the fourth specifies
!        whether the values are to be treated as replacing the core file
!        information or if they are to be treated as additive or multiplicative
!        perturbations.
!
!     Types of random elements are as follows (can be extended as needed):
!
!        KSTOCH(4*NSTELM-3) = 1 for stochastic entry in an A-matrix block
!        KSTOCH(4*NSTELM-3) = 2 for stochastic cost
!        KSTOCH(4*NSTELM-3) = 3 for stochastic RHS
!        KSTOCH(4*NSTELM-3) = 4 for stochastic range
!        KSTOCH(4*NSTELM-3) = 5 for stochastic lower bound
!        KSTOCH(4*NSTELM-3) = 6 for stochastic upper bound
!        KSTOCH(4*NSTELM-3) = 7 for stochastic entry in a Q-matrix block
!        KSTOCH(4*NSTELM-3) = 11 for stochastic upper and lower bounds ('FX' type)
!        KSTOCH(4*NSTELM-3) = 12 for stochastic demand (network option)
!        KSTOCH(4*NSTELM-3) = 13 for stochastic supply (network option)
!
!     ---------------------------------------------------------------
!            This version dated 22 April 2001. Gus Gassmann.
!     ---------------------------------------------------------------
!
      USE SAMPUTIL
      USE SMPS_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      implicit real*8 (a-h,o,p,r-z), integer*4(i-n), character*1(q)
!
! ----------------------------------------------------------------------
!
!     Start with some initializations
!
      IERR = 0
!
      DO I=1,NPER
         KREF(I) = I
      END DO
!
      NCTR = NODE0
      JPER = IPER1 - 1
  120 CONTINUE
         IF (JPER .GT. 0) THEN
            KREF(JPER) = NCTR
            NCTR = NODEINFO(NCTR)%IANCTR
            JPER = JPER - 1
            GOTO 120
         ENDIF
!
!     Go through all the random variables and simulate the ones that become
!     known in the period covered by this scenario.
!
         KDIM=0
         DO J=1,NRVAR
            IF (MDIST(8*J) .GE. IPER1) THEN
!
!     First prepare the parameters referenced by location (if there are any)
!
               IF (MDIST(8*J-1) .GT. 0) THEN
                  IF (MDIST(8*J-1) .GT. LLPAR) THEN
                     IF (.NOT. EXTMEM_LPAR(MDIST(8*J-1),IERR)) THEN
                        WRITE (IOLOG, 1100)
                        IERR = 1
                        RETURN
                     ENDIF
                  ENDIF
                  DO K=1,MDIST(8*J-1)
                     LROFF = 3*MDIST(8*J-2) + 3*(K-1)
                     XRNDPAR(K) = XTRACT(LRNDPAR(LROFF+1),      &
     &                                   LRNDPAR(LROFF+2),      &
     &                                   LRNDPAR(LROFF+3), KREF,      &
     &                                   HIST, IPER1-1, J-1, IERR)
                     IF (IERR .GT. 0) RETURN
                  END DO
               ENDIF
!
!     Next call the universal interface for all distributions
!
               CALL CALLDIST(SDIST(J),      &
     &                   IRNDPAR(MDIST(8*J-6)+1), MDIST(8*J-5),      &
     &                   RANDPAR(MDIST(8*J-4)+1), MDIST(8*J-3),      &
     &                   XRNDPAR(1),              MDIST(8*J-1),      &
     &                   HIST(KDIM+1),            MDIST(8*J-7),0,IERR)
               IF (NECHO1 .GE. 2) WRITE (IOLOG, 1000)      &
     &            J,MDIST(8*J),(HIST(JV),JV=KDIM+1,KDIM+MDIST(8*J-7))
!
            ENDIF
            IF (IERR .GT. 0) RETURN
            KDIM = KDIM + MDIST(8*J-7)
         END DO
!
!     Now store all the information into the data structure
!
      CALL SNGLBRCH(IPER1,NODE0,IERR)
      DO I=IPER1,NPER
         KREF(I) = I
      END DO
!
      RETURN
!
 1000 FORMAT(' Random variable',I4,' -- Stage',I4,1(/,10X,F12.6))
 1100 FORMAT('XXX -  FATAL  - Memory allocation failed in SAMPLE1.')
      END SUBROUTINE
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      REAL*8 FUNCTION XTRACT ( MTYPE,  MBLOCK, MLOC, KREF,      &
     &                         YPRIOR, LPRIOR, NRVP, IERR )
!
!     This function extracts one parameter referenced by location and
!     returns its current value. The type of parameter is described by
!     MTYPE, its location in the relevant array by MLOC. Possible types
!     are:
!
!      MTYPE = 1 for stochastic entry in one of the A-matrix blocks
!      MTYPE = 2 for stochastic cost
!      MTYPE = 3 for stochastic RHS
!      MTYPE = 5 for stochastic lower bound
!      MTYPE = 6 for stochastic upper bound
!      MTYPE = 7 for stochastic entry in one of the Q-matrix blocks
!      MTYPE = 11 for fixed stochastic bound (affects upper and lower)
!      MTYPE = 12 for stochastic demand (network option)
!      MTYPE = 13 for stochastic supply (network option)
!      MTYPE = 18 for reference to a primal decision variable
!      MTYPE = 19 for reference to a dual variable
!
!     Array KREF contains the reference nodes along the current scenario.
!
      USE SMPS_DATA
      USE IO_HANDLING
!
      REAL*8 YPRIOR(*)
!
!     Right-hand sides, costs, bounds and variable references are immediate
!
      INTEGER*4 KREF(*)
!
      IERR = 0
      IF (MTYPE .NE. 1 .AND. MTYPE .NE. 7) THEN
         IF (MBLOCK .LE. LPRIOR) THEN
            IF     (MTYPE .EQ.  2) THEN
               XTRACT = COST(NODEINFO(KREF(MBLOCK))%KCOST +MLOC)
            ELSEIF (MTYPE .EQ.  3) THEN
               XTRACT =  RHS(NODEINFO(KREF(MBLOCK))%KRHS  +MLOC)
            ELSEIF (MTYPE .EQ.  5) THEN
               XTRACT =  XLB(NODEINFO(KREF(MBLOCK))%KBOUND+MLOC)
            ELSEIF (MTYPE .EQ.  6) THEN
               XTRACT =  XUB(NODEINFO(KREF(MBLOCK))%KBOUND+MLOC)
            ELSEIF (MTYPE .EQ. 12) THEN
               XTRACT =  RHS(NODEINFO(KREF(MBLOCK))%KRHS  +MLOC)
            ELSEIF (MTYPE .EQ. 13) THEN
               XTRACT =  RHS(NODEINFO(KREF(MBLOCK))%KRHS  +MLOC)
            ELSEIF (MTYPE .EQ. 18) THEN
               XTRACT = 0.D0
               WRITE (IOLOG, 1100)
               IERR = 1
            ELSEIF (MTYPE .EQ. 19) THEN
               XTRACT = 0.D0
               WRITE (IOLOG, 1100)
               IERR = 1
            ENDIF
!
!     Here the reference is not part of the data structure yet.
!     This only works for proper data, not for decision variables.
!
         ELSE
            IF (MTYPE .EQ. 18 .OR. MTYPE .EQ. 19) GOTO 999
               GOTO 200
         ENDIF
!
!     Matrix elements require a little more work
!
      ELSEIF (MTYPE .EQ.  1) THEN
         DO I=1,NPER
            IF (KDATA(I+1) .GT. MBLOCK) GOTO 110
         END DO
         GOTO 999
!
  110    CONTINUE
            IF (I .GT. LPRIOR) GOTO 200
            IAREF  = KDATA(KREF(I)) + MBLOCK - KDATA(I)
            XTRACT = A(KELMA(IAREF)+MLOC)
      ELSEIF (MTYPE .EQ.  7) THEN
         DO I=1,NPER
            IF (KDATQ(I+1) .GT. MBLOCK) GOTO 130
         END DO
         GOTO 999
!
  130    CONTINUE
            IF (I .GT. LPRIOR) GOTO 200
            IQREF  = KDATQ(KREF(I)) + MBLOCK - KDATQ(I)
            XTRACT = AQMTX(IQOFF(IQREF)+MLOC)
      ELSE
         GOTO 999
      ENDIF
      RETURN
!
!     Here the reference is to a stochastic element that does not belong
!     to the previously observed history. In theory it can depend on some
!     of the random variables generated but not yet completely processed
!     (and in the Fleten, Wallace and Ziemba paper, it does).
!
!     We must go through all the stochastic elements to establish its
!     identity, then find from the D-matrix all the random variables that
!     impact on it and check that their values are known. In essence,
!     we compute V(i) = D(i,.)*Yprior and check that this equals D(i,.)*Y.
!
  200 CONTINUE
         DO ISTELM=1,NSTELM
            IF ((KSTOCH(4*ISTELM-3) .EQ. MTYPE ) .AND.      &
     &          (KSTOCH(4*ISTELM-2) .EQ. MBLOCK) .AND.      &
     &          (KSTOCH(4*ISTELM-1) .EQ. MLOC  )       ) GOTO 220
         END DO
         GOTO 999
!
  220    CONTINUE
            XTRACT = 0.D0
            KDIM = 0
            DO J=1,NRVP
               DO K=1,MDIST(8*J-7)
                  DO L=LDMTX(KDIM+K),LDMTX(KDIM+K+1)-1
                     IF (IDMTX(L) .EQ. ISTELM)      &
     &                  XTRACT = XTRACT + DMTX(L)*YPRIOR(KDIM+K)
                  END DO
               END DO
               KDIM = KDIM + MDIST(8*J-7)
            END DO
!
            DO J=NRVP+1,NRVAR
               DO K=1,MDIST(8*J-7)
                  DO L=LDMTX(KDIM+K),LDMTX(KDIM+K+1)-1
                     IF (IDMTX(L) .EQ. ISTELM) GOTO 999
                  END DO
               END DO
               KDIM = KDIM + MDIST(8*J-7)
            END DO
            RETURN
!
!     Error condition
!
  999 CONTINUE
      WRITE (IOLOG, 1000)
      XTRACT = 0.D0
      IERR = 2
      RETURN
!
 1000 FORMAT(' XXX -  FATAL  - Illegal parameter reference.',      &
     &       ' Zero returned.')
 1100 FORMAT(' XXX - WARNING - This feature not yet implemented.', &
     &       ' Zero returned.')
      END FUNCTION
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE SETPRB (IERR)
!
!     This routine sets the probabilities for the currently active sampled
!     tree. (It should also be called if the tree was only updated.)
!
!     There are two possibilities:
!        If ISMPRB = 0 (the default) then all paths have equal probability
!        If ISMPRB = 1 then all descendants of the same ancestor have equal
!                      conditional probability
!
      USE SMPS_DATA
      USE SAMPLING_DATA
      USE IO_HANDLING
!
      IERR = 0
      IF (ISMPRB .EQ. 0) THEN
         N_LEAF = 0
         DO I=NPER,1,-1
            IND = IRNGE0(I)
  100       CONTINUE
               IF (NODEINFO(IND)%NDESC .EQ. 0) THEN
                  N_LEAF = N_LEAF + 1
                  NODEINFO(IND)%PROB = 1.D0
               ELSE
                  JND = NODEINFO(IND)%IDESC
                  NODEINFO(IND)%PROB = 0.D0
  110             CONTINUE
                     NODEINFO(IND)%PROB = NODEINFO(IND)%PROB      &
     &                                  + NODEINFO(JND)%PROB
                     JND = NODEINFO(JND)%IBROTH
                     IF (JND .GT. 0) GOTO 110
               ENDIF
               IND = IABS(NODEINFO(IND)%IBROTH)
               IF (IND .GT. 0) GOTO 100
         END DO
!
         PP = 1.D0/N_LEAF
         DO I=NPER-1,1,-1
            IND = IRNGE0(I)
  130       CONTINUE
               IF (NODEINFO(IND)%NDESC .GT. 0) THEN
                  JND = NODEINFO(IND)%IDESC
  140             CONTINUE
                     NODEINFO(JND)%PROB = NODEINFO(JND)%PROB      &
     &                                  / NODEINFO(IND)%PROB
                     JND = NODEINFO(JND)%IBROTH
                     IF (JND .GT. 0) GOTO 140
               ENDIF
               IND = IABS(NODEINFO(IND)%IBROTH)
               IF (IND .GT. 0) GOTO 130
         END DO
         NODEINFO(IRNGE0(1))%PROB = 1.D0
!
      ELSE
         DO I=1,NPER-1
            IND = IRNGE0(I)
  180       CONTINUE
               IF (NODEINFO(IND)%NDESC .GT. 0) THEN
                  PP = 1.D0/NODEINFO(IND)%NDESC
                  JND = NODEINFO(IND)%IDESC
  190             CONTINUE
                     NODEINFO(JND)%PROB = PP
                     JND = NODEINFO(JND)%IBROTH
                     IF (JND .GT. 0) GOTO 190
               ENDIF
               IND = IABS(NODEINFO(IND)%IBROTH)
               IF (IND .GT. 0) GOTO 180
         END DO
         NODEINFO(IRNGE0(1))%PROB = 1.D0
      ENDIF
!
      RETURN
      END SUBROUTINE
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE TSAMPL(NSAMPL,IERR)
!
!     This subroutine samples NSAMPL scenarios from an existing scenario tree.
!     Leaf nodes in the tree are identified and are marked, depending on whether
!     they will form part of the sample or not. The conditional probabilities
!     are adjusted. The pointer array IANCTR is not altered, so as to make it
!     possible to reconstruct the original event tree later.
!
!     ----------------------------------------------------------------------
!          This version dated 7 June 1998. Written by Gus Gassmann.
!     ----------------------------------------------------------------------
!
      USE SAMPUTIL, ONLY: RNDGET
      USE SMPS_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      implicit real*8 (a-h,o,p,r-z), integer*4(i-n), character*1(q)
!
      REAL*8    PPROB,RVEC(1),CPROB
!
!     This routine should only be called if the tree is given explicitly.
!     If a sample is currently pending (TREEFMT = 1), call TRESET first.
!
      IF (NPER    .EQ. 1) GOTO 500
      IF (STFILE  .GT. 2) GOTO 510
      IF (TREEFMT .EQ. 1) CALL TRESET
!
!     We cannot select a scenario more than once, so NSAMPL must be
!     less than the number of leaf nodes.
!
      IF (NSAMPL .GT. LASTLF) THEN
         WRITE (IOLOG, 1100) LASTLF
         GOTO 999
      ELSEIF (NSAMPL .EQ. LASTLF) THEN
         GOTO 999
!
!     If more than half of the existing scenarios should be selected,
!     randomly determine the ones that will be left out instead.
!
      ELSEIF (NSAMPL .GT. LASTLF/2) THEN
         MARK = 0
         LEAVE = 1
         NMARK = LASTLF - NSAMPL
      ELSE
         MARK = 1
         LEAVE = 0
         NMARK = NSAMPL
      ENDIF
!
!     Allocate temporary arrays
!
      ALLOCATE (QPROB(IZNODE), ITAKE(IZNODE), STAT=IERR)
      IF (IERR .NE. 0) GOTO 960
!
!     Initial the selection array and compute the path probabilities.
!
      DO I=1,LASTLF
         IND = LEAFND(I)
         ITAKE(IND) = LEAVE
         PPROB = NODEINFO(IND)%PROB
  100    CONTINUE
            IND = NODEINFO(IND)%IANCTR
            IF (IND .GT. 0) THEN
               ITAKE(IND) = LEAVE
               PPROB = PPROB * NODEINFO(IND)%PROB
               GOTO 100
            ENDIF
!
         IF (I .EQ. 1) THEN
            QPROB(I) = PPROB
         ELSE
            QPROB(I) = PPROB + QPROB(I-1)
         ENDIF
      END DO
!
!     Now select the scenarios one at a time.
!
      I = 1
  120 CONTINUE
         CALL RNDGET(RVEC,1,IERR)
         TARGET = RVEC(1)
!
         DO J=1,LASTLF
            IF (TARGET .LT. QPROB(J)) GOTO 140
         END DO
!
  140    CONTINUE
            IND = LEAFND(J)
            IF (ITAKE(IND) .EQ. MARK) GOTO 120
  150    CONTINUE
            ITAKE(IND) = MARK
            IF (MARK .EQ. 1) THEN
               IND = NODEINFO(IND)%IANCTR
               IF (IND .GT. 0) GOTO 150
            ENDIF
            I = I + 1
            IF (I .LE. NMARK) GOTO 120
  160 CONTINUE
!
!     If we deselect the scenarios, we have to be careful when marking nodes
!
      IF (MARK .EQ. 0) THEN
         DO 190 I=NPER,1,-1
            ND = IRNGE0(I)
  170       CONTINUE
            IF (NODEINFO(ND)%NDESC .GT. 0) THEN
               ID = NODEINFO(IND)%IDESC
  180          CONTINUE
                  IF (ITAKE(ID) .EQ. 0) THEN
                     ID = NODEINFO(ID)%IBROTH
                     IF (ID .GT. 0) GOTO 180
                  ELSE
                     ND = ABS(NODEINFO(ND)%IBROTH)
                     IF (ND .GT. 0) GOTO 170
                        GOTO 190
                  ENDIF
                  ITAKE(ND) = 0
            ENDIF
            ND = ABS(NODEINFO(IND)%IBROTH)
            IF (ND .GT. 0) GOTO 170
  190    CONTINUE
      ENDIF
!
!     Set all the descendant pointers to zero. (We will be able
!     to reconstruct them later from the ancestor pointers.)
!
      DO I=1,LASTLF
         IND = LEAFND(I)
  200    CONTINUE
            IND = NODEINFO(IND)%IANCTR
            IF (IND .GT. 0) THEN
               NODEINFO(IND)%IDESC = 0
               NODEINFO(IND)%NDESC = 0
               GOTO 200
            ENDIF
      END DO
!
!     Now adjust the tree links and the conditional probabilities for the
!     nodes that lie on the selected paths.
!
      DO I=1,NPER
         LAST = 0
         IND = IRNGE0(I)
         IRNGE0(I) = 0
  220    CONTINUE
            IF (ITAKE(IND) .EQ. 1) THEN
               IF (IRNGE0(I) .EQ. 0) THEN
                   IRNGE0(I) = IND
                   CPROB = NODEINFO(IND)%PROB
                   ND1 = IND
               ELSE
                   IF (LAST .GT. 0) THEN
                      IF (NODEINFO(LAST)%IANCTR .EQ.      &
     &                    NODEINFO(IND)%IANCTR) THEN
                         NODEINFO(LAST)%IBROTH = IND
                         CPROB = CPROB + NODEINFO(IND)%PROB
                      ELSE
                         NODEINFO(LAST)%IBROTH = -IND
  230                    CONTINUE
                         NODEINFO(ND1)%PROB = NODEINFO(ND1)%PROB/CPROB
                         ND1 = NODEINFO(ND1)%IBROTH
                         IF (ND1 .GT. 0) GOTO 230
!
                         CPROB = NODEINFO(IND)%PROB
                         ND1 = IND
                      ENDIF
                  ENDIF
               ENDIF
               LAST = IND
               NCTR = NODEINFO(IND)%IANCTR
               IF (NCTR .NE. 0) THEN
                  IF (NODEINFO(NCTR)%IDESC .EQ. 0)      &
     &                NODEINFO(NCTR)%IDESC = IND
                  NODEINFO(NCTR)%NDESC = NODEINFO(NCTR)%NDESC + 1
               ENDIF
            ENDIF
            IND = ABS(NODEINFO(IND)%IBROTH)
            IF (IND .GT. 0) GOTO 220
               NODEINFO(LAST)%IBROTH = 0
  240          CONTINUE
                  NODEINFO(ND1)%PROB = NODEINFO(ND1)%PROB / CPROB
                  ND1 = NODEINFO(ND1)%IBROTH
                  IF (ND1 .GT. 0) GOTO 240
      END DO
      TREEFMT = 1
      GOTO 999
!
  500 CONTINUE
         WRITE (IOLOG, 1000)
         GOTO 999
!
  510 CONTINUE
         WRITE (IOLOG, 1200)
         IERR = 1
         GOTO 999
!
  960 CONTINUE
         WRITE (IOLOG, 1300)
         IERR = 1
  999 CONTINUE
      DEALLOCATE (QPROB, ITAKE, STAT=IER)
      IF (IER .NE. 0) IERR = 1
      RETURN
!
 1000 FORMAT(' This feature still under construction...')
 1100 FORMAT(' ERROR - The scenario tree contains only',I4,' paths')
 1200 FORMAT(' ERROR - Only use this routine for explicit event trees.')
 1300 FORMAT(' XXX -  FATAL  - Memory handling error in TSAMPL.')
 1510 FORMAT('(''Smpl'',I',i1,')')
 1520 FORMAT('(''Smpl'',I',i2,')')
      END SUBROUTINE
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      SUBROUTINE TRESET
!
!     This subroutine throws away information from a sampled event tree.
!     If the problem was originally specified with implicit distributions,
!     then this amounts to simply truncating the tree after the first
!     NPER nodes. If the event tree was specified explicitly (via SCENARIOS
!     or NODES), then we have a little more work to do: some leaf nodes
!     that may have been ignored before will have to be reinstated, and some
!     conditional probabilities may have to be adjusted.
!
      USE SMPS_DATA
      USE IO_HANDLING
      USE MEMORY_EXTENDER
      implicit real*8 (a-h,o,p,r-z), integer*4(i-n), character*1(q)
!
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: INUSE, JPATH, LASTND
!
      ALLOCATE (QPROB(0:NODES), INUSE(NODES), JPATH(NPER),      &
     &          LASTND(NPER),   STAT=IERR)
      IF (IERR .NE. 0) GOTO 960
!
      IF (TREEFMT .NE. 1) GOTO 999
!
!     Here we have implicit event trees
!
      IF (STFILE .GT. 2) THEN
         DO I=1,NPER
            NODEINFO(I)%IANCTR = I - 1
            NODEINFO(I)%IBROTH = 0
            IF (I .LT. NPER) THEN
               NODEINFO(I)%IDESC = I + 1
               NODEINFO(I)%NDESC = 1
            ELSE
               NODEINFO(I)%IDESC = 0
               NODEINFO(I)%NDESC = 0
            ENDIF
            IRNGE0(I) = I
         END DO
         NODES = NPER
!
         MAXROW = NODEINFO(NODES)%KROW   + NODEINFO(NODES)%NROW
         MAXCOL = NODEINFO(NODES)%KCOL   + NODEINFO(NODES)%NCOL + 1
         LASTC  = NODEINFO(NODES)%KCOST  + NODEINFO(NODES)%NCOL      &
     &                                   - NODEINFO(NODES)%NROW
         LASTR  = NODEINFO(NODES)%KRHS   + NODEINFO(NODES)%NROW
         LASTBD = NODEINFO(NODES)%KBOUND + NODEINFO(NODES)%NCOL + 1
         LASTLF = 0
!
         IF (LCHANCE .GT. 0)      &
     &       LCHANCE = KCHANCE(NODES) + NCHANCE(NODES)
         IF (LICC .GT. 0)      &
     &       LICC    = KICC(NODES)    + NICC(NODES)
         IF (LASTPLQ .GT. 0)      &
     &       LASTPLQ = LASTC
!
         IF (MARKOV) THEN
            NEWB = 2*NPER - 1
         ELSE
            NEWB = NPER*(NPER+1)/2
         ENDIF
         IF (HAVE_GC) NEWB = NEWB + NPER*(NPER-1)/2
         NEWA = 0
         NEWC = 0
         DO I=1,NEWB
            IF (NEWA .LT. KELMA(I)) NEWA = KELMA(I)
            IF (NEWC .LT. KCOLA(I)) NEWC = KCOLA(I)
         END DO
         DO I=NEWB+1,LASTBA
            IF (KELMA(I) .GT. NEWA .AND. KELMA(I) .LT. LASTA)      &
     &         LASTA = KELMA(I)
            IF (KCOLA(I) .GT. NEWC .AND. KCOLA(I) .LT. LASTCA)      &
     &         LASTCA = KCOLA(I)
         END DO
         LASTBA = NEWB
!
         IF (NQPROF .GT. 0) THEN
            NEWB = NQPROF*NPER - NQPROF*(NQPROF-1)/2
            NEWQ = 0
            NEWC = 0
            DO I=1,NEWB
               IF (NEWQ .LT. IQOFF(I)) NEWQ = IQOFF(I)
               IF (NEWC .LT. LQOFF(I)) NEWC = LQOFF(I)
            END DO
            DO I=NEWB+1,LASTQB
               IF (IQOFF(I) .GT. NEWQ .AND. IQOFF(I) .LT. LASTQ)      &
     &            LASTQ = IQOFF(I)
               IF (LQOFF(I) .GT. NEWC .AND. LQOFF(I) .LT. LASTQC)      &
     &            LASTQC = LQOFF(I)
            END DO
            LASTQB = NEWB
         ENDIF
!
!     Here the event tree is given in the stoch file in explicit form
!
      ELSEIF (STFILE .GT. 0) THEN
         QPROB(0:NODES) = 0.D0
         INUSE(1:NODES) = 0
!
         DO I=1,NPER
            ND0 = IRNGE0(I)
            IRNGE0(I) = 0
            LASTND(I) = 0
  210       CONTINUE
               INUSE(ND0) = 1
               NODEINFO(ND0)%IDESC = 0
               NODEINFO(ND0)%NDESC = 0
               ND0 = ABS(NODEINFO(ND0)%IBROTH)
               IF (ND0 .GT. 0) GOTO 210
         END DO
!
         DO I=1,LASTLF
            ND = LEAFND(I)
  230       CONTINUE
               IF (INUSE(ND) .EQ. 0) THEN
                   INUSE(ND) = 2
                   NC = NODEINFO(ND)%IANCTR
                   QPROB(NC) = QPROB(NC) + NODEINFO(ND)%PROB
                   IF (NC .GT. 0) THEN
                      ND = NC
                      GOTO 230
                   ENDIF
               ENDIF
         END DO
!
         DO I=1,LASTLF
            ND = LEAFND(I)
            IP = MXTPER
  260       CONTINUE
               IF (INUSE(ND) .EQ. 1)  NODEINFO(ND)%PROB =       &
     &             NODEINFO(ND)%PROB *(1.D0-QPROB(NODEINFO(ND)%IANCTR))
               IAN = NODEINFO(ND)%IANCTR
               JPATH(IP) = ND
               IP = IP - 1
               IF (IAN .GT. 0) THEN
                  IF (INUSE(ND) .GT. 0) THEN
                     INUSE(ND) = 0
                     IF (NODEINFO(IAN)%IDESC .EQ. 0)      &
     &                   NODEINFO(IAN)%IDESC = ND
                     NODEINFO(IAN)%NDESC = NODEINFO(IAN)%NDESC + 1
                  ENDIF
                  ND = IAN

                  GOTO 260
               ENDIF
!
               DO J=1,MXTPER-IP
                  IF (JPATH(IP+J) .NE. LASTND(J)) THEN
                     IF (LASTND(J) .EQ. 0) THEN
                        IRNGE0(J) = JPATH(IP+J)
                     ELSE
                        IF (NODEINFO(JPATH(IP+J))%IANCTR .EQ.       &
     &                      NODEINFO(LASTND(J))%IANCTR) THEN
                            NODEINFO(LASTND(J))%IBROTH =  JPATH(IP+J)
                        ELSE
                            NODEINFO(LASTND(J))%IBROTH = -JPATH(IP+J)
                        ENDIF
                     ENDIF
                     LASTND(J) = JPATH(IP+J)
                  ENDIF
               END DO
         END DO
      ENDIF
      GOTO 999
!
  960 CONTINUE
         WRITE (IOLOG, 1000)
         IERR = 1
  999 CONTINUE
      DEALLOCATE (QPROB, INUSE, JPATH, LASTND, STAT=IER)
      IF (IER .NE. 0) IERR = 1
      RETURN
!
 1000 FORMAT(' XXX -  FATAL  - TRESET: Memory handling error.')
      END SUBROUTINE
!
      END MODULE
