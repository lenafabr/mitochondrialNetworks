SUBROUTINE READKEY
  ! this subroutine reads in keywords from a parameter file
  ! it sets the various global variables defined in KEYS module
  ! name of the parameter file is param.* where * is a keyword argument
  ! if no keyword argument is supplied, the default is just a file called param
  ! The EXTRAPARAMFILES keyword will allow extra parameter files to be
  ! read in as well

  USE KEYS
  USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
  USE GENUTIL

  IMPLICIT NONE

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  INTEGER :: NUMARG ! number of command line arguments
  INTEGER :: NITEMS ! number of items on the line in the parameter file
  INTEGER :: PF ! input file unit
  LOGICAL :: FILEEND=.FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------
  INTEGER, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  INTEGER :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing random number generator
  INTEGER :: TIMEVAL(8), SEED
  ! ---------------- temporary variables ---------------
  INTEGER :: I, SC !DUMI,TMPI, DUMI1, DUMI2, DUMI3,
  CHARACTER*100 :: DUMSTR
  LOGICAL :: LDUM

  LOGICAL :: ALPHA1SET, ALPHA2SET, ALPHAPLANESET, DIFFSET, DECAYSET, FISSRADSET

  ! ------------------------
  ! set variable defaults
  ! ------------------------
  ACTION = 'NONE'
  RNGSEED = 0
  VERBOSE = .FALSE.

  ! coordinate dimension
  DIM = 3

  ! geometry and energy parameters
  MITOLEN = 5D-1; ! mitochondrial segment length
  ESTR = 1D3; ! stretching energy
  ! confinement (radius, prefactor, power)
  CELLRAD1 = -1D0 ! negative = no confinement
  ECONF = 1D4
  ! bending parameters (desired angle, energy prefactors)
  THETA0 = PI/3
  BENDMOD1 = 2D0
  BENDMOD2 = 6D0

  ! steric parameters
  STERICMOD = 1D4
  STERICRAD = 1D-1
  ! grid spacing
  USEGRID = .TRUE.
  USEEDGEGRID = .FALSE.
  USENODEGRID = .FALSE.
  GRIDSPACES = 10 ! number of grid partitions along each spatial dimension

  ! yeast model
  DOYEAST = .FALSE.
  YEASTCONF = 1D2
  YEASTONRATE = 1D0
  YEASTOFFRATE = 1D0
  YEASTBINDRANGE = 1D0

  ! fission / fusion
  ! fission rate per time per fissable node
  FISSRATE = 0D0
  ! multiplier for D3 fissions (relative to D2 fission rate)
  D3FISSMULT = 15D-1
  ! fission separation radius
  FISSRAD = 1D-1
  !fusion rate per time per node pair in range
  FUSERATE1 = 0D0
  FUSERATE2 = 0D0
  
  ! fusion bending persistence lengths
  ALPHA1 = 4D0
  ALPHA2 = 12D0

  !contact radius assumed for spherical-tipped mitochondria
  CONTACTRAD = 15D-2
  !recharge fusion machinery rate
  RECHARGERATE = -1D0

  ! using out of plane model?
  USEOUTOFPLANE = .FALSE.
  BENDMODPLANE = 0.6
  ALPHAPLANE = 1.2D0

  
  ! max number of nodes and edges
  MAXNNODE = 500
  MAXNEDGE = 250

  ! number of species & diffusivities
  NSPECIES = 0
  DIFF = 0D0
  NSUBSTEPS = 1

  ! production rates
  PRODRATE = -1D0
  NUMTRIGGERUNITS = 0
  DECAYRATE = -1D0
  FIXEDGES = .FALSE.
  UNIQUETRIGGERS = .FALSE.

  ! input/output
  USERANDOMFRAGMENTEDNETWORK = .FALSE.
  NETFILE = '*.net'
  PLANEFILE = '*.pl'
  OUTFILE = '*.out'
  !DUMPSNAPSHOTS = .FALSE. ! periodically dump chain snapshots
  SNAPSHOTEVERY = 1 ! how often to dump snapshots
  SNAPSHOTFILE = '*.snap.out' ! snapshot file
  REMODELINGFILE = '*.ffevents.out' ! fiss/fuse events file
  APPENDSNAPSHOTS = .FALSE. ! append snapshots to file rather than replacing
  PRINTEVERY = 1 ! how often to print output (to screen and OUTFILE)
  STARTSNAPSTEP =1 ! start saving snapshots at and beyond this step
  STARTDIFFSTEP = 1 ! start propagating diffusion at and beyond this step

  ! brownian dynamics
  DELT = 1D-4 ! time step
  FRICT = 1D0 ! friction coefficient
  KT = 1D0 ! temperature
  BDSTEPS = 1000 ! number of brownian steps to run
  DOBROWN = .TRUE. ! include brownian forces

  ALPHA1SET = .FALSE.
  ALPHA2SET = .FALSE.
  ALPHAPLANESET = .FALSE.
  DIFFSET = .FALSE.
  DECAYSET = .FALSE.
  FISSRADSET = .FALSE.
  
  
  ! -------------------------
  ! Read in all parameter files, starting with the ones specified on command line
  ! --------------------------

  PF = 55 ! i/o unit number to be used for parameter files

  ! get input parameter files from command line
  NPARAMFILES = 0
  NUMARG = COMMAND_ARGUMENT_COUNT()
  IF (NUMARG==0) THEN
     NPARAMFILES = 1
     PARAMFILES(1) = 'param'
     ARG = ''
  ELSE
     DO I = 1,NUMARG
        CALL GETARG(I, ARG)
        NPARAMFILES = NPARAMFILES + 1
        WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))
        PARAMFILES(NPARAMFILES) = DUMSTR
     ENDDO
     ! reset arg to its original value
     IF (NUMARG.GT.1) CALL GETARG(1,ARG)
  ENDIF

  NPARAMREAD = 0 ! keep track of how many files have been read
  DO WHILE (NPARAMREAD.LT.NPARAMFILES)
     NPARAMREAD = NPARAMREAD + 1

     PRINT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
     INQUIRE(FILE=PARAMFILES(NPARAMREAD),EXIST=LDUM)
     IF (.NOT.LDUM) THEN
        PRINT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
        STOP 1
     ENDIF
     OPEN(UNIT=PF, FILE=PARAMFILES(NPARAMREAD), STATUS='OLD')

     ! read in the keywords one line at a time
     DO
        CALL READLINE(PF,FILEEND,NITEMS)
        IF (FILEEND.and.nitems.eq.0) EXIT

        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE

        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)

        ! Skip any empty lines or any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE

        SELECT CASE(WORD) ! pick which keyword
        CASE('ACTION')
           CALL READA(ACTION, CASESET=1)
        CASE('BDSTEPS')
           CALL READI(BDSTEPS)
        CASE('CONFSPHERE')
           CALL READF(CELLRAD1)
           IF (NITEMS.GT.2) CALL READF(ECONF)
        CASE('ECONF')
           CALL READF(ECONF)
        CASE('DELT')
           CALL READF(DELT)
        CASE('DIM')
           CALL READI(DIM)
        CASE('NSPECIES')
           CALL READI(NSPECIES)
        CASE('DIFF')
           IF(NITEMS.EQ.2) THEN
              CALL READF(DIFF(1))
              DIFF = DIFF(1)
           ELSE
              DO SC = 1,NITEMS-1
                CALL READF(DIFF(SC))
              ENDDO
           ENDIF
           DIFFSET = .TRUE.
        CASE('DIFF1')
           IF(.NOT.DIFFSET) CALL READF(DIFF(1))
        CASE('DIFF2')
           IF(.NOT.DIFFSET) THEN
              CALL READF(DIFF(2))
              DIFF(2:NSPECIES) = DIFF(2)
           ENDIF
        CASE('NSUBSTEPS')
           CALL READI(NSUBSTEPS)
        CASE('ESTR')
           CALL READF(ESTR)
        CASE('BENDMOD1')
           CALL READF(BENDMOD1)
        CASE('BENDMOD2')
           CALL READF(BENDMOD2)
        CASE('BENDMODPLANE')
           CALL READF(BENDMODPLANE)
        CASE('USEOUTOFPLANE')
           CALL READO(USEOUTOFPLANE)
        CASE('THETA0')
           CALL READF(THETA0)
        CASE('STERICMOD')
           CALL READF(STERICMOD)
        CASE('STERICRAD')
           CALL READF(STERICRAD)
         CASE('FISSRAD')
            CALL READF(FISSRAD)
            FISSRADSET = .TRUE.
         CASE('D3FISSMULT')
            CALL READF(D3FISSMULT)
        CASE('FISSRATE')
           CALL READF(FISSRATE)
        CASE('FUSERATE')
           CALL READF(FUSERATE1)
           FUSERATE2 = FUSERATE1
           IF (NITEMS.GT.2) CALL READF(FUSERATE2)  
        CASE('FUSERATE1')
           CALL READF(FUSERATE1)
        CASE('FUSERATE2')
           CALL READF(FUSERATE2)
         CASE('ALPHA1')
            CALL READF(ALPHA1)
            ALPHA1SET = .TRUE.
         CASE('ALPHA2')
            CALL READF(ALPHA2)
            ALPHA2SET = .TRUE.
         CASE('ALPHAPLANE')
            CALL READF(ALPHAPLANE)
            ALPHAPLANESET = .TRUE.
        CASE('CONTACTRAD')
           CALL READF(CONTACTRAD)
        CASE('RECHARGERATE')
           CALL READF(RECHARGERATE)
        CASE('USEGRID')
           CALL READO(USEGRID)
         CASE('USEEDGEGRID')
            CALL READO(USEEDGEGRID)
         CASE('USENODEGRID')
            CALL READO(USENODEGRID)
        CASE('GRIDSPACES')
           IF(NITEMS-1.EQ.DIM) THEN
              DO SC = 1,NITEMS-1
                 CALL READI(GRIDSPACES(SC))
              ENDDO
           ELSE
              CALL READI(GRIDSPACES(1))
              GRIDSPACES = GRIDSPACES(1)
           ENDIF
        CASE('DOYEAST')
           CALL READO(DOYEAST)
        CASE('YEASTCONF')
           CALL READF(YEASTCONF)
        CASE('YEASTONRATE')
           CALL READF(YEASTONRATE)
        CASE('YEASTOFFRATE')
           CALL READF(YEASTOFFRATE)
        CASE('YEASTBINDRANGE')
           CALL READF(YEASTBINDRANGE)      
        CASE('FRICT')
           CALL READF(FRICT)
        CASE('KT')
           CALL READF(KT)
        CASE('DOBROWN')
           CALL READO(DOBROWN)
        CASE('MAXNEDGE')
           CALL READI(MAXNEDGE)
        CASE('MAXNNODE')
           CALL READI(MAXNNODE)
        CASE('MITOLEN')
           CALL READF(MITOLEN)
        CASE('USERANDOMFRAGMENTEDNETWORK')
           CALL READO(USERANDOMFRAGMENTEDNETWORK)
        CASE('NETFILE')
           CALL READA(NETFILE)
        CASE('PLANEFILE')
           CALL READA(PLANEFILE)   
        CASE('NOBROWN')
           DOBROWN = .FALSE.
        CASE('OUTFILE')
           CALL READA(OUTFILE)
        CASE('PRINTEVERY')
           CALL READI(PRINTEVERY)
        CASE('PRODRATE')
           IF (NITEMS-1.EQ.MAXSPECIES) THEN
             DO I = 1,NITEMS-1
                IF (I.GT.MAXSPECIES) EXIT
                CALL READF(PRODRATE(I))
             ENDDO
           ELSE
              CALL READF(PRODRATE(1))
              PRODRATE(2:MAXSPECIES) = PRODRATE(1)
           ENDIF
        CASE('DECAYRATE')
           IF (NITEMS-1.EQ.MAXSPECIES) THEN
              DO I = 1,NITEMS-1
                 IF(I.GT.MAXSPECIES) EXIT 
                 CALL READF(DECAYRATE(I))
              ENDDO
           ELSE
              CALL READF(DECAYRATE(1))
              DECAYRATE(2:MAXSPECIES) = DECAYRATE(1)
           ENDIF
           DECAYSET = .TRUE.
         CASE('DECAYRATE1')
            IF(.NOT.DECAYSET) CALL READF(DECAYRATE(1))
         CASE('DECAYRATE2')
            IF(.NOT.DECAYSET) THEN
               CALL READF(DECAYRATE(2))
               DECAYRATE(2:MAXSPECIES) = DECAYRATE(2)
            ENDIF
        CASE('NUMTRIGGERUNITS')
           CALL READI(NUMTRIGGERUNITS)
        CASE('FIXEDGES')
           CALL READO(FIXEDGES)
        CASE('UNIQUETRIGGERS')
           CALL READO(UNIQUETRIGGERS)
        CASE('RNGSEED')
           CALL READI(RNGSEED)
        CASE('SNAPSHOTFILE')
           CALL READA (SNAPSHOTFILE)
        CASE('REMODELINGFILE')
           CALL READA(REMODELINGFILE)
        CASE('SNAPSHOTS')
           DUMPSNAPSHOTS = .TRUE.
           IF (NITEMS.GT.1) CALL READI(SNAPSHOTEVERY)
           IF (NITEMS.GT.2) CALL READA(SNAPSHOTFILE)
           IF (NITEMS.GT.3) CALL READO(APPENDSNAPSHOTS)
        CASE('STARTDIFFSTEP')
           CALL READI(STARTDIFFSTEP)
        CASE('STARTSNAPSTEP')
           CALL READI(STARTSNAPSTEP)
        CASE('VERBOSE')
           CALL READO(VERBOSE)
        CASE DEFAULT
           print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
        END SELECT
     ENDDO
     CLOSE(PF)
  ENDDO


  ! -----------------
  ! check validity of some values, raise errors or adjust as necessary
  ! -----------------
  IF(NSPECIES.GT.MAXSPECIES) THEN
     PRINT*,'TOO MANY SPECIES. MAX IS 9', NSPECIES
     STOP 1
  ENDIF

  IF(.NOT.ALPHA1SET) ALPHA1 = BENDMOD1/(KT*MITOLEN)
  IF(.NOT.ALPHA2SET) ALPHA2 = BENDMOD2/(KT*MITOLEN)
  IF(.NOT.ALPHAPLANESET) ALPHAPLANE = BENDMODPLANE/(KT*MITOLEN)

  IF(.NOT.FISSRADSET) FISSRAD = STERICRAD

  IF (MITOLEN.LT.0) THEN
     PRINT*, 'ERROR IN MITOLEN VALUE', MITOLEN
     STOP 1
  ENDIF
  IF (ESTR.LT.0) THEN
     PRINT*, 'ERROR IN ESTR VALUE', ESTR
     STOP 1
  ENDIF
  IF (BENDMOD1.LT.0) THEN
     PRINT*, 'ERROR IN BENDMOD1 VALUE', BENDMOD1
     STOP 1
  ENDIF
  IF (BENDMOD2.LT.0) THEN
     PRINT*, 'ERROR IN BENDMOD2 VALUE', BENDMOD2
     STOP 1
  ENDIF
  IF (STERICMOD.LT.0) THEN
     PRINT*, 'ERROR IN STERICMOD VALUE', STERICMOD
     STOP 1
  ENDIF
  IF (MAXNNODE.LT.2.OR.MAXNEDGE.LT.1) THEN
     PRINT*, 'ERROR IN MAXNNODE OR MAXNEDGE VALUE', MAXNNODE,MAXNEDGE
     STOP 1
  ENDIF
  IF((USEEDGEGRID.OR.USENODEGRID).AND.(.NOT.USEGRID)) THEN
     PRINT*,"Must use base grid in order to use edge and/or node grids"
     STOP 1
  ENDIF

  ! ----------- fix file names -----------
  CALL REPLACESUBSTR(OUTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(NETFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(PLANEFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(SNAPSHOTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(REMODELINGFILE,'*',TRIM(ADJUSTL(ARG)))
  ! ---------------------------

  ! Initiate random number generator
  IF (RNGSEED.EQ.0) THEN
     ! use the current time of day in milliseconds
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = TIMEVAL(5)*3600*1000 + TIMEVAL(6)*60*1000 + TIMEVAL(7)*1000 + TIMEVAL(8)
  ELSEIF (RNGSEED.EQ.-1) THEN
     ! use the last 5 characters in the command-line argument
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)))
  ELSEIF (RNGSEED.EQ.-2) THEN
     ! use the last 4 characters in the command-line argument
     ! and additionally the millisecond time
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)),TIMEVAL(8))
  ELSE
     ! use this seed directly
     SEED = RNGSEED
  ENDIF

  print*, 'Initiating Mersenne twister random number generator with seed:', SEED
  CALL SGRND(SEED)

  print*, '------------Parameter values : -------------------'
  print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
  print*, 'Output file: ', TRIM(OUTFILE)
  IF (DUMPSNAPSHOTS) THEN
     PRINT*, 'Dumping snapshot every', SNAPSHOTEVERY,'steps. In file:', TRIM(ADJUSTL(SNAPSHOTFILE))
  ENDIF

  print*, 'Max allowed number of nodes and edges:', MAXNNODE,MAXNEDGE
  PRINT*,"confined to sphere of size, conf strength",CELLRAD1,ECONF
  print*, 'MITOLEN, ESTR, BENDMOD1, BENDMOD2, THETA0:', MITOLEN, ESTR, BENDMOD1, BENDMOD2, THETA0
  PRINT*,"use out-of-plane bending penalty?, Bending modulus",USEOUTOFPLANE,BENDMODPLANE
  print*, 'STERICMOD, STERIC RADIUS:', STERICMOD, STERICRAD
  print*, 'USEGRID?, USENODEGRID?, USEEDGEGRID? GRIDSPACES (X,Y,Z):',USEGRID, USENODEGRID, USEEDGEGRID, GRIDSPACES(1:DIM)

  PRINT*,'Starting diffusion at timestep, NUMBER OF DIFFUSING SPECIES, D COEFF:',STARTDIFFSTEP,NSPECIES,DIFF(1:NSPECIES)
  PRINT*,"number of substeps to take when doing FVM particle diffusion",NSUBSTEPS
  PRINT*,'number of producing edges, Fix edges at constant concentration?',NUMTRIGGERUNITS,FIXEDGES
  PRINT*,'production rate of each species',PRODRATE
  PRINT*,'decay rates of each species:',DECAYRATE
  PRINT*,'Do Yeast Model? Confinement strength:',DOYEAST,YEASTCONF
  PRINT*,'Yeast on rate, off rate, binding range',YEASTONRATE,YEASTOFFRATE,YEASTBINDRANGE
  PRINT*,'FISSION RATE, FUSION RATES (tip-tip,tip-side)',FISSRATE, FUSERATE1, FUSERATE2
  PRINT*,'FISSION multiplier for D3, Fission Separation Radius, Contact Radius, recharge fusion rate', D3FISSMULT, FISSRAD, &
  & CONTACTRAD, RECHARGERATE
  PRINT*,'Fusion Sensitivity parameters (alpha_i)',ALPHA1,ALPHA2, ALPHAPLANE
  PRINT*, 'Friction coefficient, KT:', FRICT, KT
  PRINT*, 'Time step, number steps, DOBROWN?:', DELT, BDSTEPS, DOBROWN
  print*, '----------------------------------------------------'


END SUBROUTINE READKEY
