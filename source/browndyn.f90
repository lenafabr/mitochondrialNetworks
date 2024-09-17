MODULE BROWNDYN
  ! Utilities for running brownian dynamics simulations
  USE KEYS, ONLY : VERBOSE
  IMPLICIT NONE
  
  TYPE DYNPARAM
     ! structure listing various parameters for the dynamic simulations
     ! (that don't necessarily belong to a particular network)
     CHARACTER(LEN=100) :: OUTFILE,SNAPSHOTFILE,REMODELINGFILE
     INTEGER :: PRINTEVERY, SNAPSHOTEVERY
     ! When to start propagating diffusion and saving snapshots
     INTEGER :: STARTDIFFSTEP, STARTSNAPSTEP
     LOGICAL :: APPENDSNAPSHOTS, DOBROWN
     ! number of substeps for diffusion
     INTEGER :: NSUBSTEPS

     LOGICAL :: TRACKENERGY
  END type DYNPARAM

CONTAINS
  SUBROUTINE RUNBROWNDYNSIM(NETP, NSTEP, DELT, DP, STARTTIME, CLOCKRATE)
    ! run a brownian dynamics simulation
    ! NSTEP: number of steps to simulate
    ! DELT: timestep
    ! OUTFILE: output file for printing out info
    ! (currently prints mechanical energy, but can switch to more useful output)
    ! PRINTEVERY: how often to print info into OUTFILE
    ! SNAPSHOTFILE: file for dumping chain configuration snapshots
    ! SNAPSHOTEVERY: how often to dump snapshots
    ! APPENDSNAPSHOTS: append snapshots to existing files (otherwise, start from scratch and append as we go)
    ! DOBROWN: include brownian forces
    ! parameters for the dynamic simulation
    USE DYNNETWORKUTIL, ONLY : DYNNETWORK, GETMECHFORCES, DYNNETSNAPSHOT, GETNEIGHBNODES, OUTSIDEPLANES
    USE DYNAMICPROCESSES, ONLY : FISSIONEVENT, FUSIONEVENT, DIFFUSEONEDGES, PRODUCEANDDECAYSPECIES, GETFUSIONRATE
    USE MT19937, ONLY : GRND
    USE GENUTIL, ONLY : RANDOMSPHEREPOINT
    IMPLICIT NONE
    TYPE(DYNNETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: NSTEP, STARTTIME
    INTEGER, INTENT(INOUT) :: CLOCKRATE
    DOUBLE PRECISION, INTENT(IN) :: DELT
    TYPE(DYNPARAM), POINTER :: DP

    DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793D0

    INTEGER :: STEP, NC, DEG1, N2, DEG2, NEIGHBS1(4),NEIGHBS2(4), EC, DC, TOP, INRANGE
    DOUBLE PRECISION :: ENERGY, ENERGY0, FORCES(NETP%DIM,NETP%NNODE), MASSLEN
    DOUBLE PRECISION :: INFO(1), CURTIME, EDGEVALS(NETP%NEDGE,1+NETP%NSPECIES)
    INTEGER, DIMENSION(7,NETP%NNODE*(DP%SNAPSHOTEVERY+1)) :: EVENTS !type, time, n1, n2, e1, e2, e3
    DOUBLE PRECISION, DIMENSION(NETP%DIM,NETP%NNODE*(DP%SNAPSHOTEVERY+1)) :: EVENTLOCS
    DOUBLE PRECISION, DIMENSION(4,NETP%NNODE*(DP%SNAPSHOTEVERY+1)) :: EVENTRATES
    DOUBLE PRECISION :: NODEVALS(NETP%NNODE,1)
    LOGICAL :: DIDFISS, CLOSENUFF, ALLOWFUSE
    DOUBLE PRECISION :: FISSPROB2, FISSPROB3, RECHARGEPROB
    DOUBLE PRECISION :: YEASTOFFPROB, YEASTONPROB
    LOGICAL :: FOUND, CANFUSENOW(NETP%NNODE), FISSABLENODE(NETP%NNODE)
    INTEGER :: CC, ACTIVE, CURRCLOCKTIME
    INTEGER :: SPECIESPRODUCER(NETP%NSPECIES,NETP%NUMTRIGGERUNITS), PC
    INTEGER :: LEFT(NETP%DIM), RIGHT(NETP%DIM), NINCURR, GC1, GC2, GC3, IDX, &
    & GRIDNODES(NETP%NNODE+1,PRODUCT(NETP%GRIDSPACES))
    DOUBLE PRECISION :: MINPOS(NETP%DIM), MAXPOS(NETP%DIM), SHIFTED(NETP%DIM), PARTSZ(NETP%DIM), DELTA
    INTEGER :: SUBSTEP

    EVENTRATES = 0D0
    
    CALL GETMECHFORCES(NETP,ENERGY0,FORCES,.FALSE.)
    PRINT*, 'STEP, ENERGY:', 0, ENERGY0

    EDGEVALS(:,1) = NETP%EDGELEN
    NODEVALS = 0D0
    EDGEVALS(:,2:NETP%NSPECIES+1) = NETP%CONC(:,1:NETP%NSPECIES)

    TOP = 0; EVENTS = 0; EVENTLOCS = 0D0

    NEIGHBS1 = 0D0; NEIGHBS2 = 0D0; CLOSENUFF = .FALSE.; ALLOWFUSE = .FALSE.

    ! Initial snapshot
    INFO = 0D0
    IF(NETP%DOYEAST) THEN
       CALL DYNNETSNAPSHOT(NETP,DP%SNAPSHOTFILE,DP%APPENDSNAPSHOTS,INFO,1,NODEVALS,1+NETP%NSPECIES,EDGEVALS)
    ELSE
       CALL DYNNETSNAPSHOT(NETP,DP%SNAPSHOTFILE,DP%APPENDSNAPSHOTS,INFO,0,NODEVALS,1+NETP%NSPECIES,EDGEVALS)
    ENDIF
    CALL REMODELINGSNAPSHOT(EVENTS(:,1:TOP),EVENTLOCS(:,1:TOP),EVENTRATES(:,1:TOP),0,DP%REMODELINGFILE,.FALSE.,NETP%DIM)

    CURTIME = 0D0

    ! fission probability for each step for each node
    FISSPROB2 = 1-EXP(-DELT*NETP%FISSRATE)
    FISSPROB3 = 1-EXP(-DELT*NETP%D3FISSMULT*NETP%FISSRATE)
    ! probability to return to a fuse-active state at each step following fusion
    RECHARGEPROB = 1-EXP(-DELT*NETP%RECHARGERATE)
    ! probabilities to tether to the cell boundary for yeast Model
    YEASTONPROB = 1-EXP(-DELT*NETP%YEASTONRATE)
    YEASTOFFPROB = 1-EXP(-DELT*NETP%YEASTOFFRATE)

    ! set the producer/fixed conc edges
    IF(NETP%FIXEDGES.OR.ANY(NETP%PRODRATE.GT.0D0)) THEN
      ACTIVE = COUNT(NETP%EDGEACT)
      SPECIESPRODUCER = 0
      IF(NETP%UNIQUETRIGGERS) THEN ! use different trigger units for each species
         DO CC = 1,NETP%NSPECIES
            DO PC = 1,NETP%NUMTRIGGERUNITS 
               FOUND = .FALSE.
               DO WHILE (.NOT.FOUND)
                  EC = CEILING(GRND()*ACTIVE)
                  IF(NETP%EDGEACT(EC).AND.ALL(SPECIESPRODUCER(CC,:).NE.EC)) THEN
                     SPECIESPRODUCER(CC,PC) = EC
                     IF(NETP%FIXEDGES) NETP%CONC(EC,CC) = 1D0
                     FOUND = .TRUE.
                     PRINT*,"Added mito ____ as trigger number ___ for species ___:",EC,PC,CC
                  ENDIF
               ENDDO
            ENDDO
            PRINT*,"TOTAL CONCENTRATION IN SPECIES ___ IS ___",CC,SUM(NETP%CONC(:,CC))
         ENDDO
      ELSE ! use the same trigger units for each species
         DO PC = 1,NETP%NUMTRIGGERUNITS
            FOUND = .FALSE.
            SEARCH: DO WHILE (.NOT.FOUND)
               EC = CEILING(GRND()*ACTIVE)
               IF(NETP%EDGEACT(EC).AND.ALL(SPECIESPRODUCER(1,:).NE.EC)) THEN
                  SPECIESPRODUCER(:,PC) = EC
                  IF(NETP%FIXEDGES) NETP%CONC(EC,:) = 1D0
                  FOUND = .TRUE.
                  PRINT*,"Added mito ____ as trigger number:",EC,PC
               ENDIF
            ENDDO SEARCH
         ENDDO
         PRINT*,"TOTAL CONCENTRATION IN EACH SPECIES ___ IS ___",SUM(NETP%CONC,1)
         PRINT*,"SPECIESPRODUCER ARRAY",SPECIESPRODUCER
       ENDIF
    ENDIF

    ! main time loop for the simulation
    DO STEP = 1,NSTEP
       CURTIME = CURTIME+DELT

       !undergo diffusion & movement from mechanical forces
       CALL LANGEVINSTEPRK4(NETP,DELT,ENERGY,DP%DOBROWN)

       !if there are diffusing species, propagate them along the network
       IF(NETP%NSPECIES.GT.0.AND.STEP.GE.DP%STARTDIFFSTEP) THEN
          ! take additional timesteps to ensure convergence of FVM solution for diffusion on the edges
          DO SUBSTEP = 1,DP%NSUBSTEPS
             CALL DIFFUSEONEDGES(NETP,DELT/DP%NSUBSTEPS,SPECIESPRODUCER)
             CALL PRODUCEANDDECAYSPECIES(NETP,DELT/DP%NSUBSTEPS,SPECIESPRODUCER)
          ENDDO
       ENDIF

       ! if we are doing the yeast model, loop thru the nodes to check if they should tether/untether to the cell membrane
       IF(NETP%DOYEAST) THEN
          ACTIONLOOP: DO NC = 1,NETP%NNODE
             IF (.NOT.NETP%NODEACT(NC)) CYCLE
             IF(.NOT.NETP%TETHERED(NC)) THEN
                ! if not tethered, check if in range and try to tether
                INRANGE = 0
                IF(SUM(NETP%NODEPOS(:,NC)**2).GT.(NETP%CELLRAD1-NETP%YEASTBINDRANGE)**2) INRANGE = 1
                IF(INRANGE.GT.0.AND.GRND().LT.YEASTONPROB) THEN
                   NETP%TETHERED(NC) = .TRUE.
                   IF(VERBOSE) PRINT*,'Getting on wall',NC,STEP
                ENDIF
             ELSEIF(GRND().LT.YEASTOFFPROB) THEN
                NETP%TETHERED(NC) = .FALSE.
                IF(VERBOSE) PRINT*,'Getting off wall',NC,STEP
             ENDIF
          ENDDO ACTIONLOOP
       ENDIF

       IF(NETP%USENODEGRID) THEN ! get some necessary info
          DO DC = 1,NETP%DIM
             MINPOS(DC) = MINVAL(NETP%NODEPOS(DC,:),NETP%NODEACT)
             MAXPOS(DC) = MAXVAL(NETP%NODEPOS(DC,:),NETP%NODEACT)
             PARTSZ(DC) = MAX((MAXPOS(DC)-MINPOS(DC))/NETP%GRIDSPACES(DC),2*NETP%CONTACTRAD)
          ENDDO
          GRIDNODES = 0
       ELSE
          CANFUSENOW = .FALSE.
       ENDIF

       FISSABLENODE = NETP%NODEACT ! don't try to fiss new nodes which were created by a fission event
       ! go thru the nodes and try to fiss them. Also, figure out which ones should be added to the fusable nodes
       FISSIONLOOP: DO NC = 1,NETP%NNODE
          ! decide whether to fission for each node
          ! only do fission if degree > 1
          ! if using node grid set that up
          IF (.NOT.NETP%NODEACT(NC)) CYCLE FISSIONLOOP
          IF(.NOT.FISSABLENODE(NC)) CYCLE FISSIONLOOP

          IF(NETP%RECHARGERATE.GT.0D0.AND.(.NOT.NETP%FUSEACTIVE(NC))) THEN
             IF(GRND().LT.RECHARGEPROB) THEN
                NETP%FUSEACTIVE(NC) = .TRUE.
             ENDIF
          ENDIF
          
          DIDFISS = .FALSE.
          DEG1 = NETP%NODEDEG(NC)
          IF(((DEG1.EQ.2).AND.(GRND().LT.FISSPROB2)).OR.((DEG1.EQ.3).AND.(GRND().LT.FISSPROB3))) THEN
             ! carry out fission
             CALL FISSIONEVENT(NETP,NC,DEG1,DIDFISS,N2)
             IF (VERBOSE) PRINT*,'Fissing Nodes at step: ',NC,N2,STEP

             IF(DIDFISS) THEN
                ! move nodes to the inactive state following fission
                IF(NETP%RECHARGERATE.GT.0D0) THEN
                   NETP%FUSEACTIVE(NC) = .FALSE.
                   NETP%FUSEACTIVE(N2) = .FALSE.
                ENDIF
                ! record some information about the fission event
                TOP=TOP+1
                IF(DEG1.EQ.2) THEN
                   ! both nodes are now degree 1, write each of their edges
                   EVENTS(:,TOP) = (/ DEG1, NC, N2, STEP, NETP%NODEEDGE(NC,1), NETP%NODEEDGE(N2,1), -1 /)
                ELSEIF(DEG1.EQ.3) THEN
                   ! fissionevent always turns nc into degree 2, write its edges first
                   EVENTS(:,TOP) = (/ DEG1, NC, N2, STEP, NETP%NODEEDGE(NC,1), NETP%NODEEDGE(NC,2), NETP%NODEEDGE(N2,1) /)
                ENDIF
                EVENTLOCS(:,TOP) = NETP%NODEPOS(:,NC)
                
                ! Keep track of rates associated with this fission event
                IF (DEG1.EQ.2) THEN
                   EVENTRATES(1,TOP) = NETP%FISSRATE
                ELSEIF (DEG1.EQ.3) THEN
                   EVENTRATES(1,TOP) = NETP%D3FISSMULT*NETP%FISSRATE
                ENDIF

                ! calculate a fusion rate for right after this fission
                CALL GETFUSIONRATE(NETP,NC,N2,DEG1,EVENTRATES(2,TOP))
             ENDIF
          ENDIF

          IF(DEG1.LT.3) THEN ! Determine which nodes can fuse. Only degrees < 3 can possibly fuse
            ! don't allow multiple fission/fusion events per node per step
            ! don't allow nodes in a fuse-inactive state to fuse
            IF(DIDFISS.OR.(NETP%RECHARGERATE.GT.0D0.AND.(.NOT.NETP%FUSEACTIVE(NC)))) CYCLE FISSIONLOOP

             IF(NETP%USENODEGRID) THEN ! set up the node grid for fusion.
                SHIFTED = NETP%NODEPOS(:,NC)-MINPOS
                LEFT = 0; RIGHT = 0
                DO DC = 1,NETP%DIM
                   LEFT(DC) = MIN(MAX(CEILING(SHIFTED(DC)/PARTSZ(DC)),1),NETP%GRIDSPACES(DC))
                   RIGHT(DC) = LEFT(DC)
                   DELTA = MOD(SHIFTED(DC),PARTSZ(DC))
                   IF(LEFT(DC).GT.1.AND.DELTA.LT.NETP%CONTACTRAD) THEN
                      LEFT(DC) = LEFT(DC) -1
                   ELSEIF(RIGHT(DC).LT.NETP%GRIDSPACES(DC).AND. &
                   & DELTA.GT.(PARTSZ(DC)-NETP%CONTACTRAD)) THEN
                      RIGHT(DC) = RIGHT(DC) +1
                   ENDIF
                ENDDO
                DO GC1 = LEFT(1),RIGHT(1)
                   DO GC2 = LEFT(2),RIGHT(2)
                      DO GC3 = LEFT(3),RIGHT(3)
                         IDX = 1 + (GC3-1) + (GC2-1)*NETP%GRIDSPACES(3) + (GC1-1)*NETP%GRIDSPACES(3)*NETP%GRIDSPACES(2)
                         NINCURR = GRIDNODES(1,IDX) + 1
                         GRIDNODES(1,IDX) = NINCURR
                         GRIDNODES(NINCURR+1,IDX) = NC
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                CANFUSENOW(NC) = .TRUE.
             ENDIF
          ENDIF

       ENDDO FISSIONLOOP

       ! do the fusion based on a spatial grid. It will speed up the code for large networks
       ! Assume the nodegrid was set above
       IF(NETP%USENODEGRID) THEN
          CALL MANAGEFUSIONSWITHGRID(NETP,GRIDNODES,EVENTS(:,TOP+1:TOP+NETP%NNODE),&
               & EVENTLOCS(:,TOP+1:TOP+NETP%NNODE),EVENTRATES(:,TOP+1:TOP+NETP%NNODE),TOP,STEP,DELT)
       ELSE ! don't use node grid     
          CALL MANAGEFUSIONSBASIC(NETP,EVENTS(:,TOP+1:TOP+NETP%NNODE),EVENTLOCS(:,TOP+1:TOP+NETP%NNODE),&
               & EVENTRATES(:,TOP+1:TOP+NETP%NNODE),TOP,CANFUSENOW,STEP,DELT)
       ENDIF

       IF (MOD(STEP,DP%PRINTEVERY).EQ.0) THEN
          ! Output information on chain status
          MASSLEN = 0D0
          DO EC = 1,NETP%NEDGE
             IF(.NOT.NETP%EDGEACT(EC)) CYCLE
             DEG1 = NETP%NODEDEG(NETP%EDGENODE(EC,1))
             DEG2 = NETP%NODEDEG(NETP%EDGENODE(EC,2))
             MASSLEN = MASSLEN + NETP%EDGELEN(EC)
             IF(DEG1.EQ.1) MASSLEN = MASSLEN + NETP%STERICRAD
             IF(DEG2.EQ.1) MASSLEN = MASSLEN + NETP%STERICRAD
          ENDDO
          CALL SYSTEM_CLOCK(CURRCLOCKTIME,CLOCKRATE)
          PRINT*, 'STEP, ENERGY, MASS, ELAPSEDTIME:', STEP, ENERGY, MASSLEN,(CURRCLOCKTIME-STARTTIME)/1000
       ENDIF

       IF (MOD(STEP,DP%SNAPSHOTEVERY).EQ.0) THEN
          IF(STEP.GE.DP%STARTSNAPSTEP-1) THEN
             ! output snapshot
             INFO = CURTIME
             EDGEVALS(:,1) = NETP%EDGELEN
             EDGEVALS(:,2:NETP%NSPECIES+1) = NETP%CONC(:,1:NETP%NSPECIES)
             IF(NETP%DOYEAST) THEN
                NODEVALS(:,1) = 0D0
                DO NC = 1,NETP%NNODE
                   IF(.NOT.NETP%NODEACT(NC)) CYCLE
                   IF(NETP%TETHERED(NC)) NODEVALS(NC,1) = 1D0
                ENDDO
                CALL DYNNETSNAPSHOT(NETP,DP%SNAPSHOTFILE,.TRUE.,INFO,1,NODEVALS,1+NETP%NSPECIES,EDGEVALS)
             ELSE
                CALL DYNNETSNAPSHOT(NETP,DP%SNAPSHOTFILE,.TRUE.,INFO,0,NODEVALS,1+NETP%NSPECIES,EDGEVALS)
             ENDIF
          ENDIF
          ! output fusion/fission events
          CALL REMODELINGSNAPSHOT(EVENTS(:,1:TOP),EVENTLOCS(:,1:TOP),EVENTRATES(:,1:TOP),TOP,DP%REMODELINGFILE,.TRUE.,NETP%DIM)
          EVENTS = 0; EVENTLOCS = 0D0; EVENTRATES = 0D0; TOP = 0
       ENDIF
    ENDDO

    PRINT*,'Active Nodes, Active Edges',SIZE(PACK(NETP%NODEACT,NETP%NODEACT)),SIZE(PACK(NETP%EDGEACT,NETP%EDGEACT))

  END SUBROUTINE RUNBROWNDYNSIM


  SUBROUTINE MANAGEFUSIONSBASIC(NETP,EVENTS,EVENTLOCS,EVENTRATES,TOP,CANFUSENOW,STEP,DELT)
    USE DYNNETWORKUTIL, ONLY : dynNETWORK, GETNEIGHBNODES
 
    IMPLICIT NONE
    TYPE(DYNNETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: STEP
    DOUBLE PRECISION, INTENT(IN) :: DELT
    LOGICAL, INTENT(IN) :: CANFUSENOW(NETP%NNODE)
    INTEGER, INTENT(INOUT) :: EVENTS(7,NETP%NNODE), TOP
    DOUBLE PRECISION, INTENT(INOUT) :: EVENTLOCS(NETP%DIM,NETP%NNODE),EVENTRATES(4,NETP%NNODE)
    LOGICAL :: DIDFUSENODES(NETP%NNODE), DIDFUSE, CLOSENUFF
    INTEGER :: NC, N2, DEG1
    INTEGER :: NEIGHBS1(3), NEIGHBS2(3), DEG2, NCN, DC, TOPOFFSET
 
    DIDFUSENODES = .FALSE.; TOPOFFSET = 0;
 
    N1LOOP: DO NC = 1,NETP%NNODE
       IF(.NOT.NETP%NODEACT(NC)) CYCLE N1LOOP
       IF((.NOT.CANFUSENOW(NC)).OR.DIDFUSENODES(NC)) CYCLE N1LOOP
       DEG1 = NETP%NODEDEG(NC)
       DIDFUSE = .FALSE.
       CALL GETNEIGHBNODES(NETP,NC,NEIGHBS1(1:DEG1))
       N2LOOP: DO N2 = NC+1,NETP%NNODE
          IF (.NOT.NETP%NODEACT(N2)) CYCLE N2LOOP
          IF((.NOT.CANFUSENOW(N2)).OR.DIDFUSENODES(N2)) CYCLE N2LOOP
          IF(ANY(NEIGHBS1(1:DEG1).EQ.N2)) CYCLE N2LOOP
          CLOSENUFF = .FALSE.
          IF(NETP%USEGRID) THEN
             DO DC = 1,NETP%DIM
                IF(NETP%GRID(DC,N2).LT.(NETP%GRID(DC,NC)-1).OR.NETP%GRID(DC,N2).GT.(NETP%GRID(DC,NC)+1)) THEN
                   CLOSENUFF = .FALSE.
                   CYCLE N2LOOP
                ELSE
                   CLOSENUFF = .TRUE.
                ENDIF
             ENDDO
          ENDIF
          IF(CLOSENUFF.OR.(.NOT.NETP%USEGRID)) THEN
             DEG2 = NETP%NODEDEG(N2)
             IF((DEG2+DEG1).GT.3) CYCLE N2LOOP
             CALL GETNEIGHBNODES(NETP,N2,NEIGHBS2(1:DEG2))
             ! shouldn't let neighbors or neighbor-neighbors fuse
             IF(ANY(NEIGHBS2(1:DEG2).EQ.NC)) CYCLE N2LOOP
             DO NCN = 1,DEG1
                IF(ANY(NEIGHBS2(1:DEG2).EQ.NEIGHBS1(NCN))) CYCLE N2LOOP
             ENDDO
             CALL TRYTOFUSE(NETP,NC,N2,DEG1,DEG2,STEP,DELT,DIDFUSE,EVENTS,EVENTLOCS,&
                  & EVENTRATES,TOP,TOPOFFSET,DIDFUSENODES)
             IF(DIDFUSE) EXIT N2LOOP
           ENDIF
       ENDDO N2LOOP
    ENDDO N1LOOP

  END SUBROUTINE MANAGEFUSIONSBASIC


  SUBROUTINE MANAGEFUSIONSWITHGRID(NETP,GRIDNODES,EVENTS,EVENTLOCS,EVENTRATES,TOP,STEP,DELT)
    USE DYNNETWORKUTIL, ONLY : dynNETWORK, GETNEIGHBNODES

    IMPLICIT NONE
    TYPE(DYNNETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: GRIDNODES(NETP%NNODE+1,PRODUCT(NETP%GRIDSPACES)),STEP
    DOUBLE PRECISION, INTENT(IN) :: DELT
    INTEGER, INTENT(INOUT) :: EVENTS(7,NETP%NNODE), TOP
    DOUBLE PRECISION, INTENT(INOUT) :: EVENTLOCS(NETP%DIM,NETP%NNODE), EVENTRATES(4,NETP%NNODE)
    LOGICAL :: INTERACTIONMATRIX(NETP%NNODE,NETP%NNODE), DIDFUSENODES(NETP%NNODE), DIDFUSE
    INTEGER :: GC, NINCURR, GNC1, NC, GNC2, N2, DEG1
    INTEGER :: NEIGHBS1(3), NEIGHBS2(3), DEG2, NCN, TOPOFFSET

    INTERACTIONMATRIX = .FALSE. ! keeps track of which node pairs have been checked (to avoid repeats)
    DIDFUSENODES = .FALSE.; TOPOFFSET = 0; NINCURR = 0
    DO GC = 1,PRODUCT(NETP%GRIDSPACES)
       NINCURR = GRIDNODES(1,GC) ! how many nodes are inside the current grid cube?
       IF(NINCURR.LT.2) CYCLE
       N1LOOP: DO GNC1 = 2,NINCURR+1
          NC = GRIDNODES(GNC1,GC)
          IF(DIDFUSENODES(NC)) CYCLE N1LOOP
          DEG1 = NETP%NODEDEG(NC)
          CALL GETNEIGHBNODES(NETP,NC,NEIGHBS1(1:DEG1))
          DIDFUSE = .FALSE.
          N2LOOP: DO GNC2 = GNC1+1,NINCURR+1
             N2 = GRIDNODES(GNC2,GC)
             ! make sure the current pair has not already been fused or checked
             IF((.NOT.DIDFUSENODES(N2)).AND.(NC.NE.N2).AND.(.NOT.INTERACTIONMATRIX(NC,N2)).AND. &
             & (.NOT.INTERACTIONMATRIX(N2,NC))) THEN
               ! remember that we have checked this pair
                INTERACTIONMATRIX(NC,N2) = .TRUE.; INTERACTIONMATRIX(N2,NC) = .TRUE.
                ! make sure the pair are not neighbors
                IF(.NOT.(ANY(NEIGHBS1(1:DEG1).EQ.N2))) THEN 
                   DEG2 = NETP%NODEDEG(N2)
                   IF(DEG1+DEG2.LT.4) THEN
                      CALL GETNEIGHBNODES(NETP,N2,NEIGHBS2(1:DEG2))
                      ! make sure the pair do not share a neighbor (we don't want duplicate edges)
                      IF(.NOT.(ANY(NEIGHBS2(1:DEG2).EQ.NC))) THEN
                         DO NCN = 1,DEG1
                            IF(ANY(NEIGHBS2(1:DEG2).EQ.NEIGHBS1(NCN))) CYCLE N2LOOP
                         ENDDO
                         CALL TRYTOFUSE(NETP,NC,N2,DEG1,DEG2,STEP,DELT,DIDFUSE,EVENTS,EVENTLOCS,EVENTRATES,TOP,TOPOFFSET,&
                         & DIDFUSENODES)
                         IF(DIDFUSE) EXIT N2LOOP
                      ENDIF
                   ENDIF
                ENDIF
             ENDIF
          ENDDO N2LOOP
       ENDDO N1LOOP
    ENDDO
  
  END SUBROUTINE MANAGEFUSIONSWITHGRID

  SUBROUTINE TRYTOFUSE(NETP,NC,N2,DEG1,DEG2,STEP,DELT,DIDFUSE,EVENTS,EVENTLOCS,EVENTRATES,&
      & TOP,TOPOFFSET,DIDFUSENODES)
    USE DYNNETWORKUTIL, ONLY : dynNETWORK, GETMECHFORCES
    USE DYNAMICPROCESSES, ONLY : GETFUSIONRATE, FUSIONEVENT
    USE MT19937, ONLY : GRND

    IMPLICIT NONE
    TYPE(DYNNETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: NC,N2,DEG1,DEG2,STEP
    DOUBLE PRECISION, INTENT(IN) :: DELT
    INTEGER, INTENT(INOUT) :: EVENTS(7,NETP%NNODE), TOP, TOPOFFSET
    LOGICAL, INTENT(INOUT) :: DIDFUSENODES(NETP%NNODE)  
    DOUBLE PRECISION, INTENT(INOUT) :: EVENTLOCS(NETP%DIM,NETP%NNODE), EVENTRATES(4,NETP%NNODE)
    LOGICAL, INTENT(OUT) :: DIDFUSE
    DOUBLE PRECISION :: FUSEPROB, DISTSQR, FUSERATE
    
    DIDFUSE = .FALSE.
    DISTSQR = SUM((NETP%NODEPOS(:,NC)-NETP%NODEPOS(:,N2))**2)
    ! possible fusion only if inside contact range
    IF(DISTSQR.LT.(2*NETP%CONTACTRAD)**2) THEN
       FUSEPROB = 0D0; FUSERATE = 0D0
       CALL GETFUSIONRATE(NETP,NC,N2,DEG1+DEG2,FUSERATE)
       FUSEPROB = 1D0 - EXP(-DELT*FUSERATE)
       IF(GRND().LT.FUSEPROB) THEN
          ! carry out fusion
          IF (VERBOSE) PRINT*,'Fusing nodes at step: ', NC,N2,STEP
          !NC will die for deg 2 fusion, the deg 1 node dies for deg 3 fusion
          CALL FUSIONEVENT(NETP,NC,N2,DEG2+DEG1,DIDFUSE)
          IF(DIDFUSE) THEN
             DIDFUSENODES(NC) = .TRUE.; DIDFUSENODES(N2) = .TRUE.
             IF(NETP%DOYEAST) THEN
                IF(NETP%TETHERED(NC).OR.NETP%TETHERED(N2)) THEN ! if either of fusing nodes was tethered, keep new one tethered
                   NETP%TETHERED(NC) = .TRUE.
                   NETP%TETHERED(N2) = .TRUE.
                ENDIF
                IF(NETP%NODEACT(NC)) THEN
                   NETP%TETHERED(N2) = .FALSE.
                ELSE
                   NETP%TETHERED(NC) = .FALSE.
                ENDIF
             ENDIF
             TOP=TOP+1; TOPOFFSET = TOPOFFSET+1
             IF(NETP%NODEACT(NC)) THEN
                EVENTS(:,TOPOFFSET) = (/ -(DEG1+DEG2), NC, N2, STEP, NETP%NODEEDGE(NC,1), NETP%NODEEDGE(NC,2), -1 /)
                IF(DEG1+DEG2.EQ.3) THEN
                   EVENTS(7,TOPOFFSET) = NETP%NODEEDGE(NC,3)
                ENDIF
                EVENTLOCS(:,TOPOFFSET) = NETP%NODEPOS(:,NC)                                        
             ELSE
                EVENTS(:,TOPOFFSET) = (/ -(DEG1+DEG2), N2, NC, STEP, NETP%NODEEDGE(N2,1), NETP%NODEEDGE(N2,2), -1 /)
                IF(DEG1+DEG2.EQ.3) THEN
                   EVENTS(7,TOPOFFSET) = NETP%NODEEDGE(N2,3)
                ENDIF
                EVENTLOCS(:,TOPOFFSET) = NETP%NODEPOS(:,N2)
             ENDIF
             ! Rate for this fusion event to have occurred
             EVENTRATES(2,TOPOFFSET)= FUSERATE
             ! rate of subsequent fission
             IF (DEG1+DEG2.EQ.2) THEN
                EVENTRATES(1,TOPOFFSET) = NETP%FISSRATE
             ELSE IF (DEG1+DEG2.EQ.3) THEN
                EVENTRATES(1,TOPOFFSET) = NETP%D3FISSMULT * NETP%FISSRATE
             ENDIF
          ELSE
             PRINT*,"DIDN't fuse when it said it would"
             STOP 1
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE TRYTOFUSE


  SUBROUTINE LANGEVINSTEPRK4(NETP,DELT,ENERGY,DOBROWN)
    ! propagate network forward in time, using a fourth-order Runge-Kutta method
    ! DELT: timestep
    ! DOBROWN: toggle whether to include brownian forces
    USE MT19937, ONLY : RNORM, GRND
    USE DYNNETWORKUTIL, ONLY : dynNETWORK, GETMECHFORCES

    IMPLICIT NONE
    TYPE(DYNNETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(IN) :: DELT
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY
    LOGICAL, INTENT(IN) :: DOBROWN
    DOUBLE PRECISION :: S2DT
    DOUBLE PRECISION, DIMENSION(NETP%DIM,NETP%NNODE) :: POS0, K1POS, K2POS, K3POS, &
         & K4POS, FORCES, FBROWN
    INTEGER :: B, I

    ENERGY = 0D0; FORCES = 0D0

    S2DT = SQRT(2*NETP%KT*NETP%FRICT/DELT)

    POS0 = NETP%NODEPOS

    ! get the brownian forces
    IF (DOBROWN) THEN
       DO B = 1,NETP%NNODE
          DO I = 1,3
             FBROWN(I,B) = RNORM()*S2DT
          ENDDO
       END DO
    ELSE
       FBROWN = 0D0
    ENDIF

    ! --------- 1ST RK STEP---------------
    CALL GETMECHFORCES(NETP,ENERGY,FORCES,.TRUE.)
    K1POS = (FBROWN + FORCES)/NETP%FRICT

    NETP%NODEPOS = POS0 + DELT/2*K1POS

    ! --------- 2ND RK STEP---------------
    CALL GETMECHFORCES(NETP,ENERGY,FORCES,.TRUE.)
    K2POS = (FBROWN + FORCES)/NETP%FRICT

    NETP%NODEPOS = POS0 + DELT/2*K2POS

    ! --------- 3RD RK STEP---------------
    CALL GETMECHFORCES(NETP,ENERGY,FORCES,.TRUE.)
    K3POS = (FBROWN + FORCES)/NETP%FRICT

    NETP%NODEPOS = POS0 + DELT*K3POS

    ! --------- 4TH RK STEP---------------
    CALL GETMECHFORCES(NETP,ENERGY,FORCES,.TRUE.)
    K4POS = (FBROWN + FORCES)/NETP%FRICT

    NETP%NODEPOS = POS0 + DELT/6*(K1POS + 2*K2POS + 2*K3POS + K4POS)

    ! calculate final energy and set the grid correctly
    CALL GETMECHFORCES(NETP,ENERGY,FORCES,.FALSE.)
  END SUBROUTINE LANGEVINSTEPRK4

  SUBROUTINE SETDYNPARAM(DP)
    ! set parameters from keyword global arguments
    USE KEYS, ONLY : OUTFILE, SNAPSHOTFILE, PRINTEVERY, SNAPSHOTEVERY, &
         &    APPENDSNAPSHOTS, DOBROWN, REMODELINGFILE, STARTSNAPSTEP,STARTDIFFSTEP

    IMPLICIT NONE
    TYPE(DYNPARAM), POINTER :: DP

    DP%OUTFILE = OUTFILE
    DP%SNAPSHOTFILE = SNAPSHOTFILE
    DP%REMODELINGFILE = REMODELINGFILE
    DP%PRINTEVERY = PRINTEVERY
    DP%SNAPSHOTEVERY = SNAPSHOTEVERY
    DP%APPENDSNAPSHOTS = APPENDSNAPSHOTS
    DP%DOBROWN = DOBROWN
    DP%STARTSNAPSTEP = STARTSNAPSTEP
    DP%STARTDIFFSTEP = STARTDIFFSTEP

  END SUBROUTINE SETDYNPARAM

  SUBROUTINE REMODELINGSNAPSHOT(EVENTS,EVENTLOCS,EVENTRATES,NEVENTS,REMODELINGFILE,APPEND,DIM)
    ! output a snapshot of the remodeling events
    ! in an efficient format (few lines per network)
    ! APPEND: append to an existing file?
    ! TOTTIME: how many total steps in the simulation
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NEVENTS
    INTEGER, INTENT(IN) :: EVENTS(7,NEVENTS)
    DOUBLE PRECISION, INTENT(IN) :: EVENTLOCS(DIM,NEVENTS), EVENTRATES(4,NEVENTS)
    CHARACTER (LEN=*), INTENT(IN) :: REMODELINGFILE
    LOGICAL, INTENT(IN) :: APPEND
    INTEGER, INTENT(IN) :: DIM
    INTEGER, PARAMETER :: FU=51 ! File output unit
    CHARACTER (LEN=100) :: FMTINFO
    INTEGER :: IC

    IF (APPEND) THEN
      OPEN(UNIT=FU,FILE=REMODELINGFILE,STATUS='UNKNOWN',ACCESS='APPEND')
      DO IC = 1,NEVENTS
         !PRINT*,'WRITING: ',EVENTS(1,IC),EVENTS(2,IC),EVENTS(3,IC),EVENTS(4,IC)

         ! this is a cleaner way to deal with multiple dimensions
         WRITE(FMTINFO,'(A,I1,A)') '(I3, I5, I5, I12, I5, I5, I5, ', DIM, 'F20.10, 4F20.10)'
         WRITE(FU,FMTINFO) EVENTS(1,IC),EVENTS(2,IC),EVENTS(3,IC),EVENTS(4,IC),&
                  &  EVENTS(5,IC),EVENTS(6,IC),EVENTS(7,IC), EVENTLOCS(1:DIM,IC), EVENTRATES(:,IC)
       
      ENDDO
    ELSE
      OPEN(UNIT=FU,FILE=REMODELINGFILE,STATUS='UNKNOWN')
      WRITE(FU,'(A)') '#IGNORE THIS LINE'
   ENDIF

   ! info line format
   !WRITE(FMTINFO,'(I8)')
   ! write the number of events
   !WRITE(FU,'(I2)') NEVENTS

   ! write each event. one line per event
   !WRITE(EFMT,'(A)') '(I8 I8 I8 F4.3)'


   CLOSE(FU)
  END SUBROUTINE REMODELINGSNAPSHOT
END MODULE BROWNDYN