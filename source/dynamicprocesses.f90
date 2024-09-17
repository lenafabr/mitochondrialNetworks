MODULE DYNAMICPROCESSES
  ! subroutines for dynamic processes on the network
  USE STACKUTIL, ONLY : STACK, PUSH, POP
  USE DYNNETWORKUTIL, ONLY : DYNNETWORK, GETNEIGHBNODES
  USE GENUTIL, ONLY : NORM, CROSS_PRODUCT

  IMPLICIT NONE  
CONTAINS
   SUBROUTINE FISSIONEVENT(NETP,NC,DEG,DIDFISS,NNEW)
    ! carry out fission event at node index NC
    ! do nothing if node degree is 1
    USE MT19937, ONLY : GRND

    IMPLICIT NONE
    TYPE(DYNNETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: NC, DEG
    LOGICAL, INTENT(OUT) :: DIDFISS
    INTEGER, INTENT(OUT) :: NNEW
    INTEGER :: EC2, EC1, EC1I, EC2I, NC2, ECC, ECC2, EC
    DOUBLE PRECISION, DIMENSION(NETP%DIM) :: DR1(NETP%DIM), DR2(NETP%DIM)

    IF (.NOT.NETP%NODEACT(NC)) THEN
       PRINT*, 'ERROR IN FISSIONEVENT: cannot fiss at inactive node', NC
       STOP 1
    ENDIF

    DIDFISS = .FALSE.; NNEW = -1

    IF (DEG.EQ.1) THEN
       ! degree 1 node. Do nothing
       ! the way the code is currently set up this should never happen
       PRINT*,"ERROR  in FISSIONEVENT, cannot fiss a degree-1 node", NC
       STOP 1
    ELSEIF (DEG.EQ.2) THEN
       ! degree 2 node. Separate the two adjacent edges

       ! which edge to break off and attach to new node
       ECC2 = 2 ! index among edges for this node
       EC2 = NETP%NODEEDGE(NC,2)
       EC2I = NETP%NODEEDGEIND(NC,2) ! For this edge, it's the 1st or 2nd node?
       ! Need to shorten the other edge as well since it is becoming a degree 1 branch
       EC1 = NETP%NODEEDGE(NC,1)
       EC1I = NETP%NODEEDGEIND(NC,1)
       NETP%EDGELEN(EC1) = NETP%EDGELEN(EC1) - NETP%STERICRAD

       !DR2 points from new node to the other node at end of its edge, DR1 points from old node to remaining neighbor node
       DR2 = NETP%NODEPOS(:,NETP%EDGENODE(EC2,MOD(EC2I,2)+1)) - NETP%NODEPOS(:,NC)
       DR2 = DR2/NORM(DR2)
       DR1 = NETP%NODEPOS(:,NETP%EDGENODE(EC1,MOD(EC1I,2)+1)) - NETP%NODEPOS(:,NC)
       DR1 = DR1/NORM(DR1)

    ELSEIF (DEG.EQ.3) THEN
       ECC2 = CEILING(3*GRND())
       EC2 = NETP%NODEEDGE(NC,ECC2)
       EC2I = NETP%NODEEDGEIND(NC,ECC2)
       DR2 = NETP%NODEPOS(:,NETP%EDGENODE(EC2,MOD(EC2I,2)+1)) - NETP%NODEPOS(:,NC)
       DR2 = DR2/NORM(DR2)
       DR1 = -DR2
    ELSE
       PRINT*, 'ERROR IN FISSIONEVENT: not currently set up to deal with degree > 3'
       STOP 1
    ENDIF

    ! Make duplicate node
    CALL POP(NETP%NODESTACK,NC2)
    NETP%NODEACT(NC2) = .TRUE.
    NETP%NODEPOS(:,NC2) = NETP%NODEPOS(:,NC)
    NNEW = NC2

    ! Broken off edge now attached to new node
    NETP%EDGENODE(EC2,EC2I) = NC2
    NETP%EDGENODEIND(EC2,EC2I) = 1
    NETP%NODEEDGE(NC2,1) = EC2
    NETP%NODEEDGEIND(NC2,1) = EC2I

    !offset the new node and change the edgelength (since it is now degree 1)
    NETP%NODEPOS(:,NC2) = NETP%NODEPOS(:,NC2) + NETP%FISSRAD*DR2
    NETP%EDGELEN(EC2) = NETP%EDGELEN(EC2) - NETP%STERICRAD
    !offset the old node as well
    NETP%NODEPOS(:,NC) = NETP%NODEPOS(:,NC) + NETP%FISSRAD*DR1

    ! shift up list of edges for the remaining old node
    DO ECC = ECC2+1,DEG
       EC = NETP%NODEEDGE(NC,ECC)
       NETP%EDGENODEIND(EC,NETP%NODEEDGEIND(NC,ECC)) = ECC-1
    ENDDO
    ! update the list of edges for the original node
    NETP%NODEEDGE(NC,ECC2:DEG-1) = NETP%NODEEDGE(NC,ECC2+1:DEG)
    NETP%NODEEDGEIND(NC,ECC2:DEG-1) = NETP%NODEEDGEIND(NC,ECC2+1:DEG)

    NETP%NODEDEG(NC) = DEG-1
    NETP%NODEDEG(NC2) = 1

    DIDFISS = .TRUE.

  END SUBROUTINE FISSIONEVENT


  SUBROUTINE FUSIONEVENT(NETP,N1,N2,DEGSUM,DIDFUSE)
    USE GENUTIL, ONLY : RAYRAYINTERSECT
    ! carry out fusion event between nodes N1, N2
    ! assuming the user only passes in nodes with degree sum = 2 or 3
    IMPLICIT NONE
    TYPE(DYNNETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: N1,N2,DEGSUM
    LOGICAL, INTENT(OUT) :: DIDFUSE
    INTEGER :: NDIE, NLIVE, EC, E1
    INTEGER :: NEIGHBS(3)
    DOUBLE PRECISION :: V1(3), V2(3)

    IF ((.NOT.NETP%NODEACT(N1)).OR.(.NOT.NETP%NODEACT(N2))) THEN
       PRINT*, 'ERROR IN FUSIONEVENT: cannot fuse at inactive node',N1,N2
       STOP 1
    ENDIF

    DIDFUSE = .FALSE.

    ! Need to modify the edgelength of degree 1 nodes, destroy one of the 2 fusing nodes.
    ! For the destroyed node, the edge that connected to it formerly should now connect
    ! to the surviving node. Also, the surviving node should connect to this edge. Need
    ! to update the degree of the surviving node, etc.

    IF(DEGSUM.EQ.2) THEN !We have 1-1 fusion
       NDIE = N1; NLIVE = N2
       EC = NETP%NODEEDGE(NDIE,1) !the broken off edge that needs to reattach
       E1 = NETP%NODEEDGE(NLIVE,1) !the intact edge. Just need to change its length

       ! make new fused point halfway between fusing tips
       CALL GETNEIGHBNODES(NETP,N1,NEIGHBS(1:1))
       CALL GETNEIGHBNODES(NETP,N2,NEIGHBS(2:2))
       V1 = NETP%NODEPOS(:,N1) - NETP%NODEPOS(:,NEIGHBS(1))
       V2 = NETP%NODEPOS(:,N2) - NETP%NODEPOS(:,NEIGHBS(2))
       NETP%NODEPOS(:,NLIVE) = 0.5*(NETP%NODEPOS(:,N1)+NETP%STERICRAD*V1/SQRT(SUM(V1**2))+ &
       & NETP%NODEPOS(:,N2)+NETP%STERICRAD*V2/SQRT(SUM(V2**2)))
       
      NETP%NODEDEG(NLIVE) = 2
      NETP%NODEEDGE(NLIVE,2) = EC
      NETP%NODEACT(NDIE) = .FALSE.
      CALL PUSH(NETP%NODESTACK,NDIE)
      NETP%EDGELEN(E1) = NETP%EDGELEN(E1) + NETP%STERICRAD
      NETP%EDGELEN(EC) = NETP%EDGELEN(EC) + NETP%STERICRAD
      ! Updating of arrays depends on whether the dying node was
      ! at the start or end of its edge...
      IF(NETP%EDGENODE(EC,1).EQ.NDIE) THEN
         NETP%EDGENODE(EC,1) = NLIVE
         NETP%NODEEDGEIND(NLIVE,2) = 1
         NETP%EDGENODEIND(EC,1) = 2
      ELSEIF(NETP%EDGENODE(EC,2).EQ.NDIE) THEN
         NETP%EDGENODE(EC,2) = NLIVE
         NETP%NODEEDGEIND(NLIVE,2) = 2
         NETP%EDGENODEIND(EC,2) = 2
      ELSE
         PRINT*,'PROBLEM IN DYNAMICPROCESSES WITH 1+1 FUSION'
      ENDIF

    ELSEIF(DEGSUM.EQ.3) THEN !We have 1-2 fusion
       IF(NETP%NODEDEG(N1).EQ.1) THEN
           NDIE = N1; NLIVE = N2
       ELSE
           NDIE = N2; NLIVE = N1
       ENDIF
       EC = NETP%NODEEDGE(NDIE,1) !the broken off edge which needs to reattach
       
       CALL GETNEIGHBNODES(NETP,NDIE,NEIGHBS(1:1))
       V1 = NETP%NODEPOS(:,NDIE) - NETP%NODEPOS(:,NEIGHBS(1))
       V1 = V1/SQRT(SUM(V1**2))
       V2 = NETP%NODEPOS(:,NDIE) + NETP%STERICRAD*V1 - NETP%NODEPOS(:,NLIVE)
       NETP%NODEPOS(:,NLIVE) = 0.5*(NETP%NODEPOS(:,NDIE) + NETP%STERICRAD*V1 + &
       & NETP%NODEPOS(:,NLIVE) + NETP%STERICRAD*V2/SQRT(SUM(V2**2)))

       NETP%NODEDEG(NLIVE) = 3
       NETP%NODEEDGE(NLIVE,3) = EC
       NETP%NODEACT(NDIE) = .FALSE.
       CALL PUSH(NETP%NODESTACK, NDIE)
       NETP%EDGELEN(EC) = NETP%EDGELEN(EC) + NETP%STERICRAD

       IF(NETP%EDGENODE(EC,1).EQ.NDIE) THEN
          NETP%EDGENODE(EC,1) = NLIVE
          NETP%NODEEDGEIND(NLIVE,3) = 1
          NETP%EDGENODEIND(EC,1) = 3
       ELSEIF(NETP%EDGENODE(EC,2).EQ.NDIE) THEN
          NETP%EDGENODE(EC,2) = NLIVE
          NETP%NODEEDGEIND(NLIVE,3) = 2
          NETP%EDGENODEIND(EC,2) = 3
       ELSE
          PRINT*,'PROBLEM IN DYNAMICPROCESSES WITH 1+2 FUSION'
       ENDIF

    ELSE
       PRINT*,'PROBLEM WITH DEGREE SUM != 2 or 3. Fusion only handles 1+1, 1+2 nodes'
       STOP 1
    ENDIF

    DIDFUSE = .TRUE.

  END SUBROUTINE FUSIONEVENT
  

  SUBROUTINE DIFFUSEONEDGES(NETP,DELT,SPECIESPRODUCER)
     IMPLICIT NONE
     TYPE(DYNNETWORK), POINTER :: NETP
     DOUBLE PRECISION , INTENT(IN) :: DELT
     INTEGER, INTENT(IN) :: SPECIESPRODUCER(NETP%NSPECIES,NETP%NUMTRIGGERUNITS)
     DOUBLE PRECISION :: FLUX21, FLUX31, FLUX32, EL1, EL2, EL3
     INTEGER :: E1, E2, E3, DEG, NC, SC, PC
     DOUBLE PRECISION :: FLUXINTOEDGE(NETP%NEDGE,NETP%NSPECIES)

     FLUXINTOEDGE = 0D0

     DO NC = 1,NETP%NNODE
        IF(.NOT.NETP%NODEACT(NC)) CYCLE
        DEG = NETP%NODEDEG(NC)
        IF(DEG.EQ.1) THEN
           CYCLE
        ELSEIF(DEG.EQ.2) THEN
           E1 = NETP%NODEEDGE(NC,1)
           E2 = NETP%NODEEDGE(NC,2)
           ! need to set the edge lengths of degree 1 tips to their full length
           IF(NETP%NODEDEG(NETP%EDGENODE(E1,1))+NETP%NODEDEG(NETP%EDGENODE(E1,2)).LT.4) THEN
             EL1 = NETP%EDGELEN(E1)+NETP%STERICRAD
           ELSE
             EL1 = NETP%EDGELEN(E1)
           ENDIF
           IF(NETP%NODEDEG(NETP%EDGENODE(E2,1))+NETP%NODEDEG(NETP%EDGENODE(E2,2)).LT.4) THEN
             EL2 = NETP%EDGELEN(E2)+NETP%STERICRAD
           ELSE
             EL2 = NETP%EDGELEN(E2)
           ENDIF
           DO SC = 1,NETP%NSPECIES
              FLUX21 = NETP%DIFF(SC)*(NETP%CONC(E2,SC)-NETP%CONC(E1,SC))/(0.5*(EL1+EL2))
              FLUXINTOEDGE(E1,SC) = FLUXINTOEDGE(E1,SC) + FLUX21/EL1
              FLUXINTOEDGE(E2,SC) = FLUXINTOEDGE(E2,SC) - FLUX21/EL2
           ENDDO
        ELSEIF(DEG.EQ.3) THEN
           E1 = NETP%NODEEDGE(NC,1)
           E2 = NETP%NODEEDGE(NC,2)
           E3 = NETP%NODEEDGE(NC,3)
           ! if the other node of the edge is deg 1, increase its length for the calculation
           IF(NETP%NODEDEG(NETP%EDGENODE(E1,1))+NETP%NODEDEG(NETP%EDGENODE(E1,2)).LT.5) THEN
             EL1 = NETP%EDGELEN(E1)+NETP%STERICRAD
           ELSE
             EL1 = NETP%EDGELEN(E1)
           ENDIF
           IF(NETP%NODEDEG(NETP%EDGENODE(E2,1))+NETP%NODEDEG(NETP%EDGENODE(E2,2)).LT.5) THEN
             EL2 = NETP%EDGELEN(E2)+NETP%STERICRAD
           ELSE
             EL2 = NETP%EDGELEN(E2)
           ENDIF
           IF(NETP%NODEDEG(NETP%EDGENODE(E3,1))+NETP%NODEDEG(NETP%EDGENODE(E3,2)).LT.5) THEN
             EL3 = NETP%EDGELEN(E3)+NETP%STERICRAD
           ELSE
             EL3 = NETP%EDGELEN(E3)
           ENDIF
           DO SC = 1,NETP%NSPECIES
              FLUX21 = NETP%DIFF(SC)*(NETP%CONC(E2,SC)-NETP%CONC(E1,SC))/(0.5*(EL1+EL2))
              FLUX31 = NETP%DIFF(SC)*(NETP%CONC(E3,SC)-NETP%CONC(E1,SC))/(0.5*(EL1+EL3))
              FLUX32 = NETP%DIFF(SC)*(NETP%CONC(E3,SC)-NETP%CONC(E2,SC))/(0.5*(EL2+EL3))
              FLUXINTOEDGE(E1,SC) = FLUXINTOEDGE(E1,SC) + FLUX21/EL1 + FLUX31/EL1
              FLUXINTOEDGE(E2,SC) = FLUXINTOEDGE(E2,SC) + FLUX32/EL2 - FLUX21/EL2
              FLUXINTOEDGE(E3,SC) = FLUXINTOEDGE(E3,SC) - FLUX32/EL3 - FLUX31/EL3
           ENDDO
        ELSE
           PRINT*,'Problem with degree > 3 in DIFFUSEONEDGES'
        ENDIF
     ENDDO
     ! update all the concentrations according to the FVM calculation
     NETP%CONC = NETP%CONC + FLUXINTOEDGE*DELT
     ! if some edges should be held at constant concentration, reset them to their original values
     IF(NETP%FIXEDGES) THEN
        DO SC = 1,NETP%NSPECIES
           DO PC = 1,NETP%NUMTRIGGERUNITS
              E1 = SPECIESPRODUCER(SC,PC)
              NETP%CONC(E1,SC) = NETP%CONC(E1,SC) - FLUXINTOEDGE(E1,SC)*DELT
           ENDDO
        ENDDO
     ENDIF
  END SUBROUTINE DIFFUSEONEDGES


  SUBROUTINE PRODUCEANDDECAYSPECIES(NETP,DELT,SPECIESPRODUCER)
    ! produce more species in each edge at a constant rate
    ! decay species in each edge at a constant rate
    TYPE(DYNNETWORK), POINTER :: NETP
    DOUBLE PRECISION,INTENT(IN) :: DELT
    INTEGER, INTENT(IN) :: SPECIESPRODUCER(NETP%NSPECIES,NETP%NUMTRIGGERUNITS)
    DOUBLE PRECISION :: CELLR1SQR,ECL,OLDCONC(NETP%NEDGE,NETP%NSPECIES)
    INTEGER :: EC, SC, PC

    CELLR1SQR = NETP%CELLRAD1**2

    IF(ANY(NETP%PRODRATE.GT.0D0)) THEN
       IF(NETP%UNIQUETRIGGERS) THEN ! only produce on certain edges
          DO SC = 1,NETP%NSPECIES
             DO PC = 1,NETP%NUMTRIGGERUNITS
                EC = SPECIESPRODUCER(SC,PC)
                IF(.NOT.NETP%EDGEACT(EC)) CYCLE
                ECL = NETP%EDGELEN(EC)
                IF(NETP%NODEDEG(NETP%EDGENODE(EC,1)).EQ.1) ECL = ECL + NETP%STERICRAD
                IF(NETP%NODEDEG(NETP%EDGENODE(EC,2)).EQ.1) ECL = ECL + NETP%STERICRAD
                NETP%CONC(EC,SC) = NETP%CONC(EC,SC) + (NETP%MITOLEN/ECL)*NETP%PRODRATE(SC)*DELT
             ENDDO
          ENDDO
       ELSE ! same rate for all
          DO EC= 1,NETP%NEDGE
             IF(.NOT.NETP%EDGEACT(EC)) CYCLE
             NETP%CONC(EC,:) = NETP%CONC(EC,:) + (NETP%MITOLEN/ECL)*NETP%PRODRATE*DELT
          ENDDO
       ENDIF
    ENDIF

    IF(ANY(NETP%DECAYRATE.GT.0D0)) THEN
       IF(NETP%FIXEDGES) THEN
          OLDCONC = NETP%CONC
          DO EC = 1,NETP%NEDGE
             IF(.NOT.NETP%EDGEACT(EC)) CYCLE
             NETP%CONC(EC,:) = NETP%CONC(EC,:)*(1-NETP%DECAYRATE*DELT) ! decay rate must not go above 1/delt
          ENDDO
          DO SC = 1,NETP%NSPECIES
             DO PC = 1,NETP%NUMTRIGGERUNITS
                NETP%CONC(SPECIESPRODUCER(SC,PC),SC) = OLDCONC(SPECIESPRODUCER(SC,PC),SC)
             ENDDO
          ENDDO
       ELSE
          DO EC = 1,NETP%NEDGE
             IF(.NOT.NETP%EDGEACT(EC)) CYCLE
             NETP%CONC(EC,:) = NETP%CONC(EC,:)*(1-NETP%DECAYRATE*DELT) ! decay rate must not go above 1/delt
          ENDDO     
       ENDIF
    ENDIF

  END SUBROUTINE PRODUCEANDDECAYSPECIES


  SUBROUTINE GETFUSIONRATE(NETP,N1,N2,DEGSUM,FUSERATE)
    USE GENUTIL, ONLY : RAYRAYINTERSECT
    IMPLICIT NONE
     TYPE(DYNNETWORK), POINTER :: NETP
     INTEGER, INTENT(IN) :: N1,N2,DEGSUM
     DOUBLE PRECISION, INTENT(OUT) :: FUSERATE
     
     INTEGER :: E1, E2, NEIGHBS(3), NEIGHBSTMP(1), D1NODE, D2NODE, E3IDX
     DOUBLE PRECISION, DIMENSION(NETP%DIM) :: Q1, Q2, Q3, Q4, NORMALVEC
     DOUBLE PRECISION :: THETA, COSTHETA, R31, R12, LEN123, LENS(3), R(3), ENERGY2, ENERGY3, FUSEFACTOR
     DOUBLE PRECISION :: NMAG, R41, COSSQRTHETA, LEN1234, EOUTOFPLANE
     DOUBLE PRECISION :: V1(3), V2(3)
 
     IF ((.NOT.NETP%NODEACT(N1)).OR.(.NOT.NETP%NODEACT(N2))) THEN
        PRINT*, 'ERROR IN CHECKBENDCOST: inactive node(s) ',N1,N2
        STOP 1
     ENDIF
 
     IF(DEGSUM.EQ.2) THEN !We have 1-1 fusion
        CALL GETNEIGHBNODES(NETP,N1,NEIGHBSTMP)
        NEIGHBS(2) = NEIGHBSTMP(1)
        CALL GETNEIGHBNODES(NETP,N2,NEIGHBSTMP)
        NEIGHBS(3) = NEIGHBSTMP(1)
        !Assume N1 and N2 will occupy the same position in space
        Q2 = NETP%NODEPOS(:,NEIGHBS(2))
        Q3 = NETP%NODEPOS(:,NEIGHBS(3))

        V1 = NETP%NODEPOS(:,N1) - NETP%NODEPOS(:,NEIGHBS(2))
        V2 = NETP%NODEPOS(:,N2) - NETP%NODEPOS(:,NEIGHBS(3))
        Q1 = 0.5*(NETP%NODEPOS(:,N1)+NETP%STERICRAD*V1/SQRT(SUM(V1**2))+ &
        & NETP%NODEPOS(:,N2)+NETP%STERICRAD*V2/SQRT(SUM(V2**2)))
                
        R31 = NORM(Q3-Q1); R12 = NORM(Q1-Q2)
        LEN123 = NETP%EDGELEN(NETP%NODEEDGE(N1,1)) + NETP%EDGELEN(NETP%NODEEDGE(N2,1))
        LEN123 = LEN123 + 2*NETP%STERICRAD ! get full length including end caps of fusing nodes
        IF(NETP%NODEDEG(NEIGHBS(2)).EQ.1) LEN123 = LEN123+NETP%STERICRAD !get the full edge length including end cap of neighbor
        IF(NETP%NODEDEG(NEIGHBS(3)).EQ.1) LEN123 = LEN123+NETP%STERICRAD !get the full edge length including end cap of neighbor
        LEN123 = 0.5*LEN123 !average of the 2 edge lengths
        ENERGY2 = NETP%ALPHA1*NETP%MITOLEN/LEN123 * (1 - DOT_PRODUCT(Q3-Q1,Q1-Q2)/(R31*R12))
        
        FUSEFACTOR = EXP(-ENERGY2)
        FUSERATE = NETP%FUSERATE1*FUSEFACTOR
     ELSE !We have 1-2 fusion
        IF(NETP%NODEDEG(N1).EQ.1) THEN
           D1NODE = N1; D2NODE = N2
        ELSE
           D1NODE = N2; D2NODE = N1
        ENDIF
        CALL GETNEIGHBNODES(NETP,D1NODE,NEIGHBS(1))
        CALL GETNEIGHBNODES(NETP,D2NODE,NEIGHBS(2:3))

        LENS(1) = NETP%EDGELEN(NETP%NODEEDGE(D1NODE,1)) + NETP%STERICRAD ! add end cap to edge length of d1 node edge
        IF(NETP%NODEDEG(NEIGHBS(1)).EQ.1) LENS(1) = LENS(1)+NETP%STERICRAD !get the full edge length including end cap
        LENS(2) = NETP%EDGELEN(NETP%NODEEDGE(D2NODE,1))
        IF(NETP%NODEDEG(NEIGHBS(2)).EQ.1) LENS(2) = LENS(2)+NETP%STERICRAD
        LENS(3) = NETP%EDGELEN(NETP%NODEEDGE(D2NODE,2))
        IF(NETP%NODEDEG(NEIGHBS(3)).EQ.1) LENS(3) = LENS(3)+NETP%STERICRAD

        V1 = NETP%NODEPOS(:,D1NODE) - NETP%NODEPOS(:,NEIGHBS(1))
        V1 = V1/SQRT(SUM(V1**2))
        V2 = NETP%NODEPOS(:,D1NODE) + NETP%STERICRAD*V1 - NETP%NODEPOS(:,D2NODE)

        ! location of new D3 node
        Q1 = 0.5*(NETP%NODEPOS(:,D2NODE)+NETP%STERICRAD*V2/SQRT(SUM(V2**2)) + NETP%NODEPOS(:,D1NODE)+NETP%STERICRAD*V1)
        R(1) = NORM(Q1-NETP%NODEPOS(:,NEIGHBS(1)))
        R(2) = NORM(Q1-NETP%NODEPOS(:,NEIGHBS(2)))
        R(3) = NORM(Q1-NETP%NODEPOS(:,NEIGHBS(3)))

        ENERGY3 = 0D0; EOUTOFPLANE = 0D0

        Q2 = NETP%NODEPOS(:,NEIGHBS(1))
        DO E2 = 2,3
           Q3 = NETP%NODEPOS(:,NEIGHBS(E2))
           LEN123 = 0.5*(LENS(1) + LENS(E2))
           COSTHETA = DOT_PRODUCT(Q3-Q1,Q1-Q2)/(R(1)*R(E2))
           IF(ABS(COSTHETA).GE.1D0.AND.ABS(COSTHETA).LT.1.01D0) THEN
              COSTHETA = NINT(COSTHETA/ABS(COSTHETA))*1D0
           ENDIF
           THETA = ACOS(COSTHETA)
           ENERGY3 = ENERGY3 + NETP%ALPHA2*NETP%MITOLEN/LEN123 * (1 - COS(THETA - NETP%THETA0))
        ENDDO
        FUSEFACTOR = EXP(-(ENERGY3))

        IF(NETP%USEOUTOFPLANE) THEN
           DO E1 = 1,3
              DO E2 = E1+1,3
                 Q2 = NETP%NODEPOS(:,NEIGHBS(E1)); Q3 = NETP%NODEPOS(:,NEIGHBS(E2))
                 LEN1234 = SUM(LENS)/3
                 E3IDX = 6 - E1 - E2
                 Q4 = NETP%NODEPOS(:,NEIGHBS(E3IDX))
                 R41 = SQRT(SUM((Q4-Q1)**2))
                 NORMALVEC = 0D0
                 CALL CROSS_PRODUCT(Q2-Q1,Q3-Q1,NORMALVEC)
                 NMAG = SQRT(SUM(NORMALVEC**2))
                 COSSQRTHETA = 1 - (DOT_PRODUCT(NORMALVEC,Q4-Q1)/(NMAG*R41))**2
                 EOUTOFPLANE = EOUTOFPLANE + NETP%ALPHAPLANE*NETP%MITOLEN/LEN1234 * (1 - SQRT(COSSQRTHETA))
              ENDDO
           ENDDO
           FUSEFACTOR = EXP(-(EOUTOFPLANE+ENERGY3))
        ENDIF
        FUSERATE = NETP%FUSERATE2*FUSEFACTOR
     ENDIF
   END SUBROUTINE GETFUSIONRATE

END MODULE DYNAMICPROCESSES