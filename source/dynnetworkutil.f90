MODULE DYNNETWORKUTIL
  ! module definiting a mitochondrial network object and basic operations using it

  USE STACKUTIL, ONLY : STACK, PUSH, POP
  USE KEYS, ONLY : VERBOSE
  USE GENUTIL, ONLY : NORM, CROSS_PRODUCT

  TYPE DYNNETWORK
     ! define a dynamic network structure (connectivity and geometry)
     ! where each edge represents (eg) a mitochondrial unit and has a preferred length
     ! nodes are junctions between fused units (eg: mitochondria)

     ! dimension of the space the network is embedded in
     ! all edges are assumed to be straight lines btwn connected nodes

     INTEGER :: DIM

     ! ----------------------
     ! information on network nodes
     ! ----------------------
     INTEGER :: NNODE ! number of nodes
     ! degree (number of branches) of each node
     INTEGER, POINTER :: NODEDEG(:)
     ! list of branch indices each node connects to
     ! NODEEDGE(i,j) = for i-th node, what is the jth edge connecting to it
     INTEGER, POINTER :: NODEEDGE(:,:)
     ! NODEEDGEIND(i,j) = for i-th node, what is its index in the jth edge EDGENODE array
     INTEGER, POINTER :: NODEEDGEIND(:,:)
     ! spatial location of node
     DOUBLE PRECISION, POINTER :: NODEPOS(:,:)
     ! grid positions of nodes
     INTEGER, POINTER :: GRID(:,:)
     ! Which nodes are active, which are fixed
     LOGICAL, POINTER :: NODEACT(:)
     ! Stack listing empty node and edge indices
     TYPE(STACK), POINTER :: NODESTACK, EDGESTACK

     ! ------------------
     ! information on network branches
     ! ------------------
     INTEGER :: NEDGE ! Number of branches
     ! nodes at the start and end of each branch
     INTEGER, POINTER :: EDGENODE(:,:)
     ! for start/end node of this edge, what is the index of the edge in the NODEEDGE array?
     INTEGER, POINTER :: EDGENODEIND(:,:)
     ! EDGELEN = ground state length for each edge
     DOUBLE PRECISION, POINTER :: EDGELEN(:)
     ! which edges are active
     LOGICAL, POINTER :: EDGEACT(:)
     ! length of non-terminal edges
     DOUBLE PRECISION :: MITOLEN

     ! ------------------
     ! information on biomolecule diffusion
     ! ------------------
     INTEGER :: NSPECIES !number of types of diffusing species
     ! diffusion coefficient of each species
     DOUBLE PRECISION, POINTER :: DIFF(:)
     ! concentration on each edge for each species
     DOUBLE PRECISION, POINTER :: CONC(:,:)
     ! production rate for each species
     DOUBLE PRECISION, POINTER :: PRODRATE(:)
     ! number of edges which produce or are held at fixed concentration
     INTEGER :: NUMTRIGGERUNITS
     ! decay rate for each species
     DOUBLE PRECISION, POINTER :: DECAYRATE(:)
     ! should source edges have fixed value?
     LOGICAL :: FIXEDGES
     ! should each diffusing species be produced by a unique trigger edge?
     LOGICAL :: UNIQUETRIGGERS


     ! ------------------
     ! yeast model
     ! ------------------
     ! whether to use yeast model
     LOGICAL :: DOYEAST  
     ! Forcing Strength of the yeast tethering
     DOUBLE PRECISION :: YEASTCONF
     ! Rates to tether onto and off of the cell boundary, range at which tethering can occur
     DOUBLE PRECISION :: YEASTONRATE, YEASTOFFRATE, YEASTBINDRANGE
     ! tells if each node is tethered or not
     LOGICAL, POINTER :: TETHERED(:)


     ! --------------
     ! enclosing planes boundary information
     ! --------------
     ! number of enclosing planes
     INTEGER :: NPLANE
     ! vertex for each plane
     DOUBLE PRECISION, POINTER :: PLANEVERTS(:,:)
     ! edges for each plane
     DOUBLE PRECISION, POINTER :: PLANEEDGES(:,:,:)
     ! normal vector for each plane
     DOUBLE PRECISION, POINTER :: PLANENORMALS(:,:)
     ! side length of each plane along the directions given by planeedges
     DOUBLE PRECISION, POINTER :: PLANESIZES(:,:)

     ! --------------
     ! fusion reactivation
     ! --------------
     ! status of node. active or no?
     LOGICAL, POINTER :: FUSEACTIVE(:)
     ! rate of returning to fuse-active state
     DOUBLE PRECISION :: RECHARGERATE


     ! --------------
     ! parameters
     ! --------------
     ! stretch modulus
     DOUBLE PRECISION :: ESTR
     ! bend modulus for degree 2 and 3 nodes, desired bending angle
     DOUBLE PRECISION :: BENDMOD1, BENDMOD2, THETA0
     ! whether to use a third bending energy (out of plane angle)
     LOGICAL :: USEOUTOFPLANE
     ! bending modulus for out-of-plane coordinate
     DOUBLE PRECISION :: BENDMODPLANE
     ! cutoff to start calculating steric forces
     DOUBLE PRECISION :: STERICRAD, STERICMOD
     ! termal energy, friction
     DOUBLE PRECISION :: KT, FRICT
     ! fission rate (per time per fissable node) and multiplier for D3 fissions
     DOUBLE PRECISION :: FISSRATE, D3FISSMULT
     ! fission radius (how far to pull back nodes on fission)
     DOUBLE PRECISION :: FISSRAD
     ! confinement (radius, prefactor, power)
     DOUBLE PRECISION :: CELLRAD1, ECONF
     ! cap distance
     DOUBLE PRECISION :: CONTACTRAD
     ! fusion: rate, bending threshold, etc
     DOUBLE PRECISION :: FUSERATE1, FUSERATE2, ALPHA1, ALPHA2, ALPHAPLANE
     ! whether to use the GRID, edgegrid, and nodegrid
     LOGICAL :: USEGRID, USEEDGEGRID, USENODEGRID
     ! spacing for spatial grid
     INTEGER, POINTER :: GRIDSPACES(:)

     
     ! arrays have been set? topology of network has been set?
     LOGICAL :: ARRAYSET=.FALSE., STRUCTURESET=.FALSE.
  END TYPE DYNNETWORK

CONTAINS

  SUBROUTINE GETMECHFORCES(NETP,ENERGY,FORCES,GETFORCES)
    ! get current mechanical (eg: elastic, steric) energy
    ! and corresponding forces on the network

    IMPLICIT NONE
    TYPE(DYNNETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY
    DOUBLE PRECISION, INTENT(OUT) :: FORCES(NETP%DIM,NETP%NNODE)
    LOGICAL, INTENT(IN) :: GETFORCES
    DOUBLE PRECISION :: ESTRETCH, EBEND, ECONF, ESTERIC, ETUBE, EOUTOFPLANE
    DOUBLE PRECISION, DIMENSION(NETP%DIM,NETP%NNODE) :: FSTRETCH, FBEND, FCONF, FSTERIC, FTUBE, FOUTOFPLANE
    INTEGER :: EC, N1, N2, NC, DC
    DOUBLE PRECISION, DIMENSION(NETP%DIM) :: Q1, Q2, MINPOS, MAXPOS, PARTSZ
    DOUBLE PRECISION :: DELTA
    INTEGER :: BOUNDS(NETP%DIM,NETP%NEDGE,2), GRIDALT(2,NETP%DIM,NETP%NNODE)
    INTEGER :: GC1, GC2, GC3, IDX, GRIDCOUNTS(PRODUCT(NETP%GRIDSPACES)) 

    ESTRETCH = 0D0; ECONF = 0D0; EBEND = 0D0; ESTERIC = 0D0; ETUBE = 0D0; EOUTOFPLANE = 0D0; ENERGY = 0D0
    FSTRETCH = 0D0; FCONF = 0D0; FBEND = 0D0; FSTERIC = 0D0; FTUBE = 0D0; FOUTOFPLANE = 0D0; FORCES = 0D0

    DO DC = 1,NETP%DIM
      MINPOS(DC) = MINVAL(NETP%NODEPOS(DC,:),NETP%NODEACT)
      MAXPOS(DC) = MAXVAL(NETP%NODEPOS(DC,:),NETP%NODEACT)
      ! the grid box size shouldn't be smaller than this
      PARTSZ(DC) = MAX((MAXPOS(DC)-MINPOS(DC))/NETP%GRIDSPACES(DC), 2*NETP%STERICRAD) 
    ENDDO
    GRIDALT = 0
    GRIDCOUNTS = 0
    BOUNDS = 0

    ! Find confinement & bending energy/forces on the nodes
    DO NC = 1,NETP%NNODE
       IF (.NOT.NETP%NODEACT(NC)) CYCLE

       IF(NETP%CELLRAD1.GT.0) THEN 
          CALL GETSPHERECONFINEMENTENERGY(NETP,NC,ECONF,FCONF,GETFORCES)
       ELSEIF(NETP%NPLANE.GT.0) THEN
          CALL GETPLANECONFINEMENT(NETP,NC,ECONF,FCONF,GETFORCES)
       ENDIF

       ! Bending forces
       CALL GETBENDINGFORCES(NETP,NC,EBEND,EOUTOFPLANE,FBEND,FOUTOFPLANE,GETFORCES)

       ! set up nodes for use with grid
       IF(NETP%USEGRID) THEN
          NETP%GRID(:,NC) = MIN(MAX(CEILING((NETP%NODEPOS(:,NC)-MINPOS)/PARTSZ),1),NETP%GRIDSPACES)
          IF(NETP%USEEDGEGRID) THEN
             DO DC = 1,NETP%DIM
                DELTA = NETP%NODEPOS(DC,NC) - (MINPOS(DC) + PARTSZ(DC)*(NETP%GRID(DC,NC)-1))
                IF(NETP%GRID(DC,NC).GT.1.AND.DELTA.LT.NETP%STERICRAD) THEN
                   GRIDALT(1,DC,NC) = -1
                ELSEIF(NETP%GRID(DC,NC).LT.NETP%GRIDSPACES(DC).AND. &
                & DELTA.GT.(PARTSZ(DC)-NETP%STERICRAD)) THEN
                   GRIDALT(2,DC,NC) = 1
                ENDIF
             ENDDO
          ENDIF
       ENDIF
    ENDDO

    IF(NETP%USEGRID.AND.NETP%USEEDGEGRID) THEN ! use edge grid
       DO EC = 1,NETP%NEDGE
          IF (.NOT.NETP%EDGEACT(EC)) CYCLE
          N1 = NETP%EDGENODE(EC,1); N2 = NETP%EDGENODE(EC,2)
          Q1 = NETP%NODEPOS(:,N1); Q2 = NETP%NODEPOS(:,N2)
   
          CALL GETSTRETCHINGONEDGE(NETP,EC,N1,N2,Q1,Q2,ESTRETCH,FSTRETCH,GETFORCES)

          BOUNDS(:,EC,1) = MIN(NETP%GRID(:,N1)+GRIDALT(1,:,N1),NETP%GRID(:,N2)+GRIDALT(1,:,N2))
          BOUNDS(:,EC,2) = MAX(NETP%GRID(:,N1)+GRIDALT(2,:,N1),NETP%GRID(:,N2)+GRIDALT(2,:,N2))
          DO GC1 = BOUNDS(1,EC,1),BOUNDS(1,EC,2)
             DO GC2 = BOUNDS(2,EC,1),BOUNDS(2,EC,2)
                DO GC3 = BOUNDS(3,EC,1),BOUNDS(3,EC,2)
                   ! add the edge to this grid box
                   ! idx = z + y * ZSize + x * ZSize * YSize
                   IDX = 1 + (GC3-1) + (GC2-1)*NETP%GRIDSPACES(3) + (GC1-1)*NETP%GRIDSPACES(3)*NETP%GRIDSPACES(2)
                   GRIDCOUNTS(IDX) = GRIDCOUNTS(IDX) + 1 ! count the number of edges added to each grid cube
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       ! use edge grid to find steric overlaps
       CALL GETSTERICSWITHEDGEGRID(NETP,BOUNDS,ESTERIC,FSTERIC,GETFORCES)

    ELSE ! do it without the edge grid
       DO EC = 1,NETP%NEDGE
          IF (.NOT.NETP%EDGEACT(EC)) CYCLE
          N1 = NETP%EDGENODE(EC,1); N2 = NETP%EDGENODE(EC,2)
          Q1 = NETP%NODEPOS(:,N1); Q2 = NETP%NODEPOS(:,N2)
 
          CALL GETSTRETCHINGONEDGE(NETP,EC,N1,N2,Q1,Q2,ESTRETCH,FSTRETCH,GETFORCES)
 
          IF(NETP%USEGRID) THEN
             BOUNDS(:,EC,1) = MIN(NETP%GRID(:,N1),NETP%GRID(:,N2))-1
             BOUNDS(:,EC,2) = MAX(NETP%GRID(:,N1),NETP%GRID(:,N2))+1             
          ENDIF
 
          ! Steric Forces on overlapping edges which are not neighbors
          CALL GETSTERICSBASIC(NETP,EC,N1,N2,Q1,Q2,BOUNDS,ESTERIC,FSTERIC,GETFORCES) ! use node grid or check every pair
       ENDDO
    ENDIF
       
    FORCES = -FSTRETCH-FCONF-FBEND-FSTERIC-FTUBE-FOUTOFPLANE
    ENERGY = ESTRETCH+ECONF+EBEND+ESTERIC+ETUBE+EOUTOFPLANE
  END SUBROUTINE GETMECHFORCES


  SUBROUTINE GETSTERICSWITHEDGEGRID(NETP,BOUNDS,ESTERIC,FSTERIC,GETFORCES)
     IMPLICIT NONE
     TYPE(DYNNETWORK), POINTER :: NETP
     INTEGER, INTENT(IN) :: BOUNDS(NETP%DIM,NETP%NEDGE,2)
     DOUBLE PRECISION, INTENT(OUT) :: ESTERIC, FSTERIC(NETP%DIM,NETP%NNODE)
     LOGICAL, INTENT(IN) :: GETFORCES
     INTEGER :: GRIDEDGES(NETP%NEDGE+1,PRODUCT(NETP%GRIDSPACES))
     LOGICAL :: INTERACTIONMATRIX(NETP%NEDGE, NETP%NEDGE)
     INTEGER :: EC, E2, N1, N2, N3, N4, GC1, GC2, GC3, IDX, NINCURR, GEC, GEC2
     DOUBLE PRECISION, DIMENSION(NETP%DIM) :: Q1, Q2

     GRIDEDGES = 0; INTERACTIONMATRIX = .FALSE.

     DO EC = 1,NETP%NEDGE
        IF(.NOT.NETP%EDGEACT(EC)) CYCLE
        N1 = NETP%EDGENODE(EC,1); N2 = NETP%EDGENODE(EC,2)
        DO GC1 = BOUNDS(1,EC,1),BOUNDS(1,EC,2)
           DO GC2 = BOUNDS(2,EC,1),BOUNDS(2,EC,2)
              DO GC3 = BOUNDS(3,EC,1),BOUNDS(3,EC,2)
                 ! add the edge to this grid box
                 !i = z + y * ZSize + x * ZSize * YSize
                 IDX = 1 + (GC3-1) + (GC2-1)*NETP%GRIDSPACES(3) + (GC1-1)*NETP%GRIDSPACES(3)*NETP%GRIDSPACES(2)
                 NINCURR = GRIDEDGES(1,IDX) + 1
                 GRIDEDGES(1,IDX) = NINCURR
                 GRIDEDGES(NINCURR+1,IDX) = EC                   
              ENDDO
           ENDDO
        ENDDO
     ENDDO 
     DO IDX = 1,PRODUCT(NETP%GRIDSPACES)
         NINCURR = GRIDEDGES(1,IDX)
         IF(NINCURR.LT.2) CYCLE
         DO GEC = 2,NINCURR+1
            EC = GRIDEDGES(GEC,IDX)
            N1 = NETP%EDGENODE(EC,1); N2 = NETP%EDGENODE(EC,2)
            Q1 = NETP%NODEPOS(:,N1); Q2 = NETP%NODEPOS(:,N2)
            DO GEC2 = GEC+1,NINCURR+1
               E2 = (GRIDEDGES(GEC2,IDX))
               IF(EC.EQ.E2.OR.INTERACTIONMATRIX(EC,E2).OR.INTERACTIONMATRIX(E2,EC)) CYCLE
               INTERACTIONMATRIX(EC,E2) = .TRUE.; INTERACTIONMATRIX(E2,EC) = .TRUE.
               N3 = NETP%EDGENODE(E2,1); N4 = NETP%EDGENODE(E2,2)
               IF((N1.NE.N3).AND.(N1.NE.N4).AND.(N2.NE.N3).AND.(N2.NE.N4)) THEN
                  CALL CHECKPAIROVERLAPANDGETFORCES(NETP,N1,N2,N3,N4,Q1,Q2,FSTERIC,ESTERIC,GETFORCES)
               ENDIF
            ENDDO
         ENDDO
     ENDDO  
  END SUBROUTINE GETSTERICSWITHEDGEGRID


  SUBROUTINE GETSTERICSBASIC(NETP,EC,N1,N2,Q1,Q2,BOUNDS,ESTERIC,FSTERIC,GETFORCES)
     IMPLICIT NONE
     TYPE(DYNNETWORK), POINTER :: NETP
     DOUBLE PRECISION, INTENT(OUT) :: ESTERIC, FSTERIC(NETP%DIM,NETP%NNODE)
     LOGICAL, INTENT(IN) :: GETFORCES
     INTEGER, INTENT(IN) :: EC,N1,N2,BOUNDS(NETP%DIM,NETP%NEDGE,2)
     DOUBLE PRECISION, DIMENSION(NETP%DIM), INTENT(IN) :: Q1,Q2
     INTEGER :: E2, N3, N4
     LOGICAL :: CLOSENUFF

     DO E2 = EC+1,NETP%NEDGE
        IF(.NOT.NETP%EDGEACT(E2)) CYCLE
        N3 = NETP%EDGENODE(E2,1); N4 = NETP%EDGENODE(E2,2)
        IF((N1.NE.N3).AND.(N1.NE.N4).AND.(N2.NE.N3).AND.(N2.NE.N4)) THEN
           IF(NETP%USEGRID) THEN
              CLOSENUFF = .TRUE.
              IF((ANY(NETP%GRID(:,N3).LT.BOUNDS(:,EC,1)).OR.ANY(NETP%GRID(:,N3).GT.BOUNDS(:,EC,2))).AND. &
              & (ANY(NETP%GRID(:,N4).LT.BOUNDS(:,EC,1)).OR.ANY(NETP%GRID(:,N4).GT.BOUNDS(:,EC,2)))) THEN
                 CLOSENUFF = .FALSE.
              ENDIF
           ENDIF
           IF((.NOT.NETP%USEGRID).OR.CLOSENUFF) THEN
              CALL CHECKPAIROVERLAPANDGETFORCES(NETP,N1,N2,N3,N4,Q1,Q2,FSTERIC,ESTERIC,GETFORCES)
           ENDIF
        ENDIF
     ENDDO
  END SUBROUTINE GETSTERICSBASIC


  SUBROUTINE CHECKPAIROVERLAPANDGETFORCES(NETP,N1,N2,N3,N4,Q1,Q2,FSTERIC,ESTERIC,GETFORCES)
     IMPLICIT NONE
     TYPE(DYNNETWORK), POINTER :: NETP
     DOUBLE PRECISION, INTENT(OUT) :: ESTERIC, FSTERIC(NETP%DIM,NETP%NNODE)
     LOGICAL, INTENT(IN) :: GETFORCES
     INTEGER, INTENT(IN) :: N1,N2,N3,N4
     DOUBLE PRECISION, DIMENSION(NETP%DIM), INTENT(IN) :: Q1,Q2
     INTEGER :: NC, I, J
     DOUBLE PRECISION, DIMENSION(NETP%DIM) :: Q3, Q4, DR
     DOUBLE PRECISION :: R, R12, R34, RSQR, MUA, MUB
     DOUBLE PRECISION :: MUBYR(NETP%DIM,2,4), DF(NETP%DIM,4)

     Q3 = NETP%NODEPOS(:,N3); Q4 = NETP%NODEPOS(:,N4)
     DR = 0.5*(Q1+Q2-Q3-Q4)
     RSQR = DOT_PRODUCT(DR,DR)
     R12 = SQRT(DOT_PRODUCT(Q1-Q2,Q1-Q2))
     R34 = SQRT(DOT_PRODUCT(Q3-Q4,Q3-Q4))
     IF(SQRT(RSQR).LT.(R12/2+R34/2+2*NETP%STERICRAD)) THEN
        ! Assign MUA,MUB as weights to determine the intersection points of closest approach
        ! PointA = Q1 + MUA*(Q2-Q1), PointB = Q3 + MUB*(Q4-Q3)
        CALL SHORTESTPATHPOINTS(Q1,Q2,Q3,Q4,NETP%DIM,MUA,MUB,MUBYR)
        DR = (1-MUA)*Q1 + MUA*Q2 - (1-MUB)*Q3 - MUB*Q4 !vector pointing from B to A
        RSQR = DOT_PRODUCT(DR,DR)
        IF(RSQR.LT.((2*NETP%STERICRAD)**2)) THEN !steric forces act inside set distance
           R = SQRT(RSQR)
           ESTERIC = ESTERIC + NETP%STERICMOD/2 * (R-2*NETP%STERICRAD)**2
           IF(GETFORCES) THEN !find the associated forces
              DF(:,1) = (1-MUA)*DR; DF(:,2) = MUA*DR
              DF(:,3) = (MUB-1)*DR; DF(:,4) = -MUB*DR
              DO NC = 1,4
                 DO I = 1,NETP%DIM
                    DO J = 1,NETP%DIM
                       DF(I,NC) = DF(I,NC) + DR(J)*(MUBYR(I,1,NC)*(Q2(J)-Q1(J)) + MUBYR(I,2,NC)*(Q3(J)-Q4(J)))
                    ENDDO
                 ENDDO
              ENDDO
              DF = DF*NETP%STERICMOD*(1-2*NETP%STERICRAD/R)
              FSTERIC(:,N1) = FSTERIC(:,N1) + DF(:,1)
              FSTERIC(:,N2) = FSTERIC(:,N2) + DF(:,2)
              FSTERIC(:,N3) = FSTERIC(:,N3) + DF(:,3)
              FSTERIC(:,N4) = FSTERIC(:,N4) + DF(:,4)
           ENDIF
        ENDIF
     ENDIF

  END SUBROUTINE CHECKPAIROVERLAPANDGETFORCES


  SUBROUTINE GETSTRETCHINGONEDGE(NETP,EC,N1,N2,Q1,Q2,ESTRETCH,FSTRETCH,GETFORCES)
     IMPLICIT NONE
     TYPE(DYNNETWORK), POINTER :: NETP
     INTEGER, INTENT(IN) :: EC,N1,N2
     DOUBLE PRECISION, INTENT(IN) :: Q1(NETP%DIM), Q2(NETP%DIM)
     DOUBLE PRECISION, INTENT(OUT) :: ESTRETCH, FSTRETCH(NETP%DIM,NETP%NNODE)
     LOGICAL, INTENT(IN) :: GETFORCES
     DOUBLE PRECISION :: DR(NETP%DIM), FTMP(NETP%DIM), NDR, LEN, COEFF, TMP 

     DR = Q2-Q1
     NDR = SQRT(DOT_PRODUCT(DR,DR)) ! current segment length
  
     LEN = NETP%EDGELEN(EC)
     IF(NETP%NODEDEG(NETP%EDGENODE(EC,1)).EQ.1) LEN = LEN + NETP%STERICRAD
     IF(NETP%NODEDEG(NETP%EDGENODE(EC,2)).EQ.1) LEN = LEN + NETP%STERICRAD
     COEFF = NETP%ESTR/(2*LEN)
     TMP = NDR-NETP%EDGELEN(EC) ! difference from ground-state length
  
     ! stretching energy
     ESTRETCH = ESTRETCH + COEFF*TMP**2
  
     IF (GETFORCES) THEN ! get the corresponding forces
        FTMP = 2*COEFF*TMP/NDR*DR
        FSTRETCH(:,N2)= FSTRETCH(:,N2)+FTMP
        FSTRETCH(:,N1) = FSTRETCH(:,N1)-FTMP
     ENDIF

  END SUBROUTINE GETSTRETCHINGONEDGE


  SUBROUTINE OUTSIDEPLANES(NETP,POS,OUTSIDE,POSZ)
     IMPLICIT NONE
     TYPE(DYNNETWORK), POINTER :: NETP
     DOUBLE PRECISION, INTENT(IN) :: POS(NETP%DIM)
     LOGICAL, INTENT(OUT) :: OUTSIDE
     DOUBLE PRECISION, INTENT(OUT) :: POSZ(NETP%NPLANE)
     LOGICAL :: INSIDE(NETP%DIM,2)
     DOUBLE PRECISION :: DENOM, POSX, POSY, MU, XSZ, YSZ
     DOUBLE PRECISION, DIMENSION(NETP%DIM) :: PROJ, XDIR, YDIR, ORIGIN, NORMAL
     INTEGER :: PC, DC
   
     POSZ = 0D0
     INSIDE = .FALSE.

     DO PC = 1,NETP%NPLANE
        POSZ(PC) = DOT_PRODUCT(POS-NETP%PLANEVERTS(:,PC),NETP%PLANENORMALS(:,PC))
     ENDDO

     ! first check to see if the node is contained by the planes
     OUTER: DO DC = 1,3
        DO PC = 1,NETP%NPLANE
           IF(POSZ(PC).LT.0) CYCLE
           NORMAL = NETP%PLANENORMALS(:,PC)
           DENOM = NORMAL(DC)
           IF(DENOM.NE.0D0) THEN
              MU = -POSZ(PC)/DENOM
              PROJ = POS
              PROJ(DC) = PROJ(DC) + MU
              ORIGIN = NETP%PLANEVERTS(:,PC)
              PROJ = PROJ-ORIGIN
              XDIR = NETP%PLANEEDGES(:,1,PC)
              XSZ = NETP%PLANESIZES(1,PC)
              YDIR = NETP%PLANEEDGES(:,2,PC)
              YSZ = NETP%PLANESIZES(2,PC)
              POSX = DOT_PRODUCT(PROJ,XDIR)
              POSY = DOT_PRODUCT(PROJ,YDIR)
              IF(0.LE.POSX.AND.XSZ.GE.POSX.AND.0.LE.POSY.AND.YSZ.GE.POSY) THEN
                 IF(MU.GT.0) THEN
                    INSIDE(DC,2) = .TRUE.
                 ELSE
                    INSIDE(DC,1) = .TRUE.
                 ENDIF
                 IF(INSIDE(DC,1).AND.INSIDE(DC,2)) CYCLE OUTER
              ELSE
                 ! intersection outside plane segment
              ENDIF
           ELSE
              ! no intersection, parallel
           ENDIF
        ENDDO
        IF((.NOT.INSIDE(DC,1)).OR.(.NOT.INSIDE(DC,2))) EXIT OUTER ! not contained in 1 of the directions
     ENDDO OUTER

     IF(ALL(INSIDE)) THEN
       OUTSIDE = .FALSE.
     ELSE
       OUTSIDE = .TRUE.
     ENDIF

  END SUBROUTINE OUTSIDEPLANES

  SUBROUTINE GETPLANECONFINEMENT(NETP,NC,ECONF,FCONF,GETFORCES)
     IMPLICIT NONE
     TYPE(DYNNETWORK), POINTER :: NETP
     DOUBLE PRECISION, INTENT(OUT) :: ECONF
     DOUBLE PRECISION, INTENT(OUT) :: FCONF(NETP%DIM,NETP%NNODE)
     LOGICAL, INTENT(IN) :: GETFORCES
     INTEGER, INTENT(IN) :: NC
     DOUBLE PRECISION, DIMENSION(NETP%DIM) :: POS, DR, PROJ, XDIR, YDIR, ORIGIN, NORMAL, NEWDR
     LOGICAL :: OUTSIDE
     DOUBLE PRECISION :: MINDISTSQR, POSX, POSY, POSZ(NETP%NPLANE), XSZ, YSZ
     INTEGER :: PC
     
     POSZ = 0D0; OUTSIDE = .FALSE.

     POS = NETP%NODEPOS(:,NC)
     CALL OUTSIDEPLANES(NETP,POS,OUTSIDE,POSZ)

     IF(OUTSIDE) THEN ! find the energy/force
       MINDISTSQR = 1e15
       DO PC = 1,NETP%NPLANE
          IF(POSZ(PC).GT.0) CYCLE

          NORMAL = NETP%PLANENORMALS(:,PC)
          ORIGIN = NETP%PLANEVERTS(:,PC)

          XDIR = NETP%PLANEEDGES(:,1,PC)
          XSZ = NETP%PLANESIZES(1,PC)
          YDIR = NETP%PLANEEDGES(:,2,PC)
          YSZ = NETP%PLANESIZES(2,PC)

          PROJ = POS - POSZ(PC)*NORMAL
          POSX = DOT_PRODUCT(PROJ-ORIGIN,XDIR)
          POSY = DOT_PRODUCT(PROJ-ORIGIN,YDIR)

          IF(0.LE.POSX.AND.XSZ.GE.POSX.AND.0.LE.POSY.AND.YSZ.GE.POSY) THEN
             ! intersection point is on the surface
             NEWDR = POSZ(PC)*NORMAL
          ELSEIF(0.LE.POSX.AND.XSZ.GE.POSX) THEN! intersection on "y" side
             IF(POSY.LT.0) THEN ! "bottom side"
                NEWDR = POS - (ORIGIN + POSX*XDIR)
             ELSE ! "top side"
                NEWDR = POS - (ORIGIN + YSZ*YDIR + POSX*XDIR)
             ENDIF
          ELSEIF(0.LE.POSY.AND.YSZ.GE.POSY) THEN ! intersection on "x" side
             IF(POSX.LT.0) THEN ! "left side"
                NEWDR = POS - (ORIGIN + POSY*YDIR)
             ELSE ! "right side"
                NEWDR = POS - (ORIGIN + XSZ*XDIR + POSY*YDIR)
             ENDIF
          ELSEIF(POSY.LT.0) THEN
             IF(POSX.LT.0) THEN !bottom left corner
                NEWDR = POS - ORIGIN
             ELSE ! bottom right corner
                NEWDR = POS - (ORIGIN + XSZ*XDIR)
             ENDIF
          ELSE
             IF(POSX.LT.0) THEN !top left corner
                NEWDR = POS - (ORIGIN + YSZ*YDIR)
             ELSE ! top right corner
                NEWDR = POS - (ORIGIN + YSZ*YDIR + XSZ*XDIR)
             ENDIF
         ENDIF
         IF(DOT_PRODUCT(NEWDR,NEWDR).LT.MINDISTSQR) THEN ! intersection is new closest
            DR = NEWDR
            MINDISTSQR = DOT_PRODUCT(NEWDR,NEWDR)
         ENDIF
       ENDDO

       ECONF = ECONF + NETP%ECONF/2*MINDISTSQR
       IF (GETFORCES) THEN
          FCONF(:,NC) = FCONF(:,NC) + NETP%ECONF*DR
       ENDIF
    ELSE
       ! no confinement for this node
    ENDIF

  END SUBROUTINE GETPLANECONFINEMENT


  SUBROUTINE GETSPHERECONFINEMENTENERGY(NETP,NC,ECONF,FCONF,GETFORCES)
     IMPLICIT NONE
     TYPE(DYNNETWORK), POINTER :: NETP
     DOUBLE PRECISION, INTENT(OUT) :: ECONF
     DOUBLE PRECISION, INTENT(OUT) :: FCONF(NETP%DIM,NETP%NNODE)
     LOGICAL, INTENT(IN) :: GETFORCES
     INTEGER, INTENT(IN) :: NC
     LOGICAL :: OUTSIDE
     DOUBLE PRECISION, DIMENSION(NETP%DIM) :: POS, DR
     DOUBLE PRECISION :: R, R1SQR, R2SQR, CELLR, CONF, R1, R2, CELLR1SQR

     OUTSIDE = .FALSE.
     R1SQR = 0D0; R2SQR = 0D0; R = 0D0; CELLR = 0D0; CONF = 0D0; R1 = 0D0; R2 = 0D0
     DR = 0D0

     POS = NETP%NODEPOS(:,NC)
     CELLR1SQR = NETP%CELLRAD1**2
     R1SQR = SUM(POS**2)
     OUTSIDE = R1SQR.GT.CELLR1SQR
     IF (OUTSIDE.OR.(NETP%DOYEAST.AND.NETP%TETHERED(NC))) THEN
        !find the confinement force
        ! in the Yeast model, if node is tethered to edge of cell then hold it there
        IF(OUTSIDE) THEN
           CONF = NETP%ECONF
        ELSE
           CONF = NETP%YEASTCONF
        ENDIF
        R = SQRT(R1SQR)
        ECONF = ECONF + CONF/2*(R-NETP%CELLRAD1)**2
        IF (GETFORCES) THEN
           FCONF(:,NC) = FCONF(:,NC) + CONF*(R-NETP%CELLRAD1)*POS/R
        ENDIF
     ENDIF

  END SUBROUTINE GETSPHERECONFINEMENTENERGY

! Takes the coordinates of 2 nodes connected by an edge, 2 more nodes
! connected by another edge, the dimension of the network space (2D, 3D, etc.)
! and fills in the values of MUA and MUB, which denote the closest points
! on the first and second edges, respectively, as weights of the original coords.
  SUBROUTINE SHORTESTPATHPOINTS(P1,P2,P3,P4,DIM,MUA,MUB,MUBYR)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: DIM
    DOUBLE PRECISION, INTENT(IN) :: P1(DIM),P2(DIM),P3(DIM),P4(DIM)
    DOUBLE PRECISION, INTENT(OUT) :: MUA, MUB, MUBYR(DIM,2,4)
    DOUBLE PRECISION :: D1343, D4321, D1321, D4343, D2121, DENOM
    DOUBLE PRECISION, DIMENSION(DIM) :: R21,R13,R43
    LOGICAL :: AOUT,BOUT

    AOUT=.FALSE.; BOUT=.FALSE.

    R21 = P2-P1
    R13 = P1-P3
    R43 = P4-P3
    D4321 = DOT_PRODUCT(R43,R21)
    D4343 = DOT_PRODUCT(R43,R43)
    D2121 = DOT_PRODUCT(R21,R21)
    D1343 = DOT_PRODUCT(R13,R43)
    D1321 = DOT_PRODUCT(R13,R21)
    DENOM = D2121*D4343 - D4321*D4321

    IF(DENOM.NE.0) THEN    
       MUA = (D1343*D4321 - D1321*D4343)/DENOM
       IF(MUA.GT.1D0) THEN
         AOUT = .TRUE.; MUA = 1D0
       ELSEIF(MUA.LT.0D0) THEN
         AOUT = .TRUE.; MUA = 0D0
       ENDIF
       MUB = (D1343 + MUA*D4321)/D4343
      
       IF(MUB.GT.1D0) THEN
          BOUT = .TRUE.; MUB = 1D0
       ELSEIF(MUB.LT.0D0) THEN
          BOUT = .TRUE.; MUB = 0D0
       ENDIF

       IF(BOUT) THEN
          MUA = (MUB*D4321 - D1321)/D2121
          IF(MUA.GT.1D0) THEN
            AOUT = .TRUE.; MUA = 1D0
          ELSEIF(MUA.LT.0D0) THEN
            AOUT = .TRUE.; MUA = 0D0
          ELSE
            AOUT = .FALSE.
          ENDIF
       ENDIF

    ELSE ! parallel lines
       MUA = 0D0
       MUB = D1343/D4343
       IF(MUB.GT.1D0) THEN
          BOUT = .TRUE.; MUB = 1D0
       ELSEIF(MUB.LT.0D0) THEN
          BOUT = .TRUE.; MUB = 0D0
       ENDIF

       IF(BOUT) THEN
          MUA = (MUB*D4321-D1321)/D2121
          IF(MUA.GT.1D0) THEN
            MUA = 1D0
          ELSEIF(MUA.LT.0D0) THEN
            MUA = 0D0
          ENDIF
       ENDIF
       ! for parallel lines the derivatives are undefined so best to set them equal to zero
       AOUT = .TRUE.; BOUT = .TRUE.
    ENDIF
    
    ! now find the derivatives
    IF(AOUT.AND.BOUT) THEN ! calculate the 4 point-point distances and 4 projection distances
       MUBYR = 0D0 ! no dependence. Both mua,mub are constants
    ELSEIF(AOUT) THEN ! Do single point Find on MUB
       !PRINT*,'MUA IS OUT'
       MUBYR(:,1,:) = 0D0
       IF(MUA.GT.5D-1) THEN
          MUBYR(:,2,1) = 0D0
          MUBYR(:,2,2) = R43/D4343
          MUBYR(:,2,3) = (2*P3-P2-P4 + 2*MUB*R43)/D4343
          MUBYR(:,2,4) = (P2-P3-2*MUB*R43)/D4343
       ELSE
          MUBYR(:,2,1) = R43/D4343
          MUBYR(:,2,2) = 0D0
          MUBYR(:,2,3) = (2*P3-P1-P4 + 2*MUB*R43)/D4343
          MUBYR(:,2,4) = (P1-P3-2*MUB*R43)/D4343
       ENDIF
    ELSEIF(BOUT) THEN ! Do single point find on MUA
       !PRINT*,'MUB IS OUT'
       MUBYR(:,2,:) = 0D0
       IF(MUB.GT.5D-1) THEN
          MUBYR(:,1,1) = (2*P1-P4-P2 + 2*MUA*R21)/D2121
          MUBYR(:,1,2) = (P4-P1-2*MUA*R21)/D2121
          MUBYR(:,1,3) = 0D0
          MUBYR(:,1,4) = R21/D2121
       ELSE
          MUBYR(:,1,1) = (2*P1-P3-P2 + 2*MUA*R21)/D2121
          MUBYR(:,1,2) = (P3-P1-2*MUA*R21)/D2121
          MUBYR(:,1,3) = R21/D2121
          MUBYR(:,1,4) = 0D0
       ENDIF
    ELSE
       ! Use the original method - no end points are used
       MUBYR(:,1,1) = (R43*(D4321-D1343) - D4343*(R21-R13) - 2*MUA*(D4321*R43 - D4343*R21)) / DENOM
       MUBYR(:,1,2) = ((D1343*R43 - D4343*R13) - 2*MUA*(D4343*R21 - D4321*R43)) / DENOM
       MUBYR(:,1,3) = (R21*(D4343-D1343) + 2*D1321*R43 - D4321*(R43+R13) - 2*MUA*(D4321*R21-D2121*R43)) / DENOM
       MUBYR(:,1,4) = (R13*D4321 + D1343*R21 - 2*D1321*R43 - 2*MUA*(D2121*R43 - D4321*R21)) / DENOM
       MUBYR(:,2,1) = ((1-MUA)*R43 + D4321*MUBYR(:,1,1)) / D4343
       MUBYR(:,2,2) = (D4321*MUBYR(:,1,2) + MUA*R43) / D4343
       MUBYR(:,2,3) = (-R43-R13 - MUA*R21 + MUBYR(:,1,3)*D4321 + 2*MUB*R43) / D4343
       MUBYR(:,2,4) = (R13 + MUA*R21 + MUBYR(:,1,4)*D4321 - 2*MUB*R43) / D4343
    ENDIF

  END SUBROUTINE SHORTESTPATHPOINTS

  SUBROUTINE GETBENDINGFORCES(NETP,NC,EBEND,EOUTOFPLANE,FBEND,FOUTOFPLANE,GETFORCES)
     IMPLICIT NONE
     TYPE(DYNNETWORK), POINTER :: NETP
     DOUBLE PRECISION, INTENT(OUT) :: EBEND, EOUTOFPLANE
     DOUBLE PRECISION, INTENT(OUT) :: FBEND(NETP%DIM,NETP%NNODE), FOUTOFPLANE(NETP%DIM,NETP%NNODE)
     LOGICAL, INTENT(IN) :: GETFORCES
     INTEGER, INTENT(IN) :: NC
     INTEGER :: N2IDX, N3IDX, DEG, NEIGHBS3(3), NEIGHBS2(2), N4IDX
     DOUBLE PRECISION, DIMENSION(NETP%DIM) :: Q1, Q2, Q3, Q4, NORMALVEC, CROSS1, CROSS2
     DOUBLE PRECISION :: THETA, COSTHETA, R31, R12, LEN123, SINEFACTOR, LEN1234, R41, NMAG, COSSQRTHETA, FORCEFACTOR, FLIP

     DEG = NETP%NODEDEG(NC)
     IF (DEG.EQ.1) THEN
        ! Force is zero. Leave energy at zero. Only bend on degree 2-3
     ELSEIF (DEG.EQ.2) THEN
        ! Node is degree 2. Find vectors to the connected nodes and find energy
        Q1 = NETP%NODEPOS(:,NC)
        CALL GETNEIGHBNODES(NETP,NC,NEIGHBS2)
        Q2 = NETP%NODEPOS(:,NEIGHBS2(1)); Q3 = NETP%NODEPOS(:,NEIGHBS2(2))
        R31 = SQRT(SUM((Q3-Q1)**2)); R12 = SQRT(SUM((Q1-Q2)**2))
        LEN123 = NETP%EDGELEN(NETP%NODEEDGE(NC,1)) + NETP%EDGELEN(NETP%NODEEDGE(NC,2))
        IF(NETP%NODEDEG(NEIGHBS2(1)).EQ.1) LEN123 = LEN123+NETP%STERICRAD !get the full edge length including end cap
        IF(NETP%NODEDEG(NEIGHBS2(2)).EQ.1) LEN123 = LEN123+NETP%STERICRAD !get the full edge length including end cap
        LEN123 = 0.5*LEN123 !average of the 2 edge lengths
        EBEND = EBEND + NETP%BENDMOD1/LEN123 * (1 - DOT_PRODUCT(Q3-Q1,Q1-Q2)/(R31*R12))
        IF (GETFORCES) THEN
          FBEND(:,NC) = FBEND(:,NC) - NETP%BENDMOD1/(LEN123*R31*R12) * &
             ((Q3+Q2-2*Q1) + DOT_PRODUCT(Q3-Q1,Q1-Q2)*((Q3-Q1)/(R31**2) -(Q1-Q2)/(R12**2)))
          FBEND(:,NEIGHBS2(1)) = FBEND(:,NEIGHBS2(1)) - NETP%BENDMOD1/(LEN123*R31*R12) * &
             ((-Q3+Q1) + DOT_PRODUCT(Q3-Q1,Q1-Q2)*(Q1-Q2)/(R12**2))
          FBEND(:,NEIGHBS2(2)) = FBEND(:,NEIGHBS2(2)) - NETP%BENDMOD1/(LEN123*R31*R12) * &
             ((Q1-Q2) - DOT_PRODUCT(Q3-Q1,Q1-Q2)*(Q3-Q1)/(R31**2))
        ENDIF
     ELSEIF (DEG.EQ.3) THEN
        ! Node is degree 3. Find vectors to the connected nodes and the resulting energy
        Q1 = NETP%NODEPOS(:,NC)
        CALL GETNEIGHBNODES(NETP,NC,NEIGHBS3)
        ! Need to calculate for each pair of edges
        DO N2IDX = 1,3
           DO N3IDX = N2IDX+1,3
              Q2 = NETP%NODEPOS(:,NEIGHBS3(N2IDX)); Q3 = NETP%NODEPOS(:,NEIGHBS3(N3IDX))
              LEN123 = NETP%EDGELEN(NETP%NODEEDGE(NC,N2IDX)) + NETP%EDGELEN(NETP%NODEEDGE(NC,N3IDX))
              IF(NETP%NODEDEG(NEIGHBS3(N2IDX)).EQ.1) LEN123 = LEN123+NETP%STERICRAD !get the full edge length including end cap
              IF(NETP%NODEDEG(NEIGHBS3(N3IDX)).EQ.1) LEN123 = LEN123+NETP%STERICRAD !get the full edge length including end cap
              LEN123 = 0.5*LEN123 !average of the 2 edge lengths
              R31 = SQRT(SUM((Q3-Q1)**2)); R12 = SQRT(SUM((Q1-Q2)**2))
              COSTHETA = DOT_PRODUCT(Q3-Q1,Q1-Q2)/(R31*R12)
              IF(ABS(COSTHETA).GE.1D0.AND.ABS(COSTHETA).LT.101D-2) THEN
                 PRINT*,'COSTHETA WAS LARGER THAN 1 IN GETBENDINGFORCES',COSTHETA
                 IF(COSTHETA.GT.0D0) THEN; THETA = 0; ELSE; THETA = ACOS(-1D0); ENDIF
              ELSE
                 THETA = ACOS(COSTHETA)
              ENDIF
              EBEND = EBEND + NETP%BENDMOD2/LEN123 * (1 - COS(THETA - NETP%THETA0))
              ! get the force
              IF(GETFORCES.AND.THETA.NE.0D0) THEN
                 SINEFACTOR = SIN(THETA-NETP%THETA0)/SIN(THETA)
                 FBEND(:,NC) = FBEND(:,NC) - SINEFACTOR * NETP%BENDMOD2/(LEN123*R31*R12) * &
                    ((Q3+Q2-2*Q1) + DOT_PRODUCT(Q3-Q1,Q1-Q2)*((Q3-Q1)/(R31**2) -(Q1-Q2)/(R12**2)))
                 FBEND(:,NEIGHBS3(N2IDX)) = FBEND(:,NEIGHBS3(N2IDX)) - SINEFACTOR * NETP%BENDMOD2/(LEN123*R31*R12) * &
                    ((-Q3+Q1) + DOT_PRODUCT(Q3-Q1,Q1-Q2)*(Q1-Q2)/(R12**2))
                 FBEND(:,NEIGHBS3(N3IDX)) = FBEND(:,NEIGHBS3(N3IDX)) - SINEFACTOR * NETP%BENDMOD2/(LEN123*R31*R12) * &
                    ((Q1-Q2) - DOT_PRODUCT(Q3-Q1,Q1-Q2)*(Q3-Q1)/(R31**2))
              ENDIF

              IF(NETP%USEOUTOFPLANE) THEN
                 N4IDX = 6 - N3IDX - N2IDX
                 LEN1234 = NETP%EDGELEN(NETP%NODEEDGE(NC,N4IDX))
                 IF(NETP%NODEDEG(NEIGHBS3(N4IDX)).EQ.1) LEN1234 = LEN1234+NETP%STERICRAD
                 LEN1234 = (LEN1234 + 2*LEN123)/3 
                 Q4 = NETP%NODEPOS(:,NEIGHBS3(N4IDX))
                 R41 = SQRT(SUM((Q4-Q1)**2))
                 NORMALVEC = 0D0
                 CALL CROSS_PRODUCT(Q2-Q1,Q3-Q1,NORMALVEC)
                 NMAG = SQRT(SUM(NORMALVEC**2))
                 FLIP = 1D0
                 IF(DOT_PRODUCT(NORMALVEC,Q4-Q1).LT.0D0) THEN
                    FLIP = -1D0
                 ENDIF
                 COSSQRTHETA = 1 - (DOT_PRODUCT(NORMALVEC,Q4-Q1)/(NMAG*R41))**2
                 IF(COSSQRTHETA.LT.0D0) COSSQRTHETA = 0D0
                 EOUTOFPLANE = EOUTOFPLANE + NETP%BENDMODPLANE/LEN1234 * (1 - SQRT(COSSQRTHETA))
                 
                 IF(GETFORCES.AND.COSSQRTHETA.GT.0D0.AND.COSSQRTHETA.LT.1D0) THEN
                    FORCEFACTOR = NETP%BENDMODPLANE/LEN1234*FLIP*SQRT(1/COSSQRTHETA-1)/(NMAG*R41)**2
                    FOUTOFPLANE(:,NEIGHBS3(N4IDX)) = FOUTOFPLANE(:,NEIGHBS3(N4IDX)) + FORCEFACTOR* &
                    & (NMAG*R41*NORMALVEC - DOT_PRODUCT(NORMALVEC,Q4-Q1)*NMAG*(Q4-Q1)/R41)

                    CROSS1 = 0D0; CROSS2 = 0D0;
                    CALL CROSS_PRODUCT(Q1-Q2,Q4-Q1,CROSS1)
                    CALL CROSS_PRODUCT(Q1-Q2,NORMALVEC,CROSS2)
                    FOUTOFPLANE(:,NEIGHBS3(N3IDX)) = FOUTOFPLANE(:,NEIGHBS3(N3IDX)) + FORCEFACTOR* &
                    & (NMAG*R41*CROSS1 - DOT_PRODUCT(NORMALVEC,Q4-Q1)*R41*CROSS2/NMAG)
                    
                    CROSS1 = 0D0; CROSS2 = 0D0;
                    CALL CROSS_PRODUCT(Q3-Q1,Q4-Q1,CROSS1)
                    CALL CROSS_PRODUCT(Q3-Q1,NORMALVEC,CROSS2)
                    FOUTOFPLANE(:,NEIGHBS3(N2IDX)) = FOUTOFPLANE(:,NEIGHBS3(N2IDX)) + FORCEFACTOR* &
                    & (NMAG*R41*CROSS1 - DOT_PRODUCT(NORMALVEC,Q4-Q1)*R41*CROSS2/NMAG)

                    CROSS1 = 0D0; CROSS2 = 0D0;
                    CALL CROSS_PRODUCT(Q2-Q3,Q4-Q1,CROSS1)
                    CALL CROSS_PRODUCT(Q2-Q3,NORMALVEC,CROSS2)
                    FOUTOFPLANE(:,NC) = FOUTOFPLANE(:,NC) + FORCEFACTOR* &
                    & (NMAG*R41*(CROSS1-NORMALVEC) - DOT_PRODUCT(NORMALVEC,Q4-Q1)*(-(Q4-Q1)*NMAG/R41 + CROSS2*R41/NMAG))
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
     ELSE ! Can't have degree > 3. Should throw some error
        PRINT*, 'Node is somehow higher degree than 3. THis is bad. Degree,NC = ',DEG,NC
        STOP 1
     ENDIF
  END SUBROUTINE GETBENDINGFORCES


  SUBROUTINE DYNNETSNAPSHOT(NETP,SNAPSHOTFILE,APPEND,INFO,NNV,NODEVALS,NEV,EDGEVALS)
    ! output a snapshot of a dynamic network configuration
    ! in an efficient format (few lines per network)
    ! APPEND: append to an existing file?
    ! INFO: float array giving info for the snapshot (eg: time, other)
    ! NNV, NEV = number of node values, edge values
    ! NODEVALS(NNODE,NNV) = array of values for each node
    ! EDGEVALS(NEDGE,NEV) = array of values for each edge
    ! INFO line at top of snapshot contains: NNODEACT, NEDGEACT, NNV, NEV, NINFO,INFO
    ! max number of active nodes / edges is 10^8-1 (7 digits allowed)
    ! only active nodes and edges are dumped out
    IMPLICIT NONE

    TYPE(DYNNETWORK), POINTER :: NETP
    CHARACTER (LEN=*), INTENT(IN) :: SNAPSHOTFILE
    LOGICAL, INTENT(IN) :: APPEND
    DOUBLE PRECISION, INTENT(IN) :: INFO(:)
    INTEGER, INTENT(IN) :: NNV, NEV
    DOUBLE PRECISION, INTENT(IN) :: NODEVALS(NETP%NNODE,NNV), EDGEVALS(NETP%NEDGE,NEV)
    INTEGER, PARAMETER :: FU=51 ! File output unit
    CHARACTER (LEN=100) :: FMTINFO, NFMT, EFMT
    INTEGER :: NINFO, NNODEACT, NEDGEACT, IC
    INTEGER :: MAPALL2ACT(NETP%NNODE), MAPACT2ALL(NETP%NNODE)

    NINFO = SIZE(INFO)
    IF (NINFO.GT.9.OR.NINFO.LT.1) THEN
       PRINT*, 'ERROR IN DYNNETSNAPSHOT: must have 1 to 9 info values', NINFO
       STOP 1
    ENDIF

    NNODEACT = COUNT(NETP%NODEACT)
    NEDGEACT = COUNT(NETP%EDGEACT)
    IF (NNODEACT.GT.1D7.OR.NEDGEACT.GT.1D7) THEN
       PRINT*, 'ERROR IN DYNNETSNAPSHOT: too many nodes or edges!', NNODEACT, NEDGEACT
       STOP 1
    ENDIF

    IF (APPEND) THEN
      OPEN(UNIT=FU,FILE=SNAPSHOTFILE,STATUS='UNKNOWN',ACCESS='APPEND')
    ELSE
      OPEN(UNIT=FU,FILE=SNAPSHOTFILE,STATUS='UNKNOWN')
   ENDIF

   ! for each active node, its index in the full node list
   MAPACT2ALL = 0
   MAPACT2ALL(1:NNODEACT) = PACK((/(IC,IC=1,NETP%NNODE)/),NETP%NODEACT)
   ! for all nodes, their index in the active list
   MAPALL2ACT(MAPACT2ALL(1:NNODEACT)) = (/(IC,IC=1,NNODEACT)/)

   ! Format strings

   ! info line format
   WRITE(FMTINFO,'(A,I1,A)') '(6I8,',NINFO,'F20.10)'
   ! write info line
   WRITE(FU,FMTINFO) NETP%DIM,NNODEACT, NEDGEACT, NNV, NEV, NINFO,INFO

   ! write position coordinates for each node (one line per dim)
   WRITE(NFMT,'(A,I8,A)') '(',NNODEACT,'F20.10)'
   DO IC = 1,NETP%DIM
      WRITE(FU,NFMT) PACK(NETP%NODEPOS(IC,:),NETP%NODEACT)
   ENDDO
   ! Write additional values for each node
   DO IC = 1,NNV
      WRITE(FU,NFMT) PACK(NODEVALS(:,IC),NETP%NODEACT)
   END DO

   ! write end nodes for each edge (2 lines)
   WRITE(EFMT,'(A,I8,A)') '(',NEDGEACT,'I8)'

   DO IC = 1,2
      WRITE(FU,EFMT) MAPALL2ACT(PACK(NETP%EDGENODE(:,IC),NETP%EDGEACT))
   ENDDO
   ! write additional values for each edge
   WRITE(EFMT,'(A,I8,A)') '(',NEDGEACT,'F20.10)'
   DO IC = 1,NEV
      WRITE(FU,EFMT) PACK(EDGEVALS(:,IC),NETP%EDGEACT)
   ENDDO
   CLOSE(FU)

  END SUBROUTINE DYNNETSNAPSHOT


  SUBROUTINE DYNNETFROMFILE(NETP,NETFILE)
    ! extract structure of dynamic network from a .net file

    USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI
    USE KEYS, ONLY : MAXDEG, DIM, MAXNNODE, MAXNEDGE, NSPECIES, NPLANE
    USE GENUTIL, ONLY : NORMALIZE
    ! Set up (allocate) a network structure, reading in connectivity and
    ! geometry from an input file
    ! Reservoir labels for nodes are a holdover from signal-spreading sims
    ! In principle, nodes with the same reservoir label will maintain
    ! identical values of some quantity (eg: protein concentration)

    IMPLICIT NONE
    TYPE(DYNNETWORK), POINTER :: NETP
    CHARACTER(LEN=*), INTENT(IN) :: NETFILE
    LOGICAL :: LDUM !, CASESET
    LOGICAL :: FILEEND=.FALSE.
    CHARACTER(LEN=100) :: WORD !, LBL, STR
    INTEGER :: NITEMS, NNODE, NEDGE, NODE1, NODE2,EID, NID !, NE
    !INTEGER :: LC, WHICHLOOP,I, CC
    !INTEGER :: TMPARRAY(MAXDEG)
    LOGICAL :: EDGELENSET(MAXNNODE)
    INTEGER :: DC, SC !NC, RC,
    INTEGER, PARAMETER :: NF = 55 ! input file unit number

    ! deallocate any previously set arrays
    IF (NETP%ARRAYSET) CALL CLEANUPDYNNETWORK(NETP)
    

    ! go through file and count nodes and branches (total and from each node)
    PRINT*, 'Reading network structure file: ', NETFILE
    INQUIRE(FILE=NETFILE,EXIST=LDUM)
    IF (.NOT.LDUM) THEN
       PRINT*, 'ERROR in NETFROMFILE: network file ', TRIM(ADJUSTL(NETFILE)), ' does not exist.'
       STOP 1
    ENDIF
    OPEN(UNIT=NF, FILE=NETFILE, STATUS='OLD')

    ! ! Go through and count number of nodes and edges in the network file
    NNODE = 0; NEDGE = 0

    DO
       CALL READLINE(NF,FILEEND,NITEMS)

       IF (FILEEND.and.nitems.eq.0) EXIT
       ! skip empty lines
       IF (NITEMS.EQ.0) CYCLE
       ! Read in the keyword for this line
       CALL READA(WORD,CASESET=1)

       IF (WORD.EQ.'NODE') THEN
          NNODE = NNODE+1
       ELSEIF (WORD.EQ.'EDGE') THEN
          NEDGE = NEDGE + 1
       ENDIF
    ENDDO
    CLOSE(NF)

    PRINT*, 'Input number of nodes, edges: ', NNODE, NEDGE
    IF(NEDGE.GT.MAXNEDGE) THEN
       PRINT*,"TOO MANY INPUT EDGES"
       STOP 1
    ENDIF
    IF(NNODE.GT.MAXNNODE) THEN
      PRINT*,"TOO MANY INPUT NODES"
      STOP 1
    ENDIF

    ! allocate arrays, starting with no initial nodes or edges
    CALL SETUPDYNNETWORK(NETP,NNODE,NEDGE,DIM,MAXDEG,MAXNNODE,MAXNEDGE,NSPECIES,NPLANE)

    ! for each edge: edge length read from file?
    EDGELENSET = .FALSE.

    PRINT*, 'Reading in node positions and edge connectivity...'
    OPEN(UNIT=NF, FILE=NETFILE, STATUS='OLD')
    DO
       CALL READLINE(NF,FILEEND,NITEMS)
       IF (FILEEND.and.nitems.eq.0) EXIT ! stop reading
       ! skip empty lines
       IF (NITEMS.EQ.0) CYCLE
       ! Read in the keyword for this line
       CALL READA(WORD,CASESET=1)
       ! Skip any comment lines
       IF (WORD(1:1).EQ.'#') CYCLE

       IF (WORD.EQ.'NODE') THEN ! read in node info
          CALL READI(NID)
          IF (NID.LT.1.OR.NID.GT.NNODE) THEN
             PRINT*, 'ERROR IN NETFROMFILE: invalid  node index', NNODE, NID
             STOP 1
          ENDIF
          DO DC = 1,DIM
             CALL READF(NETP%NODEPOS(DC,NID))
          END DO
       ELSEIF (WORD.EQ.'EDGE') THEN
          CALL READI(EID) ! edge id
          CALL READI(NODE1)
          CALL READI(NODE2)
          IF (NODE1.LT.1.OR.NODE2.LT.1.OR.EID.LT.1&
               & .OR.NODE1.GT.NNODE.OR.NODE2.GT.NNODE.OR.EID.GT.NEDGE) THEN
             PRINT*, 'ERROR IN MINNET FROM FILE: &
                  & bad edge or node indices while reading edge.', &
                  ' EID, NODE1, NODE2, NNODE, NEDGE:', &
                  & EID, NODE1, NODE2, NNODE, NEDGE
             STOP 1
          ENDIF
          IF (NITEMS.GT.4) THEN ! read edge length from file
             CALL READF(NETP%EDGELEN(EID))
             EDGELENSET(EID) = .TRUE.
             DO SC = 1,MIN(NSPECIES,NITEMS-5)
                CALL READF(NETP%CONC(EID,SC))
             ENDDO
          ENDIF
          ! nodes connected to this edge
          NETP%EDGENODE(EID,:) = (/NODE1,NODE2/)
       ENDIF
    END DO
    CLOSE(NF)
    DO EID = 1,NEDGE
       NODE1 = NETP%EDGENODE(EID,1)
       NODE2 = NETP%EDGENODE(EID,2)

       IF (.NOT.EDGELENSET(EID)) THEN
          NETP%EDGELEN(EID) = -1D0
       ENDIF


       ! increment degrees of the nodes
       NETP%NODEDEG(NODE1) = NETP%NODEDEG(NODE1)+1
       NETP%NODEDEG(NODE2) = NETP%NODEDEG(NODE2)+1

       IF ( MAX(NETP%NODEDEG(NODE1),NETP%NODEDEG(NODE2)).GT.MAXDEG) THEN
          PRINT*, 'ERROR IN MINNETFROMFILE: node degree exceeds maximum.',&
               & NODE1, NODE2, MAXDEG, NETP%NODEDEG(NODE1), NETP%NODEDEG(NODE2)
          STOP 1
       ENDIF

       ! edges connected to each node
       NETP%NODEEDGE(NODE1,NETP%NODEDEG(NODE1)) = EID
       NETP%NODEEDGEIND(NODE1,NETP%NODEDEG(NODE1)) = 1
       NETP%NODEEDGE(NODE2,NETP%NODEDEG(NODE2)) = EID
       NETP%NODEEDGEIND(NODE2,NETP%NODEDEG(NODE2)) = 2

       ! where does this edge fall among the list of edges for each node?
       NETP%EDGENODEIND(EID,1) = NETP%NODEDEG(NODE1)
       NETP%EDGENODEIND(EID,2) = NETP%NODEDEG(NODE2)
    END DO
    ! non-empty structure set for the network?
    NETP%STRUCTURESET = (NNODE.GT.0)

  END SUBROUTINE DYNNETFROMFILE

  SUBROUTINE DYNNETFROMRANDOMFRAGMENTS(NETP)
   ! set up a dynamic network of randomly scattered edges in a sphere

   USE KEYS, ONLY : MAXDEG, DIM, MAXNNODE, MAXNEDGE, NSPECIES, NPLANE, CELLRAD1, MITOLEN, STERICRAD
   USE GENUTIL, ONLY : NORMALIZE, RANDOMSPHEREPOINT
   ! Set up (allocate) a network structure of randomly scattered edges in a sphere
   ! Reservoir labels for nodes are a holdover from signal-spreading sims
   ! In principle, nodes with the same reservoir label will maintain
   ! identical values of some quantity (eg: protein concentration)

   IMPLICIT NONE
   TYPE(DYNNETWORK), POINTER :: NETP
   INTEGER :: NNODE, NEDGE, EC, NC
   DOUBLE PRECISION, DIMENSION(DIM) :: POS, DIR

   ! deallocate any previously set arrays
   IF (NETP%ARRAYSET) CALL CLEANUPDYNNETWORK(NETP)
   
   ! Set the number of nodes and edges based on MAX values
   NNODE = MAXNNODE; NEDGE = MAXNEDGE
   IF(NNODE.NE.2*NEDGE) THEN
      PRINT*,"INCOMPATIBLE MAXNEDGE, MAXNNODE in network creation",MAXNEDGE,MAXNNODE
      STOP 1
   ENDIF

   PRINT*,"Creating a random fragmented network"
   PRINT*, 'Input number of nodes, edges: ', NNODE, NEDGE

   ! allocate arrays, starting with no initial nodes or edges
   CALL SETUPDYNNETWORK(NETP,NNODE,NEDGE,DIM,MAXDEG,MAXNNODE,MAXNEDGE,NSPECIES,NPLANE)

   PRINT*, 'Setting node positions and edge connectivity...'
   
   ! shorten all edge lengths to the inset value
   NETP%EDGELEN = MITOLEN-2*STERICRAD
   
   ! start node counter
   NC = 1
   
   ! loop through the edges
   DO EC = 1,NEDGE
      ! all edges will be individual units so all nodes are degree-1
      NETP%NODEDEG(NC) = 1
      NETP%NODEDEG(NC+1) = 1
      ! place the edge randomly in space, uniform through the spherical volume
      CALL RANDOMSPHEREPOINT(CELLRAD1,.FALSE.,POS)
      CALL RANDOMSPHEREPOINT(1D0,.TRUE.,DIR)
      NETP%NODEPOS(:,NC) = POS + NETP%EDGELEN(EC)*DIR/2
      NETP%NODEPOS(:,NC+1) = POS - NETP%EDGELEN(EC)*DIR/2
      ! nodes connected to this edge
      NETP%EDGENODE(EC,1) = NC
      NETP%EDGENODE(EC,2) = NC+1
      ! edge which is connected to each node
      NETP%NODEEDGE(NC,1) = EC
      NETP%NODEEDGE(NC+1,1) = EC
      ! edges connected to each node
      NETP%NODEEDGEIND(NC,1) = 1
      NETP%NODEEDGEIND(NC+1,1) = 2
      ! where does this edge fall among the list of edges for each node?
      NETP%EDGENODEIND(EC,1) = 1
      NETP%EDGENODEIND(EC,2) = 1

      ! increment the node counter
      NC = NC+2
   ENDDO

   ! non-empty structure set for the network?
   NETP%STRUCTURESET = (NNODE.GT.0)

 END SUBROUTINE DYNNETFROMRANDOMFRAGMENTS

  SUBROUTINE SETPLANESFROMFILE(NETP,PLANEFILE)
     USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI
     USE GENUTIL, ONLY : CROSS_PRODUCT

     IMPLICIT NONE
     TYPE(DYNNETWORK), POINTER :: NETP
     CHARACTER(LEN=*), INTENT(IN) :: PLANEFILE
     LOGICAL :: LDUM
     LOGICAL :: FILEEND=.FALSE.
     CHARACTER(LEN=100) :: WORD
     INTEGER :: NITEMS
     INTEGER :: DC, PID
     INTEGER, PARAMETER :: PF = 55 ! input file unit number

     PRINT*, 'Reading plane file: ', PLANEFILE
     INQUIRE(FILE=PLANEFILE,EXIST=LDUM)
     IF (.NOT.LDUM) THEN
        PRINT*, 'ERROR in SETPLANESFROMFILE: plane file ', TRIM(ADJUSTL(PLANEFILE)), ' does not exist.'
        STOP 1
     ENDIF
     OPEN(UNIT=PF, FILE=PLANEFILE, STATUS='OLD')

     DO
        CALL READLINE(PF,FILEEND,NITEMS)

        IF (FILEEND.AND.NITEMS.EQ.0) EXIT
        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE
        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)
        ! Skip any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE

        CALL READI(PID)
        IF (PID.LT.1.OR.PID.GT.NETP%NPLANE) THEN
           PRINT*, 'ERROR IN SETPLANESFROMFILE: invalid  plane index', NETP%NPLANE, PID
           STOP 1
        ENDIF

        IF (WORD.EQ.'VERTEX') THEN
           DO DC = 1,NETP%DIM
              CALL READF(NETP%PLANEVERTS(DC,PID))
           END DO
        ELSEIF (WORD.EQ.'DIR1') THEN
           DO DC = 1,NETP%DIM
              CALL READF(NETP%PLANEEDGES(DC,1,PID))
           END DO
        ELSEIF (WORD.EQ.'DIR2') THEN
           DO DC = 1,NETP%DIM
              CALL READF(NETP%PLANEEDGES(DC,2,PID))
           END DO
        ELSE
           PRINT*,"ERROR in SETPLANESFROMFILE, invalid keyword ",WORD
           STOP 1
        ENDIF
     ENDDO
     CLOSE(PF)

     DO PID = 1,NETP%NPLANE
      ! calculate side lengths and normal from edge vectors
      CALL CROSS_PRODUCT(NETP%PLANEEDGES(:,1,PID),NETP%PLANEEDGES(:,2,PID),NETP%PLANENORMALS(:,PID))
      NETP%PLANESIZES(1,PID) = NORM(NETP%PLANEEDGES(:,1,PID))
      NETP%PLANESIZES(2,PID) = NORM(NETP%PLANEEDGES(:,2,PID))
      ! now make the direction vectors into unit length
      NETP%PLANEEDGES(:,1,PID) = NETP%PLANEEDGES(:,1,PID)/NETP%PLANESIZES(1,PID)
      NETP%PLANEEDGES(:,2,PID) = NETP%PLANEEDGES(:,2,PID)/NETP%PLANESIZES(2,PID)
      NETP%PLANENORMALS(:,PID) = NETP%PLANENORMALS(:,PID)/NORM(NETP%PLANENORMALS(:,PID))
      PRINT*,"PLANE NORMAL",NETP%PLANENORMALS(:,PID)
      PRINT*,"PLANE VERTEX",NETP%PLANEVERTS(:,PID)
      PRINT*,"PLANEEDGE1",NETP%PLANESIZES(1,PID), NETP%PLANEEDGES(:,1,PID)
      PRINT*,"PLANEEDGE2",NETP%PLANESIZES(2,PID), NETP%PLANEEDGES(:,2,PID)
     ENDDO

  ENDSUBROUTINE SETPLANESFROMFILE

  SUBROUTINE GETNEIGHBNODES(NETP,NC,NEIGHBS)
    ! list of neighbor nodes for a given node
    ! in same order as NODEEDGE
    IMPLICIT NONE
    TYPE(DYNNETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: NC
    INTEGER, INTENT(OUT) :: NEIGHBS(NETP%NODEDEG(NC))
    INTEGER :: EI, EI2, EC, ECC

    DO ECC= 1,NETP%NODEDEG(NC)
       EI = NETP%NODEEDGEIND(NC,ECC) ! this is 1st or 2nd node for this edge?
       EI2 = 3 - EI ! is neighbor node 1st or 2nd for this edge?

       ! edge index
       EC = NETP%NODEEDGE(NC,ECC)
       ! other node
       NEIGHBS(ECC) = NETP%EDGENODE(EC,EI2)
    ENDDO

  END SUBROUTINE GETNEIGHBNODES

  SUBROUTINE SETNETPARAM(NETP)
    ! set up parameters for network from global keys values
    USE KEYS, ONLY : ESTR, BENDMOD1, BENDMOD2, THETA0, STERICMOD, STERICRAD, FISSRAD, D3FISSMULT, &
      KT, FRICT, FISSRATE, FUSERATE1, FUSERATE2, CONTACTRAD, ECONF, CELLRAD1, &
      GRIDSPACES, USEGRID, NSPECIES, DIFF, DOYEAST, YEASTCONF, YEASTONRATE, YEASTOFFRATE, &
      YEASTBINDRANGE, MITOLEN, PRODRATE, ALPHA1, ALPHA2, NUMTRIGGERUNITS, &
      RECHARGERATE, USEEDGEGRID, USENODEGRID, DIM, DECAYRATE, FIXEDGES, BENDMODPLANE, USEOUTOFPLANE, &
      ALPHAPLANE, UNIQUETRIGGERS

    IMPLICIT NONE
    INTEGER :: EC
    TYPE(DYNNETWORK), POINTER :: NETP

    NETP%ESTR = ESTR
    NETP%BENDMOD1 = BENDMOD1
    NETP%BENDMOD2 = BENDMOD2
    NETP%BENDMODPLANE = BENDMODPLANE
    NETP%THETA0 = THETA0
    NETP%STERICRAD = STERICRAD
    NETP%STERICMOD = STERICMOD
    NETP%KT = KT
    NETP%FRICT = FRICT
    NETP%FISSRATE = FISSRATE
    NETP%FUSERATE1 = FUSERATE1
    NETP%FUSERATE2 = FUSERATE2
    NETP%CONTACTRAD = CONTACTRAD
    NETP%CELLRAD1 = CELLRAD1
    NETP%ECONF = ECONF
    NETP%GRIDSPACES = GRIDSPACES(1:DIM)
    NETP%USEGRID = USEGRID
    NETP%USEEDGEGRID = USEEDGEGRID
    NETP%USENODEGRID = USENODEGRID
    NETP%NSPECIES = NSPECIES
    NETP%DIFF = DIFF(1:NSPECIES)
    NETP%PRODRATE = PRODRATE(1:NSPECIES)
    NETP%DECAYRATE = DECAYRATE(1:NSPECIES)
    NETP%DOYEAST = DOYEAST
    NETP%YEASTCONF = YEASTCONF
    NETP%YEASTONRATE = YEASTONRATE
    NETP%YEASTOFFRATE = YEASTOFFRATE
    NETP%YEASTBINDRANGE = YEASTBINDRANGE
    NETP%MITOLEN = MITOLEN
    NETP%FISSRAD = FISSRAD
    NETP%D3FISSMULT = D3FISSMULT
    NETP%ALPHA1 = ALPHA1
    NETP%ALPHA2 = ALPHA2
    NETP%NUMTRIGGERUNITS = NUMTRIGGERUNITS
    NETP%RECHARGERATE = RECHARGERATE
    NETP%FIXEDGES = FIXEDGES
    NETP%USEOUTOFPLANE = USEOUTOFPLANE
    NETP%ALPHAPLANE = ALPHAPLANE
    NETP%UNIQUETRIGGERS = UNIQUETRIGGERS

    DO EC = 1,NETP%NEDGE
      IF(.NOT.NETP%EDGEACT(EC)) CYCLE
      IF(NETP%EDGELEN(EC).LT.0D0) THEN
         NETP%EDGELEN(EC) = MITOLEN
         IF(NETP%NODEDEG(NETP%EDGENODE(EC,1)).EQ.1) NETP%EDGELEN(EC) = NETP%EDGELEN(EC) - STERICRAD
         IF(NETP%NODEDEG(NETP%EDGENODE(EC,2)).EQ.1) NETP%EDGELEN(EC) = NETP%EDGELEN(EC) - STERICRAD
      ENDIF
    ENDDO

  END SUBROUTINE SETNETPARAM


  SUBROUTINE SETUPDYNNETWORK(NETP,NNODE,NEDGE,DIM,MAXDEG,MAXNNODE,MAXNEDGE,NSPECIES,NPLANE)
    ! set up a mitochondrial network by allocating arrays
    ! MITOP: Pointer to a network
    ! NNODE: number of nodes
    ! NEDGE: number of branches
    ! DIM: spatial dimensionality where network resides
    ! MAXDEG: maximum allowed degree of a node
    ! MAXNNODE: max allowed nodes
    ! MAXNEDGE: max allowed edges
    ! PREFILLSTACK: also set up the stack saving indices for empty slots that can be filled with new nodes / edges
    USE STACKUTIL, ONLY : INITIALIZESTACK
    IMPLICIT NONE
    TYPE(DYNNETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: NNODE, NEDGE, DIM, MAXDEG, MAXNNODE,MAXNEDGE,NSPECIES,NPLANE
    INTEGER :: LN(MAXNNODE-NNODE), LE(MAXNEDGE-NEDGE), I

    ! all network variables initialized with space for extra nodes/edges
    NETP%NNODE = MAXNNODE
    NETP%NEDGE = MAXNEDGE
    NETP%DIM = DIM

    ! allocate node data
    ALLOCATE(NETP%NODEEDGE(MAXNNODE,MAXDEG), NETP%NODEEDGEIND(MAXNNODE,MAXDEG))
    ALLOCATE(NETP%NODEPOS(DIM,MAXNNODE), NETP%NODEDEG(MAXNNODE),&
         & NETP%NODEACT(MAXNNODE), NETP%GRID(DIM,MAXNNODE), NETP%GRIDSPACES(DIM))

    ! allocate branch data
    ALLOCATE(NETP%EDGENODE(MAXNEDGE,2), NETP%EDGENODEIND(MAXNEDGE,2), NETP%EDGELEN(MAXNEDGE),&
         & NETP%EDGEACT(MAXNEDGE))

    !allocate diffusion data
    ALLOCATE(NETP%DIFF(NSPECIES), NETP%CONC(MAXNEDGE,NSPECIES), &
         & NETP%PRODRATE(NSPECIES),NETP%DECAYRATE(NSPECIES))

    !allocate yeast data
    ALLOCATE(NETP%TETHERED(MAXNNODE))

    !allocate planes data
    ALLOCATE(NETP%PLANEVERTS(DIM,NPLANE),NETP%PLANEEDGES(DIM,2,NPLANE),NETP%PLANENORMALS(DIM,NPLANE),NETP%PLANESIZES(2,NPLANE))

    ! allocate recharging state
    ALLOCATE(NETP%FUSEACTIVE(MAXNNODE))

    NETP%ARRAYSET = .TRUE.
    NETP%NODEDEG = 0
    NETP%NODEEDGE = 0; NETP%EDGENODE = 0
    NETP%NODEEDGEIND = 0; NETP%EDGENODEIND = 0
    NETP%NODEACT = .FALSE.
    NETP%EDGEACT = .FALSE.
    NETP%GRID = 0
    NETP%CONC = 0D0
    NETP%TETHERED = .FALSE.
    NETP%PLANEVERTS = 0D0; NETP%PLANEEDGES = 0D0; NETP%PLANENORMALS = 0D0; NETP%PLANESIZES = 0D0
    NETP%NPLANE = NPLANE
    NETP%FUSEACTIVE = .TRUE.

    NETP%NODEACT(1:NNODE) = .TRUE.
    NETP%EDGEACT(1:NEDGE) = .TRUE.

    ! allocate the pointers for stacks
    ALLOCATE(NETP%NODESTACK, NETP%EDGESTACK)

    ! set up the stacks, prefilled with unused nodes/edges
    LN = [(MAXNNODE + NNODE + 1 - I, I=(NNODE+1), MAXNNODE)]
    CALL INITIALIZESTACK(NETP%NODESTACK,MAXNNODE,LN)
    LE = [(MAXNEDGE + NEDGE + 1 - I, I=(NEDGE+1), MAXNEDGE)]
    CALL INITIALIZESTACK(NETP%EDGESTACK,MAXNEDGE,LE)

  END SUBROUTINE SETUPDYNNETWORK

  SUBROUTINE CLEANUPDYNNETWORK(NETP)
    ! deallocate arrays for the network structure
    USE STACKUTIL, ONLY : DEALLOCATESTACK
    IMPLICIT NONE
    TYPE(DYNNETWORK), POINTER :: NETP

    IF (NETP%ARRAYSET) THEN
       DEALLOCATE(NETP%NODEEDGEIND, NETP%NODEEDGE, NETP%NODEPOS, &
            & NETP%NODEDEG, NETP%NODEACT, NETP%GRID, NETP%GRIDSPACES)
       DEALLOCATE(NETP%EDGENODE, NETP%EDGENODEIND,NETP%EDGEACT, NETP%EDGELEN)
       DEALLOCATE(NETP%DIFF, NETP%CONC)
       DEALLOCATE(NETP%TETHERED, NETP%PRODRATE, NETP%DECAYRATE)
       DEALLOCATE(NETP%PLANEVERTS,NETP%PLANEEDGES,NETP%PLANENORMALS,NETP%PLANESIZES)
       DEALLOCATE(NETP%FUSEACTIVE)
       ! deallocate stack variables
       CALL DEALLOCATESTACK(NETP%NODESTACK)
       CALL DEALLOCATESTACK(NETP%EDGESTACK)
       DEALLOCATE(NETP%NODESTACK, NETP%EDGESTACK)
    ENDIF

    NETP%ARRAYSET = .FALSE.
    NETP%STRUCTURESET = .FALSE.
  END SUBROUTINE CLEANUPDYNNETWORK

END MODULE DYNNETWORKUTIL
