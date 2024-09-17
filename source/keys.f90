MODULE KEYS
  ! keyword parameters that are globally used in many different places in the code
  IMPLICIT NONE

   ! -------- General program control ---------------
  CHARACTER*100 :: ACTION
  INTEGER :: RNGSEED
  LOGICAL :: VERBOSE

  ! ----------------------
  ! mito network geometry, mechanics, diffusive spread/transport
  ! ---------------------
  INTEGER, PARAMETER :: MAXDEG = 4 ! max allowed node degree
  INTEGER, PARAMETER :: MAXSPECIES = 9 ! maximum number of diffusing species
  INTEGER :: MAXNNODE,MAXNEDGE,NSPECIES
  ! network dimension
  INTEGER :: DIM
  ! for each bond, ESTR: stretch energy prefactor, LS: ground-state length
  ! LP: persistence length (bend resistance)
  DOUBLE PRECISION :: MITOLEN, ESTR
  ! confinement (radius, prefactor)
  DOUBLE PRECISION :: CELLRAD1, ECONF
  ! number of planes if confining using planes
  INTEGER :: NPLANE
  ! Bending params: THETA0 = desired bending angle for degree 3 intersections
  ! BENDMOD1,2 = Bending energy prefactors for degree 2,3 nodes
  DOUBLE PRECISION :: THETA0, BENDMOD1, BENDMOD2

  ! Whether to include an explicit out-of-plane energy. Associated bending energy prefactor
  LOGICAL :: USEOUTOFPLANE
  DOUBLE PRECISION :: BENDMODPLANE

  ! ----------------------
  ! Steric forcing parameters. 
  ! ---------------------
  ! STERICMOD gives the energy prefactor (like spring constant) and STERICRAD gives 1/2 the 
  ! distance at which we begin caring about the interaction (close range force only)
  DOUBLE PRECISION :: STERICMOD, STERICRAD
  !whether or not to use the GRID, edgegrid, nodegrid
  LOGICAL :: USEGRID, USEEDGEGRID, USENODEGRID
  ! spacing for grid
  INTEGER :: GRIDSPACES(MAXDEG)


  ! ----------------------
  ! material spread, ER contact sites
  ! ---------------------
  ! diffusion coefficients
  DOUBLE PRECISION :: DIFF(MAXSPECIES)
  ! production rates
  DOUBLE PRECISION :: PRODRATE(MAXSPECIES)
  ! decay rates
  DOUBLE PRECISION :: DECAYRATE(MAXSPECIES)
  ! number of substeps to use for FVM diffusion calculation
  INTEGER :: NSUBSTEPS
  ! number of edges which produce secondary species
  INTEGER :: NUMTRIGGERUNITS
  ! whether to fix the trigger units at their starting concentration value
  LOGICAL :: FIXEDGES
  ! whether to use a unique edge to produce the signal for each species
  LOGICAL :: UNIQUETRIGGERS


  ! ------------------
  ! yeast model
  ! ------------------
  LOGICAL :: DOYEAST  ! whether to use yeast model, budding, polar tethering
  ! Forcing Strength of the yeast tethering
  DOUBLE PRECISION :: YEASTCONF
  ! Rates to tether onto and off of the cell boundary, range at which tethering can occur
  DOUBLE PRECISION :: YEASTONRATE, YEASTOFFRATE, YEASTBINDRANGE
  

  ! -----------------
  ! fission and fusion
  ! ----------------
  ! fission rate (per time per fissable node) and radius to separate the nodes upon fission (total separation becomes 2x this)
  ! also multiplier for D3 fissions
  DOUBLE PRECISION :: FISSRATE, FISSRAD, D3FISSMULT
  ! fusion rate for deg 2, 3 junctions respectively
  DOUBLE PRECISION :: FUSERATE1, FUSERATE2
  ! angular sensitivity for fusion
  DOUBLE PRECISION :: ALPHA1, ALPHA2, ALPHAPLANE
  ! spherical cap radius size for mitochondria tips
  DOUBLE PRECISION :: CONTACTRAD
  ! rate to recharge fusion machinery after fission
  DOUBLE PRECISION :: RECHARGERATE
  
  ! ----------------------
  ! Output / input
  ! -----------------------
  CHARACTER*100 :: OUTFILE, SNAPSHOTFILE, NETFILE, PLANEFILE, REMODELINGFILE
  LOGICAL :: DUMPSNAPSHOTS, RESTART, APPENDSNAPSHOTS
  INTEGER :: SNAPSHOTEVERY, PRINTEVERY
  LOGICAL :: USERANDOMFRAGMENTEDNETWORK ! whether to start from a random fragmented network

  ! When to start propagating diffusion and saving snapshots
  INTEGER :: STARTSNAPSTEP, STARTDIFFSTEP
  
  ! -----------------
  ! brownian dynamics
  ! -----------------
  DOUBLE PRECISION :: FRICT, DELT
  INTEGER :: BDSTEPS
  LOGICAL :: DOBROWN
  DOUBLE PRECISION :: KT


  
END MODULE KEYS
