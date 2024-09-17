PROGRAM MAIN
  ! subroutines for testing parts of the code
  USE KEYS, ONLY : ACTION, DELT, BDSTEPS, NETFILE, PLANEFILE, NPLANE, USERANDOMFRAGMENTEDNETWORK
  USE DYNNETWORKUTIL, ONLY : DYNNETWORK, DYNNETFROMFILE,SETNETPARAM, CLEANUPDYNNETWORK, SETPLANESFROMFILE, DYNNETFROMRANDOMFRAGMENTS
  USE BROWNDYN, ONLY : RUNBROWNDYNSIM, SETDYNPARAM, DYNPARAM

  IMPLICIT NONE
  TYPE(DYNNETWORK), TARGET :: NET
  TYPE(DYNNETWORK), POINTER :: NETP
  TYPE(DYNPARAM), TARGET :: DYNOBJ
  TYPE(DYNPARAM), POINTER :: DP
  INTEGER :: START, FINISH, CLOCKRATE

  CLOCKRATE = 1
  CALL SYSTEM_CLOCK(START,CLOCKRATE)

  NETP=>NET
  NETP%ARRAYSET = .FALSE.
  
  ! read in all keywork parameters
  CALL READKEY

  ! set up initial network (from file) and associated parameters
  IF(USERANDOMFRAGMENTEDNETWORK) THEN
     CALL DYNNETFROMRANDOMFRAGMENTS(NETP)
  ELSE
     CALL DYNNETFROMFILE(NETP,NETFILE)
  ENDIF
  IF(NPLANE.GT.0) CALL SETPLANESFROMFILE(NETP,PLANEFILE)
  CALL SETNETPARAM(NETP)
  ! set parameters for running dynamic sims
  DP=>DYNOBJ
  CALL SETDYNPARAM(DP)
  
  
  IF (ACTION.EQ.'RUNDYNAMICS') THEN
     CALL RUNBROWNDYNSIM(NETP, BDSTEPS, DELT, DP, START, CLOCKRATE)
  ELSE
     PRINT*, 'ERROR IN MAIN: unidentified action', ACTION
     STOP 1
  ENDIF

  CALL CLEANUPDYNNETWORK(NETP)

  CALL SYSTEM_CLOCK(FINISH,CLOCKRATE)
  PRINT*,"ELAPSED TIME (SECONDS): ",(FINISH-START)/1000

END PROGRAM MAIN
