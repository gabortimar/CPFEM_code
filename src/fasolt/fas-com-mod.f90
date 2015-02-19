!***********************************************************************************!
!																					!
!								MODULE FAS_COM										!
!-----------------------------------------------------------------------------------!
!	Module containing the type and array declarations for the main FASOLT variables !
!-----------------------------------------------------------------------------------!
!																					!
!***********************************************************************************!
MODULE FAS_COM

	IMPLICIT NONE

	INTEGER, PARAMETER:: maxnod= 65000, maxdof= 180000,  maxels= 8100 
	INTEGER, PARAMETER:: maxgps= 65000, maxinc= 10000,   maxsvs= 63
	INTEGER, PARAMETER:: maxofd= 1000

	INTEGER(4)	nlmnt, nnod, ngps, ninc, ninc0, ninc1, nsv, inc, nxs, nfs, maxit, tflag

	INTEGER(4), DIMENSION( 1:maxels)	   :: le, lt, nd, np
	INTEGER(4), DIMENSION( 1:20, 1:maxels) :: ln
	INTEGER(4), DIMENSION(  1:8, 1:maxels) :: lg
	INTEGER(4), DIMENSION( 1:maxnod)	   :: ib	
	REAL(4),	DIMENSION( 1:maxdof)	   :: x, dx, fc, ddx
	REAL(4),	DIMENSION( 1:maxgps)	   :: eps, dep
	REAL(4),	DIMENSION( 1:6, 1:maxgps)  :: st, ep, ds, de
	REAL(4),	DIMENSION( 1:3, 1:maxgps)  :: rt, dr
	REAL(4),	DIMENSION( 1:maxsvs, 1:maxgps)	:: sv, dsv

	REAL(4),	DIMENSION( 1:maxinc)	:: tiem
	REAL(4),	DIMENSION( 1:15, 1:99, 1:maxinc)	:: x_spec
	REAL(4),	DIMENSION( 1: 3, 1:99, 1:maxinc)	:: f_spec

	INTEGER(4),	DIMENSION( 1: maxofd, 1:maxdof)	:: is
	REAL(4),    DIMENSION( 1: maxofd, 1:maxdof)	:: ss

!-------------------------------------surface data
	REAL(4)	 	amu, anu, aga

!-------------------------------------numerical control data
	REAL(4)		toler, alpha, efix, smalp, smalf, smald, smale, smalr, afactor

                                                                                                        
END MODULE FAS_COM

!***********************************************************************************!
!																					!
!								MODULE FAS_CPL										!
!-----------------------------------------------------------------------------------!
!			Module containing type and array declarations for crystal				!
!				plasticity variables and material property data						!
!-----------------------------------------------------------------------------------!
!																					!
!***********************************************************************************!
MODULE FAS_CPL

	IMPLICIT NONE
	INTEGER(4), PARAMETER :: maxss=60, maxmds= 12, mgps=65000  ! this is the same as maxgps, only did this to avoid conflict

	INTEGER(4)	nmpar, nmatl, mat_code, nsa, twinlim
	REAL(4)		twinprob

	INTEGER(4), DIMENSION( 1:maxss)	::  isa
	REAL(4),    DIMENSION( 1:maxss)	::  tau, gam, wks, s0, s1
	REAL(4),    DIMENSION( 1:6, 1:maxss)  :: b		! slip rate to strain rate matrix
	REAL(4),	DIMENSION( 1:3, 1:maxss)  :: rot	! slip rate to global spin
	REAL(4),	DIMENSION( 1:maxss, 1:maxss):: ah
	
	INTEGER(4), DIMENSION( 1: maxmds):: nss
	REAL(4),	DIMENSION( 1:3, 1:maxss, maxmds)::		ab, an	! slip direction and slip plane normal vector components
	REAL(4),	DIMENSION( 1:maxss, 1:maxss , maxmds):: smatrx

	REAL(4),    	DIMENSION( 1:9, 1:maxmds)::				e
	INTEGER(4),	DIMENSION(1:mgps)::		timestwinned	! number of times IP has twinned
	INTEGER(4), 	DIMENSION(1:mgps):: 		twinsys		! indicate system for twinning
	REAL(4), 	DIMENSION( 1:maxss, 1:maxmds):: ssang	! angle of reorientation for twin system (0. if slip)
	REAL(4),	DIMENSION(1:maxss, 1:mgps):: srate	! slip rates to be written to .slp file
	REAL(4),	DIMENSION(1:maxss, 1:mgps):: cumslip	! cumulative slip activity to be written to .slp file
	REAL(4),	DIMENSION(1:10, 1:maxss, 1:maxmds)::	ampar	

END MODULE FAS_CPL
