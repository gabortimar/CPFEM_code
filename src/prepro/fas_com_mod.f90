!***********************************************************************************!
!																					!
!								MODULE FAS_COM										!
!-----------------------------------------------------------------------------------!
!	Module containing the type and array declarations for the main FASOLT variables !
!-----------------------------------------------------------------------------------!
! maxnod = maximum number of nodes (at least 8*maxels)								!
! maxels = maximum number of elements												!
! maxgps = maximum number of gauss points (st least 8*maxels)						!
! maxinc = maximum number of increments												!
! maxsvs = maximum number saved variables											!
! maxdof = ?																		!
!																					!
!***********************************************************************************!
MODULE FAS_COM

	!INTEGER, PARAMETER:: maxnod= 50000, maxdof=150000,  maxels= 6000 
	INTEGER, PARAMETER:: maxnod= 65000, maxdof= 200000,  maxels= 8500 
	INTEGER, PARAMETER:: maxgps= 65000, maxinc= 2000,   maxsvs= 63


	INTEGER(4)	nlmnt, nnod, ngps, ninc, ninc0, ninc1, nsv, inc, nxs, nfs, maxit

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
	REAL(4),	DIMENSION( 1:3,  1:99, 1:maxinc)	:: f_spec

!------------------------------------------------------------Material data
	INTEGER(4), PARAMETER :: maxss=60, maxmds= 12

	INTEGER(4)	nmpar, nmatl, mat_code, twinlim
	REAL(4)		twinprob

	INTEGER(4), DIMENSION( 1:maxss)	::  isa
	REAL(4),    DIMENSION( 1:maxss)	::  tau, gam, wks, s0, s1
	REAL(4),    DIMENSION( 1:6, 1:maxss)  :: b
	REAL(4),	DIMENSION( 1:3, 1:maxss)  :: rot
	REAL(4),	DIMENSION( 1:maxss, 1:maxss):: ah
	
	INTEGER(4), DIMENSION( 1: maxmds):: nss
	REAL(4),	DIMENSION( 1:3, 1:maxss, maxmds)::		ab, an
	REAL(4),	DIMENSION( 1:maxss, 1:maxss , maxmds):: smatrx

	REAL(4),    DIMENSION( 1:9, 1:maxmds)::				e
	REAL(4), DIMENSION( 1:maxss, 1:maxmds):: 	ssang
	REAL(4),	DIMENSION(1:10, 1:maxss, 1:maxmds)::	ampar	

!--------------------------------------surface data
	REAL(4)		 amu, anu, aga


END MODULE FAS_COM
