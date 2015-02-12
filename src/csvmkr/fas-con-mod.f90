!***********************************************************************************!
!																					!
!								MODULE FAS_COM										!
!-----------------------------------------------------------------------------------!
!	Module containing the type and array declarations for the main FASOLT variables !
!-----------------------------------------------------------------------------------!
!																					!
!***********************************************************************************!
MODULE FAS_COM


	INTEGER, PARAMETER:: maxnod= 50000, maxdof=150000,  maxels= 8000 
	INTEGER, PARAMETER:: maxgps=64000, maxinc= 2000,   maxsvs= 63
	INTEGER, PARAMETER:: maxofd= 1000

	INTEGER(4)	nlmnt, nnod, ngps, ninc, ninc0, ninc1, nsv, inc, nxs, nfs, nmpar, maxit

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

!------------------------------------------elastic contants
	INTEGER(4)  nmatl
	REAL(4), DIMENSION(1:9,1:20)::e

!------------------------------------------plot entities
	INTEGER(4)	n_point, n_ent
	INTEGER(4), DIMENSION(0:3, 1:10*maxnod)::	entlist 
	REAL(4),	DIMENSION(1:3, 1: 2*maxnod)::	point

!-------------------------------------------------restrictions
	INTEGER(4)	nelclass
	LOGICAL(4),	DIMENSION(1:maxels):: elrestrict
	LOGICAL(4),	DIMENSION(1:maxnod):: ndrestrict
	INTEGER(4),	DIMENSION(1:maxels):: elclass
 
!--------------------------------------------'contouring' 
	CHARACTER	clable*30
	INTEGER(4)	l_contour
	REAL(4)		gpv_zero, gpv_scale
	REAL(4),	DIMENSION(1:maxgps):: gpvalue 
	REAL(4),	DIMENSION(1:maxnod)::nod_cv


END MODULE FAS_COM
