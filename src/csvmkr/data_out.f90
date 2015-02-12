!*******************************************************************************!
!																				!
!   SUBROUTINE GPSTATE															!
! 		Outputs a list (.CSV) of weights and values at integration points.		!
!																				!
!*******************************************************************************!
SUBROUTINE GPSTATE(l_dummy)

	USE FAS_COM

	IMPLICIT NONE
	LOGICAL(4)	l_dummy, l_file
	INTEGER(4)	ilmnt, iout, iel, itp, nn, ng, in, ip, i, j, ks, ig
	REAL(4)		xl, weight, xg, vm_strs, VONM_STRESS

	REAL(4), DIMENSION(1:6):: stress, strain

	DIMENSION	in(1:20), ip(1:8), xl(1:3, 1:20), xg(1:3)

	PARAMETER( iout=53)

!----open file here---------------!!!!!!!
        
        
    CALL OUTFILE(iout,l_file)


!-------------------------------------write header
		WRITE( iout,'(''FAFNER3 Gauss point state values'')')
		WRITE( iout,'(''int. point,volume, x, y, z, &
			&stress(11),stress(22),stress(33),stress(23),stress(31),stress(12),&
			&inc. s(11),inc. s(22),inc. s(33),inc. s(23),inc. s(31),inc. s(12),&
			&strain(11),strain(22),strain(33),strain(23),strain(31),strain(12),&
			&inc. e(11),inc. e(22),inc. e(33),inc. e(23),inc. e(31),inc. e(12),&
			&rotation(1),rotation(2),rotation(3),inc. rot.(1),inc. rot.(2),inc. rot.(3),&
			&phi1,Phi,phi2,inc. phi1,inc. Phi,inc. phi22,&
			&el_eps(11), el_eps(22), el_eps(33), el_eps(23), el_eps(31), el_eps(12),&
			&vM eps,vM'')')

		DO  ilmnt=1,nlmnt						!	Loop on elements, with any restriction
!			IF ( elrestrict(ilmnt) ) THEN

!-------Get element type,nodes and Gauss point numbers
				iel= le(ilmnt)
				itp= lt(ilmnt)
				nn = nd(ilmnt)			! number of nodes of ilmnt-th element
				ng = np(ilmnt)			! number of gauss points of ilmnt-th element
				in(1:nn)=ln( 1:nn, ilmnt)	! list indices of nodes of ilmnt-th element
				ip(1:ng)=lg( 1:ng, ilmnt)	! list indices of gauss points of ilmnt-th element

!--------------------------------Get coordinates 
				DO i= 1, nn
					ks= 3*in(i) - 3
					DO j= 1, 3
						xl(j,i) = x(ks+j)
					END DO
				END DO

				DO ig= 1, ng			!	Loop on gauss points of ilmnt-th element

!--------------------------------Choose weight function 
					SELECT CASE( itp )

						CASE(12)
							CALL WEIGHT_4( ig, xl, xg, weight )

						CASE(21)
							CALL WEIGHT_1( ig, xl, xg, weight )

						CASE(22)
							CALL WEIGHT_2( ig, xl, xg, weight )

						CASE(23)
							CALL WEIGHT_3( ig, xl, xg, weight )

					END SELECT

					stress(1:6)=st(1:6,ip(ig))
					CALL ELASEPS(iel, sv(1,ip(ig)),sv(2,ip(ig)),sv(3,ip(ig)), stress, strain)

					vm_strs=VONM_STRESS(st(1,ip(ig)),st(2,ip(ig)),st(3,ip(ig)),&
										&	st(4,ip(ig)),st(5,i),st(6,ip(ig)))


!--------------------------------------Write weight and values
					WRITE(iout,'(I5, 48('','',G13.5))')							&
					&		ip(ig), weight,	xg(1:3),							&
					&		(st(i,ip(ig)),i=1,6),(ds(i,ip(ig)),i=1,6),			&
					&		(ep(i,ip(ig)),i=1,6),(de(i,ip(ig)),i=1,6),			&
					&		(rt(i,ip(ig)),i=1,3),(dr(i,ip(ig)),i=1,3),			&
					&		(sv(i,ip(ig)),i=1,3),(dsv(i,ip(ig)),i=1,3),			&
					&		(strain(i),i=1,6),eps(ip(ig)), vm_strs
				
				END DO
			!END IF
		END DO

	CLOSE( iout)

END SUBROUTINE






!*******************************************************************************!
!																				!
!   SUBROUTINE ELSTATE															!
! 		Outputs a list (.CSV) of weights and mean element values.				!
!																				!
!*******************************************************************************!
SUBROUTINE ELSTATE(l_dummy)

	USE FAS_COM

	IMPLICIT NONE
	LOGICAL(4)	l_dummy, l_file
	INTEGER(4)	ilmnt, iout, iel, itp, nn, ng, in, ip, i, j, ks, ig
	REAL(4)		xl, xg, weight, ph1, phi, ph2, dph1, dphi, dph2, norm
	REAL(4)		xcrd, ycrd, zcrd, accum, q_mean, dq_mean, xq
	REAL(4)		vm_eps, vm_strs, VONM_STRESS

	DIMENSION	in(1:20), ip(1:8), xl(1:3, 1:20), xg(1:3)
	DIMENSION	accum(1:31), q_mean(0:3), dq_mean(0:3), xq(0:3)

	REAL(4), DIMENSION(1:6):: stress, strain

	PARAMETER( iout=53)

!----open file here---------------!!!!!!!
        
        
        CALL OUTFILE(iout,l_file)

		WRITE( iout,'(''FAFNER3 element mean state values'')')
		WRITE( iout,'(''element, volume, x, y, z,&
			&stress(11),stress(22),stress(33),stress(23),stress(31),stress(12),&
			&inc. s(11),inc. s(22),inc. s(33),inc. s(23),inc. s(31),inc. s(12),&
			&strain(11),strain(22),strain(33),strain(23),strain(31),strain(12),&
			&inc. e(11),inc. e(22),inc. e(33),inc. e(23),inc. e(31),inc. e(12),&
			&rotation(1),rotation(2),rotation(3),inc. rot.(1),inc. rot.(2),inc. rot.(3),&
			&phi1,Phi,phi2,inc. phi1,inc. Phi,inc. phi22,&
			&el_eps(11), el_eps(22), el_eps(33), el_eps(23), el_eps(31), el_eps(12),&
			&vM eps,vM sig, phase'')')

		DO  ilmnt=1,nlmnt						!	Loop on elements, with any restriction
			!IF ( elrestrict(ilmnt) .EQ. .TRUE. ) THEN

!-------Get element type,nodes and Gauss point numbers
				iel= le(ilmnt)
				itp= lt(ilmnt)
				nn = nd(ilmnt)		! number of nodes of ilmnt-th element
				ng = np(ilmnt)		! number of gauss points of ilmnt-th element
				in(1:nn)=ln( 1:nn, ilmnt)	! list indices of nodes of ilmnt-th element
				ip(1:ng)=lg( 1:ng, ilmnt)	! list indices of gauss points of ilmnt-th element

!--------------------------------Get coordinates 
				DO i= 1, nn
					ks= 3*in(i) - 3
					DO j= 1, 3
						xl(j,i) = x(ks+j)
					END DO
				END DO

				xcrd=SUM(xl(1,1:nn)); ycrd=SUM(xl(2,1:nn)); zcrd=SUM(xl(3,1:nn)) 

!--------------------------------initialise accumulators
				accum= 0.
				q_mean=0.
				dq_mean=0.
				vm_eps=0.

				DO ig= 1, ng			!	Loop on gauss points

!--------------------------------Choose weight function 
					SELECT CASE( itp )

						CASE(12)
							CALL WEIGHT_4( ig, xl, xg, weight )

						CASE(21)
							CALL WEIGHT_1( ig, xl, xg, weight )

						CASE(22)
							CALL WEIGHT_2( ig, xl, xg, weight )

						CASE(23)
							CALL WEIGHT_3( ig, xl, xg, weight )

					END SELECT

!--------------------------------------sum weight and values
					accum( 1)= accum(1)+weight
					accum( 2: 7)= accum( 2: 7)+st(1:6,ip(ig))
					accum( 8:13)= accum( 8:13)+ds(1:6,ip(ig))
					accum(14:19)= accum(14:19)+ep(1:6,ip(ig))
					accum(20:25)= accum(20:25)+de(1:6,ip(ig))
					accum(26:28)= accum(26:28)+rt(1:3,ip(ig))
					accum(29:31)= accum(29:31)+dr(1:3,ip(ig))

					CALL QTNEUL(sv(1,ip(ig)),sv(2,ip(ig)),sv(3,ip(ig)),xq)
					q_mean= q_mean + xq

					CALL QTNEUL(dsv(1,ip(ig)),dsv(2,ip(ig)),dsv(3,ip(ig)),xq)
					dq_mean= dq_mean + xq
				
					vm_eps= vm_eps+ eps(ip(ig))

				END DO

!--------------------------------------Write values
				CALL QTNNRM(q_mean, norm)
				CALL EULQTN(q_mean,ph1,phi,ph2)
				
				CALL QTNNRM(dq_mean, norm)
				CALL EULQTN(dq_mean,dph1,dphi,dph2)

				stress(1:6)=accum(2:7)/REAL(ng, 4)
				CALL ELASEPS(iel, ph1, phi, ph2, stress, strain)
				
				                vm_strs=VONM_STRESS(stress(1),stress(2),stress(3),stress(4),stress(5),stress(6))


				WRITE(iout,'(I5, 48('','',G13.5),'','',I5)') ilmnt, accum(1),				&
						&	xcrd/REAL(nn,4),ycrd/REAL(nn,4),zcrd/REAL(nn,4),		&
						&	accum(2:31)/REAL(ng,4),ph1,phi,ph2,				&
						&	dph1,dphi,dph2,strain(1:6),vm_eps/REAL(ng,4), vm_strs,iel
			!END IF
		END DO


	CLOSE(iout)

END SUBROUTINE

!***********************************************************************!
!																		!
!	SUBROUTINE POWER													!
!  Reports the power f.dx/dt for the last increment of the whole domain	!
!																		!
!***********************************************************************!
SUBROUTINE POWER

	USE FAS_COM

	IMPLICIT NONE
	CHARACTER	text*13
	INTEGER(4)  message
	
	WRITE( text,'(G12.5)') DOT_PRODUCT( fc(1:3*nnod),dx(1:3*nnod))/tiem(ninc1)


END SUBROUTINE

!***********************************************************************!
!																		!
!   SUBROUTINE ELASEPS													!
!-----------------------------------------------------------------------!
!        Returns the elastic strains given the stress state             !
!																		!
!***********************************************************************!
SUBROUTINE ELASEPS(iel, ph1, phi, ph2, stress, strain)

	USE FAS_COM

	IMPLICIT NONE
	INTEGER(4)  iel, iset, jset, i, j, k, l
	REAL(4)		c, s, c1, c2,s1, s2, c1c2, c1s2, s1c2, s1s2
	REAL(4)		ph1, phi, ph2, mult, term1, term2, term3

	REAL(4), DIMENSION(1:6)		:: stress, strain
	REAL(4), DIMENSION(1:3,1:3) :: a
	REAL(4), DIMENSION(1:6,1:6) :: el_comp
  


!-----------------------------Form trig function terms
	c1= COS(ph1)
	s1= SIN(ph1)
	c2= COS(ph2)
	s2= SIN(ph2)
	c=  COS(phi)
	s=  SIN(phi)
	c1c2=c1*c2
	s1s2=s1*s2
	c1s2=c1*s2
	s1c2=s1*c2

!---------------------------Form the (transposed) transform matrix
	a(1,1)= c1c2 - s1s2*c
	a(1,2)=-c1s2 - s1c2*c
	a(1,3)= s1*s
	a(2,1)= s1c2 + c1s2*c
	a(2,2)=-s1s2 + c1c2*c
	a(2,3)=-c1*s
	a(3,1)= s2*s
	a(3,2)= c2*s
	a(3,3)= c

!-----------------------loop on reduced tensor terms, expand and calculate compliance

	DO iset=1,6
		SELECT CASE(iset)
			CASE(1,2,3)
				i=iset
				j=iset
				mult=1.	    ! multiplier unity if no shear involved 
			CASE(4,5,6)
				i=1+MOD(iset-3,3)
				j=1+MOD(i,3)
				mult=2.		! double if shear involved ( tensor to 'engineering')
		END SELECT

		DO jset=iset,6
			SELECT CASE(jset)
				CASE(1,2,3)
					k=jset
					l=jset
				CASE(4,5,6)
					k=1+MOD(jset-3,3)
					l=1+MOD(k,3)
			END SELECT
			
		term1=   a(i,1)*a(j,1)*a(k,1)*a(l,1)*e(1,iel)	&
			&	+a(i,2)*a(j,2)*a(k,2)*a(l,2)*e(2,iel)	&
			&	+a(i,3)*a(j,3)*a(k,3)*a(l,3)*e(3,iel)

		term2=	 a(i,2)*a(j,2)*a(k,3)*a(l,3)*e(4,iel)	&
			&	+a(i,3)*a(j,3)*a(k,2)*a(l,2)*e(4,iel)	&
			&	+a(i,3)*a(j,3)*a(k,1)*a(l,1)*e(5,iel)	&
			&	+a(i,1)*a(j,1)*a(k,3)*a(l,3)*e(5,iel)	&
			&	+a(i,1)*a(j,1)*a(k,2)*a(l,2)*e(6,iel)	&
			&	+a(i,2)*a(j,2)*a(k,1)*a(l,1)*e(6,iel)

		term3=	 a(i,2)*a(j,3)*a(k,2)*a(l,3)*e(7,iel)	&
			&	+a(i,3)*a(j,2)*a(k,3)*a(l,2)*e(7,iel)	&
			&	+a(i,3)*a(j,1)*a(k,3)*a(l,1)*e(8,iel)	&
			&	+a(i,1)*a(j,3)*a(k,1)*a(l,3)*e(8,iel)	&
			&	+a(i,1)*a(j,2)*a(k,1)*a(l,2)*e(9,iel)	&
			&	+a(i,2)*a(j,1)*a(k,2)*a(l,1)*e(9,iel)	

	el_comp(iset,jset)= mult*(term1 + term2 + 0.5*term3)
		END DO
	END DO

!-------------------------------------------------complete symmetric matrix
	DO iset=1,5
		DO jset= iset+1,6
			el_comp(jset,iset)= el_comp(iset,jset)
		END DO
	END DO

!----------------------------------------------------get the elastic strain

	strain= MATMUL(el_comp, stress)




END SUBROUTINE	

!***********************************************************************!
!																		!
!                       FASPLOT3   -WEIGHTS								!
!                        Pete Bate 2000									!
!-----------------------------------------------------------------------!
!            Gauss point volume routines for 3-d FEM					!
!-----------------------------------------------------------------------!
!																		!
!  In all routines, arguments are:										!
!        gpnumber    INTEGER(4)         The local Gauss point number	!
!        xnode       REAL(4) array(3,n) The global node coordinates		!
!        weight		 REAL				The volume weighting			!
!-----------------------------------------------------------------------!
!    NB  function routines have FAILED if weight zero					!
!																		!
!***********************************************************************!
!
!
!***********************************************************************!
!																		!
!          WEIGHT_1    10 node, 4 Gauss point tetrahedron				!
!-----------------------------------------------------------------------!
!   Corner nodes 1,3,5,10 in RH sense:									!
!				looking down from 10; 1,3,5 anticlockwise				!
!	Intermediate nodes are between:										!
!             2 (1,3) 4(3,5), 6(1,5), 7(1,10), 8(3,10), 9(5,10)			!
!   Gauss points are at volume coordinates of type						!
!                             {0.585, 0.138, 0.138, 0.138}				!
!    -with 1 nearest node 1, 2 nearest node 3 etc.						!
!-----------------------------------------------------------------------!
!  Tested uniform displacement field etc. 14/2/00						!
!																		!
!***********************************************************************!
SUBROUTINE WEIGHT_1( gpnumber, xnode, xg, weight )

	IMPLICIT NONE

!--------------------------------Argument type declarations
	INTEGER(4)	gpnumber
	REAL(4)     xnode, xg, weight

	DIMENSION	xnode( 1:3, 1:10), xg( 1:3)
!------------------------------------Local type declarations
	INTEGER(4)	i, j, k
	REAL(4)		shapef, dndn, jac_mat, sum
	DIMENSION	shapef( 1:10, 1:4), dndn( 1:3, 1:10, 1:4), jac_mat( 1:3, 1:3)

!-------------Data for shape functions and local basis shape function derivatives-----!

	DATA (shapef(i, 1), i=1,10)/   0.100000, 0.323607,-0.100000, 0.076393,-0.100000, &
								&  0.323607, 0.323607, 0.076393, 0.076393,-0.100000	 /
	DATA (shapef(i, 2), i=1,10)/  -0.100000, 0.323607, 0.100000, 0.323607,-0.100000, &
								&  0.076393, 0.076393, 0.323607, 0.076393,-0.100000  /
	DATA (shapef(i, 3), i=1,10)/  -0.100000, 0.076393,-0.100000, 0.323607, 0.100000, &
								&  0.323607, 0.076393, 0.076393, 0.323607,-0.100000  /
	DATA (shapef(i, 4), i=1,10)/  -0.100000, 0.076393,-0.100000, 0.076393,-0.100000, &
								&  0.076393, 0.323607, 0.323607, 0.323607, 0.100000  /


	DATA (dndn(1, i, 1),i=1,10)/   0.000000, 0.000000, 0.000000,-0.552786, 0.447214, &
								& -2.341641, 2.341641, 0.552786, 0.000000,-0.447214  /
	DATA (dndn(1, i, 2),i=1,10)/   0.000000, 0.000000, 0.000000,-2.341641, 0.447214, &
								& -0.552786, 0.552786, 2.341641, 0.000000,-0.447214  /
	DATA (dndn(1, i, 3),i=1,10)/   0.000000, 0.000000, 0.000000,-0.552786,-1.341641, &
								& -0.552786, 0.552786, 0.552786, 1.788854,-0.447214  /
	DATA (dndn(1, i, 4),i=1,10)/   0.000000, 0.000000, 0.000000,-0.552786, 0.447214, &
								& -0.552786, 0.552786, 0.552786,-1.788854, 1.341641  /	

	DATA (dndn(2, i, 1),i=1,10)/   0.000000, 2.341641,-0.447214, 0.000000, 0.447214, &
								& -2.341641, 0.000000, 0.552786,-0.552786, 0.000000  /
	DATA (dndn(2, i, 2),i=1,10)/   0.000000, 0.552786, 1.341641,-1.788854, 0.447214, &
								& -0.552786, 0.000000, 0.552786,-0.552786, 0.000000  /
	DATA (dndn(2, i, 3),i=1,10)/   0.000000, 0.552786,-0.447214, 1.788854,-1.341641, &
								& -0.552786, 0.000000, 0.552786,-0.552786, 0.000000  /
	DATA (dndn(2, i, 4),i=1,10)/   0.000000, 0.552786,-0.447214, 0.000000, 0.447214, &
								& -0.552786, 0.000000, 2.341641,-2.341641, 0.000000  /

	DATA (dndn(3, i, 1),i=1,10)/   1.341641, 0.552786, 0.000000,-0.552786, 0.447214, &
								& -1.788854, 0.552786, 0.000000,-0.552786, 0.000000  /
	DATA (dndn(3, i, 2),i=1,10)/  -0.447214, 2.341641, 0.000000,-2.341641, 0.447214, &
								&  0.000000, 0.552786, 0.000000,-0.552786, 0.000000  /
	DATA (dndn(3, i, 3),i=1,10)/  -0.447214, 0.552786, 0.000000,-0.552786,-1.341641, &
								&  1.788854, 0.552786, 0.000000,-0.552786, 0.000000  /
	DATA (dndn(3, i, 4),i=1,10)/  -0.447214, 0.552786, 0.000000,-0.552786, 0.447214, &
								&  0.000000, 2.341641, 0.000000,-2.341641, 0.000000  /

!---------------------------------------------Gauss point coordinates
	xg(1)= DOT_PRODUCT( shapef(1:10,gpnumber), xnode(1,1:10) )
	xg(2)= DOT_PRODUCT( shapef(1:10,gpnumber), xnode(2,1:10) )
	xg(3)= DOT_PRODUCT( shapef(1:10,gpnumber), xnode(3,1:10) )

!---------------------------------------------Sum dndn.x for Jacobian
	DO k= 1, 3				!Global coordinate
		DO j= 1, 3			!Local coordinate
			sum= 0.
			DO i=1,10		!Node number
				sum= sum + dndn( j, i, gpnumber )* xnode(k, i )
			END DO
			jac_mat( j, k)= sum 
		END DO
	END DO

!-----------------------------------Get determinant
	CALL DETERM_3( jac_mat, weight)

	weight= - 0.0416667 * weight		!Gauss weighting for all

END SUBROUTINE 

!***********************************************************************!
!																		!
!          WEIGHT_2    20 node, 8 Gauss point hexahedron				!
!-----------------------------------------------------------------------!
!   Bottom face (Z = -1),												!
!				anticlockwise view from top: corner nodes 1,3,5,7		!
!										 mid-sides    2,4,6,8			!
!   Top face (Z = 1),													!
!				anticlockwise view from top: corner nodes 13,15,17,19	!
!										 mid-sides    14,16,18,20   	!
!	Intermediate nodes are between:										!
!					9 (1,13) 10(3,15), 11(5,17), 12(7,19)				!
!   Gauss points are at volume coordinates of type						!
!                             {+-0.5773 x,y,z}							!
!   -with 1..4 anticlock. above bottom face, 6..8 below top face		!
!-----------------------------------------------------------------------!
!  Tested uniform displacement field etc. 18/2/00						!
!																		!
!***********************************************************************!
SUBROUTINE WEIGHT_2( gpnumber, xnode, xg, weight )

	IMPLICIT NONE
!--------------------------------Argument type declarations
	INTEGER(4)	gpnumber
	REAL(4)     xnode, xg, weight

	DIMENSION	xnode( 1:3, 1:20), xg( 1:3)

!------------------------------------Local type declarations
	INTEGER(4)	i, j, k
	REAL(4)		shapef, dndn, jac_mat, sum
	DIMENSION	shapef( 1:20, 1:8), dndn( 1:3, 1:20, 1:8), jac_mat( 1:3, 1:3)

!----------Data for shape functions and local basis shape function derivatives--------!

	DATA (shapef(i, 1), i=1,20)/  -0.131446, 0.414672,-0.187001, 0.111111,-0.090776, &
								&  0.111111,-0.187001, 0.414672, 0.414672, 0.111111, &
								&  0.029772, 0.111111,-0.187001, 0.111111,-0.090776, &
								&  0.029772,-0.035221, 0.029772,-0.090776, 0.111111  /
	DATA (shapef(i, 2), i=1,20)/  -0.187001, 0.414672,-0.131446, 0.414672,-0.187001, &
								&  0.111111,-0.090776, 0.111111, 0.111111, 0.414672, &
								&  0.111111, 0.029772,-0.090776, 0.111111,-0.187001, &
								&  0.111111,-0.090776, 0.029772,-0.035221, 0.029772  /
	DATA (shapef(i, 3), i=1,20)/  -0.090776, 0.111111,-0.187001, 0.414672,-0.131446, &
								&  0.414672,-0.187001, 0.111111, 0.029772, 0.111111, &
								&  0.414672, 0.111111,-0.035221, 0.029772,-0.090776, &
								&  0.111111,-0.187001, 0.111111,-0.090776, 0.029772  /
	DATA (shapef(i, 4), i=1,20)/  -0.187001, 0.111111,-0.090776, 0.111111,-0.187001, &
								&  0.414672,-0.131446, 0.414672, 0.111111, 0.029772, &
								&  0.111111, 0.414672,-0.090776, 0.029772,-0.035221, &
								&  0.029772,-0.090776, 0.111111,-0.187001, 0.111111  /
	DATA (shapef(i, 5), i=1,20)/  -0.187001, 0.111111,-0.090776, 0.029772,-0.035221, &
								&  0.029772,-0.090776, 0.111111, 0.414672, 0.111111, &
								&  0.029772, 0.111111,-0.131446, 0.414672,-0.187001, &
								&  0.111111,-0.090776, 0.111111,-0.187001, 0.414672  /
	DATA (shapef(i, 6), i=1,20)/  -0.090776, 0.111111,-0.187001, 0.111111,-0.090776, &
								&  0.029772,-0.035221, 0.029772, 0.111111, 0.414672, &
								&  0.111111, 0.029772,-0.187001, 0.414672,-0.131446, &
								&  0.414672,-0.187001, 0.111111,-0.090776, 0.111111  /
	DATA (shapef(i, 7), i=1,20)/  -0.035221, 0.029772,-0.090776, 0.111111,-0.187001, &
								&  0.111111,-0.090776, 0.029772, 0.029772, 0.111111, &
								&  0.414672, 0.111111,-0.090776, 0.111111,-0.187001, &
								&  0.414672,-0.131446, 0.414672,-0.187001, 0.111111  /
	DATA (shapef(i, 8), i=1,20)/  -0.090776, 0.029772,-0.035221, 0.029772,-0.090776, &
								&  0.111111,-0.187001, 0.111111, 0.111111, 0.029772, &
								&  0.111111, 0.414672,-0.187001, 0.111111,-0.090776, &
								&  0.111111,-0.187001, 0.414672,-0.131446, 0.414672  /

	DATA (((dndn(j, i, k),i=1,20), j=1, 3), k= 1, 8)/ &
&-0.407229, 0.718233,-0.311004, 0.262892,-0.179558, 0.192450,-0.012892,-0.262892,-0.262892, 0.262892,&
& 0.070442,-0.070442,-0.012892, 0.192450,-0.179558, 0.070442,-0.073896, 0.051567, 0.022329,-0.070442,&
&-0.407229,-0.262892,-0.012892, 0.192450,-0.179558, 0.262892,-0.311004, 0.718233,-0.262892,-0.070442,&
& 0.070442, 0.262892,-0.012892,-0.070442, 0.022329, 0.051567,-0.073896, 0.070442,-0.179558, 0.192450,&
&-0.407229,-0.262892,-0.012892,-0.070442, 0.022329,-0.070442,-0.012892,-0.262892, 0.718233, 0.192450,&
& 0.051567, 0.192450,-0.311004, 0.262892,-0.179558, 0.070442,-0.073896, 0.070442,-0.179558, 0.262892,&
& 0.311004,-0.718233, 0.407229, 0.262892, 0.012892,-0.192450, 0.179558,-0.262892,-0.262892, 0.262892,&
& 0.070442,-0.070442, 0.179558,-0.192450, 0.012892, 0.070442,-0.022329,-0.051567, 0.073896,-0.070442,&
&-0.012892,-0.262892,-0.407229, 0.718233,-0.311004, 0.262892,-0.179558, 0.192450,-0.070442,-0.262892,&
& 0.262892, 0.070442, 0.022329,-0.070442,-0.012892, 0.192450,-0.179558, 0.070442,-0.073896, 0.051567,&
&-0.012892,-0.262892,-0.407229,-0.262892,-0.012892,-0.070442, 0.022329,-0.070442, 0.192450, 0.718233,&
& 0.192450, 0.051567,-0.179558, 0.262892,-0.311004, 0.262892,-0.179558, 0.070442,-0.073896, 0.070442,&
& 0.179558,-0.192450, 0.012892, 0.262892, 0.407229,-0.718233, 0.311004,-0.262892,-0.070442, 0.070442,&
& 0.262892,-0.262892, 0.073896,-0.051567,-0.022329, 0.070442, 0.012892,-0.192450, 0.179558,-0.070442,&
& 0.179558,-0.262892, 0.311004,-0.718233, 0.407229, 0.262892, 0.012892,-0.192450,-0.070442,-0.262892,&
& 0.262892, 0.070442, 0.073896,-0.070442, 0.179558,-0.192450, 0.012892, 0.070442,-0.022329,-0.051567,&
& 0.022329,-0.070442,-0.012892,-0.262892,-0.407229,-0.262892,-0.012892,-0.070442, 0.051567, 0.192450,&
& 0.718233, 0.192450,-0.073896, 0.070442,-0.179558, 0.262892,-0.311004, 0.262892,-0.179558, 0.070442,&
&-0.012892, 0.192450,-0.179558, 0.262892,-0.311004, 0.718233,-0.407229,-0.262892,-0.070442, 0.070442,&
& 0.262892,-0.262892, 0.022329, 0.051567,-0.073896, 0.070442,-0.179558, 0.192450,-0.012892,-0.070442,&
& 0.311004,-0.262892, 0.179558,-0.192450, 0.012892, 0.262892, 0.407229,-0.718233,-0.262892,-0.070442,&
& 0.070442, 0.262892, 0.179558,-0.070442, 0.073896,-0.051567,-0.022329, 0.070442, 0.012892,-0.192450,&
&-0.012892,-0.070442, 0.022329,-0.070442,-0.012892,-0.262892,-0.407229,-0.262892, 0.192450, 0.051567,&
& 0.192450, 0.718233,-0.179558, 0.070442,-0.073896, 0.070442,-0.179558, 0.262892,-0.311004, 0.262892,&
&-0.012892, 0.192450,-0.179558, 0.070442,-0.073896, 0.051567, 0.022329,-0.070442,-0.262892, 0.262892,&
& 0.070442,-0.070442,-0.407229, 0.718233,-0.311004, 0.262892,-0.179558, 0.192450,-0.012892,-0.262892,&
&-0.012892,-0.070442, 0.022329, 0.051567,-0.073896, 0.070442,-0.179558, 0.192450,-0.262892,-0.070442,&
& 0.070442, 0.262892,-0.407229,-0.262892,-0.012892, 0.192450,-0.179558, 0.262892,-0.311004, 0.718233,&
& 0.311004,-0.262892, 0.179558,-0.070442, 0.073896,-0.070442, 0.179558,-0.262892,-0.718233,-0.192450,&
&-0.051567,-0.192450, 0.407229, 0.262892, 0.012892, 0.070442,-0.022329, 0.070442, 0.012892, 0.262892,&
& 0.179558,-0.192450, 0.012892, 0.070442,-0.022329,-0.051567, 0.073896,-0.070442,-0.262892, 0.262892,&
& 0.070442,-0.070442, 0.311004,-0.718233, 0.407229, 0.262892, 0.012892,-0.192450, 0.179558,-0.262892,&
& 0.022329,-0.070442,-0.012892, 0.192450,-0.179558, 0.070442,-0.073896, 0.051567,-0.070442,-0.262892,&
& 0.262892, 0.070442,-0.012892,-0.262892,-0.407229, 0.718233,-0.311004, 0.262892,-0.179558, 0.192450,&
& 0.179558,-0.262892, 0.311004,-0.262892, 0.179558,-0.070442, 0.073896,-0.070442,-0.192450,-0.718233,&
&-0.192450,-0.051567, 0.012892, 0.262892, 0.407229, 0.262892, 0.012892, 0.070442,-0.022329, 0.070442,&
& 0.073896,-0.051567,-0.022329, 0.070442, 0.012892,-0.192450, 0.179558,-0.070442,-0.070442, 0.070442,&
& 0.262892,-0.262892, 0.179558,-0.192450, 0.012892, 0.262892, 0.407229,-0.718233, 0.311004,-0.262892,&
& 0.073896,-0.070442, 0.179558,-0.192450, 0.012892, 0.070442,-0.022329,-0.051567,-0.070442,-0.262892,&
& 0.262892, 0.070442, 0.179558,-0.262892, 0.311004,-0.718233, 0.407229, 0.262892, 0.012892,-0.192450,&
& 0.073896,-0.070442, 0.179558,-0.262892, 0.311004,-0.262892, 0.179558,-0.070442,-0.051567,-0.192450,&
&-0.718233,-0.192450,-0.022329, 0.070442, 0.012892, 0.262892, 0.407229, 0.262892, 0.012892, 0.070442,&
& 0.022329, 0.051567,-0.073896, 0.070442,-0.179558, 0.192450,-0.012892,-0.070442,-0.070442, 0.070442,&
& 0.262892,-0.262892,-0.012892, 0.192450,-0.179558, 0.262892,-0.311004, 0.718233,-0.407229,-0.262892,&
& 0.179558,-0.070442, 0.073896,-0.051567,-0.022329, 0.070442, 0.012892,-0.192450,-0.262892,-0.070442,&
& 0.070442, 0.262892, 0.311004,-0.262892, 0.179558,-0.192450, 0.012892, 0.262892, 0.407229,-0.718233,&
& 0.179558,-0.070442, 0.073896,-0.070442, 0.179558,-0.262892, 0.311004,-0.262892,-0.192450,-0.051567,&
&-0.192450,-0.718233, 0.012892, 0.070442,-0.022329, 0.070442, 0.012892, 0.262892, 0.407229, 0.262892 /

!---------------------------------------------Gauss point coordinates
	xg(1)= DOT_PRODUCT( shapef(1:20,gpnumber), xnode(1,1:20) )
	xg(2)= DOT_PRODUCT( shapef(1:20,gpnumber), xnode(2,1:20) )
	xg(3)= DOT_PRODUCT( shapef(1:20,gpnumber), xnode(3,1:20) )

!---------------------------------------------Sum dndn.x for Jacobian
	DO k= 1, 3				!Global coordinate
		DO j= 1, 3			!Local coordinate
			sum= 0.
			DO i=1,20		!Node number
				sum= sum + dndn( j, i, gpnumber )* xnode(k, i )
			END DO
			jac_mat( j, k)= sum 
		END DO
	END DO

!-----------------------------------Get determinant
	CALL DETERM_3( jac_mat, weight)

	weight=  1.0 *weight		!Gauss weighting for all

END SUBROUTINE 

!***********************************************************************!
!																		!
!          WEIGHT_3    15 node, 6 Gauss point pentahedron				!
!-----------------------------------------------------------------------!
!   Corner nodes 1,3,5 base, 10,12,14 top								!
!				 anticlockwise	height +ve definition					!
!	intermediate nodes are (between):									!
!				2(1,3) 4(3,5), 6(5,1), 7(1,10), 8(3,12),				!
!				9(5,14), 11(10,12),13(12,14),15(14,10)					!
!   Gauss points are at area+z coordinates of type						!
!                          {0.66, 0.16, 0.16, 0.577}					!
!    -with 1 nearest node 1, 2 nearest node 3 etc.						!
!-----------------------------------------------------------------------!
!  Tested uniform displacement field etc. 21/2/00						!
!																		!
!***********************************************************************!
SUBROUTINE WEIGHT_3( gpnumber, xnode, xg, weight )

	IMPLICIT NONE

!--------------------------------Argument type declarations
	INTEGER(4)	gpnumber
	REAL(4)     xnode, xg, weight

	DIMENSION	xnode( 1:3, 1:15), xg( 1:3)

!------------------------------------Local type declarations
	INTEGER(4)	i, j, k
	REAL(4)		shapef, dndn, jac_mat, sum
	DIMENSION	shapef( 1:15, 1:6), dndn( 1:3, 1:15, 1:6), jac_mat( 1:3, 1:3)

!---------Data for shape functions and local basis shape function derivatives-----!
	DATA (shapef(i, 1), i=1,15)/  -0.046961, 0.350522,-0.143186, 0.087631,-0.143186, &
								&  0.350522, 0.444444, 0.111111, 0.111111,-0.175261, &
								&  0.093922,-0.079036, 0.023481,-0.079036, 0.093922  /
	DATA (shapef(i, 2), i=1,15)/  -0.143186, 0.350522,-0.046961, 0.350522,-0.143186, &
								&  0.087631, 0.111111, 0.444444, 0.111111,-0.079036, &
								&  0.093922,-0.175261, 0.093922,-0.079036, 0.023481  /
	DATA (shapef(i, 3), i=1,15)/  -0.143186, 0.087631,-0.143186, 0.350522,-0.046961, &
								&  0.350522, 0.111111, 0.111111, 0.444444,-0.079036, &
								&  0.023481,-0.079036, 0.093922,-0.175261, 0.093922  /
	DATA (shapef(i, 4), i=1,15)/  -0.175261, 0.093922,-0.079036, 0.023481,-0.079036, &
								&  0.093922, 0.444444, 0.111111, 0.111111,-0.046961, &
								&  0.350522,-0.143186, 0.087631,-0.143186, 0.350522  /
	DATA (shapef(i, 5), i=1,15)/  -0.079036, 0.093922,-0.175261, 0.093922,-0.079036, &
								&  0.023481, 0.111111, 0.444444, 0.111111,-0.143186, &
								&  0.350522,-0.046961, 0.350522,-0.143186, 0.087631  /
	DATA (shapef(i, 6), i=1,15)/  -0.079036, 0.023481,-0.079036, 0.093922,-0.175261, &
								&  0.093922, 0.111111, 0.111111, 0.444444,-0.143186, &
								&  0.087631,-0.143186, 0.350522,-0.046961, 0.350522  /


	DATA (dndn(1, i, 1),i=1,15)/   0.981124, 0.525783, 0.000000,-0.525784, 0.596225, &
								& -1.577352, 0.666666, 0.000000,-0.666667, 0.018875, &
								&  0.140883, 0.000000,-0.140883, 0.403775,-0.422650  /
	DATA (dndn(1, i, 2),i=1,15)/  -0.596225, 2.103135, 0.000000,-2.103135, 0.596225, &
								&  0.000000, 0.666667, 0.000000,-0.666667,-0.403775, &
								&  0.563533, 0.000000,-0.563533, 0.403775, 0.000000  /
	DATA (dndn(1, i, 3),i=1,15)/  -0.596225, 0.525784, 0.000000,-0.525783,-0.981124, &
								&  1.577352, 0.666667, 0.000000,-0.666666,-0.403775, &
								&  0.140883, 0.000000,-0.140883,-0.018875, 0.422650 /
	DATA (dndn(1, i, 4),i=1,15)/   0.018875, 0.140883, 0.000000,-0.140883, 0.403775, &
								& -0.422650, 0.666666, 0.000000,-0.666667, 0.981124, &
								&  0.525783, 0.000000,-0.525784, 0.596225,-1.577352  /
	DATA (dndn(1, i, 5),i=1,15)/  -0.403775, 0.563533, 0.000000,-0.563533, 0.403775, &
								&  0.000000, 0.666667, 0.000000,-0.666667,-0.596225, &
								&  2.103135, 0.000000,-2.103135, 0.596225, 0.000000  /
	DATA (dndn(1, i, 6),i=1,15)/  -0.403775, 0.140883, 0.000000,-0.140883,-0.018875, &
								&  0.422650, 0.666667, 0.000000,-0.666666,-0.596225, &
								&  0.525784, 0.000000,-0.525783,-0.981124, 1.577352  /

	DATA (dndn(2, i, 1),i=1,15)/   0.000000, 2.103135,-0.596225, 0.000000, 0.596225, &
								& -2.103135, 0.000000, 0.666667,-0.666667, 0.000000, &
								&  0.563533,-0.403775, 0.000000, 0.403775,-0.563533  /
	DATA (dndn(2, i, 2),i=1,15)/   0.000000, 0.525783, 0.981124,-1.577352, 0.596225, &
								& -0.525784, 0.000000, 0.666666,-0.666667, 0.000000, &
								&  0.140883, 0.018875,-0.422650, 0.403775,-0.140883  /
	DATA (dndn(2, i, 3),i=1,15)/   0.000000, 0.525784,-0.596225, 1.577352,-0.981124, &
								& -0.525783, 0.000000, 0.666667,-0.666666, 0.000000, &
								&  0.140883,-0.403775, 0.422650,-0.018875,-0.140883  /
	DATA (dndn(2, i, 4),i=1,15)/   0.000000, 0.563533,-0.403775, 0.000000, 0.403775, &
								& -0.563533, 0.000000, 0.666667,-0.666667, 0.000000, &
								&  2.103135,-0.596225, 0.000000, 0.596225,-2.103135  /
	DATA (dndn(2, i, 5),i=1,15)/   0.000000, 0.140883, 0.018875,-0.422650, 0.403775, &
								& -0.140883, 0.000000, 0.666666,-0.666667, 0.000000, &
								&  0.525783, 0.981124,-1.577352, 0.596225,-0.525784  /
	DATA (dndn(2, i, 6),i=1,15)/   0.000000, 0.140883,-0.403775, 0.422650,-0.018875, &
								& -0.140883, 0.000000, 0.666667,-0.666666, 0.000000, &
								&  0.525784,-0.596225, 1.577352,-0.981124,-0.525783  /

	DATA (dndn(3, i, 1),i=1,15)/  -0.496011,-0.222222,-0.040669,-0.055555,-0.040669, &
								& -0.222222, 0.769800, 0.192450, 0.192450,-0.273789, &
								&  0.222222,-0.151780, 0.055555,-0.151780, 0.222222  /
	DATA (dndn(3, i, 2),i=1,15)/  -0.040669,-0.222222,-0.496011,-0.222222,-0.040669, &
								& -0.055555, 0.192450, 0.769800, 0.192450,-0.151780, &
								&  0.222222,-0.273789, 0.222222,-0.151780, 0.055555  /
	DATA (dndn(3, i, 3),i=1,15)/  -0.040669,-0.055556,-0.040669,-0.222222,-0.496011, &
								& -0.222222, 0.192450, 0.192450, 0.769800,-0.151780, &
								&  0.055556,-0.151780, 0.222222,-0.273789, 0.222222  /
	DATA (dndn(3, i, 4),i=1,15)/   0.273789,-0.222222, 0.151780,-0.055555, 0.151780, &
								& -0.222222,-0.769800,-0.192450,-0.192450, 0.496011, &
								&  0.222222, 0.040669, 0.055555, 0.040669, 0.222222  /
	DATA (dndn(3, i, 5),i=1,15)/   0.151780,-0.222222, 0.273789,-0.222222, 0.151780, &
								& -0.055555,-0.192450,-0.769800,-0.192450, 0.040669, &
								&  0.222222, 0.496011, 0.222222, 0.040669, 0.055555  /
	DATA (dndn(3, i, 6),i=1,15)/   0.151780,-0.055556, 0.151780,-0.222222, 0.273789, &
								& -0.222222,-0.192450,-0.192450,-0.769800, 0.040669, &
								&  0.055556, 0.040669, 0.222222, 0.496011, 0.222222  /

!---------------------------------------------Gauss point coordinates
	xg(1)= DOT_PRODUCT( shapef(1:15,gpnumber), xnode(1,1:15) )
	xg(2)= DOT_PRODUCT( shapef(1:15,gpnumber), xnode(2,1:15) )
	xg(3)= DOT_PRODUCT( shapef(1:15,gpnumber), xnode(3,1:15) )

!---------------------------------------------Sum dndn.x for Jacobian
	DO k= 1, 3				!Global coordinate
		DO j= 1, 3			!Local coordinate
			sum= 0.
			DO i=1,15		!Node number
				sum= sum + dndn( j, i, gpnumber )* xnode(k, i )
			END DO
			jac_mat( j, k)= sum 
		END DO
	END DO

!-----------------------------------Invert Jacobian, getting determinant
	CALL DETERM_3( jac_mat, weight)

	weight= 0.16666667*weight		!Gauss weighting	

END SUBROUTINE 
!***********************************************************************!
!																		!
!          WEIGHT_4    8 node, 1 Gauss point hexahedron					!
!-----------------------------------------------------------------------!
!   Nodes 1,2,3,4 bottom face, 5,6,7,8 top face (anticlockwise,			!
!									             looking down from top)	!
!   Gauss point is at coordinates  {0, 0, 0}							!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE WEIGHT_4( gpnumber, xnode, xg, weight )

	IMPLICIT NONE
!--------------------------------Argument type declarations
	INTEGER(4)	gpnumber
	REAL(4)     xnode, xg, weight

	DIMENSION	xnode( 1:3, 1:8), xg( 1:3)

!------------------------------------Local type declarations
	INTEGER(4)	i, j, k
	REAL(4)		shapef, dndn, jac_mat, sum
	DIMENSION	shapef( 1:8, 1:1), dndn( 1:3, 1:8, 1:1), jac_mat( 1:3, 1:3)

!---------Data for shape functions and local basis shape function derivatives-----!

	DATA (shapef(i, 1), i=1,8)/  8* 0.125 	 /

	DATA (dndn(1, i, 1),i=1,8)/  -0.125, 0.125, 0.125,-0.125,-0.125, 0.125, 0.125,-0.125 /
	DATA (dndn(2, i, 1),i=1,8)/  -0.125,-0.125, 0.125, 0.125,-0.125,-0.125, 0.125, 0.125 /
	DATA (dndn(3, i, 1),i=1,8)/  -0.125,-0.125,-0.125,-0.125, 0.125, 0.125, 0.125, 0.125 /

!---------------------------------------------Gauss point coordinates
	xg(1)= DOT_PRODUCT( shapef(1:8,gpnumber), xnode(1,1:8) )
	xg(2)= DOT_PRODUCT( shapef(1:8,gpnumber), xnode(2,1:8) )
	xg(3)= DOT_PRODUCT( shapef(1:8,gpnumber), xnode(3,1:8) )

!---------------------------------------------Sum dndn.x for Jacobian
	DO k= 1, 3				!Global coordinate
		DO j= 1, 3			!Local coordinate
			sum= 0.
			DO i=1,8		!Node number
				sum= sum + dndn( j, i, gpnumber )* xnode(k, i )
			END DO
			jac_mat( j, k)= sum 
		END DO
	END DO

!-----------------------------------Invert Jacobian, getting determinant
	CALL DETERM_3( jac_mat, weight)

	weight= 8.*weight		!Gauss weighting for all

END SUBROUTINE WEIGHT_4

!***********************************************************************!
!																		!
!     SUBROUTINE DETERM_3												!
!						 Pete Bate										!
!-----------------------------------------------------------------------!
!		   Routine for determinant of a 3x3 square matrix				!
!-----------------------------------------------------------------------!
!  Arguments:   matrix		REAL ARRAY(1:3, 1:3) the matrix				!
!				determinant REAL, matrix determinant					!
!																		!
!***********************************************************************!
SUBROUTINE DETERM_3( matrix, determinant)

	IMPLICIT NONE

	REAL(4)		matrix, determinant
	DIMENSION   matrix( 1:3, 1:3 )

	INTEGER(4)	i, j, JM, ia, ib, ja, jb
	REAL(4)		local
	DIMENSION	local( 1:3, 1:3)

	JM(i) = 1 + MOD( i, 3)  !statement function for cyclic next in 1,2,3,1,..

	local = matrix			!put input matrix to local array

	DO	i= 1, 3
		ia= JM(i)
		ib= JM(ia)

		DO j= 1, 3
			ja= JM( j)
			jb= JM(ja)

			matrix(j,i)= local(ia,ja)*local(ib,jb) - local(ia,jb)*local(ib,ja)
		END DO
	END DO
	
	determinant= matrix(1,1)*local(1,1) +matrix(1,2)*local(2,1) + matrix(1,3)*local(3,1)

END SUBROUTINE

!***********************************************************************************!
!																					!
!	REAL(4) FUNCTION VONM_STRESS													!
!-----------------------------------------------------------------------------------!
!			Von Mises effective stress.	 						 					!
!																					!
!***********************************************************************************!
REAL(4) FUNCTION VONM_STRESS(s1,s2,s3,s4,s5,s6)

	IMPLICIT NONE
	REAL(4)		s1, s2, s3, s4, s5, s6, d1, d2, d3, hydro

    hydro=( s1+s2+s3 )/3.		!hydrostatic

    d1= s1 - hydro
    d2= s2 - hydro
    d3= s3 - hydro

    VONM_STRESS= SQRT( 1.5*(d1*d1 + d2*d2 + d3*d3 + 2.*(s4*s4 + s5*s5 + s6*s6) ) ) 

END FUNCTION

!***********************************************************************!
!																		!
!		SUBROUTINE OUTFILE												!
!-----------------------------------------------------------------------!
!	Routine to write out .csv file									!
!																		!
!-----------------------------------------------------------------------!
!Arguments:IN	INTEGER;	READ channel number 						!
!L_FILE		LOGICAL;    TRUE if file opened ok							!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!



SUBROUTINE 	OUTFILE(iout,l_file)
IMPLICIT 	none
LOGICAL(4) :: 		l_file
INTEGER ::			iout, indir, ierr
CHARACTER(100) ::	f_name,dir,dir1,dir2,dir3, dir4
CHARACTER(4) ::		ext

indir=59
OPEN(indir,FILE='dirname',STATUS='OLD', IOSTAT=ierr)
	IF ( ierr .EQ. 0 ) THEN
		READ(indir,*) dir1,dir2,dir3,dir4
		dir=dir3
	ELSE
		WRITE(*,'("No dirname file found. Path is current directory.")')
		dir=''
	END IF
CLOSE(indir)

WRITE(*,'("Name of .csv file (no extension, please): ")')
READ(*,*) f_name
ext=".csv"
f_name=TRIM(dir)//TRIM(f_name)//TRIM(ext)

OPEN( UNIT= iout, FILE= f_name, FORM= 'FORMATTED',		&
	&	ACCESS= 'SEQUENTIAL', STATUS= 'UNKNOWN',IOSTAT= IERR)
	IF ( ierr .EQ. 0 ) THEN
		l_file = .TRUE.
		ELSE
		WRITE(*,'("Could not open input file! ")')
		l_file = .FALSE.
	END IF
END SUBROUTINE OUTFILE

