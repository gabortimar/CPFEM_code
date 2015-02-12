!***********************************************************************!
!																		!
!       SUBROUTINE ELVALS												!
!-----------------------------------------------------------------------!
!	Evaluation and assembly of elemental quantities for deformation of	!
!	 elasto-viscoplastic material for use with Newton-Raphson scheme	!
!																		!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE ELVALS

	USE FAS_COM			! FASOLT main variables declaration

	IMPLICIT NONE

	INTEGER(4)	in, in2, ip, i, j, ii, jj, k, l, iel, itp 
	INTEGER(4)	nn, ng, ks, ig, ilmnt
	REAL(4)		sf, d, r, esm, gaussm, defl, weight, small

	REAL(4), DIMENSION (1:3, 1:20)	:: xl, dxl, sfd, fl
	REAL(4), DIMENSION ( 1:6)		:: stl, epl, dsl, del, s

	DIMENSION	r( 1:3), sf( 1:20), in( 1:20), d(1:6, 1:6), in2(1:60),	&
					&		ip(1:8), esm( 1:60, 1:60), gaussm( 1:3, 1:3)



	DATA small/ 1.e-7/			! check on element weighting ( volume)

    DO  ilmnt=1,nlmnt						!	Loop on elements

!------------Get element type, nodes and local pointer

		iel= le(ilmnt)
		itp= lt(ilmnt)
		nn = nd(ilmnt)
		ng = np(ilmnt)
		in(1:nn)=ln( 1:nn, ilmnt)

		DO i= 1, nn
			DO j= -2, 0 
				in2( 3*i + j) = 3*in(i) + j
			END DO
		END DO

!----------Initialise relevant arrays

		fl(1:3,1:nn)=0.
        esm( 1:3*nn, 1:3*nn)= 0.

!-------Get coordinates and displacements
		DO i= 1, nn
			ks= 3*in(i) - 3
			DO j= 1, 3
				xl(j,i) = x(ks+j)
				dxl(j,i)= dx(ks+j)
			END DO
		END DO


		IF ( itp .LT. 9) THEN						! 'Tie' ( 0 Gauss points) Element

			SELECT CASE ( itp)
				
				CASE(1,2,3)
					esm(itp,itp)=  efix ; esm(itp,itp+3)= -efix 

				CASE(4,5,6)
					ii= 1 + MOD(itp-3,3)
					esm(ii,ii)=    efix ; esm(ii,ii+3)=   -efix 

					ii= 1 + MOD(ii,3)
					esm(ii,ii)=    efix ; esm(ii,ii+3)=   -efix


			END SELECT

		ELSE										! 'Real' Element (>0 Gauss points)

			ip(1:ng)=lg( 1:ng, ilmnt)	!	set Gauss point numbers

			DO ig= 1, ng		!	Numerically integrate over element


!--------Choose shape function 
				SELECT CASE( itp )

					CASE(12)
						CALL SHAPE_12( ig, xl, sf, sfd, weight )

					CASE(21)
						CALL SHAPE_21( ig, xl, sf, sfd, weight )

					CASE(22)
						CALL SHAPE_22( ig, xl, sf, sfd, weight )

					CASE(23)
						CALL SHAPE_23( ig, xl, sf, sfd, weight )

				END SELECT

				IF ( weight .LT. small )	CALL FAFERR(920)

!--Get values for total strain and stress at the Gauss point
				j = ip(ig)
				epl( 1:6)= ep( 1:6, j)
				stl( 1:6)= st( 1:6, j)
				dsl( 1:6)= ds( 1:6, j)

!--------Get deformation updated strain & stress increments
				CALL DSTRNS(nn, sfd, dxl, del, r, defl)
				CALL DSTRSS(iel,j,defl,epl,del,stl,dsl,s,r,d)

!------Renew global vectors of these, and rotation, quantities
				dep(j)= defl
				de(1:6,j)= del(1:6)
				ds(1:6,j)= dsl(1:6)
				dr(1:3,j)= r(1:3)

!---------Calculate nodal force contribution and add in
				DO i= 1, nn
					fl(1, i)= fl(1, i) +weight*( s(1)*sfd(1,i) +s(6)*sfd(2,i) +s(5)*sfd(3,i) )
					fl(2, i)= fl(2, i) +weight*( s(2)*sfd(2,i) +s(4)*sfd(3,i) +s(6)*sfd(1,i) )
					fl(3, i)= fl(3, i) +weight*( s(3)*sfd(3,i) +s(5)*sfd(1,i) +s(4)*sfd(2,i) )
				END DO

!---------Calculate mechanical stiffness contribution and add in
				DO i= 1, nn
					DO  j= 1, nn

						gaussm(1,1)= ( sfd(1,i)*sfd(1,j)*(   d(1,1) -  s(1)         ) +	&
							&		   sfd(2,i)*sfd(2,j)*(   d(6,6) - (s(1)-s(2))/2.) +	&
							&          sfd(3,i)*sfd(3,j)*(   d(5,5) - (s(1)-s(3))/2.) +	&
							&          sfd(2,i)*sfd(3,j)*(d(5,6) +  s(4)/2.         ) +	&
							&          sfd(3,i)*sfd(1,j)*(d(1,5)				    ) +	&
							&	       sfd(1,i)*sfd(2,j)*(d(1,6)				    ) + &
							&          sfd(2,j)*sfd(3,i)*(d(5,6) +  s(4)/2.         ) +	&
							&          sfd(3,j)*sfd(1,i)*(d(1,5)				    ) +	&
							&	       sfd(1,j)*sfd(2,i)*(d(1,6)				    )   )

						gaussm(2,1)= ( sfd(1,i)*sfd(1,j)*(   d(1,6) -  s(6)			) +	&
							&          sfd(2,i)*sfd(2,j)*(   d(2,6) -  s(6)			) +	&
							&          sfd(3,i)*sfd(3,j)*(   d(4,5) -  s(6)/2.		) +	&
							&		   sfd(2,i)*sfd(3,j)*( d(2,5) -  s(5)/4.	    ) +	&
							&		   sfd(3,i)*sfd(1,j)*( d(1,4) -  s(4)/4.		) +	&
							&		   sfd(1,i)*sfd(2,j)*( d(6,6) - (s(1)+s(2))/4.	) + &
							&		   sfd(2,j)*sfd(3,i)*( d(4,6) -  s(5)/4.	    ) +	&
							&		   sfd(3,j)*sfd(1,i)*( d(5,6) -  s(4)/4.		) +	&
							&		   sfd(1,j)*sfd(2,i)*( d(1,2) - (s(1)+s(2))/4.  )	)

						gaussm(3,1)= ( sfd(1,i)*sfd(1,j)*(   d(1,5) -  s(5)			) +	&
							&		   sfd(2,i)*sfd(2,j)*(   d(4,6) -  s(5)/2.		) +	&
							&		   sfd(3,i)*sfd(3,j)*(   d(3,5) -  s(5)		    ) +	&
							&		   sfd(2,i)*sfd(3,j)*( d(4,5) -  s(6)/4.		) +	&
							&		   sfd(3,i)*sfd(1,j)*( d(1,3) - (s(1)+s(3))/4.  ) + &
							&		   sfd(1,i)*sfd(2,j)*( d(5,6) -  s(4)/4.	    ) + &
							&		   sfd(2,j)*sfd(3,i)*( d(3,6) -  s(6)/4.		) +	&
							&		   sfd(3,j)*sfd(1,i)*( d(5,5) - (s(1)+s(3))/4.	) + &
							&		   sfd(1,j)*sfd(2,i)*( d(1,4) -  s(4)/4.		)	)

						gaussm(1,2)= ( sfd(1,j)*sfd(1,i)*(   d(1,6) -  s(6)			) +	&
							&          sfd(2,j)*sfd(2,i)*(   d(2,6) -  s(6)			) +	&
							&          sfd(3,j)*sfd(3,i)*(   d(4,5) -  s(6)/2.		) +	&
							&		   sfd(2,j)*sfd(3,i)*( d(2,5) -  s(5)/4.	    ) +	&
							&		   sfd(3,j)*sfd(1,i)*( d(1,4) -  s(4)/4.		) +	&
							&		   sfd(1,j)*sfd(2,i)*( d(6,6) - (s(1)+s(2))/4.	) + &
							&		   sfd(2,i)*sfd(3,j)*( d(4,6) -  s(5)/4.	    ) +	&
							&		   sfd(3,i)*sfd(1,j)*( d(5,6) -  s(4)/4.		) +	&
							&		   sfd(1,i)*sfd(2,j)*( d(1,2) - (s(1)+s(2))/4.  )	)

						gaussm(2,2)= ( sfd(1,i)*sfd(1,j)*( d(6,6) - (s(2)-s(1))/2.	) +	&
							&		   sfd(2,i)*sfd(2,j)*( d(2,2) -  s(2)			) +	&
							&		   sfd(3,i)*sfd(3,j)*( d(4,4) - (s(2)-s(3))/2.	) +	&
							&		   sfd(2,i)*sfd(3,j)*( d(2,4)					) + &
							&		   sfd(3,i)*sfd(1,j)*( d(4,6) +  s(5)/2.		) +	&
							&		   sfd(1,i)*sfd(2,j)*( d(2,6)					) + &
							&		   sfd(2,j)*sfd(3,i)*( d(2,4)					) +	&
							&		   sfd(3,j)*sfd(1,i)*( d(4,6) +  s(5)/2.		) +	&
							&		   sfd(1,j)*sfd(2,i)*( d(2,6)					)   )

						gaussm(3,2)= ( sfd(1,i)*sfd(1,j)*( d(5,6) -  s(4)/2.		) +	&
							&		   sfd(2,i)*sfd(2,j)*( d(2,4) -  s(4)			) +	&
							&		   sfd(3,i)*sfd(3,j)*( d(3,4) -  s(4)			) +	&
							&		   sfd(2,i)*sfd(3,j)*( d(4,4) - (s(2)+s(3))/4.	) +	&
							&		   sfd(3,i)*sfd(1,j)*( d(3,6) -  s(6)/4.		) +	&
							&		   sfd(1,i)*sfd(2,j)*( d(2,5) -  s(5)/4.		) + &
							&		   sfd(2,j)*sfd(3,i)*( d(2,3) - (s(2)+s(3))/4.	) +	&
							&		   sfd(3,j)*sfd(1,i)*( d(4,5) -  s(6)/4.		) +	&
							&		   sfd(1,j)*sfd(2,i)*( d(4,6) -  s(5)/4.		)	)

						gaussm(1,3)= ( sfd(1,j)*sfd(1,i)*(   d(1,5) -  s(5)			) +	&
							&		   sfd(2,j)*sfd(2,i)*(   d(4,6) -  s(5)/2.		) +	&
							&		   sfd(3,j)*sfd(3,i)*(   d(3,5) -  s(5)		    ) +	&
							&		   sfd(2,j)*sfd(3,i)*( d(4,5) -  s(6)/4.		) +	&
							&		   sfd(3,j)*sfd(1,i)*( d(1,3) - (s(1)+s(3))/4.  ) + &
							&		   sfd(1,j)*sfd(2,i)*( d(5,6) -  s(4)/4.	    ) + &
							&		   sfd(2,i)*sfd(3,j)*( d(3,6) -  s(6)/4.		) +	&
							&		   sfd(3,i)*sfd(1,j)*( d(5,5) - (s(1)+s(3))/4.	) + &
							&		   sfd(1,i)*sfd(2,j)*( d(1,4) -  s(4)/4.		)	)

						gaussm(2,3)= ( sfd(1,j)*sfd(1,i)*( d(5,6) -  s(4)/2.		) +	&
							&		   sfd(2,j)*sfd(2,i)*( d(2,4) -  s(4)			) +	&
							&		   sfd(3,j)*sfd(3,i)*( d(3,4) -  s(4)			) +	&
							&		   sfd(2,j)*sfd(3,i)*( d(4,4) - (s(2)+s(3))/4.	) +	&
							&		   sfd(3,j)*sfd(1,i)*( d(3,6) -  s(6)/4.		) +	&
							&		   sfd(1,j)*sfd(2,i)*( d(2,5) -  s(5)/4.		) + &
							&		   sfd(2,i)*sfd(3,j)*( d(2,3) - (s(2)+s(3))/4.	) +	&
							&		   sfd(3,i)*sfd(1,j)*( d(4,5) -  s(6)/4.		) +	&
							&		   sfd(1,i)*sfd(2,j)*( d(4,6) -  s(5)/4.		)	)

						gaussm(3,3) =( sfd(1,i)*sfd(1,j)*(   d(5,5) - (s(3)-s(1))/2.) +	&
							&		   sfd(2,i)*sfd(2,j)*(   d(4,4) - (s(3)-s(2))/2.) +	&
							&		   sfd(3,i)*sfd(3,j)*(   d(3,3) -  s(3)         ) +	&
							&		   sfd(2,i)*sfd(3,j)*(d(3,4)                    ) +	&
							&		   sfd(3,i)*sfd(1,j)*(d(3,5) 		            ) +	&
							&		   sfd(1,i)*sfd(2,j)*(d(4,5) +  s(6)/2.			) + &
							&		   sfd(2,j)*sfd(3,i)*(d(3,4)                    ) +	&
							&		   sfd(3,j)*sfd(1,i)*(d(3,5) 		            ) +	&
							&		   sfd(1,j)*sfd(2,i)*(d(4,5) +  s(6)/2.			)   )


						gaussm= weight*gaussm		

						esm( 3*i-2, 3*j-2)= esm( 3*i-2, 3*j-2) + gaussm(1,1)
						esm( 3*i-2, 3*j-1)= esm( 3*i-2, 3*j-1) + gaussm(1,2)
						esm( 3*i-2, 3*j  )= esm( 3*i-2, 3*j  ) + gaussm(1,3)
						esm( 3*i-1, 3*j-2)= esm( 3*i-1, 3*j-2) + gaussm(2,1)
						esm( 3*i-1, 3*j-1)= esm( 3*i-1, 3*j-1) + gaussm(2,2)
						esm( 3*i-1, 3*j  )= esm( 3*i-1, 3*j  ) + gaussm(2,3)
						esm( 3*i,   3*j-2)= esm( 3*i,   3*j-2) + gaussm(3,1)
						esm( 3*i,   3*j-1)= esm( 3*i,   3*j-1) + gaussm(3,2)
						esm( 3*i,   3*j  )= esm( 3*i,   3*j  ) + gaussm(3,3)
					 
					END DO

				END DO

			END DO

		END IF

!-----------------------Assemble mechanical stiffness into global array
		DO i= 1, 3*nn
			ii= in2(i)		! row

			DO j= 1, 3*nn 
				jj= in2(j)		! column number

				IF ( jj .NE. ii ) THEN	! column location
					k=2
					DO WHILE( jj .NE. is(k, ii) )
							k=k+1
					END DO
				ELSE
					k=1 
				END IF

				ss(k,ii)= ss(k,ii) + esm(i,j)	

			END DO
		END DO

!-----------------------------Return nodal force imbalance values
		DO i= 1, nn
			k= 3*in(i) - 3
			DO j= 1, 3
				l= k + j
				fc(l)= fc(l) + fl(j,i)
			END DO
		END DO

	END DO

END SUBROUTINE ELVALS
