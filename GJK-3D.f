       subroutine collision_check(c1,v_1,v_2,n1,n2,flag)

!-----------------------------------------collision detection code (gjk algorithm 3D)------------------------------------------------!
!!c1-centroids of two particles
!v_1(d,n1)- coordinates of particle 1
!v_2(d,n2)- coordinates of particle 2
!dimension (d)=3
!n1- number of vertices for particle 1
!n2- number of vertices for particle 2
        
        integer  n1, n2 !number of vertices 
        integer  c, flag, l
        logical  a
        double precision  ra(3,n1), rb(3,n2), c1(3,2), v_1(3,n1)      
        double precision  v(3), p(3), denom, m(3,4), BCD(3,3)
        double precision  proj, ab(3), ao(3), PERP(3,2)
        double precision  ac(3), proj1, proj2, cross1(3)
        double precision  cross_value(3), v_2(3,n2), proj3
        double precision  ad(3), TETR(3,3), TRI(3,4)
                      
            ra(:,:) = v_1(:,:)
            rb(:,:) = v_2(:,:)
            v(:) = c1(:,1) - c1(:,2) 
            c = 0
            denom = (v(1))**2 +(v(2))**2 +(v(3))**2
            v(:) = v(:)/sqrt(denom)
             call support(ra, rb,v, p, c, n1,n2)
            m(:,c) = p(:)
            a = .true.
            v(:) = - p(:) !here we choose a direction in the minkowski difference that points towards the origin

        
        do while(a)

99             call support(ra, rb, v , p, c, n1, n2)
c here c is always 2
                 m(:,c) = p(:) 
                proj = (p(1) * v(1)) + (p(2) * v(2)) 
     & + (p(3) * v(3)) 
     
                if ( proj.lt.0.0) then 
                    flag = 1 	
                    return		
                else
           
                   if(c.eq.2)then
                    ab(:) = (m(:,c-1)-m(:,c)) !ab = b - a
                    ao(:) = - m(:,c) ! ao = o-a vector pointing the origin
                    proj = (ab(1)* ao(1)) + (ab(2) * ao(2))  
     & + (ab(3)* ao(3)) 

                        if(proj.ge.0.0)then 
                            call cross(ab, ao, cross_value)
                            cross1(:)= cross_value(:)
                            call cross(cross1, ab, cross_value)
                            PERP(:,1) = cross_value(:) 
                            v(:) = PERP(:,1)
                    if (PERP(1,1).eq.0.0d0.and.PERP(2,1).eq.0.0d0 
     & .and.PERP(3,1).eq.0.0d0) then
                            flag = 0                        
                         return
                    end if
                        
                        else
                            v(:) =  - m(:,c) 
                            m(:,c-1) = m(:,c) 
                            c = 1  
                                go to 99  
                        end if
			
                else !check for triangle
                        ab(:) = m(:,c-1)- m(:,c)
                        ac(:) = m(:,c-2)- m(:,c)
                        ao(:) = -m(:,c)
					
                        call cross(ab, ac, cross_value) !this gives normal above triangle
                        cross1(:) = cross_value(:)
                        call cross(ab,cross1, cross_value) !this gives perp to ab going away from triangle
                        PERP(:,1) = cross_value(:)
					
                        call cross(ab, ac, cross_value) !this gives normal above triangle
                        cross1(:) = cross_value(:)
                        call cross(cross1, ac, cross_value)!this gives perp to ac going away from triangle
                        PERP(:,2) = cross_value(:) 
							
                        proj1 = (PERP(1,1)*ao(1)) 
     & + (PERP(2,1)*ao(2)) + (PERP(3,1) * ao(3)) 

                        proj2 = (PERP(1,2) * ao(1)) 
     & + (PERP(2,2) * ao(2)) + (PERP(3,2) * ao(3))


						if(proj1.gt.0.0)then 
							proj3 = (ab(1)*ao(1))
     & + (ab(2)* ao(2)) + (ab(3)*ao(3))
     
							if(proj3.gt.0.0)then
                               m(:,1) = m(:,c-1) !c= b
                               m(:,2) = m(:,c) !points are updated as the distant point c is eliminated; b=a
                              call cross(ab,ao,cross_value)
                              cross1(:) = cross_value(:)
                              call cross(cross1,ab,cross_value)
                              v(:) = cross_value(:)
                              m(:,3) = 0.0 
                              m(:,4) = 0.0
                              c = 2
							else
                               v(:) = ao(:)
                               m(:,1) = m(:,c)
                               c = 1
                               go to 99 
							end if
						elseif(proj2.gt.0.0) then 
							proj3 = ( ac(1)*ao(1) ) 
     & + (ac(2)* ao(2)) + (ac(3)* ao(3))

							if (proj3.gt.0.0)then
                                m(:,2) = m(:,c) !b = a
                                m(:,c) = 0.0
                                call cross(ac,ao,cross_value)
                                cross1(:)= cross_value(:)
                                call cross(cross1,ac,cross_value)
                                v(:) = cross_value(:)
                                c = 2
							else
                                v(:) = ao(:)
                                m(:,1) = m(:,c)
                                c = 1
                                go to 99
								
							end if
						
						else
c-----------------------Tetrahedral check------------------------------------------									
c If all the above conditions are satisfied then check above or below triangle 
                TRI(:,1) = m(:,c)
                TRI(:,2) = m(:,c-1)
                TRI(:,3) = m(:,c-2)
                ab(:) = TRI(:,2) - TRI(:,1)
                ac(:) = TRI(:,3) - TRI(:,1)
                call cross(ab,ac, cross_value)
                TETR(:,1) = cross_value(:)
                ao(:) = - m(:,c)
                proj = (TETR(1,1)*ao(1))+ (TETR(2,1) * ao(2)) 
     & + (TETR(3,1)* ao(3))   
            if(proj.gt.0.0) then !above triangle
                v(:) = TETR(:,1)
                call support(ra, rb, v, p, c, n1,n2)
c at this point c = 4
                m(:,c) = p(:) !td(:,1)
                TRI(:,4) = p(:)
                proj= (p(1)*v(1))+ (p(2)*v(2))
     & +(p(3)*v(3))
            
                m(:,1) = TRI(:,3) !d
                m(:,2) = TRI(:,2) !c
                m(:,3) = TRI(:,1) !b
                m(:,4) = TRI(:,4) !a
            if(proj.lt.0.0)then
                flag = 1 !corrected this
                return
            else
                m(:,1) = TRI(:,3) !d
                m(:,2) = TRI(:,2) !c
                m(:,3) = TRI(:,1) !b
                m(:,4) = TRI(:,4) !a
                l = 0
                a = .true. 
                  do while(a)
                   l = l+ 1 
                   ab(:) = m(:,3) - m(:,4)
                   ao(:) = - m(:,4) !new point added
                   ac(:) = m(:,2) - m(:,4)
                   ad(:) = m(:,1) - m(:,4)
				   call cross(ab,ac, cross_value)
                   TETR(:,1) = cross_value(:)
				   proj = (TETR(1,1) * ao(1)) 
     & + (TETR(2,1)* ao(2)) + (TETR(3,1) *  ao(3))
                      if(proj.gt.0.0)then !it is above abc
						c = 3
						m(:,1) = m(:,4)
						v(:) = TETR(:,1)
						call support(ra, rb, v, p, c, n1, n2)
					   m(:,4) = p(:)
					   proj = (p(1)*v(1)) 
     & + (p(2)*v(2))+ (p(3)*v(3)) 
     						if(proj.lt.0.0)then
								flag = 1 !corrected this
								return
							end if
								
                      else
						call cross(ac, ad, cross_value)!check above acd
						TETR(:,2) = cross_value(:)
						proj = (TETR(1,2)* ao(1)) 
     & + (TETR(2,2)* ao(2)) + (TETR(3,2)* ao(3))
                           if(proj.gt.0.0)then !it is above acd
							
							  c = 3
							  m(:,3) = m(:,4)
							  v(:) = TETR(:,2)
							  call support(ra, rb, v, p, c, n1, n2)
                              m(:,4) = p(:) !new vertex generated
							  proj = (p(1)* v(1)) 
     & + (p(2)* v(2)) + (p(3)*v(3))
                                 if(proj.lt.0.0)then
                                     flag = 1
                                     return
                                end if
									
                            else 
							    call cross(ad, ab, cross_value)
								TETR(:,3) = cross_value(:)
								proj = (TETR(1,3)*ao(1)) 
     & + (TETR(2,3)*ao(2)) + (TETR(3,3)*ao(3))
                                 if(proj.gt.0.0)then
                                    c = 3
                                    m(:,2) = m(:,4)
                                    v(:) = TETR(:,3)
                              call support(ra, rb, v, p, c, n1, n2)
                                    m(:,4) = p(:) 
                                    proj = (p(1)* v(1)) 
     & + (p(2)* v(2))+ (p(3)*v(3))
										if(proj.lt.0.0) then
											flag = 1
											return
										end if
                                else
											flag = 0
											return
                                end if
							end if
					  end if
					
                    	if(l.eq.20)then
							flag = 0
							return
						end if                
			
			        end do
                    end if
        else !below triangle
            v(:) = -TETR(:,1)
            call support(ra, rb, v, p, c, n1, n2)
            m(:,c) = p(:) !td(:,1)
            TRI(:,4) = p(:)
            proj= (p(1)*v(1))+ (p(2)*v(2)) 
     & +(p(3)*v(3))
             m(:,1) = TRI(:,3) !c
             m(:,2) = TRI(:,2) !d
             m(:,3) = TRI(:,1) !b
             m(:,4) = TRI(:,4) !a
               if(proj.lt.0.0)then
                    flag = 1
                    return
               else
                    m(:,1) = TRI(:,3) !c
                    m(:,2) = TRI(:,2) !d
                    m(:,3) = TRI(:,1) !b
                    m(:,4) = TRI(:,4) !a
                    a= .true.
                    l = 0			
                 do while(a)
                    l= l+1 !iteration
                    ab(:) = m(:,3) - m(:,4)
                    ao(:) = - m(:,4) !new point added
                    ac(:) = m(:,1) - m(:,4)
                    ad(:) = m(:,2) - m(:,4)
                    call cross(ab,ac, cross_value)
                    TETR(:,1) = cross_value(:)
                    proj = (TETR(1,1) * ao(1)) 
     & + (TETR(2,1)* ao(2))+ (TETR(3,1) *ao(3)) 

                    if(proj.gt.0.0)then !above abc 
                       m(:,2) = m(:,4)
					   c = 3
					   v(:) = TETR(:,1)
                       call support(ra, rb, v, p, c, n1, n2)

						m(:,4) = p(:)
                     proj = (p(1)*v(1))
     & + (p(2)*v(2)) + (p(3)*v(3)) 

                          if(proj.lt.0.0)then
    							flag = 1
								return
     						end if
						
                     else
						call cross(ac, ad, cross_value)!check above acd
						TETR(:,2) = cross_value(:)
						proj = (TETR(1,2)* ao(1)) 
     & + (TETR(2,2)* ao(2))+ (TETR(3,2)* ao(3))
                   	 
        					if(proj.gt.0.0)then 
							   c = 3
                               m(:,3) = m(:,4)
							   v(:) = TETR(:,2)
							   call support(ra, rb, v, p, c, n1, n2)
       							m(:,4) = p(:) 
								proj = (p(1)* v(1)) 
     & + (p(2)* v(2))+ (p(3)*v(3))
								if(proj.lt.0.0)then
									flag = 1
									return
								end if
									
							else 
								call cross(ad, ab, cross_value) 
								TETR(:,3) = cross_value(:)
								proj = (TETR(1,3)*ao(1)) 
     & + (TETR(2,3)*ao(2)) + (TETR(3,3)*ao(3))
								if(proj.gt.0.0)then
                                   c = 3
                                   m(:,1) = m(:,4)
                                   v(:) = TETR(:,3)
                                call support(ra, rb, v, p, c, n1, n2)
    
    								m(:,4) = p(:) 
									proj = (p(1)* v(1)) 
     & + (p(2)* v(2)) + (p(3)*v(3))

                                      if(proj.lt.0.0) then
											flag = 1
											return
									   end if
                                else
										flag = 0 
										return
								end if
									
						  end if
							
					  end if
                          if(l.eq.20)then
                             flag = 0
                             return
                          end if
			
                    end do
                end if
         end if
	
  					end if
					
                end if
                    a =.true.
c further test for collision
		
                end if

           end do

        return
        end subroutine
                
        subroutine support(ra, rb, v, p, c, n1, n2)
c this part of the code was checked, it gives the maximum difference and points in the minkowski diff
        integer n1,n2 
        double precision  ra(3,n1), rb(3,n2), v(3), p(3), v1(3)
        double precision  M1, M2, denom
        integer  j, l,i
        integer  k, m, c
        double precision  s1(n1) , s2(n2)

        do i=1, n1
            s1(i) = v(1)*ra(1,i)+v(2)*ra(2,i)
     &             + v(3)*ra(3,i)
        end do

        M1 = s1(1)
        j = 1
        do i = 2, n1
              if (s1(i).gt.M1)then
                 M1 = s1(i)
                 j = i
              end if
        end do  
        k = j
        v1(:) = - v(:)
        do i=1 , n2
            s2(i)= v1(1)*rb(1,i)+v1(2)*rb(2,i)
     &            +v1(3)*rb(3,i)
        end do

          M2 = s2(1)
          l = 1
          do i = 2, n2
                if (s2(i).gt.M2)then
                    M2 = s2(i)
                    l = i
                end if
          end do
           m = l
           p(:) = ra(:,k) - rb(:,m)
           denom = (p(1))**2.0 + (p(2))**2.0 + (p(3))**2.0
           p(:) = p(:)/sqrt(denom)
           c = c + 1
c Unit vectors are chosen to minimize errors (numerical instability)
        return
        end subroutine


       subroutine cross(u,v,cross_value)
       double precision :: cross_value(3), u(3), v(3)
         cross_value(1)= (u(2)*v(3)) - (u(3)*v(2))
         cross_value(2)= (u(3)*v(1)) - (u(1)*v(3))
         cross_value(3)= (u(1)*v(2)) - (u(2)*v(1))
        return
        end subroutine  

