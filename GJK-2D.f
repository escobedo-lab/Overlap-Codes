	subroutine collision_check(c1,v_1,v_2,n1,n2,flag)
!-----------------------------------------collision detection code (
!GJK algorithm for 2D)------------------------------------------------!
!c1-centroids of two particles
!v_1(d,n1)- coordinates of particle 1
!v_2(d,n2)- coordinates of particle 2
!dimension (d)=2
!n1- number of vertices for particle 1
!n2-number of vertices for particle 2


        integer  n1, n2 !number of vertices 
        integer  c, flag, l, t, k
        logical  a
        double precision  ra(2,n1), rb(2,n2), c1(2,2), v_1(2,n1)
        double precision  v(2), p(2), denom, m(2,4)
        double precision  proj, ab(2), ao(2), PERP(2,2)
        double precision dot1, dot2
        double precision  ac(2), proj1, proj2
        double precision   v_2(2,n2), proj3
!            write(*,*) v_1(:,:), 'v1'
!            write(*,*) v_2(:,:), 'v2'
!            write(*,*) c1, 'c1'
          do k = 1,n1
            ra(:,k) = v_1(:,k)
         end do

          do k = 1, n2
            rb(:,k) = v_2(:,k)
          end do
            v(:) = c1(:,1) - c1(:,2)
            c = 0
            denom = (v(1))**2 +(v(2))**2
            v(:) = v(:)/sqrt(denom)
             call support(ra, rb,v, p, c, n1,n2)
            m(:,c) = p(:)
            a = .true.
            v(:) = - p(:) !here we choose a direction in the minkowskidifference that points towards the origin

        t = 0
        do while(a)

99             call support(ra, rb, v , p, c, n1, n2)
                t = t + 1
c here c is always 2
                 m(:,c) = p(:)
                proj = (p(1) * v(1)) + (p(2) * v(2))

                if ( proj.lt.0.0) then
                    flag = 1
                    return
                else

                   if(c.eq.2)then !We have a straight line (two points)
                    ab(:) = (m(:,c-1)-m(:,c)) !ab = b - a
                    ao(:) = - m(:,c) ! ao = o-a vector pointing theorigin
                    proj = (ab(1)* ao(1)) + (ab(2) * ao(2))


                        if(proj.ge.0.0)then
                          dot1 = (ab(1)*ab(1)) + (ab(2)*ab(2))
!Dot product of C.dot(A) 
                                            !write(*,*) dot1
                          dot2 = (ab(1)*ao(1))+ (ab(2)*ao(2))
                                                !write(*,*) dot2
                         PERP(1,1) = (ao(1)*(dot1) )
     &                               - ( ab(1)* dot2 )
                         PERP(2,1) = ( ao(2)*dot1  )
     &- ( ab(2)* dot2 )
!here we evaluate cross product to get the direction of the origin 
                            v(:) = PERP(:,1)

                    if (PERP(1,1).eq.0.0d0.and.PERP(2,1).eq.0.0d0) then
!                          write(*,*) 'perp condition'
                            flag = 0
                         return
                    end if

                        else
                            v(:) =  - m(:,c)
                            m(:,c-1) = m(:,c)
!                            write(*,*) c
                            c = 1
!                            write(*,*) proj, v
!                            write(*,*) 'stuck 1',t
                                go to 99
                        end if

                else !check for triangle
                        ab(:) = m(:,c-1)- m(:,c)
                        ac(:) = m(:,c-2)- m(:,c)
                        ao(:) = -m(:,c)


                        dot1 = (ab(1)*ac(1)) + (ab(2)*ac(2))
                        dot2 = (ab(1)*ab(1)) +(ab(2)*ab(2))
                        PERP(1,1) = (ab(1)*dot1) - (ac(1)* dot2)
                        PERP(2,1) = (ab(2)*dot1) - (ac(2)* dot2)



                       dot1 = (ac(1) * ab(1)) + (ac(2) * ab(2))
                       dot2 = (ac(1) * ac(1)) + (ac(2) * ac(2))
                       PERP(1,2) = (ac(1)*dot1) - (ab(1) * dot2)
                       PERP(2,2) = (ac(2)*dot1) - (ab(2) * dot2)


                        proj1 = (PERP(1,1)*ao(1))
     & + (PERP(2,1)*ao(2))

                        proj2 = (PERP(1,2) * ao(1))
     & + (PERP(2,2) * ao(2))


                           if(proj1.gt.0.0)then
                                         proj3 =(ab(1)*ao(1))
     & + (ab(2)* ao(2))

                              if(proj3.gt.0.0)then
                               m(:,1) = m(:,c-1) !c= b
                               m(:,2) = m(:,c) !points are updated asthe distant point c is eliminated; b=a
                          dot1 = (ab(1)*ab(1)) + (ab(2)*ab(2))
!Dot product of C.dot(A) 
                                            !write(*,*) dot1
                          dot2 = (ab(1)*ao(1))+ (ab(2)*ao(2))
                                                !write(*,*) dot2
                         v(1) = ( ao(1)*(dot1) )
     &                        - ( ab(1)* dot2  )
                         v(2) = ( ao(2)* dot1  )
     &                        - ( ab(2)* dot2  )  
                              m(:,3) = 0.0
                              m(:,4) = 0.0
                              c = 2
                                                        else
                               v(:) = ao(:)
                               m(:,1) = m(:,c)
                               c = 1
                              go to 99
                                                        end if
                                                elseif(proj2.gt.0.0)then
                                                 proj3 = (ac(1)*ao(1) )
     & + (ac(2)* ao(2))

                                   if(proj3.gt.0.0)then
                                         m(:,2) = m(:,c) !b = a
                                         m(:,c) = 0.0

                              dot1 = (ac(1)*ac(1)) + (ac(2)*ac(2))
!Dot product of C.dot(A) 
                                            !write(*,*) dot1
                         dot2 = (ac(1)*ao(1))+ (ac(2)*ao(2))
                                                !write(*,*) dot2
                         v(1) = (ao(1)*(dot1) )
     &                        - ( ac(1)* dot2 )
                         v(2) = ( ao(2)*dot1  )
     &                        - ( ac(2)* dot2 )
                          c = 2
                                                        else
                                v(:) = ao(:)
                                m(:,1) = m(:,c)
                                c = 1
!                                write(*,*) 'stuck 2',t
                                go to 99

                                                        end if

                                                else
!                                   write(*,*) 'staisfied'
                                   flag = 0
                                   return
         end if

                                        end if

!               
                    a =.true.
c further test for collision

                end if

           end do

        return
        end subroutine


        subroutine support(ra, rb, v, p, c, n1, n2)
c this part of the code was checked, it gives the maximum difference and points in the minkowski diff
        integer n1,n2
        double precision  ra(2,n1), rb(2,n2), v(2), p(2), v1(2)
        double precision  M1, M2, denom
        integer  j, l,i
        integer  k, m, c
        double precision  s1(n1) , s2(n2)

        do i=1, n1
            s1(i) = v(1)*ra(1,i)+v(2)*ra(2,i)
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
          denom = (p(1))**2.0 + (p(2))**2.0
          if(p(1).eq.0.0d0.and.p(2).eq.0.0d0)then
                STOP 'ZERO, Vertex equal'

          end if
           p(:) = p(:)/sqrt(denom)
           c = c + 1
c Unit vectors are chosen to minimize errors (numerical instability)
        return
        end subroutine

