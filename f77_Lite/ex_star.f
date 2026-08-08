cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE get_array_ave_std(ARRAY, N, MEAN, STD)
        IMPLICIT NONE
      INTEGER N, I
      REAL ARRAY(N), MEAN, STD, SUM, SUMSQ
      
      SUM = 0.0
      SUMSQ = 0.0
      
      do I = 1, N
          SUM = SUM + ARRAY(I)
          SUMSQ = SUMSQ + ARRAY(I)*ARRAY(I)
      enddo
   
      MEAN = SUM / FLOAT(N)
      STD = (SUMSQ - (SUM*SUM)/FLOAT(N)) / FLOAT(N-1)
      STD = SQRT(STD)
      
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ana_chi2_simple(n,map1,map2,p)
      implicit none

      integer n
      real map1(n,n),map2(n,n),p
      integer i,j,n1,n2

      n1=n/4
      n2=(n/4)*3
      p=0.
      do i=n1,n2
        do j=n1,n2
          p=p+(map1(i,j)-map2(i,j))**2
        enddo
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc