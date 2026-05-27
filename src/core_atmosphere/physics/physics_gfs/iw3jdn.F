      INTEGER FUNCTION IW3JDN(IY, IM, ID)

      IMPLICIT NONE

      INTEGER IY, IM, ID
      INTEGER A, Y, M

      A = (14 - IM) / 12
      Y = IY + 4800 - A
      M = IM + 12*A - 3

      IW3JDN = ID + (153*M + 2)/5
     &        + 365*Y + Y/4
     &        - Y/100 + Y/400
     &        - 32045

      RETURN
      END
