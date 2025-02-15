      W1=C1*C3*C10
      W2=C2*C3*C11
      DO 112 K=1,3
         W5=W1*H411(2,2,K)
         W6=W2*H411(2,3,K)
         DO 111  J=1,3
            DO 110  I=1,3
               H123(I,K,J)=-W5*H21(I,J)+W3*H421(I,J,K)
               H432(I,K,J)=-W6*H43(I,J)+W4*H421(I,J,K)
 110        CONTINUE
 111     CONTINUE
 112  CONTINUE
