      SUBROUTINE MWRAP(N,WT,NMATCH,PRCN)

      INTEGER BASIS(N),MEM(N),KA(N),PRCN
      INTEGER KB(N),ZFW,NMATCH(N),CC(N*(N-1))
      INTEGER SM(N),TMA(N),TMB(N),INDEX(N+1)
      INTEGER NBL(N*(N-1)), WT(N*N)
      INTEGER M1,NN,X,Y,N,M
      DOUBLE PRECISION  Y1(N),Y2(N),DMINUS(N),DPLUS(N)

      M=N*(N-1)/2
      M1=2*M
      NN=0
      DO 25 X=1,N
      DO 35 Y=1,N
      IF(Y.Eq. X) go to 35
      NN=NN+1
      NBL(NN)=Y
      CC(NN)=WT((X-1)*N+Y)
 35   CONTINUE
      INDEX(X)=(X-1)*(N-1)+1
 25   CONTINUE
      INDEX(N+1)=N*(N-1)+1
      EPS=10.**(-38)
      IF(PRCN.lt.6) THEN
      SUP=40000000
      ELSE
      SUP=4*10.**(PRCN+1)
      ENDIF
      CALL SAP(N,M,CC,NBL,INDEX,ZFW,
     *NMATCH,BASIS,MEM,KA,KB,SM,
     *TMA,TMB,Y1,Y2,DPLUS,DMINUS,
     *SUP,EPS)

      END SUBROUTINE
