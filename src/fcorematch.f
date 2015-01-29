 
C     **DELETE THE FIRST PART OF MONIQUEB.F TO MAKE IT **
C     **A PURE SUBROUTINE                              **

      SUBROUTINE SAP (N,M,CC,NBL,INDEX,ZFW,NMATCH,BASIS,MEM,KA,KB,
     F               SM,TMA,TMB,Y1,Y2,DPLUS,DMINUS,SUP,EPS)
C ** *****************************************************************
C     *                                                              *
C     *			FORTRAN-CODE FOR SOLVING		     *
C     *		   THE MIN-COST-PERFECT	MATCHING-PROBLEM	     *
C     *								     *
C ** *****************************************************************
C     *								     *
C     *	 1. CALL:						     *
C     *	    CALL SAP (N,M,CC,NBL,INDEX,ZFW,NMATCH,BASIS,MEM,KA,KB,   *
C     *		      SM,TMA,TMB,Y1,Y2,DPLUS,DMINUS,SUP,EPS)	     *
C     *								     *
C     *	 2. LANGUAGE:						     *
C     *	    FORTRAN IV						     *
C     *								     *
C     *	 3. METHOD:						     *
C     *	    SHORTEST AUGMENTING	PATH METHOD			     *
C     *								     *
C     *	 4. PARAMETER:						     *
C     *								     *
C     *	    INPUT:						     *
C     *	       N	  NUMBER OF NODES (EVEN)		     *
C     *	       M	  NUMBER OF EDGES			     *
C     *	       EPS	  MACHINE ACCURACY			     *
C     *	       SUP	  SUFFICIENTLY LARGE REAL NUMBER	     *
C     *	       NBL(2*M)	  LIST OF NEIGHBOURS			     *
C     *	       CC(2*M)	  COST OF EDGES	ACCORDING TO LIST OF NEIGH-  *
C     *			    BOURS				     *
C     *	       INDEX(I)	  NBL(INDEX(I))	START OF NEIGHBOURLIST OF    *
C     *			    VERTEX I (I=1,..,N)	WITH		     *
C     *			    INDEX(N+1)=2*M+1			     *
C     *								     *
C     *	    OUTPUT:						     *
C     *	       ZFW	 COST OF THE OPTIMAL MATCHING		     *
C     *	       NMATCH(N) OPTIMAL MATCHING			     *
C     *	    INTEGER ARRAYS OF LENGTH N:				     *
C     *	       BASIS,MEM,KA,KB,SM,TMA,TMB			     *
C     *	    REAL*8 ARRAYS OF LENGTH N:				     *
C     *	       Y1,Y2,DPLUS,DMINUS				     *
C     *	    INTEGER ARRAY OF LENGTH N+1:			     *
C     *	       INDEX						     *
C     *	    INTEGER ARRAY OF LENGTH 2*M:			     *
C     *	       CC						     *
C     *	    INTEGER*2 ARRAY OF LENGTH 2*M:			     *
C     *	       NBL						     *
C     *								     *
C     *	 5. EXTERNAL SUBROUTINES :				     *
C     *		 AUGMNT						     *
C     *		 EXPAND						     *
C     *		 GROW						     *
C     *		 OGRAPH						     *
C     *		 SCAN1						     *
C     *		 SCAN2						     *
C     *		 SHRINK						     *
C     *		 START						     *
C     *								     *
C ** *****************************************************************
C
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
      INTEGER    BASIS(N),MEM(N),KA(N),KB(N),ZFW,TOP,NMATCH(N)
      INTEGER    CC(2*M),SM(N),TMA(N),TMB(N),INDEX(N+1)
      INTEGER  NBL(2*M)
      DOUBLE PRECISION Y1(N),Y2(N),DMINUS(N),DPLUS(N),C0,D,DBEST,Y1B,
     *                 Y2B,YB,DFLOAT
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
C *** START
C
          TOP =N+2
          N2  =N/2 
      CALL START(N,NCARD,TOP,CC,NBL,INDEX,NMATCH,Y1)
C     IF (NCARD.EQ.N2)		       GOTO 700
      DO 100      N1=1,N
          BASIS(N1) =N1
          MEM(N1)   =N1
          Y2(N1)    =0.
          SM(N1)    =TOP
          TMA(N1)   =TOP
          TMB(N1)   =TOP
          DPLUS(N1) =SUP
          DMINUS(N1)=SUP
          KA(N1)    =0
          KB(N1)    =N1
  100 CONTINUE
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
C *** INITIALIZATION
C
          NN        =0
      DO 110      NI=1,N
      IF (NMATCH(NI).NE.TOP)            GOTO 110
          NN        =NN+1
          SM(NI)    =0
          DPLUS(NI) =0.
 110  CONTINUE
      IF (NN.LE.1)                      GOTO 700
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
C *** DETERMINATION OF THE NEW DMINUS-VALUES
C
 120  DO 140      N1=1,N
          NB1       =BASIS(N1)
      IF (SM(NB1).NE.0)                 GOTO 140
          Y1B       =Y1(NB1)
          Y2B       =Y2(N1)
          I1        =INDEX(N1)
          I2        =INDEX(N1+1)-1
      DO 130      I3=I1,I2
          N2        =NBL(I3)
          NB2       =BASIS(N2)
      IF (NB1.EQ.NB2)                   GOTO 130
          NC        =CC(I3)
          C0        =DFLOAT(NC)-Y1B-Y2B
          C0        =C0-Y1(NB2)-Y2(N2)
      IF (C0.GE.DMINUS(NB2))            GOTO 130
          KA(NB2)   =N1
          KB(NB2)   =N2
          DMINUS(NB2)=C0
 130  CONTINUE
 140  CONTINUE
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
C *** CONTROL-ROUTINE OF THE PROCEDURE
C
 200   CONTINUE
          DBEST     =SUP
      DO 215      NB=1,N
      IF (BASIS(NB).NE.NB)              GOTO 215
          D         =DMINUS(NB)
      IF (SM(NB).GE.TOP)                GOTO 205
          D         =.5*(D+DPLUS(NB))
      IF (D.GT.DBEST)                   GOTO 215
          NBEST     =NB
          DBEST     =D
      G	O T O 215
 205  IF (TMA(NB).GE.TOP)               GOTO 210
      IF (MEM(NB).EQ.NB)                GOTO 215
          D         =D+Y1(NB)
      IF (D.GE.DBEST)                   GOTO 215
          NBEST     =NB
          DBEST     =D
      G	O T O 215
  210 CONTINUE
      IF (D.GE.DBEST)                   GOTO 215
          NBEST     =NB
          DBEST     =D
  215 CONTINUE
      IF (TMA(NBEST).GE.TOP) GOTO 217
             CALL EXPAND(N,M,TOP,NMATCH,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,
     F                Y1,Y2,DPLUS,DMINUS,SUP,EPS,NBL,INDEX,NBEST,DBEST)
      GOTO 200
  217 IF (SM (NBEST).LT.TOP) GOTO 218
                 CALL GROW(N,TOP,NMATCH,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,
     F                Y1,Y2,DPLUS,DMINUS,SUP,EPS,NBL,INDEX,NBEST,DBEST)
      GOTO 200
  218     NKA       =KA(NBEST)
          NKB       =KB(NBEST)
          N1        =NBEST
          NB1       =N1
          N2        =BASIS(NKA)
          NB2       =N2
  220     TMA(NB1)  =NB2
          NK        =SM(NB1)
      IF (NK.EQ.0)                      GOTO 225
          NB2       =BASIS(NK)
          NB1       =TMA(NB2)
          NB1       =BASIS(NB1)
      G	O T O 220
  225     NB        =NB1
          NB1       =N2
          NB2       =N1
  230 IF (TMA(NB1).LT.TOP)              GOTO 235
          TMA(NB1)  =NB2
          NK        =SM(NB1)
      IF (NK.NE.0) GOTO 232
               CALL AUGMNT(N,TOP,NMATCH,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,
     F                Y1,Y2,DPLUS,DMINUS,SUP,EPS,NBL,INDEX,DBEST,N1,N2,
     F                NKA,NKB,NCARD,J700)
      IF ( J700 .EQ. 1 ) GOTO 700
      GOTO 120
  232     NB2       =BASIS(NK)
          NB1       =TMA(NB2)
          NB1       =BASIS(NB1)
      G	O T O 230
  235 IF (NB1.NE.NB) GOTO 240
               CALL SHRINK(N,TOP,NMATCH,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,
     F                Y1,Y2,DPLUS,DMINUS,SUP,EPS,NBL,INDEX,NBEST,DBEST,
     F                NB,N1,N2,NB2,NKA,NKB)
      GOTO 200
  240     NK        =TMA(NB)
          TMA(NB)   =TOP
          NM        =NMATCH(NK)
          NB        =BASIS(NM)
      G	O T O 235
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
C *** GENERATION OF THE	ORIGINAL GRAPH BY EXPANSION OF ALL
C     SHRUNKEN BLOSSOMS
C
  700 CALL OGRAPH(N,ZFW,EPS,INDEX,NBL,CC,SM,TMA,TMB,NMATCH,MEM,
     F                  BASIS,KA,KB,DPLUS,DMINUS,Y1,Y2)
      RETURN
      END
      SUBROUTINE AUGMNT(N,TOP,NMATCH,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,
     F                Y1,Y2,DPLUS,DMINUS,SUP,EPS,NBL,INDEX,DBEST,N1,N2,
     F                NKA,NKB,NCARD,JRET1)
      INTEGER    BASIS(N),MEM(N),KA(N),KB(N),ZFW
      INTEGER    CC(N*(N-1)),SM(N),TMA(N),TMB(N),INDEX(N+1)
      INTEGER    TOP,NMATCH(N)
      INTEGER  NBL(N*(N-1))
      DOUBLE PRECISION Y1(N),Y2(N),DMINUS(N),DPLUS(N),C0,D,DBEST,Y1B,
     *                 Y2B,YB,DFLOAT
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
C *** AUGMENTATION OF THE MATCHING
C     EXCHANGE OF THE MATCHING-	AND NON-MATCHING-EDGES ALONG
C     THE AUGMENTING PATH
C
          JRET1 = 0
          NB        =N1
          NK        =NKA
  605     NB1       =NB
  606     NMATCH(NB1)=NK
          NK        =SM(NB1)
          TMA(NB1)  =TOP
      IF (NK.EQ.0)                      GOTO 607
          NB2       =BASIS(NK)
          NK1       =TMA(NB2)
          NK        =TMB(NB2)
          NB1       =BASIS(NK1)
          NMATCH(NB2)=NK1
      G	O T O 606
  607 IF (NB.NE.N1)                     GOTO 608
          NB        =N2
          NK        =NKB
      G	O T O 605
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
C *** REMOVING ALL LABELS ON NON-EXPOSED BASE NODES
C
  608 CONTINUE
      DO 620      NB=1,N
      IF (BASIS(NB).NE.NB)              GOTO 620
      IF (SM(NB).GE.TOP)                GOTO 610
          D         =DBEST-DPLUS(NB)
          Y1(NB)    =Y1(NB)+D
          SM(NB)    =TOP
      IF (NMATCH(NB).NE.TOP)            GOTO 615
          SM(NB)    =0
          DPLUS(NB) =0.
      G	O T O 616
  610 IF (TMA(NB).GE.TOP)               GOTO 615
          D         =DMINUS(NB)-DBEST
          Y1(NB)    =Y1(NB)+D
          TMA(NB)   =TOP
          TMB(NB)   =TOP
  615     DPLUS(NB) =SUP
  616     DMINUS(NB)=SUP
  620 CONTINUE
      NCARD=NCARD+1
      NDIFF=N-2*NCARD
      IF(NDIFF.GT.1) RETURN
      JRET1=1
      RETURN
      END
      SUBROUTINE EXPAND(N,M,TOP,NMATCH,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,Y1,
     F                  Y2,DPLUS,DMINUS,SUP,EPS,NBL,INDEX,NBEST,DBEST)
      INTEGER    BASIS(N),MEM(N),KA(N),KB(N),ZFW
      INTEGER    CC(N*(N-1)),SM(N),TMA(N),TMB(N),INDEX(N+1)
      INTEGER    TOP,NMATCH(N)
      INTEGER  NBL(N*(N-1))
      DOUBLE PRECISION Y1(N),Y2(N),DMINUS(N),DPLUS(N),C0,D,DBEST,Y1B,
     *                 Y2B,YB,DFLOAT
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
C *** EXPANSION	OF A T-LABELED BLOSSOM
C
  500 CONTINUE
          N1        =MEM(NBEST)
          NB3       =N1
          NKA       =KA(N1)
          NK2       =N1
  505     NK1       =NK2
          NKB       =KB(NK1)
          Y1B       =Y1(NK1)
  510     BASIS(NK2)=NK1
          Y2(NK2)   =Y2(NK2)-Y1B
      IF (NK2.EQ.NKB)                   GOTO 515
          NK2       =MEM(NK2)
      GOTO 510
  515     NK2       =MEM(NKB)
          MEM(NKB)  =NK1
      IF (NK2.NE.NKA)                   GOTO 505
          Y1B       =DPLUS(N1)
          Y1(NBEST) =Y1B
          MEM(NBEST)=NKA
          NK2       =NKA
  520     Y2(NK2)   =Y2(NK2)-Y1B
      IF (NK2.EQ.NBEST)                 GOTO 525
          NK2       =MEM(NK2)
      GOTO 520
  525 CONTINUE
          NK1       =NMATCH(NBEST)
          NB1       =BASIS(NK1)
          NK2       =SM(NB1)
          NB        =BASIS(NK2)
      IF (NB.EQ.NBEST)                  GOTO 545
          NB2       =NB
  530     NK        =TMA(NB2)
          NB1       =BASIS(NK)
      IF (NB1.EQ.NBEST)                 GOTO 535
          NB2       =SM(NB1)
          NB2       =BASIS(NB2)
      GOTO 530
  535     TMA(NB)   =TMA(NBEST)
          TMA(NBEST)=TMB(NB2)
          TMB(NB)   =TMB(NBEST)
          TMB(NBEST)=NK
          NK3       =SM(NB)
          NB3       =BASIS(NK3)
          NK4       =SM(NB3)
          SM(NB)    =TOP
          NMATCH(NB)=NK1
          NB1       =NB3
  540     NK1       =TMA(NB1)
          NK2       =TMB(NB1)
          TMA(NB1)  =NK4
          TMB(NB1)  =NK3
          SM(NB1)   =NK1
          NMATCH(NB1)=NK1
          NB2       =BASIS(NK1)
          NMATCH(NB2)=NK2
          NK3       =SM(NB2)
          SM(NB2)   =NK2
      IF (NB2.EQ.NBEST)                 GOTO 545
          NB1       =BASIS(NK3)
          NK4       =SM(NB1)
          TMA(NB2)  =NK3
          TMB(NB2)  =NK4
      GOTO 540
  545 CONTINUE
          NK2       =TMB(NB)
          NB1       =BASIS(NK2)
          DMINUS(NB1)=DBEST
          N1        =0
      IF (NB1.EQ.NB)                    GOTO 555
          NK1       =TMA(NB1)
          NB3       =BASIS(NK1)
          TMA(NB1)  =TMA(NB)
          TMB(NB1)  =NK2
  550     NK        =SM(NB1)
          SM(NB1)   =TOP
          NB2       =BASIS(NK)
          NK        =TMA(NB2)
          TMA(NB2)  =TOP
          N2        =TMB(NB2)
          TMB(NB2)  =N1
          N1        =NB2
          DPLUS(NB2)=DBEST
          NB1       =BASIS(NK)
          DMINUS(NB1)=DBEST
      IF (NB1.NE.NB)                    GOTO 550
          TMA(NB)   =N2
          TMB(NB)   =NK
          SM(NB)    =TOP
      IF (NB3.EQ.NB)                    GOTO 570
  555     NB1       =0
          NB2       =NB3
  560     NK        =SM(NB2)
          SM(NB2)   =TOP
          TMA(NB2)  =TOP
          TMB(NB2)  =NB1
          NB1       =BASIS(NK)
          NK        =TMA(NB1)
          SM(NB1)   =TOP
          TMA(NB1)  =TOP
          TMB(NB1)  =NB2
          NB2       =BASIS(NK)
      IF (NB2.NE.NB)                    GOTO 560
      CALL SCAN2(NB1,N,SUP,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,
     *            Y1,Y2,DPLUS,DMINUS,NBL,INDEX)
  570 CONTINUE
  575 IF (N1.EQ.0)                      RETURN
          NB        =N1
      CALL SCAN1(NB,N,SUP,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,
     *            Y1,Y2,DPLUS,DMINUS,NBL,INDEX)
          N1        =TMB(NB)
          TMB(NB)   =TOP
      G	O T O 575
      END
      SUBROUTINE GROW (N,TOP,NMATCH,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,Y1,
     F                Y2,DPLUS,DMINUS,SUP,EPS,NBL,INDEX,NBEST,DBEST)
C
      INTEGER    BASIS(N),MEM(N),KA(N),KB(N),ZFW
      INTEGER    CC(N*(N-1)),SM(N),TMA(N),TMB(N),INDEX(N+1)
      INTEGER    TOP,NMATCH(N)
      INTEGER  NBL(N*(N-1))
      DOUBLE PRECISION Y1(N),Y2(N),DMINUS(N),DPLUS(N),DBEST,DFLOAT
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
C *** GROWING AN ALTERNATING TREE BY ADDING TWO	EDGES
C
          TMA(NBEST)=KA(NBEST)
          TMB(NBEST)=KB(NBEST)
          NM        =NMATCH(NBEST)
          NMB       =BASIS(NM)
          DPLUS(NMB)=DBEST
          SM(NMB)   =NMATCH(NMB)
      CALL SCAN1(NMB,N,SUP,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,
     *            Y1,Y2,DPLUS,DMINUS,NBL,INDEX)
      RETURN
      END
      SUBROUTINE OGRAPH(N,ZFW,EPS,INDEX,NBL,CC,SM,TMA,TMB,NMATCH,MEM,
     1                  BASIS,KA,KB,DPLUS,DMINUS,Y1,Y2)
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
      INTEGER    INDEX(N+1),CC(N*(N-1)),SM(N),TMA(N),TMB(N),N,ZFW,
     F           NMATCH(N),MEM(N),BASIS(N),KA(N),KB(N)
      INTEGER  NBL(N*(N-1))
      DOUBLE PRECISION DPLUS(N),DMINUS(N),Y1(N),Y2(N),C0,D,DBEST,YB,
     F                 Y1B,DFLOAT
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
C *** GENERATION OF THE	ORIGINAL GRAPH BY EXPANSION OF ALL
C     SHRUNKEN BLOSSOMS
C
          ZFW       =0
      DO 702     NB1=1,N
      IF (BASIS(NB1).NE.NB1)            GOTO 702
      IF (SM(NB1).LT.0)                 GOTO 702
          N2        =NMATCH(NB1)
          NB2       =BASIS(N2)
          N1        =NMATCH(NB2)
          SM(NB1)   =-1
          SM(NB2)   =-1
          I1        =INDEX(N1)
          I2        =INDEX(N1+1)-1
          DO 9000 I3=I1,I2
          IF(NBL(I3).EQ.N2) GO TO 9001
 9000     CONTINUE
 9001     NC        =CC(I3)
          D         =DFLOAT(NC)-Y1(NB1)-Y1(NB2)
          D         =D-Y2(N1)-Y2(N2)
      IF (DABS(D).GT.EPS)     CALL RPRINT(N1,N2,D)
          ZFW       =ZFW+NC
  702 CONTINUE
      DO 750      N1=1,N
  705     NB        =BASIS(N1)
      IF (NB.EQ.N1)                     GOTO 750
          NK2       =MEM(NB)
          NKA       =KA(NK2)
          NB3       =NK2
          YB        =DPLUS(NK2)
  710     NK1       =NK2
          NKB       =KB(NK1)
          Y1B       =Y1(NK1)
  715     BASIS(NK2)=NK1
          Y2(NK2)   =Y2(NK2)-Y1B
      IF (NK2.EQ.NKB)                   GOTO 720
          NK2       =MEM(NK2)
      GOTO 715
  720     NK2       =MEM(NKB)
          MEM(NKB)  =NK1
      IF (NK2.NE.NKA)                   GOTO 710
          Y1(NB)    =YB
          MEM(NB)   =NKA
          NK2       =NKA
  725     Y2(NK2)   =Y2(NK2)-YB
      IF (NK2.EQ.NB)                    GOTO 730
          NK2       =MEM(NK2)
      GOTO 725
  730     NK        =NMATCH(NB)
          NK1       =BASIS(NK)
          NK1       =NMATCH(NK1)
          NB1       =BASIS(NK1)
      IF (NB.EQ.NB1)                    GOTO 745
          NMATCH(NB1)=NK
          NB3       =TMA(NB1)
          NB3       =BASIS(NB3)
  735     NK3       =SM(NB1)
          NB2       =BASIS(NK3)
          NK1       =TMA(NB2)
          NK2       =TMB(NB2)
          NB1       =BASIS(NK1)
          NMATCH(NB1)=NK2
          NMATCH(NB2)=NK1
          I1        =INDEX(NK1)
          I2        =INDEX(NK1+1)-1
          DO 9002 I3=I1,I2
          IF(NBL(I3).EQ.NK2) GO TO 9003
 9002     CONTINUE
 9003     NC        =CC(I3)
          D         =DFLOAT(NC)-Y1(NB1)-Y1(NB2)
          D         =D-Y2(NK1)-Y2(NK2)
      IF (DABS(D).GT.EPS)     CALL RPRINT(NK1,NK2,D)
          ZFW       =ZFW+NC
      IF (NB1.NE.NB)                    GOTO 735
  740 IF (NB3.EQ.NB)                    GOTO 705
  745     N2        =SM(NB3)
          NB2       =BASIS(N2)
          N3        =SM(NB2)
          I1        =INDEX(N2)
          I2        =INDEX(N2+1)-1
          DO 9004 I3=I1,I2
          IF(NBL(I3).EQ.N3) GO TO 9005
 9004     CONTINUE
 9005     NC        =CC(I3)
          D         =DFLOAT(NC)-Y1(NB2)-Y1(NB3)
          D         =D-Y2(N2)-Y2(N3)
      IF (DABS(D).GT.EPS)     CALL RPRINT(N2,N3,D)
          ZFW       =ZFW+NC
          N3        =TMA(NB2)
          NB3       =BASIS(N3)
      GOTO 740
  750 CONTINUE
      RETURN
C7000 FORMAT(47H0ERROR , THE OPTIMALITY	CONDITIONS ARE VIOLATED,
C    *	     10X,6HEDGE	(,I3,1H,,I3,1H),5X,F15.4)
C
C     THE FOLLOWING FORMAT "7000" CORRECTS THE ORIGINAL FORMAT
C     ABOVE TO RESTRICT THE OUTPUT TO LESS THAN 80 CHARACTERS
C     PER LINE.
C
 7000 FORMAT(47H0ERROR , THE OPTIMALITY CONDITIONS ARE VIOLATED
     *   /1H ,8X,6HEDGE (,I3,1H,,I3,1H),5X,F15.4)
C
      END
C
      SUBROUTINE RPRINT(A,B,D)
      CHARACTER(55) OUTPUT
      INTEGER A,B
      DOUBLE PRECISION D
      OUTPUT='0ERROR , THE OPTIMALITY CONDITIONS ARE VIOLATED AT EDGE'
      CALL INTPR(OUTPUT, 55, [A,B], 2)
      CALL DBLEPR(OUTPUT, 0, D, 1)
      END
C
      SUBROUTINE SCAN1(NB1,N,SUP,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,
     *                  Y1,Y2,DPLUS,DMINUS,NBL,INDEX)
      INTEGER    N,TOP
      INTEGER    BASIS(N),MEM(N),KA(N),KB(N)
      INTEGER    CC(N*(N-1)),SM(N),TMA(N),TMB(N),INDEX(N+1)
      INTEGER  NBL(N*(N-1))
      DOUBLE PRECISION     Y1(N),Y2(N),DPLUS(N),DMINUS(N),D1,D2,C0,
     F                     DFLOAT
C
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
C *** SCANNING OF NODE NB1
C
          TOP       =N+2
          D1        =DPLUS(NB1)-Y1(NB1)
          DMINUS(NB1)=SUP
          D2        =D1-Y2(NB1)
          TMA(NB1)  =0
          IA        =INDEX(NB1)
          IB        =INDEX(NB1+1)-1
      DO 300      IC=IA,IB
          N2        =NBL(IC)
          NB2       =BASIS(N2)
      IF (TMA(NB2).LT.TOP)              GOTO 300
          NC        =CC(IC)
          C0        =DFLOAT(NC)+D2
          C0        =C0-Y1(NB2)-Y2(N2)
      IF (C0.GE.DMINUS(NB2))            GOTO 300
          KA(NB2)   =NB1
          KB(NB2)   =N2
          DMINUS(NB2)=C0
  300 CONTINUE
          N1        =NB1
      G	O T O 315
  305     D2        =D1-Y2(N1)
          I1        =INDEX(N1)
          I2        =INDEX(N1+1)-1
      DO 310      I3=I1,I2
          N2        =NBL(I3)
          NB2       =BASIS(N2)
      IF (TMA(NB2).LT.TOP) GO TO 310
          NC        =CC(I3)
          C0        =DFLOAT(NC)+D2
          C0        =C0-Y1(NB2)-Y2(N2)
      IF (C0.GE.DMINUS(NB2))            GOTO 310
          KA(NB2)   =N1
          KB(NB2)   =N2
          DMINUS(NB2)=C0
  310 CONTINUE
  315     N1        =MEM(N1)
      IF (N1.NE.NB1)                    GOTO 305
          TMA(NB1)  =TOP
      RETURN
      END
      SUBROUTINE SCAN2(NB,N,SUP,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,
     *                  Y1,Y2,DPLUS,DMINUS,NBL,INDEX)
      INTEGER    N,TOP
      INTEGER    BASIS(N),MEM(N),KA(N),KB(N)
      INTEGER    CC(N*(N-1)),SM(N),TMA(N),TMB(N),INDEX(N+1)
      INTEGER  NBL(N*(N-1))
      DOUBLE PRECISION     Y1(N),Y2(N),DMINUS(N),DPLUS(N),D,C0,Y1B,
     F                     Y2B,DFLOAT
C
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
C *** SCANNING OF NODE NB
C
          TOP       =N+2
  300     NB1       =NB
          NB        =TMB(NB1)
          TMB(NB1)  =TOP
          D         =SUP
          NKA       =0
          NKB       =0
          N1        =NB1
          Y1B       =Y1(NB1)
  315 CONTINUE
          Y2B       =Y2(N1)
          I1        =INDEX(N1)
          I2        =INDEX(N1+1)-1
      DO 320      I3=I1,I2
          N2        =NBL(I3)
          NB2       =BASIS(N2)
      IF (SM(NB2).GE.TOP)               GOTO 320
          NC        =CC(I3)
          C0        =DFLOAT(NC)-Y1B-Y2B
          C0        =C0-Y1(NB2)-Y2(N2)
          C0        =C0+DPLUS(NB2)
      IF (C0.GE.D)                      GOTO 320
          NKA       =N2
          NKB       =N1
          D         =C0
  320 CONTINUE
          N1        =MEM(N1)
      IF (N1.NE.NB1)                    GOTO 315
          KA(NB1)   =NKA
          KB(NB1)   =NKB
          DMINUS(NB1)=D
      IF (NB.NE.0)                      GOTO 300
      RETURN
      END
      SUBROUTINE SHRINK (N,TOP,NMATCH,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,
     F                Y1,Y2,DPLUS,DMINUS,SUP,EPS,NBL,INDEX,NBEST,DBEST,
     F                NB,N1,N2,NB2,NKA,NKB)
      INTEGER    BASIS(N),MEM(N),KA(N),KB(N),ZFW
      INTEGER    CC(N*(N-1)),SM(N),TMA(N),TMB(N),INDEX(N+1)
      INTEGER    TOP,NMATCH(N)
      INTEGER  NBL(N*(N-1))
      DOUBLE PRECISION Y1(N),Y2(N),DMINUS(N),DPLUS(N),D,DBEST,Y1B,Y2B,
     F                 YB,DFLOAT
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
C *** SHRINKING	A BLOSSOM
C
  400 CONTINUE
          YB        =Y1(NB)+DBEST-DPLUS(NB)
          Y1(NB)    =0.
          NK1       =NB
  430     Y2(NK1)   =Y2(NK1)+YB
          NK1       =MEM(NK1)
      IF (NK1.NE.NB)                    GOTO 430
          NK        =MEM(NB)
      IF (NB.NE.N2)                     GOTO 436
  435     N2        =N1
          NB2       =TMA(NB)
  436     MEM(NK1)  =NB2
          NM        =NMATCH(NB2)
          SM(NB2)   =NM
          Y1B       =Y1(NB2)+DMINUS(NB2)-DBEST
          NK1       =NB2
  440     NK2       =NK1
          Y2(NK2)   =Y2(NK2)+Y1B
          BASIS(NK2)=NB
          NK1       =MEM(NK2)
      IF (NK1.NE.NB2)                   GOTO 440
          KB(NB2)   =NK2
          Y1(NB2)   =Y1B
          NB1       =BASIS(NM)
          MEM(NK2)  =NB1
          Y1B       =Y1(NB1)+DBEST-DPLUS(NB1)
          NK2       =NB1
  445     NK1       =NK2
          Y2(NK1)   =Y2(NK1)+Y1B
          BASIS(NK1)=NB
          NK2       =MEM(NK1)
      IF (NK2.NE.NB1)                   GOTO 445
          KB(NB1)   =NK1
          Y1(NB1)   =Y1B
      IF (N2.EQ.NB1)                    GOTO 450
          NB2       =TMA(NB1)
          TMA(NB1)  =TMB(NB2)
          TMB(NB1)  =TMA(NB2)
      GOTO 436
  450 IF (N2.EQ.NBEST)                  GOTO 455
          TMA(N2)   =NKB
          TMB(N2)   =NKA
      IF (NB.NE.NBEST)                  GOTO 435
      GOTO 460
  455     TMA(NBEST)=NKA
          TMB(NBEST)=NKB
  460     MEM(NK1)  =NK
          N1        =MEM(NB)
          KA(N1)    =NK
          DPLUS(N1) =YB
          TMA(NB)   =TOP
          DPLUS(NB) =DBEST
      CALL SCAN1(NB,N,SUP,CC,BASIS,MEM,KA,KB,SM,TMA,TMB,
     *            Y1,Y2,DPLUS,DMINUS,NBL,INDEX)
      RETURN
      END
      SUBROUTINE START(N,NCARD,TOP,CE,NB,INDEX,NMATCH,Y1)
C *** ****************************************************************
C     *								     *
C     *								     *
C     *	     DETERMINATION OF AN INITIAL PARTIAL MATCHING	     *
C     *	     AND A DUAL	SOLUTION FOR STARTING THE SHORTEST	     *
C     *	     AUGMENTING	PATH CODE				     *
C     *								     *
C     *								     *
C *** ****************************************************************
C     *	 1.  CALL:						     *
C     *	     CALL START(N,NCARD,TOP,CE,NB,INDEX,NMATCH,Y1)	     *
C     *								     *
C     *	 2.  COMPUTER CODE:					     *
C     *	     FORTRAN IV						     *
C     *								     *
C     *	 3.  METHOD:						     *
C     *	     - 1-SATURATED MATCHING VIA	GREEDY			     *
C     *								     *
C     *	 4.  PARAMETERS:					     *
C     *	     INPUT:						     *
C     *		N	    NUMBER OF NODES			     *
C     *		TOP	    = N+2				     *
C     *		NB(.)	    LIST OF NEIGHBOURS			     *
C     *		CE(.)	    COSTS OF EDGES ACCORDING TO	LIST OF	     *
C     *			      NEIGHBOURS			     *
C     *		INDEX(I)    NB(INDEX(I)) FIRST NEIGHBOUR OF	     *
C     *			      VERTEX I				     *
C     *		INDEX(N+1)  = N*(N-1)+1	  ( FOR	COMPLETE GRAPHS	)    *
C     *								     *
C     *	     OUTPUT:						     *
C     *		NMATCH(.)   INITIAL PARTIAL MATCHING		     *
C     *		Y1(.)	    INITIAL DUAL SOLUTION		     *
C     *		NCARD	    CARDINALITY	OF PARTIAL MATCHING	     *
C     *								     *
C     *	     INTEGER ARRAY OF LENGTH N :			     *
C     *		NMATCH						     *
C     *								     *
C     *	     REAL*8 ARRAY OF LENGTH N :				     *
C     *		Y1						     *
C     *								     *
C     *	     INTEGER ARRAY OF LENGTH N+1 :			     *
C     *		INDEX						     *
C     *								     *
C     *	     INTEGER*2 ARRAY OF	LENGTH N*(N-1) :		     *
C     *		NB						     *
C     *								     *
C     *	     INTEGER ARRAY OF LENGTH N*(N-1) :			     *
C     *		CE						     *
C     *								     *
C     *	 5.  EXTERNAL SUBROUTINES:				     *
C     *		NONE						     *
C     *								     *
C *** ****************************************************************
C
C - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - - -	- - -
      INTEGER    CE(N*(N-1)),CQ,INDEX(N+1),NMATCH(N),TOP
      INTEGER  NB(N*(N-1))
      DOUBLE PRECISION     D,DD,Y1(N),DFLOAT
      DO  10  I=1,N
   10 NMATCH(I)=TOP
      JJCE=INDEX(1)
      CQ=CE(JJCE)
      N1=INDEX(N)-1
      DO 100  J=1,N1
         IF(CQ.GT.CE(J))               CQ=CE(J)
  100 CONTINUE
      D=DFLOAT(CQ)/2.
      DO 110  I=1,N
  110 Y1(I)=D
      NCARD=0
      DO 150  I=1,N
         IF(NMATCH(I).LT.TOP)          GOTO 150
            N1=0
            N2=INDEX(I)
            N3=INDEX(I+1)-1
            JJNB=NB(N2)
            D=DFLOAT(CE(N2))-Y1(JJNB)
            DO 130  IK=N2,N3
               J=NB(IK)
               DD=DFLOAT(CE(IK))-Y1(J)
               IF(DD.GE.D)             GOTO 120
                  N1=J
                  D=DD
                  GOTO 130
  120          IF(DD.GT.D)             GOTO 130
                  IF(NMATCH(J).LT.TOP) GOTO 130
                     N1=J
  130       CONTINUE
            IF(N1.EQ.0)                GOTO 140
               IF(NMATCH(N1).LT.TOP)   GOTO 140
                  NMATCH(I) =N1
                  NMATCH(N1)=I
                  NCARD=NCARD+1
  140    Y1(I)=D
  150 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DFLOAT(IARG)
C FUNCTION TO CONVERT INTEGER ARGUMENT IARG INTO DOUBLE	PRECISION
      DFLOAT=IARG
      RETURN
      END
