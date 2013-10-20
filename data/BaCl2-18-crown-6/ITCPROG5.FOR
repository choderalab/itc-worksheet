C    Program ITC-PLANNER by Joel Tellinguisen, Vanderbilt University.
C    (September 16, 2011)
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 FORM(80),ANS
      DIMENSION  C(10),SDC(10),AA(10,10),XY(50,5),WT(50),NU(10)
      COMMON /FROG/XY,WT,C,SDC,AA,NPAR,NADJ,NIT,NTR,NU 
      COMMON /PARS/V0,X0,XM0,NDH 
C 
C  Concentrations are in mM, DH in cal/mol, q in microcal, K(working) in 1/mM.
C
      DATA NPAR,NIT,SIGQ,SIGP,SIGV/3,11,.8,.0024,.015/
      DATA XM0NEW,V0,VTR,NTR,NFIX/0.0,1.4,25.0,10,0/
C      OPEN (UNIT=21,FILE='ITC-PRED.CDR')
      OPEN (UNIT=29,FILE='ITC-PRED.LPT')
      OPEN (UNIT=30,FILE='ITC-PRED.DAT',ACCESS='APPEND')
      PRINT 200
  200 FORMAT(' Program for simulating ITC.'/' Enter each requested item 
     1 with <10 keystrokes, and separate multiple entries by COMMAS.'//
     2' A REPEATED PROMPT MEANS PREVIOUS ENTRY WAS IN ERROR.')
      PRINT 201, SIGQ,SIGP,SIGV
  201 FORMAT(/' Var fn params, sq(ucal), sp, sv(ul) = ',F5.2,F7.4,F6.3)
      PRINT *, ' Change ? (Y or N)'
      READ 101, ANS
        IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') THEN
    2     PRINT *,' Enter new sq, sp, sv'
          READ (5,103,ERR=2) SIGQ,SIGP,SIGV
          PRINT 201, SIGQ,SIGP,SIGV
          IF (SIGQ.EQ.0.0.AND.SIGP.EQ.0.0.AND.SIGV.EQ.0.0) THEN
            PRINT *,' At least one of these must be non-zero!'
            GO TO 2
          ENDIF
        ENDIF
      PRINT 202, V0, NTR, VTR
  202 FORMAT(' Cell V0(ml) = ',F5.2,/I5,' Injections of v(ul) = ',F7.2)
      PRINT *, ' Change any of these ? (Y or N)'
      READ 101, ANS
        IF (ANS.NE.'Y'.AND.ANS.NE.'y') GO TO 58
    4 PRINT *,' Enter new V0'
          READ (5,103,ERR=4) V0
C            PRINT *, V0
            IF (V0.EQ.0.0) GO TO 4
    5     PRINT *, ' Enter # injections (max 50)'
          READ (5,103,ERR=5) XNTR
          NTR = XNTR + .1
            IF (NTR.EQ.0.OR.NTR.GT.50) GO TO 5
C      PRINT *, NTR
      PRINT *, ' Stick with constant injection volume ? (Y/N)'
      READ 101, ANS
        IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') THEN
   55   PRINT *,' Enter injection v(ul)'
        READ (5,103,ERR=55) VTR
          IF (VTR.EQ.0.0) GO TO 55
        PRINT 202, V0, NTR, VTR
        ELSE
          IVV = 1
   56   PRINT *,' Variable vi.  Give total titrant volume (ul)'
        READ (5,103,ERR=56) TITVOL
          IF (TITVOL.EQ.0.0) GO TO 56
        PRINT 203, V0,NTR,TITVOL
  203 FORMAT(' V0(ml) = ',F5.2,/I5,' Injections, total v(ul) = ',F7.2)
        ENDIF
   58 PRINT *,' The default is to fit 3 parameters.  Freeze n ? (Y/N)'
      READ 101, ANS
        IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') NFIX = 1
          NADJ = NPAR-NFIX
        IF (NFIX.NE.1) GO TO 68
      PRINT *,' Use variable-v scheme w/ exponential v distribution ?'
      READ 101, ANS
        IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') THEN
          IF (NTR.GT.10) THEN
    6     PRINT *,' Enter # injections < 11.'
          READ (5,103,ERR=6) XNTR
          NTR = XNTR + .1
            IF (NTR.EQ.0) GO TO 6
          ENDIF
          NVV = 1
C   62     PRINT *,' Enter total titrant volume (ul).'
C          READ (5,103,ERR=62) TITVOL
          BC = 4.5/NTR - .27
          VTOT = 0.0
        DO 65 I=1,NTR
   65     VTOT = VTOT + EXP(BC*I)
          CVOL = TITVOL/VTOT
          GO TO 7
        ENDIF
        GO TO 7
   68 PRINT *,' Fit n x DeltaH instead of DeltaH ? (Y/N)'
      READ 101, ANS
        IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') NDH = 1
    7 DO 8  I = 1,NADJ   
    8     NU(I) = I 
      PRINT *, ' Enter title (max 80 characters).'
      READ 101, FORM
      PRINT 101, FORM
      WRITE (29,101) FORM
      WRITE (29,110) SIGQ,SIGP,SIGV
      WRITE (29,111) V0
   11 PRINT *,' Enter concentration (mM) of titrand M '
      READ (5,103,ERR=11) XM0
        IF (XM0.EQ.0.0) GO TO 11
   12 PRINT *, ' Enter values of n, K(M^-1), DeltaH (cal/mol) '
      READ (5,103,ERR=12) SN,EK,DH
        IF (SN.EQ.0.0.OR.EK.EQ.0.0.OR.DH.EQ.0.0) GO TO 12
          EK = EK/1000.
          ECR = VTR/V0/1000.
        IF (NVV.EQ.1.OR.IVV.EQ.1) ECR = (TITVOL/NTR)/V0/1000.
          FCR = 1. - ECR
          TT = FCR**NTR
   13     ECK = EK*XM0
          TR = 13./ECK + 6.4/ECK**0.2
        IF (TR.LT.1.1) TR=1.1
          X0 = TR*XM0*TT/(1.-TT)
      PRINT 107, ECK,TR,X0
      READ 101, ANS
        IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') THEN
   14     PRINT *, ' Enter new [X]0,0 (mM)'
          READ (5,103,ERR=14) X0
        ENDIF
      WRITE (29,106) X0,XM0
        IF (XM0NEW.GT.0.001) GO TO 18
          VSUM = 0.0   
      DO 16 I=1,NTR  
          XY(I,1) = VTR
        IF (NVV.EQ.1) THEN
          XY(I,1) = CVOL*EXP(BC*I)
          GO TO 152
        ENDIF
        IF (IVV.EQ.1) THEN
   15     PRINT *,' Enter v for injection # ',I
          READ (5,103,ERR=15) XY(I,1)
          VSUM = VSUM + XY(I,1)
        IF (VSUM.GT.TITVOL) THEN
          PRINT 151, VSUM,TITVOL
  151 FORMAT(' Total titrant vol of ',F7.2,'  exceeds stated value,',
     * F7.2,'  Continue ? (Y/N) ')
          READ 101, ANS
          IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') GO TO 152
          STOP
          ENDIF
        ENDIF
  152     XY(I,4) = 0.0 
          WT(I) = 1.0    
C          NTR = I   
   16 CONTINUE
        GO TO 18
   17 PRINT *, ' Enter new [M]0 and [X]0 (mM)'
      READ (5,103,ERR=17) XM0NEW,X0
      WRITE (29,106) X0,XM0NEW
        IF (ABS(XM0NEW-XM0).LT.0.00001) GO TO 18
      PRINT *, ' [M]0 has changed; recalculate Rm ? (Y/N)'
      READ 101, ANS
        IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') THEN
          XM0 = XM0NEW
          GO TO 13
        ENDIF   
   18     C(1) = EK 
          C(2) = DH
        IF (NDH.EQ.1) THEN
C          C(1) = -EK
          C(2) = DH*SN
        ENDIF 
          C(3) = SN 
C          SIGR = SIGV/VTR
          QTOT = 0.0
      PRINT *,'  #     v        q       sig       wt'
      DO 20 I=1,NTR
          XY(I,4) = 0.0
          VTR = XY(I,1)
          YTT = FUNC(NPAR,C,XY,I)
C      PRINT *, VTR,YTT 
          SIGR = SIGV/VTR
          VAR = SIGQ**2 + (SIGR*YTT)**2 + (SIGP*YTT)**2
          SIGI = SQRT(VAR)
C      PRINT *, VAR
          WT(I) = 1.0/VAR
      PRINT 105, i,XY(I,1),ytt,SIGI,WT(I)
          XY(I,4) = YTT
          QTOT = QTOT + YTT
   20 continue
      CALL FIT1(NDH,VAR)
          TRX = XY(NTR,2)/XY(NTR,3)
      PRINT 108, ECK,TRX,QTOT
      WRITE (29,108) ECK,TRX,QTOT
      WRITE (30,109) NTR,TRX,XM0,ECK,(C(I),SDC(I),I=1,3)
      PRINT *, ' Repeat but change concentrations ? (Y/N)'
      READ 101, ANS
        IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') GO TO 17
      PRINT *, ' Repeat for new parameters ? (Y or N)'
      READ 101, ANS
        IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') GO TO 11
      STOP
  101 FORMAT(80A1)  
  102 FORMAT(16I5)
  103 FORMAT(8F10.0)
  104 FORMAT(I5,3F10.5)
  105 FORMAT(I4,F8.2,3F9.2,F9.2)   
  106 FORMAT(/1X,'Initial concs (mM) of titrant and macrom:  '1P2E11.3/)
  107 FORMAT(' c (= K[M]0) = ',G12.3,'  giving Rm = ',G12.3,
     1    ' and [X]0,0 = ',G13.4/' Change [X]0,0 ? (Y or N)')
  108 FORMAT(/' c (= Kx[M]0) = ',G10.3/
     1 ' Actual titration range (Rm) = ',F8.2/
     2 ' Total q(ucal) = ',F8.2)
  109 FORMAT(I5,F10.2,2G10.3,3(5X,G11.4,G10.3))  
  110 FORMAT(5X,' sq(ucal), sp, sv(ul) = ',F5.2,F7.4,F6.3)  
  111 FORMAT(5X,' Active volume V0(ml) = ',F6.2)
  112 FORMAT(5A1,I3,F7.3,2F6.3,2(F8.1,F6.1),2F7.3,F8.4,F7.4,F8.3)
  116 FORMAT(/5X,' # divergences = ', I5//)  
  117 FORMAT(F5.0,F15.5,F10.5)
      END 

      SUBROUTINE INVERT(N,A,AIN)   
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(10,10),AIN(10,10),STOR(10,10)   
      DO 10 I=1,N   
      DO 10 J=1,N   
   10 STOR(I,J)=A(I,J)   
      DO 11 I=1,N   
      DO 11 J=1,N   
   11 AIN(I,J)=0.0  
      DO 12 I=1,N   
   12 AIN(I,I)=1.0  
      DO 21 I=1,N   
      ATEM=STOR(I,I)
      DO 22 J=1,N   
      STOR(I,J)=STOR(I,J)/ATEM
   22 AIN(I,J)=AIN(I,J)/ATEM  
      DO 23 K=1,N   
      IF((K-I).EQ.0) GO TO 23 
      BTEM=STOR(K,I)
      DO 24 J=1,N   
      STOR(K,J)=STOR(K,J)-BTEM*STOR(I,J)
   24 AIN(K,J)=AIN(K,J)-BTEM*AIN(I,J)   
   23 CONTINUE 
   21 CONTINUE 
      RETURN   
      END 

      DOUBLE PRECISION FUNCTION FUNC(N,C,X,KK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  C(N),X(50,5) 
      COMMON /PARS/V0,X00,Y00,NDH 
C          VI = X(KK,1)/1000. 
C          V1 = (X(1,1) + C(4))/1000.
          V1 = X(1,1)/1000.     
          ECR = V1/V0    
          FCR = 1. - ECR 
          X0 = X00*ECR   
          X(1,2) = X0    
          Y0 = Y00*FCR   
          X(1,3) = Y0    
          Q = X(KK,4)
C        IF (C(1).LT.0.0) THEN
C          AK = -C(1)
C          DH = C(2)/C(3)
C        ELSE    
          AK = C(1) 
          DH = C(2)
        IF (NDH.EQ.1) DH = C(2)/C(3)
C        ENDIF 
        IF (KK.NE.1) GO TO 5  
          Y0 = Y0*C(3)   
          TSM = X0 + Y0 + 1./AK    
          TMP = TSM**2 - 4.*Y0*X0  
        IF (TMP.LT.0.0) THEN  
C      PRINT *, 'Neg TMP for MM = ',MM
          TMP = -TMP
        ENDIF  
          XC = .5*(TSM - SQRT(TMP))
C      FUNC = XC*V0*DH - Q + C(5)
      FUNC = XC*V0*DH - Q
          X(1,5) = XC
      RETURN   
    5 DO 10 I=2,KK  
          VTR = X(I,1)/1000.  
          ECR = VTR/V0   
          FCR = 1. - ECR 
          X(I,2) = X(I-1,2)*FCR + X00*ECR    
          X(I,3) = X(I-1,3)*FCR    
   10 CONTINUE 
          KM = KK-1 
          XM = X(KM,2)   
          YM = X(KM,3)   
          YM = YM*C(3)   
          TSM = XM + YM + 1./AK    
          TMP = TSM**2 - 4.*YM*XM  
        IF (TMP.LT.0.0) THEN  
C      PRINT *, 'Neg TMP for MM = ',MM
          TMP = -TMP
        ENDIF  
          XCM = .5*(TSM - SQRT(TMP))    
          VI = X(KK,1)/1000.  
          XK = X(KK,2)   
          YK = X(KK,3)   
          YK = YK*C(3)   
          TSM = XK + YK + 1./AK    
          TMP = TSM**2 - 4.*YK*XK  
        IF (TMP.LT.0.0) THEN  
C      PRINT *, 'Neg TMP for MM = ',MM
          TMP = -TMP
        ENDIF  
          XCK = .5*(TSM - SQRT(TMP))    
C          CALC = (XCK*V0 - XCM*(V0-VI))*DH + C(5)
          CALC = (XCK*V0 - XCM*(V0-VI))*DH
      FUNC = CALC - Q 
          X(KK,5) = XCK   
      RETURN   
      END 

      SUBROUTINE FIT1(NDH,VAR)    
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  AA(10,10),VEC(50,11),XY(50,5),WT(50),DD(10),    
     *    DIFF(50),C(10),SDC(10),DC(10),CS(10),NJ(10),COR(10,10)
      COMMON /FROG/XY,WT,C,SDC,AA,NPAR,NADJ,NIT,NPT,NJ 
C    
C    Iterate NIT times over Loop 70.    
C
c        IF (C(1).GT.0.0) WRITE (29,110)
C        IF (C(1).LT.0.0) WRITE (29,109)
        IF (NDH.EQ.1) THEN
          WRITE (29,109)
        ELSE
          WRITE (29,110)
        ENDIF
      PRINT *,'    K(1/mM)     DH(cal/mol)      n '    
          SSQ = 0.0 
      DO 70  KI = 1,NIT  
          KIT = KI  
      DO 30  I = 1,NPT   
          IQ = I    
        IF (KI.NE.1) GO TO 17 
          T = FUNC(NPAR,C,XY,IQ)
          VEC(I,NADJ+1) = T   
C      PRINT *, I,T 
          SSQ = SSQ + WT(I)*T*T    
   17 CONTINUE 
      DO 20  K = 1,NADJ  
          J = NJ(K) 
          CJ = C(J) 
          C(J) = .999995D0*CJ   
        IF (CJ.EQ.0.0) C(J) = -.000005  
          F1 = FUNC(NPAR,C,XY,IQ)
          C(J) = 1.000005D0*CJ  
        IF (CJ.EQ.0.0) C(J) = .000005   
          F2 = FUNC(NPAR,C,XY,IQ)
          C(J) = CJ 
          CD = .00001*CJ 
        IF (CJ.EQ.0.0) CD=.00001   
   20     VEC(I,K) = (F2-F1)/CD    
   30 CONTINUE 
      DO 35  I = 1,NADJ  
      DO 35  J = 1,NADJ  
   35     AA(I,J) = 0.0  
      DO 45  I = 1,NADJ  
      DO 45  J = 1,NADJ  
        IF (J.LT.I)  GO TO 44 
          T = 0.0   
      DO 40  K = 1,NPT   
   40     T = T + VEC(K,I)*VEC(K,J)*WT(K)    
          AA(I,J) = T    
        GO TO 45    
   44     AA(I,J) = AA(J,I)   
   45 CONTINUE 
      DO 50  I = 1,NADJ  
          T = 0.0   
      DO 48  K = 1,NPT   
   48     T = T + VEC(K,I)*VEC(K,NADJ+1)*WT(K)    
   50     DC(I) = T 
      CALL INVERT(NADJ,AA,AA) 
      DO 51  K = 1,NADJ  
          I = NJ(K) 
          CS(I) = C(I)   
   51     DD(K) = DC(K)  
      DO 54  I = 1,NADJ  
          T = 0.0   
      DO 52  J = 1,NADJ  
   52     T = T + AA(I,J)*DD(J)    
   54     DC(I) = T 
          T1 = SSQ  
          FACT = 1.0
  541 DO 55  K = 1,NADJ  
          I = NJ(K) 
   55     C(I) = CS(I) - DC(K)*FACT
          SSQ = 0.0 
      DO 60  I = 1,NPT   
          IQ = I    
          T = FUNC(NPAR,C,XY,IQ)
          DIFF(I) = T    
          VEC(I,NADJ+1) = T   
   60     SSQ = SSQ + WT(I)*T*T
        IF (SSQ.LT.1.D-7) GO TO 65    
        IF (SSQ.GT.T1) THEN   
          FACT = FACT/2. 
          GO TO 541 
        ENDIF  
   64     TEST = ABS((T1-SSQ)/SSQ) 
   65 PRINT 104,  (C(I),I=1,NPAR)   
        IF (SSQ.LT.1.D-7) GO TO 72    
        IF (TEST.LT.1.E-7)  GO TO 72
   70 CONTINUE 
      PRINT *, 'No convergence after NIT iterations.'
   72     VAR = 1.
          NOF = NPT-NADJ
        IF (NOF.GT.0)  VAR = SSQ/NOF
      DO 76  K = 1,NADJ  
          I = NJ(K) 
   76     SDC(I) = SQRT(AA(K,K))
      PRINT *,'   Standard Errors'    
      PRINT 104,  (SDC(I),I=1,NPAR)
        IF (C(1).LT.0.0) THEN
          AA(1,2) = -AA(1,2)
          AA(2,1) = -AA(2,1)
          AA(1,3) = -AA(1,3)
          AA(3,1) = -AA(3,1)
        ENDIF
          C(1) = ABS(C(1))*1000.
          SDC(1) = SDC(1)*1000.
      WRITE (29,104)  (C(I),I=1,NPAR)
      WRITE (29,111)                                             
      WRITE (29,104)  (SDC(I),I=1,NPAR)
      WRITE (29,114)
      DO 78 I=1,NADJ
   78 WRITE (29,104) (AA(I,J),J=1,NADJ)
      WRITE (29,115)
      DO 80 I=1,NADJ
      DO 80 J=1,NADJ
        IF (J.EQ.I)  GO TO 80 
          COR(I,J) = AA(I,J)/SQRT(AA(I,I)*AA(J,J))
   80 CONTINUE
      DO 82 I=1,NADJ
   82     COR(I,I) = SDC(I)/ABS(C(I))
      DO 84 I=1,NADJ
   84 WRITE (29,104) (COR(I,J),J=1,NADJ) 
      WRITE (29,112)      
      DO 85 I=1,NPT
          FRI = XY(I,5)/XY(I,3)/C(3)
          SIGI = 1./SQRT(WT(I)) 
   85 WRITE (29,105) I,(XY(I,J),J=1,3),XY(I,5),FRI,XY(I,4),SIGI,WT(I)   
      RETURN   
  104 FORMAT(1P8E14.5)  
  105 FORMAT(0P,I4,F7.2,1P3E11.3,0PF7.3,F9.2,F6.2,1PG12.3)
  109 FORMAT('     K(M^-1) n x DeltaH(cal/mol)   n   ')
  110 FORMAT('     K(M^-1)  Delta H (cal/mol)    n   ')
  111 FORMAT(' Standard errors')
  112 FORMAT(/' Calculated concs(mM), f, qs & sigs(ucal)',/'   #  v(ul) 
     1    [X]0       [M]0      [MX]     f(conv)    q     sig     wt')
  114 FORMAT(/' Variance-Covariance matrix (K in 1/mM)')    
  115 FORMAT(/' Correlation coefficients and relative errors') 
      END 

