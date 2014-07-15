      PROGRAM JOUKOWSKI

      INCLUDE 'xcilgetto.com'

      CALL OPFILE(10,'OLD','Nome file di input ?')
c      CALL OPFILEB(25,'NEW','Nome file finale binario?')
c      WRITE(6,*)' Nome file finale binario? '
c      READ(5,*)NAMEBIN
c      WRITE(6,*)NAMEBIN
C     Scritto in Subroutine FORCES
      OPEN(26,FILE='forze')
C     Scritti in Subroutine RETIC            
      OPEN(34,FILE='profilo')       
C     OPEN(35,FILE='retic')
      OPEN(36,FILE='distr_punti')
C     Scritto in Subroutine PUNTOV
C     OPEN(56,FILE='pulsvel')       
C                                    unit=71 etc
C     Scrive in Subroutine OUTPUTINT unit=61 'psiomg-tempo'                   
C                                    unit=62 'veloci-tempo' if ISTAMPAVEL=1
C     Scrive in Subroutine OUTPUTP   unit=45 'pressione-tempo'
C                                    unit=46 'vorticita-tempo'
C     Scrive in Subroutine OUTCMEDIO unit=91 'sciamedia' if ICMEDIO=1
C     Scrive in Subroutine VEL       unit=59 OMG al bordo d'uscita
C                                    unit=58 OMG al bordo d'ingresso
C     Scrive in Subroutine OUTPUTSTR unit 66-70
 
      IUSC1=0
      IUSC2=0
      ITOUT=0
      NIND_PO=1
      NIND_PRES=1

      call cost
      CALL INPUT
      call retic

      CALL INIT

      write(6,*)' dis = ',DIS
      write(6,*)' e = ',E
      write(6,*)' small = ',small
      write(6,*)' r_max = ',Rmax, ' r_int = ',Rint
      write(6,*)' a = ',a
      write(6,*)' NTETA,NZETA = ',nteta,nzeta
      write(6,*)' u_tilde = ',U_TILDE,' Strouhal = ',strouhal
      write(6,*)' beta = ',BETA
      write(6,*)' Alfa_amp = ',Alphaamp*180.d+0/pi
      write(6,*)' Alfa_phi = ',Alphaphi*180.d+0/pi
      write(6,*)' Re = ',Reyn,' RE = ',re
      write(6,*)' Kc = ',Kc
      write(6,*)' Nt ciclo e t fin = ',NTciclo,Tfin
      write(6,*)' passo temporale dt = ',dt 
      write(6,*)' N usc 1, 2, t = ',NUSC1,NUSC2,itusc
      if (icmedio.eq.1) then
        write(6,*)' x-min, x-max, nx ',xhat_min,xhat_max,nxhat
        write(6,*)' y-min, y-max, ny ',yhat_min,yhat_max,nyhat
      endif
c      stop

      l=lfun(t,Kc)
      alpha=alphafun(t,alphaamp,alphaphi)
      dl=dlfun(t,Kc)
      dalpha=dalphafun(t,alphaamp,alphaphi)
      d2l=d2lfun(t,Kc)
      d2alpha=d2alphafun(t,alphaamp,alphaphi) 

      CALL VEL

10    CONTINUE

      ITOUT=ITOUT+1     
      IUSC1=IUSC1+1
      IUSC2=IUSC2+1
 
      if (ITOUT.EQ.ITUSC) then
         write(6,*)'tempo =',t,' omg(nteta/2+1,1) ='
     $        ,omg(nteta/2+1,1,1)
         itout=0
      endif

      IF (T.LE.TFIN) THEN

          CALL VORTICITY   

          T=T+DT
          l=lfun(t,Kc)          
          alpha=alphafun(t,alphaamp,alphaphi)
          dl=dlfun(t,Kc)
          dalpha=dalphafun(t,alphaamp,alphaphi)
          d2l=d2lfun(t,Kc)
          d2alpha=d2alphafun(t,alphaamp,alphaphi)       

          CALL POISSON   
c          CALL BOUNDOMG          ! per la sol. esplicita
          CALL VEL

c          IF (ICMEDIO.EQ.1)
          IF (ICMEDIO.EQ.1.AND.T.GE.2.D+0*PI.AND.T.LT.6.D+0*PI)
     $         CALL CMEDIO

          IF (IUSC1.EQ.NUSC1) THEN
c             CALL OUTPUTstr
             CALL PUNTOV
             CALL FORCES 
             IUSC1=0
          ENDIF
          IF (IUSC2.EQ.NUSC2)THEN
c            CALL OUTPUTstr
C     inserire qui eventuale rewind file finale e sua scrittura
             REWIND 25
             CALL OUTPUT
             IUSC2=0
          ENDIF

          IF (ABS(VET_PO(NIND_PO)-T).LT.DT) then
             CALL OUTPUTINT 
             NIND_PO=NIND_PO+1
          endif
          IF (ABS(VET_PRES(NIND_PRES)-T).LT.DT) then
             call outputp
             NIND_PRES=NIND_PRES+1
          endif
          GOTO 10      
      
      ENDIF

      REWIND 25 

c     call outputint
      CALL OUTPUT
c     CALL OUTPUTstr
      CALL FORCES
c      CALL OUTPUTP

C     Stampa dei risultati sul campo medio
      CALL OUTCMEDIO

C     Chiude i file Subroutine PUNTOV
      IF (IPUNTIVEL.EQ.1) THEN 
         DO J=1,NUM_PUNTI
            CLOSE(UNIT=70+J)
         ENDDO
      ENDIF

c     Calcolo dell'efficienza
      CLOSE (26)                ! forze 
      OPEN(26,FILE='forze')
      OPEN(27,FILE='efficienza')
      CALL EFFIC
      CLOSE (26)                
      CLOSE (27)

c 100  write(6,*)'Problemi di apertura files! STOP '
c      STOP

      CLOSE (10)           ! in.dat
c      CLOSE (25)           ! output    OMG,PSI
      CLOSE (34)           ! profilo
c     CLOSE (35)           ! retic 
      CLOSE (36)           ! distr p di calcolo lungo il raggio
C     CLOSE (56)           ! 'pulsvel'
            
      END

C ***********************************************************************

      SUBROUTINE cost
      INCLUDE 'xcilgetto.com'

      pi=3.141592654d0
      c0=(0.D0,0.D0)
      ci=(0.D0,1.D0)

      return
      end

C ***********************************************************************

      SUBROUTINE INPUT
      INCLUDE 'xcilgetto.com'
      
      READ(10,*)IOPT
      READ(10,*)ICIL        ! 0 - se profilo joukowski; 1 - se cilindro
      READ(10,*)DIS
      READ(10,*)E
      READ(10,*)small
      READ(10,*)STROUHAL
      READ(10,*)BETA
      READ(10,*)Alphaamp
      READ(10,*)Alphaphi
      READ(10,*)RE
      READ(10,*)Kc
      READ(10,*)Rmax
      READ(10,*)Rint
      READ(10,*)NTETA,NZETA
      READ(10,*)NTciclo,Tfin
      READ(10,*)NUSC1,NUSC2,ITUSC
      READ(10,*)NUSC_TETA,NUSC_ZETA
      READ(10,*)ISTAMPAVEL ! 0 non stampa campi di velocita; 1 stampa
      READ(10,*)NUM_PO
      DO 20 J=1,NUM_PO
         READ(10,*)VET_PO(J)
 20   CONTINUE  
      READ(10,*)NUM_PRES
      DO 21 J=1,NUM_PRES
         READ(10,*)VET_PRES(J)
 21   CONTINUE      
      READ(10,*)ICMEDIO ! 1 calcola campo medio su 1 regione; 0 no
      IF (ICMEDIO.EQ.1) THEN 
         READ(10,*)XHAT_MIN,XHAT_MAX,NXHAT
         READ(10,*)YHAT_MIN,YHAT_MAX,NYHAT
      ENDIF
      READ(10,*)IPUNTIVEL ! 1 stampa velocita in alcuni punti
      IF(IPUNTIVEL.EQ.1) THEN 
         READ(10,*)NUM_PUNTI
         DO 22 J=1,NUM_PUNTI
            READ(10,*)VET_Z(J)
            CALL OPFILEPVEL(J,VET_Z(J))
 22      CONTINUE
      ENDIF 
      READ(10,*)NUMJ
      DO 23 J=1,NUMJ
         READ(10,*)INDJ(J)
 23   CONTINUE       

      U_TILDE=1./(PI*STROUHAL)
      REYN=RE*STROUHAL*PI/4.
      ALPHAAMP=ALPHAAMP*PI/180.
      ALPHAPHI=ALPHAPHI*PI/180.

      RETURN
      END

C *****************************************************************************

      SUBROUTINE RETIC
      INCLUDE 'xcilgetto.com'

      real(8)r0,aux1,aux2,dcsidteta,detadeta
      real(8)DXHAT,dyhat

      DTETA=2.D0*PI/FLOAT(NTETA)
      DO I=1,NTETA+1
         TETA(I)=DTETA*FLOAT(I-1)
      ENDDO

      IF (ICIL.EQ.1) THEN 
         E=0.D0
         SMALL=0.D0
         DIS=0.D0
      ENDIF
      
      r0=1.D0+e+small
      a=(rint**2-rmax*r0)/(rmax+r0-2.*rint)
      r(1)=R0
      r2(1)=r(1)**2
      ZETA(1)=DLOG(R(1)+A)
      ez(1)=dexp(zeta(1))
      ez2(1)=ez(1)**2
      ZMAX=DLOG(RMAX+A)
      DZETA=(ZMAX-ZETA(1))/FLOAT(NZETA-1)
 
      DO J=2,NZETA
         ZETA(J)=ZETA(1)+DZETA*FLOAT(J-1)
         ez(j)=exp(zeta(j))
         ez2(j)=ez(j)**2
         r(J)=EXP(ZETA(J))-A
         r2(j)=r(j)**2        
      ENDDO

C     stampa file 'distr_punti'
      do J=1,nzeta          
         write(36,37)j,zeta(j),r(j),dzeta,dteta*r(j)
      enddo
 37   format(i3,4(1x,e13.7))

      IF(ICIL.EQ.0) THEN

      DO I=1,NTETA+1
         DO J=1,NZETA
            CHI(I,J)=CMPLX(R(J)*COS(TETA(I)),R(J)*SIN(TETA(I)))
            call conf(chi(i,j),z(i,j))
            aux1=(real(chi(i,j))-e)**2
            aux2=(aimag(chi(i,j)))**2
            jacob(i,j)=1.D0+(1.D0-2.D0*(aux1-aux2))/(aux1+aux2)**2
            jacobsq(i,j)=sqrt(jacob(i,j))
            dxdcsi(i,j)=1.D0-(aux1-aux2)/(aux1+aux2)**2 
            dxdeta(i,j)=-2.D0*( real(chi(i,j))-e )*( aimag(chi(i,j)) )
     *                         /(aux1+aux2)**2 
            dydcsi(i,j)=-dxdeta(i,j) 
            dydeta(i,j)=dxdcsi(i,j) 
         ENDDO
         dcsidteta=-r(1)*sin(teta(i))        
         detadteta=r(1)*cos(teta(i))        
         aux1=real(chi(i,1))-e
         aux2=aimag(chi(i,1))
         d2xdcsi2(i)=-2.*aux1*(3.*aux2**2-aux1**2)/(aux1**2+aux2**2)**3
         d2xdcsieta(i)=2.*aux2*(3.*aux1**2-aux2**2)/(aux1**2+aux2**2)**3
         d2xdeta2=-d2xdcsi2(i)
         d2xdcsiteta(i)=d2xdcsi2(i)*dcsidteta+d2xdcsieta(i)*detadteta
         d2xdetateta(i)=d2xdcsieta(i)*dcsidteta+d2xdeta2*detadteta
         d2ydcsiteta(i)=-d2xdcsieta(i)*dcsidteta-d2xdeta2*detadteta
         d2ydetateta(i)=-d2xdeta2*dcsidteta+d2xdcsieta(i)*detadteta
         dxdteta(i)=r(1)*(-dxdcsi(i,1)*sin(teta(i))+
     $        dxdeta(i,1)*cos(teta(i)))
         dydteta(i)=r(1)*(-dydcsi(i,1)*sin(teta(i))+
     $        dydeta(i,1)*cos(teta(i)))
      ENDDO

      ELSE
         
      DO I=1,NTETA+1
         DO J=1,NZETA
            CHI(I,J)=CMPLX(R(J)*COS(TETA(I)),R(J)*SIN(TETA(I)))
            Z(I,J)=CHI(I,J)
            JACOB(I,J)=1.d0
            jacobsq(i,j)=1.d0
            dxdcsi(i,j)=1.d0 
            dxdeta(i,j)=0.d0
            dydcsi(i,j)=0.d0
            dydeta(i,j)=1.d0 
         ENDDO
         dcsidteta=-r(1)*sin(teta(i))        
         detadteta=r(1)*cos(teta(i))        
         d2xdcsi2(i)=0.d0
         d2xdcsieta(i)=0.d0
         d2xdeta2=0.d0
         d2xdcsiteta(i)=0.d0
         d2xdetateta(i)=0.d0
         d2ydcsiteta(i)=0.d0
         d2ydetateta(i)=0.d0
         dxdteta(i)=-r(1)*sin(teta(i)) 
         dydteta(i)=r(1)*cos(teta(i))   
      ENDDO

      ENDIF

C     stampa file 'profilo'
      DO 12 I=1,NTETA
         WRITE(34,*)REAL(Z(I,1)),AIMAG(Z(I,1))
 12   CONTINUE

      IF (ICMEDIO.EQ.1) THEN
         DXHAT=(XHAT_MAX-XHAT_MIN)/NXHAT
         DYHAT=(YHAT_MAX-YHAT_MIN)/NYHAT
         XHAT(1)=XHAT_MIN
         YHAT(1)=YHAT_MIN
         DO I=2,NXHAT+1
            XHAT(I)=XHAT(1)+DXHAT*FLOAT(I-1)
         ENDDO
         DO J=2,NYHAT+1
            YHAT(J)=YHAT(1)+DYHAT*FLOAT(J-1)
         ENDDO
         DT_GETTO=0.D+0
         DO I=1,NXHAT+1
            DO J=1,NYHAT+1
               U_HAT(I,J)=0.D+0
               V_HAT(I,J)=0.D+0
            ENDDO
         ENDDO
      ENDIF

c      stampa file 'retic'   
c      DO 10 I=1,NTETA+1
c         DO 20 J=1,NZETA
c           WRITE(35,22)real(z(i,j)),aimag(z(i,j))
c 20      CONTINUE
c         write(35,*)
c 10   CONTINUE
C      DO 21 J=1,NZETA
C            WRITE(35,22)real(z(1,j)),aimag(z(1,j))
C 21      CONTINUE
c 22   FORMAT(E13.7,1X,E13.7)

      RETURN
      END

C ***********************************************************************

      SUBROUTINE CONF(chix,zx)
      
      INCLUDE 'xcilgetto.com'
      
      COMPLEX(16) chix,zx
 
      zx=chix-e+1.d0/(chix-e)+dis

      RETURN
      END

C ***********************************************************************

      SUBROUTINE INIT
      INCLUDE 'xcilgetto.com'
      
      DT=2.d0*PI/float(NTciclo)
c      DT=1.d0/FLOAT(NTCICLO)

      IF (IOPT.EQ.0) THEN
         T=0.d0
         DO I=0,NTETA+1
            DO J=1,NZETA
               PSI(I,J)=0.d0
               OMG(I,J,1)=0.d0
               OMG(I,J,2)=0.d0
! 10      CONTINUE
            ENDDO
         ENDDO 
      ELSE
         CALL OPFILEB(20,'OLD','Nome file iniziale ? ')      
c        OPEN(UNIT=20,FILE='out',STATUS='NEW',err=100.,FORM='UNFORMATTED')
         READ(20)T
c         DO 11 I=1,NTETA+1
         DO 11 I=1,NTETA
            DO 12 J=1,NZETA
               READ(20)OMG(I,J,1),PSI(I,J)
               OMG(I,J,2)=OMG(I,J,1)
 12         CONTINUE
 11      CONTINUE
         CLOSE(20)
         do j=1,nzeta
            omg(0,j,1)=omg(nteta,j,1)
            omg(nteta+1,j,1)=omg(1,j,1)
            psi(0,j)=psi(nteta,j)
            psi(nteta+1,j)=psi(1,j)
         enddo
      ENDIF

      write(6,*)' t iniziale = ',t

      RETURN
      END

C ***********************************************************************

      SUBROUTINE VORTICITY
      
      INCLUDE 'xcilgetto.com'

      REAL(8) MATCONV,MATCONV1
      COMPLEX(16) XX,YY,PP
      DIMENSION XX(ITETA),YY(ITETA),PP(ITETA),MATCONV(ITETA,IZETA),
     $     MATCONV1(ITETA,IZETA)
      DIMENSION AM(ITETA+1),BM(ITETA+1),CM(ITETA+1),DM(ITETA+1),
     *           GAM(ITETA+1)
      DIMENSION AN(IZETA),BN(IZETA),CN(IZETA),DN(IZETA),
     *           WN(IZETA)

C
C    !!!!!!!!!!!!!!!! Soluzione ADI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Procedura di de-aliasing commentata
      
C     primo passo temporale implicito in teta
      
      DO J=2,NZETA-1
C         DO I=1,NTETA
C            XX(I)=CMPLX(V_R(I,J),0.)
C            YY(I)=CMPLX(0.5D0*(OMG(I,J+1,1)-OMG(I,J-1,1))/DZETA,0.)
C         ENDDO
C         CALL CFFT1(XX,NTETA,1)
C         CALL CFFT1(YY,NTETA,1)
C         CALL CONV(NTETA,XX,YY,PP)
C         CALL CFFT1(PP,NTETA,-1)
C
         DO I=1,NTETA
c         if (iusc1.eq.nusc1) 
c     $        write (70,*)teta(i),(V_R(I,J)*0.5D0*(OMG(I,J+1,1)-
c     $        OMG(I,J-1,1))/DZETA),REAL(PP(I))
         
            OX1=jacobsq(i,j)*R(J)*REYN*V_TETA(I,J)*DTETA*0.5D0
            OX2=jacob(i,j)*r2(j)*(reyn/kc)*2.D0*dteta*dteta/dt
            AM(I)=-(1.D0+ox1)
            BM(I)=2.D0+ox2
            CM(I)=-(1.D0-ox1)
            ox3=0.5*a*dzeta/r(j)
            ox4=jacobsq(i,j)*ez(j)*reyn*v_r(i,j)*dzeta*0.5
            ox5=-(dteta*r(j)/(dzeta*ez(j)))**2
            DM(I)=-omg(i,j-1,1)*ox5*(1.d0-ox3+ox4)
     *           +omg(i,j,1)*ox5*(2.d0-
     *           jacob(i,j)*ez2(j)*(reyn/kc)*2.d0*dzeta*dzeta/dt)
     *           -omg(i,j+1,1)*ox5*(1.d0+ox3-ox4) 
C     DM(I)=-omg(i,j-1,1)*ox5*(1.d0-ox3)
C     *            +omg(i,j,1)*ox5*(2.d0-
C     *               jacob(i,j)*ez2(j)*(reyn/kc)*2.d0*dzeta*dzeta/dt)
C     *            -omg(i,j+1,1)*ox5*(1.d0+ox3)
C     *            -jacobsq(i,j)*r2(j)*dteta*dteta*reyn*REAL(pp(i))/ez(j)
         ENDDO
      
         CALL TRIPER(AM,BM,CM,DM,NTETA-1,GAM)
      
         DO I=1,NTETA
            OMG(I,J,2)=DM(I)
         ENDDO
      ENDDO
      
      DO J=2,NZETA-1
         DO I=1,NTETA
            OMG(I,J,1)=omg(i,j,2)
         ENDDO
      ENDDO
     
      do j=2,nzeta-1
          OMG(0,J,1)=OMG(nteta,J,2)
          OMG(nteta+1,J,1)=OMG(1,J,2)
      enddo
      
C     secondo passo temporale implicito in zeta

C      DO J=2,NZETA-1
C         DO I=1,NTETA
C            XX(I)=CMPLX(V_TETA(I,J),0.)
C            YY(I)=CMPLX(0.5D0*(OMG(I+1,J,1)-OMG(I-1,J,1))/DTETA,0.)
C         ENDDO
C         CALL CFFT1(XX,NTETA,1)
C         CALL CFFT1(YY,NTETA,1)
C         CALL CONV(NTETA,XX,YY,PP)
C         CALL CFFT1(PP,NTETA,-1) 
C         DO I=1,NTETA
C            MATCONV(I,J)=REAL(PP(I))
C            if (iusc1.eq.nusc1) 
C     $          write (71,*)teta(i),(V_TETA(I,J)*0.5D0*(OMG(I+1,J,1)
C     $           -OMG(I-1,J,1))/DTETA),matconv(i,j)
C            
C         ENDDO
C      ENDDO
      DO I=1,NTETA
         DO J=2,NZETA-1 
            PX1=jacobsq(i,j)*Ez(j)*REYN*V_R(I,J)*DZETA*0.5D0
            PX2=0.5D0*A*DZETA/R(J)
            AN(J)=(1.D0-PX1+PX2)
            BN(J)=2.D0+jacob(i,j)*ez2(j)*(reyn/kc)*2.D0*DZETA*DZETA/DT
            CN(J)=(1.D0+PX1-PX2)
            px4=-(dzeta*ez(j)/(dteta*r(j)))**2
            px5=jacobsq(i,j)*r(j)*reyn*v_teta(i,j)*dteta*0.5d0
            dn(j)=-omg(i-1,j,1)*px4*(1.d0+px5)
     *           +omg(i,j,1)*px4*
     *           (2.d0-jacob(i,j)*r2(j)*(reyn/kc)*2.d0*dteta*dteta/dt)
     *           -omg(i+1,j,1)*px4*(1.d0-px5)
C            d1=-jacobsq(i,j)*ez2(j)*reyn*dzeta*dzeta*matconv(i,j)/r(j)
C            d2=-omg(i-1,j,1)*px4*(1.d0)
C     *            +omg(i,j,1)*px4*
C     *            (2.d0-jacob(i,j)*r2(j)*(reyn/kc)*2.d0*dteta*dteta/dt)
C     *            -omg(i+1,j,1)*px4*(1.d0)
C            dn(j)=-omg(i-1,j,1)*px4*(1.d0)
C     *            +omg(i,j,1)*px4*
C     *            (2.d0-jacob(i,j)*r2(j)*(reyn/kc)*2.d0*dteta*dteta/dt)
C     *            -omg(i+1,j,1)*px4*(1.d0)
C     *           -jacobsq(i,j)*ez2(j)*reyn*dzeta*dzeta*matconv(i,j)/r(j)  
C            if (iusc1.eq.nusc1) write(72,*)teta(i),zeta(j),d1,d2
         ENDDO  
       
c         CALL BOUNDOMG(A1,I)
	 CALL BOUNDOMG(OMG(I,1,2),I)
c         CALL TRIG(AN,BN,CN,DN,WN,NZETA,A1)
         CALL TRIG(AN,BN,CN,DN,WN,NZETA,OMG(I,1,2))
         DO J=1,NZETA
            OMG(I,J,2)=WN(J)
         ENDDO
      ENDDO

C
C      !!!!!!!!!!! Soluzione esplicita !!!!!!!!!!!!!!!!!!
C

C     Procedura di de-aliasing

c      DO J=2,NZETA-1
c         DO I=1,NTETA
c            XX(I)=CMPLX(V_R(I,J),0.d0)
c            YY(I)=CMPLX(0.5D0*(OMG(I,J+1,1)-OMG(I,J-1,1))/DZETA,0.)
c         ENDDO
c         CALL CFFT1(XX,NTETA,1)
c         CALL CFFT1(YY,NTETA,1)
c         CALL CONV(NTETA,XX,YY,PP)
c         CALL CFFT1(PP,NTETA,-1)
c         DO I=1,NTETA
c            MATCONV(I,J)=REAL(PP(I))
cc            if (iusc1.eq.nusc1) 
cc     $          write (72,*)teta(i),(V_R(I,J)*0.5D0*(OMG(I,J+1,1)
cc     $           -OMG(I,J-1,1))/DZETA),matconv(i,j)
c         ENDDO
c         DO I=1,NTETA
c            XX(I)=CMPLX(V_TETA(I,J),0.d0)
c            YY(I)=CMPLX(0.5D0*(OMG(I+1,J,1)-OMG(I-1,J,1))/DTETA,0.d0)
c         ENDDO
c         CALL CFFT1(XX,NTETA,1)
c         CALL CFFT1(YY,NTETA,1)
c         CALL CONV(NTETA,XX,YY,PP)
c         CALL CFFT1(PP,NTETA,-1) 
c         DO I=1,NTETA
c            MATCONV1(I,J)=REAL(PP(I))
cc            if (iusc1.eq.nusc1) 
cc     $          write (71,*)teta(i),(V_TETA(I,J)*0.5D0*(OMG(I+1,J,1)
cc     $           -OMG(I-1,J,1))/DTETA),matconv1(i,j)
c         ENDDO
c      ENDDO
C      DO I=1,NTETA
C         DO J=2,NZETA-1
C            COST=-(V_R(I,J)*(OMG(I,J+1,1)-OMG(I,J-1,1))/
C     *              (2.d0*DZETA*EZ(J))+
C     *             V_TETA(I,J)*(OMG(I+1,J,1)-OMG(I-1,J,1))/
C     *              (2.d0*DTETA*R(J)))*
C     *            KC/JACOBSQ(I,J)
C     *           +((OMG(I,J+1,1)-OMG(I,J-1,1))*A/
C     *              (2.d0*DZETA*EZ2(J)*R(J))+
C     *             (OMG(I,J+1,1)-2.d0*OMG(I,J,1)+OMG(I,J-1,1))/
C     *              (DZETA*DZETA*EZ2(J))+
C     *             (OMG(I+1,J,1)-2.d0*OMG(I,J,1)+OMG(I-1,J,1))/
C     *              (DTETA*DTETA*R2(J)))*
C     *            KC/(REYN*JACOB(I,J)) 
Cc            COST=-(MATCONV(I,J)/EZ(J)+
Cc     *             MATCONV1(I,J)/R(J))*KC/JACOBSQ(I,J)
Cc     *           +((OMG(I,J+1,1)-OMG(I,J-1,1))*A/
Cc     *              (2.d0*DZETA*EZ2(J)*R(J))+
Cc     *             (OMG(I,J+1,1)-2.d0*OMG(I,J,1)+OMG(I,J-1,1))/
Cc     *              (DZETA*DZETA*EZ2(J))+
cC     *             (OMG(I+1,J,1)-2.d0*OMG(I,J,1)+OMG(I-1,J,1))/
Cc     *              (DTETA*DTETA*R2(J)))*
Cc     *            KC/(REYN*JACOB(I,J))  
C           OMG(I,J,2)=OMG(I,J,1)+COST*DT
C        ENDDO
C      ENDDO
      
      do j=2,nzeta-1
          OMG(0,J,2)=OMG(nteta,J,2)
          OMG(nteta+1,J,2)=OMG(1,J,2)
      enddo

      DO I=0,NTETA+1
C         DO J=2,NZETA-1
         DO J=1,NZETA
            OMG(I,J,1)=OMG(I,J,2)
         ENDDO
      ENDDO

c      do i=1,nteta                      !?
c         do j=1,nzeta                   !?
c            omg(i,j,1)=omg(i,j,2)       !?
c         enddo                          !?
c      enddo                             !?
c      do J=1,nzeta                      !?
c         omg(0,j,1)=omg(nteta,j,1)      !?
c         omg(nteta+1,j,1)=omg(1,j,2)    !?
c      ENDDO                             !?

      RETURN
      END
  
C ***************************************************************************

      SUBROUTINE BOUNDOMG(W,I)
c
c     Soluzione esplicita
c      SUBROUTINE BOUNDOMG      

      INCLUDE 'xcilgetto.com'
C 
C     calcola  omega sulla parete all'i-esimo teta  NO
C     calcola omega alla parete e all'inf per ogni teta nell sol. esplicita
C

c      DO I=1,NTETA              

      V1=( -dl*SIN(ALPHA)+DALPHA*AIMAG(Z(I,1)) )/Kc
      V2=( -dl*COS(ALPHA)-DALPHA*REAL(Z(I,1)) )/Kc

      VA1=-dxdcsi(i,1)*sin(TETA(I))+DXDeta(i,1)*cos(TETA(I))
      VA2=-dydcsi(i,1)*sin(TETA(I))+dydeta(i,1)*cos(TETA(I))

      W=dalpha/kc
     *  -2.d0*(PSI(I,2)-PSI(I,1))/( jacob(i,1)*DZETA*DZETA*Ez2(1) )
     *  +( 2.d0/DZETA-A/R(1) )*(V1*VA1+V2*VA2)/( Jacob(I,1)*Ez(1) )
     *  +v1*(-sin(TETA(I))*dxdcsi(i,1)+cos(TETA(I))*dxdeta(i,1)
     *       +cos(TETA(I))*d2xdcsiteta(i)+sin(TETA(I))*d2xdetateta(i))
     *      /( R(1)*Jacob(i,1) )
     *  +v2*(-sin(TETA(I))*dydcsi(i,1)+cos(TETA(I))*dydeta(i,1)
     *       +cos(TETA(I))*d2ydcsiteta(i)+sin(TETA(I))*d2ydetateta(i))
     *      /( R(1)*Jacob(i,1) )  
      
c      OMG(I,1,1)=W                              
c      OMG(I,NZETA,1)=0.d0
c
c      ENDDO
c
c      OMG(0,1,1)=OMG(NTETA,1,1)
c      OMG(NTETA+1,1,1)=OMG(1,1,1)
c      OMG(0,NZETA,1)=OMG(NTETA,NZETA,1)
c      OMG(NTETA+1,NZETA,1)=OMG(1,NZETA,1)
   
      RETURN
      END

C ***********************************************************************

      SUBROUTINE POISSON
      
      INCLUDE 'xcilgetto.com'

      COMPLEX(16) A1,AM,DN,WN,MATPSI1,Y,DY1,DYM,MATPSIM
      
      DIMENSION AN(IZETA),BN(IZETA),CN(IZETA),DN(IZETA), WN(IZETA)
       
      DIMENSION MATPSI1(ITETA),Y(ITETA),DY1(ITETA)
      DIMENSION MATPSIM(ITETA),DYM(ITETA)

C !!!!!!!!!!!!!!!!!!!!!!!!!!!      
C      DO J=1,NZETA
      DO J=2,NZETA-1
         DO I=1,NTETA 
            Y(I)=CMPLX(jacob(i,j)*OMG(I,J,1),0.d0)           
         ENDDO
         CALL CFFT1(Y,NTETA,1)
         DO I=1,NTETA/2+1
            MATOMG(I,J)=Y(I)
         ENDDO
      ENDDO
      
      DO I=1,NTETA
         CALL BOUNDPSI(DY,DM,I)
         DY1(I)=CMPLX(DY,0.d0)
         DYM(I)=CMPLX(DM,0.d0)
      ENDDO
         
      CALL CFFT1(DY1,NTETA,1)
      CALL CFFT1(DYM,NTETA,1)

C     RISOLVO LE ARMONICHE NEGATIVE DA -1 A -NTETA/2
      MATPSI1(1)=C0
      MATPSIM(1)=C0
C !!!!!!!!!!!!!!!!!!!!!!!!!
C      DO I=2,NTETA/2+1
C !!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=2,NTETA/2        
         MATPSI1(I)=DY1(I)/CMPLX(0.d0,-FLOAT(I-1))
         MATPSIM(I)=DYM(I)/CMPLX(0.d0,-FLOAT(I-1))
      ENDDO

C !!!!!!!!!!!!!!!!!!!!!!!            
C      DO I=1,NTETA/2+1
C !!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,NTETA/2
         DO J=2,NZETA-1 
            AN(J)=1.d0+A*DZETA/(2.d0*R(J))
            BN(J)=2.d0+(EXP(ZETA(J))*DZETA*(-FLOAT(I-1))/R(J))**2
            CN(J)=1.d0-A*DZETA/(2.d0*R(J))
            DN(J)=MATOMG(I,J)*(DZETA*EXP(ZETA(J)))**2
         ENDDO
         A1=MATPSI1(I)
         AM=MATPSIM(I)
         CALL TRIGCMPLX(AN,BN,CN,DN,WN,NZETA,A1,AM)
         DO J=1,NZETA
            MATPSI(I,J)=WN(J)
         ENDDO
      ENDDO

      DO J=1,NZETA
         Y(1)=MATPSI(1,J)
         DO I=1,NTETA/2-1
            Y(I+1)=MATPSI(I+1,J)
            Y(NTETA+1-I)=CONJG(Y(I+1))
         ENDDO
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C         Y(NTETA/2+1)=MATPSI(NTETA/2+1,J)     
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         Y(NTETA/2+1)=C0     
         CALL CFFT1(Y,NTETA,-1)
         DO I=1,NTETA
            PSI(I,J)=REAL(Y(I))
            ENDDO
      ENDDO
    
      DO J=1,NZETA
         PSI(0,J)=PSI(NTETA,J)
         PSI(NTETA+1,J)=PSI(1,J)
      ENDDO
      
      RETURN
      END

C *******************************************************************
      
      SUBROUTINE BOUNDPSI(W,WW,I)
      
      INCLUDE 'xcilgetto.com'
 
      V1=(-DL*SIN(ALPHA)+DALPHA*AIMAG(Z(I,1)))/KC
      V2=(-DL*COS(ALPHA)-DALPHA*REAL(Z(I,1)))/KC

      VA1=DXDcsi(i,1)*COS(TETA(I))+DXDETA(i,1)*SIN(TETA(I))
      VA2=DYDcsi(i,1)*COS(TETA(I))+DYDeta(i,1)*SIN(TETA(I))
      W=-R(1)*(V1*VA1+V2*VA2)

      VA1=DXDcsi(I,NZETA)*COS(TETA(I))+DXDETA(I,NZETA)*SIN(TETA(I))
      VA2=DYDcsi(I,NZETA)*COS(TETA(I))+DYDeta(I,NZETA)*SIN(TETA(I))
      WW=R(NZETA)*U_TILDE*(VA1*COS(BETA-ALPHA)+VA2*SIN(BETA-ALPHA))

      RETURN
      END

C ***********************************************************************

      SUBROUTINE VEL

      INCLUDE 'xcilgetto.com'
 
C     calcolo delle componenti di velocita' al passo n
 
      DO I=1,NTETA
         DO J=2,NZETA-1
            V1=(-DL*SIN(ALPHA)+DALPHA*AIMAG(Z(I,J)))/KC
            V2=(-DL*COS(ALPHA)-DALPHA*REAL(Z(I,J)))/KC
            VA1=DXDcsi(i,j)*COS(TETA(I))+DXDeta(i,j)*SIN(TETA(I))
            VA2=DYDcsi(i,j)*COS(TETA(I))+DYDeta(i,j)*SIN(TETA(I))
            VA3=-DXDcsi(i,j)*SIN(TETA(I))+DXDeta(i,j)*COS(TETA(I))
            VA4=-DYDcsi(i,j)*SIN(TETA(I))+DYDeta(i,j)*COS(TETA(I))
            V_R(I,J)=( (PSI(I+1,J)-PSI(I-1,J))
     *                            /(2.d0*DTETA*R(J))
     *                                   +V1*VA1+V2*VA2 )/Jacobsq(i,j)
            V_TETA(I,J)=( -(PSI(I,J+1)-PSI(I,J-1))
     *                            /(2.d0*DZETA*ez(j))
     *                                   +V1*VA3+V2*VA4 )/Jacobsq(i,j)
         ENDDO
      ENDDO

      prod=kc*reyn*dt
      rapp=dt*kc/reyn
c      write(58,27)t,(v_r(1,200)*v_r(1,200)*2.d0*r(200)/
c     $     (dzeta+2.d0*r(200))
c     $     +v_teta(1,200)*v_teta(1,200))*prod,
c     $     (v_r(1,50)*v_r(1,50)*2.d0*r(50)/(dzeta+2.d0*r(50))
c     $     +v_teta(1,50)*v_teta(1,50))*prod,
c     $     (v_r(1,100)*v_r(1,100)*2.d0*r(100)/(dzeta+2.d0*r(100))
c     $     +v_teta(1,100)*v_teta(1,100))*prod,
c     $     (v_r(1,2)*v_r(1,2)*2.d0*r(2)/(dzeta+2.d0*r(2))
c     $     +v_teta(1,2)*v_teta(1,2))*prod
c      write(57,27)t,rapp*(1.d0/(ez2(2)*dzeta*dzeta)+1.d0/(2.d0*r(2)
c     $     *ez2(2)*dzeta)+1.d0/(r2(2)*dteta*dteta))/jacob(1,2),
c     $     rapp*(1.d0/(ez2(4)*dzeta*dzeta)+1.d0/(2.d0*r(4)
c     $     *ez2(4)*dzeta)+1.d0/(r2(4)*dteta*dteta))/jacob(1,4),
c     $     rapp*(1.d0/(ez2(7)*dzeta*dzeta)+1.d0/(2.d0*r(7)
c     $     *ez2(7)*dzeta)+1.d0/(r2(7)*dteta*dteta))/jacob(1,7),
c     $     rapp*(1.d0/(ez2(10)*dzeta*dzeta)+1.d0/(2.d0*r(10)
c     $     *ez2(10)*dzeta)+1.d0/(r2(10)*dteta*dteta))/jacob(1,10)
      write(59,27)t,omg(1,1,2),omg(1,2,2),omg(1,4,2),omg(1,10,2)
      write(58,27)t,omg(nteta/2+1,1,2),omg(nteta/2+1,2,2),
     $     omg(nteta/2+1,4,2),omg(nteta/2+1,10,2)
 27   FORMAT(5(E13.7,1X))

      DO I=1,NTETA
         V_R(I,1)=0.d0
         V_TETA(I,1)=0.d0
      ENDDO

      DO I=1,NTETA
         Ugr=U_TILDE*COS(BETA)*COS(ALPHA)+(U_TILDE*SIN(BETA)-DL)*
     $        SIN(ALPHA)+AIMAG(Z(I,NZETA))*DALPHA
         Vgr=-U_TILDE*COS(BETA)*SIN(ALPHA)+(U_TILDE*SIN(BETA)-DL)*
     $        COS(ALPHA)-REAL(Z(I,NZETA))*DALPHA
         V_CSI=(DXDCSI(I,NZETA)*Ugr+DYDCSI(I,NZETA)*Vgr)
     $        /jacobsq(i,nzeta)
         V_ETA=(DXDETA(I,NZETA)*Ugr+DYDETA(I,NZETA)*Vgr)
     $        /jacobsq(i,nzeta)
         V_R(I,NZETA)=V_CSI*COS(TETA(I))+V_ETA*SIN(TETA(I))
         V_TETA(I,NZETA)=-V_CSI*SIN(TETA(I))+V_ETA*COS(TETA(I))
      ENDDO
      DO J=1,NZETA
         V_R(NTETA+1,J)=V_R(1,J)
         V_TETA(NTETA+1,J)=V_TETA(1,J)     
      ENDDO

      RETURN
      END

C ***********************************************************************

      SUBROUTINE PUNTOV

      INCLUDE 'xcilgetto.com'

      COMPLEX(16) Zhat_P,Z_P,CHI_P
C     calcolo delle componenti di velocita' al punto (x,y)
      
      DO K=1,NUM_PUNTI

         ZHAT_P=VET_Z(K)
         Z_P=CMPLX(REAL(Zhat_P)*COS(ALPHA)+
     $        (AIMAG(Zhat_P)-L)*SIN(ALPHA),
     $        (AIMAG(Zhat_P)-L)*COS(ALPHA)-REAL(Zhat_P)*SIN(ALPHA))
         
         IF ((REAL(Z_P)-DIS).LT.0.) THEN
            CHI_P=.5d0*( (Z_P-DIS)-SQRT((Z_P-DIS)**2-4.d0) )+E
         ELSE
            CHI_P=.5d0*( (Z_P-DIS)+SQRT((Z_P-DIS)**2-4.d0) )+E
         ENDIF

         R_P=(REAL(CHI_P)*REAL(CHI_P)+AIMAG(CHI_P)*AIMAG(CHI_P))**0.5
         IF (REAL(CHI_P).GT.0.) THEN
            IF (AIMAG(CHI_P).GE.0.) 
     $           TETA_P=ATAN(AIMAG(CHI_P)/REAL(CHI_P))
            IF (AIMAG(CHI_P).LT.0.) TETA_P=ATAN(AIMAG(CHI_P)
     $           /REAL(CHI_P))+2.d0*PI
         ELSEIF (REAL(CHI_P).LT.0.) THEN
            TETA_P=ATAN(AIMAG(CHI_P)/REAL(CHI_P))+PI
         ELSEIF (REAL(CHI_P).EQ.0.) THEN
            IF (AIMAG(CHI_P).GE.0.) TETA_P=0.5d0*PI
            IF (AIMAG(CHI_P).LT.0.) TETA_P=1.5d0*PI
         ENDIF
         ZETA_P=DLOG(R_P+A)
         JINF=IDINT((ZETA_P-ZETA(1))/DZETA)+1
         JSUP=JINF+1
         IINF=IDINT(TETA_P/DTETA)+1
         ISUP=IINF+1
         
         DO I=IINF,ISUP
            DO J=JINF,JSUP

               v_csi=v_r(i,j)*cos(teta(i))-v_teta(i,j)*sin(teta(i))
               v_eta=v_r(i,j)*sin(teta(i))+v_teta(i,j)*cos(teta(i))
               IF (I.EQ.IINF) THEN
                  IF (J.EQ.JINF) THEN
                     U_GR(1,1)=(V_CSI*DYDETA(I,J)-V_ETA*DYDCSI(I,J))
     $                    /JACOBSQ(I,J)
                     V_GR(1,1)=(-V_CSI*DXDETA(I,J)+V_ETA*DXDCSI(I,J))
     $                    /JACOBSQ(I,J)
                  ELSE
                     U_GR(1,2)=(V_CSI*DYDETA(I,J)-V_ETA*DYDCSI(I,J))
     $                    /JACOBSQ(I,J)
                     V_GR(1,2)=(-V_CSI*DXDETA(I,J)+V_ETA*DXDCSI(I,J))
     $                    /JACOBSQ(I,J)
                  ENDIF
               ELSE 
                  IF (J.EQ.JINF) THEN
                     U_GR(2,1)=(V_CSI*DYDETA(I,J)-V_ETA*DYDCSI(I,J))
     $                    /JACOBSQ(I,J)
                     V_GR(2,1)=(-V_CSI*DXDETA(I,J)+V_ETA*DXDCSI(I,J))
     $                    /JACOBSQ(I,J)
                  ELSE
                     U_GR(2,2)=(V_CSI*DYDETA(I,J)-V_ETA*DYDCSI(I,J))
     $                    /JACOBSQ(I,J)
                     V_GR(2,2)=(-V_CSI*DXDETA(I,J)+V_ETA*DXDCSI(I,J))
     $                    /JACOBSQ(I,J)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         
         TT=(TETA_P-TETA(IINF))/(TETA(ISUP)-TETA(IINF))
         SS=(ZETA_P-ZETA(JINF))/(ZETA(JSUP)-ZETA(JINF))

         Ugr=(1.d0-TT)*(1.d0-SS)*U_GR(1,1)+
     $        TT*(1.d0-SS)*U_GR(2,1)+
     $        TT*SS*U_GR(2,2)+
     $        (1.d0-TT)*SS*U_GR(1,2)
         Vgr=(1.d0-TT)*(1.d0-SS)*V_GR(1,1)+
     $        TT*(1.d0-SS)*V_GR(2,1)+
     $        TT*SS*V_GR(2,2)+
     $        (1.d0-TT)*SS*V_GR(1,2)
c         Uhat=Ugr*COS(ALPHA)-Vgr*SIN(ALPHA)-DALPHA*
c     $        (AIMAG(Z_P)*COS(ALPHA)+REAL(Z_P)*SIN(ALPHA))
c         Vhat=DL+Ugr*SIN(ALPHA)+Vgr*COS(ALPHA)+DALPHA*
c     $        (-AIMAG(Z_P)*SIN(ALPHA)+REAL(Z_P)*COS(ALPHA))
         Uhat=Ugr*COS(ALPHA)-Vgr*SIN(ALPHA)-DALPHA*
     $        (AIMAG(ZHAT_P)-L)/KC
         Vhat=Ugr*SIN(ALPHA)+Vgr*COS(ALPHA)+(DL+DALPHA*
     $        REAL(ZHAT_P))/KC         

         WRITE(70+K,11)T,Uhat,Vhat,Ugr,Vgr,REAL(Z_P),AIMAG(Z_P)
      ENDDO
 11   FORMAT(7(E13.7,1X))

      RETURN
      END

C ***********************************************************************

      SUBROUTINE CMEDIO

      INCLUDE 'xcilgetto.com'
      
      COMPLEX(16) Z_P,CHI_P

C
C     calcolo della velocita' medio in una griglia zhat( , )
C      
     
      DT_GETTO=DT_GETTO+DT
      
      DO I=1,NXHAT+1
         DO J=1,NYHAT+1

            Z_P=CMPLX(XHAT(I)*COS(ALPHA)+(YHAT(J)-L)*SIN(ALPHA),
     $           (YHAT(J)-L)*COS(ALPHA)-XHAT(I)*SIN(ALPHA))
            IF ((REAL(Z_P)-DIS).LT.0.) THEN
               CHI_P=.5d0*( (Z_P-DIS)-SQRT((Z_P-DIS)**2-4.d0) )+E
            ELSE
               CHI_P=.5d0*( (Z_P-DIS)+SQRT((Z_P-DIS)**2-4.d0) )+E
            ENDIF

            R_P=(REAL(CHI_P)*REAL(CHI_P)
     $           +AIMAG(CHI_P)*AIMAG(CHI_P))**0.5
            IF (REAL(CHI_P).GT.0.) THEN
               IF (AIMAG(CHI_P).GE.0.) TETA_P=ATAN(AIMAG(CHI_P)
     $              /REAL(CHI_P))
               IF (AIMAG(CHI_P).LT.0.) TETA_P=ATAN(AIMAG(CHI_P)
     $              /REAL(CHI_P))+2.d0*PI
            ELSEIF (REAL(CHI_P).LT.0.) THEN
               TETA_P=ATAN(AIMAG(CHI_P)/REAL(CHI_P))+PI
            ELSEIF (REAL(CHI_P).EQ.0.) THEN
               IF (AIMAG(CHI_P).GE.0.) TETA_P=0.5d0*PI
               IF (AIMAG(CHI_P).LT.0.) TETA_P=1.5d0*PI
            ENDIF
            ZETA_P=DLOG(R_P+A)
            JINF=IDINT((ZETA_P-ZETA(1))/DZETA)+1
            JSUP=JINF+1
            IINF=IDINT(TETA_P/DTETA)+1
            ISUP=IINF+1

            DO II=IINF,ISUP
               DO JJ=JINF,JSUP

                  V_CSI=V_R(II,JJ)*COS(TETA(II))-
     $                 V_TETA(II,JJ)*SIN(TETA(II))
                  V_ETA=V_R(II,JJ)*SIN(TETA(II))+
     $                 V_TETA(II,JJ)*COS(TETA(II))
                  IF (II.EQ.IINF) THEN
                     IF (JJ.EQ.JINF) THEN
                        U_GR(1,1)=(V_CSI*DYDETA(II,JJ)-V_ETA*
     $                       DYDCSI(II,JJ))/JACOBSQ(II,JJ)
                        V_GR(1,1)=(-V_CSI*DXDETA(II,JJ)+V_ETA*
     $                       DXDCSI(II,JJ))/JACOBSQ(II,JJ)
                     ELSE
                        U_GR(1,2)=(V_CSI*DYDETA(II,JJ)-V_ETA*
     $                       DYDCSI(II,JJ))/JACOBSQ(II,JJ)
                        V_GR(1,2)=(-V_CSI*DXDETA(II,JJ)+V_ETA*
     $                       DXDCSI(II,JJ))/JACOBSQ(II,JJ)
                     ENDIF
                  ELSE 
                     IF (JJ.EQ.JINF) THEN
                        U_GR(2,1)=(V_CSI*DYDETA(II,JJ)-V_ETA*
     $                       DYDCSI(II,JJ))/JACOBSQ(II,JJ)
                        V_GR(2,1)=(-V_CSI*DXDETA(II,JJ)+V_ETA*
     $                       DXDCSI(II,JJ))/JACOBSQ(II,JJ)
                     ELSE
                        U_GR(2,2)=(V_CSI*DYDETA(II,JJ)-V_ETA*
     $                       DYDCSI(II,JJ))/JACOBSQ(II,JJ)
                        V_GR(2,2)=(-V_CSI*DXDETA(II,JJ)+V_ETA*
     $                       DXDCSI(II,JJ))/JACOBSQ(II,JJ)
                     ENDIF
                  ENDIF

               ENDDO
            ENDDO


            TT=(TETA_P-TETA(IINF))/(TETA(ISUP)-TETA(IINF))
            SS=(ZETA_P-ZETA(JINF))/(ZETA(JSUP)-ZETA(JINF))

            Ugr=(1.d0-TT)*(1.d0-SS)*U_GR(1,1)+
     $           TT*(1.d0-SS)*U_GR(2,1)+
     $           TT*SS*U_GR(2,2)+
     $           (1.d0-TT)*SS*U_GR(1,2)
            Vgr=(1.d0-TT)*(1.d0-SS)*V_GR(1,1)+
     $           TT*(1.d0-SS)*V_GR(2,1)+
     $           TT*SS*V_GR(2,2)+
     $           (1.d0-TT)*SS*V_GR(1,2)

            Uhat=Ugr*COS(ALPHA)-Vgr*SIN(ALPHA)-DALPHA*
     $           (YHAT(J)-L)/KC
            Vhat=Ugr*SIN(ALPHA)+Vgr*COS(ALPHA)+(DL+DALPHA*
     $           XHAT(I))/KC
            U_HAT(I,J)=U_HAT(I,J)+UHAT*DT
            V_HAT(I,J)=V_HAT(I,J)+VHAT*DT

         ENDDO
      ENDDO

      RETURN
      END

C ***********************************************************************

      SUBROUTINE OUTCMEDIO
      INCLUDE 'xcilgetto.com'

      OPEN(91,FILE='sciamedia')

      WRITE(91,*)' dt_getto = ',DT_GETTO
      DO I=1,NXHAT+1
         DO J=1,NYHAT+1
            U_HAT(I,J)=U_HAT(I,J)/DT_GETTO
            V_HAT(I,J)=V_HAT(I,J)/DT_GETTO
            WRITE(91,92)XHAT(I),YHAT(J),U_HAT(I,J),V_HAT(I,J)
         ENDDO
      ENDDO

      CLOSE(91)
 92   FORMAT(4(E13.7,1X))
      
      END
      
C ***********************************************************************

      SUBROUTINE FORCES
      INCLUDE 'xcilgetto.com'

      REAL(8) LAPL_U,LAPL_V,kcinv
      COMPLEX(16) Y,YY
      DIMENSION Y(ITETA),YY(ITETA)

c      write(90,*)' tempo = ',t
c      write(91,*)' tempo = ',t
c      write(92,*)' tempo = ',t

      kcinv=1.d+0/kc
      reyninv=1.d+0/reyn

      P0=0.d0
c      P(1)=P0

      PNEW(1)=P0
      PNEW(NTETA+1)=P0
      POD(NTETA/2+1)=P0
      DO I=1,NTETA
c
c     MODIFICATO !!!!!!
c
C         DJDCSI=2.*(DXDCSI(I,1)*d2xdcsi2(i)-dydcsi(i,1)*d2xdcsieta(I))
C         DJDETA=2.*(DXDCSI(I,1)*D2XDCSIETA(I)+DYDCSI(I,1)*D2XDCSI2(I))
         
C         LAPL=(2.*PSI(I,1)-5.*PSI(I,2)+4.*PSI(I,3)-PSI(I,4))
C     *                        /(EZ2(1)*DZETA*DZETA)+
C     *        (-3.*PSI(I,1)+4.*PSI(I,2)-PSI(I,3))*A
C     *                        /(EZ2(1)*R(1)*2.*DZETA)+
C     *        (PSI(I+1,1)-2.*PSI(I,1)+PSI(I-1,1))
C     *                        /(R(1)*R(1)*DTETA*DTETA)

C         DLAPLDZ=(-3.*PSI(I,5)+14.*PSI(I,4)-24.*PSI(I,3)+18.*PSI(I,2)-
C     *                        5.*PSI(I,1))/(2.*EZ2(1)*DZETA**3)+
C     *        (2.*PSI(I,1)-5.*PSI(I,2)+4.*PSI(I,3)-PSI(I,4))*
C     *                        (A-2.*R(1))/(R(1)*EZ2(1)*DZETA**2)-
C     *        (-3.*PSI(I,1)+4.*PSI(I,2)-PSI(I,3))*A*(2.*R(1)+EZ(1))
C     *                        /(EZ2(1)*R(1)*R(1)*2.*DZETA)+
C     *        (-3.*PSI(I+1,1)+4.*PSI(I+1,2)-PSI(I+1,3)+6.*PSI(I,1)-8.*
C     *     PSI(I,2)+2.*PSI(I,3)-3.*PSI(I-1,1)+4.*PSI(I-1,2)-PSI(I-1,3))
C     *                        /(2.*R(1)*R(1)*DZETA*DTETA**2)-
C     *        (PSI(I+1,1)-2.*PSI(I,1)+PSI(I-1,1))*2.*EZ(1)
C     *                        /(DTETA*DTETA*R(1)**3)
C         IF (I.EQ.1) THEN 
C         DLAPLDTETA=(2.*PSI(I+1,1)-2.*PSI(I-1,1)-5.*PSI(I+1,2)+5.*
C     *     PSI(I-1,2)+4.*PSI(I+1,3)-4.*PSI(I-1,3)-PSI(I+1,4)+PSI(I-1,4))
C     *                        /(2.*EZ2(1)*DTETA*DZETA**2)+
C     *        (-3.*PSI(I+1,1)+3.*PSI(I-1,1)+4.*PSI(I+1,2)-4.*PSI(I-1,2)-
C     *            PSI(I+1,3)+PSI(I-1,3))*A/(4.*EZ2(1)*R(1)*DZETA*DTETA)+
C     *        (PSI(I+2,1)-2.*PSI(I+1,1)+2.*PSI(I-1,1)-PSI(NTETA-1,1))
C     *                        /(2.*R(1)*R(1)*DTETA**3)
C         ELSEIF (I.EQ.NTETA) THEN
C         DLAPLDTETA=(2.*PSI(I+1,1)-2.*PSI(I-1,1)-5.*PSI(I+1,2)+5.*
C     *     PSI(I-1,2)+4.*PSI(I+1,3)-4.*PSI(I-1,3)-PSI(I+1,4)+PSI(I-1,4))
C     *                        /(2.*EZ2(1)*DTETA*DZETA**2)+
C     *        (-3.*PSI(I+1,1)+3.*PSI(I-1,1)+4.*PSI(I+1,2)-4.*PSI(I-1,2)-
C     *            PSI(I+1,3)+PSI(I-1,3))*A/(4.*EZ2(1)*R(1)*DZETA*DTETA)+
C     *        (PSI(2,1)-2.*PSI(I+1,1)+2.*PSI(I-1,1)-PSI(I-2,1))
C     *                        /(2.*R(1)*R(1)*DTETA**3)
C         ELSE
C         DLAPLDTETA=(2.*PSI(I+1,1)-2.*PSI(I-1,1)-5.*PSI(I+1,2)+5.*
C     *     PSI(I-1,2)+4.*PSI(I+1,3)-4.*PSI(I-1,3)-PSI(I+1,4)+PSI(I-1,4))
C     *                        /(2.*EZ2(1)*DTETA*DZETA**2)+
C     *        (-3.*PSI(I+1,1)+3.*PSI(I-1,1)+4.*PSI(I+1,2)-4.*PSI(I-1,2)-
C     *            PSI(I+1,3)+PSI(I-1,3))*A/(4.*EZ2(1)*R(1)*DZETA*DTETA)+
C     *        (PSI(I+2,1)-2.*PSI(I+1,1)+2.*PSI(I-1,1)-PSI(I-2,1))
C     *                        /(2.*R(1)*R(1)*DTETA**3)
C         ENDIF
C
C         DPDX(I)=(D2ALPHA*AIMAG(Z(I,1))+REAL(Z(I,1))*DALPHA**2-
C     *                            D2L*SIN(ALPHA))/KC+
C     *    ((DXDETA(I,1)*DJDCSI-DXDCSI(I,1)*DJDETA)*LAPL/JACOB(I,1)**3
C     *     +((-DXDETA(I,1)*COS(TETA(I))+DXDCSI(I,1)*SIN(TETA(I)))*
C     *                            DLAPLDZ/(R(1)+A)+
C     *    (DXDCSI(I,1)*COS(TETA(I))+DXDETA(I,1)*SIN(TETA(I)))*
C     *                        DLAPLDTETA/R(1))/JACOB(I,1)**2)*KC/REYN     
C
C         DPDY(I)=(-D2ALPHA*REAL(Z(I,1))+AIMAG(Z(I,1))*DALPHA**2-
C     *                            D2L*COS(ALPHA))/KC+
C     *    (-(-DYDETA(I,1)*DJDCSI+DYDCSI(I,1)*DJDETA)*LAPL/JACOB(I,1)**3
C     *     +((-DYDETA(I,1)*COS(TETA(I))+DYDCSI(I,1)*SIN(TETA(I)))*
C     *                            DLAPLDZ/(R(1)+A)+
C     *       (DYDETA(I,1)*SIN(TETA(I))+DYDCSI(I,1)*COS(TETA(I)))*
C     *                        DLAPLDTETA/R(1))/JACOB(I,1)**2)*KC/REYN      
         LAPL_U=(-(SIN(TETA(I))*DXDCSI(I,1)-COS(TETA(I))*DXDETA(I,1))*
     *            (-3.d0*OMG(I,1,1)+4.d0*OMG(I,2,1)-OMG(I,3,1))
     *                                    /(2.d0*DZETA*EZ(1))
     *           -(COS(TETA(I))*DXDCSI(I,1)+SIN(TETA(I))*DXDETA(I,1))*
     *            (OMG(I+1,1,1)-OMG(I-1,1,1))
     *                       /(2.d0*DTETA*R(1)))/JACOB(I,1)  
  
         LAPL_V=( (COS(TETA(I))*DYDETA(I,1)-SIN(TETA(I))*DYDCSI(I,1))*
     *            (-3.d0*OMG(I,1,1)+4.d0*OMG(I,2,1)-OMG(I,3,1))
     *                                    /(2.d0*DZETA*EZ(1))
     *           -(SIN(TETA(I))*DYDETA(I,1)+COS(TETA(I))*DYDCSI(I,1))*
     *            (OMG(I+1,1,1)-OMG(I-1,1,1))
     *                       /(2.d0*DTETA*R(1)))/JACOB(I,1)

         DPDX(I)=(D2ALPHA*AIMAG(Z(I,1))+REAL(Z(I,1))*DALPHA**2-
     *                            D2L*SIN(ALPHA))*KCinv*kcinv
         DPDY(I)=(-D2ALPHA*REAL(Z(I,1))+AIMAG(Z(I,1))*DALPHA**2-
     *                            D2L*COS(ALPHA))*KCinv*kcinv
         DPDX(I)=DPDX(I)+LAPL_U*REYNinv 
         DPDY(I)=DPDY(I)+LAPL_V*REYNinv 
c         write(90,22)teta(i),lapl_u,lapl_v,dpdx(i),dpdy(i)
 22   FORMAT(5(E13.7,1X))         
      ENDDO
      
      DPDX(NTETA+1)=DPDX(1)
      DPDY(NTETA+1)=DPDY(1)

C     Calcolo pressione integrando da 0 a 2*pi
C      DO I=1,NTETA
C         P(I+1)=P(I)+((DPDX(I)*DXDTETA(I)+DPDY(I)*DYDTETA(I))+
C     $        (DPDX(I+1)*DXDTETA(I+1)+DPDY(I+1)*DYDTETA(I+1)))
C     $        *DTETA*0.5D0
C      ENDDO

C
C     CALCOLO PRESSIONE COME U.B.MEHTA
C
      DO I=1,NTETA/2-1
         K=NTETA+2-I
         PNEW(I+1)=PNEW(I)+((DPDX(I)*DXDTETA(I)+DPDY(I)*DYDTETA(I))+
     $        (DPDX(I+1)*DXDTETA(I+1)+DPDY(I+1)*DYDTETA(I+1)))
     $        *DTETA*0.5D0
         PNEW(K-1)=PNEW(K)-((DPDX(K)*DXDTETA(K)+DPDY(K)*DYDTETA(K))+
     $        (DPDX(K-1)*DXDTETA(K-1)+DPDY(K-1)*DYDTETA(K-1)))
     $        *DTETA*0.5D0
      ENDDO
      P1=PNEW(NTETA/2)+((DPDX(NTETA/2)*DXDTETA(NTETA/2)+
     $     DPDY(NTETA/2)*DYDTETA(NTETA/2))+(DPDX(NTETA/2+1)*
     $     DXDTETA(NTETA/2+1)+DPDY(NTETA/2+1)*DYDTETA(NTETA/2+1)))
     $        *DTETA*0.5D0
      P2=PNEW(NTETA/2+2)-((DPDX(NTETA/2+2)*DXDTETA(NTETA/2+2)+
     $     DPDY(NTETA/2+2)*DYDTETA(NTETA/2+2))+(DPDX(NTETA/2+1)*
     $     DXDTETA(NTETA/2+1)+DPDY(NTETA/2+1)*DYDTETA(NTETA/2+1)))
     $        *DTETA*0.5D0
      PNEW(NTETA/2+1)=(P1+P2)*0.5D0     
C     Calcolo integrando da pi a 0 e da pi a 2*pi
      DO I=NTETA/2+1,3,-1
         POD(I-1)=POD(I)-((DPDX(I)*DXDTETA(I)+DPDY(I)*DYDTETA(I))+
     $        (DPDX(I-1)*DXDTETA(I-1)+DPDY(I-1)*DYDTETA(I-1)))
     $        *DTETA*0.5D0
      ENDDO
      DO I=NTETA/2+1,NTETA-1
         POD(I+1)=POD(I)+((DPDX(I)*DXDTETA(I)+DPDY(I)*DYDTETA(I))+
     $        (DPDX(I+1)*DXDTETA(I+1)+DPDY(I+1)*DYDTETA(I+1)))
     $        *DTETA*0.5D0
      ENDDO
      POD1=POD(2)-((DPDX(2)*DXDTETA(2)+DPDY(2)*DYDTETA(2))+
     $     (DPDX(1)*DXDTETA(1)+DPDY(1)*DYDTETA(1)))
     $     *DTETA*0.5D0
      POD2=POD(NTETA)+((DPDX(NTETA)*DXDTETA(NTETA)+
     $     DPDY(NTETA)*DYDTETA(NTETA))+
     $     (DPDX(NTETA+1)*DXDTETA(NTETA+1)+
     $     DPDY(NTETA+1)*DYDTETA(NTETA+1)))
     $     *DTETA*0.5D0
      POD(1)=(POD1+POD2)*0.5D0
      POD(NTETA+1)=POD(1)
C
C     CALCOLO DELLA DISTRIBUZIONE DI PRESSIONE TRAMITE FFT
C
      DO I=1,NTETA
         DPDTETA(i)=DPDX(I)*DXDTETA(I)+DPDY(I)*DYDTETA(I)
cc         WRITE(91,*)TETA(I),DPDTETA
         Y(I)=CMPLX(DPDTETA(i),0.d0)
      ENDDO
      CALL CFFT1(Y,NTETA,1)
      YY(1)=C0
      DO I=1,NTETA/2-1
         YY(I+1)=Y(I+1)/CMPLX(0.d0,-FLOAT(I))
         YY(NTETA+1-I)=CONJG(YY(I+1))
      ENDDO
      YY(NTETA/2+1)=C0
cc      WRITE(92,*)' ARMONICHE DPDTETA, P '
cc      DO I=1,NTETA
cc         WRITE(92,*)Y(I),YY(I)
cc      ENDDO
      CALL CFFT1(YY,NTETA,-1)
      DO I=1,NTETA
         Pff(I)=REAL(YY(I))
      ENDDO
      Pff(NTETA+1)=Pff(1)
cC
C     CALCOLO DELLE FORZE
C
      F_X=0.D0
      F_Xf=0.D0
      F_Xa=0.D0
      F_Y=0.D0
      F_Yf=0.D0
      F_Ya=0.D0
      M_Z=0.D0
      DO I=1,NTETA
         C2=COS(2.d0*TETA(I))
         S2=SIN(2.d0*TETA(I))
         CX1=DXDCSI(I,1)*DYDCSI(I,1)
         CX2=DYDETA(I,1)*DXDCSI(I,1)+DYDCSI(I,1)*DXDETA(I,1)
         DP1=((2.d0*PSI(I,1)-5.d0*PSI(I,2)+4.d0*PSI(I,3)-PSI(I,4))
     *               /(DZETA*DZETA)-
     *        (-3.d0*PSI(I,1)+4.d0*PSI(I,2)-PSI(I,3))
     *               /(2.d0*DZETA))/EZ2(1)
         DP2=(PSI(I+1,1)-2.d0*PSI(I,1)+PSI(I-1,1))/(DTETA*R(1))**2
         DP3=(PSI(I+1,1)-PSI(I-1,1))/(2.d0*DTETA*R(1)*R(1))
         DP4=(-3.d0*PSI(I+1,1)+3.d0*PSI(I-1,1)+4.d0*PSI(I+1,2)-
     $        4.d0*PSI(I-1,2)-PSI(I+1,3)+PSI(I-1,3))
     $        /(4.d0*R(1)*EZ(1)*DZETA*DTETA)
         DP5=(-3.d0*PSI(I,1)+4.d0*PSI(I,2)-PSI(I,3))
     $        /(2.d0*EZ(1)*R(1)*DZETA)

         D11=((C2*CX1+0.5d0*S2*CX2)*(DP1-DP2-DP5)+
     *        (2.d0*S2*CX1-C2*CX2)*(DP3-DP4)
     *                                    )/JACOB(I,1)**2
         D22=-D11
         D12=0.5d0*((-C2*CX2+2.d0*S2*CX1)*(DP1-DP2-DP5)-
     *            2.d0*(S2*CX2+2.d0*C2*CX1)*(DP3-DP4)
     *                                           )/JACOB(I,1)**2 
         
         F_Xf = F_Xf-PNEW(I)*DYDTETA(I)*dteta
         F_Xa = F_Xa+2.d0*(DYDTETA(I)*D11-DXDTETA(I)*D12)/REYN*DTETA         
         F_X = F_X+(-PNEW(I)*DYDTETA(I)+
     *               2.d0*(DYDTETA(I)*D11-DXDTETA(I)*D12)/REYN)*DTETA
         F_Yf = F_Yf+PNEW(I)*DXDTETA(I)*DTETA
         F_Ya = F_Ya+2.d0*(DYDTETA(I)*D12-DXDTETA(I)*D22)/REYN*DTETA
         F_Y = F_Y+(PNEW(I)*DXDTETA(I)+
     *               2.d0*(DYDTETA(I)*D12-DXDTETA(I)*D22)/REYN)*DTETA
         M_Z = M_Z+PNEW(I)*(DXDTETA(I)*REAL(Z(I,1))+
     *                 DYDTETA(I)*AIMAG(Z(I,1)))*DTETA+
     *           (2.D0*DTETA/REYN)*(-REAL(Z(I,1))*DXDTETA(I)*D22+
     *          (REAL(Z(I,1))*DYDTETA(I)+AIMAG(Z(I,1))*DXDTETA(I))*D12-
     *           AIMAG(Z(I,1))*DYDTETA(I)*D11)           
      ENDDO

      F_XHATf = F_Xf*COS(ALPHA)-F_Yf*SIN(ALPHA)
      F_YHATf = F_Xf*SIN(ALPHA)+F_Yf*COS(ALPHA)
      F_XHATa = F_Xa*COS(ALPHA)-F_Ya*SIN(ALPHA)
      F_YHATa = F_Xa*SIN(ALPHA)+F_Ya*COS(ALPHA)
      F_XHAT = F_X*COS(ALPHA)-F_Y*SIN(ALPHA)
      F_YHAT = F_X*SIN(ALPHA)+F_Y*COS(ALPHA)   

C     CALCOLO DELLE COMPONENTI DI VELOCITÀ NEL PIANO (X,Y)

C      DO I=1,NTETA
C         DO J=1,NZETA
C            v_csi=v_r(i,j)*cos(teta(i))-v_teta(i,j)*sin(teta(i))
C            v_eta=v_r(i,j)*sin(teta(i))+v_teta(i,j)*cos(teta(i))
C            U_GR(I,J)=(V_CSI*DYDETA(I,J)-V_ETA*DYDCSI(I,J))/JACOBSQ(I,J)
C            V_GR(I,J)=(-V_CSI*DXDETA(I,J)+V_ETA*DXDCSI(I,J))
C     *                             /JACOBSQ(I,J)
C         ENDDO
C      ENDDO
      
c      IF(IUSC1.EQ.NUSC1)THEN
C26genn         WRITE(26,24)T,F_XHAT,F_YHAT,F_X,F_Y,M_Z
      WRITE(26,25)T,F_XHAT,F_YHAT,M_Z,F_XHATf,F_YHATf,F_XHATa,F_YHATa
C         IUSC1=0
c      ENDIF
 24   FORMAT(6(E13.7,1X)) 
 25   FORMAT(8(E13.7,1X)) 
      
      RETURN
      END

C ***********************************************************************

      SUBROUTINE EFFIC
      INCLUDE 'xcilgetto.com'

      REAL(8) M_ZOLD,Mm

      ntot=25000      
C      f=0.D0
      fm_x = 0.D0
      fm_xf = 0.D0
      fm_xa = 0.D0
      fm_y = 0.D0
      fm_yf = 0.D0
      fm_ya = 0.D0
      Mm=0.D0
      delta_t=0.D0
      POT=0.D0
      POTY=0.D0
      POTM=0.D0

C26genn        10   read(26,24)T_act,F_XHAT,F_YHAT,F_X,F_Y,M_Z
 10   read(26,25)T_act,F_XHAT,F_YHAT,M_Z,F_XHATf,F_YHATf,F_XHATa,F_YHATa
      if (t_act.lt.pi*2.d+0) goto 10

      t_old=t_act
      f_xold=f_xhat
      f_yold=f_yhat
      f_xoldf=f_xhatf
      f_yoldf=f_yhatf
      f_xolda=f_xhata
      f_yolda=f_yhata
      m_zold=m_z
      dl_old=dlfun(t_old,Kc)
      dalpha_old=dalphafun(t_old,alphaamp,alphaphi)
      
      do i=1,ntot         
C26genn          read(26,24,end=20,err=20)t_act,F_Xhat,F_YHAT,F_X,F_Y,M_Z
         read(26,25,end=20,err=20)T_act,F_XHAT,F_YHAT,M_Z,F_XHATf,
     $        F_YHATf,F_XHATa,F_YHATa
         alpha=alphafun(t_act,alphaamp,alphaphi)
         
         dtloc=t_act-t_old
         dl=dlfun(t_act,Kc)
         dalpha=dalphafun(t_act,alphaamp,alphaphi)
c         f=f+0.5d0*(f_xold+f_xhat)*dtloc
         fm_x  = fm_x+0.5d0*(f_xold+f_xhat)*dtloc
         fm_xf = fm_xf+0.5d0*(f_xoldf+f_xhatf)*dtloc
         fm_xa = fm_xa+0.5d0*(f_xolda+f_xhata)*dtloc
         fm_y  = fm_y+0.5d0*(f_yold+f_yhat)*dtloc
         fm_yf = fm_yf+0.5d0*(f_yoldf+f_yhatf)*dtloc
         fm_ya = fm_ya+0.5d0*(f_xolda+f_xhata)*dtloc
         Mm = Mm+0.5d0*(M_zold+M_z)*dtloc
         delta_t = delta_t+dtloc

         POTY=POTY+0.5d0*(f_yold*dl_old+f_yhat*dl)*dtloc
         pOTM=pOTM+0.5d0*(m_zold*dalpha_old+m_z*dalpha)*dtloc
         POT=POTY+POTM
         
         t_old=t_act
         f_xold=f_xHAT
         f_yold=f_yHAT
         f_xoldf=f_xHATf
         f_yoldf=f_yHATf
         f_xolda=f_xHATa
         f_yolda=f_yHATa
         m_zold=m_z
         DL_OLD=DL
         dalpha_old=dalpha
         
      enddo
 20   continue
      
c      f=f/delta_t
      fm_x = fm_x/delta_t
      fm_xf = fm_xf/delta_t
      fm_xa = fm_xa/delta_t
      fm_y = fm_y/delta_t
      fm_yf = fm_yf/delta_t
      fm_ya = fm_ya/delta_t
      Mm = Mm/delta_t
C      ct=0.5d+0*f/(u_tilde*u_tilde)
      pot=pot/delta_t
      potY=potY/delta_t
      potM=potM/delta_t
      cT = 0.5d+0*fm_x/(u_tilde*u_tilde)
      cTf = 0.5d+0*fm_xf/(u_tilde*u_tilde)
      cTa = 0.5d+0*fm_xa/(u_tilde*u_tilde)
      cL = 0.5d+0*fm_y/(u_tilde*u_tilde)
      cLf = 0.5d+0*fm_yf/(u_tilde*u_tilde)
      cLa = 0.5d+0*fm_ya/(u_tilde*u_tilde)
      cM = Mm/(u_tilde*u_tilde*8.d0)
      cP = 0.5d+0*Pot/(u_tilde*u_tilde*u_tilde*kc)
      eff = fm_x*u_tilde*kc/pot

      WRITE(27,*)' DELTA_T = ',DELTA_T
      write(27,*)' fm_x, fm_xf, fm_xa ',fm_x,fm_xf,fm_xa
      write(27,*)' cT, cTf, cTa ',cT,cTf,cTa
      write(27,*)' fm_y, fm_yf, fm_ya ',fm_y,fm_yf,fm_ya
      write(27,*)' cL, cLf, cLa ',cL,cLf,cLa
      WRITE(27,*)' Mm, cM', Mm, cM
      WRITE(27,*)' pot, potY, potM',pot, potY, potM
      WRITE(27,*)' cP,  eff',  cP,  eff
c      WRITE(27,*)' F = ',f   
c      WRITE(27,*)' ct = ',ct 
c      WRITE(27,*)' PY = ',POTY,' PM = ',POTM,' p/delta_t = ',POT
c      WRITE(27,*)' EFF = ',EFF
c      WRITE(27,*)' FX  Ct  EFF '
c      WRITE(27,*)F,CT,EFF
 23   FORMAT(4(E13.7,1X))
 24   FORMAT(6(E13.7,1X))
 25   FORMAT(8(E13.7,1X))
      
      RETURN
      END

C ***********************************************************************

      SUBROUTINE CONV(NT,X,Y,P)

      IMPLICIT REAL(8) (A-H,O-Z)

      PARAMETER (ITETA=1024)
      COMPLEX(16) X,Y,P
      DIMENSION X(ITETA),Y(ITETA),P(ITETA)

      P(1)=X(1)*Y(1)
      DO I=1,NT/2-1
         P(1)=P(1)+X(nt+1-i)*Y(i+1)+Y(nt+1-i)*X(i+1)
      enddo

      DO N=1,NT/2-1
         P(nt+1-n)=X(1)*Y(nt+1-n)
         DO I=1,NT/2-1
            J=N-i
            if (j.gt.0) then 
               p(nT+1-n)=p(nt+1-n)+X(nt+1-i)*Y(nt+1-j)
            else
               p(nT+1-n)=p(nt+1-n)+x(nt+1-i)*Y(1-j)
            endif

            j=n+i
            if (j.lt.nt/2.and.j.gt.-nt/2) 
     $           p(nt+1-n)=p(nt+1-n)+x(1+i)*y(nt+1-j)
         enddo
         p(n+1)=conjg(p(nt+1-n))     
      enddo
      p(nt/2+1)=CMPLX(0.,0.)

      RETURN
      END

C ***********************************************************************

      real(8) function lfun(t,amp)

      real(8) t,amp

      lfun=amp*dsin(t)
c    LFUN=0.

      return
      end

C ***********************************************************************

      real(8) function dlfun(t,amp)

      real(8) t,amp

      dlfun=amp*dcos(t)
c      DLFUN=0.

      return
      end


C ***********************************************************************

      real(8) function d2lfun(t,amp)

      real(8) t,amp

      d2lfun=-amp*dsin(t)
c      D2LFUN=0.

      return
      end

C ***********************************************************************

      real(8) function alphafun(t,amp,fi)

      real(8) t,amp,fi

      alphafun=amp*dsin(t+fi)
c      ALPHAFUN=0.

      return
      end

C ***********************************************************************

      real(8) function dalphafun(t,amp,fi)

      real(8) t,amp,fi

      dalphafun=amp*dcos(t+fi)
c      DALPHAFUN=0.

      return
      end

C ***********************************************************************

      real(8) function d2alphafun(t,amp,fi)

      real(8) t,amp,fi

      d2alphafun=-amp*dsin(t+fi)
c      D2ALPHAFUN=0.

      return
      end

C **********************************************************************  
      
      SUBROUTINE TRIPER(AA,BB,CC,FF,N,GAM2)
c      INCLUDE 'xcilgetto.com'
      IMPLICIT REAL(8) (A-H,O-Z)
C
C ***************************************************************************
C     vedi Hisch I, pag.507
C     soluzione di un sistema tridiagonale di equazioni periodiche in N1 e N+2,
C     (i=1 e i=NTETA+1) 
C     
C     AA(K)*X(K-1)+BB(K)*X(K)+CC(K)*X(K+1)=FF(K)      K=N1,....,N+1
C
C     l'elemento nell'angolo superiore destro e' in A(N1)
C     l'elemento nell'angolo inferiore sinistro e' C(N+1)
C     tutti i vettori hanno dimensione N+2
C     GAM2 e' un vettore ausiliario
C     la soluzione e' messa in FF
C ***************************************************************************
C
      PARAMETER (IZETA=600)

      DIMENSION AA(IZETA),BB(IZETA),CC(IZETA),FF(IZETA),GAM2(IZETA)

      N1=1

      BB(N1)=1.D0/BB(N1)
      GAM2(N1)=-AA(N1)*BB(N1)
      AA(N1)=FF(N1)*BB(N1)
      N2=N1+1
      N1N=N1+N
      DO 10 K=N2,N
         K1=K-1
         CC(K1)=CC(K1)*BB(K1)
         BB(K)=BB(K)-AA(K)*CC(K1)
         BB(K)=1.D0/BB(K)
         GAM2(K)=-AA(K)*GAM2(K1)*BB(K)
         AA(K)=(FF(K)-AA(K)*AA(K1))*BB(K)
 10   CONTINUE
      GAM2(N)=GAM2(N)-CC(N)*BB(N)
C
C     BACK SUBSTITUTION
C
      FF(N)=AA(N)
      BB(N)=GAM2(N)
      DO 20 K1=N2,N
         K=N1N-K1
         K2=K+1
         FF(K)=AA(K)-CC(K)*FF(K2)
         BB(K)=GAM2(K)-CC(K)*BB(K2)
 20   CONTINUE

      K1=N+1
      ZAA=FF(K1)-CC(K1)*FF(N1)-AA(K1)*FF(N)
      ZAA=ZAA/(BB(K1)+AA(K1)*BB(N)+CC(K1)*BB(N1))
      FF(K1)=ZAA
      DO 30 K=N1,N
         FF(K)=FF(K)+BB(K)*ZAA
 30   CONTINUE

      FF(N+2)=FF(N1)
      RETURN
      END

C **************************************************************************

      SUBROUTINE TRIG(A,B,C,D,W,N,A1)
c      INCLUDE 'xcilgetto.com'      
      IMPLICIT REAL(8) (A-H,O-Z)
C
C ****************************************************************************
C     VEDI ROACHE APP.A
C     risolve il sistema tridiagonale per condizione di Dirichlet alla parete, 
C     dove e' nota la omega W(1) e condizione di Neumann all'infinito
C     
C      -A(K)*X(K+1)+B(K)*X(K)-C(K)*X(K-1)=D(K)
C ****************************************************************************
C
      PARAMETER (IZETA=600)

      DIMENSION A(IZETA),B(IZETA),C(IZETA),D(IZETA),W(IZETA)
      DIMENSION E(IZETA),F(IZETA)

      E(1)=0.D0
      F(1)=A1
      W(1)=F(1)
      DO M=2,N-1
         DEN=B(M)-C(M)*E(M-1)
         E(M)=A(M)/DEN
         F(M)=(D(M)+C(M)*F(M-1))/DEN
      ENDDO
      W(N)=F(N-1)/(1.d0-E(N-1))
      DO M=N-1,1,-1
         W(M)=E(M)*W(M+1)+F(M)
      ENDDO

      RETURN
      END

C *********************************************************************
C *********************************************************************

      SUBROUTINE TRIGCMPLX(A,B,C,D,W,N,A1,AM)
c      INCLUDE 'xcilgetto.com'
C ****************************************************************************
C     VEDI ROACHE APP.A
C     risolve il sistema tridiagonale per condizione di Dirichlet 
C     alla parete e all'infinito
C ****************************************************************************
      
      PARAMETER (IZETA=600)
      
      IMPLICIT REAL(8) (A-H,O-Z)
 
      COMPLEX(16) W,D,A1,AM,E,F,DEN
      
      DIMENSION A(1),B(1),C(1),D(1),W(1),E(IZETA),F(IZETA)

      E(1)=(0.,0.)
      F(1)=A1
      DO M=2,N-1
         DEN=B(M)-C(M)*E(M-1)
         E(M)=A(M)/DEN
         F(M)=(D(M)+C(M)*F(M-1))/DEN
      ENDDO

C     cond. al contorno Neumann per rmax
C     W(N)=(F(N-1)+AM)/(1.-E(N-1))
C     cond. al contorno di Dirichlet
      W(N)=AM
      DO M=N-1,1,-1
         W(M)=E(M)*W(M+1)+F(M)
      ENDDO

      RETURN
      END

C **************************************************************************

      SUBROUTINE OUTPUTINT
      
      INCLUDE 'xcilgetto.com'

      CHARACTER(14) NOMEFILE 
      CHARACTER(4) ITIME

      WRITE(ITIME,'(I4.4)') NINT(VET_PO(NIND_PO)*100)
      NOMEFILE='psiomg'//itime//'.dat'//char(0)
      OPEN(61,FILE=NOMEFILE)
      IF (ISTAMPAVEL.EQ.1) THEN 
         NOMEFILE='veloci'//itime//'.dat'//char(0)
         OPEN(62,FILE=NOMEFILE)
         WRITE(62,*)' tempo = ',t
      ENDIF
      WRITE(61,*)' tempo = ',t
      DO 10 I=1,NTETA+1,NUSC_TETA
         DO 20 J=1,NZETA,NUSC_ZETA
            COORD(I,J)=CMPLX(REAL(Z(I,J))*COS(ALPHA)-
     $           AIMAG(Z(I,J))*SIN(ALPHA),L+REAL(Z(I,J))*SIN(ALPHA)+
     $           AIMAG(Z(I,J))*COS(ALPHA))
            WRITE(61,21)REAL(COORD(I,J)),AIMAG(COORD(I,J)),
     $           OMG(I,J,1),PSI(I,J)
            IF (ISTAMPAVEL.EQ.1) THEN
               v_csi=v_r(i,j)*cos(teta(i))-v_teta(i,j)*sin(teta(i))
               v_eta=v_r(i,j)*sin(teta(i))+v_teta(i,j)*cos(teta(i))
               UGR=(V_CSI*DYDETA(I,J)-V_ETA*DYDCSI(I,J))/JACOBSQ(I,J)
               VGR=(-V_CSI*DXDETA(I,J)+V_ETA*DXDCSI(I,J))
     *              /JACOBSQ(I,J)
               
               UHAT=UGR*COS(ALPHA)-VGR*SIN(ALPHA)-DALPHA*
     $            (AIMAG(Z(I,J))*COS(ALPHA)+REAL(Z(I,J))*SIN(ALPHA))/kc
               VHAT=DL/kc+UGR*SIN(ALPHA)+VGR*COS(ALPHA)+DALPHA*
     $          (-AIMAG(Z(I,J))*SIN(ALPHA)+REAL(Z(I,J))*COS(ALPHA))/kc

               WRITE(62,21)REAL(COORD(I,J)),AIMAG(COORD(I,J)),
     $              UHAT,VHAT
            ENDIF
 20      CONTINUE

c	 write(61,*)
 10   CONTINUE
c      write(81,*)
      close(61)
      if  (ISTAMPAVEL.EQ.1) close(62)

 21   FORMAT(4(E13.7,1X))

      RETURN
      END

C *************************************************************************

      SUBROUTINE OUTPUT
      
      INCLUDE 'xcilgetto.com'

      OPEN(UNIT=25,FILE='out',FORM='UNFORMATTED')
c      OPEN(25,FILE=NAMEBIN,FORM='UNFORMATTED')
      write(25)t

c      DO 10 I=1,NTETA+1
      DO 10 I=1,NTETA
         DO 20 J=1,NZETA
            WRITE(25)OMG(I,J,1),PSI(I,J)
 20      CONTINUE
 10   CONTINUE
      
      CLOSE(25)

      RETURN
      END

C ***************************************************************************

      SUBROUTINE OUTPUTP

      INCLUDE 'xcilgetto.com'

      CHARACTER(17) NOMEFILE 
      CHARACTER(4) ITIME

      WRITE(ITIME,'(I4.4)') NINT(VET_PRES(NIND_PRES)*100)
      NOMEFILE='pressione'//itime//'.dat'//char(0)
      OPEN(45,FILE=NOMEFILE)
      NOMEFILE='vorticita'//itime//'.dat'//char(0)
      OPEN(46,FILE=NOMEFILE)
      write(45,*)' tempo = ',t  
      write(46,*)' tempo = ',t  
  
      do I=1,nteta+1
         the=DTETA*FLOAT(I-1)
         f=(DPDX(I)*DXDTETA(I)+DPDY(I)*DYDTETA(I))
c         WRITE(45,22)the,real(z(i,1)),p(i),DPDX(I),DPDY(I)
         WRITE(45,22)the,real(z(i,1)),pnew(i),POD(I),pff(i),dpdteta(i)
         WRITE(46,23)THE,real(z(i,1)),OMG(I,1,1),OMG(I,2,1),OMG(I,4,1)
      enddo
c      write(45,*)    
c      write(46,*)   
      close(45)
      close(46)

 22   FORMAT(6(E13.7,1X))
 23   FORMAT(5(E13.7,1X))
   
      RETURN
      END
  
C ***************************************************************************

      SUBROUTINE OUTPUTstr

      INCLUDE 'xcilgetto.com'
     
      write(66,*)' teta = ',teta(nteta/12+1) 
      write(67,*)' teta = ',teta(nteta/6+1) 
      write(68,*)' teta = ',teta(nteta/4+1) 
      write(69,*)' teta = ',teta(nteta/3+1) 
      write(70,*)' teta = ',teta(nteta*5/12+1) 

      write(66,*)' tempo = ',t 
      write(67,*)' tempo = ',t 
      write(68,*)' tempo = ',t 
      write(69,*)' tempo = ',t 
      write(70,*)' tempo = ',t 

      DO J=1,NZETA/2
         WRITE(66,*)R(J),REAL(Z(NTETA/12+1,J)),AIMAG(Z(NTETA/12+1,J)),
     *                  V_TETA(NTETA/12+1,J)
         WRITE(67,*)R(J),REAL(Z(NTETA/6+1,J)),AIMAG(Z(NTETA/6+1,J)),
     *                  V_TETA(NTETA/6+1,J)
         WRITE(68,*)R(J),REAL(Z(NTETA/4+1,J)),AIMAG(Z(NTETA/4+1,J)),
     *                  V_TETA(NTETA/4+1,J)
         WRITE(69,*)R(J),REAL(Z(NTETA/3+1,J)),AIMAG(Z(NTETA/3+1,J)),
     *                  V_TETA(NTETA/3+1,J)
         WRITE(70,*)R(J),REAL(Z(NTETA*5/12+1,J)),
     *              AIMAG(Z(NTETA*5/12+1,J)),V_TETA(NTETA*5/12+1,J)
      ENDDO   
      WRITE(66,*)
      WRITE(67,*)
      WRITE(68,*)
      WRITE(69,*)
      WRITE(70,*)

      RETURN
      END
    
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  
      SUBROUTINE cfft1(y,nx,iopt)

      IMPLICIT REAL(8) (A-H,O-Z)

c     complex vaux

      COMPLEX(16) Y
      DIMENSION Y(1)

      riopt=iopt

      N2X=anint(ALOG(FLOAT(NX))/ALOG(2.))

      IF (RIOPT.EQ.-1) Y(NX/2+1)=(0.,0.)

      CALL FFT(N2X,Y,riopt)

      IF (RIOPT.EQ.1) Y(NX/2+1)=(0.,0.)

      RETURN
      END
  
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  
      SUBROUTINE FFT(N,X,SIGN)
      IMPLICIT REAL(8) (A-H,O-Z)      
      
      DIMENSION M(1024),X(1)
      
      COMPLEX(16) X,WK,HOLD,Q

c      IF(N.GT.8)WRITE(6,*)' ATT. N NELLA FFT GT DEL TOLLERATO ! '
      
      LX=2**N
      DO 1 I=1,N
          M(I)=2**(N-I)
1     CONTINUE
      DO 5 L=1,N
          NBLOCK=2**(L-1)
          LBLOCK=LX/NBLOCK
          LBHALF=LBLOCK/2
          K=0
          DO 5 IBLOCK=1,NBLOCK
              FK=FLOAT(K)
              FLX=FLOAT(LX)
c              V=SIGN*6.2831853*FK/FLX
              V=SIGN*6.2831853d0*FK/FLX
              WK=CMPLX(COS(V),SIN(V))
              ISTART=LBLOCK*(IBLOCK-1)
              DO 2 I=1,LBHALF
                  J=ISTART+I
                  JH=J+LBHALF 
                  Q=X(JH)*WK
                  X(JH)=X(J)-Q
                  X(J)=X(J)+Q 
2             CONTINUE
              DO 3 I=2,N
                  II=I
                  IF(K.LT.M(I)) GO TO 4 
                  K=K-M(I)
3             CONTINUE
4             K=K+M(II) 
5     CONTINUE
      K=0 
      DO 8 J=1,LX
      IF(K.GE.J) THEN
          HOLD=X(J)
          X(J)=X(K+1)
          X(K+1)=HOLD
      ENDIF
      DO 6 I=1,N
      II=I
      IF(K.LT.M(I)) GO TO 7
      K=K-M(I)
6     CONTINUE
7     K=K+M(II)
8     CONTINUE
      IF(SIGN.LT.0.) RETURN
      DO 9 I=1,LX
          X(I)=X(I)/FLX
9     CONTINUE
      
      RETURN
      END

C *********************************************************************

      SUBROUTINE OPFILE (LU,STATUS,MESS)
C
C     LU      logic unit da assegnare al file
C     STATUS  stato richeisto 'old' o 'new'
C     MESS    stringa che contiene messaggio per la richiesta del nome
C
      
      CHARACTER*(*) STATUS,MESS
      CHARACTER(64) NAME

      IOS=-1
      DO WHILE (IOS.NE.0)
         WRITE(6,*) MESS
         READ(5,'(A)') NAME
c         NAME='st425.dat'
         IF (NAME.EQ.'/E') STOP
         OPEN(UNIT=LU,FILE=NAME,STATUS=STATUS,IOSTAT=IOS,
     *                           FORM='FORMATTED')
         IF (IOS.NE.0) THEN
            WRITE(6,*)'Errore ',IOS,' File : ',NAME
         ENDIF
      ENDDO

      END

C ****************************************************************************

      SUBROUTINE OPFILEB (LU,STATUS,MESS)
C
C     LU      logic unit da assegnare al file
C     STATUS  stato richeisto 'old' o 'new'
C     MESS    stringa che contiene messaggio per la richiesta del nome
C
      
      CHARACTER*(*) STATUS,MESS
      CHARACTER*64 NAME

      IOS=-1
      DO WHILE (IOS.NE.0)
         WRITE(6,*) MESS
         READ(5,'(A)') NAME
c         NAME='out'
         IF (NAME.EQ.'/E') STOP
         OPEN(UNIT=LU,FILE=NAME,STATUS=STATUS,IOSTAT=IOS,
     *                           FORM='UNFORMATTED')
         IF (IOS.NE.0) THEN
            WRITE(6,*)'Errore ',IOS,' File : ',NAME
         ENDIF
      ENDDO

      END

C ****************************************************************************

      SUBROUTINE OPFILEPVEL(I,ZZ)

      CHARACTER(13) NOMEFILE 
      CHARACTER(4) ITIME
      COMPLEX(16) ZZ

      NU=70+I
      WRITE(ITIME,'(I4.4)') NINT(ABS(AIMAG(ZZ)*10))
      IF (AIMAG(ZZ).LT.0.D+0) THEN
         NOMEFILE='pvel-'//itime//'.dat'//char(0)
      ELSE
         NOMEFILE='pvel+'//itime//'.dat'//char(0)
      ENDIF 
      
      OPEN(UNIT=NU,FILE=NOMEFILE)

      END

C ****************************************************************************
C ****************************************************************************
