      PROGRAM JOUKOWSKI

      INCLUDE 'sensvar.com'

     

      OPEN(1,FILE='ina.dat')
      OPEN(3,FILE='profilo')
      OPEN(4,FILE='prodistr_punti')
      open(9,file='foresults')

     
      call cost
      CALL INPUT
      CALL RETIC
 
  

      write(6,*)' dis = ',DIS
      write(6,*)' e = ',E
      write(6,*)' small = ',small
      write(6,*)' r_max = ',Rmax, ' r_int = ',Rint
      write(6,*)' a = ',a
      write(6,*)' Nmodes = ', Nmodes
      write(6,*)' NTETA,NZETA = ',nteta,nzeta
      write(6,*)' NUSC1,ITUSC = ',NUSC1,ITUSC 
      write(6,*)' u_tilde = ',U_TILDE,' f_reduced = ',freduced
      write(6,*)' Strouhal= ', Strouhal 
      write(6,*)' beta = ',BETA*180/pi 
      do i=1,Nmodes
      write(6,*)' Kc(',i,'         ) = ',Kc(i)
      write(6,*)' Alfa_tau(',i,'         ) = ',Alphatau(i)*180.d0/pi
      write(6,*)' Alfa_amp(',i,'         ) = ',Alphaamp(i)*180.d0/pi
      write(6,*)' Alfa_phi(',i,'         ) = ',Alphaphi(i)*180.d0/pi
      enddo
      write(6,*)' Re = ',Reyn,' REc = ',reyn*4/freduced
      write(6,*)' Nt ciclo e t fin = ',NTciclo,Tfin
      dt=2.d0*PI/float(NTciclo)
      write(6,*)' passo temporale dt = ',dt 
      write(6,*)' weightmu, weightnu = ',weightmu,weightnu
      write(6,*)' weightgamma, weightdelta = ',weightgamma,weightdelta
      write(6,*) 't iniziale = ', t
     
      indgdt=1
      gradientak(1)=100.d0
      savegdta=0.d0


  
     

C      do while (abs(gradientak(1)).GT.0.05)
    
      do while (indgdt.LT.2)

     
      indmod=1
     
      indforce=100.d0*indgdt+10.d0+indmod            ! indmod a enlever??
      indfa=100.d0*indgdt+20.d0+indmod
      indea=100.d0*indgdt+60.d0+indmod
      indpf=100.d0*indgdt+7


      CALL INIT

   

      IUSC1=0
      ITOUT=0


      l=lfun(t,Kc,alphatau,Nmodes)
      alpha=alphafun(t,alphaamp,alphaphi,Nmodes)
      dl=dlfun(t,Kc,alphatau,Nmodes)
      dalpha=dalphafun(t,alphaamp,alphaphi,Nmodes)
      d2l=d2lfun(t,Kc,alphatau,Nmodes)
      d2alpha=d2alphafun(t,alphaamp,alphaphi,Nmodes)

   

  

     
      CALL VEL


      OPEN(indforce)              !proforze
      OPEN(indfa)                 !proforcesensak


 
10    CONTINUE

   
      ITOUT=ITOUT+1     
      IUSC1=IUSC1+1

 
      if (ITOUT.EQ.ITUSC) then
        write(6,*)'tempo =',t,' omg(nteta/2+1,1) ='
     $        ,omg(nteta/2+1,1,1)
         write(6,*)'tempo =',t,' omgga(nteta/2+1,1) ='
     $        ,omgga(nteta/2+1,1,1,indmod)
        write(6,*) '                                                   '
C        write(6,*)'tempo =',t,' psi(nteta/2+1,1) ='
C     $        ,psi(nteta/2+1,1)
C        write(6,*)'tempo =',t,' psiga(nteta/2+1,1) ='
C     $        ,psiga(nteta/2+1,1,indmod)
        write(6,*) '***************************************************'
        write(6,*) '***************************************************'
        itout=0
      endif

      IF (T.LE.TFIN) THEN



          CALL VORTICITY
          call vortsensak


 
       
        T=T+DT
        l=lfun(t,Kc,alphatau,Nmodes)
        alpha=alphafun(t,alphaamp,alphaphi,Nmodes)
        dl=dlfun(t,Kc,alphatau,Nmodes)
        dalpha=dalphafun(t,alphaamp,alphaphi,Nmodes)
        d2l=d2lfun(t,Kc,alphatau,Nmodes)
        d2alpha=d2alphafun(t,alphaamp,alphaphi,Nmodes)
    
          
          CALL POISSON
          call psisensak



     
           CALL VEL
        
          IF (IUSC1.EQ.NUSC1) THEN
          CALL forces
          call forcesensak


          IUSC1=0
          ENDIF
       




          GOTO 10
      
      ENDIF
      





    
C     Calcolo dell'efficienza

      CLOSE (indforce)               
      close (indfa)


     

      OPEN(indforce)                      !proforze
      OPEN(indfa)                         !proforcesensak



      OPEN(indpf)                              !proefficienza
      OPEN(indea)                              !proeffisensak


      
      CALL EFFIC
      rewind(indforce)
      call efficsensak(sumak)
      CLOSE (indforce)
      close (indfa)


      gradientak(indmod)=sumak


      CLOSE (indpf)
     
      write(9,67)indgdt,indmod,beta*180/pi,Kc(indmod)
     &     ,alphaamp(indmod)*180/pi,alphatau(indmod)*180/pi,
     &           alphaphi(indmod)*180/pi,costfct,gradientak(indmod),
     &           CLnew,CPnew,CTnew,effinew 

 
      

      call magnitudeak(gradientak(indmod),savegdta,reponseak)



      savea=Kc(indmod)
      savegdta=gradientak(indmod)
      Kc(indmod)= Kc(indmod)-reponseak
   
   


      Strouhal=Kc(1)*freduced/pi
      write(indea,*)'*********************************************'
      write(indea,*)'*********************************************'
      write(indea,*)'NEW a(indmod)',Kc(indmod)
      write(indea,*)'NEW Strouhal,indmod=',Strouhal,indmod

      close (indea)

      indgdt=indgdt+1
      enddo




      CLOSE (1)           ! in.dat
      CLOSE (3)           ! profilo
      CLOSE (4)           ! distribuzione punti di calcolo lungo il raggio
      close(9)
 67    FORMAT(i3,1x,i3,1x,11(E13.7,1X))            
      END

C ***********************************************************************

      SUBROUTINE cost
      INCLUDE 'sensvar.com'

      pi=3.141592654d0
      c0=(0.D0,0.D0)
      ci=(0.D0,1.D0)

      return
      end

C ***********************************************************************

      SUBROUTINE INPUT
      INCLUDE 'sensvar.com'
      
      READ(1,*)IOPT
      READ(1,*)ICIL        ! 0 - se profilo joukowski; 1 - se cilindro
      READ(1,*)DIS
      READ(1,*)E
      READ(1,*)small
      READ(1,*)freduced
      READ(1,*)Nmodes
      READ(1,*)BETA
      do i=1,Nmodes
      READ(1,*)Kc(i)
      enddo
      do i=1,Nmodes
      READ(1,*)alphatau(i)
      enddo
      do i=1,Nmodes
      READ(1,*)Alphaamp(i)
      enddo
      do i=1,Nmodes
      READ(1,*)Alphaphi(i)
      enddo
      READ(1,*)RE
      READ(1,*)Rmax
      READ(1,*)Rint
      READ(1,*)NTETA,NZETA
      READ(1,*)NTciclo,Tfin
      READ(1,*)NUSC1,ITUSC
      READ(1,*)weightmu,weightnu
      READ(1,*)weightgamma,weightdelta


      
      Strouhal=(Kc(1)*freduced)/pi
      U_TILDE=1./(freduced)
      REYN=RE     
      do k=1,Nmodes
      ALPHAtau(k)=ALPHAtau(k)*PI/180.d0
      ALPHAAMP(k)=ALPHAAMP(k)*PI/180.d0
      ALPHAPHI(k)=ALPHAPHI(k)*PI/180.d0
      enddo     
      beta=beta*pi/180

      RETURN
      END

C *****************************************************************************
      SUBROUTINE RETIC
      INCLUDE 'sensvar.com'

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
      ZMAX=LOG(RMAX+A)
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
         write(4,37)j,zeta(j),r(j),dzeta,dteta*r(j)
      enddo
 37   format(i3,4(1x,e13.7))

      IF(ICIL.EQ.0) THEN

      DO I=1,NTETA+1
         DO J=1,NZETA
            CHI(I,J)=CMPLX(R(J)*COS(TETA(I)),R(J)*SIN(TETA(I)))
            call conf(chi(i,j),z(i,j))
            aux1=(dreal(chi(i,j))-e)**2
            aux2=(dimag(chi(i,j)))**2
	    jacob(i,j)=1.D0+(1.D0-2.D0*(aux1-aux2))/(aux1+aux2)**2
	    jacobsq(i,j)=sqrt(jacob(i,j))
            dxdcsi(i,j)=1.D0-(aux1-aux2)/(aux1+aux2)**2 
            dxdeta(i,j)=-2.D0*( dreal(chi(i,j))-e )*( dimag(chi(i,j)) )
     *                         /(aux1+aux2)**2 
            dydcsi(i,j)=-dxdeta(i,j) 
            dydeta(i,j)=dxdcsi(i,j) 
	 ENDDO
         dcsidteta=-r(1)*sin(teta(i))        
         detadteta=r(1)*cos(teta(i))        
         aux1=dreal(chi(i,1))-e
         aux2=imag(chi(i,1))
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
         WRITE(3,22)DREAL(Z(I,1)),DIMAG(Z(I,1))
         close(3)
         open (3,FILE='profilo', status='old', position='append')
 12   CONTINUE

 22   FORMAT(E13.7,1X,E13.7)

      RETURN
      END
C ***********************************************************************

      SUBROUTINE CONF(chix,zx)
      
      INCLUDE 'sensvar.com'
      
      COMPLEX*16 chix,zx
 
      zx=chix-e+1.d0/(chix-e)+dis

      RETURN
      END

C ***********************************************************************

      SUBROUTINE INIT
      INCLUDE 'sensvar.com'
      
   
C      DT=1.d0/FLOAT(NTCICLO)

      IF (IOPT.EQ.0) THEN
         T=0.d0
         DO 10 I=0,NTETA+1
            DO 10 J=1,NZETA
               PSI(I,J)=0.d0
               OMG(I,J,1)=0.d0
               OMG(I,J,2)=0.d0
               omgga(I,J,1,indmod)=0.d0
               omgga(I,J,2,indmod)=0.d0
               psiga(i,j,indmod)=0.d0



 10      CONTINUE
      ELSE
         CALL OPFILEB(8,'OLD','Nome file iniziale ? ')      
         READ(8)T

         DO 11 I=1,NTETA
            DO 12 J=1,NZETA
               READ(8) OMG(I,J,1),PSI(I,J), omgga(i,j,1,indmod),
     $          psiga(i,j,indmod)

               OMG(I,J,2)=OMG(I,J,1)
               omgga(i,j,2,indmod)=omgga(i,j,1,indmod)
  

 12         CONTINUE
 11      CONTINUE
         CLOSE(8)
         do j=1,nzeta
            omg(0,j,1)=omg(nteta,j,1)
            omg(nteta+1,j,1)=omg(1,j,1)
            omgga(0,j,1,indmod)=omgga(nteta,j,1,indmod)
            omgga(nteta+1,j,1,indmod)=omgga(1,j,1,indmod)



            psi(0,j)=psi(nteta,j)
            psi(nteta+1,j)=psi(1,j)
            psiga(0,j,indmod)=psiga(nteta,j,indmod)
            psiga(nteta+1,j,indmod)=psiga(1,j,indmod)

           


         enddo
      ENDIF

      

      RETURN
      END

C ***********************************************************************

      SUBROUTINE VORTICITY
      
      INCLUDE 'sensvar.com'

      DIMENSION AM(ITETA+1),BM(ITETA+1),CM(ITETA+1),DM(ITETA+1),
     *           GAM(ITETA+1)
      DIMENSION AN(IZETA),BN(IZETA),CN(IZETA),DN(IZETA),
     *           WN(IZETA)

C
C    !!!!!!!!!!!!!!!! Soluzione ADI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
C     primo passo temporale implicito in teta
      
      DO J=2,NZETA-1

C
         DO I=1,NTETA

            OX1=jacobsq(i,j)*R(J)*REYN*V_TETA(I,J)*DTETA*0.5D0
            OX2=jacob(i,j)*r2(j)*(reyn)*2.D0*dteta*dteta/dt
            AM(I)=-(1.D0+ox1)
            BM(I)=2.D0+ox2
            CM(I)=-(1.D0-ox1)
            ox3=0.5*a*dzeta/r(j)
            ox4=jacobsq(i,j)*ez(j)*reyn*v_r(i,j)*dzeta*0.5
            ox5=-(dteta*r(j)/(dzeta*ez(j)))**2
            DM(I)=-omg(i,j-1,1)*ox5*(1.d0-ox3+ox4)
     *           +omg(i,j,1)*ox5*(2.d0-
     *           jacob(i,j)*ez2(j)*(reyn)*2.d0*dzeta*dzeta/dt)
     *           -omg(i,j+1,1)*ox5*(1.d0+ox3-ox4)              
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

      DO I=1,NTETA
         DO J=2,NZETA-1 
            PX1=jacobsq(i,j)*Ez(j)*REYN*V_R(I,J)*DZETA*0.5D0
            PX2=0.5D0*A*DZETA/R(J)
            AN(J)=(1.D0-PX1+PX2)
            BN(J)=2.D0+jacob(i,j)*ez2(j)*(reyn)*2.D0*DZETA*DZETA/DT
            CN(J)=(1.D0+PX1-PX2)
            px4=-(dzeta*ez(j)/(dteta*r(j)))**2
            px5=jacobsq(i,j)*r(j)*reyn*v_teta(i,j)*dteta*0.5d0
            dn(j)=-omg(i-1,j,1)*px4*(1.d0+px5)
     *           +omg(i,j,1)*px4*
     *           (2.d0-jacob(i,j)*r2(j)*(reyn)*2.d0*dteta*dteta/dt)
     *           -omg(i+1,j,1)*px4*(1.d0-px5)

         ENDDO  
       

	 CALL BOUNDOMG(OMG(I,1,2),I)
         CALL TRIG(AN,BN,CN,DN,WN,NZETA,OMG(I,1,2))
         DO J=1,NZETA
            OMG(I,J,2)=WN(J)
         ENDDO
      ENDDO


      
      do j=2,nzeta-1
          OMG(0,J,2)=OMG(nteta,J,2)
          OMG(nteta+1,J,2)=OMG(1,J,2)
      enddo

      DO I=0,NTETA+1

         DO J=1,NZETA
            OMG(I,J,1)=OMG(I,J,2)
         ENDDO
      ENDDO



      RETURN
      END
  
C ***************************************************************************

      SUBROUTINE BOUNDOMG(W,I)


      INCLUDE 'sensvar.com'


      V1=( -dl*SIN(ALPHA)+DALPHA*DIMAG(Z(I,1)) )
      V2=( -dl*COS(ALPHA)-DALPHA*DREAL(Z(I,1)) )

      VA1=-dxdcsi(i,1)*sin(TETA(I))+DXDeta(i,1)*cos(TETA(I))
      VA2=-dydcsi(i,1)*sin(TETA(I))+dydeta(i,1)*cos(TETA(I))

      W=dalpha
     *  -2.d0*(PSI(I,2)-PSI(I,1))/( jacob(i,1)*DZETA*DZETA*Ez2(1) )
     *  +( 2.d0/DZETA-A/R(1) )*(V1*VA1+V2*VA2)/( Jacob(I,1)*Ez(1) )
     *  +v1*(-sin(TETA(I))*dxdcsi(i,1)+cos(TETA(I))*dxdeta(i,1)
     *       +cos(TETA(I))*d2xdcsiteta(i)+sin(TETA(I))*d2xdetateta(i))
     *      /( R(1)*Jacob(i,1) )
     *  +v2*(-sin(TETA(I))*dydcsi(i,1)+cos(TETA(I))*dydeta(i,1)
     *       +cos(TETA(I))*d2ydcsiteta(i)+sin(TETA(I))*d2ydetateta(i))
     *      /( R(1)*Jacob(i,1) )  
      

   
      RETURN
      END




C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************

      
      
      SUBROUTINE vortsensak
      
      INCLUDE 'sensvar.com'

      DIMENSION AM(ITETA+1),BM(ITETA+1),CM(ITETA+1),DM(ITETA+1),
     *          GAM(ITETA+1)
      DIMENSION AN(IZETA),BN(IZETA),CN(IZETA),DN(IZETA),
     *          WN(IZETA)




C    !!!!!!!!!!!!!!!! Soluzione ADI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
C     primo passo temporale implicito in teta
      
      DO J=2,NZETA-1

C
         DO I=1,NTETA
 
            call termsourceak(t,i,j,indmod,termSteta,termSzeta)
          
            OX1=jacobsq(i,j)*R(J)*REYN*V_TETA(I,J)*DTETA*0.5D0
            OX2=jacob(i,j)*r2(j)*reyn*2.D0*dteta*dteta/dt
            
            AM(I)=-(1.D0+ox1)
            BM(I)=2.D0+ox2
            CM(I)=-(1.D0-ox1)
            ox3=0.5*a*dzeta/r(j)
            ox4=jacobsq(i,j)*ez(j)*reyn*v_r(i,j)*dzeta*0.5
            ox5=-(dteta*r(j)/(dzeta*ez(j)))**2
            ox6=jacob(i,j)*ez2(j)*(reyn)*2.d0*dzeta*dzeta/dt          
            DM(I)=-omgga(i,j-1,1,indmod)*ox5*(1.d0-ox3+ox4)
     *            +omgga(i,j,1,indmod)*ox5*(2.d0-ox6)
     *            -omgga(i,j+1,1,indmod)*ox5*(1.d0+ox3-ox4)
     *            -r2(j)*reyn*dteta*dteta*termSteta
           

         ENDDO
      
         CALL TRIPER(AM,BM,CM,DM,NTETA-1,GAM)
      
         DO I=1,NTETA
            omgga(I,J,2,indmod)=DM(I)
             
         ENDDO
      ENDDO
      
      DO J=2,NZETA-1
	 DO I=1,NTETA
            omgga(I,J,1,indmod)=omgga(i,j,2,indmod)
         ENDDO
      ENDDO
     
    
      do j=2,nzeta-1
          omgga(0,J,1,indmod)=omgga(nteta,J,2,indmod)
          omgga(nteta+1,J,1,indmod)=omgga(1,J,2,indmod)
      enddo

      DO I=1,NTETA
         DO J=2,NZETA-1 


            call termsourceak(t,i,j,indmod,termSteta,termSzeta)
            PX1=jacobsq(i,j)*Ez(j)*REYN*V_R(I,J)*DZETA*0.5D0
            PX2=0.5D0*A*DZETA/R(J)
            AN(J)=(1.D0-PX1+PX2)
            BN(J)=2.D0+jacob(i,j)*ez2(j)*(reyn)*2.D0*DZETA*DZETA/DT
            CN(J)=(1.D0+PX1-PX2)
            px4=-(dzeta*ez(j)/(dteta*r(j)))**2
            px5=jacobsq(i,j)*r(j)*reyn*v_teta(i,j)*dteta*0.5d0
            dn(j)=-omgga(i-1,j,1,indmod)*px4*(1.d0+px5)
     *           +omgga(i,j,1,indmod)*px4*
     *           (2.d0-jacob(i,j)*r2(j)*(reyn)*2.d0*dteta*dteta/dt)
     *           -omgga(i+1,j,1,indmod)*px4*(1.d0-px5)
     *           +(0.5d0*reyn*ez(j)*dzeta)*termSzeta

        ENDDO  
       

	 CALL Boundomggak(omgga(I,1,2,indmod),I)

         CALL TRIG(AN,BN,CN,DN,WN,NZETA,omgga(I,1,2,indmod))
         DO J=1,NZETA
            omgga(I,J,2,indmod)=WN(J)
        ENDDO
      ENDDO


      
      do j=2,nzeta-1
          Omgga(0,J,2,indmod)=Omgga(nteta,J,2,indmod)
          Omgga(nteta+1,J,2,indmod)=Omgga(1,J,2,indmod)
      enddo

      DO I=0,NTETA+1

         DO J=1,NZETA
            Omgga(I,J,1,indmod)=Omgga(I,J,2,indmod)
         ENDDO
      ENDDO



      RETURN
      END





      real*8 function dhcosaxdak(time,k,alpha,tau,Nmod)
      real*8 time,alpha
      integer*8 k,Nmod
      real*8 tau(Nmod)
      dhcosaxdak=k*dcos(k*time+tau(k))*cos(alpha)
      return
      end

      real*8 function dhsinaydak(time,k,alpha,tau,Nmod)
      real*8 time,alpha
      real*8 tau(Nmod)
      integer*8 k,Nmod
      dhsinaydak=k*dcos(k*time+tau(k))*sin(alpha)     
      return 
      end


      
      real*8 function dhdtdak(time,k,tau,Nmod)
      real*8 time
      integer*8 k,Nmod
      real*8 tau(Nmod)
      dhdtdak=k*dcos(k*time+tau(k))
      return
      end
      
      real*8 function d2hdt2dak(time,k,tau,Nmod)
      real*8 time
      integer*8 k,Nmod
      real*8 tau(Nmod)
      d2hdt2dak=-k*k*dsin(k*time+tau(k))
      return
      end
      
      
      real*8 function dhsqdak(time,k,amp,tau,Nmod)
      real*8 time
      real*8 amp(Nmod),tau(Nmod)
      integer*8 k,Nmod
      aux=0
      do m=1,Nmod
      aux=aux+2.d0*amp(m)*dsin(m*time+tau(m))
      enddo
      dhsqdak=aux*dsin(k*time+tau(k))
      return
      end





     




      subroutine termsourceak(time,i,j,k,termSteta,termSzeta)
      INCLUDE 'sensvar.com'

      real*8  xchixip,xchixim,ychixip,ychixim
      real*8  dvrdak,dvtetadak,termS,souci
      integer*8 k

      k=indmod
         


      xchixip=dxdcsi(i,j)*cos(teta(i))+dxdeta(i,j)*sin(teta(i))
      xchixim=dxdeta(i,j)*cos(teta(i))-dxdcsi(i,j)*sin(teta(i))
      ychixip=dydcsi(i,j)*cos(teta(i))+dydeta(i,j)*sin(teta(i))
      ychixim=dydeta(i,j)*cos(teta(i))-dydcsi(i,j)*sin(teta(i))
      



      dvrdak=-dhsinaydak(time,k,alpha,alphatau,Nmodes)*xchixip
     *     -dhcosaxdak(time,k,alpha,alphatau,Nmodes)
     *    *ychixip+(1.d0/(2.d0*r(j)*dteta))
     *    *(psiga(i+1,j,indmod)-psiga(i-1,j,indmod))


      dvtetadak=-dhsinaydak(time,k,alpha,alphatau,Nmodes)*xchixim-
     *     dhcosaxdak(time,k,alpha,alphatau,Nmodes)*ychixim
     *     -(1.d0/(2.d0*ez(j)*dzeta))
     *     *(psiga(i,j+1,indmod)-psiga(i,j-1,indmod))






      termSteta=(1/(2*dzeta*ez(j)))*(omg(i,j+1,1)-omg(i,j-1,1))*dvrdak+
     *      (1/(2*dteta*r(j)))* (omg(i+1,j,1)-omg(i-1,j,1))*dvtetadak

      termSzeta=-dvrdak*   (omg(i,j+1,1)-omg(i,j-1,1))
     *          -dvtetadak*(omg(i+1,j,1)-omg(i-1,j,1))
     *          *(ez(j)*dzeta)/(r(j)*dteta)
 

      return
      end




      SUBROUTINE Boundomggak(WW,I)

      INCLUDE 'sensvar.com'


      V1=dhsinaydak(t,indmod,alpha,alphatau,Nmodes)
      V2=dhcosaxdak(t,indmod,alpha,alphatau,Nmodes)

      VA1=-dxdcsi(i,1)*sin(TETA(I))+DXDeta(i,1)*cos(TETA(I))
      VA2=-dydcsi(i,1)*sin(TETA(I))+dydeta(i,1)*cos(TETA(I))

      WW=-2.d0*(psiga(I,2,indmod)-psiga(I,1,indmod))/
     *   ( jacob(i,1)*DZETA*DZETA*Ez2(1) )
     *  +(A/R(1)-2.d0/DZETA )*(V1*VA1+V2*VA2)/( Jacob(I,1)*Ez(1) )
     *  +v1*(sin(TETA(I))*dxdcsi(i,1)-cos(TETA(I))*dxdeta(i,1)
     *       -cos(TETA(I))*d2xdcsiteta(i)-sin(TETA(I))*d2xdetateta(i))
     *      /( R(1)*Jacob(i,1) )
     *  +v2*(sin(TETA(I))*dydcsi(i,1)-cos(TETA(I))*dydeta(i,1)
     *       -cos(TETA(I))*d2ydcsiteta(i)-sin(TETA(I))*d2ydetateta(i))
     *      /( R(1)*Jacob(i,1) )  
      

   
      RETURN
      END





      SUBROUTINE psisensak
      
      INCLUDE 'sensvar.com'

      COMPLEX*16 A1,AM,DN,WN,MATPSI1,Y,DY1,DYM,MATPSIM
      
      DIMENSION AN(IZETA),BN(IZETA),CN(IZETA),DN(IZETA), WN(IZETA)
       
      DIMENSION MATPSI1(ITETA),Y(ITETA),DY1(ITETA)
      DIMENSION MATPSIM(ITETA),DYM(ITETA)

      DO J=2,NZETA-1
         DO I=1,NTETA 
            Y(I)=CMPLX(jacob(i,j)*Omgga(I,J,1,indmod),0.d0)           
         ENDDO
         CALL CFFT1(Y,NTETA,1)
         DO I=1,NTETA/2+1
            MATOMgga(I,J,indmod)=Y(I)
         ENDDO
      ENDDO
      
      DO I=1,NTETA
         CALL BOUNDPSIgak(DY,DM,I)
         DY1(I)=CMPLX(DY,0.d0)
         DYM(I)=CMPLX(DM,0.d0)
      ENDDO
         
      CALL CFFT1(DY1,NTETA,1)
      CALL CFFT1(DYM,NTETA,1)

C     RISOLVO LE ARMONICHE NEGATIVE DA -1 A -NTETA/2
      MATPSI1(1)=C0
      MATPSIM(1)=C0
!
      DO I=2,NTETA/2
         MATPSI1(I)=DY1(I)/CMPLX(0.d0,-FLOAT(I-1))
         MATPSIM(I)=DYM(I)/CMPLX(0.d0,-FLOAT(I-1))
      ENDDO
!
      DO I=1,NTETA/2
         DO J=2,NZETA-1 
            AN(J)=1.d0+A*DZETA/(2.d0*R(J))
            BN(J)=2.d0+(EXP(ZETA(J))*DZETA*(-FLOAT(I-1))/R(J))**2
            CN(J)=1.d0-A*DZETA/(2.d0*R(J))
            DN(J)=MATOMGga(I,J,indmod)*(DZETA*EXP(ZETA(J)))**2
         ENDDO
         A1=MATPSI1(I)
         AM=MATPSIM(I)
         CALL TRIGCMPLX(AN,BN,CN,DN,WN,NZETA,A1,AM)
         DO J=1,NZETA
            MATPSIga(I,J,indmod)=WN(J)
         ENDDO
      ENDDO

      DO J=1,NZETA
         Y(1)=MATPSIga(1,J,indmod)
         DO I=1,NTETA/2-1
            Y(I+1)=MATPSIga(I+1,J,indmod)
            Y(NTETA+1-I)=CONJG(Y(I+1))
         ENDDO
         Y(NTETA/2+1)=C0     
         CALL CFFT1(Y,NTETA,-1)
         DO I=1,NTETA
            PSIga(I,J,indmod)=DREAL(Y(I))
            ENDDO
      ENDDO
    
      DO J=1,NZETA
         PSIga(0,J,indmod)=PSIga(NTETA,J,indmod)
         PSIga(NTETA+1,J,indmod)=PSIga(1,J,indmod)
      ENDDO
      
      RETURN
      END

C *******************************************************************

      SUBROUTINE Boundpsigak(W,WW,I)
      INCLUDE 'sensvar.com'
      integer*8 k
      k=indmod
  
 
      V1=dhsinaydak(t,k,alpha,alphatau,Nmodes)
      V2=dhcosaxdak(t,k,alpha,alphatau,Nmodes)

      VA1=DXDcsi(i,1)*COS(TETA(I))+DXDETA(i,1)*SIN(TETA(I))
      VA2=DYDcsi(i,1)*COS(TETA(I))+DYDeta(i,1)*SIN(TETA(I))
      W=R(1)*(V1*VA1+V2*VA2)

       WW=0.d0
      RETURN
      END


      SUBROUTINE Forcesensak
      INCLUDE 'sensvar.com'
      COMPLEX*16 Q,QQ,S,SS
      DIMENSION Q(iteta),QQ(iteta),S(iteta),SS(iteta)



      DO I=1,NTETA
      dpdtetadak(i)=-d2hdt2dak(t,indmod,alphatau,Nmodes)*
     $    (sin(alpha)*DXDTETA(I)+cos(alpha)*DYDTETA(I))

         var1=-DXDETA(I,1)*cos(TETA(I))+DXDCSI(I,1)*sin(TETA(I))
         var2=DXDETA(I,1)*SIN(TETA(I))+DXDCSI(I,1)*COS(TETA(I))
         var3=DYDETA(I,1)*cos(TETA(I))-DYDCSI(I,1)*Sin(TETA(I))
         var4=DYDETA(I,1)*SIN(TETA(I))+DYDCSI(I,1)*COS(TETA(I))
         var5=3.d0*omgga(i,1,1,indmod)-4.d0*omgga(i,2,1,indmod)
     $       +omgga(i,3,1,indmod)
         var6=4.d0*omgga(i,2,1,indmod)-3.d0*omgga(i,1,1,indmod)
     $       -omgga(i,3,1,indmod)
        dpdwdakdteta(i)=(1.d0/(Reyn*jacob(i,1)))*(
     $  dxdteta(i)*var1*
     $   (var5)/(2.d0*dzeta*ez(1))
     $  +(1.d0/r(1))*dxdteta(i)*var2*
     $  (omgga(i-1,1,1,indmod)-omgga(i+1,1,1,indmod))/(2*dteta)
     $  +dydteta(i)*var3*
     $   (var6)/(2.d0*dzeta*ez(1))
     $  -(1.d0/r(1))*dydteta(i)*var4*
     $  (omgga(i+1,1,1,indmod)-omgga(i-1,1,1,indmod))/(2*dteta)
     $   )
     
      Q(I)=CMPLX(DPDTETAdak(i),0.d0)
      S(i)=cmplx(dpdwdakdteta(i),0.d0)
      enddo


      CALL CFFT1(Q,NTETA,1)
      CALL CFFT1(S,NTETA,1)
      
      QQ(1)=C0
      SS(1)=C0

      DO I=1,NTETA/2-1

         QQ(I+1)=Q(I+1)/CMPLX(0.d0,-FLOAT(I))
         QQ(NTETA+1-I)=CONJG(QQ(I+1))
         SS(I+1)=S(I+1)/CMPLX(0.d0,-FLOAT(I))
         SS(NTETA+1-I)=CONJG(SS(I+1))

      ENDDO

      QQ(NTETA/2+1)=C0
      SS(NTETA/2+1)=C0
      
      Call CFFT1(QQ,nteta,-1)
      Call CFFT1(SS,nteta,-1)


      DO I=1,NTETA
      dpdak(i)=Dreal(QQ(I))
      dpdwdak(i)=Dreal(SS(I))
      
      ENDDO
      dpdak(nteta+1)=dpdak(1)
      dpdwdak(nteta+1)=dpdwdak(1)

      
      dfxdpsiga=0.d0
      dfydpsiga=0.d0
      dmzdpsiga=0.d0
      Xintegdpdw=0.d0
      Yintegdpdw=0.d0
      Zintegdpdw=0.d0
      Xintegdp=0.d0
      Yintegdp=0.d0
      Zintegdp=0.d0
      
      
      
      do i=1,Nteta
      var1=cos(2.d0*teta(i))*dxdcsi(i,1)*dydcsi(i,1)
      var2=0.5d0*sin(2.d0*teta(i))*
     &    (dydeta(i,1)*dxdcsi(i,1)+ dydcsi(i,1)*dxdeta(i,1))
      var3=cos(2.d0*teta(i))*((dxdeta(i,1))**2.d0-(dxdcsi(i,1))**2.d0)
      var4=2.d0*sin(2.d0*teta(i))*dxdeta(i,1)*dxdcsi(i,1)
      d2psigadz2=(2.d0*psiga(i,1,indmod)-5.d0*psiga(i,2,indmod)
     & +4.d0*psiga(i,3,indmod)-psiga(i,4,indmod))/(dzeta*dzeta)
      dpsigadz=(-3.d0*psiga(i,1,indmod)+4.d0*psiga(i,2,indmod)
     & -psiga(i,3,indmod))/(2.d0*dzeta)
      d2psigadr2=(1.d0/ez2(1))*(d2psigadz2-dpsigadz)





      dfxdpsiga=dfxdpsiga+(1.d0/(reyn*jacob(i,1)*jacob(i,1)))*d2psigadr2
     $*(2.d0*dydteta(i)*(var1+var2)-dxdteta(i)*(var3-var4))*dteta

      dfydpsiga=dfydpsiga+(1.d0/(reyn*jacob(i,1)*jacob(i,1)))*d2psigadr2
     $*(2.d0*dxdteta(i)*(var1+var2)+dydteta(i)*(var3-var4))*dteta
     
      dmzdpsiga=dmzdpsiga+(1.d0/(reyn*jacob(i,1)*jacob(i,1)))*d2psigadr2
     $*(2.d0*(dreal(z(i,1))*dxdteta(i)-dimag(z(i,1))*dydteta(i))*
     $(var1+var2)+
     $(dreal(z(i,1))*dydteta(i)+dimag(z(i,1))*dxdteta(i))*
     $(var3-var4)   )*dteta
     
      Xintegdpdw=Xintegdpdw+dpdwdak(i)*dxdteta(i)*dteta
      
      Yintegdpdw=Yintegdpdw+dpdwdak(i)*dydteta(i)*dteta

      Zintegdpdw=Zintegdpdw+(dreal(z(i,1))*dxdteta(i)
     $+dimag(z(i,1))*dydteta(i))*dpdwdak(i)*dteta
     
      Xintegdp=Xintegdp+dpdak(i)*dxdteta(i)*dteta

      Yintegdp=Yintegdp+dpdak(i)*dydteta(i)*dteta

      Zintegdp=Zintegdp+(dreal(z(i,1))*dxdteta(i)
     $+dimag(z(i,1))*dydteta(i))*dpdak(i)*dteta
     

      enddo
      

      
      WRITE(indfa,55)T,dfxdpsiga,dfydpsiga,dmzdpsiga,Xintegdpdw,
     &Yintegdpdw,Zintegdpdw,Xintegdp,Yintegdp,Zintegdp
      close (indfa)
      open (indfa, status='old', position='append')
55    FORMAT(10(E13.7,1X))


      return
      end




      SUBROUTINE efficsensak(sum)
      INCLUDE 'sensvar.com'

      REAL*8 Mz_OLD,Mm,l_old,l_new
      integer*8 k


      delta_t=0.D0
      ntot=25000
      k=indmod

      term1 =  0.D0
      term2 =  0.D0
      term3 =  0.D0
      term4 =  0.D0
      term5 =  0.D0
      term6 =  0.D0
      term7 =  0.D0
      term8 =  0.D0
      term9 =  0.D0
      term10=  0.D0
      term11=  0.D0
      term12 = 0.D0
      term13 = 0.D0
      term14 = 0.D0
      term15 = 0.D0
      term16 = 0.D0
      term17 = 0.D0
      term18 = 0.D0


      termPbar   = 0.D0
      termFbar   = 0.D0
      termLbar   = 0.D0
      termalfbar = 0.D0
      termhbar   = 0.D0





 10   read(indfa,55)T_act,dfxdpsiga,dfydpsiga,dmzdpsiga,Xintegdpdw,
     &Yintegdpdw,Zintegdpdw,Xintegdp,Yintegdp,Zintegdp

      read(indforce,25)T_ac,F_X,F_Y,M_Z,F_XHAT,
     &F_YHAT,F_XHATa,F_YHATa,F_XHATf,F_YHATf,F_Xnew, F_Ynew




      if (t_act.lt.pi*1.99) goto 10


      t_old=t_act
      fx_old=F_X
      fy_old=F_Y
      Mz_old=M_Z
      FXbeta=F_Xnew
      FYbeta=F_Ynew
      dfx_old=dfxdpsiga
      dfy_old=dfydpsiga
      dmz_old=dmzdpsiga
      xpw_old=Xintegdpdw
      ypw_old=Yintegdpdw
      zpw_old=Zintegdpdw
      xp_old=Xintegdp
      yp_old=Yintegdp
      zp_old=Zintegdp
      l_old=lfun(t_old,Kc,alphatau,Nmodes)
      dl_old=dlfun(t_old,Kc,alphatau,Nmodes)
      alpha_old=alphafun(t_old,alphaamp,alphaphi,Nmodes)
      dalpha_old=dalphafun(t_old,alphaamp,alphaphi,Nmodes)
      dhdt_old=dhdtdak(t_old,k,alphatau,Nmodes)
      dh2_old=dhsqdak(t_old,k,Kc,alphatau,Nmodes)



      do i=1,ntot

      read(indfa,55,end=20,err=20)T_act,dfxdpsiga,dfydpsiga,
     &dmzdpsiga,Xintegdpdw,Yintegdpdw,Zintegdpdw,Xintegdp,
     &Yintegdp,Zintegdp

      read(indforce,25,end=20,err=20)T_act,F_X,F_Y,M_Z,F_XHAT,
     $     F_YHAT,F_XHATa,F_YHATa,F_XHATf,F_YHATf,F_Xnew,F_Ynew

         dtloc=t_act-t_old
         alpha_new=alphafun(t_act,alphaamp,alphaphi,Nmodes)
         dhdt_new=dhdtdak(t_act,k,alphatau,Nmodes)
         dh2_new=dhsqdak(t_act,k,Kc,alphatau,Nmodes)
         l_new=lfun(t_act,Kc,alphatau,Nmodes)
         dl_new=dlfun(t_act,Kc,alphatau,Nmodes)
         dalpha_new=dalphafun(t_act,alphaamp,alphaphi,Nmodes)




        term1= term1+0.5d0*dtloc*
     & (ypw_old*sin(alpha_old-beta)*dl_old+
     &  Yintegdpdw*sin(alpha_new-beta)*dl_new)

        term2 = term2-0.5d0*dtloc*
     & (xpw_old*cos(alpha_old-beta)*dl_old+
     &  Xintegdpdw*cos(alpha_new-beta)*dl_new)

        term3 = term3-0.5d0*dtloc*
     &     (zpw_old*dalpha_old+Zintegdpdw*dalpha_new)

        term4 = term4-0.5d0*dtloc*
     &     (ypw_old*cos(alpha_old-beta)+Yintegdpdw*cos(alpha_new-beta))

        term5= term5-0.5d0*dtloc*
     &     (xpw_old*sin(alpha_old-beta)+Xintegdpdw*sin(alpha_new-beta))

        term6 = term6-0.5d0*dtloc*
     &   (dfx_old*sin(alpha_old-beta)*dl_old+
     &    dfxdpsiga*sin(alpha_new-beta)*dl_new)

        term7 = term7-0.5d0*dtloc*
     &   (dfy_old*cos(alpha_old-beta)*dl_old+
     &    dfydpsiga*cos(alpha_new-beta)*dl_new)

        term8 = term8-0.5d0*dtloc*
     &   (dmz_old*dalpha_old+dmzdpsiga*dalpha_new)

        term9 = term9+0.5d0*dtloc*
     &   (dfx_old*cos(alpha_old-beta)+dfxdpsiga*cos(alpha_new-beta))

        term10= term10-0.5d0*dtloc*
     &   (dfy_old*sin(alpha_old-beta)+dfydpsiga*sin(alpha_new-beta))

        term11= term11-0.5d0*dtloc*
     &  (fx_old*sin(alpha_old-beta)*dhdt_old+
     &   F_X*sin(alpha_new-beta)*dhdt_new)

        term12= term12-0.5d0*dtloc*
     &  (fy_old*cos(alpha_old-beta)*dhdt_old+
     &   F_Y*cos(alpha_new-beta)*dhdt_new)

        term13= term13+0.5d0*dtloc*
     &  (yp_old*sin(alpha_old-beta)*dl_old+
     &   Yintegdp*sin(alpha_new-beta)*dl_new)

        term14= term14-0.5d0*dtloc*
     &  (xp_old*cos(alpha_old-beta)*dl_old+
     &   Xintegdp*cos(alpha_new-beta)*dl_new)

        term15= term15-0.5d0*dtloc*
     &  (zp_old*dalpha_old+Zintegdp*dalpha_new)

        term16= term16-0.5d0*dtloc*
     &  (yp_old*cos(alpha_old-beta)+Yintegdp*cos(alpha_new-beta))

        term17= term17-0.5d0*dtloc*
     &  (xp_old*sin(alpha_old-beta)+Xintegdp*sin(alpha_new-beta))

        term18= term18+0.5d0*dtloc*
     &  (dh2_old+dh2_new)


        termalfbar= termalfbar+0.5d0*dtloc*
     &  (alpha_old*alpha_old+alpha_new*alpha_new)




        termhbar= termhbar+0.5d0*dtloc*
     &  (l_old*l_old+l_new*l_new)



        termFbar= termFbar+0.5d0*dtloc*
     &  (FXbeta+F_Xnew)

        termLbar=termLbar+0.5d0*dtloc*
     &  (FYbeta+F_Ynew)
    


        termPbar= termPbar+0.5d0*dtloc*
     & (FYbeta*dl_old+F_Ynew*dl_new
     & +Mz_old*dalpha_old+M_z*dalpha_new)
 



         delta_t = delta_t+dtloc



         t_old=t_act
         fx_old=F_X
         fy_old=F_Y
         Mz_old=M_Z
         FXbeta=F_Xnew
         FYbeta=F_Ynew
         dfx_old=dfxdpsiga
         dfy_old=dfydpsiga
         dmz_old=dmzdpsiga
         xpw_old=Xintegdpdw
         ypw_old=Yintegdpdw
         zpw_old=Zintegdpdw
         xp_old=Xintegdp
         yp_old=Yintegdp
         zp_old=Zintegdp
         l_old=lfun(t_old,Kc,alphatau,Nmodes)
         dl_old=dlfun(t_old,Kc,alphatau,Nmodes)
         alpha_old=alphafun(t_old,alphaamp,alphaphi,Nmodes)
         dalpha_old=dalphafun(t_old,alphaamp,alphaphi,Nmodes)
         dhdt_old=dhdtdak(t_old,k,alphatau,Nmodes)
         dh2_old=dhsqdak(t_old,k,Kc,alphatau,Nmodes)


      enddo
 20   continue



      term1  = term1/delta_t
      term2  = term2/delta_t
      term3  = term3/delta_t
      term4  = term4/delta_t
      term5  = term5/delta_t
      term6  = term6/delta_t
      term7  = term7/delta_t
      term8  = term8/delta_t
      term9  = term9/delta_t
      term10 = term10/delta_t
      term11 = term11/delta_t
      term12 = term12/delta_t
      term13 = term13/delta_t
      term14 = term14/delta_t
      term15 = term15/delta_t
      term16 = term16/delta_t
      term17 = term17/delta_t
      term18 = term18/delta_t
     
      termalfbar= termalfbar/delta_t
      termhbar=   termhbar/delta_t
      termPbar=-  termPbar/delta_t
      termFbar=   termFbar/delta_t
      termLbar=   termLbar/delta_t



      summu= term1+term2+term3+term6+term7+term8+term11+term12+term13+
     &      term14+term15
      sumnu=term4+term5+term9+term10+term16+term17


      sumgamma=term18


      sumdelta=0.d0

      termmu=termPbar

      termnu=termFbar

      termgamma=termhbar

      termdelta=termalfbar

      costfct=weightmu*termmu+weightnu*termnu*U_tilde
     &       +weightgamma*termgamma+weightdelta*termdelta
     
      sum= weightmu*summu+weightnu*sumnu*U_tilde
     &       +weightgamma*sumgamma+weightdelta*sumdelta


      if (termnu.LT.0) then
         efficiency=-(termnu*U_tilde)/termmu
         else
          efficiency=0
          endif


      WRITE(indea,*)' DELTA_T = ',DELTA_T
      write(indea,*)' term1, term2, term3 ',term1, term2, term3
      write(indea,*)' term4, term5, term6' ,term4, term5, term6
      write(indea,*)' term7, term8, term9 ',term7, term8, term9
      write(indea,*)' term10, term11, term12',term10,term11,term12
      write(indea,*)' term13, term14, term15',term13,term14,term15
      write(indea,*)' term16, term17, term18',term16,term17,term18
      write(indea,*)'*********************************************'
      write(indea,*)'*********************************************'
      write(indea,*)' cost fct ', costfct,    '*** gradient', sum
      write(indea,*)' Power    ',termmu,   '*** powera', summu
      write(indea,*)' Thrust   ',termnu,   '*** thrusta', sumnu
      write(indea,*)' heaveamp ',termgamma,'** heavea ',sumgamma
      write(indea,*)' pitchamp ',termdelta,'*** pitcha', sumdelta
      write(indea,*)' Lift     ',termLbar,'*** efficiency', efficiency
      write(indea,*)'*********************************************'
      write(indea,*)'  dt, indgdt,a(indmod) ', dt, indgdt, Kc(indmod)
      write(indea,*)' b(indmod) , phi(indmod), tau(indmod),beta', 
     & alphaamp(indmod)*180/pi, alphaphi(indmod)*180/pi,
     & alphatau(indmod)*180/pi, beta*180/pi



 25    FORMAT(12(E13.7,1X))
 55    FORMAT(10(E13.7,1X))
      RETURN
      END



      subroutine magnitudeak(chiffre,chiffreold,order)
      INCLUDE 'sensvar.com'

      if (indgdt.EQ.1) then
        order=0.25*(chiffre/abs(chiffre))
        else
        order=((Kc(indmod)-savea)/(chiffre-chiffreold))*chiffre
      endif

     
      return
      end

 
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
C ***********************************************************************
























C ***********************************************************************

      SUBROUTINE POISSON
      
      INCLUDE 'sensvar.com'

      COMPLEX*16 A1,AM,DN,WN,MATPSI1,Y,DY1,DYM,MATPSIM
      
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
            PSI(I,J)=DREAL(Y(I))
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
      
      INCLUDE 'sensvar.com'
 
      V1=-DL*SIN(ALPHA)+DALPHA*DIMAG(Z(I,1))
      V2=-DL*COS(ALPHA)-DALPHA*DREAL(Z(I,1))

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

      INCLUDE 'sensvar.com'
 
C     calcolo delle componenti di velocita' al passo n
 
      DO I=1,NTETA
         DO J=2,NZETA-1
            V1=-DL*SIN(ALPHA)+DALPHA*DIMAG(Z(I,J))
            V2=-DL*COS(ALPHA)-DALPHA*DREAL(Z(I,J))
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


 27   FORMAT(5(E13.7,1X))

      DO I=1,NTETA
         V_R(I,1)=0.d0
         V_TETA(I,1)=0.d0
      ENDDO

      DO I=1,NTETA
         Ugr=U_TILDE*COS(BETA)*COS(ALPHA)+(U_TILDE*SIN(BETA)-DL)*
     $        SIN(ALPHA)+DIMAG(Z(I,NZETA))*DALPHA
         Vgr=-U_TILDE*COS(BETA)*SIN(ALPHA)+(U_TILDE*SIN(BETA)-DL)*
     $        COS(ALPHA)-DREAL(Z(I,NZETA))*DALPHA
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
C ***********************************************************************


      SUBROUTINE FORCES
      INCLUDE 'sensvar.com'

      REAL*8 LAPL_U,LAPL_V
      COMPLEX*16 Y
      DIMENSION Y(ITETA)


      
      reyninv=1.d+0/reyn

      P0=0.d0
C      P(1)=P0

      PNEW(1)=P0
      PNEW(NTETA+1)=P0

      DO I=1,NTETA
    
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

         DPDX(I)=(D2ALPHA*DIMAG(Z(I,1))+DREAL(Z(I,1))*DALPHA**2-
     *                   d2lfun(t,Kc,alphatau,Nmodes)*SIN(ALPHA))
         DPDY(I)=(-D2ALPHA*DREAL(Z(I,1))+DIMAG(Z(I,1))*DALPHA**2-
     *                  d2lfun(t,Kc,alphatau,Nmodes)*COS(ALPHA))
         DPDX(I)=DPDX(I)+LAPL_U*REYNinv 
         DPDY(I)=DPDY(I)+LAPL_V*REYNinv 
C         write(90,22)teta(i),lapl_u,lapl_v,dpdx(i),dpdy(i)
 22   FORMAT(5(E13.7,1X))         
      ENDDO
      
      DPDX(NTETA+1)=DPDX(1)
      DPDY(NTETA+1)=DPDY(1)



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


C     CALCOLO DELLE FORZE
C
      F_X=0.D0
      F_Xf=0.D0
      F_Xa=0.D0
      F_Y=0.D0
      F_Yf=0.D0
      F_Ya=0.D0
      M_Z=0.D0
      F_Xnew=0.D0
      F_Ynew=0.D0

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
         D11N=(C2*CX1+0.5d0*S2*CX2)*(DP1)/JACOB(I,1)**2
         D22N=-D11N
         D12N=0.5d0*(-C2*CX2+2.d0*S2*CX1)*(DP1)/JACOB(I,1)**2

         F_Xf = F_Xf-PNEW(I)*DYDTETA(I)*dteta
         F_Xa = F_Xa+2.d0*(DYDTETA(I)*D11-DXDTETA(I)*D12)/REYN*DTETA         
         F_X = F_X+(-PNEW(I)*DYDTETA(I)+
     *               2.d0*(DYDTETA(I)*D11-DXDTETA(I)*D12)/REYN)*DTETA
         F_Yf = F_Yf+PNEW(I)*DXDTETA(I)*DTETA
         F_Ya = F_Ya+2.d0*(DYDTETA(I)*D12-DXDTETA(I)*D22)/REYN*DTETA
         F_Y = F_Y+(PNEW(I)*DXDTETA(I)+
     *               2.d0*(DYDTETA(I)*D12-DXDTETA(I)*D22)/REYN)*DTETA
         M_Z = M_Z+PNEW(I)*(DXDTETA(I)*DREAL(Z(I,1))+
     *                 DYDTETA(I)*DIMAG(Z(I,1)))*DTETA+
     *           (2.D0*DTETA/REYN)*(-DREAL(Z(I,1))*DXDTETA(I)*D22+
     *          (DREAL(Z(I,1))*DYDTETA(I)+DIMAG(Z(I,1))*DXDTETA(I))*D12-
     *            DIMAG(Z(I,1))*DYDTETA(I)*D11)
      ENDDO

      F_XHATf = F_Xf*COS(ALPHA)-F_Yf*SIN(ALPHA)
      F_YHATf = F_Xf*SIN(ALPHA)+F_Yf*COS(ALPHA)
      F_XHATa = F_Xa*COS(ALPHA)-F_Ya*SIN(ALPHA)
      F_YHATa = F_Xa*SIN(ALPHA)+F_Ya*COS(ALPHA)
      F_XHAT = F_X*COS(ALPHA)-F_Y*SIN(ALPHA)
      F_YHAT = F_X*SIN(ALPHA)+F_Y*COS(ALPHA) 
      F_Xnew= F_Xhat*cos(beta)+F_Yhat*sin(beta)
      F_Ynew=-F_Xhat*sin(beta)+F_Yhat*cos(beta)

      WRITE(indforce,25)T,F_X,F_Y,M_Z,F_XHAT,F_YHAT,
     &F_XHATa,F_YHATa,F_XHATf,F_YHATf,F_Xnew, F_Ynew
      close (indforce)
      open (indforce, status='old', position='append')
 24   FORMAT(6(E13.7,1X)) 
 25   FORMAT(12(E13.7,1X)) 
      
      RETURN
      END

C ***********************************************************************

      SUBROUTINE EFFIC
      INCLUDE 'sensvar.com'

      REAL*8 M_ZOLD,Mm

      ntot=25000      

      fm_x = 0.D0
      fm_xf = 0.D0
      fm_xa = 0.D0
      fm_y = 0.D0
      fm_yf = 0.D0
      fm_ya = 0.D0
      Mm=0.D0
      delta_t=0.D0
      POT=0.D0
      Potnew=0.d0
      POTY=0.D0
      Potynew=0.d0
      POTM=0.D0
      thrust=0.d0
      dlift=0.d0


 10   read(indforce,25)T_act,F_X,F_Y,M_Z,F_XHAT,
     $ F_YHAT,F_XHATa,F_YHATa,F_XHATf,F_YHATf,F_Xnew,F_Ynew
      if (t_act.lt.pi*1.99) goto 10

      t_old=t_act
      f_xold=f_xhat
      f_yold=f_yhat
      f_xoldf=f_xhatf
      f_yoldf=f_yhatf
      f_xolda=f_xhata
      f_yolda=f_yhata
      m_zold=m_z
      Fxbeta=F_Xnew
      Fybeta=F_Ynew
      dl_old=dlfun(t_old,Kc,alphatau,Nmodes)
      dalpha_old=dalphafun(t_old,alphaamp,alphaphi,Nmodes)
      
      do i=1,ntot
         

        read(indforce,25,end=20,err=20)T_act,F_X,F_Y,M_Z,
     $  F_XHAT,F_YHAT,F_XHATa,F_YHATa,F_XHATf,F_YHATf,F_Xnew,F_Ynew
         alpha=alphafun(t_act,alphaamp,alphaphi,Nmodes)
         
         dtloc=t_act-t_old
         dl=dlfun(t_act,Kc,alphatau,Nmodes)
         dalpha=dalphafun(t_act,alphaamp,alphaphi,Nmodes)
         fm_x  = fm_x+0.5d0*(f_xold+f_xhat)*dtloc
         fm_xf = fm_xf+0.5d0*(f_xoldf+f_xhatf)*dtloc
         fm_xa = fm_xa+0.5d0*(f_xolda+f_xhata)*dtloc
         fm_y  = fm_y+0.5d0*(f_yold+f_yhat)*dtloc
         fm_yf = fm_yf+0.5d0*(f_yoldf+f_yhatf)*dtloc
         fm_ya = fm_ya+0.5d0*(f_xolda+f_xhata)*dtloc
         thrust=thrust+0.5d0*(fxbeta+f_xnew)*dtloc
         dlift=dlift+0.5d0*(fybeta+f_ynew)*dtloc
         Mm = Mm+0.5d0*(M_zold+M_z)*dtloc
         delta_t = delta_t+dtloc

         POTY=POTY+0.5d0*(f_yold*dl_old+f_yhat*dl)*dtloc
         potynew=potynew+0.5d0*(fybeta*dl_old+f_ynew*dl)*dtloc
         pOTM=pOTM+0.5d0*(m_zold*dalpha_old+m_z*dalpha)*dtloc
         POT=POTY+POTM
         potnew=potynew+potM
         potnew=-potnew
         
         t_old=t_act
         f_xold=f_xHAT
         f_yold=f_yHAT
         f_xoldf=f_xHATf
         f_yoldf=f_yHATf
         f_xolda=f_xHATa
         f_yolda=f_yHATa
         Fxbeta=F_Xnew
         Fybeta=F_Ynew
         m_zold=m_z
         DL_OLD=DL
         dalpha_old=dalpha

         
      enddo
 20   continue
      

      fm_x = fm_x/delta_t
      fm_xf = fm_xf/delta_t
      fm_xa = fm_xa/delta_t
      fm_y = fm_y/delta_t
      fm_yf = fm_yf/delta_t
      fm_ya = fm_ya/delta_t
      thrust=thrust/delta_t
      dlift=dlift/delta_t
      Mm = Mm/delta_t
      pot=pot/delta_t
      potnew=potnew/delta_t
      potY=potY/delta_t
      potynew=potynew/delta_t
      potM=potM/delta_t
      
      cT = 0.5d+0*fm_x/(u_tilde*u_tilde)
      cTnew = 0.5d+0*thrust/(u_tilde*u_tilde)
      cTf = 0.5d+0*fm_xf/(u_tilde*u_tilde)
      cTa = 0.5d+0*fm_xa/(u_tilde*u_tilde)
      cL = 0.5d+0*fm_y/(u_tilde*u_tilde)
      Clnew=0.5d+0*dlift/(u_tilde*u_tilde)
      cLf = 0.5d+0*fm_yf/(u_tilde*u_tilde)
      cLa = 0.5d+0*fm_ya/(u_tilde*u_tilde)
      cM = Mm/(u_tilde*u_tilde*8.d0)
      cP = 0.5d+0*Pot/(u_tilde*u_tilde*u_tilde)
      cPnew = 0.5d+0*Potnew/(u_tilde*u_tilde*u_tilde)
      eff = fm_x*u_tilde/pot
      effinew=-(thrust*U_tilde)/potnew


      WRITE(indpf,*)' DELTA_T = ',DELTA_T
      write(indpf,*)' fm_x, fm_xf, fm_xa ',fm_x,fm_xf,fm_xa
      write(indpf,*)' cT, cTf, cTa ',cT,cTf,cTa
      write(indpf,*)' fm_y, fm_yf, fm_ya ',fm_y,fm_yf,fm_ya
      write(indpf,*)' cL, cLf, cLa ',cL,cLf,cLa
      WRITE(indpf,*)' Mm, cM', Mm, cM
      WRITE(indpf,*)' pot, potY, potM',pot, potY, potM
      WRITE(indpf,*)' cP,  eff',  cP,  eff
      WRITE(indpf,*)' beta en degres',  beta*180/pi
      WRITE(indpf,*)' *********************************************'
      WRITE(indpf,*)' *********************************************'
      
      write(indpf,*) 'new thurst, new lift', thrust, dlift
      write(indpf,*) 'new CT, new Cl', CTnew, CLnew
      write(indpf,*) 'new CP, new CM', CPnew, CM
      write(indpf,*) 'newpot , new effic', potnew, effinew
 
 23   FORMAT(4(E13.7,1X))
 24   FORMAT(6(E13.7,1X))
 25   FORMAT(12(E13.7,1X))
     
      return
      end

C ***********************************************************************

      real*8 function lfun(t,amp,tau,Nmodes)
      real*8 t
      real*8 amp(Nmodes),tau(Nmodes)
     
      aux=0
      do k=1,Nmodes
      aux=aux+amp(k)*dsin(k*t+tau(k))
      enddo
      lfun=aux
      return
      end

C ***********************************************************************

      real*8 function dlfun(t,amp,tau,Nmodes)
      real*8 t,amp(Nmodes),tau(Nmodes)

     
      aux=0   
      do k=1,Nmodes
      aux=aux+k*amp(k)*dcos(k*t+tau(k))
      enddo 
      dlfun=aux
      return
      end


C ***********************************************************************

      real*8 function d2lfun(t,amp,tau,Nmodes)
      real*8 t,amp(Nmodes),tau(Nmodes)
 

      aux =0
      do k=1,Nmodes
      aux=aux-k**2*amp(k)*dsin(k*t+tau(k))
      enddo
      d2lfun=aux
      return
      end

C ***********************************************************************


  
      real*8 function alphafun(t,amp,fi,Nmodes)
      real*8 t
      real*8 amp(Nmodes), fi(Nmodes)
      integer*8 Nmodes

      aux=0
      do k=1,Nmodes
      aux=aux+amp(k)*dsin(k*t+fi(k))
      enddo
      alphafun=aux
      return
      end

C ***********************************************************************

      real*8 function dalphafun(t,amp,fi,Nmodes)
      real*8 t
      real*8 amp(Nmodes), fi(Nmodes)
      integer*8 Nmodes

      aux=0
      do k=1,Nmodes
      aux=aux+k*amp(k)*dcos(k*t+fi(k))
      enddo
      dalphafun=aux
      return
      end

C ***********************************************************************

      real*8 function d2alphafun(t,amp,fi,Nmodes)
      real*8 t
      real*8 amp(Nmodes),fi(Nmodes)
      integer*8 Nmodes
     
      aux=0
      do k=1,Nmodes
      aux=aux-k**2*amp(k)*dsin(k*t+fi(k))
      enddo
      d2alphafun=aux
      return
      end




 
C **********************************************************************  
      
      SUBROUTINE TRIPER(AA,BB,CC,FF,N,GAM2)
C      INCLUDE 'xcilgetto.com'
      IMPLICIT REAL*8 (A-H,O-Z)
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
      DIMENSION AA(1),BB(1),CC(1),FF(1),GAM2(1)

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
C      INCLUDE 'xcilgetto.com'
      IMPLICIT REAL*8 (A-H,O-Z)
C
C ****************************************************************************
C     VEDI ROACHE APP.A
C     risolve il sistema tridiagonale per condizione di Dirichlet alla parete, 
C     dove e' nota la omega W(1) e condizione di Neumann all'infinit6o
C     
C      -A(K)*X(K+1)+B(K)*X(K)-C(K)*X(K-1)=D(K)
C ****************************************************************************
C
      PARAMETER (IZETA=800)

      DIMENSION A(1),B(1),C(1),D(1),W(1),E(IZETA),F(IZETA)

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
C      INCLUDE 'xcilgetto.com'
C ****************************************************************************
C     VEDI ROACHE APP.A
C     risolve il sistema tridiagonale per condizione di Dirichlet 
C     alla parete e all'infinito
C ****************************************************************************
      
      PARAMETER (IZETA=800)
      
      IMPLICIT REAL*8 (A-H,O-Z)
 
      COMPLEX*16 W,D,A1,AM,E,F,DEN
      
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



C *************************************************************************

      SUBROUTINE OUTPUT
      
      INCLUDE 'sensvar.com'

      write(2)t

C      DO 10 I=1,NTETA+1
      DO 10 I=1,NTETA
         DO 20 J=1,NZETA
            WRITE(2)OMG(I,J,1),PSI(I,J),omgga(i,j,1,indmod),
     $           psiga(i,j,indmod)
 20      CONTINUE
 10   CONTINUE

      RETURN
      END
C *************************************************************************

     
    
  
      SUBROUTINE cfft1(y,nx,iopt)

      IMPLICIT REAL*8 (A-H,O-Z)

C     complex vaux

      COMPLEX*16 Y
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
      IMPLICIT REAL*8 (A-H,O-Z)      
      
      DIMENSION M(1024),X(1)
      
      COMPLEX*16 X,WK,HOLD,Q

C     IF(N.GT.8)WRITE(6,*)' ATT. N NELLA FFT GT DEL TOLLERATO ! '
      
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
C              V=SIGN*6.2831853*FK/FLX
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
      CHARACTER*64 NAME

      IOS=-1
      DO WHILE (IOS.NE.0)
         WRITE(6,*) MESS
         READ(5,'(A)') NAME
C         NAME='st425.dat'
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
C         NAME='out'
         IF (NAME.EQ.'/E') STOP
         OPEN(UNIT=LU,FILE=NAME,STATUS=STATUS,IOSTAT=IOS,
     *                           FORM='UNFORMATTED')
         IF (IOS.NE.0) THEN
            WRITE(6,*)'Errore ',IOS,' File : ',NAME
         ENDIF
      ENDDO

      END

C ****************************************************************************

