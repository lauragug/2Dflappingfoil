c     Common per il programma YOUKOWSKI
      
      parameter (iteta=512,izeta=500,idim=1)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      
      integer*8 ntciclo
      

      real*8 m_z
      real*8 Kc

      real*8 Jacob,Jacobsq

      real*8 l,lfun
      
      complex*16 c0,ci

      complex*16 chi,z,coord
      
      complex*16 vet_z
      complex*16 matomg,matomgga
      complex*16 matpsi,matpsiga

      common/costcom/    
     *       pi,             ! p greco
     *       c0,             ! zero complesso
     *       ci              ! i complesso
     
      
      common/run/
     *       iopt,           ! 0 - se è il primo run; 1 - in caso contrario
     *       icil,           ! 0 - se profilo joukowski; 1 - se cilindro
     *	     indforce,       ! index for the forces files
     *	     indfa,          ! index for ak for forces files
     *	     indea,          ! index for ak for effic files
     *	     indpf          ! index for proefficienza files




   
      common/subinit/
     *       re,            !    
     *       Kc(idim),            ! Heaving oscillations amplitude    
     *       reyn           ! numero di reynold del moto


      common/geomprofilo/
     *       dis,            ! distanza c. di rotazione-apice del profilo
     *       e,              ! parametro e del prof Youko (e/lambda)
     *       small           ! quantita' piccola

      common/reticolo/
     *       rmax,           ! rapporto tra raggio estremo e lambda
     *       rint,           ! raggio entro il quale si vogliono metà punti 
     *       a,              ! parametro a della trasformazione asse r
     *       nteta,          ! numero punti lungo teta
     *       nzeta,          ! numero punti lungo z
     *       dteta,          ! dimensione cella in teta
     *       dzeta           ! dimensione cella in zeta

      common/reticolodimension/
     *       teta(iteta+1),        ! vettore coordinata teta
     *       zeta(izeta),        ! vettore coordinata zeta
     *       r(izeta),           ! vettore coordinata r adimens.
     *       r2(izeta),          ! vettore coordinata r adimens. al quadrato
     *       ez(izeta),          ! exp coordinata adimens. verticale
     *       ez2(izeta),         ! exp**2 coordinata adimens. verticale
     *       chi(iteta+1,izeta),   ! coordinata nel piano csi,eta
     *       z(iteta+1,izeta),     ! coordinata adimens. nel piano x,y
     *       coord(iteta+1,izeta), ! coordinata adimens. nel piano x,y cappello fisso
     *       Jacob(iteta+1,izeta),     ! Jacobiano nel piano x,y
     *       Jacobsq(iteta+1,izeta),     ! Radice dello Jacobiano nel piano x,y
     *       dxdcsi(iteta+1,izeta),     ! derivata di X rispetto a csi nel piano x,y
     *       dxdeta(iteta+1,izeta),     ! derivata di X rispetto a eta nel piano x,y
     *       dydcsi(iteta+1,izeta),     ! derivata di Y rispetto a csi nel piano x,y
     *       dydeta(iteta+1,izeta),     ! derivata di Y rispetto a eta nel piano x,y
     *       dxdteta(iteta+1),         ! derivata di X rispetto a teta nel piano X,Y
     *       dydteta(iteta+1),         ! derivata di Y rispetto a teta nel piano X,Y
     *       d2xdcsi2(iteta+1),         ! derivata sec. di X rispetto a csi
     *       d2xdcsieta(iteta+1),      ! derivata sec. di X rispetto a csi e eta
     *       d2xdcsiteta(iteta+1),     ! derivata di X rispetto a csi e teta lungo il profilo
     *       d2xdetateta(iteta+1),     ! derivata di X rispetto a csi e teta lungo il profilo
     *       d2ydcsiteta(iteta+1),     ! derivata di X rispetto a csi e teta lungo il profilo
     *       d2ydetateta(iteta+1)     ! derivata di X rispetto a csi e teta lungo il profilo

      common/moto/
     *       u_tilde,       ! modulo della velocità del flusso
     *       freduced,      ! reduced frequency  substitutes the strouhal number
     *       Strouhal,      ! The strouhal number
     *       beta,          ! angolo di incidenza del flusso  
     *       alphacost,     ! angolo alpha costante
     *       alphatau(idim),      ! ampiezza della rotazione oscillante
     *       alphaamp(idim),      ! ampiezza della rotazione oscillante
     *       alphaphi(idim),      ! ampiezza della rotazione oscillante
     *       l,             ! spostamento del moto oscillante
     *       dl,            ! derivata spostamento del moto oscillante
     *       d2l,           ! derivata seconda spostamento
     *       alpha,         ! rotazione
     *       dalpha,        ! derivata della rotazione
     *       d2alpha,       ! derivata seconda rotazione
     *       Nmodes         ! Number of modes in heaving and pitching motions 


 
      common/tempo/	
     *       t,                  ! variabile tempo
     *       dt,                 ! passo temporale
     *       tfin,               ! istante finale
     *       ntciclo             ! numero passi in un ciclo

      common/stamparis/
     *       nusc1,              ! numero steps temporali uscite forze
     *       iusc1,              ! indice su nusc1
     *       itout,              ! 
     *       itusc,              ! intervallo stampa video tempo 
     *       indgdt,
     *       indmod              ! Modes index

      common/grad/
     *       weightmu,
     *       weightnu,
     *       weightgamma,
     *       weightdelta,
     *       gradient(idim),
     *       gradientak(idim),
     *       costfct,
     *       dlift,
     *       thrust,
     *       CTnew,
     *       CLnew,
     *       CPnew,
     *       effinew,	
     *       savea,
     *       savegdta



      common/soluzione/
     *	     v_r(iteta+1,izeta),   	      	! componente di velocita' v_r
     *       v_teta(iteta+1,izeta),	      	! componente di velocita' v_teta 
     *       omg(0:iteta+1,izeta,2), 		! vorticita'agli istanti n e n+1 
     *       omgga(0:iteta+1,izeta,2,idim), 	! sensitivity of the vorticity 
     *       psi(0:iteta+1,izeta),   		! funzione di corrente
     *       psiga(0:iteta+1,izeta,idim)   	! sensitivity of the stream function



      common/matrices/
     *       matomg(iteta,izeta),   		! trasformate di fourier di omega
     *       matomgga(iteta,izeta,idim),       	! FT of the vorticity sensitivity
     *       matpsi(iteta,izeta),         	! trasformate di fourier di psi
     *       matpsiga(iteta,izeta,idim)       	! FT of the stream function sensitivity 



      common/risultanti/
     *       dpdx(iteta+1),             ! derivata della pressione lungo X
     *       dpdy(iteta+1),             ! derivata della pressione lungo Y
     *       p(iteta+1),                ! pressione sul profilo
     *       dpdtetadak(iteta+1),
     *	     dpdwdak(iteta+1),
     *       dpdwdakdteta(iteta+1),
     *       p_s(iteta/2+1),            ! pressione sul profilo superiore
     *       p_i(iteta/2+1),            ! pressione sul profilo inferiore
     *       pnew(iteta+1),             ! pressione sul profilo integrando da 0
     *       dpdak(iteta+1),   
     *       f_x,			! componente risultante su X
     *       f_xf,			! coeff resistenza di forma
     *       f_xa,			! coeff resistenza di attrito
     *       f_y,			! componente risultante su Y
     *       m_z,                       ! momento intorno all'asse z
     *       f_xnew,                    ! componente risultante su X
     *       f_ynew,                    ! componente risultante su Y
     *       f_xhat,			! componente risultante su x fisso
     *       f_yhat			! componente risultante su y fisso
