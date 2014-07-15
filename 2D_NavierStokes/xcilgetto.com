c     Common per il programma YOUKOWSKI
      
      parameter (iteta=1024,izeta=600,idim=20,ixhat=500,iyhat=500)
      
      IMPLICIT REAL(8) (A-H,O-Z)

      character(12)nomebin      
      integer(8) ntciclo
      
      real(8) Kc,m_z
      real(8) Jacob,Jacobsq

      real(8) l,lfun
      
      complex(16) c0,ci

      complex(16) chi,z,coord
      
      complex(16) vet_z
      complex(16) matomg,matpsi

      common/costcom/    
     *       pi             ! p greco

      common/costcmplx/  
     *       c0,             ! zero complesso
     *       ci              ! i complesso
      
      common/run/
     *       iopt,           ! 0 - primo run; 1 - in caso contrario
     *       icil            ! 0 - joukowski; 1 - cilindro
   
      common/subinit/
     *       re,             !    
     *       reyn,           ! numero di reynold del moto
     *       kc              ! numero di keulegan-carperter

      common/geomprofilo/
     *       dis,            ! distanza c. di rot-apice del profilo
     *       e,              ! parametro e del prof Youko (e/lambda)
     *       small           ! quantita' piccola

      common/reticolo/
     *       rmax,           ! rapporto tra raggio estremo e lambda
     *       rint,           ! raggio entro il quale metà punti 
     *       a,              ! parametro a della trasformazione asse r
     *       nteta,          ! numero punti lungo teta
     *       nzeta,          ! numero punti lungo z
     *       dteta,          ! dimensione cella in teta
     *       dzeta           ! dimensione cella in zeta

      common/reticolodimension/
     *       teta(iteta+1),        ! vettore coordinata teta
     *       zeta(izeta),        ! vettore coordinata zeta
     *       r(izeta),           ! vettore coordinata r adimens.
     *       r2(izeta),          ! v coordinata r adimens. al quadrato
     *       ez(izeta),          ! exp coordinata adimens. verticale
     *       ez2(izeta),         ! exp**2 coordinata adimens. vert
     *       Jacob(iteta+1,izeta),     ! Jacobiano nel piano x,y
     *       Jacobsq(iteta+1,izeta),     ! Radice dello Jacobiano x,y
     *       dxdcsi(iteta+1,izeta)     ! d di X rispetto a csi x,y

      common/reticolocmplx/
     *       chi(iteta+1,izeta),   ! coordinata nel piano csi,eta
     *       z(iteta+1,izeta),     ! coordinata adimens. nel piano x,y
     *       coord(iteta+1,izeta) ! c adimens. nel piano x,y fisso

      common/retidimension/
     *       dxdeta(iteta+1,izeta),     ! d di X rispetto a eta x,y
     *       dydcsi(iteta+1,izeta),     ! d di Y rispetto a csi x,y
     *       dydeta(iteta+1,izeta),     ! d di Y rispetto a eta x,y
     *       dxdteta(iteta+1),         ! d di X rispetto a teta x,y
     *       dydteta(iteta+1),         ! d di Y rispetto a teta x,y
     *       d2xdcsi2(iteta+1),         ! da sec. di X rispetto a csi
     *       d2xdcsieta(iteta+1),      ! d sec. di X ris a csi e eta
     *       d2xdcsiteta(iteta+1),     ! d di X r a csi e teta l prof
     *       d2xdetateta(iteta+1),     ! d di X r a csi e teta l prof
     *       d2ydcsiteta(iteta+1),     ! d di X r a csi e teta l prof
     *       d2ydetateta(iteta+1)     ! d di X ris a csi e teta l prof

      common/moto/
     *       u_tilde,       ! modulo della velocità del flusso
     *       strouhal,      ! numero di strouhal
     *       beta,          ! angolo di incidenza del flusso  
     *       alphaamp,      ! ampiezza della rotazione oscillante
     *       alphaphi,      ! ampiezza della rotazione oscillante
     *       l,             ! spostamento del moto oscillante
     *       dl,            ! derivata spostamento del moto oscillante
     *       d2l,           ! derivata seconda spostamento
     *       alpha,         ! rotazione
     *       dalpha,        ! derivata della rotazione
     *       d2alpha        ! derivata seconda rotazione 

      common/gettoint/
     *       icmedio,            ! 1 per il cal del campo medio, 0 no
     *       nxhat,              !
     *       nyhat               !

      common/getto/
     *       xhat_min,           !
     *       xhat_max,           !
     *       yhat_min,           !
     *       yhat_max,           !
     *       dt_getto,           !
     *       u_gr(2,2),          ! comp. u nel sist.rif. (X,Y)
     *       v_gr(2,2),          ! comp. v nel sist.rif. (X,Y)
     *       xhat(ixhat),        !
     *       yhat(iyhat),        ! 
     *       u_hat(ixhat,iyhat), !
     *       v_hat(ixhat,iyhat)  !

      common/tempo/
     *       t,                  ! variabile tempo
     *       dt,                 ! passo temporale
     *       tfin,               ! istante finale
     *       ntciclo             ! numero passi in un ciclo

      common/stamparisint/
     *       istampavel,          ! se 0 non stampa campi vel, 1 stampa
     *       ipuntivel,          ! se 1 da vel nel tempo pti, 0 no 
     *       nusc1,              ! numero steps temporali uscite
     *       nusc2,              ! numero steps temporali uscite
     *       nfusc2,             ! numero file uscita 2
     *       iusc1,              ! indice su nusc1
     *       iusc2,              ! indice su nusc2	
     *       itout,              ! 
     *       itusc,              ! intervallo stampa video tempo 
     *       numj,               ! n punti zeta su cui uscite
     *       num_po,             ! n istanti t in cui uscite campi
     *       num_pres,           ! n istanti t in cui uscite pressioni
     *       num_punti,          ! n punti in cui stampa vel nel tempo
     *       nind_po,            ! indice stampa campi
     *       nind_pres          ! indice stampa pressioni 

      common/stamparis/
     *       indj(idim),         ! vettore punti uscite lungo zeta
     *       vet_po(idim),       ! vettore istanti uscite campo
     *       vet_pres(idim),      ! vettore istanti uscite campo
     *       vet_z(idim)         ! vettore punti in cui stampa vel

      common/soluzione/
     *       v_r(iteta+1,izeta),        ! c di velocita' v_r
     *       v_teta(iteta+1,izeta),     ! c di velocita' v_teta 
     *       omg(0:iteta+1,izeta,2),    ! vort agli istanti n e n+1 
     *       psi(0:iteta+1,izeta),      ! funzione di corrente
     *       matomg(iteta,izeta),       ! trasf di fourier di omega
     *       matpsi(iteta,izeta)        ! trasf di fourier di psi 

      common/risultanti/
     *       dpdx(iteta+1),             ! d della pressione lungo X
     *       dpdy(iteta+1),             ! d della pressione lungo Y
     *       dpdteta(iteta+1),          ! d della pressione lungo teta
     *       p(iteta+1),                ! pressione sul profilo
     *       p_s(iteta/2+1),            ! p sul profilo superiore
     *       p_i(iteta/2+1),            ! p sul profilo inferiore
     *       pnew(iteta+1),             ! p sul profilo integrando da 0
     *       pod(iteta+1),              ! p sul profilo integr da pi/2
     *       pff(iteta+1)              ! p sul profilo tramite fft

      common/risulta/
     *       f_x,                       ! componente risultante su X
     *       f_xf,                      ! coeff resistenza di forma
     *       f_xa,                      ! coeff resistenza di attrito
     *       f_y,                       ! componente risultante su Y
     *       f_yf,                      ! coeff resistenza di forma
     *       f_ya,                      ! coeff resistenza di attrito
     *       m_z,                       ! momento intorno all'asse z
     *       f_xhat,                    ! c risultante su x fisso
     *       f_yhat,                    ! c risultante su y fisso
     *       f_xhatf,                   ! c risultante su x fisso
     *       f_yhatf,                   ! c risultante su y fisso
     *       f_xhata,                   ! c risultante su x fisso
     *       f_yhata                    ! c risultante su y fisso
