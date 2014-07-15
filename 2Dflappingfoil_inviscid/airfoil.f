      program airfoil

      include 'airfoil_com.f'
      
      call input

      indder=0

      tcorr = .000d+0
      
      u0P=0.d+0
      call alfacorr
      call alfapcorr
      call alfappcorr
      call rlcorr
      call rlpcorr
      call rlppcorr
     
      call init

c per fare queste uscite sarebbe necessario aggiornare il tempo    
c      call output(tcorr)
      call output1(tcorr)

      ncount = 1

c  il tempo qui di seguito aggiornato e' quello iniziale

      tcorr=tcorr+dt

      u0P=0.d+0
      call alfacorr
      call alfapcorr
      call alfappcorr
      call rlcorr
      call rlpcorr
      call rlppcorr
     
      do while (tcorr.lt.tfin)

c la soluzione iniziale viene aggiornata

      do n=-2,nvort
        vort_zita_v(n)=vort_zita(n)
        vort_eta_v(n)=vort_eta(n)
        vort_gamma_v(n)=vort_gamma(n)
      enddo

c     write(6,*)' tcorr = ',tcorr

C 1 passo del metodo di Runge-Kutta
c vortici liberi: il risultato viene messo in vort_zita_n e vort_eta_n

      do i=1,nvort
        call vort_mov1(i)
      enddo

c vortici emessi dalle singolarita'
c nella subroutine il tempo viene aggiornato

      call new_vort1

c il tempo viene aggiornato. La posizione dei vortici
c e' gia' stata aggiornata

      tcorr=tcorr+dt

      u0P=0.d+0
      call alfacorr
      call alfapcorr
      call alfappcorr
      call rlcorr
      call rlpcorr
      call rlppcorr

c     write(6,*)' primo passo '
c     write(6,*)' 1 ',tcorr,vort_zita(-1),vort_eta(-1),vort_gamma(-1)
c     write(6,*)' 2 ',tcorr,vort_zita(-2),vort_eta(-2),vort_gamma(-2)

C 2 passo del metodo di Runge-Kutta
c vortici liberi

      do i=1,nvort
        call vort_mov2(i)
      enddo

c vortici emessi dalle singolarita'
c nella subroutine il tempo viene aggiornato

      call new_vort2

c il tempo viene aggiornato. La posizione dei vortici
c e' gia' stata aggiornata

      tcorr=tcorr+dt
      
      u0P=0.d+0
      call alfacorr
      call alfapcorr
      call alfappcorr
      call rlcorr
      call rlpcorr
      call rlppcorr
      
c     write(6,*)' secondo passo '
c     write(6,*)' 1 ',tcorr,vort_zita(-1),vort_eta(-1),vort_gamma(-1)
c     write(6,*)' 2 ',tcorr,vort_zita(-2),vort_eta(-2),vort_gamma(-2)

c      call force
      
      if (ncount.eq.ndtusc) then
c        call output(tcorr)
        call output1(tcorr)
        ncount = 1
      else
        ncount = ncount + 1
      endif

      enddo !while

      write(6,*)'rend = ',rendn/rendd,' util , diss = ',rendn,rendd

      end

C ***********************************************************************************
C ***********************************************************************************
C ***
C ***        subroutine input
C ***
C ***     lettura parametri da input
C ***
C ***********************************************************************************
C ***********************************************************************************

      subroutine input

      include 'airfoil_com.f'

      open(20,file='dati7',status='old')
      read(20,*)dt,tfin,ndtusc
      write(6,*)' dt tfin  ndtusc =',dt,tfin,ndtusc
      read(20,*)lambda,dp
      write(6,*)' lambda,dp =',lambda,dp
      read(20,*)u0,beta
      write(6,*)' u0,beta =',u0,beta
      read(20,*)rlh,omega
      write(6,*)' rlh,omega =',rlh,omega
      read(20,*)alfah,phi_a
      write(6,*)' alfah,phi_a =',alfah,phi_a
      read(20,*)perc,dec,enne
      write(6,*)' perc,dec,enne =',perc,dec,enne
      read(20,*)t1,t2
      write(6,*)' t1 ,t2 =',t1,t2
      
      close(20)

      return
      end

C ***********************************************************************************
C ***********************************************************************************
C ***
C ***        subroutine init
C ***
C ***     inizializzaione variabili e vortici emessi dalla piastra
C ***
C ***********************************************************************************
C ***********************************************************************************

      subroutine init

      include 'airfoil_com.f'

      real*8 tetac

c  inizializzazione costanti

      pi = 2.d+0*asin(1.d+0)
      ci = (0.d+0,1.d+0)
      cg = (pi/4.d+0)**(2.d+0/3.d+0)
      rho=1000.d+0
      gamma_max = 0.d+0

c inizializzazione vettori e quantità

      do n=-2,nv
	vort_zita(n)=0.d+0
	vort_eta(n)=0.d+0
	vort_gamma(n)=0.d+0
	vort_zita_n(n)=0.d+0
	vort_eta_n(n)=0.d+0
	vort_gamma_n(n)=0.d+0
      enddo !n

      nvort=0

c   inizializzazione vortici emessi dal bordo della piastra

      call init_vort(1,1)

      do i=-2,-1
        vort_zita(i)=vort_zita_n(i)
	vort_eta(i)=vort_eta_n(i)
	vort_gamma(i)=vort_gamma_n(i)
      enddo

	rendn=0.d+0
	rendd=0.d+0
	iren=0

      return 
      end

C ***********************************************************************************
C ***********************************************************************************
C ***
C ***        subroutine new_vort1
C ***
C ***     calcola la posizione dei vortici attaccati alla piastra ed 
C ***       eventualmente li stacca
C *** 
C ***********************************************************************************
C ***********************************************************************************

      subroutine new_vort1

      
      include 'airfoil_com.f'

      real*8 rJg,Xg,Yg,dxgdzita,dygdzita,dxgdeta,dygdeta
      real*8 dpsi_t_dzita,dpsi_t_deta,uvor,vvor

      real*8 zitac,etac,gammac,dzita,deta
      
      real*8 zita_n,eta_n,gamma_n

      dimension eta_n(-2:-1),zita_n(-2:-1),gamma_n(-2:-1),ind(-2:-1)

C conti pagina b quaderno arancio a meta
      
      call coeff   
      
c     write(6,*)' ************* '
c     
c     write(6,*)'1 dzit1,2,dt = ',dzit1dt,dzit2dt
c     write(6,*)'1 det1,2,dt = ',det1dt,det2dt
c     write(6,*)'1 dgam1,2,dt = ',dgam1dt,dgam2dt
c   
      vort_zita(-2) = vort_zita_v(-2)+dt*dzit2dt
      vort_zita(-1) = vort_zita_v(-1)+dt*dzit1dt
      vort_eta(-2) = vort_eta_v(-2)+dt*det2dt
      vort_eta(-1) = vort_eta_v(-1)+dt*det1dt
      vort_gamma(-2) = vort_gamma_v(-2)+dt*dgam2dt 
      vort_gamma(-1) = vort_gamma_v(-1)+dt*dgam1dt 

C aggiorno la posizione anche degli altri vortici

      do i=1,nvort
        vort_zita(i)=vort_zita_n(i)
        vort_eta(i)=vort_eta_n(i)
      enddo

      return 
      end

C ***********************************************************************************
C ***********************************************************************************
C ***********************************************************************************
C ***
C ***        subroutine new_vort2
C ***
C ***     calcola la posizione dei vortici attaccati alla piastra ed 
C ***       eventualmente li stacca
C *** 
C ***********************************************************************************
C ***********************************************************************************

      subroutine new_vort2

      include 'airfoil_com.f'

      real*8 rJg,Xg,Yg,dxgdzita,dygdzita,dxgdeta,dygdeta
      real*8 dpsi_t_dzita,dpsi_t_deta,uvor,vvor

      real*8 raux
    
      real*8 zitac,etac,gammac,dzita,deta
      
      real*8 zita_n,eta_n,gamma_n
      real*8 dgam

      dimension eta_n(-2:-1),zita_n(-2:-1),gamma_n(-2:-1),ind(-2:-1)
      dimension dgam(-2:-1)

C conti pagina b quaderno arancio a meta
      
      do ii=-2,-1
        if (abs(vort_gamma_v(ii)).gt.gamma_max) 
     .                               gamma_max=abs(vort_gamma_v(ii))
      enddo
      
      call coeff   
   
c     write(6,*)' ************* '
      
c     write(6,*)'2 dzit1,2,dt = ',dzit1dt,dzit2dt
c     write(6,*)'2 det1,2,dt = ',det1dt,det2dt
c     write(6,*)'2 dgam1,2,dt = ',dgam1dt,dgam2dt
    
      zita_n(-2) = vort_zita_v(-2)+2.d+0*dt*dzit2dt
      zita_n(-1) = vort_zita_v(-1)+2.d+0*dt*dzit1dt
      eta_n(-2) = vort_eta_v(-2)+2.d+0*dt*det2dt
      eta_n(-1) = vort_eta_v(-1)+2.d+0*dt*det1dt
      gamma_n(-2) = vort_gamma_v(-2)+2.d+0*dt*dgam2dt 
      gamma_n(-1) = vort_gamma_v(-1)+2.d+0*dt*dgam1dt 
      
      dgam(-2)=dgam2dt
      dgam(-1)=dgam1dt

      raux=( vort_zita_v(-2)**2+vort_eta_v(-2)**2 )**2
      dZncapdt(-2)= (
     .   ( raux-lambda**2*(vort_zita_v(-2)**2-vort_eta_v(-2)**2) )/raux
     .  +ci*2.d+0*lambda**2*vort_zita_v(-2)*vort_eta_v(-2)/raux
     .              )*( dzit2dt+ci*det2dt )
      raux=( vort_zita_v(-1)**2+vort_eta_v(-1)**2 )**2
      dZncapdt(-1)= (
     .   ( raux-lambda**2*(vort_zita_v(-1)**2-vort_eta_v(-1)**2) )/raux
     .  +ci*2.d+0*lambda**2*vort_zita_v(-1)*vort_eta_v(-1)/raux
     .              )*( dzit1dt+ci*det1dt )

C aggiorno la posizione degli altri vortici

      do i=1,nvort
        vort_zita(i)=vort_zita_n(i)
        vort_eta(i)=vort_eta_n(i)
      enddo

C assegnazioni

      do ii=-2,-1
        ind(ii)=0
      enddo

      do ii = -2,-1,1

        if ( gamma_n(ii)*dgam(ii).lt.0..
     .           and.abs(dgam(ii)).gt.perc) then
       
C stacco                      

          ind(ii)=1
          nvort = nvort + 1
          vort_zita(nvort)  = zita_n(ii)
          vort_eta(nvort)   = eta_n(ii)
          vort_gamma(nvort) = vort_gamma_v(ii)
	  t_cut(nvort)=tcorr 
	  vort_gamma(ii)=0.d+0
	  indcut(nvort)=ii
          dZncapdt(nvort)=(0.d+0,0.d+0)
        else

C il vortice sta crescendo

          vort_zita(ii) = zita_n(ii)
          vort_eta(ii) = eta_n(ii)
	  vort_gamma(ii) = gamma_n(ii) 
        
	endif
      enddo  !ii
      
      call init_vort(ind(-2),ind(-1))

      do ii=-2,-1
        if (ind(ii).eq.1) then 
          vort_eta(ii)=vort_eta_n(ii)
	  vort_zita(ii)=vort_zita_n(ii)
	  vort_gamma(ii)=vort_gamma_n(ii)
        endif
      enddo
      
      return 
      end

C ***********************************************************************************
C ***********************************************************************************
C ***
C ***        subroutine vort_mov1
C ***
C ***   conti a pag 18 - arancio con d_gamma/dt=0.
C *** 
C ***********************************************************************************
C ***********************************************************************************

      subroutine vort_mov1(n)

      include 'airfoil_com.f'

      real*8 rJg,dxgdzita,dygdzita,dxgdeta,dygdeta
      real*8 dpsi_t_dzita,dpsi_t_deta,uvor,vvor
      
      real*8 zitac,etac,dzita,deta

      real*8 raux,tetaaux

      zitac = vort_zita(n)
      etac = vort_eta(n)

      dzita = 
     .      dpsi_t_deta(zitac,etac,n)  
     .   +  (    uvor(n)*dygdeta(zitac,etac) -
     .           vvor(n) *dxgdeta(zitac,etac)    )

      deta = 
     .      - dpsi_t_dzita(zitac,etac,n) + 
     .      ( - uvor(n)*dygdzita(zitac,etac) +
     .          vvor(n)*dxgdzita(zitac,etac)  )  

      dzita = dzita/rJg(zitac,etac)
      deta = deta/rJg(zitac,etac)

      vort_zita_n(n) = vort_zita(n)+dt*dzita 
      vort_eta_n(n) = vort_eta(n)+dt*deta   

      raux=sqrt(vort_zita_n(n)**2.d+0+vort_eta_n(n)**2.d+0)
      if(raux.lt.lambda)then
        if(vort_zita_n(n).eq.0.d+0)then
          tetaaux=pi/2.d+0
	  if(vort_eta_n(n).lt.0.d+0)tetaaux=-pi/2.d+0
	endif
        tetaaux=atan(vort_eta_n(n)/vort_zita_n(n))
        if(vort_zita_n(n).lt.0.d+0)tetaaux=tetaaux+pi
        raux=2.d+0*lambda-raux
        vort_zita_n(n)=raux*cos(tetaaux)
        vort_eta_n(n)=raux*sin(tetaaux)
      endif

      return
      end

C ***********************************************************************************
C ***********************************************************************************
C ***********************************************************************************
C ***
C ***        subroutine vort_mov2
C ***
C ***   conti a pag 18 - arancio con d_gamma/dt=0.
C *** 
C ***********************************************************************************2C ***********************************************************************************

      subroutine vort_mov2(n)

      include 'airfoil_com.f'

      real*8 rJg,dxgdzita,dygdzita,dxgdeta,dygdeta
      real*8 dpsi_t_dzita,dpsi_t_deta,uvor,vvor
      
      real*8 zitac,etac,dzita,deta

      real*8 raux,tetaaux

      if (n.lt.1.or.n.gt.nvort) then
        write(6,*)' errore indice vort_mov',n
        stop
      endif

      zitac = vort_zita(n)
      etac = vort_eta(n)

      dzita = 
     .      dpsi_t_deta(zitac,etac,n)  
     .   +  (    uvor(n)*dygdeta(zitac,etac) -
     .           vvor(n) *dxgdeta(zitac,etac)    )

      deta = 
     .      - dpsi_t_dzita(zitac,etac,n) + 
     .      ( - uvor(n)*dygdzita(zitac,etac) +
     .          vvor(n)*dxgdzita(zitac,etac)  )  

      dzita = dzita/rJg(zitac,etac)
      deta = deta/rJg(zitac,etac)

      raux=( vort_zita_v(n)**2+vort_eta_v(n)**2 )**2
      dZncapdt(n)= (
     .   ( raux-lambda**2*(vort_zita_v(n)**2-vort_eta_v(n)**2) )/raux
     .  +ci*2.d+0*lambda**2*vort_zita_v(n)*vort_eta_v(n)/raux
     .              )*( dzita+ci*deta )

      vort_zita_n(n)=vort_zita_v(n)+2.d+0*dt*dzita 
      vort_eta_n(n)=vort_eta_v(n)+2.d+0*dt*deta   

      raux=sqrt(vort_zita_n(n)**2.d+0+vort_eta_n(n)**2.d+0)
      if(raux.lt.lambda)then
        if(vort_zita_n(n).eq.0.d+0)then
          tetaaux=pi/2.d+0
	  if(vort_eta_n(n).lt.0.d+0)tetaaux=-pi/2.d+0
	endif
        tetaaux=atan(vort_eta_n(n)/vort_zita_n(n))
        if(vort_zita_n(n).lt.0.d+0)tetaaux=tetaaux+pi
        raux=2.d+0*lambda-raux
        vort_zita_n(n)=raux*cos(tetaaux)
        vort_eta_n(n)=raux*sin(tetaaux)
      endif

      return
      end

C ***********************************************************************************
C ***********************************************************************************
C ***
C ***     subroutine init_vort
C ***
C ***   inizializza i vortici emessi dalla piastra con Graham 1977
C ***      -1 = vortice destro , -2 = vortice sinistro conti pag.31,32
C ***
C ***********************************************************************************
C ***********************************************************************************
      
      subroutine init_vort(n2,n1)
      
      include 'airfoil_com.f'
      
      real*8 pos,ds
      real*8 v_eta
      real*8 dzita_em
      
      real*8 z1v,e1v,z2v,e2v
      real*8 zit1,et1,zit2,et2
      
      real*8 a11,a12,a21,a22
      real*8 d1,d2
      
      real*8 dpsi_u_dzita,dpsi_h_dzita,dpsi_p_dzita,dpsi_n_dzita
      
      if (n2.eq.1)then
        pos=-lambda
        ds = cg*(abs(v_eta(pos,0.d+0,n))*2.d+0*dt)**(2.d+0/3.d+0)
        dzita_em=sqrt(lambda*ds/2.d+0)
c       dzita_em=.01                 
        vort_zita_n(-2)=-lambda-dzita_em
        vort_eta_n(-2)=0.d+0
      endif
      if (n1.eq.1)then
        pos=lambda
        ds = cg*(abs(v_eta(pos,0.d+0,n))*2.d+0*dt)**(2.d+0/3.d+0)
        dzita_em=sqrt(lambda*ds/2.d+0)
c       dzita_em=.01                 
        vort_zita_n(-1)=lambda+dzita_em
        vort_eta_n(-1)=0.d+0
      endif

      zit1=lambda
      et1=0.d+0
      zit2=-lambda
      et2=0.d+0

      z1v=vort_zita_n(-1)
      e1v=vort_eta_n(-1)
      z2v=vort_zita_n(-2)
      e2v=vort_eta_n(-2)

      if(n2.eq.1.and.n1.eq.1)then
        a11=( lambda**2-z1v**2-e1v**2 )
     .   /( (z1v-lambda)**2+e1v**2 )
        a12=( lambda**2-z2v**2-e2v**2 )
     .   /( (z2v-lambda)**2+e2v**2 )
        a21=( -lambda**2+z1v**2+e1v**2 )
     .   /( (z1v+lambda)**2+e1v**2 )
        a22=( -lambda**2+z2v**2+e2v**2 )
     .   /( (z2v+lambda)**2+e2v**2 )
        d1=2.d+0*pi*lambda*(
     .   dpsi_u_dzita(zit1,et1)
     .      +dpsi_h_dzita(zit1,et1)
     .          +dpsi_p_dzita(zit1,et1)
     .                   )
        do i=1,nvort
          d1=d1+2.d+0*pi*lambda*dpsi_n_dzita(zit1,et1,i)
        enddo
        d2=2.d+0*pi*lambda*(
     .   dpsi_u_dzita(zit2,et2)
     .      +dpsi_h_dzita(zit2,et2)
     .          +dpsi_p_dzita(zit2,et2)
     .                   )
        do i=1,nvort
          d2=d2+2.d+0*pi*lambda*dpsi_n_dzita(zit2,et2,i)
        enddo
        vort_gamma_n(-1)=(d1*a22-d2*a12)/(a11*a22-a21*a12)
        vort_gamma_n(-2)=(d2*a11-d1*a21)/(a11*a22-a21*a12)
        if(abs(vort_gamma_n(-1)).lt..00005)vort_gamma_n(-1)=
     .      .00005*( vort_gamma_n(-1)/abs(vort_gamma_n(-1)) )
        if(abs(vort_gamma_n(-2)).lt..00005)vort_gamma_n(-2)=
     .      .00005*( vort_gamma_n(-2)/abs(vort_gamma_n(-2)) )
	write(6,*)tcorr,' gamma 1 e 2 = ',vort_gamma_n(-1),vort_gamma_n(-2)
c       read(5,*)io
      endif
      if(n2.eq.0.and.n1.eq.1)then
        a11=( lambda**2-z1v**2-e1v**2 )
     .   /( (z1v-lambda)**2+e1v**2 )
        d1=2.d+0*pi*lambda*(
     .   dpsi_u_dzita(zit1,et1)
     .      +dpsi_h_dzita(zit1,et1)
     .          +dpsi_p_dzita(zit1,et1)
     .                   )
	d1=d1+2.d+0*pi*lambda*dpsi_n_dzita(zit1,et1,-2)
        do i=1,nvort
          d1=d1+2.d+0*pi*lambda*dpsi_n_dzita(zit1,et1,i)
        enddo
        vort_gamma_n(-1) = d1/a11	 
	write(6,*)tcorr,' gamma 1 = ',vort_gamma_n(-1)
        if(abs(vort_gamma_n(-1)).lt..00005)vort_gamma_n(-1)=
     .      .00005*( vort_gamma_n(-1)/abs(vort_gamma_n(-1)) )
      endif
      if(n2.eq.1.and.n1.eq.0)then
        a22=( -lambda**2+z2v**2+e2v**2 )
     .   /( (z2v+lambda)**2+e2v**2 )
        d2=2.d+0*pi*lambda*(
     .   dpsi_u_dzita(zit2,et2)
     .      +dpsi_h_dzita(zit2,et2)
     .          +dpsi_p_dzita(zit2,et2)
     .                   )
        d2=d2+2.d+0*pi*lambda*dpsi_n_dzita(zit2,et2,-1)
        do i=1,nvort
          d2=d2+2.d+0*pi*lambda*dpsi_n_dzita(zit2,et2,i)
        enddo
        vort_gamma_n(-2) = d2/a22	 
	write(6,*)tcorr,' gamma 2 = ',vort_gamma_n(-2)
        if(abs(vort_gamma_n(-2)).lt..00005)vort_gamma_n(-2)=
     .      .00005*( vort_gamma_n(-2)/abs(vort_gamma_n(-2)) )
      endif
      
      return
      end

C ***********************************************************************************
C ***********************************************************************************
C ***
C ***        subroutine force 
C ***
C ***  calcola forze a momento conti da pag.23    quaderno arancio
C ***
C ***********************************************************************************
C ************************************************************************************

      subroutine force

      include 'airfoil_com.f'

      real*8 Xg,Yg

      real*8 Gam
      real*8 GamS
      real*8 GamD
      real*8 S1,S2
      real*8 ucap,vcap,ucapP,vcapP
      real*8 Mom

      complex*16 zit
      complex*16 Z1,Z2
      complex*16 For1,For2
    
C queste devono essere le componenti di velocita' dell'airfoil rispetto
C ad un sistema X,Y fermo istantaneamente coincidente con il sistema X,Y 
C solidale con l'airfoil
     
      ucap=-u0*cos(alfa-beta)+rlp*sin(alfa)
      vcap= u0*sin(alfa-beta)+rlp*cos(alfa)+alfap*dp
      ucapP=u0*sin(alfa-beta)*alfap
     .  +rlpp*sin(alfa)+rlp*cos(alfa)*alfap
      vcapP=u0*cos(alfa-beta)*alfap
     .  +rlpp*cos(alfa)-rlp*sin(alfa)*alfap+alfapp*dp

      GamS=0.d+0
      GamD=0.d+0
      Gam=0.d+0
      Z1=(0.d+0,0.d+0)
      Z2=(0.d+0,0.d+0)
      Z3=(0.d+0,0.d+0)
      Z4=(0.d+0,0.d+0)
      Z5=(0.d+0,0.d+0)
      S1=0.d+0
      S2=0.d+0

      t_cut(-2)=tcorr
      t_cut(-1)=tcorr
     
      do n=-2,nvort
        if(n.ne.0)then
 	  Zncap(n)=cmplx( 
     .	Xg(vort_zita(n),vort_eta(n))-dp,Yg(vort_zita(n),vort_eta(n))
     .                  )
          zit=cmplx(vort_zita(n),vort_eta(n))
          if(n.eq.-2)then	
            GamS=GamS+vort_gamma(n)
     .      *exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne)         
     .                                                  /(2.*pi)
          endif
          if(n.eq.-1)then	
            GamD=GamD+vort_gamma(n)
     .      *exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne)          
     .                                                  /(2.*pi)
          endif
c         if(n.gt.0)then
c           if(indcut(n).eq.-2)then
c             GamS=GamS+vort_gamma(n)
c    .      *exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne)         
c    .                                                  /(2.*pi)
c           endif
c           if(indcut(n).eq.-1)then
c             GamD=GamD+vort_gamma(n)
c    .      *exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne)          
c    .                                                  /(2.*pi)
c           endif
c         endif

          Gam=Gam+vort_gamma(n)
     .      *exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne)/(2.*pi)
          Z1=Z1+( vort_gamma(n)
     .      *exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne) )
     .       *lambda**2*( 1.d+0/zit + 1.d+0/conjg(zit) )
          Z2=Z2+( vort_gamma(n)
     .      *exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne) )
     .       *( dZncapdt(n)+cmplx(ucap,vcap)+ci*alfap*Zncap(n) )
          S1=S1+( vort_gamma(n)
     .      *exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne) )
     .       *dreal( lambda**4/zit**2 )
          S2=S2+( vort_gamma(n)
     .      *exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne) )
     .       *dreal( 
     .conjg( dZncapdt(n)+cmplx(ucap,vcap)+ci*alfap*Zncap(n) )*Zncap(n) 
     .             )

        endif
      enddo

c     do n=-2,nvort
c       if(n.ne.0)then
c	  Zncap(n)=cmplx( 
c    .	Xg(vort_zita(n),vort_eta(n))-dp,Yg(vort_zita(n),vort_eta(n))
c    .                  )
c         zit=cmplx(vort_zita(n),vort_eta(n))
c         if(n.eq.-2)then	
c           GamS=GamS+vort_gamma(n)
c    .                                                  /(2.*pi)
c         endif
c         if(n.eq.-1)then	
c           GamD=GamD+vort_gamma(n)
c    .                                                  /(2.*pi)
c         endif
c         if(n.gt.0)then
c           if(indcut(n).eq.-2)then
c             GamS=GamS+vort_gamma(n)
c    .                                                  /(2.*pi)
c           endif
c           if(indcut(n).eq.-1)then
c             GamD=GamD+vort_gamma(n)
c    .                                                  /(2.*pi)
c           endif
c         endif

c         Gam=Gam+vort_gamma(n)
c    .      /(2.*pi)
c         Z1=Z1+ vort_gamma(n)
c    .       *lambda**2*( 1.d+0/zit + 1.d+0/conjg(zit) )
c         Z2=Z2+ vort_gamma(n)
c    .       *( dZncapdt(n)+cmplx(ucap,vcap)+ci*alfap*Zncap(n) )
c         S1=S1+ vort_gamma(n)
c    .       *dreal( lambda**4/zit**2 )
c         S2=S2+ vort_gamma(n)
c    .       *dreal( 
c    .conjg( dZncapdt(n)+cmplx(ucap,vcap)+ci*alfap*Zncap(n) )*Zncap(n) 
c    .             )

c       endif
c     enddo

      if(indder.ne.0)then
      
      For1 =
     .      -4.d+0*pi*lambda**2*ci*vcapP
     .      -ci*(Z1-Z1_o)/(2.d+0*dt)
     .      +ci*2.d+0*lambda*(GamD-GamD_o)/(2.d+0*dt)
     .      -ci*2.d+0*lambda*(GamS-GamS_o)/(2.d+0*dt)
     .      +ci*Z2
c Nel termine seguente il valore di 
c Gam e' la somma della circolazione dei vortici emessi infatti
c l'espressione corretta sarebbe con il segno + e la circolazione
c dei vortici immagine
     .      -ci*cmplx(ucap,vcap)*Gam
     .      +4.d+0*pi*lambda**2*alfap*vcap
     .      +alfap*Z1
   
      Mom = -4.d+0*pi*lambda**4*alfapp
     .      -4.d+0*pi*lambda**2*ucap*vcap
     .      -(S1-S1_o)/(2.d+0*dt)
     .      +lambda**2*(Gam-Gam_o)/(2.d+0*dt)
     .      +S2
     .      -dreal( cmplx(ucap,-vcap)*Z1 )

C calcolo delle potenze

      FgX=dreal(For1)
      FgY= imag(For1)
      
      Fgxpic=FgX*cos(alfa)-FgY*sin(alfa)
      Fgypic=FgX*sin(alfa)+FgY*cos(alfa)

      if(rlh.ne.0.)then
        FgXadim=FgX/( (rlh*omega)**2*4.*lambda )
        FgYadim=FgY/( (rlh*omega)**2*4.*lambda )
        rMomadim=Mom/( (rlh*omega)**2*16.*lambda**2 )
      else
        FgXadim=FgX/( u0**2*4.*lambda )
        FgYadim=FgY/( u0**2*4.*lambda )
        rMomadim=Mom/( u0*16.*lambda**2 )
      endif
     
      Fgxpicadim=FgXadim*cos(alfa)-FgYadim*sin(alfa)
      Fgypicadim=FgXadim*sin(alfa)+FgYadim*cos(alfa)
     
c     write(14,10)tcorr,FgX,FgY,Mom
c     write(14,10)tcorr,FgXadim,FgYadim,rMomadim
10    format(7(1x,e11.5))

      pots=-Fgypic*rlp - Mom*alfap
      potu=Fgxpic*u0*cos(beta)
      
      if (tcorr.ge.t1.and.tcorr.le.t2) then
        if (iren.eq.0) then
          rendn1=potu
          rendd1=pots
          iren=1
        else
          rendn2=potu
          rendd2=pots
          rendn=rendn+0.5d+0*(2.d+0*dt)*(rendn1+rendn2)
          rendd=rendd+0.5d+0*(2.d+0*dt)*(rendd1+rendd2)
          rendn1=rendn2
          rendd1=rendd2
        endif
      endif
      
c     write(16,10)tcorr,rendn,rendd

      endif

      indder=1
      
      Z1_o=Z1
      S1_o=S1
      GamS_o=GamS
      GamD_o=GamD
      Gam_o=Gam

      return
      end

C ***********************************************************************************
C ***********************************************************************************
C ***
C ***        subroutine output
C ***
C ***********************************************************************************
C ***********************************************************************************

      subroutine output(time)

      include 'airfoil_com.f'
      character*70 filpri
      character*6 prti
     
      real*8 Xg,Yg,time
      
      itim=nint(time*1000.)
      write(prti,77)itim
 77   format(i6.6)      

      filpri = 'fld'//prti//'.dat'
      open(13,file=filpri,form='formatted')
      rewind(13)
     
      do i=-2,-1
       write(13,*)Xg(vort_zita(i),vort_eta(i)),
     .            Yg(vort_zita(i),vort_eta(i)),
     .            vort_zita(i),vort_eta(i),vort_gamma(i)
      enddo
       write(20,*)time,Xg(vort_zita(-1),vort_eta(-1)),
     .          Yg(vort_zita(-1),vort_eta(-1)),vort_gamma(-1)
       write(21,*)time,Xg(vort_zita(-2),vort_eta(-2)),
     .          Yg(vort_zita(-2),vort_eta(-2)),vort_gamma(-2)
      
      do i=1,nvort
       write(13,*)Xg(vort_zita(i),vort_eta(i)),
     .            Yg(vort_zita(i),vort_eta(i)),
     .            vort_zita(i),vort_eta(i),vort_gamma(i)
      enddo

      
      close(13)
      
      return
      end

c *****************************************************

      subroutine output1(time)

      include 'airfoil_com.f'
      
      real*8 Xg,Yg,time
      
      write(31,10)tcorr,Xg(vort_zita(-1),vort_eta(-1)),
     .            Yg(vort_zita(-1),vort_eta(-1)),
     .            vort_gamma(-1),dgam1dt
      write(32,10)tcorr,Xg(vort_zita(-2),vort_eta(-2)),
     .            Yg(vort_zita(-2),vort_eta(-2)),
     .            vort_gamma(-2),dgam2dt

      do i=1,nvort
        index = 40+i
	write(index,10)tcorr,Xg(vort_zita(i),vort_eta(i)),
     .            Yg(vort_zita(i),vort_eta(i)),
     .     vort_gamma(i)*exp(-dec*(omega*(tcorr-t_cut(i))/pi)**enne)
	
      enddo
 10   format(5(1x,e12.6))

      return
      end



