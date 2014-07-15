
C***************************************************************
C***************************************************************
C***
C***    subroutines per il calcolo del campo di moto
C***
C***************************************************************
C***************************************************************
     
      real*8 function dpsi_u_dzita(zita,eta)

C  calcolo contributo d_psi_dzita associato al moto base
C  pagina d degli appunti piano zita-eta

      include 'airfoil_com.f'

      real*8 zita,eta,aux
	
      aux = zita**2 + eta**2

      dpsi_u_dzita = u0*(
     .    2.d+0*zita*lambda**2*
     .         ( eta*cos(alfa-beta)+zita*sin(alfa-beta) )/aux**2 +
     .         (1.d+0-lambda**2.d+0/aux)*sin(alfa-beta) 
     .                 )

      return
      end

C***************************************************************
     
      real*8 function dpsi_u_dzitaP(zita,eta)

C  calcolo contributo d_psi_dzita associato al moto base
C  pagina d degli appunti piano zita-eta

      include 'airfoil_com.f'

      real*8 zita,eta,aux
	
      aux = zita**2 + eta**2

      dpsi_u_dzitaP = u0P*(
     .    2.d+0*zita*lambda**2*
     .         ( eta*cos(alfa-beta)+zita*sin(alfa-beta) )/aux**2 +
     .         (1.d+0-lambda**2.d+0/aux)*sin(alfa-beta) 
     .                    )
     .              + u0*(
     .    2.d+0*zita*lambda**2*
     .   alfap*( -eta*sin(alfa-beta)+zita*cos(alfa-beta) )/aux**2 +
     .   alfap*(1.d+0-lambda**2.d+0/aux)*cos(alfa-beta) 
     .                   )

      return
      end

C***************************************************************

      real*8 function dpsi_u_deta(zita,eta)

C  calcolo contributo d_psi_deta associato al moto base
C  pagina d degli appunti piano zita-eta

      include 'airfoil_com.f'
       
      real*8 zita,eta,aux

      aux = zita**2 + eta**2

      dpsi_u_deta = u0*(
     .      2.d+0*eta*lambda**2*
     .         ( eta*cos(alfa-beta)+zita*sin(alfa-beta) )/aux**2 +
     .         (1.-lambda**2/aux)*cos(alfa-beta)
     .                   )

      return
      end

C***************************************************************

      real*8 function dpsi_h_dzita(zita,eta)

C    calcolo contributo d_psi_dzita associato al moto di oscillazione in 
C    direzione verticale
C    pagina d degli appunti piano zita-eta

      include 'airfoil_com.f'
      
      real*8 zita,eta,aux

      aux = zita**2 + eta**2

      dpsi_h_dzita = rlp*(  2.*zita*lambda**2*
     .              ( zita*cos(alfa)-eta*sin(alfa) )/aux**2 +
     .               (1.-lambda**2/aux)*cos(alfa)  )

      return
      end

C***************************************************************

      real*8 function dpsi_h_dzitaP(zita,eta)

C    calcolo contributo d_psi_dzita associato al moto di oscillazione in 
C    direzione verticale
C    pagina d degli appunti piano zita-eta

      include 'airfoil_com.f'
      
      real*8 zita,eta,aux

      aux = zita**2 + eta**2

      dpsi_h_dzitaP = rlpp*(  2.*zita*lambda**2*
     .              ( zita*cos(alfa)-eta*sin(alfa) )/aux**2 +
     .               (1.-lambda**2/aux)*cos(alfa)  )
     .               + rlp*(  2.*zita*lambda**2*
     .        alfap*( -zita*sin(alfa)-eta*cos(alfa) )/aux**2 +
     .         alfap*(1.-lambda**2/aux)*sin(alfa)  )

      return
      end

C***************************************************************

      real*8 function dpsi_h_deta(zita,eta)

C    calcolo contributo d_psi_deta associato al moto di oscillazione in 
C    direzione verticale
C    pagina d degli appunti piano zita-eta

      include 'airfoil_com.f'
      
      real*8 zita,eta,aux

      aux = zita**2 + eta**2

      dpsi_h_deta = rlp*(  2.*eta*lambda**2*
     .              ( zita*cos(alfa)-eta*sin(alfa) )/aux**2 -
     .               (1.-lambda**2/aux)*sin(alfa)  )

      return
      end

C***************************************************************

      real*8 function dpsi_p_dzita(zita,eta)

C    calcolo contributo d_psi_deta associato al moto di rotazione
C    pagina d degli appunti piano zita-eta

      include 'airfoil_com.f'
      
      real*8 zita,eta,aux

      aux = zita**2 + eta**2

      dpsi_p_dzita = alfap*(  zita - zita*lambda**4/aux**2 +
     .                lambda**2*(zita**2-eta**2)*
     .                (4.*zita*lambda**2/aux**3 - 2.*zita/aux**2 ) +
     .                 2.*zita*lambda**2*(1/aux-lambda**2/aux**2) +
     .                 2.*dp*zita**2*lambda**2/aux**2 +
     .                 dp*(1.-lambda**2/aux)    )

      return
      end

C***************************************************************

      real*8 function dpsi_p_dzitaP(zita,eta)

C    calcolo contributo d_psi_deta associato al moto di rotazione
C    pagina d degli appunti piano zita-eta

      include 'airfoil_com.f'
      
      real*8 zita,eta,aux

      aux = zita**2 + eta**2

      dpsi_p_dzitaP = alfapp*(  zita - zita*lambda**4/aux**2 +
     .                lambda**2*(zita**2-eta**2)*
     .                (4.*zita*lambda**2/aux**3 - 2.*zita/aux**2 ) +
     .                 2.*zita*lambda**2*(1/aux-lambda**2/aux**2) +
     .                 2.*dp*zita**2*lambda**2/aux**2 +
     .                 dp*(1.-lambda**2/aux)    )

      return
      end

C***************************************************************

      real*8 function dpsi_p_deta(zita,eta)

C    calcolo contributo d_psi_deta associato al moto di rotazione
C    pagina d degli appunti piano zita-eta

      include 'airfoil_com.f'
      
      real*8 zita,eta,aux

      aux = zita**2 + eta**2

      dpsi_p_deta = alfap*(  eta - eta*lambda**4/aux**2 +
     .                lambda**2*(zita**2-eta**2)*
     .                (4.*eta*lambda**2/aux**3 - 2.*eta/aux**2 ) -
     .                 2.*eta*lambda**2*(1/aux-lambda**2/aux**2) +
     .                 2.*dp*zita*eta*lambda**2/aux**2  )

      return
      end

C***************************************************************

      real*8 function dpsi_n_dzita(zita,eta,n)

C    calcolo contributo d_psi_dzita associato al moto indotto dal 
C    vortice n-esimo
C    pagina e degli appunti piano zita-eta

      include 'airfoil_com.f'
      
      real*8 zita,eta

      if(n.eq.-1.or.n.eq.-2)t_cut(n)=tcorr

      dpsi_n_dzita=-vort_gamma(n)
     .      *exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne)/(2.*pi)*
     .     (
     .           (zita-vort_zita(n))/
     .           ((zita-vort_zita(n))**2 + (eta-vort_eta(n))**2) -
     .       (  
     .        
     .     (zita*vort_zita(n)+eta*vort_eta(n)-lambda**2)*vort_zita(n) -
     .        (eta*vort_zita(n)-zita*vort_eta(n) )*vort_eta(n) 
     .        ) /
     .       ( (zita*vort_zita(n)+eta*vort_eta(n)-lambda**2)**2 +
     .        (eta*vort_zita(n)-zita*vort_eta(n))**2  )
     .     )

      return
      end

C***************************************************************

      real*8 function dpsi_n_dzitaP(zita,eta,n)

C    calcolo contributo d_psi_dzita associato al moto indotto dal 
C    vortice n-esimo
C    pagina e degli appunti piano zita-eta

      include 'airfoil_com.f'
      
      real*8 zita,eta
      real*8 zitac,etac
      real*8 znP,enP

      real*8 dxgdzita,dygdzita,dxgdeta,dygdeta,rJg
      real*8 dpsi_t_dzita,dpsi_t_deta,uvor,vvor

      if(n.eq.-1.or.n.eq.-2)t_cut(n)=tcorr

      dpsi_n_dzitaP=vort_gamma(n)
     .      *dec*enne*(omega/pi)*(omega*(tcorr-t_cut(n))/pi)**(enne-1.)
     .      *(exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne)/(2.*pi))*
     .     (
     .           (zita-vort_zita(n))/
     .           ((zita-vort_zita(n))**2 + (eta-vort_eta(n))**2) -
     .       (  
     .        
     .     (zita*vort_zita(n)+eta*vort_eta(n)-lambda**2)*vort_zita(n) -
     .        (eta*vort_zita(n)-zita*vort_eta(n) )*vort_eta(n) 
     .        ) /
     .       ( (zita*vort_zita(n)+eta*vort_eta(n)-lambda**2)**2 +
     .        (eta*vort_zita(n)-zita*vort_eta(n))**2  )
     .     )

      zitac = vort_zita(n)
      etac = vort_eta(n)

      znP = 
     .      dpsi_t_deta(zitac,etac,n)  
     .   +  (    uvor(n)*dygdeta(zitac,etac) -
     .           vvor(n) *dxgdeta(zitac,etac)    )

      enP = 
     .      - dpsi_t_dzita(zitac,etac,n) + 
     .      ( - uvor(n)*dygdzita(zitac,etac) +
     .          vvor(n)*dxgdzita(zitac,etac)  )  

      znP = znP/rJg(zitac,etac)
      enP = enP/rJg(zitac,etac)

      dpsi_n_dzitaP=dpsi_n_dzitaP
     .              -vort_gamma(n)
     .      *(exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne)/(2.*pi))*
     .     (
     .       -znP/
     .        ( (zita-vort_zita(n))**2+(eta-vort_eta(n))**2 )
     .       +2.d+0*(zita-vort_zita(n))*
     .             ( (zita-vort_zita(n))*znP+(eta-vort_eta(n))*enP )/
     .        ( (zita-vort_zita(n))**2+(eta-vort_eta(n))**2 )**2
     .-(  
     .     (zita*znP+eta*enP)*vort_zita(n)
     .    +(zita*vort_zita(n)+eta*vort_eta(n)-lambda**2)*znP
     .    -(eta*znP-zita*enP)*vort_eta(n) 
     .    -(eta*vort_zita(n)-zita*vort_eta(n))*enP 
     . ) /( (zita*vort_zita(n)+eta*vort_eta(n)-lambda**2)**2 +
     .             (eta*vort_zita(n)-zita*vort_eta(n))**2  )
     .+(    
     .   ( (zita*vort_zita(n)+eta*vort_eta(n)-lambda**2)*vort_zita(n)
     .    -(eta*vort_zita(n)-zita*vort_eta(n))*vort_eta(n) ) 
     .  *2.d+0*(    
     .           (zita*vort_zita(n)+eta*vort_eta(n)-lambda**2)
     .            *(zita*znP+eta*enP)
     .          +(eta*vort_zita(n)-zita*vort_eta(n))
     .            *(eta*znP-zita*enP) 
     .         )
     . ) /( (zita*vort_zita(n)+eta*vort_eta(n)-lambda**2)**2 +
     .             (eta*vort_zita(n)-zita*vort_eta(n))**2  )**2

     .     )
      
      return
      end

C***************************************************************

      real*8 function dpsi_n_deta(zita,eta,n)

C    calcolo contributo d_psi_deta associato al moto indotto dal 
C    vortice n-esimo
C    pagina e degli appunti piano zita-eta

      include 'airfoil_com.f'
      
      real*8 zita,eta
      
      if(n.eq.-1.or.n.eq.-2)t_cut(n)=tcorr

      dpsi_n_deta=-vort_gamma(n)
     .         *exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne)/(2.*pi)*
     .     (
     .           (eta-vort_eta(n))/
     .           ((zita-vort_zita(n))**2 + (eta-vort_eta(n))**2) -
     .       (  
     .        
     .     (zita*vort_zita(n)+eta*vort_eta(n)-lambda**2)*vort_eta(n) +
     .        (eta*vort_zita(n)-zita*vort_eta(n) )*vort_zita(n) 
     .        ) /
     .       ( (zita*vort_zita(n)+eta*vort_eta(n)-lambda**2)**2 +
     .        (eta*vort_zita(n)-zita*vort_eta(n))**2  )
     .     )

      return
      end

C***************************************************************

      subroutine coeff

C    calcolo dei diversi coefficienti che compaino nelle equazioni
C    originate dalla condizione di Kutta e dalle equazioni di Brown
C    and Michael

      include 'airfoil_com.f'
     
      real*8 z1v,e1v,z2v,e2v
      real*8 zit1,et1,zit2,et2
 
      real*8 a11,a12,a21,a22
      real*8 b11,b12,b21,b22
      real*8 c11,c12,c21,c22
      real*8 h11,h12,h21,h22
      real*8 l1,l2
      real*8 Ag11,Bg11,Ag12,Bg12,Ag21,Bg21,Ag22,Bg22
      real*8 d1P,d2P

      real*8 Xg,Yg,dxgdzita,dxgdeta,dygdzita,dygdeta,rJg
      real*8 dpsi_u_dzitaP,dpsi_h_dzitaP,dpsi_p_dzitaP,dpsi_n_dzitaP
      real*8 dpsi_t_dzita,dpsi_t_deta,uvor,vvor

      z1v=vort_zita(-1)
      e1v=vort_eta(-1)
      z2v=vort_zita(-2)
      e2v=vort_eta(-2)

      zit1=lambda
      zit2=-lambda
      et1=0.d+0
      et2=0.d+0
      
      a11=( lambda**2-vort_zita(-1)**2-vort_eta(-1)**2 )
     .   /( (vort_zita(-1)-lambda)**2+vort_eta(-1)**2 )
      a12=( lambda**2-vort_zita(-2)**2-vort_eta(-2)**2 )
     .   /( (vort_zita(-2)-lambda)**2+vort_eta(-2)**2 )
      a21=( -lambda**2+vort_zita(-1)**2+vort_eta(-1)**2 )
     .   /( (vort_zita(-1)+lambda)**2+vort_eta(-1)**2 )
      a22=( -lambda**2+vort_zita(-2)**2+vort_eta(-2)**2 )
     .   /( (vort_zita(-2)+lambda)**2+vort_eta(-2)**2 )

      Ag11=2.d+0*lambda
     .    *( (vort_zita(-1)-lambda)**2-vort_eta(-1)**2 )
     .    /( (vort_zita(-1)-lambda)**2+vort_eta(-1)**2 )**2

      Bg11=-4.d+0*lambda
     .    *vort_eta(-1)*( lambda-vort_zita(-1) )
     .    /( (vort_zita(-1)-lambda)**2+vort_eta(-1)**2 )**2

      Ag12=2.d+0*lambda
     .    *( (vort_zita(-2)-lambda)**2-vort_eta(-2)**2 )
     .    /( (vort_zita(-2)-lambda)**2+vort_eta(-2)**2 )**2

      Bg12=-4.d+0*lambda
     .    *vort_eta(-2)*( lambda-vort_zita(-2) )
     .    /( (vort_zita(-2)-lambda)**2+vort_eta(-2)**2 )**2

      Ag21=2.d+0*lambda
     .    *( (vort_zita(-1)+lambda)**2-vort_eta(-1)**2 )
     .    /( (vort_zita(-1)+lambda)**2+vort_eta(-1)**2 )**2

      Bg21=4.d+0*lambda
     .    *vort_eta(-1)*( lambda+vort_zita(-1) )
     .    /( (vort_zita(-1)+lambda)**2+vort_eta(-1)**2 )**2

      Ag22=2.d+0*lambda
     .    *( (vort_zita(-2)+lambda)**2-vort_eta(-2)**2 )
     .    /( (vort_zita(-2)+lambda)**2+vort_eta(-2)**2 )**2

      Bg22=4.d+0*lambda
     .    *vort_eta(-2)*( lambda+vort_zita(-2) )
     .    /( (vort_zita(-2)+lambda)**2+vort_eta(-2)**2 )**2

      b11=(   (Xg(z1v,e1v)-(2.d+0*lambda+dp))*dygdeta(z1v,e1v)
     .       -Yg(z1v,e1v)*dxgdeta(z1v,e1v)  )
     .  /(vort_gamma(-1)*rJg(z1v,e1v)) 
      b12=(  -(Xg(z1v,e1v)-(2.d+0*lambda+dp))*dygdzita(z1v,e1v)
     .       +Yg(z1v,e1v)*dxgdzita(z1v,e1v)  )
     .  /(vort_gamma(-1)*rJg(z1v,e1v)) 
      b21=(   (Xg(z2v,e2v)+2.d+0*lambda-dp)*dygdeta(z2v,e2v)
     .       -Yg(z2v,e2v)*dxgdeta(z2v,e2v)  )
     .  /(vort_gamma(-2)*rJg(z2v,e2v)) 
      b22=(  -(Xg(z2v,e2v)+2.d+0*lambda-dp)*dygdzita(z2v,e2v)
     .       +Yg(z2v,e2v)*dxgdzita(z2v,e2v)  )
     .  /(vort_gamma(-2)*rJg(z2v,e2v)) 
      
      c11=(   dpsi_t_deta(z1v,e1v,-1)+uvor(-1)*dygdeta(z1v,e1v)
     .            -vvor(-1)*dxgdeta(z1v,e1v)  )/rJg(z1v,e1v)
      c12=(  -dpsi_t_dzita(z1v,e1v,-1)-uvor(-1)*dygdzita(z1v,e1v)
     .            +vvor(-1)*dxgdzita(z1v,e1v)  )/rJg(z1v,e1v)

      c21=(   dpsi_t_deta(z2v,e2v,-2)+uvor(-2)*dygdeta(z2v,e2v)
     .            -vvor(-2)*dxgdeta(z2v,e2v)  )/rJg(z2v,e2v)
      c22=(  -dpsi_t_dzita(z2v,e2v,-2)-uvor(-2)*dygdzita(z2v,e2v)
     .            +vvor(-2)*dxgdzita(z2v,e2v)  )/rJg(z2v,e2v)

      h11=a11-vort_gamma(-1)*(Ag11*b11+Bg11*b12)
      h12=a12-vort_gamma(-2)*(Ag12*b21+Bg12*b22)
      h21=a21-vort_gamma(-1)*(Ag21*b11+Bg21*b12)
      h22=a22-vort_gamma(-2)*(Ag22*b21+Bg22*b22)
      
      d1P=2.d+0*pi*lambda*(
     .   dpsi_u_dzitaP(zit1,et1)
     .      +dpsi_h_dzitaP(zit1,et1)
     .          +dpsi_p_dzitaP(zit1,et1)
     .                   )
     
      do i=1,nvort
        d1P=d1P+2.d+0*pi*lambda*dpsi_n_dzitaP(zit1,et1,i)
      enddo
      
      d2P=2.d+0*pi*lambda*(
     .   dpsi_u_dzitaP(zit2,et2)
     .      +dpsi_h_dzitaP(zit2,et2)
     .          +dpsi_p_dzitaP(zit2,et2)
     .                   )
      
      do i=1,nvort
        d2P=d2P+2.d+0*pi*lambda*dpsi_n_dzitaP(zit2,et2,i)
      enddo

c     write(6,*)' **** ' 
c     write(6,*)' d1P = ',d1P
c     write(6,*)' d2P = ',d2P
c     write(6,*)' Gamma1 = ',vort_gamma(-1)
c     write(6,*)' Gamma2 = ',vort_gamma(-2)
c     write(6,*)' zita1 = ',zit1
c     write(6,*)' zita2 = ',zit2
c     write(6,*)' eta1 = ',et1
c     write(6,*)' eta2 = ',et2
c     write(6,*)' zita1 V = ',z1v
c     write(6,*)' zita2 V = ',z2v
c     write(6,*)' eta1 V = ',e1v
c     write(6,*)' eta2 V = ',e2v
c    
      l1=d1P-vort_gamma(-1)*(Ag11*c11+Bg11*c12)
     .      -vort_gamma(-2)*(Ag12*c21+Bg12*c22)
      l2=d2P-vort_gamma(-1)*(Ag21*c11+Bg21*c12)
     .      -vort_gamma(-2)*(Ag22*c21+Bg22*c22)

      dgam1dt=(l1*h22-l2*h12)/(h11*h22-h21*h12)
      dgam2dt=(l2*h11-l1*h21)/(h11*h22-h21*h12)

c     dzit1dt=-b11*( (l1*h22-l2*h12)/(h11*h22-h21*h12) )+c11
c     det1dt=-b12*( (l1*h22-l2*h12)/(h11*h22-h21*h12) )+c12
c     dzit2dt=-b21*( (l2*h11-l1*h21)/(h11*h22-h21*h12) )+c21
c     det2dt=-b22*( (l2*h11-l1*h21)/(h11*h22-h21*h12) )+c22
      
      dzit1dt1=-b11*dgam1dt
      dzit1dt2=c11
      dzit2dt1=-b21*dgam2dt
      dzit2dt2=c21
      
      dzit1dt=-b11*dgam1dt+c11
      det1dt =-b12*dgam1dt+c12
      dzit2dt=-b21*dgam2dt+c21
      det2dt =-b22*dgam2dt+c22
     
c     write(6,*)' **** ' 
c     write(6,*)' b11 = ',b11
c     write(6,*)' b21 = ',b21
c     write(6,*)' c11 = ',c11
c     write(6,*)' c12 = ',c12
c     write(6,*)' c21 = ',c21
c     write(6,*)' c22 = ',c22
c     write(6,*)' h11 = ',h11
c     write(6,*)' h12 = ',h12
c     write(6,*)' h21 = ',h21
c     write(6,*)' h22 = ',h22
c     write(6,*)' l1 = ',l1
c     write(6,*)' l2 = ',l2
c     write(6,*)' dzit1dt1 = ',dzit1dt1
c     write(6,*)' dzit1dt2 = ',dzit1dt2
c     write(6,*)' dzit2dt1 = ',dzit2dt1
c     write(6,*)' dzit2dt2 = ',dzit2dt2
c     write(6,*)' det1dt = ',det1dt
c     write(6,*)' det2dt = ',det2dt
c     write(6,*)' dgam1dt = ',dgam1dt
c     write(6,*)' dgam2dt = ',dgam2dt
c     write(6,*)' **** ' 
      
      return
      end

C***************************************************************

      real*8 function rJg(zita,eta)

C  calcolo lo jacobiano nel piano zita-eta
C     pagina 8 quaderno viola con e=0.

      include 'airfoil_com.f'

      real*8 zita,eta

      rJg = 1. + lambda**2*
     .     (lambda**2 - 2.*(zita**2-eta**2))/
     .     (zita**2+eta**2)**2

      return
      end

C***************************************************************

      real*8 function dxgdzita(zita,eta)

C  calcolo di d-xgranda su dzita nel piano zita-eta
C     pagina 7 quaderno viola con e=0.

      include 'airfoil_com.f'
      
      real*8 zita,eta

      dxgdzita = 1. - lambda**2*
     .     (zita**2-eta**2)/
     .     (zita**2+eta**2)**2

      return
      end


C***************************************************************

      real*8 function dxgdeta(zita,eta)

C  calcolo di d-xgranda su deta nel piano zita-eta
C     pagina 7 quaderno viola con e=0.

      include 'airfoil_com.f'
      
      real*8 zita,eta

      dxgdeta = - 2.*lambda**2*
     .      zita*eta/(zita**2+eta**2)**2

      return
      end

C***************************************************************

      real*8 function dygdzita(zita,eta)

C  calcolo di d-ygrande su dzita nel piano zita-eta
C  pagina 7 quaderno viola con e=0.

      include 'airfoil_com.f'
      
      real*8 zita,eta

      dygdzita =  2.*lambda**2*
     .      zita*eta/(zita**2+eta**2)**2

      return
      end

C***************************************************************

      real*8 function dygdeta(zita,eta)

C  calcolo di d-ygrande su deta nel piano zita-eta
C  pagina 7 quaderno viola con e=0.

      include 'airfoil_com.f'
      
      real*8 zita,eta

      dygdeta = 1. - lambda**2*
     .      (zita**2-eta**2)/(zita**2+eta**2)**2

      return
      end

C***************************************************************

      real*8 function uvor(n)

C calcola la componente u velocità del vortice n-esimo
C   piano zita-eta
C pagina c -arancio - attenzione segno cambiato

     
      include 'airfoil_com.f'

      complex*16 aux,aux1,aux2

      aux1 = vort_zita(n)-ci*vort_eta(n)
      aux2 = vort_zita(n)+ci*vort_eta(n)
      
      if(n.eq.-1.or.n.eq.-2)t_cut(n)=tcorr

      aux = ci*vort_gamma(n)
     .     *exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne)
     .	*aux1 / (2.*pi* (aux1**2-lambda**2))*
     .      ( - aux1*aux2/(aux1*aux2-lambda**2)-
     .          lambda**2/(aux1**2-lambda**2)   )
      uvor = dreal(aux)

      return
      end

C***************************************************************

      real*8 function vvor(n)

C calcola la componente u velocità del vortice n-esimo
C piano zita-eta
C pagina c -arancio - attenzione segno cambiato

      include 'airfoil_com.f'
	
      complex*16  aux,aux1,aux2

      aux1 = vort_zita(n)-ci*vort_eta(n)
      aux2 = vort_zita(n)+ci*vort_eta(n)
      
      if(n.eq.-1.or.n.eq.-2)t_cut(n)=tcorr

      aux = ci*vort_gamma(n)*
     .     exp(-dec*(omega*(tcorr-t_cut(n))/pi)**enne)
     .  *aux1 / (2.*pi* (aux1**2-lambda**2))*
     .      ( - aux1*aux2/(aux1*aux2-lambda**2)-
     .          lambda**2/(aux1**2-lambda**2)   )
      vvor = imag(aux)

      return
      end

C***************************************************************

      real*8 function v_eta_aux(zita,eta)

C calcola la componente eta di v per la init - calcoli retro p.16

      include 'airfoil_com.f'
     
      real*8 zita,eta
      real*8 dpsidz
      real*8 dpsi_u_dzita,dpsi_h_dzita,dpsi_p_dzita,dpsi_n_dzita
      
      dpsidz = dpsi_u_dzita(zita,eta) + dpsi_h_dzita(zita,eta) +
     .         dpsi_p_dzita(zita,eta)
      
      ii=-2
      dpsidz=dpsidz+dpsi_n_dzita(zita,eta,ii)

      ii=-1
      dpsidz=dpsidz+dpsi_n_dzita(zita,eta,ii)

      if (nvort.ge.1) then
        do nn=1,nvort
          dpsidz = dpsidz + dpsi_n_dzita(zita,eta,nn)
        enddo
      endif
      
      v_eta_aux = -dpsidz

      return
      end

C***************************************************************

      real*8 function v_eta(zita,eta,n)

C calcola la componente eta di v per la init - calcoli retro p.16

      include 'airfoil_com.f'
     
      real*8 zita,eta
      real*8 dpsidz
      real*8 dpsi_u_dzita,dpsi_h_dzita,dpsi_p_dzita,dpsi_n_dzita
      
      dpsidz = dpsi_u_dzita(zita,eta) + dpsi_h_dzita(zita,eta) +
     .         dpsi_p_dzita(zita,eta)
      
      if(n.ne.-1)then 
        ii=-1 
        dpsidz=dpsidz+dpsi_n_dzita(zita,eta,ii)
      endif

      if(n.ne.-2)then
        ii=-2
        dpsidz=dpsidz+dpsi_n_dzita(zita,eta,ii)
      endif

      if (nvort.ge.1) then
        do nn=1,nvort
          dpsidz = dpsidz + dpsi_n_dzita(zita,eta,nn)
        enddo
      endif
      
      v_eta = -dpsidz

      return
      end

C***************************************************************

      real*8 function Xg(zita,eta)
	
C calcola la trasformazione pag a retro

      include 'airfoil_com.f'
      
      real*8 zita,eta

      Xg = dp + zita*(1.+lambda**2/(zita**2+eta**2))
      
      return
      end

C***************************************************************

      real*8 function Yg(zita,eta)
	
C calcola la trasformazione pag a retro

      include 'airfoil_com.f'
        
      real*8 zita,eta

      Yg = eta*(1.-lambda**2/(zita**2+eta**2))
	
      return
      end

C***************************************************************

      real*8 function dpsi_t_dzita(zita,eta,n)
       
C n = -1 vortice dx, -2 sin
       
      include 'airfoil_com.f'
	
      real*8 zita,eta,aux
      real*8 dpsi_u_dzita,dpsi_h_dzita,dpsi_p_dzita,dpsi_n_dzita
       
      aux = 0.
       
      if (n.ne.-1) aux = aux + dpsi_n_dzita(zita,eta,-1)
      if (n.ne.-2) aux = aux + dpsi_n_dzita(zita,eta,-2)
       
      do i=1,nvort
        if (i.ne.n) aux = aux + dpsi_n_dzita(zita,eta,i)
      enddo
       
      dpsi_t_dzita = dpsi_u_dzita(zita,eta)
     .              +dpsi_h_dzita(zita,eta)
     .              +dpsi_p_dzita(zita,eta)
     .              +aux

      return
      end

C***************************************************************       

      real*8 function dpsi_t_deta(zita,eta,n)
       
C n = -1 vortice dx, -2 sin
       
      include 'airfoil_com.f'
	
      real*8 zita,eta,aux
      real*8 dpsi_u_deta,dpsi_h_deta,dpsi_p_deta,dpsi_n_deta
       
      aux = 0.d+0
       
      if (n.ne.-1) aux = aux + dpsi_n_deta(zita,eta,-1)
      if (n.ne.-2) aux = aux + dpsi_n_deta(zita,eta,-2)
       
      do i=1,nvort
        if (i.ne.n) aux = aux + dpsi_n_deta(zita,eta,i)
      enddo
       
      dpsi_t_deta = dpsi_u_deta(zita,eta) +
     .               dpsi_h_deta(zita,eta) +
     .               dpsi_p_deta(zita,eta) +
     .               aux

      return
      end

C***************************************************************       

      subroutine rlcorr
       
      include 'airfoil_com.f'
      
      rl = rlh*sin(omega*tcorr)
       
      return
      end

C***************************************************************       

      subroutine alfacorr
       
      include 'airfoil_com.f'
       
      alfa = alfah*sin(omega*tcorr+phi_a)
       
      return
      end

C***************************************************************       

      subroutine rlpcorr
       
      include 'airfoil_com.f'
       
      rlp = rlh*omega*cos(omega*tcorr)
      
       
      return
      end

C***************************************************************       

      subroutine alfapcorr
       
      include 'airfoil_com.f'
       
      alfap = alfah*omega*cos(omega*tcorr+phi_a)
       
      return
      end

C***************************************************************       

      subroutine rlppcorr
       
      include 'airfoil_com.f'
       
      rlpp = -rlh*omega**2*sin(omega*tcorr)
       
      return
      end

C***************************************************************       

      subroutine alfappcorr
       
      include 'airfoil_com.f'
       
      alfapp = -alfah*omega**2*sin(omega*tcorr+phi_a)
       
      return
      end

C***************************************************************       

