cS non c'e' piu' vort_o_gamma
      parameter (NP=205,nv=400,num=5000)

      IMPLICIT REAL*8 (A-H,O-Z)

      integer enne

      real*8 pi,cg,rho
      real*8 dt,tfin,tcorr,gamma_max,perc,t1,t2
cS u0P
      real*8 lambda,dp,u0,u0P,beta,rlh,omega,alfah,phi_a
      real*8 rl,rlp,rlpp,alfa,alfap,alfapp
      real*8 rendn1,rendn2,rendn,rendd1,rendd2,rendd  
      real*8 Gam_o,GamS_o,GamD_o,S1_o
      real*8 r0,cost0,cost1

      complex*16 ci
      complex*16 Zncap,dZncapdt
      complex*16 Z1_o
      complex*16 VET_Z

      real*8 vort_zita,vort_eta,vort_gamma,
     .       t_cut,dec, 
     .       vort_zita_n,vort_eta_n,vort_gamma_n,
     .       vort_zita_v,vort_eta_v,vort_gamma_v 
            
      real*8 dzit1dt,dzit2dt,det1dt,det2dt,dgam1dt,dgam2dt
	
      dimension vort_zita(-2:nv),vort_eta(-2:nv),vort_gamma(-2:nv),
     .       t_cut(-2:nv),
     .       vort_zita_n(-2:nv),vort_eta_n(-2:nv),vort_gamma_n(-2:nv),
     .       vort_zita_v(-2:nv),vort_eta_v(-2:nv),vort_gamma_v(-2:nv),
     .      indcut(nv)
	
      common/gen/pi,ci,cg,rho
     
      common/run/dzit1dt,dzit2dt,det1dt,det2dt,dgam1dt,dgam2dt
      
      common/vort/Zncap(-2:nv),dZncapdt(-2:nv),
     .            vort_zita,vort_eta,vort_gamma,
     .            t_cut,dec,
     .            vort_zita_n,vort_eta_n,vort_gamma_n,
     .            vort_zita_v,vort_eta_v,vort_gamma_v,
     .            Gam_o,GamS_o,GamD_o,S1_o,Z1_o,indder,enne,nvort
      
      common/intgr/dt,tfin,tcorr,gamma_max,perc,t1,t2,ndtusc
      
      common/geom/lambda,dp,
     .  	u0,u0P,beta,rlh,omega,alfah,phi_a,
     .    	rl,rlp,rlpp,alfa,alfap,alfapp
      
      common/pot/rendn1,rendn2,rendn,rendd1,rendd2,rendd,iren

      COMMON/PUNTIV/VET_Z(NP),NUM_PUNTI,IPUNTIVEL
      
      COMMON/CAMPOMED/U_HAT(NP,NP),V_HAT(NP,NP),YHAT(NP),XHAT(NP),
     $     DTGETTO,YHAT_MAX,YHAT_MIN,XHAT_MAX,XHAT_MIN,
     $     PIN,PFIN,NYHAT,NXHAT,ICMEDIO
   
      COMMON/CORE/R0,COST0,COST1
