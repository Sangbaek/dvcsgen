    	program main
      implicit real*8(a-h,k-m,o-z)
#include "dvcs.inc"
	    real*8 ebeam,s,q2,x,t,phi,egamma,thetag,phig,stot,delta
        integer*4 ipol,ikeygene,ichnnel,iddd
	integer*4 iappr,hel
        parameter(nxq2poi=1)
        parameter(ntpoi=4)
        parameter(nphipoi=10)		
        parameter(ndelta=6)		
        real*8 arx(nxq2poi)
        real*8 arq2(nxq2poi)
        real*8 art(ntpoi)
		real*8 ardelta(ndelta)
        data arx/0.4d0/     
        data arq2/1.8d0/      
        data art/-4.045d0,-2.774d0,-1.504d0,-0.234d0/     
        data ardelta/1d-6,1d-5,1d-4,1d-3,1d-2,1d-1/

		 
         ebeam=10.604

         ipol=0 !target polarization switch:0-unpolarized, 1-long, 2-trans along x-axis, 3-trans along y-axis
         iappr=2   !1-BH only, 2-BMK 
         Ivar = 3
         Ich = -1
         IGPD = 3
         Ipn = 1
         Ed = ebeam

       if(iappr.eq.2)then
         hel=1     !             hel=1(0) to (not) include polarized part   
       endif

c             vv2cut=.3d0    ! missing mass cut        
            vv2cut=999999999.3d0    ! no missing mass cut used

           open(8,file='xykin.dat')
           do i=1,100000  


           read(8,*,end=111)x,q2,tpl,phi		   

				        t=-tpl
 
				  
          
c		  print *,' --> 1'
		  
         do iii=4,4  
         ikeygene=iii

c		 do iddd=1,ndelta
		 do iddd=6,6
         
		 delta=ardelta(iddd)
		 
        call bhradgen(ebeam,x,q2,t,phi,vv2cut,delta,egamma,thetag,phig,ichnnel,ipol,iappr,hel,nevent,ikeygene,stot)

		enddo
		
         enddo


         enddo
111		 continue
		 
		 
        end



       	  subroutine bhradgen(ebeam,xff,q2ff,tff,phiff,vv2cut,delta,egamma,thetag,
     +  phig,ichnnel,ipolff,iapprff,helff,nevent,ikeygene,stot)
!
!         lab system; OZ is along \vec q, OXZ -- scattering plane
!
!          ebeam -- beam energy
!          xff (or x), q2ff (or q2), tff (ot r) -- kinematic invariants   
!          phiff (phi) -- phi of final proton
!          egamma, thetag, phig -- variables of simulated (second) photon, if ichnnel>1
!          ichnnel -- channel of scattering: 1-no second photon; 2-second photon along initial electron, 3-second photon along final electron
!          ipolff (or ipol) -- target polarization parameter (0-unpolarized, 1-longitudinally polarized alog \vec q)  
!          iapprff (or iappr) --  1-BH only, 2-BMK
!          helff (or hel) --  hel=1(0) to (not) include polarized part   
!          nevent -- the required number of simulated events for given kinematical point (i.e., for fixed x,q2,t,phi)  
!
        implicit none
        integer*4 nevent,ikeygene,iacc
        integer*4 nzd,nzdphi,ichnnel,ipolff,iapprff,helff,iaddcontr_p,iaddcontr_s
	real*8 ebeam,s,q2,x,t,phi,sx,xx,mp,mp2,ml,ml2,pi,alpha,barn,
     +	sborn,siborn,sig0,sirad,siradtest,xi,z1m,z2m,z2cur,z1,z2,lll,be,bb,tcol,
     +  tmin,stot,vv2cut,vmax,vv2calc,vv2min_s,vv2min_p
     +  ,sitottest,sitottest_rad,sitottest_non
	real*8 tsp,tpp,zspeak,zppeak,ran2
	real*8 onmz1,onmz2,dz1,dz2,w2,t1lim,t2lim,lay,lat,layz,az,cz,
     +	sz,sqdz,cpr,spr,sbr,cphbr,sphbr,xsp
        real*8 epr,epst,egcut,z1mc,z2mc,sig0s,sig0p,sum1,sum2,vacpol,
     +	egsim,u0,w0,lauw,eta,etaz(3),phis,sum1delta,sum2delta
        real*8 sirad0,sirad000,siradsin,siradsin2,siradcos,siradcos2,
     +	siradcos3,sibor0,sibor000,siborsin,siborsin2,siborcos,siborcos2,
     +  siborcos3,siradadd,llog,sum1phi,sum2phi,sum1test,sum2test
        real*8 sum1add,sum2add,sum1tde,sum2tde,egamma,random1,random2,
     +	thetag,phig
	 real*8 fracacc_s,naccev_s,ntotev_s,fracacc_p,naccev_p,ntotev_p

        real*8 xff,q2ff,phiff,tff
        parameter(nzd=55)
        parameter(nzdphi=5)
        real*8 deltaz1,deltaz2,ikeymc,delta,fd1(0:nzd),fd2(0:nzd),fd1phi(0:nzdphi),fd2phi(0:nzdphi),
     +	zd1(0:nzd),zd2(0:nzd),zd1phi(0:nzdphi),zd2phi(0:nzdphi),probn,probs,probp
	real*4 urand
	integer*4 ipg,in,i,nev,iy,itkey,itpoi,iepoi,ixq2poi,nepoi,nxq2poi,
     +	ntpoi,nphipoi,
     +	ikeyphiint,isamp
	integer*4 ikeydd,ikeyfo,ipol,iappr,hel,ntreg1,ntreg2,ntreg3,
     +	nttot,nphipoi1,nphipoi2,izd
        real*8 v2p,v2m,a1,a2s,a2p,v2s1,v2s2,v2p1,v2p2,dds,ddp,sp  
     
	common/const/alpha,barn,mp,mp2,ml2,ml,pi
        common/kinpoi/s,q2,t,phi,x,xx,sborn,eta(3),ipol
        common/iappr/iappr,hel
        data alpha/0.729735d-2/,barn/0.389379d9/,ml2/0.261112d-6/,
     +	mp/0.938272d0/,iy/12/

        ikeyfo=ikeygene

        x=xff
        q2=q2ff
        t=tff
        phi=phiff      
        ipol=ipolff 
        iappr=iapprff
        hel=helff
        ml=0.000511D0!sqrt(ml2)
        mp2=mp**2
	pi=3.1415926536D0!atan(1d0)*4d0

	 s=2d0*mp*ebeam

         lll=log(q2/ml2)
         be=2d0*alpha/pi*(LLL-1d0)
         bb=be/2d0
	 epr=ebeam-q2/2d0/mp/x
	 sx=q2/x
	 xx=s-sx
	 sp=s+xx
	 w2=sx-q2+mp2
  	 lay=sx**2+4d0*mp2*q2
	 t2lim=-0.5d0*((sx-q2)*(sx+sqrt(lay))+2d0*mp2*q2)/w2
	 t1lim=-0.5d0*((sx-q2)*(sx-sqrt(lay))+2d0*mp2*q2)/w2
         tsp=-Q2*xx/(s-q2)
         tpp=-Q2*s/(xx+q2)

 	 lat=t*(t-4d0*mp2)
	 xi=sqrt(lat/lay)

         if(ipol.eq.0)then
           eta(1)=0d0
           eta(2)=0d0
           eta(3)=0d0
         endif
         if(ipol.eq.1)then
           eta(1)=0d0
           eta(2)=0d0
           eta(3)=1d0
         endif
         if(ipol.ge.2)then
           phis=(ipol-2)*pi/2d0
           eta(1)=cos(phis)
           eta(2)=sin(phis)
           eta(3)=0d0
         endif

	 sborn=siborn(s,q2,x,t,cos(phi),sin(phi),1,ipol,eta)

	  vmax=(sqrt(lay*lat)+sx*t)/2d0/mp2-q2+t

      vv2calc=min(vv2cut,vmax)



        a1=4.*mp*cos(phi)*sqrt(q2*(s*xx-mp2*q2))
       v2p=(sx*t+sqrt(lay*lat))/2d0/mp2-q2+t
       v2m=(sx*t-sqrt(lay*lat))/2d0/mp2-q2+t
       a2s=q2*(sp*(sx+2.*t)-lay)-t*(lay+sp*sx)
       a2p=q2*(sp*(sx+2.*t)+lay)+t*(lay-sp*sx)
        dds=a1**2*v2p*v2m+a2s**2
        ddp=a1**2*v2p*v2m+a2p**2
      zspeak=1d0-4.*v2p*v2m*lay/((v2m-v2p)*sqrt(dds)+a2s*(v2p+v2m)+2*v2p*v2m*(lay+sp*sx)) 
      zppeak=1d0-4.*v2p*v2m*lay/((v2m-v2p)*sqrt(ddp)+a2p*(v2p+v2m)+2*v2p*v2m*(lay+sp*sx)) 
        v2s1=v2p*dds/(a1**2*v2p**2+a2s**2)
        v2s2=2d0/((1./v2m+1./v2p)+a2s*(1./v2p-1./v2m)/sqrt(dds))
        v2p1=v2p*ddp/(a1**2*v2p**2+a2p**2)
        v2p2=2d0/((1./v2m+1./v2p)+a2p*(1./v2p-1./v2m)/sqrt(ddp))
        z1m=1d0-2d0*vv2calc*lay/(a1*sqrt((v2p-vv2calc)*(vv2calc-v2m))+a2s+vv2calc*(sp*sx+lay))
        z2m=1d0-2d0*vv2calc*lay/(a1*sqrt((v2p-vv2calc)*(vv2calc-v2m))+a2p+vv2calc*(sp*sx+lay))

       if(vv2calc.le.v2s2.or.cos(phi).ge.0d0)then
	  iaddcontr_s=0
      z1mc=z1m
      else
	  iaddcontr_s=1
      z1mc=zspeak
      endif
	  
      if(vv2calc.le.v2p2.or.cos(phi).ge.0d0)then
	  iaddcontr_p=0
      z2mc=z2m
      else
	  iaddcontr_p=1
      z2mc=zppeak
      endif


         deltaz1=min(delta/ebeam,1d0-z1m)  
         deltaz2=min(delta/(epr+delta),1d0-z2m)

		 
         do izd=0,nzd
           zd1(izd)=z1mc+(1d0-deltaz1-z1mc)*izd/nzd
           zd2(izd)=z2mc+(1d0-deltaz2-z2mc)*izd/nzd
           fd1(izd)=0d0
           fd2(izd)=0d0
         enddo 

         do izd=0,nzdphi
           zd1phi(izd)=z1mc+(z1m-z1mc)*izd/nzdphi
           zd2phi(izd)=z2mc+(z1m-z2mc)*izd/nzdphi
           fd1phi(izd)=0d0
           fd2phi(izd)=0d0
         enddo 

          i=ikeyfo

	  nev=2**i
	  sum1=0d0
      sum2=0d0
	  sum1delta=0d0
      sum2delta=0d0
	  sum1phi=0d0
      sum2phi=0d0
            sborn=siborn(s,q2,x,t,cos(phi),sin(phi),1,ipol,eta) 
	 do in=1,nev
c 	  z1=z1mc+(1d0-deltaz1-z1mc)*urand(iy)
 	  z1=z1mc+(1d0-1e-16-z1mc)*urand(iy)
          z2=1.d0
	    layz=(z1*s-xx/z2)**2+4d0*mp2*z1*q2/z2
	    az=-((z1*s-xx/z2)*t+2d0*mp2*(t-z1*q2/z2))/sqrt(lat*layz)
	    cz=(sx*(z1*s-xx/z2)+2d0*(1d0/z2+z1)*mp2*q2)/sqrt(lay*layz)
	    sz=2d0*(1d0/z2-z1)*mp*sqrt(q2*(s*xx-mp2*q2))/sqrt(lay*layz)
	    sqdz=sqrt(cz**2+sz**2*cos(phi)**2-az**2)
	    cpr=(az*cz+sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    spr=(sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    sbr=sqrt(1d0-(cpr*cz-spr*sz*cos(phi))**2)	
	    cphbr=(cz*spr*cos(phi)+sz*cpr)/sbr
            sphbr=spr*sin(phi)/sbr 
	    xsp=z1*q2/(z1*z2*s-xx)
            etaz(1)=cz*eta(1)+sz*eta(3)  
            etaz(2)=eta(2)
            etaz(3)=-sz*eta(1)+cz*eta(3)
	    xsp=z1*q2/(z1*z2*s-xx)
c	    xb1=z1/(1d0/xb-s*(1d0-z1)/Q2)
	    sig0s=siborn(z1*s,z1*q2,xsp,t,cphbr,sphbr,1,ipol,etaz)
		if(z1.le.1d0-deltaz1)then
	    sum1add=(1d0-z1mc)*(1d0+z1**2)
     +	    *(spr/sqdz*(xsp/x)**2*sig0s   )/(1d0-z1)		
ccc     +	    *(spr/sqdz*(xsp/x)**2*sig0s-sborn)/(1d0-z1)		
     	 sum1=sum1+sum1add/dble(nev)
            do izd=1,nzd
             if(z1.le.zd1(izd))fd1(izd)=fd1(izd)+sum1add/dble(nev)
            enddo
        else
	    sum1add=(1d0-z1mc)*(1d0+z1**2)
     +	    *(spr/sqdz*(xsp/x)**2*sig0s-sborn)/(1d0-z1)		
	    sum1delta=sum1delta+sum1add/dble(nev)		
		endif	

		if(z1.lt.z1m.and.iaddcontr_s.eq.1)then 
	    cpr=(az*cz-sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    spr=(-sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    sbr=sqrt(1d0-(cpr*cz-spr*sz*cos(phi))**2)	
	    cphbr=(cz*spr*cos(phi)+sz*cpr)/sbr
        sphbr=spr*sin(phi)/sbr 
		sig0s=siborn(z1*s,z1*q2,xsp,t,cphbr,sphbr,1,ipol,etaz)
	    sum1add=(1d0-z1mc)*(1d0+z1**2)*(spr/sqdz*(xsp/x)**2*sig0s    )/(1d0-z1)		
	    sum1phi=sum1phi+sum1add/dble(nev) 
            do izd=1,nzdphi
             if(z1.le.zd1phi(izd))fd1phi(izd)=fd1phi(izd)+sum1add/dble(nev)
            enddo
		endif
		
		
          z1=1.d0
c	  z2=z2mc+(1d0-deltaz2-z2mc)*urand(iy)
	  z2=z2mc+(1d0-1e-16-z2mc)*urand(iy)
	    layz=(z1*s-xx/z2)**2+4d0*mp2*z1*q2/z2
	    az=-((z1*s-xx/z2)*t+2d0*mp2*(t-z1*q2/z2))/sqrt(lat*layz)
	    cz=(sx*(z1*s-xx/z2)+2d0*(1d0/z2+z1)*mp2*q2)/sqrt(lay*layz)
	    sz=2d0*(1d0/z2-z1)*mp*sqrt(q2*(s*xx-mp2*q2))/sqrt(lay*layz)
	    sqdz=sqrt(cz**2+sz**2*cos(phi)**2-az**2)
	    cpr=(az*cz+sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    spr=(sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    sbr=sqrt(1d0-(cpr*cz-spr*sz*cos(phi))**2)	
	    cphbr=(cz*spr*cos(phi)+sz*cpr)/sbr
            sphbr=spr*sin(phi)/sbr 
            etaz(1)=cz*eta(1)+sz*eta(3)  
            etaz(2)=eta(2)
            etaz(3)=-sz*eta(1)+cz*eta(3)
	    xsp=z1*q2/(z1*z2*s-xx)
c	    xb2=1d0/(1d0/xb-(1d0-z2)*s/Q2)
            sig0p=siborn(s,q2/z2,xsp,t,cphbr,sphbr,1,ipol,etaz)
		if(z2.le.1d0-deltaz2)then
            sum2add=(1d0-z2mc)*(1d0+z2**2)*(spr/sqdz
     +	    *(xsp/x)**2*sig0p/z2      )/(1d0-z2)		
ccc     +	    *(xsp/x)**2*sig0p/z2-sborn)/(1d0-z2)		
	    sum2=sum2+sum2add/dble(nev)
            do izd=1,nzd
             if(z2.le.zd2(izd))fd2(izd)=fd2(izd)+sum2add/dble(nev)
            enddo
		else
            sum2add=(1d0-z2mc)*(1d0+z2**2)*(spr/sqdz
     +	    *(xsp/x)**2*sig0p/z2-sborn)/(1d0-z2)		
	    sum2delta=sum2delta+sum2add/dble(nev)
        endif		
		
        if(z2.lt.z2m.and.iaddcontr_p.eq.1)then 
   	      cpr=(az*cz-sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	      spr=(-sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	      sbr=sqrt(1d0-(cpr*cz-spr*sz*cos(phi))**2)	
	      cphbr=(cz*spr*cos(phi)+sz*cpr)/sbr
          sphbr=spr*sin(phi)/sbr 
          sig0p=siborn(s,q2/z2,xsp,t,cphbr,sphbr,1,ipol,etaz)
          sum2add=(1d0-z2mc)*(1d0+z2**2)*(spr/sqdz*(xsp/x)**2*sig0p/z2      )/(1d0-z2)		
   	      sum2phi=sum2phi+sum2add/dble(nev)
            do izd=1,nzdphi
             if(z2.le.zd2phi(izd))fd2phi(izd)=fd2phi(izd)+sum2add/dble(nev)
            enddo         
		 endif

		
	 enddo	


		   

          LLog= log(q2/ml2)-1d0
 		  

	 
          sirad=alpha    /Pi*vacpol(-t)*sborn
     &  	+alpha/2d0/pi*LLog*(sum1+sum2+sum1delta+sum2delta+sum1phi+sum2phi+
     &   (3d0+2.d0*log(deltaz1)+2d0*log(deltaz2)-2d0*deltaz1+0.5d0*deltaz1**2-2d0*deltaz2+0.5d0*deltaz2**2)*sborn)

	 
c          sirad= 
c     &  	+alpha/2d0/pi*(log(q2/ml2)    )
c     &          *(sum1+sum1delta+sum2+sum2delta+(z1m*(z1m+2d0)/2d0+z2m*(z2m+2d0)/2d0
c     &          +2d0*log(1d0-z1m)+2d0*log(1d0-z2m))*sborn)
	 
	 
	 
c	            write(*,'(4f8.3,4g17.8)')sirad/sborn-1.,sum1,sum1delta,sum1+sum1delta,sum2,sum2delta,sum2+sum2delta
	 
          probn=((1d0+alpha/Pi*vacpol(-t))*sborn
     &	  +alpha/2d0/pi*LLog*(sum1delta+sum2delta+(3d0+2.d0*log(deltaz1)+2d0*log(deltaz2)-2d0*deltaz1+0.5d0*deltaz1**2-2d0*deltaz2+0.5d0*deltaz2**2)*sborn))/(sborn+sirad)

          probs=
     &  	alpha/2d0/pi*LLog
     &          *(sum1+sum1phi)/(sborn+sirad)

          probp=
     &  	alpha/2d0/pi*LLog
     &          *(sum2+sum2phi)/(sborn+sirad)

c          write(47,'(g8.1,4f7.3,2g12.4,5f8.2)')delta,x,q2,t,phi,
c     1	  alpha,pi,vacpol(-t),1d2*alpha/Pi*vacpol(-t)
c	 2   ,1d2*(1.d0/(1d0-0.5*alpha/Pi*vacpol(-t))**2-1d0)
	 
	 
c	      write(*,'(6g12.4)')sum1,sum2,sum1delta,sum2delta,sum1phi,sum2phi
	 
	 
          write(66,'(g8.1,4f7.3,3f8.2,2f10.2,2f8.2)')delta,x,q2,t,phi,
     1	  1d2*alpha/2d0/pi*LLog*3d0
     2	 ,1d2*alpha/Pi*vacpol(-t)
     3   ,1d2*alpha/2d0/pi*LLog*(2.d0*log(deltaz1)+2d0*log(deltaz2)-2d0*deltaz1+0.5d0*deltaz1**2-2d0*deltaz2+0.5d0*deltaz2**2)
     4   ,1d2*alpha/2d0/pi*LLog*(sum1+sum2+sum1phi+sum2phi)/sborn
     5   ,1d2*alpha/2d0/pi*LLog*(sum1delta+sum2delta)/sborn
     6   ,1d2*sirad/sborn
     6   ,1d2*siradtest/sborn
		  

c          write(51,'(5f7.3,4f8.2,3g12.4)')delta*1d3,x,q2,t,phi
          write(*,'(5f7.3,3g12.4)')delta*1d3,x,q2,t,phi, probn, probp, probs
c	 6   ,1d2*siradtest/sborn
c     6   ,1d2*(sitottest/sborn-1)
cc	 6   ,1d2*(sitottest_rad/sborn)
cc	 6   ,1d2*(sitottest_non/sborn-1)
c	 7   ,1d2*(exp(alpha/pi*LLog*log(deltaz1*deltaz2))-1d0-alpha/pi*LLog*log(deltaz1*deltaz2))


		  
		  
		  ntotev_s=0.
		  ntotev_p=0.
		  naccev_s=0.
		  naccev_p=0.
          random1=urand(iy)          
          if(random1.le.probn)then 
             ichnnel=1
             egamma=0d0 
             thetag=0d0
             phig=0d0 
          elseif(random1.le.probn+probs)then
             ichnnel=2
             ran2=(random1-probn)/probs
			 if(ran2.le.sum1/(sum1+sum1phi))then 
			 random2=ran2*(sum1+sum1phi)   !*sum1/sum1
             do izd=1,nzd
                if(random2.lt.fd1(izd).and.random2.ge.fd1(izd-1))then
                   z1=zd1(izd-1)+(zd1(izd)-zd1(izd-1))*(random2
     +		   -fd1(izd-1))/(fd1(izd)-fd1(izd-1))
                   Egamma=(1.-z1)*Ebeam 
                endif   
             enddo
             else
			 random2=ran2*(sum1+sum1phi)-sum1   !(ran2-sum1/(sum1+sum1phi))*(sum1+sum1phi)/sum1phi*sum1phi
             do izd=1,nzdphi
                if(random2.lt.fd1phi(izd).and.random2.ge.fd1phi(izd-1))then
                   z1=zd1phi(izd-1)+(zd1phi(izd)-zd1phi(izd-1))*(random2
     +		   -fd1phi(izd-1))/(fd1phi(izd)-fd1phi(izd-1))
                   Egamma=(1.-z1)*Ebeam 
                endif   
             enddo		 
			 endif
             phig=0d0
             thetag=acos((s*sx+2*mp2*q2)/sqrt(lay)/s)
		  else
             ichnnel=3
             ran2=(random1-probn-probs)/probs
			 if(ran2.le.sum2/(sum2+sum2phi))then 
			 random2=ran2*(sum2+sum2phi)   !*sum2/sum2
             do izd=1,nzd
                if(random2.lt.fd2(izd).and.random2.ge.fd2(izd-1))then
                   z2=zd2(izd-1)+(zd2(izd)-zd2(izd-1))
     +		   *(random2-fd2(izd-1))/(fd2(izd)-fd2(izd-1))
                   Egamma=(1.-z2)/z2*Epr 
                endif   
             enddo
             else
			 random2=ran2*(sum2+sum2phi)-sum2   !(ran2-sum2/(sum2+sum2phi))*(sum2+sum2phi)/sum2phi*sum2phi
             do izd=1,nzdphi
                if(random2.lt.fd2phi(izd).and.random2.ge.fd2phi(izd-1))then
                   z2=zd2phi(izd-1)+(zd2phi(izd)-zd2phi(izd-1))
     +		   *(random2-fd2phi(izd-1))/(fd2phi(izd)-fd2phi(izd-1))
                   Egamma=(1.-z2)/z2*Epr 
                endif   
             enddo		 
			 endif
             phig=0d0
             thetag=acos((xx*sx-2*mp2*q2)/sqrt(lay)/xx)
          endif   
c           write(71,*)ichnnel,Egamma,thetag

          if(ichnnel.gt.1)then
		  call cutacc(s,xx,q2,phi,t,ichnnel,egamma,iacc)
		  if(ichnnel.eq.2)then
		  if(iacc.eq.1)naccev_s=naccev_s+1. 
		  ntotev_s=ntotev_s+1.
		  endif
		  if(ichnnel.eq.3)then
		  if(iacc.eq.1)naccev_p=naccev_p+1. 
		  ntotev_p=ntotev_p+1.
          endif

		  endif
      print *, ichnnel, thetag, egamma
		  




      fracacc_s=naccev_s/ntotev_s
      fracacc_p=naccev_p/ntotev_p
	stot=(sborn+sirad)*(probn+fracacc_s*probs+fracacc_p*probp)
	  
	end
      
	  
	  subroutine cutacc(s,xx,q2,phi,t,ichnnel,egamma,iacc)
       implicit none
	   real*8 s,xx,q2,phi,t,egamma,alpha,barn,mp,mp2,ml2,ml,pi	
	   real*8 missingmass2,complanarity,complanarity0,ptfin,ang_m_c,ptfin2
	   real*8 sx,aly,sqly,cos1,sin1,cos2,sin2,e1,e2,epr,ppr,costp,sintp,scalarpr
	   Integer*4 ichnnel,iacc,i,iacctt
	   	  common/const/alpha,barn,mp,mp2,ml2,ml,pi	
	   real*8 k1(0:3),k2(0:3),k(0:3),p1(0:3),p2(0:3),k_3kin(0:3),k_4kin(0:3),kgene(0:3),fin(0:3)

	   sx=s-xx
	   aly=sx**2+4.*mp2*q2
	   sqly=sqrt(aly)

	   cos1=(s*sx+2.*mp2*q2)/(s*sqly)
	   sin1=sqrt(1.-cos1**2)
	   e1=s/(2.*mp)	   
       k1(0)=e1   
       k1(1)=e1*sin1
       k1(2)=0.
       k1(3)=e1*cos1
	   
	   cos2=(xx*sx-2.*mp2*q2)/(xx*sqly)
	   sin2=sqrt(1.-cos2**2)
	   e2=xx/(2.*mp)	   
       k2(0)=e2   
       k2(1)=e2*sin2
       k2(2)=0.
       k2(3)=e2*cos2
       
	   p1(0)=mp
	   p1(1)=0.
	   p1(2)=0.
	   p1(3)=0.

	   epr=mp-t/(2.*mp)
       ppr=sqrt(epr**2-mp2)
	   costp=(2.*mp2*q2-2.*mp2*t-t*sx)/(sqly*sqrt(t*(t-4.*mp2)))
	   sintp=sqrt(1.-costp**2)
	   p2(0)=epr
	   p2(1)=ppr*sintp*cos(phi)
	   p2(2)=ppr*sintp*sin(phi)
	   p2(3)=ppr*costp
	   
	   kgene(0)=egamma
	   if(ichnnel.eq.2)then
       kgene(1)=egamma*sin1
       kgene(2)=0.
       kgene(3)=egamma*cos1
	   endif

	   if(ichnnel.eq.3)then
       kgene(1)=egamma*sin2
       kgene(2)=0.
       kgene(3)=egamma*cos2
	   endif
	   
	   do i=0,3
	   k_3kin(i)=k1(i)+p1(i)-k2(i)-p2(i)
	   k_4kin(i)=k1(i)+p1(i)-k2(i)-p2(i)-kgene(i)
	   fin(i)=k2(i)+p2(i)+k_4kin(i)
	   enddo
       
c	   write(*,*)' test me2 ',scalarpr(k1,k1),scalarpr(k2,k2),ml2
c	   write(*,*)' test mp2 ',scalarpr(p1,p1),scalarpr(p2,p2),mp2
c	   write(*,*)' test lam ',scalarpr(k_3kin,k_3kin),scalarpr(k_4kin,k_4kin)
c	   write(*,*)' test q2  ',2.*scalarpr(k1,k2),q2
c	   write(*,*)' test s   ',2.*scalarpr(k1,p1),s
c	   write(*,*)' test xx  ',2.*scalarpr(k2,p1),xx
c	   write(*,*)' test t   ',2.*mp2-2*scalarpr(p1,p2),t

       missingmass2=2.*scalarpr(kgene,p2) 
	   complanarity=abs(pi-acos((p2(1)*k_4kin(1)+p2(2)*k_4kin(2))/(sqrt(p2(1)**2+p2(2)**2)*sqrt(k_4kin(1)**2+k_4kin(2)**2))))*180./pi
	   complanarity0=abs(pi-acos((p2(1)*k_3kin(1)+p2(2)*k_3kin(2))/(sqrt(p2(1)**2+p2(2)**2)*sqrt(k_3kin(1)**2+k_3kin(2)**2))))*180./pi
       ptfin2=(fin(1)**2+fin(2)**2+fin(3)**2-((fin(1)*k1(1)+fin(2)*k1(2)+fin(3)*k1(3))/e1)**2)
       ang_m_c=acos((k_3kin(1)*k_4kin(1)+k_3kin(2)*k_4kin(2)+k_3kin(3)*k_4kin(3))/sqrt(k_3kin(1)**2+k_3kin(2)**2+k_3kin(3)**2)/sqrt(k_4kin(1)**2+k_4kin(2)**2+k_4kin(3)**2))*180./pi

	   
	   iacc=1
	   iacctt=10000
	   if(complanarity.gt.1.47313d0)iacctt=iacctt+100 
	   if(ptfin2.gt.0.0921189d0**2)iacctt=iacctt+10 
	   if(ang_m_c.gt.1.00133d0)iacctt=iacctt+1 
	   if(iacctt.gt.10000)iacc=0
	   
c      write(*,'(i4,6f8.4,i7)')ichnnel,egamma,missingmass2,complanarity,complanarity0,1000*ptfin2,ang_m_c,iacctt	   

	   
*here are the experimental cuts used in my analysis.
*1) Missing mass squared of the (e,p) system :
*-0.226248 < MM2_ep < 0.198154 GeV2.
*2) Coplanarity, angle between the (virtual photon, proton) and (virtual photon, photon) planes : -1.47313 *< delta_phi < 1.26223 degrees.
*3) Missing transverse momentum of the (e,p,gamma) system : sqrt(p_x*p_x+p_y*p_y) < 0.0921189 GeV.
*4) Angle between predicted and measured photon : theta_gamma_X < 1.00133 degrees.
*Let me know if you need any additional information.
	   
	   end

	   real*8 function scalarpr(aaa,bbb)
	   real*8 aaa(0:3),bbb(0:3)
	   scalarpr=aaa(0)*bbb(0)-aaa(1)*bbb(1)-aaa(2)*bbb(2)-aaa(3)*bbb(3)
	   end
	   
	  
      subroutine dzsub(vv,q2,onmz,ipo,dzfun)
      implicit real*8(a-h,k-m,o-z)
	common/const/alpha,barn,mp,mp2,ml2,ml,pi

      lll=log(q2/ml2)
      be=2d0*alpha/pi*(LLL-1d0)
      z=1d0-onmz


      cc0=(1d0+3d0/8d0*be-be**2/48d0*(lll/3d0+pi**2-47d0/8d0))
      dga=0.5d0*be   !*(1d0-z)**(be/2d0-1d0)
     . *cc0
       if(onmz.gt.0d0)then
      zbzb=onmz**(1d0-be/2d0)
       dga=dga+zbzb*(
     . -be/4d0*(1d0+z)+be**2/32d0*(-4d0*(1d0+z)*log(onmz)
     . -(1d0+3d0*z**2)/onmz*log(z)-5d0-z) )



      endif

       if(1d0-z-4d0*mp*ml/vv.gt.0d0)then
	 zer=max(1d0-z-4d0*mp*ml/vv,0d0)
	 lll1=lll+log(onmz**2)
	 ll53=lll1-5d0/3d0
	 deen=zbzb*alpha**2/pi**2/12d0*zer**(be/2d0)/onmz
     .	  *ll53**2*(1d0+z**2+be/6d0*ll53)
       else
	 deen=0.d0
       endif
      if(ipo.eq.1)then
      dees=zbzb*alpha**2/4d0/pi**2*lll**2*(2d0*(1d0-z**3)/z+.5d0*(1d0-z)
     . +(1d0+z)*log(z))
      elseif(ipo.eq.2)then
      dees=zbzb*alpha**2/4d0/pi**2*lll**2*(5d0*(1d0-z)/2d0
     . +(1d0+z)*log(z))
      endif

      dzfun=dga+deen+dees
***	 dzfun=dzfun*(1d0-z)**(1d0-be/2d0)

      end


	double precision function siborn(s,q2,x,t,cphi,sphi,ia,ipol,eta)
        implicit none
	real*8 s,q2,x,t,cphi,sphi,vacpol,alphar,alpha,barn,mp2,ml,mp,
     +	ml2,pi,f1,f2,u0,w0,tt1,tt2,lauw,lay,xx,sx,sp,w2,tmin,tmax,
     +  factorvac
	real*8 tt3,tt4,f3,f4,phis,etak1,etak2,etap2,eta(3),sqly,sqlsxq
        real*8 tt10,tt20,tt30,tt40,tt1m,tt2m,tt3m,tt4m
        real*8 phiphi
	integer*4 ia,ipol
	integer*4 iappr,hel
	common/const/alpha,barn,mp,mp2,ml2,ml,pi
        common/iappr/iappr,hel

          if(iappr.eq.2)then
            if(sphi.ge.0d0)phiphi=acos(min(1d0,max(-1d0,cphi)))
            if(sphi.lt.0d0)phiphi=2.*pi-acos(min(1d0,max(-1d0,cphi)))
            call bmkxsec(3, 3, 1, -1, 1,hel, s/2./mp, x, Q2, t,  pi-phiphi,siborn,eta,ipol,3)
c            call bmkxsec2(x, Q2, t, pi-phiphi,pi-phiphi,siborn)
            return   
          endif
        
        
	return
	end
      

        double precision function peak(z)
        implicit none
        real*8 z,lay,layz,az,cz,sz,sqdz,cpr,spr,sbr,cphbr,sphbr,xsp,sigp,z1,z2,lat
        real*8 siborn,sigsp,s,q2,t,phi,x,xx,sborn,sx,alpha,barn,mp,mp2,ml2,ml,pi,eta,etaz(3)
        integer*4 isp,ipol  ,isubs
 	common/const/alpha,barn,mp,mp2,ml2,ml,pi
        common/kinpoi/s,q2,t,phi,x,xx,sborn,eta(3),ipol
        common/test/isp,isubs
         if(isp.eq.1)then
          z1=z
          z2=1d0
         endif
         if(isp.eq.2)then
          z1=1d0
          z2=z
         endif
            sx=s-xx 
 	    lat=t*(t-4d0*mp2)
            lay=sx**2+4d0*mp2*q2  
	    layz=(z1*s-xx/z2)**2+4d0*mp2*z1*q2/z2
	    az=-((z1*s-xx/z2)*t+2d0*mp2*(t-z1*q2/z2))/sqrt(lat*layz)
	    cz=(sx*(z1*s-xx/z2)+2d0*(1d0/z2+z1)*mp2*q2)/sqrt(lay*layz)
	    sz=2d0*(1d0/z2-z1)*mp*sqrt(q2*(s*xx-mp2*q2))/sqrt(lay*layz)
	    if(cz**2+sz**2*cos(phi)**2-az**2.le.0d0) then
	    peak=0d0
	    return
	    endif
	    sqdz=sqrt(cz**2+sz**2*cos(phi)**2-az**2)
c	    print*,cz**2+sz**2*cos(phi)**2-az**2,zc
c	    pause
	    cpr=(az*cz+sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    spr=(sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    sbr=sqrt(1d0-(cpr*cz-spr*sz*cos(phi))**2)	
	    cphbr=(cz*spr*cos(phi)+sz*cpr)/sbr
            sphbr=spr*sin(phi)/sbr 
            etaz(1)=cz*eta(1)+sz*eta(3)  
            etaz(2)=eta(2)
            etaz(3)=-sz*eta(1)+cz*eta(3)
	    xsp=z1*q2/(z1*z2*s-xx)
	    sigsp=siborn(z1*s,z1*q2/z2,xsp,t,cphbr,sphbr,1,ipol,etaz)
	    peak=(1d0+z**2)*(spr/sqdz*(xsp/x)**2*sigsp/z2-isubs*sborn)/(1d0-z)	
        end



	

****************** vacpol *************************************

      double precision function vacpol(t)
c contribution from vacuum polarization by leptons (suml) and hadrons (sumh)
      implicit real*8(a-h,l,m,o-z)
	common/const/alpha,barn,mp,mp2,ml2,ml,pi
c      common/cmp/pi,alpha,amp,amp2,aml,aml2,barn
c      include 'const.inc'
      dimension am2(3)
c
c    am2 : squared masses of charge leptons
c
      data am2/.26110d-6,.111637d-1,3.18301d0/

      suml=0.
      do 10 i=1,3
	 a2=2.*am2(i)
	 sqlmi=dsqrt(t*t+2.*a2*t)
	 allmi=dlog((sqlmi+t)/(sqlmi-t))/sqlmi
  10  suml=suml+2.*(t+a2)*allmi/3.-10./9.+4.*a2*(1.-a2*allmi)/3./t
      if(t.lt.1.d0)then
	aaa = -1.345d-9
	bbb = -2.302d-3
	ccc = 4.091
      elseif(t.lt.64d0)then
	aaa = -1.512d-3
	bbb =  -2.822d-3
	ccc = 1.218
      else
	aaa = -1.1344d-3
	bbb = -3.0680d-3
	ccc = 9.9992d-1
      endif
      sumh = -(aaa+bbb*log(1.+ccc*t)) *2*pi/alpha

      vacpol=suml+sumh
c        vacpol=0.
c      print *,t,vacpol,suml,sumh,-(aaa+bbb*log(1.+ccc*t)) *2*pi/alpha,aaa,bbb,ccc,pi,alpha  

      end

      subroutine bmkxsec(Ivar, IGPD, Ipn, Ich, uel, hel,   
c     &                 E, xb, Q2, del2, Phi_e, Phi_g,dsigma)
     &                 E, xb, Q2, del2, Phi_g,dsigma,eta,ipol,iborn)
        include "dvcsmom.inc"
c      double precision E, xb, Q2, del2,Phi_e,Phi_g,Phi_s,Phi_gb,dsigma
      double precision E, xb, Q2, del2,Phi_g,Phi_gb,dsigma
      double precision nu,W2,W,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
      double precision Mp, mele, pi
      real*8 alpha000,barn000,mp000,mp2000,ml2000,ml000,pi000
      common/const/alpha000,barn000,mp000,mp2000,ml2000,ml000,pi000
      common/myconst/ Mp, mele, pi
c      data Mp/0.938D0/, mele/0.000511D0/, pi/3.1415926536D0/
      integer istat,Ivar, IGPD, Ipn, Ich, hel, uel
         integer*4 ipol
         real*8 eta(3)
      real corraul,rP1,rP2
c      print *, pi000, ml000, mp000
      pi=pi000
      mele=ml000
      mp=mp000
c
c
      istat=0                    ! kinematic range OK
      xmin1 = Q2/(2D0*Mp*E)
      xmax1 = 1D0
      nu  = Q2/(2D0*Mp*xb)
      W2  = Mp**2 + 2D0*Mp*nu - Q2
      W   = sqrt(W2)
      qmod = sqrt(nu**2 + Q2)
c
      E1cm = Mp*(Mp + nu)/W
      P1cm = Mp*qmod/W
      E2cm = (W2 + Mp**2)/(2D0*W)
      P2cm = (W2 - Mp**2)/(2D0*W)
      del2max = 2D0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
      del2min = 2D0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)
c
      if( xb.le.xmin1 .or. xb.gt.xmax1 ) istat=1           !    x  out of range
      if( del2.ge.del2min .or. del2.le.del2max ) istat=2   ! delta out of range

      if (istat.eq.0) then
      call resetmom()
!      call getphoton(E,  xb,  Q2,  del2,  Phi_e,Phi_g)
      call getphoton(E,  xb,  Q2,  del2,  Phi_g)
      Phi_gb=pi - Phi_g
c      Phi_s=Phi_e
c      call bhdvcs(E,  xb,  Q2,  del2,  Phi_gb) 
      call bhdvcs(E, xb,  Q2,  del2,  Phi_s,Phi_gb,rP1,rP2) 
c
c      dsBH =hc0BH +hc1BH*cos(Phi_gb)+hc2BH*cos(2D0*Phi_gb)
c      dsBHlp =hel*help*(hc0BHlp +hc1BHlp*cos(Phi_gb))
c      dsBHtp =hel*help*cos(Phi_s)*(hc0BHtpcos +hc1BHtpcos*cos(Phi_gb)) 
c     6       +hel*hs1BHtpsin*sin(Phi_s)*sin(Phi_gb)
c      dsDVCS=hc0dvcs+hc1dvcs*cos(Phi_gb)+hs1dvcs*sin(Phi_gb)+hs2dvcs*sin(2*Phi_gb)
c      dsIunp=hc0Iunp +hel*hs1Iunp*sin(Phi_gb)+hc1Iunp*cos(Phi_gb) 
c     6     +hel*hs2Iunp*sin(2*Phi_gb)+hc2Iunp*cos(2*Phi_gb)
c      dsIlp=hel*help*hc0Ilp+hel*help*hc1Ilp*cos(Phi_gb)+help*hs1Ilp*sin(Phi_gb)
c      dsItp=hel*hc0Itpcos*cos(Phi_s)+hc0Itpsin*sin(Phi_s)
c     6     +hc1Itpsin*cos(Phi_gb)*sin(Phi_s)+hel*hc1Itpcos*cos(Phi_gb)*cos(Phi_s)
c     6     +hel*hs1Itpsin*sin(Phi_gb)*sin(Phi_s)+hs1Itpcos*sin(Phi_gb)*cos(Phi_s)


       cosphis=eta(1)*cos(Phi_gb)+eta(2)*sin(Phi_gb)
       sinphis=eta(2)*cos(Phi_gb)-eta(1)*sin(Phi_gb)
       eta3bmk=-eta(3)

      dsBH =uel*(hc0BH +hc1BH*cos(Phi_gb)+hc2BH*cos(2D0*Phi_gb))
      dsBHlp =hel*eta3bmk*(hc0BHlp +hc1BHlp*cos(Phi_gb))
      dsBHtp =hel*cosphis*(hc0BHtpcos +hc1BHtpcos*cos(Phi_gb)) 
     6       +hel*sinphis*hs1BHtpsin*sin(Phi_gb)
c      dsDVCS=hc0dvcs+hc1dvcs*cos(Phi_gb)+hs1dvcs*sin(Phi_gb)+hs2dvcs*sin(2*Phi_gb)
      dsIunp=uel*hc0Iunp +hel*hs1Iunp*sin(Phi_gb)+uel*hc1Iunp*cos(Phi_gb) 
     6     +hel*hs2Iunp*sin(2*Phi_gb)+uel*hc2Iunp*cos(2*Phi_gb)
      dsIlp=eta3bmk*(hel*hc0Ilp+hel*hc1Ilp*cos(Phi_gb)+uel*hs1Ilp*sin(Phi_gb))
      dsItp=hel*hc0Itpcos*cosphis+uel*hc0Itpsin*sinphis
     6     +uel*hc1Itpsin*cos(Phi_gb)*sinphis+hel*hc1Itpcos*cos(Phi_gb)*cosphis
     6     +hel*hs1Itpsin*sin(Phi_gb)*sinphis+uel*hs1Itpcos*sin(Phi_gb)*cosphis

*	     write(61,'(8g11.3)')phi_gb,hc0Iunp,hc1Iunp,hc2Iunp,hs1Iunp,hs2Iunp
c	     write(*,'(8g11.3)')xb,q2,del2,Phi_g ,hc0Iunp,hs1Iunp,hc1Iunp,hs2Iunp,hc2Iunp
c		 stop

	 
c      
CC       dsigma=dsBH +dsDVCS-Ich*dsIunp !+dsBHlp -Ich*dsIlp !-Ich*dsItp

      if(iborn.eq.1.and.ipol.eq.0) dsigma=2.*pi*dsBH
      if(iborn.eq.1.and.ipol.gt.0) dsigma=2.*pi*(dsBHlp + dsBHtp)
      if(iborn.eq.2.and.ipol.eq.0) dsigma=2.*pi*(-Ich*dsIunp)
      if(iborn.eq.2.and.ipol.gt.0) dsigma=2.*pi*(-Ich*dsIlp-Ich*dsItp)
      if(iborn.eq.3.and.ipol.eq.0) dsigma=2.*pi*(dsBH-Ich*dsIunp)
      if(iborn.eq.3.and.ipol.gt.0) dsigma=2.*pi*(dsBHlp + dsBHtp-Ich*dsIlp-Ich*dsItp)

c      write(*,*)dsigma,eta(1),eta(2),eta(3),dsBHlp,dsBHtp,cos(Phi_gb),sin(Phi_gb)
c      stop

cc       dsigma=2.*pi*dsBH !+dsDVCS-Ich*dsIunp !+dsBHlp -Ich*dsIlp !-Ich*dsItp
c       dsigma=2.*pi*dsIunp
c      dsigma=2.*pi*dsBH
       
c       write(*,'(6f8.3,g16.5)')E, xb, Q2, del2, Phi_e, Phi_g,dsigma


c       write(*,11) dsBH,hc0BH,hc1BH,cos(Phi_gb),hc2BH,cos(2D0*Phi_gb)
!xb,del2,Phi_g ,hc0BH,hc1BH/hc0BH,hc0BHlp/hc0BH,hc1BHlp/hc0BH,dsIunp
 11   format(3h+++ ,15G11.3)
      else
       dsigma=0
c      print *,'out of limits ',xb,xmin1,xmax1,del2,del2min,del2max,istat      
      endif
      return
      end
c

C        Ich  -  positron(electron) = +1(-1)
C int   hel  -  lepton beam polarization
C int   help -  target polarization
C double  dsigma   d^5\sig / dx dQ^2 d|\Delta^2| d\phi_e d\phi_\gamma
C                  (pb GeV^{-4})
c
      subroutine bmkxsec2(xb, Q2, del2, Phi_e, Phi_g,dsigma)
      implicit none
#include "dvcsmom.inc"
#include "dvcs.inc"
#include "ntupgdvcs.inc"
      double precision E, xb, Q2, del2,Phi_e,Phi_g,Phi_s,Phi_gb,dsigma
      double precision nu,W2,W,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
      double precision yb,ymax,ycol
      double precision  xmin1,xmax1
      double precision dsIunp,dsBHlp, dsIlp,dsBHtp,dsItp
      real corraul,rP1,rP2
C       Ivar = 3
C       Ich = -1
C       IGPD = 3
C       Ipn = 1
c
c
      istatus=0                    ! kinematic range OK
      Ed=cl_be
      E=Ed
      xmin1 = Q2/(2D0*Mp*E)
      xmax1 = 1D0
      nu  = Q2/(2D0*Mp*xb)
      W2  = Mp**2 + 2D0*Mp*nu - Q2
      W   = sqrt(W2)
      qmod = sqrt(nu**2 + Q2)
c
      E1cm = Mp*(Mp + nu)/W
      P1cm = Mp*qmod/W
      E2cm = (W2 + Mp**2)/(2D0*W)
      P2cm = (W2 - Mp**2)/(2D0*W)
      del2max = 2D0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
      del2min = 2D0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)
c
      if( xb.le.xmin1 .or. xb.gt.xmax1 ) istatus=1           !    x  out of range
      if( del2.ge.del2min .or. del2.le.del2max ) istatus=2   ! delta out of range
      yb=nu/Ed
      if(yb.gt.cl_ymax) istatus=3                               ! y<ymax
      call dvcsycol(del2,xb,Q2,ycol)

      Phi_gb=pi - Phi_g
      Phi_s=Phi_e
      call bhdvcs(xb,  Q2,  del2,  Phi_s,Phi_gb,rP1,rP2) 
      call bhdvcs(3,3,Ipn,E,  xb,  Q2,  del2, Phi_gb) 
c
      if(abs(rP1).le.cl_ycol ) istatus=4                   ! y-too big
c      if((ycol-yb).le.cl_ycol ) istatus=4                   ! y-too big

      if (istatus.eq.0) then
c
      dsBH =hc0BH +hc1BH*cos(Phi_gb)+hc2BH*cos(2D0*Phi_gb)
      dsDVCS=hc0dvcs+hc1dvcs*cos(Phi_gb)+hs1dvcs*sin(Phi_gb)+hs2dvcs*sin(2*Phi_gb)
      dsIunp=hc0Iunp +heli*hs1Iunp*sin(Phi_gb)+hc1Iunp*cos(Phi_gb) 
     6     +heli*hs2Iunp*sin(2*Phi_gb)+hc2Iunp*cos(2*Phi_gb)
c      
c
c      scale the sin\phi moment
       dsIunp=dsIunp*cl_scale
c
       if(cl_ktcor) then
         dsIunp=dsIunp*corraul(del2,Q2)     ! correct~20% of A_LU at small t
       endif
c
       dsigma=dsBH
       if(cl_bh.gt.1)  dsigma=dsigma+dsDVCS-Ich*dsIunp 

c
c      L-POL
c
       if(cl_pol.eq.1) then 
         dsBHlp =heli*helpi*(hc0BHlp +hc1BHlp*cos(Phi_gb))
         dsIlp=heli*helpi*hc0Ilp+heli*helpi*hc1Ilp*cos(Phi_gb)
     +         +helpi*hs1Ilp*sin(Phi_gb)
         dsigma=dsigma+dsBHlp -Ich*dsIlp
       endif
c
c      T-POL
c 
       if(cl_pol.eq.2) then 
         dsBHtp =heli*helpi*cos(Phi_s)*(hc0BHtpcos +hc1BHtpcos*cos(Phi_gb)) 
     6       +heli*hs1BHtpsin*sin(Phi_s)*sin(Phi_gb)
         dsItp=heli*hc0Itpcos*cos(Phi_s)+hc0Itpsin*sin(Phi_s)
     6     +hc1Itpsin*cos(Phi_gb)*sin(Phi_s)+heli*hc1Itpcos*cos(Phi_gb)*cos(Phi_s)
     6     +heli*hs1Itpsin*sin(Phi_gb)*sin(Phi_s)+hs1Itpcos*sin(Phi_gb)*cos(Phi_s)
         dsigma=dsigma+dsBHtp -Ich*dsItp
       endif
c
c
c
      else
        dsigma=0
c       print *,'out of limits ',xb,xmin1,xmax1,del2,del2min,del2max,istat      
      endif
       gwbh=dsBH
       gwdvcs=dsDVCS
       gvint=dsIunp
c       gsin=hs1Iunp/(hc0BH+hc0dvcs+hc0Iunp)
       gsin=hs1Iunp/hc0BH
       gsin2=dsigma
        gproh=helpi
        ghp1=hcp1
        ghp2=hcp2
        ghccb=hccb
        ghcci=hcci
        gh0bh=hc0BH
        gh1bh=hc1BH
        gh2bh=hc2BH
        gh0dvcs=hc0dvcs
        gh1dvcs=hc1dvcs
        ghs1dvcs=hs1dvcs
        gh0iunp=hc0Iunp
        gh1iunp=hc1Iunp
        ghs1iunp=hs1Iunp
        ghs2iunp=hs2Iunp
       if(cl_pol.eq.1) then
        gh0bhlp=hc0bhlp
        gh1bhlp=hc1bhlp
        ghs1ilp=hs1ilp
        gh1ilp=hc1ilp
        gh0ilp=hc0ilp
       else if(cl_pol.eq.2) then
        gh0bhtpcos = hc0BHtpcos
        gh1bhtpcos = hc1BHtpcos
        ghs1bhtpsin= hs1BHtpsin
        ghs1itpcos = hs1Itpcos
        ghs1itpsin = hs1Itpsin
        gh1itpcos  = hc1Itpcos
        gh1itpsin  = hc1Itpsin
        gh0itpcos  = hc0Itpcos
        gh0itpsin  = hc0Itpsin
      endif
c      print *,'bhdvcs rP1= ',rP1,rP2,dsigma,cl_ycol,ycol-yb,xb,Q2,del2
      return
      end
c

      subroutine bhdvcs(E, xbd, Q2d, del2d, phield, phigd,rP1,rP2)
      implicit double precision (A-H,O-Z)
#include "dvcsmom.inc"
#include "dvcs.inc"
c
C BH+DVCS+INT by V. Korotkov
C Tpol and Lpol added by  H.Avakian
c
C   Diff. cross-section for Bethe-Heitler, DVCS and interf.term 
C
C         d^5\sig / dx dQ^2 d|\Delta^2| d\phi_e d\phi_\gamma
C              (pb GeV^{-4})
C INPUT:
C        Ivar -  1 (BH only), 2 (DVCS + int only), 3 (BH + DVCS + int)
C        IGPD -  GPD variant
C        Ipn  -  proton(neutron) target = 1(2)
C real   E(GeV), xb, Q2(GeV^2), del2(GeV^2)(negative) 
C real   phiel(rad) - scatt. electron azimuthal angle
C real   phig(rad)  - photon azimuthal angle around momentum transfer vector q
C OUTPUT:
C real   hs*,hc* sin and cos moments of 5-fold diff. cross-section for
C                            BH, DVCS and interf. terms
*

      parameter ( alpha = 1D0/137.036D0, hc2  = 0.38938D0)
      parameter ( coeff = 1D+9*hc2*alpha**3 )
c
      real rP1,rP2
      double precision nu, k1pl, k2pl, Kfac, Jfac
      common/todvcs/ xb, yb, Q2, del2, del2min,phip,phipel,P1,P2,Kfac,Jfac,ds
*
      common/formfac/ F1pn(2), F2pn(2)
*
      xb    = xbd
      Q2    = Q2d
      del2  = del2d
*
      nu  = Q2d/(2D0*Mp*xbd)
      qmod = sqrt(nu**2 + Q2d)
      yb = nu/E
      Esc = E - nu

c
      eps = 2D0*xbd*Mp/sqrt(Q2)
      eps2=eps*eps
      qeps2=1D0 + eps2
      sqeps2=sqrt(qeps2)
      ds = coeff*(xbd*yb**2/(16D0*pi**2*Q2**2))/sqeps2
c
      del2min=-Q2*(2D0*(1D0-xb)*(1D0-sqeps2)+eps2)
      del2min=del2min/(4D0*xb*(1D0-xb)+eps2)
      tau  = del2d/(4D0*Mp**2)
      taum1=1D0-tau
      xtau=xbd*xbd/tau
      del2q2=del2d/Q2
      del2q4=del2q2*del2q2
      del2q2m1=1D0-del2q2
*
*
      phip = phigd 
      phipel = phield 
*
      call nuclFF( del2d )
*
      y1eps=1D0 - yb - yb*yb*eps2/4D0
      sqy1eps=sqrt(y1eps)
      Kfac = sqrt((-del2q2)*(1D0 - xbd)*y1eps*
     *       (1D0 - del2min/del2d)*(sqrt(1D0 + eps2) + 
     *  ((4D0*xbd*(1D0 - xbd) + eps2)/(4D0*(1D0 - xbd)))*
     *                                       ((del2d - del2min)/Q2)))
      Jfac = (1D0 - yb - yb*eps2/2D0)*(1D0 + del2q2) - 
     *       (1D0 - xbd)*(2D0 - yb)*del2q2
      P1 = -(Jfac + 2D0*Kfac*cos(phip))/(yb*(1D0 + eps2))
      P2 = 1D0 + del2q2 - P1
      rP1=P1
      rP2=P2
c
        if( Ivar.eq.1 .or. Ivar.eq.3 ) then   
       F1 = F1pn(Ipn)
       F2 = F2pn(Ipn)    
       F1_M_F = F1 + tau*F2
       F1_M_F2 = F1**2 - tau*F2**2
       F1_P_F = F1 + F2
       F1_P_F2 = F1_P_F*F1_P_F
       c01_BH = 8D0*Kfac**2*((2D0 + 3D0*eps2)*Q2*F1_M_F2/del2d +
     *                                          2D0*xbd**2*F1_P_F2)
       c02_BH = (2D0 - yb)**2*((2D0 + eps2)*F1_M_F2*
     *       ((2D0*xbd*Mp)**2*(1D0 + del2q2)**2/del2d + 
     *                         4D0*(1D0 - xbd)*(1D0 + xbd*del2q2)) +
     *     4D0*xbd**2*F1_P_F2*(xbd + 
     *                   (1D0 - xbd + eps2/2D0)*(del2q2m1)**2 -
     *                    xbd*(1D0 - 2D0*xbd)*(del2q2)**2))
       c03_BH = 8D0*(1D0 + eps2)*(1D0 - yb - yb*yb*eps2/4D0)*
     *       (2D0*eps2*(1D0 - del2d/(4D0*Mp**2))*F1_M_F2 - 
     *           xbd**2*(del2q2m1)**2*F1_P_F2)
       c0_BH = c01_BH + c02_BH + c03_BH
       c1_BH = 8D0*Kfac*(2D0 - yb)*(
     *      F1_M_F2*(4D0*(xbd*Mp)**2/del2d - 2D0*xbd - eps2) +
     *      F1_P_F2*2D0*xbd**2*(1D0 - (1D0 - 2D0*xbd)*del2q2))
       c2_BH = 8D0*(xbd*Kfac)**2*(F1_M_F2*4D0*Mp**2/del2d + 2D0*F1_P_F2)
c
c     BH-lpol part
c
      bhkin1=8D0*xbd*yb*sqeps2/taum1*F1_P_F
      c01_BHlp=(xbd/2D0*del2q2m1-tau)/2D0
      c02_BHlp=2D0-xbd-2D0*(1D0-xbd)**2*del2q2+eps2*del2q2m1
     6         -xbd*(1D0-2D0*xbd)*del2q4
      c03_BHlp=1D0-(1D0-xbd)*del2q2
      c04_BHlp=xtau/4D0*(1D0+del2q2)**2+(1D0-xbd)*(1D0+xbd*del2q2)
      c0_BHlp=c01_BHlp*c02_BHlp*F1_P_F+c03_BHlp*c04_BHlp*F1_M_F
      c0_BHlp=c0_BHlp*(2D0-yb)*bhkin1
c
      c11_BHlp=(2D0*tau-xbd*del2q2m1)*(1D0-xbd+xbd*del2q2)  
      c12_BHlp=1D0+xbd-(3D0-2D0*xbd)*(1D0+xbd*del2q2)-xtau*(1D0+del2q4)
      c1_BHlp=c11_BHlp*F1_P_F+c12_BHlp*F1_M_F
      c1_BHlp=-c1_BHlp*Kfac*bhkin1
c
c
c     BH-Tpol part
c
      c01_BHtpcos=-8D0*(2D0-yb)*yb*sqrt(q2d)/Mp*sqeps2*Kfac/sqy1eps
c
      c02_BHtpcos=xbd*(xbd*Mp)**2/q2d*del2q2m1*F1_P_F
      c03_BHtpcos=1D0-(1D0-xbd)*del2q2
      c04_BHtpcos=xtau/4D0*del2q2m1*F1+xbd/2D0*F2
      c0_BHtpcos=c01_BHtpcos*F1_P_F*(c02_BHtpcos+c03_BHtpcos*c04_BHtpcos)
c
      c11_BHtpcos=-16D0*xbd*yb*sqy1eps*Mp/sqrt(q2d)
      c11_BHtpcos= c11_BHtpcos*sqeps2*F1_P_F
      c12_BHtpcos=2D0*Kfac*Kfac*q2d/del2/y1eps
      c13_BHtpcos=xbd*del2q2m1*F1+tau*F2
      c14_BHtpcos=qeps2*xbd*del2q2m1*F1_M_F
      c1_BHtpcos=c11_BHtpcos*(c12_BHtpcos*c13_BHtpcos+c14_BHtpcos)
c
      s11_BHtpsin=16D0*yb*xbd*xbd*sqy1eps/sqrt(q2d)*Mp
      s12_BHtpsin=qeps2*sqeps2*del2q2m1*F1_P_F*F1_M_F
      s1_BHtpsin=s11_BHtpsin*s12_BHtpsin
c
       BHfact=ds/((xbd*yb*(1D0 + eps2))**2*del2d*P1*P2)
       hcP1=P1
       hcP2=P2
       hccb=BHfact
       hc0BH=c0_BH*BHfact
       hc1BH=c1_BH*BHfact
       hc2BH=c2_BH*BHfact
c
       ac0bh=c0_BH
       ac1bh=c1_BH
       ac2bh=c2_BH
c
c
       hc0BHlp=c0_BHlp*BHfact
       hc1BHlp=c1_BHlp*BHfact
c
       hc0BHtpcos=c0_BHtpcos*BHfact
       hc1BHtpcos=c1_BHtpcos*BHfact
       hs1BHtpsin=s1_BHtpsin*BHfact
c
        endif
c
        if( Ivar.eq.2 .or. Ivar.eq.3 ) then       
         call dvcsfun()
        endif
c
      return
      end


      subroutine dvcsfun()
      implicit double precision (A-H,O-Z)
#include "dvcs.inc"
#include "dvcsmom.inc"
      double precision Kfac,Jfac,Intfac
      common/todvcs/ x, y,Q2,del2,del2min,phip,phipel,P1,P2,Kfac,Jfac,ds
      common/formfac/ F1pn(2), F2pn(2)
c
      skew = x/(2D0 - x)
      if(IGPD.lt.100) then
      call amptab(skew, del2, 
     &             H1_RE, H1_IM, H1T_RE, H1T_IM,
     &             E1_RE, E1_IM, E1T_RE, E1T_IM )
      else
      call amptabgag(skew, del2, 
     &             H1_RE, H1_IM, H1T_RE, H1T_IM,
     &             E1_RE, E1_IM, E1T_RE, E1T_IM )
      endif
*
C         proton/neutron
      F1 = F1pn(Ipn)
      F2 = F2pn(Ipn)
*
      deldel    = 1D0 - del2min/del2
      deldel_sq = sqrt(deldel)
      del2m2    = -del2/Mp**2
      del2m4    = -del2m2/4D0
      delm2_sq  = sqrt(del2m2)
      cy2   = 2D0 - 2D0*y + y**2
      Intfac=ds/(x*y**3*P1*P2*(-del2))
      hcci=Intfac
*
C  DVCS
C     
      a1 = H1_RE**2 + H1_IM**2 + H1T_RE**2 + H1T_IM**2
      a2 = 2.*( H1_RE*E1_RE  +  H1_IM*E1_IM +
     &         H1T_RE*E1T_RE + H1T_IM*E1T_IM )
      a3 =  E1_RE**2 +  E1_IM**2
      a4 = E1T_RE**2 + E1T_IM**2 
      C_DVCS = ( 4D0*(1D0-x)*a1 - a2*x**2 - (x**2 + (2D0-x)**2*del2m4)*a3 
     &           - x**2*del2m4*a4 )/(2D0 - x)**2
      C_DVCS_eff = -x*C_DVCS
      c0_DVCS = 2D0*cy2*C_DVCS
      c1_DVCS = 8D0*((2D0 - y)/(2D0 - x))*C_DVCS_eff
      T_DVCS  = (c0_DVCS + Kfac*c1_DVCS*cos(phip))/(y**2*Q2)
      DVCSfac=ds/(y**2*Q2)
      hc0dvcs=c0_DVCS*DVCSfac
      hc1dvcs=Kfac*c1_DVCS*DVCSfac
      ac0dvcs=c0_DVCS
      ac1dvcs=Kfac*c1_DVCS
C
C  INTERF
C
      C_I_re = F1*H1_RE + x/(2D0-x)*(F1+F2)*H1T_RE - del2m4*F2*E1_RE
      C_I_im = F1*H1_IM + x/(2D0-x)*(F1+F2)*H1T_IM - del2m4*F2*E1_IM
      RE2    = x/(2D0-x)*(H1_RE + E1_RE) + H1T_RE
      C_I_re_eff = -x*C_I_re
      C_I_im_eff = -x*C_I_im
      b1 = (2D0 - x)*(1D0 - y) 
      b1= b1 + (2D0-y)**2/(1D0-y)*Kfac*Kfac/del2*Q2 ! old - (1D0 - x)*(2D0 - y)**2*deldel
      b2 = (1D0 - y)*x*(F1 + F2)
      c0_I = -8D0*(2D0 - y)*( b1*C_I_re - b2*RE2 )
      c1_I = -8D0*cy2*C_I_re
      s1_I =  8D0*y*(2D0 - y)*C_I_im
      c2_I = -16D0*((2D0 - y)/(2D0 - x))*C_I_re_eff
      s2_I =  16D0*(y/(2D0 - x))*C_I_im_eff
c
c
C     moments
c
      hs2Iunp=Kfac*Kfac*s2_I*Intfac
      hs1Iunp=Kfac*s1_I*Intfac
      hc2Iunp=Kfac*Kfac*c2_I*Intfac
      hc1Iunp=Kfac*c1_I*Intfac
      hc0Iunp=del2/Q2*c0_I*Intfac
c
      ac0int=del2/Q2*c0_I
      ac1int=Kfac*c1_I
      as1int=Kfac*s1_I
c
c      print *,'mysl-dvcs',hs1Iunp,hc1Iunp,hc0Iunp
C
C  LPOL
C
      C_LP_re = (F1+F2)*skew*(H1_RE+x/2D0*E1_RE)+
     6          F1*H1T_RE-skew*(x/2D0*F1+del2m4*F2)*E1T_RE
c
c
      C_LP_im = (F1+F2)*skew*(H1_IM+x/2D0*E1_IM)+
     6          F1*H1T_IM-skew*(x/2D0*F1+del2m4*F2)*E1T_IM
c
      DC_LP_re=-skew*(F1+F2)*(H1_RE+x/2D0*E1_RE+skew
     6*(H1T_RE+x/2D0*E1T_RE))
      DC_LP_im=-skew*(F1+F2)*(H1_IM+x/2D0*E1_IM+skew
     6*(H1T_IM+x/2D0*E1T_IM))
c
      yf2=(2D0-y)**2/(1D0-y)+2
      c0lp_I =-8*y*(Kfac**2*yf2*C_LP_re+(1D0-y)*(2D0-x)*del2/Q2
     6        *(DC_LP_re+C_LP_re))
      s1lp_I = 8*Kfac*cy2*C_LP_im
      c1lp_I = -8*Kfac*y*(2D0-y)*C_LP_re

C
C  TPOL
C
      xb2=x*skew
      C_TPP_re_s = (F1+F2)*(xb2*(H1_RE+x/2D0*E1_RE)+x*del2m4*E1_RE)
     6          -xb2*F1*(H1T_RE+x/2D0*E1T_RE)
      C_TPP_re_b1 =del2m4*4D0*(1D0-x)/(2D0-x)*F2*H1T_RE
      C_TPP_re_b21 =-del2m4*x*F1*E1T_RE
      C_TPP_re_b22 =-del2m4*xb2*F2*E1T_RE
      C_TPP_re_b =C_TPP_re_b1+C_TPP_re_b21+C_TPP_re_b22
c
      C_TPP_im_s = (F1+F2)*(xb2*(H1_IM+x/2D0*E1_IM)+x*del2m4*E1_IM)
     6          -xb2*F1*(H1T_IM+x/2D0*E1T_IM)
      C_TPP_im_b =del2m4*(4D0*(1D0-x)/(2D0-x)*F2*H1T_IM
     6            -(x*F1+xb2*F2)*E1T_IM)
c
c
c
      C_TPM_re_s = 1D0/(2D0-x)*(x*x*F1-(1D0-x)*4D0*del2m4*F2)*H1_RE
     6          -xb2*(F1+F2)*(H1T_RE+del2m4*E1T_RE)
      C_TPM_re_b = (del2m4*((2D0-x)*F1+xb2*F2) +xb2*F1)*E1_RE
c
      C_TPM_im_s = 1D0/(2D0-x)*(x*x*F1-(1D0-x)*4D0*del2m4*F2)*H1_IM
     6          -xb2*(F1+F2)*(H1T_IM+del2m4*E1T_IM)
      C_TPM_im_b = (del2m4*((2D0-x)*F1+xb2*F2) +xb2*F1)*E1_IM
c
c
      C_TPM_re=C_TPM_re_s+C_TPM_re_b
      C_TPP_re=C_TPP_re_s+C_TPP_re_b
      C_TPM_im=C_TPM_im_s+C_TPM_im_b
      C_TPP_im=C_TPP_im_s+C_TPP_im_b
c
c
c
      DC_TPP_re=-4D0*del2m4*(F2*H1T_RE-x/(2D0-x)*(F1+x/2D0*F2)*E1T_RE)
      DC_TPM_re=4D0*del2m4*(F2*H1_RE-F1*E1_RE)
c
      DC_TPP_im=-4D0*del2m4*(F2*H1T_IM-x/(2D0-x)*(F1+x/2D0*F2)*E1T_IM)
      DC_TPM_im=4D0*del2m4*(F2*H1_IM-F1*E1_IM)
c
c
c
      qm8=8D0*Mp*sqrt(1D0-y)/sqrt(Q2)
      c0tpcos_I = -qm8*Kfac*y*(((2D0-y)**2/(1D0-y)+2D0)*C_TPP_re+DC_TPP_re)
      c0tpsin_I = qm8*Kfac*(2D0-y)*((2D0-y)**2/(1D0-y)*C_TPM_im+DC_TPM_im)
c
      c1tpcos_I = -qm8*y*(2D0-y)*C_TPP_re
      c1tpsin_I = qm8*cy2*C_TPM_im
c
      s1tpcos_I = qm8*cy2*C_TPP_im
      s1tpsin_I = -qm8*y*(2D0-y)*C_TPM_re

      hc0Itpcos = c0tpcos_I*Intfac
      hc0Itpsin = c0tpsin_I*Intfac
      hc1Itpcos = c1tpcos_I*Intfac
      hc1Itpsin = c1tpsin_I*Intfac
      hs1Itpcos = s1tpcos_I*Intfac
      hs1Itpsin = s1tpsin_I*Intfac
c
c 
      hs1Ilp = s1lp_I*Intfac
      hc1Ilp = c1lp_I*Intfac
      hc0Ilp = c0lp_I*Intfac
c
      return
      end
c
C       implicit double precision (A-H,O-Z)
C #include "dvcs.inc"
C #include "dvcsmom.inc"
C       double precision Kfac,Jfac,Intfac
C       common/todvcs/ x, y,Q2,del2,del2min,phip,phipel,P1,P2,Kfac,Jfac,ds
C       common/formfac/ F1pn(2), F2pn(2)
C c
C       real zskew,zF1,zF2,zHp,zEp,zHtp,zEtp,zCIim
C       common /zgpds/zskew,zF1,zF2,zHp,zEp,zHtp,zEtp,zCIim
C c
C       skew = x/(2D0 - x)
C       call amptab( IGPD, Ipn, skew, del2, 
C      &             H1_RE, H1_IM, H1T_RE, H1T_IM,
C      &             E1_RE, E1_IM, E1T_RE, E1T_IM )
C *
C c
C c
C C         proton/neutron
C       F1 = F1pn(Ipn)
C       F2 = F2pn(Ipn)
C c
C c
C c      print *,skew,F1,F2,del2
C       zskew=skew
C       zF1=F1
C       zF2=F2
C       zHp=H1_IM
C       zEp=E1_IM
C       zHtp=H1T_IM
C       zEtp=E1T_IM
C c
C *
C       deldel    = 1D0 - del2min/del2
C       deldel_sq = sqrt(deldel)
C       del2m2    = -del2/Mp**2
C       del2m4    = -del2m2/4D0
C       delm2_sq  = sqrt(del2m2)
C       cy2   = 2D0 - 2D0*y + y**2
C       Intfac=ds/(x*y**3*P1*P2*(-del2))
C c
C c
C C  DVCS
C C     
C       a1 = H1_RE**2 + H1_IM**2 + H1T_RE**2 + H1T_IM**2
C       a2 = 2.*( H1_RE*E1_RE  +  H1_IM*E1_IM +
C      &         H1T_RE*E1T_RE + H1T_IM*E1T_IM )
C       a3 =  E1_RE**2 +  E1_IM**2
C       a4 = E1T_RE**2 + E1T_IM**2 
C       C_DVCS = ( 4D0*(1D0-x)*a1 - a2*x**2 - (x**2 + (2D0-x)**2*del2m4)*a3 
C      &           - x**2*del2m4*a4 )/(2D0 - x)**2
C       C_DVCS_eff = -x*C_DVCS
C       c0_DVCS = 2D0*cy2*C_DVCS
C       c1_DVCS = 8D0*((2D0 - y)/(2D0 - x))*C_DVCS_eff
C       T_DVCS  = (c0_DVCS + Kfac*c1_DVCS*cos(phip))/(y**2*Q2)
C       DVCSfac=ds/(y**2*Q2)
C       hc0dvcs=c0_DVCS*DVCSfac
C       hc1dvcs=Kfac*c1_DVCS*DVCSfac
C C
C C  INTERF
C C
C c
C c
C       C_I_re = F1*H1_RE + x/(2D0-x)*(F1+F2)*H1T_RE - del2m4*F2*E1_RE
C       C_I_im1 = F1*H1_IM
C       C_I_im2 = x/(2D0-x)*(F1+F2)*H1T_IM
C       C_I_im3  = -del2m4*F2*E1_IM
      
C       C_I_im = C_I_im1 + C_I_im2 + C_I_im3
C       zCIim=-C_I_im/F1 
C       RE2    = x/(2D0-x)*(H1_RE + E1_RE) + H1T_RE
C       C_I_re_eff = -x*C_I_re
C       C_I_im_eff = -x*C_I_im
C       b1 = (2D0 - x)*(1D0 - y) 
C       b1= b1 + (2D0-y)**2/(1D0-y)*Kfac*Kfac/del2*Q2 ! old - (1D0 - x)*(2D0 - y)**2*deldel
C       b2 = (1D0 - y)*x*(F1 + F2)
C       c0_I = -8D0*(2D0 - y)*( b1*C_I_re - b2*RE2 )
C       c1_I = -8D0*cy2*C_I_re
C       s1_I =  8D0*y*(2D0 - y)*C_I_im
C       c2_I = -16D0*((2D0 - y)/(2D0 - x))*C_I_re_eff
C       s2_I =  16D0*(y/(2D0 - x))*C_I_im_eff
C c
C       hs2Iunp=Kfac*Kfac*s2_I*Intfac
C       hs1Iunp=Kfac*s1_I*Intfac
C       hc2Iunp=Kfac*Kfac*c2_I*Intfac
C       hc1Iunp=Kfac*c1_I*Intfac
C       hc0Iunp=del2/Q2*c0_I*Intfac
C c
C c      print *, 'vvv ',hs1Iunp/hc0BH,C_I_im2/C_I_im1,C_I_im3/C_I_im1
C c
C c      print *,'mysl-dvcs',hs1Iunp,hc1Iunp,hc0Iunp
C C
C C  LPOL
C C
C       C_LP_re = (F1+F2)*skew*(H1_RE+x/2D0*E1_RE)+
C      6          F1*H1T_RE-skew*(x/2D0*F1+del2m4*F2)*E1T_RE
C c
C c
C       C_LP_im = (F1+F2)*skew*(H1_IM+x/2D0*E1_IM)+
C      6          F1*H1T_IM-skew*(x/2D0*F1+del2m4*F2)*E1T_IM
C c
C       DC_LP_re=-skew*(F1+F2)*(H1_RE+x/2D0*E1_RE+skew
C      6*(H1T_RE+x/2D0*E1T_RE))
C       DC_LP_im=-skew*(F1+F2)*(H1_IM+x/2D0*E1_IM+skew
C      6*(H1T_IM+x/2D0*E1T_IM))
C c
C       yf2=(2D0-y)**2/(1D0-y)+2
C       c0lp_I =-8*y*(Kfac**2*yf2*C_LP_re+(1D0-y)*(2D0-x)*del2/Q2
C      6        *(DC_LP_re+C_LP_re))
C       s1lp_I = 8*Kfac*cy2*C_LP_im
C       c1lp_I = -8*Kfac*y*(2D0-y)*C_LP_re

C C
C C  TPOL
C C
C       xb2=x*skew
C       C_TPP_re_s = (F1+F2)*(xb2*(H1_RE+x/2D0*E1_RE)+x*del2m4*E1_RE)
C      6          -xb2*F1*(H1T_RE+x/2D0*E1T_RE)
C       C_TPP_re_b1 =del2m4*4D0*(1D0-x)/(2D0-x)*F2*H1T_RE
C       C_TPP_re_b21 =-del2m4*x*F1*E1T_RE
C       C_TPP_re_b22 =-del2m4*xb2*F2*E1T_RE
C       C_TPP_re_b =C_TPP_re_b1+C_TPP_re_b21+C_TPP_re_b22
C c
C       C_TPP_im_s = (F1+F2)*(xb2*(H1_IM+x/2D0*E1_IM)+x*del2m4*E1_IM)
C      6          -xb2*F1*(H1T_IM+x/2D0*E1T_IM)
C       C_TPP_im_b =del2m4*(4D0*(1D0-x)/(2D0-x)*F2*H1T_IM
C      6            -(x*F1+xb2*F2)*E1T_IM)
C c
C c
C c
C       C_TPM_re_s = 1D0/(2D0-x)*(x*x*F1-(1D0-x)*4D0*del2m4*F2)*H1_RE
C      6          -xb2*(F1+F2)*(H1T_RE+del2m4*E1T_RE)
C       C_TPM_re_b = (del2m4*((2D0-x)*F1+xb2*F2) +xb2*F1)*E1_RE
C c
C       C_TPM_im_s = 1D0/(2D0-x)*(x*x*F1-(1D0-x)*4D0*del2m4*F2)*H1_IM
C      6          -xb2*(F1+F2)*(H1T_IM+del2m4*E1T_IM)
C       C_TPM_im_b = (del2m4*((2D0-x)*F1+xb2*F2) +xb2*F1)*E1_IM
C c
C c
C       C_TPM_re=C_TPM_re_s+C_TPM_re_b
C       C_TPP_re=C_TPP_re_s+C_TPP_re_b
C       C_TPM_im=C_TPM_im_s+C_TPM_im_b
C       C_TPP_im=C_TPP_im_s+C_TPP_im_b
C c
C c
C c
C       DC_TPP_re=-4D0*del2m4*(F2*H1T_RE-x/(2D0-x)*(F1+x/2D0*F2)*E1T_RE)
C       DC_TPM_re=4D0*del2m4*(F2*H1_RE-F1*E1_RE)
C c
C       DC_TPP_im=-4D0*del2m4*(F2*H1T_IM-x/(2D0-x)*(F1+x/2D0*F2)*E1T_IM)
C       DC_TPM_im=4D0*del2m4*(F2*H1_IM-F1*E1_IM)
C c
C c
C c
C       qm8=8D0*Mp*sqrt(1D0-y)/sqrt(Q2)
C       c0tpcos_I = -qm8*Kfac*y*(((2D0-y)**2/(1D0-y)+2D0)*C_TPP_re+DC_TPP_re)
C       c0tpsin_I = qm8*Kfac*(2D0-y)*((2D0-y)**2/(1D0-y)*C_TPM_im+DC_TPM_im)
C c
C       c1tpcos_I = -qm8*y*(2D0-y)*C_TPP_re
C       c1tpsin_I = qm8*cy2*C_TPM_im
C c
C       s1tpcos_I = qm8*cy2*C_TPP_im
C       s1tpsin_I = -qm8*y*(2D0-y)*C_TPM_re

C       hc0Itpcos = c0tpcos_I*Intfac
C       hc0Itpsin = c0tpsin_I*Intfac
C       hc1Itpcos = c1tpcos_I*Intfac
C       hc1Itpsin = c1tpsin_I*Intfac
C       hs1Itpcos = s1tpcos_I*Intfac
C       hs1Itpsin = s1tpsin_I*Intfac
C c
C c 
C       hs1Ilp = s1lp_I*Intfac
C       hc1Ilp = c1lp_I*Intfac
C       hc0Ilp = c0lp_I*Intfac
C c
C       return
C       end

      subroutine amptabgag(skew, del2, 
     &                   H1_RE, H1_IM, H1T_RE, H1T_IM,
     &                   E1_RE, E1_IM, E1T_RE, E1T_IM )
#include "dvcs.inc"
c
      REAL xi
      REAL mt
      REAL cff_re
      REAL cff_im
      REAL v_re,v_im
c
      double precision skew, del2, H1_RE, H1_IM, H1T_RE, H1T_IM,
     &                             E1_RE, E1_IM, E1T_RE, E1T_IM 
      double precision F1pn(2), F2pn(2)
      common/formfac/ F1pn, F2pn
      real mpi/0.1396/
      data init/1/
c
      SAVE
c
      if( init .eq. 1 ) then
        init = 0
        call GPD_Init(IGPD-100)
      endif
*
      F1u = 2D0*F1pn(1) + F1pn(2)
      F1d = 2D0*F1pn(2) + F1pn(1)
      F2u = 2D0*F2pn(1) + F2pn(2)
      F2d = 2D0*F2pn(2) + F2pn(1)
c
      xi=skew
      mt=-del2
      v_re = GPD_GET_CFF(xi,mt,0)
      v_im = GPD_GET_CFF(xi,mt,3)
      H1_RE=-v_re
      H1_IM=-v_im
      v_re = GPD_GET_CFF(xi,mt,1)
      v_im = GPD_GET_CFF(xi,mt,4)
      H1T_RE=-v_re
      H1T_IM=-v_im
      v_re = GPD_GET_CFF(xi,mt,2)
      v_im = GPD_GET_CFF(xi,mt,5)
      E1_RE=-v_re
      E1_IM=-v_im
      E1T_IM=0
      E1T_RE=0
c
      return
      end
c
c
c
      subroutine amptab(skew, del2, 
     &                   H1_RE, H1_IM, H1T_RE, H1T_IM,
     &                   E1_RE, E1_IM, E1T_RE, E1T_IM )
#include "dvcs.inc"
      double precision skew, del2, H1_RE, H1_IM, H1T_RE, H1T_IM,
     &                             E1_RE, E1_IM, E1T_RE, E1T_IM 
      double precision F1pn(2), F2pn(2)
      common/formfac/ F1pn, F2pn
*
      common/retbl/ vh1ure(51,5),  vh1dre(51,5),
     &              vh1ture(51,5), vh1tdre(51,5),
     &              ve1ure(51,5),  ve1dre(51,5),
     &              ve1ture(51,5), ve1tdre(51,5)
      common/imtbl/ vh1uim(51,5),  vh1dim(51,5),
     &              vh1tuim(51,5), vh1tdim(51,5),
     &              ve1uim(51,5),  ve1dim(51,5),
     &              ve1tuim(51,5), ve1tdim(51,5)
*
      common/retbl1/ vh1urenf(51,21,2), vh1drenf(51,21,2) 
      common/imtbl1/ vh1uimnf(51,21,2), vh1dimnf(51,21,2)
*
      real mpi/0.1396/
      data init/1/
      real  h1uimag,h1dimag,h1tuimag,h1tdimag
      common /h1imag/h1uimag,h1dimag,h1tuimag,h1tdimag
*
      SAVE
*
      if( init .eq. 1 ) then
        init = 0
        call rtable
        skewmin  = 0.01D0
        skewmax  = 1.00D0
        skewminl = log10(skewmin)
        skewmaxl = log10(skewmax)
        dskewl   = (skewmaxl - skewminl)/51D0
        dlmin    = 0.01D0
        dlmax    = 1.00D0
        dlminl   = log10(dlmin)
        dlmaxl   = log10(dlmax)
        ddll     = (dlmaxl - dlminl)/20D0
      endif
*
      F1u = 2D0*F1pn(1) + F1pn(2)
      F1d = 2D0*F1pn(2) + F1pn(1)
      F2u = 2D0*F2pn(1) + F2pn(2)
      F2d = 2D0*F2pn(2) + F2pn(1)
*
      DS   = (log10(skew) - skewminl)/dskewl + 1
      IS   = int(DS)
      FDS  = DS - real(IS)
      FDS1 = 1. - FDS
*
      DT   = (log10(-del2) - dlminl)/ddll + 1
      IT   = int(DT)
      FDT  = DT - real(IT)
      FDT1 = 1. - FDT
*
      if( IGPD .ge. 6 ) then
        IGPD1 = IGPD - 4
        IGPD2 = IGPD - 5
      else
        IGPD1 = IGPD
      endif
*
        if( IGPD .le. 5 ) then
      h1ure =   vh1ure(IS,IGPD1)*FDS1 + vh1ure(IS+1,IGPD1)*FDS
      h1dre =   vh1dre(IS,IGPD1)*FDS1 + vh1dre(IS+1,IGPD1)*FDS
*
      h1uim =   vh1uim(IS,IGPD1)*FDS1 + vh1uim(IS+1,IGPD1)*FDS
      h1dim =   vh1dim(IS,IGPD1)*FDS1 + vh1dim(IS+1,IGPD1)*FDS
        else
      h1ure = (vh1urenf(IS,IT,IGPD2)*FDS1 + 
     &           vh1urenf(IS+1,IT,IGPD2)*FDS)*FDT1 +
     &        (vh1urenf(IS,IT+1,IGPD2)*FDS1 +
     &           vh1urenf(IS+1,IT+1,IGPD2)*FDS)*FDT
      h1dre = (vh1drenf(IS,IT,IGPD2)*FDS1 + 
     &           vh1drenf(IS+1,IT,IGPD2)*FDS)*FDT1 +
     &        (vh1drenf(IS,IT+1,IGPD2)*FDS1 +
     &           vh1drenf(IS+1,IT+1,IGPD2)*FDS)*FDT
*
      h1uim = (vh1uimnf(IS,IT,IGPD2)*FDS1 + 
     &           vh1uimnf(IS+1,IT,IGPD2)*FDS)*FDT1 +
     &        (vh1uimnf(IS,IT+1,IGPD2)*FDS1 +
     &           vh1uimnf(IS+1,IT+1,IGPD2)*FDS)*FDT
      h1dim = (vh1dimnf(IS,IT,IGPD2)*FDS1 + 
     &           vh1dimnf(IS+1,IT,IGPD2)*FDS)*FDT1 +
     &        (vh1dimnf(IS,IT+1,IGPD2)*FDS1 +
     &           vh1dimnf(IS+1,IT+1,IGPD2)*FDS)*FDT
        endif
*
      h1ture = vh1ture(IS,IGPD1)*FDS1 + vh1ture(IS+1,IGPD1)*FDS
      h1tdre = vh1tdre(IS,IGPD1)*FDS1 + vh1tdre(IS+1,IGPD1)*FDS
*
      e1ure =   ve1ure(IS,IGPD1)*FDS1 + ve1ure(IS+1,IGPD1)*FDS
      e1dre =   ve1dre(IS,IGPD1)*FDS1 + ve1dre(IS+1,IGPD1)*FDS
*
      e1ture = ve1ture(IS,IGPD1)*FDS1 + ve1ture(IS+1,IGPD1)*FDS
      e1tdre = ve1tdre(IS,IGPD1)*FDS1 + ve1tdre(IS+1,IGPD1)*FDS
*
      h1tuim = vh1tuim(IS,IGPD1)*FDS1 + vh1tuim(IS+1,IGPD1)*FDS
      h1tdim = vh1tdim(IS,IGPD1)*FDS1 + vh1tdim(IS+1,IGPD1)*FDS
*
      e1uim =   ve1uim(IS,IGPD1)*FDS1 + ve1uim(IS+1,IGPD1)*FDS
      e1dim =   ve1dim(IS,IGPD1)*FDS1 + ve1dim(IS+1,IGPD1)*FDS
*
      e1tuim = ve1tuim(IS,IGPD1)*FDS1 + ve1tuim(IS+1,IGPD1)*FDS
      e1tdim = ve1tdim(IS,IGPD1)*FDS1 + ve1tdim(IS+1,IGPD1)*FDS
*
      gA  = 1.267D0/(1D0 - del2/0.84D0)**2
      gA0 = 0.6D0*gA
      gAu = 0.5D0*( gA + gA0)/( 0.8D0*1.267D0)
      gAd = 0.5D0*(-gA + gA0)/(-0.2D0*1.267D0)
      ha = 4D0*Mp**2*1.267D0/(mpi**2-del2)
*
        if( Ipn .eq. 1 ) then
          if( IGPD .le. 5 ) then
      h1uimag=(F1u/2D0)*h1uim
      h1dimag=(F1d/2D0)*h1dim
      H1_RE = (4D0*(F1u/2D0)*h1ure + F1d*h1dre)/9D0
      H1_IM = (4D0*(F1u/2D0)*h1uim + F1d*h1dim)/9D0
          else
      h1uimag=h1uim
      h1dimag=h1dim
      H1_RE = (4D0*h1ure + h1dre)/9D0
      H1_IM = (4D0*h1uim + h1dim)/9D0
          endif
      h1tuimag=gAu*h1tuim
      h1tdimag=gAd*h1tdim
      H1T_RE = (4D0*gAu*h1ture + gAd*h1tdre)/9D0
      H1T_IM = (4D0*gAu*h1tuim + gAd*h1tdim)/9D0
      E1_RE = (4D0*(F2u/2D0)*e1ure + F2d*e1dre)/9D0
      E1_IM = (4D0*(F2u/2D0)*e1uim + F2d*e1dim)/9D0
      E1T_RE = ha*(4D0*e1ture + e1tdre)/9D0
      E1T_IM = ha*(4D0*e1tuim + e1tdim)/9D0
        elseif( Ipn .eq. 2 ) then
      H1_RE = (4D0*(F1d/2D0)*h1dre + F1u*h1ure)/9D0
      H1_IM = (4D0*(F1d/2D0)*h1dim + F1u*h1uim)/9D0
      H1T_RE = (4D0*gAd*h1tdre + gAu*h1ture)/9D0
      H1T_IM = (4D0*gAd*h1tdim + gAu*h1tuim)/9D0
      E1_RE = (4D0*(F2d/2D0)*e1dre + F2u*e1ure)/9D0
      E1_IM = (4D0*(F2d/2D0)*e1dim + F2u*e1uim)/9D0
      E1T_RE = ha*(4D0*e1tdre + e1ture)/9D0
      E1T_IM = ha*(4D0*e1tdim + e1tuim)/9D0
        endif
*
      return
      end

      subroutine rtable
*
*      A - 1, B - 2, C - 3, D - 4, E - 5 
*
      common/retbl/ vh1ure(51,5),  vh1dre(51,5),
     &              vh1ture(51,5), vh1tdre(51,5),
     &              ve1ure(51,5),  ve1dre(51,5),
     &              ve1ture(51,5), ve1tdre(51,5)
      common/imtbl/ vh1uim(51,5),  vh1dim(51,5),
     &              vh1tuim(51,5), vh1tdim(51,5),
     &              ve1uim(51,5),  ve1dim(51,5),
     &              ve1tuim(51,5), ve1tdim(51,5)
*
*         6, 7
*
      common/retbl1/ vh1urenf(51,21,2), vh1drenf(51,21,2) 
      common/imtbl1/ vh1uimnf(51,21,2), vh1dimnf(51,21,2)
* 
      open(unit=11,file='gpd.dat',status='old')
      do i = 1,51
        read(11,101) vh1ure(i,1),vh1dre(i,1),vh1ture(i,1),vh1tdre(i,1),
     &               ve1ture(i,1),ve1tdre(i,1) 
      enddo
      call ucopy(vh1ure(1,1),ve1ure(1,1),51)
      call ucopy(vh1dre(1,1),ve1dre(1,1),51)
      do j=2,5
        call ucopy(ve1ture(1,1),ve1ture(1,j),51)
        call ucopy(ve1tdre(1,1),ve1tdre(1,j),51)       
      enddo
*
      do j = 2,3
       do i = 1,51
        read(11,102) vh1ure(i,j),vh1dre(i,j),vh1ture(i,j),vh1tdre(i,j)
       enddo
      enddo
      call ucopy(vh1ure(1,2),ve1ure(1,2),102)
      call ucopy(vh1dre(1,2),ve1dre(1,2),102)
*
      do j = 4,5
       do i = 1,51
        read(11,101) vh1ure(i,j),vh1dre(i,j),vh1ture(i,j),vh1tdre(i,j),
     &               ve1ure(i,j),ve1dre(i,j)
       enddo
      enddo
*
      do i = 1,51
        read(11,101) vh1uim(i,1),vh1dim(i,1),vh1tuim(i,1),vh1tdim(i,1),
     &               ve1tuim(i,1),ve1tdim(i,1) 
      enddo
      call ucopy(vh1uim(1,1),ve1uim(1,1),51)
      call ucopy(vh1dim(1,1),ve1dim(1,1),51)
      do j=2,5
        call ucopy(ve1tuim(1,1),ve1tuim(1,j),51)
        call ucopy(ve1tdim(1,1),ve1tdim(1,j),51)       
      enddo
*
      do j = 2,3
       do i = 1,51
        read(11,102) vh1uim(i,j),vh1dim(i,j),vh1tuim(i,j),vh1tdim(i,j)
       enddo
      enddo
      call ucopy(vh1uim(1,2),ve1uim(1,2),102)
      call ucopy(vh1dim(1,2),ve1dim(1,2),102)
*
      do j = 4,5
       do i = 1,51
        read(11,101) vh1uim(i,j),vh1dim(i,j),vh1tuim(i,j),vh1tdim(i,j),
     &               ve1uim(i,j),ve1dim(i,j)
       enddo
      enddo
*
      do IG = 1,2
        do IT = 1,21
          do IS = 1,50
            read(11,103) vh1urenf(IS,IT,IG), vh1drenf(IS,IT,IG)
          enddo
        enddo
      enddo
*
      do IG = 1,2
        do IT = 1,21
          do IS = 1,50
            read(11,103) vh1uimnf(IS,IT,IG), vh1dimnf(IS,IT,IG)
          enddo
        enddo
      enddo
      close( 11 )
*
 101  format(6e12.4)
 102  format(4e12.4)
 103  format(2e12.4)
      end

	   SUBROUTINE UCOPY (A,B,N)
C
C CERN PROGLIB# V301    UCOPY           .VERSION KERNAPO  1.24  920511
C ORIG. 01/03/85  R.BRUN
C
      Dimension      A(*), B(*)

	  do i=1,n
       b(i)=a(i)
	  enddo 
      END

	  
	  
      subroutine nuclFF( del2 )
C
C  Elastic nucleon's formfactors
C
      implicit double precision (A-H,O-Z)
#include "dvcs.inc"
      common/formfac/ F1pn(2), F2pn(2)
      double precision Mv, kp, kn
      parameter (Mv = 0.843D0, kp = 1.79285D0, kn = -1.91D0)
C
      dipol = 1D0/(1D0 - del2/Mv**2)**2
*
      GE_p = dipol
      GE_n = 0D0
      GM_p = (1D0 + kp)*dipol
      GM_n =        kn*dipol
*
      delm = del2/(2D0*Mp)**2
*
      F1pn(1) = (GE_p - delm*GM_p)/(1D0-delm)
      F1pn(2) = (GE_n - delm*GM_n)/(1D0-delm)
      F2pn(1) = (GM_p - GE_p)/(1D0-delm)         
      F2pn(2) = (GM_n - GE_n)/(1D0-delm)         
      return
      end

      subroutine V3subd( A, B, C)
      implicit double precision (A-H,O-Z)
      dimension A(3), B(3), C(3)
      do i = 1,3
        C(i) = A(i) - B(i)
      enddo
      end

      double precision function V3dotd( A, B)
      implicit double precision (A-H,O-Z)
      dimension A(3), B(3)
      S = 0D0
      do i = 1,3
        S = S + A(i)*B(i)
      enddo
      V3dotd = S
      end
c
      subroutine resetmom()
c
C     set to 0 all moments
      implicit none
        include "dvcsmom.inc"
c
c  BH
c
       hc0BH=0D0
       hc1BH=0D0
       hc2BH=0D0
       hc0BHlp=0D0
       hc1BHlp=0D0
       hc0BHtpcos=0D0
       hc1BHtpcos=0D0
       hs1BHtpsin=0D0
c
c  DVCS
c      
       hc0dvcs=0D0
       hc1dvcs=0D0
       hs1dvcs=0D0
       hs2dvcs=0D0
c
c  INTERF
c       
       hs1Iunp=0D0
       hc1Iunp=0D0
       hc0Iunp=0D0
       hs2Iunp=0D0
       hc2Iunp=0D0
c       
       hs1Ilp=0D0
       hc1Ilp=0D0
       hc0Ilp=0D0
c
       hs1Itpcos=0D0
       hs1Itpsin=0D0
       hc1Itpcos=0D0
       hc1Itpsin=0D0
       hc0Itpcos=0D0
       hc0Itpsin=0D0
c            
         return
         end
c
!      subroutine getphoton(Ed,xbd,Q2d,del2d,phield,phigd)
      subroutine getphoton(Ed,xbd,Q2d,del2d,phigd)
c
C     set to 0 all moments
      implicit none
*
      double precision Mp, mele, pi
      common/myconst/ Mp, mele, pi
      dimension V3k1(3), V3k2(3), V3gam(3), V3q(3)
      dimension  V3p1(3), V3p2(3), V3plus(3)
      double precision V3k1,V3k2,V3gam,V3q,V3p1, V3p2,V3plus
!      double precision Ed,xbd,Q2d,del2d,phield,phigd
      double precision Ed,xbd,Q2d,del2d,phigd
      double precision yb,nu,Esc,sintel,costel,costVq,sintVq,qmod
      double precision Ep,Egam,costgg,sintgg,Vgx,Vgy,Vgz
      double precision costeg,teteg,coste1g,tete1g,V3dotd
c
      nu  = Q2d/(2D0*Mp*xbd)
      qmod = sqrt(nu**2 + Q2d)
      yb = nu/Ed
      Esc = Ed - nu
      costel = 1D0 - Q2d/(2D0*Ed*Esc)
      sintel = sqrt(1D0 - costel**2)
c
      Ep   = Mp - del2d/(2D0*Mp)
      Egam = nu + del2d/(2D0*Mp)
c
      V3k1(1) = 0D0
      V3k1(2) = 0D0
      V3k1(3) = Ed
      V3k2(1) = Esc*sintel
      V3k2(2) = 0D0
      V3k2(3) = Esc*costel
*
      call V3subd( V3k1, V3k2, V3q)
      costVq = V3q(3)/qmod
      sintVq = sqrt(1. - costVq**2)
*
      costgg = (2D0*Egam*(Mp + nu) + Q2d - 2D0*Mp*nu)/(2D0*Egam*qmod)
      sintgg = sqrt(1D0 - costgg**2)
      Vgx = Egam*sintgg*cos(phigd)
      Vgy = Egam*sintgg*sin(phigd)
      Vgz = Egam*costgg
*
      V3gam(1) =  Vgx*costVq - Vgz*sintVq
      V3gam(2) =  Vgy
      V3gam(3) =  Vgx*sintVq + Vgz*costVq
*
      call V3subd( V3q, V3gam, V3p2)
      costeg = V3dotd( V3k1, V3gam)/(Ed*Egam)
      teteg  = acos( costeg )
      coste1g = V3dotd( V3k2, V3gam)/(Esc*Egam)
      tete1g  = acos( coste1g )
      return
      end
      
      
      
      
      
      
c$nodebug
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C *						      *
      FUNCTION URAND(IY)
C *						      *
C *   This is a standard pseudo-random generator      *
C *   that work on IBM-370 and IBM-PC. We don't       *
C *   know does it work on SUN? 		      *
C *						      *
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
      REAL*4  URAND,S
      INTEGER*4 IY
      INTEGER*4 A,C,MASK
      PARAMETER (A  = 843314861)
      PARAMETER (C  = 453816693)
      PARAMETER (S  = 4.6566128E-10)

      IY=IAND(A*IY+C,Z'7FFFFFFF')
      URAND=FLOAT(IY)*S
      END

C********************************************************************* 
c
       subroutine vsumm(a,b,c,n)
       real a(n),b(n),c(n)
       integer i,n
        do i=1,n
          c(i)=a(i)+b(i)
        enddo
       return
       end
c
       subroutine vdifm(a,b,c,n)
       real a(n),b(n),c(n)
       integer i,n
        do i=1,n
          c(i)=a(i)-b(i)
        enddo
       return
       end
c
c
       real function vdotm(a,b,n)
       real a(n),b(n),s
       integer i,n
       s=0.0
       do i=1,3
         s=s+a(i)*b(i)
       enddo
       if(n.eq.4) s=s-a(n)*b(n)
       vdotm=s
       return
       end
c   
       real function vangle(a,b,c,d)
       real a(3),b(3),c(3),d(3),xm,ym,vcos
       real x(3),y(3),pi
       pi=acos(-1.0)
       call crossm(a,b,x)
       call crossm(c,d,y)
       xm=vdotm(x,x,3)
       ym=vdotm(y,y,3)
       if(xm.gt.0.0 .and. ym.gt.0.0) then
         vcos=vdotm(x,y,3)/sqrt(xm)/sqrt(ym)
         if(abs(vcos).lt.1.0) then
            vangle=acos(vcos)
         else
            if(vcos.ge.1.0)  vangle=0
            if(vcos.le.-1.0)  vangle=pi
         endif 
       else
         vangle=0
       endif
       return
       end
c
       subroutine crossm(a,b,c)
       real a(3),b(3),c(3)
       c(1)=a(2)*b(3)-a(3)*b(2)
       c(2)=a(3)*b(1)-a(1)*b(3)
       c(3)=a(1)*b(2)-a(2)*b(1)
       return
       end
c
       real function corraul(delt,dQ2)
       double precision delt,dQ2
       real t,Q2,qfac,tfac
       t=-delt
       Q2=dQ2
       if(Q2.lt.1.5) Q2=1.5
       if(Q2.gt.2.6) Q2=2.6
       qfac=(0.69*t-1.04*t*t)*(q2-1.5)/1.1+1.0
       tfac=(0.908-alog(t)/9.5)
       corraul=1.0/qfac*tfac
       return 
       end
c
