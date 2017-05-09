! en este programa se establece el potencial químico cos(wt) en la frontera
! w=0.01, movilidades en 10
         program modelob

        implicit none

         integer nx,ny,distinic,vecinos,veci,nxe,nye,np,cl,cc
!cl es de ciclos largos: numero de periodos de oscilacion
!cc es el ciclo de una oscilacion
         double precision mui,dt,dx,e,phiii,mov1,mov2,tri,w
!mui es el valor de mu impuesto en la frontera
!phiii es el valor equivalete impuesto en phi (resolver la cubica)
         integer conosin,nr,dm
!conosin debe ser 1 si las condiciones en mu serán impuestas por mui, de lo contrario deberá ser cero

         parameter (distinic=20,mui=0.08,conosin=1,mov1=1.,mov2=1.)
         parameter (nx=50,ny=30,vecinos=2,np=nx*ny,veci=2*vecinos)
         parameter (dx=1.0,e=1.,tri=0.1)
	 

	     parameter (w=0.001)
         		 
         integer muestreo,nump

!cl de ciclos largos es el número de iteraciones que hará el programa el ciclo interno. Es el número de archivos de salida que se quiere
!cc de ciclo corto es el número de iteraciones que se hace de la subrutina modelo b antes de imprimr la salida en un fort.nr

         parameter (nye=ny+veci,nxe=nx+veci)

         integer ix,iy
         double precision phi(nye,nxe),mu(nye,nxe),lap(nye,nxe)
	     double precision phi1(nye,nxe)
         double precision M(nye,nxe),ayi(nye,nxe),coordx(np),coordy(np)
         double precision axd(nye,nxe),ayd(nye,nxe),axi(nye,nxe),tiempo
         integer dx2,i,ncoord,mood
         character oi

       dt=0.1
	   cc=2.*3.1416/(w*dt)
	   cl=100
	   nump=100
	
	   muestreo=100
 
        call condinic(phi,nye,nxe,distinic,phi1)

	  call coordenadas(phi,coordx,coordy,ncoord,nye,nxe,np)
	  dm=coordy(15)

        open(unit=10,file="posicioninterfaz.txt",status="unknown")
        write(*,*)"mui=",mui,"  distancia inicial=",distinic
        write(*,*)"tama�o del sistema ",nx," X",ny
       	write(*,*)"el modulo de impresión es",muestreo 
        write(*,*)"el n�mero de puntos por periodo es",nump
	write(*,*)"la posicion inicial de la interfaz medida es",dm," se corre por:"
&,cl,"ciclos largos de",cc, " ciclos cortos cada uno"
	write(*,*)""
	write(*,*)"frecuencia",w



        do ix=1,nxe
        do iy=1,nye
        M(iy,ix)=1
        enddo
        enddo


        do nr=12,cl+12



        do i=1,cc

       tiempo=((nr-12)*cc+i)*dt

        call coordenadas(phi,coordx,coordy,ncoord,nye,nxe,np)
!        call movilidad(M,nxe,nye,mov1,mov2,coordx,coordy,np)
        call laplaciano(phi,lap,ny,nx,vecinos,dx,nxe,nye)
        call modelo(mu,phi,nxe,nye,e,dx,dt,lap,M,mui,nx,ny,vecinos
&,conosin,phiii,phi1,tri,tiempo,w,coordy,dm,np)
        call condfrontera(phi,nxe,nye,vecinos,veci,mu)
       mood=MOD(i,muestreo)
       if(nr.gt.12+70) then
       if(mood.eq.0) then
        call coordprint(w,coordx,coordy,np,tiempo,mu,phi,nxe,nye)
       endif
       endif

       enddo

       call coordenadas(phi,coordx,coordy,ncoord,nye,nxe,np)

       enddo

       call guardar(phi1,phi,mu,M,nye,nxe,tiempo,mov1,mov2
&,w,tri)
       end


       subroutine movilidad(M,nxe,nye,mov1,mov2,coordx,coordy,np)
       integer ix,iy,np,nxe,nye,ixi
       double precision coordx(np),coordy(np),M(nye,nxe),mov1,mov2
       double precision ixip

!!        open(unit=99,file="movilidad",status="unknown")

        do iy=1,nye

        ixi=nint(coordy(iy+1))
!!        write(99,*)ixi,coordy(iy+1)
        do ix=ixi+2,nxe
        M(iy,ix)=mov2
        enddo
        
        do ix=1,ixi+1
        M(iy,ix)=mov1
        enddo

        enddo

!!        close(99)
        end


        subroutine modelo(mu,phi,nxe,nye,e,dx,dt,lap,M,mui,nx,
&ny,vecinos,conosin,phiii,phi1,tri,tiempo,w,coordy,dm,np)

        integer nx,ny,ix,iy,nxe,nye,vecinos,conosin,dm,np
        double precision mu(nye,nxe),phi(nye,nxe),axd(nye,nxe)
        double precision axi(nye,nxe),M(nye,nxe),w,coordy(np)
        double precision ayd(nye,nxe),ayi(nye,nxe),lap(nye,nxe)
        double precision e,dx,dt,mui,phiii,dt2,zz,yy,tiempo
	double precision LD,phi1(nye,nxe),tri,phiaux(nye,nxe)

! la idea es imponer mu en las primeras dos filas en x
!conosin vale 1 para condiciones en mu dadas
!        conosin=1

         dt2=dt*dt

        do ix=vecinos+1,nxe-1
        do iy=vecinos,nye-1
        mu(iy,ix)=((phi(iy,ix))**3)-phi(iy,ix)-(e**2)*lap(iy,ix)
        enddo
        enddo

        call phiauxiliar(phi,phiaux,nxe,nye)
        call condmu(mu,phi,nxe,nye,mui,vecinos,phiii,tiempo,w
&,coordy,dm,np)

! phi se va a reescribir en el interior de la celda, por lo que esto será desde x=2, pues el cálculo requiere una celda extra para el lapaciano

        dx2=dx*dx
        do ix=vecinos+1,nxe-1
        do iy=vecinos,nye-1

        axd(iy,ix)=M(iy,ix+1)*(mu(iy,ix+1)-mu(iy,ix))
&-M(iy,ix)*(mu(iy,ix)-mu(iy,ix-1))
        ayd(iy,ix)=M(iy+1,ix)*(mu(iy+1,ix)-mu(iy,ix))
&-M(iy,ix)*(mu(iy,ix)-mu(iy-1,ix))
        axi(iy,ix)=M(iy,ix)*(mu(iy,ix+1)-mu(iy,ix))
&-M(iy,ix-1)*(mu(iy,ix)-mu(iy,ix-1))
        ayi(iy,ix)=M(iy,ix)*(mu(iy+1,ix)-mu(iy,ix))
&-M(iy-1,ix)*(mu(iy,ix)-mu(iy-1,ix))
        LD=0.5*(axd(iy,ix)+ayd(iy,ix)+axi(iy,ix)+ayi(iy,ix))/dx2

        yy=(2*tri-dt)
        zz=1/(2*tri+dt)



        phi(iy,ix)=LD*2*dt2*zz+4*tri*phi(iy,ix)*zz
        phi(iy,ix)=phi(iy,ix)-phi1(iy,ix)*yy*zz
		
        enddo
        enddo

        call recuperarphi(phiaux,phi1,nye,nxe)

        end

********************************************

        subroutine phiauxiliar(phi,phiaux,nxe,nye)

         integer ix,iy,nxe,nye
         double precision phi(nye,nxe),phiaux(nye,nxe)
         do ix=1,nxe
         do iy=1,nye

         phiaux(iy,ix)=phi(iy,ix)

         enddo
         enddo
         end

********************************************

        subroutine recuperarphi(phiaux,phi1,nye,nxe)
 
        integer ix,iy
        double precision phi1(nye,nxe),phiaux(nye,nxe)

	    do ix=1,nxe
        do iy=1,nye

        phi1(iy,ix)=phiaux(iy,ix)

        enddo
        enddo

        end

********************************************

        subroutine condfrontera(phi,nxe,nye,vecinos,veci,mu)

        integer vecinos,veci,nxe,nye,ix,iy
        double precision phi(nye,nxe),mu(nye,nxe)
!       condiciones de espejo
        do ix=1,vecinos
        do iy=1,nye
!       phi(iy,ix)=phi(iy,veci+1-ix)
        phi(iy,nxe-ix+1)=mu(iy,nxe-veci+ix)
        mu(iy,nxe-ix+1)=mu(iy,nxe-veci+ix)
        enddo
        enddo

        do ix=1,nxe
        do iy=1,vecinos
        mu(iy,ix)=mu(nye-veci+iy,ix)
        mu(nye-vecinos+iy,ix)=mu(iy+vecinos,ix)

        phi(iy,ix)=phi(nye-veci+iy,ix)
        phi(nye-vecinos+iy,ix)=phi(iy+vecinos,ix)
        enddo
        enddo

        end

***************************************************

        subroutine export(phi,mu,nye,nxe,nr,tiempo,phi1)

        integer ey,ex,ny,nxe
        double precision phi(nye,nxe),mu(nye,nxe),tiempo
	double precision phi1(nye,nxe)

        open(unit=3,file="recover",status="unknown")
	write(3,*) tiempo
        do ex=1,nxe
        do ey=1,nye
	write(3,*) phi(ey,ex),phi1(ey,ex)
        enddo
        enddo
        close(3)
        end

******************************************************

        subroutine condinic(phi,nye,nxe,distinic,phi1)

        integer nye,nxe,distinic
        double precision phi(nye,nxe),phi1(nye,nxe),a
        integer ix,iy

        a=sqrt(2.0)

        do ix=1,nxe
        do iy=1,nye

        phi(iy,ix)=tanh((ix-101.5)/a)
        phi1(iy,ix)=tanh((ix-101.5)/a)

        enddo
        enddo

        end

*********************************************************

        subroutine recover(phi,nye,nxe,tiempo)

        integer nye,nxe,ix,iy,tiempo
        double precision phi(nye,nxe),w

        open(98,file="recover",status="old")

	    read(98,*) tiempo
        do ix=1,nxe
        do iy=1,nye
        read(98,*)w
        phi(iy,ix)=w
        enddo
        enddo

        close(98)
        return

        end

****************************************************

        subroutine laplaciano(phi,lap,ny,nx,vecinos,dx,nxe,nye)

        integer ny,nx,nxe,nye
        double precision phi(nye,nxe),lap(nye,nxe),dx
        integer ix,iy,xi1,xi2,xd1,xd2,yi1,yi2,yd1,yd2,vecinos

        dx2=1/(dx*dx)
        do ix=vecinos,nxe-1
        do iy=vecinos,nye-1
        xi1=ix-1
        xi2=ix-2
        yi1=iy-1
        yi2=iy-2
        xd1=ix+1
        xd2=ix+2
        yd1=iy+1
        yd2=iy+2

        lap(iy,ix)=dx2*((phi(iy,xd1)+phi(iy,xi1))+(phi(yi1,ix)+
&phi(yd1,ix)-4.*phi(iy,ix)))

        enddo
        enddo
        continue

        end

***********************************************

        subroutine condmu(mu,phi,nxe,nye,mui,vecinos,phiii,tiempo,w
&,coordy,dm,np)

        integer ix,iy,nxe,nye,vecinos,np,dm
        double precision mu(nye,nxe),phi(nye,nxe),w,tiempo
        double precision mui,fron,phiii,muii,coordy(np)

       muii=mui*sin(w*tiempo)
        call phii(muii,fron)	

        do iy=1,nye
	    phi(iy,1)=-1
	    phi(iy,2)=fron
        mu(iy,1)=muii.
	
        enddo

        end

*************************************************
!PHII

        subroutine phii(muii,fron)
        double precision fron,muii,b,c,PI,a,phiii
        a=1
        b = a/3
        c = -muii/2
        fron=-2*sqrt(b)*cos((acos(c/(sqrt((b**3)))))/3)
        end

        subroutine perfil(phi,mu,nye,nxe,nr)

        integer iy,ix,nye,nxe,nr
        double precision phi(nye,nxe),mu(nye,nxe)
        
        iy=10
        
        do ix=1,nxe
        write(nr,*)ix,phi(iy,ix),mu(iy,ix)
!        write(nr,*)ix,mu(iy,ix)
        enddo
 
	close(nr)
	end

        subroutine perfilmu(phi,mu,nye,nxe,nr)

        integer iy,ix,nye,nxe,nr
        double precision phi(nye,nxe),mu(nye,nxe)
        
        iy=10
        
        do ix=1,nxe
!        write(nr,*)ix,phi(iy,ix)
        write(nr,*)ix,mu(iy,ix)
        enddo
        
       
        close(nr)

        end


		subroutine guardar(phi1,phi,mu,M,nye,nxe,tiempo,mov1,mov2
&,w,tri)
		
		
		integer nye,nxe,ix,iy
		double precision phi1(nye,nxe),mu(nye,nxe),phi(nye,nxe)
		double precision tiempo,M(nye,nxe),mov1,mov2,w,tri
        open(unit=1,file="recover",status="unknown")
		write(1,*)tiempo,mov1,mov2,w,tri
		do ix=1,nxe
		do iy=1,nye
		write(1,*)phi1(iy,ix),phi(iy,ix),mu(iy,ix),M(iy,ix)
		enddo
		enddo
		


		end

        subroutine coordprint(w,coordx,coordy,np,tiempo,mu,phi
&,nxe,nye)
      
        integer np,wk,nk,nr,nye,nxe
        double precision coordx(np),coordy(np),tiempo,w
        double precision mu(nye,nxe),phi(nye,nxe)

        wk=15
        
        write(10,*)tiempo,coordy(wk),mu(wk,nint(coordy(wk)))-
&mu(wk,nint(coordy(wk))-1),phi(wk,coordy(wk))

        end

      subroutine coordenadas(p,x,y,npun,lx,ly,ninter)
      implicit double precision (a-h,o-z)
      integer lx,ly,ninter
      DIMENSION p(lx,ly)
      DIMENSION X(ninter),Y(ninter)
      DIMENSION KN(14,0:4),IN(4,lx),JN(4,ly)
      DIMENSION NRX(lx),NLX(lx)
      DIMENSION NRY(ly),NLY(ly)
      DIMENSION IP(lx,ly),IS(lx,ly)

      DATA KN/2,1,1,1,1,2,3,3,2,1,1,1,1,2,                              FRA00110
     &        2,0,1,0,2,0,4,4,0,4,0,1,0,2,                              FRA00120
     &        0,0,0,1,1,2,3,3,2,3,1,0,0,0,                              FRA00130
     &        0,2,3,4,4,0,0,0,0,2,4,3,2,0,                              FRA00140
     &        3,1,0,0,3,4,0,0,4,1,0,0,1,3/                              FRA00150


       imin=1
       imax=ly
       npun=0
                                                                        FRA00170
       NLx(1)=lx                                                        FRA00180
       NRx(1)=2                                                         FRA00190
       DO I=2,lx-1
         NLx(I)=I-1                                                     FRA00210
         NRx(I)=I+1                                                     FRA00220
       enddo
       NLx(lx)=lx-1                                                     FRA00230
       NRx(lx)=1

       NLy(1)=ly                                                        FRA00180
       NRy(1)=2                                                         FRA00190
       DO I=2,ly-1
         NLy(I)=I-1                                                     FRA00210
         NRy(I)=I+1                                                     FRA00220
       enddo
       NLy(ly)=ly-1                                                     FRA00230
       NRy(ly)=1                                                        FRA00120

       DO I=1,Lx                                                        FRA00250
        IN(1,I)=NRx(I)                                                  FRA00260
        IN(2,I)=I                                                       FRA00270
        IN(3,I)=NLx(I)                                                  FRA00280
        IN(4,I)=I                                                       FRA00290
      enddo
       DO I=1,Ly                                                        FRA00250
        JN(1,I)=I                                                       FRA00300
        JN(2,I)=NLy(I)                                                  FRA00310
        JN(3,I)=I
        JN(4,I)=NRy(I)                                                  FRA00330
      enddo

      N=Lx*Ly                                                           FRA00360
      RMESH=1.0
                                                                        FRA00370
      DO 11 J=1,LY
      DO 11 I=1,LX
         IP(I,J)=0
         IF(P(I,J).GT.0.)IP(I,J)=1
 11   CONTINUE                                                          FRA00510
      DO 1 J=1,LY                                                       FRA00520
      DO 1 I=1,LX                                                       FRA00530
        IJ=IP(I,J)+2*IP(NRx(I),J)+4*IP(NRx(I),NRy(J))
     & +8*IP(I,NRy(J))
        IF(IJ.EQ.15)IJ=0                                                FRA00550
   1    IS(I,J)=IJ

      DO 2 LN=1,NINTER                                                  FRA00580

c         DO 4 J0=1,LY                                                  FRA00590
          DO 4 J0=imin,imax
          DO 4 I0=1,LX
             IF(IS(I0,J0).NE.0)GOTO 5                                   FRA00610
  4       CONTINUE                                                      FRA00620
          GOTO 99                                                       FRA00630
  5       CONTINUE                                                      FRA00640
          CALL SMOOTH(Lx,ly,P,NRX,NRY,I0,J0,RI,RJ)
          X(1)=RI                                                       FRA00660
          Y(1)=RJ                                                       FRA00670
          K=KN(IS(I0,J0),0)                                             FRA00680
           I=IN(K,I0)
           J=JN(K,J0)
      IF(IS(I0,J0).NE.5.AND.IS(I0,J0).NE.10)IS(I0,J0)=0                 FRA00710
         IF(IS(I0,J0).EQ.5.AND.(K.EQ.1.OR.K.EQ.4))IS(I0,J0)=1           FRA00720
         IF(IS(I0,J0).EQ.5.AND.(K.EQ.2.OR.K.EQ.3))IS(I0,J0)=4           FRA00730
         IF(IS(I0,J0).EQ.10.AND.(K.EQ.3.OR.K.EQ.4))IS(I0,J0)=2          FRA00740
         IF(IS(I0,J0).EQ.10.AND.(K.EQ.1.OR.K.EQ.2))IS(I0,J0)=8          FRA00750
        CALL SMOOTH(LX,LY,P,NRX,NRY,I,J,RI,RJ)
           X(2)=RI                                                      FRA00770
           Y(2)=RJ                                                      FRA00780
       DO 3 IJ=3,lx*ly
c      DO 3 IJ=imin+2,imax*imax
           K=KN(IS(I,J),K)                                              FRA00800
         IF(IS(I,J).NE.5.AND.IS(I,J).NE.10)IS(I,J)=0                    FRA00810
         IF(IS(I,J).EQ.5.AND.(K.EQ.1.OR.K.EQ.4))IS(I,J)=1               FRA00820
         IF(IS(I,J).EQ.5.AND.(K.EQ.2.OR.K.EQ.3))IS(I,J)=4               FRA00830
         IF(IS(I,J).EQ.10.AND.(K.EQ.3.OR.K.EQ.4))IS(I,J)=2              FRA00840
         IF(IS(I,J).EQ.10.AND.(K.EQ.1.OR.K.EQ.2))IS(I,J)=8              FRA00850
           I=IN(K,I)                                                    FRA00860
           J=JN(K,J)                                                    FRA00870
         CALL SMOOTH(LX,LY,P,NRX,NRY,I,J,RI,RJ)
            X(IJ)=RI                                                    FRA00890
           Y(IJ)=RJ                                                     FRA00900
           IF((I.EQ.I0).AND.(J.EQ.J0))GOTO 6                            FRA00910
  3     CONTINUE                                                        FRA00120
        GOTO 99                                                         FRA00930
  6     CONTINUE                                                        FRA00940

                                                                        FRA00950
        IF(IJ.GE.32)THEN

        npun=npun+ij
                                                                        FRA00970
        goto 99

        ENDIF                                                           FRA01330
                                                                        FRA01340
  2   CONTINUE
 99   CONTINUE                                                          FRA01380
      END

      SUBROUTINE SMOOTH(LX,LY,P,NRX,NRY,I,J,RI,RJ)
                                                                        FRA01430
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                FRA01440
      DIMENSION p(lx,ly),nrx(lx),nry(ly)
      rlx=dfloat(lx)
      rly=dfloat(ly)
      RI0=I-1                                                           FRA01460
      RJ0=J-1                                                           FRA01470
C CCCCCCCC                                                              FRA01480
      IF(I.EQ.1)RI0=LX                                                  FRA01500
      IF(J.EQ.1)RJ0=LY                                                  FRA01500
C CCCCCCCC       

      P1=P(I,J)                                                         FRA01520
      P2=P(NRx(I),J)                                                    FRA01530
      P3=P(NRx(I),NRy(J))                                               FRA01540
      P4=P(I,NRy(J))                                                    FRA01550
      P0=0.25*(P1+P2+P3+P4)                                             FRA01560
      PX=0.5*(P3-P4+P2-P1)                                              FRA01570
      PY=0.5*(P4-P1+P3-P2)                                              FRA01580
      P2=PX*PX+PY*PY                                                    FRA01590
      RI=RI0-P0*PX/P2                                                   FRA01600
      RJ=RJ0-P0*PY/P2                                                   FRA01610
      IF(RI.GT.RLx)RI=RI-RLx                                            FRA01620
      IF(RI.LT.0.)RI=RI+RLx                                             FRA01630
      IF(RJ.GT.RLy)RJ=RJ-RLy                                            FRA01640
      IF(RJ.LT.0.)RJ=RJ+RLy                                             FRA01650
      RETURN                                                            FRA01660
      END

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               REA00080
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

