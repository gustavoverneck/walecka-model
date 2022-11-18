!==================================================================
      !legenda
      !rmb = vetor com 8 componentes - fermions - r real mb massa barions
      !rml = vetor com 2 componentes - leptons - r real mb massa barions
      !qb = vetor c 8 compont - cargas
      !i, j, k, l, m, n = contadores inteiros
      !35 99130-6304  ==  luiz zap
!==================================================================	  
!para rodar os codigos:	  

! gfortran wal.f90 -w -o executavel
!==================================================================
      
program bomeos
!     calculates the eos for the baryon octet model (nlwm) at t=0
!     uses glendenning couplig constants
!     by chiapparini and modified by menezes
use nr
USE nrtype

implicit none

logical check
REAL(SP), DIMENSION(5):: x,fvec			!x sao 5 variaveis, fvec equacoes
REAL(SP), DIMENSION(8):: rmb,qb,ri3,gs,gv,gr,dr,sb	
REAL(SP), DIMENSION(2):: rnl,rml
REAL(SP):: rlam,csi,rnbt,schargeb,xls,xlw,xlr,xs,xw,xr,rk,rmns,rmue,rmun,rmpu,rmpummev,rmpumfm
REAL(SP):: rnbsup,rnbinf,rm,rms,rmv,rmr,dnb,dener,ener,press,rkfe,rkfn,vsigma,vomega,vrho,rnbtd,rne0,rnn0,rphi,rv0
REAL(SP), parameter:: pi2=pi*pi
integer:: nh, i, j, numberh, npoint,gm


! FIXME Ian: Coloquei esse rnb aqui mas n sei o tamanho que deveria ser
REAL(SP), dimension(8) :: rnb
! END FIXME

!dados fisicos das particulas  poderia ser definido direto em :  real rmb(8)/ xxx,xxx, ... /
rmb = [939.d0,939.d0,1116.d0,1193.d0,1193.d0,1193.d0,1318.d0,1318.d0]  ! definindo vetor de 8 componentes
rml = [0.511d0,105.66d0]
qb = [0.0d0,1.0d0,0.0d0,-1.0d0,0.0d0,1.0d0,-1.0d0,0.0d0]
ri3 = [-0.5d0,0.5d0,0.0d0,-1.0d0,0.0d0,1.0d0,-0.5d0,0.5d0]
sb = [0.0d0,0.0d0,-1.d0,-1.d0,-1.d0,-1.d0,-2.d0,-2.d0]
		
!indices das especies das particulas nucleares
!     1=neutron,2=proton,3=lambda,4=sigma-,5=sigma0,6=sigma+,7=cascade-,8=cascade0
!     1=electron,2=muon
!massas
!mb(1) = 939.d0, mb(2) = 939.d0, mb(3) = 1116.d0, mb(4) = 1193.d0, mb(5) = 1193.d0,  mb(6) = 1193.d0
!mb(7) = 1318.d0,  mb(8) = 1318.d0
!ml(1) = 0.511d0, ml(2) = 105.66d0
		
!cargas
!qb(1)=0.0d0, qb(2)=1.0d0, qb(3)=0.0d0, qb(4)=-1.0d0, qb(5)0.0d0, qb(6)=1.0d0, qb(7)=-1.0d0, qb(8)=0.0d0
!ql(1)=1.0d0, ql(2)=1.0d0,
		
!isospin barions
!i3(1)=-0.5d0, i3(2)=0.5d0, i3(3)=0.0d0, i3(4)=-1.0d0, i3(5)=0.0d0, i3(6)=1.0d0, i3(7)=-0.5d0, i3(8)=0.5d0
		
		
write(6,*)'2 for nucleons and 8 for nucleons + hyperons'
read(5,*) nh 
write(6,*) 'vc está rodando matéria estelar'

!data input

open(newunit=gm,file='gm1.inp')		!lendo os dados do arquivo gm1.inp
read(gm,*) rm				!939.d0	nucleon mass		  proton e neutron
read(gm,*) rms				!512.d0	scalar meson mass	   m_sigma
read(gm,*) rmv				!783.d0	vector meson mass	   m_omega
read(gm,*) rmr				!770.d0	vector-isovector meson mass	   m_rho
read(gm,*) xls				!11.785d0 	gsn^2/ms^2		xls 
read(gm,*) xlw				!7.1480d0  	gwn^2/mw^2		xlw
read(gm,*) xlr				!4.939d0  	grn^2/mr^2		xlr
read(gm,*) rk				!0.005894d0	k/m			rk	termo cubico n linear
read(gm,*) rlam				!-0.00642d0   	lambda			rlam	termo quartico n linear
read(gm,*) csi				!0.0		csi			csi
read(gm,*) xs				!0.7		gs/gn			xs sigma-hiperon coupling constant
read(gm,*) xw				!0.783       	gw/gn 			xw omega-hiperon coupling constant
read(gm,*) xr				!0.783       	gr/gn 			xr rho-hiperon coupling constant
read(gm,*) rnbinf			!0.03d0	nbinf [fm^-3]		rnbinf	range limite inferior
read(gm,*) rnbsup			!1.23d0	nbsup [fm^-3]		rnbsup	range limite superior
read(gm,*) npoint			!400		number of points	npoint

close(gm)

rms=rms/rm 				!admensionalidade massa do meson escalar m_sigma
rmv=rmv/rm  				!admensionalidade massa do meson escalar m_omega     
rmr=rmr/rm 				!admensionalidade massa do meson escalar m_rho
rnbinf=rnbinf/(rm/197.32d0)**3		!admensionalidade limite inferior
rnbsup=rnbsup/(rm/197.32d0)**3		!admensionalidade limite inferior
dnb=(rnbsup-rnbinf)/(npoint-1)		!numero de passos

do i=3,8
	rmb(i)=rmb(i)/rm		!admensionalidade barions/hiperions
enddo

do i=1,2				!admensionalidade barions/nucleons e dos leptons
	rmb(i)=1.d0
	rml(i)=rml(i)/rm
enddo
		
xls=xls*(rm/197.33d0)**2		!admensionalidade da const acoplamento meson scalar sigma
xlw=xlw*(rm/197.33d0)**2		!admensionalidade da const acoplamento meson scalar omega
xlr=xlr*(rm/197.33d0)**2		!admensionalidade da const acoplamento meson scalar rho

do i=1,2
	gs(i)=sqrt(xls)*rms			!const acoplamento
	gv(i)=sqrt(xlw)*rmv
	gr(i)=sqrt(xlr)*rmr
enddo

do i=3,8
	gs(i)=xs*gs(1)    		!sigma-hiperon coupling constant
	gv(i)=xw*gv(1)    		!omega-hiperon coupling constant
	gr(i)=xr*gr(1)    		!rho-hiperon coupling constant
enddo

rk=rk*gs(1)**3				!constante k do termo sigma cubico do lagrangeano
rlam=rlam*gs(1)**4			!constante lambda do termo sigma quartico do lagrangeano

!     initializing variables		!condicoes iniciais das variaveis

rne0=0.1d0*rnbinf				
rnn0=0.9d0*rnbinf      
x(1)=(3.d0*pi2*rne0)**0.1666666667d0
x(2)=(3.d0*pi2*rnn0)**0.1666666667d0
x(3)=dasin(sqrt(0.005d0/0.999999d0))     	!     x(4)=gvn/rmv*sqrt(rnbinf)
x(4)=sqrt(0.22731d-05)			!     x(5)=(grn/rmr)**2*(ri3(1)*rnn0+ri3(2)*rne0)
x(5)=-0.57641d-06

!main loop 	!loop principal que calcula as equações de estado

do i=1,npoint
	rnbt=rnbinf+(i-1)*dnb

            
	call broydn(x,check)
		if(check)pause 'there is no root in broydn...'

	call mapping(x,rkfe,rkfn,vsigma,vomega,vrho)

	call eos(rkfe,rkfn,vsigma,vomega,vrho,ener,press,rv0,rphi)    !the eos
 
	rnbtd=rnbt*(rm/197.32d0)**3
	rmns=(rmb(1)-vsigma)

	rmun=vomega+ri3(1)*vrho + sqrt((3.d0*pi2*rnb(1))**0.6666666667d0+rmns**2)   !neutron chemical potential
	rmue=sqrt((3.d0*pi2*rnl(1))**0.6666666667d0+rml(1)**2) !electron chemical potential
         
	dener=ener/rnbt-1.d0   		!dens de energia de ligação - massa (valor 1)
	dener=dener*rm			!multiplicado por rm para ficar em mev


!     output in fm^-4 for the program that calculates stars properties; 
!     if you choose this output, comment the above one
	press=press*(rm/197.32)**4			!pressao em fm^-4
	ener=ener*(rm/197.32)**4			!energia em fm^-4


	write(22,10) rnbtd, (rnl(j)/rnbt,j=1,2), (rnb(j)/rnbt,j=1,nh)
	
	write(100,*) rnbtd,ener,press,rmun*rm
  
	rmpu = rmb(1)*rnb(1)+rmb(2)*rnb(2) +rmb(3)*rnb(3) + rmb(4)*rnb(4) &
	+ rmb(5)*rnb(5)+rmb(6)*rnb(6)+rmb(7)*rnb(7)  + rmb(8)*rnb(8)
             
	rmpummev = ((rmpu/(rnbtd)))*(939.d0/197.33d0)**3*939.d0             !massa média por partícula
	rmpumfm =  rmpummev/197.33d0
     
enddo

10 	format(11(1x,e12.5)) 
      
stop

end program

!--------------------------------------------------------------------------
subroutine mapping(x,rkfe,rkfn,vsigma,vomega,vrho)
!     this subroutine sets the potential fields between the physical values    

	use nrtype
	implicit none 
	 
	REAL(SP), DIMENSION(5):: x			!x sao 5 variaveis, fvec equacoes
	REAL(SP):: rkfe,rkfn,vsigma,vomega,vrho
   
	rkfe=x(1)**2                            	!to assure kfe greater than zero
	rkfn=x(2)**2                            	!to assure kfn greater than zero
	vsigma=0.999999d0*sin(x(3))**2   		!mapping vsigma in (0,mb(1))
	vomega=x(4)**2                       		!to assure vomega greater than zero
	vrho=x(5)                              	!no mapping on vrho
      
	return
	end   
         
!-----------------------------------------------------------------------
subroutine skfb2(i,rkfe,rkfn,vsigma,vomega,vrho,rkfb2)
!     this subroutine calculates the fermi momentum squared of the ith baryon
!     using the constrains on the chemical potential

	use nrtype
	implicit none  
	REAL(SP), DIMENSION(8):: rmb,rnb,qb,ri3,gs,gv,gr,dr,sb
	REAL(SP):: rkfe,rkfn,vsigma,vomega,vrho,rkfe2,rkfn2,rkfb,rkfb2,rmns,rmbs,efn,rmue,rmun,rmub
	integer:: nh, i, j
	

	REAL(SP), DIMENSION(2):: rnl,rml

	
	rkfe2=rkfe**2					!quadrado de  rkfe=x(1)**2
	rkfn2=rkfn**2					!quadrado de  rkfn=x(2)**2

	rmns=rmb(1)-vsigma            		!neutron effective mass  eq 4.4    rmns=rmb(1)-vsigma*gs(1)/gs(1)
	rmbs=rmb(i)-vsigma*gs(i)/gs(1)          	!baryon effective mass  eq 4.4

	efn=sqrt(rkfn2+rmns**2)                     	!neutron fermi energy 3.7
	rmue=sqrt(rkfe2+rml(1)**2)                  	!electron chemical potential  3.7 mas para o lepton eletron
	rmun=vomega+ri3(1)*vrho+efn                 	!neutron chemical potential 5.7
	rmub=rmun-qb(i)*rmue                         	!baryon chemical potential eq 5.33
	rkfb2=(rmub-gv(i)/gv(1)*vomega-gr(i)/gr(1)*ri3(i)*vrho)**2-rmbs**2  !baryon fermi mom squared
	!kfb^2 = (mu_b - g_omega*omega - g_rho*i3*rho)^2 - mb_sigma*^2   eq 5.7 isolando kfb^2

	return
	end      

!-----------------------------------------------------------------------
subroutine eos(rkfe,rkfn,vsigma,vomega,vrho,ener,press,rphi,rv0)
!	esta subrotina calcula o eos desta matéria estelar

	use nrtype
	implicit none  
	REAL(SP), DIMENSION(5):: x,fvec			!x sao 5 variaveis, fvec equacoes
	REAL(SP), DIMENSION(8):: rmb,qb,ri3,gs,gv,gr,dr,sb	!
	REAL(SP), DIMENSION(2):: rnl,rml
	REAL(SP):: rlam,csi,rnbt,rnb,schargeb,rma,rnu,rmu,rk,kfmu2,rkfmu,rkfmu2
	REAL(SP):: rkfe,rkfn,vsigma,vomega,vrho,rkfb2,rkfb,rphi,rv0,b0
	REAL(SP):: enerf,ener,rmv,rmr,rms,rmbs,efn,rmue,rmun,press,enerbar,pressbar,enerlep,presslep,pressf
	REAL(SP):: re3,re4,re3le,re4le,re3lmu,re4lmu,f3,f4
	integer:: nh, i, j, numberh

	REAL(SP) :: rmns

	

	rphi=vsigma/gs(1)
	rv0=vomega/gv(1)
	b0=-vrho/gr(1)

	enerf=rmv**2*rv0**2/2.d0+csi*gv(1)**4*rv0**4/8.d0 &	!dens de energia no momento de fermi 5.36
	+ rmr**2*b0**2/2.d0 &
	+ rms**2*rphi**2/2.d0 + rk*rphi**3/6.d0 &
	+ rlam*rphi**4/24.d0

	pressf=rmv**2*rv0**2/2.d0+csi*gv(1)**4*rv0**4/24.d0 &  	!pressao no momento de fermi 5.37
	+rmr**2*b0**2/2.d0 & 
	-rms**2*rphi**2/2.d0 - rk*rphi**3/6.d0 &
	-rlam*rphi**4/24.d0

	enerbar=0.d0		!energia barionica inicial
	pressbar=0.d0		!pressao barionica inicial

	do i=1,nh
		call skfb2(i,rkfe,rkfn,vsigma,vomega,vrho,rkfb2)         
		if(rkfb2.gt.0.d0)then                                  
			rkfb=sqrt(rkfb2)     
		else
			rkfb=0.d0
		endif     

		rmbs=rmb(i)-vsigma*gs(i)/gs(1)               
		rma=rmbs
		rmns=rmb(1)-vsigma   
		efn=sqrt(rkfn**2+rmns**2)                    
		rmue=sqrt(rkfe**2+rml(1)**2)                  
		rmun=vomega+ri3(1)*vrho+efn                 
		rmu=rmun-qb(i)*rmue
		rnu=rmu-gv(i)/gv(1)*vomega-gr(i)/gr(1)*ri3(i)*vrho
	
		call gauss(f3,0.d0,rkfb,10,re3)
		enerbar=enerbar+re3
         
		call gauss(f4,0.d0,rkfb,10,re4)
		pressbar=pressbar+re4
	enddo

	enerlep=0.d0
	presslep=0.d0
	rkfmu=0.d0                         
	kfmu2=rkfe**2+rml(1)**2-rml(2)**2                           
              
	if(rkfmu2.gt.0.d0)rkfmu=sqrt(rkfmu2)      

		rmue=sqrt(rkfe**2+rml(1)**2)   
		rnu=rmue
		rma=rml(1)
		
	call gauss(f3,0.d0,rkfe,10,re3le)
	enerlep=enerlep+re3le
		
	call gauss(f4,0.d0,rkfe,10,re4le)
	presslep=presslep+re4le
	rma=rml(2)
		
	call gauss(f3,0.d0,rkfmu,10,re3lmu)
	enerlep=enerlep+re3lmu
		
	call gauss(f4,0.d0,rkfmu,10,re4lmu)		
	presslep=presslep+re4lmu
	ener=enerf+enerbar+enerlep		!soma de todos os termos de energia
	press=pressf+pressbar+presslep		!soma de todos os termos de pressao

	return
end subroutine
	
	!------------------------------------------------------------------
	REAL(SP) pure function f3(x)
	use nrtype
	implicit none
        REAL(SP), intent(in) :: x
	REAL(SP):: e,rma
	REAL(SP):: pi2 

	e=sqrt(x*x+rma*rma)	
	f3=x*x*e/pi2		!dens energia 5.36 sem mesons

	return
	end function

	!--------------------------------------------------------
	REAL(SP) pure function f4(x)
	use nrtype
	implicit none
        REAL(SP), intent(in) :: x
	REAL(SP):: e,rma
	REAL(SP):: pi2 

	e=sqrt(x*x+rma*rma)		!
	f4=x**4/e/(3.d0*pi2)		!pressao 5.37

	return
	end function

!---------------------------------------------------------------
subroutine gauss(f,ug,og,nn,fxint)		!programa integracao
	use nrtype
	implicit none
	REAL(SP):: x, o, u, ug, ou, og, fxint, ri, hh,f
	REAL(SP), DIMENSION(10):: wg,zg
	integer nh, i, kk, nn

	
	wg = [0.6667134430869d-01,0.1494513491506d00,0.2190863625160d00,0.2692667193100d00,0.2955242247148d00,&
	0.2955242247148d00,0.2692667193100d00,0.2190863625160d00,0.1494513491506d00,0.6667134430869d-01]
	zg = [-0.9739065285172d00,-0.8650633666890d00,-0.6794095682990d00,-0.4333953941292d00,-0.1488743389816d00,&
	+0.1488743389816d00,+0.4333953941292d00,+0.6794095682990d00,+0.8650633666890d00,+0.9739065285172d00]

	fxint=0.d0
	hh=(og-ug)/dble(float(nn))
	u=ug
	o=u+hh
	kk=1
24	ou=o+u
	ri=0.d0
	
	do 26 i=1,10
		x=0.5_sp*(zg(i)*hh+ou)
26		ri=ri+wg(i)*f(x)
		fxint=ri*hh/2.0_sp+fxint
		kk=kk+1
		
		if(kk-nn)28,28,9999
28			u=o
			o=o+hh
		go to 24
9999    return
end subroutine
!------------------------------------------
!-----------------------------------------------------------------------
subroutine funcv(n,x,fvec)
!     esta subrotinas é utilizada pelo metodo de broydn e da as funcoes para ser zero
! 	  no vetor fvec

	use nrtype
	implicit none 
	logical check
	REAL(SP), DIMENSION(5):: x,fvec			!x sao 5 variaveis, fvec equacoes
	REAL(SP), DIMENSION(8):: rmb,ri3,gs,gv,gr,dr,sb,qb,rnb	!
	REAL(SP), DIMENSION(2):: rnl,rml
	REAL(SP):: rlam,csi,rnbt,schargeb,rma,rnu,rmun,rmue,rmu,rkfb,rkfe,rkfn,rkfb2,dense,densb,densmu,rkfmu,rkfmu2,rk
	REAL(SP):: fsigma,fomega,frho,charge,rnumber,fsigmanl,fomeganl,rmbs,efb,efn,rint,rmns,rmr,rms,rmv
	REAL(SP):: re1,vomega,vrho,vsigma
	REAL(SP), parameter:: pi2=pi*pi
	integer:: nh, n, i, j, numberh
	
	call mapping(x,rkfe,rkfn,vsigma,vomega,vrho)

	fsigma=0.d0
	fomega=0.d0
	frho=0.d0
	charge=0.d0                           !electric charge density
	rnumber=0.d0                           !baryonic number density 
	schargeb=0.0

	fsigmanl=-rk/(2.d0*rms**2)*vsigma**2/gs(1) -rlam/(6.d0*rms**2)*vsigma**3/gs(1)**2  	
	!funcao em termos nao lineares sigma nl eq 5.39
      
	fomeganl=-csi*gv(1)**4/(6.d0*rmv**2)*vomega**3/gv(1)**2  	
	!funcao em termos nao lineares omega nl 5.38

	do i=1,nh
		call skfb2(i,rkfe,rkfn,vsigma,vomega,vrho,rkfb2)         

		if(rkfb2.gt.0.d0)then                                  
			rkfb=sqrt(rkfb2)     
		else
                 	rkfb=0.d0
            	endif                   

		rmbs=rmb(i)-vsigma*gs(i)/gs(1)           !baryon effective mass

		efb=sqrt(rkfb2+rmbs**2)     !baryon fermi energy

		rint=0.5d0*(rkfb*efb-rmbs**2*log((rkfb+efb)/rmbs))  
 
		rma=rmbs
		rmns=rmb(1)-vsigma   
		efn=sqrt(rkfn**2+rmns**2)                    
		rmue=sqrt(rkfe**2+rml(1)**2)                  
		rmun=vomega+ri3(1)*vrho+efn                 
		rmu=rmun-qb(i)*rmue
		rnu=rmu-gv(i)/gv(1)*vomega-gr(i)/gr(1)*ri3(i)*vrho
               
		call gauss(f1,0.d0,rkfb,10,re1)
               
		fsigma=fsigma+gs(i)*gs(1)/rms**2*re1 

		densb=rkfb**3/3.d0/pi2

		fomega=fomega+gv(i)*gv(1)/rmv**2*densb    		!omega equation of motion 5.38
		frho=frho+ri3(i)*gr(i)*gr(1)/rmr**2*densb  		!rho equation of motion 5.40
                  
		charge=charge+qb(i)*densb      
		rnumber=rnumber+densb
		schargeb=schargeb+abs(sb(i))*densb/3.d0                  
		rnb(i)=densb
          
	enddo      

	dense=rkfe**3/3.d0/pi2
	charge=charge-dense 
	rkfmu2=rkfe**2+rml(1)**2-rml(2)**2
                                           
	if(rkfmu2.gt.0.d0)then                                  
		rkfmu=sqrt(rkfmu2)     
	else
		rkfmu=0.d0
	endif
	     
	densmu=rkfmu**3/3.d0/pi2

	charge=charge-densmu

	rnl(1)=dense
	rnl(2)=densmu


    	fsigma=fsigma+fsigmanl
    	fomega=fomega+fomeganl

    	fvec(1)=fsigma-vsigma
    	fvec(2)=fomega-vomega
    	fvec(3)=frho-vrho
    	fvec(4)=charge
    	fvec(5)=rnumber-rnbt

	return
	end subroutine

	!-----------------------------------------------------------------------
	REAL(SP) pure function f1(x)
	use nrtype
	implicit none
        REAL(SP), intent(in) :: x
	REAL(SP):: e,rma
	REAL(SP):: pi2 

	e=sqrt(x*x+rma*rma)
	f1=x*x/e*rma/pi2

	return
	end function

