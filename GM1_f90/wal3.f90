!==================================================================
      !legenda
      !rmb = vetor com 8 componentes - fermions - r real mb massa barions
      !rml = vetor com 2 componentes - leptons - r real mb massa barions
      !qb = vetor c 8 compont - cargas
      !i, j, k, l, m, n = contadores inteiros
      !35 99130-6304  ==  luiz zap
!==================================================================	  
!para rodar os codigos:

!indices das especies das particulas nucleares
!     1=neutron,2=proton,3=lambda,4=sigma-,5=sigma0,6=sigma+,7=cascade-,8=cascade0
!     1=electron,2=muon
!massas
!mb(1) = 939.0, mb(2) = 939.0, mb(3) = 1116.0, mb(4) = 1193.0, mb(5) = 1193.0,  mb(6) = 1193.0
!mb(7) = 1318.0,  mb(8) = 1318.0
!ml(1) = 0.5110, ml(2) = 105.660
		
!cargas
!qb(1)=0.00, qb(2)=1.00, qb(3)=0.00, qb(4)=-1.00, qb(5)0.00, qb(6)=1.00, qb(7)=-1.00, qb(8)=0.00
!ql(1)=1.00, ql(2)=1.00,
		
!isospin barions
!i3(1)=-0.50, i3(2)=0.50, i3(3)=0.00, i3(4)=-1.00, i3(5)=0.00, i3(6)=1.00, i3(7)=-0.50, i3(8)=0.50
	
	  

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
REAL(SP), dimension(8) :: rnb

!dados fisicos das particulas  poderia ser definido direto em :  real rmb(8)/ xxx,xxx, ... /
rmb = [939.0,939.0,1116.0,1193.0,1193.0,1193.0,1318.0,1318.0]  ! definindo vetor de 8 componentes
rml = [0.5110,105.660]
qb = [0.00,1.00,0.00,-1.00,0.00,1.00,-1.00,0.00]
ri3 = [-0.50,0.50,0.00,-1.00,0.00,1.00,-0.50,0.50]
sb = [0.00,0.00,-1.0,-1.0,-1.0,-1.0,-2.0,-2.0]
	
		
write(6,*)'2 for nucleons and 8 for nucleons + hyperons'
read(5,*) nh 
write(6,*) 'vc está rodando matéria estelar'

!data input

open(1,file='gm1.inp')		!lendo os dados do arquivo gm1.inp
read(1,*) rm				!939.0	nucleon mass		  proton e neutron
read(1,*) rms				!512.0	scalar meson mass	   m_sigma
read(1,*) rmv				!783.0	vector meson mass	   m_omega
read(1,*) rmr				!770.0	vector-isovector meson mass	   m_rho
read(1,*) xls				!11.7850 	gsn^2/ms^2		xls 
read(1,*) xlw				!7.14800  	gwn^2/mw^2		xlw
read(1,*) xlr				!4.9390  	grn^2/mr^2		xlr
read(1,*) rk				!0.0058940	k/m			rk	termo cubico n linear
read(1,*) rlam				!-0.006420   	lambda			rlam	termo quartico n linear
read(1,*) csi				!0.0		csi			csi
read(1,*) xs				!0.7		gs/gn			xs sigma-hiperon coupling constant
read(1,*) xw				!0.783       	gw/gn 			xw omega-hiperon coupling constant
read(1,*) xr				!0.783       	gr/gn 			xr rho-hiperon coupling constant
read(1,*) rnbinf			!0.030	nbinf [fm^-3]		rnbinf	range limite inferior
read(1,*) rnbsup			!1.230	nbsup [fm^-3]		rnbsup	range limite superior
read(1,*) npoint			!400		number of points	npoint

close(1)

rms=rms/rm 				!admensionalidade massa do meson escalar m_sigma
rmv=rmv/rm  				!admensionalidade massa do meson escalar m_omega     
rmr=rmr/rm 				!admensionalidade massa do meson escalar m_rho
rnbinf=rnbinf/(rm/197.320)**3		!admensionalidade limite inferior
rnbsup=rnbsup/(rm/197.320)**3		!admensionalidade limite inferior
dnb=(rnbsup-rnbinf)/(npoint-1)		!numero de passos

do i=3,8
	rmb(i)=rmb(i)/rm		!admensionalidade barions/hiperions
enddo

do i=1,2				!admensionalidade barions/nucleons e dos leptons
	rmb(i)=1.0
	rml(i)=rml(i)/rm
enddo
		
xls=xls*(rm/197.330)**2		!admensionalidade da const acoplamento meson scalar sigma
xlw=xlw*(rm/197.330)**2		!admensionalidade da const acoplamento meson scalar omega
xlr=xlr*(rm/197.330)**2		!admensionalidade da const acoplamento meson scalar rho

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

rne0=0.10*rnbinf				
rnn0=0.90*rnbinf      
x(1)=(3.0*pi2*rne0)**0.16666666670
x(2)=(3.0*pi2*rnn0)**0.16666666670
x(3)=asin(sqrt(0.0050/0.9999990))     	!     x(4)=gvn/rmv*sqrt(rnbinf)
x(4)=sqrt(0.22731d-05)			!     x(5)=(grn/rmr)**2*(ri3(1)*rnn0+ri3(2)*rne0)
x(5)=-0.57641d-06

!main loop 	!loop principal que calcula as equações de estado

do i=1,npoint
	rnbt=rnbinf+(i-1)*dnb

            
	call broydn(x,check)
		if(check)pause 'there is no root in broydn...'

	call mapping(x,rkfe,rkfn,vsigma,vomega,vrho)

	call eos(rkfe,rkfn,vsigma,vomega,vrho,ener,press,rv0,rphi)    !the eos
 
	rnbtd=rnbt*(rm/197.320)**3
	rmns=(rmb(1)-vsigma)

	rmun=vomega+ri3(1)*vrho + sqrt((3.0*pi2*rnb(1))**0.66666666670+rmns**2)   !neutron chemical potential
	rmue=sqrt((3.0*pi2*rnl(1))**0.66666666670+rml(1)**2) !electron chemical potential
         
	dener=ener/rnbt-1.0   		!dens de energia de ligação - massa (valor 1)
	dener=dener*rm			!multiplicado por rm para ficar em mev


!     output in fm^-4 for the program that calculates stars properties; 
!     if you choose this output, comment the above one
	press=press*(rm/197.32)**4			!pressao em fm^-4
	ener=ener*(rm/197.32)**4			!energia em fm^-4


	write(22,10) rnbtd, (rnl(j)/rnbt,j=1,2), (rnb(j)/rnbt,j=1,nh)
	
	write(100,*) rnbtd,ener,press,rmun*rm
  
	rmpu = rmb(1)*rnb(1)+rmb(2)*rnb(2) +rmb(3)*rnb(3) + rmb(4)*rnb(4) &
	+ rmb(5)*rnb(5)+rmb(6)*rnb(6)+rmb(7)*rnb(7)  + rmb(8)*rnb(8)
             
	rmpummev = ((rmpu/(rnbtd)))*(939.0/197.330)**3*939.0             !massa média por partícula
	rmpumfm =  rmpummev/197.330
     
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
	vsigma=0.9999990*sin(x(3))**2   		!mapping vsigma in (0,mb(1))
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
	external f3,f4
	REAL(SP) :: rmns

	

	rphi=vsigma/gs(1)
	rv0=vomega/gv(1)
	b0=-vrho/gr(1)

	enerf=rmv**2*rv0**2/2.0+csi*gv(1)**4*rv0**4/8.0 &	!dens de energia no momento de fermi 5.36
	+ rmr**2*b0**2/2.0 &
	+ rms**2*rphi**2/2.0 + rk*rphi**3/6.0 &
	+ rlam*rphi**4/24.0

	pressf=rmv**2*rv0**2/2.0+csi*gv(1)**4*rv0**4/24.0 &  	!pressao no momento de fermi 5.37
	+rmr**2*b0**2/2.0 & 
	-rms**2*rphi**2/2.0 - rk*rphi**3/6.0 &
	-rlam*rphi**4/24.0

	enerbar=0.0		!energia barionica inicial
	pressbar=0.0		!pressao barionica inicial

	do i=1,nh
		call skfb2(i,rkfe,rkfn,vsigma,vomega,vrho,rkfb2)         
		if(rkfb2.gt.0.0)then                                  
			rkfb=sqrt(rkfb2)     
		else
			rkfb=0.0
		endif     

		rmbs=rmb(i)-vsigma*gs(i)/gs(1)               
		rma=rmbs
		rmns=rmb(1)-vsigma   
		efn=sqrt(rkfn**2+rmns**2)                    
		rmue=sqrt(rkfe**2+rml(1)**2)                  
		rmun=vomega+ri3(1)*vrho+efn                 
		rmu=rmun-qb(i)*rmue
		rnu=rmu-gv(i)/gv(1)*vomega-gr(i)/gr(1)*ri3(i)*vrho
	
		call gauss(f3,0.0,rkfb,10,re3)
		enerbar=enerbar+re3
         
		call gauss(f4,0.0,rkfb,10,re4)
		pressbar=pressbar+re4
	enddo

	enerlep=0.0
	presslep=0.0
	rkfmu=0.0                         
	kfmu2=rkfe**2+rml(1)**2-rml(2)**2                           
              
	if(rkfmu2.gt.0.0)rkfmu=sqrt(rkfmu2)      

		rmue=sqrt(rkfe**2+rml(1)**2)   
		rnu=rmue
		rma=rml(1)
		
	call gauss(f3,0.0,rkfe,10,re3le)
	enerlep=enerlep+re3le
		
	call gauss(f4,0.0,rkfe,10,re4le)
	presslep=presslep+re4le
	rma=rml(2)
		
	call gauss(f3,0.0,rkfmu,10,re3lmu)
	enerlep=enerlep+re3lmu
		
	call gauss(f4,0.0,rkfmu,10,re4lmu)		
	presslep=presslep+re4lmu
	ener=enerf+enerbar+enerlep		!soma de todos os termos de energia
	press=pressf+pressbar+presslep		!soma de todos os termos de pressao

	return
end subroutine
	
!------------------------------------------------------------------
REAL(SP) pure function f3(sk)
	use nrtype
	implicit none
    	REAL(SP), intent(in) :: sk
	REAL(SP):: e,rma
	REAL(SP):: pi2 

	e=sqrt(sk*sk+rma*rma)	
	f3=sk*sk*e/pi2		!dens energia 5.36 sem mesons

	return
end function

!--------------------------------------------------------
REAL(SP) pure function f4(sk)
	use nrtype
	implicit none
    	REAL(SP), intent(in) :: sk
	REAL(SP):: e,rma
	REAL(SP):: pi2 

	e=sqrt(sk*sk+rma*rma)		!
	f4=sk**4/e/(3.0*pi2)		!pressao 5.37

	return
end function

!------------------------------------------
SUBROUTINE gauss(f,a,b,n,int_g)
	use nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(10):: w, t
	REAL(SP):: a,b,f,sk,h,int_g
	INTEGER:: i,n

	w = [0.6667134430869e-01,0.149451349150600,0.219086362516000,0.269266719310000,0.295524224714800,0.295524224714800,&
	0.269266719310000,0.219086362516000,0.149451349150600,0.6667134430869e-01]
	
	t = [-0.973906528517200, -0.865063366689000,-0.679409568299000,-0.433395394129200,-0.148874338981600,+0.148874338981600,&
	+0.433395394129200,+0.679409568299000,+0.865063366689000,+0.973906528517200]

	h = (b - a)/2

	int_g = 0.0

	DO i=1, 10
		sk = 0.50*(t(I)*(b - a) + (b + a))
		int_g = int_g + h*w(I)*f(sk)
	END DO

!	RETURN
END SUBROUTINE
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
	REAL(SP):: re1,vomega,vrho,vsigma,f1
	REAL(SP), parameter:: pi2=pi*pi
	integer:: nh, n, i, j, numberh
	external f1
	
	call mapping(x,rkfe,rkfn,vsigma,vomega,vrho)

	fsigma=0.0
	fomega=0.0
	frho=0.0
	charge=0.0                           !electric charge density
	rnumber=0.0                           !baryonic number density 
	schargeb=0.0

	fsigmanl=-rk/(2.0*rms**2)*vsigma**2/gs(1) -rlam/(6.0*rms**2)*vsigma**3/gs(1)**2  	
	!funcao em termos nao lineares sigma nl eq 5.39
      
	fomeganl=-csi*gv(1)**4/(6.0*rmv**2)*vomega**3/gv(1)**2  	
	!funcao em termos nao lineares omega nl 5.38

	do i=1,nh
		call skfb2(i,rkfe,rkfn,vsigma,vomega,vrho,rkfb2)         

		if(rkfb2.gt.0.0)then                                  
			rkfb=sqrt(rkfb2)     
		else
                 	rkfb=0.0
            	endif                   

		rmbs=rmb(i)-vsigma*gs(i)/gs(1)           !baryon effective mass

		efb=sqrt(rkfb2+rmbs**2)     !baryon fermi energy

		rint=0.50*(rkfb*efb-rmbs**2*log((rkfb+efb)/rmbs))  
 
		rma=rmbs
		rmns=rmb(1)-vsigma   
		efn=sqrt(rkfn**2+rmns**2)                    
		rmue=sqrt(rkfe**2+rml(1)**2)                  
		rmun=vomega+ri3(1)*vrho+efn                 
		rmu=rmun-qb(i)*rmue
		rnu=rmu-gv(i)/gv(1)*vomega-gr(i)/gr(1)*ri3(i)*vrho
               
		call gauss(f1,0.0,rkfb,10,re1)
               
		fsigma=fsigma+gs(i)*gs(1)/rms**2*re1 

		densb=rkfb**3/3.0/pi2

		fomega=fomega+gv(i)*gv(1)/rmv**2*densb    		!omega equation of motion 5.38
		frho=frho+ri3(i)*gr(i)*gr(1)/rmr**2*densb  		!rho equation of motion 5.40
                  
		charge=charge+qb(i)*densb      
		rnumber=rnumber+densb
		schargeb=schargeb+abs(sb(i))*densb/3.0                  
		rnb(i)=densb
          
	enddo      

	dense=rkfe**3/3.0/pi2
	charge=charge-dense 
	rkfmu2=rkfe**2+rml(1)**2-rml(2)**2
                                           
	if(rkfmu2.gt.0.0)then                                  
		rkfmu=sqrt(rkfmu2)     
	else
		rkfmu=0.0
	endif
	     
	densmu=rkfmu**3/3.0/pi2

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
REAL(SP) pure function f1(sk)
	use nrtype
	implicit none
    	REAL(SP), intent(in) :: sk
	REAL(SP):: e,rma
	REAL(SP):: pi2 

	e=sqrt(sk*sk+rma*rma)
	f1=sk*sk/e*rma/pi2

	return
end function

