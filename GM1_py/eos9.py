'''
#==================================================================
      #legenda
      #rmb = vetor com 8 componentes - fermions - r real mb massa barions
      #rml = vetor com 2 componentes - leptons - r real mb massa barions
      #qb = vetor c 8 compont - cargas
      #i, j, k, l, m, n = contadores inteiros
      #35 99130-6304  ==  luiz zap
#==================================================================	  
#para rodar os codigos:	  

# gfortran wal.f90 -w -o executavel
#==================================================================
#indices das especies das particulas nucleares
#     1=neutron,2=proton,3=lambda,4=sigma-,5=sigma0,6=sigma+,7=cascade-,8=cascade0
#     1=electron,2=muon
#massas

#mb(1) = 939.0, mb(2) = 939.0, mb(3) = 1116.0, mb(4) = 1193.0, mb(5) = 1193.0,  mb(6) = 1193.0
#mb(7) = 1318.0,  mb(8) = 1318.0
#ml(1) = 0.5110, ml(2) = 105.660

#cargas
#qb(1)=0.00, qb(2)=1.00, qb(3)=0.00, qb(4)=-1.00, qb(5)0.00, qb(6)=1.00, qb(7)=-1.00, qb(8)=0.00
#ql(1)=1.00, ql(2)=1.00,

#isospin barions
#i3(1)=-0.50, i3(2)=0.50, i3(3)=0.00, i3(4)=-1.00, i3(5)=0.00, i3(6)=1.00, i3(7)=-0.50, i3(8)=0.50

'''
#--------------------------------------------------------------------------
# initial data

import numpy as np

#--------------------------------------------------------------------------
def mapping(x):		#esta sub-rotina define os campos potenciais entre os valores físicos   

	x = np.empty(n)	
	x = np.array([x[0],x[1],x[2],x[3],x[4]])
	
	kfe = x[0]**2                            	#to assure kfe greater than zero
	kfn = x[1]**2                            	#to assure kfn greater than zero
	vsigma = 0.9999990*np.sin(x[2])**2   		#mapping vsigma in (0,mb(1))
	vomega = x[3]**2                       		#to assure vomega greater than zero
	vrho = x[4]                              	#no mapping on vrho
      
	return kfe,kfn,vsigma,vomega,vrho 
         
#-----------------------------------------------------------------------
def skfb2(i,kfe,kfn,vsigma,vomega,vrho,kfb2):	# bug é chamado o kfb2 na função mas ele é também chamado globalmente
#esta sub-rotina calcula o momento fermi ao quadrado do i-ésimo barion usando as restrições do potencial químico

	global mb;ml;qb;i3;sb;gs;gv;gr;nh;pi2;k;ms;lam;csi
	
	
	kfe2 = kfe**2								#quadrado de  kfe=x(1)**2
	kfn2 = kfn**2								#quadrado de  kfn=x(2)**2

	mns = mb[0] - vsigma            			#neutron effective mass  eq 4.4    mns=rmb(1)-vsigma*gs(1)/gs(1)
	mbs = mb[i] - vsigma*gs[i]/gs[0]          	#baryon effective mass  eq 4.4

	efn = np.sqrt(kfn2 + mns**2)                #neutron fermi energy 3.7
	mue = np.sqrt(kfe2 + ml[0]**2)                  	#electron chemical potential  3.7 mas para o lepton eletron
	mum = vomega + i3[0]*vrho + efn                 	#neutron chemical potential 5.7
	mub = mum - qb[i]*mue                         	#baryon chemical potential eq 5.33
	kfb2 = (mub-gv[i]/gv[0]*vomega - gr[i]/gr[0]*i3[i]*vrho)**2 - mbs**2  #baryon fermi mom squared
	  #kfb^2 = (mu_b - g_omega*omega - g_rho*i3*rho)^2 - mb_sigma*^2   eq 5.7 isolando kfb^2

	return i,efn,mue,mum,mub,vrho,kfb2
   
#-----------------------------------------------------------------------
def fvec(x,n):
#     esta subrotinas é utilizada pelo metodo de broyden e da as funcoes para ser zero no vetor fvec
	global mb;ml;qb;i3;sb;gs;gv;gr;nh;pi2;k;ms;lam;csi

	x = np.empty(n) #retorna um tensor preenchido com dados não inicializados
	x = np.array([x[0],x[1],x[2],x[3],x[4]])
	fv = np.empty(n)
	
	n = len(x)
	
	kfe,kfn,vsigma,vomega,vrho = mapping(x)
	
	fsigma = 0.001
	fomega =  0.001
	#fsigmanl =  0.001
	#fomeganl = 0.001	
	frho = 0.001
	charge = 0.001                           #electric charge density
	rnumber =0.001                           #baryonic number density 
	
	fsigmanl = -k/(2.0*ms**2)*vsigma**2/gs[0] - lam/(6.0*ms**2)*vsigma**3/gs[0]**2  			#funcao em termos nao lineares sigma nl eq 5.39
	fomeganl = -csi*gv[0]**4/(6.0*mv**2)*vomega**3/gv[0]**2  	#funcao em termos nao lineares omega nl 5.38

	global kfb2

	global nb
	nb = np.empty(nh)
	
	global nl
	nl = np.empty(2)
	
	for i in range(0,nh):
		
		i,efn,mue,mum,mub,vrho,kfb2 = skfb2(i,kfe,kfn,vsigma,vomega,vrho,kfb2)         
		
		if kfb2 > 0.0:                                  
			kfb = np.sqrt(kfb2)     
		else:
			kfb = 0.0
		
		mbs = mb[i] - vsigma*gs[i]/gs[0]           	#baryon effective mass
		efb = np.sqrt(kfb2 + mbs**2)     		#baryon fermi energy
#		rint = 0.50*(kfb*efb - mbs**2*np.log((kfb + efb)/mbs))  
		rma = mbs
		mns = mb[0] - vsigma   
		efn = np.sqrt(kfn**2 + mns**2)                    
		mue = np.sqrt(kfe**2 + ml[0]**2)                  
		mum = vomega + i3[0]*vrho + efn                 
#		rmu = mum - qb[i]*mue
#		rnu = rmu - gv[i]/gv[0]*vomega - gr[i]/gr[0]*i3[i]*vrho
		
		re1 = gauss_legendre(f1(kfb, rma),0.0,kfb,10)      ##call gauss(f1,0.0,kfb,10,re1,ii) bug alterei o f1 pela função, e adicionei x e rma, pode estar errado
		
		fsigma = fsigma + gs[i]*gs[0]/ms**2*re1 
		densb = kfb**3/3.0/pi2
		fomega = fomega + gv[i]*gv[0]/mv**2*densb    		#omega equation of motion 5.38
		frho = frho + i3[i]*gr[i]*gr[0]/mr**2*densb  		#rho equation of motion 5.40
		charge = charge + qb[i]*densb      
		rnumber = rnumber + densb
		nb[i] = densb

	dense = kfe**3/3.0/pi2
	charge = charge - dense
	kfmu2 = kfe**2 + ml[0]**2 - ml[1]**2
		
	if kfmu2 > 0.0:                                  
		kfmu = np.sqrt(kfmu2)     
	else:
		kfmu = 0.0
		
	densmu = kfmu**3/3.0/pi2
	charge = charge - densmu
	nl[0] = dense
	nl[1] = densmu
			
	fsigma = fsigma + fsigmanl
	fomega = fomega + fomeganl
	
	fv[0] = fsigma - vsigma 
	fv[1] = fomega - vomega
	fv[2] = frho  - vrho 
	fv[3] = charge 
	fv[4] = rnumber - nbt 
	
	return np.array([fv[0],fv[1],fv[2],fv[3],fv[4]])
			
#-----------------------------------------------------------------------
def eos(kfe,kfn,vsigma,vomega,vrho,ener,press,phi,v0):
#	esta subrotina calcula o eos desta matéria estelar

	global mb; ml; qb; i3; sb; gs; gv; gr; x; nh; n; pi2; k; ms; lam; csi

	phi = vsigma/gs[0]
	v0 = vomega/gv[0]
	b0 = -vrho/gr[0]

	#dens de energia no momento de fermi 5.36
	enerf = mv**2*v0**2/2.0 + csi*gv[0]**4*v0**4/8.0 + mr**2*b0**2/2.0 + ms**2*phi**2/2.0 + k*phi**3/6.0 + lam*phi**4/24.0

	#pressao no momento de fermi 5.37
	pressf = mv**2*v0**2/2.0 + csi*gv[0]**4*v0**4/24.0 + mr**2*b0**2/2.0 - ms**2*phi**2/2.0 - k*phi**3/6.0 - lam*phi**4/24.0

	enerbar = 0.01		#energia barionica inicial
	pressbar = 0.01		#pressao barionica inicial
	
	global kfb2
	
	for i in range(0,nh):
	
		i,efn,mue,mum,mub,vrho,kfb2 = skfb2(i,kfe,kfn,vsigma,vomega,vrho,kfb2)          
		
		if kfb2 > 0.0:                                  
			kfb=np.sqrt(kfb2)     	# estava dsqrt(kfb2)
		else:
			kfb=0.0
  
		mbs = mb[i] - vsigma*gs[i]/gs[0]               
		rma = mbs
		mns = mb[0] - vsigma   
		efn = np.sqrt(kfn**2 + mns**2)                    
		mue = np.sqrt(kfe**2 + ml[0]**2)                  
		mum = vomega + i3[0]*vrho + efn                 
		rmu = mum - qb[i]*mue
		rnu = rmu - gv[i]/gv[0]*vomega - gr[i]/gr[0]*i3[i]*vrho

		re3 = gauss_legendre(f3(kfb, rma),0.0,kfb,10)	#call gauss(f3,0.0,kfb,10,re3,ii)
		enerbar = enerbar + re3
 		
		re4 = gauss_legendre(f4(kfb, rma),0.0,kfb,10) 	#call gauss(f4,0.0,kfb,10,re4,ii)
		pressbar = pressbar + re4

	enerlep = 0.001
	presslep = 0.001
	kfmu = 0.001                         
	kfmu2 = kfe**2 + ml[0]**2 - ml[1]**2                           
                                           
    
	if kfmu2 > 0.0:
		kfmu = np.sqrt(kfmu2)      
		mue = np.sqrt(kfe**2 + ml[0]**2)   
		rnu = mue
		rma = ml[0]
	
	re3le = gauss_legendre(f3(kfb, rma),0.0,kfe,10)	#call gauss(f3,0.0,kfe,10,re3le,ii)
	print(re3le)
	
	enerlep = enerlep + re3le
	
	re4le = gauss_legendre(f4(kfb, rma),0.0,kfe,10)	#call gauss(f4,0.0,kfe,10,re4le,ii)
	
	presslep = presslep + re4le
	rma = ml[1]
	
	re3lmu = gauss_legendre(f3(kfb, rma),0.0,kfmu,10)	#call gauss(f3,0.0,kfmu,10,re3lmu,ii)
	
	enerlep = enerlep + re3lmu
	
	re4lmu = gauss_legendre(f4(kfb, rma),0.0,kfmu,10)	#call gauss(f4,0.0,kfmu,10,re4lmu,ii)		
	
	presslep = presslep + re4lmu
	ener = enerf + enerbar + enerlep			#soma de todos os termos de energia
	press = pressf + pressbar + presslep		#soma de todos os termos de pressao

	return kfe,kfn,vsigma,vomega,vrho,ener,press,presslep,enerlep

#-----------------------------------------------------------------------
def f1(kr, rma):
	global pi2
	f1 = kr*kr*rma/(pi2*np.sqrt(kr*kr + rma*rma))
	return f1

#------------------------------------------------------------------
def f3(kr, rma):

	global pi2
	
	f3 = kr*kr*np.sqrt(kr*kr + rma*rma)*pi2		#dens energia 5.36 sem mesons

	return f3

#--------------------------------------------------------
def f4(kr, rma):

	global pi2

	f4 = kr**4/(np.sqrt(kr*kr + rma*rma)*(3.0*pi2))		#pressao 5.37

	return f4

#---------------------------------------------------------------
#função que resolve integrais via método de gauss_legendre

def gauss_legendre(f,a,b,n_int):		#programa integracao usando o método de gauss-legendre
	w = [0.6667134430869e-01,0.149451349150600,0.219086362516000,0.269266719310000,0.295524224714800,0.295524224714800,0.269266719310000,0.219086362516000,0.149451349150600,0.6667134430869e-01]
	
	t = [-0.973906528517200, -0.865063366689000,-0.679409568299000,-0.433395394129200,-0.148874338981600,+0.148874338981600,+0.433395394129200,+0.679409568299000,+0.865063366689000,+0.973906528517200]

	int = 0.0			#valor inicial a ser integrado
	h = (b - a)/2
	
	for i in range(0,10):
		x = 0.50*(t[i]*(b - a) + (b + a))
		int = int + h*w[i]*f	# bug a função f(x) não existe?
		
	return int


#from numpy import dot

#------------------------------------------
#função que resolve sistema não-linear via método de broyden
def gauss_jordan(a, b):
    n = len(b)
    for k in range(0,n-1): 
        for i in range(k+1,n): 
            if a[i][k] != 0.0: 
                m = a[i][k]/a[k][k] 
                for j in range(k + 1, n):
                    a[i][j] = a[i][j] - m*a[k][j] 
                b[i] = b[i] - m*b[k] 
    for k in range(n-1,-1,-1):   #Começa no último "n-1", até o zero "-1" e de sima para baixo "-1"
        b[k] = (b[k] - np.dot(a[k][k:n],b[k:n]))/a[k][k]
    return b

# bug estava com um erro no np.abs() - não entendi e refiz manualmente

#import scipy.linalg as la

#print ('via scipy   X =', la.solve(C,d))

#função que resolve o método de broyden para sist não-lineares	
def broyden(xk, n, f, B, eps, k_max):
	
	k = 0
	fk = f(xk,n)
	B = B(xk)
#Função é usada para calcular uma das oito normas de matrizes diferentes ou uma das normas vetoriais.
	f_norma = np.linalg.norm(fk, ord=2)   

	while abs(f_norma) > eps and k < k_max:
		
		s = gauss_jordan(B, -1*fk)
		
		#s = sla.solve(B, -1*fk)
        
		xk = xk + s
		f_novo = f(xk,n)
		
		df = f_novo - fk
        
		B = B + (np.outer((df - np.dot(B,s)),s)) / (np.dot(s,s))
        
		f_norma = np.linalg.norm(fk, ord=2)
		
		fk = f_novo
        
		k = k + 1
		
	return xk

#testar matriz identidade
def B(x):
	return np.array([[1.0,0.0,0.0, 0.0,0.0],[0.0,1.0,0.0, 0.0,0.0],[0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,1.0]])

#----------------------------------------------------------------------------

print("estamos rodando matéria estelar")

file = open('eos.dat', 'w')

nh = 8

mb = np.array([939.0,939.0,1116.0,1193.0,1193.0,1193.0,1318.0,1318.0])  # definindo vetor de 8 componentes
ml = np.array([0.5110,105.660])
qb = np.array([0.00,1.00,0.00,-1.00,0.00,1.00,-1.00,0.00])
i3 = np.array([-0.50,0.50,0.00,-1.00,0.00,1.00,-0.50,0.50])
sb = np.array([0.00,0.00,-1.0,-1.0,-1.0,-1.0,-2.0,-2.0])


# bug defini o kfb2 apenas para não dar erro
kfb2 = 0.001
kfmu2 = 0.001
ener = 0.001
press = 0.001
kfe = 0.001
kfn = 0.001
vsigma = 0.001
vomega =0.001
vrho = 0.001
phi = 0.001
v0 = 0.001

gs = np.empty(nh)
gv = np.empty(nh)
gr = np.empty(nh)

n = 5
x = np.empty(n)
x = np.array([x[0],x[1],x[2],x[3],x[4]])

eps = 0.000001
k_max = 1000

pi = np.pi; pi2 = pi**2   	#pi2/9.869604410090

# data input

rm = 939.0				#nucleon mass		  proton e neutron
ms = 512.0				#scalar meson mass	   m_sigma
mv = 783.0				#vector meson mass	   m_omega
mr = 770.0				#vector-isovector meson mass	   m_rho
xls = 11.7850 				#gsn^2/ms^2		xls 
xlw = 7.14800  			# gwn^2/mw^2		xlw
xlr = 4.9390  				# grn^2/mr^2		xlr
k = 0.0058940				# k/m			k	termo cubico n linear
lam = -0.006420   			#lambda			lam	termo quartico n linear
csi = 0.0000001				#csi			csi
xs = 0.7				#gs/gn			xs sigma-hiperon coupling constant
xw = 0.783      			# gw/gn 			xw omega-hiperon coupling constant
xr = 0.783     			# gr/gn 			xr rho-hiperon coupling constant

nbinf = 0.030				#nbinf [fm^-3]		nbinf	range limite inferior
nbsup = 1.230				#nbsup [fm^-3]		nbsup	range limite superior
npoint = 400				#number of points	npoint

ms = ms/rm 				#admensionalidade massa do meson escalar m_sigma
mv = mv/rm  				#admensionalidade massa do meson escalar m_omega     
mr = mr/rm 				#admensionalidade massa do meson escalar m_rho
nbinf = nbinf/(rm/197.320)**3		#admensionalidade limite inferior
nbsup = nbsup/(rm/197.320)**3		#admensionalidade limite inferior
dnb = (nbsup-nbinf)/(npoint-1)		#numero de passos

for i in range(2,8):
	mb[i] = mb[i]/rm		#admensionalidade barions/hiperions

for i in range(0,2):			#admensionalidade barions/nucleons e dos leptons
	mb[i] = 1.0
	ml[i] = ml[i]/rm
		
xls = xls*(rm/197.330)**2		#admensionalidade da const acoplamento meson scalar sigma
xlw = xlw*(rm/197.330)**2		#admensionalidade da const acoplamento meson scalar omega
xlr = xlr*(rm/197.330)**2		#admensionalidade da const acoplamento meson scalar rho

for i in range(0,2):
	gs[i] = np.sqrt(xls)*ms		#const acoplamento
	gv[i] = np.sqrt(xlw)*mv
	gr[i] = np.sqrt(xlr)*mr

for i in range(2,8):
	gs[i] = xs*gs[0]    		#sigma-hiperon coupling constant
	gv[i] = xw*gv[0]    		#omega-hiperon coupling constant
	gr[i] = xr*gr[0]    		#rho-hiperon coupling constant

k = k*gs[0]**3			#constante k do termo sigma cubico do lagrangeano
lam = lam*gs[0]**4			#constante lambda do termo sigma quartico do lagrangeano

# initializing variables		#condicoes iniciais das variaveis

ne0 = 0.10*nbinf				
nn0 = 0.90*nbinf

x[0] = (3.0*pi2*ne0)**0.16666666670
x[1] = (3.0*pi2*nn0)**0.16666666670
x[2] = np.sin(np.sqrt(0.0050/0.9999990))     	#     x(4)=gvn/mv*dsqrt(nbinf)
x[3] = np.sqrt(0.22731e-5)			#     x(5)=(grn/mr)**2*(i3(1)*nn0+i3(2)*ne0)
x[4] = -0.57641e-6

x = np.array([x[0],x[1],x[2],x[3],x[4]])

#x = np.array([1.0,1.0,1.0,1.0,1.0])
#--------------------------------------------------------------------------

# loop principal que calcula as equações de estado
for i in range(0,npoint):
	nbt = nbinf + (i - 1)*dnb
	
	n = len(x)
	
	x = broyden(x, n, fvec, B, eps, k_max)
		
	kfe,kfn,vsigma,vomega,vrho = mapping(x)
	
	kfe,kfn,vsigma,vomega,vrho,ener,press,presslep,enerlep = eos(kfe,kfn,vsigma,vomega,vrho,ener,press,v0,phi)    #the eos

	nbtd = nbt*(rm/197.320)**3
	mns = (mb[0] - vsigma)

	mum = vomega + i3[0]*vrho + np.sqrt((3.0*pi2*nb[0])**0.66666666670 + mns**2)   	#neutron chemical potential
	mue = np.sqrt((3.0*pi2*nl[0])**0.66666666670 + ml[0]**2) 				#electron chemical potential
         
	dener = ener/nbt - 1.0   		#dens de energia de ligação - massa (valor 1)
	dener = dener*rm			#multiplicado por rm para ficar em mev

#     output in fm^-4 for the program that calculates stars properties; 
#     if you choose this output, comment the above one
	press = press*(rm/197.32)**4			#pressao em fm^-4
	ener = ener*(rm/197.32)**4			#energia em fm^-4

#	file.write("%.10f %.10f %.10f \n" %(nbtd,ener,press,mum*rm))
#	write(22,10) nbtd, (nl(j)/nbt,j=1,2), (nb(j)/nbt,j=1,nh) xxxxxxxxxxx acertar xxxxxxxxxxxxx

	file.write("%.10f %.10f %.10f %.10f\n" %(nbtd,ener,press,mum*rm))	
  
	rmpu = mb[0]*nb[0] + mb[1]*nb[1] + mb[2]*nb[2] + mb[3]*nb[3]+ mb[4]*nb[4] + mb[5]*nb[5] + mb[6]*nb[6] + mb[7]*nb[7]

	rmpum = ((rmpu/(nbtd)))*(939.0/197.330)**3
	rmpummev = rmpum*939.0                                #massa média por partícula
	rmpumfm =  rmpummev/197.330

	
