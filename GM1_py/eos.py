'''
#==================================================================
      #Legenda
      #RMB = vetor com 8 componentes - fermions - R real MB massa barions
      #RML = vetor com 2 componentes - leptons - R real MB massa barions
      #QB = vetor c 8 compont - cargas
      #I, J, K, L, M, N = contadores inteiros
      #35 99130-6304  ==  Luiz zap
#==================================================================	  
#Para rodar os codigos:	  

# gfortran wal.f90 -w -o executavel
#==================================================================
#Indices das especies das particulas nucleares
#     1=neutron,2=proton,3=lambda,4=sigma-,5=sigma0,6=sigma+,7=cascade-,8=cascade0
#     1=electron,2=muon
#massas

#MB(1) = 939.0, MB(2) = 939.0, MB(3) = 1116.0, MB(4) = 1193.0, MB(5) = 1193.0,  MB(6) = 1193.0
#MB(7) = 1318.0,  MB(8) = 1318.0
#ML(1) = 0.5110, ML(2) = 105.660
		
#Cargas
#QB(1)=0.00, QB(2)=1.00, QB(3)=0.00, QB(4)=-1.00, QB(5)0.00, QB(6)=1.00, QB(7)=-1.00, QB(8)=0.00
#QL(1)=1.00, QL(2)=1.00,
		
#Isospin barions
#I3(1)=-0.50, I3(2)=0.50, I3(3)=0.00, I3(4)=-1.00, I3(5)=0.00, I3(6)=1.00, I3(7)=-0.50, I3(8)=0.50

   
#PROGRAM BOMEOS
##     CALCULATES THE EOS FOR THE BARYON OCTET MODEL (NLWM) AT T=0
##     uses Glendenning couplig constants
##     BY CHIAPPARINI AND MODIFIED BY MENEZES

#IMPLICIT NONE

#LOGICAL check
#REAL(16), DIMENSION X(5),FVEC(5),sb(8)			#X sao 5 variaveis, FVEC equacoes
#REAL(16), DIMENSION RMB(8),RML(2),QB(8),RI3(8),GS(8),GV(8),GR(8), dr(8)	#
#REAL(16), DIMENSION RNB(8),RNL(2)
#INTEGER NH, I, J, K, L

#COMMON/CDATA/RMB,RML,RMS,RMV,RMR,RK,RLAM,CSI, GS,GV,GR,QB,RI3,sb
#COMMON/CNBT/RNBT,RNB,RNL
#COMMON/CFVEC/FVEC
#COMMON/SCH/SCHARGEB
#COMMON/LOOP/I
#COMMON/NUMBERH/NH

##dados fisicos das particulas  poderia ser definido direto em :  REAL RMB(8)/ xxx,xxx, ... /
'''

import numpy as np
import tqdm
import time

start = time.time()

print("Rodando Matéria Estelar")

file = open('eos.dat', 'w')

NH = int(input("Insert\n2 - FOR NUCLEONS\n8 - FOR NUCLEONS + HYPERONS\nYour choice: \n")) 

RMB = np.array([939.0,939.0,1116.0,1193.0,1193.0,1193.0,1318.0,1318.0])  # definindo vetor de 8 componentes
RML = np.array([0.5110,105.660])
QB = np.array([0.00,1.00,0.00,-1.00,0.00,1.00,-1.00,0.00])
RI3 = np.array([-0.50,0.50,0.00,-1.00,0.00,1.00,-0.50,0.50])
SB = np.array([0.00,0.00,-1.0,-1.0,-1.0,-1.0,-2.0,-2.0])

GS = np.zeros(NH)
GV = np.zeros(NH)
GR = np.zeros(NH)

pi = np.pi; PI2 = pi**2   	#PI2/9.869604410090

#     Data input

RM = 939.0					#NUCLEON MASS		  proton e neutron
RMS = 512.0					#SCALAR MESON MASS	   M_sigma
RMV = 783.0					#VECTOR MESON MASS	   M_omega
RMR = 770.0					#VECTOR-ISOVECTOR MESON MASS	   M_rho
XLS = 11.7850 				#gsN^2/ms^2		XLS 
XLW = 7.14800  				# gwN^2/mw^2		XLW
XLR = 4.9390  				# grN^2/mr^2		XLR
RK = 0.0058940				# K/M			RK	termo cubico n linear
RLAM = -0.006420   			#LAMBDA			RLAM	termo quartico n linear
CSI = 0.0					#CSI			CSI
XS = 0.7					#gs/gN			XS sigma-hiperon coupling constant
XW = 0.783      			# gw/gN 			XW omega-hiperon coupling constant
XR = 0.783     				# gr/gN 			XR rho-hiperon coupling constant

RNBINF = 0.030				#nbinf [fm^-3]		RNBINF	range limite inferior
RNBSUP = 1.230				#nbsup [fm^-3]		RNBSUP	range limite superior
NPOINT = 400				#number of points	NPOINT

RMS = RMS/RM 				#admensionalidade massa do meson escalar M_sigma
RMV = RMV/RM  				#admensionalidade massa do meson escalar M_omega     
RMR = RMR/RM 				#admensionalidade massa do meson escalar M_rho
RNBINF = RNBINF/(RM/197.320)**3		#admensionalidade limite inferior
RNBSUP = RNBSUP/(RM/197.320)**3		#admensionalidade limite inferior
DNB = (RNBSUP-RNBINF)/(NPOINT-1)		#numero de passos

for I in range(3,8):
	RMB[I] = RMB[I]/RM		#admensionalidade barions/hiperions

for I in range(1,2):			#admensionalidade barions/nucleons e dos leptons
	RMB[I] = 1.0
	RML[I] = RML[I]/RM
		
XLS = XLS*(RM/197.330)**2		#admensionalidade da const acoplamento meson scalar sigma
XLW = XLW*(RM/197.330)**2		#admensionalidade da const acoplamento meson scalar omega
XLR = XLR*(RM/197.330)**2		#admensionalidade da const acoplamento meson scalar rho

for I in tqdm.tqdm(range(1,2)):
	GS[I] = np.sqrt(XLS)*RMS		#const acoplamento
	GV[I] = np.sqrt(XLW)*RMV
	GR[I] = np.sqrt(XLR)*RMR

for I in tqdm.tqdm(range(3,8)):
	GS[I] = XS*GS[1]    		#sigma-hiperon coupling constant
	GV[I] = XW*GV[1]    		#omega-hiperon coupling constant
	GR[I] = XR*GR[1]    		#rho-hiperon coupling constant

RK = RK*GS[1]**3			#constante k do termo sigma cubico do lagrangeano
RLAM = RLAM*GS[1]**4			#constante lambda do termo sigma quartico do lagrangeano

#     Initializing variables		#condicoes iniciais das variaveis

N = 5
X = np.zeros(N)

RNE0 = 0.10*RNBINF				
RNN0 = 0.90*RNBINF
      
X[0] = (3.0*PI2*RNE0)**0.16666666670
X[1] = (3.0*PI2*RNN0)**0.16666666670
X[2] = np.sin(np.sqrt(0.0050/0.9999990))     	#     X(4)=GVN/RMV*DSQRT(RNBINF)
X[3] = np.sqrt(0.22731e-5)			#     X(5)=(GRN/RMR)**2*(RI3(1)*RNN0+RI3(2)*RNE0)
X[4] = -0.57641e-6


##############################################################################
def FUNCV(N,X,FVEC):
#     ESTA SUBROTINAS É UTILIZADA PELO METODO DE BROYDEN E DA AS FUNCOES PARA SER ZERO
# 	  NO VETOR FVEC
	
	X = np.zeros(N)
	FVEC = np.zeros(N)
	
	X, RKFE, RKFN, VSIGMA, VOMEGA, VRHO = MAPPING(X)

	FSIGMA = 0.0
	FOMEGA = 0.0
	FRHO = 0.0
	CHARGE = 0.0                           #electric charge density
	RNUMBER = 0.0                           #baryonic number density 
	SCHARGEB = 0.0

	FSIGMANL = -RK/(2.0*RMS**2)*VSIGMA**2/GS[1] - RLAM/(6.0*RMS**2)*VSIGMA**3/GS[1]**2  			#Funcao em termos nao lineares sigma NL eq 5.39
      
	FOMEGANL = -CSI*GV[1]**4/(6.0*RMV**2)*VOMEGA**3/GV[1]**2  	#Funcao em termos nao lineares omega NL 5.38
	RNB = np.zeros(NH)
	for I in range(1,NH):
		RKFE, RKFN, VSIGMA, VOMEGA, VRHO, RKFB2 = SKFB2(I,RKFE,RKFN,VSIGMA,VOMEGA,VRHO)         

		if RKFB2 > 0.0:                                  
			RKFB = np.sqrt(RKFB2)     
		else:
                 	RKFB = 0.0
                

		RMBS = RMB[I] - VSIGMA*GS[I]/GS[1]           #baryon effective mass

		EFB = np.sqrt(RKFB2 + RMBS**2)     #baryon fermi energy

		RINT = 0.50*(RKFB*EFB - RMBS**2*np.log((RKFB + EFB)/RMBS))  
 
		RMA = RMBS
		RMNS = RMB[1] - VSIGMA   
		EFN = np.sqrt(RKFN**2 + RMNS**2)                    
		RMUE = np.sqrt(RKFE**2 + RML[1]**2)                  
		RMUN = VOMEGA + RI3[1]*VRHO + EFN                 
		RMU = RMUN - QB[I]*RMUE
		RNU = RMU - GV[I]/GV[1]*VOMEGA - GR[I]/GR[1]*RI3[I]*VRHO
		RE1 = Gauss_Legendre(F1,0.0,RKFB,10, RMA)      ##CALL GAUSS(F1,0.0,RKFB,10,RE1,II)
               
		FSIGMA = FSIGMA + GS[I]*GS[1]/RMS**2*RE1 

		DENSB = RKFB**3/3.0/PI2

		FOMEGA = FOMEGA + GV[I]*GV[1]/RMV**2*DENSB    		#omega equation of motion 5.38
		FRHO = FRHO + RI3[I]*GR[I]*GR[1]/RMR**2*DENSB  		#rho equation of motion 5.40
                  
		CHARGE = CHARGE + QB[I]*DENSB      
		RNUMBER = RNUMBER + DENSB
		SCHARGEB = SCHARGEB + abs(SB[I])*DENSB/3.0                  
		RNB[I] = DENSB
          

	DENSE = RKFE**3/3.0/PI2
	CHARGE = CHARGE - DENSE 
	RKFMU2 = RKFE**2 + RML[0]**2 - RML[1]**2
                                           
	if RKFMU2 > 0.0:                                  
		RKFMU = np.sqrt(RKFMU2)     
	else:
		RKFMU = 0.0
	
	DENSMU = RKFMU**3/3.0/PI2
	CHARGE -= DENSMU
	RNL[0] = DENSE
	RNL[1] = DENSMU
	FSIGMA += FSIGMANL
	FOMEGA += FOMEGANL
    
	FVEC[0] = FSIGMA - VSIGMA
	FVEC[1] = FOMEGA - VOMEGA
	FVEC[2] = FRHO - VRHO
	FVEC[3] = CHARGE
	FVEC[4] = RNUMBER - RNBT
	return FVEC
##############################################################################

#-----------------------------------------------------------------------
def SKFB2(I,RKFE,RKFN,VSIGMA,VOMEGA,VRHO):
#Esta SUB-ROTINA calcula o momento fermi ao quadrado do i-ésimo barion usando as restrições do potencial químico

	RKFE2 = RKFE**2					#quadrado de  RKFE=X(1)**2
	RKFN2 = RKFN**2					#quadrado de  RKFN=X(2)**2

	RMNS = RMB[1] - VSIGMA            		#neutron effective mass  eq 4.4    RMNS=RMB(1)-VSIGMA*GS(1)/GS(1)
	RMBS = RMB[I] - VSIGMA*GS[I]/GS[1]          	#baryon effective mass  eq 4.4

	EFN = np.sqrt(RKFN2 + RMNS**2)                     	#neutron fermi energy 3.7
	RMUE = np.sqrt(RKFE2 + RML[1]**2)                  	#electron chemical potential  3.7 mas para o lepton eletron
	RMUN = VOMEGA + RI3[1]*VRHO + EFN                 	#neutron chemical potential 5.7
	RMUB = RMUN - QB[I]*RMUE                         	#baryon chemical potential eq 5.33
	RKFB2 = (RMUB-GV[I]/GV[1]*VOMEGA - GR[I]/GR[1]*RI3[I]*VRHO)**2 - RMBS**2  #baryon fermi mom squared
	  #kfB^2 = (mu_B - g_omega*omega - g_rho*I3*rho)^2 - mB_sigma*^2   eq 5.7 isolando kfB^2

	return RKFE, RKFN, VSIGMA, VOMEGA, VRHO, RKFB2
   
#-----------------------------------------------------------------------

	     
#-----------------------------------------------------------------------
def EOS(RKFE,RKFN,VSIGMA,VOMEGA,VRHO,ENER,PRESS,RPHI,RV0):
#	Esta SUBROTINA calcula o EOS desta matéria estelar

	RPHI = VSIGMA/GS[1]
	RV0 = VOMEGA/GV[1]
	B0 = -VRHO/GR[1]

#dens de energia no momento de fermi 5.36
	ENERF = RMV**2*RV0**2/2.0 + CSI*GV[1]**4*RV0**4/8.0 + RMR**2*B0**2/2.0 + RMS**2*RPHI**2/2.0 + RK*RPHI**3/6.0 + RLAM*RPHI**4/24.0

#pressao no momento de fermi 5.37
	PRESSF = RMV**2*RV0**2/2.0 + CSI*GV[1]**4*RV0**4/24.0 + RMR**2*B0**2/2.0 - RMS**2*RPHI**2/2.0 - RK*RPHI**3/6.0 - RLAM*RPHI**4/24.0

	ENERBAR = 0.0		#energia barionica inicial
	PRESSBAR = 0.0		#pressao barionica inicial

	for I in range(1,NH):
	
		RKFE, RKFN, VSIGMA, VOMEGA, VRHO, RKFB2 = SKFB2(I,RKFE,RKFN,VSIGMA,VOMEGA,VRHO,RKFB2)          
		
		if RKFB2 > 0.0:                                  
			RKFB=DSQRT(RKFB2)     
		else:
			RKFB=0.0
  

		RMBS = RMB[I] - VSIGMA*GS[I]/GS[1]               
		RMA = RMBS
		RMNS = RMB[1] - VSIGMA   
		EFN = np.sqrt(RKFN**2 + RMNS**2)                    
		RMUE = np.sqrt(RKFE**2 + RML[1]**2)                  
		RMUN = VOMEGA + RI3[1]*VRHO + EFN                 
		RMU = RMUN - QB[I]*RMUE
		RNU = RMU - GV[I]/GV[1]*VOMEGA - GR[I]/GR[1]*RI3[I]*VRHO

		RE3 = Gauss_Legendre(F3,0.0,RKFB,10)	#CALL GAUSS(F3,0.0,RKFB,10,RE3,II)
		ENERBAR = ENERBAR + RE3
 		
		RE4 = Gauss_Legendre(F4,0.0,RKFB,10) 	#CALL GAUSS(F4,0.0,RKFB,10,RE4,II)
		PRESSBAR = PRESSBAR + RE4

	ENERLEP = 0.0
	PRESSLEP = 0.0
	RKFMU = 0.0                         
	KFMU2 = RKFE**2 + RML[1]**2 - RML[2]**2                           
              
	if RKFMU2 > 0.0:
		RKFMU = np.sqrt(RKFMU2)      
		RMUE = np.sqrt(RKFE**2 + RML[1]**2)   
		RNU = RMUE
		RMA = RML[1]
	
	RE3LE = Gauss_Legendre(F3,0.0,RKFE,10)	#CALL GAUSS(F3,0.0,RKFE,10,RE3LE,II)
	
	ENERLEP = ENERLEP + RE3LE
	
	RE4LE = Gauss_Legendre(F4,0.0,RKFE,10)	#CALL GAUSS(F4,0.0,RKFE,10,RE4LE,II)
	
	PRESSLEP = PRESSLEP + RE4LE
	RMA = RML[2]
	
	RE3LMU = Gauss_Legendre(F3,0.0,RKFMU,10)	#CALL GAUSS(F3,0.0,RKFMU,10,RE3LMU,II)
	
	ENERLEP = ENERLEP + RE3LMU
	
	RE4LMU = Gauss_Legendre(F4,0.0,RKFMU,10)	#CALL GAUSS(F4,0.0,RKFMU,10,RE4LMU,II)		
	
	PRESSLEP = PRESSLEP + RE4LMU
	ENER = ENERF + ENERBAR + ENERLEP			#soma de todos os termos de energia
	PRESS = PRESSF + PRESSBAR + PRESSLEP		#soma de todos os termos de pressao

	return RKFE, RKFN, VSIGMA, VOMEGA, VRHO, ENER, PRESS, RPHI, RV0

#-----------------------------------------------------------------------
def F1(X, RMA):
	F1 = X*X*RMA/(PI2*np.sqrt(X*X + RMA*RMA))
	return F1

#------------------------------------------------------------------
def F3(X, RMA):
	F3 = X*X*np.sqrt(X*X + RMA*RMA)*PI2		#dens energia 5.36 sem mesons

	return F3

#--------------------------------------------------------
def F4(X, RMA):

	F4 = X**4/(np.sqrt(X*X + RMA*RMA)*(3.0*PI2))		#pressao 5.37

	return F4

#---------------------------------------------------------------
#Função que resolve integrais via método de Gauss_Legendre

def Gauss_Legendre(F,a,b,N_int, RMA):		#programa integracao usando o método de Gauss-Legendre
	w = [0.6667134430869e-01,0.149451349150600,0.219086362516000,0.269266719310000,0.295524224714800,0.295524224714800,0.269266719310000,0.219086362516000,0.149451349150600,0.6667134430869e-01]
	
	t = [-0.973906528517200, -0.865063366689000,-0.679409568299000,-0.433395394129200,-0.148874338981600,+0.148874338981600,+0.433395394129200,+0.679409568299000,+0.865063366689000,+0.973906528517200]

	Int = 0.0			#Valor inicial a ser integrado
	h = (b - a)/2
	
	for I in range(0,10):
		x = 0.50*(t[I]*(b - a) + (b + a))
		Int = Int + h*w[I]*F(x, RMA)
		
	return Int

#------------------------------------------
#Função que resolve sistema não-linear via método de Broyden

def gauss_jordan(a,b): 
    n = len(b) 
    for k in range(0,n-1): 
        for i in range(k+1,n): 
            if a[i][k] != 0.0: 
                m = a [i][k]/a[k][k] 
                a[i][k+1:n] = a[i][k+1:n] - m*a[k][k+1:n] 
                b[i] = b[i] - m*b[k] 
    for k in range(n-1,-1,-1): 
        b[k] = (b[k] - np.dot(a[k][k+1:n],b[k+1:n]))/a[k][k] 
    return b

#função que resolve o método de Broyden para sist não-lineares	
def Broyden(F, B, x, N, eps, k_max):
	k = 0
	F_v = F(x)
	B = B(x)
	F_norma = np.linalg.norm(F_v, ord=2)
	N = len(F_v)

	while abs(F_norma) > eps and k < k_max:
		
		s = gauss_jordan(B, -F_v)
        
		x = x + s
		F_novo = F(x)
		
		dF = F_novo - F_v
        
		B = B + (np.outer((dF - np.dot(B,s)),s))/(np.dot(s,s))
        
		F_norma = np.linalg.norm(F_v, ord=2)
		
		F_v = F_novo
        
		k = k + 1

	return x

#testar matriz identidade
def B(x):
	return np.array([[1.0, 0.0, 0.0, 0.0, 0.0],[2.0, 1.0, 0.0, 0.1, 0.1],[0.0,0.1,1.0,0.1,0.1],[0.0,0.1,1.0,0.1,0.1],[0.0,0.1,1.0,0.1,0.1]])


#--------------------------------------------------------------------------
def MAPPING(X):
#Esta SUB-ROTINA define os campos potenciais entre os valores físicos   

	RKFE = X[0]**2                            	#to assure kfe greater than zero
	RKFN = X[1]**2                            	#to assure kfn greater than zero
	VSIGMA = 0.9999990*np.sin(X[2])**2   		#mapping vsigma in (0,MB(1))
	VOMEGA = X[3]**2                       		#to assure vomega greater than zero
	VRHO = X[4]                              	#no mapping on vrho
      
	return X, RKFE, RKFN, VSIGMA, VOMEGA, VRHO


#Main loop 	#loop principal que calcula as equações de estado

print("Inside Main Loop")
for I in tqdm.tqdm(range(1,NPOINT)):
	RNBT = RNBINF + (I - 1)*DNB
    
	eps = 0.0000001
	k_max = 1000
	FVEC = np.zeros(N)
	FVEC = FUNCV(N,X,FVEC)						#função eqs de movimento das não-lineares
	X = Broyden(FVEC, B, x, N, eps, k_max)		#função resolve as eqs n lineares


	X, RKFE, RKFN, VSIGMA, VOMEGA, VRHO = MAPPING(X)

	RKFE, RKFN, VSIGMA, VOMEGA, VRHO, ENER, PRESS, RPHI, RV0 = EOS(RKFE,RKFN,VSIGMA,VOMEGA,VRHO,ENER,PRESS,RV0,RPHI)    #the EOS
 
	RNBTD = RNBT*(RM/197.320)**3
	RMNS = (RMB[1] - VSIGMA)

	RMUN = VOMEGA + RI3[1]*VRHO + np.sqrt((3.0*PI2*RNB[1])**0.66666666670 + RMNS**2)   	#neutron chemical potential
	RMUE = DSQRT((3.0*PI2*RNL[1])**0.66666666670 + RML[1]**2) 				#electron chemical potential
         
	DENER = ENER/RNBT - 1.0   		#dens de energia de ligação - massa (valor 1)
	DENER = DENER*RM			#multiplicado por RM para ficar em MeV


#     OUTPUT IN fm^-4 FOR THE PROGRAM THAT CALCULATES STARS PROPERTIES; 
#     IF YOU CHOOSE THIS OUTPUT, COMMENT THE ABOVE ONE
	PRESS = PRESS*(RM/197.32)**4			#pressao em fm^-4
	ENER = ENER*(RM/197.32)**4			#energia em fm^-4


#	file.write("%.10f %.10f %.10f \n" %(RNBTD,ENER,PRESS,RMUN*RM))
#	WRITE(22,10) RNBTD, (RNL(j)/RNBT,j=1,2), (RNB(j)/RNBT,j=1,NH) xxxxxxxxxxx acertar xxxxxxxxxxxxx

	file.write("%.10f %.10f %.10f \n" %(RNBTD,ENER,PRESS,RMUN*RM))	
#	WRITE(100,*) RNBTD,ENER,PRESS,RMUN*RM xxxxxxxxxxxxxxx acertar xxxxxxxxxxxxxxx
  
	RMPU = RMB(1)*RNB(1) + RMB(2)*RNB(2) + RMB(3)*RNB(3) + RMB(4)*RNB(4)+ RMB(5)*RNB(5) + RMB(6)*RNB(6) + RMB(7)*RNB(7) + RMB(8)*RNB(8)
             
	RMPUM = ((RMPU/(RNBTD)))*(939.0/197.330)**3
	RMPUMMEV = RMPUM*939.0                                #massa média por partícula
	RMPUMFM =  RMPUMMEV/197.330
     

print("Total execution time: ", start-time.time()," seconds")