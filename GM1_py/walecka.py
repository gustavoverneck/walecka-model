# Libraries
import numpy as np
import scipy.linalg as sla
import tqdm
import time

class BaryonOctetModel:
    def __init__(self):
        print("<----- Walecka Model for Stellar Matter ----->")
        self.modelChoice()
        self.getInput()
        self.main()

    def modelChoice(self):
        # -------- Option of Model --------- #
        chose = False
        print("\nChoose your model:")
        while not chose:
            self.nh = int(input("Insert:\n\t2 - for Nucleons\n\t8 - for Nucleons + Hyperons\n\nYour choice: "))
            nh = self.nh
            if nh == 2 or nh == 8:
                chose = True
            else:
                print("Invalid choice. Please try again.")
        if nh == 2:
            print("You chose: Nucleon\n")
        elif nh == 8:
            print("You choice: Nucleon + Hyperons\n")

    def getInput(self):
        # -------- GET INPUT DATA FROM FILE-------- #
        with open("/home/gverneck/codes/tcc/GM1_py/input.txt", 'r') as f:
            lines = f.readlines()
            self.rm = float(lines[0].split("\t")[0]) # rm
            self.ms = float(lines[1].split("\t")[0]) / self.rm # ms
            self.mv = float(lines[2].split("\t")[0]) / self.rm # mv
            self.mr = float(lines[3].split("\t")[0]) / self.rm # mr
            self.xls = float(lines[4].split("\t")[0])
            self.xlw = float(lines[5].split("\t")[0])
            self.xlr = float(lines[6].split("\t")[0])
            self.k = float(lines[7].split("\t")[0])
            self.lam = float(lines[8].split("\t")[0])
            self.csi = float(lines[9].split("\t")[0])
            self.xs = float(lines[10].split("\t")[0])
            self.xw = float(lines[11].split("\t")[0])
            self.xr = float(lines[12].split("\t")[0])
            self.nbnif = float(lines[13].split("\t")[0]) / (self.rm/197.32)**3
            self.nbsup = float(lines[14].split("\t")[0]) / (self.rm/197.32)**3
            self.n_points = int(lines[15].split("\t")[0])

            self.dnb = (self.nbsup-self.nbnif)/(self.n_points-1)

        # -------- MANUAL INPUT DATA -------- #
            self.mb = np.array([939.0, 939.0, 1116.0, 1193.0, 1193.0, 1193.0, 1318.0, 1318.0]) #definindo vetor de 8 componentes
            self.ml = np.array([0.5110, 105.660])
            self.qb = np.array([0.00,1.00,0.00,-1.00,0.00,1.00,-1.00,0.00])
            self.i3 = np.array([-0.50,0.50,0.00,-1.00,0.00,1.00,-0.50,0.50])
            self.sb = np.array([0.00,0.00,-1.0,-1.0,-1.0,-1.0,-2.0,-2.0])
            self.gs = np.zeros(self.nh)
            self.gv = np.zeros(self.nh)
            self.gr = np.zeros(self.nh)

        for i in range(2,8):
            self.mb[i] /= self.rm # Adimensionalidade barions/hiperions
        
        for i in range(0,2):
            #self.mb = 1.0  # DUVIDA - POR QUE TORNAR 1.0? O for anterior a define..
            self.ml[i] /= self.rm

        self.xls = self.xls*(self.rm/197.330)**2
        self.xlw = self.xlw*(self.rm/197.330)**2
        self.xlr = self.xlr*(self.rm/197.330)**2

        for i in range(0,2):
            self.gs[i] = np.sqrt(self.xls)*self.ms  # Constante de acoplamento
            self.gv[i] = np.sqrt(self.xlw)*self.mv
            self.gr[i] = np.sqrt(self.xlr)*self.mr

        for i in range(2,8):    # DUVIDA - é gs[0] mesmo?
            self.gs[i] = self.gs[0] * self.xs   # sigma-hiperion Coup. Const.
            self.gv[i] = self.gs[0] * self.xw   # omega-hiperion coup. const.
            self.gv[i] = self.gs[0] * self.xr   # rho-hiperion coup. const
        
        self.k *= self.gs[1]**3    # k do termo sigma cubico do Lagrang.
        self.lam *= self.gs[1]**4  # lambda do termo sigma quadratico do Lagrang.

        # Initializing Variables
        self.n = 5
        self.x = np.zeros(self.n)

        self.ne0 = 0.10 * self.nbnif
        self.nn0 = 0.90 * self.nbnif

        self.x[0] = (3.0*np.pi**2*self.ne0)**0.16666666670
        self.x[1] = (3.0*np.pi**2*self.nn0)**0.16666666670
        self.x[2] = np.sin(np.sqrt(0.0050/0.9999990))
        self.x[3] = np.sqrt(0.22731e-5)
        self.x[4] = -0.57641e-6

        self.eps = 10e-8
        self.k_max = 1000

    def main(self):
        file = open("eos.dat", "w+")
        for i in tqdm.tqdm(range(0, self.n_points)):
            self.nbt = self.nbnif + (i+1)*self.dnb
            x = self.Broyden(self.x, self.fvec())
            kfe, kfn, vsigma, vomega, vrho = self.mapping()
            
            ener, press, presslep, enerlep = self.eos(kfe, kfn, vsigma, vomega, vrho)

            self.nbtd = self.nbt*(self.rm/197.320)**3
            self.mns = self.mb[0] - vsigma
            mum = vomega + self.i3[0]*vrho + np.sqrt((3.0*np.pi**2*self.nb[0])**0.66666666670 + self.mns**2)   	#neutron chemical potential
            mue = np.sqrt((3.0*np.pi**2*self.nl[0])**0.66666666670 + self.ml[0]**2)
            dener = ener/self.nbt - 1.0     # Densidade de energia
            dener *= self.rm                # Mult. por rm para ficar em MeV
            file.write("%.10f %.10f %.10f %.10f\n" %(self.nbtd,ener,press,mum*self.rm))
            rmpu = self.mb[0]*self.nb[0] + self.mb[1]*self.nb[1] + self.mb[2]*self.nb[2] + self.mb[3]*self.nb[3]+ self.mb[4]*self.nb[4] + self.mb[5]*self.nb[5] + self.mb[6]*self.nb[6] + self.mb[7]*self.nb[7]
            rmpum = ((rmpu/(self.nbtd)))*(939.0/197.330)**3
            rmpummev = rmpum*939.0
            rmpumfm =  rmpummev/197.330
            #DUVIDA - O CÓDIGO GERA OUTPUT E CALCULA ESSAS ÚLTIMAS COISAS
        
    def fvec(self):     # returns fv - Duvida, reseta as variaveis
        fv = np.zeros(self.n)
        x = np.zeros(self.n)
        self.nb = np.zeros(self.nh)
        self.nl = np.zeros(2)
        kfe, kfn, vsigma, vomega, vrho = self.mapping()

        fsigma = 0.001
        fomega = 0.001
        frho = 0.0
        charge = 0.0
        rnumber = 0.0
        scahrge = 0.0

        fsigmanl = -self.k/(2.0*self.ms**2)*vsigma**2/self.gs[0] - self.lam/(6.0*self.ms**2)*vsigma**3/self.gs[0]**2
        fomeganl = -self.csi*self.gv[0]**4/(6.0*self.mv**2)*vomega**3/self.gv[0]**2
        
        for i in range(1, self.nh):
            self.kfb2 = self.skfb2(i, kfe, kfn, vsigma, vomega, vrho)
            kfb2 = self.kfb2
            if kfb2 > 0.0:
                kfb = np.sqrt(kfb2)
            else:
                kfb = 0.0
        
            mbs = self.mb[i] - vsigma*self.gs[i]/self.gs[0]
            efb = np.sqrt(kfb2 + mbs**2)
            #rint=0.50*(kfb*efb - mbs**2*np.log((kfb + efb)/mbs))
            rma = mbs
            mns = self.mb[0] - vsigma
            efn = np.sqrt(kfn**2 + mns**2)
            mue = np.sqrt(kfe**2 +self.ml[0]**2)
            mum = vomega + self.i3[0]*vrho + efn
            #rmu = mum - qb[i]*mue
            #rnu = rmu - self.gv[i]/self.gv[0]*vomega - self.gr[i]/self.gr[0]*i3[i]*vrho

            re1 = self.GaussLegendre(self.f1(kfb, rma), 0.0, kfb, 10)
            fsigma += self.gs[i]*self.gs[0]/self.ms**2*re1
            densb = kfb**3/3.0/np.pi**2
            fomega += self.gv[i]*self.gv[0]/self.mv**2*densb
            charge += self.qb[i]*densb
            rnumber += densb
            self.nb[i] = densb
        
        dense = kfe**3/3.0/np.pi**2
        charge -= dense
        kfmu2 = kfe**2 + self.ml[0]**2 - self.ml[1]**2

        if kfmu2 > 0.0:
            kfmu = np.sqrt(kfmu2)
        else:
            kfmu = 0.0
        
        densmu = kfmu**3/3.0/np.pi**2
        charge == densmu
        self.nl[0] = dense
        self.nl[1] = densmu

        fsigma += fsigmanl
        fomega += fomeganl

        fv[0] = fsigma + fsigmanl
        fv[1] = fomega - vomega
        fv[2] = frho - vrho
        fv[3] = charge
        fv[4] = rnumber - self.nbt

        return np.array([fv[0], fv[1], fv[2], fv[3], fv[4]])
    
    def mapping(self):  # returns kfe, kfn, vsigma, vomega, vrho
        kfe = self.x[0]**2
        kfn = self.x[1]**2
        vsigma = 0.9999990*np.sin(self.x[2])**2 # DUVIDA - pq multiplicar por 0.9...
        vomega = self.x[3]**2
        vrho = self.x[4]

        return kfe, kfn, vsigma, vomega, vrho

    def skfb2(self, i, kfe, kfn, vsigma, vomega, vrho): # returns kfb2
        kfe2 = kfe**2
        kfn2 = kfn**2
        mns = self.mb[1] - vsigma
        mbs = self.mb[i] - vsigma*self.gs[i]/self.gs[1]
        efn = np.sqrt(kfn2 + mns**2)
        mue = np.sqrt(kfe2 + mns**2)
        mun = vomega + self.i3[i]*vrho + efn
        mub = mun - self.qb[i]*mue
        kfb2 = (mub - self.gv[i]/self.gv[1]*vomega - self.gv[i]/self.gv[1]*self.i3[i]*vrho)**2 - mbs**2
        return kfb2

    def GaussLegendre(self, f, a, b, n_int):
        w = [0.6667134430869e-01,0.149451349150600,0.219086362516000,0.269266719310000,0.295524224714800,0.295524224714800,0.269266719310000,0.219086362516000,0.149451349150600,0.6667134430869e-01]
        t = [-0.973906528517200, -0.865063366689000,-0.679409568299000,-0.433395394129200,-0.148874338981600,+0.148874338981600,+0.433395394129200,+0.679409568299000,+0.865063366689000,+0.973906528517200]
        it = 0.0
        h = (b-a)/2
        for i in range(0,10):
            x = 0.50*(t[i]*(b - a) + (b + a))
            it += h*w[i]*f
        return it

    def GaussJordan(self, a, b):
        n = len(b)
        for k in range(0, n-1):
            for i in range(k+1, n):
                if a[i][k] != 0.0:
                    m = a[i][k]/a[k][k]
                    for j in range(k+1, n):
                        a[i][j] = a[i][j] - m*a[k][j]
                    b[i] = b[i] - m*b[k]
        for k in range(n-1, -1, -1):
            b[k] = (b[k] - np.dot(a[k][k:n], b[k:n]))/a[k][k]
        return b

    def f1(self, kr, rma):
        f1 = kr*kr*rma/(np.pi**2*np.sqrt(kr*kr + rma*rma))
        return f1

    def f3(self, kr, rma):
        f3 = kr*kr*np.sqrt(kr*kr + rma*rma)*np.pi**2
        return f3

    def f4(self, kr, rma):
        f4 = kr**4 / (np.sqrt(kr*kr + rma*rma)*(3.0*np.pi**2))
        return f4

    def Broyden(self, xk, f):
        k = 0
        n = len(xk)
        fk = f
        B = self.B()
        f_norma = np.linalg.norm(fk, ord=2)

        while abs(f_norma) > self.eps and k < self.k_max:
            s = self.GaussJordan(B, -1*fk)
            xk += s
            f_novo = f
            df = f_novo - fk
            B = B + (np.outer((df - np.dot(B, s)), s)) / (np.dot(s, s)) # ERRO NA DIVISÃO - GERANDO INFINITOS

            f_norma = np.linalg.norm(fk, ord=2)
            fk = f_novo
            k += 1
        
        return xk

    def B(self):
        return np.array([[1.0,0.0,0.0, 0.0,0.0],[0.0,1.0,0.0, 0.0,0.0],[0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,1.0]])

    def eos(self, kfe, kfn, vsigma, vomega, vrho):  # returns ener, press, presslep, enerlep
        phi = vsigma / self.gs[0]
        v0 = vomega / self.gv[0]
        b0 = -vrho / self.gr[0]

        # Densidade de Energia no Momento de Fermi - 5.36
        enerf = self.mv**2*v0**2/2.0 + self.csi*self.gv[0]**4*v0**4/24.0 + self.mr**2*b0**2/2.0 - self.ms**2*phi**2/2.0 - self.k*phi**3/6.0 - self.lam*phi**4/24.0
        pressf = self.mv**2*v0**2/2.0 + self.csi*self.gv[0]**4*v0**4/24.0 + self.mr**2*b0**2/2.0 - self.ms**2*phi**2/2.0 - self.k*phi**3/6.0 - self.lam*phi**4/24.0
        enerbar = 0.01      # energia barionica inicial
        pressbar = 0.01     # pressão barionica inicial

        for i in range(0, self.nh):
            self.kfb2 = self.skfb2(i, kfe, kfn, vsigma, vomega, vrho)
            if self.kfb2 > 0.0:
                kfb = np.sqrt(kfb2)
            else:
                kfb = 0.0

            mbs = self.mb[i] - vsigma*self.gs[i]/self.gs[0]
            rma = mbs
            mns = self.mb[0] - vsigma
            efn = np.sqrt(kfn**2 + mns**2)
            mue = np.sqrt(kfe**2 + self.ml[0]**2)
            mum = vomega + self.i3[0]*vrho + efn
            rmu = mum + self.qb[i]*mue
            rnu = rmu - self.gv[i]/self.gv[0]*vomega - self.gr[i]/self.gr[0]*self.i3[i]*vrho

            re3 = self.GaussLegendre(self.f3(kfb, rma), 0.0, kfb, 10)
            enerbar += re3

            re4 = self.GaussLegendre(self.f4(kfb, rma), 0.0, kfb, 10)
            pressbar += re4
        
        enerlep = 0.001
        presslep = 0.001
        kfmu = 0.001
        kfmu2 = kfe**2 + self.ml[0]**2 - self.ml[1]**2

        if kfmu2 > 0.0:
            kfmu = np.sqrt(kfmu2)
            mue = np.sqrt(kfe**2 + self.ml[0]**2)
            rnu = mue
            rma = self.ml[0]

        re3le = self.GaussLegendre(self.f3(kfb, rma), 0.0, kfe, 10)
        enerlep += re3le

        re4le = self.GaussLegendre(self.f4(kfb, rma), 0.0, kfe, 10)
        presslep += re4le
        rma = self.ml[1]

        re3lmu = self.GaussLegendre(self.f3(kfb, rma), 0.0, kfmu, 10)
        enerlep += re3lmu

        re4lmu = self.GaussLegendre(self.f4(kfb, rma), 0.0, kfmu, 10)

        presslep += re4lmu
        ener = enerf + enerbar + enerlep        # Soma de todos os termos de energia
        press = pressf + pressbar + presslep    # Soma de todos os termos de pressão

        return ener, press, presslep, enerlep


# <------------------------------------------>

# Just run the code in a exquisite manner
if __name__ == "__main__":
    start = time.time()
    BaryonOctetModel()
    print("Total execution time: ", round(time.time()-start, 3)," seconds")