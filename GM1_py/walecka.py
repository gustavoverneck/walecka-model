# Reescrita do código de F90 para Python

# --- Lista de Funções a serem implementadas ---
#   FUNCV(N, X, FVEC)
#   SKFB2(I, RKFE, RKFN, VSIGMA, VOMEGA, VRHO)
#   EOS(RKFE, RKFN, VSIGMA, VOMEGA, VRHO, ENER, PRESS, RPHI, RV0)
#   F1, F3, F4
#   GAUSS-LEGENDRE
#   GAUSS-JORDAN
#   BROYDEN
#   B
#   MAPPING
#   


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
            nh = int(input("Insert:\n\t2 - for Nucleons\n\t8 - for Nucleons + Hyperons\n\nYour choice: "))
            if nh == 2 or nh == 8:
                chose = True
            else:
                print("Invalid choice. Please try again.")
        if nh == 2:
            print("You chose: Nucleon")
        elif nh == 8:
            print("You choice: Nucleon + Hyperons")
        self.NH = nh

    def getInput(self):
        # -------- GET INPUT DATA FROM FILE-------- #
        with open("/home/gverneck/codes/tcc/GM1_py/input.txt", 'r') as f:
            lines = f.readlines()
            self.RM = float(lines[0].split("\t")[0]) # RM
            self.RMS = float(lines[1].split("\t")[0]) / self.RM # RMS
            self.RMV = float(lines[2].split("\t")[0]) / self.RM # RMV
            self.RMR = float(lines[3].split("\t")[0]) / self.RM # RMR
            self.XLS = float(lines[4].split("\t")[0])
            self.XLW = float(lines[5].split("\t")[0])
            self.RK = float(lines[6].split("\t")[0])
            self. RLAM = float(lines[7].split("\t")[0])
            self.CSI = float(lines[8].split("\t")[0])
            self.XS = float(lines[9].split("\t")[0])
            self. XW = float(lines[10].split("\t")[0])
            self.XR = float(lines[11].split("\t")[0])
            self.RNBINF = float(lines[12].split("\t")[0]) / (self.RM/197.32)**3
            self.RNBSUP = float(lines[13].split("\t")[0]) / (self.RM/197.32)**3
            self.n_points = float(lines[14].split("\t")[0])
            self.DNB = (self.RNBSUP-self.RNBINF)/(self.n_points-1)

        # -------- MANUAL INPUT DATA -------- #
            self.RMB = np.array([939.0, 939.0, 1116.0, 1193.0, 1193.0, 1193.0, 1318.0, 1318.0]) #definindo vetor de 8 componentes
            self.QB = np.array([0.00,1.00,0.00,-1.00,0.00,1.00,-1.00,0.00])
            self.RI3 = np.array([-0.50,0.50,0.00,-1.00,0.00,1.00,-0.50,0.50])
            self.SB = np.array([0.00,0.00,-1.0,-1.0,-1.0,-1.0,-2.0,-2.0])
            self.GS = np.zeros(self.NH)
            self.GV = np.zeros(self.NH)
            self.GR = np.zeros(self.NH)

        for i in range(3,8):
            self.RMB[i] /= self.RM # Adimensionalidade barions/hiperions
        
        for i in range(1,2):
            #self.RMB = 1.0  # DUVIDA - POR QUE TORNAR 1.0?
            self.RMB[i] /= self.RM

        self.XLS = self.XLS*(self.RM/197.330)**2
        self.XLW = self.XLW*(self.RM/197.330)**2
        self.XLR = self.XLR*(self.RM/197.330)**2

        for i in range(1,2):
            self.GS[i] = np.sqrt(self.XLS)*self.RMS  # Constante de acoplamento
            self.GV[i] = np.sqrt(self.XLW)*self.RMV
            self.GR[i] = np.sqrt(self.XLR)*self.RMR

        for i in range(3,8):
            self.GS[i] *= self.XS   # sigma-hiperion Coup. Const.
            self.GV[i] *= self.XW   # omega-hiperion coup. const.
            self.GR[i] *= self.XR   # rho-hiperion coup. const
        
        self.RK *= self.GS[1]**3    # k do termo sigma cubico do Lagrang.
        self.RLAM *= self.GS[1]**4  # lambda do termo sigma quadratico do Lagrang.

        # Initializing Variables
        self.N = 5
        self.X = np.zeros(N)

        self.RNE0 = 0.10 * self.RNBINF
        self.RNN0 = 0.90 * self.RNBINF

        X[0] = (3.0*np.pi**2*self.RNE0)**0.16666666670
        X[1] = (3.0*np.pi**2*self.RNN0)**0.16666666670
        X[2] = np.sin(np.sqrt(0.0050/0.9999990))
        X[3] = np.sqrt(0.22731e-5)
        X[4] = -0.57641e-6

    def main(self):
        for i in tqdm.tqdm(range()):
            self.RNBT = self.RNBINF + (i+1)*self.DNB

            eps = 10e-8
            k_max = 1000
            self.FVEC = np.zeros(self.N)
            self.FVEC = ...



# <------------------------------------------>

# Just run the code in a exquisite manner
if __name__ == "__main__":
    start = time.time()
    BaryonOctetModel()
    print("Total execution time: ", start-time.time()," seconds")