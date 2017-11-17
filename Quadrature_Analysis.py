################################################################################
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import Gauss_Legendre as Leg
################################################################################

ASCFILE = "relativeError.asc"

def Repeated_Quadrature(XMLFILE,n1,n2):
    error = [0]*(abs(n2-n1)+1)
    n = [i for i in range(n1,n2+1)]
    with open(ASCFILE,'w') as errorFile:
        for i in range(n1,n2+1):
            [numeric, analytic] = Leg.Return_Quadrature(XMLFILE,i)
            error[i-n1] = abs(numeric-analytic)/analytic
            errorFile.write(str(error[i-n1])+"\t"+str(i)+"\n")
    
def Convergence_Graph(ASCFILE,n1,n2):
    with open(ASCFILE,'r') as plotFile:
        error = [0]*abs(n2-n1)
        n = [0]*abs(n2-n1)     
        for i in range(0,abs(n2-n1)):
            temp = plotFile.readline().split()
            error[i]    = float(temp[0])
            n[i]        = int(temp[1])
    plt.plot(n,error)
    plt.show()
    
def main():
    os.chdir(sys.path[0])
    XMLFILE = sys.argv[1]
    PROGRAM = int(sys.argv[2])
    n1      = int(sys.argv[3])
    n2      = int(sys.argv[4])
    if PROGRAM == 1:
        Repeated_Quadrature(XMLFILE,n1,n2)
    if PROGRAM == 2:
        Convergence_Graph(ASCFILE,n1,n2)
        
if __name__ == '__main__':
    main()