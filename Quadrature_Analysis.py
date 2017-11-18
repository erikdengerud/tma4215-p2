################################################################################
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import Gauss_Legendre as Leg
################################################################################


##Repated Gauss-Legendre quadrature for given nodes {n1, n1+1,...,n2-1,n2}. Save everything to an .asc-file
def Repeated_Quadrature(ASCFILE,XMLFILE,n1,n2):
    with open('ASCfiles/'+ASCFILE,'w') as errorFile:
        error = [0]*(abs(n2-n1)+1)
        n = [i for i in range(n1,n2+1)]
        for i in range(n1,n2+1):
            [numeric, analytic] = Leg.Return_Quadrature('XMLfiles/' + XMLFILE,i)
            error[i-n1] = abs((numeric-analytic)/analytic)
            errorFile.write(str(error[i-n1])+"\t"+str(i)+"\n")
        print("File successfully written.")
        
##Plot the convergence for the performed quadrature. Save plot to .png-file    
def Convergence_Graph(ASCFILE,n1,n2):
    with open('ASCfiles/'+ASCFILE,'r') as plotFile:
        error = [0]*(abs(n2-n1)+1)
        n = [0]*(abs(n2-n1)+1)     
        for i in range(0,abs(n2-n1)+1):
            temp = plotFile.readline().split()
            error[i]    = float(temp[0])
            n[i]        = int(temp[1])
        print("File successfully read. Close plot window to exit.")
        
    plt.plot(n,error)
    plt.xlabel(r'# of integration nodes $n$')
    plt.ylabel(r'Absolute error $e(n)$')
    plt.xticks(np.arange(n1, n2+1, 1.0))
    pngName = ASCFILE.split('.')
    pngName = pngName[0]+"Figure.png"
    plt.savefig('ErrorPlots/'+pngName)
    
    #plt.show()
   
    
def main():
    os.chdir(sys.path[0])
    XMLFILE = sys.argv[1]
    PROGRAM = int(sys.argv[2])
    n1      = int(sys.argv[3])
    n2      = int(sys.argv[4])
    
    ASCFILE = XMLFILE.split('.')
    ASCFILE = ASCFILE[0]+"RelativeError"+str(n1)+"to"+str(n2)+".asc"
     
    if PROGRAM == 1:
        Repeated_Quadrature(ASCFILE,XMLFILE,n1,n2)
    if PROGRAM == 2:
        Convergence_Graph(ASCFILE,n1,n2)
        
if __name__ == '__main__':
    main()