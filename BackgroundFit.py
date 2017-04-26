##This code implements the reaction plane method proposed by C. Nattrass et al. in Phys. Rev. C 93, 044915.
##It takes as input histograms in a ROOT file. It returns the result of the fit in a MIGRAD object 

## Version: v0, in progress
## Author: Miguel Arratia. 
import ROOT
import numpy as np
from iminuit import Minuit, describe, Struct
from math import cos
from math import sin
from math import pi
from ROOT import gROOT
#
from FitModel import Signal, Background, TotalPDF
from PlotFitResult import Plot
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors as mcolors

dataALL = {}
dataBKG = {}
phi_s = {}
c  = {}

def FromTH1toArray(histo):
  # Calculate how many y bins we have
  xaxis = histo.GetXaxis()
  nx = xaxis.GetNbins()
  xbins = range(1, nx+1) # Removes under/overflow
  d = {}
  d["x_center"] = np.array([xaxis.GetBinCenter(i) for i in xbins])
  d["x_low"]    = np.array([xaxis.GetBinLowEdge(i) for i in xbins])
  d["x_up"]     = np.array([xaxis.GetBinUpEdge(i) for i in xbins])
  d["y"]        = np.array([histo.GetBinContent(i) for i in xbins])
  d['dy']       = np.array([histo.GetBinError(i) for i in xbins])
  return d


   

def Chi2(A3,A4,s3,s4,C2,
         A5,A6,s5,s6,C3,
         A7,A8,s7,s8,C4,
         B, v2_t, v2_a, v4_t, v4_a, V1, V3):
    params_Signal = {}
    params_BKG = {}

    params_Signal["Midplane"]  = {}
    params_Signal["Inplane"]  = {}
    params_Signal["Outplane"] = {}
 
    params_Signal["Outplane"]["A1"] = A3
    params_Signal["Outplane"]["A2"] = A4
    params_Signal["Outplane"]["s1"] = s3
    params_Signal["Outplane"]["s2"] = s4
    params_Signal["Outplane"]["C1"] = C2
    
    params_Signal["Midplane"]["A1"] = A5
    params_Signal["Midplane"]["A2"] = A6
    params_Signal["Midplane"]["s1"] = s5
    params_Signal["Midplane"]["s2"] = s6
    params_Signal["Midplane"]["C1"] = C3

    params_Signal["Inplane"]["A1"] = A7
    params_Signal["Inplane"]["A2"] = A8
    params_Signal["Inplane"]["s1"] = s7
    params_Signal["Inplane"]["s2"] = s8
    params_Signal["Inplane"]["C1"] = C4

    
    params_BKG["B"] = B
    params_BKG["v2_t"] = v2_t
    params_BKG["v2_a"] = v2_a
    params_BKG["v4_t"] = v4_t
    params_BKG["v4_a"] = v4_a
    params_BKG["V1"] = V1
    params_BKG["V3"] = V3

    return Chi2_All(params_Signal, params_BKG) + Chi2_BKG(params_BKG)

def Chi2_All(params_Signal, params_BKG):
    #### FIT TO SIGNAL + BACKGROUND
    total_chi2 = 0
    for key in dataBKG.keys():
        x = dataALL[key]["x_center"]
        y = dataALL[key]["y"]
        sigma = dataALL[key]["dy"]
        #print x, y
        total_chi2 = total_chi2 + np.sum( np.power( y - Signal( x , params_Signal[key] ) - Background( x , params_BKG, phi_s[key], c[key] ) , 2.0)/ np.power(sigma,2.0) )
    return total_chi2

def Chi2_BKG(params_BKG):
    total_chi2 = 0
    ### BACKGROUND ONLY
    for key in dataBKG.keys():
        x = dataBKG[key]["x_center"]
        y = dataBKG[key]["y"]
        sigma = dataBKG[key]["dy"]
        length = len(x)
        x = x[:int(length/2)]
        y = y[:int(length/2)] #this only considers half the histogram (i.e, "Dphi < pi/2")
        sigma = sigma[:int(length/2)]
        #print x, y
        total_chi2 = total_chi2+ np.sum( np.power( y - Background( x , params_BKG, phi_s[key], c[key] )  , 2.0)/ np.power(sigma,2.0))
    return total_chi2

def PerformFitTotal():
    print ' About to start MINUIT'
    sigma_init = 0.07
    sigma_init_away = 0.2
    limit_sigmaup = 0.35
    limit_sigmado = 0.025
    m = Minuit(Chi2,
               s3=sigma_init, limit_s3=(limit_sigmado,limit_sigmaup), error_s3=0.001,
               s4=sigma_init_away, limit_s4=(limit_sigmado,limit_sigmaup), error_s4=0.001,
               s5=sigma_init, limit_s5=(limit_sigmado, limit_sigmaup), error_s5=0.001, 
               s6=sigma_init_away, limit_s6=(limit_sigmado, limit_sigmaup), error_s6=0.001,
               s7=sigma_init, limit_s7=(limit_sigmado, limit_sigmaup), error_s7=0.001, 
               s8=sigma_init_away,  limit_s8=(limit_sigmado, limit_sigmaup), error_s8=0.001,
               C2=0.0, fix_C2=True, 
               C3=0.0, fix_C3=True,
               C4=0.0, fix_C4=True,
               #A3 =0.55, limit_A3=(0.1,1.0), error_A3=0.01, 
               #A4=0.05, limit_A4=(0.01, 1.0), error_A4=0.0001,
               #A5 =0.55, limit_A5=(0.1,1.0), error_A5=0.01, 
               #A6=0.05, limit_A6=(0.01, 1.0), error_A6=0.0010,
               #A7 =0.55, limit_A7=(0.1,1.0), error_A7=0.01, 
               #A8=0.05, limit_A8=(0.01, 1.0), error_A8=0.0001,
               v2_t=0.02, limit_v2_t =(0,0.50), error_v2_t=0.001, 
               v2_a=0.02, limit_v2_a =(0,0.50), error_v2_a=0.001,
               v4_t=0.01, limit_v4_t =(0,0.50), error_v4_t=0.001, 
               v4_a=0.01, limit_v4_a =(0,0.50), error_v4_a=0.001,
               V3=0 , limit_V3 = (-0.1, 0.5), error_V3 =0.001,
               V1=0.0, fix_V1=True) 
               #V1=0.0, limit_V1 = (-0.1, 0.5), error_V1=0.0001)
               #B=2.0 , limit_B = (0.1, 10.0), error_B = 0.01)
    m.migrad()
    #m.minos()  #more sophisticated error estimation
    m.print_matrix() #correlation
    plt.show()
    return m
    

filename = "fInput"


finput = ROOT.TFile("%s.root" %(filename),"READ")
llaves_plane = {"Inplane", "Midplane", "Outplane"}

##Read data from root files
for key in llaves_plane:
    dataALL[key] = FromTH1toArray( finput.Get("SignalRegion_%s" %(key)))
    dataBKG[key] = FromTH1toArray( finput.Get("BackgroundRegion_%s" %(key)))
    #dataBKG[key] = FromTH1toArray( finput.Get("BKG_BackgroundRegion_%s" %(key)))
## Define globals
#######################################
phi_s["Inplane"]  = 0
phi_s["Midplane"] =pi/4.0
phi_s["Outplane"] =pi/2.0

c["Inplane"]  = pi/6.0
c["Midplane"] = pi/12.0
c["Outplane"] = pi/6.0
########################################

#Minuit = PerformFitTotal()
Minuit = PerformFitTotal()
#Minuit.draw_contour('B','V3', bound=5, show_sigma=True);

def PlotMinosContour(var1,var2,SigmaMax=3, npoints=10):
    for Sigma in range(1,SigmaMax+1):
        errx, erry, line = Minuit.mncontour(var1, var2, numpoints=npoints, sigma=Sigma)
        #print line
        x = [pair[0] for pair in line]
        y = [pair[1] for pair in line]
        plt.plot(x,y,'-o',label='%i sigma' %int(Sigma))
    plt.xlabel(var1)
    plt.ylabel(var2)
    plt.legend(loc='best')
    plt.savefig('MinosContour_%s_%s.png' %(var1,var2))
    

PlotMinosContour('B','V3', SigmaMax=3, npoints=20)
    
Plot(Minuit, dataALL, dataBKG, phi_s, c)



