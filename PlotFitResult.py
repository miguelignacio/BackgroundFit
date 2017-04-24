
import numpy as np
from iminuit import Minuit, describe, Struct
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors as mcolors
from FitModel import Background, Signal, TotalPDF



mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['axes.labelsize']  = 14
#mpl.rcParams['legend.fontsize'] = 0.5 
#mpl.rcParams['legend.linewidth'] = 0.8
mpl.rcParams['lines.markersize'] = 3.8
mpl.rcParams['legend.numpoints'] =1
mpl.rcParams['legend.fontsize'] = 'small'
from math import cos
from math import sin
from math import pi


def Plot(m, dataALL, dataBKG, phi_s, c):

    fitresult_BKG = {}
    fitresult_BKG["B"]    = m.values["B"]
    fitresult_BKG["v2_t"] = m.values["v2_t"]
    fitresult_BKG["v4_t"] = m.values["v4_t"]
    fitresult_BKG["v2_a"] = m.values["v2_a"]
    fitresult_BKG["v4_a"] = m.values["v4_a"]
    fitresult_BKG["V1"]   = m.values["V1"]
    fitresult_BKG["V3"]   = m.values["V3"]
 
    fitresult_Signal = {}

    fitresult_Signal["Outplane"] = {}
    fitresult_Signal["Outplane"]["A1"] = m.values["A3"]
    fitresult_Signal["Outplane"]["A2"] = m.values["A4"]
    fitresult_Signal["Outplane"]["s1"] = m.values["s3"]
    fitresult_Signal["Outplane"]["s2"] = m.values["s4"]
    fitresult_Signal["Outplane"]["C1"] = m.values["C2"]

    fitresult_Signal["Midplane"] = {}
    fitresult_Signal["Midplane"]["A1"] = m.values["A5"]
    fitresult_Signal["Midplane"]["A2"] = m.values["A6"]
    fitresult_Signal["Midplane"]["s1"] = m.values["s5"]
    fitresult_Signal["Midplane"]["s2"] = m.values["s6"]
    fitresult_Signal["Midplane"]["C1"] = m.values["C3"]

    fitresult_Signal["Inplane"] = {}
    fitresult_Signal["Inplane"]["A1"] = m.values["A7"]
    fitresult_Signal["Inplane"]["A2"] = m.values["A8"]
    fitresult_Signal["Inplane"]["s1"] = m.values["s7"]
    fitresult_Signal["Inplane"]["s2"] = m.values["s8"]
    fitresult_Signal["Inplane"]["C1"] = m.values["C4"]


    f, axes = plt.subplots(1,3, sharex=True, sharey=True)


    xfit = np.linspace(-0.5, 0.5, num=150, endpoint=True)
    xfit_all = np.linspace(-0.5, 1.5, num=150, endpoint=True)
    ##FIT RESULT
    for j, key in enumerate(dataALL.keys()):
        model = Background(xfit, fitresult_BKG, phi_s[key], c[key])
        axes[j].plot(xfit, model,'-', color='dodgerblue') ##plot fit function
        model = TotalPDF(xfit_all, fitresult_Signal[key], fitresult_BKG, phi_s[key], c[key])
        axes[j].plot(xfit_all, model,'-', color='darkorange') ##plot fit function
        x = dataALL[key]["x_center"]
        y = dataALL[key]["y"]
        dy = dataALL[key]["dy"]
        axes[j].errorbar(x,y,yerr=dy, fmt='o', label='All, \n $|\Delta\eta|<0.5$', color = 'darkorange')
        x = dataBKG[key]["x_center"]#
        y = dataBKG[key]["y"]
        dy = dataBKG[key]["dy"]
        axes[j].errorbar(x,y,yerr=dy, fmt='o', label='All, \n $1.0<|\Delta\eta|<1.5$', color='dodgerblue')
        axes[j].set_title(key)
        axes[j].set_xlabel(r' $\Delta\phi/\pi$')

        axes[j].locator_params(axis='x', nticks=3)
    axes[0].set_ylabel(r' $dN/\Delta\phi$')
    axes[0].legend(loc='best', borderaxespad=0., frameon=False)
    f.subplots_adjust(hspace=0, wspace=0)
    plt.subplots_adjust(hspace=0, wspace=0)
    f.savefig('fitresult_A.png')
    
    
    for j, key in enumerate(dataALL.keys()):
       modelSignal = Signal(xfit_all, fitresult_Signal[key])
       axes[j].plot(xfit_all, modelSignal, '-', color='tomato')
    f.savefig('fitresult_B.png')