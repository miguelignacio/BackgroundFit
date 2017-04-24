import numpy as np
from math import cos
from math import sin
from math import pi
from math import sqrt

def TotalPDF(x, params_Signal, params_BKG , phi,c):
    return Signal(x,params_Signal) + Background(x,params_BKG, phi, c)

def Signal(x, params):

    A1 = params["A1"]
    A2 = params["A2"]
    s1 = params["s1"]
    s2 = params["s2"]
    C1 = params["C1"]
    return A1*1/(s1*sqrt(2*pi))*np.exp(-(x-0.0)**2/(2*s1**2) ) + A2*1/(s2*sqrt(2*pi))*np.exp(-(x-1.0)**2/(2*s2**2) ) + C1

def Background(x, params, phi , c ):
    B    = params["B"]
    v2_t = params["v2_t"]
    v2_a = params["v2_a"]
    v4_t = params["v4_t"]
    v4_a = params["v4_a"]
    V1   = params["V1"]
    V3   = params["V3"]
    num = v2_t + cos(2*phi)*sin(2*c)/(2*c) +  v4_t*cos(2*phi)*sin(2*c)/(2*c) +  v2_t*cos(4*phi)*sin(4*c)/(4*c) + v4_t*cos(6*phi)*sin(6*c)/(6*c)
    den =  1 + 2*v2_t*cos(2*phi)*sin(2*c)/(2*c) + 2*v4_t*cos(4*phi)*sin(4*c)/(4*c)
    v2R = num/den
    num2 = v4_t + cos(4*phi)*sin(4*c)/(4*c) + v2_t*cos(2*phi)*sin(2*c)/(2*c)+ v2_t*cos(6*phi)*sin(6*c)/(6*c) + v4_t*cos(8*phi)*sin(8*c)/(8*c)
    v4R = num2/den
    BR = B*den*c*2/np.pi
    factor = 1.0
    if(c==np.pi/12.0): factor=2.0
    BR = BR*factor
    return BR*(1 + 2*V1*np.cos(pi*x) + 2*v2R*v2_a*np.cos(2*pi*x) + 2*V3*np.cos(3*pi*x) + 2*v4R*v4_a*np.cos(4*pi*x))