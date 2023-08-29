#Parameters
difm = 0.015206
mum = 60*24*3*10**-6
cm = 0.1
chi = 0.3
barm = 350.0
baro = 400.0
rm = 60*24*6*10**-7

def initialConditionMicroglia(L, hx, xSize, microglia):
    tot = 0
    for k in range(xSize * xSize):
        i = int(k/xSize)
        j = k%xSize
        if (i - int(xSize/2))**2 + (j - int(xSize/2))**2 < 5 / (hx**2):
            microglia[i][j] = barm/3    
            tot = tot + barm/3
    print(tot)
    return microglia

from sympy import Derivative, Max, Min
from skfdiff import Model, Simulation
import pylab as pl
import numpy as np
from scipy.signal.windows import gaussian
from scipy import integrate
model = Model(["difm*(dxxM + dyyM) + mum*M*(barm - M) - cm*M- chi*(upwind(dxO, M/(barm + M), x, 1) + upwind(dyO, M/(barm + M), y, 1))",
               "rm*(M*M/(barm + M))*(baro - O)"],
               ["M(x, y)", "O(x, y)"],
               parameters=["difm", "mum", "cm", "chi", "barm", "baro", "rm"],
               boundary_conditions="noflux")
tmax = 28
dt =0.1
hx = 0.5
L = 20
xSize = int(L/hx)
x = y = np.linspace(-L / 2, L / 2, xSize)
xx, yy = np.meshgrid(x, y, indexing="ij")

M = np.zeros((int(L/hx), int(L/hx)))
M = initialConditionMicroglia(L, hx, xSize, M)
O = np.zeros_like(M)



initial_fields = model.Fields(x=x, y=y, M=M, O=O,
                              difm = difm, cm = cm, chi = chi, mum = mum, barm = barm, baro = baro, rm = rm)

simulation = Simulation(model, initial_fields, dt=dt, tmax=tmax, time_stepping=False)
container = simulation.attach_container()
tmax, final_fields = simulation.run()

container.data.O[-1].values.tolist()