import cmath
import numpy as np
import matplotlib.pyplot as plt


#Saddle points


Ip= 0.57915
Up = 0.658307
w0 = 0.057
w0p = 2*np.pi*300/800
npoints = 1000

phirret1 = []
phirret2 = []
harmonicnumber = []


harmonics = np.arange(11,49,2)

phireturnreal = np.linspace(0.25*2*np.pi,1*2*np.pi,npoints)
phireturnima = np.linspace(-0.1*2*np.pi,0.05*2*np.pi,npoints)
[phirreturn1,phirreturn2] = np.meshgrid(phireturnreal,phireturnima)

for i in harmonics:
    gamman = np.sqrt(((i * w0) - Ip) / (2 * Up))
    p1 = np.cosh(phirreturn2) * np.sin(phirreturn1) + gamman
    p2 = np.sinh(phirreturn2) * np.cos(phirreturn1)
    gammabar = np.sqrt(Ip/2*Up) + p2
    P = p1 ** 2 + gammabar ** 2 + 1
    D = np.sqrt(P ** 2 - 4 * (p1 ** 2))
    phiion1 = np.arcsin(np.sqrt((P - D) / 2))
    phiion2 = np.arccosh(np.sqrt((P + D) / 2))
    F1 = p1 * (phirreturn1 - phiion1) - p2 * (phirreturn2 - phiion2) - np.cos(phiion1) * np.cosh(phiion2) + np.cosh(
        phirreturn2) * np.cos(phirreturn1)
    F2 = p1 * (phirreturn2 - phiion2) + p2 * (phirreturn1 - phiion1) + np.sin(phiion1) * np.sinh(phiion2) - np.sinh(
        phirreturn2) * np.sin(phirreturn1)
    F = F1 ** 2 + F2 ** 2
    [m, n] = np.where(F == np.min(F))
    m = int(m)
    n = int(n)
    harmonicnumber.append(i)
    phirret2.append(float(phireturnima[[m]]))
    phirret1.append(float(phireturnreal[[n]]))
    # Empty the cells
    F[m,:] = None
    F[:,n] = None
    # Second miminum
    [a, b] = np.where(F == np.nanmin(F))
    a = int(a)
    b = int(b)
    harmonicnumber.append(i)
    phirret2.append(float(phireturnima[[a]]))
    phirret1.append(float(phireturnreal[[b]]))


#Ionization phase





rango = np.arange(0,len(harmonicnumber),1)
phiion1again = []
phiion2again = []

for i in rango:
    gamman = np.sqrt(((float(harmonicnumber[i]) * w0) - Ip) / 2 * Up)
    p1 = (np.cosh(phirret2[i]) * np.sin(phirret1[i])) + gamman
    p2 = np.sinh(phirret2[i]) * np.cos(phirret1[i])
    gammabar = np.sqrt(Ip/2*Up) + p2
    P = p1 ** 2 + gammabar ** 2 + 1
    D = np.sqrt((P ** 2 - 4 * p1**2))
    phiion1again.append(np.arcsin(np.sqrt((P-D)/2))/w0p)
    phiion2again.append(np.arccosh(np.sqrt((P+D)/2)))

phirret1[:] = [x / w0p for x in phirret1]


#Plotting


plt.figure(1)
plt.plot(phirret1,harmonicnumber,'*',color = 'r',label ="Saddle points")
plt.plot(phiion1again,harmonicnumber,'*',color = 'r')
plt.grid()
plt.xlabel('t (fs)')
plt.ylabel('Harmonic order')
plt.suptitle('Saddle point method for Argon')
plt.title('Ar at 800nm with I = 3*10^14 W/cm^2')
plt.legend()
plt.show()


