import cmath
import numpy as np
import matplotlib.pyplot as plt



#Classical model

ti = np.arange(0,0.65,0.01)
tr = np.arange(0.65,2.6,0.01)
w0p = 2*np.pi*300/800
Ipnou = 21.5645
Upnou = 23.8848
hw = (6.582*10**-16)*(w0p/(1*10**-15))
vi = (np.sin((np.pi/2) - 3*np.arcsin((2/np.pi)*w0p*ti-1))-np.sin(w0p*ti))
vr = (np.sin(w0p*tr) -np.cos((np.pi/2)*np.sin((1/3)*(-(np.pi/2)+w0p*tr))))
E1 = (Ipnou+2*Upnou*vi**2)/hw
E2 = (Ipnou+2*Upnou*vr**2)/hw





#Saddle points


Ip= 0.792
Up = 0.876
w0 = 0.057
npoints = 500

phirret1 = []
phirret2 = []
harmonicnumber = []


harmonics = np.arange(15,83,2)

phireturnreal = np.linspace(0.25*2*np.pi,1*2*np.pi,npoints) #limits based on the plot from the graph in the book
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
plt.plot(ti,E1,'o',color = 'b',label ="Classical ")
plt.plot(tr,E2,'o',color = 'b')
plt.grid()
plt.xlabel('t (fs)')
plt.ylabel('Harmonic order')
plt.suptitle('Comparison between Saddle point and classical methods')
plt.title('Ne at 800nm with I = 4*10^14 W/cm^2')
plt.legend()
plt.show()


stack = np.stack((harmonicnumber,phirret1,phiion1again),axis=1)

np.savetxt('test.txt',stack,delimiter=',',newline='\r\n',fmt='%8.8f')


'''
Coments regarding comparsion: 

It can be seen that the shape of the saddle points method reaches higher harmonics meaning that the approximation 
for the energy in the classical method is not that accurate since it does not take into account the extra value 
of the kinetic energy than the returning electron accumulates between the origin and the exit point of the barrier.
'''