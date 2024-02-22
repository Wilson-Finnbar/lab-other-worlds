# %%
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as cn
import uncertainties as uc
import uncertainties.umath as um
np.set_printoptions(formatter={'float': lambda x: format(x, '6.2E')})

'''
star main has coordiantes x = 424.97 and y = 520.39
star 1 has coordinates x = 614.35 and y = 806.13
star 2 has coordinates x = 318.62 and y = 331.07
star 3 has coordiantes x = 355.5 and y = 211.32
'''
image_time = ['2011-08-23 04:32:20','2011-08-23 04:43:09','2011-08-23 04:55:25','2011-08-23 05:07:43','2011-08-23 05:21:05','2011-08-23 05:33:25','2011-08-23 05:45:00','2011-08-23 05:56:35','2011-08-23 06:08:11','2011-08-23 06:19:44','2011-08-23 06:31:20','2011-08-23 06:42:55']
time = ['04:32:20','04:43:09','04:55:25','05:07:43','05:21:05','05:33:25','05:45:00','05:56:35','06:08:11','06:19:44','06:31:20','06:42:55']
starmain = np.array([0.16949e+06,0.16348e+06,0.16628e+06,0.16089e+06,0.16490e+06,0.16781e+06,0.16973e+06,0.17165e+06,0.17244e+06,0.17564e+06,0.17711e+06,0.17804e+06])
star1 = np.array([99143,0.10254E+06,98864,95639,0.10180E+06,0.10334E+06,0.10423E+06,0.10638E+06,0.10495E+06,0.10522E+06,0.10577E+06,0.10538E+06])
star2 = np.array([0.10370E+06,0.11788E+06,0.11839E+06,0.12670E+06,0.13055E+06,0.13137E+06,0.13214E+06,0.13307E+06,0.13041E+06,0.13083E+06,0.13286E+06,0.13103E+06])
star3 = np.array([0.10750E+06,0.10256E+06,0.10108E+06,0.10733E+06,0.10950E+06,0.11302E+06,0.11248E+06,0.11271E+06,0.11356E+06,0.11508E+06,0.11351E+06,0.11419E+06])

flux1 = np.zeros(12)
flux2 = np.zeros(12)
flux3 = np.zeros(12)

for i in range(0,12):
    flux1[i] = starmain[i] / star1[i]
    flux2[i] = starmain[i] / star2[i]
    flux3[i] = starmain[i] / star3[i]

plt.plot(time,flux1, '.-', label='Flux 1')
plt.plot(time,flux2, '.-', label='Flux 2')
plt.plot(time,flux3, '.-', label='Flux 3')
plt.xticks(rotation=90)
plt.legend()
plt.xlabel('Time')
plt.ylabel('Normalised flux')
plt.xlim([min(time),max(time)])
plt.ylim(1.2,1.8)

import tikzplotlib
tikzplotlib.save("qatartransit.pgf", axis_height='8cm', axis_width='13cm')

dt = 6450
R = uc.ufloat(0.823,0.025) * cn.R_sun.value
M = uc.ufloat(0.85,0.03) * cn.M_sun.value
G = cn.G.value
maxf = uc.ufloat(1.68,0.01)
minf = uc.ufloat(1.64,0.01)

a = (dt**2 * G * M)/(4*R**2)
P = (dt*np.pi*a)/(R)
Rsq = (maxf-minf)/(maxf)
R_p = um.sqrt((maxf-minf)/(maxf)) *R
print(f"{a/cn.au.value}\n{P/86400}\n{R_p / cn.R_jup.value}\n{Rsq}")

'''
from 4:55:25 to 6:42:55
error in sum in apature is sqrt(count) which is 1000 this has minimal impact on the values
'''
# %%
from uncertainties import unumpy

phase = np.array([0.2014,0.8320,0.5648,0.2350,0.6610,0.0581,0.7582,0.4653,0.1646])
velocity = np.array([-37.9647,-37.4835,-37.5939,-37.9013,-37.5739,-37.8350,-37.4732,-37.7496,-37.8753])
velocity = velocity*1000
uncertainty = np.array([9.2622,4.6936,14.1425,8.3571,6.8672,7.0690,7.0483,25.3662,7.2845])
joined = unumpy.uarray(velocity,uncertainty)

vel_w = np.zeros(len(velocity))
mean_v = np.average(velocity)
for i in range(len(velocity)):
    vel_w[i] = (velocity[i] - mean_v)

x_values = np.arange(0,1,0.01)

def sinfunc(x):
    return (0.25*1000) * np.sin( 2*np.pi*x+np.pi)

plt.plot(x_values,sinfunc(x_values))
plt.plot(phase,vel_w, '.')
plt.xlim(0,1)
plt.ylim(-300,300)
plt.xlabel('Phase')
plt.ylabel('Relative Velocity (m/s)')
print(sum(joined)/len(joined))
tikzplotlib.save("qatarradial.pgf", axis_height='8cm', axis_width='13cm')

M_psini = (250 * M**(2/3) * 122688**(1/3) ) / ((2 * np.pi * G)**(1/3))
print(M_psini / cn.M_jup.value)
print((M_psini / (4/3 * np.pi * (R_p)**3)))
# %%
print(f'{starmain}')
# %%
import pandas as pd
#pd.set_option('display.float_format', lambda x: '%6.2E' % x)
data = {
    'Qatar1-b': starmain,
    'Star 1': star1,
    'Star 2': star2,
    'Star 3': star3,
    'Flux 1': flux1,
    'Flux 2': flux2,
    'Flux 3': flux3
}
df = pd.DataFrame(data, index = time)
#print(df.to_latex())

radial = {
    'Phase': phase,
    'Measured Velocity (m/s)': velocity,
    'Relative Velocity (m/s)': vel_w,
    'Uncertainty (m/s)': uncertainty
}
dfr = pd.DataFrame(radial)
print(dfr.to_latex())

# %%
import uncertainties as uc

Ar = uc.ufloat([186,46])

print(0.04*0.03*(197/Ar))
# %%
