import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import string
# Minimize the error over the whole lambda(eps) curve
# data from 10.1016/j.jcomc.2020.100005
filename = "1mmpermin.txt" # Experimental data (%F/Fmax,lambda) 
# Then, transform it (eps,lambda) and save in epsLambdaExperim.txt
fnamebase = filename.removesuffix(".txt").replace(".",",")
data = np.loadtxt(filename, delimiter=",", skiprows=1)
# *** DDM results from KEY=1 *** for KEY=-2 make sure you just run FILE5.inp with KEY=1 !!!
filename2 = "FILE7.out"

fnamebase = filename2.removesuffix(".out").replace(".",",")
data2 = np.loadtxt(filename2, skiprows=1)
# Hexply F6376 data for reference only
Xt = 2274 # MPa Turbin(2024) Table 1
Xc = 1849  # MPa
Yt = 102 # Mpa
Yc = 255 # MPa
S12 = 63 # MPa
S23 = 35 # Mpa
E1 = 142000 # MPa
E2 = 9100 # MPa
G12 = 5200 # MPa
v12 = 0.27
v23 = 0.3 #(G23 = 3500 MPa)
tt = 0.8 # mm, transition thickness estimate for C/Ep
a0 = tt/4 # crack size, 4.18 ICMD3e
L22 = E2*E1**3/2/(E1**3-v12**2*E2**3) # (4.74) in ICMD3e
GIc = 1.12**2*np.pi*a0*Yt**2/L22 # (4.76)
GIIc = np.pi*a0*S12**2/G12 # (4.79)
print("(4.76)GIc=",round(GIc,3),"(4.79)GIIc=",round(GIIc,3))
# LSS: [0_2/90_4]_s
UTSarray = np.array([707,685,691,673,699])# MPa, Table 3
LoadRate = np.array([.019,.19,1.9,19,1])# kN/sec,Table 3
UTS = np.average(UTSarray); print("UTS=",UTS)
# carefull: using percentage strain
x = data[:, 0] # %F/Fmax
y = data[:, 1] # crack density [1/mm]
x = x * UTS/E1 # Fig.3: %F/Fmax*(A0/A0)=sigma/XT, mult. by XT/E1 = sigma/E1 = %eps1
# Transformed to (eps,lambda) in epsLambdaExperim.txt, to be read with DDM KEY=-2
xydata = np.column_stack((x/100, y))# build table with actual strain and lambda values
np.savetxt("epsLambdaExperim.txt", xydata, fmt="%.6f", header="epsx.... lambda", comments="") # delimiter=", "
# Fit a line
coeff = np.polyfit(x, y, 1)
a, b = coeff # a=slope, b=y_intercept, y = b + a*x
epsx0 = -b/a # x_intercept & crack initiation at lambda=0.02
print("slope=",round(a,3), "y_intercept=",round(b,3), "x_intercept=",round(epsx0,3))
# Evaluate the fitted line
y_fit = np.polyval(coeff, x)
# Fit a quadratic 
# coeff = np.polyfit(x, y, 2)
# a, b, c = coeff
# print(a, b, c)

# get DDM results 
x2 = data2[:,1] # strain
y2 = data2[:,3]
# find lowest x2 when y2>0.02
idx = np.where(y2 > 0.02)[0]
if len(idx) > 0:
    xmin = x2[idx[0]]
else:
    xmin = None

print(f"first crack at epsx = {xmin:.3g}") # formatting-string
x2 = data2[:,1]*100 # transform to %strain for visualization
fig, ax = plt.subplots()
# Minor ticks 
# ax.xaxis.set_minor_locator(MultipleLocator(.01)) # spaced this much ()
# ax.grid(which='minor', linestyle=':', linewidth=0.5) #turn minor grid on

xlabel = "$\\epsilon_x=(F/F_{max})*UTS/E_1$ [%]"
ylabel = "crack density [1/mm]"

ax.plot(x, y, label="load-rate 1 mm/min Fig. 3", linestyle='none', marker='o', color="grey", mfc='none')
ax.text(0., .18, "[0_2/90_4]_s Hexply F6376 UTS=691", fontsize=9, color="black")

ax.plot(x2[::3], y2[::3], label="DDM6 GIc=0.157", linestyle='none', marker='+', color="blue", ms=10, mew=1)
ax.axhline(0, color='black', linewidth=1, zorder=0)#zorder=10 puts it in front of the data

ax.minorticks_on()
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend()
# plt.title("Crack density vs. Strain")
plt.savefig(fnamebase)
plt.show()
