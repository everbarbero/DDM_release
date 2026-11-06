import numpy as np
import matplotlib.pyplot as plt

density = 11 #enough density to a) display some LCF data w/lambda<0.02, and
#								b) write NcyclesLambdaFit.inp
# data from 10.1016/j.ijfatigue.2022.106811 
filename = "F8-1Lc.txt" # Figure 8, row 1, L:left plot, c:circle symbol
fnamebase = filename.removesuffix("c.txt").replace(".",",")
data1Lc = np.loadtxt(filename, skiprows=0)
filename2 = "F8-1Ls.txt" # L:left, s:square
data1Ls = np.loadtxt(filename2, skiprows=0)
filename3 = "F8-1Rs.txt" # R:right, s:square
data1Rs = np.loadtxt(filename3, skiprows=0)
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
GIc = 1.6 # previous study
GIIc = 2*GIc # our guess
Tg	= 190
Tgwet = 150
TCure = 175
CST	= 130 #continuous service temperature
# 130 < SFT < 190 range SFT
# Tmin = -Tmax : Ever's guess
x1 = data1Lc[:, 0]# cycles plot row 1 L(left) in Fig. 
y1 = data1Lc[:, 1]# crack density
x2 = data1Ls[:, 0]# cycles plot row1 L(left)
y2 = data1Ls[:, 1]# crack density
x3 = data1Rs[:, 0]# cycles plot row1 R(left)
y3 = data1Rs[:, 1]# crack density
# concatenate
x = np.concatenate((x1, x2, x3))
y = np.concatenate((y1, y2, y3))

# Fit with model equation to smooth out the experimental data
coeff = np.polyfit(np.log10(x), y, 2) # y = a (log10​x)^2 + b (log10​x) + c
a, b, c = coeff # a=quadratic coeff, b=slope, c=y-intercept
params = np.array([a,b,c]).round(5)
print("a=quadratic coeff, b=slope, c=y-intercept:",params)

# generate points
logN = np.linspace(np.log10(501), np.log10(30501), density) 
yN = np.polyval(coeff, logN)

# save data to file for DDM6 KEY=-2 or you could write NcyclesLambdaFit.inp with experimental data
xy = np.column_stack((10**logN, yN))
sorted_indices = np.argsort(xy[:, 0])
xy_sorted = xy[sorted_indices]
np.savetxt('NcyclesLambdaFit.txt', xy_sorted, fmt='%12.5e', header=' Ncycles      lambda_fit_F8-row-1', comments='')

fig, ax = plt.subplots()
xlabel = "N cycles"
ylabel = "crack density $\\lambda$ (1/mm)"
ax.semilogx(x1, y1, label="74% UTS Group 1", linestyle='none', marker='o', color="blue",fillstyle='none',)
ax.semilogx(x2, y2, label="74% UTS Group 1", linestyle='none', marker='s', color="black",fillstyle='none',)
ax.semilogx(x3, y3, label="74% UTS Group 2", linestyle='none', marker='^', color="brown",fillstyle='none',)
ax.text(2000, 0, "[0_2/90_4]_s Hexply F6376 Fig. 8", fontsize=10, color="black")
ax.semilogx(10**logN, yN, linestyle='-', marker='.',label="quadratic fit")
ax.annotate(f"quadratic={a:.3g}", (1362, .18), horizontalalignment='right',color="black",fontsize=10)
ax.annotate(f"    slope={b:.3g}", (1362, .16), horizontalalignment='right',color="black",fontsize=10)
ax.annotate(f"intercept={c:.3g}", (1362, .14), horizontalalignment='right',color="black",fontsize=10)

plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend()
# plt.title("Crack density vs. N cycles")
plt.savefig(fnamebase)
plt.show()
