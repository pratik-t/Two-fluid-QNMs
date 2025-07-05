import matplotlib.pyplot as plt
import numpy as np

def DM_graph(en, eos_name, mass_frac, DM_mass, cs, cv, ls):
    
    file_name = f"{eos_name}_{mass_frac:.2f}_{str(cs).zfill(3)}_{str(cv).zfill(3)}_data.dat"
    
    file_path = f"./DATA/DDATA_{DM_mass}_{en}/{file_name}"
    print(file_path)
    a  = 0.9

    try:

        en, ed, rn, rd, mn, md, yR, g, f, p = np.loadtxt(file_path, usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], skiprows=1, unpack=True)
        
        if len(mn)!=0:
            ws = np.ones(10)/10
            max_index= np.argmax(mn+md)
            
            en = np.convolve(en, ws, mode='valid')[:max_index + 1]
            ed = np.convolve(ed, ws, mode='valid')[:max_index + 1]
            rn = np.convolve(rn, ws, mode='valid')[:max_index + 1]
            rd = np.convolve(rd, ws, mode='valid')[:max_index + 1]
            mn = np.convolve(mn, ws, mode='valid')[:max_index + 1]
            md = np.convolve(md, ws, mode='valid')[:max_index + 1]
            yR = np.convolve(yR, ws, mode='valid')[:max_index + 1]
            g  = np.convolve(g, ws, mode='valid')[:max_index + 1]
            f  = np.convolve(f, ws, mode='valid')[:max_index + 1]
            p  = np.convolve(p, ws, mode='valid')[:max_index + 1]
            
            halo=0
            for i in range(len(rd)-1, -1, -1):
                if (rd[i] < rn[i]):
                    continue
                else:
                    halo = i
                    break
            core = 0
            for i in range(len(rd)-1, -1, -1):
                if (rd[i] > rn[i]):
                    continue
                else:
                    core = i
                    break
                    # Every mass below 'X' is DM Halo
                    # Every mass above 'X' is DM core
                    
            R_tide = np.maximum(rn, rd)
            Mass = mn+md
            b = Mass*1.4766/R_tide
            k2 = ((8/5)*(b**5)*((1-2*b)**2)*(2-yR+2*b*(yR-1)))/(2*b*(6-3*yR+3*b*(5*yR-8)) +
                4*(b**3)*(13-11*yR+b*(3*yR-2)+2*b**2*(1+yR)) + 3*((1-2*b)**2)*(2-yR+2*b*(yR-1))*np.log(1-2*b))
            tide = (2./3.)*k2*(R_tide/(Mass*1.4766))**5
            
            plt.plot(rn, mn+md, ls=ls)

    except FileNotFoundError:
        print(file_name, 'not found')


mf = [0.04]
dm = [1,1.0]
cs = [30]
cv = [50]

for i in range(len(mf)):
    for j in range(len(dm)):
        for k in range(len(cs)):
            for l in range(len(cv)):
                # DM_graph('lso5e5', 'DDME2', mf[i], dm[j], cs[k], cv[l], '-')
                # DM_graph('lso1e5', 'DDME2', mf[i], dm[j], cs[k], cv[l], '--')
                # DM_graph('1e6', 'DDME2', mf[i], dm[j], cs[k], cv[l], ':')
                DM_graph('GeV', 'DDME2', mf[i], dm[j], cs[k], cv[l], '-.')

# plt.ylim([0,1000])



# import numpy as np
# import matplotlib.pyplot as plt

# r,p1,p2,m1,m2,phi = np.loadtxt('fort.10', usecols=[0,1,2,3,4,5], unpack=True)
# plt.plot(r, p1)

# r, p1, p2, m1,m2,phi = np.loadtxt('fort.11', usecols=[0,1,2,3,4,5], unpack=True)
# plt.plot(r, p1, ':')

# r, y = np.loadtxt('jj.dat', usecols=[0,1], unpack=True)
# plt.plot(r, y, '.')

# r, y = np.loadtxt('ji.dat', usecols=[0, 1], unpack=True)
# plt.plot(r, y, '.')

# r, y = np.loadtxt('ni.dat', usecols=[0, 1], unpack=True)
# plt.plot(r, y, '.')
        

# plt.axhline(0,c='k')
# plt.yscale('log')
plt.show()