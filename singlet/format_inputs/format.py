#!/usr/bin/env 
import pandas as pd
import numpy as np
import os.path

import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.ROOT)

import hist
from hist import Hist

folder_results = "singlet/13TteV/original_logs"

sinTheta = range(0,10, 1)
lam112   = range(-300, 400, 100)
masses = [300, 600, 1000]

#missing_points = "%s/missing_points.txt" % folder_results

df = pd.DataFrame(
    columns=(
        'mass', 
        'sinTheta', 
        'lam112',
        'kl', 
        'CX(pb)', 
        )
    )

#df_missing = pd.DataFrame(
#    columns=(
#        'mass', 
#        'sinTheta', 
#        'lam112',
#        'kl', 
#        )
#    )

count_entry = 0
count_entry_missing = 0
for mass in masses: 
    for ss in sinTheta:
        for ll in lam112:
            sinTheta_file = str(ss/10).replace(".","p") if ss != 0 else "0"
            lam112_file   = str(ll).replace("-","m")
            file_name = "%s/log_M%s_ST%s_L%s_K1.log" % (folder_results, str(mass), sinTheta_file, lam112_file)
            check_file = os.path.isfile(file_name)
            cross_section = -0.1

            if check_file:
                with open(file_name) as ff:
                    for line in ff:
                        if "Cross-section :" in line:
                            cross_section = line.split()[2]
            else:
                #df_missing[count_entry_missing] = [mass, ss/10, ll, 1]
                count_entry_missing += 1

            #print(file_name, mass, ss/10, ll, check_file, cross_section)
            def get_value(x):
                if type(x) == str:
                    return eval(x)
                else:
                    return x
            if cross_section == -0.1:
                print("ERROR", file_name, mass, ss/10, ll, check_file, cross_section)
            df.loc[count_entry] = [get_value(mass), get_value(ss)/10, get_value(ll), 1,get_value(cross_section)]
            count_entry += 1 

print (df[df["mass"] == 300])

df.to_csv("singlet/13TteV/XS.csv", index=False)
print(count_entry)
#print (df_missing)
print(count_entry_missing)

print(    df['CX(pb)'].values,  #z, 
    df['lam112'].values, # x, 
    df['sinTheta'].values, #y, 
    )
zlabel = "XS (pp -> X -> HH) (pb)"
fig = plt.figure(figsize=(8, 8))
plt.xlabel('$\\lambda_{112}$') # lam112
plt.ylabel('$\\sin\\theta$') # sinTheta

# test_hist = Hist(hist.axis.Regular(10, -1, 1, name="sinTheta", label="sinTheta"),hist.axis.Regular(10, -300, 300, name="lam112", label="lam112"), hist.axis.Regular(10, 0, 1000, name="CX(pb)", label="CX(pb)"))
# test_hist.fill(df['sinTheta'].values, df['lam112'].values, df['CX(pb)'].values)

# test_hist.project("CX(pb)").plot()

# cbar = hep.hist2dplot(
#     df['CX(pb)'].values,  #z, 
#     df['lam112'].values, # x, 
#     df['sinTheta'].values, #y, 
#     flow=None)
# cbar.cbar.ax.set_ylabel(zlabel, rotation=0, labelpad=labelpad)
# cbar.cbar.ax.tick_params(axis='y', labelrotation=0)

print("test",df[df['mass']==300])

CX = df['CX(pb)'][df['mass']==300].values.reshape(10,7)
# print("lambda",df['lam112'].values.reshape(30,7))

plt.imshow(CX, cmap='hot', interpolation='nearest')
plt.xticks(np.arange(0, 7, 1), df['lam112'][df['mass']==300].values.reshape(10,7)[0],rotation=45)
plt.yticks(np.arange(0, 10, 1), df['sinTheta'][df['mass']==300].values.reshape(10,7)[:,0])
plt.colorbar()


for ext in ('.png', '.pdf'):
    plt.savefig("singlet/13TteV/XS_MX_300_" + ext)

"""
fig, ax = plt.subplots()
cc = ax.pcolormesh(
    df['sinTheta'].values, 
    df['lam112'].values, 
    df['CX(pb)'].values, 
    cmap='RdBu', 
    #vmin=z_min, vmax=z_max
    )
ax.set_title('MX = %s, kl = 1' % str(2))
# set the limits of the plot to the limits of the data
#ax.axis([x.min(), x.max(), y.min(), y.max()])
fig.colorbar(cc, ax=ax)
"""

# https://stackoverflow.com/questions/39727040/matplotlib-2d-plot-from-x-y-z-values
# keep on naming convention of file kl1 (is that the 'K) ?
# extens the lam112 to 1000
# include sinTheta = 1

# draw 2D as histogram and as interpolated colormap