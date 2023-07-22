#!/usr/bin/env 
import pandas as pd
import numpy as np
import os.path

import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.ROOT)

from scipy.interpolate import interp2d

import hist
from hist import Hist

folder_results = "singlet/13TteV/original_logs/full"

sinTheta = [ss/10 for ss in range(10,0, -1)]
lam112   = range(-300, 400, 100)
masses = [300, 600, 1000]
print(len(sinTheta), len(lam112))

size_x = len(lam112)
size_y = len(sinTheta)
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
last_value = 0.01408 # 300, 0
for mass in masses: 
    for ss in sinTheta:
        for ll in lam112:
            sinTheta_file = str(ss).replace(".","p") if ss != 0 else "0"
            lam112_file   = str(ll).replace("-","m")
            file_name = "%s/log_M%s_ST%s_L%s_K1.log" % (folder_results, str(mass), sinTheta_file, lam112_file)
            check_file = os.path.isfile(file_name)
            cross_section = -0.1

            if check_file:
                with open(file_name) as ff:
                    for line in ff:
                        if "Cross-section :" in line:
                            cross_section = line.split()[2]
                            last_value = cross_section
            else:
                #df_missing[count_entry_missing] = [mass, ss/10, ll, 1]
                count_entry_missing += 1
                cross_section = last_value
                

            #print(file_name, mass, ss/10, ll, check_file, cross_section)
            def get_value(x):
                if type(x) == str:
                    return eval(x)
                else:
                    return x
            #if cross_section == -0.1:
            #    print("ERROR", file_name, mass, ss/10, ll, check_file, cross_section)
            df.loc[count_entry] = [get_value(mass), get_value(ss), get_value(ll), 1,get_value(cross_section)]
            count_entry += 1 

## fill with quadratic with lam112 beyond the computation
#lam112   = range(400, 1000, 100)
#for mass in masses: 
#    for ss in sinTheta:
#        for ll in lam112:
#            print(df[(df[]== &)])

df.to_csv("singlet/13TteV/XS.csv", index=False)
print("total entries expected", count_entry)
print("missing entries", count_entry_missing)

zlabel = "XS (pp -> X -> HH) (pb)"

# that did not worked well, but for docs
# https://github.com/bfonta/genproductions/blob/master/bin/MadGraph5_aMCatNLO/plots/plot_xsec.py#L131
# cbar = hep.hist2dplot(
#     df['CX(pb)'].values,  #z, 
#     df['lam112'].values, # x, 
#     df['sinTheta'].values, #y, 
#     flow=None)
# cbar.cbar.ax.set_ylabel(zlabel, rotation=0, labelpad=labelpad)
# cbar.cbar.ax.tick_params(axis='y', labelrotation=0)

masses_with_predictions = {
    300 : 0.273, 
    600 : 0.03115234375, 
    1000 : 0.004833984375
    }

for mass in masses:
    fig = plt.figure(figsize=(8, 8))
    plt.xlabel('$\\lambda_{112}$') # lam112
    plt.ylabel('$\\sin\\theta$') # sinTheta

    CX = df['CX(pb)'][df['mass']==mass].values.reshape(10,7)
    # print("lambda",df['lam112'].values.reshape(30,7))
    print("==================== limit", masses_with_predictions[mass])
    try:
        print("==================== limit", masses_with_predictions[mass])
        CS1 = plt.contour(lam112, sinTheta, CX)
        plt.clabel(CS1, levels=[masses_with_predictions[mass]], inline=True, fontsize=10)
    except:
        pass

    plt.imshow(CX, cmap='hot', interpolation='nearest')
    plt.xticks(np.arange(0, size_x, 1), df['lam112'][df['mass']==mass].values.reshape(size_y,size_x)[0],rotation=45)
    plt.yticks(np.arange(0, size_y, 1), df['sinTheta'][df['mass']==mass].values.reshape(size_y,size_x)[:,0])
    plt.colorbar()


    for ext in ('.png', '.pdf'):
        plt.savefig("singlet/13TteV/XS_MX_%s%s" % (mass, ext))

#"""
for mass in masses:
    # https://stackoverflow.com/questions/39727040/matplotlib-2d-plot-from-x-y-z-values

    ff = interp2d(
        df['lam112'].values,   #x_list,
        df['sinTheta'].values, #y_list,
        df['CX(pb)'].values,   #z_list,
        kind="quintic")

    lam112_scan = np.arange(-300, 1000, 10)
    sinTheta_scan = np.arange(0, 1, 0.01)

    fig = plt.figure(figsize=(8, 8))
    for ll in lam112:
        plt.plot(
            df[(df['lam112']==ll) & (df['mass']==mass) ]['sinTheta'].values,
            df[(df['lam112']==ll) & (df['mass']==mass) ]['CX(pb)'].values,
            label="lam112 = %s" % ll,
            marker='o',
        )
        #interpolated_CX = [ff(ll, sinTheta) for sinTheta in sinTheta_scan]
        #plt.plot(
        #    sinTheta_scan,
        #    interpolated_CX,
        #    linestyle='dashed', marker='o',
        #)
    #plt.set_title("MX = %s GeV" % str(mass))
    plt.legend(loc="upper right", title="MX = %s GeV" % str(mass))
    #for ext in ('.pdf'):
    plt.savefig("singlet/13TteV/XS_intepolated_MX_%s_sTheta.pdf" % (mass))

    fig = plt.figure(figsize=(8, 8))
    for ll in sinTheta:
        plt.plot(
            df[(df['sinTheta']==ll) & (df['mass']==mass) ]['lam112'].values,
            df[(df['sinTheta']==ll) & (df['mass']==mass) ]['CX(pb)'].values,
            label="sinTheta = %s" % ll
        )
    #plt.set_title("MX = %s GeV" % str(mass))
    plt.legend(loc="upper right", title="MX = %s GeV" % str(mass))
    #for ext in ('.pdf'):
    plt.savefig("singlet/13TteV/XS_intepolated_MX_%s_lam112.pdf" % (mass))


    fig = plt.figure(figsize=(8, 8))
    CX = df['CX(pb)'][df['mass']==mass].values.reshape(10,7)
    #contours = plt.contourf(lam112, sinTheta, CX, 3)
    CS = plt.contour(lam112, sinTheta, CX)

    plt.clabel(CS, CS.levels, inline=True, fontsize=10)
    plt.savefig("singlet/13TteV/XS_intepolated_MX_%s.pdf" % (mass))
 #"""

# https://stackoverflow.com/questions/39727040/matplotlib-2d-plot-from-x-y-z-values
# keep on naming convention of file kl1 (is that the 'K) ?
# extens the lam112 to 1000
# include sinTheta = 1

# draw 2D as histogram and as interpolated colormap