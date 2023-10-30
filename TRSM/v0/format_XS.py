#!/usr/bin/env 
import pandas as pd
import numpy as np
import time
from array import array
import seaborn as sns
# based on https://gitlab.cern.ch/cms-b2g/diboson-combination/combination-2016/-/blob/master/hvt.py

import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

name_br = "heavy_SM_higgs/Heavy_Higgs_BR.csv"
name_cx = "TRSM/v0/XS_original.csv"

Hbb = 5.77E-01
Hgg = 2.28E-03

df_cx = pd.read_csv(name_cx)
df_br = pd.read_csv(name_br)

# The Heavy Higgs is the Y-particle here
df_br.rename({
    'MX[GeV]': 'MS', 
    }, axis=1, inplace=True)


df_cx = df_cx.merge(df_br, on=['MS']).fillna(method='ffill')

print(df_cx[(df_cx["MH3"] == 1000)])
print(df_cx[(df_cx["MS"] == 125)])

print(df_cx.columns)
df_cx["X_to_HY_from_4b(fb)"] = df_cx["Xbbbb(fb)"]/df_cx["H_to_bb"]
df_cx["X_to_HY_from_bbgg(fb)"] = df_cx["Xbbgg+ggbb(fb)"]/(Hgg*df_cx["H_to_bb"] + Hbb*df_cx["H_to_gg"])
for column in ["X_to_HY_from_4b(fb)", "X_to_HY_from_bbgg(fb)"] :
    df_cx[column] = np.where(df_cx[column] < 0, -1, df_cx[column])

#df_cx["X_to_HY(fb)_1"] = df_cx[["X_to_HY_from_4b(fb)", "X_to_HY_from_bbgg(fb)"]].min(axis=1)

df_cx["X_to_HY(fb)"] = -1

df_cx["X_to_HY(fb)"] = np.where(df_cx["X_to_HY_from_4b(fb)"] < 0, df_cx["X_to_HY_from_bbgg(fb)"], df_cx["X_to_HY(fb)"])
df_cx["X_to_HY(fb)"] = np.where(df_cx["X_to_HY_from_bbgg(fb)"] < 0, df_cx["X_to_HY_from_4b(fb)"], df_cx["X_to_HY(fb)"])
df_cx["X_to_HY(fb)"] = np.where(
    (df_cx["X_to_HY_from_4b(fb)"] >0) & (df_cx["X_to_HY_from_bbgg(fb)"] > 0), 
    df_cx["X_to_HY_from_4b(fb)"], 
    df_cx["X_to_HY(fb)"])



print(df_cx[["MH3", "MS", "X_to_HY_from_4b(fb)", "X_to_HY_from_bbgg(fb)", "X_to_HY(fb)"]])

# Print validation
colors = ["r", "b", "g", "g", "y", "c"]
fig = plt.figure()
ax = fig.add_subplot()
for cc, column in enumerate(["X_to_HY_from_4b(fb)", "X_to_HY_from_bbgg(fb)", "X_to_HY(fb)"]) :
    # remove negative
    if column == "X_to_HY_from_4b(fb)" :
        linestyle = "-"
        label = "From (4b boosted)"
    elif column == "X_to_HY_from_bbgg(fb)" :
        linestyle = "-"
        label = r'From (bb$\gamma\gamma$)'
    else:
        label = "Used in reinterpretation"
        linestyle = "-."
    df_local = df_cx[(df_cx[column]> 0) & (df_cx["MS"]== 125)]

    plt.plot(
        df_local["MH3"].values, 
        df_local[column].values, 
        linestyle,
        label=label,
        color=colors[cc] )


plt.xlabel("MX (GeV)")
plt.ylabel(r'$\sigma(pp \rightarrow X \rightarrow$ HY) (fb)')
 
plt.yscale('log')
#plt.scale('log')
plt.title("TRSM cross section (MY = 125 GeV)")
plt.legend(loc='lower left', ncol=1) #, bbox_to_anchor=(-0.1,1.15),  frameon=True, edgecolor='black',framealpha=1,fancybox=False, ncol=4) 
nameout = "TRSM/v0/compare_versions"
for typePlot in ["pdf", "png"]:
    plt.savefig("%s.%s" % (nameout, typePlot))
    print("Saved XCheck figure %s.%s" % (nameout, typePlot))
plt.clf()

nameout = "TRSM/v0/XS_formatted"
df_cx.to_csv("%s.csv" % nameout, index=False)
print("saved in tabular format %s" % nameout)