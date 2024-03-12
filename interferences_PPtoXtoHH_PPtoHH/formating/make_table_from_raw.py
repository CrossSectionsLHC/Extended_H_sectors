#!/usr/bin/env 
import pandas as pd
import numpy as np
import time
from array import array
import seaborn as sns
# based on https://gitlab.cern.ch/cms-b2g/diboson-combination/combination-2016/-/blob/master/hvt.py

import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
import glob
import os.path

masses = [260, 270, 300, 320, 340, 350, 360, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000]

# Draw G/M of LHC TWIKI resonance only
name = 'heavy_SM_higgs/Heavy_Higgs_Total_Width.csv'
print(name)
df_width = pd.read_csv(name)
df_width["GoM"] = df_width["GammaH[GeV]"]/df_width["MH[GeV]"]
print(df_width)

plt.plot(df_width["MH[GeV]"], df_width["GoM"], marker=".") 
plt.xlabel("Mass [GeV]")
plt.ylabel("G/M")
plt.savefig('width_o_mass_lhc_twiki.png')
plt.clf()

name_XS = 'heavy_SM_higgs/Heavy_Higgs_13TeV_NNLO_NNLL_columb.csv'
df_cx_twiki = pd.read_csv(name_XS)
print(df_cx_twiki)
df_cx_twiki = df_cx_twiki.rename(columns={ "sigma[pb]" : "sigma[pb]_twiki"}) # "Mass[GeV]": "Mass[Gev]",
df_cx_twiki = df_cx_twiki[df_cx_twiki["Mass[Gev]"].isin(masses)]
#df_cx_twiki.reset_index()
#df_cx_twiki["Mass[Gev]"] = df_cx_twiki["Mass[Gev]"].astype(float)
#print(df_cx_twiki)

# TODO: will pass these to this git repo unde "v0"
eos_folder = "/eos/user/a/acarvalh/interferences_HH/raw_cross_sections/"
types_cx = ["SChan_eta0", "SChan_h_SChan_eta0", "BOX_SChan_eta0"] # 
#["vary_masses_result_width_0.001", "vary_masses_result_width_0.05"]


###################
# Inputs to compute non-resonant k factors
# to use when we have MX and width dependent k-factors 
###################
BOX_LO = 0.0305
SChan_h_LO	= 0.0055

BOX_NLO = 0.0704
BOX_SChan_h_NLO = 0.0111

###################
# Formatting cross sections
###################

g_o_ms = [0.001, 0.05, 0.01, 0.1]
df_cx_total = pd.DataFrame()

for g_o_m in g_o_ms:
    file_read = "CX_LHC13_decay_%s.txt" % g_o_m
    widths_read = "vary_masses_result_width_%s" % g_o_m
    df_cx = pd.DataFrame()

    for tt, type_cx in enumerate(types_cx):
        file_xs = "%s/%s/%s/%s" % (eos_folder,type_cx,widths_read,file_read)
        if not os.path.isfile(file_xs):
            print("DOES NOT EXIST (yet) %s " % file_xs)
            continue
        print("read %s " % file_xs)
        if tt == 0:
            df_cx = pd.read_csv(file_xs)
            #df_cx.reset_index()
            try:
                df_cx = df_cx.drop(['units'], axis=1) 
            except:
                pass

            df_cx = df_cx.rename(columns={"sigma[pb]": "sigma[pb]_%s" % (type_cx)})
            df_cx = df_cx.merge(df_cx_twiki, on=["Mass[Gev]"]  , how='outer')
        else:
            df_local = pd.read_csv(file_xs)
            try:
                df_local = df_local.drop(['units'], axis=1)
            except:
                pass

            df_local = df_local.rename(columns={"sigma[pb]": "sigma[pb]_%s" % (type_cx)})
            df_local = df_local.merge(df_cx_twiki, on=['Mass[Gev]']  , how='outer')

            df_cx = df_cx.merge(df_local, on=['Mass[Gev]', "sigma[pb]_twiki", 'scale_pos[%]',  "scale_neg[%]", "pdf_as[%]", "EW[%]", "aS[%]"]  , how='outer')

        df_cx["GoM"] = float(g_o_m)
    
    df_cx_total = pd.concat([df_cx_total, df_cx])

df_cx_total.to_csv("interferences_PPtoXtoHH_PPtoHH/cross_sections_for_mass_width.csv", index=False)
#print(df_cx_total)
#print(df_cx_total.keys)

df_cx_total["width_HH" ] = 31.803878252**2*(1-4*(125/df_cx_total['Mass[Gev]'])**2)**0.5/(8*3.1415*df_cx_total['Mass[Gev]'])
df_cx_total["BR_HH" ] = df_cx_total["width_HH"]/(df_cx_total["GoM"]*df_cx_total['Mass[Gev]'])
df_cx_total["sigma_prod_pb_%s" % ("SChan_eta0")] = df_cx_total["sigma[pb]_%s" % "SChan_eta0"]/df_cx_total["BR_HH"] 
df_cx_total["kfactor_%s" % ("SChan_eta0")] = df_cx_total["sigma[pb]_twiki"]/df_cx_total["sigma_prod_pb_%s" % ("SChan_eta0")]
df_cx_total["kfactor_SChan_eta0_eff"] = 2

print("IMPORTANT: in this table we are already normalizing the interference part to kt = 1")
print("Box-resonance scale to 2")
print("Triangle-resonance scale to 2*1.5 = 3 ")
df_cx_total["sigma[pb]_SChan_h_SChan_eta0"] = 3*df_cx_total["sigma[pb]_SChan_h_SChan_eta0"]
df_cx_total["sigma[pb]_BOX_SChan_eta0"] = 2*df_cx_total["sigma[pb]_BOX_SChan_eta0"]

print("IMPORTANT: in this table we already apply k-factors to to all resonant and resonant interference parts")
print("*** By now we are assuming kfactor resonance only = 2 *****")
df_cx_total["kfactor_SChan_h_SChan_eta0_eff"] = ((BOX_SChan_h_NLO/SChan_h_LO)*df_cx_total["kfactor_SChan_eta0_eff"])/2.
df_cx_total["kfactor_BOX_SChan_eta0_eff"]     = ((BOX_NLO/BOX_LO)*df_cx_total["kfactor_SChan_eta0_eff"])/2.

#print(df_cx_total[["kfactor_SChan_h_SChan_eta0_eff", "kfactor_BOX_SChan_eta0_eff"]])
df_cx_total["sigma[pb]_SChan_eta0"]         = df_cx_total["kfactor_SChan_eta0_eff"]*df_cx_total["sigma[pb]_SChan_eta0"]
df_cx_total["sigma[pb]_SChan_h_SChan_eta0"] = df_cx_total["kfactor_SChan_h_SChan_eta0_eff"]*df_cx_total["sigma[pb]_SChan_h_SChan_eta0"]
df_cx_total["sigma[pb]_BOX_SChan_eta0"]     = df_cx_total["kfactor_BOX_SChan_eta0_eff"]*df_cx_total["sigma[pb]_BOX_SChan_eta0"]

for plot_type in ["width_HH", "BR_HH", "sigma_prod_pb_%s" % ("SChan_eta0"), "kfactor_%s" % ("SChan_eta0")]:
    for g_o_m in g_o_ms : 
        plt.plot(df_cx_total["Mass[Gev]"].loc[df_cx_total["GoM"] == g_o_m], df_cx_total["%s" % plot_type].loc[df_cx_total["GoM"] == g_o_m], marker=".", label='G/M = %s' % str(g_o_m),linewidth=3.0) 
    plt.xlabel("Mass [GeV]")
    plt.ylabel(plot_type)
    plt.legend()
    plt.savefig('%s.png' % plot_type)
    plt.clf()

#################################################

plt.plot(df_cx_total["Mass[Gev]"].loc[df_cx_total["GoM"] == 0.001], df_cx_total["sigma[pb]_SChan_eta0"].loc[df_cx_total["GoM"] == 0.001], marker=".", label='Pure resonance, G/M = 0.001') 
plt.plot(df_cx_total["Mass[Gev]"].loc[df_cx_total["GoM"] == 0.001], df_cx_total["sigma[pb]_SChan_h_SChan_eta0"].loc[df_cx_total["GoM"] == 0.001], marker=".", label='Triangle X resonance, G/M = 0.001') 
plt.plot(df_cx_total["Mass[Gev]"].loc[df_cx_total["GoM"] == 0.001], df_cx_total["sigma[pb]_BOX_SChan_eta0"].loc[df_cx_total["GoM"] == 0.001], marker=".", label='Box X resonance, G/M = 0.001') 
plt.legend()
#plt.title("G/M = 0.1%")
plt.xlabel("Mass [GeV]")
plt.ylabel("sigma(pp to hh)[pb] ")
plt.savefig('cross_section_resonance_GoM_0p001.png')
plt.clf()

#################################################

for g_o_m in g_o_ms : 
    plt.plot(df_cx_total["Mass[Gev]"].loc[df_cx_total["GoM"] == g_o_m], df_cx_total["sigma[pb]_SChan_h_SChan_eta0"].loc[df_cx_total["GoM"] == g_o_m], marker=".", label='Triangle X resonance, G/M = %s'  % str(g_o_m)) 
for g_o_m in g_o_ms : 
    plt.plot(df_cx_total["Mass[Gev]"].loc[df_cx_total["GoM"] == g_o_m], df_cx_total["sigma[pb]_BOX_SChan_eta0"].loc[df_cx_total["GoM"] == g_o_m], marker=".", label='Box X resonance, G/M = %s' % str(g_o_m)) 

plt.xlabel("Mass [GeV]")
plt.ylabel("sigma(pp to hh)[pb]")
plt.legend()
plt.savefig('cross_section_interference.png')
plt.clf()

#################################################
# Draw full CX
#################################################
#def full_cx(box, triangle, nonres_int, resonant, res_box, res_trangle, kt, ktH, kap_hhh, kap_hhH):
#    # for a given mass