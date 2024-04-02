import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

#################################################
# Draw full CX in kappas and angles
# kl , sin_alpha , MX
#################################################
def cross_sections_parts(GoM, MX, sin_alpha, kl, table):

    cos_alpha = (1 - sin_alpha**2)**0.5
    SM_lam = 31.8
    mh = 125.0

    kt = cos_alpha
    ktH = sin_alpha
    kap_hhh = SM_lam * kl
    kap_hhH = SM_lam * (2*mh**2 + MX**2)*(cos_alpha**2*sin_alpha + (cos_alpha**3-kl)*cos_alpha/sin_alpha)/mh**2

    #NNLO values of LHC HH TWIKI (pb)
    xs_box_nnlo = 0.0703874
    xs_Sh_nnlo = 0.0110595
    xs_box_Sh_int_nnlo = -0.0504111

    #print(table)
    #print(table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)][["sigma[pb]_SChan_eta0", "sigma[pb]_BOX_SChan_eta0", "sigma[pb]_SChan_h_SChan_eta0"]])
    # taken from our table with all correction factors on (pb)
    resonant = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_SChan_eta0"].values[0] 
    res_box = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_BOX_SChan_eta0"].values[0]
    res_trangle = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_SChan_h_SChan_eta0"].values[0]

    nonresonant_XS = xs_box_nnlo * kt**4 + xs_Sh_nnlo * kt**2 * (kap_hhh)**2 + xs_box_Sh_int_nnlo * kt**3 * kap_hhh
    resonant_XS = resonant * ktH**2 * kap_hhH**2
    res_nonres_int_XS = kap_hhH * (res_box * kt**2 * ktH  + res_trangle * kap_hhh * ktH * kt )

    width_HH = (1-4*mh**2/MX**2)**0.5/(8*3.1415*MX)# kap_hhH**2 #*
    width_SM = sin_alpha**2*table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["GammaH_twiki[GeV]"].values[0]

    return {
        "nonresonant_XS" : nonresonant_XS,
        "resonant_XS" : resonant_XS,
        "res_nonres_int_XS" : res_nonres_int_XS,
        "full" : nonresonant_XS + resonant_XS + res_nonres_int_XS,
        "ratio" : (res_nonres_int_XS)/(nonresonant_XS + resonant_XS),
        "BH_HH" : width_HH/(width_HH + width_SM), ### this needs to multiply the resonant only part
        "width_HH" : width_HH,
        "width_SM" : width_SM,
        "total_width" : width_HH + width_SM,
        "k_hhH" : kap_hhH/SM_lam
    }

table = pd.read_csv("interferences_PPtoXtoHH_PPtoHH/cross_sections_for_mass_width.csv")

masses = np.unique(table["Mass[Gev]"].values) 
GoM = 0.001
print(masses)


# draw width parts X sin_alpha to kl = 1 , 10 to MX = 600 GeV
MX = 600
kl = 1  # from -15 , 15
sin_alpha_scan = [*range(1, 10, 1)] + [*range(10, 100, 10)]
width_HH_scan_kl1_MX_600 = [ cross_sections_parts(GoM, MX, sin_alpha_loc/100, kl, table)["total_width"] for sin_alpha_loc in sin_alpha_scan]
print(sin_alpha_scan)
print(width_HH_scan_kl1_MX_600)

columns = [
    'Mass[Gev]', 
    'sin_alpha',
    'kl', 
    'GoM',
    "BR_HH",
    "ratio",
    "ratio_nonres",
    "resonant_XS",
    "k_hhH"
    ] # do first to kl = 1


for gg, gom in enumerate([0.001, 0.01, 0.05, 0.1]):
    print("G/M = %s" % str(gom))
    df = pd.DataFrame(columns=columns)
    klam = [-5, 0, 1, 2.45, 5] #[*range(-5, 5, 1)]
    counter = 0
    for mass in masses:
        for kl in klam:
            for sin_alpha in sin_alpha_scan: 
                sin_alpha_loc = sin_alpha/100
                computations = cross_sections_parts(gom, mass, sin_alpha_loc, kl, table)
                df.loc[counter] = [
                    mass/1000,
                    sin_alpha_loc, 
                    kl,
                    computations["total_width"]/mass,
                    computations["BH_HH"],
                    computations["ratio"],
                    (computations["resonant_XS"]-computations["nonresonant_XS"])/(computations["resonant_XS"]+computations["nonresonant_XS"]),
                    computations["resonant_XS"],
                    computations["k_hhH"]
                ]
                counter += 1 

    print(df)
    #print(df_kl_1.loc[(df['Mass[Gev]'] == 600)])
    # GoM

    
    for kl in klam:
        variables = ['GoM', "BR_HH", "ratio", "ratio_nonres"] if gg == 0 and kl == 1 else ["ratio", "ratio_nonres",  "resonant_XS", "k_hhH"] 

        df_kl_1 = df.loc[(df['kl'] == kl)]
        # https://stackoverflow.com/questions/70714388/plot-a-2d-colormap-heatmap-in-matplotlib-with-x-y-z-data-from-a-pandas-dataframe
        #fig, ax = plt.subplots()
        fig = plt.figure(figsize=(9,9))
        ax = fig.add_subplot(111)
        for var in variables:
            df_pivoted = df_kl_1.pivot(columns='Mass[Gev]', index='sin_alpha', values=var)
            ax = sns.heatmap(data=df_pivoted, annot=False, fmt='f', cmap='RdYlGn', cbar=True, cbar_kws={'label': var}, square=True) #  annot=True,
            ax.tick_params(labelrotation=0)
            plt.title("kl = %s, raw GoM = %s" % (str(kl), str(gom)))
            #plt.tight_layout(pad=1)
            outfile = 'interferences_PPtoXtoHH_PPtoHH/scan_singlet_Z2/%s_scan_singlet_Z2_raw_GoM_%s_kl_%s.png' % (var, str(gom).replace(".", "o"), str(kl).replace(".", "o"))
            plt.savefig(outfile)
            print("Saved %s" % outfile)
            plt.clf()
            #########################################


