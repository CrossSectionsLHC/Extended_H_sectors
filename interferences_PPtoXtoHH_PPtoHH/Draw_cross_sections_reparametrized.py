import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
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
    vev = 246

    kt = cos_alpha
    ktH = sin_alpha
    kap_hhh =  kl # SM_lam *
    kap_hhH = (2*mh**2 + MX**2)*(cos_alpha**2*sin_alpha + (cos_alpha**3-kl)*cos_alpha/sin_alpha)/(2 * vev)

    #NNLO values of LHC HH TWIKI (pb)
    xs_box_nnlo = 0.0703874
    xs_Sh_nnlo = 0.0110595
    xs_box_Sh_int_nnlo = -0.0504111

    #print(table)
    #print(table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)][["sigma[pb]_SChan_eta0", "sigma[pb]_BOX_SChan_eta0", "sigma[pb]_SChan_h_SChan_eta0"]])
    # taken from our table with all correction factors on (pb)
    resonant = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_SChan_eta0"].values[0] 
    #resonant = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_SChan_eta0_no_HH_decay"].values[0] 
    res_box = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_BOX_SChan_eta0"].values[0]
    res_trangle = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_SChan_h_SChan_eta0"].values[0]

    nonresonant_XS = xs_box_nnlo * kt**4 + xs_Sh_nnlo * kt**2 * (kap_hhh)**2 + xs_box_Sh_int_nnlo * kt**3 * kap_hhh
    res_nonres_int_XS = kap_hhH * (res_box * kt**2 * ktH  + res_trangle * kap_hhh * ktH * kt )

    width_HH = kap_hhH**2*(1-4*mh**2/MX**2)**0.5/(8*3.1415*MX) #*
    width_SM = sin_alpha**2*table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["GammaH_twiki[GeV]"].values[0]
    BR_HH = width_HH/(width_HH + width_SM)

    #resonant_XS = (resonant * ktH**2 )*BR_HH
    resonant_XS = (resonant * ktH**2 * kap_hhH**2)
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
        "k_hhH" : kap_hhH #/SM_lam
    }

def cross_sections_parts_singlet_no_Z2(GoM, MX, sin_alpha, kl, khhH, table):

    cos_alpha = (1 - sin_alpha**2)**0.5
    SM_lam = 31.8
    mh = 125.0

    kt = cos_alpha
    ktH = sin_alpha
    kap_hhh = SM_lam * kl #* cos_alpha
    kap_hhH = SM_lam * khhH #* sin_alpha #(2*mh**2 + MX**2)*(cos_alpha**2*sin_alpha + (cos_alpha**3-kl)*cos_alpha/sin_alpha)/mh**2

    #NNLO values of LHC HH TWIKI (pb)
    xs_box_nnlo = 0.0703874
    xs_Sh_nnlo = 0.0110595
    xs_box_Sh_int_nnlo = -0.0504111

    #print(table)
    #print(table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)][["sigma[pb]_SChan_eta0", "sigma[pb]_BOX_SChan_eta0", "sigma[pb]_SChan_h_SChan_eta0"]])
    # taken from our table with all correction factors on (pb)
    resonant = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_SChan_eta0"].values[0] 
    #resonant = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_SChan_eta0_no_HH_decay"].values[0] 
    res_box = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_BOX_SChan_eta0"].values[0]
    res_trangle = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_SChan_h_SChan_eta0"].values[0]

    nonresonant_XS = xs_box_nnlo * kt**4 + xs_Sh_nnlo * kt**2 * (kap_hhh)**2 + xs_box_Sh_int_nnlo * kt**3 * kap_hhh
    res_nonres_int_XS = kap_hhH * (res_box * kt**2 * ktH  + res_trangle * kap_hhh * ktH * kt )

    width_HH = (1-4*mh**2/MX**2)**0.5/(8*3.1415*MX)*kap_hhH**2 #*
    width_SM = sin_alpha**2*table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["GammaH_twiki[GeV]"].values[0]
    BR_HH = width_HH/(width_HH + width_SM)

    resonant_XS = (resonant * ktH**2 * kap_hhH**2)#*BR_HH
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
        "k_hhH" : khhH # kap_hhH/SM_lam
    }


table = pd.read_csv("interferences_PPtoXtoHH_PPtoHH/cross_sections_for_mass_width.csv")

masses = np.unique(table["Mass[Gev]"].values) 
GoM = 0.001
print(masses)

z2_symmetry = True 
########################################
# Draw for singlet with Z2 symmetry
#########################################
sin_alpha_scan = [*range(1, 10, 1)] + [*range(10, 100, 10)]

if z2_symmetry:
    # draw width parts X sin_alpha to kl = 1 , 10 to MX = 600 GeV
    MX = 600
    kl = 0.58  # from -15 , 15
    width_HH_scan_kl1_MX_600 = [ cross_sections_parts(GoM, MX, sin_alpha_loc/100, kl, table)["total_width"] for sin_alpha_loc in sin_alpha_scan]
    print(sin_alpha_scan)
    print(width_HH_scan_kl1_MX_600)
    colors = ["r", "b", "m", "y", "c"]
    # draw k_hhH with kl to few masses
    #scan_kap_hhh =  [*range(-120, -20, 10)] + [*range(-10, 10, 1)] + [*range(20, 120, 10)] #+ [*range(0, 41, 1)] 
    scan_kap_hhh =  [*range(80, 120, 1)] 
    scan_kap_hhh_real = [kap_hhh/100 for kap_hhh in scan_kap_hhh]
    k_hhH = 1
    for ss, sinA in enumerate([0.2, 0.16, 0.1]) :
        cx_kl_1_kt_1_600 = [cross_sections_parts(GoM, 600, sinA, kls/100, table)["k_hhH"] for kls in scan_kap_hhh]
        plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".", label="sinA = %s" % sinA, linestyle="--",linewidth=3.0, color=colors[ss]) 

        #cx_kl_1_kt_1_600 = [cross_sections_parts(GoM, 300, 0.16, kls/100, table)["k_hhH"] for kls in scan_kap_hhh]
        #plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".", label="MX = 300, sinA = 0.16",linewidth=3.0) 

        cx_kl_1_kt_1_600 = [cross_sections_parts(GoM, 300, sinA, kls/100, table)["k_hhH"] for kls in scan_kap_hhh]
        plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".",linewidth=3.0, color=colors[ss]) 

        #cx_kl_1_kt_1_600 = [cross_sections_parts(GoM, 300, 0.1, kls/100, table)["k_hhH"] for kls in scan_kap_hhh]
        #plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".", label="MX = 300, sinA = 0.1",linewidth=3.0) 


    plot_type = "kap_hhH_kl_scan"
    plt.xlabel("kl")
    #plt.title("sinA = 0.17")
    plt.ylabel("kX")
    #plt.yscale("log")
    #plt.xscale("log")
    plt.legend()
    figsave = 'interferences_PPtoXtoHH_PPtoHH/scan_singlet_Z2/%s.png' % plot_type
    plt.savefig(figsave)
    plt.clf()
    print("saved", figsave)

    # draw k_hhH with kl to few masses
    scan_kap_hhh =  [*range(90, 110, 1)] #+ [*range(20, 220, 10)] #+ [*range(0, 41, 1)] 
    scan_kap_hhh_real = [kap_hhh/100 for kap_hhh in scan_kap_hhh]

    k_hhH = 1
    MX = 300
    for ss, sinA in enumerate([0.3, 0.2, 0.1, 0.01]) :
        cx_kl_1_kt_1_600 = [cross_sections_parts(GoM, MX, sinA, kls/100, table)["resonant_XS"] for kls in scan_kap_hhh]
        plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".", label="sinA = %s" % (sinA), linestyle="--",linewidth=3.0, color=colors[ss]) 

        #for sinA in [0.2, 0.16, 0.1, 0.05, 0.01] :
        #    cx_kl_1_kt_1_600 = [cross_sections_parts(GoM, 300, sinA, kls/100, table)["nonresonant_XS"] for kls in scan_kap_hhh]
        #    plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".", label="MX = 300, sinA = %s" % sinA,linewidth=3.0, linestyle="--") 

        # for ss, sinA in [0.2, 0.16, 0.1, 0.05, 0.01] :
        cx_kl_1_kt_1_600 = [cross_sections_parts(GoM, MX, sinA, kls/100, table)["full"] for kls in scan_kap_hhh]
        plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".", linewidth=3.0,  color=colors[ss]) # label="MX = 300, sinA = %s" % sinA,

        cx_kl_1_kt_1_600 = [abs(cross_sections_parts(GoM, MX, sinA, kls/100, table)["res_nonres_int_XS"]) for kls in scan_kap_hhh]
        plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".", linewidth=2.0, linestyle=":", color=colors[ss]) # label="MX = 300, sinA = %s" % sinA,


    #cx_kl_1_kt_1_600 = [cross_sections_parts(GoM, 600, 0.1, kls/100, table)["resonant_XS"] for kls in scan_kap_hhh]
    #plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".", label="MX = 600, sinA = 0.11",linewidth=3.0) 

    #cx_kl_1_kt_1_600 = [cross_sections_parts(GoM, 300, 0.1, kls/100, table)["resonant_XS"] for kls in scan_kap_hhh]
    #plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".", label="MX = 300, sinA = 0.11",linewidth=3.0) 

    #cx_kl_1_kt_1_600 = [cross_sections_parts(GoM, 600, 0.05, kls/100, table)["resonant_XS"] for kls in scan_kap_hhh]
    #plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".", label="MX = 600, sinA = 0.05",linewidth=3.0) 

    #cx_kl_1_kt_1_600 = [cross_sections_parts(GoM, 300, 0.05, kls/100, table)["resonant_XS"] for kls in scan_kap_hhh]
    #plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".", label="MX = 300, sinA = 0.05",linewidth=3.0) 

    #cx_kl_1_kt_1_600 = [cross_sections_parts(GoM, 600, 0.01, kls/100, table)["resonant_XS"] for kls in scan_kap_hhh]
    #plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".", label="MX = 600, sinA = 0.01",linewidth=3.0) 

    #cx_kl_1_kt_1_600 = [cross_sections_parts(GoM, 300, 0.01, kls/100, table)["resonant_XS"] for kls in scan_kap_hhh]
    #plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_600, marker=".", label="MX = 300, sinA = 0.01",linewidth=3.0) 

    plot_type = "resXS_kl_scan"
    plt.xlabel("kl")
    plt.title("MX = %s GeV" % MX)
    plt.ylabel("sigma (pp > HH) [pb]")
    plt.yscale("log")
    #plt.xscale("log")
    plt.legend()
    figsave = 'interferences_PPtoXtoHH_PPtoHH/scan_singlet_Z2/%s.png' % plot_type
    plt.savefig(figsave)
    plt.clf()
    print("saved", figsave)

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
        klam = [-5, -0.58, 0, 0.58, 1, 2.45, 5] #[*range(-5, 5, 1)]
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
            #variables = ['GoM', "BR_HH", "ratio", "ratio_nonres"] if gg == 0 and kl == 1 else ["ratio", "ratio_nonres",  "resonant_XS", "k_hhH"] 
            variables = ['GoM', "BR_HH", "ratio", "ratio_nonres",  "resonant_XS", "k_hhH"] 

            df_kl_1 = df.loc[(df['kl'] == kl)]
            # https://stackoverflow.com/questions/70714388/plot-a-2d-colormap-heatmap-in-matplotlib-with-x-y-z-data-from-a-pandas-dataframe
            #fig, ax = plt.subplots()
            fig = plt.figure(figsize=(9,9))
            ax = fig.add_subplot(111)
            for var in variables:
                df_pivoted = df_kl_1.pivot(columns='Mass[Gev]', index='sin_alpha', values=var)
                if var in ["resonant_XS", 'GoM', "BR_HH"] :
                    ax = sns.heatmap(data=df_pivoted, annot=False, fmt='f', cmap='RdYlGn', cbar=True, cbar_kws={'label': var}, square=True, norm=LogNorm()) #  annot=True,
                else:
                    ax = sns.heatmap(data=df_pivoted, annot=False, fmt='f', cmap='RdYlGn', cbar=True, cbar_kws={'label': var}, square=True) 
                ax.tick_params(labelrotation=0)
                plt.title("kl = %s, raw GoM = %s" % (str(kl), str(gom)))
                #plt.tight_layout(pad=1)
                outfile = 'interferences_PPtoXtoHH_PPtoHH/scan_singlet_Z2/%s_scan_singlet_Z2_raw_GoM_%s_kl_%s.png' % (var, str(gom).replace(".", "o"), str(kl).replace(".", "o"))
                plt.savefig(outfile)
                

                print("Saved %s" % outfile)
                plt.clf()
                #########################################
else:
    # plot the same plots as the DIB review
    m_to_to = [270, 400, 500, 600, 700, 800]
    # [ 260  270  300  320  340  350  360  400  450  500  550  600  650  700  750  800  900 1000]
    lam_HHX_scan = [*range(-20, 20, 2)]
    kl = 1

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
        klam = [1] #[*range(-5, 5, 1)]
        counter = 0
        for lam_HHX_value in lam_HHX_scan:
            for mass in m_to_to:
                for kl in klam:
                    for sin_alpha in sin_alpha_scan: 
                        sin_alpha_loc = sin_alpha/100
                        computations = cross_sections_parts_singlet_no_Z2(gom, mass, sin_alpha_loc, kl, lam_HHX_value, table)
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
        for mass in m_to_to:
            # draw lam_HHX X sinA
            variables = ['GoM', "BR_HH", "ratio", "ratio_nonres"] if gg == 0 and kl == 1 else ["ratio", "ratio_nonres",  "resonant_XS"] 

            df_kl_1 = df.loc[(df['Mass[Gev]'] == mass/1000)]
            print(df_kl_1)
            # https://stackoverflow.com/questions/70714388/plot-a-2d-colormap-heatmap-in-matplotlib-with-x-y-z-data-from-a-pandas-dataframe
            #fig, ax = plt.subplots()
            fig = plt.figure(figsize=(9,9))
            ax = fig.add_subplot(111)
            for var in variables:
                df_pivoted = df_kl_1.pivot(columns='sin_alpha', index="k_hhH", values=var)
                if var in ["resonant_XS", 'GoM', "BR_HH"] :
                    ax = sns.heatmap(data=df_pivoted, annot=False, fmt='f', cmap='RdYlGn', cbar=True, cbar_kws={'label': var}, square=True, norm=LogNorm()) #  annot=True,
                else:
                    ax = sns.heatmap(data=df_pivoted, annot=False, fmt='f', cmap='RdYlGn', cbar=True, cbar_kws={'label': var}, square=True) 
                ax.tick_params(labelrotation=0)
                plt.title("kl = %s, MX = %s GeV, raw GoM = %s" % (str(kl), mass, str(gom)))
                #plt.tight_layout(pad=1)
                outfile = 'interferences_PPtoXtoHH_PPtoHH/scan_singlet_noZ2/%s_scan_singlet_noZ2_raw_GoM_%s_mass_%s.png' % (var, str(gom).replace(".", "o"), str(mass).replace(".", "o"))
                plt.savefig(outfile)
                

                print("Saved %s" % outfile)
                plt.clf()            
