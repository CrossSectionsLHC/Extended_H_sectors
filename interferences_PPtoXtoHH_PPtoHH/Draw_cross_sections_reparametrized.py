import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

#################################################
# Draw full CX in kappas and angles
# kl , sin_alpha , MH
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

    width_HH = kap_hhH**2*(1-4*mh**2/MX**2)**0.5/(8*3.1415*MX)
    width_sum = sin_alpha**2*table["GammaH_twiki[GeV]"]

    return {
        "nonresonant_XS" : nonresonant_XS,
        "resonant_XS" : resonant_XS,
        "res_nonres_int_XS" : res_nonres_int_XS,
        "full" : nonresonant_XS + resonant_XS + res_nonres_int_XS
        "BH_HH" : width_HH/(width_HH + width_sum)
        "width_HH" : width_HH,
        "width_sum" : width_sum
    }


