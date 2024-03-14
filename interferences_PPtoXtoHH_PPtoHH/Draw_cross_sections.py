import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

#################################################
# Draw full CX
#################################################
def cross_sections_parts(GoM, MX, kt, ktH, k_hhh, kap_hhH, table):
    #NNLO values of LHC HH TWIKI (pb)
    xs_box_nnlo = 0.0703874
    xs_Sh_nnlo = 0.0110595
    xs_box_Sh_int_nnlo = -0.0504111

    SM_kl = 31.8
    kap_hhh = k_hhh * SM_kl

    #print(table)
    #print(table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)][["sigma[pb]_SChan_eta0", "sigma[pb]_BOX_SChan_eta0", "sigma[pb]_SChan_h_SChan_eta0"]])
    # taken from our table with all correction factors on (pb)
    resonant = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_SChan_eta0"].values[0] 
    res_box = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_BOX_SChan_eta0"].values[0]
    res_trangle = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_SChan_h_SChan_eta0"].values[0]

    nonresonant_XS = xs_box_nnlo * kt**4 + xs_Sh_nnlo * kt**2 * (kap_hhh)**2 + xs_box_Sh_int_nnlo * kt**3 * kap_hhh
    resonant_XS = resonant * ktH**2 * kap_hhH**2
    res_nonres_int_XS = kap_hhH * (res_box * kt**2 * ktH  + res_trangle * kap_hhh * ktH * kt )

    return {
        "nonresonant_XS" : nonresonant_XS,
        "resonant_XS" : resonant_XS,
        "res_nonres_int_XS" : res_nonres_int_XS,
        "full" : nonresonant_XS + resonant_XS + res_nonres_int_XS
    }

print("ahhhhh")
kt = 1 
ktH = 1
k_hhh = 1
kap_hhH = 31.8
GoM = 0.001
MX = 300
table = pd.read_csv("interferences_PPtoXtoHH_PPtoHH/cross_sections_for_mass_width.csv")
print(cross_sections_parts(GoM, MX, kt, ktH, k_hhh, kap_hhH, table)) 

# draw a line scanning kap_hhH from 0 - 100
scan_kap_hhH = [*range(0, 40, 1)] + [*range(50, 110, 10)] 
print(scan_kap_hhH)
cx_kl_1_kt_1_full = [cross_sections_parts(GoM, MX, kt, ktH, k_hhh, kap_hhH_scan, table)["full"] for kap_hhH_scan in scan_kap_hhH]
cx_kl_1_kt_1_res = [cross_sections_parts(GoM, MX, kt, ktH, k_hhh, kap_hhH_scan, table)["resonant_XS"] for kap_hhH_scan in scan_kap_hhH]
cx_kl_1_kt_1_int = [cross_sections_parts(GoM, MX, kt, ktH, k_hhh, kap_hhH_scan, table)["res_nonres_int_XS"] for kap_hhH_scan in scan_kap_hhH]

plt.plot(scan_kap_hhH, cx_kl_1_kt_1_full, marker=".", label="full",linewidth=3.0) 
plt.plot(scan_kap_hhH, cx_kl_1_kt_1_res, marker=".", label="resonant-only",linewidth=3.0) 
plt.plot(scan_kap_hhH, cx_kl_1_kt_1_int, marker=".", label="interference",linewidth=3.0) 

plot_type = "kap_hhH_scan_others_SM"
plt.xlabel("kappa hhH")
plt.title("MX = %s GeV | G/M = %s | kl = kt = 1" % (str(MX), str(GoM)))
plt.ylabel("sigma ( pp to (X to) HH)[pb]")
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.savefig('interferences_PPtoXtoHH_PPtoHH/%s.png' % plot_type)
plt.clf()

k_hhh = 10

cx_kl_1_kt_1_full = [cross_sections_parts(GoM, MX, kt, ktH, k_hhh, kap_hhH_scan, table)["full"] for kap_hhH_scan in scan_kap_hhH]
cx_kl_1_kt_1_res = [cross_sections_parts(GoM, MX, kt, ktH, k_hhh, kap_hhH_scan, table)["resonant_XS"] for kap_hhH_scan in scan_kap_hhH]
cx_kl_1_kt_1_int = [cross_sections_parts(GoM, MX, kt, ktH, k_hhh, kap_hhH_scan, table)["res_nonres_int_XS"] for kap_hhH_scan in scan_kap_hhH]

plt.plot(scan_kap_hhH, cx_kl_1_kt_1_full, marker=".", label="full",linewidth=3.0) 
plt.plot(scan_kap_hhH, cx_kl_1_kt_1_res, marker=".", label="resonant-only",linewidth=3.0) 
plt.plot(scan_kap_hhH, cx_kl_1_kt_1_int, marker=".", label="interference",linewidth=3.0) 

plot_type = "kap_hhH_scan_kl_10_kt_1"
plt.xlabel("kappa hhH")
plt.title("MX = %s GeV | G/M = %s | kl = 10 | kt = 1" % (str(MX), str(GoM)))
plt.ylabel("sigma ( pp to (X to) HH)[pb]")
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.savefig('interferences_PPtoXtoHH_PPtoHH/%s.png' % plot_type)
plt.clf()


