import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

#################################################
# Draw full CX in couplings
#################################################
def cross_sections_parts_original(GoM, MX, kt, ktH, k_hhh, k_hhH, table):
    #NNLO values of LHC HH TWIKI (pb)
    xs_box_nnlo = 0.0703874
    xs_Sh_nnlo = 0.0110595
    xs_box_Sh_int_nnlo = -0.0504111

    SM_kl = 31.8
    mh = 125.0
    kap_hhh = k_hhh #* SM_kl
    kap_hhH = k_hhH #* SM_kl

    #print(table)
    #print(table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)][["sigma[pb]_SChan_eta0", "sigma[pb]_BOX_SChan_eta0", "sigma[pb]_SChan_h_SChan_eta0"]])
    # taken from our table with all correction factors on (pb)
    resonant = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_SChan_eta0"].values[0] 
    #resonant = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_SChan_eta0_no_HH_decay"].values[0] 
    res_box = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_BOX_SChan_eta0"].values[0]
    res_trangle = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_SChan_h_SChan_eta0"].values[0]
    resonant_twiki = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["sigma[pb]_twiki"].values[0] 

    nonresonant_XS = xs_box_nnlo * kt**4 + xs_Sh_nnlo * kt**2 * (kap_hhh)**2 + xs_box_Sh_int_nnlo * kt**3 * kap_hhh
    #resonant_XS = resonant * ktH**2 * kap_hhH**2
    res_nonres_int_XS = kap_hhH * (res_box * kt**2 * ktH  + res_trangle * kap_hhh * ktH * kt )

    width_HH = kap_hhH**2*(1-4*mh**2/MX**2)**0.5/(8*3.1415*MX)
    width_SM = table.loc[(table["GoM"] == GoM) & (table["Mass[Gev]"] == MX)]["GammaH_twiki[GeV]"].values[0] # sin_alpha**2* # not changing ktop
    BR_HH = width_HH/(width_HH + width_SM)

    resonant_XS = (resonant * ktH**2 * kap_hhH**2) #*0.01 #*BR_HH

    return {
        "nonresonant_XS" : nonresonant_XS,
        "resonant_XS" : resonant_XS,
        "res_nonres_int_XS" : res_nonres_int_XS,
        "full" : nonresonant_XS + resonant_XS + res_nonres_int_XS,
        "ratio" : (res_nonres_int_XS)/(nonresonant_XS + resonant_XS),
        "ratio_res" : (resonant_XS)/(nonresonant_XS + resonant_XS),
        "BH_HH" : width_HH/(width_HH + width_SM),
        "width_HH" : width_HH,
        "width_SM" : width_SM,
        "GoM" : (width_HH + width_SM)/MX,
        "kfact" : resonant_twiki/resonant
    }


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
        "full" : nonresonant_XS + resonant_XS + res_nonres_int_XS,
        "BH_HH" : width_HH/(width_HH + width_sum),
        "width_HH" : width_HH,
        "width_sum" : width_sum
    }


print("ahhhhh")
kt = 1 
ktH = 1
k_hhh = 1
kap_hhH = 31.8
k_hhH = 1
GoM = 0.001
MX = 400
table = pd.read_csv("interferences_PPtoXtoHH_PPtoHH/cross_sections_for_mass_width.csv")
print(cross_sections_parts_original(GoM, MX, kt, ktH, k_hhh, kap_hhH, table)) 

# [*range(-1100, 200, 100)] + [*range(-110, -50, 10)] + [*range(50, 110, 10)]   + [*range(200, 1100, 100)]
# draw a line scanning kap_hhH from 0 - 100
#print(range(-40, 40))
#scan_kap_hhH =  [*range(-40, -20, 10)] + [*range(-21, 20, 1)] + [*range(20, 50, 10)] #+ [*range(0, 41, 1)] 
scan_kap_hhH =   [*range(-100, 100, 1)] 
scan_kap_hhH_real = [kap_hhH/100 for kap_hhH in scan_kap_hhH]

print(scan_kap_hhH)
cx_kl_1_kt_1_full = [cross_sections_parts_original(GoM, MX, kt, ktH, k_hhh, kap_hhH/100, table)["full"] for kap_hhH in scan_kap_hhH]
cx_kl_1_kt_1_res = [cross_sections_parts_original(GoM, MX, kt, ktH, k_hhh, kap_hhH/100, table)["resonant_XS"] for kap_hhH in scan_kap_hhH]
cx_kl_1_kt_1_int = [abs(cross_sections_parts_original(GoM, MX, kt, ktH, k_hhh, kap_hhH/100, table)["res_nonres_int_XS"]) for kap_hhH in scan_kap_hhH]
# (GoM, MX, kt, ktH, k_hhh, k_hhH, table)

plt.plot(scan_kap_hhH_real, cx_kl_1_kt_1_full, marker=".", label="full",linewidth=3.0) 
plt.plot(scan_kap_hhH_real, cx_kl_1_kt_1_res, marker=".", label="resonant-only",linewidth=3.0) 
plt.plot(scan_kap_hhH_real, cx_kl_1_kt_1_int, marker=".", label="|interference|",linewidth=3.0) 

plot_type = "kap_hhH_scan_others_SM"
plt.xlabel("kappa hhH")
plt.title("MX = %s GeV | G/M = %s | kl = kt = = ktX = 1" % (str(MX), str(GoM)))
plt.ylabel("sigma ( pp to (X to) HH)[pb]")
plt.yscale("log")
#plt.xscale("log")
plt.legend()
figsave = 'interferences_PPtoXtoHH_PPtoHH/%s.png' % plot_type
plt.savefig(figsave)
plt.clf()
print("saved", figsave)

###############

scan_kap_hhH =  [*range(-40, -20, 10)] + [*range(-10, 10, 1)] + [*range(20, 50, 10)] #+ [*range(0, 41, 1)] 
scan_kap_hhH_real = [kap_hhH/10 for kap_hhH in scan_kap_hhH]
print(scan_kap_hhH)
cx_kl_1_kt_1_full = [cross_sections_parts_original(GoM, MX, kt, ktH, k_hhh, kap_hhH/10, table)["BH_HH"] for kap_hhH in scan_kap_hhH]
plt.plot(scan_kap_hhH_real, cx_kl_1_kt_1_full, marker=".", label="BR(X -> HH)",linewidth=3.0) 

plot_type = "BR_kap_hhH_scan_others_SM"
plt.xlabel("kappa hhH")
plt.ylabel("BR(X -> HH)")
plt.title("MX = %s GeV | kl = kt = = ktX = 1" % (str(MX)))
plt.ylabel("BR(X -> HH)")
plt.yscale("log")
#plt.xscale("log")
plt.legend()
figsave = 'interferences_PPtoXtoHH_PPtoHH/%s.png' % plot_type
plt.savefig(figsave)
plt.clf()
print("saved", figsave)

cx_kl_1_kt_1_full = [100*cross_sections_parts_original(GoM, MX, kt, ktH, k_hhh, kap_hhH/10, table)["GoM"] for kap_hhH in scan_kap_hhH]
plt.plot(scan_kap_hhH, cx_kl_1_kt_1_full, marker=".", label="MX = 400 GeV",linewidth=3.0) 

#cx_kl_1_kt_1_full_1tev = [100*cross_sections_parts_original(GoM, 500, kt, ktH, k_hhh, kap_hhH/10, table)["GoM"] for kap_hhH in scan_kap_hhH]
#plt.plot(scan_kap_hhH, cx_kl_1_kt_1_full_1tev, marker=".", label="MX = 500 GeV",linewidth=3.0) 

plot_type = "GoM__kap_hhH_scan_others_SM"
plt.xlabel("kappa hhH")
plt.ylabel("\\Gamma/M (%)")
plt.title("MX = %s GeV | kl = kt = = ktX = 1" % (str(MX)))
#plt.ylabel("\\Gamma/M")
#plt.yscale("log")
#plt.xscale("log")
plt.legend()
figsave = 'interferences_PPtoXtoHH_PPtoHH/%s.png' % plot_type
plt.savefig(figsave)
plt.clf()
print("saved", figsave)

k_hhh = 10

cx_kl_1_kt_1_full = [cross_sections_parts_original(GoM, MX, kt, ktH, k_hhh, kap_hhH_scan, table)["full"] for kap_hhH_scan in scan_kap_hhH]
cx_kl_1_kt_1_res = [cross_sections_parts_original(GoM, MX, kt, ktH, k_hhh, kap_hhH_scan, table)["resonant_XS"] for kap_hhH_scan in scan_kap_hhH]
cx_kl_1_kt_1_int = [cross_sections_parts_original(GoM, MX, kt, ktH, k_hhh, kap_hhH_scan, table)["res_nonres_int_XS"] for kap_hhH_scan in scan_kap_hhH]

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


scan_kap_hhh =  [*range(-40, -20, 10)] + [*range(-10, 10, 1)] + [*range(20, 50, 10)] #+ [*range(0, 41, 1)] 
scan_kap_hhh_real = [kap_hhh/10 for kap_hhh in scan_kap_hhh]
k_hhH = 1
print(scan_kap_hhH)
cx_kl_1_kt_1_full = [cross_sections_parts_original(GoM, MX, kt, ktH, k_hhh_scan/10, k_hhH, table)["full"] for k_hhh_scan in scan_kap_hhh]
cx_kl_1_kt_1_res = [cross_sections_parts_original(GoM, MX, kt, ktH, k_hhh_scan/10, k_hhH, table)["resonant_XS"] for k_hhh_scan in scan_kap_hhh]
cx_kl_1_kt_1_int = [abs(cross_sections_parts_original(GoM, MX, kt, ktH, k_hhh_scan/10, k_hhH, table)["res_nonres_int_XS"]) for k_hhh_scan in scan_kap_hhh]

plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_full, marker=".", label="full",linewidth=3.0) 
plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_res, marker=".", label="resonant-only",linewidth=3.0) 
plt.plot(scan_kap_hhh_real, cx_kl_1_kt_1_int, marker=".", label="interference",linewidth=3.0) 

plot_type = "kap_hhh_scan_others_SM"
plt.xlabel("k#lambda")
plt.title("MX = %s GeV | kX = kt = = ktX = 1" % (str(MX)))
plt.ylabel("sigma ( pp to (X to) HH)[pb]")
plt.yscale("log")
#plt.xscale("log")
plt.legend()
figsave = 'interferences_PPtoXtoHH_PPtoHH/%s.png' % plot_type
plt.savefig(figsave)
plt.clf()
print("saved", figsave)

plot_type = "GoM_scan_kl_1_kap_hhH_31p8_kt_1"
GoM = [0.001, 0.01, 0.05, 0.1]

k_hhh = 1
kap_hhH = 31.8

cx_kl_1_kt_1_full = [cross_sections_parts_original(GoMs, MX, kt, ktH, k_hhh, k_hhH, table)["full"] for GoMs in GoM]
cx_kl_1_kt_1_res = [cross_sections_parts_original(GoMs, MX, kt, ktH, k_hhh, k_hhH, table)["resonant_XS"] for GoMs in GoM]
cx_kl_1_kt_1_int = [cross_sections_parts_original(GoMs, MX, kt, ktH, k_hhh, k_hhH, table)["res_nonres_int_XS"] for GoMs in GoM]

plt.plot(GoM, cx_kl_1_kt_1_full, marker=".", label="full",linewidth=3.0) 
plt.plot(GoM, cx_kl_1_kt_1_res, marker=".", label="resonant-only",linewidth=3.0) 
plt.plot(GoM, cx_kl_1_kt_1_int, marker=".", label="interference",linewidth=3.0) 

plt.xlabel("kappa hhH")
plt.title("MX = %s GeV | kap_hhH = 31.8 | kl = 1 | kt = 1" % (str(MX)))
plt.ylabel("sigma ( pp to (X to) HH)[pb]")
plt.xlabel("G/M")
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.savefig('interferences_PPtoXtoHH_PPtoHH/%s.png' % plot_type)
plt.clf()

plot_type = "kfact_scan_kl_1_kap_hhH_31p8_kt_1"
# kfactor with MX
m_to_to = [ 260,  270,  300,  320,  340,  350,  360,  400,  450,  500,  550,  600,  650,  700,  750,  800,  900, 1000]
for GoMs in GoM:
    cx_kl_1_kt_1_full = [cross_sections_parts_original(GoMs, mass, kt, ktH, k_hhh, kap_hhH, table)["kfact"] for mass in m_to_to]
    plt.plot(m_to_to, cx_kl_1_kt_1_full, marker=".", label="simulation G/M = %s" % str(GoMs),linewidth=3.0) 
plt.xlabel("MX (GeV)")
#plt.title("MX = %s GeV | kap_hhH = 31.8 | kl = 10 | kt = 1" % (str(MX)))
plt.ylabel("k-factor")
plt.xlabel("G/M")
#plt.yscale("log")
#plt.xscale("log")
plt.legend()
plt.savefig('interferences_PPtoXtoHH_PPtoHH/%s.png' % plot_type)
plt.clf()

# make kX ktX for MX = 400

plot_type = "ratio_scan_kl_1_kap_hhH_31p8_kt_1"
for ktX_val in [1, 0.5, 0.1 ]:
    cx_kl_1_kt_1_full = [100*cross_sections_parts_original(GoMs, MX, kt, ktX_val, k_hhh, k_hhH_val, table)["ratio"] for k_hhH_val in scan_kap_hhH]
    plt.plot(scan_kap_hhH, cx_kl_1_kt_1_full, marker=".", label="ktX = %s" % str(ktX_val),linewidth=3.0) 
plt.xlabel("kX")
plt.title("MX = %s GeV | kl = kt = 1" % (str(MX)))
plt.ylabel("interference/(res+nonres) (%)")
#plt.yscale("log")
#plt.xscale("log")
plt.legend()
plt.savefig('interferences_PPtoXtoHH_PPtoHH/%s.png' % plot_type)
plt.clf()