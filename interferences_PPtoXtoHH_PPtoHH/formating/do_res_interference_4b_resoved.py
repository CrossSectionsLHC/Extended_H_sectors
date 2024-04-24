from __future__ import absolute_import
from __future__ import print_function
import CombineHarvester.CombineTools.ch as ch
#print("CombineHarvester", ch.__version__)
#import chtools.utils as utils
#import ROOT
#import time
#import datetime
#import sys
#import argparse
import re
#import pprint
import os
#import gzip,shutil
#import numpy as np
import pandas as pd
#from collections import defaultdict
#import six
#import glob 

# slc7
#  cd /afs/cern.ch/work/a/acarvalh/HHH/CMSSW_11_3_4/src/ ; cmsenv ; cd -
# It will  duplicate and scale the resonant process

#datacards_folder_bbgg = "/eos/cms/store/group/phys_b2g/XToHH_HY_combination/bbgg/spin0_HH/v1/"
#datacards_folder_bbgg = "/afs/cern.ch/work/a/acarvalh/interferences_HH/tests_card/"
datacards_folder_bbgg = "/eos/cms/store/group/phys_b2g/XToHH_HY_combination/bbbb/spin0_HH/v2/"
br_HH_to_factor = 0.339
decay = "hbbhbb"
masses = [400] # [260, 270, 280, 300, 320, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000]
outfolder = "/afs/cern.ch/work/a/acarvalh/interferences_HH/tests_card/"

file_xsec = "Extended_H_sectors/interferences_PPtoXtoHH_PPtoHH/cross_sections_for_mass_width.csv"
table_res_xsec = pd.read_csv(file_xsec)
era = ["2016", "2017", "2018"]

MX = 400
GoM = 0.001 
resonant     = table_res_xsec.loc[(table_res_xsec["GoM"] == GoM) & (table_res_xsec["Mass[Gev]"] == MX)]["sigma[pb]_SChan_eta0_no_HH_decay"].values[0]  #BR_HH === to be treated on the physics model  
res_box      = table_res_xsec.loc[(table_res_xsec["GoM"] == GoM) & (table_res_xsec["Mass[Gev]"] == MX)]["sigma[pb]_BOX_SChan_eta0"].values[0]
res_triangle = table_res_xsec.loc[(table_res_xsec["GoM"] == GoM) & (table_res_xsec["Mass[Gev]"] == MX)]["sigma[pb]_SChan_h_SChan_eta0"].values[0]

print("XSsec to res, res-box, res-triangle", resonant, res_box, res_triangle)

duplicate = {
    "sig_NMSSM_bbbb_MX_400_MY_125" : {
    #"Radionhh500_2018" : {
        #"ggHH_kl_0_kt_1_kX_0_ktX_0_restriangle_0_resbox_0_nonresonly_1" : 0.069725,
        "ggHH_kl_1_kt_1_kX_0_ktX_0_restriangle_0_resbox_0_nonresonly_1" : 0.031047*br_HH_to_factor,
        "ggHH_kl_2p45_kt_1_kX_0_ktX_0_restriangle_0_resbox_0_nonresonly_1" : 0.013124*br_HH_to_factor,
        "ggHH_kl_5_kt_1_kX_0_ktX_0_restriangle_0_resbox_0_nonresonly_1" : 0.091172*br_HH_to_factor,
        "ggHH_kl_0_kt_0_kX_1_ktX_1_restriangle_0_resbox_0_nonresonly_0" : resonant*0.01*br_HH_to_factor,
        "ggHH_kl_0_kt_1_kX_1_ktX_1_restriangle_0_resbox_1_nonresonly_0" : res_box*br_HH_to_factor,
        "ggHH_kl_1_kt_1_kX_1_ktX_1_restriangle_1_resbox_0_nonresonly_0" : res_triangle*br_HH_to_factor  
    }
}

"""
ggHH_kl_()_kt_()_kX_()_ktX_()_restriangle_()_resbox_()_nonresonly_()
add_ggf_sample(kl=0.0, kt=1.0, kX=0.0, ktX=0.0, resbox=0, restriangle=0, nonresonly=1, xs=0.069725, label="ggHH_kl_0_kt_1_kX_0_ktX_0_restriangle_0_resbox_0_nonresonly_1") 
add_ggf_sample(kl=1.0, kt=1.0, kX=0.0, ktX=0.0, resbox=0, restriangle=0, nonresonly=1, xs=0.031047, label="ggHH_kl_1_kt_1_kX_0_ktX_0_restriangle_0_resbox_0_nonresonly_1")
add_ggf_sample(kl=2.45, kt=1.0, kX=0.0, ktX=0.0, resbox=0, restriangle=0, nonresonly=1, xs=0.013124, label="ggHH_kl_2p45_kt_1_kX_0_ktX_0_restriangle_0_resbox_0_nonresonly_1") 
add_ggf_sample(kl=5.0, kt=1.0, kX=0.0, ktX=0.0, resbox=0, restriangle=0, nonresonly=1, xs=0.091172, label="ggHH_kl_5_kt_1_kX_0_ktX_0_restriangle_0_resbox_0_nonresonly_1") 

add_ggf_sample(kl=0.0, kt=0.0, kX=1.0, ktX=1.0, resbox=0, restriangle=0, nonresonly=0, xs=resonant, label="ggHH_kl_0_kt_0_kX_1_ktX_1_restriangle_0_resbox_0_nonresonly_0") 
add_ggf_sample(kl=0.0, kt=1.0, kX=1.0, ktX=1.0, resbox=1, restriangle=0, nonresonly=0, xs=res_box, label="ggHH_kl_0_kt_1_kX_1_ktX_1_restriangle_0_resbox_1_nonresonly_0") 
add_ggf_sample(kl=1.0, kt=1.0, kX=1.0, ktX=1.0, resbox=0, restriangle=1, nonresonly=0, xs=res_triangle, label="ggHH_kl_1_kt_1_kX_1_ktX_1_restriangle_1_resbox_0_nonresonly_0") 
"""
for mass in masses :
    card_to_import = "%s/datacard_mass_%s.txt" % (datacards_folder_bbgg, str(mass))
    datacard_out   = "%s/datacard_mass_%s.txt" % (outfolder, str(mass))

    cb = ch.CombineHarvester()
    print("Making fake coupling version of: %s" % card_to_import)

    cb.ParseDatacard(card_to_import, analysis='comb', channel="Hbbgg", mass='')
    print('>> Finish Parsing cards')
    cb.PrintObs()

    catList=cb.bin_set()

    for cat in catList:
        print("Processing cat : ",cat)
        procList = cb.cp().bin([cat]).process_set()
        for prc in procList:
            proc_= prc #"----" #utl.getC3D4ProcName(prc)
            print("Looking at ",prc,"  -> ",proc_ )
            # ["Radionhh260_2016"]
            procs_to_clone = ["sig_NMSSM_bbbb_MX_400_MY_125"]
            if prc in procs_to_clone:
                print("  > Got proc : ",prc," to be cloned", duplicate["sig_NMSSM_bbbb_MX_400_MY_125"])
                #k3,k4=utl.getKlK4(proc_)
                #src_proc_xSection =  utl.getXSec(k3,k4)
                for target_prc in duplicate[proc_]:
                    #print(target_prc)
                    target_proc = "%s_%s" % (target_prc, decay) 
                    #print( target_proc, duplicate["sig_NMSSM_bbbb_MX_500_MY_125"][target_prc])
                    #k3,k4      = utl.getKlK4(target_prc)
                    #targetXSec_scale = utl.getXSec(k3,k4)/src_proc_xSection
                    #lumiScale= float(utl.getLumiForProc(prc))*1e3
                    # The original card is in 1 pb
                    target_rate      = duplicate[proc_][target_prc] #targetXSec_scale*lumiScale
                    target_prc_name = prc.replace(proc_, target_proc)
                    print("     > Cloning ",prc," for ",target_prc,"  | Xsec  scale : ",target_rate," to ",target_prc_name)
                    ch.CloneProcsAndSysts(
                        cb.cp().bin([cat]).process([prc]),
                        cb,
                        #lambda x: x.set_process(target_prc_name )
                        lambda x: x.set_process(x.process().replace("sig_NMSSM_bbbb_MX_400_MY_125", target_prc_name))
                    )
                    #    #cb.deep().bin([cat]).process([prc]),

                    def updateRate(x):
                        if x.process()==target_prc_name:
                            rate_now = x.rate()
                            print("rate", rate_now, target_rate)
                            x.set_rate(target_rate*rate_now)
                    cb.cp().bin([cat]).ForEachProc( updateRate )
                    print("     >  Rates for ",target_prc_name,"  set as ",cb.cp().bin([cat]).process([target_prc_name]).GetRate())
            #
            #for procs in procs_to_clone:
            #    #cb.FilterAll(lambda x: re.match(procs,x.process())!=None )
            #cb.FilterAll(procs_to_clone)
            #cb.FilterAll(lambda x: x.process() != "sig_NMSSM_bbbb_MX_500_MY_125" )
            #def remove_proc(p):
            #    if p.process() in procs_to_clone:
            #        p.remove(True)
            #        print("Processes being removed : ",procs)
            #        #set(procList).remove(procs)
            #cb.ForEachProc(remove_proc)
            

                

            #cb.process(procs_to_clone,False)

    if outfolder:
        if not os.path.exists(outfolder):
            os.system('mkdir -p '+outfolder)
        ofname=card_to_import.replace('.txt','_4b_interferences.txt')
        cb.FilterAll(lambda x: re.match("sig_*",x.process())!=None )
        # Filter scale_sig_bbbb         rateParam  *          ggHH_kl_1_kt_1_kX_0_ktX_0_restriangle_0_resbox_0_nonresonly_1_hbbhbb 0.339
        #cb.FilterAll(lambda x: x.process() == "sig_NMSSM_bbbb_MX_500_MY_125" )
        cb_cat= cb.cp() #.bin([cat])
        ofname=outfolder+'/'+ofname.split('/')[-1]
        print("Exporting card ",ofname)
        cb_cat.WriteDatacard( ofname , ofname.replace('.txt',".root"))


#cb.PrintProcs();
