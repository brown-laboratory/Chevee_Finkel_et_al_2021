###############################################################################
# Generate_Figures_and_Results_CheveeFinkeletal2021.py
# The following files are necessary and available upon request:
#   - All_units_OptoMetrics_df.pkl
#   - master_log (11 pkl parts, A through K)
#   - S1master_log.pkl
#   - XCorr_df.pkl
#   - S1XCorr_df.pkl
#   - master_DREADD.pkl
# Note: You will need to change the paths when loading data depending on where 
# you downloaded it
###############################################################################
import numpy as np # '1.16.2'
import pandas as pd # '0.24.2'
import scipy as sp # '1.2.1'
from scipy import signal
from scipy import stats
from scipy.io import loadmat
import matplotlib.pyplot as plt # '3.0.3'
import seaborn as sns # '0.9.0'
import os
import sys
import math
from sklearn import metrics # '0.20.3'
from sklearn import linear_model
import pickle # '4.0'
import h5py # '2.10.0'
sys.path.insert(0,'C:\\Users\Brown Lab\Documents\Maxime_tools')

# %% Load master_log
from pathlib import Path
data_folder = Path('F:/Claustrum/Maxime_ClaustrumRevision/test')

file_to_open = data_folder / "master_log-partA.pkl"
with open(file_to_open, 'rb') as f:
   master1a = pickle.load(f, encoding='latin1')

file_to_open = data_folder / "master_log-partB.pkl"
with open(file_to_open, 'rb') as f:
   master1b = pickle.load(f, encoding='latin1')

file_to_open = data_folder / "master_log-partC.pkl"
with open(file_to_open, 'rb') as f:
   master1c = pickle.load(f, encoding='latin1')

file_to_open = data_folder / "master_log-partD.pkl"
with open(file_to_open, 'rb') as f:
   master1d = pickle.load(f, encoding='latin1')
   
file_to_open = data_folder / "master_log-partE.pkl"
with open(file_to_open, 'rb') as f:
   master2a = pickle.load(f, encoding='latin1')

file_to_open = data_folder / "master_log-partF.pkl"
with open(file_to_open, 'rb') as f:
   master2b = pickle.load(f, encoding='latin1')

file_to_open = data_folder / "master_log-partG.pkl"
with open(file_to_open, 'rb') as f:
   master2c = pickle.load(f, encoding='latin1')

file_to_open = data_folder / "master_log-partH.pkl"
with open(file_to_open, 'rb') as f:
   master2d = pickle.load(f, encoding='latin1')

file_to_open = data_folder / "master_log-partI.pkl"
with open(file_to_open, 'rb') as f:
   master3a = pickle.load(f, encoding='latin1')

file_to_open = data_folder / "master_log-partJ.pkl"
with open(file_to_open, 'rb') as f:
   master3b = pickle.load(f, encoding='latin1')

file_to_open = data_folder / "master_log-partK.pkl"
with open(file_to_open, 'rb') as f:
   master3c = pickle.load(f, encoding='latin1')   
   
master_log=pd.concat([master1a,master1b,master1c,master1d,master2a,master2b,master2c, master2d, master3a,master3b, master3c], axis=0)
del master1a,master1b,master1c,master1d,master2a,master2b,master2c, master2d, master3a, master3b, master3c

# %%
##############################################################################
# Figure 1
##############################################################################

# Figure 1A
# Diagram

# Figure 1B
# Geneated in Excel file EFtask performance.xlsx. Numbers were generated using Generate_Log_IntanData_forBehavior.m in MATLAB

# Figure 1C
from Figure1_functions_20210912 import Fraction_of_trials
# Behavior data:Fraction of trials
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37', 'Claustrum18', 'Claustrum23', 'Claustrum25'] 
HITs, HITsSEM, CRs, CRsSEM, MISSs, MISSsSEM, FAs, FAsSEM = Fraction_of_trials(mice, master_log)

# Figure 1D
# Diagram

# Figure 1E
# Images

# Figure 1F
# Images

# Figure 1G
data_folder = Path('C:\\Users\Brown Lab\Desktop\master_log')
file_to_open = data_folder / 'All_units_OptoMetrics_df.pkl'
with open(file_to_open, 'rb') as f:
  All_units_OptoMetrics_df = pickle.load(f, encoding='latin1')

from Figure1_functions_20210912 import Opto_plots
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37','Claustrum18', 'Claustrum23', 'Claustrum25'] 
Opto, Salt, OptoSameTT, SaltSameTT, InBetween = Opto_plots(mice, master_log, All_units_OptoMetrics_df) #Need to go inside the function to change which examples are plotted
# outputs number of units in each category
# 73 / 142 / 91 / 117 / 122 / 

# %%
##############################################################################
# Figure S1
############################################################################## 

# Figure S1A
# Images

# Figure S1B
# Generated in MATLAB using Make_plot_pvalues_SALT.m

# Figure S1C, D, E (Same function as Figure 1G)
data_folder = Path('C:\\Users\Brown Lab\Desktop\master_log')
file_to_open = data_folder / 'All_units_OptoMetrics_df.pkl'
with open(file_to_open, 'rb') as f:
  All_units_OptoMetrics_df = pickle.load(f, encoding='latin1')

from Figure1_functions_20210912 import Opto_plots
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37','Claustrum18', 'Claustrum23', 'Claustrum25'] 
Opto, Salt, OptoSameTT, SaltSameTT, InBetween = Opto_plots(mice, master_log, All_units_OptoMetrics_df) #Need to go inside the function to change which examples are plotted
# outputs number of units in each category
# 73 / 142 / 91 / 117 / 122 / 

# Figure S1F
# Diagram

# Figure S1G-N
# Images

# %%
##############################################################################
# Figure 2
##############################################################################

# Figure 2A
from Figure2_functions_20210912 import example_claustrum
example_name='Cl4_06-05-17_TT4clst1'
example_claustrum( unit_log=master_log.loc[np.equal(master_log['unit_name'],example_name)] , yaxis_range=[0,50] )

# Figure 2B
from Figure2_functions_20210912 import Heatmap_AllTrialTypes
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37'] #n=XXX neurons
Heatmap_AllTrialTypes(mice, master_log)

# Figure 2C-E
from Figure2_functions_20210912 import Responsiveness
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37','Claustrum18', 'Claustrum23', 'Claustrum25'] 
for trial_type in ['SomCR', 'VisCR']:
    Number_Resp, Number_Total= Responsiveness(trial_type, master_log, mice)
    print(Number_Resp)
    print(Number_Resp/Number_Total)

# %%
##############################################################################
# Figure S2
############################################################################## 

# Figure S2A
from FigureS2_functions_20210912 import example_claustrum_lickaligned
example_name='Cl4_06-05-17_TT4clst1'
example_claustrum_lickaligned( unit_log=master_log.loc[np.equal(master_log['unit_name'],example_name)] , yaxis_range=[0,50] )

# Figure S2B
from FigureS2_functions_20210912 import Heatmap_AllTrialTypes_lickaligned
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37'] 
Heatmap_AllTrialTypes_lickaligned(mice, master_log)

# Figure S2C
#####################################
#  need to update pandas from version '0.24.2' to version 1.1.3 to download and read S1master_log
#####################################
# Import S1 data frame 
from pathlib import Path
import pickle
data_folder = Path('D:/Claustrum/master_log/FINAL')
file_to_open = data_folder / "S1master_logC.pkl"
with open(file_to_open, 'rb') as f:
   S1master_log = pickle.load(f, encoding='latin1')   

from FigureS2_functions_20210912 import S1_CRresponse
mice=['EF0074', 'EF0076','EF0077','EF0079']
Mean_Hit, SEM_Hit, Mean_CR, SEM_CR= S1_CRresponse (S1master_log, mice)
# (13.282718095725109,
#  0.5827184254832429,
#  12.092587607553773,
#  0.5200816997742062)   
# Wilcoxon p=7.73e-5
# N=4 mice; n=754 neurons

# Figure S2D,E
opto_log=master_log[master_log['Category']=='OptoTag']
from Figure2_functions_20210912 import Heatmap_AllTrialTypes
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37'] #n=XXX neurons
Heatmap_AllTrialTypes(mice, opto_log)
mice=['Claustrum18', 'Claustrum23', 'Claustrum25'] #n=XXX neurons
Heatmap_AllTrialTypes(mice, opto_log)

# Figure S2F 
from FigureS2_functions_20210912 import Stim_Lick_aligned
trial_type='SomHit'
examples=['Cl4_05-23-17_TT7clst4',
'Cl5_05-25-17_TT8clst1',
'Cl5_06-04-17_TT6clst4',
'Cl5_06-09-17_TT2clst2',
'Claustrum23_20190724_TT6clst1',
'Claustrum23_20190724_TT6clst2',
'Claustrum23_20190725_TT3clst3',
'Claustrum25_20190808_TT1clst2',
'Claustrum25_20190812_TT5clst5',
'Claustrum37_20191125_TT8clst1']
examples=np.unique(master_log['unit_name'])
for example in examples:
    Stim_Lick_aligned(master_log,trial_type, example)

# %%
##############################################################################
# Figure 3
##############################################################################

# Figure 3A-D
from Figure2_functions_20210912 import Responsiveness
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37','Claustrum18', 'Claustrum23', 'Claustrum25'] 
for trial_type in ['SomHit', 'VisHit']:
    Number_Resp, Number_Total= Responsiveness(trial_type, master_log, mice)
    print(Number_Resp)
    print(Number_Resp/Number_Total)

# Figure 3E-H
from Figure3_functions_20210912 import LDA_figure
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37','Claustrum18', 'Claustrum23', 'Claustrum25'] 
meanSessionSize, minSessionSize,  maxSessionSize,  semSessionSize, meanScoreLick, semScoreLick, meanScoreStim, semScoreStim, meanScoreBlock, semScoreBlock, p_anova, p_lick_vs_stim, p_lick_vs_block, p_stim_vs_block = LDA_figure(master_log, mice)
#(7.465753424657534,
# 1,
# 20,
# 0.5063241232288771,
# 0.8141083694691659,
# 0.009847095085922236,
# 0.6025373358084632,
# 0.007423418864006554,
# 0.6519113144534906,
# 0.009196302816454602,
# 3.7953742291780915e-42,
# 1.3368529388047567e-13,
# 1.1063294062047036e-12,
# 3.927201982727906e-06)

# Figure 3I
# Diagram

# Figure 3J,K
from Figure3_functions_20210912 import AUC_OptoTag_Modality
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37','Claustrum18', 'Claustrum23', 'Claustrum25'] 
Touch_mean, Touch_sem, Visual_mean, Visual_sem = AUC_OptoTag_Modality(mice, master_log, BinNumber=6)
#(0.549608553985441,
# 0.01218985608498346,
# 0.5723787689067167,
# 0.010951122502857662)

# %%
##############################################################################
# Figure S3
##############################################################################

# Figure S3A
from FigureS3_functions_20210912 import Hit_vs_FA
MeanHit, SEMHit, MeanFA, SEMFA, p = Hit_vs_FA(master_log, mice)
#(9.372464255450652,
# 0.47748253288006265,
# 8.507185975275663,
# 0.4202464827433787,
# 5.648828541173874e-08)

# FIgure S3B
from FigureS3_functions_20210912 import TouchHit_vs_VisualHit_OPTO
MeanHitTouch, SEMHitTouch, MeanHitVis, SEMHitVis, p = TouchHit_vs_VisualHit_OPTO(master_log, mice)
#(19.648242919851857,
# 1.7247426889756312,
# 20.545703349893145,
# 1.8456247960209429,
# 0.12575225348464605)

# Figure S3C
from FigureS3_functions_20210912 import TouchHit_vs_VisualHit
MeanHitTouch, SEMHitTouch, MeanHitVis, SEMHitVis, p = TouchHit_vs_VisualHit(master_log, mice)
#(9.403198944953399,
# 0.4921036705851308,
# 9.356472674619248,
# 0.4810439750771753,
# 0.7942600025528961)

# Figure S3D
from FigureS3_functions_20210912 import meanAUC_category
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37']
Category='OptoTag' 
BinNumber=6
Som_mean, Som_sem, Vis_mean, Vis_sem, p = meanAUC_category(master_log, mice, Category, BinNumber)
#(0.5454028763881733,
# 0.02196788475346742,
# 0.5700347453961705,
# 0.019300787221593556,
# 0.7343252914391281)

# Figure S3E
mice=['Claustrum18', 'Claustrum23', 'Claustrum25']
Category='OptoTag' 
BinNumber=6
Som_mean, Som_sem, Vis_mean, Vis_sem, p = meanAUC_category(master_log, mice, Category, BinNumber)
#(0.5527127445929481,
# 0.01361792073940051,
# 0.574108881497834,
# 0.012616875211977059,
# 0.24303547341118348)

# Figure S3F
from FigureS3_functions_20210912 import All_AUC_modality
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37','Claustrum18', 'Claustrum23', 'Claustrum25'] 
BinNumber=6
Som_mean, Som_sem, Vis_mean, Vis_sem,p, X, X, X = All_AUC_modality(master_log, mice, BinNumber)

# Figure S3G
from FigureS3_functions_20210912 import Onset_cumfreq
Som_mean, Som_sem, Vis_mean, Vis_sem, x, y, temp_units1, temp_units2, temp_median1, temp_median2, p = Onset_cumfreq(master_log, mice)

# Figure S3H
from FigureS3_functions_20210912 import Onset_scatter_bar
p = Onset_scatter_bar(master_log, mice,  x, y, temp_units1, temp_units2, temp_median1, temp_median2)

# %%
##############################################################################
# Figure 4
##############################################################################

# Figure 4A
# Diagram

# Figure 4B,C
from Figure4_functions_20210912 import AUC_All_LickDirection
example1=232
example2=258
BinNumber=6
Right_mean, Right_sem, Left_mean, Left_sem, p_all, p_sig = AUC_All_LickDirection(mice, master_log, BinNumber, example1, example2)
#(0.5383366135032738,
# 0.0034059694149615096,
# 0.5331840733668348,
# 0.003411665094204253,
# 0.06058894053473373,
# 0.06767012459376)

# Figure 4D
# Diagram

# Figure 4G
from Figure4_functions_20210912 import ITNB_scatter
peak_Right_mean, peak_Right_sem, peak_Left_mean, peak_Left_sem, testR, testL,  temp_units1, temp_units2, temp_units1_totallicks, temp_units2_totallicks, temp_unit_names = ITNB_scatter(mice, master_log, example1, example2)

# Figure 4E,F
# Need to run Figure 4G first
from Figure4_functions_20210912 import example_ITBNandAUC
for example in  [example1, example2]:
    example_ITBNandAUC(mice, master_log, example, temp_units1, temp_units2, temp_units1_totallicks, temp_units2_totallicks, temp_unit_names)

# Figure 4H
from Figure4_functions_20210912 import R2_correlation
R2_correlation(mice, master_log, testR, testL, example1, example2)

from Figure4_functions_20210912 import Loglikelihood_decoder
# 6 mice contribute (3 and 3, 16 and 15 each)
Score_mean, Score_sem, Control_mean, Control_sem, p0, p1, p2,Number_of_units = Loglikelihood_decoder(mice, master_log)
#(0.6891110881800185,
# 0.024343526208324456,
# 0.5029912613664099,
# 0.006850576235667066,
# 1.037423329720682e-05,
# 8.646299427329582e-06,
# 0.6203880213542354)

# %%
##############################################################################
# Figure S4
##############################################################################
#Images only

# %%
##############################################################################
# Figure 5
##############################################################################

##########################
#LOAD XCorr_df
##########################
from pathlib import Path
data_folder = Path('C:\\Users\Brown Lab\Desktop\master_log\FINAL')

file_to_open = data_folder / "XCorr_df_20200731.pkl"
with open(file_to_open, 'rb') as f:
   XCorr_df = pickle.load(f, encoding='latin1')

from pathlib import Path
data_folder = Path('C:\\Users\Brown Lab\Desktop\master_log')

file_to_open = data_folder / "S1XCorr_df_20200730.pkl"
with open(file_to_open, 'rb') as f:
   S1XCorr_df = pickle.load(f, encoding='latin1')

# Figure 5A
from Figure5_functions_20210912 import Example_XCorr_pair
pair=2091
Example_XCorr_pair(master_log, XCorr_df, pair)

# Figure 5B
from Figure5_functions_20210912 import XCorr_example
pairs=[11,23]
for pair in pairs:
    XCorr_example( XCorr_df,pair)

# Figure 5C
from Figure5_functions_20210912 import Claustrum_v_S1_XCorr
p, Sig, S1Sig, LickPref1_Sig, LickPref2_Sig, FRSig, S1FRSig = Claustrum_v_S1_XCorr(XCorr_df, S1XCorr_df)
#Chi square p: 0.00900365040481692

# Figure 5D
from Figure5_functions_20210912 import CumFreq_Claustrum_v_S1
Sig_mean, Sig_sem, S1Sig_mean, S1Sig_sem = CumFreq_Claustrum_v_S1(Sig, S1Sig)
#(0.0717712305784988,
# 0.01154616200165663,
# 0.06325272572791209,
# 0.00872277782933359)

# Figure 5E
from Figure5_functions_20210912 import Quandrant_plot
All, FRAll, LickPref1_All, LickPref2_All = Quandrant_plot(XCorr_df, Sig, LickPref1_Sig, LickPref2_Sig, cutoff=0.1)

from Figure5_functions_20210912 import PiePlot_LickPref
p = PiePlot_LickPref(Sig, LickPref1_Sig, LickPref2_Sig, LickPref1_All, LickPref2_All)
#p=3.38e-6

# %%
##############################################################################
# Figure S5
##############################################################################

# Figure S5B,C
from FigureS5_functions_20210912 import ITNB_onset
mean_right, sem_right, mean_left, sem_left, temp_unit_names, temp_units1_totallicks, temp_units2_totallicks, temp_units1, temp_units2 = ITNB_onset(master_log, mice)

# Figure S5A
from FigureS5_functions_20210912 import S5_examples
directory='C:\\Users\\Brown Lab\\Desktop\\Figure 2020-07\\TO distribute in Figure folders\\'
S5_examples(master_log, mice, temp_units1, temp_units2, temp_unit_names,  temp_units1_totallicks, temp_units2_totallicks, directory)

# %%
##############################################################################
# Figure 6
##############################################################################

# Figure 6A
from Figure6_functions_20210912 import LickPrefCutoff
Right_pref, Left_pref = LickPrefCutoff(mice, master_log)

# Figure 6B
from Figure6_functions_20210912 import example_baseline
example_baseline(master_log, 'Claustrum31_20191107_TT7clst1') #'Claustrum31_20191107_TT7clst1'   'Claustrum31_20191114_TT3clst1' 
example_baseline(master_log, 'Claustrum31_20191114_TT3clst1' ) #'Claustrum31_20191107_TT7clst1'   'Claustrum31_20191114_TT3clst1' 

# Figure 6C
from Figure6_functions_20210912 import Intention_index
Right_mean, Right_sem, Left_mean, Left_sem, p1, p2, p3 = Intention_index(master_log, mice)

# Figure 6D
# Image

# Figure 6E
#####################################
#  need to update pandas from version '0.24.2' to version 1.1.3 to download and read S1master_log
#####################################file_to_open = 'F:\\Claustrum\DREADD\master_DREADD.pkl'
with open(file_to_open, 'rb') as f:
   master_DREADD = pickle.load(f, encoding='latin1')
   
mice=['ClaustrumO','ClaustrumP','ClaustrumQ','ClaustrumR','ClaustrumS','ClaustrumT','ClaustrumU','ClaustrumV','ClaustrumW','ClaustrumX', 'ClaustrumY']
from Figure6_functions_20210912 import DREADD_anova_wrapper
test_function='Hit_rate'
p_dreaddSALINE_dreaddAGONIST, p_controlSALINE_controlAGONIST, p_dreaddSALINE_controlSALINE,p_dreaddAGONIST_controlAGONIST = DREADD_anova_wrapper(mice, master_DREADD, test_function)   
# =============
# ANOVA SUMMARY
# =============

# Source          SS    DF1    DF2     MS      F    p-unc    np2      eps
# -----------  -----  -----  -----  -----  -----  -------  -----  -------
# virus        0.000      1      9  0.000  0.001    0.978  0.000  nan
# injection    0.000      1      9  0.000  0.448    0.520  0.047    1.000
# Interaction  0.002      1      9  0.002  1.500    0.252  0.143  nan
# (0.30247600783821565,
#  0.5902347260874201,
#  0.7522844631738326,
#  0.711648774599416)

# Figure 6F
test='FA_rate_impulsivity'
p_dreaddSALINE_dreaddAGONIST, p_controlSALINE_controlAGONIST, p_dreaddSALINE_controlSALINE,p_dreaddAGONIST_controlAGONIST = DREADD_anova_wrapper(mice, master_DREADD, test) 
# =============
# ANOVA SUMMARY
# =============

# Source          SS    DF1    DF2     MS       F    p-unc    np2      eps
# -----------  -----  -----  -----  -----  ------  -------  -----  -------
# virus        0.000      1      9  0.000   0.075    0.790  0.008  nan
# injection    0.002      1      9  0.002  10.146    0.011  0.530    1.000
# Interaction  0.001      1      9  0.001   6.047    0.036  0.402  nan
# (0.004398341775302218,
#  0.7946945523578082,
#  0.40165355383593326,
#  0.7436523609367314)

# Figure 6G
test_function='FA_rate_rule_exclusive'
p_dreaddSALINE_dreaddAGONIST, p_controlSALINE_controlAGONIST, p_dreaddSALINE_controlSALINE,p_dreaddAGONIST_controlAGONIST = DREADD_anova_wrapper(mice, master_DREADD, test_function)   
# =============
# ANOVA SUMMARY
# =============

# Source          SS    DF1    DF2     MS      F    p-unc    np2      eps
# -----------  -----  -----  -----  -----  -----  -------  -----  -------
# virus        0.000      1      9  0.000  0.482    0.505  0.051  nan
# injection    0.000      1      9  0.000  0.916    0.363  0.092    1.000
# Interaction  0.000      1      9  0.000  0.000    0.992  0.000  nan
   
# (0.5527097665840065,
#  0.48890215081316185,
#  0.6053260662882047,
#  0.4678664639038731)

# %%
##############################################################################
# Figure S6
##############################################################################

from FigureS6_functions_20210912 import Claustrum_v_S1_XCorr_poststim
Sig_mean, Sig_sem, S1Sig_mean, S1Sig_sem, p = Claustrum_v_S1_XCorr_poststim(XCorr_df, S1XCorr_df)

# Figure S6C-H
from FigureS6_functions_20210912 import Compare_XCorrWidth_Claustrum_vs_S1
mean_w, sem_w, mean_S1w, sem_S1w, p = Compare_XCorrWidth_Claustrum_vs_S1(XCorr_df, S1XCorr_df)
   
# %%
##############################################################################
# Figure S7
##############################################################################
   
# Figure S7A,B
from FigureS7_functions_20210912 import clustering
clustering(mice, master_log)

# Figure S7C
from FigureS7_functions_20210912 import Clusters_bargraph
Num_per_group, Mice_in_each_group = Clusters_bargraph(mice, master_log)
   
# %%
##############################################################################
# Figure Reviewers
##############################################################################   

# Figure 1: Geometrical relationship
from FigureReviewers_functions_20210912 import Tetrode_map_LickPref
Tetrode_map_LickPref(mice, master_log)
   
# Figure 2: does intention index figure depend on outlyers?
from FigureReviewers_functions_20210912 import Intention_noOutlyers
Right_mean, Right_sem, Left_mean, Left_sem, p1, p2, p3 = Intention_noOutlyers(mice, master_log)
   
# %%
##############################################################################
# Number reported in text
##############################################################################
   
# Number of responsive OptoTag neurons  
##############################################################################
from Figure2_functions_20210912 import Responsiveness
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37','Claustrum18', 'Claustrum23', 'Claustrum25'] 
Opto_log=master_log[master_log['Category']=='OptoTag']
for trial_type in ['SomCR', 'VisCR']:
    Number_Resp, Number_Total= Responsiveness(trial_type, Opto_log, mice)
    print(Number_Resp)
    print(Number_Resp/Number_Total)
#SomCR: 4/73=0.0548
#VisCR:2/73=0.0274
   
# Number of responsive S1 neurons
##############################################################################
from TextResults_functions_20210912 import Responsiveness_S1
trial_type='VisCR'
mice=['EF0074', 'EF0076','EF0077','EF0079']
numer_resp, numer_total=Responsiveness_S1(trial_type, S1master_log, mice) #Function from Figure1_functions
#Results: 143/754 (19.0%)

trial_type='SomCR'
mice=['EF0074', 'EF0076','EF0077','EF0079']
numer_resp, numer_total=Responsiveness_S1(trial_type, S1master_log, mice) #Function from Figure1_functions
#Results: 4/754 (0.5%)  
   
# Number of responsive S1 neurons
##############################################################################
from TextResults_functions_20210912 import SomCR_vs_VisCR 
Num = SomCR_vs_VisCR(master_log, mice, BinNumber=6)   
#Num=34 
  
# Number of licks produced during Hit trials and FA trials
##############################################################################
from TextResults_functions_20210912 import Num_of_licks
meanHit, semHit, meanFA, semFA, nHit, nFA, nFA, p = Num_of_licks(mice, master_log)
# meanHit:7.20 ,semHit:0.0233 , meanFA:6.092 , semFA:0.0389 , nHit=71454 trials, nFA=24521 trials, p=3.0e-244    

# Correlation between number of licks produced and FR
##############################################################################
from TextResults_functions_20210912 import Corr_lick_FR_AND_HitvsFApreLick
Hitmean_exc, Hitmsem_exc,FAmean_exc,FAmsem_exc,Hitmean_inh, Hitmsem_inh, FAmean_inh, FAmsem_inh = Corr_lick_FR_AND_HitvsFApreLick(master_log, mice)
#  75(60 activated (25pos_corr, 35 neg_corr) +15 inhibited (7pos_cor, 8neg_cor)) out of 347 responsive (21.6%)
#For inhibited neurons: p=4.576e-4, meanHit=4.645 ,semHit=0.527 , meanFA=5.387 ,semFA=0.585 , n=87
    #For activated neurons: p=2.550e-15, meanHit=15.50 ,semHit=0.926 , meanFA=13.328 ,semFA=0.8246 , n=260
   
 # Split different groups to test AUC Touch and Vis
##############################################################################   
from FigureS3_functions_20210912 import All_AUC_modality
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37','Claustrum18', 'Claustrum23', 'Claustrum25'] 
BinNumber=6
X,X,X,X,X,X, Num_SIG, Num_EXC, NUM_INH = All_AUC_modality(master_log, mice, BinNumber)    
#    ([0.561513715225209,
#  0.007943313424492816,
#  0.5656273376017031,
#  0.007627456823079985,
#  0.7755736928815845],
# [0.6261701853005067,
#  0.0044910256402111055,
#  0.621197818369402,
#  0.004856299429509752,
#  0.11906036605960235],
# [0.4329782052372923,
#  0.00451772286917138,
#  0.4336270357206496,
#  0.004853021954263879,
#  0.6259107217420403])

 # R2 for optotag only
############################################################################## 
# Need to have generated testR and testL from ITNB_scatter (Figure 4)
from TextResults_functions_20210912 import R2_optoOnly
r,p=R2_optoOnly(mice, master_log, testR, testL)
#(0.9135758802849533, 1.8575011822866597e-29)


 # Compared the proportions of correlated neurons across block type and response types
############################################################################## 
from TextResults_functions_20210912 import  XCorr_VisualBlock
count_Vis, total_count = XCorr_VisualBlock(XCorr_df)
from TextResults_functions_20210912 import  XCorr_TouchBlock
count_Som, total_count = XCorr_TouchBlock(XCorr_df)
from scipy.stats import chi2_contingency
#TOUCH: 55 pairs
#Vision: 58 pairs
g, p, dof, expctd=chi2_contingency([[total_count,count_Vis],[total_count,count_Som]])

from TextResults_functions_20210912 import  XCorr_CorrectTrial
count_Correct, total_count = XCorr_CorrectTrial(XCorr_df)
from TextResults_functions_20210912 import  XCorr_IncorrectTrial
count_Incorrect, total_count = XCorr_IncorrectTrial(XCorr_df)
#CORRECT: 65 pairs
#INCORRECT: 48 pairs
g, p, dof, expctd=chi2_contingency([[total_count,count_Correct],[total_count,count_Incorrect]])

 # Compared baseline FR between Claustrum and S1
############################################################################## 
from TextResults_functions_20210912 import FRpre_Claustrum_vs_S1
Mean, SEM, S1Mean, S1SEM, n, nS1, p = FRpre_Claustrum_vs_S1(XCorr_df, S1XCorr_df)
#Results: meanFR:8.75 , semFR:0.674 , MeanFRS1:12.00 , semFRS1:0.639 , n=75 , nS1=129, p=0.00118

 # One more XCorr
##############################################################################
from TextResults_functions_20210912 import XCorr_strength_byLickSide
meanRIGHT, semRIGHT, meanLEFT, semLEFT, p = XCorr_strength_byLickSide(XCorr_df)
#(0.10510813249912059,
# 0.02867523223876897,
# 0.09668915499613566,
# 0.029266555881769702,
# 0.13865474831098074)


 # ILI rate DREADD
##############################################################################
#####################################
#  need to update pandas from version '0.24.2' to version 1.1.3 to download and read S1master_log
#####################################file_to_open = 'F:\\Claustrum\DREADD\master_DREADD.pkl'
with open(file_to_open, 'rb') as f:
   master_DREADD = pickle.load(f, encoding='latin1')
   
mice=['ClaustrumO','ClaustrumP','ClaustrumQ','ClaustrumR','ClaustrumS','ClaustrumT','ClaustrumU','ClaustrumV','ClaustrumW','ClaustrumX', 'ClaustrumY']
from TextResults_functions_20210912 import ILIrate_DREADD
p_dreaddSALINE_dreaddAGONIST, p_controlSALINE_controlAGONIST, p_dreaddSALINE_controlSALINE,p_dreaddAGONIST_controlAGONIST = ILIrate_DREADD(mice, master_DREADD)   
# =============
# ANOVA SUMMARY
# =============

# Source          SS    DF1    DF2     MS      F    p-unc    np2      eps
# -----------  -----  -----  -----  -----  -----  -------  -----  -------
# virus        0.018      1      9  0.018  0.947    0.356  0.095  nan
# injection    0.001      1      9  0.001  0.337    0.576  0.036    1.000
# Interaction  0.004      1      9  0.004  1.753    0.218  0.163  nan

# (0.25253535756073153,
#  0.3971490277008885,
#  0.21524824961691802,
#  0.6215237166067376)  
   
