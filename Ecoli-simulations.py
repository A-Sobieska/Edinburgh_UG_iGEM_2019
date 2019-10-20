"""
Flux balance analysis on Escherichia coli

Aleksandra Sobieska
University of Edinburgh Undergraduate iGEM team 2019
"""
from cobra import Reaction, Metabolite
import cobra.test
from cobra.flux_analysis import flux_variability_analysis
import matplotlib.pyplot as plt
#%% load iJO1366 reconstruction of E. coli from the imported library
WTmodel = cobra.test.create_test_model("ecoli") 
#close glucose uptake 
WTmodel.reactions.get_by_id("EX_glc_e").lower_bound = 0.
#set E. coli's objective to biomass growth rate
WTmodel.objective = 'Ec_biomass_iJO1366_WT_53p95M'

#list of considered carbon sources and their respective names
carbons = ["EX_glc_e", "EX_gal_e", "EX_lac__D_e"] 
carbon_names = ["glucose", "galactose", "lactate"]

#list of wild type and other gentically modified E. coli strains
strains = ["WT", "SHI", "HydA", "FDH"]
#%% simulating anaerobic conditions

#close oxygen uptake and production
WTmodel.reactions.get_by_id("EX_o2_e").lower_bound = 0.
WTmodel.reactions.get_by_id("EX_o2_e").upper_bound = 0.
#implement metabolic modifications according to
#Seppälä JJ et al. Prospecting hydrogen production of Escherichia coli by metabolic network modeling. 
#International Journal of Hydrogen Energy. 2013;28(27):11780-11789. doi:10.1016/j.ijhydene.2013.07.002.

#open the reaction of formate hydrogen lyase
WTmodel.reactions.get_by_id("FHL").lower_bound = -1000.0
WTmodel.reactions.get_by_id("FHL").upper_bound = 1000.0 
#deactivate pyruvate dehydrogenase
WTmodel.reactions.get_by_id("PDH").lower_bound = 0.0 
WTmodel.reactions.get_by_id("PDH").upper_bound = 0.0 
#deactivate 2-oxoglutarate dehydrogenase 
WTmodel.reactions.get_by_id("AKGDH").lower_bound = 0.0 
WTmodel.reactions.get_by_id("AKGDH").upper_bound = 0.0 
#make isocitrate dehydrogenase irreversible
WTmodel.reactions.get_by_id("ICDHyr").lower_bound = 0.0
WTmodel.reactions.get_by_id("ICDHyr").upper_bound = 1000.0
#%% create SHI strain
#make a copy of WT reconstruction as the metabolic basis for SHI strain
SHI_model = WTmodel.copy()

SHI_rxn = Reaction('SHI')
SHI_rxn.name = 'SHI hydrogenase'
SHI_rxn.subsystem = 'Hydrogen production'
SHI_rxn.lower_bound = -1000.0
SHI_rxn.upper_bound = 1000.0

#retrieve these metabolites from the reconstruction
NADPH_c = SHI_model.metabolites.get_by_id("nadph_c")
NADP_c = SHI_model.metabolites.get_by_id("nadp_c")
H_c = SHI_model.metabolites.get_by_id("h_c")
H2_c = SHI_model.metabolites.get_by_id("h2_c")

#make the reaction
SHI_rxn.add_metabolites({
    NADPH_c: -1.0,
    H_c: -1.0,
    NADP_c: 1.0,
    H2_c: 1.0
})

#add the reaction to the reconstruction
SHI_model.add_reactions([SHI_rxn])
#%% make HydA strain - introducing ferredoxin
HydA_model = WTmodel.copy()

FDX_rxn = Reaction('FDX')
FDX_rxn.name = 'NADH-ferredoxin reductase'
FDX_rxn.subsystem = 'Hydrogen production'
FDX_rxn.lower_bound = -1000.0
FDX_rxn.upper_bound = 1000.0

FNR_rxn = Reaction('FNR')
FNR_rxn.name = 'NADPH-ferredoxin reductase'
FNR_rxn.subsystem = 'Hydrogen production'
FNR_rxn.lower_bound = -1000.0
FNR_rxn.upper_bound = 1000.0


HydA_rxn = Reaction('HydA')
HydA_rxn.name = 'FeFe hydrogenase'
HydA_rxn.subsystem = 'Hydrogen production'
HydA_rxn.lower_bound = -1000.0
HydA_rxn.upper_bound = 1000.0

#create a new metabolite from scratch
rfer_c = Metabolite(
    'rfer_c',
    name='reduced 2Fe-2S ferredoxin',
    compartment='c') #cytosol

oxfer_c = Metabolite(
    'oxfer_c',
    name='oxidised 2Fe-2S ferredoxin',
    compartment='c) #cytosol

#retrieve these metabolites from the reconstruction
NADH_c = HydA_model.metabolites.get_by_id("nadh_c")
NAD_c = HydA_model.metabolites.get_by_id("nad_c")
NADPH_c = HydA_model.metabolites.get_by_id("nadph_c")
NADP_c = HydA_model.metabolites.get_by_id("nadp_c")
H_c = HydA_model.metabolites.get_by_id("h_c")
H2_c = HydA_model.metabolites.get_by_id("h2_c")

FDX_rxn.add_metabolites({
    oxfer_c: -2.0,
    NADH_c: -1.0,
    rfer_c: 2.0,
    H_c: 1.0,
    NAD_c: 1.0,
})

FNR_rxn.add_metabolites({
    oxfer_c: -2.0,
    NADPH_c: -1.0,
    rfer_c: 2.0,
    H_c: 1.0,
    NADP_c: 1.0,
})

HydA_rxn.add_metabolites({
    rfer_c: -2.0,
    H_c: -2.0,
    oxfer_c: 2.0,
    H2_c: 1.0
})

HydA_model.add_reactions([FDX_rxn]) 
HydA_model.add_reactions([FNR_rxn])
HydA_model.add_reactions([HydA_rxn])
#%% create FDH strain
FDH_model = WTmodel.copy()

FDH_rxn = Reaction('FDH')
FDH_rxn.name = 'Formate dehydrogenase'
FDH_rxn.subsystem = 'Hydrogen production'
FDH_rxn.lower_bound = 0.0
FDH_rxn.upper_bound = 1000.0

CO2_c = FDH_model.metabolites.get_by_id("co2_c")
formate_c = FDH_model.metabolites.get_by_id("for_c")
NADH_c = FDH_model.metabolites.get_by_id("nadh_c")
NAD_c = FDH_model.metabolites.get_by_id("nad_c")
H_c = FDH_model.metabolites.get_by_id("h_c")

FDH_rxn.add_metabolites({
    CO2_c: -1.0,
    NADH_c: -1.0,
    formate_c: 1.0,
    NAD_c: 1.0
})

FDH_model.add_reactions([FDH_rxn])

#%% plot standard FBA results (biomass growth set as objective)

#Each model copy is for one carbon source for a given strain
WT_models = [WTmodel.copy(), WTmodel.copy(), WTmodel.copy()]
SHI_models = [SHI_model.copy(), SHI_model.copy(), SHI_model.copy()]
HydA_models = [HydA_model.copy(), HydA_model.copy(), HydA_model.copy()]
FDH_models = [FDH_model.copy(), FDH_model.copy(), FDH_model.copy()]

#each list element will correspond to a different carbon source for a given strain
WT_fluxes = []
SHI_fluxes = []
HydA_fluxes = []
FDH_fluxes = []

for i in range(0, 3):
    #set uptake rate for one carbon source for each strain
    WT_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    SHI_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    HydA_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    FDH_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    
    #run FBA optimisation (set to find the maximum possible value of the objective)
    result = WT_models[i].optimize()
    #append flux result of H2 exchange to the list for plotting
    WT_fluxes.append(result.fluxes["EX_h2_e"])
    
    result = SHI_models[i].optimize()
    SHI_fluxes.append(result.fluxes["EX_h2_e"])
    
    result = HydA_models[i].optimize()
    HydA_fluxes.append(result.fluxes["EX_h2_e"])
    
    result = FDH_models[i].optimize()
    FDH_fluxes.append(result.fluxes["EX_h2_e"])

#plot the hydrogen flux results for each strain and carbon source
labels = []
labels.append(plt.scatter(carbon_names, WT_fluxes, c = "blue", s = 100, label = strains[0]))

labels.append(plt.scatter(carbon_names, SHI_fluxes, c = "red", s = 50, label = strains[1]))

labels.append(plt.scatter(carbon_names, HydA_fluxes, c = "green", s = 100, label = strains[2]))

labels.append(plt.scatter(carbon_names, FDH_fluxes, c = "orange", s = 50, label = strains[3]))

plt.ylabel(r'$H_2$ flux [mmol gD$W^{-1}h^{-1}]$')
plt.legend(labels, strains)
plt.show()

#%% FVA (biomass growth set as objective)

WT_models = [WTmodel.copy(), WTmodel.copy(), WTmodel.copy()]
SHI_models = [SHI_model.copy(), SHI_model.copy(), SHI_model.copy()]
HydA_models = [HydA_model.copy(), HydA_model.copy(), HydA_model.copy()]
FDH_models = [FDH_model.copy(), FDH_model.copy(), FDH_model.copy()]

WT_h2_max = [] #stores maximum H2 flux value for WT strain; each list element corresponds to a different carbon source
WT_h2_min = [] #stores minimum possible H2 flux value for WT strain
SHI_h2_max = []
SHI_h2_min = []
HydA_h2_max = []
HydA_h2_min = []
FDH_h2_max = []
FDH_h2_min = []

for i in range(0, 3):
    #set uptake rate for one carbon source for each strain
    WT_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    SHI_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    HydA_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    FDH_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    
    #run FVA for H2 exchange reaction
    #loopless=True does not allow loops in the optimisation to make the results thermodynamically feasible
    #fraction_of_optimum=1 makes sure FVA corresponds to the standard FBA results with the highest objective value
    result = flux_variability_analysis(WT_models[i], reaction_list=WT_models[i].reactions.get_by_id('EX_h2_e'), loopless=True, fraction_of_optimum=1)
    WT_h2_max.append(result["maximum"][0]) #append maximum H2 flux to the list 
    WT_h2_min.append(result["minimum"][0])
    
    result = flux_variability_analysis(SHI_models[i], reaction_list=SHI_models[i].reactions.get_by_id('EX_h2_e'), loopless=True, fraction_of_optimum=1)
    SHI_h2_max.append(result["maximum"][0])
    SHI_h2_min.append(result["minimum"][0])
    
    result = flux_variability_analysis(HydA_models[i], reaction_list=HydA_models[i].reactions.get_by_id('EX_h2_e'), loopless=True, fraction_of_optimum=1)
    HydA_h2_max.append(result["maximum"][0])
    HydA_h2_min.append(result["minimum"][0])
    
    result = flux_variability_analysis(FDH_models[i], reaction_list=FDH_models[i].reactions.get_by_id('EX_h2_e'), loopless=True, fraction_of_optimum=1)
    FDH_h2_max.append(result["maximum"][0])
    FDH_h2_min.append(result["minimum"][0])

##plot FVA maximum and minimum hydrogen flux results for each strain and carbon source    
labels = []

labels.append(plt.scatter(carbon_names, WT_h2_max, c = "blue", s = 100, label = strains[0]))
plt.scatter(carbon_names, WT_h2_min, c = "blue", s = 100, label = strains[0])

labels.append(plt.scatter(carbon_names, SHI_h2_max, c = "red", s = 50, label = strains[1]))
plt.scatter(carbon_names, SHI_h2_min, c = "red", s = 50, label = strains[1])

labels.append(plt.scatter(carbon_names, HydA_h2_max, c = "green", s = 150, label = strains[2]))
plt.scatter(carbon_names, HydA_h2_min, c = "green", s = 150, label = strains[2])

labels.append(plt.scatter(carbon_names, FDH_h2_max, c = "orange", s = 100, label = strains[3]))
plt.scatter(carbon_names, FDH_h2_min, c = "orange", s = 100, label = strains[3])

plt.ylabel(r'$H_2$ flux [mmol gD$W^{-1}h^{-1}]$')
plt.legend(labels, strains)
plt.show()

#%% make copies of strain models for FBA/FVA optimization with hydrogen exchange as an objective
WT_H2 = WTmodel.copy()
WT_H2.objective = "EX_h2_e"

SHI_H2 = SHI_model.copy()
SHI_H2.objective = "EX_h2_e"

HydA_H2 = HydA_model.copy()
HydA_H2.objective = "EX_h2_e"

FDH_H2 = FDH_model.copy()
FDH_H2.objective = "EX_h2_e"
#%% plotting basic FBA results (H2 exchange as objective)

WT_models = [WT_H2.copy(), WT_H2.copy(), WT_H2.copy()]
SHI_models = [SHI_H2.copy(), SHI_H2.copy(), SHI_H2.copy()]
HydA_models = [HydA_H2.copy(), HydA_H2.copy(), HydA_H2.copy()]
FDH_models = [FDH_H2.copy(), FDH_H2.copy(), FDH_H2.copy()]

#each list element will correspond to a different carbon source
WT_fluxes = []
SHI_fluxes = []
HydA_fluxes = []
FDH_fluxes = []

for i in range(0, 3):
    WT_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    SHI_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    HydA_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    FDH_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    
    result = WT_models[i].optimize()
    WT_fluxes.append(result.fluxes["EX_h2_e"])
    
    result = SHI_models[i].optimize()
    SHI_fluxes.append(result.fluxes["EX_h2_e"])
    
    result = HydA_models[i].optimize()
    HydA_fluxes.append(result.fluxes["EX_h2_e"])
    
    result = FDH_models[i].optimize()
    FDH_fluxes.append(result.fluxes["EX_h2_e"])


labels = []
labels.append(plt.scatter(carbon_names, WT_fluxes, c = "blue", s = 100, label = strains[0]))

labels.append(plt.scatter(carbon_names, SHI_fluxes, c = "red", s = 50, label = strains[1]))

labels.append(plt.scatter(carbon_names, HydA_fluxes, c = "green", s = 100, label = strains[2]))

labels.append(plt.scatter(carbon_names, FDH_fluxes, c = "orange", s = 50, label = strains[3]))

plt.ylabel(r'$H_2$ flux [mmol gD$W^{-1}h^{-1}]$')
plt.legend(labels, strains)
plt.show()

#%% FVA (H2 exchange as objective)
WT_models = [WT_H2.copy(), WT_H2.copy(), WT_H2.copy()]
SHI_models = [SHI_H2.copy(), SHI_H2.copy(), SHI_H2.copy()]
HydA_models = [HydA_H2.copy(), HydA_H2.copy(), HydA_H2.copy()]
FDH_models = [FDH_H2.copy(), FDH_H2.copy(), FDH_H2.copy()]

WT_h2_max = []
WT_h2_min = []
SHI_h2_max = []
SHI_h2_min = []
HydA_h2_max = []
HydA_h2_min = []
FDH_h2_max = []
FDH_h2_min = []

for i in range(0, 3):
    WT_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    SHI_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    HydA_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    FDH_models[i].reactions.get_by_id(carbons[i]).lower_bound = -20.
    
    result = flux_variability_analysis(WT_models[i], reaction_list=WT_models[i].reactions.get_by_id('EX_h2_e'), loopless=True, fraction_of_optimum=1)
    WT_h2_max.append(result["maximum"][0])
    WT_h2_min.append(result["minimum"][0])
    
    result = flux_variability_analysis(SHI_models[i], reaction_list=SHI_models[i].reactions.get_by_id('EX_h2_e'), loopless=True, fraction_of_optimum=1)
    SHI_h2_max.append(result["maximum"][0])
    SHI_h2_min.append(result["minimum"][0])
    
    result = flux_variability_analysis(HydA_models[i], reaction_list=HydA_models[i].reactions.get_by_id('EX_h2_e'), loopless=True, fraction_of_optimum=1)
    HydA_h2_max.append(result["maximum"][0])
    HydA_h2_min.append(result["minimum"][0])
    
    result = flux_variability_analysis(FDH_models[i], reaction_list=FDH_models[i].reactions.get_by_id('EX_h2_e'), loopless=True, fraction_of_optimum=1)
    FDH_h2_max.append(result["maximum"][0])
    FDH_h2_min.append(result["minimum"][0])
       
labels = []
labels.append(plt.scatter(carbon_names, WT_h2_max, c = "blue", s = 100, label = strains[0]))
plt.scatter(carbon_names, WT_h2_min, c = "blue", s = 100, label = strains[0])

labels.append(plt.scatter(carbon_names, SHI_h2_max, c = "red", s = 50, label = strains[1]))
plt.scatter(carbon_names, SHI_h2_min, c = "red", s = 50, label = strains[1])

labels.append(plt.scatter(carbon_names, HydA_h2_max, c = "green", s = 100, label = strains[2]))
plt.scatter(carbon_names, HydA_h2_min, c = "green", s = 100, label = strains[2])

labels.append(plt.scatter(carbon_names, FDH_h2_max, c = "orange", s = 50, label = strains[3]))
plt.scatter(carbon_names, FDH_h2_min, c = "orange", s = 50, label = strains[3])

plt.ylabel(r'$H_2$ flux [mmol gD$W^{-1}h^{-1}]$')
plt.legend(labels, strains)
plt.show()
