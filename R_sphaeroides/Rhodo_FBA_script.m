%% Flux balance analysis on Rhodobacter sphearoides
%Aleksandra Sobieska
%University of Edinburgh Undergraduate iGEM team 2019

%if needed, initialise Cobra Toolbox and choose a solver
initCobraToolbox
changeCobraSolver('gurobi', 'all');
%%
%load Rhodobacter sphearoides reconstruction 
%(assuming it is in the same folder as this script)
modelRho = readCbModel('iRsp1140.xml'); 
%set H2 exchange as an objective
modelRho = changeObjective(modelRho,'RXN1205'); 
%% find all exchange reactions
index = find(contains(modelRho.rxnNames, 'exchange'));
%convert index numbers of the exchange reactions to their names
exchangeRxns = modelRho.rxns(index);
%close flux in all exchange reactions 
modelRho = changeRxnBounds(modelRho, exchangeRxns, 0, 'l');
%% set anaerobic conditions
modelRho = changeRxnBounds(modelRho, 'RXN1331', -1000, 'l'); %open photon exchange
modelRho = changeRxnBounds(modelRho, 'RXN0223', 0, 'l'); %remove oxygen uptake
%% incorporate metabolic modifications by changing reaction bounds of these reactions
% which portray anaerobic conditions
% from the supplemental data file S1 found in: 
% Burger BT et al. Combining Genome-Scale Experimental and Computational Methods 
% To Identify Essential Genes in Rhodobacter sphaeroides. mSystems. 2017;2(3). 
% doi:10.1128/mSystems.00015-17.

%"_c" added to metabolite's name means that the metabolite is found in cytosol

%set flux of O2 transport via diffusion to zero
modelRho = changeRxnBounds(modelRho, 'RXN0949', 0, 'b'); 
%close o2_c + fmnh2_c --> e4p_c + h2o_c + dmbzid_c
modelRho = changeRxnBounds(modelRho, 'RXN1318', 0, 'b');
%close flux of 3 o2_c + pppg9_c --> 3 h2o2_c + ppp9_c
modelRho = changeRxnBounds(modelRho, 'RXN0661', 0, 'b');
%close the reaction of coproporphyrinogen III oxidase
modelRho = changeRxnBounds(modelRho, 'RXN0630', 0, 'b');
%close the flux of the reaction of magnesium-protoporphyrin IX monomethyl ester cyclase
modelRho = changeRxnBounds(modelRho, 'RXN0663', 0, 'b');
%close the reaction fluxes of pyridoxamine 5-phosphate oxidases
modelRho = changeRxnBounds(modelRho, 'RXN1063', 0, 'b');
modelRho = changeRxnBounds(modelRho, 'RXN1064', 0, 'b');
%close the reaction flux of cytochrome c oxidase
modelRho = changeRxnBounds(modelRho, 'RXN0515', 0, 'b');
%close the reaction flux of 2-octaprenylphenol hydroxylase
modelRho = changeRxnBounds(modelRho, 'RXN1020', 0, 'b');
%close the reaction flux of fdxrd_c + o2_c + hedacp_c + 2.0 h_c <=> cll_428_c + fdxox_c + 2.0 h2o_c
modelRho = changeRxnBounds(modelRho, 'RXN1345', 0, 'b');
%make cobyrinic acid a,c-diamide synthase irreversible
modelRho = changeRxnBounds(modelRho, 'RXN0648', 0, 'l');
%% set minimal medium composition - open lower bounds to allow uptake of substances
%taken from the already mentioned supplemental data file.

modelRho = changeRxnBounds(modelRho, 'RXN0213', -1000, 'l'); %phosphate
modelRho = changeRxnBounds(modelRho, 'RXN0224', -3, 'l'); %succinate
modelRho = changeRxnBounds(modelRho, 'RXN0214', -1, 'l'); %glutamate
modelRho = changeRxnBounds(modelRho, 'RXN0196', -1000, 'l'); %sulfate
modelRho = changeRxnBounds(modelRho, 'RXN1326', -1000, 'l'); %biotin
modelRho = changeRxnBounds(modelRho, 'RXN1329', -1000, 'l'); %nicotinate
modelRho = changeRxnBounds(modelRho, 'RXN1158', -1000, 'l'); %magnesium
modelRho = changeRxnBounds(modelRho, 'RXN1354', -1000, 'l'); %thiamine
modelRho = changeRxnBounds(modelRho, 'RXN0192', -1000, 'l'); %iron
modelRho = changeRxnBounds(modelRho, 'RXN0222', -1, 'l'); %CO2
modelRho = changeRxnBounds(modelRho, 'RXN0217', -1000, 'l'); %water
modelRho = changeRxnBounds(modelRho, 'RXN0187', -1000, 'l'); %cobalamin
modelRho = changeRxnBounds(modelRho, 'RXN0220', -1000, 'l'); %sodium
modelRho = changeRxnBounds(modelRho, 'RXN1167', -1000, 'l'); %potassium
modelRho = changeRxnBounds(modelRho, 'RXN0221', -1000, 'l'); %calcium
modelRho = changeRxnBounds(modelRho, 'RXN0191', -1000, 'l'); %manganese
modelRho = changeRxnBounds(modelRho, 'RXN0190', -1000, 'l'); %zinc
modelRho = changeRxnBounds(modelRho, 'RXN0188', -1000, 'l'); %calcium
modelRho = changeRxnBounds(modelRho, 'RXN0197', -1000, 'l'); %molybdate
modelRho = changeRxnBounds(modelRho, 'RXN0193', -1000, 'l'); %cobalt
modelRho = changeRxnBounds(modelRho, 'RXN1091', -1000, 'l'); %nickel


%% run standard FBA optimisation on wild type R. sphaeroides strain
% 'max' setting in optimizeCbModel allows finding the highest objective flux (in this case H2
% production rate) possible
% 'false' setting - does not allow loops in the optimisation to make it
% thermodynamically feasible
FBAsol = optimizeCbModel(modelRho, 'max', false);
%print flux of H2 exchange reaction
disp('WT strain - H2 exchange flux')
disp(FBAsol.x(1076))


%% create HydA strain of R. spharoides
% by adding 3 reactions to WT strain

% all metabolites assumed to be found in cytosol

%2 reduced ferredoxin + 2 H+ <-> 2 oxidised ferredoxin + H2
modelHydA = addReaction(modelRho, 'HydA hydrogenase',...
'metaboliteList', {'CPD0651', 'CPD0063', 'CPD0652', 'CPD0627'},...
'stoichCoeffList', [-2; -2; 2; 1]);

%2 oxidised ferredoxin + NADPH <-> 2 reduced ferredoxin + H+ + NADP+
modelHydA = addReaction(modelHydA, 'NADPH ferredoxin reductase',...
'metaboliteList', {'CPD0652', 'CPD0005', 'CPD0651', 'CPD0063', 'CPD0006'},...
'stoichCoeffList', [-2; -1; 2; 1; 1]);

%2 oxidised ferredoxin + NADH <-> 2 reduced ferredoxin + H+ + NAD+
modelHydA = addReaction(modelHydA, 'NADH ferredoxin reductase',...
'metaboliteList', {'CPD0652', 'CPD0004', 'CPD0651', 'CPD0063', 'CPD0003'},...
'stoichCoeffList', [-2; -1; 2; 1; 1]);
%% run FBA optimisation on HydA strain
FBAsolHydA = optimizeCbModel(modelHydA, 'max', false);
disp('HydA strain - H2 exchange flux')
disp(FBAsolHydA.x(1076))


%% create SHI strain with the following reaction: 
modelSHI = addReaction(modelRho, 'SHI',...
'metaboliteList', {'CPD0005', 'CPD0063', 'CPD0006', 'CPD0627'},...
'stoichCoeffList', [-1; -1; 1; 1]);
% all metabolites assumed to be found in cytosol

%% run FBA optimisation on SHI strain
FBAsolSHI = optimizeCbModel(modelSHI, 'max', false);
disp('SHI strain - H2 exchange flux')
disp(FBAsolSHI.x(1076))
