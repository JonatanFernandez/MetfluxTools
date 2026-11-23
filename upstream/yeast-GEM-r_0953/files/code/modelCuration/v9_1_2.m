% This scripts applies curations to be applied on yeast-GEM release 9.1.1,
% to get to yeast-GEM release 9.1.2.
%
% Curations in this release:
%   - Removal of the thermodynamically infeasible pyridoxal oxidase r_0953
%     (Issue #NNN)

%% Load yeast-GEM 9.1.1 (requires local yeast-GEM git repository)
cd ..
codeDir=pwd();
model = getEarlierModelVersion('9.1.1');
model.id='yeastGEM_develop';
model.version='';
% dataDir=fullfile(pwd(),'..','data','modelCuration','v9.1.2'); % No dataDir required for these curations
cd modelCuration

%% ========================================================================
% Remove pyridoxal oxidase r_0953 (Issue #NNN)
%
% r_0953 (BiGG PYDXO, MetaNetX MNXR103365, ModelSEED rxn09566; no KEGG or
% Rhea reaction) reads
%   pyridoxal + ammonium + 2 H2O + 0.5 O2 => pyridoxamine + 2 H2O2
% It is mass and charge balanced, but it reduces pyridoxal to pyridoxamine
% with electrons taken from water, which is oxidised to hydrogen peroxide.
% Its Gibbs energy is +375 +/- 7 kJ/mol with all reactants at 1 mM
% (eQuilibrator 0.8.1, pH 7.2, I = 0.25 M, pMg 3), of which +262 kJ/mol is
% the formation of 3 H2O2 from water and O2 alone (NBS tables), so no
% intracellular concentrations can make it run in the direction allowed by
% its bounds. It is one of the 22 alarms of the reversibility check in
% Issue #428.
% Correcting it to the oxidative deamination of pyridoxamine (the
% non-phosphorylated analogue of r_0954) is not supported either: EC
% 1.4.3.5, Rhea and UniProt (P38075, PDX3_YEAST) list only the
% 5'-phosphates as substrates, and the E. coli orthologue PdxH did not
% accept non-phosphorylated substrates (doi:10.1016/0167-4838(94)00235-9).
% The reactions of the PLP salvage route from pyridoxamine, pyridoxamine
% kinase (r_2026) and pyridoxamine 5'-phosphate oxidase (r_0954), are kept.
% Growth and gene essentiality are unchanged. r_0954 becomes blocked:
% pyridoxamine and its 5'-phosphate had no other source in the model (the
% corrected form would block it as well). See Issue #NNN for the analysis.
model = removeReactions(model,'r_0953',true,true,true);

%% ========================================================================

%% DO NOT CHANGE OR REMOVE THE CODE BELOW THIS LINE.
% Show some metrics:
cd(fullfile(codeDir,'modelTests'))
disp('Run gene essentiality analysis')
[new.accuracy,new.tp,new.tn,new.fn,new.fp] = essentialGenes(model);
fprintf('Genes in model: %d\n',numel(model.genes));
fprintf('Gene essentiality accuracy: %.4f\n', new.accuracy);
fprintf('True non-essential genes: %d\n', numel(new.tp));
fprintf('True essential genes: %d\n', numel(new.tn));
fprintf('False non-essential genes: %d\n', numel(new.fp));
fprintf('False essential genes: %d\n', numel(new.fn));
fprintf('\nRun growth analysis\n')
R2=growth(model);
fprintf('R2 of growth prediction: %.4f\n', R2);

% Save model:
cd(codeDir)
saveYeastYaml(model);
cd modelCuration
