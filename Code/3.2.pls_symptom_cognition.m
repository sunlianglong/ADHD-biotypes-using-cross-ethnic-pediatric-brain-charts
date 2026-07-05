%% PLS analysis for brain_symptoms
% setting dir
clear all
root_dir = '/Users/buxuan/BNU_Helab/data/nm_subtyping_adhd/analysis_v2/CT_subtypes_v2/pls';
data_dir = fullfile(root_dir,'input/');
out_dir = fullfile(root_dir,'output_symptoms/');

addpath(fullfile(root_dir,'code','PLS_MIPlab'));
addpath(fullfile(root_dir,'code','replication','code'));

% Options
normalization_img = 1; % normalization options for brain data (X)
normalization_behav = 1; % normalization options for behavior data (Y)
% 0 = no normalization
% 1 = zscore across all subjects

%% Load data and PLS analysis
brain = csvread([data_dir,'z1_symp.csv'],1,5); 
clinic = readtable([data_dir,'symp1.csv']);
symptom = clinic(:,["inattention","hyperactive_impulsive","ADHD_index"]);
symptom = table2array(symptom);

X = brain;
Y = symptom;
[U,S,V,Lx,Ly,explCovLC,LC_behav_loadings,LC_brain_loadings] = ...
    myPLS_analysis(X,Y,normalization_img,normalization_behav);

%% Permutation testing
nPerms = 5000;
nSubj = size(X,1);
diagnosis_grouping = ones(nSubj,1);
pvals_LC = myPLS_permut(X,Y,U,S,nPerms,diagnosis_grouping,normalization_img,normalization_behav,[]);
myStats=[explCovLC'; pvals_LC'];

%% Save results
csvwrite([out_dir,'stats_subtype1.csv'],myStats);

csvwrite([out_dir,'brain_scores1.csv'],Lx(:,1));
csvwrite([out_dir,'brain_weights1.csv'],V(:,1));
csvwrite([out_dir,'brain_loadings1.csv'],LC_brain_loadings(:,1));

csvwrite([out_dir,'symptom_scores1.csv'],Ly(:,1));
csvwrite([out_dir,'symptom_weights1.csv'],U(:,1));
csvwrite([out_dir,'symptom_loadings1.csv'],LC_behav_loadings(:,1));

%% Correlation between significant latent components (Lx and Ly)
% true correlation
[r_corr,p_corr_self] = corr(Lx(:,1),Ly(:,1));

% permutation 5000 times for correlation (permuted the behavioral data)
perm_index = zeros(nPerms,size(Y,1));
for i=1:nPerms
    perm_index(i,:) = randperm(size(Y,1));
end

corr_perm = zeros(1,nPerms);
for j = 1:nPerms
    disp(j);
    Yp=Ly(perm_index(j,:),1);
    corr_perm(1,j) = corr(Lx(:,1),Yp);
end
clear Yp
p_corr = length(find(corr_perm(1,:)>r_corr))/(nPerms+1);
myCorr = [r_corr,p_corr];

csvwrite([out_dir,'PLS_corr1.csv'],myCorr);

%% Bootstrapping (5000 times)
signif_LC = [1];
nBootstraps = 5000;
grouping = ones(nSubj,1);

% run bootstrapping
[LC_brain_loadings_boot,LC_behav_loadings_boot,all_boot_orders] = CBIG_VK2019_bootstrap_loadings...
    (X,Y,U,signif_LC,nBootstraps,grouping,normalization_img,normalization_behav,[]);

% stats bootstrapping (sd, z-score, p-val)
% behavior loadings
nBehav = size(LC_behav_loadings,1);
for iter_lc = 1:size(signif_LC,1)   
    for iter_behav = 1:nBehav

        this_lc = signif_LC(iter_lc);
        
        % Std across samples
        std_behav_boot(iter_behav,iter_lc) = std(LC_behav_loadings_boot(iter_behav,iter_lc,:));
               
        % Z-score (original loading / std loading across samples)
        zscore_behav_boot(iter_behav,iter_lc) = ...
            LC_behav_loadings(iter_behav,this_lc) / std_behav_boot(iter_behav,iter_lc);
        
        % P-values for z-scores
        if zscore_behav_boot(iter_behav,iter_lc) >= 0
            pvals_behav_boot(iter_behav,iter_lc) = 1-cdf('norm',zscore_behav_boot(iter_behav,iter_lc),0,1);
        elseif zscore_behav_boot(iter_behav,iter_lc) < 0
            pvals_behav_boot(iter_behav,iter_lc) = cdf('norm',zscore_behav_boot(iter_behav,iter_lc),0,1);
        end
        
    end    
end

% brain loadings
nRois = size(LC_brain_loadings,1);
for iter_lc = 1:size(signif_LC,1)   
    for iter_brain = 1:nRois

        this_lc = signif_LC(iter_lc);
        
        % Std across samples
        std_brain_boot(iter_brain,iter_lc) = std(LC_brain_loadings_boot(iter_brain,iter_lc,:));
               
        % Z-score (original loading / std loading across samples)
        zscore_brain_boot(iter_brain,iter_lc) = ...
            LC_brain_loadings(iter_brain,this_lc) / std_brain_boot(iter_brain,iter_lc);
        
        % P-values for z-scores
        if zscore_brain_boot(iter_brain,iter_lc) >= 0
            pvals_brain_boot(iter_brain,iter_lc) = 1-cdf('norm',zscore_brain_boot(iter_brain,iter_lc),0,1);
        elseif zscore_brain_boot(iter_brain,iter_lc) < 0
            pvals_brain_boot(iter_brain,iter_lc) = cdf('norm',zscore_brain_boot(iter_brain,iter_lc),0,1);
        end
        
    end    
end

% save
csvwrite([out_dir,'sym_boot_z1.csv'],zscore_behav_boot);
csvwrite([out_dir,'sym_boot_p1.csv'],pvals_behav_boot);
csvwrite([out_dir,'brain_boot_z1.csv'],zscore_brain_boot);
csvwrite([out_dir,'brain_boot_p1.csv'],pvals_brain_boot);

% save all results
clear i iter_brain iter_behav iter_lc j this_lc;
save(fullfile(out_dir,'PLSresults_symp_subtype1.mat'));
clear
clc
