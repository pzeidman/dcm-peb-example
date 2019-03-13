%% Load PEB prerequisites

% Load design matrix
dm = load('../design_matrix.mat');
X        = dm.X;
X_labels = dm.labels;

% Import downloaded GCM file if needed
if ~exist('../analyses/GCM_full.mat','file')
    copyfile('../analyses/GCM_full_pre_estimated.mat', ...
        '../analyses/GCM_full.mat');
end

% Load GCM
GCM=load('../analyses/GCM_full.mat');
GCM=GCM.GCM;

% PEB settings
M = struct();
M.Q      = 'all';
M.X      = X;
M.Xnames = X_labels;
M.maxit  = 256;
%% Build PEB (using B parameters)
[PEB_B,RCM_B] = spm_dcm_peb(GCM,M,{'B'});
save('../analyses/PEB_B.mat','PEB_B','RCM_B');
%% Automatic search
BMA_B = spm_dcm_peb_bmc(PEB_B);
save('../analyses/BMA_search_B.mat','BMA_B');     
%% Hypothesis-based analysis (B)

% Load estimated PEB
load('../analyses/PEB_B.mat');

% Load template models
templates = load('../analyses/GCM_templates.mat');

% Run model comparison
[BMA,BMR] = spm_dcm_peb_bmc(PEB_B, templates.GCM);

% Show connections in winning model 4
BMA.Kname(BMA.K(4,:)==1)

% Show connections in winning model 15
BMA.Kname(BMA.K(15,:)==1)

save('../analyses/BMA_B_28models.mat','BMA','BMR');

%% Family analysis

% Load the result from the comparison of 28 reduced models
load('../analyses/BMA_B_28models.mat');

% Compare families
[BMA_fam_task,fam_task] = spm_dcm_peb_bmc_fam(BMA, BMR, templates.task_family, 'ALL');

[BMA_fam_b_dv,fam_b_dv] = spm_dcm_peb_bmc_fam(BMA, BMR, templates.b_dv_family, 'NONE');

[BMA_fam_b_lr,fam_b_lr] = spm_dcm_peb_bmc_fam(BMA, BMR, templates.b_lr_family, 'NONE');

save('../analyses/BMA_fam_task.mat','BMA_fam_task','fam_task');
save('../analyses/BMA_fam_b_dv.mat','BMA_fam_b_dv','fam_b_dv');
save('../analyses/BMA_fam_b_lr.mat','BMA_fam_b_lr','fam_b_lr');
%% LOO
[qE,qC,Q] = spm_dcm_loo(GCM,M,{'B(4,4,3)'});
save('../analyses/LOO_rdF_words.mat','qE','qC','Q');
%% Correlate rdF
B = cellfun(@(x)x.Ep.B(4,4,3),GCM(:,1));
LI = X(:,2);
figure;scatter(LI,B);
lsline;
[R,P] = corrcoef(LI,B);
%% Build PEB (A)
[PEB_A,RCM_A] = spm_dcm_peb(GCM(:,1),M,{'A'});
save('../analyses/PEB_A.mat','PEB_A','RCM_A');  
%% Search-based analysis (A)
load('../analyses/PEB_A.mat');
BMA_A = spm_dcm_peb_bmc(PEB_A);
save('../analyses/BMA_search_A.mat','BMA_A');
spm_dcm_peb_review(BMA_A,GCM);