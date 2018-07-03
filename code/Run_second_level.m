%% Load PEB prerequisites

% Load design matrix
dm = load('../design_matrix.mat');
X        = dm.X;
X_labels = dm.labels;

% Load GCM
GCM=load('../analyses/GCM_full.mat');
GCM=GCM.GCM;

M = struct();
M.alpha  = 1;
M.beta   = 16;
M.hE     = 0;
M.hC     = 1/16;
M.Q      = 'all';
M.X      = X;
M.Xnames = X_labels;
%% Build PEB (B,C)
[PEB_BC,RCM_BC] = spm_dcm_peb(GCM(:,1),M,{'B','C'});
save('../analyses/PEB_BC.mat','PEB_BC','RCM_BC');
%% Automatic search (B,C)
BMA_BC = spm_dcm_peb_bmc(PEB_BC);
save('../analyses/BMA_search_B.mat','BMA_BC');          
%% Hypothesis-based analysis (B,C)
load('../analyses/PEB_BC.mat');
templates = load('../analyses/GCM_templates.mat');

% Run model comparison
[BMA,BMR] = spm_dcm_peb_bmc(PEB_BC, templates.GCM);

% Show connections in winning model 31
BMA.Kname(BMA.K(31,:)==1)

% Show connections in winning model 69
BMA.Kname(BMA.K(69,:)==1)

save('../analyses/BMA_BC_112models.mat','BMA','BMR');

%% Family analysis
[BMA_fam_task,fam_task] = spm_dcm_peb_bmc_fam(BMA,BMR,task_family,'NONE');

[BMA_fam_c,fam_c]       = spm_dcm_peb_bmc_fam(BMA,BMR,c_family,'NONE');

[BMA_fam_b_dv,fam_b_dv] = spm_dcm_peb_bmc_fam(BMA,BMR,b_dv_family,'NONE');

[BMA_fam_b_lr,fam_b_lr] = spm_dcm_peb_bmc_fam(BMA,BMR,b_lr_family,'NONE');

save('../analyses/BMA_fam_task.mat','BMA_fam_task','fam_task');
save('../analyses/BMA_fam_c.mat','BMA_fam_c','fam_c');
save('../analyses/BMA_fam_b_dv.mat','BMA_fam_b_dv','fam_b_dv');
save('../analyses/BMA_fam_b_lr.mat','BMA_fam_b_lr','fam_b_lr');
%% LOO
[qE,qC,Q] = spm_dcm_loo(GCM(:,1),M,{'B(4,4,3)'});
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