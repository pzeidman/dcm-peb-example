% The script runs the step-by-step GUI example from start to end using the
% batch complete_gui_batch.mat.
% -------------------------------------------------------------------------

% Where to store analysis results
dir_out = '../analyses/';

% Where to find manually created template DCMs
dcm_full = '../GLM/sub-01/DCM_full.mat';
dcm_alt  = {'../GLM/sub-01/DCM_no_ldF_modulation.mat'};

% Find SPM.mat files for all subjects (timing)
spms = cellstr(spm_select('FPListRec','../GLM/','SPM.mat'));
assert(length(spms) == 60);

% Find timeseries for all subjects
lvF = cellstr(spm_select('FPListRec','../GLM/','VOI_lvF'));
ldF = cellstr(spm_select('FPListRec','../GLM/','VOI_ldF'));
rvF = cellstr(spm_select('FPListRec','../GLM/','VOI_rvF'));
rdF = cellstr(spm_select('FPListRec','../GLM/','VOI_rdF'));
assert(length(lvF) == 60);
assert(length(ldF) == 60);
assert(length(rvF) == 60);
assert(length(rdF) == 60);

% Prepare the batch
load('complete_gui_batch.mat');
matlabbatch{1}.spm.dcm.spec.fmri.group.output.dir       = cellstr(dir_out);
matlabbatch{1}.spm.dcm.spec.fmri.group.template.fulldcm = cellstr(dcm_full);
matlabbatch{1}.spm.dcm.spec.fmri.group.template.altdcm  = dcm_alt;
matlabbatch{1}.spm.dcm.spec.fmri.group.data.spmmats     = spms;
matlabbatch{1}.spm.dcm.spec.fmri.group.data.region = {
                                                      lvF
                                                      ldF
                                                      rvF
                                                      rdF
                                                      }';

% Run
spm_jobman('run',matlabbatch);