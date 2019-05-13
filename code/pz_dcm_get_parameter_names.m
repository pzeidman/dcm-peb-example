function names = pz_dcm_get_parameter_names(DCM)
% Gets the names of the parameters in a DCM
if isfield(DCM,'Ep')
    P = DCM.Ep;
elseif isfield(DCM,'M') && isfield(DCM.M,'pE')
    P = DCM.M.pE;
else
    error('Sorry this only works for estimated DCMs');
end

np = length(spm_vec(P));

names = spm_fieldindices(P,1:np);