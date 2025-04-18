function [] = makeTable(DataStruc, VarList)

% All runs within DataStruc share the same following parameters:
Wo = DataStruc.(VarList{1,1}).Params_phys.Wo_paper;
Ar = DataStruc.(VarList{1,1}).Params_phys.Ar_paper;
Sigma = DataStruc.(VarList{1,1}).Params_phys.sigma_dimal;
StrainMeasure = DataStruc.(VarList{1,1}).Params_phys.strainmeasure;
Strain = DataStruc.(VarList{1,1}).Params_phys.frac;

Title = strcat('Results/Wo=', num2str(Wo), '_Ar=', num2str(Ar), ...
    '_sigma0=', num2str(Sigma), '_StrainM=', StrainMeasure, '_Strain=', num2str(Strain));

ColumnNames = {'K - forw. (mN/m)', 'G - forw. (mN/m)', 'K - Tensiom. (mN/m)', 'G - Tensiom. (mN/m)', 'K - Circ. (mN/m)', ...
    'Tau_apex - forw. (mN/m)', 'Tau_apex - Tensiom. (mN/m)', 'Tau_apex - Circ. (mN/m)', ...
    'Dil.Str_apex - forw ( )', 'Dil.Str_apex - Tensiom. ( )', 'Dil.Str_apex - Circ. ( )'};

idx = 1;
for i = 1:size(VarList, 1)
    for j = 1:size(VarList, 2)
        K_forw_vals(idx,1) = DataStruc.(VarList{i,j}).Params_phys.Kmod_dimal;
        G_forw_vals(idx,1) = DataStruc.(VarList{i,j}).Params_phys.Gmod_dimal;
        K_Tensiom_vals(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv.Kmod_dimal;
        G_Tensiom_vals(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv.Gmod_dimal;
        K_Circle_vals(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv.KmodCircleFit_dimal;
        Tau_forw_vals(idx,1) = DataStruc.(VarList{i,j}).DefSt_fwd.sigmas_dimal(1,1);
        Tau_Tensiom_vals(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv.sigmas_dimal(1,1);
        Tau_Circle_vals(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv.SigmaCircleFit_dimal(1,1);
        DilStr_forw_vals(idx,1) = DataStruc.(VarList{i,j}).DefSt_fwd.lams(1,1)*DataStruc.(VarList{i,j}).DefSt_fwd.lamp(1,1);
        DilStr_Tensiom_vals(idx,1) = DataStruc.(VarList{i,j}).DefSt_inv.lams(1,1)*DataStruc.(VarList{i,j}).DefSt_inv.lamp(1,1);
        DilStr_Circle_vals(idx,1) = DataStruc.(VarList{i,j}).Params_phys.frac;
        idx = idx + 1;
    end
end

T = table(K_forw_vals, G_forw_vals, K_Tensiom_vals, G_Tensiom_vals, K_Circle_vals, ...
    Tau_forw_vals, Tau_Tensiom_vals, Tau_Circle_vals, DilStr_forw_vals, DilStr_Tensiom_vals, DilStr_Circle_vals, ...
    'VariableNames', ColumnNames);
disp(sprintf('<strong> %100s </strong>',Title));
disp(T);

end