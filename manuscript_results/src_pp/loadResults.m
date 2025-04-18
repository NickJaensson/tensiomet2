function [DataStruc] = loadResults(path)

RefSt_fwd = load(strcat(path, '/vars_sol_ref.mat'));
RefSt_fwd = RefSt_fwd.vars_sol_ref;
RefSt_num = load(strcat(path, '/vars_num_ref.mat'));
RefSt_fwd.s = RefSt_num.vars_num_ref.s;

DefSt_fwd = load(strcat(path, '/vars_sol.mat'));
DefSt_fwd = DefSt_fwd.vars_sol;
DefSt_num = load(strcat(path, '/vars_num.mat'));
DefSt_fwd.s = DefSt_num.vars_num.s;

RefSt_inv = load(strcat(path, '/vars_sol_ref_fit.mat'));
RefSt_inv = RefSt_inv.vars_sol_ref_fit;
RefSt_num = load(strcat(path, '/vars_num_ref_fit.mat'));
RefSt_inv.s = RefSt_num.vars_num_ref_fit.s; 

DefSt_inv = load(strcat(path, '/vars_sol_fit.mat'));
DefSt_inv = DefSt_inv.vars_sol_fit;
DefSt_num = load(strcat(path, '/vars_num_fit.mat'));
DefSt_inv.s = DefSt_num.vars_num_fit.s;

Params_phys = load(strcat(path, '/params_phys.mat'));
Params_phys = Params_phys.params_phys;
Params_num = load(strcat(path, '/params_num.mat'));
Params_num = Params_num.params_num;

DataStruc.RefSt_fwd = RefSt_fwd;
DataStruc.DefSt_fwd = DefSt_fwd;
DataStruc.RefSt_inv = RefSt_inv;
DataStruc.DefSt_inv = DefSt_inv;
DataStruc.Params_phys = Params_phys;
DataStruc.Params_num = Params_num;

end