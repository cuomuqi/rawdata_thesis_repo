function [par, metaPar, txtPar] = pars_init_Ostrea_edulis(metaData)

metaPar.model = 'asj'; 

% reference parameter (not to be changed)
par.T_ref = C2K(20); free.T_ref = 0;          units.T_ref = 'K';              label.T_ref = 'Reference temperature';

%% core primary parameters
par.z = 0.7033;         free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 0.839;         free.F_m   = 1;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;        free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;        free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.01132;         free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.9148;       free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;       free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 13.33;         free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;            free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.001964;     free.k_J   = 1;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 2349;         free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 0.0001877;   free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hr  = 0.0003829;  free.E_Hr  = 1;   units.E_Hr  = 'J';        label.E_Hr  = 'maturity at release';
par.E_Hv = 0.002047;   free.E_Hv = 1;    units.E_Hv = 'J';         label.E_Hv = 'maturity at start metamorphosis velum';
par.E_Hs = 0.002766;    free.E_Hs  = 1;   units.E_Hs = 'J';         label.E_Hs = 'maturity at settlement'; 
par.E_Hj = 2.13;       free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hp = 177.1;       free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 1.213e-09;    free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 4.966e-05;    free.s_G   = 1;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% Temperature parameters
par.T_A = 4080;         free.T_A   = 1;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 

%% auxiliary parameters
par.del_M = 0.2847;     free.del_M = 1;   units.del_M = '-';        label.del_M = 'shape coefficient';
par.del_Mb= 0.6773;     free.del_Mb =1;   units.del_Mb = '-';       label.del_Mb = 'shape coefficient at birth';

%% environmental parameters (temperatures are in auxData)
par.f = 1.0;            free.f     = 0;   units.f = '-';            label.f    = 'scaled functional response for 0-var data';
par.f_Saur = 0.4693;    free.f_Saur = 1;  units.f_Saur = '-';       label.f_Saur = 'scaled functional response for tL_adult1 data'; 
par.f_Labarta = 1;      free.f_Labarta= 0;units.f_Labarta = '-';    label.f_Labarta = 'food level for Labarta';
par.L_0Saur = 1.445;    free.L_0Saur = 1; units.L_0Saur = 'cm';     label.L_0Saur = 'initial physical length in tL_Saur data'; 

par.kx = 5.978e6;       free.kx = 1;      units.kx = 'cm';          label.kx = 'half saturation food constant';
par.e_r = 1;       free.e_r = 1;     units.e_r = '-';          label.e_r = 'energy density at release';
par.MC2 = 1.961e+05;        free.MC2= 1;      units.MC2 = '-';          label.MC2 = 'the accleration factor to make the Mcof decreasing faster';
par.E_X = 8.639e-08;    free.E_X = 1;     units.E_X = 'J/cell';     label.E_X = 'energy content of the feed (cell mix)';
%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class);

%% Pack output:
txtPar.units = units; txtPar.label = label; par.free = free; 

