for i=1:10
    
close all; 
global pets

pets = {'Ostrea_edulis'};
check_my_pet(pets); 

estim_options('default'); 
estim_options('max_step_number',5000); 
estim_options('max_fun_evals',5000);  

estim_options('pars_init_method', 1);
estim_options('results_output', 2);
estim_options('method', 'nm');

estim_pars; 

end
