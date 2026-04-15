rho_num = 1024;
theta_num = 1024;

variable_name =  ["more_R", "more_Z", "more_B", "more_lambda","more_J", "more_J_B", ...
    "more_gcon_rr", "more_gcon_tt", "more_gcon_zz", ...
    "more_gcon_rt", "more_gcon_rz", "more_gcon_tz", ...
    "more_gcov_rr", "more_gcov_tt", "more_gcov_zz", ...
    "more_gcov_rt", "more_gcov_rz", "more_gcov_tz", ...
    "more_kappa_n", "more_kappa_g", "more_nu", "more_dnu_dt", "more_jp"];

variable_num = length(variable_name);

data_2d = zeros(rho_num,theta_num);

for variable_id = 1:variable_num

    variable_id;

    for rho_id = 1:rho_num
        data_2d(rho_id,:) = eval(variable_name(variable_id)+"_"+string(rho_id)+";");
    end

    eval([char(variable_name(variable_id)), '= data_2d;']);
    
    if(variable_id==1)
        save(strcat(inputPath,'more_equilibrium.mat'), variable_name(variable_id));
    else
        save(strcat(inputPath,'more_equilibrium.mat'), variable_name(variable_id), '-append');
    end

end