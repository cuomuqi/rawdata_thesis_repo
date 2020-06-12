function deLRH = dget_eLRH(t, eLRH, Food, Temp, kx, v, L_T, l_v, l_s, l_j, L_m, g, kap, E_m, k_J, T_ref, Trange, MC2, E_Hp)
%set temp
    if size(Temp,2)==1                              %zero_variate_model, temperature is constant
        Tcorr = tempcorr(Temp, T_ref, Trange);      %Calculate temperature correction factor
    else                                            % uni_variate_model, temperature is variable
        T = spline1(t,Temp);                        %find temperature at certain time t
        Tcorr = tempcorr(C2K(T), T_ref, Trange);         % Calculate temperature correction factor
    end
    
%set food
    if size(Food,2)==1       %  zero_variate_model, food is constant
        f= Food;             % if food is constant, f is needed
    else                    % uni_variate_model, food is variable
      x = spline1(t, Food); % -, food density at t
      %X = 1e-3 * d_X * x/ w_X * celltomg;               % C-mol/l, food concentration
      f = x /(kx + x);                                   %-, functional response (between 0-1) 
    end
    
% find initial conditions 
  ee = eLRH(1); L = eLRH(2); ER = eLRH(3); EH = eLRH(4); %set state variables at time 0
  
% Calculate parameters to be used in diff equations
%run get_ts
  L_v = l_v * L_m; L_s = L_m * l_s; L_j = L_m * l_j;    %Calculate structural length at s, j and m
  L_m_j = L_m * L_j/ L_s;                               %Calculates max structural length (L_j/L_s = acceleration factor)
  s_M = min(L_j/ L_s, max(1, L/ L_s));                  % -, find acceleration factor (1, before s; L/L_s between s and j; L_j/L_s after j) 

% Change parameters to temp and acceleration  
  vT = v * Tcorr * s_M;                 % v is a rate, and should therefore be adapted to temperature and acceleration factor
  kTJ = Tcorr * k_J;                    % Energy spend on maturity maintenance
  
% calculate r and mobilisation
  r = vT * (ee/L - (1 + L_T/ L)/ L_m_j)/ (ee + g); % 1/d, spec growth rate
  p_C = (vT/ L - r) * ee * E_m * L^3;              % mobilisation

        Mcof = min (1, (MC2*Tcorr.* (L-L_v).^2+0.2)) ; % so when L-L_s =0, this coefficient is 1, and MC2 is a factor to make it to increase faster becuase in all L matrix there is not the case it so outside th L_v and L_s
  % Differential equations                        
  de = (f*Mcof - ee) * vT/ L; % 1/d, change in scaled reserve density e
  dL = L * r/ 3; % cm/d, change in structural length L 

  if EH < E_Hp
      dER = 0;
      dEH = (1 - kap) * p_C - kTJ * EH; % change in maturity
  else
      dER = (1 - kap) * p_C - kTJ * E_Hp; % J/d, change in reproduction buffer
      dEH = 0;
  end
  
% Pack to output
  %deLRH = [de; dL; dER; dEH; s_M; Mcof; vT; kTJ; r; p_C ]; % pack output
  deLRH = [de; dL; dER; dEH]; % pack output
  
end
