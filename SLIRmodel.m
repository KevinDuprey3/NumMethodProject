%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SLIR function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLIRmodel combines all the SLIR equations into one system.  
% the dependent variables are stored in y as:
% y(1) = S (susceptible population)
% y(2) = L (latent population)
% y(3) = I (infectious population)
% y(4) = R (recovered/removed population)
%
% with parameters
% p(1) = beta (rate of new infections)
% p(2) = mu_l (inverse length of latent period in days)
% p(3) = mu_i (inverse length of the infectious period in days)
% p(4) = k (population growth rate)
% p(5) = e (rate of new infections from external sources)

function [dydt] = SLIRmodel(t,y,p,T)
    dt = 1/24;
        % calculate mu_L
        t = floor(t*24)+1; % round value for now to not throw error
%         Nsteps = length(T);
%         mu_L_target = 6;
%    
%        mu_l_inv = latentperiod(t,dt,Nsteps,mu_L_target,...
%           zeros(size(T)),T);
        % temperature

mu_L_target = p(2);
istart=1; 
PT   = zeros(size(T));
mu_L = zeros(size(T));
for i=1:length(T)
    PT(i) = Sall_temp_effect(T(i));
    mu_L(i) = sum(PT(istart:i));
    while(mu_L(i)>mu_L_target)
        istart=istart+1; %increment until we are at the min value timestep
        mu_L(i) = sum(PT(istart:i));
    end
end
mu_l_inv = 1./mu_L; %we want the inverse of the latent period length


        if T(t) > 0 && T(t) < 35
            Tb = 0.000214 * T(t)^(2.06737) * (35 - T(t))^(0.72859);
        else
            Tb = 0;
        end

        Te = -0.35968 + 0.10789*T(t) - 0.00214*(T(t)^2); %STOLL CHANGE
    % day
    tday = t*dt + 30; % find current day based on index %STOLL CHANGE

    %assign parameters
    beta = Tb * p(1);
    mu_i = p(3);
    k    = p(4);
    e    = p(5);
    Ap   = p(6);
    dpldt = (1.33 * tday) * Te;
    mu_i_inv = 1/mu_i;
    

    %assign variables
    S = y(1);
    L = y(2);
    I = y(3);
    R = y(4);
    P = y(5);
    Pb = y(6);

    dydt(6) = (0.1724*Pb - 0.0000212*(Pb^2)) * Te; % Pb
    dydt(5) = dydt(6) + dpldt;          % P
    dydt(1) = -beta*S*I + dydt(5) * (1/Ap); % S
    dydt(2) = beta*S*I - mu_l_inv(t)*L + e;         % L %STOLL CHANGE
    dydt(3) = mu_l_inv(t)*L-mu_i_inv*I;            % I
    dydt(4) = mu_i_inv*I;                   % R
    
end