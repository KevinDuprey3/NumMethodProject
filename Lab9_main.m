%%%%%%%%%%%%%
% Kevin Duprey
% u1406068
% MEEN 2450 - 04
% 11/30/2023
%%%%%%%%%%%%%%

clear all, clc, close all

load('EnvironmentalForcing.mat'); % load given arrays


% 
% % initial conditions
% Bmax = 1; % (1/day)
% uL_min = 10; % min number days latent
% uI = 10; % number days infected
% e = 0.001; % rate of introduction from external sources
% Ap = 5000; % normalization factor
% 
% 
% % initial primary variables into vector
% params = zeros(6,1);
% params(5) = 1.33 * 30 * (-0.35968 + 0.10789*15 - 0.00214*15*15); % P0
% params(1) = params(6)/Ap;     % S0
% params(2) = 0.01 * params(1); % L0
% params(3) = 0;                % I0
% params(4) = uI * params(3);   % R0
% 
% % initialize integrated vectors with initial conditions
% S = zeros(length(tspan),1);
% S(1) = params(1);
% L =  zeros(length(tspan),1);
% L(1) = params(2);
% I =  zeros(length(tspan),1);
% I(1) = params(3);
% R =  zeros(length(tspan),1);
% R(1) = params(4);
% P =  zeros(length(tspan),1);
% P(1) = params(6);
% 
% vec = SLIRP_ddt(Bmax, )
% %

%% actual

% initial conditions
Tb = zeros(1,length(tspan));
for i = 1:length(tspan)
    if T(i) > 0 || T(i) < 35
        Tb(i) = 0.000214 * T(i)^(2.06737) * (35 - T(i))^(0.72859);
    else
        Tb(i) = 0;
    end
end
Te(:) = -0.35968 + 0.10789 .* T(:) + 0.00214*T(:).^2;
beta = 2;    %rate infection increases (1/day)
mu_L = sum(Tb);  %rate latent period ends (inverse of number of days latent)
mu_I = 0.25; %rate infection clears (inverse of number of days infectious)
days = 60;   %length of simulation (day)
dt   = 0.05; %timestep (fraction of a day)
Ap   = 5000; %normalization factor
P_i = 1.33*30*(-0.35968 + 0.10789 *15 + 0.00214*15*15)*30; % initial population
S_i  = P_i/Ap;    %initial size of the population (normalized)
L_i  = S_i*0.01; %initial fraction of population that is latent
I_i  = 0;    %initial fraction of population that is infectious
R_i  = I_i*mu_I;  %initial fraction of population that is recovered
Pb_i = 1; % initial berry population fraction
k    = 0.01;  %population growth rate (fraction per day)
e    = 0.005; %rate of introduction from external sources

[S,L,I,R,P,Pb,time] = PathogenGrowth_0D(S_i,L_i,I_i,R_i,P_i,Pb_i,beta,mu_L,mu_I,k,e,Ap,T,days,dt)








