clear all
close all

rng default
%% Run Monolix fit first - then load the estimated parameter distribution

%beta p rho phi delta t00 TT I_rho (these are from the population estimates)
%can also be loaded directly from monolix bestfit.
Pop_est = [8.06, 2.65, 0.17, 6.49, 1.85, 13.1, 2.06, 0.67];
SD_RE   = [0.52, 0.035, 0.99, 0.077, 0.28, 0.61, 0.39, 0.95];
logit_index = [1, 3, 4, 6];

NUM_PAR = length(Pop_est);

%beta_exp is logit between 7.5 and 9
%rho is logit between 0 and 1
%phi_exp is logit between 5 and 12
%t00 is logit between 7 and 28

%% Set up random number to pick parameters
NUM   = 50000; %number of random samples
RANDOM_PARS = nan(NUM,NUM_PAR);
% note: a and b are the bounded range of the parameters (logit).

% beta
    b = 9; a = 7.5;
    logit_beta = log((Pop_est(1) - a)/(b - (Pop_est(1))));
    m_beta = normrnd(logit_beta,SD_RE(1)^2, [1,NUM]);

    RANDOM_PARS(:,1) = (b.*exp(m_beta) + a)./(1+exp(m_beta));

% p
    eta_p = normrnd(0,SD_RE(2)^2, [1,NUM]);

    RANDOM_PARS(:,2) = Pop_est(2).*exp(eta_p);

% rho
    b = 1; a = 0;
    logit_rho = log((Pop_est(3) - a)/(b - (Pop_est(3))));
    m_rho = normrnd(logit_rho,SD_RE(3)^2, [1,NUM]);

    RANDOM_PARS(:,3) = (b.*exp(m_rho) + a)./(1+exp(m_rho));

% phi
    b = 12; a = 5;
    logit_phi = log((Pop_est(4) - a)/(b - (Pop_est(4))));
    m_phi = normrnd(logit_phi,SD_RE(4)^2, [1,NUM]);

    RANDOM_PARS(:,4) = (b.*exp(m_phi) + a)./(1+exp(m_phi));

% delta
    eta_delta = normrnd(0,SD_RE(5)^2, [1,NUM]);

    RANDOM_PARS(:,5) = Pop_est(5).*exp(eta_delta);

% t00
    b = 28; a = 7;
    logit_t00 = log((Pop_est(6) - a)/(b - (Pop_est(6))));
    m_t00 = normrnd(logit_t00,SD_RE(6)^2, [1,NUM]);

    RANDOM_PARS(:,6) = (b.*exp(m_t00) + a)./(1+exp(m_t00));

% TT
    eta_TT = normrnd(0,SD_RE(7)^2, [1,NUM]);

    RANDOM_PARS(:,7) = Pop_est(7).*exp(eta_TT);

% I_rho
    eta_Irho = normrnd(0,SD_RE(8)^2, [1,NUM]);

    RANDOM_PARS(:,8) = Pop_est(8).*exp(eta_Irho);

%% Simulation 
T0 = 80000000;
R0 = 0;
E0 = 1;
I0 = 0;
V0 = 0;

ICs  = [T0 R0 E0 I0 V0];

tf = 28;
separation = 100; %every 100 index = 1 day
tspan = linspace(0,tf,tf*separation+1);

%% LHS
Admissible_params_rebound = nan(NUM,NUM_PAR);
Admissible_params_non_rebound = nan(NUM,NUM_PAR);

V_cutoff     = 6; 
T_cutoff_min = 2;
T_cutoff_max = 7;

for ii=1:NUM
params    = RANDOM_PARS(ii,:);
TT        = params(7);
params(7) = []; %remove the TT from the input parameters
params(6) = 13; %fixed timing of adaptive immune response at 12 day post infection
params(8) = 100; %no treatment is indicated with 100

        % Initialization
            T0 = 8*10^7;
            R0 = 0;
            E0 = 1;
            I0 = 0;
            V0 = 0;

            X_int = [T0, R0, E0, I0, V0];

        % Time vector
            t_f = 28; %day post infection
            n = t_f*100+1;
            t_span = linspace(0,t_f,n);

        % Solve ODE
            opts  = odeset('RelTol',1e-10,'AbsTol',1e-10);
            [t y] = ode45(@model,t_span,X_int,opts,params);

        % Assign variables
            T = y(:,1);
            R = y(:,2);
            E = y(:,3);
            I = y(:,4);
            V = log10(y(:,5));

%Note: every 100 index is one day
[V_max,T_max] = max(y(:,5));
    
    y(:,4) = real(y(:,4));
    TF = islocalmin(y(:,4));
    first_local_min_index = find(TF>0,1);

[V_max_2,T_max_2] = max(y(first_local_min_index:end,4));
T_max_2 = T_max + T_max_2;

if log10(V_max) > 6 %bound range for peak (VL)

    if T_max > 200 && T_max < 701 %bound range for peak (time)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(TF) == 0 %if local min is not empty

            if log10(y(first_local_min_index,4)) < 4 %must dip below 3 log10

                if log10(V_max_2) > 4 %must rise above 3 log10 (second peak)

                    if log10(y(end,4)) < 2 %rebound but must be cleared by day 28

                        Admissible_params_rebound(ii,:) = params;

                    end

                else %rebound but less than 3 log10

                    if log10(y(end,4)) < 2 %must be below 2 log10 by day 28

                        Admissible_params_non_rebound(ii,:) = params;

                    end

                end
            end 
        else %if local min is empty (no rebound at all)

            if log10(y(end,4)) < 2 %must be below 2 log10 by day 28 if continuously decreasing
            
                    Admissible_params_non_rebound(ii,:) = params;

            end

        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

end

Admissible_params_rebound = rmmissing(Admissible_params_rebound);
Admissible_params_non_rebound = rmmissing(Admissible_params_non_rebound);


%% Rebound
figure(1);

t1 = tiledlayout(2,4); 

nexttile; hold on; box on;
histogram((Admissible_params_rebound(:,1)),25);
xlabel('-log10 \beta')

nexttile; hold on; box on;
histogram((Admissible_params_rebound(:,2)),25);
xlabel('log10 p')

nexttile; hold on; box on;
histogram((Admissible_params_rebound(:,3)),25);
xlabel('\rho')

nexttile; hold on; box on;
histogram((Admissible_params_rebound(:,4)),25);
xlabel('- log10 phi')

nexttile; hold on; box on;
histogram(Admissible_params_rebound(:,5),25);
xlabel('\delta_0')

nexttile; hold on; box on;
histogram(Admissible_params_rebound(:,6),25);
xlabel('t^*')

nexttile; hold on; box on;
histogram(Admissible_params_rebound(:,7),25);
xlabel('log10 I_{\rho}')

nexttile; hold on; box on;
histogram(Admissible_params_rebound(:,8),25);
xlabel('n')

title(t1, 'Rebound')

f1 = gcf;
exportgraphics(f1,'distribution_admissible_pars_rebound.png','Resolution',300)

%% NON-Rebound
figure(2); 

t2 = tiledlayout(2,4); 
title(t2, 'Non-rebound')

nexttile; hold on; box on;
histogram((Admissible_params_non_rebound(:,1)),25);
xlabel('-log10 \beta')

nexttile; hold on; box on;
histogram((Admissible_params_non_rebound(:,2)),25);
xlabel('log10 p')

nexttile; hold on; box on;
histogram((Admissible_params_non_rebound(:,3)),25);
xlabel('\rho')

nexttile; hold on; box on;
histogram((Admissible_params_non_rebound(:,4)),25);
xlabel('- log10 phi')

nexttile; hold on; box on;
histogram(Admissible_params_non_rebound(:,5),25);
xlabel('\delta_0')

nexttile; hold on; box on;
histogram(Admissible_params_non_rebound(:,6),25);
xlabel('t^*')

nexttile; hold on; box on;
histogram(Admissible_params_non_rebound(:,7),25);
xlabel('log10 I_{\rho}')

nexttile; hold on; box on;
histogram(Admissible_params_non_rebound(:,8),25);
xlabel('n')

f2 = gcf;
exportgraphics(f2,'distribution_admissible_pars_non_rebound.png','Resolution',300)

%% plot a sample (Rebound)
figure(3); hold on; box on;
title('Rebound')
sample_size = 5;
random_indexes = randi(length(Admissible_params_rebound(:,1)),sample_size);

for jj = 1:sample_size
params = Admissible_params_rebound(random_indexes(jj),:);
[t,y] = ode45(@model,tspan,ICs,[],params);
plot(t,log10(y(:,4)),'-k','LineWidth',2)

end

xline(T_cutoff_min,'--g')
xline(T_cutoff_max,'--g')
yline(V_cutoff,'--r')
% yline(4,'-.r')

ylim([2 inf])
xlim([0 tf])
ylabel('Log10 viral load (/mL)')
xlabel('Days post infection')


f3 = gcf;
exportgraphics(f3,'samples_accepted_trajectories_rebound.png','Resolution',300)

%% plot a sample (Non-rebound)
figure(4); hold on; box on;
title('Non-rebound')
sample_size = 5; 
random_indexes = randi(length(Admissible_params_non_rebound(:,1)),sample_size);

for jj = 1:sample_size
params = Admissible_params_non_rebound(random_indexes(jj),:);
[t,y] = ode45(@model,tspan,ICs,[],params);
plot(t,log10(y(:,4)),'-k','LineWidth',2)

end

xline(T_cutoff_min,'--g')
xline(T_cutoff_max,'--g')
yline(V_cutoff,'--r')

ylim([2 inf])
xlim([0 tf])
ylabel('Log10 viral load (/mL)')
xlabel('Days post infection')


f4 = gcf;
exportgraphics(f4,'samples_accepted_trajectories_non_rebound.png','Resolution',300)

%% save
save('Admissible_params')

%%
function dY = model(t,y,param)
    beta_exp   = param(1); beta  = 10^(-beta_exp);
    p_exp      = param(2); p     = 10^p_exp;
    rho        = param(3);
    phi_exp    = param(4); phi   = 10^(-phi_exp);
    delta0     = param(5);

    t00        = param(6);
    I_rho_exp  = param(7); I_rho = 10^I_rho_exp;
    n          = 1;

    ti         = param(8);
    te         = ti+4.5;

    k          = 4;
    c          = 10;
    sigma      = 0.5;
    EC50       = 62;
    eps        = 1;

    %adaptive immunity

    delta_max = 20*delta0;

    if t<t00
	    delta = delta0;
    else
        delta = delta_max - (delta_max - delta0)*exp(-sigma*(t-t00));
    end

    param_A = [ti,te];
    if ti < 100
        C = mAb(t,param_A);
        eff = eps*C/(C+EC50);
    else
        eff = 0;
    end
    
    T = y(1);  
    R = y(2);
    E = y(3);      
    I = y(4);
    V = y(5);

    dT = - beta*T*V - phi*I*T + rho*(1-I^n/(I^n+I_rho^n))*R;

    dR = phi*I*T - rho*(1-I^n/(I^n+I_rho^n))*R;

    dE = beta*T*V - k*E;
                        
    dI = k*E - delta*I;

    dV = (1-eff)*p*I - c*V;

    dY = [dT; dR; dE; dI; dV];
end

%% sub-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = mAb(t,p)
Id = 1/2; %drug taken once every 12 hours
ka   = 17.53;   %callculated from tmax and ke
ke   = 2.77;    %based on elimination half-life of 6.05 hours
Ch   = 6.25e3;
ti = p(1);
te = p(2);
Ce = 2.4774e3;

kae = ka/(ke-ka);
% te  = te-1/2;

    if t<ti
        
        C = 0;

    elseif t>ti && t<=te

        eke = exp(-ke*t);
        Nd  = floor(t/Id) + 1;

            C = Ch*kae*(eke/(exp(ka*Id)-1))*(1 - exp((ke-ka)*t)*(1-exp(Nd*ka*Id)) ...
            + (exp(ke*Id) - exp(ka*Id))*((exp((Nd-1)*ke*Id) -1)/(exp(ke*Id)-1)) ...
            - exp(((Nd-1)*ke+ka)*Id));
        
    elseif t>te

        t = t-te;
        eke = exp(-ke*t);
        eka = exp(-ka*t);
        C = Ce*exp(-ke*(t+1/2)) + Ch*kae*(eka-eke);
    
    end
end