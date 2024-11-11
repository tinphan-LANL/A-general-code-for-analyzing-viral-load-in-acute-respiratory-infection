clear all
close all
rng default

%% Set up
cd ..
load('Admissible_params.mat')
cd 'Treatment_on_day_6_8d'

Admissible_params = [Admissible_params_non_rebound; Admissible_params_rebound];
close all;
%%
    rebound_threshold = 4; %log10
    
    t_i   = 6; %day of treatment
    t_e   = t_i+8;

% total population served by DITP
T0 = 80000000;
R0 = 0;
E0 = 1;
I0 = 0;
V0 = 0;

ICs  = [T0 R0 E0 I0 V0];

% tf = 28;
% tspan = linspace(0,tf,tf*100+1);

%% Take some samples
sample_size = 100; %number of patients per test
test_size = 20; %number of tests
param_size = length(Admissible_params(1,:)); %number of params
REBOUND_HOLDER = nan(1,test_size);
INDEX_HOLDER = zeros(test_size,sample_size);

PARAM_HOLDER = nan(param_size,sample_size);
PARAM_SUPER_HOLDER = cell(1,test_size); % this holds PARAM_HOLDER

for kk = 1:test_size

figure(2+kk); hold on; box on;
NUM_REBOUND = 0;

ymax = 11;
area([t_i t_e],[ymax ymax])
colororder([0.5 0.5 0.5])

random_indexes = randi(length(Admissible_params),sample_size);

for jj = 1:sample_size

params    = Admissible_params(random_indexes(jj),:);
params(8) = t_i;

PARAM_HOLDER(:,jj) = params';

[t,y] = ode45(@model,tspan,ICs,[],params);
plot(t,log10(y(:,5)),'-k','LineWidth',2)


[V_max,T_max] = max(y(t_e*100:end,5));

if log10(V_max) > rebound_threshold
    
    NUM_REBOUND = NUM_REBOUND + 1;
    INDEX_HOLDER(kk,jj) = 1; 

end

end

yline(log10(rebound_threshold),'--g')

xline(t_i)
xline(t_e)

ylim([0 inf])
xlim([0 tf])
ylabel('Log10 viral load (/mL)')
xlabel('Days post infection')

REBOUND_HOLDER(kk) = NUM_REBOUND;
PARAM_SUPER_HOLDER{kk} = PARAM_HOLDER;

end
% [beta_exp, p_exp, rho, phi_exp, delta_0, t00, I_rho_exp, n, t_i]
Rebound_index = find(INDEX_HOLDER'>0);
%%%BETA
REBOUND_BETA = nan(1,length(Rebound_index));
NON_REBOUND_BETA = nan(1,test_size*sample_size - length(Rebound_index));
%%%PI
REBOUND_PI = nan(1,length(Rebound_index));
NON_REBOUND_PI = nan(1,test_size*sample_size - length(Rebound_index));
%%%RHO
REBOUND_RHO = nan(1,length(Rebound_index));
NON_REBOUND_RHO = nan(1,test_size*sample_size - length(Rebound_index));
%%%PHI
REBOUND_PHI = nan(1,length(Rebound_index));
NON_REBOUND_PHI = nan(1,test_size*sample_size - length(Rebound_index));
%%%DELTA
REBOUND_DELTA = nan(1,length(Rebound_index));
NON_REBOUND_DELTA = nan(1,test_size*sample_size - length(Rebound_index));
%%%TT0
REBOUND_TT0 = nan(1,length(Rebound_index));
NON_REBOUND_TT0 = nan(1,test_size*sample_size - length(Rebound_index));
%%%I_RHO
REBOUND_I_RHO = nan(1,length(Rebound_index));
NON_REBOUND_I_RHO = nan(1,test_size*sample_size - length(Rebound_index));
%%%N
REBOUND_N = nan(1,length(Rebound_index));
NON_REBOUND_N = nan(1,test_size*sample_size - length(Rebound_index));


for mm = 1:test_size

    if mm == 1
        Rebound_index_holder = Rebound_index(mm:REBOUND_HOLDER(mm));
    else
        Rebound_index_holder = Rebound_index(sum(REBOUND_HOLDER(1:mm-1)) + 1 : sum(REBOUND_HOLDER(1:mm)));
        Rebound_index_holder = Rebound_index_holder - (mm-1)*sample_size;
    end

PARAM_HOLDER = PARAM_SUPER_HOLDER{mm};

%% Beta
BETA_HOLDER = PARAM_HOLDER(1,:);

    if mm == 1
        REBOUND_BETA(mm:REBOUND_HOLDER(mm)) = BETA_HOLDER(Rebound_index_holder);
            BETA_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_BETA(mm:sample_size - REBOUND_HOLDER(mm)) = BETA_HOLDER;
    else
        REBOUND_BETA(sum(REBOUND_HOLDER(1:mm-1)) + 1 : sum(REBOUND_HOLDER(1:mm))) = BETA_HOLDER(Rebound_index_holder);
            BETA_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_BETA((mm-1)*sample_size - sum(REBOUND_HOLDER(1:mm-1)) + 1: mm*sample_size - sum(REBOUND_HOLDER(1:mm))) = BETA_HOLDER;
    end

%% Pi
PI_HOLDER = PARAM_HOLDER(2,:);

    if mm == 1
        REBOUND_PI(mm:REBOUND_HOLDER(mm)) = PI_HOLDER(Rebound_index_holder);
            PI_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_PI(mm:sample_size - REBOUND_HOLDER(mm)) = PI_HOLDER;
    else
        REBOUND_PI(sum(REBOUND_HOLDER(1:mm-1)) + 1 : sum(REBOUND_HOLDER(1:mm))) = PI_HOLDER(Rebound_index_holder);
            PI_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_PI((mm-1)*sample_size - sum(REBOUND_HOLDER(1:mm-1)) + 1: mm*sample_size - sum(REBOUND_HOLDER(1:mm))) = PI_HOLDER;
    end

%% Rho
RHO_HOLDER = PARAM_HOLDER(3,:);

    if mm == 1
        REBOUND_RHO(mm:REBOUND_HOLDER(mm)) = RHO_HOLDER(Rebound_index_holder);
            RHO_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_RHO(mm:sample_size - REBOUND_HOLDER(mm)) = RHO_HOLDER;
    else
        REBOUND_RHO(sum(REBOUND_HOLDER(1:mm-1)) + 1 : sum(REBOUND_HOLDER(1:mm))) = RHO_HOLDER(Rebound_index_holder);
            RHO_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_RHO((mm-1)*sample_size - sum(REBOUND_HOLDER(1:mm-1)) + 1: mm*sample_size - sum(REBOUND_HOLDER(1:mm))) = RHO_HOLDER;
    end

%% Phi
PHI_HOLDER = PARAM_HOLDER(4,:);

    if mm == 1
        REBOUND_PHI(mm:REBOUND_HOLDER(mm)) = PHI_HOLDER(Rebound_index_holder);
            PHI_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_PHI(mm:sample_size - REBOUND_HOLDER(mm)) = PHI_HOLDER;
    else
        REBOUND_PHI(sum(REBOUND_HOLDER(1:mm-1)) + 1 : sum(REBOUND_HOLDER(1:mm))) = PHI_HOLDER(Rebound_index_holder);
            PHI_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_PHI((mm-1)*sample_size - sum(REBOUND_HOLDER(1:mm-1)) + 1: mm*sample_size - sum(REBOUND_HOLDER(1:mm))) = PHI_HOLDER;
    end

%% Delta
DELTA_HOLDER = PARAM_HOLDER(5,:);

    if mm == 1
        REBOUND_DELTA(mm:REBOUND_HOLDER(mm)) = DELTA_HOLDER(Rebound_index_holder);
            DELTA_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_DELTA(mm:sample_size - REBOUND_HOLDER(mm)) = DELTA_HOLDER;
    else
        REBOUND_DELTA(sum(REBOUND_HOLDER(1:mm-1)) + 1 : sum(REBOUND_HOLDER(1:mm))) = DELTA_HOLDER(Rebound_index_holder);
            DELTA_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_DELTA((mm-1)*sample_size - sum(REBOUND_HOLDER(1:mm-1)) + 1: mm*sample_size - sum(REBOUND_HOLDER(1:mm))) = DELTA_HOLDER;
    end

%% tt0
TT0_HOLDER = PARAM_HOLDER(6,:);

    if mm == 1
        REBOUND_TT0(mm:REBOUND_HOLDER(mm)) = TT0_HOLDER(Rebound_index_holder);
           TT0_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_TT0(mm:sample_size - REBOUND_HOLDER(mm)) = TT0_HOLDER;
    else
        REBOUND_TT0(sum(REBOUND_HOLDER(1:mm-1)) + 1 : sum(REBOUND_HOLDER(1:mm))) = TT0_HOLDER(Rebound_index_holder);
            TT0_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_TT0((mm-1)*sample_size - sum(REBOUND_HOLDER(1:mm-1)) + 1: mm*sample_size - sum(REBOUND_HOLDER(1:mm))) = TT0_HOLDER;
    end

%% I_RHO
I_RHO_HOLDER = PARAM_HOLDER(7,:);

    if mm == 1
        REBOUND_I_RHO(mm:REBOUND_HOLDER(mm)) = I_RHO_HOLDER(Rebound_index_holder);
            I_RHO_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_I_RHO(mm:sample_size - REBOUND_HOLDER(mm)) = I_RHO_HOLDER;
    else
        REBOUND_I_RHO(sum(REBOUND_HOLDER(1:mm-1)) + 1 : sum(REBOUND_HOLDER(1:mm))) = I_RHO_HOLDER(Rebound_index_holder);
            I_RHO_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_I_RHO((mm-1)*sample_size - sum(REBOUND_HOLDER(1:mm-1)) + 1: mm*sample_size - sum(REBOUND_HOLDER(1:mm))) = I_RHO_HOLDER;
    end

%% N
N_HOLDER = PARAM_HOLDER(8,:);

    if mm == 1
        REBOUND_N(mm:REBOUND_HOLDER(mm)) = N_HOLDER(Rebound_index_holder);
           N_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_N(mm:sample_size - REBOUND_HOLDER(mm)) = N_HOLDER;
    else
        REBOUND_N(sum(REBOUND_HOLDER(1:mm-1)) + 1 : sum(REBOUND_HOLDER(1:mm))) = N_HOLDER(Rebound_index_holder);
            N_HOLDER(Rebound_index_holder) = []; 
        NON_REBOUND_N((mm-1)*sample_size - sum(REBOUND_HOLDER(1:mm-1)) + 1: mm*sample_size - sum(REBOUND_HOLDER(1:mm))) = N_HOLDER;
    end

end

%% parameter distribution 1
ONE   = ones(size(REBOUND_BETA));
TWO   = 2*ones(size(NON_REBOUND_BETA));

if isempty(REBOUND_BETA) ~= 1
g1 = figure(100);
g1.WindowState = 'maximized';
tt1 = tiledlayout(2,4);

%%
nexttile

hh1 = boxplot(([REBOUND_BETA, NON_REBOUND_BETA]),[ONE, TWO], ...
        'Labels',{'Rebound', 'Non-rebound'}, ...
        'BoxStyle','outline','Colors', [1 0 0; 0 0 1]);

set(gca,'FontSize',16)
title('- log_{10} (\beta)', 'FontSize', 18)
% axis tight
xlim([0.5 2.5])
% ylim([7.7 8])
%%
nexttile

boxplot(([REBOUND_PI, NON_REBOUND_PI]),[ONE, TWO], ...
        'Labels',{'Rebound', 'Non-rebound'}, ...
        'BoxStyle','outline','Colors', [1 0 0; 0 0 1]);

set(gca,'FontSize',16)
title('log_{10} (\pi)', 'FontSize', 18)
% axis tight
xlim([0.5 2.5])
% ylim([2.4 3])

%%
nexttile

boxplot(([REBOUND_RHO, NON_REBOUND_RHO]),[ONE, TWO], ...
        'Labels',{'Rebound', 'Non-rebound'}, ...
        'BoxStyle','outline','Colors', [1 0 0; 0 0 1]);

set(gca,'FontSize',16)
title('\rho', 'FontSize', 18)
% axis tight
xlim([0.5 2.5])
% ylim([0.1 0.7])

%%
nexttile

boxplot(([REBOUND_PHI, NON_REBOUND_PHI]),[ONE, TWO], ...
        'Labels',{'Rebound', 'Non-rebound'}, ...
        'BoxStyle','outline','Colors', [1 0 0; 0 0 1]);

set(gca,'FontSize',16)
title('-log_{10} (\phi)', 'FontSize', 18)
% axis tight
xlim([0.5 2.5])
% ylim([6.6 6.8])

% [beta_exp, p_exp, rho, phi_exp, delta_0, t00, I_rho_exp, n, t_i]

%%
nexttile

boxplot(([REBOUND_DELTA, NON_REBOUND_DELTA]),[ONE, TWO], ...
        'Labels',{'Rebound', 'Non-rebound'}, ...
        'BoxStyle','outline','Colors', [1 0 0; 0 0 1]);

set(gca,'FontSize',16)
title('\delta_0', 'FontSize', 18)
% axis tight
xlim([0.5 2.5])
% ylim([1 4.2])

%%
nexttile

boxplot(([REBOUND_TT0, NON_REBOUND_TT0]),[ONE, TWO], ...
        'Labels',{'Rebound', 'Non-rebound'}, ...
        'BoxStyle','outline','Colors', [1 0 0; 0 0 1]);

set(gca,'FontSize',16)
title('t^*', 'FontSize', 18)
% axis tight
xlim([0.5 2.5])
% ylim([10 22])
%%
nexttile

boxplot(([REBOUND_I_RHO, NON_REBOUND_I_RHO]),[ONE, TWO], ...
        'Labels',{'Rebound', 'Non-rebound'}, ...
        'BoxStyle','outline','Colors', [1 0 0; 0 0 1]);

set(gca,'FontSize',16)
title('I_{\rho}', 'FontSize', 18)
% axis tight
xlim([0.5 2.5])
% ylim([1 2.5])

%%
nexttile

boxplot(([REBOUND_N, NON_REBOUND_N]),[ONE, TWO], ...
        'Labels',{'Rebound', 'Non-rebound'}, ...
        'BoxStyle','outline','Colors', [1 0 0; 0 0 1]);

set(gca,'FontSize',16)
title('n', 'FontSize', 18)
% axis tight
xlim([0.5 2.5])
% ylim([0.3 1.5])

%%

p_beta  = ranksum(REBOUND_BETA,NON_REBOUND_BETA);
p_pi    = ranksum(REBOUND_PI,NON_REBOUND_PI);
p_rho   = ranksum(REBOUND_RHO,NON_REBOUND_RHO);
p_phi   = ranksum(REBOUND_PHI,NON_REBOUND_PHI);
p_delta = ranksum(REBOUND_DELTA,NON_REBOUND_DELTA);
p_tt0   = ranksum(REBOUND_TT0,NON_REBOUND_TT0);
p_I_rho = ranksum(REBOUND_I_RHO,NON_REBOUND_I_RHO);
p_n     = ranksum(REBOUND_N,NON_REBOUND_N);

p_values = [p_beta p_pi p_rho p_phi p_delta p_tt0 p_I_rho p_n];

% ax = gca;
% ax.FontSize = 16;

% title('Treatment on day 2', 'FontSize',16)

f1 = gcf;
exportgraphics(f1,'Params_distribution_treatment_day_6.png','Resolution',300)
end

save('REBOUND_COUNT_6','REBOUND_HOLDER')

%% parameter distribution 2
% g2 = figure(200); box on; hold on;
% g2.WindowState = 'maximized';
% 
% X = [REBOUND_BETA', REBOUND_PI', REBOUND_PHI', REBOUND_RHO', REBOUND_TT0', REBOUND_I_RHO'];
% [~,ax] = plotmatrix(X,'.k');
% 
% ax(1,1).YLabel.String = '\beta';
% ax(2,1).YLabel.String = '\pi';
% ax(3,1).YLabel.String = '\phi';
% ax(4,1).YLabel.String = '\rho';
% ax(5,1).YLabel.String = 't^0';
% ax(6,1).YLabel.String = '\sigma';
% 
% 
% ax(6,1).XLabel.String = '\beta';
% ax(6,2).XLabel.String = '\pi';
% ax(6,3).YLabel.String = '\phi';
% ax(6,4).YLabel.String = '\rho';
% ax(6,5).XLabel.String = 't^0';
% ax(6,6).XLabel.String = '\sigma';
% 
% ax(1,1).FontSize = 16;
% ax(2,1).FontSize = 16;
% ax(3,1).FontSize = 16;
% ax(4,1).FontSize = 16;
% ax(5,1).FontSize = 16;
% ax(6,1).FontSize = 16;
% 
% ax(6,1).FontSize = 16;
% ax(6,2).FontSize = 16;
% ax(6,3).FontSize = 16;
% ax(6,4).FontSize = 16;
% ax(6,5).FontSize = 16;
% ax(6,6).FontSize = 16;
% 
% title('Rebound - Treatment on day 2','FontSize',16)
% 
% %    
% f2 = gcf;
% exportgraphics(f2,'Rebound_treatment_day_2.png','Resolution',300)

%% parameter distribution 2
% g3 = figure(300); box on; hold on;
% g3.WindowState = 'maximized';
% 
% X = [NON_REBOUND_BETA', NON_REBOUND_PI', NON_REBOUND_PHI', NON_REBOUND_RHO', NON_REBOUND_TT0', NON_REBOUND_I_RHO'];
% [~,ax] = plotmatrix(X,'.k');
% 
% ax(1,1).YLabel.String = '\beta';
% ax(2,1).YLabel.String = '\pi';
% ax(3,1).YLabel.String = '\phi';
% ax(4,1).YLabel.String = '\rho';
% ax(5,1).YLabel.String = 't^0';
% ax(6,1).YLabel.String = '\sigma';
% 
% 
% ax(6,1).XLabel.String = '\beta';
% ax(6,2).XLabel.String = '\pi';
% ax(6,3).YLabel.String = '\phi';
% ax(6,4).YLabel.String = '\rho';
% ax(6,5).XLabel.String = 't^0';
% ax(6,6).XLabel.String = '\sigma';
% 
% ax(1,1).FontSize = 16;
% ax(2,1).FontSize = 16;
% ax(3,1).FontSize = 16;
% ax(4,1).FontSize = 16;
% ax(5,1).FontSize = 16;
% ax(6,1).FontSize = 16;
% 
% ax(6,1).FontSize = 16;
% ax(6,2).FontSize = 16;
% ax(6,3).FontSize = 16;
% ax(6,4).FontSize = 16;
% ax(6,5).FontSize = 16;
% ax(6,6).FontSize = 16;
% 
% title('Non-rebound - Treatment on day 2','FontSize',16)
% 
% %    
% f3 = gcf;
% exportgraphics(f3,'Non_rebound_treatment_day_2.png','Resolution',300)


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function dY = model(t,y,param)
    beta_exp   = param(1); beta  = 10^(-beta_exp);
    p_exp      = param(2); p     = 10^p_exp;
    rho        = param(3);
    phi_exp    = param(4); phi   = 10^(-phi_exp);
    delta0     = param(5);

    t00        = param(6);
    I_rho_exp  = param(7); I_rho = 10^I_rho_exp;
    n          = 1; %param(8);

    ti         = param(8);
    te         = ti+7.5;

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