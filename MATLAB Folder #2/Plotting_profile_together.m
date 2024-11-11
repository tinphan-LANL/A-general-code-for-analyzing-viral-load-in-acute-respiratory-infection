clear all
close all
%%
rng('default')
load('Admissible_params.mat')
close all

Admissible_params = [Admissible_params_non_rebound; Admissible_params_rebound];
%% Rebound & non-rebound
figure(1); 
binsizeZ = 25;
sizeT = 10;

t1 = tiledlayout(2,3); 
% [beta_exp, p_exp, rho, phi_exp, delta_0, t00, I_rho_exp, n]
nexttile; hold on; box on;
h1 = histogram((Admissible_params(:,1)),binsizeZ);%,'Normalization','probability');
xlabel('-log_{10} \beta','FontSize',sizeT)
ylim([0 10000])

nexttile; hold on; box on;
h2 = histogram((Admissible_params(:,2)),binsizeZ);%,'Normalization','probability');
xlabel('log_{10} p','FontSize',sizeT)
ylim([0 10000])

nexttile; hold on; box on;
h3 = histogram((Admissible_params(:,3)),binsizeZ);%,'Normalization','probability');
xlabel('\rho','FontSize',sizeT)
ylim([0 10000])

nexttile; hold on; box on;
h4 = histogram((Admissible_params(:,4)),binsizeZ);%,'Normalization','probability');
xlabel('- log_{10} \phi','FontSize',sizeT)
ylim([0 10000])

nexttile; hold on; box on;
h5 = histogram(Admissible_params(:,5),binsizeZ);%,'Normalization','probability');
xlabel('\delta_0','FontSize',sizeT)
ylim([0 10000])

nexttile; hold on; box on;
h6 = histogram(Admissible_params(:,7),binsizeZ);%,'Normalization','probability');
xlabel('log_{10} K_{\rho}','FontSize',sizeT)
ylim([0 10000])
color1 = [0.5350 0.5350 0.5350];
h1.FaceColor = color1;
h2.FaceColor = color1;
h3.FaceColor = color1;
h4.FaceColor = color1;
h5.FaceColor = color1;
h6.FaceColor = color1;



% title(t1, 'Rebound', 'FontSize',16)
t1.TileSpacing = 'tight';
t1.Padding = 'compact';
ylabel(t1,'Number of in silico individuals','FontSize',sizeT)

f1 = gcf;
f1.Position = [100 100 500 300];
exportgraphics(f1,'distribution_admissible_pars_rebound_and_nonrebound.png','Resolution',300)
exportgraphics(f1,'distribution_admissible_pars_rebound_and_nonrebound.pdf','Resolution',300)


%% plot a sample (Rebound)
figure(2); 
t2 = tiledlayout(1,1); 
nexttile; hold on; box on;
% title('Rebound')
sample_size = 100;%min(20,floor(0.05*length(Admissible_params)));
random_indexes = randi(length(Admissible_params(:,1)),sample_size);
% 
% xline(T_cutoff_min,'Color',[0.7 0.7 0.7],'LineWidth',4)
% xline(T_cutoff_max,'Color',[0.7 0.7 0.7],'LineWidth',4)

% xline(13,'--k','LineWidth',1.5) %adaptive immune

for jj = 1:sample_size
params = Admissible_params(random_indexes(jj),:);
[t,y] = ode45(@model,tspan,ICs,[],params);
plot(t,log10(y(:,5)),'Color',[0.7 0.7 0.7],'LineWidth',0.5)

end

% xline(T_cutoff_min,'Color',[0.7 0.7 0.7],'LineWidth',4)
% xline(T_cutoff_max,'Color',[0.7 0.7 0.7],'LineWidth',4)
% yline(V_cutoff,'--k','LineWidth',2)
% yline(4,'-.r')

yline(4,':','Color',[0.75 0.75 0.75 0.8],'LineWidth',2)

ylim([0 10])
% yticks([2 3 4 5 6 7 8 9 10])
xlim([0 21])
xticks([0 7 14 21])% 28])
ylabel({'Log10 viral copies/mL';'nasal swab'})
xlabel('Days post infection')

set(gca,'FontSize',9)

t2.TileSpacing = 'tight';
t2.Padding = 'compact';

f2 = gcf;
f2.Position = [100 100 300 240];
exportgraphics(f2,'samples_accepted_trajectories_rebound_and_nonrebound.png','Resolution',300)
exportgraphics(f2,'samples_accepted_trajectories_rebound_and_nonrebound.pdf','Resolution',300)

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