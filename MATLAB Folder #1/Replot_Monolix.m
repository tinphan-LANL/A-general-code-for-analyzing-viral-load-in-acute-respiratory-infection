% This one runs the code for fixed n=1
clear all 
close all 
%% Define some constants
LOD = 1;

%% Data
cd .. ; cd 'Monolix fit'
T = readtable('paxlovid_processed_data_nov_29_2023.csv');
cd ..; cd 'MATLAB Folder #1'

VL_all   = T{:,3};  %viral load log10/mL
Time_all = T{:,2};  %time post symptom onset
Censor   = T{:,4};  %censoring
Pax_all  = T{:,5};  %pax time post symptom
PID      = T{:,1};  %patient ID
Group    = T{:,6};  %Group 1-4
DEF      = T{:,7};  %clinical classification of rebound

[C, ia, ic] = unique(PID,'rows','stable'); %ia is the index of unique patient id
Num_patients = length(C);

%% Best fit parameters
cd ..; cd 'Monolix fit/Model_fit/IndividualParameters'
P = readtable('estimatedIndividualParameters.txt');
cd ..; cd ..; cd ..; cd 'MATLAB Folder #1'

P_all    = P{1:end,18:25};

P_best_fit = nan(size(P_all));

Num_fitted_patients = min(Num_patients, length(P_all(:,1)));

V_PAX = nan(2,Num_fitted_patients); % store VL at initiation and end of Paxlovid.

%% Loop for plots
plot_width  = 9; %number of figures in a row
plot_height = 6; %number of figures in a columb
plot_size   = plot_width*plot_height;
num_plots   = ceil(Num_fitted_patients/plot_size);

for kk = 1:num_plots

    figure(kk);
    t1 = tiledlayout(plot_height,plot_width);

    for ii=1+plot_size*(kk-1): min(plot_size*kk,Num_fitted_patients)

        if ii<Num_patients
            VL   = VL_all(ia(ii):ia(ii+1)-1);
            Time = Time_all(ia(ii):ia(ii+1)-1);
            pid  = PID(ia);
            def  = DEF(ia);
            pax  = Pax_all(ia(ii):ia(ii+1)-1);
            pax  = rmmissing(pax); %remove nan

        else
            VL   = VL_all(ia(ii):end);
            Time = Time_all(ia(ii):end);
            pid  = PID(ia);
            def  = DEF(ia);
            pax  = Pax_all(ia(ii):end);
            pax  = rmmissing(pax); %remove nan
        end
    
        % Collect patient-specific best fit parameters
            beta_exp    = P_all(ii,1);
                beta = 10^(-beta_exp);
            p_exp       = P_all(ii,2);
                p = 10^p_exp;
            rho         = P_all(ii,3); 
            phi_exp     = P_all(ii,4);
                phi = 10^(-phi_exp);
            delta0      = P_all(ii,5);
            t00         = P_all(ii,6);
            TT          = P_all(ii,7);
            I_rho_exp   = P_all(ii,8);
                I_rho = 10^(I_rho_exp);
            n           = 1; %P_all(ii,9);

        % Shift time of infection according to TT
            Time = Time + TT;
            ti   = pax(1) + TT; %time starting Pax
            t00  = t00+TT;

        % Parameter vector
            param = [beta, p, rho, phi, delta0, ti, t00, I_rho, n];

        % Time vector
            x_axis_upper_limit = 28; %day post infection
            t_f = max(Time(end),x_axis_upper_limit);
            n = t_f*100+1;
            t_span = linspace(0,t_f,n);

        % Initialization
            T0 = 8*10^7;
            R0 = 0;
            E0 = 1;
            I0 = 0;
            V0 = 0;

            X_int = [T0, R0, E0, I0, V0];

        % Solve ODE
            opts  = odeset('RelTol',1e-10,'AbsTol',1e-10);
            [t y] = ode45(@model,t_span,X_int,opts,param);

        % Assign variables
            T = y(:,1);
            R = y(:,2);
            E = y(:,3);
            I = y(:,4);
            V = log10(y(:,5));

        % Plotting data
            nexttile; hold on; box on;
        
            below_indexes = find(VL == 1);
            above_indexes = find(VL >  1);

            if strcmp(def{ii}, {'No Viral Rebound'}) ~= 1  %rebound
                scatter(Time(below_indexes), VL(below_indexes),60,([228 147 179]-30)/255,'LineWidth',1.5)
                scatter(Time(above_indexes), VL(above_indexes),80,([228 147 179]-30)/255,'filled','LineWidth',1.5)
                yline(LOD,':','Color',([228 147 179]-35)/255,'LineWidth',1.5)
            else
                scatter(Time(below_indexes), VL(below_indexes),60,[0 0.447 0.741],'LineWidth',1.5)
                scatter(Time(above_indexes), VL(above_indexes),80,[0 0.447 0.741],'filled','LineWidth',1.5)
                yline(LOD,':','Color',[0 0.4 0.7],'LineWidth',1.5)
            end

            if ti<100
                x_points = [ti, ti, ti+5, ti+5];
                y_points = [0, 10, 10, 0];
                color_sh = [0.75 0.75 0.75];

                a = fill(x_points, y_points, color_sh,'LineStyle','none');
                a.FaceAlpha = 0.25;
            end

                xline(TT,':','Color',[0 0 0],'LineWidth',1.5)

        % Plotting model solution
            if strcmp(def{ii}, {'No Viral Rebound'}) ~= 1  %rebound
                plot(t,V,'Color',[228 147 179]/255,'LineWidth',2.5)
            else
                plot(t,V,'Color',[0 0.5 0.7],'LineWidth',2.5)
            end
            set(gca,'FontSize',12)

        % title and axis label
        if strcmp(def{ii}, {'No Viral Rebound'}) == 1  %no rebound

                title(strcat('PID',{' '},num2str(ii)),'FontSize',10,'Color',[0 0.5 0.7]) %
        else %rebound
                title(strcat('PID',{' '},num2str(ii)),'FontSize',10,'Color',([228 147 179]-35)/255)
        end

        xlim([0 x_axis_upper_limit])
        xticks([0 7 14 21 28])
        ylim([0 10])
%         aaaa = 1;

    end

    xlabel(t1,'Days post infection','FontSize',16)
    ylabel(t1,{'Log10 viral copies/mL';'nasal swab'},'FontSize',16)

    f = gcf;
    f.WindowState = 'maximize';
    exportgraphics(f,strcat('Fig_2a.png'),'Resolution',300)
    exportgraphics(f,strcat('Fig_2a.pdf'),'Resolution',600)

end

% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boxplot
% Recall there is 20 in group 1, 31 in group 2.
% [beta, p, rho, phi, delta0, TT, t00, I_rho, n]
rng('default')
Group1 = P_all(1:20,:);
Group2 = P_all(21:21+30,:);

Pax_G1 = Pax_all(ia(1:20));
Pax_G2 = Pax_all(ia(21:51));

Group1(:,9) = Pax_G1 + Group1(:,7); %pax from time of infection
Group2(:,9) = Pax_G2 + Group2(:,7);

Group1(:,10) = Pax_G1; %pax from time of symptom onset
Group2(:,10) = Pax_G2;

Holder1 = Group1(:,8);
Holder2 = Group2(:,8);

Group1(:,8) = Group1(:,7);
Group2(:,8) = Group2(:,7);

Group1(:,7) = Holder1;
Group2(:,7) = Holder2;

namess = {'- log_{10} \beta', 'log_{10} p', '\rho', '\phi', '\delta_0', 't^* (dpi)',...
          'log_{10} K_\rho', {'Days from symptom onset';'to infection time'},...
          {'Days from infection'; 'to N-R initiation'}, {'Days from symptom onset'; 'to N-R initiation'},...
          {'Days from end of N-R to onset';'of adaptive immune response'}};

Group1(:,4) = 10.^(-Group1(:,4));
Group2(:,4) = 10.^(-Group2(:,4));

Group1(:,11) = Group1(:,6)-Group1(:,9)-4.5;
Group2(:,11) = Group2(:,6)-Group2(:,9)-4.5;

figure(kk+1);
    t2 = tiledlayout(2,6);

for ii=1:length(Group1(1,:))

nexttile; hold on; box on;

X1 = Group1(:,ii)';
X2 = Group2(:,ii)';

X   = [X1 X2];
grp = [0.5.*ones(1,length(X1)), 1.5.*ones(1,length(X2))];

h = boxchart(grp,X);
h.BoxFaceColor = [0.85 0.325 0.098];
h.MarkerStyle = 'none'; %remove outliers plot
h.BoxFaceAlpha = 0.1;
h.LineWidth = 1;
h.WhiskerLineColor = [0.85 0.325 0.098];
axis tight
xlim([0 2])

% plot data 1
xCenter = 0.5;
spread  = 0.3;
scatter(rand(size(X1))*spread - (spread/2) + xCenter,X1,15,[0.714 0.273 0.08232],'LineWidth',0.5);

% plot data 2
xCenter = 1.5;
spread  = 0.3;
scatter(rand(size(X2))*spread - (spread/2) + xCenter,X2,15,[0.714 0.273 0.08232],'LineWidth',0.5);

xticks([0.5 1.5])

set(gca,'XTickLabel',{'Rebound','Nonrebound'},'FontSize',9)

if ii<8
    ylabel(namess{ii},'FontSize',12)
else
    ylabel(namess{ii},'FontSize',9)
end

p_score = ranksum(X1,X2); %calculate p-score
p_score = round(p_score,4);

yt = get(gca, 'YTick');
xt = get(gca, 'XTick');

% title(strcat('p-score',{' '},num2str(p_score)),'FontSize',10)
title(strcat('p-value',{' '},num2str(p_score)),'FontSize',10)

end

f = gcf; 
% f.WindowState = 'maximize';
f.Position = [50 50 1500 500];
exportgraphics(f,'Fig_2b.png','Resolution',300)
exportgraphics(f,'Fig_2b.pdf','Resolution',600)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function dY = model(t,y,param)
    beta   = param(1);
    p      = param(2);
    rho    = param(3);
    phi    = param(4);
    delta0 = param(5);

    ti     = param(6);
    te     = ti+4.5;

    t00    = param(7);
    I_rho  = param(8);
    n      = param(9);

    k      = 4;
    c      = 10;
    sigma  = 0.5;
    EC50   = 62;
    eps    = 1;
    
    %adaptive immunity

    delta_max = 10*delta0;

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