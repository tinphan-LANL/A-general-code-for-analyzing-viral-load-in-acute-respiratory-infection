close all
clear all
rng('default')
%% 5-day treatment

cd 'Treatment_on_day_4'
load('REBOUND_COUNT_4.mat')
cd ..
REBOUND_4 = REBOUND_HOLDER;

cd 'Treatment_on_day_5'
load('REBOUND_COUNT_5.mat')
cd ..
REBOUND_5 = REBOUND_HOLDER;

cd 'Treatment_on_day_6'
load('REBOUND_COUNT_6.mat')
cd ..
REBOUND_6 = REBOUND_HOLDER;

cd 'Treatment_on_day_7'
load('REBOUND_COUNT_7.mat')
cd ..
REBOUND_7 = REBOUND_HOLDER;

%% 6-day treatment

cd 'Treatment_on_day_4_6d'
load('REBOUND_COUNT_4.mat')
cd ..
REBOUND_4_6d = REBOUND_HOLDER;

cd 'Treatment_on_day_5_6d'
load('REBOUND_COUNT_5.mat')
cd ..
REBOUND_5_6d = REBOUND_HOLDER;

cd 'Treatment_on_day_6_6d'
load('REBOUND_COUNT_6.mat')
cd ..
REBOUND_6_6d = REBOUND_HOLDER;

cd 'Treatment_on_day_7_6d'
load('REBOUND_COUNT_7.mat')
cd ..
REBOUND_7_6d = REBOUND_HOLDER;

%% 7-day treatment

cd 'Treatment_on_day_4_7d'
load('REBOUND_COUNT_4.mat')
cd ..
REBOUND_4_7d = REBOUND_HOLDER;

cd 'Treatment_on_day_5_7d'
load('REBOUND_COUNT_5.mat')
cd ..
REBOUND_5_7d = REBOUND_HOLDER;

cd 'Treatment_on_day_6_7d'
load('REBOUND_COUNT_6.mat')
cd ..
REBOUND_6_7d = REBOUND_HOLDER;

cd 'Treatment_on_day_7_7d'
load('REBOUND_COUNT_7.mat')
cd ..
REBOUND_7_7d = REBOUND_HOLDER;

%% 8-day treatment

cd 'Treatment_on_day_4_8d'
load('REBOUND_COUNT_4.mat')
cd ..
REBOUND_4_8d = REBOUND_HOLDER;

cd 'Treatment_on_day_5_8d'
load('REBOUND_COUNT_5.mat')
cd ..
REBOUND_5_8d = REBOUND_HOLDER;

cd 'Treatment_on_day_6_8d'
load('REBOUND_COUNT_6.mat')
cd ..
REBOUND_6_8d = REBOUND_HOLDER;

cd 'Treatment_on_day_7_8d'
load('REBOUND_COUNT_7.mat')
cd ..
REBOUND_7_8d = REBOUND_HOLDER;

%% 10-day treatment

cd 'Treatment_on_day_4_10d'
load('REBOUND_COUNT_4.mat')
cd ..
REBOUND_4_10d = REBOUND_HOLDER;

cd 'Treatment_on_day_5_10d'
load('REBOUND_COUNT_5.mat')
cd ..
REBOUND_5_10d = REBOUND_HOLDER;

cd 'Treatment_on_day_6_10d'
load('REBOUND_COUNT_6.mat')
cd ..
REBOUND_6_10d = REBOUND_HOLDER;

cd 'Treatment_on_day_7_10d'
load('REBOUND_COUNT_7.mat')
cd ..
REBOUND_7_10d = REBOUND_HOLDER;

%%
figure(1); hold on; box on;

stepZ  = 0.15;
stepZZ = stepZ;

ONE     = ones(size(REBOUND_4));
ONEH    = (1+stepZ)*ONE;
ONEHH   = (1+stepZ+stepZZ)*ONE;
ONEHHH   = (1+stepZ+stepZZ+stepZZ)*ONE;
ONEHHHH   = (1+stepZ+stepZZ+stepZZ+stepZZ)*ONE;

TWO     = 2*ONE;
TWOH    = (2+stepZ)*ONE;
TWOHH   = (2+stepZ+stepZZ)*ONE;
TWOHHH   = (2+stepZ+stepZZ+stepZZ)*ONE;
TWOHHHH   = (2+stepZ+stepZZ+stepZZ+stepZZ)*ONE;

THREE   = 3*ONE;
THREEH  = (3+stepZ)*ONE;
THREEHH = (3+stepZ+stepZZ)*ONE;
THREEHHH = (3+stepZ+stepZZ+stepZZ)*ONE;
THREEHHHH = (3+stepZ+stepZZ+stepZZ+stepZZ)*ONE;

FOUR    = 4*ONE;
FOURH   = (4+stepZ)*ONE;
FOURHH  = (4+stepZ+stepZZ)*ONE;
FOURHHH  = (4+stepZ+stepZZ+stepZZ)*ONE;
FOURHHHH  = (4+stepZ+stepZZ+stepZZ+stepZZ)*ONE;

FIVE   = 5*ONE;
SIX    = 6*ONE;
SEVEN  = 7*ONE;

ymax = 25;
ymin = 0;
ylim([ymin ymax])

x_points = [0.5, 0.5, 2.7, 2.7];
y_points = [ymin, ymax, ymax, ymin];
color_sh = [1 1 1]; %the shading

a = fill(x_points, y_points, color_sh,'LineStyle','none');
a.FaceAlpha = 0.25;

h = boxchart([ONE, TWO, THREE,  FOUR],...
             [REBOUND_4, REBOUND_5,  REBOUND_6,  REBOUND_7]);
k = boxchart([ONEH, TWOH, THREEH, FOURH],...
             [REBOUND_4_6d, REBOUND_5_6d, REBOUND_6_6d, REBOUND_7_6d]);
l = boxchart([ONEHH, TWOHH, THREEHH, FOURHH],...
             [REBOUND_4_7d, REBOUND_5_7d, REBOUND_6_7d, REBOUND_7_7d]);
u = boxchart([ONEHHH, TWOHHH, THREEHHH, FOURHHH],...
             [REBOUND_4_8d, REBOUND_5_8d, REBOUND_6_8d, REBOUND_7_8d]);
w = boxchart([ONEHHHH, TWOHHHH, THREEHHHH, FOURHHHH],...
             [REBOUND_4_10d, REBOUND_5_10d, REBOUND_6_10d, REBOUND_7_10d]);

size_box = 0.1;

h.BoxFaceColor = [0.85 0.325 0.098];
h.MarkerStyle = 'none'; 
h.BoxFaceAlpha = 0.1;
h.LineWidth = 1.5;
h.BoxWidth = size_box;
h.WhiskerLineColor = [0.85 0.325 0.098];

k.BoxFaceColor = [86 70 160]/255;
k.MarkerStyle = 'none'; 
k.BoxFaceAlpha = 0.1;
k.LineWidth = 1.5;
k.BoxWidth = size_box;
k.WhiskerLineColor = [86 70 160]/255;

l.BoxFaceColor = [62 147 153]/255;
l.MarkerStyle = 'none'; 
l.BoxFaceAlpha = 0.1;
l.LineWidth = 1.5;
l.BoxWidth = size_box;
l.WhiskerLineColor = [62 147 153]/255;

u.BoxFaceColor = [150 150 150]/255;
u.MarkerStyle = 'none'; 
u.BoxFaceAlpha = 0.1;
u.LineWidth = 1.5;
u.BoxWidth = size_box;
u.WhiskerLineColor = [150 150 150]/255;

w.BoxFaceColor = [50 50 50]/255;
w.MarkerStyle = 'none'; 
w.BoxFaceAlpha = 0.1;
w.LineWidth = 1.5;
w.BoxWidth = size_box;
w.WhiskerLineColor = [50 50 50]/255;

axis tight

%% Data scatter 5-day
% plot treated 4
xCenter = 1;
spread  = 0.05;
scatter(rand(size(REBOUND_4))*spread - (spread/2) + xCenter,REBOUND_4,15,[0.714 0.273 0.08232],'LineWidth',0.5);

% plot treated 5
xCenter = 2;
scatter(rand(size(REBOUND_5))*spread - (spread/2) + xCenter,REBOUND_5,15,[0.714 0.273 0.08232],'LineWidth',0.5);

% plot treated 6
xCenter = 3;
scatter(rand(size(REBOUND_6))*spread - (spread/2) + xCenter,REBOUND_6,15,[0.714 0.273 0.08232],'LineWidth',0.5);

% plot treated 7
xCenter = 4;
scatter(rand(size(REBOUND_7))*spread - (spread/2) + xCenter,REBOUND_7,15,[0.714 0.273 0.08232],'LineWidth',0.5);

%% untreated with data
% plot untreated
% xCenter = 5;
% scatter(rand(size(REBOUND_100))*spread - (spread/2) + xCenter,REBOUND_100,15,[73 87 98]/255,'LineWidth',0.5);

% data
% 1. From this data set (extended/ongoing - Edelstein et al.), 3 rebound out of 44. (6.8%)
% 2. From Deo et al. Annals of Internal Med. "Symptom and viral rebound in untreated..."
%    13% (34) out of 261. Rebound (chosen def of high viral rebound - due to presencence of infectious): increase of >0.5log to >=5log
% 3. From Wong et al. Lancet Infectious Diseases. "Viral burden rebound in hospitalized patients"
%    4.5% (170 of 3787). Rebound defined as a reduction in >=3 ct.
% 4. Epic-HR trial (Anderson et al. nejm - nirmatrelvir-ritonavir and viral rebound in covid-19)
%    17 (1.7% of 980). Hidden in supplementary of Anderson et al. >= 0.5 log increase to >= 2.7 log10. Only use 3 time points (5, 10, 14).
%    The Edelstein et al. noted that using the same definition, only 3 out of 124 (2.4%) has VR detected and 13 of 16 (81.2%) rebound events were not captured.
% 5. Wong et al. Jama network open. "Incidence of viral rebound after..."
%    68 (0.6% of 11688). Def: Ct>40, then Ct<40.
% 6. Pandit et al. Clinical Infectious Diseases. "The coronavirus diseases 2019 rebound study..." - no access.
%    9.3% (out of 43). Definition unclear to due to no access.
% 7. Edelstein et al. (original)
%    1.8% (1 out of 55). Def: two folds - but for untreated - it's infectious culture detection following negative sample.
% 8. Dai et al. MedRxiv "Viral kinetics of severe acute..."
%    1 out of 25 (or 4%). Def: At least 2 neg Ct>=25, then at least 2 pos Ct<=25.
% 9. Hay et al. Elife. "Quantifying the impact of immune history ..."
%    Different criterias are compared: range (0-3%) of 749 to 1334 cohort.
%    Take 3% (40 out of 1334) - def 2+ VL >=30, then 2+ VL<30. Rationale: short/transient rebound is often observed in the other study. Ct 30 is lower than Ct25 ~7log.
% X_vect = xCenter+spread+linspace(0.05,0.5,9);
% Untreated_rebound_data = [6.8, 13.0, 4.5, 1.7, 0.6, 9.3, 1.8, 4, 3];
% x_loc = (5+2*stepZ)*ones(size(Untreated_rebound_data));

% uD = boxchart(x_loc,Untreated_rebound_data);

% uD.BoxFaceColor = [86 86 86]/255;
% uD.MarkerStyle = 'none'; 
% uD.BoxFaceAlpha = 0.1;
% uD.LineWidth = 1.5;
% uD.BoxWidth = stepZ - 0.05;
% uD.WhiskerLineColor = [86 86 86]/255;

% xCenter2 = xCenter+2*stepZ;
% scatter(rand(size(Untreated_rebound_data))*spread - (spread/2) + xCenter2,Untreated_rebound_data,15,[46 46 46]/255,'v','LineWidth',0.5);
% 
% % add p-score via ranksum test
% p_score = ranksum(REBOUND_100,Untreated_rebound_data); %calculate p-score
% p_score = round(p_score,2);
% 
% plot([xCenter xCenter2], [15 15],'k','LineWidth',1.5);
% plot([xCenter xCenter], [14.5 15],'k','LineWidth',1.5);
% plot([xCenter2 xCenter2], [14.5 15],'k','LineWidth',1.5);
% text(xCenter+0.075,16,strcat('P =',{' '},num2str(p_score)));

%% Data scatter 6-day
% plot treated 4
xCenter = 1+stepZ;
scatter(rand(size(REBOUND_4_6d))*spread - (spread/2) + xCenter,REBOUND_4_6d,15,[71 63 158]/255,'LineWidth',0.5);

% plot treated 5
xCenter = 2+stepZ;
scatter(rand(size(REBOUND_5_6d))*spread - (spread/2) + xCenter,REBOUND_5_6d,15,[71 63 158]/255,'LineWidth',0.5);

% plot treated 6
xCenter = 3+stepZ;
scatter(rand(size(REBOUND_6_6d))*spread - (spread/2) + xCenter,REBOUND_6_6d,15,[71 63 158]/255,'LineWidth',0.5);

% plot treated 7
xCenter = 4+stepZ;
scatter(rand(size(REBOUND_7_6d))*spread - (spread/2) + xCenter,REBOUND_7_6d,15,[71 63 158]/255,'LineWidth',0.5);

%% Data scatter 7-day
% plot treated 4
xCenter = 1+stepZ+stepZZ;
scatter(rand(size(REBOUND_4_7d))*spread - (spread/2) + xCenter,REBOUND_4_7d,15,[88 150 81]/255,'LineWidth',0.5);

% plot treated 5
xCenter = 2+stepZ+stepZZ;
scatter(rand(size(REBOUND_5_7d))*spread - (spread/2) + xCenter,REBOUND_5_7d,15,[88 150 81]/255,'LineWidth',0.5);

% plot treated 6
xCenter = 3+stepZ+stepZZ;
scatter(rand(size(REBOUND_6_7d))*spread - (spread/2) + xCenter,REBOUND_6_7d,15,[88 150 81]/255,'LineWidth',0.5);

% plot treated 7
xCenter = 4+stepZ+stepZZ;
scatter(rand(size(REBOUND_7_7d))*spread - (spread/2) + xCenter,REBOUND_7_7d,15,[88 150 81]/255,'LineWidth',0.5);

%% Data scatter 8-day
% plot treated 4
xCenter = 1+stepZ+stepZZ+stepZZ;
scatter(rand(size(REBOUND_4_8d))*spread - (spread/2) + xCenter,REBOUND_4_8d,15,[135 135 135]/255,'LineWidth',0.5);

% plot treated 5
xCenter = 2+stepZ+stepZZ+stepZZ;
scatter(rand(size(REBOUND_5_8d))*spread - (spread/2) + xCenter,REBOUND_5_8d,15,[135 135 135]/255,'LineWidth',0.5);

% plot treated 6
xCenter = 3+stepZ+stepZZ+stepZZ;
scatter(rand(size(REBOUND_6_8d))*spread - (spread/2) + xCenter,REBOUND_6_8d,15,[135 135 135]/255,'LineWidth',0.5);

% plot treated 7
xCenter = 4+stepZ+stepZZ+stepZZ;
scatter(rand(size(REBOUND_7_8d))*spread - (spread/2) + xCenter,REBOUND_7_8d,15,[135 135 135]/255,'LineWidth',0.5);

%% Data scatter 10-day
% plot treated 4
xCenter = 1+stepZ+stepZZ+stepZZ+stepZZ;
scatter(rand(size(REBOUND_4_10d))*spread - (spread/2) + xCenter,REBOUND_4_10d,15,[35 35 35]/255,'LineWidth',0.5);

% plot treated 5
xCenter = 2+stepZ+stepZZ+stepZZ+stepZZ;
scatter(rand(size(REBOUND_5_10d))*spread - (spread/2) + xCenter,REBOUND_5_10d,15,[35 35 35]/255,'LineWidth',0.5);

% plot treated 6
xCenter = 3+stepZ+stepZZ+stepZZ+stepZZ;
scatter(rand(size(REBOUND_6_10d))*spread - (spread/2) + xCenter,REBOUND_6_10d,15,[35 35 35]/255,'LineWidth',0.5);

% plot treated 7
xCenter = 4+stepZ+stepZZ+stepZZ+stepZZ;
scatter(rand(size(REBOUND_7_10d))*spread - (spread/2) + xCenter,REBOUND_7_10d,15,[35 35 35]/255,'LineWidth',0.5);

%%
scaleZ = 0.2;
xlim([0.5+scaleZ 5-scaleZ])
xticks([1 2 3 4])
set(gca,'XTickLabel',{'Day 1','Day 2','Day 3','Day 4'})

ylim([ymin ymax])
% set(b,'linew',1.5)
ylabel('Percent rebound')
xlabel('Day treatment initiated relative to symptom onset')

set(gca,'Fontsize',14)
set(gca,'LineWidth',1)

f = gcf;
f.Position = [100 100 1750 450];
exportgraphics(f,'percent_of_rebound_so3.png','Resolution',300)
exportgraphics(f,'percent_of_rebound_so3.pdf','Resolution',300)

