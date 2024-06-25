function SV_Behavioural_Analysis_Rewardwise(participant_id)

% ---------------------------------------------------------------------------------------------------------------------------------------%
% INPUT ARGUMENTS:
% participant_id: 'PSPDXXX'
%
% OUTPUTS:
% No putputs. Only makes plots. Can be modified to output variables of
% choice
% ---------------------------------------------------------------------------------------------------------------------------------------%

% Path to the data

% participant_id = 'PSHC95999' % hard code for now
datadir = 'M:\Data_Masterfile\H20-00572_All-Dressed\AllDressed_WorkOnData\All-Dressed_Second_Visit';
behavdir = sprintf('%s\\%s',datadir,participant_id);
gvslab = {'Pink','Sham','Best GVS'};
sqtval = 0.3; % Set as this for now 
lowreward = 1;
highreward =5;

% If analyzing each run independently

% Load data for first run
filename = sprintf('%s_SV_Run_1.mat',participant_id);
if exist(filename,'file')
    input_ds = load(filename);
    disp('Data loaded')
end
res = create_ds_valid_v1_sv(input_ds,sqtval);
% Sort behaviour data by GVS [GVS1 (Sham), GVS2,...,GVS9]
resall = sortbyGVS_valid_sv(res);


%% Analyse MID task data (Effect of GVS)

ngvs = 3; % number of GVS
reward_level = [1, 5];
nreward = length(reward_level);


% Initialize variables
meanrewardPT = zeros(nreward, ngvs);
meanrewardPres = zeros(nreward, ngvs);
meanrewardVigour = zeros(nreward, ngvs);

serewardPT = zeros(nreward, ngvs);
serewardForce = zeros(nreward, ngvs);
serewardVigour = zeros(nreward, ngvs);

res = resall;

% For each GVS...
for g = 1:ngvs

    % For eah reward...
    for i = 1:nreward

        % Trials to keep: all non-negative reward trials with successful
        % squeeze
        Kp = (res.reward(:,g) == reward_level(i) & res.goodtrials(:,g) == 1);

        % Peak time per GVS
        pt = res.peak_time(:,g);
        pt = pt(Kp);

        % Peak pressure per GVS
        pres = res.peak_pressure(:,g);
        pres = pres(Kp);

        % Remove trials where peak time is nan
        pres(isnan(pt)) = [];
        pt(isnan(pt)) = [];
        vigour = pres./pt;

        % Mean and standard error of peak times
        meanrewardPT(i,g) = mean(pt);
        serewardPT(i,g) = std(pt)/sqrt(length(pt));

        % Mean and standard error of peak pressure
        meanrewardPres(i,g) = mean(pres);
        serewardForce(i,g) = std(pres)/sqrt(length(pres));

        % Mean and standard error of vigour
        meanrewardVigour(i,g) = mean(vigour);
        serewardVigour(i,g) = std(vigour)/sqrt(length(vigour));
    end
end

% Compute min and max values for setting axis limits
Mpos1Vig = meanrewardVigour(1,:); % reward $1 or $0
SEpos1Vig = serewardVigour(1,:); % reward $1 or $0

Mpos5Vig = meanrewardVigour(2,:); % reward $5
SEpos5Vig = serewardVigour(2,:); % reward $5

vmin = min([Mpos1Vig-SEpos1Vig,Mpos5Vig-SEpos5Vig]);
vmax = max([Mpos1Vig+SEpos1Vig,Mpos5Vig+SEpos5Vig]);

Mpos1PT = meanrewardPT(1,:); % reward $1 or $0
SEpos1PT = serewardPT(1,:); % reward $1 or $0

Mpos5PT = meanrewardPT(2,:); % reward $5
SEpos5PT = serewardPT(2,:); % reward $5

pmin = min([Mpos1PT-SEpos1PT,Mpos5PT-SEpos5PT]);
pmax = max([Mpos1PT+SEpos1PT,Mpos5PT+SEpos5PT]);

Mpos1Pres = meanrewardPres(1,:); % reward $1 or $0
SEpos1Pres = serewardForce(1,:); % reward $1 or $0

Mpos5Pres = meanrewardPres(2,:); % reward $5
SEpos5Pres = serewardForce(2,:); % reward $5

fmin = min([Mpos1Pres-SEpos1Pres,Mpos5Pres-SEpos5Pres]);
fmax = max([Mpos1Pres+SEpos1Pres,Mpos5Pres+SEpos5Pres]);


%% Lets make some plots

% Initialize colormap
nr = bluewhitered(400);
nb = nr;
nb(:,3) = nr(:,1); nb(:,1) = nr(:,2);
nb = nb(80:40:400,:);

figure(100);
tiledlayout(nreward,1);
set(gcf, 'Position', [0 0 500 700]);

% Plot vigour data

% Plot GVS and sham behaviour for $0/$1 reward first
nexttile;

bar(1, Mpos1Vig(1),'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on;
errorbar(1,Mpos1Vig(1),SEpos1Vig(1),'Color','k','Linewidth',1);
for gvs = 2:ngvs
    bar(gvs,Mpos1Vig(gvs),'FaceColor',nb(gvs-1,:),'EdgeColor',nb(gvs-1,:)); hold on;
    errorbar(gvs,Mpos1Vig(gvs),SEpos1Vig(gvs),'Color','k','Linewidth',1);
end

xlim([0,ngvs+1]); ylim([vmin vmax]);
ylabel('motor vigour (V/s)','Fontweight','bold');
set(gca,'XTick',1:ngvs,'XTickLabel',gvslab);
box off;
title(sprintf('Reward $%d',lowreward),'FontSize',12);

% Plot GVS and sham behaviour for $5 reward next
nexttile;

bar(1, Mpos5Vig(1),'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on;
errorbar(1,Mpos5Vig(1),SEpos5Vig(1),'Color','k','Linewidth',1);
for gvs = 2:ngvs
    bar(gvs,Mpos5Vig(gvs),'FaceColor',nb(gvs-1,:),'EdgeColor',nb(gvs-1,:)); hold on;
    errorbar(gvs,Mpos5Vig(gvs),SEpos5Vig(gvs),'Color','k','Linewidth',1);
end

xlim([0,ngvs+1]); ylim([vmin vmax]);
ylabel('motor vigour (V/s)','Fontweight','bold');
set(gca,'XTick',1:ngvs,'XTickLabel',gvslab);
box off;
title(sprintf('Reward $%d',highreward),'FontSize',12);

sgtitle(sprintf('%s',participant_id),'FontSize',16,'Fontweight','bold');

% save figure
figname = sprintf('%s/motor_vigour_rewardwise_%s', behavdir, participant_id);
print(gcf,figname,'-dpng','-r300');
% close;


% Plot peak time data
figure(200);
tiledlayout(nreward,1);
set(gcf, 'Position', [0 0 500 700]);

% plot GVS and sham behaviour for $0/$1 reward first
nexttile;

bar(1, Mpos1PT(1),'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on;
errorbar(1,Mpos1PT(1),SEpos1PT(1),'Color','k','Linewidth',1);
for gvs = 2:ngvs
    bar(gvs,Mpos1PT(gvs),'FaceColor',nb(gvs-1,:),'EdgeColor',nb(gvs-1,:)); hold on;
    errorbar(gvs,Mpos1PT(gvs),SEpos1PT(gvs),'Color','k','Linewidth',1);
end

xlim([0,ngvs+1]); ylim([pmin pmax]);
ylabel('peak time (s)','Fontweight','bold');
set(gca,'XTick',1:ngvs,'XTickLabel',gvslab);
box off;
title(sprintf('Reward $%d',lowreward),'FontSize',12);

% plot GVS and sham behaviour for $5 reward next
nexttile;

bar(1, Mpos5PT(1),'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on;
errorbar(1,Mpos5PT(1),SEpos5PT(1),'Color','k','Linewidth',1);
for gvs = 2:ngvs
    bar(gvs,Mpos5PT(gvs),'FaceColor',nb(gvs-1,:),'EdgeColor',nb(gvs-1,:)); hold on;
    errorbar(gvs,Mpos5PT(gvs),SEpos5PT(gvs),'Color','k','Linewidth',1);
end

xlim([0,ngvs+1]); ylim([pmin pmax]);
ylabel('peak time (s)','Fontweight','bold');
set(gca,'XTick',1:ngvs,'XTickLabel',gvslab);
box off;
title(sprintf('Reward $%d',highreward),'FontSize',12);

sgtitle(sprintf('%s',participant_id),'FontSize',16,'Fontweight','bold');

% save figure
figname = sprintf('%s/pt_rewardwise_%s', behavdir,participant_id);
print(gcf,figname,'-dpng','-r300');
% close;


% Plot peak pressure data
figure(300);
tiledlayout(nreward,1);
set(gcf, 'Position', [0 0 500 700]);

% plot GVS and sham behaviour for $0/$1 reward first
nexttile;

bar(1, Mpos1Pres(1),'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on;
errorbar(1,Mpos1Pres(1),SEpos1Pres(1),'Color','k','Linewidth',1);
for gvs = 2:ngvs
    bar(gvs,Mpos1Pres(gvs),'FaceColor',nb(gvs-1,:),'EdgeColor',nb(gvs-1,:)); hold on;
    errorbar(gvs,Mpos1Pres(gvs),SEpos1Pres(gvs),'Color','k','Linewidth',1);
end

xlim([0,ngvs+1]); ylim([fmin fmax]);
ylabel('peak pressure (V)','Fontweight','bold');
set(gca,'XTick',1:ngvs,'XTickLabel',gvslab);
box off;
title(sprintf('Reward $%d',lowreward),'FontSize',12);

% plot GVS and sham behaviour for $5 reward next
nexttile;

bar(1, Mpos5Pres(1),'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on;
errorbar(1,Mpos5Pres(1),SEpos5Pres(1),'Color','k','Linewidth',1);
for gvs = 2:ngvs
    bar(gvs,Mpos5Pres(gvs),'FaceColor',nb(gvs-1,:),'EdgeColor',nb(gvs-1,:)); hold on;
    errorbar(gvs,Mpos5Pres(gvs),SEpos5Pres(gvs),'Color','k','Linewidth',1);
end
xlim([0,ngvs+1]); ylim([fmin fmax]);
ylabel('peak pressure (V)','Fontweight','bold');
set(gca,'XTick',1:ngvs,'XTickLabel',gvslab);
box off;
title(sprintf('Reward $%d',highreward),'FontSize',12);

sgtitle(sprintf('%s',participant_id),'FontSize',16,'Fontweight','bold');

% save figure
figname = sprintf('%s/peak_pressure_rewardwise_%s', behavdir,participant_id);
print(gcf,figname,'-dpng','-r300');
% close;
end