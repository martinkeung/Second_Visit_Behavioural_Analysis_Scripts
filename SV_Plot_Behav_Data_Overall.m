function SV_Plot_Behav_Data_Overall

% ---------------------------------------------------------------------------------------------------------------------------------------%
% INPUT ARGUMENTS:
% No inputs
%
% OUTPUTS:
% No outputs. Only makes plots. Can be modified to output variables of
% choice
% ---------------------------------------------------------------------------------------------------------------------------------------%

% Path to all the data
datadir = 'M:\Data_Masterfile\H20-00572_All-Dressed\AllDressed_WorkOnData\All-Dressed_Second_Visit';
addpath(genpath('M:\Data_Masterfile\H20-00572_All-Dressed\AllDressed_WorkOnData\All-Dressed_Second_Visit'));
addpath(genpath('M:\H20-00572_All-Dressed\Second Visit\Analysis_Scripts'));

cd 'M:\Data_Masterfile\H20-00572_All-Dressed\AllDressed_WorkOnData\All-Dressed_Second_Visit'

gvslab = {'Pink','Sham','Best GVS'};

d = dir('PSHC*');

for ii = 1:size(d,1)
    dirflags(ii,1) = d(ii).isdir;
end
d(~dirflags) = [];
nsubj = size(d,1);

% Initialize variables
ngvs = 3; % number of GVS

Mvigour = zeros(nsubj,ngvs);
SEvigour = zeros(nsubj,ngvs);

Mpt = zeros(nsubj,ngvs);
SEpt = zeros(nsubj,ngvs);

Mpres = zeros(nsubj,ngvs);
SEforce = zeros(nsubj,ngvs);
sqtval = 0.3;

for s = 1:size(d,1)
    
    participant_id = d(s).name;
    disp(participant_id);
    behavdir = sprintf('%s/%s',datadir,participant_id);
    
    % Load data for specific run
    filename = sprintf('%s_SV_Run_1.mat',participant_id);
    if exist(filename,'file')
        input_ds = load(filename);
    end
    res = create_ds_valid_v1_sv(input_ds,sqtval);
    
    % Sort behaviour data by GVS [GVS1 (Sham), GVS2,...,GVS9]
    res = sortbyGVS_valid_sv(res);
    
    
    %% Analyse MID task data (Effect of GVS)
    
    for g = 1:ngvs
    
        % Trials to keep: all non-negative reward trials with successful
        % squeeze
        Kp = (res.reward(:,g) >=0 & res.goodtrials(:,g) == 1);
    
        % Pak time per GVS
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
        Mpt(s,g) = nanmean(pt);
        SEpt(s,g) = nanstd(pt)/sqrt(length(pt));
    
        % Mean and standard error of peak pressure
        Mpres(s,g) = nanmean(pres);
        SEforce(s,g) = nanstd(pres)/sqrt(length(pres));
    
        % Mean and standard error of vigour
        Mvigour(s,g) = nanmean(vigour);
        SEvigour(s,g) = nanstd(vigour)/sqrt(length(vigour));
    
    end
    
    % Compute min and max values for setting axis limits
    vmin(s,:) = min(Mvigour(s,:) - SEvigour(s,:));
    vmax(s,:) = max(Mvigour(s,:) + SEvigour(s,:));
    vmin(isnan(vmin)) = 0; vmax(isnan(vmax)) = 1;
    
    pmin(s,:) = min(Mpt(s,:) - SEpt(s,:));
    pmax(s,:) = max(Mpt(s,:) + SEpt(s,:));
    pmin(isnan(pmin)) = 0; pmax(isnan(pmax)) = 1;
    
    fmin(s,:) = min(Mpres(s,:) - SEforce(s,:));
    fmax(s,:) = max(Mpres(s,:) + SEforce(s,:));
    fmin(isnan(fmin)) = 0; fmax(isnan(fmax)) = 1;

end

%% Let's make some plots

% Note: tiledlayout parameters have been initialized to accomodate 12
% subjects. Edit as needed.

% Initialize colormap
nr = bluewhitered(400);
nb = nr;
nb(:,3) = nr(:,1); nb(:,1) = nr(:,2);
nb = nb(80:40:400,:);
close;

% Plot vigour data
figure;
tiledlayout(4,5); % adjust as needed, this will fit 20 participants
set(gcf, 'Position', get(0,'Screensize'));

for s = 1:size(d,1)
    nexttile
    bar(1, Mvigour(s,1),'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on;
    errorbar(1,Mvigour(s,1),SEvigour(s,1),'Color','k','Linewidth',1);
    for gvs = 2:ngvs
        bar(gvs,Mvigour(s,gvs),'FaceColor',nb(gvs-1,:),'EdgeColor',nb(gvs-1,:)); hold on;
        errorbar(gvs,Mvigour(s,gvs),SEvigour(s,gvs),'Color','k','Linewidth',1);
    end

    xlim([0,ngvs+1]);
    ylim([vmin(s,:) vmax(s,:)]);
    ylabel('motor vigour (V/s)','Fontweight','bold');
    set(gca,'XTick',1:ngvs,'XTickLabel',gvslab);
    box off;
    title(sprintf('%s',d(s).name),'FontSize',12);

end

sgtitle(sprintf('Motor vigour'),'FontSize',16,'Fontweight','bold');

figname = sprintf('%s/motor_vigour_allsubj', datadir);
print(gcf,figname,'-dpng','-r300');
% close;


% Plot peak time data
figure;
tiledlayout(4,5);
set(gcf, 'Position', get(0,'Screensize'));

for s = 1:size(d,1)

    % Plot GVS and sham behaviour
    nexttile;
    if strcmp(d(s).name,'PSPD004') && strcmp(med_status,'ON')
        set(gca,'Visible','off');
        continue;
    end

    bar(1, Mpt(s,1),'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on;
    errorbar(1,Mpt(s,1),SEpt(s,1),'Color','k','Linewidth',1);
    for gvs = 2:ngvs
        bar(gvs,Mpt(s,gvs),'FaceColor',nb(gvs-1,:),'EdgeColor',nb(gvs-1,:)); hold on;
        errorbar(gvs,Mpt(s,gvs),SEpt(s,gvs),'Color','k','Linewidth',1);
    end

    xlim([0,ngvs+1]);
    ylim([pmin(s,:) pmax(s,:)]);
    ylabel('peak time (s)','Fontweight','bold');
    set(gca,'XTick',1:ngvs,'XTickLabel',gvslab);
    box off;
    title(sprintf('%s',d(s).name),'FontSize',12);

end

sgtitle(sprintf('Peak time'),'FontSize',16,'Fontweight','bold');

figname = sprintf('%s/pt_allsubj', datadir);
print(gcf,figname,'-dpng','-r300');
% close;

% plot pressure data
figure;
tiledlayout(4,5);
set(gcf, 'Position', get(0,'Screensize'));

for s = 1:size(d,1)

    nexttile
    bar(1, Mpres(s,1),'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on;
    errorbar(1,Mpres(s,1),SEforce(s,1),'Color','k','Linewidth',1);
    for gvs = 2:ngvs
        bar(gvs,Mpres(s,gvs),'FaceColor',nb(gvs-1,:),'EdgeColor',nb(gvs-1,:)); hold on;
        errorbar(gvs,Mpres(s,gvs),SEforce(s,gvs),'Color','k','Linewidth',1);
    end

    xlim([0,ngvs+1]);
    ylim([fmin(s,:) fmax(s,:)]);
    ylabel('peak pressure (V)','Fontweight','bold');
    set(gca,'XTick',1:ngvs,'XTickLabel',gvslab);
    box off;
    title(sprintf('%s',d(s).name),'FontSize',12);

end

sgtitle(sprintf('Peak pressure'),'FontSize',16,'Fontweight','bold');

figname = sprintf('%s/peak_pressure_allsubj', datadir);
print(gcf,figname,'-dpng','-r300');
% close;


end
