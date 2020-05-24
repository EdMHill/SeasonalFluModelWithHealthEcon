%Purpose:
%Function to plot temporal profiles (of influenza infection prevalence)
% from en masse model simulations

%Used with age-structured models!
%--------------------------------------------------------------------------

function InfTemporalProfile_FMAltMultiRuns(InputFileName,SaveFileName,SimnRunType,varargin)
%Inputs:
%   InputFileName - (string) File containing required model simulation outputs
%   SaveFileName - (string) Specify location to save plots
%   SimnRunType - (scalar) Specifies whether run was inference (fitting to historical data), or forward simulation
%                    1 -> historical; 2 -> forward simn
%   varargin  - input to contain number of additional seasons
load(InputFileName,'SimnData')

%Load cells containing infected timeseries
InfTimeProfileAllRuns = SimnData(:,2);
PopnDistProfileAllRuns = SimnData(:,3);

%% Amass runs into single array

%Initialise storage array
%  -> Row per run
%  -> column per timestep
RunNum = numel(InfTimeProfileAllRuns);
LengthT = length(InfTimeProfileAllRuns{1});
InfTimeProfileAllRunsArray = zeros(RunNum,LengthT,101);
PopnDistAllRunsArray = zeros(RunNum,LengthT,101);

%Iterate through each individual run, add timeseries info to InfTimeProfileAllRunsArray
for ii = 1:RunNum
    InfTimeProfileAllRunsArray(ii,:,:) = InfTimeProfileAllRuns{ii};
    PopnDistAllRunsArray(ii,:,:) = PopnDistProfileAllRuns{ii};
end

%Get proportion OF EACH AGE GROUP that is infected at each time step.
InfAsPropnAgeGrpAllRunsArray = InfTimeProfileAllRunsArray./PopnDistAllRunsArray;

%% Plot all timeseries on single plot
figure('Color',[1 1 1]);
clf
position = [100, 100, 2*550, 450];
set(0, 'DefaultFigurePosition', position);
hold on

%Plot temporal profiles for each simulation replicate
plot(1:LengthT,squeeze(InfAsPropnAgeGrpAllRunsArray(:,:,18)),'Color',[0.5,0.5,0.5],'LineWidth',0.5);

if SimnRunType == 1  %Historical seasons only
    TotalYrs = varargin{1};

    %Specify x ticks
    xticks(0:365:TotalYrs*365)

    %Construct x tick labels
    XVals = cell(TotalYrs+1,1);
    XVals{1} = '0';
    for ii = 1:TotalYrs
        XVals{ii+1} = num2str(ii);
        xticklabels(XVals)
    end

    %Add lines denoting mid-season mark
    MidSeasonMark = (0:366:(TotalYrs-1)*366) + (366/2);
    for ii = 1:TotalYrs
        plot([MidSeasonMark(ii) MidSeasonMark(ii)],[0 1],'--','Color',[0.8,0,0],'LineWidth',1)
    end

    %Specify x tick labels
    %xticks([0 1*365 2*365 3*365 4*365 5*365 6*365 7*365 8*365])
    xticks(0:365:TotalYrs*365)
    xticklabels({'0','1','2','3','4','5','6','7','8','9'})

    %Add lines denoting mid-season mark
    MidSeasonMark = (0:366:(TotalYrs-1)*366) + (366/2);
    for ii = 1:9
        plot([MidSeasonMark(ii) MidSeasonMark(ii)],[0 1],'--','Color',[0.8,0,0],'LineWidth',1)
    end

elseif SimnRunType == 2 %Forward simulation, includes additional seasons
    AdditionalSeasonNum = varargin{1};

    TotalYrs = 8 + AdditionalSeasonNum;

    %Specify x ticks
    xticks(0:365:TotalYrs*365)

    %Construct x tick labels
    XVals = cell(TotalYrs+1,1);
    XVals{1} = '0';
    for ii = 1:TotalYrs
        XVals{ii+1} = num2str(ii);
        xticklabels(XVals)
    end

    %Add lines denoting mid-season mark
    MidSeasonMark = (0:366:(TotalYrs-1)*366) + (366/2);
    for ii = 1:TotalYrs
        plot([MidSeasonMark(ii) MidSeasonMark(ii)],[0 1],'--','Color',[0.8,0,0],'LineWidth',1)
    end

else
    error('Incorrect SimnRunType value (must be 1 or 2)');
end

%Specify axes labels
xlabel('Time (years)')
ylabel('Proportion of 17-18yrs infected')

%Specify ylim
ylim([0 max(InfAsPropnAgeGrpAllRunsArray(:))+0.01])

%Specify plot properties
set(gca,'FontSize',16)
set(gca,'LineWidth',1)

box on


%Save to file
export_fig(SaveFileName,'-pdf','-transparent','-r1200')
