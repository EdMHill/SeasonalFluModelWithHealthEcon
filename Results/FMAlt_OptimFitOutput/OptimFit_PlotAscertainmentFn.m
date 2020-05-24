% Purpose:
% Plot linear ascertainment function
% Include median and 95% credible interval
% 
% Author: Ed Hill
% Date: Ocotber 2019
%--------------------------------------------------------------------------

clear variables

%--------------------------------------------------------------------------
% Import data
%--------------------------------------------------------------------------
SingleAscertainParamsFileName = '../FMAlt_OptimFitOutputFiles/OptimiserParamsTrace#14_Emp_TransContactMatrix_JobArrayCombined.txt';
SingleAscertainModelParamsTemp = dlmread(SingleAscertainParamsFileName);

%--------------------------------------------------------------------------
% Extract relevant portions of data
%--------------------------------------------------------------------------

%Reference ascertainment values (ages 100+)
SingleAscertainModel_AscertainmentRefParams = SingleAscertainModelParamsTemp(:,12:17);

%Relative ascertainment values, scaling factors at knot ages
SingleAscertainModel_AscertainmentScaleParams = SingleAscertainModelParamsTemp(:,18:22);

%--------------------------------------------------------------------------
% Compute quantiles for 95% credible interval
%--------------------------------------------------------------------------
AscertainmentRefParams_PrctVals = prctile(SingleAscertainModel_AscertainmentRefParams,[2.5,50,97.5]);
AscertainmentScaleParams_PrctVals = prctile(SingleAscertainModel_AscertainmentScaleParams,[2.5,50,97.5]);

%--------------------------------------------------------------------------
% Median & 95% credible interval profiles
%--------------------------------------------------------------------------
SeasonNum = size(AscertainmentRefParams_PrctVals,2);
KnotVals = size(AscertainmentScaleParams_PrctVals,2) + 1;

%Initialise storage array for profiles
AscertainmentProfileBySeason_Median = zeros(SeasonNum,KnotVals);
AscertainmentProfileBySeason_LB = zeros(SeasonNum,KnotVals);
AscertainmentProfileBySeason_UB = zeros(SeasonNum,KnotVals);

for SeasonIdx = 1:SeasonNum
    
    %Assign value to 100+ age group
    AscertainmentProfileBySeason_Median(SeasonIdx,end) = AscertainmentRefParams_PrctVals(2,SeasonIdx);
    AscertainmentProfileBySeason_LB(SeasonIdx,end) = AscertainmentRefParams_PrctVals(1,SeasonIdx);
    AscertainmentProfileBySeason_UB(SeasonIdx,end) = AscertainmentRefParams_PrctVals(3,SeasonIdx);

    %Scale values at knot points
    AscertainmentProfileBySeason_Median(SeasonIdx,1:(end-1)) = AscertainmentRefParams_PrctVals(2,SeasonIdx).*AscertainmentScaleParams_PrctVals(2,:);
    AscertainmentProfileBySeason_LB(SeasonIdx,1:(end-1)) = AscertainmentRefParams_PrctVals(1,SeasonIdx).*AscertainmentScaleParams_PrctVals(1,:);
    AscertainmentProfileBySeason_UB(SeasonIdx,1:(end-1)) = AscertainmentRefParams_PrctVals(3,SeasonIdx).*AscertainmentScaleParams_PrctVals(3,:);
end

%%
%--------------------------------------------------------------------------
% Construct plot (panel per season)
%--------------------------------------------------------------------------

%Set up x-axis plotting points
xVal = [0 2 18 65 85 100];

%Get limits for y-axis
MaxYval = max(max(AscertainmentProfileBySeason_UB)) + 0.005;

%Define figure title vector
TitleVec = {'2012/13','2013/14','2014/15','2015/16','2016/17','2017/18'};

%Plot initialisation 
fig = figure(); 
clf;
set(fig,'Color', [1 1 1])
position = [10, 10, 3.5*550, 2.2*450];
set(0, 'DefaultFigurePosition', position);

%Iterate through each season
for SeasonIdx = 1:SeasonNum
    subplot(2,3,SeasonIdx)
    hold on
    
    
    %Fill area between the two credible interval boundary curves
    p2 = patch([xVal fliplr(xVal)], [AscertainmentProfileBySeason_LB(SeasonIdx,:) fliplr(AscertainmentProfileBySeason_UB(SeasonIdx,:))],...
                                [0 0 1],... %patch colour
                                'FaceAlpha',0.3,...
                                'DisplayName','95% credible region');
                            
    %Plot median line
    p1 = plot(xVal,AscertainmentProfileBySeason_Median(SeasonIdx,:),'x--',...
        'Color',[0 0 1],'LineWidth',1.5,'MarkerSize',10,...
        'DisplayName','Median');
    
    %Plot 95% credible interval
    plot(xVal,AscertainmentProfileBySeason_LB(SeasonIdx,:),'--','Color',[0.5 0.5 0.5],'LineWidth',1.5)
    plot(xVal,AscertainmentProfileBySeason_UB(SeasonIdx,:),'--','Color',[0.5 0.5 0.5],'LineWidth',1.5)
    
    %Set axes limits
    ylim([0 MaxYval])
    
    %Set axes tick labels
    xticks([-1 2 18 65 85 100])
    xticklabels({'0','2','18','65','85','100+'})
    
    %Add title
    title(TitleVec(SeasonIdx))
    
    %Set axes labels
    if SeasonIdx > 3
        xlabel('Age (years)')
    end
    
    if mod(SeasonIdx,3) == 1
        ylabel('Ascertainment probability')
    end
    
    %Add legend
    if SeasonIdx == 3
        legend([p1 p2])
    end
    
    %Specify general axis properties
    set(gca,'FontSize',20)
    set(gca,'LineWidth',1)
    box on
    
end

%Save file
FileName='FMAlt_OptimFitParam_LinearAscertainProfiles';
export_fig(FileName,'-pdf','-transparent','-r1200','-painters')
