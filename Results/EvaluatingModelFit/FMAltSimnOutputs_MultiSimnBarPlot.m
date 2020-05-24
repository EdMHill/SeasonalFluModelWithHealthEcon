% Purpose:
% Back-to-back bar plots for FMAlt simulation output
% Plot multiple replicates against data
% Per run, GP consultations for flu per strain (per 100,000 population)
%
% Panel per season. Within each season, horizontal stacked bar per age band
%--------------------------------------------------------------------------
clear variables

%% SET FLAG VARIABLES
EmpDataFlag = 1; %Binary indicator. 1 for empirical data. 0 for synthetic data
JuliaFlag = 1; %Binary indicator. 1 for Julia. 0 otherwise.

if EmpDataFlag ~= 0 && EmpDataFlag ~= 1
    error('EmpDataFlag must take value 0 or 1.')
end

if JuliaFlag ~= 0 && JuliaFlag ~= 1
    error('JuliaFlag must take value 0 or 1.')
end

%% SET simnID
simnID = 'TransContactArray_HistStrat';

%% LOAD DATA

if EmpDataFlag == 1
    %Empirical data
    Truth=load('../EmpiricalData_FMAlt/FMAlt_FluPosData_2009to2018.csv');
    TruthRetained = Truth(4:end,:);
else
    %Synthetic data
    Truth=load('../SynthData_FMAlt/FMAlt_Simn#4B_FluPosData.csv');
    TruthRetained = Truth(1:end,:);
end

%Put into 3D array
SeasonsToPlot = 6;
M = 101; %Number of age classes
NumOfStrains = 4;
TruthRetainedByStrain = reshape(TruthRetained,[SeasonsToPlot,M,NumOfStrains]);

%% LOAD MODEL SIMULATION OUTPUTS

%Simulated data
if JuliaFlag == 1
    load(strcat('FMAlt_SimnOutputFiles_Julia/FMAlt_SimnJuliaV1ModelRun_#',simnID,'.mat'))
else
    load(['FMAlt_SimnOutputFiles_Matlab/FMAlt_SimnMatlabModelRun_#',num2str(simnID),'.mat'])
end
FMAltSimn_CellOutput = SimnData(:,1); %From data file, get array outputs from each season

%Put into 4D array (seasons,age classes, strains, # of simns)
SimnNum = size(FMAltSimn_CellOutput,1);
FMAltSimn_GPFluConsultOutput = zeros(SeasonsToPlot,M,NumOfStrains,SimnNum);

for ii=1:SimnNum
    FMAltSimn_GPFluConsultOutput(:,:,:,ii) = FMAltSimn_CellOutput{ii};
end

%% AGGREGATE DATA & MODEL SIMULATION OUTPUT TO CONFORM TO USER-SPECIFIED AGE BANDS
AgeBandLowerBounds = [0 2 11 18 41 65 85];
AgeBandUpperBounds = [1 10 17 40 64 84 100];
AgeBandNum = length(AgeBandUpperBounds);

TruthRetainedAgg = zeros(SeasonsToPlot,AgeBandNum,NumOfStrains);
FMAltSimn_GPFluConsultOutputAgg = zeros(SeasonsToPlot,AgeBandNum,NumOfStrains,SimnNum);

for ii = 1:AgeBandNum
    StartIdx = AgeBandLowerBounds(ii) + 1;
    EndIdx = AgeBandUpperBounds(ii) + 1;

    %Mean across age groups into user-specified age bands
    TruthRetainedAgg(:,ii,:) = mean(TruthRetainedByStrain(:,StartIdx:EndIdx,:),2);
    FMAltSimn_GPFluConsultOutputAgg(:,ii,:,:) = mean(FMAltSimn_GPFluConsultOutput(:,StartIdx:EndIdx,:,:),2);
end

%% CONSTRUCT PLOT - Panel per season

%Define figure title vector
TitleVec = {'2012/13','2013/14','2014/15','2015/16','2016/17','2017/18'};
YAxisLabels = {'0-1','2-10','11-17','18-40','41-64','65-84','85+'};

%Plot initialisation
fig = figure();
clf;
set(fig,'Color', [1 1 1])
position = [10, 10, 3.5*550, 2.2*450];
set(0, 'DefaultFigurePosition', position);

YLimMax = AgeBandNum+1.5;

%Iterate through each season
for Year = 1:SeasonsToPlot
    subplot(2,3,Year)
    hold on
    plot([0 0],[0.8 YLimMax],'-k'); %Have central dividing line

    for AgeGrp = 1:AgeBandNum
        clear R;
        R(1:NumOfStrains,1:SimnNum) = FMAltSimn_GPFluConsultOutputAgg(Year,AgeGrp,:,:);
        AH1=R(1,:); AH3=R(2,:); BYam=R(3,:); BVic=R(4,:);
        AH1_Truth = TruthRetainedAgg(Year,AgeGrp,1); AH3_Truth = TruthRetainedAgg(Year,AgeGrp,2); %Empirical type A data
        BYam_Truth = TruthRetainedAgg(Year,AgeGrp,3); BVic_Truth = TruthRetainedAgg(Year,AgeGrp,4); %Empirical type B data
        y = 0.1 + 0.7*[0:SimnNum-1]/(SimnNum-1);
        o1=[1:SimnNum]; o2=[1:SimnNum];
        % if you include the next line the two sides of the plot are independent,
        % if you comment it out the figure is messier, but a more acurate
        % representation.
        %[x,o]=sort(BYam+BVic); o2=o(q);
%
%         %Strain order: A(H3N2), A(H1N1)pdm09, B/Yamagata, B/Victoria
%         h1=fill([0 -AH1_Truth-AH3_Truth -AH1_Truth-AH3_Truth 0],AgeGrp+[0.85 0.85 0.95 0.95],'r','DisplayName','A(H3N2)'); set(h1,'FaceColor',[1 0.6 0]);
%         h2=fill([0 -AH1_Truth -AH1_Truth 0],AgeGrp+[0.85 0.85 0.95 0.95],'r','DisplayName','A(H1N1)pdm09'); set(h2,'FaceColor',[1 0 0]);
%         h3=fill([0 BYam_Truth+BVic_Truth BYam_Truth+BVic_Truth 0],AgeGrp+[0.85 0.85 0.95 0.95],'b','DisplayName','B/Victoria'); set(h3,'FaceColor','c');
%         h4=fill([0 BYam_Truth BYam_Truth 0],AgeGrp+[0.85 0.85 0.95 0.95],'b','DisplayName','B/Yamagata'); set(h4,'FaceColor',[0 0 1]);
%
%         h=fill([0 -AH1(o1)-AH3(o1) 0],AgeGrp+[0.1 y 0.8],'r'); set(h,'FaceColor',[1 0.7 0]); set(h,'EdgeColor',[0.99 0.99 0.99],'LineStyle','none');
%         h=fill([0 -AH1(o1) 0],AgeGrp+[0.1 y 0.8],'r'); set(h,'FaceColor',[1 0 0],'LineStyle','none');
%         h=fill([0 BYam(o2)+BVic(o2) 0],AgeGrp+[0.1 y 0.8],'c'); set(h,'LineStyle','none');
%         h=fill([0 BYam(o2) 0],AgeGrp+[0.1 y 0.8],'b'); set(h,'LineStyle','none'); set(h,'FaceColor',[0.4 0.4 1])
%
        %Strain order: A(H1N1)pdm09, A(H3N2), B/Victoria, B/Yamagata
        h1=fill([0 -AH1_Truth-AH3_Truth -AH1_Truth-AH3_Truth 0],AgeGrp+[0.85 0.85 0.95 0.95],'r','DisplayName','A(H1N1)pdm09'); set(h1,'FaceColor',[1 0 0]);
        h2=fill([0 -AH3_Truth -AH3_Truth 0],AgeGrp+[0.85 0.85 0.95 0.95],'r','DisplayName','A(H3N2)'); set(h2,'FaceColor',[1 0.6 0]);
        h3=fill([0 BYam_Truth+BVic_Truth BYam_Truth+BVic_Truth 0],AgeGrp+[0.85 0.85 0.95 0.95],'b','DisplayName','B/Yamagata'); set(h3,'FaceColor',[0 0 1]);
        h4=fill([0 BVic_Truth BVic_Truth 0],AgeGrp+[0.85 0.85 0.95 0.95],'b','DisplayName','B/Victoria'); set(h4,'FaceColor','c');

        h=fill([0 -AH1(o1)-AH3(o1) 0],AgeGrp+[0.1 y 0.8],'r'); set(h,'FaceColor',[1 0 0]); set(h,'EdgeColor',[0.99 0.99 0.99],'LineStyle','none');
        h=fill([0 -AH3(o1) 0],AgeGrp+[0.1 y 0.8],'r'); set(h,'FaceColor',[1 0.7 0],'LineStyle','none');
        h=fill([0 BYam(o2)+BVic(o2) 0],AgeGrp+[0.1 y 0.8],'b'); set(h,'LineStyle','none'); set(h,'FaceColor',[0.4 0.4 1]);
        h=fill([0 BVic(o2) 0],AgeGrp+[0.1 y 0.8],'c'); set(h,'LineStyle','none');

    end
    set(gca,'YLim',[0.8 YLimMax]);


    %Set y-axis tick length to zero
    ax = gca;
    ax.YAxis.TickLength = [0 0];
    yticklabels(YAxisLabels);

    %Set y-axis limits and label, based on age groups
    if AgeBandNum == 7
        yticks([1.5 2.5 3.5 4.5 5.5 6.5 7.5])
    end

    %Set x-axis limits and tick labels
    xlim([-200 200])
    xticks([-200,-100,0,100,200])
    xticklabels({'200','100','0','100','200'})

%     %Add influenza type descriptive text x-position
%     TextLabelXPos1 = -750;
%     TextLabelXPos2 = 750;

%     %Add influenza type descriptive text
%     TextLabelYPos = YLimMax - 0.3;
%     txt1 = 'Type A';
%     text(TextLabelXPos1,TextLabelYPos,txt1,'FontWeight','bold','FontSize',15)
%
%     txt2 = 'Type B';
%     text(TextLabelXPos2,TextLabelYPos,txt2,'FontWeight','bold','FontSize',15)
%
%
    %Axes labels
    if Year == 5
        xlabel('GP consultations attributable to influenza (per 100,000)')
    end

    if mod(Year,3) == 1
        ylabel('Age group (yrs)')
    end

    %Add title
    title(TitleVec(Year))

    %Add legend
    if Year == 3
        leg1 = legend([h1;h2;h4;h3],'Location','northeast','FontSize',15,...
            'Position',[0.244848487000445 0.834175086703761 0.0987012965493387 0.0904040378753584]);
        %leg1=legend('A - H1N12009','A - H3','B - Yamagata','B - Victoria','B - All undet.');
    end

    %Specify general axis properties
    set(gca,'FontSize',18)
    set(gca,'LineWidth',1)
    %set(leg1,'FontSize',12)
    box on
end

% %%% Add labels to panel %%
% annotation(fig,'textbox',...
%     [0.25 0.03 0.50 0.05],...
%     'String',{'GP consultations attributable to influenza (per 100,000)'},...
%     'LineStyle','none','FontSize',18,'FontWeight','normal');
%%
% Output to file
if JuliaFlag == 1
    FileName=['MATLABfigs/FMAlt_BarPlots/FMAlt_MKPlotStyle_JuliaSimnNum', num2str(simnID)];
else
    FileName=['MATLABfigs/FMAlt_BarPlots/FMAlt_MKPlotStyle_MatlabSimnNum', num2str(simnID)];
end

export_fig(FileName,'-pdf','-transparent','-r1200')
