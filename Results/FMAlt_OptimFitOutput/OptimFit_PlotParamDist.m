%Purpose:
%Plot parameter distributions for particles obtained from gradient
%free optimisation scheme
%--------------------------------------------------------------------------

clear variables

%--------------------------------------------------------------------------
%Add required directories to path
%--------------------------------------------------------------------------
addpath('../../../ModelAgnosticScripts')

%--------------------------------------------------------------------------
%%% Declare if want plots to indicate the median
%--------------------------------------------------------------------------
MedianFlag = 0;

%--------------------------------------------------------------------------
% Specify files to be accessed
%--------------------------------------------------------------------------
SingleAscertainParamsFileName = '../FMAlt_OptimFitOutputFiles/OptimiserParamsTrace#14_Emp_TransContactMatrix_JobArrayCombined.txt';
DualAscertainParamsFileName = '../FMAlt_OptimFitOutputFiles/OptimiserParamsTrace#15_Emp_TransContactMatrix_JobArrayCombined.txt';

%--------------------------------------------------------------------------
% Load data
%--------------------------------------------------------------------------
SingleAscertainModelParamsTemp = dlmread(SingleAscertainParamsFileName);
%DualAscertainModelParams = dlmread(DualAscertainParamsFileName);

%--------------------------------------------------------------------------
% Process data
%--------------------------------------------------------------------------

%Swap B/Yamagata and B/Victoria transmissibility param value columns.
SingleAscertainModelParams = SingleAscertainModelParamsTemp;
SingleAscertainModelParams(:,3) = SingleAscertainModelParamsTemp(:,4); %B/Victoria, was column 4, move to column 3
SingleAscertainModelParams(:,4) = SingleAscertainModelParamsTemp(:,3); %B/Yamagata, was column 3, move to column 4

% DualAscertainModelParams = DualAscertainModelParamsTemp;
% DualAscertainModelParams(:,3) = DualAscertainModelParamsTemp(:,4); %B/Victoria, was column 4, move to column 3
% DualAscertainModelParams(:,4) = DualAscertainModelParamsTemp(:,3); %B/Yamagata, was column 3, move to column 4


%--------------------------------------------------------------------------
% Use error files to split parameter sets based on closeness to empirical data
%--------------------------------------------------------------------------

%Load file
SingleAscertainErrorFileName = '../FMAlt_OptimFitOutputFiles/OptimiserErrors#14_Emp_TransContactMatrix_JobArrayCombined.txt';

% Load data
SingleAscertainModelErrors = dlmread(SingleAscertainErrorFileName);
SingleAscertainModelErrors = SingleAscertainModelErrors;

HighErrorThreshold = 30000;
SingleAscertainModel_LowErrorIdx = SingleAscertainModelErrors<=HighErrorThreshold;
SingleAscertainModel_HighErrorIdx = SingleAscertainModelErrors>HighErrorThreshold;

%Split parameters into distinct variables
SingleAscertainModelParams_LowError = SingleAscertainModelParams(SingleAscertainModel_LowErrorIdx,:);
SingleAscertainModelParams_HighError = SingleAscertainModelParams(SingleAscertainModel_HighErrorIdx,:);

%--------------------------------------------------------------------------
% Plot label set up
%--------------------------------------------------------------------------

%Declare parameter labels
SingleAscertainModelLabels =...
    {'$R_{0_{A(H1N1)pdm09}}$','$R_{0_{A(H3N2)}}$','$R_{0_{B/Victoria}}$','$R_{0_{B/Yamagata}}$',...
            'a (Nat. Inf. Mod. Sus.)','b (Inf. B cross-reactivity)','$\xi$',...
            '$\sigma_{0-17}$','$\sigma_{18-64}$','$\sigma_{65-84}$','$\sigma_{85+}$',...
            '$\epsilon_{2012/13}$','$\epsilon_{2013/14}$','$\epsilon_{2014/15}$',...
            '$\epsilon_{2015/16}$','$\epsilon_{2016/17}$','$\epsilon_{2017/18}$',...
            '$\tilde{\epsilon}_{0yrs}$','$\tilde{\epsilon}_{2yrs}$','$\tilde{\epsilon}_{18yrs}$','$\tilde{\epsilon}_{65yrs}$','$\tilde{\epsilon}_{85yrs}$'};

DualAscertainModelLabels =...
    {'R_{0_{A(H1N1)pdm09}}','R_{0_{A(H3N2)}}','R_{0_{B/Victoria}}','R_{0_{B/Yamagata}}',...
            'a (Nat. Inf. Mod. Sus.)','b (Inf. B cross-reactivity)','\xi',...
            '\sigma_{0-17}','\sigma_{18-64}','\sigma_{65+}',...
            '\epsilon_{2012/13}^{A}','\epsilon_{2013/14}^{A}','\epsilon_{2014/15}^{A}',...
            '\epsilon_{2015/16}^{A}','\epsilon_{2016/17}^{A}','\epsilon_{2017/18}^{A}',...
            '\epsilon_{0yrs}^{A}','\epsilon_{2yrs}^{A}','\epsilon_{18yrs}^{A}','\epsilon_{65yrs}^{A}','\epsilon_{85yrs}^{A}',...
             '\epsilon_{2012/13}^{B}','\epsilon_{2013/14}^{B}','\epsilon_{2014/15}^{B}',...
            '\epsilon_{2015/16}^{B}','\epsilon_{2016/17}^{B}','\epsilon_{2017/18}^{B}',...
            '\epsilon_{0yrs}^{B}','\epsilon_{2yrs}^{B}','\epsilon_{18yrs}^{B}','\epsilon_{65yrs}^B}','\epsilon_{85yrs}^{B}'};

x_label = SingleAscertainModelLabels;

%--------------------------------------------------------------------------
%%% SET UP PRIOR DISTRIBUTION PLOTTING VARIABLES
%--------------------------------------------------------------------------
PriorDataVals = 1;
PriorDataFlag = 1; %Flag variable: 0 (1) -> Do not (Do) include prior distribution on plots
%PriorPdfSampleVals = [];
%PriorPdfVals = [];

ParamNum = numel(x_label);

%Sepcify number of prior pdf points to plot
PriorPdfSamplePts = 100; %Number of points to compute value of prior pdf at
PriorPdfSampleVals = zeros(ParamNum,PriorPdfSamplePts); %Storage array for prior val sample points
PriorPdfVals = zeros(ParamNum,PriorPdfSamplePts); %Storage array for prior values

if PriorDataVals == 1

    %Specify prior distribution plot range
    R0_lb = 1.5*ones(4,1); R0_ub = 2.5*ones(4,1);
    ExpHist_lb = [0;0;0]; ExpHist_ub = [1;1;1];

    Suscep_lb = 0*ones(4,1); Suscep_ub = 1*ones(4,1);

    AscertainProb_lb = 0*ones(6,1); AscertainProb_ub = 0.1*ones(6,1);
    AscertainProbMod_lb = 0*ones(5,1); AscertainProbMod_ub = 1*ones(5,1);

    AllParam_lb = [R0_lb;ExpHist_lb;Suscep_lb;AscertainProb_lb;AscertainProbMod_lb];
    AllParam_ub = [R0_ub;ExpHist_ub;Suscep_ub;AscertainProb_ub;AscertainProbMod_ub];

    %Evaluate prior pdf at requested number of points
    for ii = 1:ParamNum
        PriorPdfSampleVals(ii,:) = linspace(AllParam_lb(ii),AllParam_ub(ii),PriorPdfSamplePts);
        PriorPdfVals(ii,:) = unifpdf(PriorPdfSampleVals(ii,:),AllParam_lb(ii),AllParam_ub(ii));
    end
end

%Aggregate prior distribution information into a single cell array
PriorDataVar = {PriorDataFlag,PriorPdfSampleVals,PriorPdfVals};

%--------------------------------------------------------------------------
%%% SET UP SYNTH DATA CELL
%--------------------------------------------------------------------------
SynthDataFlag = 0;
% ParamTrueVal = [];

ParamTrueVal_R0 = [1.65,1.45,1.6,1.5]; %Note, have flipped B/Vic and B/Yam values!
ParamTrueVal_ExpHist = [0.3,0.5,0.3];
ParamTrueVal_AscertainProb = [0.005,0.01,0.03,0.01,0.0075,0.01];

if PriorDataVals == 1
    ParamTrueVal_Suscep = ones(1,6);
    ParamTrueVal = [ParamTrueVal_R0 ParamTrueVal_ExpHist ParamTrueVal_Suscep ParamTrueVal_AscertainProb];
elseif PriorDataVals == 2
    ParamTrueVal_Suscep = ones(1,6);
    ParamTrueVal_AscertainScale = [1,1,1];
    ParamTrueVal = [ParamTrueVal_R0 ParamTrueVal_ExpHist ParamTrueVal_Suscep ParamTrueVal_AscertainProb,ParamTrueVal_AscertainScale];
elseif PriorDataVals == 3
    ParamTrueVal_Suscep = ones(1,7);
    ParamTrueVal_AscertainScale = [1,1,1];
    ParamTrueVal = [ParamTrueVal_R0 ParamTrueVal_ExpHist ParamTrueVal_Suscep ParamTrueVal_AscertainProb,ParamTrueVal_AscertainScale];
elseif PriorDataVals == 4
    ParamTrueVal_Suscep = ones(1,3);
    ParamTrueVal_AscertainScale = [1,1,1];
    ParamTrueVal = [ParamTrueVal_R0 ParamTrueVal_ExpHist ParamTrueVal_Suscep ParamTrueVal_AscertainProb,ParamTrueVal_AscertainScale];
end

SynthDataVar = {SynthDataFlag,ParamTrueVal};

%--------------------------------------------------------------------------
%%% SPECIFY PLOTTING VARIABLES
%--------------------------------------------------------------------------
SubPlotRowNum = 5;
SubPlotColNum = 6;
FigSizeWidthScale = 4.5;
FigSizeHeightScale = 3.5;

%Set plotting panels to be used
PlotPanelsUsed = [1:4 7:9 13:16 19:24 25:29];

%--------------------------------------------------------------------------
%%% DECLARE HISTOGRAM NORMALISATION OPTION
%--------------------------------------------------------------------------
NormOption = 1; %1 - pdf; otherwise, plots relative probability

%--------------------------------------------------------------------------
%%% STATE SAVE FILENAME
%--------------------------------------------------------------------------
SaveFileName_SingleAscertainModel_Histogram = 'FMAlt_OptimFitParamHistograms_SingleAscertainModel_TransContactMatrix_March2020';
SaveFileName_SingleAscertainModel_ScaledDensity = 'FMAlt_OptimFitParamScaledDensity_SingleAscertainModel_TransContactMatrix_March2020';
SaveFileName_SingleAscertainModel_Overlayed = 'FMAlt_OptimFitParamHistogramsNoPriorOverlayed_SingleAscertainModel_TransContactMatrix_March2020';

SaveFileName_DualAscertainModel_Histogram = 'FMAlt_OptimFitParamHistograms_DualAscertainModel';
SaveFileName_DualAscertainModel_ScaledDensity = 'FMAlt_OptimFitParamScaledDensity_DualAscertainModel';
SaveFileName_DualAscertainModel_Overlayed = 'FMAlt_OptimFitParamHistogramsNoPriorOverlayed_DualAscertainModel';

% %Plot histograms with priors
% PostDistHistograms(SingleAscertainModelParams_LowError,x_label,SynthDataVar,PriorDataVar,...
%                                 SubPlotRowNum,SubPlotColNum,...
%                                 FigSizeWidthScale,FigSizeHeightScale,...
%                                 PlotPanelsUsed,...
%                                 MedianFlag,NormOption,SaveFileName_SingleAscertainModel_ScaledDensity)


% %Plot histograms omitting priors
NormOption = 1; %1 - pdf; otherwise, plots relative probability
PriorDataFlag = 0;
PostDistHistograms(SingleAscertainModelParams_LowError,x_label,SynthDataVar,PriorDataVar,...
                                SubPlotRowNum,SubPlotColNum,...
                                FigSizeWidthScale,FigSizeHeightScale,...
                                PlotPanelsUsed,...
                                MedianFlag,NormOption,SaveFileName_SingleAscertainModel_Histogram)


% OverlayedPostDistHistograms(SingleAscertainModelParams(1:72,:),...
%                                 SingleAscertainModelParams_LowError,SingleAscertainModelParams_HighError,...
%                                 x_label,SynthDataVar,PriorDataVar,...
%                                 SubPlotRowNum,SubPlotColNum,...
%                                 FigSizeWidthScale,FigSizeHeightScale,...
%                                 MedianFlag,NormOption,SaveFileName_SingleAscertainModel_Overlayed)
