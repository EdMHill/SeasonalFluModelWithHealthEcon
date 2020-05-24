%Purpose:
%Construct arrays storing weekly Influenza positivity values
%monitored through the RCGP sentinel swabbing scheme in England
%Save to text/MAT files

%Method:
%Values obtained from Public Health England weekly & annual national influenza reports
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%RCGP sentinel swabbing
%Data from these schemes are used to calculate the proportion of ILI cases consulting
%in primary care who are swabbed and who are positive for influenza each week.

%Subset of general practices in the RCGP Weekly Returns Service
%submits respiratory samples for virological testing from patients presenting in primary
%care with an ILI.
%--------------------------------------------------------------------------

%% Arrays containing data

%Initialise arrays
% Row per season (Row 1 - 2009/2010, Row 2 - 2010/2011, ... , Final row - 2016/2017)
% Column per week (Col 1 - week 40, Col 2 - week 41, ... , Final Col - week 20)

RCGP_OverallPropnFluPositive = cell(1,9);


%Season 2009/2010
%Data taken from 2010/11 PHE report, has curve for entire 2009/10 year
RCGP_OverallPropnFluPositive{1} =...
    [8.0,13.7,11.7,11.4,... %Weeks 36-39
    17.2,23.3,35.5,41.4,40.5,32.5,... %Weeks 40-45
    33.0,27.5,22.9,22.4,22.9,17.8,21.7,27.0,... %Weeks 46-53
    9.4,9.6,10.3,8.7,4.3,3.2,1.8,... %Weeks 1-7
    3.7,6.4,9.2,6.9,2.7,0.0,9.8,... %Weeks 8-14
    0.0,2.9,0.0,0.0,0.0,0.0,... %Weeks 15-20
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]/100;  %Weeks 21-35


%Season 2010/2011
%Week 36-39 set to zero
%Weeks 14-20 subjectively chosen values
RCGP_OverallPropnFluPositive{2} =...
    [0,0,0,0,... %Weeks 36-39
    1.4,2.2,6.2,14.6,9.5,6.2,... %Weeks 40-45
    5.4,17.8,29.7,50.0,65.3,65.0,61.0,... %Weeks 46-52
    40.5,31.1,22.4,22.7,20.0,14.1,15.4,... %Weeks 1-7
    3.2,10.8,6.5,10.5,4.8,7.8,5.0,... %Weeks 8-14
    5.0,5.0,5.0,5.0,5.0,5.0,... %Weeks 15-20
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]/100;  %Weeks 21-35


%Season 2011/2012
%Week 36-39 set to zero
%Weeks 16-20 subjectively chosen values
RCGP_OverallPropnFluPositive{3} =...
    [0,0,0,0,... %Weeks 36-39
    0.0,1.6,1.6,1.6,0.0,0.0,... %Weeks 40-45
    0.0,0.0,1.9,2.2,0.8,3.2,10.5,... %Weeks 46-52
    6.0,4.9,9.0,13.5,16.2,32.4,34.9,... %Weeks 1-7
    34.7,31.5,30.0,36.5,37.2,15.5,15.0,... %Weeks 8-14
    19.5,15.0,10.0,6.0,4.0,2.0,... %Weeks 15-20
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]/100;  %Weeks 21-35


%Season 2012/2013 (NOTE: USES DATAMART DATA!)
%Week 36-39 set to zero
RCGP_OverallPropnFluPositive{4} =...
    [0,0,0,0,... %Weeks 36-39
    0.3,0.7,0.9,1.0,2.2,1.0,... %Weeks 40-45
    2.9,3.4,4.3,5.9,13.8,21.4,25.7,... %Weeks 46-52
    23.8,19.7,15.5,19.8,25.5,22.9,21.9,... %Weeks 1-7
    26.6,20.3,19.0,20.7,19.3,15.9,13.6,... %Weeks 8-14
    15.5,8.4,6.7,3.4,2.9,0.7... %Weeks 15-20
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]/100;  %Weeks 21-35


%Season 2013/2014
%Week 36-39 set to zero
%Weeks 16 subjectively chosen values
RCGP_OverallPropnFluPositive{5} =...
    [0,0,0,0,... %Weeks 36-39
    0,0,0,2.8,2.2,1.9,... %Weeks 40-45
    0,0,1.9,3.2,8.0,8.9,2.3,... %Weeks 46-52
    5.2,13.9,24.6,18.1,23.0,35.8,39.6,... %Weeks 1-7
    31.5,27.8,44.1,43.0,18.1,15.8,14.4,... %Weeks 8-14
    4.1,5.0,0.0,0.0,0.0,0.0... %Weeks 15-20
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]/100;  %Weeks 21-35


%Season 2014/2015
%Week 36-39 set to zero
%Weeks 16-18, 20 from weekly reports
%Weeks 19 subjectively chosen values
RCGP_OverallPropnFluPositive{6} =...
    [0,0,0,0,... %Weeks 36-39
    0,3.7,3.0,0,4.5,0,... %Weeks 40-45
    5.2,8.9,7.4,9.7,30.5,40.8,44.9,... %Weeks 46-52
    40.7,37,26.7,22.3,36.3,23.0,37.0,... %Weeks 1-7
    26.7,37.0,28.2,41.5,39.3,34.1,... %Weeks 8-13
    28.1,37.8,30.0,22.2,14.3,10.0,0.0... %Weeks 14-20
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]/100;  %Weeks 21-35

% Readings in centimetres. 1.35cm to 10%
%     [0,0.5,0.4,0,0.6,0,... %Weeks 40-45
%     0.7,1.2,1,1.3,4.1,5.5,6,... %Weeks 46-52
%     5.5,5,3.6,3,4.9,3.1,5,... %Weeks 1-7
%     3.6,5,3.8,5.6,5.3,4.6,3.8,... %Weeks 8-14
%     5.1,16,17,18,19,20]; %Weeks 15-20

%Season 2015/2016
%Week 36-39 set to zero
%Weeks 18-20 partially informed by FluNet
%Weeks 21 onward set to zero (weekly reports)
RCGP_OverallPropnFluPositive{7} =...
    [0,0,0,0,... %Weeks 36-39
    0.0,0.0,0.0,1.5,2.2,3.0,... %Weeks 40-45
    0.0,8.9,3.0,5.2,15.6,13.3,17.0,21.5,... %Weeks 46-53
    14.1,21.5,24.5,43.0,37.0,37.8,38.5,... %Weeks 1-7
    42.0,37.5,55.0,61.7,58.5,44.4,35.6,... %Weeks 8-14
    23.0,20.7,19.3,22.4,21.65,15.9... %Weeks 15-20
    0.0,0.0,0.0,0.0,0,0,0,0,0,0,0,0,0,0,0]/100;  %Weeks 21-35

% Readings in centimetres. 2.7cm to 20%
%     [0.0,0.0,0.0,0.2,0.3,0.4,... %Weeks 40-45
%     0.0,1.2,0.4,0.7,2.1,1.8,2.3,2.9,... %Weeks 46-53
%     1.9,2.9,3.3,5.8,5.0,5.1,5.2,... %Weeks 1-7
%     5.6,5.0,7.3,8.2,7.9,6.0,4.8,... %Weeks 8-14
%     3.1,2.8,2.6,18,19,20]; %Weeks 15-20

%Season 2016/2017
%Week 36-39 set to zero
%Weeks 14-20 from weekly reports (all zero, though few samples)
RCGP_OverallPropnFluPositive{8} =...
    [0,0,0,0,... %Weeks 36-39
    0,0,0,0.9,0,4.2,... %Weeks 40-45
    4.5,5.5,4.8,18.2,34.2,30.0,41.8,... %Weeks 46-52
    39.0,40.6,38.5,30.9,40.3,22.1,35.5,... %Weeks 1-7
    26.1,24.2,11.8,4.2,7.0,0.0,0.0,... %Weeks 8-14
    0.0,0.0,0.0,0.0,0.0,0.0... %Weeks 15-20
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]/100;  %Weeks 21-35

%Season 2017/2018
%Week 36-39 set to zero
%Week 16 based on weekly reports
%Weeks 17 set to WHO flunet values
%Week 18 onward all zeros, based on weekly reports
RCGP_OverallPropnFluPositive{9} =...
    [0,0,0,0,... %Weeks 36-39
    0,0,1.2,1.2,4.2,6.1,... %Weeks 40-45
    3.9,11.2,9.1,13.3,19.7,45.75,45.75,... %Weeks 46-52
    41.2,39.4,46.65,52.1,48.5,45.45,46.65,... %Weeks 1-7
    53.0,50.9,53.65,55.75,56.0,38.2,39.4,33.95,... %Weeks 8-15
    28.5,20.0,14.1,0.0,0.0,0.0... %Weeks 16-21
    0.0,0.0,0,0,0,0,0,0,0,0,0,0,0,0]/100;  %Weeks 21-35

%--------------------------------------------------------------------------

%% Save to MAT file
save('RCGPSamplePositivityDataFile_2009to2018.mat','RCGP_OverallPropnFluPositive')

%% Save to text file
fileID = fopen('RCGPSamplePositivityTxtFile_2009to2018.txt','w');

[nrows,ncols] = size(RCGP_OverallPropnFluPositive);
for ii = 1:ncols
    fprintf(fileID, '%f ',RCGP_OverallPropnFluPositive{ii});
    fprintf(fileID, '\n',RCGP_OverallPropnFluPositive{ii});
end

fclose(fileID);
