%Purpose:
%Compute posterior distribution summary statistics using weighted samples
%gathered from an inference procedure.
%--------------------------------------------------------------------------

function OutputPrctileVals = PostDistStats(PosteriorParams,ParamWeights,ReqPrctiles)
%Inputs
%PosteriorParams - (2D Array) Row per simulation run. Column per parameter.
%ParamWeights - (Column vector)
%ReqPrctiles - (1D vector)

%Outputs
%OutputPrctileVals - (2D array) Row per requested prctile value. Column per parameter.

%--------------------------------------------------------------------------
%%% VARIABLE INITIALISATION
%--------------------------------------------------------------------------

%Number of parameters (columns of StatInput)
ParamNum = size(PosteriorParams,2);

%Number of prctile evaluations requested
PrctileNum = numel(ReqPrctiles);

%Convert percentage to proportion
ReqPropn = ReqPrctiles/100;

%Initialise output array
OutputPrctileVals = zeros(PrctileNum,ParamNum);

%--------------------------------------------------------------------------
%%% NORMALISE WEIGHTS
%--------------------------------------------------------------------------
NormParamWeights = ParamWeights/sum(ParamWeights);

%--------------------------------------------------------------------------
%%% WEIGHTED PERCENTILE COMPUTATION
%--------------------------------------------------------------------------

%Sort each column of StatInput into ascending order
[SortedStatInput,SortIdx] = sort(PosteriorParams,1);

%For each parameter, compute cumulative sum of weights
CumSumWeightsByParam = cumsum(NormParamWeights(SortIdx),1);

%Iterate over requested prctile values.
for ii = 1:PrctileNum

    if ReqPropn(ii) < 0.5
        %------------------------------------------------------------------
        %%% PRCTILES BELOW 50 (use infimum. Greatest lower bound)
        %------------------------------------------------------------------

        for jj = 1:ParamNum
              SelectedParamIdx = find(CumSumWeightsByParam(:,jj) <= ReqPropn(ii),1,'last');
              OutputPrctileVals(ii,jj) = SortedStatInput(SelectedParamIdx,jj);
        end

    elseif ReqPropn(ii) == 0.5

        %------------------------------------------------------------------
        %%% MEDIAN
        %------------------------------------------------------------------
        for jj = 1:ParamNum
            OutputPrctileVals(ii,jj) = weightedMedian(PosteriorParams(:,jj),ParamWeights);
        end
    else
        %------------------------------------------------------------------
        %%% PRCTILES ABOVE 50 (use supremum. Least upper bound)
        %------------------------------------------------------------------

        for jj = 1:ParamNum
              SelectedParamIdx = find(CumSumWeightsByParam(:,jj) >= ReqPropn(ii),1,'first');
              OutputPrctileVals(ii,jj) = SortedStatInput(SelectedParamIdx,jj);
        end
    end
end
