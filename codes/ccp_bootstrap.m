function [BootMean,BootSigma,Peaks410,Peaks660,Amps410,Amps660] = ccp_bootstrap(Matrix,Weights,YAxis,BootstrapNumber,TraceNumber)
%CCP_BOOTSTRAP Bootstrap resampling and stacking
%   CCP_BOOTSTRAP(Matrix,BootstrapNumber,TraceNumber) stacks the
%   vectors in Matrix columns a BootstrapNumber of times.  This function
%   uses "ignorenan" which ignores NaN values in the matrix you are
%   bootstrapping.
%
%   10/17/2008
%   We are no longer outputting the whole bootstrap matrix (BootMatrix).
%   Instead, we are doing all processes on the matrix inside this funcion.
%   This includes finding "peaks" in the data near depths of 410 and 660.
%   It requires another input of the y-axis (the depth axis).
%
%   10/21/2008
%   I added the "nrootmean" to this function.  Comment or uncomment as
%   needed.  Edit "nrootmean" function to change N value.
%   I also added the output of the amplitudes for each depth pick.  For
%   now, this will be temporarily stored in the "P410sDepthError" and
%   "P660sDepthError" fields. 
%
%   Author: Kevin C. Eagar
%   Date Created: 06/20/2007
%   Last Updated: 10/21/2008

Range410 = [380 440];
Range660 = [620 700];

% Rows
TraceLength = size(Matrix,1);
% Columns
TotalTraces = size(Matrix,2);

% Preallocate BootMatrix
%--------------------------------------------------------------------------
BootMatrix = zeros(size(Matrix,1),BootstrapNumber);

if TotalTraces == 1
    
    BootMean = colvector(Matrix);
    BootSigma = zeros(length(BootMean),1);
    
    for n = 1:BootstrapNumber
        BootMatrix(:,n) = BootMean;
    end
    
else

    for n = 1:BootstrapNumber

        %Pick out the m number of random traces (they can repeat)
        %----------------------------------------------------------------------
        %rand('state',sum(100*clock))
        % Init the random number generator
        rng('shuffle','twister')
        BootstrapIndices = ceil(TotalTraces.*rand(TraceNumber,1));

        % Take the mean (stack) of the new matrix and store the stack as an
        % array in another new matrix
        %----------------------------------------------------------------------
        Check1 = isnan(Matrix(:,BootstrapIndices));
        Check2 = 0;
        for m = 1:size(Matrix,1)
            if (any(Check1(m,:)) && ~all(Check1(m,:))); Check2 = 1; break; end
        end
        if Check2 == 0
            %BootMatrix(:,n) = mean(Matrix(:,BootstrapIndices),2);
            %BootMatrix(:,n) = nrootmean(Matrix(:,BootstrapIndices),2);
            BootWeights = Weights(:,BootstrapIndices) .* Matrix(:,BootstrapIndices);
            BootMatrix(:,n) = sum(BootWeights,2)./sum(Weights(:,BootstrapIndices),2);
        else
            %BootMatrix(:,n) = ignoreNaN(Matrix(:,BootstrapIndices),@mean,2);
            %BootMatrix(:,n) = ignoreNaN(Matrix(:,BootstrapIndices),@nrootmean,2);
            BootWeights = Weights(:,BootstrapIndices) .* Matrix(:,BootstrapIndices);
%             BootMatrix(:,n) = ignoreNaN(BootWeights,@mean,2)/sum(Weights(~isnan(Weights(BootstrapIndices))));
            BootMatrix(:,n) = ignoreNaN(BootWeights,@sum,2)./ignoreNaN(Weights,@sum,2);
        end
    end

    % Take the mean of the new matrix of stacks (a stack of stacks) and the
    % standard deviation of the stack
    %--------------------------------------------------------------------------
    Check1 = isnan(BootMatrix);
    Check2 = 0;
    for m = 1:size(BootMatrix,1)
        if (any(Check1(m,:)) && ~all(Check1(m,:))); Check2 = 1; break; end
    end
    if Check2 == 0
        BootMean = colvector(mean(BootMatrix,2));
        BootSigma = colvector(std(BootMatrix') * 1.96);
    else
        BootMean = colvector(ignoreNaN(BootMatrix,@mean,2));
        BootSigma = colvector(ignoreNaN(BootMatrix',@std) * 1.96);
    end
    
end

% Within the ranges around the 410 and 660 depths, we are finding the
% maximum peak, determining its index in the matrix, and relating that to
% the same index in the y-axis to get a depth pick.  Then, we will output
% the array of depth picks for each of the two discontinuities.
%--------------------------------------------------------------------------

% Get the row indices for the 410 and 660 ranges
% Indices410 = intersect(find(YAxis >= Range410(1)),find(YAxis <= Range410(2)));
% Indices660 = intersect(find(YAxis >= Range660(1)),find(YAxis <= Range660(2)));

% Preallocate the vectors to make sure they are row vectors
Amps410 = nan * ones(1,size(BootMatrix,2));
Amps660 = nan * ones(1,size(BootMatrix,2));
Peaks410 = nan * ones(1,size(BootMatrix,2));
Peaks660 = nan * ones(1,size(BootMatrix,2));

% % Get the index of the maximum value down each column (each trace), only
% % for the indices within our ranges.
% for n = 1:size(BootMatrix,2)
%     if all(isnan(BootMatrix(Indices410,n)))
%         PIR410(n) = nan;
%         Amps410(n) = nan;
%     else
%         [Amps410(n),PIR410(n)] = max(BootMatrix(Indices410,n),[],1);
%     end
%     if all(isnan(BootMatrix(Indices660,n)))
%         PIR660(n) = nan;
%         Amps660(n) = nan;
%     else
%         [Amps660(n),PIR660(n)] = max(BootMatrix(Indices660,n),[],1);
%     end
% end
% 
% % Get the true index of the maximum value in each column within our ranges.
% PI410 = PIR410 + (Indices410(1) - 1);
% PI660 = PIR660 + (Indices660(1) - 1);
% 
% % Get the depth value of the indices by finding the values in the YAxis
% % (depth axis)
% for n = 1:length(PI410)
%     if ~isnan(PI410(n))
%         Peaks410(n) = YAxis(PI410(n));
%     else
%         Peaks410(n) = nan;
%     end
% end
% for n = 1:length(PI660)
%     if ~isnan(PI660(n))
%         Peaks660(n) = YAxis(PI660(n));
%     else
%         Peaks660(n) = nan;
%     end
% end