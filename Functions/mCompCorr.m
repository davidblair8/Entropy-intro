function [FDR, Bonferroni, Sidak] = mCompCorr(N, p, pval)
% MCOMPCORR implements mutliple-comparison correction with FDR, Bonferroni,
% and Dunn-Sidak thresholds.
%   INPUTS:
%       N:    number of hypothesis tests to correct
%       p:    vector of p-values from each hypothesis test.  Must be a column
%             vector of length N.
%       pval: desired p-value threshold for correction.  May be scalar
%             (single threshold) or vector (multiple thresholds).  If
%             vector, has length S.
%   OUTPUTS:
%   	FDR:        NxS logical array of FDR results.  A true value
%                   indicates that the test survives correction.
%       Sidak:      NxS logical array of Sidak threshold results.  A true
%                   value indicates that the test survives correction.
%       Bonferroni:	NxS logical array of Bonferroni threshold results.  A 
%                 	true value indicates that the test survives correction.

% Set default values
if isempty(N)
	N = size(p,1);
end

% Run FDR correction once per desired threshold
FDR = zeros(N, length(pval));
for s = 1:length(pval)
	FDR(sort(FDR_benjHoch(p, pval(s))), s) = 1;
end
FDR = logical(squeeze(FDR));

% Run Bonferroni multiple comparison correction
Bonferroni = zeros(N, length(pval));
for s = 1:length(pval)
	Bonferroni(:,s) = (p < (pval(s)/N));
end
Bonferroni = logical(squeeze(Bonferroni));

% Run Dunn-Sidak multiple comparison correction
Sidak = zeros(N, length(pval));
for s = 1:length(pval)
	alpha = 1-(1-pval(s))^(1/N);
	Sidak(:,s) = (p < alpha);
end
Sidak = logical(squeeze(Sidak));

end

