function [h, p, tstat, FDR, Bonferroni, Sidak] = robustTests(A, B, N, varargin)
% ROBUSTTESTS is a wrapper function for a statistical test paired with
% multiple comparison correction.
%	

% Set default values
if isempty(N)
	N = size(A,1);
end
pval = 0.05;
testtype = 'kstest2';
permutations = 10000;
exact = 0;

% unpack varagin
for k = 1:2:length(varargin)
	switch varargin{k}
		case 'p'
			pval = varargin{k+1};
		case 'testtype'
			testtype = varargin{k+1};
		case 'permutation'
			permutations = varargin{k+1};
		case 'exact'
			exact = varargin{k+1};
	end
end

% Preallocate storage arrays
h = nan(N, 1);
p = nan(N, 1);
tstat = nan(N, 1);

% Run tests for each row
switch testtype
	case 'kstest2'
		for r = 1:N
			[h(r), p(r), tstat(r)] = kstest2(A(r,:), B(r,:), 'Alpha',pval(1));
		end
		clear r
	case 'ttest2'
        [h, p, ~, stats] = ttest2(A', B', 'Alpha',pval(1));
        tstat = stats.tstat;
	case 'permutation'
		for r = 1:N
			[p(r), ~, tstat(r)] = permutationTest(A(r,:), B(r,:), permutations, 'exact',exact, 'sidedness','both');
			h(r) = p(r) < pval(1);
		end
	case 'ranksum'
		for r = 1:N
			[h(r), p(r), stats] = ranksum(A(r,:), B(r,:), 'Alpha',pval(1));
			tstat(r) = stats.ranksum;
		end
	otherwise
		
end

% Run mutliple comparison corrections on p-values
[FDR, Bonferroni, Sidak] = mCompCorr(N, p, pval);

end