function [h, p] = sigtest(e)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% preallocate arrays
h = nan(1,3);
p = nan(1,3);

% Student's t-test
isnorm = ~(jbtest(squeeze(e(:,1))) || jbtest(squeeze(e(:,2))));
if isnorm
    [h(1), p(1)] = ttest2(squeeze(e(:,1)), squeeze(e(:,2)));
else
    h(1) = 0;
    p(1) = NaN;
end

% Kolmogorov-Smirnov test
[h(2), p(2)] = kstest2(e(:,1), e(:,2));

% permutation test
p(3) = permutationTest(e(:,1), e(:,2), 10000, 'exact',0, 'sidedness','both');
h(3) = (p(3) < 0.05);

% convert to table
h = array2table(logical(h), 'VariableNames',["t" "ks" "perm"]);
p = array2table(p, 'VariableNames',["t" "ks" "perm"]);

end