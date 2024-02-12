function [prop] = propLeading(D)
% PROPLEADING calculates the amount of variance captured by the leading
% eigenvector(s) of a subject-time timeseries array.
%	D is a cell containing the eigenvalues of 

prop = nan(1, size(D,2));
for t = 1:size(D,2)
    prop(t) = D(1,t)/sum(D(:,t));
end

end

