function S=subsetSum(N,M)
%%SUBSETSUM computes all positive possibilities of 0,...,N-digits summing up to
% M using dynamical programming: S(N+1,M) is deduced by using optimal 
% substructure in M and M.
%   Input:
%       N (int): number of digits
%       M (int): sum of all digits
%   Output:
%       S (n x M matrix): the subset sums
%   Usage:
%       S=subsetSum(N,M)

% Initialize empty cell array for all subproblems.
subS=cell(M+1,N); % need one more row than M, because of zero digits

% N=1 equals 1 digit, meaning the only possiblity is 0,...,M itself
subS(:,1)=num2cell(0:M);

for iN=2:N
    % to determine all possiblities with iN digits, use all known possibilities
    % with iN-1 digits
    for iM=0:M
        % to determine all possibilities for iM+1, use all known
        % possibilities up to iM
        for m=0:iM
            % known possiblities up to iM-m with iN-1 digits
            h=subS{iM-m+1,iN-1};
            % new possibilities with iN digits by adding m as next digit to
            % subproblem h
            h(:,end+1)=m;
            subS{iM+1,iN}=[subS{iM+1,iN};h];
        end        
    end
end
S=subS{M+1,N};
end