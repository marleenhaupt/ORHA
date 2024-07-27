function [SignificantVariables,crit_p,adjusted_pvalues] = fdr_corrected_perm_test_MH(data, n_perm, q_value,tail)

% IMPORTANT: 
% This function is usually used in the lab under the name fdr_permutation_cluster_1sample_alld.m
% As I was always confused by the term "cluster" appearing in the file name
% while this test does not perform any cluster-related calculation, I
% renamed it. Greta already had the same idea which is why it is identical
% to her file named fdr_corrected_perm_test.m

% The only thing that needs to be adjusted in the code below is the path to
% the fdr_bh folder!

% Author: Dimitrios Pantazis, December 2015
% Modified: Siying Xie, XX XX
% Modified: Agnessa Karapetian, June 2021
% Modified: Greta Häberle, April 2022
% Renamed: Marleen Haupt, February 2023

% WHAT IT DOES
% Performs one-sided (>0) or two-sided pemutation test + FDR correction on the 'data'.
% The data array can have any number of dimensions (1D, 2D, 3D, 4D etc) in
% the order [observations x variable1 x variable2 x ...]
% Data from each observation are randomly multiplied by +-1 to create
% permutation samples. These samples are converted to pvalues and then the
% fdr_bh function is applied to identify statistically significant points.

% INPUT:
% data: observations x variable1 x variable2 x variable3 x ... (supports all dimensions)
% nperm: number of permutations
% q_value: alpha level for FDR correction
% tail: string 'right' or 'both'. Default is 'both'.
% StatMapPermPV: (optional) permutation pvalue map (see output)

% OUTPUT:
% SignificantVariables:  times with significant values
% pvalues: p-values of the data (ground truth)
% crit_p: critical (statictically decisive) p-value in FDR
% adjusted_pvalues: FDR corrected p-values of the ground truth

% add the path of the fdr_bh folder in your folder structure
addpath('/scratch/haum92/aging/eeg/code/stats/fdr_bh/');

    data = data -50;   
    N = ndims(data);
    nobservations = size(data,1);
    for n = 2:N
        nvariable(n-1) = size(data,n);
    end
    if ~exist('tail') || strcmp(tail,'right') %if two tail t-test
    func = '';
        else %if two-sided t-test
    func = 'abs';
    end
    
    cln = repmat({':'},1,N-1);
    
    %create permutation samples and convert them to pvalues (perms x variable1 x variable2 x ...)
    if ~exist('StatMapPermPV') %if pvalues have not been precomputed     
        StatMapPerm = single(zeros([n_perm nvariable]));

        %first permutation sample is original data
        StatMapPerm(1,cln{:}) = mean(data,1) ./ std(data);    

        %perform permutations for n-1 permutations as the first permutation
        %is the orginal data 
        for i = 2:n_perm %par
            if ~rem(i,100)
                disp(['Create permutation samples: ' num2str(i) ' out of ' num2str(n_perm)]);
            end
            perm = single(sign(rand(nobservations,1)-0.5));
            data_perm = repmat(perm,1,nvariable(1)) .* data; 
            StatMapPerm(i,cln{:}) = mean(data_perm,1) ./ std(data_perm);         
        end    

        %convert to pvalues
        eval([ 'StatMapPermPV = (n_perm+1 - tiedrank(' func '(StatMapPerm)))/n_perm;' ]);    
    end

    clear StatMapPerm;

    %Perform fdr
    pvalues = squeeze(StatMapPermPV(1,cln{:}));
    [SignificantVariables,crit_p,~,adjusted_pvalues] = fdr_bh(pvalues,q_value,'pdep');

end

