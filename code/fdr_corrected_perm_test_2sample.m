function [SignificantVariables,crit_p,adjusted_pvalues] = fdr_corrected_perm_test_2samples_MH(data1, data2, n_perm, q_value,tail)

% IMPORTANT: 
% This function is based of fdr_corrected_perm_test_MH.m which comes from
% the following source.

% Author: Dimitrios Pantazis, December 2015
% Modified: Siying Xie, XX XX
% Modified: Agnessa Karapetian, June 2021
% Modified: Greta Häberle, April 2022
% Renamed: Marleen Haupt, February 2023

% And it is based on
% https://github.com/behinger/permtest/blob/master/permtest.m

% INPUT:
% data1: first sammple observations x variable1 x variable2 x variable3 x ... (supports all dimensions)
% data2: first sammple observations x variable1 x variable2 x variable3 x ... (supports all dimensions)
% important: BOTH DATA NEED SAME number of dimensions
% nperm: number of permutations
% q_value: alpha level for FDR correction
% tail: string 'right' or 'both'. Default is 'both'.

% OUTPUT:
% SignificantVariables:  time points with significant values
% pvalues: p-values of the data (ground truth)
% crit_p: critical (statictically decisive) p-value in FDR
% adjusted_pvalues: FDR corrected p-values of the ground truth

% add the path of the fdr_bh folder in your folder structure
addpath('/scratch/haum92/aging/eeg/code/stats/fdr_bh/');
   
    data1 = data1 -50;  
    data2 = data2 -50; 
    data = [data1;data2];
    
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
        StatMapPerm(1,cln{:}) = mean(data1,1)-mean(data2,1);    

        %perform permutations for n-1 permutations as the first permutation
        %is the orginal data 
        for i = 2:n_perm %par
            if ~rem(i,100)
                disp(['Create permutation samples: ' num2str(i) ' out of ' num2str(n_perm)]);
            end
         
            helpvector = [ones(1, size(data1,1)), (-1*ones(1,size(data2,1)))]';
            perm = helpvector(randperm(length(helpvector)));
            data_perm1 = data(perm==1,:);
            data_perm2 = data(perm==-1,:);
            StatMapPerm(i,cln{:}) = mean(data_perm1,1)-mean(data_perm2,1);         
        end    

        %convert to pvalues
        eval([ 'StatMapPermPV = (n_perm+1 - tiedrank(' func '(StatMapPerm)))/n_perm;' ]);    
    end

    clear StatMapPerm;

    %Perform fdr
    pvalues = squeeze(StatMapPermPV(1,cln{:}));
    [SignificantVariables,crit_p,~,adjusted_pvalues] = fdr_bh(pvalues,q_value,'pdep'); %pdep

end

