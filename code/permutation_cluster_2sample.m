function [significantVarWei,significantVarMax,pValWei,pValMax,clusters] = permutation_cluster_2samples_weight_alld_MH (data1, data2, n_perm, cluster_th, significance_th, tail)
% Based on permutation_cluster_1sample_weight_alld_MH
% requires two data sets with same dimensions
%
%Performs one-sided (>0) or two-sided cluster-size/weight test on the 'data'.
% The data array can have any number of dimensions (1D, 2D, 3D, 4D etc) in
% the order [observations x variable1 x variable2 x ...]
% Data from each observation are randomly multiplied by +-1 to create
% permutation samples. These samples are converted to pvalues and then the
% cluster_th threshold (in pvalue units) is applied to identify suprathreshold clusters. The
% distribution of the size and the weights of suprathreshold clusters is used to assign statistical
% significance to the clusters (also in pvalue units) of the original data.
%
% INPUT:
%   data: observations x variable1 x variable2 x variable3 x ... (supports all dimemsions)
%   nperm: number of permutations
%   cluster_th: cluster defining threshold (in pvalue units)
%   significance_th: significance threshold (alpha value)
%   tail: string 'right' or 'both'. Default is 'right'.
%
% % for example:
%   load data;
%   nperm = 1000;
%   cluster_th = 0.05;
%   significance_th = 0.05;
%   tail = 'right';
%
% OUTPUT:
%   significantVarMax:  significant clusters (maximum cluster size)
%   significantVarWei:  significant clusters (maximum cluster weights)
%   clusters:           clusters above cluster defining threshold
%   pValMax:            p-value of maximum cluster size found in original data
%   pValWei:            p-value of maximum cluster weight found in original data
%
% Author: Dimitrios Pantazis, December 2015
% Modified: Siying Xie, August 2020
% Modified: Marleen Haupt, February 2023 (added data = data -50)

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
            data_perm1 = data(perm==1,:,:);
            data_perm2 = data(perm==-1,:,:);
            StatMapPerm(i,cln{:}) = mean(data_perm1,1)-mean(data_perm2,1);         
        end    

        %convert to pvalues
        eval([ 'StatMapPermPV = (n_perm+1 - tiedrank(' func '(StatMapPerm)))/n_perm;' ]);    
    end

%find maximum cluster size and maximum weighted cluster for all permutation samples
if strcmp(tail,'left')
    [clusterMaxSize(1), clusterMaxWei(1), clusters, clustersize, clusterweight] = find_clusters_weight_alld(squeeze(StatMapPerm(1,cln{:})*-1), squeeze(StatMapPermPV(1,cln{:})<=cluster_th));
else
    [clusterMaxSize(1), clusterMaxWei(1), clusters, clustersize, clusterweight] = find_clusters_weight_alld(abs(squeeze(StatMapPerm(1,cln{:}))), squeeze(StatMapPermPV(1,cln{:})<=cluster_th));
end

for i = 2:n_perm
    if ~rem(i,100)
        disp(['Compute maximum cluster: ' num2str(i) ' out of ' num2str(n_perm)]);
    end
    if strcmp(tail,'left')
        [clusterMaxSize(i), clusterMaxWei(i)] = find_clusters_weight_alld(squeeze(StatMapPerm(i,cln{:})*-1), squeeze(StatMapPermPV(i,cln{:})<=cluster_th));
    else
        [clusterMaxSize(i), clusterMaxWei(i)] = find_clusters_weight_alld(abs(squeeze(StatMapPerm(i,cln{:}))), squeeze(StatMapPermPV(i,cln{:})<=cluster_th));
    end
end

%find cluster threshold
clusterMaxSize_sorted = sort(clusterMaxSize, 'descend');
clusterMaxWei_sorted = sort(clusterMaxWei, 'descend');
th_max = clusterMaxSize_sorted(n_perm*significance_th );
th_wei = clusterMaxWei_sorted(n_perm*significance_th );

%find significant variables
if length(nvariable) == 1
    significantVarMax = zeros(nvariable,1);
    significantVarWei = zeros(nvariable,1);
else
    significantVarMax = zeros(nvariable);
    significantVarWei = zeros(nvariable);
end

%apply threshold on found clusters
significantVarMax([clusters{clustersize>th_max}]) = 1;
significantVarWei([clusters{clusterweight>th_wei}]) = 1;
if ~isempty (clustersize)
    pValMax = ( find(clusterMaxSize_sorted == clusterMaxSize(1), 1, 'first') ) / length(clusterMaxSize_sorted);
    pValWei = ( find(clusterMaxWei_sorted == clusterMaxWei(1), 1, 'first') ) / length(clusterMaxWei_sorted);
    disp(['pValue(max) = ', num2str(pValMax), '; pValue(weighted) = ', num2str(pValWei), '.']);
else
    pValMax = NaN;
    pValWei = NaN;
    disp('No cluster found.');
end

end
