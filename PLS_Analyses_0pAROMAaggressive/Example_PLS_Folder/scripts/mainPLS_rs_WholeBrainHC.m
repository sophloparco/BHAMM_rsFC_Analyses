%% Script to set options and run PLS analysis.
clc;
% clear all;
close all;
currDir=pwd;

addpath('~/Documents/MATLAB/spm12');
addpath('~/Documents/MATLAB/Pls');
addpath('~/Documents/MATLAB/Pls/plscmd');

% addpath(genpath('/data/laborajah2/users/Lizzy/functionalConnectivity_Lifespan/scripts/'));

%% Paths
behavPath='path/to/PLS_Analysis';
matInpath='path/to/Connectivity_Matrices';
% groups={'Post_Meno_Matrices', 'Pre-Peri_Meno_Matrices'};

%% Age group samples
% Set up subjects
allSubjectList={};

%separate subjects for groups
group1={};
group2={};
group3={};
%% 




subjectList = {allSubjectList}; %Single group
% subjectList= {group1, group2, group3}; % By groups

conditions = {'Rest'};

nsubj=cellfun('length',subjectList);
ncond=numel(conditions);

nRoi = 204; %no of ROis
behav_vector=spm_load(fullfile(behavPath, 'tab-delimited_file_with_participants_IDs_and_variables_of_interest.txt'));
behav_vector=behav_vector(:,['Columns #s for variables of interest']);
% contrast_vector=spm_load('NRBPLS_EncRetME_Contrast.txt'); % For Non-Rotated PLS

flag=0; %set this to 1 to stack Z map (Fisher R-Z transformed)

% Call script to stack datamat
% disp('Stacking datamats....');
datamat = stackPLSdatamat(subjectList, conditions, matInpath, nRoi,flag);


%% Options to set before running PLS analysis
% few options for the analysis, for more help see pls_analysis.m
% method: This option will decide which PLS method that the
% program will use:
%			1. Mean-Centering Task PLS
%			2. Non-Rotated Task PLS
%			3. Regular Behavior PLS
%			4. Regular Multiblock PLS
%			5. Non-Rotated Behavior PLS
%			6. Non-Rotated Multiblock PLS
%  num_boot - no. of bootstraps
%  num_perm - no of permutations

option.method = 3;
option.meancentering_type = 1;
option.num_boot = 500; % Number of bootstraps
option.num_perm = 500; % Number of permutations
option.stacked_behavdata=behav_vector; % Set this for BPLS
% option.stacked_designdata=contrast_vector; % Set this for Non-Rotated PLS

outFileName='PLS_WBHC_result'; % Filename for results.mat

% run pls
disp('Running PLS analysis...');
result = pls_analysis(datamat,nsubj,ncond,option);

% save result.mat to current working directory
disp('Saving result file...');
save(strcat(currDir, '/', outFileName, '.mat'), 'result');

%% Looking at results

% Plot p-values
pval = result.perm_result.sprob
nLV=numel(pval);
figure;
bar(pval,'r');
hold on;
h = zeros(nLV,1);
for i=1:nLV
    h(i)=plot(NaN,NaN, '.r');
end
legend(h,strcat('LV', num2str([1:nLV]'), {' - '} ,num2str(pval)));
title(['Permuted values greater than observed, ', num2str(option.num_perm), ' permutation tests']);
hold off;

% Plot effect sizes (% crossblock covariance)
pcov = result.s.^2 / sum(result.s.^2)
figure;
bar(pcov);
hold on;
h = zeros(nLV,1);
for i=1:nLV
    h(i)=plot(NaN,NaN, '.');
end
legend(h,strcat('LV', num2str([1:nLV]'), {' - '} ,num2str(pcov*100), '%'));
title('Percent Crossblock covariance');
hold off;

% Display Results for each LV (displays task LVs and BSR gris plots)
% disp('Displaying Figures...')
displayPLSresults(result,conditions,1,nRoi);
% displayPLSresults(result,conditions,2,nRoi);
% displayPLSresults(result,conditions,3,nRoi);

% save rawBSR matrix as a text file
disp('Saving raw BSR...')
saveOutputFiles(result,1,nRoi,outFileName);
% saveOutputFiles(result,2,nRoi,outFileName);
% saveOutputFiles(result,3,nRoi,outFileName);
% saveOutputFiles(result,4,nRoi,outFileName);


% TODO: Compute mean connectivity values. output txt file of mean R



%% Caluclate loading
stacked_datamat = vertcat(datamat{:});

u_loading=corr(stacked_datamat, result.usc);
z_u_loading=0.5*log((1+u_loading)./(1-u_loading));

%% display u_loading
lv =1;
loading_lv=u_loading(:,lv);

b= triu(ones(nRoi),1);
b(b==1)=loading_lv;
rawLoadingMat=b'+b; %raw loading matrix


dlmwrite(strcat(pwd,'/', outFileName, '_rawLoadingMat_LV',num2str(lv),'.csv'), rawLoadingMat, ',');

% lv =2;
% loading_lv=u_loading(:,lv);
% 
% b= triu(ones(nRoi),1);
% b(b==1)=loading_lv;
% rawLoadingMat=b'+b; %raw loading matrix
% 
% 
% dlmwrite(strcat(pwd,'/', outFileName, '_rawLoadingMat_LV',num2str(lv),'.csv'), rawLoadingMat, ',');

