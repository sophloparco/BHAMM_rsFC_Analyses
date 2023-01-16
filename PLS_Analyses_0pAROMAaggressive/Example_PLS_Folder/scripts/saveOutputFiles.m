function [  ] = saveOutputFiles( result,lv,nRoi,outFileName, saveThresholdedMat, pctl )
%% saveOutputFiles - writes the raw BSR values as a nRoi x nRoi matrix into a txt file (with suffix _rawBSR)
%   Inputs :
%       result - result.mat file from the PLS analysis
%       lv - which Latent Variable you want to save
%       nRoi - number of ROIs
%       outFileName - prefix of the output txt file

bootstrap_ratio_lv1=result.boot_result.compare_u(:,lv);
if nargin < 5
    saveThresholdedMat =0;
    pctl=prctile(bootstrap_ratio_lv1,95)
end
saveThresholdedMat
if saveThresholdedMat ~=1 && saveThresholdedMat ~=0
    error('saveThresholdedMat can only be either 1 or 0');
end
%Threshold BSR
% % % % % pctl=prctile(bootstrap_ratio_lv1,95)
% % % % % pctl=3.28
bootstrap_ratio_lv1=result.boot_result.compare_u(:,lv);
b= triu(ones(nRoi),1);
b(b==1)=bootstrap_ratio_lv1;
test=b'+b;
[~,outFileName,~]=fileparts(outFileName);
m=0;sig=1;
pVal=(1+erf((abs(pctl)-m)/(sqrt(2)*sig)))/2;
pVal=(1-pVal)*2
negBSR=bootstrap_ratio_lv1(bootstrap_ratio_lv1<0);
pctl_neg=prctile(negBSR,95)
posBSR=bootstrap_ratio_lv1(bootstrap_ratio_lv1>0);
pctl_pos=prctile(posBSR,95)
threshBSR=zeros(size(test));

for i = 1:nRoi
    for j=1:nRoi
        if test(i,j) < 0 && test(i,j) < -(abs(pctl))
            threshBSR(i,j)=test(i,j);
            
        elseif test(i,j) > 0 && test(i,j) > abs(pctl)
            threshBSR(i,j)=test(i,j);
            
        end
    end
end

strPCTL=strrep(num2str(pctl), '.', '');
if saveThresholdedMat == 0
    dlmwrite(strcat(pwd,'/', outFileName, '_rawBSR_LV',num2str(lv),'.csv'), test, ',');
else
    dlmwrite(strcat(pwd,'/', outFileName, '_thresh',strPCTL, 'BSR_LV',num2str(lv),'.csv'), threshBSR, ',');
end


end

