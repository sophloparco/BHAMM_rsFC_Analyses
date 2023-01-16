function [threshBSR,test ] = displayPLSresults( result, conditions, lv, nRoi, singValPlot )
%% displayPLSresults - 
%   Inputs :
%       result - result.mat file from the PLS analysis
%       conditions - conditions - 1xN containing list of conditions
%       lv - which Latent Variable you want to view
%       nRoi - number of ROIs
%       method: This option is PLS method that you used:
%			1. Mean-Centering Task PLS
%			2. Non-Rotated Task PLS
%			3. Regular Behavior PLS
%			4. Regular Multiblock PLS
%			5. Non-Rotated Behavior PLS
%			6. Non-Rotated Multiblock PLS




% load('rri_color_code');


numOfGroups=size(result.num_subj_lst,2);
num_conds=length(conditions);
method=result.method;
if  method == 3
    noOfBehavVars = size(result.stacked_behavdata,2);
end

if nargin == 4
    singValPlot=0;
end

if singValPlot ~=1 && singValPlot ~=0
    error('saveThresholdedMat can only be either 1 or 0');
end

%display singular values plot
if singValPlot == 1
    num_perm=500;
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
    title(['Permuted values greater than observed, ', num2str(num_perm), ' permutation tests']);
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
end

figure;
% display dominant contrast for lv 1
if method == 3 || method == 5
    numOfBehavVecs=size(result.stacked_behavdata,2);
    upperLim=result.boot_result.ulcorr-result.boot_result.orig_corr;
    lowerLim=result.boot_result.orig_corr-result.boot_result.llcorr;
    barResult=result.boot_result.orig_corr;
    for g=1:numOfGroups
        for k=1:num_conds
            bar_hdl = bar((g-1)*numOfBehavVecs*num_conds + [1:numOfBehavVecs] + numOfBehavVecs*(k-1), ...
                barResult((g-1)*numOfBehavVecs*num_conds + [1:numOfBehavVecs] + numOfBehavVecs*(k-1), ...
                lv)); hold on;
%             set(bar_hdl,'facecolor',color_code(k,:));
        end
    end
    errorbar(1:size(barResult,2),barResult(:,lv), lowerLim(:,lv),upperLim(:,lv), 'k.'); hold off;
    y_label='Correlations';
    plot_title='Correlations Overview ';
else
    upperLim=result.boot_result.ulusc-result.boot_result.orig_usc;
    lowerLim=result.boot_result.orig_usc-result.boot_result.llusc;
    barResult=result.boot_result.orig_usc;
    
       
    for g=1:numOfGroups
        for k=1:num_conds
            bar_hdl = bar((g-1)*num_conds + k,barResult((g-1)*num_conds + k,lv)); hold on;
%             set(bar_hdl,'facecolor',color_code(k,:));
            
        end
    end
    errorbar(1:size(barResult,2),barResult(:,lv), lowerLim(:,lv),upperLim(:,lv), 'k.'); hold off;
    y_label='Brain Scores';
    plot_title='Task PLS Brain Scores with CI ';
    
    
end

% [l_hdl, o_hdl] = legend(conditions, 0);
% [l_hdl, o_hdl] = legend(conditions); %matlab 2019 version
% legend_txt(o_hdl);
% set(l_hdl,'color',[0.9 1 0.9]);
% setappdata(gca,'LegendHdl2',[{l_hdl} {o_hdl}]);

% set(gca, 'XGrid', 'on');
set(gca, 'GridLineStyle', '--');
% 
% 
% xlabel('Groups');
% ylabel(y_label);
% if method == 1
%     set(gca,'XTick',([2:numOfGroups] - 1)*num_conds +0.5)
% elseif method ==3
%     set(gca,'XTick',([2:numOfGroups] - 1)*num_conds*noOfBehavVars +0.5)
% end
% set(gca,'XTickLabel',1:numOfGroups);
% title([plot_title, 'of LV: ', num2str(lv)]);



bootstrap_ratio_lv1=result.boot_result.compare_u(:,lv);
% z_bsr=0.5*log((1+bootstrap_ratio_lv1)./(1-bootstrap_ratio_lv1));

b= triu(ones(nRoi),1);
b(b==1)=bootstrap_ratio_lv1;
test=b'+b;

%Threshold BSR
% pctl=prctile(bootstrap_ratio_lv1,95)
pctl=3.00
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



% 
% cbar=create_colourbar(bootstrap_ratio_lv1, abs(pctl), -abs(pctl));
% cbar_raw=create_colourbar(bootstrap_ratio_lv1, 0, 0);
% if nRoi == 64
%     figure; imagesc(test); hold on; plot([0,nRoi+10],[25,25],[25,25],[0,nRoi+10],[nRoi+10,0], [nRoi+10,0], 'Color', 'k','LineWidth', 2) ; hold off; colormap(cbar_raw); colorbar
%     title('Raw BSR')
%     figure; imagesc(threshBSR); hold on; plot([0,nRoi+10],[25,25],[25,25],[0,nRoi+10],[nRoi+10,0], [nRoi+10,0], 'Color', 'k','LineWidth', 2) ; hold off;
%     title(['BSR thresholded at ' num2str(abs(pctl)) ' (p-value = ' num2str(pVal) ')'])
% else
%     figure;imagesc(test, [-5, 5]); colormap(cbar_raw); colorbar 
%      title('Raw BSR')
%     figure; imagesc(threshBSR);
%     title(['BSR thresholded at ' num2str(abs(pctl)) ' (p-value = ' num2str(pVal) ')'])
% end
% colormap(cbar)
% colorbar
end

