

clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus'
user_settings
cd(W)
savemat = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/2_pmaps/';

load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/2_pmaps/feed_pmap_code.mat')

for ii = 4:size(feed,1)
    tic
    %----------------------------------------------------------------
    % GET BASIC PARAMETERS FROM STRUCTURE
    %----------------------------------------------------------------
                        clearvars -except ii path feed savemat sel_row feedf setting slope nfish VSDI ref_movie movies reject_on outfield trial_kinds    
nfish = feed{ii,1};
    VSDI = TORus('load', nfish) ;
                    
    for ref_movi=  1:numel(feed{ii,4})
        ref_movie = feed{ii,4}{ref_movi};

        for reji=  1:numel(feed{ii,5})
            reject_on = feed{ii,5}(reji);
            
            for fieldi = 1
                outfield = 'peakmean';
                
                for condi = 1:size(feed{ii,3},1) % blank conditions needs to be included and be the first for the code to properly work
                    trial_kinds = feed{ii,3}(condi,:);
                                        
                  
                    if strcmpi(outfield, 'peakmean')
                        flagpeakmean = 1;
                    else
                        flagpeakmean = 0;
                    end

                    
                    %----------------------------------------------------------------
                    % START CODE
                    %----------------------------------------------------------------
                    % SETTINGS THAT ALLOW LOADING STRUCTURE
                    n_perm = 1000;  %@ SET
                    
                    %----------------------------------------------------------------
                    % LOAD THE STRUCTURE WITH PERMUTATION RESULTS (if exists)
                    %----------------------------------------------------------------
                    
                    name = [num2str(VSDI.ref) '_TFCE' num2str(n_perm) 'rep_reject' num2str(reject_on), 'data', ref_movie   ,'_Pmaps.mat'];
                    
     
                        load(fullfile(savemat,name))
                 
                        
                        if isfield(TFCEmaps, 'measurelist')
                            jj = numel(TFCEmaps.measurelist);
                        else
                            jj = 0;
                        end
                        TFCEmaps.measurelist{1,jj+1} = 'peakmean' outfield ;

                    save(fullfile(savemat,name),'TFCEmaps')
                    
                    % PLOT STATISTICAL DIFFERENCE BETWEEN CONDITIONS (P-THRESHOLDED)
                    
                    %     % PLOT WITH A THRESHOLD OF 0.001 IN THE second ROW
                    %     pthresh = 0.05; %to threshold out when plotting
                    % %     custom_map = colormap_loadBV();
                    %     custom_map = jet;
                    %
                    %         clims= [];
                    %
                    %     for condi= 1:length(TFCEmaps.codedef)
                    %         maps = TFCEmaps.(outfield);
                    %
                    %         ax(1) = subplot(1, 2, 1); %can plot only 3 activity
                    %         backgr = VSDI.backgr(:,:,VSDI.nonanidx(1));
                    %         alpha = TFCEmaps.(outfield).Pmap(:,:,ploti)< pthresh;
                    %         %    imagesc(imtiles(:,:,ploti))
                    % %         plot_logpmap_overlaid(maps.Pmap(:,:,condi),VSDI.backgr(:,:,VSDI.nonanidx(1)),pthresh ,0, ax(ploti), 1, flipud(parula));
                    %         plot_framesoverlaid2(TFCEmaps.(outfield).diffmap(:,:,ploti),backgr, alpha  ,0, ax(ploti), clims, 1, [], custom_map);
                    %
                    % %         idxcond = find(VSDI.condition(:,1) == newkinds(condi)); idxcond = idxcond(1);
                    % %         tit= ['c=',num2str(newkinds(condi)), '(', num2str(VSDI.condition(idxcond,4)), 'mA)'];
                    %         tit = TFCE.codedef{condi};
                    %         set(ax(ploti),'XColor', 'none','YColor','none')
                    %
                    %         ax(1).Title.String = tit;
                    %         ax(1).Visible = 'off';
                    %
                    %         ax(2) = subplot(1,2,2);
                    %         map = maps.meanmap(:,:,condi);
                    %
                    %         imagesc(map); axis image
                    %         title([outfield ' mean map'])
                    %         ax(2).Visible = 'off';
                    %
                    %     end
                    %     sgtit = [num2str(VSDI.ref) ':' num2str(pthresh) ' p-thresh ' TFCE.map_measure '(TFCE n' TFCEmaps.perm ') rej' num2str(reject_on)];
                    %     sgtitle(sgtit)
                    %
                    %     %             ax(9) = subplot(3,3,9)
                    %     %             imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)))
                    %     %             colormap(ax(9), bone)
                    %     %             colorbar('off')
                    %     %             ax(9).Visible = 'off';
                    %     %             axis image
                    %
                    %         localname = ['P_Map_thresh' num2str(VSDI.ref) 'block' num2str(condB)  ' (TFCEperm' num2str(n_perm)  'rep) of_ ' outfield  '_p' num2str(pthresh) '_' ref_movie 'reject' num2str(reject_on) '.jpg'];
                    %
                    %     %SAVE
                    %     saveas(gcf, fullfile(savein,localname ), 'jpg')
                    %     close all
                    %     %             clear maps  name pathsave moviesA moviesB idxA idxB nA nB cond_def
                    %     blob()
                end %trial_kinds
 
            end % field
        end %rejection
    end
end % nfishi
% blob(); pause(0.1); blob()
% end %P-MAPS EXTRACTION LOOP
% ----------------------------------------------
