%% Sequential Patch Deletion: Time-series analysis
%
% This compares the time series of a given perm/effect of perturbation
% over the course of the repeated experiment.

is_reset = true;

if (is_reset)
    clear;
    set(0,'DefaultAxesFontSize',15);

    %% Data collection and preparation:
    %
    % (Min/ave/max, elim/disp, random/targeted, reserve type, time-series)
    temp_noreserves = zeros(3,2,2,1,1001);
    perm_noreserves = zeros(3,2,2,1,37);
    temp_reserves = zeros(3,2,2,12,1001);
    perm_reserves = zeros(3,2,2,12,31);

    % set placeholder max values
    temp_noreserves(1,:,:,:,:) = 99999;
    perm_noreserves(1,:,:,:,:) = 99999;
    temp_reserves(1,:,:,:,:) = 99999;
    perm_reserves(1,:,:,:,:) = 99999;

    simulation_numbers = [100,101,102,103,104,105,106,107,108,109,110,111,112,115,118];
    reserve_seq_length = [1; 12];
    for sim = 1:1:15
        sim_num = simulation_numbers(sim);
        for perm = 1:1:2
            for effect = 1:1:2
                for is_reserves = 1:1:2
                    for reserve_seq = 1:1:reserve_seq_length(is_reserves)
                        for deletion_seq = 1:1:11

                            if (reserve_seq<10)
                                res_str = sprintf('0%d',reserve_seq);
                            else
                                res_str = sprintf('%d',reserve_seq);
                            end

                            if (deletion_seq<10)
                                del_str = sprintf('0%d',deletion_seq);
                            else
                                del_str = sprintf('%d',deletion_seq);
                            end

                            base_str = 'SPD Sequential Patch Deletion/Ongoing/WW_2020_01_1pt0e05_0pt6_mig1em2';
                            filename = sprintf('%s_%d_meta_robust_ongoing_%d_%d_%d_00010000_%s_%s.txt',base_str,sim_num,perm,effect,is_reserves,res_str,del_str);

                            % read in the individual data
                            individual_placeholder = importdata(filename);

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %
                            % store the needed properties

                            % (Min/ave/max, elim/disp, random/targeted, reserve type, time-series)

                            if (perm==1)
                                if (is_reserves==1)
                                    if (deletion_seq<11)
                                        temp_noreserves(1,effect,1,reserve_seq,1:1001) = min(squeeze(temp_noreserves(1,effect,1,reserve_seq,1:1001)),individual_placeholder(:,11));
                                        temp_noreserves(2,effect,1,reserve_seq,1:1001) = squeeze(temp_noreserves(2,effect,1,reserve_seq,1:1001))+individual_placeholder(:,11);
                                        temp_noreserves(3,effect,1,reserve_seq,1:1001) = max(squeeze(temp_noreserves(3,effect,1,reserve_seq,1:1001)),individual_placeholder(:,11));
                                    else
                                        temp_noreserves(1,effect,2,reserve_seq,1:1001) = min(squeeze(temp_noreserves(1,effect,2,reserve_seq,1:1001)),individual_placeholder(:,11));
                                        temp_noreserves(2,effect,2,reserve_seq,1:1001) = squeeze(temp_noreserves(2,effect,2,reserve_seq,1:1001)) + individual_placeholder(:,11);
                                        temp_noreserves(3,effect,2,reserve_seq,1:1001) = max(squeeze(temp_noreserves(3,effect,2,reserve_seq,1:1001)),individual_placeholder(:,11));
                                    end
                                else
                                    if (deletion_seq<11)
                                        temp_reserves(1,effect,1,reserve_seq,1:1001) = min(squeeze(temp_reserves(1,effect,1,reserve_seq,1:1001)),individual_placeholder(:,11));
                                        temp_reserves(2,effect,1,reserve_seq,1:1001) = squeeze(temp_reserves(2,effect,1,reserve_seq,1:1001)) + individual_placeholder(:,11);
                                        temp_reserves(3,effect,1,reserve_seq,1:1001) = max(squeeze(temp_reserves(3,effect,1,reserve_seq,1:1001)),individual_placeholder(:,11));
                                    else
                                        temp_reserves(1,effect,2,reserve_seq,1:1001) = min(squeeze(temp_reserves(1,effect,2,reserve_seq,1:1001)),individual_placeholder(:,11));
                                        temp_reserves(2,effect,2,reserve_seq,1:1001) = squeeze(temp_reserves(2,effect,2,reserve_seq,1:1001)) + individual_placeholder(:,11);
                                        temp_reserves(3,effect,2,reserve_seq,1:1001) = max(squeeze(temp_reserves(3,effect,2,reserve_seq,1:1001)), individual_placeholder(:,11));
                                    end
                                end
                            else
                                if (is_reserves==1)
                                    if (deletion_seq<11)
                                        perm_noreserves(1,effect,1,reserve_seq,1:37) = min(squeeze(perm_noreserves(1,effect,1,reserve_seq,1:37)),individual_placeholder(:,11));
                                        perm_noreserves(2,effect,1,reserve_seq,1:37) = squeeze(perm_noreserves(2,effect,1,reserve_seq,1:37)) + individual_placeholder(:,11);
                                        perm_noreserves(3,effect,1,reserve_seq,1:37) = max(squeeze(perm_noreserves(3,effect,1,reserve_seq,1:37)),individual_placeholder(:,11));
                                    else
                                        perm_noreserves(1,effect,2,reserve_seq,1:37) = min(squeeze(perm_noreserves(1,effect,2,reserve_seq,1:37)),individual_placeholder(:,11));
                                        perm_noreserves(2,effect,2,reserve_seq,1:37) = squeeze(perm_noreserves(2,effect,2,reserve_seq,1:37)) + individual_placeholder(:,11);
                                        perm_noreserves(3,effect,2,reserve_seq,1:37) = max(squeeze(perm_noreserves(3,effect,2,reserve_seq,1:37)),individual_placeholder(:,11));
                                    end
                                else
                                    if (deletion_seq<11)
                                        perm_reserves(1,effect,1,reserve_seq,1:31) = min(squeeze(perm_reserves(1,effect,1,reserve_seq,1:31)), individual_placeholder(:,11));
                                        perm_reserves(2,effect,1,reserve_seq,1:31) = squeeze(perm_reserves(2,effect,1,reserve_seq,1:31)) + individual_placeholder(:,11);
                                        perm_reserves(3,effect,1,reserve_seq,1:31) = max(squeeze(perm_reserves(3,effect,1,reserve_seq,1:31)), individual_placeholder(:,11));
                                    else
                                        perm_reserves(1,effect,2,reserve_seq,1:31) = min(squeeze(perm_reserves(1,effect,2,reserve_seq,1:31)), individual_placeholder(:,11));
                                        perm_reserves(2,effect,2,reserve_seq,1:31) = squeeze(perm_reserves(2,effect,2,reserve_seq,1:31)) + individual_placeholder(:,11);
                                        perm_reserves(3,effect,2,reserve_seq,1:31) = max(squeeze(perm_reserves(3,effect,2,reserve_seq,1:31)), individual_placeholder(:,11));
                                    end
                                end
                            end

                            %
                            %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            clear individual_placeholder;

                        end
                    end
                end
            end
        end
    end

    % Now the averages by simulation and by deletion sequence:
    % Random:
    temp_noreserves(2,:,1,:,:) = temp_noreserves(2,:,1,:,:)./150;
    perm_noreserves(2,:,1,:,:) = perm_noreserves(2,:,1,:,:)./150;
    temp_reserves(2,:,1,:,:) = temp_reserves(2,:,1,:,:)./150;
    perm_reserves(2,:,1,:,:) = perm_reserves(2,:,1,:,:)./150;
    % Targeted:
    temp_noreserves(2,:,2,:,:) = temp_noreserves(2,:,2,:,:)./15;
    perm_noreserves(2,:,2,:,:) = perm_noreserves(2,:,2,:,:)./15;
    temp_reserves(2,:,2,:,:) = temp_reserves(2,:,2,:,:)./15;
    perm_reserves(2,:,2,:,:) = perm_reserves(2,:,2,:,:)./15;

    % Discard un-needed iterators:
    clearvars -except perm_noreserves perm_reserves temp_noreserves temp_reserves;
end
    
%% Graphics:
%
% (Min/ave/max, elim/disp, random/targeted, reserve type, time-series)

% For each of the four kinds of perturbation and two kinds of targeting,
% show the min-mean-max of each of the 13 reserve types and no reserve on
% the same plot.

reserve_labels = {'No reserves', 'Remote block', 'Central block', 'Remote line', 'Central line', ...
        'Remote 2 blocks', 'Central 2 blocks', '3 blocks', 'Max-dispersal individuals',...
        'Med-dispersal individuals', 'Low-dispersal individuals', 'Highest diversity', 'Lowest range'};

% include min-max shading (1) or not (2)
for version = 1:1:2
    % Temporary or Permanent
    for perm = 1:1:2

        if (perm == 1)
            % temporary
            dummy_nores = temp_noreserves;
            dummy_res = temp_reserves;
            x_max = 1001;
            x_max_res = 1001;
            leg_loc = 'NorthEast';
        else
            % permanent
            dummy_nores = perm_noreserves;
            dummy_res = perm_reserves;
            x_max = 36;
            x_max_res = 31;
            leg_loc = 'SouthWest';
        end

        for effect = 1:1:2
            for del_type = 1:1:2

                % New plot
                hFig=figure('visible','on');
                clf(hFig);
                hold off;

                x_vec = transpose(linspace(0,x_max-1,x_max));
                x_vec_res = transpose(linspace(0,x_max_res-1,x_max_res));
                col_set = {'r',[0.3010, 0.7450, 0.9330],[200/256 200/256 200/256],[180/256 180/256 180/256],...
                    [160/256 160/256 160/256],[140/256 140/256 140/256],[120/256 120/256 120/256],[100/256 100/256 100/256],...
                    [80/256 80/256 80/256],[60/256 60/256 60/256],[40/256 40/256 40/256],[0 200/256 0],[127/256 0 255/256]};

                % No reserves
                % Ave
                plot(x_vec(1:x_max), log(squeeze(dummy_nores(2,effect,del_type,1,1:x_max))),'LineWidth',1.0,'LineStyle','-','color',col_set{1});
                hold on;
                
                if (version == 1)
                    % Min
                    plot(x_vec(1:x_max), log(squeeze(dummy_nores(1,effect,del_type,1,1:x_max))),'LineWidth',0.1,'LineStyle','-','color',col_set{1},'HandleVisibility','off');
                    % Max
                    plot(x_vec(1:x_max), log(squeeze(dummy_nores(3,effect,del_type,1,1:x_max))),'LineWidth',0.1,'LineStyle','-','color',col_set{1},'HandleVisibility','off');
                    % Fill
                    x2 = [transpose(x_vec(1:x_max)), fliplr(transpose(x_vec(1:x_max)))];
                    inBetween = [transpose(log(squeeze(dummy_nores(1,effect,del_type,1,1:x_max)))), fliplr(transpose(log(squeeze(dummy_nores(3,effect,del_type,1,1:x_max)))))];
                    fill(x2, inBetween, col_set{1},'FaceAlpha',0.3,'HandleVisibility','off');
                end
                % Then reserves of different kinds
                for res_type = 1:1:12
                    
                    % Ave
                    plot(x_vec_res(1:x_max_res), log(squeeze(dummy_res(2,effect,del_type,res_type,1:x_max_res))),'LineWidth',1.0,'LineStyle','-','color',col_set{res_type+1});
                    
                    if (version == 1)
                        % Min
                        plot(x_vec_res(1:x_max_res), log(squeeze(dummy_res(1,effect,del_type,res_type,1:x_max_res))),'LineWidth',0.1,'LineStyle','-','color',col_set{res_type+1},'HandleVisibility','off');
                        % Max
                        plot(x_vec_res(1:x_max_res), log(squeeze(dummy_res(3,effect,del_type,res_type,1:x_max_res))),'LineWidth',0.1,'LineStyle','-','color',col_set{res_type+1},'HandleVisibility','off');
                        % Fill
                        x2 = [transpose(x_vec_res(1:x_max_res)), fliplr(transpose(x_vec_res(1:x_max_res)))];
                        inBetween = [transpose(log(squeeze(dummy_res(1,effect,del_type,res_type,1:x_max_res)))), fliplr(transpose(log(squeeze(dummy_res(3,effect,del_type,res_type,1:x_max_res)))))];
                        fill(x2, inBetween, col_set{res_type+1},'FaceAlpha',0.3,'HandleVisibility','off');
                    end
                end

                % Print plot
                legend(reserve_labels,'Location',leg_loc);
                xlabel('Perturbation events');
                ylabel('log(Global diversity)');
                box on;
                xlim([0 x_max]);
                %ylim([36 550])
                printname=sprintf('SPD Sequential Patch Deletion/FIGURES/SPD_sequential_%d_%d_%d_%d',version,perm,del_type,effect);
                print(printname,'-dpng','-r400','-painters');
                close(hFig);

            end
        end
        clearvars dummy_res dummy_nores;
    end
end

        