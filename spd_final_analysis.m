%% Sequential Patch Deletion: Final Analysis
%
% This specifically compares the final meta-communities that remain after
% the SPD perturbation experiments are completed, comparing the overall
% (not the time-series) effect of different repeated perturbation types and
% reserve placements.

clear;
set(0,'DefaultAxesFontSize',15);

%% Load the data

% Loop over each of the 15 repeat simulations
count1 = 0;
individual_data_FN = cell(15,1);

for s1 = 100:1:120
    
    % Collect the data that has both the non-reserve and the targeted
    % reserve data corrected
    filename=sprintf('SPD Sequential Patch Deletion/Final/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_meta_robust_final_00010000_FIXED2.txt',s1);
    
    if (exist(filename,'file')==2)
        
        count1 = count1+1;
        
        % store the individual data in a cell array
        individual_placeholder=readmatrix(filename);
        individual_data_FN{count1}=individual_placeholder;
        clear individual_placeholder;
    end
end

clearvars s1 count1 filename;

%% Create set of dependent variables:

% change in global diversity
% change in average local diversity
% change in average (maximum) SCTL
% change in average (average) SCTL

% have already looked at bodysize effects etc. so this is really going to
% just focus on conservation of diversity against different types of
% perturbation and the role of reserves on just this outcome

% initial data is found in the first line of the file, so an additional
% file is not required

dependent_var = zeros(15,2,2,2,12,11,4); % (simulation, permanence, effect, reserves?, reserve (or only) sequence, deletion sequence, property)
dependent_prop_cols = [10, 12, 13, 16];


for sim = 1:1:15
    
    % open the cell row
    
    for perm = 1:1:2
        for effect = 1:1:2
            for is_reserves = 1:1:2
                
                reserve_seq_length = [1; 12];
                
                for reserve_seq = 1:1:reserve_seq_length(is_reserves)
                    for deletion_seq = 1:1:11
                        
                        row = 1+286*(perm-1)+143*(effect-1)+11*(is_reserves-1)+11*(reserve_seq-1)+deletion_seq;

                        for prop = 1:1:4

                            col = dependent_prop_cols(prop);
                            dependent_var(sim,perm,effect,is_reserves,reserve_seq,deletion_seq,prop) = (individual_data_FN{sim}(row,col)-individual_data_FN{sim}(1,col))/individual_data_FN{sim}(1,col);
                        end
                    end
                end
            end
        end
    end
end
clearvars row prop col deletion_seq reserve_seq_length;

%% Plot ALL



% Plot effect by No/Reserves(placement) for each disruption scenario.

perm_str = {'temporary', 'permanent'};
effect_str = {'elimination', 'displacement'};
bp_labels = {'No reserves', 'Remote block', 'Central block', 'Remote line', 'Central line', ...
        'Remote two blocks', 'Central two blocks', 'Three blocks', 'Max-dispersal individuals',...
        'Med-dispersal individuals', 'Low-dispersal individuals', 'Highest diversity', 'Lowest range'};
deletion_type_str = {'random','targeted'};
deletion_type_sel = {1:10,11};

is_boxplots = false;
if (is_boxplots)
    for deletion_type = 1:1:2
        for perm = 1:1:2
           for effect = 1:1:2

               hFig=figure('visible','on');
               clf(hFig);
               hold off;
               x_loc  = 1;

               %%%%%%%%%%%%%%%%%%%%%%%%
               % No reserves
               A = reshape( squeeze(dependent_var(:,perm,effect,1,1,deletion_type_sel{deletion_type},1)), 1, []);
               boxplot(A,'positions',x_loc);
               hold on;

               %%%%%%%%%%%%%%%%%%%%%%%%
               % Reserves by placement

               for is_reserves = 2:1:2

                    reserve_seq_length = [1; 12];

                    for reserve_seq = 1:1:12
                        x_loc = x_loc+1;
                        A = reshape( squeeze(dependent_var(:,perm,effect,is_reserves,reserve_seq,deletion_type_sel{deletion_type},1)), 1, []);
                        boxplot(A,'positions',x_loc);
                        hold on;
                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%
                title_str = sprintf('Effect of %s %s %s sequences', deletion_type_str{deletion_type}, perm_str{perm}, effect_str{effect});
                title(title_str);
                ylabel('Fractional change in global diversity')
                xlabel('Reserve placement')
                xticks([1:13]);
                xticklabels(bp_labels);
                xtickangle(45);
                xlim([0 14])
                ylim([-1.01 -0.6])
                printname=sprintf('SPD Sequential Patch Deletion/FIGURES/SPD_final_reserves_%d_%d_%d',deletion_type,perm,effect);
                print(printname,'-dpng','-r400','-painters');
                close(hFig);

            end
        end
    end
end

%% ANOVA:
%
% Use one-way ANOVA to quantify differences within these scenarios
% (between different reserve placements)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = cell(2,2,2); % contains the overall p-value for whether or not there is "a difference" between ANY of the 13 no/reserve sets
t = cell(2,2,2); % contains the statistics table output with each ANOVA test
c = cell(2,2,2); % contains the results of multcompare - column 6 has the p-value for each pair of reserve types having the same mean.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for deletion_type = 1:1:2
    for perm = 1:1:2
       for effect = 1:1:2

           %%%%%%%%%%%%%%%%%%%%%%%%
           % Begin ANOVA
           y = [];
           g = [];
           set_num = 0;
           for is_reserves = 1:1:2
                reserve_seq_length = [1; 12];
                for reserve_seq = 1:1:reserve_seq_length(is_reserves)
                    
                    set_num = set_num+1;
                    new_y = transpose(reshape( squeeze(dependent_var(:,perm,effect,is_reserves,reserve_seq,deletion_type_sel{deletion_type},1)), 1, []));
                    data_size = size(new_y);
                    data_len = data_size(1);
                    new_g = transpose(linspace(set_num,set_num,data_len));
                    
                    y = [y; new_y];
                    g = [g; new_g];
                    
                    clearvars data_len data_size new_y new_g;
                end
           end
           [p{deletion_type,perm,effect}, t{deletion_type,perm,effect}, stats] = anova1(y,g);
           
           % Save the figure:
           % title_str = sprintf('Effect of %s %s %s sequences', deletion_type_str{deletion_type}, perm_str{perm}, effect_str{effect});
           % title(title_str);
           ylabel('Fractional change in global diversity')
           xlabel('Reserve placement')
           xticks([1:13]);
           xticklabels(bp_labels);
           xtickangle(45);
           xlim([0 14])
           ylim([-1.01 -0.6])
           printname=sprintf('SPD Sequential Patch Deletion/FIGURES/SPD_final_reserves_ANOVA_%d_%d_%d',deletion_type,perm,effect);
           print(printname,'-dpng','-r400','-painters');
           close(gcf);
           
           c{deletion_type,perm,effect} = multcompare(stats);
           close(gcf);
           
           clearvars y g set_num printname title_str;
           
           %
           %%%%%%%%%%%%%%%%%%%%%%%%
           
        end
    end
end   

clearvars deletion_type perm effect;




