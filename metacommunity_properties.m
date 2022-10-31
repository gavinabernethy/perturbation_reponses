% Gather average initial properties of the 15 meta-communities
%
% Using data from:
%   - Local Initial data
%   - Robustness Correlation data (this has the true local diversity when
%   each cell is isolated).
% NOTE: there is some error when reading the robustness_correlation files:
% In almost all cases it reads an additional non-existent column at the
% start, that it fills with NaN's. To ensure that it does this with every
% file consistently, use the additional arguments in the readmatrix(). Then
% remove this column to avoid having to add an extra column index whenever
% accessing this data.

clear;
set(0,'DefaultAxesFontSize',15);
color_array={ [0.05, 0.65, 0.25], [247/256, 182/256, 2/256], [0.15, 0.8, 0.4], [61/256, 89/256, 252/256], [252/256, 88/256, 76/256] };


%% Collect data:
% Loop over each of the 15 repeat simulations

count = 0;
individual_data_LI=cell(15,1);
individual_data_RC=cell(15,6,6);

for s1 = 100:1:120
     
    filename=sprintf('Local Initial data/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_meta_robust_localinitial_00010000.txt',s1);
    
    if (exist(filename,'file')==2)
        
        count = count+1;
        
        % store the individual local initial data in a cell array
        individual_placeholder=readmatrix(filename);
        individual_data_LI{count}=individual_placeholder;
        clear individual_placeholder;
        
        % store the individual patch-specific meta-robustness data in a
        % cell array
        for s2 = 1:1:6
            for s3 = 1:1:6
                
                filename2=sprintf('Robustness Correlation data/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_robustnesscorrelation_%d_%d.txt',s1,s2,s3);
                individual_placeholder=readmatrix(filename2,'Delimiter',' ','ConsecutiveDelimitersRule','join');
                individual_data_RC{count,s2,s3}=individual_placeholder;
                clear individual_placeholder;
                
            end
        end
    end
end

% Remove the faulty extra column from the RobustnessCorrelation data:
individual_data_RC2 = cell(15,6,6);
for r1 = 1:1:15
    for r2 = 1:6
        for r3 = 1:1:6
            individual_data_RC2{r1,r2,r3} =  individual_data_RC{r1,r2,r3}(:,2:14);
        end
    end
end
clearvars individual_data_RC;
individual_data_RC = individual_data_RC2;
clearvars individual_data_RC2;

%% Average local diversity and other properties:
local_properties = zeros(6,6,9);
prop_list = {'local_diversity','max_SCTL','link_density','connectance',...
    'average_SCTL','omnivory_frac','ave_population','ave_bodysize','ave_range'};
for k1 = 1:1:15
    line = 0;
    for k2 = 1:1:6
        for k3 = 1:1:6
            line = line+1;
            for prop = 1:1:9
                local_properties(k2,k3,prop) = local_properties(k2,k3,prop) + individual_data_LI{k1}(line,13+prop);
            end
        end
    end
end
local_properties = local_properties/15.0;

% now combine the inner/middle/outer cells
local_properties_comb = zeros(6,6,9);
for prop = 1:1:9
    
    % inner
    local_properties_comb(3:4,3:4,prop) = sum(local_properties(3:4,3:4,prop), 'all')/4.0;
    % middle
    middle_ave = (sum(local_properties(2,2:5,prop), 'all') + sum(local_properties(5,2:5,prop), 'all') + sum(local_properties(3:4,2,prop), 'all') + sum(local_properties(3:4,5,prop), 'all'))/12.0;
    local_properties_comb(2,2:5,prop) = middle_ave;
    local_properties_comb(5,2:5,prop) = middle_ave;
    local_properties_comb(3:4,2,prop) = middle_ave;
    local_properties_comb(3:4,5,prop) = middle_ave;
    % outer
    outer_ave = (sum(local_properties(1,1:6,prop), 'all') + sum(local_properties(6,1:6,prop), 'all') + sum(local_properties(2:5,1,prop), 'all') + sum(local_properties(2:5,6,prop), 'all'))/20.0;
    local_properties_comb(1,1:6,prop) = outer_ave;
    local_properties_comb(6,1:6,prop) = outer_ave;
    local_properties_comb(2:5,1,prop) = outer_ave;
    local_properties_comb(2:5,6,prop) = outer_ave;
end
    
% Create the colormap limits of the minimum-maximum values (2 d.p.) for each
% property:
cmap_limits_1 = {[7 37], [1 5], [1 3.16], [0.14 0.53], [1 2.31], [0.1 0.65], [417 3069], [2.52 6.56], [1.65 4.6]};
for prop = 1:1:9
    hFig=figure('visible','on');
    clf(hFig);
    hold off;
    imagesc(local_properties_comb(:,:,prop));
    set(gca,'YDir','normal');
    colorbar;
    caxis(cmap_limits_1{prop});
    colormap(brewermap([],'Greys'));
    % caxis([0 30]);
    box on;
    xlabel('Patch $x$', 'Interpreter', 'Latex');
    ylabel('Patch $y$', 'Interpreter', 'Latex');
    printname=sprintf('Initial Metacommunity Properties/aggregated_cells_and_simulations/initial_%s',prop_list{prop});
    print(printname,'-dpng','-r400','-painters');
    close(hFig);
end    

%% True local diversity (after isolation)
true_local_diversity = zeros(6,6,3);
for k1 = 1:1:15
    for k2 = 1:1:6
        for k3 = 1:1:6
            true_local_diversity(k2,k3,1) = true_local_diversity(k2,k3,1) + individual_data_RC{k1,k2,k3}(1,3);
        end
    end
end
true_local_diversity(:,:,1) = true_local_diversity(:,:,1)./15.0;
% net increase over true local diversity
true_local_diversity(:,:,2) = local_properties(:,:,1) - true_local_diversity(:,:,1);
% fractional increase over true local diversity
true_local_diversity(:,:,3) = (local_properties(:,:,1) - true_local_diversity(:,:,1))./true_local_diversity(:,:,1);

% then convert to comb also
true_local_diversity_comb = zeros(6,6,3);
for prop = 1:1:3
    % inner
    true_local_diversity_comb(3:4,3:4,prop) = sum(true_local_diversity(3:4,3:4,prop), 'all')/4.0;
    % middle
    middle_ave = (sum(true_local_diversity(2,2:5,prop), 'all') + sum(true_local_diversity(5,2:5,prop), 'all') + sum(true_local_diversity(3:4,2,prop), 'all') + sum(true_local_diversity(3:4,5,prop), 'all'))/12.0;
    true_local_diversity_comb(2,2:5,prop) = middle_ave;
    true_local_diversity_comb(5,2:5,prop) = middle_ave;
    true_local_diversity_comb(3:4,2,prop) = middle_ave;
    true_local_diversity_comb(3:4,5,prop) = middle_ave;
    % outer
    outer_ave = (sum(true_local_diversity(1,1:6,prop), 'all') + sum(true_local_diversity(6,1:6,prop), 'all') + sum(true_local_diversity(2:5,1,prop), 'all') + sum(true_local_diversity(2:5,6,prop), 'all'))/20.0;
    true_local_diversity_comb(1,1:6,prop) = outer_ave;
    true_local_diversity_comb(6,1:6,prop) = outer_ave;
    true_local_diversity_comb(2:5,1,prop) = outer_ave;
    true_local_diversity_comb(2:5,6,prop) = outer_ave;
end

cmap_limits_2 = {[0 12], [0 2.67], [-0.73 0]};
for z=1:1:3
    hFig=figure('visible','on');
    clf(hFig);
    hold off;
    imagesc(true_local_diversity_comb(:,:,z));
    set(gca,'YDir','normal');
    colorbar;
    colormap(brewermap([],'Greys'));
    caxis(cmap_limits_2{z});
    box on;
    xlabel('Patch $x$', 'Interpreter', 'Latex');
    ylabel('Patch $y$', 'Interpreter', 'Latex');
    printname=sprintf('Initial Metacommunity Properties/aggregated_cells_and_simulations/initial_true_local_diversity_%d',z);
    print(printname,'-dpng','-r400','-painters');
    close(hFig);
end


%% Correlation and analysis of cell-level properties

% Correlation coefficients

prop_labels = {'Local Diversity','Max SCTL','Link Density','Connectance',...
    'Average SCTL','Omnivory Fraction','Average Population','Average Bodysize',...
    'Average Range','Centrality','Topology','True Local Diversity',...
    'Raw increase over true LD', 'Fractional increase over true LD', 'Fractional decrease over apparent LD'};

centrality = transpose([1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2, 3, 3, 2, 1, 1, 2, 3, 3, 2, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1]);
topology = transpose([2, 3, 3, 3, 3, 2, 3, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 3, 2, 3, 3, 3, 3, 2]);

% add these columns to each local_initial matrix
for j=1:1:15
    individual_data_LI{j}(:,23) = centrality;
    individual_data_LI{j}(:,24) = topology;
end

% add "true local diversity" as a column
for k1 = 1:1:15
    true_local_column = zeros(36,3);
    line = 0;
    for k2 = 1:1:6
        for k3 = 1:1:6
            line = line+1;
            true_local_column(line,1) = individual_data_RC{k1,k2,k3}(1,3);
            true_local_column(line,2) = individual_data_LI{k1}(line,14) - individual_data_RC{k1,k2,k3}(1,3);
            true_local_column(line,3) = (individual_data_LI{k1}(line,14) - individual_data_RC{k1,k2,k3}(1,3))/individual_data_RC{k1,k2,k3}(1,3);
            true_local_column(line,4) = (individual_data_RC{k1,k2,k3}(1,3) - individual_data_LI{k1}(line,14))/individual_data_LI{k1}(line,14);
        end
    end
    individual_data_LI{k1}(:,25) = true_local_column(:,1);
    individual_data_LI{k1}(:,26) = true_local_column(:,2);
    individual_data_LI{k1}(:,27) = true_local_column(:,3);
    individual_data_LI{k1}(:,28) = true_local_column(:,4);
end



% now loop over the 15 properties and plot and correlate them against each
% other
corr_coef_values=zeros(15,15);

type_list1=[repmat("double",1,15)];
correlation_table=table('Size',[15 15],'VariableNames',prop_labels,'RowNames',prop_labels,'VariableTypes',type_list1);
type_list2=[repmat("double",1,5)];
stats_table=table('Size',[15 5],'VariableNames',{'Minimum', 'Mean - 1SD', 'Mean', 'Mean + 1SD', 'Maximum'},'RowNames',prop_labels,'VariableTypes',type_list2);

prop_lim = {[5 40], [0 6], [0 4], [0.1 0.6], [0.5 2.5], [0 0.7], [0 3500], [2 7], [1 5], [0 4], [1 5], [0 35], [-2 15], [-0.5 3], [-1 0.5]};

for t1=1:1:15
    
    % extract the x-data into a single vector
    x_data = zeros(540,1);
    for sim = 1:1:15
        x_data(36*(sim-1)+1:36*sim) = individual_data_LI{sim}(:,t1+13);
    end
    
    for t2=1:1:15
        
        % extract the y-data into a single vector
        y_data = zeros(540,1);
        for sim = 1:1:15
            y_data(36*(sim-1)+1:36*sim) = individual_data_LI{sim}(:,t2+13);
        end

        corr_matrix=corrcoef(x_data,y_data);
        corr_coef_values(t1,t2)=corr_matrix(2,1); % Calculate the linear correlation coefficient
        correlation_table{t1,t2}=corr_matrix(2,1);

        hFig=figure('visible','on');
        clf(hFig);
        hold off;
        scatter(x_data,y_data,25,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',0.8);
        box on;
        xlabel(prop_labels{t1});
        ylabel(prop_labels{t2});
        xlim(prop_lim{t1});
        ylim(prop_lim{t2});
        legend(num2str(corr_coef_values(t1,t2),3),'Location','NorthEast');
        printname=sprintf('Initial Metacommunity Properties/cell_properties/initial_cell_correlations_%d_%d',t1,t2);
        print(printname,'-dpng','-r400','-painters');
        close(hFig);

    end
    
    % Mean, standard deviation, and range of each cell-level property across
    stats_table{t1,1} = min(x_data);
    stats_table{t1,2} = mean(x_data) - std(x_data);
    stats_table{t1,3} = mean(x_data);
    stats_table{t1,4} = mean(x_data) + std(x_data);
    stats_table{t1,5} = max(x_data);
    
    % Box-and-whisker diagrams of the distribution of each property:
    hFig=figure('visible','on');
    clf(hFig);
    hold off;
    box on;
    boxplot(x_data);
    ylabel(prop_labels{t1});
    set(gca,'xtick',[]);
    printname=sprintf('Initial Metacommunity Properties/property_distributions/initial_cell_property_distributions_%d',t1);
    print(printname,'-dpng','-r400','-painters');
    close(hFig);
end


%% Format the Correlations Table
writetable(correlation_table,'Initial Metacommunity Properties/correlation_table.xlsx','WriteRowNames',true);

%% Format the Statistics Table
writetable(stats_table,'Initial Metacommunity Properties/stats_table.xlsx','WriteRowNames',true);



%% Correlation of species-level properties
% This is done in "spatial_stability_analysis.m"








