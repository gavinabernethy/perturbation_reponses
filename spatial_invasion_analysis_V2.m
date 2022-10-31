%% Spatial Invasion (SSI) analysis
% Code to analyse the spatial_invasion output file from this subroutine.

% Version 2 includes the data from:
% TEST4 - which looks at the distribution of re-invasions when there are no
% other non-resource species.
% TEST5 - which looks at the ability of randomly-generated species to
% invade the meta-community.
% Both of these should be compared to the original experiments.

% This V2 script omits much of the other analysis, so DO NOT DELETE the
% original _analysis.m script!

clear;
set(0,'DefaultAxesFontSize',15);
print_images=true; 

color_array={ [0.05, 0.65, 0.25], [247/256, 182/256, 2/256], [0.15, 0.8, 0.4], [61/256, 89/256, 252/256], [252/256, 88/256, 76/256] };

%% Collect data:
% Loop over each of the 15 repeat simulations

count = 0;
individual_data=cell(15,1);
individual_data_SSS=cell(15,1);
total_data=zeros(0,24);
lengths_store=zeros(15,1);

for s1 = 100:1:120
    
    filename=sprintf('SSI Species Spatial Invasion/DATA/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_spatial_invasion_00010000.txt',s1);
    
    if (exist(filename,'file')==2)
        
        count = count+1;
        
        % store the individual data in a cell array
        %
        % N x 24, where N is the number of non-resource species in the
        % metacommunity (i.e. n_spec - 36)
        individual_placeholder=readmatrix(filename);
        individual_data{count}=individual_placeholder;
        
        % concatenate them into a single vector
        current_length = size(total_data,1);
        add_length = size(individual_data{count},1);
        lengths_store(count)=size(individual_data{count},1);
        new_length = current_length+add_length;
        total_data_placeholder = zeros(new_length,24);
        total_data_placeholder(1:current_length,:)=total_data;
        total_data_placeholder(current_length+1:new_length,:)=individual_placeholder;
        clearvars total_data;
        total_data = total_data_placeholder;
        clearvars total_data_placeholder;
        
        % load the individual files for spatial_stability, as we need the
        % results of this for the initial properties of the metacommunity
        % before the species is re-introduced everywhere
        filename=sprintf('SSS Spatial Stability data/DATA/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_spatial_stability_00010000.txt',s1);
        individual_placeholder=readmatrix(filename);
        individual_data_SSS{count}=individual_placeholder;
        
    end
end

clearvars -except individual_data individual_data_SSS total_data lengths_store print_images color_array;

%% %% Modification of the response variables to "fractional" rather than absolute change:
% This needs to be fractional relative to the situation AFTER the species
% was deleted, which is the resulting information obtained during the SSS
% experiments.

control_prop_num=[4,8,9,10,11,12,13]; % 7 inputs
control_prop_name={'Range','Biomass','Population','SCTL (patch average)','SCTL (population average)',...
        'SCTL (biomass average)','Bodysize'};

response_prop_num=[5,16,17,18,19,20,21,22,23,24,7]; % 11 output properties
response_prop_num_init=[5,15,16,17,18,19,20,21,22,23]; % index for the first 10 (i.e. not including num_cells
% invaded, as this is kept absolute) in the spatial_stability_files that have the initial info
response_prop_name={'Relative secondary extinctions','Local diversity response','Max SCTL response','Link density response',...
    'Connectance response','Average SCTL response','Omnivory response','Population response','Bodysize response','Range response','Number of invasions'};

% for response properties, need the "fractional change" rather than the raw
% change or the resultant value:
individual_data_frac = individual_data;
total_data_frac = total_data;

% Modify the individual data:
for s2=1:1:15 % the simulation
    individual_data_frac{s2}(:,response_prop_num(1)) = (1+individual_data_SSS{s2}(2:end,6)+individual_data_SSS{s2}(2:end,response_prop_num_init(1)) - individual_data{s2}(:,response_prop_num(1)))./(individual_data_SSS{s2}(2:end,response_prop_num_init(1))-36);
    for s3=2:1:10 % not the 11th property (num_cells_invaded)
        % fractional change in properties = (new property - original property)/original property
        individual_data_frac{s2}(:,response_prop_num(s3)) = (individual_data{s2}(:,response_prop_num(s3)) - individual_data_SSS{s2}(2:end,response_prop_num_init(s3)))./individual_data_SSS{s2}(2:end,response_prop_num_init(s3));
    end
end

% Modify the concatenated data:
start_line=1;
end_line=0;
for s4=1:1:15 % simulation
    end_line=end_line+lengths_store(s4);
        
    total_data_frac(start_line:end_line,response_prop_num(1)) = (1+individual_data_SSS{s4}(2:end,6) + individual_data_SSS{s4}(2:end,response_prop_num_init(1)) - total_data(start_line:end_line,response_prop_num(1)))./(individual_data_SSS{s4}(2:end,response_prop_num_init(1))-36);
    for s5=2:1:10
        % fractional change in other properties = (new property - original property)/original property
        total_data_frac(start_line:end_line,response_prop_num(s5)) = (total_data(start_line:end_line,response_prop_num(s5)) - individual_data_SSS{s4}(2:end,response_prop_num_init(s5)))./individual_data_SSS{s4}(2:end,response_prop_num_init(s5));
    end
                 
    start_line=start_line+lengths_store(s4);
end

clearvars start_line end_line s2 s3 s4 s5;


%% Collect and modify the TEST4 and TEST5 data!


% TEST4
count = 0;
total_test4_data=zeros(0,24);
lengths_store=zeros(15,1);

for s1 = 100:1:120
    
    filename=sprintf('SSI Species Spatial Invasion/TEST4/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_spatial_invasion_00010000.txt',s1);
    
    if (exist(filename,'file')==2)
        
        count = count+1;
        individual_placeholder=readmatrix(filename);
        % concatenate them into a single vector
        current_length = size(total_test4_data,1);
        add_length = size(individual_placeholder,1);
        lengths_store(count)=size(individual_placeholder,1);
        new_length = current_length+add_length;
        total_data_placeholder = zeros(new_length,24);
        total_data_placeholder(1:current_length,:)=total_test4_data;
        total_data_placeholder(current_length+1:new_length,:)=individual_placeholder;
        clearvars total_data individual_placeholder;
        total_test4_data = total_data_placeholder;
        clearvars total_data_placeholder;    
    end
end


% TEST5
count = 0;
total_test5_data=zeros(0,24);
lengths_store=zeros(15,1);

for s1 = 100:1:120
    
    filename=sprintf('SSI Species Spatial Invasion/TEST5/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_spatial_invasion_00010000.txt',s1);
    
    if (exist(filename,'file')==2)
        
        count = count+1;
        individual_placeholder=readmatrix(filename);
        % concatenate them into a single vector
        current_length = size(total_test5_data,1);
        add_length = size(individual_placeholder,1);
        lengths_store(count)=size(individual_placeholder,1);
        new_length = current_length+add_length;
        total_data_placeholder = zeros(new_length,24);
        total_data_placeholder(1:current_length,:)=total_test5_data;
        total_data_placeholder(current_length+1:new_length,:)=individual_placeholder;
        clearvars total_data individual_placeholder;
        total_test5_data = total_data_placeholder;
        clearvars total_data_placeholder;    
    end
end


total_test5_data_frac = total_test5_data;
% Modify the concatenated data:
start_line=1;
end_line=0;
for s4=1:1:15 % simulation
    end_line=end_line+lengths_store(s4);
    % fractional secondary extinctions = (1+original property - new)/original property
    total_test5_data_frac(start_line:end_line,response_prop_num(1)) = (1.0+individual_data_SSS{s4}(1,response_prop_num_init(1))-total_test5_data(start_line:end_line,response_prop_num(1)))./individual_data_SSS{s4}(1,response_prop_num_init(1));      
    start_line=start_line+lengths_store(s4);
end




% TEST6
count = 0;
total_test6_data=zeros(0,24);
lengths_store=zeros(15,1);

for s1 = 100:1:120
    
    filename=sprintf('SSI Species Spatial Invasion/TEST6/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_spatial_invasion_00010000.txt',s1);
    
    if (exist(filename,'file')==2)
        
        count = count+1;
        individual_placeholder=readmatrix(filename);
        % concatenate them into a single vector
        current_length = size(total_test6_data,1);
        add_length = size(individual_placeholder,1);
        lengths_store(count)=size(individual_placeholder,1);
        new_length = current_length+add_length;
        total_data_placeholder = zeros(new_length,24);
        total_data_placeholder(1:current_length,:)=total_test6_data;
        total_data_placeholder(current_length+1:new_length,:)=individual_placeholder;
        clearvars total_data individual_placeholder;
        total_test6_data = total_data_placeholder;
        clearvars total_data_placeholder;    
    end
end

clearvars start_line end_line s1 s2 s3 s4 s5;


%% Analysis - Correlations:

collect_correlations = false;
if (collect_correlations)
     
    % Correlate modified responses against controls for the whole dataset:
   
    control_prop_xlim = {[0 7], [-5000 100000], [-2000 25000], [0 5], [0 5], [0 5], [0 15]};
    
    response_frac_prop_ylim = {[-0.01 0.15], [-0.1 0.08], [-0.1 0.12], [-0.08 0.1], [-0.08 0.1], [-0.08 0.1],...
        [-0.2 0.2], [-0.2 0.1], [-0.05 0.15], [-0.1 0.8], [0 40]};
    
    corr_coef_values = zeros(4,2);
    
    % In each of the following three cases, compare the original
    % experimental data with TEST5:
    
    % Number of invaded patches vs. bodysize
    hFig=figure('visible','on');
    clf(hFig);
    hold off;
    corr_matrix=corrcoef(total_data_frac(:,control_prop_num(7)),total_data_frac(:,response_prop_num(11)));
    corr_coef_values(1,1)=corr_matrix(2,1); % Calculate the linear correlation coefficient
    scatter(total_data_frac(:,control_prop_num(7)),total_data_frac(:,response_prop_num(11)),15,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',0.5,'MarkerFaceAlpha',0.9);
    box on;    
    xlabel(control_prop_name{7});
    ylabel(response_prop_name{11});
    xlim([0.5 13]);
    ylim([-1 37]);
    str1 = ['Resident species: ' num2str(corr_coef_values(1,1),3)];
    legend(str1,'Location','SouthEast');
    printname=sprintf('SSI Species Spatial Invasion/FIGURES/V2/spatial_invasion_correlation_%d_%d_resident',7,11);
    print(printname,'-dpng','-r400','-painters');
    close(hFig);
    
    hFig=figure('visible','on');
    clf(hFig);
    hold off;
    corr_matrix=corrcoef(total_test5_data_frac(:,control_prop_num(7)),total_test5_data_frac(:,response_prop_num(11)));
    corr_coef_values(1,2)=corr_matrix(2,1); % Calculate the linear correlation coefficient
    scatter(total_test5_data_frac(:,control_prop_num(7)),total_test5_data_frac(:,response_prop_num(11)),15,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',0.5,'MarkerFaceAlpha',0.9);
    box on;    
    xlabel(control_prop_name{7});
    ylabel(response_prop_name{11});
    xlim([0.5 13]);
    ylim([-1 37]);
    str2 = ['Random species: ' num2str(corr_coef_values(1,2),3)]; 
    legend(str2,'Location','SouthEast');
    printname=sprintf('SSI Species Spatial Invasion/FIGURES/V2/spatial_invasion_correlation_%d_%d_random',7,11);
    print(printname,'-dpng','-r400','-painters');
    close(hFig);
    
    % When there are no other species?
    hFig=figure('visible','on');
    clf(hFig);
    hold off;
    corr_matrix=corrcoef(total_test4_data(:,control_prop_num(7)),total_test4_data(:,response_prop_num(11)));
    corr_coef_values(2,1)=corr_matrix(2,1); % Calculate the linear correlation coefficient
    scatter(total_test4_data(:,control_prop_num(7)),total_test4_data(:,response_prop_num(11)),15,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',0.5,'MarkerFaceAlpha',0.9);
    box on;    
    xlabel(control_prop_name{7});
    ylabel(response_prop_name{11});
    xlim([0.5 13]);
    ylim([9 37]);
    str1 = ['Resident species: ' num2str(corr_coef_values(2,1),3)];
    legend(str1,'Location','SouthEast');
    printname=sprintf('SSI Species Spatial Invasion/FIGURES/V2/spatial_invasion_correlation_%d_%d_resident_empty',7,11);
    print(printname,'-dpng','-r400','-painters');
    close(hFig);
    
    hFig=figure('visible','on');
    clf(hFig);
    hold off;
    corr_matrix=corrcoef(total_test6_data(:,control_prop_num(7)),total_test6_data(:,response_prop_num(11)));
    corr_coef_values(2,2)=corr_matrix(2,1); % Calculate the linear correlation coefficient
    scatter(total_test6_data(:,control_prop_num(7)),total_test6_data(:,response_prop_num(11)),15,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',0.5,'MarkerFaceAlpha',0.9);
    box on;    
    xlabel(control_prop_name{7});
    ylabel(response_prop_name{11});
    xlim([0.5 13]);
    ylim([9 37]);
    str2 = ['Random species: ' num2str(corr_coef_values(2,2),3)]; 
    legend(str2,'Location','SouthEast');
    printname=sprintf('SSI Species Spatial Invasion/FIGURES/V2/spatial_invasion_correlation_%d_%d_random_empty',7,11);
    print(printname,'-dpng','-r400','-painters');
    close(hFig);
    
    %%
    % Relative change in global non-resource biodiversity vs. bodysize
    hFig=figure('visible','on');
    clf(hFig);
    hold off;
    corr_matrix=corrcoef(total_data_frac(:,control_prop_num(7)),total_data_frac(:,response_prop_num(1)));
    corr_coef_values(3,1)=corr_matrix(2,1); % Calculate the linear correlation coefficient
    scatter(total_data_frac(:,control_prop_num(7)),total_data_frac(:,response_prop_num(1)),10,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',0.5,'MarkerFaceAlpha',0.9);
    box on;
    xlabel(control_prop_name{7});
    ylabel(response_prop_name{1});
    xlim([0.5 13]);
    ylim([-0.002 0.135]);
    str1 = ['Resident species: ' num2str(corr_coef_values(3,1),3)];
    legend(str1, 'Location','NorthWest');
    printname=sprintf('SSI Species Spatial Invasion/FIGURES/V2/spatial_invasion_correlation_%d_%d_resident',7,1);
    print(printname,'-dpng','-r400','-painters');
    close(hFig);
    
    hFig=figure('visible','on');
    clf(hFig);
    hold off;
    corr_matrix=corrcoef(total_test5_data_frac(:,control_prop_num(7)),total_test5_data_frac(:,response_prop_num(1)));
    corr_coef_values(3,2)=corr_matrix(2,1); % Calculate the linear correlation coefficient
    scatter(total_test5_data_frac(:,control_prop_num(7)),total_test5_data_frac(:,response_prop_num(1)),10,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',0.5,'MarkerFaceAlpha',0.9);
    box on;
    xlabel(control_prop_name{7});
    ylabel(response_prop_name{1});
    xlim([0.5 13]);
    ylim([-0.002 0.135]);
    str2 = ['Random species: ' num2str(corr_coef_values(3,2),3)]; 
    legend(str2,'Location','NorthWest');
    printname=sprintf('SSI Species Spatial Invasion/FIGURES/V2/spatial_invasion_correlation_%d_%d_random',7,1);
    print(printname,'-dpng','-r400','-painters');
    close(hFig);

    
    %
    % final element (sec. extinctions vs. number of invaded patches, both of which are RESPONSES)
    hFig=figure('visible','on');
    clf(hFig);
    hold off;
    corr_matrix=corrcoef(total_data_frac(:,response_prop_num(11)),total_data_frac(:,response_prop_num(1)));
    corr_coef_values(4,1)=corr_matrix(2,1); % Calculate the linear correlation coefficient
    scatter(total_data_frac(:,response_prop_num(11)),total_data_frac(:,response_prop_num(1)),10,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',0.5,'MarkerFaceAlpha',0.9);
    box on;
    xlabel(response_prop_name{11});
    ylabel(response_prop_name{1});
    xlim([-1 37]);
    ylim([-0.002 0.135]);
    str1 = ['Resident species: ' num2str(corr_coef_values(4,1),3)];
    legend(str1, 'Location','NorthWest');
    printname=sprintf('SSI Species Spatial Invasion/FIGURES/V2/spatial_invasion_correlation_SPECIAL_1_11_resident');
    print(printname,'-dpng','-r400','-painters');
    close(hFig);
    
    hFig=figure('visible','on');
    clf(hFig);
    hold off;
    corr_matrix=corrcoef(total_test5_data_frac(:,response_prop_num(11)),total_test5_data_frac(:,response_prop_num(1)));
    corr_coef_values(4,2)=corr_matrix(2,1); % Calculate the linear correlation coefficient
    scatter(total_test5_data_frac(:,response_prop_num(11)),total_test5_data_frac(:,response_prop_num(1)),10,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',0.5,'MarkerFaceAlpha',0.9);
    box on;
    xlabel(response_prop_name{11});
    ylabel(response_prop_name{1});
    xlim([-1 37]);
    ylim([-0.002 0.135]);
    str2 = ['Random species: ' num2str(corr_coef_values(4,2),3)]; 
    legend(str2,'Location','NorthWest');
    printname=sprintf('SSI Species Spatial Invasion/FIGURES/V2/spatial_invasion_correlation_SPECIAL_1_11_random');
    print(printname,'-dpng','-r400','-painters');
    close(hFig);
end

%% 
% RE-SAMPLING

is_resampling = true;
if (is_resampling)
    
    % Generate the data
    iterate_datasets = {total_test5_data, total_test6_data};
    iterate_num_resamples = [20000, 20000];
    iterate_num_datapoints = [800, 5000];
    hold_resamples = cell(2,1);
    
    cc_mean = zeros(2,1);
    cc_cis = cell(2,1);
    distribution_stats = cell(2,1);
    
    for p1 = 1:1:2
        
        dataset = iterate_datasets{p1};
        num_resamples = iterate_num_resamples(p1);
        num_datapoints = iterate_num_datapoints(p1);
        spatial_invasion_resampling;
        hold_resamples{p1} = set_of_samples;
    
        %%%%%%%%%%%
        % For each resample, calculate the bodysize-invasion correlation coefficent:
        cc_store = zeros(num_resamples,1);
        for p2 = 1:1:num_resamples
            fprintf('p1 = %d, part A: p2 = %d\n', p1, p2);
            corr_matrix = corrcoef(hold_resamples{p1}{p2}(:,13),hold_resamples{p1}{p2}(:,7));
            cc_store(p2) = corr_matrix(2,1); % Calculate the linear correlation coefficient
        end

        % Calculate the mean and CI of the CC
        cc_mean(p1) = mean(cc_store);
        SEM = std(cc_store)/sqrt(length(cc_store));        % Standard Error
        ts = tinv([0.025  0.975],length(cc_store)-1);      % T-Score
        cc_cis{p1} = mean(cc_store) + ts*SEM;              % Confidence Intervals
        
        %%%%%%%%%%%
        % Now calculate the distribution of invaded patches for each
        % resample
        rescaled_values_store = zeros(num_resamples,37);
        for p2 = 1:1:num_resamples
            fprintf('p1 = %d, part B: p2 = %d\n', p1, p2);
            [values, edges] = histcounts(hold_resamples{p1}{p2}(:,7),[-0.5 : 1 : 36.5]); % specify that all cases use the same set of all bins
            rescaled_values = 1/sum(values)*values;
            rescaled_values_store(p2,:) = rescaled_values;
        end
        distribution_centers = (edges(1:end-1)+edges(2:end))/2;  

        % Get the mean and CI of each datapoint in the distribution
        distribution_data = zeros(37,3);
        for num_patches = 1:1:37
            this_patch_data = rescaled_values_store(:,num_patches);
            distribution_data(num_patches,2) = mean(this_patch_data);
            SEM = std(this_patch_data)/sqrt(length(this_patch_data));           % Standard Error
            ts = tinv([0.025  0.975],length(this_patch_data)-1);                % T-Score
            patch_confidence_intervals = mean(this_patch_data) + ts*SEM;        % Confidence Intervals
            distribution_data(num_patches,1) = patch_confidence_intervals(1);
            distribution_data(num_patches,3) = patch_confidence_intervals(2);
        end
        distribution_stats{p1} = distribution_data;
    
    end
    writematrix(cc_mean,'SSI Species Spatial Invasion/cc_mean.csv');
    writematrix(cc_cis{1},'SSI Species Spatial Invasion/cc_cis_1.csv');
    writematrix(cc_cis{2},'SSI Species Spatial Invasion/cc_cis_2.csv');
end




%% Distributions:
%
% Distribution of:
% - the number of patches invaded

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

is_distributions = true;
if (is_distributions)
    
    mega_data = {total_data_frac, total_test5_data, total_test4_data, total_test6_data};
    r_squared = zeros(4,1);
    
    medians = zeros(2,3);
    for c=1:1:2
        
        hFig=figure('visible','on');
        clf(hFig);
        hold off;
        
        for c0 = 1:1:2
            d=2*(c-1)+c0;
            
                   
            clearvars rescaled_values values compare_norm

            % Number of patches invaded:
            patch_invasion_stats = zeros(4,1);
            patch_invasion_stats(1) = min(mega_data{d}(:,7));
            patch_invasion_stats(2) = mean(mega_data{d}(:,7));
            patch_invasion_stats(3) = max(mega_data{d}(:,7));
            patch_invasion_stats(4) = std(mega_data{d}(:,7));

            % Plot the distribution:
            [values, edges] = histcounts(mega_data{d}(:,7));
            centers = (edges(1:end-1)+edges(2:end))/2;

            compare_norm = normpdf([patch_invasion_stats(1):1:patch_invasion_stats(3)],patch_invasion_stats(2),patch_invasion_stats(4));

            % Determine R^2:
            rescaled_values = 1/sum(values)*values;
            ss_res = 0.0;
            ss_tot = 0.0;
            for x=1:1:size(rescaled_values,2)
                ss_res = ss_res+(rescaled_values(x)-compare_norm(x))^2.0;
                ss_tot = ss_tot + (rescaled_values(x)-mean(rescaled_values))^2.0;
            end
            r_squared(d) = 1-ss_res/ss_tot;


            if (d<2)
                plot([patch_invasion_stats(1):0.1:patch_invasion_stats(3)], normpdf([patch_invasion_stats(1):0.1:patch_invasion_stats(3)],patch_invasion_stats(2),patch_invasion_stats(4)),'LineWidth',2.5,'Color',color_array{c0});
            end
            hold on;
            scatter(centers, rescaled_values,30,'o','filled','MarkerFaceColor',color_array{c0},'MarkerEdgeColor','k','LineWidth',1.2);
            
            % locate_median
            running_cdf = 0.0;
            for x=1:1:size(rescaled_values,2)
                running_cdf = running_cdf + rescaled_values(x);
                if (running_cdf >= 0.5)
                    medians(c,c0) = centers(x);
                    break;
                end
            end
                 
        
        end
        
        % Now plot the resampled distribution data
        hold on;
        errorbar(distribution_centers,distribution_stats{c}(:,2),abs(distribution_stats{c}(:,2)-distribution_stats{c}(:,1)),abs(distribution_stats{c}(:,2)-distribution_stats{c}(:,3)),'.','Color',color_array{4},'HandleVisibility','off');
        scatter(distribution_centers, distribution_stats{c}(:,2),30,'o','filled','MarkerFaceColor',color_array{4},'MarkerEdgeColor','k','LineWidth',1.2);
        sum(distribution_stats{c}(:,2))
        
        running_cdf = 0.0;
        for x=1:1:size(distribution_stats{c}(:,2),1)
            running_cdf = running_cdf + distribution_stats{c}(x,2);
            if (running_cdf >= 0.5)
                medians(c,3) = distribution_centers(x);
                break;
            end
        end
        %
        
        box on;
        xlabel('Number of patches invaded');
        ylabel('Probability');
       
        
        str1 = sprintf('Fitted normal distribution PDF: %2.3f',r_squared(1));
        str2 = sprintf('Norm. data: resident re-invasion. Median: %2.0f', medians(1,1));
        str3 = sprintf('Norm. data: random invasion. Median: %2.0f', medians(1,2));
        str4 = sprintf('Norm. data: empty resident re-invasion. Median: %2.0f', medians(2,1));
        str5 = sprintf('Norm. data: empty random invasion. Median: %2.0f', medians(2,2));
        str6 = sprintf('Norm. data: re-sampled random invasion. Median: %2.0f', medians(1,3));
        str7 = sprintf('Norm. data: re-sampled empty random invasion. Median: %2.0f', medians(2,3));
        
        xlim([0 37]);
        if (c==1)
            legend(str1,str2,str3,str6,'Location','NorthWest')
            ylim([0 0.11]);
        else
            legend(str4,str5,str7,'Location','NorthWest')
            ylim([0 0.35]);
        end
        printname=sprintf('SSI Species Spatial Invasion/FIGURES/V2/spatial_invasion_num_patches_distribution_%d',c);
        print(printname,'-dpng','-r400','-painters');
        close(hFig);
        clearvars x temp_ss_res ss_tot rescaled_values;
    end

    clearvars para edges values centers hFig;
end

