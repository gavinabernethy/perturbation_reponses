%% Spatial Invasion (SSI) analysis
% Code to analyse the spatial_invasion output file from this subroutine.


clear;
set(0,'DefaultAxesFontSize',15);
print_images=true; % only print the final special correlation

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

% 1. k3;   2. n_spec_store;   3. g1 (id of invading species);   4. initial_range(g1);   5. n_spec;
% 6. secondary_extinctions;   7. num_cells_invaded;   8. initial_biomass(g1);   9. initial_population(g1);
% 10. trophic_levels_average(g1,1);   11. trophic_levels_average(g1,2);   12. trophic_levels_average(g1,3);
% 13. bodysize_store(g1,1).  14. bodysize_store(g1,3).  15. bodysize_store(g1,4)
% -
% 16. current_n_spec_global;   17. current_max_sc_tl_global;   18. current_link_density_global;
% 19. current_connectance_global;   20. current_ave_sc_tl_global;   21. current_omnivory_frac_global;
% 22. current_ave_population_global;   23. current_ave_bodysize_global;   24. current_ave_distribution_global;

control_prop_num=[4,8,9,10,11,12,13]; % 7 inputs
control_prop_name={'Range','Biomass','Population','SCTL (patch average)','SCTL (population average)',...
        'SCTL (biomass average)','Bodysize'};
% Note that biomass and population properties are totalled across all
% instances of the species in different patches

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
    
    % calculate the fraction of "secondary" extinctions due to the
    % invasion, assuming that everything that was driven extinct during the
    % initial deletion of the species was accidentally allowed to be
    % resurrected when the species was re-introduced, so add the
    % num_secondary_extinctions from SSS(2:end,6) to the absolute loss of
    % species as well as 1 to account for the species itself being
    % re-introduced for the invasion.
    individual_data_frac{s2}(:,response_prop_num(1)) = (1+individual_data_SSS{s2}(2:end,6)+individual_data_SSS{s2}(2:end,response_prop_num_init(1)) - individual_data{s2}(:,response_prop_num(1)))./(individual_data_SSS{s2}(2:end,response_prop_num_init(1))-36);
    
    % The first three terms here are equivalent to just the initial global
    % diversity before the initial deletion of the species, becuase we are
    % essentially assuming that all extinctions in this process (due to the
    % mistake) are only those that occur due to *invasion* specificially
    % (as the cells that the species was originally present in won't be
    % able to change anyway due to the low-population invasion and that the
    % timeframe does not allow any other species to evolve, so there should
    % not really be significant change in the ecosystem from removing and
    % re-introduing the species and any species whose populations dipped
    % under the threshold). Therefore, we are simply asking how many such
    % extinctions occuredd, AS A FRACTION OF THE SPECIES WHO COULD SURVIVE
    % IN THE ABSENCE OF THIS SPECIES (i.e. dividing by the result of the correct deletion experiment rather than that initial value).
    
    % calculate the fractional change in other existing metacommunity
    % properties as a result of the invasion
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

%% Analysis - Correlations:

collect_correlations = false;
if (collect_correlations)
    % Begin by collecting all 7x11 +  1 (sec. extinctions vs. number of invaded patches) correlations individually for each of
    % the 15 simulations and storing them in a table
    row_list = {'Sim 1','Sim 2','Sim 3','Sim 4','Sim 5','Sim 6','Sim 7','Sim 8','Sim 9','Sim 10','Sim 11','Sim 12','Sim 13','Sim 14','Sim 15','Combined'};

    column_list={
        'Range vs. Secondary extinctions response','Range vs. Local diversity response','Range vs. Max SCTL response','Range vs. Link density response','Range vs. Connectance response',...
        'Range vs. Average SCTL response','Range vs. Omnivory response','Range vs. Population response','Range vs. Bodysize response','Range vs. Range response','Range vs. Number of Invasions response',...
        'Biomass vs. Secondary extinctions response','Biomass vs. Local diversity response','Biomass vs. Max SCTL response','Biomass vs. Link density response','Biomass vs. Connectance response',...
        'Biomass vs. Average SCTL response','Biomass vs. Omnivory response','Biomass vs. Population response','Biomass vs. Bodysize response','Biomass vs. Range response','Biomass vs. Number of Invasions response',...  
        'Population vs. Secondary extinctions response','Population vs. Local diversity response','Population vs. Max SCTL response','Population vs. Link density response',...
        'Population vs. Connectance response','Population vs. Average SCTL response','Population vs. Omnivory response','Population vs. Population response','Population vs. Bodysize response','Population vs. Range response','Population vs. Number of Invasions response',...
        'SCTL (patch average) vs. Secondary extinctions response','SCTL (patch average) vs. Local diversity response','SCTL (patch average) vs. Max SCTL response','SCTL (patch average) vs. Link density response','SCTL (patch average) vs. Connectance response',...
        'SCTL (patch average) vs. Average SCTL response','SCTL (patch average) vs. Omnivory response','SCTL (patch average) vs. Population response','SCTL (patch average) vs. Bodysize response','SCTL (patch average) vs. Range response','SCTL (Patch average) vs. Number of Invasions response',...
        'SCTL (population average) vs. Secondary extinctions response','SCTL (population average) vs. Local diversity response','SCTL (population average) vs. Max SCTL response','SCTL (population average) vs. Link density response','SCTL (population average) vs. Connectance response',...
        'SCTL (population average) vs. Average SCTL response','SCTL (population average) vs. Omnivory response','SCTL (population average) vs. Population response','SCTL (population average) vs. Bodysize response','SCTL (population average) vs. Range response','SCTL (Population average) vs. Number of Invasions response',...
        'SCTL (biomass average) vs. Secondary extinctions response','SCTL (biomass average) vs. Local diversity response','SCTL (biomass average) vs. Max SCTL response','SCTL (biomass average) vs. Link density response','SCTL (biomass average) vs. Connectance response',...
        'SCTL (biomass average) vs. Average SCTL response','SCTL (biomass average) vs. Omnivory response','SCTL (biomass average) vs. Population response','SCTL (biomass average) vs. Bodysize response','SCTL (biomass average) vs. Range response','SCTL (biomass average) vs. Number of Invasions response',...
        'Bodysize vs. Secondary extinctions response','Bodysize vs. Local diversity response','Bodysize vs. Max SCTL response','Bodysize vs. Link density response','Bodysize vs. Connectance response',...
        'Bodysize vs. Average SCTL response','Bodysize vs. Omnivory response','Bodysize vs. Population response','Bodysize vs. Bodysize response','Bodysize vs. Range response','Bodysize vs. Number of Invasions response',...
        'Secondary extinctions response vs. Number of Invasions response'};

        
    type_list=[repmat("double",1,78)];
    correlation_table=table('Size',[16 78],'VariableNames',column_list,'RowNames',row_list,'VariableTypes',type_list);

    % Mostly populate it with correlations of properties within the individual
    % simulations
    % control_prop_xlim={[0 7],[-4000 90000],[-2000 25000],[0 4.5],[0 4.5],[0 4.5],[0 15]};

    % cycle simulations 1-15
    for j0=1:1:15

        % cycle control properties x
        for j1=1:1:7

            % cycle response properties y
            for j2=1:1:11
                corr_matrix=corrcoef(individual_data_frac{j0}(:,control_prop_num(j1)),individual_data_frac{j0}(:,response_prop_num(j2)));
                % Calculate the linear correlation coefficient
                correlation_table{j0,(j1-1)*11+j2}=corr_matrix(2,1);
            end
        end
        
        % final column (sec. extinctions vs. number of invaded patches)
        corr_matrix=corrcoef(individual_data_frac{j0}(:,response_prop_num(1)),individual_data_frac{j0}(:,response_prop_num(11)));
        correlation_table{j0,78}=corr_matrix(2,1);
        
    end
     
    % Correlate modified responses against controls for the whole dataset:
    corr_coef_values=zeros(7,11);
    
    control_prop_xlim = {[0 7], [-5000 100000], [-2000 25000], [0 5], [0 5], [0 5], [0 15]};
    
    response_frac_prop_ylim = {[-0.01 0.15], [-0.1 0.08], [-0.1 0.12], [-0.08 0.1], [-0.08 0.1], [-0.08 0.1],...
        [-0.2 0.2], [-0.2 0.1], [-0.05 0.15], [-0.1 0.8], [0 40]};

    for t1=1:1:7
        for t2=1:1:11

            corr_matrix=corrcoef(total_data_frac(:,control_prop_num(t1)),total_data_frac(:,response_prop_num(t2)));
            corr_coef_values(t1,t2)=corr_matrix(2,1); % Calculate the linear correlation coefficient
            
            % save in the final row of the table
            correlation_table{16,(t1-1)*11+t2}=corr_matrix(2,1);

            if (print_images==true)
                hFig=figure('visible','on');
                clf(hFig);
                hold off;
                scatter(total_data_frac(:,control_prop_num(t1)),total_data_frac(:,response_prop_num(t2)),25,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',0.8);
                box on;
                xlabel(control_prop_name{t1});
                ylabel(response_prop_name{t2});
                xlim(control_prop_xlim{t1});
                ylim(response_frac_prop_ylim{t2});
                legend(num2str(corr_coef_values(t1,t2),3),'Location','NorthEast');
                printname=sprintf('SSI Species Spatial Invasion/FIGURES/spatial_invasion_correlation_%d_%d',t1,t2);
                print(printname,'-dpng','-r400','-painters');
                close(hFig);
            end

        end
    end
    
    
    % final element (sec. extinctions vs. number of invaded patches, both of which are RESPONSES)
    corr_matrix=corrcoef(total_data_frac(:,response_prop_num(1)),total_data_frac(:,response_prop_num(11)));
    correlation_table{16,78}=corr_matrix(2,1);
    hFig=figure('visible','on');
    clf(hFig);
    hold off;
    scatter(total_data_frac(:,response_prop_num(11)),total_data_frac(:,response_prop_num(1)),25,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',0.8);
    box on;
    xlabel(response_prop_name{11});
    ylabel(response_prop_name{1});
    xlim(response_frac_prop_ylim{11});
    ylim(response_frac_prop_ylim{1});
    legend(num2str(corr_matrix(2,1),3),'Location','NorthWest');
    printname=sprintf('SSI Species Spatial Invasion/FIGURES/spatial_invasion_correlation_SPECIAL_1_11');
    print(printname,'-dpng','-r400','-painters');
    close(hFig);

    %% Format the Correlations Table
    writetable(correlation_table,'SSI Species Spatial Invasion/spatial_invasion_correlation_table_RAW.xlsx','WriteRowNames',true);

end

%% Analysis - Multiple linear regression:

%%%%%%%%
% Multiple regression of each of the 11 dependent ("response") variables
% in turn against a regression of 5 of the dependent ("control") variables
% (as always, only include one measure of SCTL).
% None of the correlation images clearly indicated a specific nonlinear
% relationship, but this will need to be re-checked when the correct data
% is run for the number of secondary extinctions - check its correlation
% images!

is_regress=false;
if (is_regress==true)
    
    mult_regress_res=cell(11,2,2,2,2,2);
    
    best_rsquareadjusted=zeros(11,1);
    best_zlist=cell(11,1);
    best_zlist(cellfun('isempty',best_zlist))={[0,0,0,0,0,0]};
    list_of_indepvar = {'Range','Biomass','Population','SCTL (patch average)','Bodysize'};
    
    
    % Combine all the independent variables into a single predictor matrix X
    X = zeros(7022,5);
    ind_prop_regress_rows = [4,8,9,10,13];
    
    for spec=1:1:7022
        for ind_prop=1:1:5
            X(spec,ind_prop) = total_data_frac(spec,ind_prop_regress_rows(ind_prop));
        end
    end
    

    % Build a series of models for each dependent variable:
    for dep_prop=1:1:11
        
        y = total_data_frac(:,response_prop_num(dep_prop));

                        
        % Now try including each of the variables
        clearvars X3;
        for z1=0:1:1
             for z2=0:1:1
                 for z3=0:1:1
                     for z4=0:1:1
                         for z5=0:1:1

                             z_list=[z1,z2,z3,z4,z5];
                             if (sum(z_list)>0)
                                 [X3,name_list] = x3_constructor(z_list,X,list_of_indepvar,response_prop_name{dep_prop});
                                 mult_regress_res{dep_prop,z1+1,z2+1,z3+1,z4+1,z5+1}=fitlm(X3,y,'VarNames',name_list);

                                 % check if this is the
                                 % best model so far?

                                 rsq_adj = mult_regress_res{dep_prop,z1+1,z2+1,z3+1,z4+1,z5+1}.Rsquared.Adjusted;
                                 if (rsq_adj > best_rsquareadjusted(dep_prop))
                                     best_rsquareadjusted(dep_prop) = rsq_adj;
                                     best_zlist{dep_prop} = z_list;
                                 end
                             end
                         end
                     end
                 end
             end
         end     
    end
end


%% Distributions:
%
% Distribution of:
% - the number of patches invaded
% - absolute number of secondary extinctions
% - fraction of secondary extinctions
%
% Also return the mean, maximum and minimum of these properties.
% Then fit normal, power and exponential functions to them?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

is_distributions = true;
if (is_distributions)

    % Number of patches invaded:
    patch_invasion_stats = zeros(4,1);
    patch_invasion_stats(1) = min(total_data_frac(:,7));
    patch_invasion_stats(2) = mean(total_data_frac(:,7));
    patch_invasion_stats(3) = max(total_data_frac(:,7));
    patch_invasion_stats(4) = std(total_data_frac(:,7));

    % Plot the distribution:
    [values, edges] = histcounts(total_data_frac(:,7));
    centers = (edges(1:end-1)+edges(2:end))/2;

    compare_norm = normpdf([1:1:36],patch_invasion_stats(2),patch_invasion_stats(4));

    % identify the best re-scaling parameter:
    best_ss_res = 999999999999.0;
    best_para = 0;
    for para = 0.0:0.000000001:0.01
        rescaled_values = para*values;
        temp_ss_res = 0.0;
        for x=1:1:36
            temp_ss_res = temp_ss_res+(rescaled_values(x)-compare_norm(x))^2.0;
        end
        if (temp_ss_res<best_ss_res)
            best_ss_res = temp_ss_res;
            best_para = para;
        end
    end

    % Determine R^2:
    ss_tot = 0.0;
    rescaled_values = best_para*values;
    for x=1:1:36
        ss_tot = ss_tot + (rescaled_values(x)-mean(rescaled_values))^2.0;
    end
    r_squared = 1-best_ss_res/ss_tot; % = 0.9853

    hFig=figure('visible','on');
    clf(hFig);
    hold off;
    plot([patch_invasion_stats(1):0.1:patch_invasion_stats(3)], normpdf([patch_invasion_stats(1):0.1:patch_invasion_stats(3)],patch_invasion_stats(2),patch_invasion_stats(4)),'LineWidth',2.5);
    hold on;
    scatter(centers, best_para*values,30,'o','filled','MarkerFaceColor','red','MarkerEdgeColor','k','LineWidth',1.2);
    box on;
    xlabel('Number of patches invaded');
    ylabel('Probability density');
    legend('Normal distribution PDF','Normalised data','Location','NorthEast')
    printname=sprintf('SSI Species Spatial Invasion/FIGURES/spatial_invasion_num_patches_distribution');
    print(printname,'-dpng','-r400','-painters');
    close(hFig);

    clearvars x temp_ss_res ss_tot rescaled_values para edges values centers hFig;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Number of secondary extinctions:
    absolute_sec_ext = zeros(7022,1);
    relative_sec_ext = zeros(7022,1);
    start_line=1;
    end_line=0;
    for s4=1:1:15 % simulation
        end_line=end_line+lengths_store(s4);

        absolute_sec_ext(start_line:end_line) = (1+individual_data_SSS{s4}(2:end,6) + individual_data_SSS{s4}(2:end,response_prop_num_init(1)) - total_data(start_line:end_line,response_prop_num(1)));
        relative_sec_ext(start_line:end_line) = (1+individual_data_SSS{s4}(2:end,6) + individual_data_SSS{s4}(2:end,response_prop_num_init(1)) - total_data(start_line:end_line,response_prop_num(1)))./(individual_data_SSS{s4}(2:end,response_prop_num_init(1))-36);

        start_line=start_line+lengths_store(s4);
    end
    clearvars s4 start_line end_line;

    % For both absolute and relative secondary extinctions
    best_para_a = zeros(2,2);
    best_para_n = zeros(2,2);
    type_str = {'Absolute','Relative'};
    for type = 1:1:2
        for law = 1:1:2

            if (type==1)
                data = absolute_sec_ext;
                x_start = 7.5; % from the center of the largest bin
                x_end = 49.5; % to the center of the last bin before height non-zero
                para_n_start = 0;
                para_n_end = 4;
                para_n_inc = 0.0001;
                if (law==1)
                    para_a_start=800000;
                    para_a_end=2500000;
                    para_a_inc=100;
                else
                    para_a_start=2300;
                    para_a_end=2500;
                    para_a_inc=1;
                end
            else
                para_n_start = 50;
                para_n_end = 100;
                para_n_inc = 0.0001;
                para_a_start=2000;
                para_a_end=5000;
                para_a_inc=1;
                if (law==1)
                    % need to shift
                    data = relative_sec_ext+1;
                    x_start = 1.0155;
                    x_end = 1.1115;
                else
                    data = relative_sec_ext;
                    x_start = 0.0135;
                    x_end = 0.1185;
                end
            end

            [values, edges] = histcounts(data);
            centers = (edges(1:end-1)+edges(2:end))/2;
            % x_start and x_end should be values in "centers"

            ind_start = find(abs(centers-x_start)<=0.000001);
            ind_end = find(abs(centers-x_end)<=0.0001);

            if (law==1)
                % fit power law y = ax^-n => log(y) = -nlog(x) + log(a)
                x_data = log(centers(ind_start:ind_end));
                y_data = log(values(ind_start:ind_end));
                if (type==1)
                    xlabel_str = sprintf('log(%s secondary extinctions)',type_str{type});
                else
                    xlabel_str = sprintf('log(%s secondary extinctions + 1)',type_str{type});
                end
            else
                % fit exponential law y = ae^-nx => log(y) = - nx + log(a)
                x_data = centers(ind_start:ind_end);
                y_data = log(values(ind_start:ind_end));
                xlabel_str = sprintf('%s secondary extinctions',type_str{type});
            end

            % identify the best parameters a and n:

            best_ss_res = 999999999999.0;

            for para_a = para_a_start : para_a_inc : para_a_end
                for para_n = para_n_start : para_n_inc : para_n_end

                    temp_ss_res = 0.0;
                    for x=1:1:size(x_data,2)
                        res = y_data(x) - (-para_n*x_data(x) + log(para_a));
                        temp_ss_res = temp_ss_res + res^2.0;
                    end
                    if (temp_ss_res<best_ss_res)
                        best_ss_res = temp_ss_res;
                        best_para_a(type,law) = para_a;
                        best_para_n(type,law) = para_n;
                    end

                end
            end

            % Determine R^2:
            ss_tot = 0.0;
            for x=1:1:size(x_data,2)
                ss_tot = ss_tot + (y_data(x)-mean(y_data))^2.0;
            end
            r_squared(type,law) = 1-best_ss_res/ss_tot;

            % Print
            hFig=figure('visible','on');
            clf(hFig);
            hold off;
            plot(x_data,-best_para_n(type,law)*x_data+log(best_para_a(type,law)),'LineWidth',2.5,'Color','red');
            hold on;
            scatter(x_data, y_data,30,'o','filled','MarkerFaceColor','blue','MarkerEdgeColor','k','LineWidth',1.2);
            box on;
            xlabel(xlabel_str);
            ylabel('log(Frequency)');
            r_str = num2str(r_squared(type,law),3);
            if (type==1)
                if (law==1)
                    xlim([1.8 4.2]);
                else
                    xlim([0 55]);
                end
            else
                if (law==1)
                    xlim([0 0.1]);
                else
                    xlim([0 0.12]);
                end
            end
            ylim([-0.5 8]);
            legend_str = sprintf('Fitted model: R^2 = %s',r_str);
            legend('Data',legend_str,'Location','SouthWest')
            printname=sprintf('SSI Species Spatial Invasion/FIGURES/spatial_invasion_sec_ext_dist_%d_%d',type,law);
            print(printname,'-dpng','-r400','-painters');
            close(hFig);
        end

        % Then create the original plot    
        if (type==1)
            data = absolute_sec_ext;
            x_start = 7.5;
            x_end = 49.5;
            x_inc = 0.0001;
        else
            data = relative_sec_ext;
            x_start = 0.015;
            x_end = 0.10;
            x_inc = 0.0000001;
        end
        x_vals = [x_start:x_inc:x_end];
        [values, edges] = histcounts(data);
        centers = (edges(1:end-1)+edges(2:end))/2;
        hFig=figure('visible','on');
        clf(hFig);
        hold off; % power law
        if (type==1)
            plot(x_vals,best_para_a(type,1)*x_vals.^-best_para_n(type,1),'LineWidth',2.5,'Color',color_array{1});
        else
            plot(x_vals,best_para_a(type,1)*(x_vals+1).^-best_para_n(type,1),'LineWidth',2.5,'Color',color_array{1});
        end
        hold on; % exponential law
        plot(x_vals,best_para_a(type,2)*exp(-x_vals*best_para_n(type,2)),'LineWidth',2.5,'Color',color_array{4});
        hold on; % data
        scatter(centers,values,30,'o','filled','MarkerFaceColor',color_array{2},'MarkerEdgeColor','k','LineWidth',1.2);
        box on;
        xlabel_str = sprintf('%s secondary extinctions',type_str{type});
        xlabel(xlabel_str);
        ylabel('Frequency');
        a_str = num2str(best_para_a(type,1),7);
        n_str = num2str(best_para_n(type,1),3);
        if (type==1)
            ylim([0 600]);
            power_str = sprintf('Power law: N(s) = %s(s)^{-%s}',a_str,n_str);
        else
            ylim([0 900]);
            power_str = sprintf('Power law: N(s) = %s(s+1)^{-%s}',a_str,n_str);
        end
        a_str = num2str(best_para_a(type,2),6);
        n_str = num2str(best_para_n(type,2),3);
        exp_str = sprintf('Exponential law: N(s) = %sexp(-%ss)',a_str,n_str);
        legend(power_str,exp_str,'Data','Location','NorthEast')
        printname=sprintf('SSI Species Spatial Invasion/FIGURES/spatial_invasion_sec_ext_dist_full_%d',type);
        print(printname,'-dpng','-r400','-painters');
        close(hFig);

    end
end

%% Function: 

function [Xmatrix, var_names] = x3_constructor(z_list,allX,master_names,res_name)
    num_cols = sum(z_list);
    Xmatrix = zeros(7022,num_cols);
    current_col = 0;
    var_names = cell(num_cols+1,1);
    for z = 1:1:5
        if (z_list(z)==1)
            current_col=current_col+1;
            Xmatrix(:,current_col) = allX(:,z);
            var_names{current_col} = master_names{z};
        end
    end
    var_names{num_cols+1} = res_name;
end







