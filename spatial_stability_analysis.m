%% Spatial Stability analysis
% Code to analyse the spatial_stability output file from this subroutine.

% Consider if new or stronger correlations (e.g. between bodysize and
% range) emerge when the different simulations are considered in isolation.
% Or, if things like 'range' are judged RELATIVE to the individual
% simulation. It might be better to consider these as 'small', 'large'
% etc. as it may be that larger bodysize corresponds to large range or
% SCTL, but the absolute values/scales of these quantities may vary across
% simulations and so the patterns could be obscured when the data is all
% mixed together.


clear;
set(0,'DefaultAxesFontSize',15);
print_images=true;
color_array={ [0.05, 0.65, 0.25], [247/256, 182/256, 2/256], [0.15, 0.8, 0.4], [61/256, 89/256, 252/256], [252/256, 88/256, 76/256] };



%% Collect data:
% Loop over each of the 15 repeat simulations

count = 0;
individual_data=cell(15,1);
total_data=zeros(0,23);
lengths_store=zeros(15,1);

for s1 = 100:1:120
    
    filename=sprintf('SSS Spatial Stability data/DATA/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_spatial_stability_00010000.txt',s1);
    
    if (exist(filename,'file')==2)
        
        count = count+1;
        
        % store the individual data in a cell array
        individual_placeholder=readmatrix(filename);
        individual_data{count}=individual_placeholder;
        
        % concatenate them (without the 'zero' line in each case) into a single vector
        current_length = size(total_data,1);
        add_length = size(individual_data{count},1)-1;
        lengths_store(count)=size(individual_data{count},1);
        new_length = current_length+add_length;
        total_data_placeholder = zeros(new_length,23);
        total_data_placeholder(1:current_length,:)=total_data;
        total_data_placeholder(current_length+1:new_length,:)=individual_placeholder(2:end,:);
        clearvars total_data;
        total_data = total_data_placeholder;
        clearvars total_data_placeholder;
        
       
        % 23 x N, where N is the number of non-resource species + 1 (first
        % row is general data on the global ecosystem before any perturbation 
        % - e.g. number of species, average local connectance, etc.)
        %
        % Columns:
        % 1 - timestep (always 10,000)
        % 2 - total number of species before deletion (will be constant)
        % 3 - species deleted in this row (0, 37, 38, 39, ...)
        % 4 - range of deleted species
        % 5 - remaining number of species
        % 6 - number of secondary extinctions
        % 7 - biomass of deleted species 
        % 8 - population of deleted species
        % 9 - average trophic level of deleted species by cell
        % 10 - average trophic level of deleted species by population
        % 11 - average trophic level of deleted species by biomass
        % 12 - bodysize of deleted species
        % 13 - fraction of non-directional links that can be traversed by
        % the deleted species [AS WE ARE NOT USING TRAIT-GATED TRAVERSAL
        % FOR THIS EXPERIMENT, THIS WILL ALWAYS BE 1]
        % 14 - average number of matched traits across all non-directional
        % links for the deleted species [AS WE ARE NOT USING TRAIT-GATED 
        % TRAVERSAL FOR THIS EXPERIMENT, THIS WILL ALWAYS BE 1]
        % 15 - average local diversity
        % 16 - average local max SCTL
        % 17 - average local link density 
        % 18 - average local connectance
        % 19 - average local average SCTL
        % 20 - average local omnivory fraction
        % 21 - average local average population
        % 22 - average local average bodysize
        % 23 - average local average range of local species
        
    end
end

clearvars -except individual_data total_data lengths_store print_images color_array;

%% Modification:

response_prop_num=[6,15,16,17,18,19,20,21,22,23];
response_prop_name={'Relative secondary extinctions','Local diversity Response','Max SCTL Response','Link density Response',...
    'Connectance Response','Average SCTL Response','Omnivory Response','Population Response','Bodysize Response','Range Response'};
response_frac_prop_ylim={[-0.001 0.022],[-0.025 0.025],[-0.05 0.05],[-0.025 0.025],[-0.035 0.035],[-0.025 0.025],[-0.07 0.07],[-0.07 0.07],[-0.02 0.02],[-0.04 0.04]};

% for response properties, need the "fractional change" rather than the raw
% change or the resultant value:
individual_data_frac = individual_data;
total_data_frac = total_data;

% Modify the individual data:
for s2=1:1:15
    for s3=1:1:10
        if (s3==1)
            % fractional secondary extinctions  = num secondary extinctions/(original num species -1-36)
            individual_data_frac{s2}(:,response_prop_num(s3)) = individual_data{s2}(:,response_prop_num(s3))./(individual_data{s2}(:,2)-1-36);
        else
            % fractional change in other properties = (new property - original property)/original property
            individual_data_frac{s2}(:,response_prop_num(s3)) = (individual_data{s2}(:,response_prop_num(s3)) - individual_data{s2}(1,response_prop_num(s3)))./individual_data{s2}(1,response_prop_num(s3));
        end
    end
end

% Modify the concatenated data:
start_line=1;
end_line=0;
for s4=1:1:15
    end_line=end_line+lengths_store(s4)-1;
        
    for s5=1:1:10
        if (s5==1)
            % fractional secondary extinctions  = num secondary extinctions/(original num species -1 -36)
            total_data_frac(start_line:end_line,response_prop_num(s5)) = total_data(start_line:end_line,response_prop_num(s5))./(total_data(start_line:end_line,2)-1-36);
        else
            % fractional change in other properties = (new property - original property)/original property
            total_data_frac(start_line:end_line,response_prop_num(s5)) = (total_data(start_line:end_line,response_prop_num(s5)) - individual_data{s4}(1,response_prop_num(s5)))./individual_data{s4}(1,response_prop_num(s5));
        end
    end
                 
    start_line=start_line+lengths_store(s4)-1;
end

clearvars start_line end_line s1 s2 s3 s4 s5 r1 r2 r3 r4 r5

%% Analysis: Correlation

control_prop_num=[4,7,8,9,10,11,12];
control_prop_name={'Range','Biomass','Population','SCTL (patch average)','SCTL (population average)',...
        'SCTL (biomass average)','Bodysize'};

collect_correlations = false;
if (collect_correlations)
    % Begin by collecting all 7x(10+7) correlations individually for each of
    % the 15 simulations and storing them in a table
    row_list = {'Sim 1','Sim 2','Sim 3','Sim 4','Sim 5','Sim 6','Sim 7','Sim 8','Sim 9','Sim 10','Sim 11','Sim 12','Sim 13','Sim 14','Sim 15','Combined'};

    column_list={'Range vs. Range','Range vs. Biomass','Range vs. Population','Range vs. SCTL (patch average)','Range vs. SCTL (population average)',...
        'Range vs. SCTL (biomass average)','Range vs. Bodysize',...
        'Range vs. Secondary extinctions response','Range vs. Local diversity response','Range vs. Max SCTL response','Range vs. Link density response',...
        'Range vs. Connectance response','Range vs. Average SCTL response','Range vs. Omnivory response','Range vs. Population response','Range vs. Bodysize response','Range vs. Range response',...
        'Biomass vs. Range','Biomass vs. Biomass','Biomass vs. Population','Biomass vs. SCTL (patch average)','Biomass vs. SCTL (population average)',...
        'Biomass vs. SCTL (biomass average)','Biomass vs. Bodysize',...
        'Biomass vs. Secondary extinctions response','Biomass vs. Local diversity response','Biomass vs. Max SCTL response','Biomass vs. Link density response',...
        'Biomass vs. Connectance response','Biomass vs. Average SCTL response','Biomass vs. Omnivory response','Biomass vs. Population response','Biomass vs. Bodysize response','Biomass vs. Range response',...
        'Population vs. Range','Population vs. Biomass','Population vs. Population','Population vs. SCTL (patch average)','Population vs. SCTL (population average)',...
        'Population vs. SCTL (biomass average)','Population vs. Bodysize',...
        'Population vs. Secondary extinctions response','Population vs. Local diversity response','Population vs. Max SCTL response','Population vs. Link density response',...
        'Population vs. Connectance response','Population vs. Average SCTL response','Population vs. Omnivory response','Population vs. Population response','Population vs. Bodysize response','Population vs. Range response',...
        'SCTL (patch average) vs. Range','SCTL (patch average) vs. Biomass','SCTL (patch average) vs. Population','SCTL (patch average) vs. SCTL (patch average)','SCTL (patch average) vs. SCTL (population average)',...
        'SCTL (patch average) vs. SCTL (biomass average)','SCTL (patch average) vs. Bodysize',...
        'SCTL (patch average) vs. Secondary extinctions response','SCTL (patch average) vs. Local diversity response','SCTL (patch average) vs. Max SCTL response','SCTL (patch average) vs. Link density response',...
        'SCTL (patch average) vs. Connectance response','SCTL (patch average) vs. Average SCTL response','SCTL (patch average) vs. Omnivory response','SCTL (patch average) vs. Population response','SCTL (patch average) vs. Bodysize response','SCTL (patch average) vs. Range response',...
        'SCTL (population average) vs. Range','SCTL (population average) vs. Biomass','SCTL (population average) vs. Population','SCTL (population average) vs. SCTL (patch average)','SCTL (population average) vs. SCTL (population average)',...
        'SCTL (population average) vs. SCTL (biomass average)','SCTL (population average) vs. Bodysize',...
        'SCTL (population average) vs. Secondary extinctions response','SCTL (population average) vs. Local diversity response','SCTL (population average) vs. Max SCTL response','SCTL (population average) vs. Link density response',...
        'SCTL (population average) vs. Connectance response','SCTL (population average) vs. Average SCTL response','SCTL (population average) vs. Omnivory response','SCTL (population average) vs. Population response','SCTL (population average) vs. Bodysize response','SCTL (population average) vs. Range response',...
        'SCTL (biomass average) vs. Range','SCTL (biomass average) vs. Biomass','SCTL (biomass average) vs. Population','SCTL (biomass average) vs. SCTL (patch average)','SCTL (biomass average) vs. SCTL (population average)',...
        'SCTL (biomass average) vs. SCTL (biomass average)','SCTL (biomass average) vs. Bodysize',...
        'SCTL (biomass average) vs. Secondary extinctions response','SCTL (biomass average) vs. Local diversity response','SCTL (biomass average) vs. Max SCTL response','SCTL (biomass average) vs. Link density response',...
        'SCTL (biomass average) vs. Connectance response','SCTL (biomass average) vs. Average SCTL response','SCTL (biomass average) vs. Omnivory response','SCTL (biomass average) vs. Population response','SCTL (biomass average) vs. Bodysize response','SCTL (biomass average) vs. Range response',...
        'Bodysize vs. Range','Bodysize vs. Biomass','Bodysize vs. Population','Bodysize vs. SCTL (patch average)','Bodysize vs. SCTL (population average)',...
        'Bodysize vs. SCTL (biomass average)','Bodysize vs. Bodysize',...
        'Bodysize vs. Secondary extinctions response','Bodysize vs. Local diversity response','Bodysize vs. Max SCTL response','Bodysize vs. Link density response',...
        'Bodysize vs. Connectance response','Bodysize vs. Average SCTL response','Bodysize vs. Omnivory response','Bodysize vs. Population response','Bodysize vs. Bodysize response','Bodysize vs. Range response'};

    type_list=[repmat("double",1,119)];
    correlation_table=table('Size',[16 119],'VariableNames',column_list,'RowNames',row_list,'VariableTypes',type_list);

    % Mostly populate it with correlations of properties within the individual
    % simulations
    control_prop_xlim={[0 7],[-4000 90000],[-2000 25000],[0 4.5],[0 4.5],[0 4.5],[0 15]};

    % cycle simulations 1-15
    for j0=1:1:15

        % cycle control properties x
        for j1=1:1:7

            % cycle control properties y
            for j2=1:1:7
                corr_matrix=corrcoef(individual_data_frac{j0}(:,control_prop_num(j1)),individual_data_frac{j0}(:,control_prop_num(j2)));
                % Calculate the linear correlation coefficient
                correlation_table{j0,(j1-1)*17+j2}=corr_matrix(2,1);
            end

            % cycle response properties y
            for j3=1:1:10
                corr_matrix=corrcoef(individual_data_frac{j0}(:,control_prop_num(j1)),individual_data_frac{j0}(:,response_prop_num(j3)));
                % Calculate the linear correlation coefficient
                correlation_table{j0,(j1-1)*17+j3+7}=corr_matrix(2,1);
            end
        end
    end


    %%        
    % Correlate these modified responses against controls:
    corr_coef_values=zeros(7,10);

    for t1=1:1:7
        for t2=1:1:10

            corr_matrix=corrcoef(total_data_frac(:,control_prop_num(t1)),total_data_frac(:,response_prop_num(t2)));
            corr_coef_values(t1,t2)=corr_matrix(2,1); % Calculate the linear correlation coefficient
            correlation_table{16,(t1-1)*17+t2+7}=corr_matrix(2,1);

            if (print_images==true)
                hFig=figure('visible','on');
                clf(hFig);
                hold off;
                scatter(total_data_frac(:,control_prop_num(t1)),total_data_frac(:,response_prop_num(t2)),50,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',1.2);
                box on;
                xlabel(control_prop_name{t1});
                ylabel(response_prop_name{t2});
                xlim(control_prop_xlim{t1});
                ylim(response_frac_prop_ylim{t2});
                legend(num2str(corr_coef_values(t1,t2),3),'Location','NorthEast');
                printname=sprintf('SSS Spatial Stability data/FIGURES/Correlation/spatial_stability_correlation_%d_%d',t1,t2);
                print(printname,'-dpng','-r400','-painters');
                close(hFig);
            end

        end
    end



    %% Control vs Control
    % ALSO COLLECT CORRELATION BETWEEN THE RANGE OF SPECIES AND THEIR OTHER
    % PROPERTIES (OMNIVIORY AND SCTL), APART FROM DELETION
    corr_coef_values_control=zeros(7,7);
    for w1=1:1:7
        for w2=1:1:7

            corr_matrix=corrcoef(total_data_frac(:,control_prop_num(w1)),total_data_frac(:,control_prop_num(w2)));
            corr_coef_values_control(w1,w2)=corr_matrix(2,1); % Calculate the linear correlation coefficient
            correlation_table{16,(w1-1)*17+w2}=corr_matrix(2,1);

            if (print_images==true)
                hFig=figure('visible','on');
                clf(hFig);
                hold off;
                scatter(total_data_frac(:,control_prop_num(w1)),total_data_frac(:,control_prop_num(w2)),50,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',1.2);
                box on;
                xlabel(control_prop_name{w1});
                ylabel(control_prop_name{w2});
                xlim(control_prop_xlim{w1});
                ylim(control_prop_xlim{w2});
                legend(num2str(corr_coef_values_control(w1,w2),3),'Location','NorthEast');
                printname=sprintf('SSS Spatial Stability data/FIGURES/Correlation_control/spatial_stability_correlation_control_%d_%d',w1,w2);
                print(printname,'-dpng','-r400','-painters');
                close(hFig);
            end

        end
    end

    %% Boxes

    % Establish the limits of the boxes
    control_prop_boxlim={{[0.5 1.5],[1.5 2.5],[2.5 3.5],[3.5 4.5],[4.5 5.5],[5.5 6.5]}
        {[0 10000],[10000 20000],[20000 30000],[30000 40000],[40000 50000],[50000 60000],[60000 70000],[70000 80000],[80000 90000]}
        {[0 2000],[2000 4000],[4000 6000],[6000 8000],[8000 10000],[10000 12000],[12000 14000],[14000 16000],[16000 18000],[18000 20000],[20000 22000],[22000 24000]}
        {[0 0.5],[0.5 1],[1 1.5],[1.5 2],[2 2.5],[2.5 3],[3 3.5],[3.5 4],[4 4.5]} % for this one, want "1" to be in "[0.5 1]" rather than "[1 1.5]"
        {[0 0.5],[0.5 1],[1 1.5],[1.5 2],[2 2.5],[2.5 3],[3 3.5],[3.5 4],[4 4.5]}
        {[0 0.5],[0.5 1],[1 1.5],[1.5 2],[2 2.5],[2.5 3],[3 3.5],[3.5 4],[4 4.5]}
        {[0 1],[1 2],[2 3],[3 4],[4 5],[5 6],[6 7],[7 8],[8 9],[9 10],[10 11],[11 12],[12 13]}};
    
    % and the corresponding tick labels
    control_prop_labels = {{'1', '2', '3', '4', '5', '6'}
        {'5000', '15000', '25000', '35000', '45000', '55000', '65000', '75000', '85000'}
        {'1000', '3000', '5000', '7000', '9000', '11000', '13000', '15000', '17000', '19000' '21000', '23000'}
        {'0.25', '0.75', '1.25', '1.75', '2.25', '2.75', '3.25', '3.75', '4.25'}
        {'0.25', '0.75', '1.25', '1.75', '2.25', '2.75', '3.25', '3.75', '4.25'}
        {'0.25', '0.75', '1.25', '1.75', '2.25', '2.75', '3.25', '3.75', '4.25'}
        {'0.5', '1.5', '2.5', '3.5', '4.5', '5.5', '6.5', '7.5', '8.5', '9.5', '10.5', '11.5', '12.5'}};

    response_prop_boxlim={{[0 0.002],[0.002 0.004],[0.004 0.006],[0.006 0.008],[0.008 0.010],[0.010 0.012],[0.012 0.014],[0.014 0.016],[0.016 0.018],[0.018 0.020],[0.020 0.022]}
        {[-0.025 -0.020],[-0.020 -0.015],[-0.015 -0.010],[-0.010 -0.005],[-0.005 0],[0 0.005],[0.005 0.010],[0.010 0.015]}
        {[-0.05 -0.04],[-0.04 -0.03],[-0.03 -0.02],[-0.02 -0.01],[-0.01 0],[0 0.01],[0.01 0.02],[0.02 0.03],[0.03 0.04]}
        {[-0.020 -0.015],[-0.015 -0.010],[-0.010 -0.005],[-0.005 0],[0 0.005],[0.005 0.010],[0.010 0.015],[0.015 0.020]}
        {[-0.03 -0.02],[-0.02 -0.01],[-0.01 0],[0 0.01],[0.01 0.02],[0.02 0.03]}
        {[-0.025 -0.020],[-0.020 -0.015],[-0.015 -0.010],[-0.010 -0.005],[-0.005 0],[0 0.005],[0.005 0.010],[0.010 0.015],[0.015 0.020]}
        {[-0.06 -0.04],[-0.04 -0.02],[-0.02 0],[0 0.02],[0.02 0.04],[0.04 0.06]}
        {[-0.05 -0.04],[-0.04 -0.03],[-0.03 -0.02],[-0.02 -0.01],[-0.01 0],[0 0.01],[0.01 0.02],[0.02 0.03],[0.03 0.04],[0.04 0.05],[0.05 0.06],[0.06 0.07]}
        {[-0.020 -0.015],[-0.015 -0.010],[-0.010 -0.005],[-0.005 0],[0 0.005],[0.005 0.010],[0.010 0.015]}
        {[-0.02 -0.01],[-0.01 0],[0 0.01],[0.01 0.02],[0.02 0.03],[0.03 0.04]}};
    
    response_prop_labels={{'0.001','0.003','0.005','0.007','0.009','0.011','0.013','0.015','0.017','0.019','0.021'}
        {'-0.0225','-0.0175','-0.0125','-0.0075','-0.0025','0.0025','0.0075','0.0125'}
        {'-0.045','-0.035','-0.025','-0.015','-0.005','0.005','0.015','0.025','0.035'}
        {'-0.0175','-0.0125','-0.0075','-0.0025','0.0025','0.0075','0.0125','0.0175'}
        {'-0.025','-0.015','-0.005','0.005','0.015','0.025'}
        {'-0.0225','-0.0175','-0.0125','-0.0075','-0.0025','0.0025','0.0075','0.0125', '0.0175'}
        {'-0.05','-0.03','-0.01','0.01','0.03','0.05'}
        {'-0.045','-0.035','-0.025','-0.015','-0.005','0.005','0.015','0.025','0.035','0.045','0.055','0.065'}
        {'-0.0175','-0.0125','-0.0075','-0.0025','0.0025','0.0075','0.0125'}
        {'-0.015','-0.005','0.005','0.015','0.025','0.035'}};
    

    % Sort the values into the boxes

    box_array=cell(7,10);
    box_norm_array=cell(7,10);

    for r1=1:1:7 % control 
        for r2=1:1:10 % response

            % Determine the dimensions of the matrix
            num_x_boxes = size(control_prop_boxlim{r1},2);
            num_y_boxes = size(response_prop_boxlim{r2},2);

            % Generate the empty matrix
            box_array{r1,r2} = zeros(num_x_boxes,num_y_boxes);
            box_norm_array{r1,r2} = zeros(num_x_boxes,num_y_boxes);

            % Now go through each element
            for r3=1:1:size(total_data_frac,1)

                % For element r3, we go through the limits for both the x and
                % y boxes
                x_box=0;
                for r4=1:1:num_x_boxes
                    if (total_data_frac(r3,control_prop_num(r1))<control_prop_boxlim{r1}{r4}(2))
                        x_box=r4;
                        break;
                    end
                end
                % error detection
                if (x_box==0)
                    disp('Error! x_box not assigned.')
                    pause;
                end

                y_box=0;
                for r5=1:1:num_y_boxes
                    if (total_data_frac(r3,response_prop_num(r2))<response_prop_boxlim{r2}{r5}(2))
                        y_box=r5;
                        break;
                    end
                end
                % error detection
                if (y_box==0)
                    disp('Error! y_box not assigned.')
                    pause;
                end

                % Now count in the chosen box
                box_array{r1,r2}(x_box,y_box)=box_array{r1,r2}(x_box,y_box)+1;

            end 

            % Normalise the boxes by column, so that each column sums to 1.
            for h1=1:1:num_x_boxes
                column_sum = sum(box_array{r1,r2}(h1,:));
                box_norm_array{r1,r2}(h1,:) = box_array{r1,r2}(h1,:)./column_sum;
            end

            if (print_images==true)
                
                % Image the raw counts
                hFig=figure('visible','on');
                clf(hFig);
                hold off;
                imagesc(transpose(box_array{r1,r2}));
                set(gca,'YDir','normal');
                colorbar;
                colormap('jet');
                % caxis([8000 55000]);
                box on;
                xlabel(control_prop_name{r1});
                ylabel(response_prop_name{r2});
                set(gca,'xticklabels',control_prop_labels{r1})
                set(gca,'yticklabels',response_prop_labels{r2})
                printname=sprintf('SSS Spatial Stability data/FIGURES/Rawbox/spatial_stability_rawbox_%d_%d',r1,r2);
                print(printname,'-dpng','-r400','-painters');
                close(hFig);


                % Image the normalised data
                hFig=figure('visible','on');
                clf(hFig);
                hold off;
                imagesc(transpose(box_norm_array{r1,r2}));
                set(gca,'YDir','normal');
                colorbar;
                colormap('jet');
                % caxis([8000 55000]);
                box on;
                xlabel(control_prop_name{r1});
                ylabel(response_prop_name{r2});
                set(gca,'xticklabels',control_prop_labels{r1})
                set(gca,'yticklabels',response_prop_labels{r2})
                printname=sprintf('SSS Spatial Stability data/FIGURES/Normbox/spatial_stability_normbox_%d_%d',r1,r2);
                print(printname,'-dpng','-r400','-painters');
                close(hFig);
            end

        end
        
        %%%% Now control vs control boxes
        for r2=1:1:7 % control again

            % Determine the dimensions of the matrix
            num_x_boxes = size(control_prop_boxlim{r1},2);
            num_y_boxes = size(control_prop_boxlim{r2},2);

            % Generate the empty matrix
            box_array{r1,r2} = zeros(num_x_boxes,num_y_boxes);
            box_norm_array{r1,r2} = zeros(num_x_boxes,num_y_boxes);

            % Now go through each element
            for r3=1:1:size(total_data_frac,1)

                % For element r3, we go through the limits for both the x and
                % y boxes
                x_box=0;
                for r4=1:1:num_x_boxes
                    if (total_data_frac(r3,control_prop_num(r1))<control_prop_boxlim{r1}{r4}(2))
                        x_box=r4;
                        break;
                    end
                end
                % error detection
                if (x_box==0)
                    disp('Error! x_box not assigned.')
                    pause;
                end

                y_box=0;
                for r5=1:1:num_y_boxes
                    if (total_data_frac(r3,control_prop_num(r2))<control_prop_boxlim{r2}{r5}(2))
                        y_box=r5;
                        break;
                    end
                end
                % error detection
                if (y_box==0)
                    disp('Error! y_box not assigned.')
                    pause;
                end

                % Now count in the chosen box
                box_array{r1,r2}(x_box,y_box)=box_array{r1,r2}(x_box,y_box)+1;

            end 

            % Normalise the boxes by column, so that each column sums to 1.
            for h1=1:1:num_x_boxes
                column_sum = sum(box_array{r1,r2}(h1,:));
                box_norm_array{r1,r2}(h1,:) = box_array{r1,r2}(h1,:)./column_sum;
            end

            if (print_images==true)
                
                % Image the raw counts
                hFig=figure('visible','on');
                clf(hFig);
                hold off;
                imagesc(transpose(box_array{r1,r2}));
                set(gca,'YDir','normal');
                colorbar;
                colormap('jet');
                % caxis([8000 55000]);
                box on;
                xlabel(control_prop_name{r1});
                ylabel(control_prop_name{r2});
                set(gca,'xticklabels',control_prop_labels{r1})
                set(gca,'yticklabels',control_prop_labels{r2})
                printname=sprintf('SSS Spatial Stability data/FIGURES/Control_boxes/spatial_stability_control_rawbox_%d_%d',r1,r2);
                print(printname,'-dpng','-r400','-painters');
                close(hFig);


                % Image the normalised data
                hFig=figure('visible','on');
                clf(hFig);
                hold off;
                imagesc(transpose(box_norm_array{r1,r2}));
                set(gca,'YDir','normal');
                colorbar;
                colormap('jet');
                % caxis([8000 55000]);
                box on;
                xlabel(control_prop_name{r1});
                ylabel(control_prop_name{r2});
                set(gca,'xticklabels',control_prop_labels{r1})
                set(gca,'yticklabels',control_prop_labels{r2})
                printname=sprintf('SSS Spatial Stability data/FIGURES/Control_boxes/spatial_stability_control_normbox_%d_%d',r1,r2);
                print(printname,'-dpng','-r400','-painters');
                close(hFig);
            end

        end
        
    end


    %% Format the Correlations Table
    writetable(correlation_table,'SSS Spatial Stability data/spatial_stability_correlation_table_RAW.xlsx','WriteRowNames',true);

end

%% Statistics table for the control properties

type_list=[repmat("double",1,4)];
stats_table=table('Size',[7 4],'VariableNames',{'Minimum', 'Mean', '1SD', 'Maximum'},'RowNames',control_prop_name,'VariableTypes',type_list);

% Mean, standard deviation, and range of each cell-level property across
for t0 = 1:1:7
    x_data = total_data(:,control_prop_num(t0));
    stats_table{t0,1} = min(x_data);
    stats_table{t0,2} = mean(x_data);
    stats_table{t0,3} = std(x_data);
    stats_table{t0,4} = max(x_data);
end    
writetable(stats_table,'SSS Spatial Stability data/species_initial_stats_table.xlsx','WriteRowNames',true);


%% Multiple Linear Regression

%%%%%%%%
% Multiple regression of each of the 10 dependent ("response") variables
% in turn against a regression of 5 of the dependent ("control") variables
% (only one version of SCTL), plus SCTL^2 as this looks like it could have
% a quadratic relationship with FractionOfSecondaryExtinctions

is_regress=false;
if (is_regress==true)
    
    mult_regress_res=cell(10,2,2,2,2,2,2);
    
    best_rsquareadjusted=zeros(10,1);
    best_zlist=cell(10,1);
    best_zlist(cellfun('isempty',best_zlist))={[0,0,0,0,0,0]};
    list_of_indepvar = {'Range','Biomass','Population','SCTL (patch average)','Bodysize','SCTL (patch average) Squared'};
    
    % Combine all the independent variables into a single predictor matrix X
    X = zeros(7022,6);
    ind_prop_regress_rows = [4,7,8,9,12];
    for spec=1:1:7022
        for ind_prop=1:1:5
            X(spec,ind_prop) = total_data_frac(spec,ind_prop_regress_rows(ind_prop));
        end
        X(spec,6) = total_data_frac(spec,9)*total_data_frac(spec,9);
    end
    

    % Build a series of models for each dependent variable:
    for dep_prop=1:1:10
        
        y = total_data_frac(:,response_prop_num(dep_prop));

                        
        % Now try including each of the variables
        clearvars X3;
        for z1=0:1:1
             for z2=0:1:1
                 for z3=0:1:1
                     for z4=0:1:1
                         for z5=0:1:1
                             for z6=0:1:1

                                 z_list=[z1,z2,z3,z4,z5,z6];
                                 if (sum(z_list)>0)
                                     [X3,name_list] = x3_constructor(z_list,X,list_of_indepvar,response_prop_name{dep_prop});
                                     mult_regress_res{dep_prop,z1+1,z2+1,z3+1,z4+1,z5+1,z6+1}=fitlm(X3,y,'VarNames',name_list);

                                     % check if this is the
                                     % best model so far?

                                     rsq_adj = mult_regress_res{dep_prop,z1+1,z2+1,z3+1,z4+1,z5+1,z6+1}.Rsquared.Adjusted;
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Number of secondary extinctions:
    absolute_sec_ext = total_data(:,6);
    relative_sec_ext = total_data_frac(:,6);


    % For both absolute and relative secondary extinctions
    best_para_a = zeros(2,2);
    best_para_n = zeros(2,2);
    type_str = {'Absolute','Relative'};
    for type = 1:1:2
        for law = 1:1:2

            if (type==1)
                data = absolute_sec_ext;
                [values, edges] = histcounts(data);
                x_start = 1; % from the center of the largest bin (but not x=0)
                x_end = 8; % to the center of the last bin before height non-zero
                para_n_start = 0;
                para_n_end = 5;
                para_n_inc = 0.00001;
                if (law==1)
                    para_a_start=3000;
                    para_a_end=5000;
                    para_a_inc=1;
                else
                    para_a_start=3000;
                    para_a_end=5000;
                    para_a_inc=1;
                end
            else
                para_n_start = 200;
                para_n_end = 600;
                para_n_inc = 0.001;
                para_a_start=1000;
                para_a_end=5000;
                para_a_inc=1;
                if (law==1)
                    % need to shift
                    data = relative_sec_ext+1;
                    [values, edges] = histcounts(data,[0.999:0.002:1.03]);
                    x_start = 1.002;
                    x_end = 1.02;
                else
                    data = relative_sec_ext;
                    [values, edges] = histcounts(data,[-0.001:0.002:0.03]);
                    x_start = 0.002;
                    x_end = 0.02;
                end
            end
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
                    xlim([-0.02 2.2]);
                else
                    xlim([0 9]);
                end
            else
                if (law==1)
                    xlim([0 0.022]);
                else
                    xlim([0 0.022]);
                end
            end
            ylim([0 9]);
            legend_str = sprintf('Fitted model: R^2 = %s',r_str);
            legend('Data',legend_str,'Location','SouthWest')
            printname=sprintf('SSS Spatial Stability data/FIGURES/Extinction_distribution/spatial_stability_sec_ext_%d_%d',type,law);
            print(printname,'-dpng','-r400','-painters');
            close(hFig);
        end

        % Then create the original plot    
        if (type==1)
            data = absolute_sec_ext;
            x_start = 1;
            x_end = 8;
            x_inc = 0.0001;
            [values, edges] = histcounts(data);
        else
            data = relative_sec_ext;
            x_start = 0.002;
            x_end = 0.02;
            x_inc = 0.0000001;
            [values, edges] = histcounts(data,[-0.001:0.002:0.03]);
        end
        x_vals = [x_start:x_inc:x_end];
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
            ylim([0 4500]);
            power_str = sprintf('Power law: N(s) = %s(s)^{-%s}',a_str,n_str);
        else
            ylim([0 4500]);
            xlim([-0.001 0.022]);
            power_str = sprintf('Power law: N(s) = %s(s+1)^{-%s}',a_str,n_str);
        end
        a_str = num2str(best_para_a(type,2),6);
        n_str = num2str(best_para_n(type,2),3);
        exp_str = sprintf('Exponential law: N(s) = %sexp(-%ss)',a_str,n_str);
        legend(power_str,exp_str,'Data','Location','NorthEast')
        printname=sprintf('SSS Spatial Stability data/FIGURES/Extinction_distribution/spatial_stability_sec_ext_full_%d',type);
        print(printname,'-dpng','-r400','-painters');
        close(hFig);

    end
end

%% 
% Bodysize Distributions

% get the average non-resource bodysize in each simulation:
ave_bodysize = zeros(15,1);
for y1 = 1:1:15
    ave_bodysize(y1) = mean(individual_data{y1}(2:end,12));
end
raw_bodysize_set = total_data(:,12);
relative_bodysize_set = zeros(7022,1);
end_line=0;
for y2=1:1:15
    start_line = end_line+1;
    end_line = start_line + size(individual_data{y2},1)-2;
    relative_bodysize_set(start_line:end_line) = individual_data{y2}(2:end,12)/ave_bodysize(y2); 
end

% Raw histogram
hFig=figure('visible','on');
clf(hFig);
hold off;
histogram(raw_bodysize_set,60,'FaceColor',color_array{1},'LineWidth',0.8,'Normalization','probability');
box on;
xlabel('Bodysize');
ylabel('Probability');
ylim([0 0.1]);
printname=sprintf('SSS Spatial Stability data/FIGURES/Bodysize_distribution/bodysize_raw_hist');
print(printname,'-dpng','-r400','-painters');
close(hFig);

[freq,bin_edges] = histcounts(raw_bodysize_set,60);
freq = freq/sum(freq);
writematrix(transpose(bin_edges),'SSS Spatial Stability data/FIGURES/Bodysize_distribution/bodysize_bin_edges.csv');
writematrix(transpose(freq),'SSS Spatial Stability data/FIGURES/Bodysize_distribution/bodysize_freq.csv');

% Relative histogram
hFig=figure('visible','on');
clf(hFig);
hold off;
histogram(relative_bodysize_set,60,'FaceColor',color_array{1},'LineWidth',0.8,'Normalization','probability');
box on;
xlabel('Relative bodysize of species to meta-community mean');
ylabel('Probability');
ylim([0 0.1]);
printname=sprintf('SSS Spatial Stability data/FIGURES/Bodysize_distribution/bodysize_relative_hist');
print(printname,'-dpng','-r400','-painters');
close(hFig);
                
                

%% Function:

function [Xmatrix, var_names] = x3_constructor(z_list,allX,master_names,res_name)
    num_cols = sum(z_list);
    Xmatrix = zeros(7022,num_cols);
    current_col = 0;
    var_names = cell(num_cols+1,1);
    for z = 1:1:6
        if (z_list(z)==1)
            current_col=current_col+1;
            Xmatrix(:,current_col) = allX(:,z);
            var_names{current_col} = master_names{z};
        end
    end
    var_names{num_cols+1} = res_name;
end




