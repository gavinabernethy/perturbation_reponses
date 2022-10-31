%% Individual Patch Deletion experiments
%
% For each simulation, we have individually deleted each of the 36 patches
% using all 2x2 combinations of methods.
% Before each set, the initial state of the system is printed.
% Therefore, there should be 4x(36+1)=148 lines in each file.
% 
% However, for the initial properties of the cell that was deleted (aside
% from its diversity and position), we would need to consult a different
% file which is the _meta_robust_initial_


% IN THIS VERSION (UPDATED OVER _OLD), WE CREATE A UNIQUE X-VECTOR AND
% Y-VECTOR FOR EACH INPUT-OUTPUT FOR PLOTS, CORRELATION, AND REGRESSION.
% THIS IS TO EXCLUDE DEFAULT VALUES OF 0 IN EITHER THE INPUT OR OUTPUT THAT
% OCCUR DUE TO AN EMPTY SET - EG. "AVERAGE BODYSIZE OF SOME SET OF SPECIES"
% IS RECORDED AS ZERO BECAUSE THERE WERE ACTUALLY NO SUCH SPECIES.

clear;
set(0,'DefaultAxesFontSize',15);
color_array={ [0.05, 0.65, 0.25], [247/256, 182/256, 2/256], [0.15, 0.8, 0.4], [61/256, 89/256, 252/256], [252/256, 88/256, 76/256] };



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collect data:
% Loop over each of the 15 repeat simulations

count = 0;
individual_data_MRC=cell(15,1);
individual_data_LI=cell(15,1);
individual_data_RC=cell(15,36);

for s1 = 100:1:120
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % _meta_robust_corr:
    filename=sprintf('IPD Individual Patch Deletion data/EXTRA/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_meta_robust_corr_00010000.txt',s1);
    if (exist(filename,'file')==2)
        
        count = count+1;
        
        % store the individual data in a cell array
        individual_placeholder=readmatrix(filename);
        individual_data_MRC{count}=individual_placeholder;
        clear individual_placeholder;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % _meta_robust_localinitial:
    filename=sprintf('Local Initial data/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_meta_robust_localinitial_00010000.txt',s1);
    if (exist(filename,'file')==2)
        % store the individual data in a cell array
        individual_placeholder=readmatrix(filename);
        individual_data_LI{count}=individual_placeholder;
        clear individual_placeholder;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % _robustnesscorrelation:
    patch_count = 0;
    for x=1:1:6
        for y=1:1:6
            patch_count = patch_count+1;
            filename=sprintf('Robustness Correlation data/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_robustnesscorrelation_%d_%d.txt',s1,x,y);
            if (exist(filename,'file')==2)
                % store the individual data in a cell array
                individual_placeholder=readmatrix(filename,'Delimiter',' ','ConsecutiveDelimitersRule','join');
                individual_data_RC{count,patch_count}=individual_placeholder;
                clear individual_placeholder;
            end
        end
    end
end

% Remove the faulty extra column from the RobustnessCorrelation data:
individual_data_RC2 = cell(15,36);
for r1 = 1:1:15
    for r2 = 1:1:36
        individual_data_RC2{r1,r2} =  individual_data_RC{r1,r2}(:,2:14);
    end
end
clearvars individual_data_RC;
individual_data_RC = individual_data_RC2;
clearvars individual_data_RC2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Individual MRC data is therefore {15}(148x43)
% 
% For each simulation, rows 1, 38, 75, 112 should be identical as the
% initial state is given before the next list of cells with a different
% type of deletion.
%
% i2, i3, k3, n_spec_store, i5, i6, x_max, y_max, n_spec, deleted_cell_topology, num_local_spec_stored(i5,i6), current_n_spec_global,
% current_max_sc_tl_global, current_link_density_global, current_connectance_global, current_ave_sc_tl_global,
% current_omnivory_frac_global, current_ave_population_global, current_ave_bodysize_global, current_ave_distribution_global
%
% Columns:
%
% 1 = i2 = disruption permanence (1 = temporary; 2 = permanent)
% 2 = i3 = disruption effect (1 = elimination; 2 = displacement)
% 3 = k3 = time step, this is always 10000.
% 4 = n_spec_store = number of initial species in the global ecosystem,
% this is always the same within each row of the simulation.
% 5 = i5 = the x-coordinate of the patch that is deleted (1-6).
% 6 = i6 = the y-coordinate of the patch that is deleted (1-6).
% 7 = x_max = always 6.
% 8 = y_max = always 6.
% 9 = n_spec = the global diversity remaining after the perturbation. The
% change between column 4 and this is the impact of the deletion on global
% biodiversity.
% 10 = deleted_cell_topology = number of patches that the deleted patch was
% connected to, INCLUDING ITSELF (3 on a side; 4 on a corner; 5 on the 
% interior).
% 11 = num_local_spec_stored(i5,i6) = the number of (not necessarily
% unique) species that were in the patch that was deleted. i.e. it's local
% diversity
% 12 = current_n_spec_global = average local diversity
% 13 = current_max_sc_tl_global = average maximum SCTL
% 14 = current_link_density_global = average link density
% 15 = current_connectance_global = average connectance
% 16 = current_ave_sc_tl_global = average of the average (across species
% within the patch) SCTL
% 17 = current_omnivory_frac_global = average fraction of the species that
% are omnivores
% 18 = current_ave_population_global = average of the average (across
% species within the patch) population size of species
% 19 = current_ave_bodysize_global = average of the average (across
% species within the patch) bodysize of species
% 20 = current_ave_distribution_global = average of the average (across
% species within the patch) range of species
%
% For columns 12-20, each of these is an average of a local property across
% the 36 (or 35 if perturbation was permanent) patches following the 
% perturbation. Compare in each row with the initial value in the initial
% rows 1, 38, 75 or 112 to see the difference.
%
%
% Following the revision in _FIX2:
%
% 21 = num_residents = the number of species that were present in the
%  perturbed patch
% 22 = num_residents_survived = the number of species that were present in the
%  perturbed patch prior to the experiment, and in either the perturbed
%  patch or any neighbouring patch following the perturbation.
% 23 = num_unique_residents = the number of species that were present in
%  the perturbed patch and NOT in ANY of the immediate neighbours prior to
%  the perturbation.
% 24 = num_unique_residents_survived  = the number of species that were present in
%  the perturbed patch and NOT in ANY of the immediate neighbours prior to
%  the perturbation, and who were present in either the perturbed patch or
%  any neighbouring patch following the perturbation.
% 25 = num_unique_residents_invaded = the number of species that were present in
%  the perturbed patch and NOT in ANY of the immediate neighbours prior to
%  the perturbation, and who were present in any neighbouring patch
%  following the perturbation.
% 26 = num_neighbours = the number of species who were present in any
%  neighbouring patch prior to the experiment.
% 27 = num_neighbours_survived = the number of species who were present in any
%  neighbouring patch prior to the experiment, and in either the perturbed
%  patch or any neighbouring patch following the perturbation.
% 28 = num_unique_neighbours = the number of species who were present in any
%  (or multiple) neighbouring patch(s) prior to the experiment but NOT in the
%  perturbed patch.
% 29 = num_unique_neighbours_survived = the number of species who were present in any
%  (or multiple) neighbouring patch(s) prior to the experiment but NOT in the
%  perturbed patch, and who were present in either the perturbed patch or
%  any neighbouring patch following the perturbation.
%
% 30 = ave_bodysize_residents
% 31 = ave_bodysize_residents_survived
% 32 = ave_bodysize_residents_died
% 33 = ave_bodysize_unique_residents
% 34 = ave_bodysize_unique_residents_survived
% 35 = ave_bodysize_unique_residents_died
% 36 = ave_bodysize_unique_residents_invaded
% 37 = ave_bodysize_unique_residents_failedtoinvade = the average bodysize
%   of all species who were residents in the perturbed patch and not in
%   neighbouring patches prior to the perturbation, and who did not
%   successfully invade a neighbouring patch (including the possibility
%   that the died altogether -  that is "failed to invade" is NOT just 
%   "survived and failed to invade")
% 38 = ave_bodysize_neighbours
% 39 = ave_bodysize_neighbours_survived
% 40 = ave_bodysize_neighbours_died
% 41 = ave_bodysize_unique_neighbours
% 42 = ave_bodysize_unique_neighbours_survived
% 43 = ave_bodysize_unique_neighbours_died




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the _localinitial_ data:
% 
% k3, y1, y2, i2, i3, i1, i7, i8, n_spec_store, x_max, y_max, n_spec, num_res_temp, current_n_spec_local(y1,y2),
% current_max_sc_tl_local(y1,y2), current_link_density_local(y1,y2), current_connectance_local(y1,y2),
% current_ave_sc_tl_local(y1,y2), current_omnivory_frac_local(y1,y2), current_ave_population_local(y1,y2),
% current_ave_bodysize_local(y1,y2), current_ave_distribution_local(y1,y2)
%
% Columns:
%
% 1 = k3 = timestep, always 10000
% 2 = y1 = x-value of the patch
% 3 = y2 = y-value of the patch
% 4 = i2 = disruption permanence (1 = temporary; 2 = permanent). Should be
% 1 as this data is only collected on the first iteration.
% 5 = i3 = disruption effect (1 = elimination; 2 = displacement). Should be
% 1 as this data is only collected on the first iteration.
% 6 = i1 = presence of reserves (1 = no reserves; 2 = reserves). Should be
% 1 as this data is only collected on the first iteration.
% 7 = i7 = which reserve set? Should be 1 as this data is only collected
% on the first iteration.
% 8 = i8 = which deletion sequence? Should be 1 as this data is only
% collected on the first iteration.
% 9 = n_spec_store = total global diversity (always the same for all
% patches within the same imulation)
% 10 = x_max = always 6
% 11 = y_max = always 6
% 12 = n_spec = total global diversity (should be the same as col. 9)
% 13 = num_res_temp = number of resources (always 36)
% 14 = current_n_spec_local(y1,y2) = local diversity
% 15 = current_max_sc_tl_local(y1,y2) = local maximum SCTL
% 16 = current_link_density_local(y1,y2) = local link density
% 17 = current_connectance_local(y1,y2) = local connectance
% 18 = current_ave_sc_tl_local(y1,y2) = local average SCTL
% 19 = current_omnivory_frac_local(y1,y2) = local omnivory fraction
% 20 = current_ave_population_local(y1,y2) = local average population
% 21 = current_ave_bodysize_local(y1,y2) = local average bodysize
% 22 = current_ave_distribution_local(y1,y2) = local average range


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the _robustnesscorrelation data:
%
% k3, n_spec_store, n_spec_temp_store, a4, n_spec, secondary_extinctions,
% sds_bodysize_store(a4,1), sds_population_store(1,a4,hj12,hj13), temp_trophic_levels(a4),
% prey_links(a4), predator_links(a4), sds_bodysize_store(a4,3), sds_bodysize_store(a4,4)

% 3 = the true local diversity of the cell for this text file.





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Re-organise the important data:

% We want a single matrix containing all of the 19 independent variables,
% one with all 25 dependent variables (as relative change where needed).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INDEPENDENT VARIABLES

% Prepare the independent variables from the MR _corr file:
% (i) local diversity, (ii) topology-1, (iii)-(x) resident/neighbour
% properties.
%
% Then the additional initial properties of patches that were collected from 
% the _localinitial_ file:
% (i) max SCTL, (ii) link density, (iii) connectance, (iv) ave SCTL, 
% (v) omnivory fraction, (vi) ave population, (vii) ave bodysize, 
% (viii) ave range
%
% Finally, the true diversity after isolation from the
% _robustnesscorrelation file.

ind_prop_label={'Local diversity','Patch topology',...
    ...
    'Residents', 'Unique residents', 'Neighbours', 'Unique neighbours', 'Ave. bodysize residents',...
    'Ave. bodysize unique residents', 'Ave. bodysize neighbours', 'Ave. bodysize unique neighbours',...
    ...
    'Maximum SCTL','Link density','Connectance', 'Average SCTL',...
    'Fraction of omnivorous species','Average population','Average bodysize','Average range',...
    ...
    'True local diversity',...
    ...
    'Residents : Neighbours', 'Unique residents : unique neighbours',...
    'Bodysize of unique residents : bodysize of unique neighbours'};

ind_prop_label_SHORT={'Local diversity','Patch topology',...
    'Residents', 'Uni. residents', 'Neighbours', 'Uni. neighbours', 'BS residents',...
    'BS uni. residents', 'BS neighbours', 'BS uni. neighbours',...
    'Max SCTL','Link density','Connectance', 'Ave SCTL',...
    'Omnivory','Population','Bodysize','Range', 'True local diversity',...
    'Residents:neighbours', 'Uni.res. : uni.neighbours', 'BS uni.res.: uni.neighbours'};


independent_var = zeros(15,36,19); % (sim,patch,prop)
independent_prop_colsMRC = [11,10,21,23,26,28,30,33,38,41]; % from _meta_robust_corr_
independent_prop_colsLI = [15,16,17,18,19,20,21,22]; % from _localinitial_
independent_prop_colsRC = [3]; % from _robustnesscorrelation


for sim = 1:1:15
    for patch = 1:1:36
        for prop = 1:1:22
            
            % _meta_robust_corr:
            if (prop<11)
                col = independent_prop_colsMRC(prop);
                row=patch+1;
                if (prop==2)
                    % Topology (modify to reduce by 1)
                    independent_var(sim,patch,prop) = individual_data_MRC{sim}(row,col)-1.0;
                else
                    independent_var(sim,patch,prop) = individual_data_MRC{sim}(row,col);
                end
            
            % _local_initial:
            elseif (prop<19)
                row=patch;
                col = independent_prop_colsLI(prop-10);
                independent_var(sim,patch,prop) = individual_data_LI{sim}(row,col);
            
            % _robustness_correlation_x_y:
            elseif (prop<20)
                col = independent_prop_colsRC(prop-18);
                independent_var(sim,patch,prop) = individual_data_RC{sim,patch}(1,col);
                
            % Now the ratios of residents:neighbours,
            % unique_residents:unique_neighbours
            % and bodysize unique_residents : bodysize unique_neighbours
            elseif (prop==20)
                independent_var(sim,patch,prop) = independent_var(sim,patch,3)/independent_var(sim,patch,5);
            elseif (prop==21)
                independent_var(sim,patch,prop) = independent_var(sim,patch,4)/independent_var(sim,patch,6);
            elseif (prop==22)
                independent_var(sim,patch,prop) = independent_var(sim,patch,8)/independent_var(sim,patch,10);
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEPENDENT VARIABLES

% From the _meta_robust_corr files:
% - (i) fractional change in global diversity
% - fractional change in average:
% (ii) local diversity, (iii) max SCTL, (iv) link density, (v) connectance,
% (vi) ave SCTL, (vii) omnivory fraction, (viii) ave population,
% (ix) ave bodysize, (x) ave range.
%
% Then the additional properties from the same file.
% 
% (xi) 22 = num_residents_survived
% (xii) 24 = num_unique_residents_survived
% (xiii) 25 = num_unique_residents_invaded
% (xiv) 27 = num_neighbours_survived
% (xv) 29 = num_unique_neighbours_survived
% (xvi) 31 = ave_bodysize_residents_survived
% (xvii) 32 = ave_bodysize_residents_died
% (xviii) 34 = ave_bodysize_unique_residents_survived
% (xix) 35 = ave_bodysize_unique_residents_died
% (xx) 36 = ave_bodysize_unique_residents_invaded
% (xxi) 37 = ave_bodysize_unique_residents_failedtoinvade
% (xxii) 39 = ave_bodysize_neighbours_survived
% (xxiii) 40 = ave_bodysize_neighbours_died
% (xxiv) 42 = ave_bodysize_unique_neighbours_survived
% (xxv) 43 = ave_bodysize_unique_neighbours_died

dependent_var = zeros(15,2,2,36,30); % (simulation, permanence, effect, cell, property)
dependent_prop_cols = [9,12,13,14,15,16,17,18,19,20,22,24,25,27,29,31,32,34,35,36,37,39,40,42,43];
dep_prop_label = {'Fractional change in global diversity', 'Fractional change in average local diversity',...
    'Fractional change in average max. SCTL', 'Fractional change in average link density', 'Fractional change in average connectance',...
    'Fractional change in average ave. SCTL', 'Fractional change in average omnivory fraction',...
    'Fractional change in average ave. population', 'Fractional change in average ave. bodysize', 'Fractional change in average ave. range',...
     ...
    'Fraction of residents survived', 'Fraction of unique residents survived', 'Fraction of unique residents invaded',...
    'Fraction of neighbours survived', 'Fraction of unique neighbours survived',... 
    ...
	{'Relative average bodysize of'; 'residents who survived'}, {'Relative average bodysize of'; 'residents who died'},...
    {'Relative average bodysize of'; 'unique residents who survived'}, {'Relative average bodysize of'; 'unique residents who died'},...
    {'Relative average bodysize of'; 'unique residents who invaded'}, {'Relative average bodysize of unique'; 'residents who failed-to-invade'},...
    {'Relative average bodysize of'; 'neighbours who survived'}, {'Relative average bodysize of'; 'neighbours who died'},...
    {'Relative average bodysize of'; 'unique neighbours who survived'}, {'Relative average bodysize of'; 'unique neighbours who died'},...
    ...
    'Fraction of residents who died', 'Fraction of unique residents who died', {'Fraction of unique residents'; 'who failed-to-invade'},...
    'Fraction of neighbours who died', 'Fraction of unique neighbours who died'};

dep_prop_label_SHORT = {'Global diversity', 'Local diversity',...
    'Max SCTL', 'Link density', 'Connectance', 'Ave SCTL', 'Omnivory',...
    'Population', 'Bodysize', 'Range',...
    'Residents survived', 'Uni. residents survived', 'Uni. residents invaded',...
    'Neighbours survived', 'Uni. neighbours survived',... 
	'BS surviving residents', 'BS dead residents',...
    'BS uni. residents survived', 'BS uni. residents died',...
    'BS uni. residents invaded', 'BS uni. residents not-invade',...
    'BS neighbours survived', 'BS neighbours died',...
    'BS uni. neighbours survived', 'BS uni. neighbours died',...
    'Residents died', 'Uni. residents died', 'Uni. residents not-invade',...
    'Neighbours died', 'Uni. neighbours died'};

for sim = 1:1:15
    for perm = 1:1:2
        for effect = 1:1:2
            for patch = 1:1:36

                row=1+74*(perm-1)+37*(effect-1)+patch;

                for prop = 1:1:30
                    
                    if (prop<26)
                        col = dependent_prop_cols(prop);
                    end
                    
                    if (prop<11)
                        % these are about how global properties are
                        % changing when this patch is perturbed
                        dependent_var(sim,perm,effect,patch,prop) = (individual_data_MRC{sim}(row,col)-individual_data_MRC{sim}(1,col))/individual_data_MRC{sim}(1,col);
                    
                    elseif (prop<26)
                        
                        % calculate resident/neighbour survivors as
                        % fractions of their group, and average bodysizes
                        % now should be relative to the average
                        % for the same group prior to the
                        % perturbation? (so that we can understand for
                        % example of the average bodysize has gone up or
                        % down as a result of the perturbation)
                        indep_var_basis = [21,23,23,26,28,30,30,33,33,33,33,38,38,41,41];
                        col2 = indep_var_basis(prop-10);
                        dependent_var(sim,perm,effect,patch,prop) = individual_data_MRC{sim}(row,col)/individual_data_MRC{sim}(row,col2);
                        % these are about how properties
                        % (resident/neighbour) that are DEPENDENT ON THE
                        % CELL BEING PERTURBED respond to its perturbation
                        %
                        % hence the final part is ..}(row,col2) rather than
                        % ..}(1,col2), as in the case of the other props.
                        %
                        % This also automatically will exclude any empty
                        % cases, as if for example there are zero unique
                        % residents, then they will have a default value of
                        % zero for their average bodysize. Then when we go
                        % to calculate the ratio by which this output has
                        % changed in such an empty case, we recieve NaN
                        % which means that MatLab will handle it by
                        % excluding that datapoint automatically.
                    
                        
                        % create the last 5 properties from the ones that
                        % already exist. Make sure to refer to the raw
                        % data, not the fractions that have already been
                        % rescaled:
                      
                                     
                    elseif (prop == 26) % - fraction residents died = (num residents - num residents survived)/num residents
                        dependent_var(sim,perm,effect,patch,prop) = abs((individual_data_MRC{sim}(row,21)-individual_data_MRC{sim}(row,22))/individual_data_MRC{sim}(row,21));
                    elseif (prop == 27) % - fraction unique residents died = (num unique residents - num unique residents survived)/num unique residents
                        dependent_var(sim,perm,effect,patch,prop) = abs((individual_data_MRC{sim}(row,23)-individual_data_MRC{sim}(row,24))/individual_data_MRC{sim}(row,23));
                    elseif (prop == 28) % - fraction unique residents failed to invade = (num unique residents - num unique residents invaded)/num unique residents
                        dependent_var(sim,perm,effect,patch,prop) = abs((individual_data_MRC{sim}(row,23)-individual_data_MRC{sim}(row,25))/individual_data_MRC{sim}(row,23));
                    elseif (prop == 29) % - fraction neighbours died = (num neighbours - num neighbours survived)/num neighbours
                        dependent_var(sim,perm,effect,patch,prop) = abs((individual_data_MRC{sim}(row,26)-individual_data_MRC{sim}(row,27))/individual_data_MRC{sim}(row,26));
                    elseif (prop == 30) % - fraction unique neighbours died = (num unique neighbours - num unique neighbours survived)/num unique neighbours
                        dependent_var(sim,perm,effect,patch,prop) = abs((individual_data_MRC{sim}(row,28)-individual_data_MRC{sim}(row,29))/individual_data_MRC{sim}(row,28));
                        
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FILTER THE DATA FOR DEFAULT NON-MEANINGFUL VALUES
%
% We create a set of X-vector and Y-vector for every (X,Y) pair, excluding
% all of the non-meaningful values.

corrected_data_set = cell(16,2,2,21,30); % 1-15 sims (16 is all), 2 permanence, 2 effect, 19 inputs, 30 outputs
% each entry in this cell array should itself be a 2x1 cell array
% containing two equal length column vectors (X and Y).


% which properties do each property need to be compared against?
indprop_check_list = [0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 4];
% do not count ave bodysize unique residents if there ARE no unique
% residents, so check this value also of independent_var;

%                     1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30          
depprop_check_list = [0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  26 12 27 13 28 0  29 0  30 0  0  0  0  0];
% We add new columns (26-30) to the dependent_var data, so that they
% can be referenced here to to check that they are non-zero:
% - 26: num residents died
% - 27: num unique residents died
% - 28: num unique residents failed to invade
% - 29: num neighbours died
% - 30: num unique neighbours died


for indprop = 1:1:22
    for depprop = 1:1:30
        for perm = 1:1:2
            for effect = 1:1:2
                
                count_total_admissable = 0;
                count_cell = 0;
                x_total = [];
                y_total = [];
                cell_list = [];
                
                for sim = 1:1:15
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % iterate through each x-y value and check that neither
                    % is inadmissible (zeros or NaN) before appending to 
                    % the x-y vector pair
                    count_ind_admissable = 0;
                    x_ind = [];
                    y_ind = [];
                    
                    for patch = 1:1:36
                        
                        count_cell = count_cell + 1;
                        
                        this_patch_succeeds = 1;
                        % check the input and output for NaN and Inf first:
                        if (isnan(independent_var(sim,patch,indprop))) || (isnan(dependent_var(sim,perm,effect,patch,depprop))) || (isinf(independent_var(sim,patch,indprop))) || (isinf(dependent_var(sim,perm,effect,patch,depprop)))
                            this_patch_succeeds = 0;
                        end
                        % check the reference properties
                        if (depprop_check_list(depprop) ~= 0)
                            % check the dependent reference
                            if (dependent_var(sim,perm,effect,patch,depprop_check_list(depprop))==0) || (isnan(dependent_var(sim,perm,effect,patch,depprop_check_list(depprop)))) || (isinf(dependent_var(sim,perm,effect,patch,depprop_check_list(depprop))))
                                this_patch_succeeds = 0;
                            end
                        end
                        if (indprop_check_list(indprop) ~= 0)
                            % check the independent reference
                            if (independent_var(sim,patch,indprop_check_list(indprop))==0) || (isnan(independent_var(sim,patch,indprop_check_list(indprop)))) || (isinf(independent_var(sim,patch,indprop_check_list(indprop))))
                                this_patch_succeeds = 0;
                            end
                        end
                        
                        if (this_patch_succeeds == 1)
                            % append the acceptable data
                            count_ind_admissable = count_ind_admissable + 1;
                            count_total_admissable = count_total_admissable + 1;
                            x_ind(count_ind_admissable) = independent_var(sim,patch,indprop);
                            y_ind(count_ind_admissable) = dependent_var(sim,perm,effect,patch,depprop);
                            x_total(count_total_admissable) = independent_var(sim,patch,indprop);
                            y_total(count_total_admissable) = dependent_var(sim,perm,effect,patch,depprop);
                            cell_list(count_total_admissable) = count_cell;
                        end
                    
                    end
                    % store the individual simulation data
                    corrected_data_set{sim, perm, effect, indprop, depprop} = {x_ind; y_ind};
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                end
                % store the concatenated data
                corrected_data_set{16, perm, effect, indprop, depprop} = {cell_list; x_total; y_total};
            end
        end
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scatter plots and Pearson correlation coefficient:

% Conduct a multivariate analysis of the 19 independent variables, against
% each of the 30 dependent variables.

% ISSUE:
% Note that often in the _specincell_ data, there are species whose SCTL is <1 or even 0. Should this be the case, and what would it mean?
% May have to make sure this information is not used in the analysis if it is unreliable.

indep_prop_lims = {
    [0 40], [1 5], [0 40], [-5 30], [0 100],...
    [0 100], [1 12], [1 12], [1 12], [1 12],...
    [0 6], [0 4], [0 0.6], [0 3], [0 1],...
    [0 3500], [1 12], [1 5], [0 40], [-0.1 1.4], [-0.1 1.4] [0 3]};
 
dep_prop_lims = {
    [-0.1 0.01], [-0.1 0.1], [-0.1 0.1], [-0.1 0.1], [-0.12 0.12],...
    [-0.06 0.06], [-0.12 0.12], [-0.5 1.0], [-0.06 0.06], [-0.1 0.2],...
    [0 1.2], [-0.05 1.05], [-0.05 1.05], [0.5 1.05], [0.5 1.05],...
    [0.2 2.4], [0.2 2.4], [0.2 2.4], [0.2 2.4], [0.2 2.4],...
    [0.2 2.4], [0.2 2.4], [0.2 2.4], [0.2 2.4], [0.2 2.4],...
    [-0.05 1.05], [-0.05 1.05], [-0.05 1.05], [-0.05 0.5], [-0.05 0.5]};

is_plots=false;
if (is_plots==true)
    % First, a plot (one for each of the 30 dependent variables):
    % x-axis: patch 1-36 for each simulation 1-15
    % y-axis: the dependent variable, different colour for each of the 4 groups

    x_label={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'};

    for dep_prop=1:1:30

        hFig=figure('visible','on');
        clf(hFig);
        hold off;

        count=0;
        for perm=1:1:2
            for effect=1:1:2
                count = count+1; % this is for colour for each type
                
                x_vec = corrected_data_set{16, perm, effect, 1, dep_prop}{1}; % the number of the cell that this data came from
                y_vec = corrected_data_set{16, perm, effect, 1, dep_prop}{3};
                % output properties may not be well-defined for all 540
                % patches, but they should not be further limited by input
                % property 1 at least.

                scatter(x_vec,y_vec,20,'o','filled','MarkerFaceColor',color_array{count+1},'MarkerEdgeColor','k','LineWidth',0.8);
                hold on;

            end
        end

        box on;
        xlabel('Simulation and patch');
        set(gca,'xtick',1:37:540)
        xticklabels(x_label);
        ylabel(dep_prop_label{dep_prop});
        xlim([0 541]);
        ylim(dep_prop_lims{dep_prop});
        legend('Temporary Elimination','Temporary Displacement','Permanent Elimination','Permanent Displacement','Location','SouthEast');
        printname=sprintf('IPD Individual Patch Deletion data/FIGURES/ind_patch_deletion_%d',dep_prop);
        print(printname,'-dpng','-r400','-painters');
        close(hFig);

    end    

    %%%%%%%%
    % Now plot the dependent variables against the independent variables,
    % rather than the patch number
    type_list = [repmat("double",1,4)];
    correlation_table_DivOnly = table('Size',[22 4],'VariableNames',{'Temporary Elimination','Temporary Displacement','Permanent Elimination','Permanent Displacement'},'RowNames',ind_prop_label_SHORT,'VariableTypes',type_list);
    
    correlation_table = cell(2,2); % separate table for each perturbation type
    type_list = [repmat("double",1,660)];
    column_list = {};
    for ind_prop=1:1:22
        for dep_prop=1:1:30
            column_list{end+1} = [dep_prop_label_SHORT{dep_prop}, ' vs ', ind_prop_label_SHORT{ind_prop}];
        end
    end
    row_list = {'Sim 1','Sim 2','Sim 3','Sim 4','Sim 5','Sim 6','Sim 7','Sim 8','Sim 9','Sim 10','Sim 11','Sim 12','Sim 13','Sim 14','Sim 15','Combined'};
    for perm=1:1:2
        for effect=1:1:2
            correlation_table{perm,effect} = table('Size',[16 660],'VariableNames',column_list,'RowNames',row_list,'VariableTypes',type_list);
        end
    end
    table_column = 0;
    for ind_prop=1:1:22

        for dep_prop=1:1:30
            
            table_column = table_column+1;
            hFig=figure('visible','on');
            clf(hFig);
            hold off;
            
            count=0;
            for perm=1:1:2
                for effect=1:1:2
                    
                    count=count+1;

                    x_vec = corrected_data_set{16, perm, effect, ind_prop, dep_prop}{2};
                    y_vec = corrected_data_set{16, perm, effect, ind_prop, dep_prop}{3};
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                    % Calculate the Pearson correlation coefficient:
                    for sim=1:1:15
                        corr_matrix=corrcoef(corrected_data_set{sim, perm, effect, ind_prop, dep_prop}{1},corrected_data_set{sim, perm, effect, ind_prop, dep_prop}{2});
                        if isnan(corr_matrix(2,1))
                            correlation_table{perm,effect}{sim,table_column}=0.0;
                        else
                            correlation_table{perm,effect}{sim,table_column}=corr_matrix(2,1);
                        end
                    end
                    % Then over all patches and simulations:
                    corr_matrix=corrcoef(x_vec,y_vec);
                    if isnan(corr_matrix(2,1))
                        overall_coef = 0.0;
                    else
                        overall_coef = corr_matrix(2,1);
                    end
                    correlation_table{perm,effect}{16,table_column} = overall_coef;
                    
                    % also store the values for global diversity against
                    % the properties of all simulations in a dedicated
                    % table
                    if (dep_prop==1)
                        correlation_table_DivOnly{ind_prop,count} = overall_coef;
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                    scatter(x_vec,y_vec,25,'o','filled','MarkerFaceColor',color_array{count+1},'MarkerEdgeColor','k','LineWidth',0.8);
                    hold on;

                end
            end

            box on;
            xlabel(ind_prop_label{ind_prop});
            ylabel(dep_prop_label{dep_prop});
            xlim(indep_prop_lims{ind_prop});
            ylim(dep_prop_lims{dep_prop});
            legend('Temporary Elimination','Temporary Displacement','Permanent Elimination','Permanent Displacement','Location','SouthEast');
            printname=sprintf('IPD Individual Patch Deletion data/FIGURES/ind_patch_deletion_%d_%d',ind_prop,dep_prop);
            print(printname,'-dpng','-r400','-painters');
            close(hFig);

        end
        
    end
    
    % Write the five Correlations Tables
    for perm=1:1:2
        for effect=1:1:2
            tablename=sprintf('IPD Individual Patch Deletion data/ipd_correlation_table_%d_%d_RAW.xlsx',perm,effect);
            writetable(correlation_table{perm,effect},tablename,'WriteRowNames',true);
        end
    end
    writetable(correlation_table_DivOnly,'IPD Individual Patch Deletion data/ipd_correlation_table_diversityOnly_RAW.xlsx','WriteRowNames',true);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Average, maximum and minimum:

% Then calculate the min, max and average of each of the 30 dependent properties over all of the 15*36 
% patches where appropriate and present this in a table against the 4 kinds of perturbation

type_list = [repmat("double",1,30)];
row_list = {'Temporary elimination', 'Temporary displacement', 'Permanent elimination', 'Permanent displacement'};
average_table = table('Size',[4 30],'VariableNames',dep_prop_label_SHORT,'RowNames',row_list,'VariableTypes',type_list);
min_table = table('Size',[4 30],'VariableNames',dep_prop_label_SHORT,'RowNames',row_list,'VariableTypes',type_list);
max_table = table('Size',[4 30],'VariableNames',dep_prop_label_SHORT,'RowNames',row_list,'VariableTypes',type_list);
table_row = 0;
for perm=1:1:2
    for effect = 1:1:2
        table_row = table_row+1;
        for dep_prop=1:1:30
            average_table{table_row, dep_prop} = mean(corrected_data_set{16, perm, effect, 1, dep_prop}{3});
            min_val = min(corrected_data_set{16, perm, effect, 1, dep_prop}{3});
            max_val = max(corrected_data_set{16, perm, effect, 1, dep_prop}{3});
            if isempty(min_val)
                min_val = 0;
            end
            if isempty(max_val)
                max_val = 0;
            end
            min_table{table_row, dep_prop}  = min_val;
            max_table{table_row, dep_prop} = max_val;
        end
    end
end
writetable(average_table,'IPD Individual Patch Deletion data/ipd_average_table_RAW.xlsx','WriteRowNames',true);
writetable(min_table,'IPD Individual Patch Deletion data/ipd_min_table_RAW.xlsx','WriteRowNames',true);
writetable(max_table,'IPD Individual Patch Deletion data/ipd_max_table_RAW.xlsx','WriteRowNames',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% When is displacement worse than elimination?

is_disp_subtract_elim_test = true;
if (is_disp_subtract_elim_test)
    temp_disp_subtract_elim = corrected_data_set{16,1,2,1,1}{3} - corrected_data_set{16,1,1,1,1}{3};
    perm_disp_subtract_elim = corrected_data_set{16,2,2,1,1}{3} - corrected_data_set{16,2,1,1,1}{3};

    type_list = [repmat("double",1,2)];
    row_list = ind_prop_label_SHORT;
    disp_subt_elim_table = table('Size',[22 2],'VariableNames',{'Temporary','Permanent'},'RowNames',row_list,'VariableTypes',type_list);
    for indprop=1:1:22

        hFig=figure('visible','on');
        clf(hFig);
        hold off;

        x_vec = [];
        temp_disp_subtract_elim_TEMP = [];
        perm_disp_subtract_elim_TEMP = [];
        count_admissable = 0;
        count = 0;
        for sim = 1:1:15
            for patch = 1:1:36     
                count = count+1;
                this_patch_succeeds = 1;
                % check the input and output for NaN and Inf first:
                if (isnan(independent_var(sim,patch,indprop))) || (isinf(independent_var(sim,patch,indprop)))
                    this_patch_succeeds = 0;
                end
                if (indprop_check_list(indprop) ~= 0)
                    % check the independent reference
                    if (independent_var(sim,patch,indprop_check_list(indprop))==0) || (isnan(independent_var(sim,patch,indprop_check_list(indprop)))) || (isinf(independent_var(sim,patch,indprop_check_list(indprop))))
                        this_patch_succeeds = 0;
                    end
                end
                if (this_patch_succeeds == 1)
                    % append the acceptable data
                    count_admissable = count_admissable + 1;
                    x_vec(count_admissable) = independent_var(sim,patch,indprop);
                    temp_disp_subtract_elim_TEMP(count_admissable) = temp_disp_subtract_elim(count);
                    perm_disp_subtract_elim_TEMP(count_admissable) = perm_disp_subtract_elim(count);
                end
            end
        end
        %yline(0,'HandleVisibility','off','LineWidth',1);
        hold on;
        % shading
        X = [indep_prop_lims{indprop}(1) indep_prop_lims{indprop}(2)];
        Y = [0.03 0; 0.03 0];
        a = area(X,Y,'HandleVisibility','off');
        a(1).FaceColor = [0.9 0.9 0.9];
        a(2).FaceColor = [1 1 1];
        scatter(x_vec,temp_disp_subtract_elim_TEMP,25,'o','filled','MarkerFaceColor',color_array{1},'MarkerEdgeColor','k','LineWidth',0.85,'MarkerFaceAlpha',1);
        hold on;
        scatter(x_vec,perm_disp_subtract_elim_TEMP,25,'o','filled','MarkerFaceColor',color_array{2},'MarkerEdgeColor','k','LineWidth',0.85,'MarkerFaceAlpha',0.85);
        xlim(indep_prop_lims{indprop});
        ylim([-0.05 0.03]);
        box on;
        xlabel(ind_prop_label{indprop});
        ylabel({'Difference between fractional change in diversity'; ' due to displacement and change due to elimination'});
        legend('Temporary','Permanent','Location','SouthEast');
        printname=sprintf('IPD Individual Patch Deletion data/FIGURES/ipd_elimination_subtract_disp_%d',indprop);
        print(printname,'-dpng','-r400','-painters');
        close(hFig);

        % Calculate the Pearson correlation coefficient:
        corr_matrix=corrcoef(x_vec,temp_disp_subtract_elim_TEMP);
        disp_subt_elim_table{indprop,1}=corr_matrix(2,1);
        corr_matrix=corrcoef(x_vec,perm_disp_subtract_elim_TEMP);
        disp_subt_elim_table{indprop,2}=corr_matrix(2,1);


    end
    writetable(disp_subt_elim_table,'IPD Individual Patch Deletion data/ipd_dispsubtelim_table_RAW.xlsx','WriteRowNames',true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regression


% Then, multiple regression within each of these 4 groups, of each of the
% 30 dependent variabes in turn against a regression of all of the 21 
% independent variables

is_regress=false;
if (is_regress==true)
    
    mult_regress_res=cell(2,2,30,2,2,2,2,2,2,2,2);
    
    best_rsquareadjusted=zeros(2,2,30);
    best_zlist=cell(2,2,30);
    best_zlist(cellfun('isempty',best_zlist))={[0,0,0,0,0,0,0,0]};
    
    % Combine all independent variables into a single predictor matrix X:
    X = zeros(540,22);
    for sim=1:1:15
        for ind_prop=1:1:22
            X(1+36*(sim-1):36*sim,ind_prop) = independent_var(sim,:,ind_prop);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %     Full set of independent variables:
    %         1. 'Local diversity', 2. 'Patch topology', 3. 'Residents', 4. 'Uni. residents'
    %         5. 'Neighbours', 6. 'Unique neighbours', 7. 'Bodysize residents'
    %         8. 'Bodysize uni. residents', 9. 'Bodysize neighbours', 10. 'Bodysize uni. neighbours'
    %         11. 'Max SCTL', 12. 'Link density', 13. 'Connectance', 14. 'Ave SCTL'
    %         15. 'Omnivory', 16. 'Population', 17. 'Bodysize', 18. 'Range', 19. 'True local diversity'
    %         20. 'Residents:neighbours', 21. 'Uni.residents:uni.neighbours'
    %         22. 'Bodysize Uni.residents : Bodysize uni.neighbours'
    
    % Then select 8 for this analysis:
    X2 = zeros(540,10);
    x2_labels = cell(10,1);
    x2_column_list = [1,2,4,21,8,22,16,18];
    for x=1:1:8
        X2(:,x)=X(:,x2_column_list(x));
        x2_labels{x} = ind_prop_label_SHORT{x2_column_list(x)};
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This is necessary to avoid naming conflict between "Local diversity"
    % in the input and ouput variables in the linear model:
    dep_prop_label_OUT = {'Output Global diversity', 'Output Local diversity',...
    'Output Max SCTL', 'Output Link density', 'Output Connectance', 'Output Ave SCTL', 'Output Omnivory',...
    'Output Population', 'Output Bodysize', 'Output Range',...
    'Residents survived', 'Uni. residents survived', 'Uni. residents invaded',...
    'Neighbours survived', 'Uni. neighbours survived',... 
	'Bodysize surviving residents', 'Bodysize dead residents',...
    'Bodysize uni. residents survived', 'Bodysize uni. residents died',...
    'Bodysize uni. residents invaded', 'Bodysize uni. residents not-invade',...
    'Bodysize neighbours survived', 'Bodysize neighbours died',...
    'Bodysize uni. neighbours survived', 'Bodysize uni. neighbours died',...
    'Residents died', 'Uni. residents died', 'Uni. residents not-invade',...
    'Neighbours died', 'Uni. neighbours died'};


    for perm=1:1:2
        for effect=1:1:2
            
            % Recall that:
            % disruption permanence (1 = temporary; 2 = permanent)
            % disruption effect (1 = elimination; 2 = displacement)

            % Build a series of models for each dependent variable:
            for dep_prop=1:1:30

                y = zeros(540,1);
                for sim=1:1:15
                    y(1+36*(sim-1):36*sim)=dependent_var(sim,perm,effect,:,dep_prop);
                end
                
                
                % Now try including each of the variables
                clearvars X3;
                for z1=0:1:1
                     for z2=0:1:1
                         for z3=0:1:1
                             for z4=0:1:1
                                 for z5=0:1:1
                                     for z6=0:1:1
                                         for z7=0:1:1
                                             for z8=0:1:1
                                                 z_list=[z1,z2,z3,z4,z5,z6,z7,z8];
                                                 if (sum(z_list)>0)
                                                     [X3,name_list] = x3_constructor(z_list,X2,x2_labels,dep_prop_label_OUT{dep_prop});
                                                     mult_regress_res{perm,effect,dep_prop,z1+1,z2+1,z3+1,z4+1,z5+1,z6+1,z7+1,z8+1}=fitlm(X3,y,'VarNames',name_list);

                                                     % check if this is the
                                                     % best model so far?

                                                     rsq_adj = mult_regress_res{perm,effect,dep_prop,z1+1,z2+1,z3+1,z4+1,z5+1,z6+1,z7+1,z8+1}.Rsquared.Adjusted;
                                                     if (rsq_adj > best_rsquareadjusted(perm,effect,dep_prop))
                                                         best_rsquareadjusted(perm,effect,dep_prop) = rsq_adj;
                                                         best_zlist{perm,effect,dep_prop} = z_list;
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
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANCOVA

% Then, an ANCOVA test to determine the role of the 4 sampled "groups" 
% (permanence and effect of the perturbation).

% We only use this concerning the extinctions/ change to global diversity

% TEST %%%%
ancova_test=false;
if (ancova_test==true)
    pert=zeros(2160,1);
    pert(1:540)=1;
    pert(541:1080)=2;
    pert(1081:1620)=3;
    pert(1621:2160)=4;

    global_diversity_change=zeros(2160,1);
    initial_local_diversity=zeros(2160,1);

    x_vec=zeros(540,1);
    for sim=1:1:15
        x_vec(1+36*(sim-1):36*sim)=independent_var(sim,:,1);
    end
    initial_local_diversity(1:540)=x_vec;
    initial_local_diversity(541:1080)=x_vec;
    initial_local_diversity(1081:1620)=x_vec;
    initial_local_diversity(1621:2160)=x_vec;

    count=0;
    for perm=1:1:2
        for effect=1:1:2

            count=count+1;

            % Prepare the y-vector
            y_vec=zeros(540,1);
            for sim=1:1:15
                y_vec(1+36*(sim-1):36*sim)=dependent_var(sim,perm,effect,:,1);
            end

            global_diversity_change(1+540*(count-1):540*count)=y_vec;

        end
    end
    [h,atab,ctab,stats] = aoctool(initial_local_diversity,global_diversity_change,pert);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function used in the multiple linear regression:

function [Xmatrix, var_names] = x3_constructor(z_list,allX,master_names,res_name)
    num_cols = sum(z_list);
    Xmatrix = zeros(540,num_cols);
    current_col = 0;
    var_names = cell(num_cols+1,1);
    for z = 1:1:8
        if (z_list(z)==1)
            current_col=current_col+1;
            Xmatrix(:,current_col) = allX(:,z);
            var_names{current_col} = master_names{z};
        end
    end
    var_names{num_cols+1} = res_name;
end







        