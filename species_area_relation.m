
%% Species-area relation
% How does the total gamma-diversity increase as we increase the number of
% patches considered?


% I think this will need the more detailed local data files for species
% stability, as these are the only ones that record which species are
% present in each specific patch and thus would allow us to determine
% beta-diversity and measures of gamma-diversity at different scales.

collect_data = false;
if (collect_data)
    
    clear;

    % Begin by loading the individual cell data
    simulation_data = cell(15,6,6);
    sim = 0;
    for sim_number = 100:1:120

        test_filename = sprintf('SAR Species Area Relation/DATA/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_bodysizeactlocal_1_1.txt',sim_number);
        if (exist(test_filename,'file')==2)

            sim = sim+1;

            for x = 1:1:6
                for y = 1:1:6
                    % store the individual patch data in a cell array
                    filename = sprintf('SAR Species Area Relation/DATA/WW_2020_01_1pt0e05_0pt6_mig1em2_%d_bodysizeactlocal_%d_%d.txt',sim_number,x,y);
                    placeholder = readmatrix(filename,'Delimiter',' ','ConsecutiveDelimitersRule','join');
                    % Now remove the extra column at the start
                    simulation_data{sim,x,y} = placeholder(:,2:end);
                    clearvars placeholder;
                end
            end
        end
    end
    


    % set the final vector 
    patch_size_vector = [1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 30, 36];
    averaged_diversity = zeros(size(patch_size_vector,2),1);
    % for each patch size
    for patch_size_iter = 1:1:18
        
        fprintf('Beginning patch size iterate: %d\n',patch_size_iter);

        % obtain the set of patch sub-cell co-ordinates
        patch_size = patch_size_vector(patch_size_iter);
        patch_set = patch_selector(patch_size);

        counter = 0;

        % Now look at each meta-community
        for sim = 1:1:15

            % now for each (multi-cell) patch element of patch_set
            for patch = 1:1:size(patch_set,2)

                % count the total number of unique species that occur across the
                % cells in this patch in this meta-community
                counter = counter + 1;

                cell_list = patch_set{patch};
                num_cells = size(cell_list,2);

                % build an empty list of all species numbers in this
                % meta-community
                num_rows = size(simulation_data{sim,1,1},1);
                num_spec = simulation_data{sim,1,1}(num_rows,4);

                species_list = zeros(num_spec,1);

                for patch_cell = 1:1:num_cells
                    % look at what species are in this cell
                    x = cell_list{patch_cell}(1);
                    y = cell_list{patch_cell}(2);
                    for k = 1:1:num_rows
                        if ((simulation_data{sim,x,y}(k,1) == 10000) && (simulation_data{sim,x,y}(k,5) == 1))
                            species_list(simulation_data{sim,x,y}(k,4)) = 1;
                        end
                    end
                end

                % count the total non-resource diversity in this (multi-cell) patch
                diversity = sum(species_list,'all') - num_cells;
                averaged_diversity(patch_size_iter) = averaged_diversity(patch_size_iter) + diversity;
            end
        end
        % average the summed-up diversity over all the patch selections and
        % meta-communities
        if (counter>0)
            averaged_diversity(patch_size_iter) = averaged_diversity(patch_size_iter)/double(counter);
        else
            averaged_diversity(patch_size_iter) = 0.0;
        end
    end
end



%% Plot the graph
set(0,'DefaultAxesFontSize',15);
color_array={ [0.05, 0.65, 0.25], [247/256, 182/256, 2/256], [0.15, 0.8, 0.4], [61/256, 89/256, 252/256], [252/256, 88/256, 76/256] };


% Determine best fit and R^2:
x_data = log(patch_size_vector);
y_data = log(averaged_diversity);

ss_tot = 0.0;
for val=1:1:size(x_data,2)
    ss_tot = ss_tot + (y_data(val)-mean(y_data))^2.0;
end

best_ss_res = 9999999.0;
parameter_c = 0.0;
parameter_m = 0.0;
for try_c = 1:1:1000
    for try_m = 1:1:1000
        model = transpose((1.5*try_m/1000)*x_data + 2+2*try_c/1000);
        ss_res=0;
        for val=1:1:size(x_data,2)
            ss_res = ss_res+(y_data(val)-model(val))^2.0;
        end
        if (ss_res < best_ss_res)
            best_ss_res = ss_res;
            parameter_c = 2+2*try_c/1000;
            parameter_m = 1.5*try_m/1000;
        end
    end
end
r_squared = 1-best_ss_res/ss_tot;


hFig=figure('visible','on');
clf(hFig);
hold off;
scatter(x_data, y_data, 50,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceAlpha',1);
hold on;
plot(x_data, parameter_m*x_data+parameter_c, 'LineWidth', 2, 'Color', color_array{2});
scatter(x_data, y_data, 50,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceAlpha',1,'HandleVisibility','off');
% this is just a repeat to get it overlaid on the plot
box on;
xlabel('log(Number of patches)');
ylabel('log(Total non-resource diversity)');
str1 = 'Data';
para_m_str = num2str(parameter_m,3);
para_c_str = num2str(parameter_c,3);
r_sq_str = num2str(r_squared,3);
str2 = sprintf('log(S) = %s log(A) + %s',para_m_str, para_c_str);
legend(str1, str2,'Location','NorthWest');
printname=sprintf('SAR Species Area Relation/species_area_log_log');
print(printname,'-dpng','-r400','-painters');
close(hFig);

% Non log graph
hFig=figure('visible','on');
clf(hFig);
hold off;
x_data = patch_size_vector;
y_data = averaged_diversity;
model = 20.086*x_data.^0.882;
scatter(x_data, y_data, 50,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceAlpha',1);
hold on;
plot(x_data, model, 'LineWidth', 2, 'Color', color_array{2});
scatter(x_data, y_data, 50,'o','filled','MarkerFaceColor',color_array{3},'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceAlpha',1,'HandleVisibility','off');
% this is just a repeat to get it overlaid on the plot
box on;
xlabel('Number of patches');
ylabel('Total non-resource diversity');
legend('Data','$S = 20.086 A^{0.882}$','Interpreter','latex','Location','NorthWest');
printname=sprintf('SAR Species Area Relation/species_area');
print(printname,'-dpng','-r400','-painters');
close(hFig);




    
%% Patch sub-cell co-ordinate selector    
function [patch_set] = patch_selector(patch_size)
    if (patch_size == 1)

        % for patch size 1 (1x1)
        patch_set = {};
        for m1=1:1:6
            for m2=1:1:6
                patch_set{end+1} = {[m1,m2]};
            end
        end

    elseif (patch_size == 2)

        % for patch size 2 (2x1)
        patch_set = {};
        for m1=1:1:5
            for m2=1:1:6
                patch_set{end+1} = {[m1,m2], [m1+1,m2]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1]};
            end
        end

    elseif (patch_size == 3)

        % for patch size 3 (3x1)
        patch_set = {};
        for m1=1:1:4
            for m2=1:1:6
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2]};
            end
        end

    elseif (patch_size == 4)

        % for patch size 4 (2x2) and (4x1)
        patch_set = {};
        for m1=1:1:3
            for m2=1:1:6
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2], [m2,m1+3]};
            end
        end
        for m1=1:1:5
            for m2=1:1:5
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1,m2+1], [m1+1,m2+1]};
            end
        end

    elseif (patch_size == 5)

        % for patch size 5 (5x1)
        patch_set = {};
        for m1=1:1:2
            for m2=1:1:6
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1+4,m2]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2], [m2,m1+3], [m2,m1+4]};
            end
        end

    elseif (patch_size == 6)

        % for patch size 6 (3x2) and (6x1)
        patch_set = {};
        for m1=1:1:1
            for m2=1:1:6
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1+4,m2], [m1+5,m2]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2], [m2,m1+3], [m2,m1+4], [m2,m1+5]};
            end
        end
        for m1=1:1:4
            for m2=1:1:5
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2], [m2+1,m1], [m2+1,m1+1], [m2+1,m1+2]};
            end
        end

    elseif (patch_size == 8)

        % for patch size 8 (4x2)
        patch_set = {};
        for m1=1:1:3
            for m2=1:1:5
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1], [m1+3,m2+1]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2], [m2,m1+3], [m2+1,m1], [m2+1,m1+1], [m2+1,m1+2], [m2+1,m1+3]};
            end
        end

    elseif (patch_size == 9)

        % for patch size 9 (3x3)
        patch_set = {};
        for m1=1:1:4
            for m2=1:1:4
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1], [m1,m2+2], [m1+1,m2+2], [m1+2,m2+2]}; 
            end
        end

    elseif (patch_size == 10)

        % for patch size 10 (5x2)
        patch_set = {};
        for m1=1:1:2
            for m2=1:1:5
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1+4,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1], [m1+3,m2+1], [m1+4,m2+1]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2], [m2,m1+3], [m2,m1+4], [m2+1,m1], [m2+1,m1+1], [m2+1,m1+2], [m2+1,m1+3], [m2+1,m1+4]};
            end
        end

    elseif (patch_size == 12)

        % for patch size 12 (4x3) and (6x2)
        patch_set = {};
        for m1=1:1:3
            for m2=1:1:4
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1], [m1+3,m2+1], [m1,m2+2], [m1+1,m2+2], [m1+2,m2+2], [m1+3,m2+2]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2], [m2,m1+3], [m2+1,m1], [m2+1,m1+1], [m2+1,m1+2], [m2+1,m1+3], [m2+2,m1], [m2+2,m1+1], [m2+2,m1+2], [m2+2,m1+3]};
            end
        end
        for m1=1:1:1
            for m2=1:1:5
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1+4,m2], [m1+5,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1], [m1+3,m2+1], [m1+4,m2+1], [m1+5,m2+1]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2], [m2,m1+3], [m2,m1+4], [m2,m1+5], [m2+1,m1], [m2+1,m1+1], [m2+1,m1+2], [m2+1,m1+3], [m2+1,m1+4], [m2+1,m1+5]};
            end
        end

    elseif (patch_size == 15)

        % for patch size 15 (5x3)
        patch_set = {};
        for m1=1:1:2
            for m2=1:1:4
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1+4,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1], [m1+3,m2+1], [m1+4,m2+1], [m1,m2+2], [m1+1,m2+2], [m1+2,m2+2], [m1+3,m2+2], [m1+4,m2+2]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2], [m2,m1+3], [m2,m1+4], [m2+1,m1], [m2+1,m1+1], [m2+1,m1+2], [m2+1,m1+3], [m2+1,m1+4], [m2+2,m1], [m2+2,m1+1], [m2+2,m1+2], [m2+2,m1+3], [m2+2,m1+4]};
            end
        end

    elseif (patch_size == 16)

        % for patch size 16 (4x4)
        patch_set = {};
        for m1=1:1:3
            for m2=1:1:3
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1], [m1+3,m2+1], [m1,m2+2], [m1+1,m2+2], [m1+2,m2+2], [m1+3,m2+2], [m1,m2+3], [m1+1,m2+3], [m1+2,m2+3], [m1+3,m2+3]};
            end
        end

    elseif (patch_size == 18)

        % for patch size 18 (6x3)
        patch_set = {};
        for m1=1:1:1
            for m2=1:1:4
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1+4,m2], [m1+5,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1], [m1+3,m2+1], [m1+4,m2+1], [m1+5,m2+1], [m1,m2+2], [m1+1,m2+2], [m1+2,m2+2], [m1+3,m2+2], [m1+4,m2+2], [m1+5,m2+2]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2], [m2,m1+3], [m2,m1+4], [m2,m1+5], [m2+1,m1], [m2+1,m1+1], [m2+1,m1+2], [m2+1,m1+3], [m2+1,m1+4], [m2+1,m1+5], [m2+2,m1], [m2+2,m1+1], [m2+2,m1+2], [m2+2,m1+3], [m2+2,m1+4], [m2+2,m1+5]};
            end
        end

    elseif (patch_size == 20)

        % for patch size 20 (5x4)
        patch_set = {};
        for m1=1:1:2
            for m2=1:1:3
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1+4,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1], [m1+3,m2+1], [m1+4,m2+1], [m1,m2+2], [m1+1,m2+2], [m1+2,m2+2], [m1+3,m2+2], [m1+4,m2+2], [m1,m2+3], [m1+1,m2+3], [m1+2,m2+3], [m1+3,m2+3], [m1+4,m2+3]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2], [m2,m1+3], [m2,m1+4], [m2+1,m1], [m2+1,m1+1], [m2+1,m1+2], [m2+1,m1+3], [m2+1,m1+4], [m2+2,m1], [m2+2,m1+1], [m2+2,m1+2], [m2+2,m1+3], [m2+2,m1+4], [m2+3,m1], [m2+3,m1+1], [m2+3,m1+2], [m2+3,m1+3], [m2+3,m1+4]};
            end
        end

    elseif (patch_size == 24)

        % for patch size 24 (6x4)
        patch_set = {};
        for m1=1:1:1
            for m2=1:1:3
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1+4,m2], [m1+5,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1], [m1+3,m2+1], [m1+4,m2+1], [m1+5,m2+1], [m1,m2+2], [m1+1,m2+2], [m1+2,m2+2], [m1+3,m2+2], [m1+4,m2+2], [m1+5,m2+2], [m1,m2+3], [m1+1,m2+3], [m1+2,m2+3], [m1+3,m2+3], [m1+4,m2+3], [m1+5,m2+3]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2], [m2,m1+3], [m2,m1+4], [m2,m1+5], [m2+1,m1], [m2+1,m1+1], [m2+1,m1+2], [m2+1,m1+3], [m2+1,m1+4], [m2+1,m1+5], [m2+2,m1], [m2+2,m1+1], [m2+2,m1+2], [m2+2,m1+3], [m2+2,m1+4], [m2+2,m1+5], [m2+3,m1], [m2+3,m1+1], [m2+3,m1+2], [m2+3,m1+3], [m2+3,m1+4], [m2+3,m1+5]};
            end
        end

    elseif (patch_size == 25)

        % for patch size 25 (5x5)
        patch_set = {};
        for m1=1:1:2
            for m2=1:1:2
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1+4,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1], [m1+3,m2+1], [m1+4,m2+1], [m1,m2+2], [m1+1,m2+2], [m1+2,m2+2], [m1+3,m2+2], [m1+4,m2+2], [m1,m2+3], [m1+1,m2+3], [m1+2,m2+3], [m1+3,m2+3], [m1+4,m2+3], [m1,m2+4], [m1+1,m2+4], [m1+2,m2+4], [m1+3,m2+4], [m1+4,m2+4]};
            end
        end

    elseif (patch_size == 30)

        % for patch size 30 (6x5)
        patch_set = {};
        for m1=1:1:1
            for m2=1:1:2
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1+4,m2], [m1+5,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1], [m1+3,m2+1], [m1+4,m2+1], [m1+5,m2+1], [m1,m2+2], [m1+1,m2+2], [m1+2,m2+2], [m1+3,m2+2], [m1+4,m2+2], [m1+5,m2+2], [m1,m2+3], [m1+1,m2+3], [m1+2,m2+3], [m1+3,m2+3], [m1+4,m2+3], [m1+5,m2+3], [m1,m2+4], [m1+1,m2+4], [m1+2,m2+4], [m1+3,m2+4], [m1+4,m2+4], [m1+5,m2+4]};
                patch_set{end+1} = {[m2,m1], [m2,m1+1], [m2,m1+2], [m2,m1+3], [m2,m1+4], [m2,m1+5], [m2+1,m1], [m2+1,m1+1], [m2+1,m1+2], [m2+1,m1+3], [m2+1,m1+4], [m2+1,m1+5], [m2+2,m1], [m2+2,m1+1], [m2+2,m1+2], [m2+2,m1+3], [m2+2,m1+4], [m2+2,m1+5], [m2+3,m1], [m2+3,m1+1], [m2+3,m1+2], [m2+3,m1+3], [m2+3,m1+4], [m2+3,m1+5], [m2+4,m1], [m2+4,m1+1], [m2+4,m1+2], [m2+4,m1+3], [m2+4,m1+4], [m2+4,m1+5]};
            end
        end

    elseif (patch_size == 36)

        % for patch size 36 (6x6)
        patch_set = {};
        for m1=1:1:1
            for m2=1:1:1
                patch_set{end+1} = {[m1,m2], [m1+1,m2], [m1+2,m2], [m1+3,m2], [m1+4,m2], [m1+5,m2], [m1,m2+1], [m1+1,m2+1], [m1+2,m2+1], [m1+3,m2+1], [m1+4,m2+1], [m1+5,m2+1], [m1,m2+2], [m1+1,m2+2], [m1+2,m2+2], [m1+3,m2+2], [m1+4,m2+2], [m1+5,m2+2], [m1,m2+3], [m1+1,m2+3], [m1+2,m2+3], [m1+3,m2+3], [m1+4,m2+3], [m1+5,m2+3], [m1,m2+4], [m1+1,m2+4], [m1+2,m2+4], [m1+3,m2+4], [m1+4,m2+4], [m1+5,m2+4], [m1,m2+5], [m1+1,m2+5], [m1+2,m2+5], [m1+3,m2+5], [m1+4,m2+5], [m1+5,m2+5]};
            end
        end

    end
end
