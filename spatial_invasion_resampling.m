%% spatial_invasion_resampling
%
% This is the function that is called to gather up the requisite number of
% resamples of specified number of datapoints each.
%
% Note that we are re-sampling with replacement as it us much easier!
%
% Make sure to declare these two variables, and the dataset in advance in
% the the spatial_invasion_analysis_V2.m script.

dataset = sortrows(dataset,13);

dataset_rows = size(dataset,1);
bin_edges = readmatrix('SSI Species Spatial Invasion/bodysize_bin_edges.csv');
bin_probability = readmatrix('SSI Species Spatial Invasion/bodysize_freq.csv'); 
datapoints_per_bin = round(bin_probability*num_datapoints); % this is how many need to be drawn for each bin, rounded

set_of_samples = cell(num_resamples,1); % this is the important output

for k = 1:1:num_resamples

    sampled_dataset = zeros(ceil(num_datapoints/2),24); % this will store our drawn samples.
    % Note that it may be altered during the process, due to rounding. So
    % make it smaller than necessary to ensure no empty rows, but expect a
    % result [0.9-1.1]*num_datapoints in number of rows.
    sampled_datapoint_counter = 0; % this indexes the rows for storage

    % then loop over the bins
    for bin = 1:1:60

        % identify the bounds of the rows in the dataset
        bin_lower_bound = bin_edges(bin);
        bin_upper_bound = bin_edges(bin+1);
        bin_start = 0;
        bin_end = dataset_rows;
        for row=1:1:dataset_rows
            if ((dataset(row,13)>bin_lower_bound) && (bin_start==0))
                bin_start = row;
            end
            if (dataset(row,13)>bin_upper_bound)
                bin_end = row;
                break
            end
        end

        % draw each of the samples for that bin
        if (datapoints_per_bin(bin) > 0)
            for bin_samp = 1:1:datapoints_per_bin(bin)

                sampled_datapoint_counter = sampled_datapoint_counter + 1;
                random_row = round(bin_start + (bin_end - bin_start)*rand);
                sampled_dataset(sampled_datapoint_counter,:) = dataset(random_row,:);

            end
        end 
    end 
    set_of_samples{k} = sampled_dataset;
    clearvars sampled_dataset;
end 