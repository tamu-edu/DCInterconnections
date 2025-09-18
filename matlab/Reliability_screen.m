

%clear all; clc;
t_start = tic;

%load case file
mpc_base      = loadcase('YOUR MATPOWER Case File');

%load annual demand, wind, solar, hydro profiles
D = load('YOUR_demand_profiles.mat');   area_load_arr = D.area_load;   clear D
W = load('YOUR_wind_profiles.mat');     wind_arr      = W.wind_MW;     wind_bus = W.wind_bus; clear W
S = load('YOUR_solar_profiles.mat');    solar_arr     = S.solar_MW;    solar_bus = S.solar_bus; clear S
H = load('YOUR_hydro_profiles.mat');    hydro_arr     = H.hydro_MW;    hydro_bus = H.hydro_bus; clear H


fprintf('===== Starting N-1 OPF Contingency Analysis =====\n');
define_constants;

% Base case
mpc_full = mpc_base;

% Get number of branches
n_branches = size(mpc_full.branch, 1);
n_gens = size(mpc_full.gen, 1);

% Generator matrix column indices
GEN_STATUS = 8;  

% Generator status
gen_status = mpc_full.gen(:, GEN_STATUS);

% Count online and offline generators
n_online = sum(gen_status == 1);
n_offline = sum(gen_status == 0);

% PMIN analysis
gen_pmin = mpc_full.gen(:, PMIN);
n_pmin_zero = sum(gen_pmin == 0);

% Filter only online generators
online_idx = find(gen_status == 1);
online_pmin = mpc_full.gen(online_idx, PMIN);
online_pmax = mpc_full.gen(online_idx, PMAX);

% Total capacity
total_online_pmin = sum(online_pmin);
total_online_pmax = sum(online_pmax);

% Display gen results
fprintf('Total generators: %d\n', length(mpc_full.gen(:,1)));
fprintf('Online generators: %d\n', n_online);
fprintf('Offline generators: %d\n', n_offline);
fprintf('Total PMAX of online generators: %.2f MW\n', total_online_pmax);



% Select target day(s) and hour(s) for reliability assessment (e.g. peak day, peak hour)
Days = 224; %peak day, add other concerned days
hrs = 16; %peak hour, add other concerned hours

%the set of concerning hours 
selected_hours = (Days-1).*24 + hrs;


mpc_hour = mpc_full;
% mpc_hour.bus(:, PD) = area_load_arr(selected_hour, :)';

% total_load = sum(mpc_hour.bus(:, PD));
% fprintf('Total load at hour %d: %.2f MW\n', selected_hour, total_load);


%map the gen list to mpc.gen idx
gen_buses = mpc_full.gen(:, 1);

%aasign the renewable buses to mpc gen buses
wind_idx = assign_gen_P(mpc_full, wind_bus);
solar_idx = assign_gen_P(mpc_full, solar_bus);
hydro_idx = assign_gen_P(mpc_full, hydro_bus);


ng = size(mpc_full.gen, 1);
all_idx = (1:ng)';
renewable_idx = unique([wind_idx; solar_idx; hydro_idx]);
conventional_idx = setdiff(all_idx, renewable_idx);

fprintf('Renewable generators: %d\n', length(renewable_idx));
fprintf('Conventional generators: %d\n', length(conventional_idx));


%N-1 contingencies
%line and transformer contingencies
%filter the n-1 cases causing islanding
from_bus = mpc_full.branch(:, 1);
to_bus = mpc_full.branch(:, 2);
n_branch = length(from_bus);

all_buses = unique([from_bus; to_bus]);
kv_threshold = 69;

% map bus numbers to consecutive indices
[~, ~, bus_idx_map] = unique([from_bus; to_bus]);
from_idx = bus_idx_map(1:n_branch);
to_idx   = bus_idx_map(n_branch+1:end);

bus_kv = mpc_full.bus(:, 10); 
kv_from = bus_kv(from_idx);  
kv_to   = bus_kv(to_idx);  
kv_min  = min(kv_from, kv_to);

% original graph
G = graph(from_idx, to_idx);

valid_branch_indices = [];

%check if islanding will be inccurred
j=0;
k=0;
for i = 1:n_branch
    if kv_min(i) < kv_threshold
        j=j+1;
        continue;
    end

    Gi = rmedge(G, from_idx(i), to_idx(i));
    bins = conncomp(Gi);
    if numel(unique(bins)) == 1
        valid_branch_indices(end+1) = i;
    else
        k=k+1;
    end
end

fprintf('Total branches: %d\n', n_branch);
fprintf('Non-islanding N-1 cases: %d\n', numel(valid_branch_indices));
n_valid_ctgs = length(valid_branch_indices);


%gen contengencies
%only keep online gens
valid_gen_indices = find(mpc_hour.gen(:, GEN_STATUS) == 1);  
n_valid_gen_ctgs = length(valid_gen_indices);



%locations to be scanned
bus_info = readtable('bus_info.csv'); %detailed bus information

is_bus = (bus_info.BusNomVolt == 115) | (bus_info.BusNomVolt == 161) | (bus_info.BusNomVolt == 230);
filtered_bus_info = bus_info(is_bus, :);

extra_load_buses_list = filtered_bus_info.BusNum + 3000000; %to match the default bus name format in the case file


%data center load levels
capacity_list = [300, 500, 640, 680, 800, 1000, 1280, 1360, 1600, 2000]; %or other capacities concerned


output_csv = sprintf('reliability_screening_results.csv');
result_per_bus = cell(length(extra_load_buses_list), 1);
results_gen_per_bus = cell(length(extra_load_buses_list), 1);


%parallel computing pool configurations
cluster = parcluster('local');
cluster.JobStorageLocation = scratch_dir;
parpool(cluster, str2double(getenv('SLURM_CPUS_PER_TASK')));
disp(gcp);


%Main loop
for j = 1:length(extra_load_buses_list)
    bus = extra_load_buses_list(j);
    fprintf('-> Bus %d \n', bus);
    
    converged_ctg_mask = true(1, n_valid_ctgs);
    result_per_bus{j} = cell(length(capacity_list), 1);

    for cap_idx = 1:length(capacity_list)
        additional_mw = capacity_list(cap_idx);
        fprintf('-> Capacity %d MW\n', additional_mw);
        
        n_selected_hours = length(selected_hours);
        shed_values = nan(n_valid_ctgs, n_selected_hours);
        shed_flags = false(n_valid_ctgs, n_selected_hours);
            
        mpc_hour = mpc_full;
        mpc_hour.bus(:, PD) = area_load_arr(selected_hour, :)' * load_factor;
        row_idx = find(mpc_hour.bus(:, 1) == bus);
        if isempty(row_idx)
            warning('Bus %d not found in the case for hour %d.', bus, selected_hour);
        else
            mpc_hour.bus(row_idx, PD) = mpc_hour.bus(row_idx, PD) + additional_mw;
        end
        mpc_hour.bus(:, QD) = mpc_hour.bus(:, PD) * tan(theta);
    
        mpc_hour.gen(wind_idx, PMAX) = wind_arr(selected_hour, :)';
        mpc_hour.gen(solar_idx, PMAX) = solar_arr(selected_hour, :)';
        mpc_hour.gen(hydro_idx, PMAX) = hydro_arr(selected_hour, :)';
    
        % N-1
        mpopt_dc = mpoption('opf.dc.solver', 'Gurobi', 'verbose', 0, 'out.all', 0);
        to_run_mask = converged_ctg_mask;
        ctg_list_to_run = valid_branch_indices(to_run_mask);
        n_run = length(ctg_list_to_run);
        fprintf('n_run: %d', n_run);
    
        local_converged = false(1, n_run);
        local_failed_mask = false(1, n_run);
    
        % Step 1: DCOPF
        parfor i_valid = 1:n_run
            ctg_idx = ctg_list_to_run(i_valid);
            mpc_ctg = mpc_hour;
            mpc_ctg.branch(ctg_idx, :) = [];
            % fprintf('ctg_idx:%d', ctg_idx);
            try
                result = rundcopf(mpc_ctg, mpopt_dc); % no slack
                % fprintf('result.success:%d', result.success);
                if result.success
                    local_converged(i_valid) = true;
                else
                    local_failed_mask(i_valid) = true;
                end
            catch ME
                fprintf('DCOPF error in ctg %d: %s\n', ctg_idx, ME.message);
                local_failed_mask(i_valid) = true;
            end
        end
    
        % % Step 2: calculate load shedding amount, if necessary (for future use)
        % if sum(local_failed_mask) < n_branch*0.1 % only consider those buses whose failed contingency cases are lower than a threshold
        %     fprintf('-> Hour %d: %d failed ctgs, entering Step 2\n', selected_hour, sum(local_failed_mask));
        % 
        %     shed_tmp = nan(1, n_run);
        %     shed_flag_tmp = false(1, n_run);
        % 
        %     parfor i_valid = 1:n_run
        %         if ~local_converged(i_valid)
        %             ctg_idx = ctg_list_to_run(i_valid);
        %             mpc_ctg = mpc_hour;
        %             mpc_ctg.branch(ctg_idx, :) = [];
        %             ng_orig = size(mpc_hour.gen, 1);
        %             mpc_ctg = add_load_shedding_slack(mpc_ctg);
        % 
        % 
        %             try
        %                 result = rundcopf(mpc_ctg, mpopt_dc);
        %                 if result.success
        %                     shed = sum(result.gen(ng_orig+1:end, PG));
        %                     shed_tmp(i_valid) = shed;
        %                     shed_flag_tmp(i_valid) = shed > 1e-3;
        %                     % fprintf('shed:%f\n', shed);
        %                 end
        %             catch ME
        %                 fprintf('DCOPF error in ctg %d: %s\n', ctg_idx, ME.message);
        %             end
        %         end
        %     end
        %     shed_values(to_run_mask, hr_idx) = shed_tmp;
        %     shed_flags(to_run_mask, hr_idx) = shed_flag_tmp;
        % else
        %     fprintf('-> Hour %d: too many failed ctgs (%d), skipping Step 2\n', selected_hour, sum(local_failed_mask));
        % end
    
        converged_ctg_mask(to_run_mask) = local_converged;

        shed_count = sum(shed_flags(:));
        avg_shed = mean(shed_values(shed_flags), 'omitnan');

        % Save this cap result
        result_per_bus{j}{cap_idx} = struct( ...
            'dc_converged_count', sum(converged_ctg_mask), ...
            'shed_count', shed_count, ...
            'shed_avg', avg_shed);

        fprintf('===== Completed N-1 OPF Contingency Analysis =====\n');
    end
end

%gen
for j = 1:length(extra_load_buses_list)
    bus = extra_load_buses_list(j);
    fprintf('-> Bus %d \n', bus);
    
    converged_ctg_mask = true(1, n_valid_gen_ctgs);
    results_gen_per_bus{j} = cell(length(capacity_list), 1);

    for cap_idx = 1:length(capacity_list)
        additional_mw = capacity_list(cap_idx);
        fprintf('-> Capacity %d MW\n', additional_mw);
        
        n_selected_hours = length(selected_hours);
        shed_values = nan(n_valid_gen_ctgs, n_selected_hours);
        shed_flags = false(n_valid_gen_ctgs, n_selected_hours);
        
        selected_hour = selected_hours(1);
    
        mpc_hour = mpc_full;
        mpc_hour.bus(:, PD) = area_load_arr(selected_hour, :)' * load_factor;
        row_idx = find(mpc_hour.bus(:, 1) == bus);
        if isempty(row_idx)
            warning('Bus %d not found in the case for hour %d.', bus, selected_hour);
        else
            mpc_hour.bus(row_idx, PD) = mpc_hour.bus(row_idx, PD) + additional_mw;
        end
        mpc_hour.bus(:, QD) = mpc_hour.bus(:, PD) * tan(theta);
    
        mpc_hour.gen(wind_idx, PMAX) = wind_arr(selected_hour, :)';
        mpc_hour.gen(solar_idx, PMAX) = solar_arr(selected_hour, :)';
        mpc_hour.gen(hydro_idx, PMAX) = hydro_arr(selected_hour, :)';
    
        mpopt_dc = mpoption('opf.dc.solver', 'Gurobi', 'verbose', 0, 'out.all', 0);
        to_run_mask = converged_ctg_mask;
        ctg_list_to_run = valid_gen_indices(to_run_mask);
        n_run = length(ctg_list_to_run);
        fprintf('n_run: %d', n_run);
    
        local_converged = false(1, n_run);
        local_failed_mask = false(1, n_run);
        
   
        % Step 1: DCOPF
        parfor i_valid = 1:n_run
            ctg_idx = ctg_list_to_run(i_valid);
            mpc_ctg = mpc_hour;
            mpc_ctg.gen(ctg_idx, PMAX) = 0;     % gen outage
            mpc_ctg.gen(ctg_idx, GEN_STATUS) = 0;

            try
                result = rundcopf(mpc_ctg, mpopt_dc); % no slack
                % fprintf('result.success:%d', result.success);
                if result.success
                    local_converged(i_valid) = true;
                else
                    local_failed_mask(i_valid) = true;
                end
            catch ME
                fprintf('DCOPF error in ctg %d: %s\n', ctg_idx, ME.message);
                local_failed_mask(i_valid) = true;
            end
        end
    
        % %Step 2: calculate load shedding amount, if necessary (for future use)
        % if sum(local_failed_mask) < n_gens*0.1 % only consider those buses whose failed contingency cases are lower than a threshold
        % 
        %     fprintf('-> Hour %d: %d failed ctgs, entering Step 2\n', selected_hour, sum(local_failed_mask));
        % 
        %     shed_tmp = nan(1, n_run);
        %     shed_flag_tmp = false(1, n_run);
        % 
        %     parfor i_valid = 1:n_run
        %         if ~local_converged(i_valid)
        %             ctg_idx = ctg_list_to_run(i_valid);
        %             mpc_ctg = mpc_hour;
        %             mpc_ctg.branch(ctg_idx, :) = [];
        %             ng_orig = size(mpc_hour.gen, 1);
        %             mpc_ctg = add_load_shedding_slack(mpc_ctg);
        % 
        %             try
        %                 result = rundcopf(mpc_ctg, mpopt_dc);
        %                 if result.success
        %                     shed = sum(result.gen(ng_orig+1:end, PG));
        %                     shed_tmp(i_valid) = shed;
        %                     shed_flag_tmp(i_valid) = shed > 1e-3;
        %                     % fprintf('shed:%f\n', shed);
        %                 end
        %             catch ME
        %                 fprintf('DCOPF error in ctg %d: %s\n', ctg_idx, ME.message);
        %             end
        %         end
        %     end
        %     shed_values(to_run_mask, hr_idx) = shed_tmp;
        %     shed_flags(to_run_mask, hr_idx) = shed_flag_tmp;
        % else
        %     fprintf('-> Hour %d: too many failed ctgs (%d), skipping Step 2\n', selected_hour, sum(local_failed_mask));
        % end
        % 
        % converged_ctg_mask(to_run_mask) = local_converged;
        % 
        % shed_count = sum(shed_flags(:));
        % avg_shed = mean(shed_values(shed_flags), 'omitnan');

        % Save this cap result
        results_gen_per_bus{j}{cap_idx} = struct( ...
            'dc_converged_count', sum(converged_ctg_mask), ...
            'shed_count', shed_count, ...
            'shed_avg', avg_shed);

        fprintf('===== Completed N-1 OPF Contingency Analysis =====\n');
    end
end



%N-k contingencies, according to user need
%Users can set the N-k contingencies of concerns and similarly calculate
%the N-k reliability results





% Create summary
summary_table = table();
summary_table.Bus = extra_load_buses_list(:);

for cap_i = 1:length(capacity_list)
    conv_rate = zeros(length(extra_load_buses_list), 1);

    for j_bus = 1:length(extra_load_buses_list)
        if isempty(result_per_bus{j_bus}) || isempty(result_per_bus{j_bus}{cap_i})
            conv_rate(j_bus) = 1;
        else
            res = result_per_bus{j_bus}{cap_i};
            conv_rate(j_bus) = dc_converged_count/n_valid_ctgs;
        end
    end
    summary_table.(sprintf('%dMW', capacity_list(cap_i))) = round(conv_rate, 3);
end

writetable(summary_table, output_csv);
fprintf('\nSaved result to %s\n', output_csv);


elapsed_time = toc(t_start);
fprintf('Elapsed time: %.4f seconds\n', elapsed_time);



%related functions
%%assign the gen buses to mpc gen buses
function assigned_idx = assign_gen_P(mpc, bus_list)
    gen_buses = mpc.gen(:, 1);
    gen_buses_available = gen_buses;  % Copy to track availability
    assigned_idx = zeros(length(bus_list), 1);

    for i = 1:length(bus_list)
        bus_i = bus_list(i);
        candidate_idxs = find(gen_buses_available == bus_i);

        if isempty(candidate_idxs)
            error('No available generator in mpc.gen for bus %d (input index %d)', bus_i, i);
        end

        assigned_idx(i) = candidate_idxs(1);
        gen_buses_available(candidate_idxs(1)) = NaN;  % Mark this one as used
    end
    % fprintf('Number of gens: %d\n', length(assigned_idx));
end

%calculate load shedding by adding slack generators to load buses
function mpc = add_load_shedding_slack(mpc)
    define_constants;

    % load_buses = mpc.bus(mpc.bus(:, PD) > 0, BUS_I);
    % existing_gen_buses = unique(mpc.gen(:, GEN_BUS));
    % buses_to_add = setdiff(load_buses, existing_gen_buses);
    % 
    % n_load_buses = length(buses_to_add);

    areas = unique(mpc.bus(:, BUS_AREA));
    n_area = length(areas);
    n_top = 120;
    n_slack = n_top*n_area;
    % fprintf('areas: %d', areas);
    gen_cols = size(mpc.gen, 2);
    n_pts = 5;

    x_points = linspace(0, 1, n_pts) * 200;       % 0, 25, 50, 75, 100 MW
    y_points = 1000000 * x_points;           %set a very high cost for slack gens

    shed_gens = [];
    shed_costs = {};
    shed_fuels = {};

    for i = 1:n_area
        area_i = areas(i);

        bus_idx = find(mpc.bus(:, BUS_AREA) == area_i & mpc.bus(:, PD) > 0);
        % fprintf('load bus in area %d: %d\n', area_i, length(bus_idx));

        [~, sorted_idx] = sort(mpc.bus(bus_idx, PD), 'descend');
        top_buses = bus_idx(sorted_idx(1:min(n_top, length(sorted_idx))));
        % fprintf('top_buses in area %d: %d\n', area_i, length(top_buses));

        for j = 1:length(top_buses)            
            bus_num = mpc.bus(top_buses(j), BUS_I);
            Pd = mpc.bus(top_buses(j), PD);
            % fprintf('bus id selected in area %d: %d\n', area_i, bus_num)
    
            % slack gen
            gen_row = zeros(1, gen_cols);
            gen_row(GEN_BUS) = bus_num;
            gen_row(PG) = 0;
            gen_row(QG) = 0;
            gen_row(QMAX) = 0;
            gen_row(QMIN) = 0;
            gen_row(VG) = 1.0;
            gen_row(MBASE) = mpc.baseMVA;
            gen_row(GEN_STATUS) = 1;
            gen_row(PMAX) = Pd;
            gen_row(PMIN) = 0;
    
            %[MODEL STARTUP SHUTDOWN N X1 Y1 X2 Y2 ... XN YN]
            cost_vec = zeros(1, 4 + 2 * n_pts);
            cost_vec(1) = 1;        % MODEL=1
            cost_vec(2) = 0;        % STARTUP
            cost_vec(3) = 0;        % SHUTDOWN
            cost_vec(4) = n_pts;    % number of points
            for k = 1:n_pts
                cost_vec(4 + 2*k - 1) = x_points(k);  % P
                cost_vec(4 + 2*k) = y_points(k);      % Cost
            end
            shed_gens = [shed_gens; gen_row];
            shed_costs = [shed_costs; cost_vec];
            shed_fuels = [shed_fuels; {'LS'}];
        end
    end

    % disp('slack gens:\n');
    % disp(slack_gens);
    mpc.gen = [mpc.gen; shed_gens];
    mpc.gencost = [mpc.gencost; cell2mat(shed_costs)];
    mpc.genfuel = [mpc.genfuel; shed_fuels];
end