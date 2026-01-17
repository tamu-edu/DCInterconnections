%% run_SCUC_1GW_sensitivity_HPRC_resume.m
% Robust, resume-able SCUC/SCED driver for HPRC with atomic bus–day-level checkpoints.
% - If a (bus,day) file exists and is valid, it is skipped.
% - Crashes/timeouts only require rerunning the remaining (bus,day) pairs.
% - Aggregation only runs once all (bus,day) pairs are valid for that bus.

clearvars; clc;

%% 0) Parallel pool / scratch setup
RESUME_ONLY = false;   
AGG_ONLY    = false;  

numWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(numWorkers) || numWorkers < 1
    numWorkers = max(1, feature('numcores')-1);
end

scratchDir  = '/proactive_siting/';
jobDir = fullfile(scratchDir, 'matlab_jobstorage', getenv('SLURM_JOB_ID'));
if ~exist(jobDir,'dir'); mkdir(jobDir); end

needPool = ~(AGG_ONLY || RESUME_ONLY);
if needPool
    c = parcluster('Processes');
    c.JobStorageLocation = jobDir;
    c.NumWorkers = numWorkers;
    if isempty(gcp('nocreate')); parpool(c, numWorkers); end
end

%% 1) Constants & result root
define_constants;
lam_p  = LAM_P; mu_sf = MU_SF; mu_st = MU_ST;

result_root = fullfile(scratchDir, 'results');
if ~exist(result_root,'dir'); mkdir(result_root); end

%% 2) Load MATPOWER & XGD cases
mpc_base = loadcase('DC_case2000');
xgd_base = loadxgendata('DC_GD_case2000', mpc_base);

%% 3) Load time-series data
D = load('/data/texas_2020_demand.mat');   area_load_arr = D.area_load;   clear D
W = load('/data/texas_2020_wind.mat');     wind_arr      = W.wind_MW;     clear W
S = load('/data/texas_2020_solar.mat');    solar_arr     = S.solar_MW;    clear S
H = load('/data/texas_2020_hydro.mat');    hydro_arr     = H.hydro_MW;    clear H

%% 4) MOST options
mpopt = mpoption( ...
    'out.all',0, 'verbose',3, ...
    'gurobi.threads',1, ...
    'gurobi.opts.TimeLimit',3600, ...
    'gurobi.opts.MIPGap',1e-2, ...
    'most.dc_model',1, ...
    'most.skip_prices',0, ...
    'most.uc.run',1);

%% 5) Days and buses
start_day = 226;
end_day   = 226;
days      = start_day:end_day;
nt        = 24;

load('reliability.mat','reliability');
testBuses = reliability;
nB = numel(testBuses);

%% 6) Scenario setup
flexMode  = 'base';
sizeDC    = 1000;  
tol_mu    = 1e-6;
peakHours = 16:19;
scenario_tag = sprintf('%s_%dMW', flexMode, sizeDC);

%% 7) Pre-create per-bus directories
busRoots = cell(nB,1); outDirs = cell(nB,1);
for jj = 1:nB
    b = testBuses(jj);
    busRoots{jj} = fullfile(result_root, sprintf('bus%04d_%s', b, scenario_tag));
    outDirs{jj}  = fullfile(busRoots{jj}, 'scuc_sced');
    if ~exist(busRoots{jj},'dir'), mkdir(busRoots{jj}); end
    if ~exist(outDirs{jj}, 'dir'), mkdir(outDirs{jj});  end
end

%% 8) Flatten all (bus,day) pairs
pairs = [];
for ii = 1:nB
    for d = days
        pairs = [pairs; ii, d]; %#ok<AGROW>
    end
end
nPairs = size(pairs,1);

%% 9) parfor over flattened pairs
nbus = size(mpc_base.bus,1);
nln  = size(mpc_base.branch,1);

parfor pp = 1:nPairs
    if AGG_ONLY, continue; end

    ii = pairs(pp,1);
    d  = pairs(pp,2);
    b  = testBuses(ii);

    this_root = busRoots{ii};
    out_dir   = outDirs{ii};
    fn        = fullfile(out_dir, sprintf('24h-SCUC-day-%d.mat', d));

    % skip if already valid
    if is_valid_dayfile(fn, nbus, nln), continue; end

    % build modified load (fixed for this bus)
    dc_24h   = build_dc_profile(flexMode, sizeDC);
    area_mod = area_load_arr;
    area_mod(:,b) = area_mod(:,b) + repmat(dc_24h,365,1);

    % run SCUC-SCED
    try
        [ok,lambda_bus,mu_branch,flow_branch,Fmax] = ...
            try_run_SCUC_day(d, mpc_base, xgd_base, ...
                              area_mod, wind_arr, solar_arr, hydro_arr, ...
                              mpopt, lam_p, mu_sf, mu_st);
    catch ME
        warning('Day %d at bus %d crashed: %s', d, b, ME.message);
        ok = false;
    end
    if ok
        S = struct('lambda_bus',lambda_bus, ...
                   'mu_branch',mu_branch, ...
                   'flow_branch',flow_branch, ...
                   'Fmax',Fmax);
        atomic_save(fn, S);
    end
end

%% 10) Tear down pool
if needPool, delete(gcp('nocreate')); end

%% 11) Aggregation (per bus)
mkPM = @(tag,mA,rA,sA,mP,rP,sP,mO,rO,sO) struct( ...
    'scenario',        char(tag), ...
    'meanLMP_all',     mA, 'range95_5_all', rA, 'stdLMP_all',  sA, ...
    'meanLMP_peak',    mP, 'range95_5_peak',rP, 'stdLMP_peak', sP, ...
    'meanLMP_off',     mO, 'range95_5_off', rO, 'stdLMP_off',  sO);

% UPDATED: include legacy metrics too (bindingHours_* and sumAbsShadow_*)
mkCM = @(tag,bA,cA,bP,cP,bO,cO, ...
            bHa,sAa,bHp,sAp,bHo,sAo) struct( ...
    'scenario',              char(tag), ...
    'bindingShare_all',      bA,   'congRent_all_Mpd',   cA, ...
    'bindingShare_peak',     bP,   'congRent_peak_Mpd',  cP, ...
    'bindingShare_off',      bO,   'congRent_off_Mpd',   cO, ...
    'bindingHours_all',      bHa,  'sumAbsShadow_all',   sAa, ...
    'bindingHours_peak',     bHp,  'sumAbsShadow_peak',  sAp, ...
    'bindingHours_off',      bHo,  'sumAbsShadow_off',   sAo );

for ii = 1:nB
    this_root = busRoots{ii};
    out_dir   = outDirs{ii};
    done_flag = fullfile(this_root, 'DONE.txt');
    if exist(done_flag,'file') && exist(fullfile(this_root,'metrics_summary.mat'),'file')
        continue;
    end

    % check all days valid
    all_ok = true;
    for d = days
        fn = fullfile(out_dir, sprintf('24h-SCUC-day-%d.mat', d));
        if ~is_valid_dayfile(fn, nbus, nln)
            all_ok = false; break;
        end
    end
    if ~all_ok, continue; end

    % aggregation
    LMP=[]; MU=[]; FLOW=[]; Fmax=[];
    for d = days
        M = load(fullfile(out_dir, sprintf('24h-SCUC-day-%d.mat', d)), ...
                 'lambda_bus','mu_branch','flow_branch','Fmax');
        LMP=[LMP; M.lambda_bus];
        MU=[MU; M.mu_branch];
        FLOW=[FLOW; M.flow_branch];
        if isempty(Fmax), Fmax = M.Fmax(:)'; end
    end

    T=size(LMP,1); hod=mod((0:T-1),24)+1;
    idx_peak=ismember(hod,peakHours); idx_off=~idx_peak;

    nodeMean=mean(LMP,2); nodeStd=std(LMP,0,2);
    p95=prctile(LMP,95,2); p05=prctile(LMP,5,2);
    nodeRange=p95-p05; bind_any=any(MU>tol_mu,2);
    rent_hour=sum(MU.*abs(FLOW),2);

    [mA,rA,sA,bA,cA]=reduce_window(true(T,1),nodeMean,nodeRange,nodeStd,bind_any,rent_hour,numel(days));
    [mP,rP,sP,bP,cP]=reduce_window(idx_peak,nodeMean,nodeRange,nodeStd,bind_any,rent_hour,numel(days));
    [mO,rO,sO,bO,cO]=reduce_window(idx_off,nodeMean,nodeRange,nodeStd,bind_any,rent_hour,numel(days));

    % === ADD: legacy-style congestion metrics (match old script semantics) ===
    % line×hour bindings counts and sum of absolute shadow prices
    bindingHours_all  = sum(MU(:) > tol_mu);                 % line-hour count
    sumAbsShadow_all  = sum(abs(MU(:)));                     % $/MW aggregated over lines×hours

    bindingHours_peak = sum(sum(MU(idx_peak,:) > tol_mu));   % line-hour count in peak window
    sumAbsShadow_peak = sum(sum(abs(MU(idx_peak,:))));

    bindingHours_off  = sum(sum(MU(idx_off,:)  > tol_mu));   % line-hour count off-peak
    sumAbsShadow_off  = sum(sum(abs(MU(idx_off,:))));

    b=testBuses(ii);
    tag=sprintf('bus%04d_%s', b, scenario_tag);
    pm=mkPM(tag,mA,rA,sA,mP,rP,sP,mO,rO,sO);
    cm=mkCM(tag,bA,cA,bP,cP,bO,cO, ...
             bindingHours_all, sumAbsShadow_all, ...
             bindingHours_peak, sumAbsShadow_peak, ...
             bindingHours_off,  sumAbsShadow_off);

    atomic_save(fullfile(this_root,'metrics_summary.mat'), struct('pm',pm,'cm',cm));
    fid=fopen(done_flag,'w'); if fid>0, fprintf(fid,'done\n'); fclose(fid); end
end

%% 12) Combine all buses
allPM={}; allCM={};
for ii=1:nB
    ms=fullfile(busRoots{ii},'metrics_summary.mat');
    if exist(ms,'file')
        S=load(ms,'pm','cm');
        allPM{end+1}=S.pm; %#ok<AGROW>
        allCM{end+1}=S.cm; %#ok<AGROW>
    end
end
if ~isempty(allPM)
    priceMetrics_struct=[allPM{:}];
    congMetrics_struct =[allCM{:}];
    atomic_save(fullfile(result_root, sprintf('all_metrics_%s.mat', scenario_tag)), ...
        struct('priceMetrics_struct',priceMetrics_struct,'congMetrics_struct',congMetrics_struct));
end

% ========================================================================
% HELPER FUNCTIONS (build_dc_profile, daily_SCUD_SCED, reduce_window,
% is_valid_dayfile, atomic_save, try_run_SCUC_day)
% ========================================================================

function p = build_dc_profile(mode,sizeDC)
% Returns a 24×1 MW profile; peak window = hours 16–19 (1-based)
nor = 0.8*sizeDC;          % “normal” 80 % load

switch lower(mode)
    case 'pause'  % 1–15 nor | 16–19 curtailed | 20–24 nor
        pausePCT = 0.90;
        p = [repmat(nor,15,1); repmat(pausePCT*nor,4,1); repmat(nor,5,1)];

    case 'shift'  % energy-neutral shift: 16–19 at 80% of nor, off-peak makeup
        shiftPCT = 0.80;                         % 80 % during peak
        curtailed = (1-shiftPCT)*nor*4;          % MWh curtailed over 4 h
        bump     = curtailed / 20;               % spread over 20 off-peak hours
        p = [repmat(nor + bump,15,1); ...        % 1–15
             repmat(shiftPCT*nor,4,1); ...       % 16–19
             repmat(nor + bump,5,1)];            % 20–24

    case 'base'   % constant 0.8 · sizeDC all day
        p = repmat(nor,24,1);

    otherwise
        error('unknown flex mode %s',mode);
end
end

function [mdo, ms] = daily_SCUD_SCED(mpc, xgd, area_load, wind_pen, solar_pen, hydro_pen, init, mpopt)

define_constants;

% ----- sizes -----
nb   = size(mpc.bus,1);
nl   = size(mpc.branch,1);
ng   = size(mpc.gen,1);
nt   = 24;          % hours
area_ind = (1:nb)'; 

% relax-bounds scaling (unchanged)
if isfield(mpopt,'relax_bounds')
  ub_scale = 1 + mpopt.relax_bounds;
  lb_scale = 1 - mpopt.relax_bounds;
else
  ub_scale = 1; lb_scale = 1;
end

% (optional) renumbering if needed — keep order stable with F/T mapping
persistent did_warn_reindex;
if isempty(did_warn_reindex)
    fprintf('Renumbering buses to consecutive IDs...\n');
    did_warn_reindex = true;
end
old_bus = mpc.bus(:,BUS_I);
new_bus = (1:nb)';

mpc_ordered = mpc;
mpc_ordered.bus(:,BUS_I) = new_bus;
for i = 1:nl
  f = mpc.branch(i,F_BUS);  t = mpc.branch(i,T_BUS);
  from_idx = find(old_bus == f, 1, 'first');
  to_idx   = find(old_bus == t, 1, 'first');
  mpc_ordered.branch(i,[F_BUS T_BUS]) = [new_bus(from_idx), new_bus(to_idx)];
end
for i = 1:ng
  g = mpc.gen(i,GEN_BUS);
  gi = find(old_bus == g, 1, 'first');
  mpc_ordered.gen(i,GEN_BUS) = new_bus(gi);
end

% find renewables (assumes genfuel exists)
iwind  = find(strcmp(mpc_ordered.genfuel,'wind'));
isolar = find(strcmp(mpc_ordered.genfuel,'solar'));
ihydro = find(strcmp(mpc_ordered.genfuel,'hydro'));
nwind  = numel(iwind); nsolar = numel(isolar); nhydro = numel(ihydro);

%% —— bus-level load profile ——
disp('creating bus-level load profiles');
busprofile = struct('type','mpcData','table',CT_TBUS,'rows',area_ind,'col',PD,'chgtype',CT_REP,'values',[]);
busprofile.values(:,1,:) = area_load;  % nt×1×nb
profiles = getprofiles(busprofile);

disp('creating wind profiles');
windprofile_max = struct('type','mpcData','table',CT_TGEN,'rows',(1:nwind)','col',PMAX,'chgtype',CT_REP,'values',[]);
windprofile_max.values(:,1,:) = wind_pen * ub_scale;
profiles = getprofiles(windprofile_max, profiles, iwind);

disp('creating solar profiles');
solarprofile_max = struct('type','mpcData','table',CT_TGEN,'rows',(1:nsolar)','col',PMAX,'chgtype',CT_REP,'values',[]);
solarprofile_max.values(:,1,:) = solar_pen * ub_scale;
profiles = getprofiles(solarprofile_max, profiles, isolar);

disp('creating hydro profiles');
hydroprofile_max = struct('type','mpcData','table',CT_TGEN,'rows',(1:nhydro)','col',PMAX,'chgtype',CT_REP,'values',[]);
hydroprofile_max.values(:,1,:) = hydro_pen * ub_scale;
profiles = getprofiles(hydroprofile_max, profiles, ihydro);

% build MOST model and run
mdi = loadmd(mpc_ordered, nt, xgd, [], [], profiles);
if ~isempty(init)
  mdi.UC.InitialState = init.commit;
  mdi.UC.InitialPg    = init.dispatch;
end

try
    mdo = most(mdi, mpopt);
    ms.exitflag = mdo.QP.exitflag;
catch
    ms.exitflag = -4;
    mdo        = struct();
    return
end
end

function [mMean, mRange, mStd, bShare, rentMpd] = reduce_window(mask, nodeMean, nodeRange, nodeStd, bind_any, rent_hour, nDays)
% mask: T×1 logical; rent_hour: $/h
H = max(1, sum(mask)); %#ok<NASGU>
mMean  = mean(nodeMean(mask));
mRange = mean(nodeRange(mask));
mStd   = mean(nodeStd(mask));
bShare = 100 * mean(bind_any(mask));                 % [% of hours]
rent_per_day = sum(rent_hour(mask)) / max(1, nDays); % $/day
rentMpd = rent_per_day / 1e6;                        % $M/day
end

function ok = is_valid_dayfile(fname, nbus, nline)

ok = false;
if ~exist(fname,'file'); return; end
try
    S = load(fname, 'lambda_bus','mu_branch','flow_branch','Fmax');
    ok = isfield(S,'lambda_bus') && isfield(S,'mu_branch') && ...
         isfield(S,'flow_branch') && isfield(S,'Fmax') && ...
         isequal(size(S.lambda_bus), [24, nbus]) && ...
         isequal(size(S.mu_branch),  [24, nline]) && ...
         isequal(size(S.flow_branch),[24, nline]) && ...
         (isequal(size(S.Fmax), [1, nline]) || isequal(size(S.Fmax), [nline, 1]) || numel(S.Fmax)==nline);
catch
    ok = false;
end
end

function atomic_save(final_path, vars_struct)
[tmpdir, tmpname, tmpext] = fileparts(final_path); %#ok<ASGLU>
if isempty(tmpdir), tmpdir = pwd; end
uuid = char(java.util.UUID.randomUUID);
tmpfile = fullfile(tmpdir, [tmpname, '.part_', uuid]);
save(tmpfile, '-struct', 'vars_struct', '-v7');  
[ok,msg] = movefile(tmpfile, final_path, 'f');
if ~ok
    warning('atomic_save movefile failed: %s', msg);
    if exist(tmpfile,'file'); delete(tmpfile); end
    error('atomic_save failed for %s', final_path);
end
end

function [ok,lambda_bus,mu_branch,flow_branch,Fmax] = try_run_SCUC_day( ...
    day, mpc_base, xgd_base, area_load, wind_arr, solar_arr, hydro_arr, ...
    mpopt, lam_p, mu_sf, mu_st)

ok = false;
lambda_bus=[]; mu_branch=[]; flow_branch=[]; Fmax=[];

define_constants;

nt   = 24;
base = (day-1)*nt;


al = area_load(base+(1:nt),:);
wp = wind_arr (base+(1:nt),:);
sp = solar_arr(base+(1:nt),:);
hp = hydro_arr(base+(1:nt),:);


[mdo, ms] = daily_SCUD_SCED(mpc_base, xgd_base, al, wp, sp, hp, [], mpopt);
if any(ms.exitflag <= 0), return; end

nbus = size(mpc_base.bus,1);
nln  = size(mpc_base.branch,1);

lambda_bus = zeros(24, nbus);
mu_branch  = zeros(24, nln);
flow_branch= zeros(24, nln);
Fmax       = mpc_base.branch(:, RATE_A)';

for t = 1:nt
    BT = mdo.flow(t).mpc.bus;
    BR = mdo.flow(t).mpc.branch;
    lambda_bus(t,:) = BT(:, lam_p).';                   % $/MWh per bus
    mu_sum          = max(BR(:, mu_sf),0) + max(BR(:, mu_st),0);  % $/MW
    mu_branch(t,:)  = mu_sum.';
    flow_branch(t,:)= abs(BR(:, PF)).';                 % MW (|PF|)
end
if all(lambda_bus(:)==0) && all(mu_branch(:)==0), return; end
ok = true;
end