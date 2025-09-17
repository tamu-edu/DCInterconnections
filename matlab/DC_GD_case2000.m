function xgd_table = DC_GD_case2000(mpc)

% get min-on (MinUp) time and min-off (MinDown) time
% in number of periods (this case delta_t = 1hour)
define_constants;
load('texas_2020_generation.mat'); % get min_on, min_off
ng = size(min_on,1);

commitkey = ones(ng,1);
commitkey(strcmp(mpc.genfuel,'nuclear')) = 2;

xgd_table.colnames = {'CommitKey','CommitSched','MinUp','MinDown'};
xgd_table.data = [commitkey, ones(ng,1), min_on, min_off];

end