% Quick dependency checks for MATPOWER/MOST/YALMIP/GUROBI
fprintf('== Checking key dependencies on MATLAB path ==\n');

% MATPOWER / MOST
req = {'define_constants','mpoption','most','loadcase','loadxgendata'};
for k = 1:numel(req)
    w = which(req{k},'-all');
    if isempty(w)
        fprintf(2,'[MISSING] %s not found on path\n', req{k});
    else
        fprintf('[OK] %s -> %s\n', req{k}, w{1});
    end
end

% GUROBI MEX
w = which('gurobi','-all');
if isempty(w)
    fprintf(2,'[MISSING] Gurobi MATLAB interface (gurobi.mexa64) not found\n');
else
    fprintf('[OK] Gurobi MATLAB interface -> %s\n', w{1});
end

% YALMIP
w = which('sdpsettings','-all');
if isempty(w)
    fprintf(2,'[MISSING] YALMIP (sdpsettings) not found\n');
else
    fprintf('[OK] YALMIP -> %s\n', w{1});
end

% MATPOWER have_fcn
w = which('have_fcn','-all');
if isempty(w)
    fprintf(2,'[WARN] have_fcn not found (MATPOWER lib)\n');
else
    fprintf('[OK] have_fcn -> %s\n', w{1});
end

fprintf('== Done ==\n');
