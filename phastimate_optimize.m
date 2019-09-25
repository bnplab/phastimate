function [optimal_parameters ga_output] = phastimate_optimize(epochs, truephase, filter_objects_by_order, bounds_filter_order, bounds_window, bounds_edge, bounds_ar_order, hilbertwindow)
%PHASTIMATE_OPTIMIZE run genetic algorithm to optimize parameters for given set of epochs
%   optimal_parameters = phastimate_optimize(epochs, truephase, filter_objects)
%
%   Input:
%     data is a time x epoch matrix
%     truephase is the true phase in radians at the end of the epoch
%     filter_objects_by_order is a cell array of digitalFilter objects indexed by order
%     bounds_NN is lower and upper bound, e.g. [160 500]
%
%   Out: 
%     optimal_parameters : parameters for the winning run

%TODO: make population_size a parameter

%TODO: check that hilbert window is a power of 2

assert(~isempty(which('ga')), 'genetic algorithm function ga.m not found, is the Global Optimization Toolbox installed?')

assert(numel(truephase) == size(epochs,2), 'number of elements in truephase vector needs to match number of epochs')
assert(iscell(filter_objects_by_order), 'filter_objects_by_order must be a cell array of digitalFilter objects')
cellfun(@(x) assert(isa(x, 'digitalFilter'), 'filter_objects_by_order cell array members within bounds_filter_order must be of type digitalFilter'), filter_objects_by_order(bounds_filter_order(1):bounds_filter_order(end)))

problem.options = optimoptions('ga', 'PlotFcn', @gaplotbestf, 'PopulationSize', 100);
problem.solver = 'ga';

%PHASTIMATE(data, D, edge, ord, hilbertwindow, [offset_correction], [iterations], [armethod])

ang_var_of_diff = @(x, y) 1-abs(mean(exp(1i*x)./exp(1i*y)));

problem.fitnessfcn = @(x) ang_var_of_diff(truephase, ...
    phastimate(epochs((end-x(1)+1):end,:), filter_objects_by_order{x(2)}, x(3), x(4), hilbertwindow) ...
    );
% x(1) = window_length
% x(2) = filter_order
% x(3) = edge
% x(4) = ar_order

problem.nvars = 4;
problem.intcon = 1:4;
problem.lb = [bounds_window(1), bounds_filter_order(1), bounds_edge(1), bounds_ar_order(1)];
problem.ub = [bounds_window(2), bounds_filter_order(2), bounds_edge(2), bounds_ar_order(2)];

problem.x0 = ceil(mean([problem.lb; problem.ub]));

% A*x ≤ b
problem.Aineq = [-1 3 0 0]; %window_length x(1) > 3 * filter_order x(2)
problem.bineq = -1;

% check if the function evaluates
feval(problem.fitnessfcn, problem.x0);

% run solver
[x,fval,exitflag,ga_output] = ga(problem);

optimal_parameters = [];
optimal_parameters.window_length = x(1);
optimal_parameters.filter_order= x(2);
optimal_parameters.edge = x(3);
optimal_parameters.ar_order = x(4);
optimal_parameters.fval = fval;
    
end