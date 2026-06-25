function [setup_time, solve_time, obj_val] = polygon_bell(n, mm_level, gamma)
    if nargin < 3
        gamma = [];
    end
    tic;
    % One operator per vertex of the polygon: vars_A = ["X1","X2",...,"Xn"]
    vars_A = arrayfun(@(i) sprintf("X%d", i), 1:n);
    % Algebraic scenario with n hermitian operators
    setting = AlgebraicScenario(vars_A);

    rules = setting.OperatorRulebook;

    % next add projector constraints
    for i = 1:length(vars_A)
        rules.MakeProjector(vars_A(i));
    end
    % next add commutation constraints between adjacent vertices:
    % X_{i} commutes with X_{i+1} (cyclic)
    for k = 1:n
        next_k = mod(k, n) + 1;
        rules.AddCommutator(vars_A(k), vars_A(next_k));
    end
    setting.Complete(100);
    % ops(i) is the operator for vertex i
    ops = setting.getAll();

    % Objective function (projector form: dichotomic = 1 - 2*projector)
    obj = 0;
    for i = 1:n
        next_i = mod(i, n) + 1;
        if ~isempty(gamma)
            c = gamma(i);
        elseif mod(i, 4) == 0
            c = -1;
        else
            c = 1;
        end
        obj = obj + c*(1-2*ops(i))*(1-2*ops(next_i));
    end

    mm = setting.MomentMatrix(mm_level);
    yalmip('clear');

    % Declare basis variables a (real)
    model_vars = setting.yalmipVars();

    % Compose moment matrix in these basis variables
    M = mm.Apply(model_vars);

    constraints = [model_vars(1) == 1];
    constraints = [constraints, M>=0];

    % Objective function (maximize)
    objective = obj.Apply(model_vars);
    setup_time=toc;

    tic
    % Solve
    optimize(constraints, objective);
    solve_time =  toc;

    obj_val = value(objective);
end

function [avgsetuptime,avgsolvetime] = avg_time(t_runs,n,mm_level)
    total_setup_time = 0;
    total_solve_time = 0;
    for run = 1:t_runs
        [setup_time,solve_time,val] = polygon_bell(n,mm_level);
        total_setup_time = total_setup_time + setup_time;
        total_solve_time = total_solve_time + solve_time;
    end
    avgsetuptime = total_setup_time/t_runs;
    avgsolvetime = total_solve_time/t_runs;
end
%[avgsetuptime,avgsolvetime] = avg_time(10,5,2);
%fprintf('Average setup time: %.4f seconds\n', avgsetuptime);
%fprintf('Average solve time: %.4f seconds\n', avgsolvetime);
[sut,sot,val] = polygon_bell(5,2)