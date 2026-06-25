function [setup_time, solve_time, obj_val] = polygon_bell(n, mm_level)
    tic;
    % to generate the matrix of Operatornames, where each row corresponds to names for variables in a vertex of polygon 
    vars_M = arrayfun(@(i,j) sprintf("X%d%d", i, j), repmat((1:n)', 1, 2), repmat(1:2, n, 1), 'UniformOutput', false);
    % Below we transpose and then flatten the matrix vars_M to get vars_A,
    % thsi transpose and flatten helps set the order of variables that
    % seems to give complete rules atleast for n=3.
    vars_MT = vars_M';
    vars_A = [vars_MT{:}];
    % Algebraic scenario with 2*n hermitian operators
    setting = AlgebraicScenario(vars_A);

    rules = setting.OperatorRulebook;

    % next add projector constraints
    for i = 1:length(vars_A)
        rules.MakeProjector(vars_A(i));
    end
    % next add commutation constraints between adjacent vertices
    for k = 1:n
        next_k = mod(k, n) + 1;
        for i = 1:2
            for j = 1:2
                rules.AddCommutator(vars_M{k,i}, vars_M{next_k,j});
            end
        end
    end
    setting.Complete(100);
    % vars_A is row-major: ops(2k-1)=row k col1, ops(2k)=row k col2
    % Here it gives warning of "Supplied ruleset was not completed."
    ops = setting.getAll();

    % Cyclic CHSH objective: CHSH between each pair of consecutive rows
    obj = 0;
    for k = 1:n
        next_k = mod(k, n) + 1;
        a1 = ops(2*k-1); a2 = ops(2*k);
        b1 = ops(2*next_k-1); b2 = ops(2*next_k);
        obj = obj + (1-2*a1)*(1-2*b1) + (1-2*a1)*(1-2*b2) + (1-2*a2)*(1-2*b1) - (1-2*a2)*(1-2*b2);
    end

    % Here it gives error of some Symbol X... not found in Symbols table.
    mm = setting.MomentMatrix(mm_level);
    yalmip('clear');

    % Declare basis variables a (real)
    model_vars = setting.yalmipVars();
    
    % Compose moment matrix in these basis variables
    M = mm.Apply(model_vars);

    constraints = [model_vars(1) == 1];
    constraints = [constraints, M>=0];

    % Objective function (maximize)
    objective = -obj.Apply(model_vars);
    options = sdpsettings('verbose',0);
    setup_time=toc;

    tic
    % Solve
    optimize(constraints, objective,options);
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
[setup_time, solve_time, obj_val] = polygon_bell(7, 2)