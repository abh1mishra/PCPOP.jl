[opt_val] = gmnl(4,4);

function [opt_val] = gmnl(n, k)
    % Initialize Locality Scenario with n parties, 2 measurements, 2 outcomes
    tic
    scenario = LocalityScenario(n, 2, 2);
    
    % Generate operator objects
    % In LocalityScenario, operators are automatically projectors and commute 
    % between different parties. With 2 outcomes, it generates 1 operator per 
    % measurement, yielding 2 operators per party (total 2*n).
    A_flat = scenario.getAll();
    
    A = cell(n, 2);
    for i = 1:n
        A{i, 1} = A_flat(2*(i-1) + 1);
        A{i, 2} = A_flat(2*(i-1) + 2);
    end
    
    I_id = scenario.id();
    
    % Objective: Product of A{i,1}
    obj = A{1, 1};
    for i = 2:n
        obj = obj * A{i, 1};
    end
    
    % Equality constraints on moments
    eq_list = cell(n+1, 1);
    for i = 1:(n-1)
        eq_list{i} = A{i, 2} * A{i+1, 1};
    end
    eq_list{n} = A{n, 2} * A{1, 1};
    
    prod_term = I_id - A{1, 2};
    for i = 2:n
        prod_term = prod_term * (I_id - A{i, 2});
    end
    eq_list{n+1} = prod_term;
        
    % Create Moment Matrix
    mm = scenario.MomentMatrix(k);
    %cvx_solver mosek
    cvx_begin sdp

        scenario.cvxVars('a');
        M = mm.Apply(a);
        objective = obj.Apply(a);
        
        a(1) == 1;
        M >= 0;
        
        % Apply moment equalities separately as YALMIP constraints
        for i = 1:length(eq_list)
             eq_list{i}.Apply(a) == 0;
        end
        
        maximize(objective)
        
        opt_val = value(objective);
    cvx_end
    toc
end
