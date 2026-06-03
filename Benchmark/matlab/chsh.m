%% Example: cvx_chsh.m 
% Solves CHSH scenario, with CVX.
% Expected answer: 2 sqrt 2 (2.828...)

%% Define scenario
% Two parties
tic
scenario = LocalityScenario(2);

Alice = scenario.Parties(1);
Bob = scenario.Parties(2);

% Each party with two measurements
A0 = Alice.AddMeasurement(2);
A1 = Alice.AddMeasurement(2);
B0 = Bob.AddMeasurement(2);
B1 = Bob.AddMeasurement(2);

%% Make moment matrix
matrix = scenario.MomentMatrix(10);


% Alternatively, make via full-correlator
CHSH_ineq = scenario.FCTensor([[0 0 0]; [0 1 1]; [0 1 -1]]);


%% Define and solve SDP
cvx_solver mosek
cvx_begin sdp

    % Declare basis variables a (real) and b (imaginary)
    scenario.cvxVars('a');
    
    % Compose moment matrix from these basis variables
    M = matrix.Apply(a);

    % Normalization
    a(1) == 1;

    % Positivity
    M >= 0;

    %CHSH inequality (maximize!)
    solve_chsh_ineq = CHSH_ineq.Apply(a);
    maximize(solve_chsh_ineq);
cvx_end
toc
%% Print out values found (should be identical!)
%chsh_max_val = CHSH_ineq.Apply(a)
