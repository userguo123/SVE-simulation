clear;
clc;

%% ============================================================
% Panel A: Quantitative verification of the n^{-H} convergence rate
%
% We estimate
%
%       E_n = E[|X_1^n - X_1|]
%
% by Monte Carlo and fit
%
%       log(E_n) = C - H log(n).
%
% Therefore:
%
%       fitted slope = -H
%
% and
%
%       fitted H = - fitted slope.
%
% No figure is generated.
%% ============================================================


%% Parameters

H_list = [0.1, 0.2, 0.3, 0.4];

% Discretization levels
n_list = [3, 5, 8, 20, 50, 100, 200, 500];

% Monte Carlo sample size
M = 20000;

% Parameter in theta_t = 1 + alpha |W_t|
alpha = 0.5;


%% Store results

mean_error = zeros(length(H_list), length(n_list));
std_error  = zeros(length(H_list), length(n_list));

% log(n)
log_n = log(n_list);

% Fitted slope
slope_list = zeros(length(H_list),1);

% Fitted H = -slope
H_fitted = zeros(length(H_list),1);


%% ============================================================
% Main simulation
%% ============================================================

for h = 1:length(H_list)

    H = H_list(h);

    GammaH = gamma(H + 0.5);


    for k = 1:length(n_list)

        n = n_list(k);

        error_MC = zeros(M,1);


        %% Monte Carlo simulation

        for m = 1:M

            % ----------------------------------------------------
            % Initial values
            % ----------------------------------------------------

            tau = 0;
            Wtau = 0;

            Xn = 0;
            Xexact = 0;


            % ----------------------------------------------------
            % Generate irregular sampling intervals
            % ----------------------------------------------------

            while tau < 1

                % Sampling intensity
                theta = 1 + alpha * abs(Wtau);

                % Proposed time step
                Delta_tau = 1 / (n * theta);

                % Avoid overshooting t = 1
                if tau + Delta_tau > 1
                    Delta_tau = 1 - tau;
                end

                tau_next = tau + Delta_tau;


                % ------------------------------------------------
                % Brownian increment
                % ------------------------------------------------

                dW = sqrt(Delta_tau) * randn;


                % ------------------------------------------------
                % Kernel at left endpoint
                % ------------------------------------------------

                Ktau = (1 - tau)^(H - 0.5) / GammaH;


                % ------------------------------------------------
                % Euler approximation
                % ------------------------------------------------

                Xn = Xn ...
                    + Ktau * Delta_tau ...
                    + Ktau * dW;


                % ------------------------------------------------
                % Exact stochastic integral on this interval
                % ------------------------------------------------

                a = tau;
                b = tau_next;
                dt = Delta_tau;


                % Integral of K(1-s) ds
                A = ( ...
                    (1-a)^(H+0.5) ...
                    - (1-b)^(H+0.5) ...
                    ) / ((H+0.5) * GammaH);


                % Integral of K(1-s)^2 ds
                Q = ( ...
                    (1-a)^(2*H) ...
                    - (1-b)^(2*H) ...
                    ) / (2*H * GammaH^2);


                % Conditional variance
                conditional_variance = ...
                    Q - A^2 / dt;

                % Numerical protection against round-off error
                conditional_variance = ...
                    max(conditional_variance, 0);


                Z = randn;

                J = (A/dt) * dW ...
                    + sqrt(conditional_variance) * Z;


                % Exact contribution
                Xexact = Xexact + A + J;


                % Update Brownian motion and time
                Wtau = Wtau + dW;
                tau = tau_next;

            end


            % ----------------------------------------------------
            % Pathwise error
            % ----------------------------------------------------

            error_MC(m) = abs(Xn - Xexact);

        end


        %% Monte Carlo estimate

        mean_error(h,k) = mean(error_MC);

        std_error(h,k) = std(error_MC) / sqrt(M);

    end


    %% ============================================================
    % Log-log regression
    %
    %       log(E_n) = intercept + slope * log(n)
    %
    % Theoretical:
    %
    %       slope = -H
    %
    % Hence:
    %
    %       fitted H = -slope
    %% ============================================================

    x = log(n_list(:));

    y = log(mean_error(h,:)');

    regression_coefficients = polyfit(x, y, 1);

    slope = regression_coefficients(1);

    intercept = regression_coefficients(2);

    slope_list(h) = slope;

    H_fitted(h) = -slope;

end


%% ============================================================
% Output 1:
% log(n) and log(E_n) for every H and every n
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('             log(n) and log(E_n) data\n');
fprintf('============================================================\n');

for h = 1:length(H_list)

    fprintf('\n');
    fprintf('H = %.1f\n', H_list(h));
    fprintf('------------------------------------------------------------\n');
    fprintf('       n              log(n)              log(E_n)\n');
    fprintf('------------------------------------------------------------\n');

    for k = 1:length(n_list)

        fprintf('%8d        %12.6f        %15.6f\n', ...
            n_list(k), ...
            log_n(k), ...
            log(mean_error(h,k)));

    end

end


%% ============================================================
% Output 2:
% Summary table
%
% H
% Fitted H
% Fitted slope
% Theoretical H
% Theoretical slope
%% ============================================================

fprintf('\n\n');
fprintf('============================================================\n');
fprintf('             Panel A: Regression Results\n');
fprintf('============================================================\n');

fprintf(['   H       Fitted H       Fitted slope       ', ...
         'Theoretical H       Theoretical slope\n']);

fprintf('------------------------------------------------------------\n');

for h = 1:length(H_list)

    fprintf(' %.1f       %.4f          %.4f             %.1f                %.4f\n', ...
        H_list(h), ...
        H_fitted(h), ...
        slope_list(h), ...
        H_list(h), ...
        -H_list(h));

end

fprintf('============================================================\n');


%% ============================================================
% Optional: create MATLAB tables
%
% These tables can also be exported to CSV if needed.
%% ============================================================


%% Detailed data table

H_column = [];
n_column = [];
logn_column = [];
logE_column = [];

for h = 1:length(H_list)

    H_column = [H_column; ...
                repmat(H_list(h), length(n_list), 1)];

    n_column = [n_column; ...
                n_list(:)];

    logn_column = [logn_column; ...
                   log(n_list(:))];

    logE_column = [logE_column; ...
                   log(mean_error(h,:)')];

end

data_table = table( ...
    H_column, ...
    n_column, ...
    logn_column, ...
    logE_column, ...
    'VariableNames', ...
    {'H','n','log_n','log_E'});


%% Summary table

summary_table = table( ...
    H_list(:), ...
    H_fitted, ...
    slope_list, ...
    H_list(:), ...
    -H_list(:), ...
    'VariableNames', ...
    {'H','Fitted_H','Fitted_Slope', ...
     'Theoretical_H','Theoretical_Slope'});


%% Display MATLAB tables

fprintf('\n');
fprintf('Detailed data table:\n');
disp(data_table);

fprintf('\n');
fprintf('Summary table:\n');
disp(summary_table);


%% ============================================================
% Save tables
%
% Uncomment the following two lines if you want CSV files.
%% ============================================================

% writetable(data_table, 'PanelA_log_data.csv');
% writetable(summary_table, 'PanelA_summary.csv');
