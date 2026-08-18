clear;
clc;
close all;

%% ============================================================
%  Panel A: Quantitative verification of the n^{-H} convergence rate
%
%  Model:
%       b(x) = 1
%       sigma(x) = 1
%
%  Sampling intensity:
%       theta_t = 1 + alpha |W_t|
%
%  We compute the pathwise error
%
%       |X_1^n - X_1|
%
%  and estimate
%
%       E[|X_1^n - X_1|]
%
%  by Monte Carlo.
%
%  The exact value of X_1 is simulated on the SAME Brownian path
%  using the exact Gaussian integral on each random sampling interval.
%
%  The theoretical prediction is
%
%       E[|X_1^n - X_1|] = O(n^{-H}).
%
%  Therefore the fitted log-log slope should be approximately -H.
%% ============================================================


%% Parameters

H_list = [0.1, 0.3, 0.4];

% Discretization levels
n_list = [3, 5, 8, 20, 50, 100, 200, 500];

% Monte Carlo sample size
M = 20000;

% Parameter in theta_t = 1 + alpha |W_t|
alpha = 0.5;


%% Store results

mean_error = zeros(length(H_list), length(n_list));
std_error  = zeros(length(H_list), length(n_list));
slope_list = zeros(length(H_list),1);


%% Figure

figure;
hold on;


%% Main simulation

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

            Xn = 0;      % Euler approximation
            Xexact = 0;  % Exact X on the same Brownian path


            % ----------------------------------------------------
            % Generate the irregular sampling intervals
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
                % Kernel at the left endpoint
                % ------------------------------------------------

                Ktau = (1 - tau)^(H - 0.5) / GammaH;


                % ------------------------------------------------
                % Euler scheme
                % ------------------------------------------------

                Xn = Xn + Ktau * Delta_tau + Ktau * dW;


                % ------------------------------------------------
                % Exact stochastic integral on this interval
                %
                % J = integral_tau^{tau_next} K(1-s) dW_s
                %
                % (J, dW) are jointly Gaussian.
                % ------------------------------------------------

                a = tau;
                b = tau_next;
                dt = Delta_tau;


                % Integral of K(1-s) ds
                A = ( ...
                    (1-a)^(H+0.5) - (1-b)^(H+0.5) ...
                    ) / ((H+0.5) * GammaH);


                % Integral of K(1-s)^2 ds
                Q = ( ...
                    (1-a)^(2*H) - (1-b)^(2*H) ...
                    ) / (2*H * GammaH^2);


                % Cov(J, dW) = A
                %
                % Var(dW) = dt
                %
                % Var(J) = Q
                %
                % Therefore:
                %
                % J = (A/dt) dW + sqrt(Q-A^2/dt) Z


                conditional_variance = Q - A^2 / dt;

                % Numerical protection against tiny negative values
                conditional_variance = max(conditional_variance, 0);


                Z = randn;

                J = (A/dt) * dW ...
                    + sqrt(conditional_variance) * Z;


                % Exact contribution to X
                Xexact = Xexact + A + J;


                % Update Brownian motion and time
                Wtau = Wtau + dW;
                tau = tau_next;

            end


            % ----------------------------------------------------
            % Pathwise discretization error
            % ----------------------------------------------------

            error_MC(m) = abs(Xn - Xexact);

        end


        %% Monte Carlo estimate

        mean_error(h,k) = mean(error_MC);

        std_error(h,k) = std(error_MC) / sqrt(M);

    end


    %% ------------------------------------------------------------
    % Log-log regression
    %
    % log(E|X^n-X|) = intercept + slope * log(n)
    %
    % Theoretical slope = -H
    %% ------------------------------------------------------------

    x = log(n_list(:));
    y = log(mean_error(h,:)');

    regression_coefficients = polyfit(x,y,1);

    slope = regression_coefficients(1);

    slope_list(h) = slope;


    %% ------------------------------------------------------------
    % Plot
    %% ------------------------------------------------------------

    loglog(n_list, mean_error(h,:), 'o-', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', ...
        ['H = ', num2str(H), ...
         ', fitted slope = ', num2str(slope,'%.3f')]);

end


%% Reference slopes

% We normalize the first point of each reference line
% to make it visually comparable with the Monte Carlo curve.

for h = 1:length(H_list)

    H = H_list(h);

    reference = mean_error(h,1) ...
        * (n_list / n_list(1)).^(-H);

    loglog(n_list, reference, '--', ...
        'LineWidth', 1.0, ...
        'HandleVisibility', 'off');

end


%% Figure formatting

xlabel('$logn$', 'Interpreter', 'latex');
ylabel('$logE[|X_1^n-X_1|]$', ...
    'Interpreter', 'latex');

title('Convergence rate of the irregular Euler scheme');

legend('Location','southwest');

grid on;
box on;


%% Print fitted slopes

fprintf('\n');
fprintf('===============================================\n');
fprintf('       Panel A: Fitted convergence rates\n');
fprintf('===============================================\n');
fprintf('   H        Fitted slope       Theoretical slope\n');
fprintf('-----------------------------------------------\n');

for h = 1:length(H_list)

    fprintf(' %.1f          %.4f              %.4f\n', ...
        H_list(h), ...
        slope_list(h), ...
        -H_list(h));

end

fprintf('===============================================\n');
fprintf('\n');
