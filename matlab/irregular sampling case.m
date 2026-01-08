
clear; clc;

H_list = [0.1, 0.3, 0.4];
n_list = [3, 5, 8, 20, 100];
M = 70000;          % Monte Carlo size
alpha = 0.5;        % parameter in theta_t = 1 + alpha |W_t|

figure; hold on;

for h = 1:length(H_list)
    H = H_list(h);
    GammaH = gamma(H + 0.5);

    % Exact solution
    EX = 1 / (GammaH * (H + 0.5));

    err = zeros(length(n_list),1);

    for k = 1:length(n_list)
        n = n_list(k);

        Xn_MC = zeros(M,1);

        for m = 1:M
            % Initialize
            tau = 0;
            Wtau = 0;
            Xn = 0;

            % Generate random sampling times
            while tau < 1
                % theta(tau)
                theta = 1 + alpha * abs(Wtau);

                % random time step
                Delta_tau = 1 / (n * theta);

                % avoid overshooting 1
                if tau + Delta_tau > 1
                    Delta_tau = 1 - tau;
                end

                % Brownian increment
                dW = sqrt(Delta_tau) * randn;

                % kernel
                Ktau = (1 - tau)^(H - 0.5) / GammaH;

                % update numerical solution
                Xn = Xn + Ktau * Delta_tau + Ktau * dW;

                % update time and Brownian motion
                tau = tau + Delta_tau;
                Wtau = Wtau + dW;
            end

            Xn_MC(m) = Xn;
        end

        % Monte Carlo approximation 
        EXn_MC = mean(Xn_MC);

        err(k) = abs(EXn_MC - EX);
    end

    plot(n_list, err, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['H = ', num2str(H)]);
end

xlabel('n');
ylabel('|Xn - X| ');
title('Error with random sampling times');
legend('Location','northeast');
grid on;



