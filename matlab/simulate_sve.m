
clear; clc;

H_list = [0.1, 0.3, 0.4];
n_list = [5, 8,100];
M = 80000;     % Monte Carlo size

figure; hold on;

for h = 1:length(H_list)
    H = H_list(h);
    GammaH = gamma(H + 0.5);

    % Exact solution
    EX = 1 / (GammaH * (H + 0.5));

    err = zeros(length(n_list),1);

    for k = 1:length(n_list)
        n = n_list(k);
        Delta = 1 / n;

        Xn_MC = zeros(M,1);

        for m = 1:M
            Xn = 0;
            for j = 1:n
                t = (j-1) / n;
                Kj = (1 - t)^(H - 0.5) / GammaH;

                dW = sqrt(Delta) * randn;

                Xn = Xn + Kj * Delta + Kj * dW;
            end
            Xn_MC(m) = Xn;
        end

        % Monte Carlo approximation 
        EXn_MC = mean(Xn_MC);

        err(k) = abs(EXn_MC - EX);
    end

    loglog(n_list, err, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', ['H = ', num2str(H)]);
end

xlabel('n');
ylabel('| X_n-X|');
title('Error of the (Xn-X)');
legend('Location','southwest');
grid on;



