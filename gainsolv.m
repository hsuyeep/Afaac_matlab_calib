function g = gainsolv(tolopt, R0, Rhat, g0)
% g = mysolvec4(tolopt, R0, Rhat, g0)
%
% Gain estimation using algorithm developed by Stef Salvini, edited by SJW
% for readability.
%
% Argumets
% tolopt : stop criterion on the change in the norm of the estimated gain
%          vectors in consecutive iterations
% R0     : Nelem x Nelem array covariance matrix mode, zero entries are
%          flagged
% Rhat   : Nelem x Nelem measured array covariance matrix, zero entries are
%          flagged
% g0     : Nelem x 1 vector with initial estimate for gains
%
% Return values
% g      : Nelem x 1 estiamted gain vector
%
% SS & SJW, March 7, 2012

% initialization
Nelem = length(g0);
g = g0;
colnorm = zeros(Nelem, 1);
Rhatnorm = Rhat;
% normalize the measured array covariance matrix (weighting)
for idx = 1:Nelem
    colnorm(idx) = Rhat(:, idx)' * Rhat(:, idx);
    Rhatnorm(:, idx) = Rhat(:, idx) / colnorm(idx);
end
% initial calibration
Rcolcal = complex(zeros(Nelem), zeros(Nelem));
for idx = 1:Nelem
    Rcolcal(:, idx) = g .* R0(:, idx);
end
% start iterating
calest = complex(zeros(Nelem, 1), zeros(Nelem, 1));
gresult = zeros(Nelem, 800);
for iter = 1:800
    % update gain estimate
    for idx = 1:Nelem
        calest(idx) = Rhatnorm(:, idx)' * Rcolcal(:, idx);
    end
    g = 1 ./ conj(calest);
    if (mod(iter, 2) > 0)
        % find new convergence curve
        gest1 = g;
    else
        % update calibration results so far and check for convergence
        gold = g;
        g = (g + gest1) / 2;
        if (iter >= 2)
            dgnorm = norm(g - gold);
            gnorm = norm(g);
            if (dgnorm / gnorm <= tolopt)
                % disp(['convergence reached after ' num2str(iter) ' iterations']);
                % disp(['relative change in norm(g) in the last iteration: ' num2str(dgnorm / gnorm)]);
                break
            end
        end
    end
    % update calibration of array covariance matrix
    for idx = 1:Nelem
        Rcolcal(:, idx) = g .* R0(:, idx);
    end
    gresult(:, iter) = g;
end
%save gainiter gresult
