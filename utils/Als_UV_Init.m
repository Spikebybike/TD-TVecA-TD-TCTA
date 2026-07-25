function [U, V] = Als_UV_Init(S, mask, e_rank, max_iter)
    [m, n] = size(S);
    if nargin < 4
        max_iter = 1;
    end

    V=zeros(n, e_rank)+0.01;

    epsilon = 1e-6;  
    for iter = 1:max_iter

        U = zeros(m, e_rank);
        for i = 1:m
            idx = find(mask(i, :) == 1);
            if length(idx) >= e_rank
                Vsub = V(idx, :);
                Ssub = S(i, idx);
                U(i, :) = (Vsub' * Vsub + epsilon * eye(e_rank)) \ (Vsub' * Ssub(:));
            else
                U(i, :) = 0.01 * randn(1, e_rank);
            end
        end


        V_new = zeros(n, e_rank);
        for j = 1:n
            idx = find(mask(:, j) == 1);
            if length(idx) >= e_rank
                Usub = U(idx, :);
                Ssub = S(idx, j);
                V_new(j, :) = (Usub' * Usub + epsilon * eye(e_rank)) \ (Usub' * Ssub(:));
            else
                V_new(j, :) = 0.01 * randn(1, e_rank);
            end
        end
        V = V_new;
    end
    fprintf('factorization inited\n');
end
