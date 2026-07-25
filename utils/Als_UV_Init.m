%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Supporting Function Als_UV_Init for
%%% "Structured Sparse Millimeter-Wave 3-D SAR Imaging via Truncated-DCT and Toeplitz
%%% Matrix Methods"
%%% Zhuohang Tan, Zeyu Chen, Yiyi Liu, Zhi Li, Yiguang Liu
%%% Developed by Zhuohang Tan, Sichuan University
%%% 2026-7-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------
% Input:  S: observed matrix to factorize, typically a TCT lifted matrix
%         mask: observation mask of S (1 for observed entries and 0 for missing entries)
%         e_rank: estimated factorization rank
%         max_iter: optional number of ALS iterations (default: 1)
% Output: U,V: submatrices satisfying S approximately equal to U*V'
%-------------------------------------------------------------------
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
