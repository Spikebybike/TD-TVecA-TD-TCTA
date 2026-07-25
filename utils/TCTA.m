%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Core Function TCTA for
%%% "Structured Sparse Millimeter-Wave 3-D SAR Imaging via Truncated-DCT and Toeplitz Matrix Methods",
%%% Zhuohang Tan, Zeyu Chen, Yiyi Liu, Zhi Li, Yiguang Liu
%%% Developed by Zhuohang Tan, Sichuan University
%%% 2026-7-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------
% Input:  S: sparsely sampled range slice; zero entries denote missing samples
%         P,Q: first- and second-level pencil parameters of the TCT lifting
%         mu: ADMM penalty parameter
%         max_iter_K: maximum number of ADMM iterations
%         e_rank: estimated rank of the lifted TCT matrix
% Output: X_completed: completed range slice in the original spatial domain
%-------------------------------------------------------------------
function X_completed = TCTA(S, P, Q, mu, max_iter_K,e_rank)

%% initialization
mask=S~=0;
mask=double(mask);     % observed-entry mask
[M, N] = size(S);
s_TCT = TCT(S, P, Q);
mask_TCT= TCT(mask, P, Q);
[U, V] = Als_UV_Init(s_TCT, mask_TCT, e_rank);   % factorization of  TCT matrix
Z = zeros(size(s_TCT));
X = S;

%% ADMM
for k = 1:max_iter_K
    X0=X;
    % update X
    T_UV = U * V' - Z/mu;
    X = mask .* S + (1 - mask) .* Inverse_TCT(T_UV, M, N, P, Q);

    % update U/V
    T_X = TCT(X, P, Q);
    temp = mu*T_X + Z;
    U = temp * V / (eye(size(V,2)) + mu * (V' * V));
    V = temp' * U / (eye(size(U,2)) + mu * (U' * U));

    % update Z
    Z = Z+mu*(T_X - U * V');

    toll=norm(X - X0, 'fro')/norm(X0, 'fro');
    if mod(k,5)==0
        fprintf('Iteration %d/%d: toll = %.6f \n',k,max_iter_K,toll);
    end
end
X_completed = X;
end
