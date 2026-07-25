%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Core Function TVNNA for
%%% "Fast Structured Sparse Millimeter-wave 3D SAR Imaging Based on Low-rank and
%%% Smooth Matrix Completion"
%%% Zhuohang Tan, Zeyu Chen, Haojie Tang, Peng Mou, Yiguang Liu
%%% Developed by Zhuohang Tan, Sichuan University
%%% 2026-7-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------
% Input:  S: sparsely sampled range slice; zero entries denote missing samples
%         lambda: weight of the second-order TV smoothness regularization
%         rho: initial ADMM penalty parameter
%         mu: multiplier used to increase rho after each iteration
%         maxIter: maximum number of ADMM iterations
%         tol: relative Frobenius-norm change used as the convergence threshold
% Output: Sc: completed range slice satisfying the observed entries of S
%-------------------------------------------------------------------
function Sc = TVNNA(S, lambda, rho, mu, maxIter, tol)
[M, N] = size(S);
mask = S ~= 0;   % observed-entry mask

%% construct TV difference matrices L (row direction) and R (column direction)

% % first-order TV difference matrix
% L = full(spdiags([-ones(m-1,1), ones(m-1,1)], [0 1], m-1, m));
% R = full(spdiags([-ones(n-1,1), ones(n-1,1)], [0 1], n, n-1));

% % second-order TV difference matrix
L = full(spdiags([ones(M-2,1), -2*ones(M-2,1), ones(M-2,1)],[0 1 2], M-2, M));
R = full(spdiags([ones(N-2,1), -2*ones(N-2,1), ones(N-2,1)],[0 1 2], N, N-2));

LHL_I=L'*L+eye(M);
RRH=R*R';

%% initialize variables
Sc = S;                % initial value of Sc
Sc0=S;
A = L*Sc;
B = Sc*R;
C = Sc;
Y1 = zeros(size(A));
Y2 = zeros(size(B));
Y3 = zeros(size(C));

%% ADMM
for k = 1:maxIter

    % update A
    A = (rho/(rho + 2*lambda))*(L*Sc - Y1/rho);

    % update B
    B = (rho/(rho + 2*lambda))*(Sc*R - Y2/rho);

    % update C
    W = Sc - Y3/rho;
    [U,G,V] = svd(W,'econ');
    S_thresh = diag(max(diag(G) - 1/rho, 0));
    C= U*S_thresh*V';
    C(mask) = S(mask);   % preserve observed entries

    % update Sc (sylvester equation)
    T1 = A + Y1/rho;
    T2 = B + Y2/rho;
    T3 = C + Y3/rho;
    C_rhs = L'*T1 + T2*R' + T3;
    Sc = sylvester(LHL_I, RRH, C_rhs);

    % update lagrange multipliers
    Y1 = Y1 + rho*(A - L*Sc);
    Y2 = Y2 + rho*(B - Sc*R);
    Y3 = Y3 + rho*(C - Sc);
    rho = min(mu * rho, 1e8);

    % check convergence
    err = norm(Sc0-Sc,'fro')/max(norm(Sc0,'fro'),1);
    err = min(err,1);
    if mod(k,5) == 0
        fprintf('Iteration %d/%d: toll = %.6f \n',k,maxIter,err);
    end
    if err < tol
        fprintf('ADMM converged: %d iterations, residual = %.6e\n', k, err);
        break;
    end
    Sc0 = Sc;
end
end
