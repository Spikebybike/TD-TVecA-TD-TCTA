function Sc = TVNNA(S, lambda, rho, mu, maxIter, tol)
% TVNNA: TV + 核范数 +ADMM 矩阵补全
% 输入:
%   S       - 缺失矩阵，缺失值为0
%   lambda  - TV惩罚参数
%   rho     - ADMM惩罚参数
%   mu     - rho增速
%   maxIter - 最大迭代次数
%   tol     - 收敛阈值（例如1e-4）
% 输出:
%   Sc      - 补全后的矩阵

[m, n] = size(S);
Omega = S ~= 0;   % 观测位置掩码

%% 构造差分矩阵 L (行方向) 和 R (列方向)
% % first-order TV
% L = full(spdiags([-ones(m-1,1), ones(m-1,1)], [0 1], m-1, m));
% R = full(spdiags([-ones(n-1,1), ones(n-1,1)], [0 1], n, n-1));

% % second-order TV
L = full(spdiags([ones(m-2,1), -2*ones(m-2,1), ones(m-2,1)],[0 1 2], m-2, m));
R = full(spdiags([ones(n-2,1), -2*ones(n-2,1), ones(n-2,1)],[0 1 2], n, n-2));

LHL_I=L'*L+eye(m);
RRH=R*R';

%% 初始化变量
Sc = S;                % X 的初值
Sc0=S;
A = L*Sc;
B = Sc*R;
C = Sc;

Y1 = zeros(size(A));
Y2 = zeros(size(B));
Y3 = zeros(size(C));

%% ADMM迭代
for k = 1:maxIter

    % 1. 更新 A
    A = (rho/(rho + 2*lambda))*(L*Sc - Y1/rho);

    % 2. 更新 B
    B = (rho/(rho + 2*lambda))*(Sc*R - Y2/rho);

    % 3. 更新 C (核范数 + 观测约束)
    W = Sc - Y3/rho;
    [U,G,V] = svd(W,'econ');
    S_thresh = diag(max(diag(G) - 1/rho, 0));
    C= U*S_thresh*V';
    C(Omega) = S(Omega);   % 强制观测位置为原值

    % 4. 更新 Sc (Sylvester方程)
    T1 = A + Y1/rho;
    T2 = B + Y2/rho;
    T3 = C + Y3/rho;
    C_rhs = L'*T1 + T2*R' + T3;
    Sc = sylvester(LHL_I, RRH, C_rhs);

    % 5. 更新拉格朗日乘子
    Y1 = Y1 + rho*(A - L*Sc);
    Y2 = Y2 + rho*(B - Sc*R);
    Y3 = Y3 + rho*(C - Sc);

    rho = min(mu * rho, 1e8);

    % 6. 收敛判断
    err = norm(Sc0-Sc,'fro')/max(norm(Sc0,'fro'),1);
    err = min(err,1);

    if mod(k,5) == 0
        fprintf('迭代 %d 次, 收敛残差 = %.6e\n', k, err);
    end

    if err < tol
        fprintf('ADMM收敛: 迭代 %d 次, 收敛残差 = %.6e\n', k, err);
        break;
    end

    Sc0 = Sc;
end

end
