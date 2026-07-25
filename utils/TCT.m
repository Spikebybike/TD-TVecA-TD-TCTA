%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Supporting Function TCT for
%%% "Structured Sparse Millimeter-Wave 3-D SAR Imaging via Truncated-DCT and Toeplitz
%%% Matrix Methods"
%%% Zhuohang Tan, Zeyu Chen, Yiyi Liu, Zhi Li, Yiguang Liu
%%% Developed by Zhuohang Tan, Sichuan University
%%% 2026-7-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------
% Input:  S: M-by-N range slice to be lifted
%         P: pencil parameter of the column-wise Toeplitz transform
%         Q: pencil parameter of the block-level Toeplitz transform
% Output: T: (M*N)-by-(P*Q) two-level Toeplitz-Column-Toeplitz matrix
%-------------------------------------------------------------------
function T= TCT(S, P, Q)
[M,N] = size(S);
T_blocks = cell(N,1);
for i=1:N
    slice=S(:,i);
    T_blocks{i}=toeplitz(slice(1:M),slice([1,M:-1:M-P+2]));
end

T_list = cell(1,Q);

for i = 1:Q
    T_list{i}=cat(1,T_blocks{1:N});
    TEMP=T_blocks{N};
    [T_blocks{2:N}]=deal(T_blocks{1:N-1});
    T_blocks{1}=TEMP;
end

T= cat(2, T_list{:});
end
