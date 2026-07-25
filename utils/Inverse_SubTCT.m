%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Supporting Function Inverse_SubTCT for
%%% "Structured Sparse Millimeter-Wave 3-D SAR Imaging via Truncated-DCT and Toeplitz
%%% Matrix Methods"
%%% Zhuohang Tan, Zeyu Chen, Yiyi Liu, Zhi Li, Yiguang Liu
%%% Developed by Zhuohang Tan, Sichuan University
%%% 2026-7-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------
% Input:  T: M-by-R Toeplitz submatrix produced by the first TCT level
% Output: V: M-by-1 source vector obtained by diagonal averaging and index unwrapping
%-------------------------------------------------------------------
function V = Inverse_SubTCT(T)
[M, R] = size(T);
L = M + R - 1;
[row, col] = ndgrid(1:M, 1:R);

diag_idx = col - row;
pos = diag_idx + M;

sums=accumarray(pos(:),T(:),[L,1],@sum,0);
counts=accumarray(pos(:),1,[L,1],@sum,0);
V = sums ./ counts;
V(1:M)=V(M:-1:1);
V(M+1:M+R-1)=V(M+R-1:-1:M+1);
V(M-R+2:M)=(V(M-R+2:M)+V(M+1:M+R-1))/2;
V=V(1:M);
end
