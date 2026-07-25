%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Supporting Function Inverse_TCT for
%%% "Structured Sparse Millimeter-Wave 3-D SAR Imaging via Truncated-DCT and Toeplitz
%%% Matrix Methods"
%%% Zhuohang Tan, Zeyu Chen, Yiyi Liu, Zhi Li, Yiguang Liu
%%% Developed by Zhuohang Tan, Sichuan University
%%% 2026-7-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------
% Input:  T: (M*N)-by-(P*Q) two-level Toeplitz lifted matrix
%         M,N: row and column counts of the original range slice
%         P,Q: first- and second-level TCT pencil parameters
% Output: S: M-by-N range slice recovered by averaging repeated Toeplitz entries
%-------------------------------------------------------------------
function S= Inverse_TCT(T, M, N, P,Q)
S = zeros(M, N);
rb= M;cb = P;
diag_blocks = zeros((N+Q-1)*rb,cb);
for diag_idx=-(N-1):Q-1
    current_block=zeros(rb,cb);
    num=0;
    if diag_idx <= 0
        start_row = 1-diag_idx;
        end_row = N;
    else
        start_row = 1;
        end_row = Q - diag_idx;
    end
    for row = start_row:min(end_row, N)
        col = row + diag_idx;
        if col >= 1 && col <= Q
            current_block = current_block + T((row-1)*rb+1:row*rb,(col-1)*cb+1:col*cb);
            num=num+1;
        end
    end
    current_block=current_block/num;
    diag_blocks((diag_idx+N-1)*rb+1:(diag_idx+N)*rb,:) = current_block;
end
blocks_cell = mat2cell(diag_blocks, repmat(rb, N+Q-1, 1), cb);
blocks_cell(1:N) = blocks_cell(N:-1:1);
blocks_cell(N+1:N+Q-1) = blocks_cell(N+Q-1:-1:N+1);

for i=1:Q-1
    T1=blocks_cell{N-Q+1+i};
    T2=blocks_cell{N+i};
    T=(T1+T2)/2;
    blocks_cell{N-Q+1+i}=T;
end

for i=1:N
    S(:,i)=Inverse_SubTCT(blocks_cell{i});
end

end
