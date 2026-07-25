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