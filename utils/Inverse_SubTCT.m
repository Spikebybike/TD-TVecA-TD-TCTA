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