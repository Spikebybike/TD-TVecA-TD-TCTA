function ITDS=ITD(S,N,L)
B=zeros(size(S,1),N-L);
ITDS=cat(2,S,B);
ITDS=idct(ITDS,[],2);
end