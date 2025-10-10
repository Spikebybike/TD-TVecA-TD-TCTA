%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function ITD for zero-padding and IDCT
%%% Developed by Zhuohang Tan, Sichuan University
%%% 2025-10-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ITDS=ITD(S,N,L)
B=zeros(size(S,1),N-L);
ITDS=cat(2,S,B);
ITDS=idct(ITDS,[],2);
end
