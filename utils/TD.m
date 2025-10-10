%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function TD for Truncated-DCT
%%% Developed by Zhuohang Tan, Sichuan University
%%% 2025-10-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TDS=TD(S,L)
TDS=dct(S,[],2);
TDS=TDS(:,1:L);
end

