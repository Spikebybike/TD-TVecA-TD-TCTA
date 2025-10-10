function TDS=TD(S,L)
TDS=dct(S,[],2);
TDS=TDS(:,1:L);
end