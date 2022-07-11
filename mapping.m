function  symb = mapping(bits,map)

card = map.nbps;
l_data = length(bits);
tmp = reshape (bits,card,l_data/card); 
symb = map.constellation(bi2de(tmp','left-msb' )+1);


% Local Variables: ***
% mode: Matlab ***
% End: ***
