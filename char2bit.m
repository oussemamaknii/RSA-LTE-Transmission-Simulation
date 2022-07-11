function [suite_bits] = char2bit(chaine_char)

% Converts a character chain into a bit sequence (0/1) 
% thanks to an 8 bit-ASCII code.
% parameter : in quotes ('example')

chaine_dec=double(chaine_char);
chaine_bin=dec2bin(chaine_dec,8);
bits_ligne=reshape(chaine_bin',1,length(chaine_dec)*8);
suite_bits=str2num(bits_ligne')';
