function [chaine_char] = bit2char(suite_bits)

% Converts a (0/1) sequence into a character chain.
% An 8-bits ASCII code is considered. 
% The original sequence is truncated if 8 does not divide its length.
% parameter: suite_bits is a vector.

dims=size(suite_bits);
if (dims(1) == 1)
    suite_bits=suite_bits';
end;
long = length(suite_bits);
if (mod(long,8) ~= 0)
    long = floor(long/8)*8;
    suite_bits = suite_bits(1:long);
end;
bits_char=num2str(suite_bits)';
chaine_bin=reshape(bits_char,8,long/8)';
chaine_dec=bin2dec(chaine_bin);
chaine_char=(char(chaine_dec))';