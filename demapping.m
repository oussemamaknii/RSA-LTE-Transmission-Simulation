function dec_b = demapping (vect , map)

for ind = 1:length(vect)
    tmp = find( abs(vect(ind) - map.constellation) < 1e-1)-1;    
    tmp_bits(:,ind) = de2bi(tmp,map.nbps,'left-msb');
end
dec_b = reshape(tmp_bits,1, length(vect)*map.nbps);