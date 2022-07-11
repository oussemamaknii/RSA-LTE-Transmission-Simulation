function [map] = mapping_def(map_type) 
%===================================================================
%=                                                                 =
%= [map] = mapping_def(map_type)                                   =
%=                                                                 =
%= This function checks if the selected mapping (map_type)         =
%= is implemented.                                                 =
%= If so, the number of bits per symbol (nbps), the normalization  =
%= factor of the constellation (norm), the constellation           =
%= (constellation), and the bit weights to build symbols (weights) = 
%= are defined.                                                    =
%= All parameters are saved in the structure map.                  =
%=                                                                 =
%===================================================================

switch map_type
    %=========================
    case 'BPSK'
        %==========================
        map.nbps = 1;
        map.norm = 1;
        map.card = log2(2);
        map.constellation = [-1 1];
        map.cmin = 0;
        map.cmax = 1;

    case 'QPSK'
        %==========================
        map.nbps = 2;
        map.norm = 2;
        map.card = log2(4);
        map.constellation = [-1-1i -1+1i 1-1i 1+1i];
        map.cmin = 0;
        map.cmax = 1;

    case '8QAM'
        %==========================
        map.nbps = 3;
        map.norm = 6;
        map.card = log2(8);
        map.constellation = [-3+1i -1+1i 3+1i 1+1i -3-1i -1-1i 3-1i 1-1i];
        map.cmin = 0;
        map.cmax = 2;

    case '16QAM'
        %==========================
        map.nbps = 4;
        map.norm = 10;
        map.card = log2(16);
        map.constellation = ...
            [-3-3i -3-1i -3+3i -3+1i -1-3i -1-1i -1+3i -1+1i ...
            3-3i 3-1i 3+3i 3+1i 1-3i 1-1i 1+3i 1+1i];
        map.cmin = 0;
        map.cmax = 3;

    case '64QAM'
        %==========================
        map.nbps = 6;
        map.norm = 42;
        map.card = log2(64);
        map.constellation = ...
            [-7-7i -7-5i -7-1i -7-3i -7+7i -7+5i -7+1i -7+3i ...
            -5-7i -5-5i -5-1i -5-3i -5+7i -5+5i -5+1i -5+3i ...
            -1-7i -1-5i -1-1i -1-3i -1+7i -1+5i -1+1i -1+3i ...
            -3-7i -3-5i -3-1i -3-3i -3+7i -3+5i -3+1i -3+3i...
            +7-7i +7-5i +7-1i +7-3i +7+7i +7+5i +7+1i +7+3i...
            +5-7i +5-5i +5-1i +5-3i +5+7i +5+5i +5+1i +5+3i...
            +1-7i +1-5i +1-1i +1-3i +1+7i +1+5i +1+1i +1+3i...
            +3-7i +3-5i +3-1i +3-3i +3+7i +3+5i +3+1i +3+3i];
        map.cmin = 0;
        map.cmax = 7;


    case '256QAM'
        %==========================
        map.nbps = 8;
        map.norm = 170;
        map.constellation = ...
            (repmat([-15i;-13i;-9i;-11i;-1i;-3i;-7i;-5i;15i;13i;9i;11i;i;3i;7i;5i],16,1)...
            +kron([-15;-13;-9;-11;-1;-3;-7;-5;15;13;9;11;1;3;7;5],ones(16,1)));
        map.constellation = map.constellation.';
        map.cmin = 0;
        map.cmax = 15;

    otherwise
        %==========================
        error('\n ERROR: Non defined constellation \n')

end
%==========================

map.type = map_type;
map.weights = 2.^(map.nbps-1:-1:0);
map.mapping = 'Z8';
%map.constellation = 1/sqrt(map.norm)*map.constellation;
map.constellation = map.constellation;



% Local Variables: ***
% mode: Matlab ***
% End: ***