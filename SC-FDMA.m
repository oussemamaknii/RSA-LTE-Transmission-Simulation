clc
clear all
close all
 

%% LTE system parameters 
% Carrier frequency 
fc = 800e6; 

% Bandwidth LTE: B = 5 MHZ and NFFT = 512; 
Bandwidth = 5e6; 
NRB = 25;
Ns = 7; 
Nc = 12;
Bsc = 15e3; 
NbofusedSubcarrier = NRB*Nc; 
Bu = NbofusedSubcarrier*Bsc;
TTI = 0.5e-3; 

% Bande totale occupée: On ne sert pas des parties de la bande qui débordent
% des 5 MHz prévus dans la norme LTE 
NFFT  = 512; 
Bsc = 15e3; 
Total_BW = Bsc*NFFT; 

% Temps symbole
Ts = 1/Total_BW; 

Fe = 7.68e6;

% Temps de garde total entre les blocs
Tg = TTI - Ns*NFFT*Ts; 

Tg2 = 36*Ts 
Tg1 = 40*Ts 

CP = [40 36 36 36 36 36 36]; 
RB_size = 512*7+CP(1)+6*CP(2)

%% Charactérisation du canal radio

% Profil de délai du canal radio 
Path_dB =[0,-1,-9,-10,-15,-20];% Puissance des trajets en (dB)
delay =[0,0.3,0.7,1.1,1.7,2.5];% délai du trajet (micro sec)

% Puissance du trajet normalisé en échelle linéaire
Path_var =sqrt(10.^(Path_dB./10))/(sqrt(sum(10.^(0.1*Path_dB))));

% Trouver les valeurs des délais calés avec la période d'échantillonnage
scaled_delay = round(delay*10^(-6)/Ts);

% Ainsi que les puissances calées sur les périodes d'échantillonnage
Ps=zeros(1,max(scaled_delay)+1);

for i=1:size(scaled_delay,2)
    Ps(scaled_delay(i)+1)=Path_var(i);
end

load Multi_Path;
L  =length(Multi_path); 
% Multi_path = zeros(1,L); 
% Multi_path(1) = 1; 
% 

tau = Ts*(0:1:L-1);

Multi_path_mod = Multi_path.*exp(-1i*2*pi*fc*tau);


% Spectre du signal transmis modulé
figure;

f = Fe*linspace(0,1,NFFT);
plot(f,20*log10(abs(fft(Multi_path,NFFT))))
hold on 
plot(f,20*log10(abs(fft(Multi_path_mod,NFFT))), 'r')
title('Spectre fréquentiel du canal multi-trajet')
grid on


%% Emetteur OFDM: Modulation multiporteuse 

load 'text.mat';
bit_stream = char2bit(text); 


map = mapping_def('QPSK');
M = map.nbps; 
NbofUnusedSubcarrier = NFFT - NRB*Nc; 
NbofQAMsymbols = length(bit_stream)/M; 
Nbofblocks = ceil(NbofQAMsymbols/NbofusedSubcarrier); 

symb_stream = mapping(bit_stream,map); 

add_zeros = zeros(1,Nbofblocks*NbofusedSubcarrier - NbofQAMsymbols);
symb_stream = [symb_stream add_zeros]; 

blocks_useful = reshape(symb_stream,length(symb_stream)/Nbofblocks,Nbofblocks).'; 
blocks_zeros = zeros(Nbofblocks,NbofUnusedSubcarrier); 

blocks = [blocks_useful blocks_zeros]; 

blocks_ifft = zeros(Nbofblocks,NFFT); 
blocks_ifft_centre = zeros(Nbofblocks,NFFT); 


index_shift = mod(NbofusedSubcarrier/2:NFFT-1+NbofusedSubcarrier/2,NFFT) + 1; 

for ind_blocks = 1:Nbofblocks
    blocks_ifft(ind_blocks,:) = ifft(blocks(ind_blocks,:));     
    blocks_ifft_centre(ind_blocks,:)  = ifft(blocks(ind_blocks,index_shift)); 
end

% Spectre du signal en bande de base
figure;
f = Fe*linspace(0,1,NFFT)/1e6;
plot(f, 20*log10(abs(blocks_ifft(1,:))))
xlabel('f (MHz)')
ylabel('Amplitude en dB');
title('Spectre fréquentiel du symbol OFDM étalé sur les 5 MHz')
grid on 

Cyclic_Prefix = blocks_ifft(:, NFFT - L + 2:NFFT); 
Cyclic_Prefix_centre = blocks_ifft_centre(:, NFFT - L + 2:NFFT); 

blocks_ifftCP = [Cyclic_Prefix blocks_ifft];
blocks_ifftCP_centre = [Cyclic_Prefix_centre blocks_ifft_centre];

stream_ifftsymbs = []; 
Before_CP1 = zeros(1,CP(1)-L + 1);
Before_CP2 = zeros(1,CP(2)-L + 1);

for i = 1:Nbofblocks
    if (i==1)
        stream_ifftsymbs = [Before_CP1 blocks_ifftCP(i,:)];
        stream_ifftsymbs_centre = [Before_CP1 blocks_ifftCP_centre(i,:)];
    else
        stream_ifftsymbs = [stream_ifftsymbs Before_CP2 blocks_ifftCP(i,:)];  %#ok<AGROW>
        stream_ifftsymbs_centre = [stream_ifftsymbs_centre Before_CP2 blocks_ifftCP_centre(i,:)];  %#ok<AGROW>
    end
end

temps = Ts*(0:length(stream_ifftsymbs)-1);


% Spectre du signal en bande de base
figure;
f = Fe*linspace(0,1,NFFT)/1e6;
plot(f,20*log10(abs(fft(stream_ifftsymbs,NFFT))))
xlabel('f (MHz)')
ylabel('Amplitude en dB');
title('Spectre fréquentiel du signal transmis étalé sur les 5 MHz')
grid on


% Spectre du signal en bande de base
figure;
f = Fe*linspace(0,1,NFFT)/1e6;
plot(f,20*log10(abs(fft(stream_ifftsymbs_centre,NFFT))))
title('Spectre fréquentiel du signal transmis centré autour de 0')
xlabel('f (MHz)')
ylabel('Amplitude en dB');
grid on



%% Mise sous porteuse 


transmitted_signal = stream_ifftsymbs_centre.*exp(2i*pi*fc*temps); 
% Spectre du signal transmis modulé
transmitted_signal_1 = stream_ifftsymbs_centre.*exp(2i*pi*fc*temps);
transmitted_signal_2 = stream_ifftsymbs_centre.*exp(-2i*pi*fc*temps);

figure;
hold on 
k = round(fc/Fe);
f = Fe*linspace(0,1,NFFT)/1e6+ (k*Fe)/1e6;
plot(f,20*log10(abs(fft(transmitted_signal,NFFT))) ,'g')
hold on 
%plot(f,20*log10(abs(fft(transmitted_signal_1,NFFT))) ,'r')
%hold on 
%plot(f,20*log10(abs(fft(transmitted_signal_2,NFFT))) ,'b')
xlabel('f (MHz)')
ylabel('Amplitude en dB');
title('Spectre fréquentiel du signal modulé étalé sur les 5 MHz')
grid on


%% Transmission

load bruit;

received_signal = filter(Multi_path,1,transmitted_signal) + bruit(1:length(transmitted_signal)) ; 

figure;
f = Fe*linspace(0,1,NFFT)/1e6; %+ (k*Fe)/1e6;
plot(f,20*log10(abs(fft(received_signal,NFFT))))
title('Spectre fréquentiel du signal reçu sur les 5 MHz')
xlabel('f (MHz)')
ylabel('Amplitude en dB');
grid on


%% Remise en bande de base et filtrage de la composante autour de 2*fc
demodulated_signal = received_signal.*exp(-2i*pi*fc*temps); 

figure;
f = Fe*linspace(0,1,NFFT)/1E6;
plot(f,20*log10(abs(fft(demodulated_signal,NFFT))))
title('Spectre fréquentiel du signal recu apès remise en bande de base')
xlabel('f (MHz)')
ylabel('Amplitude en dB');

grid on

load 'PasseBas.mat';
demodulated_filter_signal = filter(PasseBas,1,demodulated_signal); 
figure;
f = Fe*linspace(0,1,NFFT)/1E6;
plot(f,20*log10(abs(fft(demodulated_filter_signal,NFFT))))
title('Spectre fréquentiel du signal recu apès remise en bande de base et filtrage')
xlabel('f (MHz)')
ylabel('Amplitude en dB');
grid on


%% Modem OFDM

%demodulated_filter_signal = filter(Multi_path,1,stream_ifftsymbs_centre);
demod_symbs_stream = []; demod_symbs_streamp = []; 
Freq_Channel = fft(Multi_path_mod,NFFT); 
Cf = Freq_Channel; 
Index_start_block(1) = 1; 
Index_end_block(1) = CP(1) +NFFT; 

for i = 2:7
    Index_start_block(i) = Index_end_block(i-1)+1; 
    Index_end_block(i) = Index_start_block(i) + CP(i)+ NFFT - 1; 
end

for i = 1:Nbofblocks
    clear demod_blocks_signal demod_blocks_signal_CP tmp fft_demod_signal_stream demod_symbs ;
    demod_blocks_signal = demodulated_signal(Index_start_block(i):Index_end_block(i));
    demod_blocks_signal_CP = demod_blocks_signal(CP(i)+1:end);
    demod_blocks_signal_CP = demod_blocks_signal_CP;
    fft_demod_signal_stream = fft(demod_blocks_signal_CP,NFFT);
    demod_symbs = fft_demod_signal_stream./Cf;
    demod_symbs_center(index_shift) = demod_symbs;
    demod_symbs_streamp = [demod_symbs_streamp demod_symbs_center];
    demod_symbs_center = demod_symbs_center(1:NbofusedSubcarrier);
    demod_symbs_stream = [demod_symbs_stream demod_symbs_center];
end

figure;
f = Fe*linspace(0,1,NFFT)/1E6;
plot(f,20*log10(abs(fft(demod_symbs_streamp,NFFT))))
hold on 
plot(f,20*log10(abs(fft(stream_ifftsymbs,NFFT))),'r')

title('Spectre fréquentiel du signal recu sur les multiporteuses')
xlabel('f (MHz)')
ylabel('Amplitude en dB');
grid on




figure
plot(demod_symbs_stream(1:NbofQAMsymbols),'g*')
hold on 
plot(map.constellation,'b+')
grid on 


decoded_symbs_stream  = sign(real(demod_symbs_stream)) + 1i*sign(imag(demod_symbs_stream)); 

decoded_bit_stream = demapping(decoded_symbs_stream(1:NbofQAMsymbols),map);

bit2char(decoded_bit_stream)