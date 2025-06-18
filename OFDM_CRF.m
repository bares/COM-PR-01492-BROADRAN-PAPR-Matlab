% Reducción de la PAPR:
clc 
clear all
close all


% Método de Creast Factor Reduction (CRF). Se implementará sobre una
% modulación OFDM estándar, con demodulación OFDM también.
tic
% Medir la memoria antes del algoritmo
memAntes = memory;

% Parámetros de la señal OFDM
N = 1024; % Número de subportadoras
simbolos = 1024; % Número de símbolos OFDM
prefijo = 64; % Longitud del prefijo cíclico
M = 16; % Modulación (16-QAM)
SNR = 10; % Relación señal a ruido (dB)

% Se generar unos datos aleatorios y se modulan
datos = randi([0 M-1], N, simbolos); % Datos aleatorios
datosmodulados = qammod(datos, M, 'UnitAveragePower', true); % Modulación QAM

% Se realiza la modulación y demodulación OFDM normal
simbolosOFDM = ifft(datosmodulados); %IFFT
simbolos_serial = simbolosOFDM(:);

simbolos_dibujo = abs(simbolosOFDM(:)).^2;

% Cálculo PAPR:
maxima = max(simbolos_dibujo);
media = mean(simbolos_dibujo);
PAPR = 10*log10(maxima/media);
% Realizamos el clipping:

amplitud  = rms(simbolosOFDM(:))/3; % Amplitud máxima, por encima de aqui se hará el recorte
amplitud_compleja = rms(simbolosOFDM(:))/1.3 + rms(simbolosOFDM(:))/1.3*1i;
simbolos_clip = simbolos_serial;
for i = 1:length(simbolos_dibujo)

    if(real(simbolos_dibujo(i))>amplitud)
        simbolos_clip(i) = amplitud_compleja;
    elseif(imag(simbolos_dibujo(i))>amplitud)  
        simbolos_clip(i) = amplitud_compleja;
    else
        simbolos_clip(i) = simbolos_serial(i);
    end
    
end

% Señal después del clipping:
clip = abs(simbolos_clip(:).^2);

% PAPR después del recorte
maximo_clip = max(clip);
media_clip = mean(clip);
papr_clip = 10*log10(maximo_clip/media_clip);
% Filtrado, según el método del paper
fc = 0.7;  % Este parametro hay que ajustarlo 

alpha = exp(-2.0 * pi * fc); % Factor de atenuación

senal_salida = zeros(size(simbolos_clip));
senal_salida(1) = (1 - alpha) * simbolos_clip(1); % Primer punto a filtar

for i = 2:length(simbolos_clip)
     senal_salida(i) = (1 - alpha) * simbolos_clip(i) + alpha * senal_salida(i - 1);
end

filtering = senal_salida;

simbolos_filt = abs(filtering(:).^2);

% PAPR después de filtering
maxima_filt = max(simbolos_filt(:));
media_filt = mean(simbolos_filt(:));
papr_filt = 10*log10(maxima_filt/media_filt);
% Proceso con filtering:

senal_filt = reshape(filtering, [N simbolos]);

% Agregar prefijo cíclico 
prefijoOFDM_filt = [senal_filt(end-prefijo+1:end, :); senal_filt];

% Serializar la señal OFDM
ofdm_filt = prefijoOFDM_filt(:);

% Transmitir la señal a través de un canal gaussiano con ruido blanco
senal_recibida_filt = awgn(ofdm_filt, SNR, 'measured');

% Demodulación: Eliminar el prefijo cíclico
ofdm_recibida_filt = reshape(senal_recibida_filt, N + prefijo, simbolos);
ofdm_recibida_filt = ofdm_recibida_filt(prefijo+1:end, :);

% Se realiza la fft
fftodm1_filt = fft(ofdm_recibida_filt, N);

%Se realiza la demodulación
demodulacion_filt = qamdemod(fftodm1_filt, M, 'UnitAveragePower', true);


% Medir memoria después del algoritmo
memDespues = memory;

% Cálculo de la memoria consumida total
memConsumida = memDespues.MemUsedMATLAB - memAntes.MemUsedMATLAB;

disp(['Consumo durante el algoritmo CRF: ', num2str(memConsumida / 1e6), ' MB']);
toc
% Proceso sin clipping
% Agregar prefijo cíclico 
prefijoOFDM = [simbolosOFDM(end-prefijo+1:end, :); simbolosOFDM];

% Serializar la señal OFDM
ofdm = prefijoOFDM(:);

% Transmitir la señal a través de un canal gaussiano con ruido blanco
senal_recibida = awgn(ofdm, SNR, 'measured');

% Demodulación: Eliminar el prefijo cíclico
ofdm_recibida = reshape(senal_recibida, N + prefijo, simbolos);
ofdm_recibida = ofdm_recibida(prefijo+1:end, :);

% Se realiza la fft
fftodm1 = fft(ofdm_recibida, N);

%Se realiza la demodulación
demodulacion = qamdemod(fftodm1, M, 'UnitAveragePower', true);

%Se halla la BER del proceso sin correción
errores = sum(sum(datos ~= demodulacion));
bits_totales =  N * simbolos * log2(M);
BER = errores / bits_totales;
%disp(['El valor de la BER es: ' num2str(BER) ]); 
% Proceso con clipping:
senal_clip = reshape(simbolos_clip, [N simbolos]);

% Agregar prefijo cíclico 
prefijoOFDM_clip = [senal_clip(end-prefijo+1:end, :); senal_clip];

% Serializar la señal OFDM
ofdm_clip = prefijoOFDM_clip(:);

% Transmitir la señal a través de un canal gaussiano con ruido blanco
senal_recibida_clip = awgn(ofdm_clip, SNR, 'measured');

% Demodulación: Eliminar el prefijo cíclico
ofdm_recibida_clip = reshape(senal_recibida_clip, N + prefijo, simbolos);
ofdm_recibida_clip = ofdm_recibida_clip(prefijo+1:end, :);

% Se realiza la fft
fftodm1_clip = fft(ofdm_recibida_clip, N);

% Se realiza la demodulación
demodulacion_clip = qamdemod(fftodm1_clip, M, 'UnitAveragePower', true);

% Se halla la BER del proceso con filtrado
errores_filt = sum(sum(datos ~= demodulacion_filt));
BER_filt = errores_filt / bits_totales;
%disp(['El valor de la BER es: ' num2str(BER_filt)]);
% Muestras de los resultados:

disp(['La PAPR original es: ' num2str(PAPR) ' dB']);
disp(['La PAPR después de aplicar CRF es: ' num2str(papr_filt) ' dB']);
mejora_db = PAPR - papr_filt;
mejora_porcentaje = 100-((papr_filt/PAPR)*100);

disp(['La mejora en dB es: ' num2str(mejora_db) ]);
disp(['La mejora porcentual es: ' num2str(mejora_porcentaje) ' %']);
% En cuanto a la BER:

disp(['La BER original es: ' num2str(BER)]);
disp(['La BER con CRF es: ' num2str(BER_filt)]);
porcentaje_BER = 100 - ((BER/BER_filt) * 100);
disp(['La diferencia porcentual es: ' num2str(porcentaje_BER) ' %']);