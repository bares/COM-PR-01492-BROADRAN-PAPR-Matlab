% Código para implementar el método de Tone Reservation de reducción de la
% PAPR sobre una modulación OFDM. 

% Las portadoras reservadas se usan para insertar celdas que reducen la
% PAPR. 

clc; 
clear; 
close all;

tic
% Memoria antes del inicio del algoritmo
memAntes = memory;

% Parámetros
N_total = 1024;  % Número de subportadoras,ese decir, subportadoras normales y reservadas
tonos = 32;  % Número de subportadoras reservadas, se podría poner de forma relativa al Númeto total también
N = N_total - tonos; % Subportadoras que contienen los datos
nfft = 1024; % Número de puntos para la fft
M = 16;  % Orden de la modulación QAM
iteraciones = 25; % Número de iteraciones del gradiente
prefijo = 64; % Longitud del prefijo cíclico
SNR = 10; % Relación señal a ruido en dB

% Inicio  de Tone reservation

rng('shuffle'); % Random, para distintos resultados para cada iteracción
posiciones_reservadas = randperm(N_total, tonos);

% Generación de los símbolos QAM
datos = randi([0 M-1], N_total, nfft); % Datos aleatorios
datos_modulados = qammod(datos, M, 'UnitAveragePower', true); % Modulación QAM

% Se insertan los tonos reservados, con símbolos QAM
Datos2 = datos_modulados;
for i = 1:nfft
    for j = 1:N_total
        if ismember(j, posiciones_reservadas)
            % Se asigna al tono reservado un símbolo random QAM
            Datos2(j, i) = datos_modulados(randi(N,1), i);         
        end
    end
end

% Reducción del Pico del Kernel
kernel_ini = zeros(N_total, 1);
kernel_ini(posiciones_reservadas) = 1;
Pico_Kernel = sqrt(N_total)/tonos * ifft(kernel_ini);

% Algoritmo para reservar los tonos

Datos_final_TR = ifft(Datos2, N_total);
PAPRdB = zeros(1,iteraciones);
    
    for k = 1:iteraciones
        
        % Se busca el máximo
        [~, position] = max(abs(Datos_final_TR));
        pico = Datos_final_TR(position);

        % Se escala y se halla el tamaño del paso
        escalado = pico / abs(pico);
        differencia = abs(pico) - sqrt(1.2 * mean(abs(Datos_final_TR).^2));
        alpha = escalado * differencia;

        % Redución del pico
        valor_reduccion = alpha * circshift(Pico_Kernel, position-1);

        % Se elimina el pico y se actualiza la señal
        Datos_final_TR(posiciones_reservadas, :) = Datos_final_TR(posiciones_reservadas, :) - valor_reduccion;

        % Se calcula la PAPR
        PAPR = max(abs(Datos_final_TR).^2) / mean(abs(Datos_final_TR).^2);
        PAPRdB(k) = 10*log10(PAPR);
        
    end

PAPR_final_TR = min(PAPRdB);

% Se continua la modulación OFDM con la señal TR:
prefijo_TR = [Datos_final_TR(end-prefijo+1:end, :); Datos_final_TR]; % Se añade el prefijo cíclio
prefijo_TR_serial = prefijo_TR(:); % Se serializa

% Se transmite la señal a través de un canal gaussiano con ruido blanco
senal_recibida_TR = awgn(prefijo_TR_serial, SNR, 'measured');

% Demodulación
ofdm_recibida_TR = reshape(senal_recibida_TR, N_total + prefijo, nfft);
ofdm_recibida_TR = ofdm_recibida_TR(prefijo+1:end, :); %Eliminar el prefijo cíclico

% FFT
fftodm_TR = fft(ofdm_recibida_TR, N_total);

% Se eliminan de la señal las subportadoras anteriormente reservadas:
posiciones_reservadas = sort(posiciones_reservadas);
fftodm_TR(posiciones_reservadas, :) = [];
datos_subp = datos;
datos_subp(posiciones_reservadas, :) =[]; % Lo mismo con los datos, para la comparación para la BER

% Se realiza la demodulación
demodulacion_TR = qamdemod(fftodm_TR, M, 'UnitAveragePower', true);

% Se halla la BER

errores = sum(sum(datos_subp ~= demodulacion_TR));
bits_totales =  size(demodulacion_TR, 1) * size(demodulacion_TR, 2) * log2(M);
BER_TR = errores / bits_totales;
disp(['El valor de la BER usando TR es: ' num2str(BER_TR) ]); 
% Se mide memoria después del algoritmo
memDespues = memory;

% Cálculo de la memoria consumida total
memConsumida = memDespues.MemUsedMATLAB - memAntes.MemUsedMATLAB;

disp(['Consumo durante el algoritmo TR: ', num2str(memConsumida / 1e6), ' MB']);
% Se continua la modulación OFDM con la señal original y se calcula su PAPR
data_ofdm = ifft(datos_modulados, N_total, 1);

prefijo_ofdm = [data_ofdm(end-prefijo+1:end, :); data_ofdm]; % Se añade el prefijo cíclio
prefijo_serial = prefijo_ofdm(:); % Se serializa

% Se transmite la señal a través de un canal gaussiano con ruido blanco
senal_recibida = awgn(prefijo_serial, SNR, 'measured');

% Demodulación
ofdm_recibida = reshape(senal_recibida, N_total + prefijo, nfft);
ofdm_recibida = ofdm_recibida(prefijo+1:end, :); %Eliminar el prefijo cíclico

% Se realiza la fft
fftodm = fft(ofdm_recibida, N_total);

% Se realiza la demodulación
demodulacion = qamdemod(fftodm, M, 'UnitAveragePower', true);

% Se halla la BER
errores = sum(sum(datos ~= demodulacion));
bits_totales =  N_total * nfft * log2(M);
BER = errores / bits_totales;
disp(['El valor de la BER es: ' num2str(BER) ]); 
PAPR_original = max(abs(data_ofdm(:)).^2) ./ mean(abs(data_ofdm(:)).^2);
PAPRdB_original = 10*log10(PAPR_original);
% Mostrar las PAPRs y la mejora, tanto en valor absoluto como en porcentaje
disp(['La PAPR original es: ' num2str(PAPRdB_original) ' dB']);
disp(['La PAPR con TR es: ' num2str(PAPR_final_TR) ' dB']);
mejora_db = PAPRdB_original - PAPR_final_TR;
mejora_porcentaje = 100 - ((PAPR_final_TR/PAPRdB_original) * 100);

disp(['La mejora en dB es: ' num2str(mejora_db)]);
disp(['La mejora porcentual es: ' num2str(mejora_porcentaje) ' %']);
% En cuanto a la BER: 
disp(['La BER original es: ' num2str(BER)]);
disp(['La BER con TR es: ' num2str(BER_TR)]);
porcentaje_BER = 100 - ((BER/BER_TR) * 100);
disp(['La diferencia porcentual es: ' num2str(porcentaje_BER) ' %']);
toc