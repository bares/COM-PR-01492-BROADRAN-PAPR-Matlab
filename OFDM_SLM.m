clc
clear all
close all

% Reducción de la PAPR:

% Método de Selective Mapping (SLM). Se implementará sobre una
% modulación OFDM estándar, con demodulación OFDM también

tic % Contador del tiempo de ejecución

% Memoria antes de iniciar algoritmo
memAntes = memory;

% Parámetros de la señal OFDM
N = 1024; % Número de subportadoras
simbolos = 1024; % Número de símbolos OFDM
prefijo = 64; % Longitud del prefijo cíclico
M = 16; % Modulación
SNR = 10; % Relación señal a ruido (dB)


% Se generar unos datos aleatorios y se modulan
datos = randi([0 M-1], N, simbolos); % Datos aleatorios
datosmodulados = qammod(datos, M, 'UnitAveragePower', true); % Modulación QAM

% Se realiza la modulación y demodulación OFDM normal
simbolosOFDM = ifft(datosmodulados); %IFFT
simbolos_serial = simbolosOFDM(:); % Se serializan los datos

simbolos_dibujo = abs(simbolosOFDM(:)).^2; % Símbolos para el cálculo de la PAPR

% Cálculo PAPR:
maxima = max(simbolos_dibujo);
media = mean(simbolos_dibujo);
PAPR = 10*log10(maxima/media); % PAPR antes de iniciar el algoritmo
% Método de SLM:

subportadora_actual = datosmodulados;

fase = exp(1j * 2*pi * rand(N,simbolos)); % Fase con la que se va a variar a la subportadora

for i = 1:simbolos

      subportadora = subportadora_actual .* fase(:,i); % Se aplica la variación de fase en la subportadora actual
      fft_inversa = ifft(subportadora, N);
      %Se halla la PAPR de cada variación
      valor_abs = abs(fft_inversa(:)).^2;
      maximo(i) = max(valor_abs);
      media(i) = mean(valor_abs);  
      PAPR2(i) = 10*log10(maximo(i)/media(i)); 

end
PAPRminima = min(PAPR2); % Nos quedamos con la PAPR menor
% Se selecciona el desfase que provoca la menor PAPR en la
% subportadora actual

indice_minimo = min(find(PAPR2 == PAPRminima));
fase_elegida = fase(:, indice_minimo);
subportadora_elegida = subportadora_actual .* fase_elegida;

% Modulación después de realizar la SLM

% Se realiza la modulación y demodulación OFDM normal
simbolosOFDM_slm = ifft(subportadora_elegida); %IFFT

% Agregar prefijo cíclico
prefijoOFDM_slm = [simbolosOFDM_slm(end-prefijo+1:end, :); simbolosOFDM_slm];

% Proceso con SLM:

% Transmitir la señal a través de un canal gaussiano con ruido blanco
senal_recibida_slm = awgn(prefijoOFDM_slm, SNR, 'measured');

% Demodulación: Eliminar el prefijo cíclico
ofdm_recibida_slm = reshape(senal_recibida_slm, N + prefijo, simbolos);
ofdm_recibida_slm = ofdm_recibida_slm(prefijo+1:end, :);

% Se realiza la fft
fftodm1_slm = fft(ofdm_recibida_slm, N);

% Se multiplica la señal por el desfase elegido:

fftodm1_slm = fftodm1_slm.*conj(fase_elegida);

%Se realiza la demodulación
demodulacion_slm = qamdemod(fftodm1_slm, M, 'UnitAveragePower', true);

%Se halla la BER
errores_slm = sum(sum(datos ~= demodulacion_slm));
bits_totales =  N * simbolos * log2(M);
BER_slm = errores_slm / bits_totales;
disp(['El valor de la BER aplicando SLM es: ' num2str(BER_slm)]);
% Medir memoria después del algoritmo
memDespues = memory;

% Cálculo de la memoria consumida total
memConsumida = memDespues.MemUsedMATLAB - memAntes.MemUsedMATLAB;

disp(['Consumo durante el algoritmo SLM: ', num2str(memConsumida / 1e6), ' MB']);
% Proceso sin SLM
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

%Se halla la BER
errores = sum(sum(datos ~= demodulacion));
bits_totales =  N * simbolos * log2(M);
BER = errores / bits_totales;
disp(['El valor de la BER es: ' num2str(BER) ]); 
% Mostrar los resultados:

disp(['La PAPR original es: ' num2str(PAPR) ' dB']);
disp(['La PAPR con SLM es: ' num2str(PAPRminima) ' dB']);
mejora_db = PAPR - PAPRminima;
mejora_porcentaje = 100 - ((PAPRminima/PAPR)*100);

disp(['La mejora en dB es: ' num2str(mejora_db) ]);
disp(['La mejora porcentual es: ' num2str(mejora_porcentaje) ' %']);

% En cuanto a la BER:

disp(['La BER original es: ' num2str(BER)]);
disp(['La BER con SLM es: ' num2str(BER_slm)]);
porcentaje_BER = 100 - ((BER/BER_slm) * 100);
disp(['La diferencia porcentual es: ' num2str(porcentaje_BER) ' %']);
toc