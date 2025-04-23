% Mô phỏng Alamouti 2x1 với kênh TDL-D 5G, tích hợp OFDM
clear all; close all; clc;

frmLen = 1000; % chiều dài frame (số ký hiệu)
numPackets = 100; % số lượng packet
EbN0 = 0:2:20; % giá trị Eb/N0 từ 0 đến 20 dB, bước 2 dB
N = 2; % số lượng anten Tx
M = 1; % số lượng anten thu Rx

% Cấu hình OFDM (5G NR)
bandwidth = 100e6; % Băng thông 100 MHz
scs = 30e3; % Subcarrier spacing 30 kHz
numSubcarriers = 3300; % Số subcarrier (theo chuẩn 5G NR cho băng thông 100 MHz)
fftSize = 4096; % Kích thước FFT
cpLength = 288; % Cyclic prefix length (theo chuẩn 5G NR)
sampleRate = 122.88e6; % Tần số lấy mẫu (122.88 MHz)

% Cấu hình kênh TDL-D (5G)
channel = nrTDLChannel;
channel.DelayProfile = 'TDL-D'; % Kênh TDL-D (Rayleigh, không LOS)
channel.NumTransmitAntennas = N;
channel.NumReceiveAntennas = M;
channel.DelaySpread = 300e-9; % Trễ đa đường (300 ns)
channel.MaximumDopplerShift = 324; % Doppler shift (100 km/h, 3.5 GHz)
channel.SampleRate = sampleRate; % Tần số lấy mẫu

% Danh sách các loại điều chế
modTypes = {'QPSK', '16QAM', '64QAM'};
P_values = [4, 16, 64]; % Số điểm trong chòm sao
bitsPerSymbol = [2, 4, 6]; % Số bit trên mỗi ký hiệu

% Khởi tạo mảng lưu trữ BER và dung năng
BER = zeros(length(modTypes), length(EbN0)); % BER cho từng loại điều chế
capacity = zeros(length(modTypes), length(EbN0)); % Dung năng trung bình
error = zeros(length(modTypes), numPackets); % Lỗi cho từng packet

% Vòng lặp cho các loại điều chế
for modIdx = 1:length(modTypes)
    P = P_values(modIdx); % Độ dài điều chế
    bitsPerSym = bitsPerSymbol(modIdx); % Số bit trên mỗi ký hiệu
    
    for idx = 1:length(EbN0)
        capacityPerPacket = zeros(1, numPackets); % Dung năng cho mỗi packet
        for packetIdx = 1:numPackets
            % Tạo dữ liệu
            data = randi([0 P-1], frmLen, 2); % Dữ liệu (frmLen x 2)
            X = pskmod(data, P, 0); % Điều chế (QPSK, 16QAM, 64QAM)

            % Truyền tín hiệu Alamouti trên từng subcarrier
            x1 = X(:,1); % Ký hiệu x1
            x2 = X(:,2); % Ký hiệu x2

            % OFDM: Chia tín hiệu thành các subcarrier
            txSig = zeros(fftSize, 2); % Tín hiệu trên 2 anten
            txSig(1:numSubcarriers/2, 1) = x1(1:numSubcarriers/2); % Anten 1, khe 1
            txSig(1:numSubcarriers/2, 2) = x2(1:numSubcarriers/2); % Anten 2, khe 1
            txSig(end-numSubcarriers/2+1:end, 1) = -conj(x2(1:numSubcarriers/2)); % Anten 1, khe 2
            txSig(end-numSubcarriers/2+1:end, 2) = conj(x1(1:numSubcarriers/2)); % Anten 2, khe 2

            % IFFT và thêm Cyclic Prefix
            txTime = zeros(fftSize + cpLength, 2);
            for ant = 1:2
                txTime(:,ant) = [ifft(txSig(:,ant), fftSize); zeros(cpLength, 1)];
                txTime(1:cpLength, ant) = txTime(end-cpLength+1:end, ant); % Thêm CP
            end

            % Truyền qua kênh TDL-D
            rxTime = channel(txTime); % Truyền qua kênh
            rxTime = awgn(rxTime, EbN0(idx) + 10*log10(bitsPerSym), 'measured'); % Thêm nhiễu AWGN

            % Loại bỏ Cyclic Prefix và FFT
            rxSig = fft(rxTime(cpLength+1:end,:), fftSize); % FFT

            % Ước lượng kênh (giả định dùng pilot)
            Hr = (randn(1, N, M) + 1i*randn(1, N, M)) / sqrt(2); % Kênh giả định
            Hr = repmat(Hr, numSubcarriers, 1, 1);
            h1 = Hr(:,1,1); % Kênh từ anten 1
            h2 = Hr(:,2,1); % Kênh từ anten 2

            % Giải mã MRC trên từng subcarrier
            r1 = rxSig(1:numSubcarriers/2, 1); % Khe 1
            r2 = rxSig(end-numSubcarriers/2+1:end, 1); % Khe 2
            norm = abs(h1).^2 + abs(h2).^2;
            x1Hat = (conj(h1).*r1 + h2.*conj(r2))./norm;
            x2Hat = (conj(h2).*r1 - h1.*conj(r2))./norm;

            % Gộp tín hiệu
            decodedSymbols = [x1Hat; x2Hat];
            decodedData = pskdemod(decodedSymbols, P, 0);

            % Tính lỗi
            error(modIdx, packetIdx) = biterr(decodedData, data(1:numSubcarriers,:));

            % Tính dung năng
            rho_dB = EbN0(idx) + 10*log10(bitsPerSym);
            rho = 10^(rho_dB/10);
            normH = abs(h1).^2 + abs(h2).^2;
            C = log2(1 + (rho/N) * normH);
            capacityPerPacket(packetIdx) = mean(C);
        end

        % Tính BER
        BER(modIdx, idx) = sum(error(modIdx,:)) / (numPackets * numSubcarriers * bitsPerSym);
        % Tính dung năng trung bình
        capacity(modIdx, idx) = mean(capacityPerPacket);
        fprintf('Modulation: %s, Eb/N0 = %d dB, BER = %e, Capacity = %f bits/s/Hz\n', ...
            modTypes{modIdx}, EbN0(idx), BER(modIdx, idx), capacity(modIdx, idx));
        end
end

% Vẽ đồ thị BER
h1 = figure;
colors = {'r*-', 'b*-', 'g*-'};
for modIdx = 1:length(modTypes)
    semilogy(EbN0, BER(modIdx,:), colors{modIdx}, 'LineWidth', 1.5, 'DisplayName', modTypes{modIdx});
    hold on;
end
grid on;
set(gca, 'yscale', 'log', 'xlim', [EbN0(1), EbN0(end)], 'ylim', [1e-5 1]);
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('Hệ thống Alamouti 2x1 MRC với các loại điều chế trên kênh TDL-D 5G');
legend('show');
set(h1, 'color', [1 1 1]);
set(h1, 'NumberTitle', 'off');
set(h1, 'renderer', 'zbuffer');

% Vẽ đồ thị dung năng
h2 = figure;
for modIdx = 1:length(modTypes)
    plot(EbN0, capacity(modIdx,:), colors{modIdx}(1), 'LineWidth', 1.5, 'DisplayName', modTypes{modIdx});
    hold on;
end
grid on;
xlabel('Eb/N0 (dB)'); ylabel('Capacity (bits/s/Hz)');
title('Dung năng trung bình của hệ thống Alamouti 2x1 trên kênh TDL-D 5G');
legend('show');
set(h2, 'color', [1 1 1]);
set(h2, 'NumberTitle', 'off');
set(h2, 'renderer', 'zbuffer');