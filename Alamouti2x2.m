clear all; close all; clc;

frmLen = 100; % Chiều dài frame
numPackets = 1000; % Số lượng packet (tăng để BER chính xác hơn)
EbNo = 0:2:20;
N = 2; % Số anten phát Tx
M = 2; % Số anten thu Rx (không dùng trong mã này)
P = 4; % Điều chế QPSK

for idx = 1:length(EbNo)
    error22 = zeros(1, numPackets);
    for packetIdx = 1:numPackets
        % Tạo dữ liệu
        data = randi([0 P-1], frmLen, 2); % 2 ký hiệu cho Alamouti 2x2
        X = pskmod(data, P, 0); % Điều chế QPSK
        
        % Symbol
        s1 = X(:,1); % 1x1
        s2 = X(:,2); % 1x1

        % Ma trận truyền
        tx1 = [s1 s2]; % Thời điểm 1
        tx2 = [-conj(s2) conj(s1)]; % Thời điểm 2

        % Kênh truyền
        Hr = (randn(frmLen, N, 2) + 1i*randn(frmLen, N, 2)) / sqrt(2); % 2 kênh
        h1 = reshape(Hr(:,:,1), frmLen, N); % 1x2
        h2 = reshape(Hr(:,:,2), frmLen, N); % 1x2

        % Nhiễu AWGN
        Es = 1; % Công suất ký hiệu
        Eb = Es / log2(P); % Công suất bit
        N0 = Eb / (10^(EbNo(idx)/10)); % Mật độ nhiễu
        noise = sqrt(N0/2) * (randn(frmLen, N) + 1i*randn(frmLen, N));

        % Tín hiệu thu
        r1 = h1.*s1 + h2.*s2 + noise; % Thời điểm 1
        r2 = h1.*(-conj(s2)) + h2.*conj(s1) + noise; % Thời điểm 2
        
        % Giải mã Alamouti
        z1 = conj(h1).*r1 + h2.*conj(r2);
        z2 = conj(h2).*r1 - h1.*conj(r2);
        
        z(:,1) = sum(z1,2);
        z(:,2) = sum(z2,2);
        % Giải điều chế
        demod22 = pskdemod(z, P, 0);
        error22(packetIdx) = biterr(demod22, data);
    end
    BER22(idx) = sum(error22) / (numPackets * frmLen * log2(P));
    fprintf('Eb/N0 = %d dB, BER = %e\n', EbNo(idx), BER22(idx));
end

% Vẽ biểu đồ
h = figure;
semilogy(EbNo, BER22, 'r*-', 'LineWidth', 1.5, 'DisplayName', 'Alamouti 2x2');
grid on;
set(gca, 'yscale', 'log', 'xlim', [EbNo(1), EbNo(end)], 'ylim', [1e-5 1]);
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('Hệ thống MIMO 2x2 Alamouti');
legend('show');
set(h, 'color', [1 1 1]);
set(h, 'NumberTitle', 'off');
set(h, 'renderer', 'zbuffer');