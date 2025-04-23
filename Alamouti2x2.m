clear all; close all; clc;

frmLen = 1000; % Chiều dài frame
numPackets = 1000; % Số lượng packet (tăng để BER chính xác hơn)
EbN0 = 0:2:20;
N = 2; % Số anten phát Tx
M = 2; % Số anten thu Rx (không dùng trong mã này)
P = 4; % Điều chế QPSK

modTypes = {'QPSK', '16QAM', '64QAM'}; % Các loại điều chế: 'QPSK', '16QAM', '64QAM'
BER_MRC = zeros(length(modTypes), length(EbN0)); % Khởi tạo mảng lưu trữ BER cho mỗi loại điều chế

for modIdx = 1:length(modTypes)
    modType = modTypes{modIdx}; % Loại điều chế hiện tại
    for idx = 1:length(EbN0)
        error22 = zeros(1, numPackets);
        for packetIdx = 1:numPackets

            % Tạo dữ liệu điều chế tùy chọn
            switch modType
                case 'QPSK'
                    data = randi([0 P-1], frmLen, 2); % Dữ liệu QPSK (frmLen x 2)
                    X = pskmod(data, P, 0); % Điều chế QPSK
                case '16QAM'
                    P = 16; % Đổi độ dài cho 16QAM
                    data = randi([0 P-1], frmLen, 2); % Dữ liệu 16QAM
                    X = qammod(data, P); % Điều chế 16QAM
                case '64QAM'
                    P = 64; % Đổi độ dài cho 64QAM
                    data = randi([0 P-1], frmLen, 2); % Dữ liệu 64QAM
                    X = qammod(data, P); % Điều chế 64QAM
                otherwise
                    error('Modulation type not recognized!');
            end
             
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
            N0 = Eb / (10^(EbN0(idx)/10)); % Mật độ nhiễu
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
            switch modType
                case 'QPSK'
                    decodedData_MRC = pskdemod(z, P, 0); % Giải điều chế QPSK
                case '16QAM'
                    decodedData_MRC = qamdemod(z, P); % Giải điều chế 16QAM
                case '64QAM'
                    decodedData_MRC = qamdemod(z, P); % Giải điều chế 64QAM
            end
            errorMRC(packetIdx) = biterr(decodedData_MRC, data);
        end
        % Tính BER cho mỗi loại điều chế
        BER_MRC(modIdx, idx) = sum(errorMRC) / (numPackets * frmLen * log2(P));
        fprintf('ModType: %s, Eb/N0 = %d dB, BER = %e\n', modType, EbN0(idx), BER_MRC(modIdx, idx));
    end
end

% Vẽ đồ thị
h = figure;
semilogy(EbN0, BER_MRC(1,:), 'r*-', 'LineWidth', 1.5, 'DisplayName', 'QPSK');
hold on;
semilogy(EbN0, BER_MRC(2,:), 'bs-', 'LineWidth', 1.5, 'DisplayName', '16QAM');
semilogy(EbN0, BER_MRC(3,:), 'g^-', 'LineWidth', 1.5, 'DisplayName', '64QAM');
grid on;
set(gca, 'yscale', 'log', 'xlim', [EbN0(1), EbN0(end)], 'ylim', [1e-5 1]);
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('Hệ thống Alamouti 2x2 với các dạng điều chế');
legend('show');
set(h, 'color', [1 1 1]);
set(h, 'NumberTitle', 'off');
set(h, 'renderer', 'zbuffer');