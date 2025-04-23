% Mô phỏng STBC Alamouti 2x1 với MRC và ML
clear all; close all; clc;

frmLen = 1000; % chiều dài frame (số ký hiệu)
numPackets = 1000; % số lượng packet
EbN0 = 0:2:20; % giá trị Eb/N0 từ 0 đến 20 dB, bước 2 dB
N = 2; % số lượng anten Tx
M = 1; % số lượng anten thu Rx
P = 4; % độ dài điều chế QPSK (mặc định)

modTypes = {'QPSK', '16QAM', '64QAM'}; % Các loại điều chế: 'QPSK', '16QAM', '64QAM'
BER_MRC = zeros(length(modTypes), length(EbN0)); % Khởi tạo mảng lưu trữ BER cho mỗi loại điều chế

% Vòng lặp cho từng loại điều chế
for modIdx = 1:length(modTypes)
    modType = modTypes{modIdx}; % Loại điều chế hiện tại

    % Vòng lặp cho Eb/N0 điểm
    for idx = 1:length(EbN0)
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

            % Mã hóa Alamouti
            x1 = X(:,1); % Ký hiệu x1 (frmLen x 1)
            x2 = X(:,2); % Ký hiệu x2 (frmLen x 1)

            % Kênh Rayleigh fading
            Hr = (randn(frmLen, N, M) + 1i*randn(frmLen, N, M)) / sqrt(2);
            h1 = Hr(:,1,1); % Kênh từ anten 1
            h2 = Hr(:,2,1); % Kênh từ anten 2

            % Truyền tín hiệu
            s1 = h1.*x1 + h2.*x2; % Tín hiệu tại khe thời gian 1
            s2 = -h1.*conj(x2) + h2.*conj(x1); % Tín hiệu tại khe thời gian 2
            r1 = awgn(s1, EbN0(idx) + 3, 'measured'); % Thêm nhiễu (QPSK: Es/N0 = Eb/N0 + 3 dB)
            r2 = awgn(s2, EbN0(idx) + 3, 'measured');

            % Giải mã MRC
            norm = abs(h1).^2 + abs(h2).^2;
            x1Hat_MRC = (conj(h1).*r1 + h2.*conj(r2))./norm;
            x2Hat_MRC = (conj(h2).*r1 - h1.*conj(r2))./norm;

            % Gộp tín hiệu MRC
            decodedSymbols_MRC = [x1Hat_MRC x2Hat_MRC];
            
            % Giải điều chế tùy thuộc vào loại điều chế
            switch modType
                case 'QPSK'
                    decodedData_MRC = pskdemod(decodedSymbols_MRC, P, 0); % Giải điều chế QPSK
                case '16QAM'
                    decodedData_MRC = qamdemod(decodedSymbols_MRC, P); % Giải điều chế 16QAM
                case '64QAM'
                    decodedData_MRC = qamdemod(decodedSymbols_MRC, P); % Giải điều chế 64QAM
            end

            % Tính lỗi cho MRC
            errorMRC(packetIdx) = biterr(decodedData_MRC, data);
        end

        % Tính BER cho mỗi loại điều chế
        BER_MRC(modIdx, idx) = sum(errorMRC) / (numPackets * frmLen * log2(P));
        fprintf('ModType: %s, Eb/N0 = %d dB, BER = %e\n', modType, EbN0(idx), BER_MRC(modIdx, idx));
    end
end

% Vẽ đồ thị
h = figure;
semilogy(EbN0, BER_MRC(1,:), 'r*-', 'LineWidth', 1.5, 'DisplayName', 'Tx=2,Rx=1,QPSK');
hold on;
semilogy(EbN0, BER_MRC(2,:), 'bs-', 'LineWidth', 1.5, 'DisplayName', 'Tx=2,Rx=1,16QAM');
semilogy(EbN0, BER_MRC(3,:), 'g^-', 'LineWidth', 1.5, 'DisplayName', 'Tx=2,Rx=1,64QAM');
grid on;
set(gca, 'yscale', 'log', 'xlim', [EbN0(1), EbN0(end)], 'ylim', [1e-5 1]);
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('Hệ thống Alamouti 2x1 với các dạng điều chế');
legend('show');
set(h, 'color', [1 1 1]);
set(h, 'NumberTitle', 'off');
set(h, 'renderer', 'zbuffer');
