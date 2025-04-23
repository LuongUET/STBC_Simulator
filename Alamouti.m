clear all; close all; clc;

frmLen = 1000; % chiều dài frame (số ký hiệu)
numPackets = 1000; % số lượng packet
EbN0 = 0:2:20; % giá trị Eb/N0 từ 0 đến 20 dB, bước 2 dB
N = 2; % số lượng anten Tx (cho Alamouti)
M = 1; % số lượng anten thu Rx
P = 4; % độ dài điều chế QPSK

% Khởi tạo mảng lưu trữ BER
BERMRC = zeros(1, length(EbN0)); % BER cho Alamouti 2x1 MRC
BERSISO = zeros(1, length(EbN0)); % BER cho SISO 1x1
errorMRC = zeros(1, numPackets); % Lỗi cho Alamouti MRC
errorSISO = zeros(1, numPackets); % Lỗi cho SISO

% Vòng lặp cho Eb/N0 điểm
for idx = 1:length(EbN0)
    for packetIdx = 1:numPackets
        % Tạo dữ liệu QPSK
        data = randi([0 P-1], frmLen, 2); % Dữ liệu QPSK (frmLen x 2)
        X = pskmod(data, P, 0); % Điều chế QPSK

        % --- Mô phỏng SISO 1x1 ---
        % Tạo dữ liệu QPSK cho SISO (dùng cột 1 của dữ liệu gốc)
        dataSISO = data(:,1); % Dữ liệu cho SISO (frmLen x 1)
        xSISO = pskmod(dataSISO, P, 0); % Điều chế QPSK

        % Kênh Rayleigh fading (SISO: 1 Tx, 1 Rx)
        hSISO = (randn(frmLen, 1) + 1i*randn(frmLen, 1)) / sqrt(2);

        % Truyền tín hiệu (SISO)
        rSISO = hSISO.*xSISO; % Tín hiệu nhận được
        rSISO = awgn(rSISO, EbN0(idx) + 3, 'measured'); % Thêm nhiễu

        % Giải mã SISO (bù kênh đơn giản)
        xHat_SISO = conj(hSISO).*rSISO./(abs(hSISO).^2);
        decodedData_SISO = pskdemod(xHat_SISO, P, 0);

        % Tính lỗi cho SISO
        errorSISO(packetIdx) = biterr(decodedData_SISO, dataSISO);


        % --- Mô phỏng Alamouti 2x1 ---
        % Mã hóa Alamouti
        x1 = X(:,1); % Ký hiệu x1 (frmLen x 1)
        x2 = X(:,2); % Ký hiệu x2 (frmLen x 1)

        % Kênh Rayleigh fading (Alamouti)
        Hr = (randn(frmLen, N, M) + 1i*randn(frmLen, N, M)) / sqrt(2);
        h1 = Hr(:,1,1); % Kênh từ anten 1
        h2 = Hr(:,2,1); % Kênh từ anten 2

        % Truyền tín hiệu (Alamouti)
        s1 = h1.*x1 + h2.*x2; % Tín hiệu tại khe thời gian 1
        s2 = -h1.*conj(x2) + h2.*conj(x1); % Tín hiệu tại khe thời gian 2
        r1 = awgn(s1, EbN0(idx) + 3, 'measured'); % Thêm nhiễu (QPSK: Es/N0 = Eb/N0 + 3 dB)
        r2 = awgn(s2, EbN0(idx) + 3, 'measured');

        % Giải mã MRC (Alamouti)
        norm = abs(h1).^2 + abs(h2).^2;
        x1HatAla21 = (conj(h1).*r1 + h2.*conj(r2))./norm;
        x2HatAla21 = (conj(h2).*r1 - h1.*conj(r2))./norm;

        % Gộp tín hiệu MRC (Alamouti)
        decodedSymbols_Ala21 = [x1HatAla21 x2HatAla21];
        decodedData_Ala21 = pskdemod(decodedSymbols_Ala21, P, 0);

        % Tính lỗi cho Alamouti MRC
        errorAla21(packetIdx) = biterr(decodedData_Ala21, data);
        
        % --- Mô phỏng Alamouti 2x1 ---
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

            decodedData_Ala22 = pskdemod(z, P, 0); % Giải điều chế QPSK
            errorAla22(packetIdx) = biterr(decodedData_Ala22, data);
    end

    % Tính BER (QPSK: 2 bit/ký hiệu)
    BERMRC21(idx) = sum(errorAla21) / (numPackets * frmLen * 2 * 2); % Alamouti: 2 cột dữ liệu
    BERMRC22(idx) = sum(errorAla22) / (numPackets * frmLen * log2(P));
    BERSISO(idx) = sum(errorSISO) / (numPackets * frmLen * 2); % SISO: 1 cột dữ liệu
    fprintf('Eb/N0 = %d dB, BER Alamouti = %e, BER SISO = %e, BER SISO = %e\n', EbN0(idx), BERMRC21(idx), BERMRC22(idx),BERSISO(idx));
end

% Vẽ đồ thị
h = figure;
semilogy(EbN0, BERMRC21, 'r*-', 'LineWidth', 1.5, 'DisplayName', 'Alamouti 2x1'); % Đường Alamouti
hold on;
semilogy(EbN0, BERMRC22, 'bs-', 'LineWidth', 1.5, 'DisplayName', 'Alamouti 2x2'); % Đường Alamouti
hold on;
semilogy(EbN0, BERSISO, 'g^-', 'LineWidth', 1.5, 'DisplayName', 'SISO 1x1'); % Đường SISO
grid on;
set(gca, 'yscale', 'log', 'xlim', [EbN0(1), EbN0(end)], 'ylim', [1e-5 1]);
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('Mô phỏng mã hóa Alamouti');
legend('show');
set(h, 'color', [1 1 1]);
set(h, 'NumberTitle', 'off');
set(h, 'renderer', 'zbuffer');