clear all; close all; clc;

% Thông số chung
frmLen     = 1000;             % chiều dài frame
numPackets = 1000;            % số lượng packet
EbN0       = 0:2:20;          % Eb/N0 từ 0 đến 20 dB, bước 2 dB
N          = 4;               % số lượng anten Tx
M          = 1;               % số lượng anten Rx

% Các phương thức điều chế cần so sánh
modTypes = {'BPSK','QPSK','8PSK'};

% Khởi tạo ma trận lưu BER
BER44 = zeros(length(modTypes), length(EbN0));

% Vòng lặp theo từng điều chế
for modIdx = 1:length(modTypes)
    modType = modTypes{modIdx};
    % Xác định P (số điểm constellation) theo điều chế
    switch modType
        case 'BPSK'
            P = 2;
        case 'QPSK'
            P = 4;
        case '8PSK'
            P = 8;
        otherwise
            error('Modulation type not supported');
    end
    
    % Vòng lặp theo từng điểm Eb/N0
    for idx = 1:length(EbN0)
        errors = zeros(1, numPackets);
        
        for pkt = 1:numPackets
            % Tạo dữ liệu ngẫu nhiên và điều chế
            data = randi([0 P-1], frmLen, N);
            X    = pskmod(data, P, 0);  % PSK modulation
            x0 = X(:,1);
            x1 = X(:,1:4);
            x2(:,1) = -X(:,2); x2(:,2) = X(:,1); x2(:,3) = -X(:,4); x2(:,4) = X(:,3);
            x3(:,1) = -X(:,3); x3(:,2) = X(:,4); x3(:,3) = X(:,1); x3(:,4) = -X(:,2);
            x4(:,1) = -X(:,4); x4(:,2) = -X(:,3); x4(:,3) = X(:,2); x4(:,4)= X(:,1);
       
            x5 = conj(x1);
            x6 = conj(x2);
            x7 = conj(x3);
            x8 = conj(x4);

            
            % Sinh kênh Rayleigh (frmLen × Nt × 4)
            Hr = (randn(frmLen, N, 4) + 1i*randn(frmLen, N, 4))/sqrt(2);
            
            % Nhận tín hiệu tại 1 anten Rx
            for n = 1:4
                H = reshape(Hr(:,:,n), frmLen, N);
                % Tính tín hiệu nhận cho 8 bản truyền
                r1(:,n) = awgn(sum(H.*x1,2)/sqrt(N), EbN0(idx));
                r2(:,n) = awgn(sum(H.*x2,2)/sqrt(N), EbN0(idx));
                r3(:,n) = awgn(sum(H.*x3,2)/sqrt(N), EbN0(idx));
                r4(:,n) = awgn(sum(H.*x4,2)/sqrt(N), EbN0(idx));
                r5(:,n) = awgn(sum(H.*x5,2)/sqrt(N), EbN0(idx));
                r6(:,n) = awgn(sum(H.*x6,2)/sqrt(N), EbN0(idx));
                r7(:,n) = awgn(sum(H.*x7,2)/sqrt(N), EbN0(idx));
                r8(:,n) = awgn(sum(H.*x8,2)/sqrt(N), EbN0(idx));
                
                % Quy tụ tín hiệu theo Alamouti giải mã
            z1_1(:,n) = r1(:,n).*conj(H(:,1)) + r2(:,n).*conj(H(:,2)) + r3(:,n).*conj(H(:,3)) + r4(:,n).*conj(H(:,4));
            z1_2(:,n) = conj(r5(:,n)).*H(:,1) + conj(r6(:,n)).*H(:,2) + conj(r7(:,n)).*H(:,3) + conj(r8(:,n)).*H(:,4);
            z1(:,n) = z1_1(:,n) + z1_2(:,n);

            z2_1(:,n) = r1(:,n).*conj(H(:,2)) - r2(:,n).*conj(H(:,1)) - r3(:,n).*conj(H(:,4)) + r4(:,n).*conj(H(:,3));
            z2_2(:,n) = conj(r5(:,n)).*H(:,2) - conj(r6(:,n)).*H(:,1) - conj(r7(:,n)).*H(:,4) + conj(r8(:,n)).*H(:,3);
            z2(:,n) = z2_1(:,n) + z2_2(:,n);

            z3_1(:,n) = r1(:,n).*conj(H(:,3)) + r2(:,n).*conj(H(:,4)) - r3(:,n).*conj(H(:,1)) - r4(:,n).*conj(H(:,2));
            z3_2(:,n) = conj(r5(:,n)).*H(:,3) + conj(r6(:,n)).*H(:,4) - conj(r7(:,n)).*H(:,1) - conj(r8(:,n)).*H(:,2);
            z3(:,n) = z3_1(:,n) + z3_2(:,n);

            z4_1(:,n) = r1(:,n).*conj(H(:,4)) - r2(:,n).*conj(H(:,3)) + r3(:,n).*conj(H(:,2)) - r4(:,n).*conj(H(:,1));
            z4_2(:,n) = conj(r5(:,n)).*H(:,4) - conj(r6(:,n)).*H(:,3) + conj(r7(:,n)).*H(:,2) - conj(r8(:,n)).*H(:,1);
            z4(:,n) = z4_1(:,n) + z4_2(:,n);   
            end
            
            % Tổng hợp tín hiệu từ 4 khung
            r_combined = [sum(z1,2), sum(z2,2), sum(z3,2), sum(z4,2)];
            
            % Giải điều chế
            demod = pskdemod(r_combined, P, 0);
            
            % Tính số lỗi bit
            errors(pkt) = biterr(demod, data);
        end
        
        % Tính BER trung bình
        BER44(modIdx, idx) = sum(errors) / (numPackets * frmLen);
        fprintf('Mod=%s Eb/N0=%2d dB → BER=%e\n', modType, EbN0(idx), BER44(modIdx, idx));
    end
end

% Vẽ đồ thị BER
figure;
semilogy(EbN0, BER44(1,:), 'r*-','LineWidth',1.5','DisplayName','BPSK'); hold on;
semilogy(EbN0, BER44(2,:), 'bs-','LineWidth',1.5','DisplayName','QPSK');
semilogy(EbN0, BER44(3,:), 'g^-','LineWidth',1.5','DisplayName','8PSK');
grid on;
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('BER – Alamouti Tx=4, Rx=1 (BPSK/QPSK/8PSK)');
legend('show');
set(gcf,'color',[1 1 1]);
