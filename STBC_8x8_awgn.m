clc; clear;

N = 8; % Số anten phát
M = 8; % Số anten thu
P = 4; % QPSK
EbN0_dB = 0:2:20;
frmLen = 1000; % Số frame để tính BER
BER = zeros(1, length(EbN0_dB));

for idx = 1:length(EbN0_dB)
    errors = 0; total = 0;
    snr = EbN0_dB(idx);
    
    for pkt = 1:frmLen
        %Tạo dữ liệu và điều chế QPSK
        data = randi([0 P-1], N, 1); % 8 symbol
        x = pskmod(data, P, pi/4);  % QPSK
        
        %Mã hóa STBC theo G8
        X_coded = zeros(16, 8); % 16 time slot × 8 anten
        % Time slot 1–8
        X_coded(1,:)  = [x(1)  x(2)  x(3)  x(4)  x(5)  x(6)  x(7)  x(8)];
        X_coded(2,:)  = [-x(2) x(1)  x(4) -x(3)  x(6) -x(5) -x(8)  x(7)];
        X_coded(3,:)  = [-x(3) -x(4) x(1)  x(2)  x(7)  x(8) -x(5) -x(6)];
        X_coded(4,:)  = [-x(4) x(3) -x(2) x(1)  x(8) -x(7)  x(6) -x(5)];
        X_coded(5,:)  = [-x(5) -x(6) -x(7) -x(8) x(1)  x(2)  x(3)  x(4)];
        X_coded(6,:)  = [-x(6) x(5) -x(8) x(7) -x(2) x(1) -x(4)  x(3)];
        X_coded(7,:)  = [-x(7) x(8) x(5) -x(6) -x(3)  x(4) x(1) -x(2)];
        X_coded(8,:)  = [-x(8) -x(7) x(6) x(5) -x(4) -x(3)  x(2)  x(1)];
        % Time slot 9–16 (conjugate phần tử)
        for i = 1:8
            X_coded(8+i,:) = conj(X_coded(i,:));
        end
        
        % Kênh truyền MIMO Rayleigh (8x8)
        H = (randn(N, M) + 1i*randn(N, M)) / sqrt(2); % N×M

        for n = 1:16
            %ma tran phat
            H = reshape(Hr(:,:,t),frmLen,N);
            %tin hieu tai anten thu
            r1(:,n)=awgn(sum(H.*X_coded(1,:),2)/sqrt(N),EbN0(idx));
            r2(:,n)=awgn(sum(H.*X_coded(2,:),2)/sqrt(N),EbN0(idx));
            r3(:,n)=awgn(sum(H.*X_coded(3,:),2)/sqrt(N),EbN0(idx));
            r4(:,n)=awgn(sum(H.*X_coded(4,:),2)/sqrt(N),EbN0(idx));
            r5(:,n)=awgn(sum(H.*X_coded(5,:),2)/sqrt(N),EbN0(idx));
            r6(:,n)=awgn(sum(H.*X_coded(6,:),2)/sqrt(N),EbN0(idx));
            r7(:,n)=awgn(sum(H.*X_coded(7,:),2)/sqrt(N),EbN0(idx));
            r8(:,n)=awgn(sum(H.*X_coded(8,:),2)/sqrt(N),EbN0(idx));
            r9(:,n)=awgn(sum(H.*X_coded(9,:),2)/sqrt(N),EbN0(idx));
            r10(:,n)=awgn(sum(H.*X_coded(10,:),2)/sqrt(N),EbN0(idx));
            r11(:,n)=awgn(sum(H.*X_coded(11,:),2)/sqrt(N),EbN0(idx));
            r12(:,n)=awgn(sum(H.*X_coded(12,:),2)/sqrt(N),EbN0(idx));
            r13(:,n)=awgn(sum(H.*X_coded(13,:),2)/sqrt(N),EbN0(idx));
            r14(:,n)=awgn(sum(H.*X_coded(14,:),2)/sqrt(N),EbN0(idx));
            r15(:,n)=awgn(sum(H.*X_coded(15,:),2)/sqrt(N),EbN0(idx));
            r16(:,n)=awgn(sum(H.*X_coded(16,:),2)/sqrt(N),EbN0(idx));

            


        end
        % Thêm noise
        Y = awgn(Y, snr, 'measured');

        % 4. Giải mã MRC
        x_hat = zeros(N,1); % ước lượng symbol
        for k = 1:N
            temp = 0;
            for m = 1:M
                temp = temp + H(k,m)' * sum(Y(:,m));
            end
            x_hat(k) = temp;
        end
        
        % 5. Giải điều chế
        data_hat = pskdemod(x_hat, P, pi/4);
        errors = errors + sum(data ~= data_hat);
        total = total + N;
    end
    
    BER(idx) = errors / total;
    fprintf("Eb/N0 = %d dB => BER = %.4e\n", snr, BER(idx));
end

% Vẽ biểu đồ BER
semilogy(EbN0_dB, BER, '-o', 'LineWidth', 2);
xlabel('Eb/N0 (dB)'); ylabel('Bit Error Rate');
title('STBC MIMO 8x8 sử dụng MRC');
grid on;
