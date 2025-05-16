function simu_TH4_STBC4x4(frmLen,numPackets,EbN0)
        
        P = 2; 
       
        BER41 = ber_mimo_4x4(frmLen, numPackets, EbN0, 4, 1, P);
        BER42 = ber_mimo_4x4(frmLen, numPackets, EbN0, 4, 2, P);
        BER43 = ber_mimo_4x4(frmLen, numPackets, EbN0, 4, 3, P);
        BER44 = ber_mimo_4x4(frmLen, numPackets, EbN0, 4, 4, P);
        
        h = figure;
        semilogy(EbN0, BER41, 'r*-', 'LineWidth', 1, 'DisplayName', 'Tx=4, Rx=1, BPSK'); hold on;
        semilogy(EbN0, BER42, 'bs-', 'LineWidth', 1, 'DisplayName', 'Tx=4, Rx=2, BPSK'); hold on;
        semilogy(EbN0, BER43, 'g^-', 'LineWidth', 1, 'DisplayName', 'Tx=4, Rx=3, BPSK'); hold on;
        semilogy(EbN0, BER44, 'k>-', 'LineWidth', 1, 'DisplayName', 'Tx=4, Rx=4, BPSK'); hold on;
        grid on;
        set(gca, 'yscale', 'log', 'xlim', [EbN0(1), EbN0(end)], 'ylim', [1e-7 1]);
        xlabel('Eb/N0 (dB)');
        ylabel('Bit Error Rate (BER)');
        title('BER vs SNR for STBC 4x4 with BPSK');
        legend('show');
        set(h, 'Color', [1 1 1]);
        set(h, 'NumberTitle', 'off');
        set(h, 'Renderer', 'zbuffer');
end

function BER = ber_mimo_4x4(frmLen, numPackets, EbN0, N, M, P)

    BER = zeros(size(EbN0));

    for idx = 1:length(EbN0)
        errors = zeros(1, numPackets);
        for pkt = 1:numPackets
            data = randi([0 P-1], frmLen, N);
            X    = pskmod(data, P, 0);

            x1 = X(:,1:4);
            x2 = [-X(:,2),  X(:,1), -X(:,4),  X(:,3)];
            x3 = [-X(:,3),  X(:,4),  X(:,1), -X(:,2)];
            x4 = [-X(:,4), -X(:,3),  X(:,2),  X(:,1)];

            Hr = (randn(frmLen, N, M) + 1i*randn(frmLen, N, M)) / sqrt(2);
            SNR_dB = EbN0(idx) + 10*log10(log2(P));

            for n = 1:M
                H = reshape(Hr(:,:,n), frmLen, N);

                r1(:,n) = awgn(sum(H .* x1, 2) / sqrt(N), SNR_dB,'measured');
                r2(:,n) = awgn(sum(H .* x2, 2) / sqrt(N), SNR_dB,'measured');
                r3(:,n) = awgn(sum(H .* x3, 2) / sqrt(N), SNR_dB,'measured');
                r4(:,n) = awgn(sum(H .* x4, 2) / sqrt(N), SNR_dB,'measured');

                z1(:,n) = r1(:,n).*conj(H(:,1)) + r2(:,n).*conj(H(:,2)) + r3(:,n).*conj(H(:,3)) + r4(:,n).*conj(H(:,4));
                z2(:,n) = r1(:,n).*conj(H(:,2)) - r2(:,n).*conj(H(:,1)) - r3(:,n).*conj(H(:,4)) + r4(:,n).*conj(H(:,3));
                z3(:,n) = r1(:,n).*conj(H(:,3)) + r2(:,n).*conj(H(:,4)) - r3(:,n).*conj(H(:,1)) - r4(:,n).*conj(H(:,2));
                z4(:,n) = r1(:,n).*conj(H(:,4)) - r2(:,n).*conj(H(:,3)) + r3(:,n).*conj(H(:,2)) - r4(:,n).*conj(H(:,1));
            end

           
            norm_factor = sum(sum(abs(Hr).^2, 3), 2);
            r_comb = [sum(z1,2), sum(z2,2), sum(z3,2), sum(z4,2)] ./ norm_factor;

            demod = pskdemod(r_comb, P, 0);
            errors(pkt) = biterr(demod, data);
        end

        BER(idx) = sum(errors) / (numPackets * frmLen * log2(P) * 4);
        fprintf('MIMO STBC 4x4 Tx=%d, Rx=%d: Eb/N0=%.1f dB → BER=%e\n', N, M, EbN0(idx), BER(idx));
    end
end
