function simu_TH2_STBC3x4(frmLen,numPackets,EbN0)

    BER41 = ber_mimo_3x4(frmLen, numPackets, EbN0, 3, 4, 4);
    BER42 = ber_mimo_3x4(frmLen, numPackets, EbN0, 3, 4, 16);
    BER43 = ber_mimo_3x4(frmLen, numPackets, EbN0, 3, 4, 64);

    h = figure;
    semilogy(EbN0, BER41, 'r*-', 'LineWidth', 1, 'DisplayName', 'Tx=3, Rx=4, 4QAM');
    hold on;
    semilogy(EbN0, BER42, 'bs-', 'LineWidth', 1, 'DisplayName', 'Tx=3, Rx=4, 16QAM');
    hold on;
    semilogy(EbN0, BER43, 'g^-', 'LineWidth', 1, 'DisplayName', 'Tx=3, Rx=4, 64QAM');
    grid on;

    set(gca, 'yscale', 'log', 'xlim', [EbN0(1), EbN0(end)], 'ylim', [1e-8 1]);
    xlabel('Eb/N0 (dB)');
    ylabel('Bit Error Rate (BER)');
    title('BER vs SNR for STBC 3x4 with QAM');
    legend('show');
    set(h, 'Color', [1 1 1]);
    set(h, 'NumberTitle', 'off');
    set(h, 'Renderer', 'zbuffer');
end

function BER = ber_mimo_3x4(frmLen, numPackets, EbN0, N, M, P)
    BER = zeros(size(EbN0)); 

    for idx = 1:length(EbN0)
        errors = zeros(1, numPackets); 
        for pkt = 1:numPackets
            
            data = randi([0 3], frmLen, 4);
            X = qammod(data, P, 'UnitAveragePower', true);

            x1 = [ X(:,1),  X(:,2),  X(:,3)];
            x2 = [-X(:,2),  X(:,1), -X(:,4)];
            x3 = [-X(:,3),  X(:,4),  X(:,1)];
            x4 = [-X(:,4), -X(:,3),  X(:,2)];
            x5 = conj(x1);
            x6 = conj(x2);
            x7 = conj(x3);
            x8 = conj(x4);

            Hr = (randn(frmLen, N, M) + 1i * randn(frmLen, N, M)) / sqrt(2);
            SNR_dB = EbN0(idx) + 10*log10(log2(P));
          
            for n = 1:M
                H = reshape(Hr(:,:,n), frmLen, N);

                r1(:,n) = awgn(sum(H .* x1, 2) / sqrt(N), SNR_dB,'measured');
                r2(:,n) = awgn(sum(H .* x2, 2) / sqrt(N), SNR_dB,'measured');
                r3(:,n) = awgn(sum(H .* x3, 2) / sqrt(N), SNR_dB,'measured');
                r4(:,n) = awgn(sum(H .* x4, 2) / sqrt(N), SNR_dB,'measured');
                r5(:,n) = awgn(sum(H .* x5, 2) / sqrt(N), SNR_dB,'measured');
                r6(:,n) = awgn(sum(H .* x6, 2) / sqrt(N), SNR_dB,'measured');
                r7(:,n) = awgn(sum(H .* x7, 2) / sqrt(N), SNR_dB,'measured');
                r8(:,n) = awgn(sum(H .* x8, 2) / sqrt(N), SNR_dB,'measured');

                % STBC decoding
                z1_1(:,n) = r1(:,n).*conj(H(:,1)) + r2(:,n).*conj(H(:,2)) + r3(:,n).*conj(H(:,3));
                z1_2(:,n) = conj(r5(:,n)).*H(:,1) + conj(r6(:,n)).*H(:,2) + conj(r7(:,n)).*H(:,3);
                z1(:,n) = z1_1(:,n) + z1_2(:,n);

                z2_1(:,n) = r1(:,n).*conj(H(:,2)) - r2(:,n).*conj(H(:,1)) + r4(:,n).*conj(H(:,3));
                z2_2(:,n) = conj(r5(:,n)).*H(:,2) - conj(r6(:,n)).*H(:,1) + conj(r8(:,n)).*H(:,3);
                z2(:,n) = z2_1(:,n) + z2_2(:,n);

                z3_1(:,n) = r1(:,n).*conj(H(:,3)) - r3(:,n).*conj(H(:,1)) - r4(:,n).*conj(H(:,2));
                z3_2(:,n) = conj(r5(:,n)).*H(:,3) - conj(r7(:,n)).*H(:,1) - conj(r8(:,n)).*H(:,2);
                z3(:,n) = z3_1(:,n) + z3_2(:,n);

                z4_1(:,n) = -r2(:,n).*conj(H(:,3)) + r3(:,n).*conj(H(:,2)) - r4(:,n).*conj(H(:,1));
                z4_2(:,n) = -conj(r6(:,n)).*H(:,3) + conj(r7(:,n)).*H(:,2) - conj(r8(:,n)).*H(:,1);
                z4(:,n) = z4_1(:,n) + z4_2(:,n);
            end

            norm_factor = sum(sum(abs(Hr).^2, 3), 2);
            r_comb = [sum(z1,2), sum(z2,2), sum(z3,2), sum(z4,2)] ./ norm_factor;
            demod = qamdemod(r_comb, P, 'UnitAveragePower', true);
            errors(pkt) = biterr(demod, data);
        end

        BER(idx) = sum(errors) / (numPackets * frmLen * log2(P) * 4);
        fprintf('MIMO STBC 3X4 Tx=%d, Rx=%d, %d-QAM: Eb/N0=%2.1f dB → BER=%e\n', N, M,P, EbN0(idx), BER(idx));
    end
end