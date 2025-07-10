function D=spectral_similarity(c1,c2,Fs,plot_flag)

    if plot_flag
        figure; 
    end
    minLen=min([length(c1) length(c2)]);
    C={c1(1:minLen),c2(1:minLen)};
    for j=1:2
        Sig=C{j};
        Y = fft(Sig);
        L=length(Sig);
        P2 = abs(Y/L);
        P1 = P2(int32(1):int32(L/2+1));
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        Pm=movmean(P1,12); 
        if plot_flag
           plot(normalize(Pm)); hold on;
        end
        P_w(j,:)=normalize(Pm);
    end
    D = norm(P_w(2,:) - P_w(1,:));
end


