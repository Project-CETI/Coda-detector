
function Slant_Delay=estimate_slant_delays(c_full,F_ds)

    M=2; Slant_Delay=[]; 

        click_cand=c_full;
        Rxx=xcorr(click_cand,click_cand); 
        Rxx=abs(Rxx/max(Rxx));
        [p_ey,l_ey] =findpeaks(Rxx,F_ds,'MinPeakDistance',1.5e-3);
        [p_sort,I_sort]=sort(p_ey,'descend');
        l_sort=l_ey(I_sort);
        p_chosen=p_sort(1:M); l_chosen=l_sort(1:M);
        % figure; plot(Rxx); hold on; plot(l_chosen*F_ds,p_chosen,'x'); 
        Slant_Delay=1e3*abs(diff(l_chosen));
        
end