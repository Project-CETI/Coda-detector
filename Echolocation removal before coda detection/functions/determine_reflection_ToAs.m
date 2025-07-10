
function [ref_ToAs,ref_scores]=determine_reflection_ToAs(test,locs,Frame1,roi,F_ds,plot_flag)
            
    M=2; ref_ToAs=[]; ref_scores=[]; snr=[];
    for i=1:length(locs)
        if  (locs(i)+roi)*F_ds<length(test)
            click_cand=test(int32((locs(i)-5e-3)*F_ds):int32((locs(i)+roi)*F_ds));
        else
            break;
        end
        Rxx=xcorr(click_cand,click_cand); 
        Rxx=abs(Rxx/max(Rxx));
        % [p_ey,l_ey] =findpeaks(Rxx,F_ds,'MinPeakDistance',1.5e-3);
         [p_ey,l_ey] =findpeaks(Rxx,F_ds,'MinPeakDistance',2.5e-3);

        
        [p_sort,I_sort]=sort(p_ey,'descend');
        l_sort=l_ey(I_sort);
        p_chosen=p_sort(1:M); l_chosen=l_sort(1:M);
        % figure; plot(Rxx); hold on; plot(l_chosen*F_ds,p_chosen,'x'); 
        ref_ToAs(i)=1e3*abs(diff(l_chosen));
        p_chosen=p_chosen/max(p_chosen);
        ref_scores(i)=p_chosen(2);
        if plot_flag
            % plot(locs(i),r(i),'kx'); hold on;
            if ref_scores(i)>0
                plot(Frame1+locs(i),ref_ToAs(i),'kx'); hold on;
            else
                plot(Frame1+locs(i),ref_ToAs(i),'ro'); hold on;
            end
        end
    end
    % ylim([0 roi*1e3]);
     grid on;
end