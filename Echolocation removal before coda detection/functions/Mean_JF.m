function Confidence=Mean_JF(Observations,buffer,F_ds,Buffer_length,t_index)

% Observations=Detections(t_index);
     W=3e-3;
    
    for i=1:length(Observations.Pkk)       
        locs=Observations.ToAs{i}-(t_index-1)*Buffer_length;
        Pkk_vec=Observations.Pkk{i};
        ref_vec=Observations.ref{i};
        alpha_DTW=[];
        alpha_Spectral=[];
        alpha_Pkk=[];
        alpha_Slant_Delay=[];
        for j=1:length(Observations.Pkk{i})-1
                if locs(j)-W<0
                    init=1;
                else
                    init=int32((locs(j)-W)*F_ds);
                end
                if locs(j+1)+W>length(buffer)
                    final=length(buffer);
                else
                    final=int32((locs(j+1)+W)*F_ds);
                end
                Click_current=buffer(init:int32((locs(j)+W)*F_ds));
                Click_next=buffer(int32((locs(j+1)-W)*F_ds):final);
                alpha_DTW(j)=dtw(Click_current/max(Click_current),Click_next/max(Click_next));
                alpha_Spectral(j)=spectral_similarity(Click_current,Click_next,F_ds,0);
                alpha_Pkk(j)=abs(log(Pkk_vec(j+1)/Pkk_vec(j)));
                alpha_Slant_Delay(j)=abs(log(ref_vec(j+1)/ref_vec(j)));
        end
    
        ICI_class(i)=jains_fairness_index(Observations.ICI{i});
        Pkk_class(i)=jains_fairness_index(alpha_Pkk);
        ref_class(i)=jains_fairness_index(alpha_Slant_Delay);
        DTW_class(i)=jains_fairness_index(alpha_DTW);
        Spectral_class(i)=jains_fairness_index(alpha_Spectral);
    
        Confidence(i)=median([ICI_class Pkk_class ref_class DTW_class Spectral_class]);
    end

end


