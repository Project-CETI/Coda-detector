function [alpha_Slant_Delay,alpha_Pkk,DTW,Membership]=extract_click_pair_features_V3(F_weights,ICI,c1,c2,c1_full,c2_full,F_ds,All_objs,Tag_flag)
        
        Pkk_left=max(c1)-min(c1);
        Pkk_right=max(c1)-min(c2);
       
        Slant_Delay_1=estimate_slant_delays(c1_full,F_ds);
        Slant_Delay_2=estimate_slant_delays(c2_full,F_ds);
        DTW=dtw(c1/max(c1),c2/max(c2));
        % alpha_Pkk=(Pkk_right-Pkk_left)/max([Pkk_right Pkk_left]);
        % alpha_Slant_Delay=(Slant_Delay_1-Slant_Delay_2)/max([Slant_Delay_1 Slant_Delay_2]);
        alpha_Pkk=log(Pkk_right/Pkk_left);
        alpha_Slant_Delay=log(Slant_Delay_1/Slant_Delay_2);
        Spectral=spectral_similarity(c1,c2,F_ds,0);

        Features=[alpha_Slant_Delay alpha_Pkk DTW Spectral ICI];
        for q=1:4
            objA=All_objs(q).objs;
            N_val=All_objs(q).N_val;
            Mem_val(q)=pdf(objA,Features(q))/N_val;
        end
       % Membership=sum(Mem_val);
       Membership=sum(F_weights.*Mem_val);

       if Tag_flag
           Membership=sum((1/3)*Mem_val(2:4));
       end
       
end

