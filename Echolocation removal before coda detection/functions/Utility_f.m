function [U,U_pair_max,Transient]=Utility_f(L_tot,T_inds,ToAs_vec)
     
     Transient=[];
     U_pair=0;
     for i=1:length(T_inds)-1
         U_pair(i)=L_tot(T_inds(i),T_inds(i+1));
     end
     [B,TFrm] = rmoutliers(U_pair);
     if sum(TFrm)==2
         U_pair=B;
         El=find(TFrm==1);
         Transient=El(2);
         ToAs_vec(El(2))=[];
     end

     ICI_vec=diff(ToAs_vec);
     Consi=0;
    for c_idx=1:length(ICI_vec)-1
          Consi(c_idx)=Consi_Score(ICI_vec(c_idx),ICI_vec(c_idx+1));
    end

    % U=mean(U_pair)+mean(Consi);
    U=-mean(Consi);
    U_pair_max=max(U_pair);
end