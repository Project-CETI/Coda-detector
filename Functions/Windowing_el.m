function Inds_filt=Windowing_el(Inds,T_haifa)

        for j=1:length(Inds)
             ToAs=T_haifa(Inds(j),2:end);
             click=ToAs{1};
             for k=1:length(ToAs)-1
                 if ~ismissing(ToAs{k+1})
                     click(k+1)=ToAs{k+1};
                     ICIs=diff(click);                    
                 end
             end
             All_clicks(j)={click};
             All_ICIs(j)={ICIs};
        end

        Inds_el=[];

        for i=1:length(All_clicks)-1

            L_max=max([length(All_ICIs{i}) length(All_ICIs{i+1})]);
            [L_min,min_ind]=min([length(All_ICIs{i}) length(All_ICIs{i+1})]);
            ICI_1=zeros(1,L_max); ICI_2=zeros(1,L_max);
            ICI_1(1:length(All_ICIs{i}))=All_ICIs{i};
            ICI_2(1:length(All_ICIs{i+1}))=All_ICIs{i+1};
            Change=abs(ICI_1-ICI_2);

            for Q=1:L_max-1
                Change=abs(ICI_2-circshift(ICI_1, [0, Q]));
                Sum_max(Q)=sum(Change<5e-3);
            end
            Q_max=find(Sum_max==max(Sum_max));
            Change=abs(ICI_2-circshift(ICI_1, [0, Q_max]));
            
            if sum(Change<5e-3)/L_min>0.7
                if min_ind==1
                    Inds_el=[Inds_el i];
                elseif min_ind==2
                    Inds_el=[Inds_el i+1];
                end
            end
        end

        for i=1:length(All_clicks)-2

            L_max=max([length(All_ICIs{i}) length(All_ICIs{i+2})]);
            [L_min,min_ind]=min([length(All_ICIs{i}) length(All_ICIs{i+2})]);
            ICI_1=zeros(1,L_max); ICI_2=zeros(1,L_max);
            ICI_1(1:length(All_ICIs{i}))=All_ICIs{i};
            ICI_2(1:length(All_ICIs{i+2}))=All_ICIs{i+2};

            for Q=1:L_max-1
                Change=abs(ICI_2-circshift(ICI_1, [0, Q]));
                Sum_max(Q)=sum(Change<5e-3);
            end
            Q_max=find(Sum_max==max(Sum_max));
            Change=abs(ICI_2-circshift(ICI_1, [0, Q_max]));
            
            if sum(Change<5e-3)/L_min>0.7
                if min_ind==1
                    Inds_el=[Inds_el i];
                elseif min_ind==2
                    Inds_el=[Inds_el i+2];
                end
            end
        end

       Inds_filt=Inds;
       Inds_filt(Inds_el)=[];

end