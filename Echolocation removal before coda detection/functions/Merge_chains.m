function Detected_subtrains=Merge_chains(Detected_subtrains,locs)
    
    while(1)
        c=0; Merge_idx=[]; Merge_flag=0;
        if length(Detected_subtrains)>1
           comb=nchoosek(1:length(Detected_subtrains),2);
        else
            break;
        end
        for q=1:size(comb,1)
            a=locs(Detected_subtrains{comb(q,1)});
            b=locs(Detected_subtrains{comb(q,2)});
          
            cand_ICI_a=diff(sort(a));
            cand_ICI_b=diff(sort(b));
            if size(a,1)>1
                cand_ICI_merged=diff(sort([a' b']));
            else
                cand_ICI_merged=diff(sort([a b]));
            end

            % if length(rmoutliers(diff(sort([a b]))))<4
            %    cand_ICI_merged=rmoutliers(diff(sort([a b])));
            % else
            %     cand_ICI_merged=diff(sort([a b]));
            % end
            Consi_a=[]; Consi_b=[]; Consi_merged=[];
            Consi_merged_mean=[]; Consi_separated_mean=[];
            if length(cand_ICI_a)>1
                for j=1:length(cand_ICI_a)-1
                    Consi_a(j)=log(cand_ICI_a(j+1)/cand_ICI_a(j));
                end
            end
            if length(cand_ICI_b)>1       
                for j=1:length(cand_ICI_b)-1
                    Consi_b(j)=log(cand_ICI_b(j+1)/cand_ICI_b(j));
                end
            end
            for j=1:length(cand_ICI_merged)-1
                Consi_merged(j)=log(cand_ICI_merged(j+1)/cand_ICI_merged(j));
            end  
            Consi_merged_mean=mean(abs(Consi_merged));
            Consi_separated_mean=[Consi_a Consi_b];
            if max(abs(Consi_merged))<0.16
                c=c+1;
                Merge_idx(c,:)=[comb(q,:) Consi_merged_mean];
                Merge_flag=1;
            end

            % if ~isempty(Consi_separated_mean)
            %     if Consi_merged_mean<mean(abs(Consi_separated_mean)) || max(abs(Consi_merged))<0.16
            %         c=c+1;
            %         Merge_idx(c,:)=[comb(q,:) Consi_merged_mean];
            %         Merge_flag=1;
            %     end
            % elseif Consi_merged_mean<0.17
            %         c=c+1;
            %         Merge_idx(c,:)=[comb(q,:) Consi_merged_mean];
            %         Merge_flag=1;
            % end

         end
       
        if ~Merge_flag
            break;
        end
    
        [~,min_idx]=min(Merge_idx(:,3));
        Combine=Merge_idx(min_idx,1:2);
        Remain = setxor(1:length(Detected_subtrains),Combine);
        Detected_subtrains=[Detected_subtrains(Remain) {sort(cell2mat(Detected_subtrains(Combine)'))}];
    end 

end
 

