function AllCombs_mat=mat2inds(AllCombs,N)

    AllCombs_mat=zeros(size(AllCombs,1),N);
    for i=1:size(AllCombs,1)
        AllCombs_mat(i,AllCombs(i,:))=1;
    end

end


