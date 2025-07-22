function Pk=Peaks_extract(test,locs,F_ds)

    Pkk=[];
    for i=1:length(locs)
        click_cand=test(int32((locs(i)-2e-3)*F_ds):int32((locs(i)+2e-3)*F_ds));
        Pk(i)=max(click_cand);
    end

end