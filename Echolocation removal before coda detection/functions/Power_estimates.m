function Pkk=Power_estimates(test,locs,F_ds)

    Pkk=[];
    for i=1:length(locs)
        click_cand=test(int32((locs(i)-2e-3)*F_ds):int32((locs(i)+2e-3)*F_ds));
        Pkk(i)=max(click_cand)-min(click_cand);
    end

end