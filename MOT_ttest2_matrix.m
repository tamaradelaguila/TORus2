for condi = 1:4
    
    % SELECT TRIALS
for triali = makeRow(seltrials)
    for roii = 1:3
       mats.peak(triali, roii, condii) = ;
    end
end 

end

for condi = 1:4
    
    % SELECT TRIALS
for triali = makeRow(seltrials)
    for roi1 = 1:3
        for roi2 = 1:3
       values1=  mats.peak(triali, roii, condii);
       values2=  mats.peak(triali, roii, condii);
       
       tresults(roi1,roi2, condi) = ttest2(values1, values2);
       
        end
    end
end 

end

imagesc(:,:,condi);