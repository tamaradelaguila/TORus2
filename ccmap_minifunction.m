function ccmap = load_BV(m)
BVcolor = [0.0333 0.25 0;
    0.0667 0.5 0;
    0.1 0.75 0.015;
    0.97 1 0.0153;
    0.9647 0.815 0.0155;
    0.9608 0.63 0.03;
    0.957 0.445 0.0447;
    0.958 0.26 0.06;
    1 0.227 0.537]; 

l = length(BVcolor);
cindex = linspace(1,9,m);
r = interp1([1:l],BVcolor(:,1),cindex);
g = interp1([1:l],BVcolor(:,2),cindex);
b = interp1([1:l],BVcolor(:,3),cindex);
ccmap = [r' g' b'];
end