
function percCounter(i,ds,n)

%Keeps track of how much farther to go
percDone = floor(100*(i/(n)));
percDoneLast = floor(100*((i-ds)/n));

if isequal(percDone,percDoneLast) == 0
    disp(sprintf('%s',num2str(percDone),'% done...'))
end
end
