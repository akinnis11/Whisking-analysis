
function minFollicleDistance = findFollicle_b(fx,fy,faceEdgeX,faceEdgeY)

dist = [];
minDist = [];
dist = [];
for j = 1:size(faceEdgeX,2)
    dist(j) = sqrt((faceEdgeX(j) - fx)^2 + (faceEdgeY(j) - fy)^2);
end

[a,b]=min(dist);
minFollicleDistance = a;
end