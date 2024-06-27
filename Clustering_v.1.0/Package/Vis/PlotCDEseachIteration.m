function h = PlotCDEseachIteration(Data,v,label,param)
c = param.kClust;
fm = param.mFuzzy;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Visualization
lower1=min(Data(1,:));upper1=max(Data(1,:));
lower2=min(Data(2,:));upper2=max(Data(2,:));

[x,y] = meshgrid(linspace(lower1, upper1, 200), linspace(lower2, upper2, 200));
pair = [x(:) y(:)];
[pair1,~] = size(pair);

for i = 1 : c
    xv = pair - ones(pair1, 1)*v(i,:);
    d(:,i)= sum((xv*eye(2).*xv),2);%
end
d = (d+1e-10).^(-1/(fm-1));
f0 = (d ./ (sum(d,2)*ones(1,c)));

f = max(f0')';
Z = reshape(f,size(x,1),size(x,2));


h = figure;

contour(x,y,Z, "DisplayName", "Contour"); hold on; 

numCluster = max(label);
colors = lines(numCluster);  
colororder("meadow")

for i = 1:numCluster
    clusterData = Data(:,label == i);
    scatter(clusterData(1,:), clusterData(2,:), 18, colors(i,:), 'filled');  % 36 is the marker size
end

scatter(v(:,1), v(:,2), 72, colors(i,:), 'filled', "k", "DisplayName", "Centroids");  % 36 is the marker size

end