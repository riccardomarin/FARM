function plotShape(o,points,f)
if not(exist('f','var'))
    f=o.v(:,3);
end
d=0.001;

figure();
 colormap(jet);
 trisurf(o.f,o.v(:,1),o.v(:,2),o.v(:,3),f,'FaceAlpha', 0.6,'EdgeColor', 'none');
 hold on;
 for i=1:length(points)
    scatter3(o.v(points(i),1),o.v(points(i),2), o.v(points(i),3),'filled')
    text(o.v(points(i),1)+d, o.v(points(i),2)+d, o.v(points(i),3)+d,num2str(i));
 end
 hold off;
end