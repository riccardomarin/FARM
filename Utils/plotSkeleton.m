function plotSkeleton(giunti,KT,label,c,l)
d=0.01;
if nargin==4
    l=24;
end

for i=2:l
    point=[giunti(i,:);giunti(KT(i),:)];
    plot3(point(:,1),point(:,2) ,point(:,3),c );
    t=[num2str(i) label];
    text(point(1,1)+d, point(1,2)+d, point(1,3)+d,t);
hold on;
end
t=[num2str(1) label];
text(giunti(1,1)+d, giunti(1,2)+d, giunti(1,3)+d,t);
axis equal
%axis([-1.5 1.5 -3 3 -3 3]);
end