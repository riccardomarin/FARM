%%%%%
% Code for article:
% Marin, R. and Melzi, S. and Rodolà, E. and Castellani, U., FARM: Functional Automatic Registration Method for 3D Human Bodies, CGF 2019
% Github: https://github.com/riccardomarin/FARM/
% Project Page: http://profs.scienze.univr.it/~marin/farm/index.html
%%%%%

function landmarks = SymBreaker(landmarks_actual, Joints_Target, KT,print)

syms d a b c x y z
landmarks=landmarks_actual;

v1=Joints_Target(12,:)-Joints_Target(9,:);
p1=v1(1)*(x-Joints_Target(9,1))+v1(2)*(y-Joints_Target(9,2))+v1(3)*(z-Joints_Target(9,3));
n1=v1;


po1=Joints_Target(9,:);
po2=Joints_Target(6,:);
po3=Joints_Target(12,:);
vec1=po1-po2;
vec2=po3-po2;
vec_n=cross(vec1,vec2);
p2=vec_n(1)*(x-Joints_Target(9,1))+vec_n(2)*(y-Joints_Target(9,2))+vec_n(3)*(z-Joints_Target(9,3));
n2=vec_n;

c = double(fliplr(coeffs(p2)));
po1=Joints_Target(9,:);
po2=Joints_Target(6,:);
po3=Joints_Target(9,:)+c([1,2,3]);
vec1=po1-po2;
vec2=po3-po2;
vec_n=cross(vec1,vec2);
p3=vec_n(1)*(x-Joints_Target(9,1))+vec_n(2)*(y-Joints_Target(9,2))+vec_n(3)*(z-Joints_Target(9,3));
n3=vec_n;

po1=Joints_Target(3,:);
po2=Joints_Target(6,:);
po3=po2+c([1,2,3]);
vec1=po1-po2;
vec2=po3-po2;
vec_n=cross(vec1,vec2);
p4=vec_n(1)*(x-Joints_Target(6,1))+vec_n(2)*(y-Joints_Target(6,2))+vec_n(3)*(z-Joints_Target(6,3));
n4=vec_n;

c = double(fliplr(coeffs(p4)));
po1=Joints_Target(3,:);
po2=Joints_Target(6,:);
po3=Joints_Target(3,:)+c([1,2,3]);
vec1=po1-po2;
vec2=po3-po2;
vec_n=cross(vec1,vec2);
p5=vec_n(1)*(x-po1(1))+vec_n(2)*(y-po1(2))+vec_n(3)*(z-po1(3));
n5=vec_n;

c = double(fliplr(coeffs(p5)));
po1=Joints_Target(1,:);
po2=Joints_Target(3,:);
po3=Joints_Target(3,:)+c([1,2,3]);
vec1=po1-po2;
vec2=po3-po2;
vec_n=cross(vec1,vec2);
p6=vec_n(1)*(x-po1(1))+vec_n(2)*(y-po1(2))+vec_n(3)*(z-po1(3));
n6=vec_n;

p1_explicit=solve(p1,z);
p2_explicit=solve(p2,z);
p3_explicit=solve(p3,z);
p4_explicit=solve(p4,z);
p5_explicit=solve(p5,z);
p6_explicit=solve(p6,z);

c = double(fliplr(coeffs(p2)));
po1=Joints_Target(9,:);
po3=Joints_Target(6,:)+c([1,2,3]);
po2=Joints_Target(6,:);
P=[po1;po2;po3];
T=[1 2 3];
TR = triangulation(T,P);
FN = faceNormal(TR);
N=double(fliplr(coeffs(p3)));
N=N([1,2,3]);
N=N/norm(N);

po2=Joints_Target(3,:);
po3=Joints_Target(6,:)+c([1,2,3]);
po1=Joints_Target(6,:);
P=[po1;po2;po3];
T=[1 2 3];
TR = triangulation(T,P);
FN2 = faceNormal(TR);
N2=double(fliplr(coeffs(p4)));
N2=N2([1,2,3]);
N2=N/norm(N);

c = double(fliplr(coeffs(p5)));

po2=Joints_Target(3,:);
po3=Joints_Target(3,:)+c([1,2,3]);
po1=Joints_Target(1,:);
P=[po1;po2;po3];
T=[1 2 3];
TR = triangulation(T,P);
FN3 = faceNormal(TR);
N3=double(fliplr(coeffs(p4)));
N3=N3([1,2,3]);
N3=N/norm(N); 

FN=FN/10;
FN2=FN2/10;
FN3=FN3/10;


if nargin==4

FS=n1/norm(n1)/10
  subplot(2,4,[1 5])
  plotSkeleton(Joints_Target,KT,'','b');
  span=[Joints_Target(9,1)-0.1 Joints_Target(9,1)+0.1 Joints_Target(9,2)-0.1 Joints_Target(9,2)+0.1];
  fsurf(p1_explicit,span,'FaceAlpha',0.5)
  quiver3(Joints_Target(9,1),Joints_Target(9,2),Joints_Target(9,3),FS(1),FS(2),FS(3),1,'LineWidth',2)
  

subplot(2,4,[2 6])
 plotSkeleton(Joints_Target,KT,'','b');
 span=[Joints_Target(6,1)-0.1 Joints_Target(6,1)+0.1 Joints_Target(6,2)-0.1 Joints_Target(6,2)+0.1];
 fsurf(p3_explicit,span,'FaceAlpha',0.5)
 quiver3(Joints_Target(6,1),Joints_Target(6,2),Joints_Target(6,3),FN(1),FN(2),FN(3),2,'LineWidth',2)

subplot(2,4,[3 7])
 plotSkeleton(Joints_Target,KT,'','b');
 span=[Joints_Target(6,1)-0.05 Joints_Target(6,1)+0.05 Joints_Target(6,2)-0.1 Joints_Target(6,2)+0.1];
 fsurf(p4_explicit,span,'FaceAlpha',0.5)
 quiver3(Joints_Target(6,1),Joints_Target(6,2),Joints_Target(6,3),FN2(1),FN2(2),FN2(3),2,'LineWidth',2)
subplot(2,4,[4 8])
 %figure;
 plotSkeleton(Joints_Target,KT,'','b');
 span=[Joints_Target(3,1)-0.05 Joints_Target(3,1)+0.05 Joints_Target(3,2)-0.05 Joints_Target(3,2)+0.05];
 fsurf(p6_explicit,span,'FaceAlpha',0.5)
 quiver3(Joints_Target(3,1),Joints_Target(3,2),Joints_Target(3,3),FN3(1),FN3(2),FN3(3),2,'LineWidth',2)
 quiver3(Joints_Target(10,1),Joints_Target(10,2),Joints_Target(10,3),FN3(1),FN3(2),FN3(3),2,'LineWidth',2)

end

top=Joints_Target(4,:);

origin=Joints_Target(1,:)-Joints_Target(1,:);
z=top-Joints_Target(1,:);
y=FN3;

uno=Joints_Target(2,:)-Joints_Target(1,:);
due=Joints_Target(3,:)-Joints_Target(1,:);

ax_y=origin-y;
ax_z=origin-z;

ax_x=cross(ax_y,ax_z);
ax_x=ax_x/norm(ax_x);
ax_y=ax_y/norm(ax_y);
ax_z=ax_z/norm(ax_z);
uno=uno/norm(uno);
due=due/norm(due);

if(norm(ax_x-uno)>norm(ax_x-due))
   dx=Joints_Target(3,:);
else
   dx=Joints_Target(2,:);
   landmarks(3)=landmarks_actual(4);
   landmarks(4)=landmarks_actual(3);
end

%ARMS
syms d a b c x y z;
po1=Joints_Target(16,:);
po2=Joints_Target(13,:);
po3=Joints_Target(10,:);
vec1=po1-po2;
vec2=po3-po2;
vec_n=cross(vec1,vec2);
p2=vec_n(1)*(x-Joints_Target(16,1))+vec_n(2)*(y-Joints_Target(16,2))+vec_n(3)*(z-Joints_Target(16,3));
p2_explicit=solve(p2,z);


c = double(fliplr(coeffs(p2)));
po1=Joints_Target(10,:);
po2=Joints_Target(13,:);
po3=Joints_Target(10,:)+c([1,2,3]);
vec1=po1-po2;
vec2=po3-po2;
vec_n=cross(vec1,vec2);
p3=vec_n(1)*(x-Joints_Target(10,1))+vec_n(2)*(y-Joints_Target(10,2))+vec_n(3)*(z-Joints_Target(10,3));
p3_explicit=solve(p3,z);

c = double(fliplr(coeffs(p2)));

po1=Joints_Target(13,:);
po2=Joints_Target(10,:);
po3=Joints_Target(10,:)+c([1,2,3]);
P=[po1;po2;po3];
T=[1 2 3];
TR = triangulation(T,P);
FN = faceNormal(TR);

top=Joints_Target(13,:);

origin=Joints_Target(10,:)-Joints_Target(10,:);
z=top-Joints_Target(10,:);
y=FN3;
uno=Joints_Target(14,:)-Joints_Target(10,:);
due=Joints_Target(15,:)-Joints_Target(10,:);

ax_y=origin-y;
ax_z=origin-z;

ax_x=cross(ax_y,ax_z);
ax_x=ax_x/norm(ax_x);
ax_y=ax_y/norm(ax_y);
ax_z=ax_z/norm(ax_z);
uno=uno/norm(uno);
due=due/norm(due);

if(norm(ax_x-uno)>norm(ax_x-due))
   dx=Joints_Target(15,:);
else
   dx=Joints_Target(14,:);
   landmarks(2)=landmarks_actual(5);
   landmarks(5)=landmarks_actual(2);
end
end

