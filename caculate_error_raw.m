syms x1 y1 z1;    %建立符号变量
syms x2 y2 z2;
syms x3 y3 z3;
syms x4 y4 z4;
syms r1 r2 r3;
syms dr1 dr2 dr3
r1=((x4-x1)^2+(y4-y1)^2+(z4-z1)^2)^0.5;   %吉利公式
r2=((x4-x2)^2+(y4-y2)^2+(z4-z2)^2)^0.5;
r3=((x4-x3)^2+(y4-y3)^2+(z4-z3)^2)^0.5;

a= [x4 y4 z4];
f=[r1;r2;r3];

x=jacobian(f,a)%雅克比矩阵（代入坐标值前，输出结果为符号矩阵） 
X=vpa(subs(x,[x1,y1,z1;x2,y2,z2;x3,y3,z3;x4,y4,z4],[1,0,0;0,1,0;0,0,1;1,2,3]),3) %求雅克比矩阵(代入坐标值后，输出结果为数值矩阵)

dR=transpose([dr1 dr2 dr3]);  %转置
vpa(inv(X)*dR,3)  %写成线性方程，包含一次求逆使dr来表示

vpa(sum(inv(X).*inv(X),2),3)% 误差传递实际值

%inv(X)*dR
%X*A
%b=[1 2];
%y=subs(x,a,b)