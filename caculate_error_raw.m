syms x1 y1 z1;    %�������ű���
syms x2 y2 z2;
syms x3 y3 z3;
syms x4 y4 z4;
syms r1 r2 r3;
syms dr1 dr2 dr3
r1=((x4-x1)^2+(y4-y1)^2+(z4-z1)^2)^0.5;   %������ʽ
r2=((x4-x2)^2+(y4-y2)^2+(z4-z2)^2)^0.5;
r3=((x4-x3)^2+(y4-y3)^2+(z4-z3)^2)^0.5;

a= [x4 y4 z4];
f=[r1;r2;r3];

x=jacobian(f,a)%�ſ˱Ⱦ��󣨴�������ֵǰ��������Ϊ���ž��� 
X=vpa(subs(x,[x1,y1,z1;x2,y2,z2;x3,y3,z3;x4,y4,z4],[1,0,0;0,1,0;0,0,1;1,2,3]),3) %���ſ˱Ⱦ���(��������ֵ��������Ϊ��ֵ����)

dR=transpose([dr1 dr2 dr3]);  %ת��
vpa(inv(X)*dR,3)  %д�����Է��̣�����һ������ʹdr����ʾ

vpa(sum(inv(X).*inv(X),2),3)% ����ʵ��ֵ

%inv(X)*dR
%X*A
%b=[1 2];
%y=subs(x,a,b)