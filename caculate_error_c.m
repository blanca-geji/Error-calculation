function [ error_xyz  ] = caculate_error_c( input_args )%#codegen
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

X=ones(3,3);
X(1,1)=-(2*input_args(1,1) - 2*input_args(4,1))/(2*((input_args(1,1) - input_args(4,1))^2 + (input_args(1,2) - input_args(4,2))^2 + (input_args(1,3) - input_args(4,3))^2)^(1/2));
X(1,2)=-(2*input_args(1,2) - 2*input_args(4,2))/(2*((input_args(1,1) - input_args(4,1))^2 + (input_args(1,2) - input_args(4,2))^2 + (input_args(1,3) - input_args(4,3))^2)^(1/2));
X(1,3)=-(2*input_args(1,3) - 2*input_args(4,3))/(2*((input_args(1,1) - input_args(4,1))^2 + (input_args(1,2) - input_args(4,2))^2 + (input_args(1,3) - input_args(4,3))^2)^(1/2));
X(2,1)=-(2*input_args(2,1) - 2*input_args(4,1))/(2*((input_args(2,1) - input_args(4,1))^2 + (input_args(2,2) - input_args(4,2))^2 + (input_args(2,3) - input_args(4,3))^2)^(1/2));
X(2,2)=-(2*input_args(2,2) - 2*input_args(4,2))/(2*((input_args(2,1) - input_args(4,1))^2 + (input_args(2,2) - input_args(4,2))^2 + (input_args(2,3) - input_args(4,3))^2)^(1/2));
X(2,3)=-(2*input_args(2,3) - 2*input_args(4,3))/(2*((input_args(2,1) - input_args(4,1))^2 + (input_args(2,2) - input_args(4,2))^2 + (input_args(2,3) - input_args(4,3))^2)^(1/2));
X(3,1)=-(2*input_args(3,1) - 2*input_args(4,1))/(2*((input_args(3,1) - input_args(4,1))^2 + (input_args(3,2) - input_args(4,2))^2 + (input_args(3,3) - input_args(4,3))^2)^(1/2));
X(3,2)=-(2*input_args(3,2) - 2*input_args(4,2))/(2*((input_args(3,1) - input_args(4,1))^2 + (input_args(3,2) - input_args(4,2))^2 + (input_args(3,3) - input_args(4,3))^2)^(1/2));
X(3,3)=-(2*input_args(3,3) - 2*input_args(4,3))/(2*((input_args(3,1) - input_args(4,1))^2 + (input_args(3,2) - input_args(4,2))^2 + (input_args(3,3) - input_args(4,3))^2)^(1/2));

%x =[ -(2*x1 - 2*x4)/(2*((x1 - x4)^2 + (y1 - y4)^2 + (z1 - z4)^2)^(1/2)), -(2*y1 - 2*y4)/(2*((x1 - x4)^2 + (y1 - y4)^2 + (z1 - z4)^2)^(1/2)), -(2*z1 - 2*z4)/(2*((x1 - x4)^2 + (y1 - y4)^2 + (z1 - z4)^2)^(1/2));
%     -(2*x2 - 2*x4)/(2*((x2 - x4)^2 + (y2 - y4)^2 + (z2 - z4)^2)^(1/2)), -(2*y2 - 2*y4)/(2*((x2 - x4)^2 + (y2 - y4)^2 + (z2 - z4)^2)^(1/2)), -(2*z2 - 2*z4)/(2*((x2 - x4)^2 + (y2 - y4)^2 + (z2 - z4)^2)^(1/2));
%     -(2*x3 - 2*x4)/(2*((x3 - x4)^2 + (y3 - y4)^2 + (z3 - z4)^2)^(1/2)), -(2*y3 - 2*y4)/(2*((x3 - x4)^2 + (y3 - y4)^2 + (z3 - z4)^2)^(1/2)), -(2*z3 - 2*z4)/(2*((x3 - x4)^2 + (y3 - y4)^2 + (z3 - z4)^2)^(1/2))];

%X=subs(x,[x1,y1,z1;x2,y2,z2;x3,y3,z3;x4,y4,z4],input_args);

error_xyz=sum(inv(X).*inv(X),2);

end
