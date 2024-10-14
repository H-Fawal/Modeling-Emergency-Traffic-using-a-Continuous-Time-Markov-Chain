clc;
clear all;
C = 3;
Picj = [];
CR_M2M = 0;
CR_H2H=0;
CR_SOS=0;
L1=1;
L2=1;
L3=2;
%L3=1;
M1=1;
M2=1;
M3=2;
%M3=1;
disp([' --- Continuous Time Markov Chain model- Single Traffic']);
disp([' --- C = ', num2str(C)]);
syms P0c0 P1c0 P2c0 P3c0 P0c1 P0c2 P0c3 P1c1 P1c2 P2c1
% p1=P0c0=(0,0)
% p2=P1c0=(1,0)
% p3=P2c0=(2,0)
% p4=P3c0=(3,0)
% p5=P0c1=(0,1)
% p6=P0c2=(0,2)
% p7=P0c3=(0,3)
% p8=P1c1=(1,1)
% p9=P1c2=(1,2)
% p10=P2c1=(2,1)
%EQUATIONS
eqn0=P0c0 +P1c0+ P2c0 +P3c0 +P0c1 +P0c2 +P0c3+ P1c1 +P1c2 +P2c1 ==1
%EMTY
eqn1 = ((2*L1 + 2*L2 + (L1 & L2) )*P0c0 )- ((M1 * P1c0) + ( M2 * P0c1) + ( M1 * P2c0) + (M2 * P0c2) + (((M3)) * P1c1)) ==0
%OCCUPIED
eqn2 = ((L1 * P0c0) + (M1 * P2c0) + (M2 * P1c1) + (M1 * P3c0) + (M2 * P1c2) + ((M3) * P2c1)) - ((2*L1 + 2*L2 + (L3) + M1 ) * P1c0)==0
eqn3 = ((L2*(P0c0))+ (M1* P1c1)+(M2* P0c2)+ (M1* P2c1)+(M2* P0c3)+((M3)* P1c2))-((( 2*L1+ 2*L2 + ((M3))+ M2)*P0c1)) ==0
eqn4= ((L1*P1c0)+(L1*P0c0)+(M2*P2c1)+(M1*P3c0))- ((2*M1 +L1+L2)*P2c0)==0
eqn5= ((L2*P0c0)+(L2*P0c1)+(M2*P0c3)+(M1*P1c2))-(((2*M2) + L1+ L2)*P0c2) ==0
eqn6= ((L2*P1c0)+(L1*P0c1)+ ((L1 & L2)*P0c0)+(M2*P1c2)+(M1*P2c1))- (((M3)+M1+M2+L1+L2)*P1c1) ==0
%FULL
eqn7 =((L1*P0c2)+(L2*P1c1)+(L2*P1c0)+((L3)*P0c1))-((M1+2*M2+((M3)))*P1c2) ==0
eqn8 =(L1*P1c1)+(L2*P2c0)+(L1*P0c1)+((L3)*P1c0)-(((2*M1+M2 +(M3) )*P2c1))==0
eqn9 = ((L2 * P0c2) + (L2* P0c1)) -((2*M2)*P0c3)==0
% eqn10 = ((L1* P2c0) + (L1* P1c0)) - ((2*M1)*P3c0)==0
%[i,j] = equationsToMatrix(eqn9,eqn1,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn2,eqn10)
% [M,B] = equationsToMatrix(eqn0,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10)
[M,B] = equationsToMatrix(eqn0,eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9)%
%M = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
% 9, -1, -1, 0, -1, -1, 0, -2, 0, 0;
% 2, -11, 1, 1, 0, 0, 0, 1, 1, 2;
% 2, 0, 0, 0, -11, 1, 1, 1, 2, 1;
% 2, 2, -6, 1, 0, 0, 0, 0, 0, 1;
% 2, 0, 0, 0, 2, -6, 1, 0, 1, 0;
% 1, 2, 0, 0, 2, 0, 0, -8, 1, 1;
% 0, 2, 0, 0, 2, 2, 0, 2, -5, 0;
% 0, 2, 2, 0, 2, 0, 0, 2, 0, -5;
% 0, 0, 0, 0, 2, 2, -2, 0, 0, 0;]
% B = [1;0;0;0;0;0;0;0;0;0];
disp([M,B]);
Rank=rank([M,B]);
disp(['Rank:'])
disp(Rank);
Determinant=abs(det([M]));
disp(['determinant:'])
disp(Determinant);
% Determinant=det([i,j]);
% disp(['determinant:'])
% disp(Determinant);
% answer = linsolve(i,j);
answer = linsolve(M,B);
figure, plot(answer);
%hist(answer)
%answer=inv(A)*b;
disp(['Displaying the result:'])
disp(answer);
disp(eval(answer));
disp(['summation:'])
summation = sum(answer)
vpa (summation)
% normalizedAns = answer/summation
% vpa(normalizedAns)
% checkvalues = sum(normalizedAns)
% subplot(212)
% plot(normalizedAns)
xlabel('state')
ylabel('probability to be at a certain state')
%
CR_M2M = [];
for i=1:C+1
for j=1:C+1
if((i+j)-2>C)
continue;
end
CR_M2M = [CR_M2M (i-1).*M1.*answer(i)];
%strM2M = strcat('CR_M2M = CR_M2M + ', num2str(i-1), '*M1*P',num2str(i-1),'c',num2str(j-1),';')
end
end
CR_M2M = sum(CR_M2M);
disp('the SCR OF M2M')
disp(CR_M2M);
vpa (CR_M2M)
CR_H2H = [];
for i=1:C+1
for j=1:C+1
if((i+j)-2>C)
continue;
end
CR_H2H = [CR_H2H (j-1).*M1.*answer(i)];
%strM2M = strcat('CR_M2M = CR_M2M + ', num2str(i-1), '*M1*P',num2str(i-1),'c',num2str(j-1),';')
end
end
CR_H2H = sum(CR_H2H);
disp('the SCR OF H2H')
disp(CR_H2H);
vpa (CR_H2H)
% CR_SOS = [];
%
% for i=1:C+1
% for j=1:C+1
% if((i+j)-2>C)
% continue;
% end
% CR_SOS = [CR_SOS ((i-1)&(j-1)).*M3.*(answer(i)& answer(j))];
% %strM2M = strcat('CR_M2M = CR_M2M + ', num2str(i-1), '*M1*P',num2str(i-1),'c',num2str(j-1),';')
% end
% end
%
% CR_SOS = sum(CR_SOS);
% disp('the SCR OF SOS')
% disp(CR_SOS);
% vpa (CR_SOS)
% CR_H2H = []
% for i=1:C+1
% for j=1:C+1
% if((i+j)-2>C)
% continue;
% end
% CR_H2H = [CR_H2H (j-1)*M2*Picj]
%
%
%
% end
% end
%
% CR_H2H = sum(CR_H2H)
%
% disp('the SCR OF H2H')
% disp(CR_H2H);
%
% for i=1:C+1
% for j=1:C+1
% if((i+j)-2>C)
% continue;
% end
% CR_M2M = sum((i-1)*M1*Picj)
% CR_M2M = CR_M2M + (i-1).*M1.*answer(i);
% %strM2M = strcat('CR_M2M = CR_M2M + ', num2str(i-1), '*M1*P',num2str(i-1),'c',num2str(j-1),';')
% disp('the SCR OF M2M')
% disp(CR_M2M);
% M2M_SCR = sum(CR_M2M)
% X=vpa (M2M_SCR)
% disp(X)
%
%
%
% end
% end
%
%
% for i=1:C+1
% for j=1:C+1
% if((i+j)-2>C)
% continue;
% end
% CR_H2H = sum((j-1)*M2*Picj)
% CR_H2H = CR_H2H + (j-1).*M2.*answer(j);
% %strM2M = strcat('CR_M2M = CR_M2M + ', num2str(i-1), '*M1*P',num2str(i-1),'c',num2str(j-1),';')
% disp('the SCR OF H2H')
% disp(CR_H2H);
% H2H_SCR = sum(CR_H2H)
% Y=vpa (H2H_SCR)
% disp(Y)
%
%
%
% end
% end
% %
% for i=1:C+1
% for j=1:C+1
% if((i+j)-2>C)
% continue;
% end
% CR_H2H = sum((i*j)*M3*Picj)
% CR_H2H = CR_H2H + ((i-1).*(j-1).*M3.*answer(i).*answer(j));
% %strM2M = strcat('CR_M2M = CR_M2M + ', num2str(i-1), '*M1*P',num2str(i-1),'c',num2str(j-1),';')
% disp('the SCR OF SOS')
% disp(CR_SOS);
% SOS_SCR = sum(CR_SOS)
% vpa (SOS_SCR)
%
%
% end
% end
