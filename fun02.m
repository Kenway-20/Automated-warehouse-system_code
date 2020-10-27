function f2=fun02(matrix01)
L0=1;
row01=size(matrix01,1);
col01=size(matrix01,2);
m=xlsread('auto',1,'D2:D25');
for i=1:row01
    fenzi(i,1)=0;
    for j=1:col01/3   
    fenzi(i,1)=m(j,1)*matrix01(i,3*j)*L0+fenzi(i,1);
    end
end
f2=fenzi;

