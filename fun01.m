function f1=fun01(matrix01)
row01=size(matrix01,1);
col01=size(matrix01,2);
load zhouzhuanlv;
p=zhouzhuanlv;
vx=1;vy=1;vz=1;
f1=zeros(row01,1);
for i=1:row01
    adress=matrix01(i,:);
    f1(i,1)=0;
    for j=1:col01/3
    x=adress(1,3*j-2);y=adress(1,3*j-1);z=adress(1,3*j);
    %f1=f1+p(i)*(x/vx+y/vy+z/vz);
    f1(i,1)=p(j)*(x/vx+y/vy+z/vz)+f1(i,1);
    end
end
    




