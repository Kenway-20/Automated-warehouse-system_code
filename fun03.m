function f3=fun03(matrix01)
row01=size(matrix01,1);
col01=size(matrix01,2);
f3=zeros(row01,1);
load auto;
zhonglei=auto;
zhonglei=cell2mat(zhonglei);
for i=1:row01
    adress=matrix01(i,:);
    f3(i,1)=0;
    for j=1:col01/3
      if zhonglei(j,1)=='A'
        na=3;nb=3;nc=1;
      elseif zhonglei(j,1)=='B'
        na=2;nb=1;nc=2;
      else
        na=4;nb=4;nc=1;
      end
    x=adress(1,3*j-2);y=adress(1,3*j-1);z=adress(1,3*j);
    f3(i,1)=sqrt((x-na)^2+(y-nb)^2+(z-nc)^2)+f3(i,1);
    end
end
