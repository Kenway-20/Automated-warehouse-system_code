clc;clear;close all
range=4;  %立方体货架任意方向的箱子格数
goods=24;  %货物数量

NIND=40;   %种群数目
MAXGEN=200;  %遗传代数
px=0.7;  %交叉概率
pm=0.01;  %变异概率
NVAR=3*goods;  %染色体长度
PRECI=15;   %变量二进制位数
GGAP=0.9;   %代沟
trace=zeros(MAXGEN,2);  %
%tracepbj=zeros(40,MAXGEN);
%BaseV=zeros(1,NVAR)+5;
FieldD=[rep([PRECI],[1,NVAR]);rep([1;range],[1,NVAR]);rep([1;0;1;1],[1,NVAR])];  %
Chrom=crtbp(NIND,NVAR*PRECI);  %生成种群长度为NIND，个体长度为MVAR*PRECI的二进制种群
%Chrom=crtbp(NIND,NVAR,BaseV);
gen=0;
w1=0.2;w2=0.2;w3=0.6;

%ObjV=1/(w1*fun01(bs2rv(Chrom,FieldD))+w2*fun02(bs2rv(Chrom,FieldD))+w3*fun03(bs2rv(Chrom,FieldD))+1);
matrix01=bs2rv(Chrom,FieldD);
%matrix01=Chrom;
ObjV=w1*fun01(matrix01)+w2*fun02(matrix01)+w3*fun03(matrix01);
FitnV=zeros(NIND,1);
while gen<MAXGEN,
    FitnV=ranking(ObjV);  %计算适应度
    %FitnV=1/(ObjV+1);
    SelCh=select('sus',Chrom,FitnV,GGAP);  %选择高级个体
    SelCh=recombin('xovsp',SelCh,px);    %交叉重组
    SelCh=mut(SelCh,pm);   %变异
    %ObjVSel=1rv/(w1*fun01(bs2rv(SelCh,FieldD))+w2*fun02(bs2rv(SelCh,FieldD))+w3*fun03(bs2rv(SelCh,FieldD))+1);
    matrix01=bs2rv(SelCh,FieldD);
    %matrix01=SelCh;
    ObjVSel=w1*fun01(matrix01)+w2*fun02(matrix01)+w3*fun03(matrix01);
    [Chrom, ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);  %重插入，将部分子代与父代替换
    matrix01=bs2rv(Chrom,FieldD);
    gen=gen+1;
    trace(gen,1)=min(ObjV);
    trace(gen,2)=sum(ObjV)/length(ObjV);
end
%%取整
    for j=1:goods
       xyz=matrix01(1,(3*j-2):(3*j));
       x=xyz(1);y=xyz(2);z=xyz(3);
       value0=10000000;
       for ix=floor(x):floor(x)+1
           for iy=floor(y):floor(y)+1
               for iz=floor(z):floor(z)+1
                   ixyz=[ix iy iz];
                   value=w1*fun01(ixyz)+w2*fun02(ixyz)+w3*fun03(ixyz);
                   if value<value0
                       value0=value;
                       record(1,(3*j-2):(3*j))=[ix iy iz];
                   end 
               end 
           end
       end
    end
%%找出重复数据
%变成列形式
for j=1:goods
    xyz=record(1,(3*j-2):(3*j));
    record01(j,:)=xyz;
end
record02=record01;
hang=0;
for j=1:goods-1
    xyz=record02(j,:);
    num=0;
    for i=j+1:goods
       comxyz=record02(i,:);
       if comxyz==xyz
           num=num+1;  %重复的个数（除去本身）
           record02(i,1:3)=[rand rand rand];
       end
    end
    if num>0
        hang=hang+1;  %有重复则记录
        repeat(hang,1:3)=xyz;
        repeat(hang,4)=num;
    end
end
%%处理重复数据
record03=record01;
m=size(repeat,1);
final01=zeros(1,3);
for i=1:m
    rexyz=repeat(i,:);
    rex=rexyz(1);rey=rexyz(2);rez=rexyz(3);
    time=0;
    for x=rex-1:rex+1
        for y=rey-1:rey+1
            for z=rez-1:rez+1
                if x>=1 && y>=1 && z>=1 && x<=range && y<=range && z<=range
                    xyz1=[x y z];aa=0;
                    for j=1:size(record03,1)
                        buer=isequal(record03(j,:),xyz1);
                        if buer==1  %1表示相等                       
                        aa=aa+1; 
                        end   
                    end
                    if aa==0
                       time=time+1;
                       store(time,:)=[x,y,z];
                    end                   
                end
            end
        end
        
    end
    target=w1*fun01(store)+w2*fun02(store)+w3*fun03(store);
        target(:,2:4)=store;
        final=sortrows(target,1);
        final=final(1:repeat(i,4),:);
        record03=[record03;final(:,2:4)];
       
        final01=[final01;final(:,2:4)];
        final01(1,:)=[];
end
bbk=0;
for i=1:goods
    xyz01=record02(i,:);
    if xyz01(1)>=1
        bbk=bbk+1;
        rbr(bbk,1:3)=xyz01;
    end
end
finalboss=[final01;rbr];
X=finalboss(:,1);
Y=finalboss(:,2);
Z=finalboss(:,3);

figure
plot3(X,Y,Z,'o');
axis([1 range 1 range 1 range]);
grid on 

figure
plot(trace(:,1));hold on;
plot(trace(:,2),'-.');grid;
legend('解的变化','种群均值的变化');
%matrix02=floor(matrix01);
for i=1:size(X,1)
    [x,y,z]=meshgrid(X(i)-1:X(i),Y(i)-1:Y(i),Z(i)-1:Z(i));
    fun=(x+y+z)/3;
    slice(x,y,z,fun,X(i)-1:X(i),Y(i)-1:Y(i),Z(i)-1:Z(i));
    hold on
end
