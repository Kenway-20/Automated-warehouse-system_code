clc;clear;close all
range=4;  %������������ⷽ������Ӹ���
goods=24;  %��������

NIND=40;   %��Ⱥ��Ŀ
MAXGEN=200;  %�Ŵ�����
px=0.7;  %�������
pm=0.01;  %�������
NVAR=3*goods;  %Ⱦɫ�峤��
PRECI=15;   %����������λ��
GGAP=0.9;   %����
trace=zeros(MAXGEN,2);  %
%tracepbj=zeros(40,MAXGEN);
%BaseV=zeros(1,NVAR)+5;
FieldD=[rep([PRECI],[1,NVAR]);rep([1;range],[1,NVAR]);rep([1;0;1;1],[1,NVAR])];  %
Chrom=crtbp(NIND,NVAR*PRECI);  %������Ⱥ����ΪNIND�����峤��ΪMVAR*PRECI�Ķ�������Ⱥ
%Chrom=crtbp(NIND,NVAR,BaseV);
gen=0;
w1=0.2;w2=0.2;w3=0.6;

%ObjV=1/(w1*fun01(bs2rv(Chrom,FieldD))+w2*fun02(bs2rv(Chrom,FieldD))+w3*fun03(bs2rv(Chrom,FieldD))+1);
matrix01=bs2rv(Chrom,FieldD);
%matrix01=Chrom;
ObjV=w1*fun01(matrix01)+w2*fun02(matrix01)+w3*fun03(matrix01);
FitnV=zeros(NIND,1);
while gen<MAXGEN,
    FitnV=ranking(ObjV);  %������Ӧ��
    %FitnV=1/(ObjV+1);
    SelCh=select('sus',Chrom,FitnV,GGAP);  %ѡ��߼�����
    SelCh=recombin('xovsp',SelCh,px);    %��������
    SelCh=mut(SelCh,pm);   %����
    %ObjVSel=1rv/(w1*fun01(bs2rv(SelCh,FieldD))+w2*fun02(bs2rv(SelCh,FieldD))+w3*fun03(bs2rv(SelCh,FieldD))+1);
    matrix01=bs2rv(SelCh,FieldD);
    %matrix01=SelCh;
    ObjVSel=w1*fun01(matrix01)+w2*fun02(matrix01)+w3*fun03(matrix01);
    [Chrom, ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);  %�ز��룬�������Ӵ��븸���滻
    matrix01=bs2rv(Chrom,FieldD);
    gen=gen+1;
    trace(gen,1)=min(ObjV);
    trace(gen,2)=sum(ObjV)/length(ObjV);
end
%%ȡ��
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
%%�ҳ��ظ�����
%�������ʽ
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
           num=num+1;  %�ظ��ĸ�������ȥ����
           record02(i,1:3)=[rand rand rand];
       end
    end
    if num>0
        hang=hang+1;  %���ظ����¼
        repeat(hang,1:3)=xyz;
        repeat(hang,4)=num;
    end
end
%%�����ظ�����
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
                        if buer==1  %1��ʾ���                       
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
legend('��ı仯','��Ⱥ��ֵ�ı仯');
%matrix02=floor(matrix01);
for i=1:size(X,1)
    [x,y,z]=meshgrid(X(i)-1:X(i),Y(i)-1:Y(i),Z(i)-1:Z(i));
    fun=(x+y+z)/3;
    slice(x,y,z,fun,X(i)-1:X(i),Y(i)-1:Y(i),Z(i)-1:Z(i));
    hold on
end
