clc;
clear all;
N=2400, M=2000, DoYouWantACE=0;
H=ProgressiveEdgeGrowthACE(N, M, DoYouWantACE)

[k n]=size(H);

for i=1:k
    j(i)=sum(H(i,:));
end

j=transpose(j);

a=find(j==5);%5
b=find(j==7);%7

for i=1:length(a)
    aa(i,:)=find(H(a(i),:)==1);
    bb(i,:)=find(H(b(i),:)==1);
end

[sa1 sa2]=size(aa);
[sb1 sb2]=size(bb);

if sb2>sa2
    sa=sb2-sa2;
    sa22=zeros(length(a),sa);
    aa2=[aa sa22];
else
    sa=sa2-sb2;
    sa22=zeros(length(a),sa);
    aa2=[bb sa22];
end

for i=1:length(a)
    aaa(i,:)=H(a(i),:);
    bbb(i,:)=H(b(i),:);
end
fprintf('aaa \n')
aaa2=aaa;
bbb2=bbb;

for i=1:length(a)
    for j=1:N        
       
            if sb2>sa2
                if aaa(i,j)==0 && bbb(i,j)==1
                fprintf('change')
                i
                j
                aaa2(i,j)=1;
                bbb2(i,j)=0;
                break;
                end
            else
                if aaa(i,j)==1 && bbb(i,j)==0
                aaa2(i,j)=0;
                bbb2(i,j)=1;
                break;
                end
            end
    end
end
H2=H;

for i=1:length(a)
    H2(a(i),:)=aaa2(i,:);
    H2(b(i),:)=bbb2(i,:);
end

filename='PEG_2400.mat';
save(filename)
