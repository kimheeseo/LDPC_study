clc
clear all
load '9.mat'
%gallager matrix : 500 * 600
%VN degree=5, CN degree = 6
%v_degree=3, c_degree=6;
H=[1 0 1 0 1;0 0 0 1 1;0 1 0 1 0;0 1 0 0 1]
[K N]=size(H)
%K=9, N=12,M=N-K; M: H의 row의 갯수, N : H의 열의 갯수
%K=500, M=100, N=600

A=H;%불러온 H값을 A에 대입
A2=A;
if rank(A) ==K
    fprintf('full rank');
else
    fprintf('not full rank');
end
%B=rref(H) %rref만들기, B는 rref값
%B2=mod(B,2) %mod2버전 H값 B2 size : N*K

%% ref : row echelon form
[m,n] = size(A);

% Loop over the entire matrix.
i = 1;  
j = 1;

while (i <= m-1) && (j <= n)
   % Find value and index of largest element in the remainder of column j.
   k = find(A(i:m,j),m) + i - 1;
   k2=min(k);
   % Swap i-th and k-th rows.
   A([i k2],j:n) = A([k2 i],j:n);

   if isempty(k) == 1
       while isempty(k) == 1
        j=j+1;
        fprintf('error')
        k = find(A(i:m,j),m) + i - 1;
        k2=min(k);
        if length(k) ~=1
            for h=2:length(k)
            A(k(h),:)=A(k2,:)+A(k(h),:);
            A(k(h),:)=mod(A(k(h),:),2);
            end
        end
       end
   end
   i=j;
   k = find(A(i:m,j),m) + i - 1;
   k2=min(k);
   if length(k) ~=1
    for h=2:length(k)
        A(k(h),:)=A(k2,:)+A(k(h),:);
        A(k(h),:)=mod(A(k(h),:),2);
    end
   end

   i = i + 1;
   j = i;
end
A
A=rref(A)
A=mod(A,2)

%%codeword
%codeword=zeros(N,1);
codeword=ones(N,1);
%information=randi([0,1],1,N-rank(B));

%B2(rref H) : H' * codeword
for i = 1:rank(A)
    a(i)=min(find(A(i,:)==1)); %pivot 위치 파악 : pivot position - parity bits
end

c=(1:N);
d=setdiff(c,a); %non-pivot position : information bits

%information 만들기 : randi값 만들기
for i=1:length(d)
    codeword(d(i))=randi([0,1])
end

for i=1:rank(A)
    parity(i)=codeword(a(i))
end

codeword2=codeword
%parity용 B2
for i=1:K
    if mod(A(i,:)*codeword,2) ==1
        codeword2(a(i))=0
    end
end

qwe=A*codeword2 %qwe= H*codeword
qwe=mod(qwe,2)

qwe2=A2*codeword2 %qwe2= rref(H) * codeword
qwe2=mod(qwe2,2)

index=0
for i=1:K    
    if qwe(i,1)==0
        index=index+1
    end
end

if index==K
    fprintf('parity check *codeword = no error')
end

index=0
for i=1:K    
    if qwe2(i,1)==0
        index=index+1
    end
end

if index==K
    fprintf('rref parity check *codeword = no error\n')
end

%% BPSK
fprintf('bpsk전 codeword값') %BPSK 전
codeword2

%BPSK
fprintf('BPSK')
for i=1:N
    if codeword2(i)==1
        codeword2(i)=-1;
    else
        codeword2(i)=1;
    end
end
fprintf('bpsk후 codeword값')