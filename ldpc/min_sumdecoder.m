clc
clear all
%load '12.mat'
%gallager matrix : 500 * 600
%VN degree=5, CN degree = 6
%v_degree=3, c_degree=6;
%H=[1 0 0 0 1 1;0 1 0 1 1 0;1 1 1 0 0 0;0 1 0 1 0 1;1 1 0 0 1 1]
H=[0 0 1 1 0 0 0 1 0 1 0 0;1 0 0 0 1 0 0 0 1 0 1 0;0 1 0 0 0 1 1 0 0 0 1 0;0 0 1 0 1 0 1 0 0 0 0 1;1 0 0 0 0 1 0 1 0 0 0 1;0 1 0 1 0 0 0 0 1 1 0 0];
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
end %codeword -> codeword2 : A*codewrd2, A2*codeword2가 0이 되게 하기 위해서.

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

fprintf('BPSK')
for i=1:N
    if codeword2(i)==1
        codeword2(i)=-1;
    else
        codeword2(i)=1;
    end
end
fprintf('bpsk후 codeword값')

%% decoder
dB=[0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0];

for i=1:length(dB)
   no=1/(exp(dB(i)*log(10)/10));%noise
   tx=codeword2+sqrt(no/2)*randn(size(codeword2));%received signal
end

cn=0;
tx=[-9.1 4.9 -3.2 3.6 -1.4 3.1 0.3 1.6 -6.1 -2.5 -7.8 -6.8];
tx=transpose(tx);
%% Check node processing (Iteration 1)
% Variable node processing (iteration 1)은 0으로 초기화한 값이기에, 생략 
for j=1:K
    Q=find(H(j,:));
    for i=1:length(Q)
        cn(j,i)=tx(Q(i));%H의 j번째 row에서 1의 위치 파악
    end
end

A3=A2;
for j=1:K
    Q=find((H(j,:)));
    for i=1:length(find(A2(j,:)))
        A3(j,Q(i))=cn(j,i);
    end
end

A4=zeros(K,N)
for j=1:K
    for i=1:length(find(A2(j,:)))
        Q=find((A2(j,:)));
        A4(j,i)=A3(j,Q(i));
    end
end

%% min값 찾기
for i=1:K
    A5=1;
    for j=1:length(find(A4(i,:)))
        A5=1;
        Q=find(A4(i,:));
        qw=setdiff(Q,Q(j));
        min3=abs(A4(i,qw(1)));
        for k=1:length(qw)       
            if min3>abs(A4(i,qw(k)))
                min3=abs(A4(i,qw(k)));
            end
        end
        min4(i,j)=min3;
    end
end
%% sign 찾기
for i=1:K
    sign5=1;
    for j=1:length(find(A4(i,:)))
        sign5=1;
        Q=find(A4(i,:));
        qw=setdiff(Q,Q(j));
        for k=1:length(qw)
            sign5=sign5*sign(A4(i,qw(k)));
        end
        sign6(i,j)=sign5;
    end
end
sign6

QAX=min4.*sign6;
A7=A2;

for i=1:K
    scv=find(A2(i,:))
    for j=1:length(scv)
        A7(i,scv(j))=QAX(i,j)
    end
end

%% Code Estimation (Iteration 1)
tx2=tx;

for i=1:N
    if tx2(i)>0
        tx2(i)=0;
    else
        tx2(i)=1;
    end
end

syndrome=A2*tx2;
syndrome=mod(syndrome,2);

error=0;
for i=1:K
    error=syndrome(i)+error;
end

if error ~= 0
    fprintf('error')
end

%% iteration 2 : variable node update
iter=2;
iter2=1;
for iteration = 1:iter-1
aa=zeros(K,N);
for i=1:N
    bb=find(A2(:,i))
    cc=length(bb);
    for j=1:cc
        aa(i,j)=bb(j);
    end
end

[zz xx]=size(aa)

op=zeros(K,N);
for j=1:N 
    aq(j)=sum(A7(:,j))+tx(j);
end %aq는 variable node processing
%tx값을 op의 갯수만큼 빼기 수정해야함.

%% min값 찾기
for j=1:K
    Q=find(H(j,:));
    for i=1:length(Q)
        cn(j,i)=aq(Q(i));%H의 j번째 row에서 1의 위치 파악
    end
end

A3=A2;
for j=1:K
    Q=find((H(j,:)));
    for i=1:length(find(A2(j,:)))
        A3(j,Q(i))=cn(j,i);
    end
end%variable node processing값을 A2에 배치

A4=zeros(K,N)
for j=1:K
    for i=1:length(find(A2(j,:)))
        Q=find((A2(j,:)));
        A4(j,i)=A3(j,Q(i));
    end
end

%% min값 찾기
for i=1:K
    A5=1;
    for j=1:length(find(A4(i,:)))
        A5=1;
        Q=find(A4(i,:));
        qw=setdiff(Q,Q(j));
        min3=abs(A4(i,qw(1)));
        for k=1:length(qw)       
            if min3>abs(A4(i,qw(k)))
                min3=abs(A4(i,qw(k)));
            end
        end
        min4(i,j)=min3;
    end
end
%% sign 찾기
for i=1:K
    sign5=1;
    for j=1:length(find(A4(i,:)))
        sign5=1;
        Q=find(A4(i,:));
        qw=setdiff(Q,Q(j));
        for k=1:length(qw)
            sign5=sign5*sign(A4(i,qw(k)));
        end
        sign6(i,j)=sign5;
    end
end
sign6

QAX=min4.*sign6;
A8=A2;

for i=1:K
    scv=find(A2(i,:))
    for j=1:length(scv)
        A8(i,scv(j))=QAX(i,j)
    end
end

%% Code Estimation (Iteration 1)
tx3=aq;

for i=1:N
    if tx3(i)>0
        tx3(i)=0;
    else
        tx3(i)=1;
    end
end
tx3=transpose(tx3);
syndrome=A2*tx3;
syndrome=mod(syndrome,2);

error=0;
for i=1:K
    error=syndrome(i)+error;
end

if error ~= 0
    fprintf('error')
else
    fprintf('no error \n')
    break;
end
iter2=iter2+1;
end
fprintf('decoding end')
fprintf('total_Number of iterations')
iter2=iter2+1;
iter2