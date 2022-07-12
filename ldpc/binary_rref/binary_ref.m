clear all; clc;
%A2=[1 0 1 1 1;0 0 0 1 1;0 1 0 1 1;0 1 0 0 1];
%A2=[1 0 1 0 1;0 0 0 1 1;0 1 0 1 0;0 1 0 0 1];
%A2=[1 0 1 1 1;0 1 0 1 1;0 1 0 0 1;0 0 0 1 0];
%A2=[0 1 0 1 1; 0 0 1 1 1;1 0 1 1 1;1 0 0 0 1];
%A2=[1 0 0 0 1 1;0 1 0 1 1 0;1 1 1 0 0 0;0 1 0 1 0 1;1 1 0 0 1 1]
%A2=[1 0 1 1 0 0;0 1 0 1 1 0;0 1 1 0 0 1;0 1 0 1 1 1;0 1 1 0 1 1]
A2=[1 1 1 0 0 0;0 0 0  1 1 1;0 0 0 1 1 1;1 1 1 0 0 0]
A=A2;
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
