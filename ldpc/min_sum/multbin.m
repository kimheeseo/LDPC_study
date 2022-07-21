%Multiplies two matrices in binary
%Nov.2005:written by Tinoosh Mohsenin. 
function[result]=multbin(a,b)

a=[1 0 1
   1 1 0];
b=[1 1  
   1 1
   1 0]
[row_size m]=size(a);
[m col_size]=size(b);
result=zeros(row_size,col_size);
for i=1:row_size
    for j=1:col_size
        index=find(b(:,j)>0);
        for (x=1:length(index))
         result(i,j)=xor(a(i,index(x)),result(i,j));
        end
    end
end