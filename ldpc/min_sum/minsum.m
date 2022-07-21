%function [estimated_code,syndrome,numblock,errorblock, sum_errorbit,rword,tword,fer,ber] = log_SPA()

% =================================================================
% minsum example.m
% By: Tinoosh Mohsenin 2004
% =================================================================

clear all;
hmatrix=[0 0 1 1 0 0 0 1 0 1 0 0;1 0 0 0 1 0 0 0 1 0 1 0;0 1 0 0 0 1 1 0 0 0 1 0;0 0 1 0 1 0 1 0 0 0 0 1;1 0 0 0 0 1 0 1 0 0 0 1;0 1 0 1 0 0 0 0 1 1 0 0];

n=12;                                     %code length
information_length=7;
m=6;                         %number of parity checks or number of rows
row_weight=4;
column_weight=2;

row_list = zeros(m,row_weight);
column_list=zeros(column_weight,n)
for i=1:m
    row_list(i,:)=find(hmatrix(i,:)~=0);    
end
for i=1:n
    column_list(:,i)=find(hmatrix(:,i)~=0);    
end

%AWGN channel is assumed

%RATE=(n-rank)/n;

codeword= [1     0     1     0     1     0     1     0     1     1     1     1]; 
tword=1-2*codeword;
%tword=1-2*zeros(1,n); % means that if x=0 => tword(BPSK)=1, if x=1=> BPSK=-1
numblock=0;
errorblock=0;
sum_errorbit=0;
max_errorblock=1;
iter_max = 5;				% maximum number of loop
iteration=0;
flag = 0;

while (errorblock<max_errorblock & flag==0)
lambda=[-9.1    4.9   -3.2    3.6   -1.4    3.1    0.3    1.6   -6.1   -2.5   -7.8   -6.8];
        
%------------------------------------
% First do the hard decesion and find if the received data is the actual
% code
%-----------------------------------
for i=1:n
      if lambda(i)>=0 
          estimated_code(i)=0;
       else
          estimated_code(i)=1;
      end    
end  
% syndrom check
syndrome=zeros(1,m);  
syndrome=multbin(estimated_code,hmatrix');
%sum_syndrome=sum(syndrome);

%Just for printing
beta_initial=zeros(m,n);
for i=1:m
    for j=1:row_weight
        beta_initial(i,row_list(i,j))=lambda(row_list(i,j));
    end
end
beta_initial
%if (sum_syndrome>0)    
 if (find (syndrome>0)>0)
% -----------------------------------
% misum decoding algorithm
% -----------------------------------
while (iteration <= iter_max)

%------------------------------------    
%row operation
%------------------------------------
  for i = 1:m
    for j = 1:row_weight
      sum = 0.0; prod = 1.0; min=100;
      for k = 1:row_weight
    	if k ~= j                   %process all elements in the row except the element j
	    tmp = lambda(row_list(i,k)) + beta(i,row_list(i,k));  %
        prod = prod* sign_tm(tmp);
        abs_tmp=abs(tmp);
        if (abs_tmp<min)
            min=abs_tmp;
        end  %if
        end % if
      end % for k
   alpha(i,row_list(i,j)) =prod * min;
   end% for j
  end % for i  

% column operation
  for i = 1:n
    sum2 = 0.0;
    for j = 1:column_weight
      sum1 = 0.0;
      for k = 1:column_weight
	   if k ~= j
	   sum1 = sum1 +  alpha(column_list(k,i),i);
       end;
      end;
      beta(column_list(j,i),i) = sum1;
      sum2 = sum2 + alpha(column_list(j,i),i);
    end;
    posterior(i) = lambda(i) + sum2;
  end;
  beta
  posterior
  %fprintf('----------------------------- loop = %d\n',loop_counter);
  %fprintf('log posterior probability ratio:\n');
  %disp(posterior);

  for i=1:n
      if posterior(i)>=0 
          estimated_code(i)=0;
      else
          estimated_code(i)=1;
      end
      
  end
% syndrom check
syndrome=multbin(estimated_code,hmatrix');
%sum_syndrome=sum(syndrome);

%  if (sum_syndrome==0)
%      iteration=iter_max+1;
%      numblock=numblock+1;
%  else
     if (iteration < iter_max)
        iteration=iteration+1;     
     elseif (iteration == iter_max)
     % increase number of Blocks, sum up the number of errors of the estimated code to previous value 
        for i=1:n
            sum_errorbit=sum_errorbit+estimated_code(i)~=tword(i);
        end
        errorblock=errorblock+1;
        iteration=iter_max+1;
        numblock=numblock+1;
     end
  %end
 % sum_syndrome      
  %iteration
end % end of the iteration
else
   % flag=1;
    numblock=numblock+1;
end % end if hard deceision 

%numblock
%errorblock
end % end of the max_errorblocks
numblock
errorblock
sum_errorbit
fer= max_errorblock/numblock;
ber=sum_errorbit/(numblock*n);

fid=fopen('minsum-report.txt','w');
fprintf(fid,'numblock= %d  errorblock= %d  sum_errorbit= %d fer= %f ber= %f',numblock,errorblock,sum_errorbit, fer, ber);
fclose(fid)
