function [sign] = sign_tm(num)
if num>0
    sign=1;
elseif num<0
    sign=-1;
else
    sign=0;
end
    