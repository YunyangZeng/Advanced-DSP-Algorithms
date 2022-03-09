function [a,e,k] = levinson_durbin_recursion(r,n)
% Performs Levinson-Derbin recursion
P=n;
e=[r(1)];
a=zeros(P);
k=zeros(P,1);
for i=1:P
    if i<2
        k(i)=(r(i+1)-0)/e(i);
        a(i,i)=k(i);
        e(i+1)=(1-k(i)^2)*e(i);
    else
        aj_phii_j=0;
        for j=1:i-1
            
            aj_phii_j=aj_phii_j+a(i-1,j)*r(i-j+1);
        end
        k(i)=(r(i+1)-aj_phii_j)/e(i);
        a(i,i)=k(i);
        for j=1:i-1
            a(i,j)=a(i-1,j)-k(i)*a(i-1,i-j);
        end
        e(i+1)=(1-k(i)^2)*e(i);
    end
end
a=a(P,1:P);
end

