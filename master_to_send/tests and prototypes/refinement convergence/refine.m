function [v_ref] = refine(v)
%returns design vector for doubly refined grid
N= size(v,1);
d=floor(sqrt(N));
v_ref = zeros(4*N,1)
for l=1:N
i = 1+ mod(l-1,d);
j= ceil(l/d);    

m = i+(d)*(j-1);

v_ref((2*i-1)+(2*d)*(2*(j)-2))=v(l);
v_ref((2*i)+(2*d)*(2*(j)-2))=v(l);
v_ref((2*i-1)+(2*d)*(2*(j)-1))=v(l);
v_ref((2*i)+(2*d)*(2*(j)-1))=v(l);
end
end

