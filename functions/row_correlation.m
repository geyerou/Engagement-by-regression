function r = row_correlation(A,B)
A=double(A); B=double(B);
A=A-mean(A,2); B=B-mean(B,2);
num=sum(A.*B,2);
den=sqrt(sum(A.^2,2).*sum(B.^2,2));
r=single(num./max(den,eps));
end
