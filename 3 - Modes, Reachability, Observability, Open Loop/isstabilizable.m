function result = isstabilizable(A,B)

[V,J]= eig(A);
eigenvectors= V;
eigenvalues= diag(J);

% all modes stable condition
if all(real(eigenvalues)<0)
    result=true;
    return
end

% unstable modes stabilisability check
idx= find(real(eigenvalues)>=0);

R= ctrb(A,B); 
tol= max(size(R))*eps(norm(R)); % Default tolerance for Rank()
for i=1:length(idx)
    if norm(R'* eigenvectors(:,idx(i))) < tol
        result=false;
        return
    end
end

% Otherwise unstable modes are stabilizable
result=true;

end