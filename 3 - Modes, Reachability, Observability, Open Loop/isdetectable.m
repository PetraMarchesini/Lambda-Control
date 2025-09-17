function result = isdetectable(A,C)

[V,J]= eig(A);
eigenvectors= V;
eigenvalues= diag(J);

% all modes stable condition
if all(real(eigenvalues)<0)
    result=true;
    return
end

% unstable modes observability check
idx= find(real(eigenvalues)>=0);

O= obsv(A,C);
tol= max(size(O))*eps(norm(O)); % Default tolerance for Rank()
for i=1:length(idx)
    if norm(O* eigenvectors(:,idx(i))) < tol
        result=false;
        return
    end
end

% Otherwise unstable modes are detectable
result=true;

end