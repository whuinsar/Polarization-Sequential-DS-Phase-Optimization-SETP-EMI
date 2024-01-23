function gamma = Regularization(gamma,stacksize,nrep)
[~,d] = eig(double(gamma));penalty = 1e-6;reg = eye(size(gamma));
count=1;
while min(diag(d)) < 0
        temp = diag(gamma);
        gamma = penalty*gamma;
        gamma(logical(eye(size(gamma))))= temp;
        [~,d] = eig(double(gamma));
    count=count+1;
    if count>stacksize
       return
    end
end
end