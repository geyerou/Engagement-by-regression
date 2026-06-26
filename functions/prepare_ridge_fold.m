function model = prepare_ridge_fold(Xtrain,Xtest,Ytrain,lambda)
model.muX=mean(Xtrain,1);
model.sdX=std(Xtrain,0,1); model.sdX(model.sdX<eps('single'))=1;
Xtr=double((Xtrain-model.muX)./model.sdX);
model.Xtest=double((Xtest-model.muX)./model.sdX);
[model.U,S,model.V]=svd(Xtr,'econ');
s=diag(S);
model.shrink=s./(s.^2+size(Xtr,1)*lambda);
model.muY=mean(Ytrain,1);
model.sdY=std(Ytrain,0,1); model.sdY(model.sdY<eps('single'))=1;
end
