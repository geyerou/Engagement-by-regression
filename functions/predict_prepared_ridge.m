function prediction = predict_prepared_ridge(model,Ytrain,block_size)
n_vox=size(Ytrain,2);
prediction=zeros(size(model.Xtest,1),n_vox,'single');
for first=1:block_size:n_vox
    idx=first:min(first+block_size-1,n_vox);
    yz=double((Ytrain(:,idx)-model.muY(idx))./model.sdY(idx));
    B=model.V*(model.shrink.*(model.U'*yz));
    prediction(:,idx)=single((model.Xtest*B).* ...
        double(model.sdY(idx))+double(model.muY(idx)));
end
end
