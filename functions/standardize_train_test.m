function [XtrZ,XteZ,YtrZ,YteRaw,muY,sdY] = ...
    standardize_train_test(Xtr,Xte,Ytr,Yte)
muX = mean(Xtr,1);
sdX = std(Xtr,0,1);
sdX(sdX < eps('single')) = 1;
XtrZ = (Xtr-muX)./sdX;
XteZ = (Xte-muX)./sdX;

muY = mean(Ytr,1);
sdY = std(Ytr,0,1);
sdY(sdY < eps('single')) = 1;
YtrZ = (Ytr-muY)./sdY;
YteRaw = Yte;
end
