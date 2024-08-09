x1 = randn(100,1) + 3; y1 = randn(100,1) - 4;
x2 = randn(120,1) - 2; y2 = randn(120,1) + 1;
figure; plot(x1,y1,'o'); hold on; grid on; plot(x2,y2,'x'); 
xlabel('x'); ylabel('y');
SvmMdl = fitclinear([x1,y1; x2,y2], [zeros(size(x1)); ones(size(x2))]);
w = SvmMdl.Beta; b = SvmMdl.Bias; 
xb = linspace(min([x1;x2]), max([x1;x2]));
yb = (-b - w(1)*xb) / w(2);
plot(xb, yb, '.')