close all;

H = -5 : .0001 : 5;

N = normpdf(H , rand(1) , rand(1))';



tic
Integration(.0001, N, "1D")
toc

byHand = toc;

tic
trapz(H , N)
toc

byMatlab = toc;


sign(byHand-byMatlab);

plot(N);