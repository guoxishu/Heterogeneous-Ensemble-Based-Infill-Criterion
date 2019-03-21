function y=RBF2test(x, p, N)

Centers=p.Centers;Spreads=p.Spreads;
W2=p.W2;B2=p.B2;
TestDistance = dist(Centers,x');
TestSpreadsMat = repmat(Spreads,1,N);
TestHiddenUnitOut = radbas(TestDistance./TestSpreadsMat);
y= (W2*TestHiddenUnitOut+repmat(B2,1,N))'; 