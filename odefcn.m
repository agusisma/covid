function dydt = odefcn(t,y,beta,gamma)

dydt    = zeros(3,1);
dydt(1) = -beta*y(1)*y(2);
dydt(2) = beta*y(1)*y(2)-gamma*y(2);
dydt(3) = gamma*y(2);
