

nQ = [3 1 -2];
v = [2,9,4]';

for i=1:1000

  nQ(3) = i/100;
  var_vec(i) = i/100;
  [rot, diffrot] = diffqn(nQ);

  res(i,1:3) = (rot*v)';
  dres(i,1:3) = (diffrot(1:3,7:9)*v)';
endfor

clf;
plot(var_vec, res(:,3))
hold on
plot(var_vec, dres(:,3), 'g+')
plot(var_vec(1:999), diff(res(:,3))./diff(var_vec'), 'r-')

[rot, diffrot] = diffqn( [1 2 3] );

rot
diffrot