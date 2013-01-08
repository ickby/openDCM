%backend ("fltk")

run('output.m');

param = 1:12;

for i=param
  clf;
  plot(parameter(i,:),results(i,:));
  hold on;
  plot(parameter(i,:), gradient(i,:), 'g+');
  plot(parameter(i,1:999), diff(results(i,:))./diff(parameter(i,:)), 'r-');
  figure;
endfor

