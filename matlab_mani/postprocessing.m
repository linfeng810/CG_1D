% this is postprocessing script
% when using, please copy the '../results.out' file into the following 
% bracket: (whole file)
data = [];

t = data(:,1);
phi = data(:,2:end);
nnod = length(phi(1,:));
dx = 1/(nnod-1);
x = 0:dx:1;
figure(1);clf(1);
for i = 1:length(t)
    plot(x,phi(i,:));
    xlabel('x'); ylabel('p');
    title({'time:',num2str(t(i))})
    pause(0.1);
end