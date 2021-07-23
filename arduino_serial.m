clc
clear 
close all
delete(instrfind({'Port'},{'COM5'}));
s = serial('COM5', 'BaudRate', 115200);
fopen(s);

for i = 1:inf
    a=str2num(fscanf(s));
        if length(a)==5
            b(i,1:5) = a;
            plot(b(:,4),b(:,2),'LineWidth',2)
            hold on
            pause(0.000001);
        end
end