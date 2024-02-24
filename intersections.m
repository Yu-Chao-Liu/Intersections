clear
close all
clc

addpath(genpath('J:\Data\Matlab'));             
%% Load file
[spreadsheet,sheetpath] = uigetfile('*.txt','Select a txt file'); % Load the excel spreadsheet
cd(sheetpath);
t = readtable(spreadsheet); % number of segments, index of seg type(1:soma; 2:axon; 3:dend),x,y,z, dia, parent seg (-1: origin)
%% loda parameters
[parameters,par_path] = uigetfile('*.xlsx','Select a excel file'); % Load the summary excel spreadsheet
cd(par_path);
par = readtable(parameters);

%% define parameters
date = str2double(extractBefore(spreadsheet,9));
cellnumber = extractAfter(spreadsheet,9);
cellnumber = str2double(cellnumber(1:end-7));

index = find(par.(2) == date & par.(3) == cellnumber);

celltype = par.(4){index};
pi_s = par.(7)(index);
theta = pi_s*pi;
tenmicron_pixelnumber = par.(12)(index);
soma_x = par.(13)(index);
soma_y = par.(14)(index);
star_pixel = par.(16)(index);
end_pixel = par.(18)(index);

index_axon = find(t.(2) == 2);
index_somandend = find(t.(2) == 1 | t.(2) == 3);

%% analysis
count_x = 0;
count_y = 0;
i = star_pixel;
k = 1;

t.(8) = t.(3)-soma_x;    % x_translation
t.(9) = t.(4)-soma_y;    % y_translation

t.(10) = cos(theta)*t.(8) - sin(theta)*t.(9);   % x_rotation
t.(11) = sin(theta)*t.(8) + cos(theta)*t.(9);   % y_rotation
t.(12)(:) = nan; % x_intersection, 1 or 0
t.(13)(:) = nan; % y_intersection, 1 or 0

while i < end_pixel
    for j = 1:length(index_axon)       
        b1 = t.(10)(index_axon(j));
        b3 = t.(11)(index_axon(j));
                
        temp = t.(7)(index_axon(j));
        temp_index = find(t.(1) == temp);
        
        b2 = t.(10)(temp_index);
        b4 = t.(11)(temp_index);
               
        if (b1<=i && b2>i) || (b2<=i && b1>i)
            count_x = count_x+1;
            t.(12)(index_axon(j)) = 1;
        elseif t.(12)(index_axon(j)) ~= 1
            t.(12)(index_axon(j)) = 0;         
        end
        if (b3<=i && b4>i) || (b4<=i && b3>i)
            count_y = count_y+1;
            t.(13)(index_axon(j)) = 1;
        elseif t.(13)(index_axon(j)) ~= 1
            t.(13)(index_axon(j)) = 0;            
        end
    end
    k
    summary(k,1) = i;
    summary(k,2) = count_x;
    summary(k,3) = count_y;
    
    k = k+1;
    i = i+tenmicron_pixelnumber;
    count_x = 0;
    count_y = 0;   
end

summary(:,4) = -800:10:800;    % distance from soma (um)
summary(:,5) = 100*summary(:,2)/sum(summary(:,2));  % normalized tengential axonal density (%)
summary(:,6) = 100*summary(:,3)/sum(summary(:,3));  % normalized vertical axonal density (%)
intersec_x = find(t.(12) == 1);
intersec_y = find(t.(13) == 1);

%% visulization
figure    % original tracing
plot(t.(3)(index_somandend),t.(4)(index_somandend),'k.')
hold on
plot(t.(3)(index_axon),t.(4)(index_axon),'r.')
title('original tracing')
x0=700;
y0=300;
x_max = abs(max(t.(3)));
y_max = abs(max(t.(4)));
xy_max = max(x_max, y_max);
xlim([0,xy_max]);
ylim([0,xy_max]);
width = 500;
height = 500;
set(gcf,'position',[x0,y0,width,height])


figure    % translated and rotated tracing + x_intersections
g_y=[star_pixel:tenmicron_pixelnumber:end_pixel]; % user defined grid Y [start:spaces:end]
g_x=[star_pixel:tenmicron_pixelnumber:end_pixel]; % user defined grid X [start:spaces:end]
gray = 0.95;
grayColor = [gray gray gray];
for a=1:length(g_x)
   plot([g_x(a) g_x(a)],[g_y(1) g_y(end)],'color', grayColor) %y grid lines
   hold on    
end
h = zeros(1,3);
h(1) = plot(t.(10)(index_somandend),t.(11)(index_somandend),'k.','DisplayName', char({'soma n dend'}));
h(2) = plot(t.(10)(index_axon),t.(11)(index_axon),'r.','DisplayName', char({'axon'}));
h(3) = plot(t.(10)(intersec_x),t.(11)(intersec_x),'bo','DisplayName', char({'x intersects'})); % x_intersections
legend(h(1:3),'Location','Best');

title('translated and rotated tracing')
x0=700;
y0=300;
x_max = abs(max(t.(10)));
y_max = abs(max(t.(11)));
xy_max = max(x_max, y_max);
xlim([-xy_max,xy_max]);
ylim([-xy_max,xy_max]);
width = 500;
height = 500;
set(gcf,'position',[x0,y0,width,height])


figure    % translated and rotated tracing + y_intersections
for a=1:length(g_y)
   plot([g_x(1) g_x(end)],[g_y(a) g_y(a)],'color', grayColor) %x grid lines
   hold on    
end
h = zeros(1,3);
h(1) = plot(t.(10)(index_somandend),t.(11)(index_somandend),'k.','DisplayName', char({'soma n dend'}));
h(2) = plot(t.(10)(index_axon),t.(11)(index_axon),'r.','DisplayName', char({'axon'}));
h(3) = plot(t.(10)(intersec_y),t.(11)(intersec_y),'go','DisplayName', char({'y intersects'})); % y_intersections
legend(h(1:3),'Location','Best');

title('translated and rotated tracing')
xlim([-xy_max,xy_max]);
ylim([-xy_max,xy_max]);
set(gcf,'position',[x0,y0,width,height])



%% output file
cd(sheetpath);
print(1,'-dpng','-r300',strcat(erase(spreadsheet,'.txt'),' ori_tracing')) %save plot as png (looks better)
print(2,'-dpng','-r300',strcat(erase(spreadsheet,'.txt'),' x_intersec')) %save plot as png (looks better)
print(3,'-dpng','-r300',strcat(erase(spreadsheet,'.txt'),' y_intersec')) %save plot as png (looks better)
writematrix(summary,strcat(erase(spreadsheet,'.txt'),'.xlsx'))

disp('job done');
sound(sin(1:3000));