clear
close all
clc

addpath(genpath('J:\Data\Matlab'));             
%% Load file
[spreadsheet,sheetpath] = uigetfile('*.txt','Select a txt file'); % Load the excel spreadsheet
cd(sheetpath);
t = readtable(spreadsheet); % 1: number of segments, 2: index of seg type(1:soma; 2:axon; 3:dend),3:x,4:y,5:z, 6:dia, 7:parent seg (-1: origin)
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

pixel_size = par.(11)(index); % scale factor
%% visulization
count_x = 0;
count_y = 0;
i = star_pixel;
k = 1;

t.(8) = t.(3)-soma_x;    % x_translation
t.(9) = t.(4)-soma_y;    % y_translation

t.(10) = cos(theta)*t.(8) - sin(theta)*t.(9);   % x_rotation
t.(11) = sin(theta)*t.(8) + cos(theta)*t.(9);   % y_rotation
t.(10) = t.(10)*pixel_size;
t.(11) = t.(11)*pixel_size;
t.(12)(:) = nan; % x_intersection, 1 or 0
t.(13)(:) = nan; % y_intersection, 1 or 0

% while i < end_pixel
%     for j = 1:length(index_axon)       
%         b1 = t.(10)(index_axon(j));
%         b3 = t.(11)(index_axon(j));
%                 
%         temp = t.(7)(index_axon(j));
%         temp_index = find(t.(1) == temp);
%         
%         b2 = t.(10)(temp_index);
%         b4 = t.(11)(temp_index);
%                
%         if (b1<=i && b2>i) || (b2<=i && b1>i)
%             count_x = count_x+1;
%             t.(12)(index_axon(j)) = 1;
%         elseif t.(12)(index_axon(j)) ~= 1
%             t.(12)(index_axon(j)) = 0;         
%         end
%         if (b3<=i && b4>i) || (b4<=i && b3>i)
%             count_y = count_y+1;
%             t.(13)(index_axon(j)) = 1;
%         elseif t.(13)(index_axon(j)) ~= 1
%             t.(13)(index_axon(j)) = 0;            
%         end
%     end
%     k
%     summary(k,1) = i;
%     summary(k,2) = count_x;
%     summary(k,3) = count_y;
%     
%     k = k+1;
%     i = i+tenmicron_pixelnumber;
%     count_x = 0;
%     count_y = 0;   
% end
% 
% summary(:,4) = -800:10:800;    % distance from soma (um)
% summary(:,5) = 100*summary(:,2)/sum(summary(:,2));  % normalized tengential axonal density (%)
% summary(:,6) = 100*summary(:,3)/sum(summary(:,3));  % normalized vertical axonal density (%)
% intersec_x = find(t.(12) == 1);
% intersec_y = find(t.(13) == 1);

fig = figure;
for a = 1:size(t,1)
    temp = t.(1)(a);
    temp_index = find(t.(7) == temp);
    if ~isempty(temp_index)
        for aa = 1:length(temp_index)
            x = t.(10)(temp_index(aa))-t.(10)(a);
            y = t.(11)(temp_index(aa))-t.(11)(a);
            phi = atan(y/x);        
        
            x1 = cos(phi+0.5*pi)*(t.(6)(a)/2) + t.(10)(a);
            y1 = sin(phi+0.5*pi)*(t.(6)(a)/2) + t.(11)(a);
            x2 = cos(phi-0.5*pi)*(t.(6)(a)/2) + t.(10)(a);
            y2 = sin(phi-0.5*pi)*(t.(6)(a)/2) + t.(11)(a);
        
            x3 = cos(phi-0.5*pi)*(t.(6)(temp_index(aa))/2) + t.(10)(temp_index(aa));
            y3 = sin(phi-0.5*pi)*(t.(6)(temp_index(aa))/2) + t.(11)(temp_index(aa));
            x4 = cos(phi+0.5*pi)*(t.(6)(temp_index(aa))/2) + t.(10)(temp_index(aa));
            y4 = sin(phi+0.5*pi)*(t.(6)(temp_index(aa))/2) + t.(11)(temp_index(aa));
            if x1 == x2 && y1 == y2
                pgon = polyshape([x1 x3 x4],[y1 y3 y4]);
            elseif x3 == x4 && y3 == y4
                pgon = polyshape([x1 x2 x3],[y1 y2 y3]);
            else
                pgon = polyshape([x1 x2 x3 x4],[y1 y2 y3 y4]);
            end
            if t.(2)(a) == 1 || t.(2)(a) == 3
                plot(pgon,'FaceColor','k','FaceAlpha',1,'LineStyle','none') % soma n dendrites
            elseif t.(2)(a) == 2
                plot(pgon,'FaceColor','r','FaceAlpha',1,'LineStyle','none') % axons
            end
            hold on
        end
    end   
    a
end
title(strcat(erase(spreadsheet,'.txt'),' translated and rotated tracing'))
x0=700;
y0=100;
x_max = 1000; % micrometer
y_max = 1000; % micrometer
% x_max = abs(max(t.(10)));
% y_max = abs(max(t.(11)));
xy_max = max(x_max, y_max);
xlim([-xy_max,xy_max]);
ylim([-xy_max,xy_max]);
width = 800;
height = 800;
set(gcf,'position',[x0,y0,width,height])

%% output file
cd(sheetpath);
exportgraphics(fig,strcat(erase(spreadsheet,'.txt'),'.pdf'),'ContentType','vector'); %salva in pdf
spreadsheet
disp('job done');
sound(sin(1:3000));