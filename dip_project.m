%% Topic : image enhancement using energy equalisation with clipping limit
clc;
clear;
close all;
% Step1 : input RGB
image_rgb= imread('OPI (2).jpg'); %Reading an RGB input image
image_rgb= double(image_rgb); % converting data type into double
image_red= image_rgb(:,:,1); % Red Component
image_green= image_rgb(:,:,2); % Green Component
image_blue= image_rgb(:,:,3); % Blue Component
[M,N]= size(image_red); % Size of image
L= 256; % Number of gray levels
l= 0:L-1;
%% Step 2 : energy curve
B_red= zeros(M,N,L);
B_green= zeros(M,N,L);
B_blue= zeros(M,N,L);
for p= 0:L-1
 % Generation of Bg Matrix
B_red(:,:,p+1)= 2*(image_red > p)-1;
B_green(:,:,p+1)= 2*(image_green > p)-1;
B_blue(:,:,p+1)= 2*(image_blue > p)-1;
end
N_d= ones(3); % for 3*3 neighbourhood
N_d((length(N_d)+1)/2,(length(N_d)+1)/2)= 0; % For center element 0
C= ones(M,N); % Generation of C matrix which is identity matrix
E_red= zeros(1,L);
E_green= zeros(1,L);
E_blue= zeros(1,L);
for p= 1:L
% Energy Calculation at a particular intensity level
E_red(p)= sum(sum(-1*B_red(:,:,p).*conv2(B_red(:,:,p),N_d,'same') + C.*conv2(C,N_d,'same'))); % for red
E_green(p)= sum(sum(-1*B_green(:,:,p).*conv2(B_green(:,:,p),N_d,'same') + C.*conv2(C,N_d,'same'))); % for green
E_blue(p)= sum(sum(-1*B_blue(:,:,p).*conv2(B_blue(:,:,p),N_d,'same') + C.*conv2(C,N_d,'same'))); % for blue
end
%% Step 3 : clippling limit apply
C_clip_red = (median(E_red) + mean(E_red))/2; % Clipping limit of Red component
C_clip_green = (median(E_green) + mean(E_green))/2; % Clipping limit of Green component
C_clip_blue = (median(E_blue) + mean(E_blue))/2; % Clipping limit of Blue component
E_red = E_red;
E_green = E_green;
E_blue = E_blue;
E_red(E_red >= C_clip_red) = C_clip_red; % Clipped Energy curve for Red channel
E_green(E_green >= C_clip_green) = C_clip_green; % Clipped Energy curve for Green channel
E_blue(E_blue >= C_clip_blue) = C_clip_blue; % Clipped Energy curve for Blue channel
%% Step 4 : standard deviation
l_mean_red = sum(l.*E_red)/sum(E_red);
l_mean_green = sum(l.*E_green)/sum(E_green);
l_mean_blue = sum(l.*E_blue)/sum(E_blue);
SD_red = sqrt(sum(((l - l_mean_red).^2).*E_red)/sum(E_red)); % standard deviation of red component
SD_green = sqrt(sum(((l - l_mean_green).^2).*E_green)/sum(E_green)); % standard deviation of green component
SD_blue = sqrt(sum(((l - l_mean_blue).^2).*E_blue)/sum(E_blue)); % standard deviation of blue component
SD_red = round(SD_red);  % round-off value of standard deviation for red component 
SD_green = round(SD_green); % round-off value of standard deviation for green component 
SD_blue = round(SD_blue); % round-off value of standard deviation for blue component 
%%  component partitioning 

L_0_r = min(image_red,[],'all'); % Minimum intensity level for Red component
L_L_1_r = max(image_red,[],'all'); % Maximum intensity level for Red component
L_low_red = L_0_r + SD_red; % Lower limit for Red component
L_high_red = L_L_1_r - SD_red; % Upper limit for Red component

L_0_g = min(image_green,[],'all'); % Minimum intensity level for Green component
L_L_1_g = max(image_green,[],'all'); % Maximum intensity level for Green component
L_low_green = L_0_g + SD_green; % Lower limit for Green component
L_high_green = L_L_1_g - SD_green; % Upper limit for Green component
L_0_b = min(image_blue,[],'all'); % Minimum intensity level for Blue component
L_L_1_b = max(image_blue,[],'all'); % Maximum intensity level for Blue component
L_low_blue = L_0_b + SD_blue; % Lower limit for Blue component
L_high_blue = L_L_1_b - SD_blue; % Upper limit for Blue component
%% Step 6 : pdf and cdf calculation for each component
% PDF for energy curves of Red component
pdf_L_red = E_red(1:L_low_red+1)/sum(E_red(1:L_low_red+1));
pdf_M_red = E_red(L_low_red+2:L_high_red+1)/sum(E_red(L_low_red+2:L_high_red+1));
pdf_U_red = E_red(L_high_red+2:L)/sum(E_red(L_high_red+2:L));
% CDF for energy curves of Red component
cdf_L_red = cumsum(pdf_L_red);
cdf_M_red = cumsum(pdf_M_red);
cdf_U_red = cumsum(pdf_U_red);
% PDF for energy curves of Green component
pdf_L_green = E_green(1:L_low_green+1)/sum(E_green(1:L_low_green+1));
pdf_M_green = E_green(L_low_green+2:L_high_green+1)/sum(E_green(L_low_green+2:L_high_green+1));
pdf_U_green = E_green(L_high_green+2:L)/sum(E_green(L_high_green+2:L));
% CDF for energy curves of Green component
cdf_L_green = cumsum(pdf_L_green);
cdf_M_green = cumsum(pdf_M_green);
cdf_U_green = cumsum(pdf_U_green);
% PDF for energy curves of Blue component
pdf_L_blue = E_blue(1:L_low_blue+1)/sum(E_blue(1:L_low_blue+1));
pdf_M_blue = E_blue(L_low_blue+2:L_high_blue+1)/sum(E_blue(L_low_blue+2:L_high_blue+1));
pdf_U_blue = E_blue(L_high_blue+2:L)/sum(E_blue(L_high_blue+2:L));
% CDF for energy curves of Blue component
cdf_L_blue = cumsum(pdf_L_blue);
cdf_M_blue = cumsum(pdf_M_blue);
cdf_U_blue = cumsum(pdf_U_blue);
%% Transfer Function Generation for each energy curves of Red component
Tf_L_red = L_low_red*cdf_L_red;
Tf_M_red = (L_low_red +1) + (L_high_red - L_low_red +1)*cdf_M_red;
Tf_U_red = (L_high_red +1) + (L - L_high_red +1)*cdf_U_red;
% Final Transfer Function for Red component
Tf_red= [Tf_L_red, Tf_M_red, Tf_U_red];
%% Transfer Function Generation for each energy curves for Green component
Tf_L_green= L_low_green*cdf_L_green;
Tf_M_green= (L_low_green +1) + (L_high_green - L_low_green +1)*cdf_M_green;
Tf_U_green= (L_high_green +1) + (L - L_high_green +1)*cdf_U_green;
% Final Transfer Function for Green component
Tf_green= [Tf_L_green, Tf_M_green, Tf_U_green];
%% Transfer Function Generation for each energy curves for Blue component
Tf_L_blue= L_low_blue*cdf_L_blue;
Tf_M_blue= (L_low_blue +1) + (L_high_blue - L_low_blue +1)*cdf_M_blue;
Tf_U_blue= (L_high_blue +1) + (L - L_high_blue +1)*cdf_U_blue;
% Final Transfer Function for Blue component
Tf_blue= [Tf_L_blue, Tf_M_blue, Tf_U_blue];
%% Step 7: output RGB image components
image_out_red = image_red;
image_out_green = image_green;
image_out_blue = image_blue;
% Now applying transfer function on RGB
for p= 0:L-1
image_out_red(image_out_red == p) = Tf_red(p+1);
image_out_green(image_out_green == p) = Tf_green(p+1);
image_out_blue(image_out_blue == p) = Tf_blue(p+1);
end
% output energy curve for RGB
B_out_red= zeros(M,N,L);
B_out_green= zeros(M,N,L);
B_out_blue= zeros(M,N,L);
for p= 0:L-1
 %Generation of Bg Matrix for output RGB image
B_out_red(:,:,p+1)= 2*(image_out_red > p)-1;
B_out_green(:,:,p+1)= 2*(image_out_green > p)-1;
B_out_blue(:,:,p+1)= 2*(image_out_blue > p)-1;
end
E_out_red= zeros(1,L);
E_out_green= zeros(1,L);
E_out_blue= zeros(1,L);
for p= 1:L
 % Energy Calculation for RGB
E_out_red(p)= sum(sum(-1*B_out_red(:,:,p).*conv2(B_out_red(:,:,p),N_d,'same') + C.*conv2(C,N_d,'same')));
E_out_green(p)= sum(sum(-1*B_out_green(:,:,p).*conv2(B_out_green(:,:,p),N_d,'same') + C.*conv2(C,N_d,'same')));
E_out_blue(p)= sum(sum(-1*B_out_blue(:,:,p).*conv2(B_out_blue(:,:,p),N_d,'same') + C.*conv2(C,N_d,'same')));
end
%% Step 8 : final image output
image_out(:,:,1)= image_out_red;
image_out(:,:,2)= image_out_green;
image_out(:,:,3)= image_out_blue;
image_out= uint8(image_out); 
%imwrite(image_out,'enhanced_image.png');
%% Step 9 : plots for each parts
image_rgb= uint8(image_rgb); % Converting data type into uint8 
image_red= uint8(image_red); 
image_green= uint8(image_green); 
image_blue= uint8(image_blue);
image_out_red= uint8(image_out_red); 
image_out_green= uint8(image_out_green); 
image_out_blue= uint8(image_out_blue);
figure(1)
image(image_rgb);
title('Input Image');
figure(2)
image(image_out);
title('Output Image');
figure(3)
hold on
p(1)= plot(l,E_red,'-r','DisplayName','Input energy_RED');
p(2)= plot(l,E_green,'-g','DisplayName','Input energy_GREEN');
p(3)= plot(l,E_blue,'-b','DisplayName','Input energy_BLUE');
p(4)= plot(l,E_red,'--r','DisplayName','Clipped energy_RED');
p(5)= plot(l,E_green,'--g','DisplayName','Clipped energy_GREEN');
p(6)= plot(l,E_blue,'--b','DisplayName','Clipped energy_BLUE');
xlabel('Intensity');
ylabel('Energy');
title('Input and Clipped Energy curves');
xline(L_low_red,'-.r'); xline(L_high_red,'-.r');
xline(L_low_green,'-.g'); xline(L_high_green,'-.g');
xline(L_low_blue,'-.b'); xline(L_high_blue,'-.b');
legend(p(1:6));
figure(4)
hold on
plot(l,Tf_red,'-r');
plot(l,Tf_green,'-g');
plot(l,Tf_blue,'-b');
xlabel('Input Intensity values');
ylabel('Output Intensity values');
title('Transfer functions');
legend('R','G','B')
figure(5)
hold on
plot(l,E_out_red,'-r');
plot(l,E_out_green,'-g');
plot(l,E_out_blue,'-b');
xlabel('Intensity');
ylabel('Energy');
title('Output Energy curves');
legend('Output energy_R','Output energy_G','Output energy_B')
figure(6)
hold on
plot(l,histogram(image_red),'-r');
plot(l,histogram(image_green),'-g');
plot(l,histogram(image_blue),'-b');
xlabel('Intensity values');
ylabel('Number of pixels');
title('Input Histogram');
legend('R','G','B')
figure(7)
hold on
plot(l,histogram(image_out_red),'-r');
plot(l,histogram(image_out_green),'-g');
plot(l,histogram(image_out_blue),'-b');
xlabel('Intensity values');
ylabel('Number of pixels');
title('Output Histogram');
legend('R','G','B')