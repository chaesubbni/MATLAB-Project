clear; clc;
rad = 0 : 0.01 : 0.2;  % 반지름
num_cone = 20;        % 원뿔 둘레 점의 개수
num_Z = length(rad);    % 원의 개수 = Z축으로 생기는 점의 개수
[X, Y, Z] = cylinder(rad(end) - rad, num_cone);   % 원뿔 만들기



rad_br = 0 : 0.004 : 0.09;  % 반지름
num_cone_br = 20;        % 원뿔 둘레 점의 개수
num_Z_br = length(rad);    % 원의 개수 = Z축으로 생기는 점의 개수
[X_br, Y_br, Z_br] = cylinder(rad_br(end) - rad_br, num_cone);   % 원뿔 만들기

radius = 0.2; % 구의 반지름
resolution = 20; % 구의 해상도 (정점의 개수)

% 구 생성
[sx, sy, sz] = sphere(resolution); % 정규화된 구의 정점 좌표 생성
x_curcle = radius * sx; % 구의 x 좌표
y_curcle = radius * sy; % 구의 y 좌표
z_curcle = radius * sz; % 구의 z 좌표

[sx, sy, sz] = sphere(resolution); % 정규화된 구의 정점 좌표 생성
x_c = radius * sx; % 구의 x 좌표
y_c = radius * sy; % 구의 y 좌표
z_c = radius * sz; % 구의 z 좌표



% 원뿔 점들을 회전변환 시킬 수 있게 1행으로 행렬을 재배열
pre_G(1,:) = reshape(X', 1, []);
pre_G(2,:) = reshape(Y', 1, []);
pre_G(3,:) = reshape(Z', 1, []);
pre_G(4,:) = 0;

pre_G_br(1,:) = reshape(X_br', 1, []);
pre_G_br(2,:) = reshape(Y_br', 1, []);
pre_G_br(3,:) = reshape(Z_br', 1, []);
pre_G_br(4,:) = 0;

pre_G_gr(1,:) = reshape(x_curcle',1,[]);
pre_G_gr(2,:) = reshape(y_curcle',1,[]);
pre_G_gr(3,:) = reshape(z_curcle',1,[]);
pre_G_gr(4,:) = 0;

pre_G_gr_c(1,:) = reshape(x_curcle',1,[]);
pre_G_gr_c(2,:) = reshape(y_curcle',1,[]);
pre_G_gr_c(3,:) = reshape(z_curcle',1,[]);
pre_G_gr_c(4,:) = 0;

f= 0.05;                                % 주파수
dt = 1;                               % 시간 간격
t = 0:dt:15;                            % 시간
num_temp = 1;            % 기둥

r = 10*2;                                % 각 점의 각도 가중치 변화율 요소
vel_conv = 20;                           % 시간에 따른 수렴 속도

psi2 = (30/180*pi) * (sin(2*pi*f*t));   % +- 30도의 크기의 사인파 생성 (각 점의 회전 각도)

figure(1); clf
hold on



for mm = 1:num_temp
    Temp{mm} =  surf(X, Y, Z,'FaceColor', [0.5 0.35 0.05],  'edgecolor', 'none');    % 중심 기둥을 그립니다.
    grass_c{mm} = surf(x_curcle, y_curcle, z_curcle, 'FaceColor', [0,0.5, 0], 'EdgeColor', 'none'); % 중심 가지의 나뭇잎을 그립니다.

    grass{1}{mm} = surf(x_curcle, y_curcle, z_curcle, 'FaceColor', [0,0.5, 0], 'EdgeColor', 'none');
    grass{2}{mm} = surf(x_curcle, y_curcle, z_curcle, 'FaceColor', [0,0.5, 0], 'EdgeColor', 'none');
    grass{3}{mm} = surf(x_curcle, y_curcle, z_curcle, 'FaceColor', [0,0.5, 0], 'EdgeColor', 'none');
    grass{4}{mm} = surf(x_curcle, y_curcle, z_curcle, 'FaceColor', [0,0.5, 0], 'EdgeColor', 'none');
    grass{5}{mm} = surf(x_curcle, y_curcle, z_curcle, 'FaceColor', [0,0.5, 0], 'EdgeColor', 'none');
    
% 잔디의 중심 값 랜덤으로 형성
    x_loc(mm) = 0; % 본가지 x위치
    y_loc(mm) = 0; % 본가지 y위치
    z_high(mm) =5; % 본가지 z좌표 높이
    brench_num(mm) = 5; % 5개의 곁가지 생성
    
    
    Brench{1}{mm} = surf(X_br,Y_br,Z_br,'FaceColor', [0.5 0.35 0.05],  'edgecolor', 'none');    % 가지
    pos{1}{mm} = 10;                             % 가지의 위치의 Z값 설정
    the{1}{mm} = -60;
    brench_high{1}{mm} = 5;
    
    Brench{2}{mm} = surf(X_br,Y_br,Z_br,'FaceColor',[0.5 0.35 0.05],  'edgecolor', 'none');    % 가지 
    pos{2}{mm} = 8;                             % 가지의 위치의 Z값 설정
    the{2}{mm} = 60;
    brench_high{2}{mm} = 5;
    
    Brench{3}{mm} = surf(X_br,Y_br,Z_br,'FaceColor',[0.5 0.35 0.05],  'edgecolor', 'none');    % 가지
    pos{3}{mm} = 7;                             % 가지의 위치의 Z값 설정
    the{3}{mm} = 60;
    brench_high{3}{mm} = 5;
    
    Brench{4}{mm} = surf(X_br,Y_br,Z_br,'FaceColor', [0.5 0.35 0.05],  'edgecolor', 'none');    % 가지 
    pos{4}{mm} = 5;                             % 가지의 위치의 Z값
    the{4}{mm} = -60;
    brench_high{4}{mm} = 5;
    
    Brench{5}{mm} = surf(X_br,Y_br,Z_br,'FaceColor', [0.5 0.35 0.05],  'edgecolor', 'none');    % 가지 
    pos{5}{mm} = 9;                             % 가지의 위치의 Z값 설정
    the{5}{mm} = -60;
    brench_high{5}{mm} = 5;
    
    ax_wind = plot(0, 0,'-b','LineWidth',2);
   
end

box on; grid on;
xlim([-5 5]); ylim([-5 5]); zlim([0 9]);    
xlabel('x'); ylabel('y'); zlabel('z')     
view(-20, 30);

view_count = 0;
for pp = [0 30 -30 -60]       % 바람의 방향이 0도 30도 45도 85도로 for문 반복
    t_vector = dt * ones(num_temp,1);
    
    wind_x_c = [];
    wind_y_c = [];
    
    wind_direction = pp / 180 * pi;     % 바람의 방향
    
    xy_dist_old = ones(num_temp,1) * 10;
    
    % 바람의 위치와 각 나무와의 거리를 계산하여 바람이 부는 점에서 가까운쪽부터 나무가 흔들리기 시작함
    for j = 1 : length(t)
        wind_x = t(j) * 2 - 5;
        wind_y = tan(wind_direction) * wind_x;
        
        for mm = 1 : num_temp
            psi(mm) = psi2(round(t_vector(mm) / dt));
            
            xy_dist_new(mm) = sqrt((x_loc(mm) - wind_x)^2 + (y_loc(mm) - wind_y)^2);
            if xy_dist_new(mm) < xy_dist_old(mm)
                t_vector(mm) = dt;
            end
        end
        xy_dist_old = xy_dist_new;
        
        % 나무가 흔들리는 모션을 넣어주는 for문
        for j = 1 : num_Z
            weight_a = (exp((j - num_Z) / r));      % 바람 방향에 따른 회전변환
            Rz1 = [cos(wind_direction)   -sin(wind_direction)  0    0;
                   sin(wind_direction)   cos(wind_direction)   0    0;
                            0                    0             1    0;
                            0                    0             0    1];
            
                        
                        
            for mm = 1 : num_temp
                weight_b{mm} = exp(-t_vector(mm) / vel_conv);       % 흔들림에 따른 회전변환
                
                Rx{mm} = [cos(weight_b{mm} * weight_a * psi(mm))   0   sin(weight_b{mm} * weight_a * psi(mm))       0;
                                           0                       1                       0                        0;
                         -sin(weight_b{mm} * weight_a * psi(mm))   0   cos(weight_b{mm} * weight_a * psi(mm))       0;
                                           0                       0                       0                        1];
                
                    for l = 1 : brench_num(mm)
                    % 가지가 누워있게 하는 회전 행렬
                    R = [1            0              0              0;
                         0      cosd(the{l}{mm})  -sind(the{l}{mm})       0;
                         0      sind(the{l}{mm})   cosd(the{l}{mm})       0;
                         0            0              0              1];
                                           
                    R2 = [cosd(the{l}{mm})  0          sind(the{l}{mm})              0;
                          0                 1                 0                      0;
                          -sind(the{l}{mm}) 0          cosd(the{l}{mm})              0;
                          0                 0                 0                      1];
                     
                     
                     for k = 1 : num_cone
                        y1{mm}(:,num_cone*(j - 1) + k) = Rz1 * Rx{mm} *pre_G(:,num_cone*(j - 1) + k);  % 나무에 회전변환을 적용
                        Grass_c{mm}(:,num_cone*(j - 1) + k) = Rz1 * Rx{mm} *pre_G_gr(:,num_cone*(j - 1) + k);
                           
                        mo = mod(l,2);
                        if (mo == 0)
                           br{l}{mm}(:,num_cone*(j - 1) + k) =   Rx{mm}* Rz1 * R * pre_G_br(:,num_cone*(j - 1) + k);    % 가지에 회전변환을 적용
                                                   
                        elseif (l == 5)
                           br{l}{mm}(:,num_cone*(j - 1) + k) =   Rx{mm}* Rz1 * R * pre_G_br(:,num_cone*(j - 1) + k);    % 가지에 회전변환을 적용
                                          
                        else (mo ~= 0)
                           br{l}{mm}(:,num_cone*(j - 1) + k) =   Rx{mm}* Rz1 * R2 * pre_G_br(:,num_cone*(j - 1) + k);    % 가지에 회전변환을 적용
                           
                        end
                        Grass{l}{mm}(:,num_cone*(j - 1) + k) =  Rz1 * Rx{mm} * pre_G_gr(:,num_cone*(j - 1) + k);    
                        
                    end
                end
                
            end
        end
   
        
        wind_x_c = [wind_x_c wind_x];
        wind_y_c = [wind_y_c wind_y];
        
        
        for mm = 1: num_temp
            
            % 나무를 회전변환하여 생긴값을 원래 나무 그림에 업데이트 시켜줌
            Temp{mm}.XData = reshape(y1{mm}(1,:), num_cone, []) + x_loc(mm);
            Temp{mm}.YData = reshape(y1{mm}(2,:), num_cone, []) + y_loc(mm);
            Temp{mm}.ZData = reshape(y1{mm}(3,:), num_cone, []) * z_high(mm);
            grass_c{mm}.XData = reshape(Grass_c{mm}(1,:), num_cone, []) + Temp{mm}.XData(end,end);
            grass_c{mm}.YData = reshape(Grass_c{mm}(2,:), num_cone, []) + Temp{mm}.YData(end,end);
            grass_c{mm}.ZData = reshape(Grass_c{mm}(3,:), num_cone, []) + Temp{mm}.ZData(end,end);
            
            
            % 가지를 회전변환하여 생긴 값을 원래 가지 그림에 업데이트 시켜줌
            for j=1:brench_num(mm)
                Brench{j}{mm}.XData = reshape(br{j}{mm}(1,:), num_cone, []) / 2 + (Temp{mm}.XData(1, pos{j}{mm}) + Temp{mm}.XData(11, pos{j}{mm})) / 2;
                Brench{j}{mm}.YData = reshape(br{j}{mm}(2,:), num_cone, []) / 2 + (Temp{mm}.YData(1, pos{j}{mm}) + Temp{mm}.YData(11, pos{j}{mm})) / 2;
                Brench{j}{mm}.ZData = reshape(br{j}{mm}(3,:), num_cone, []) * brench_high{j}{mm} + Temp{mm}.ZData(1, pos{j}{mm});
                
                grass{j}{mm}.XData = reshape(Grass{j}{mm}(1,:), num_cone, []) + Brench{j}{mm}.XData(end,end);
                grass{j}{mm}.YData = reshape(Grass{j}{mm}(2,:), num_cone, []) + Brench{j}{mm}.YData(end,end);
                grass{j}{mm}.ZData = reshape(Grass{j}{mm}(3,:), num_cone, []) + Brench{j}{mm}.ZData(end,end);
            end
            ax_wind.XData = wind_x_c;
            ax_wind.YData = wind_y_c;
            
            t_vector(mm) = t_vector(mm) + dt;   % 시간의 증가
            drawnow limitrate                   % 그림을 좀더 부드럽게 그려지게 함
           
        end
        
        
        
        % 나무를 보는 각도가 조금씩 바뀜
        view_count = view_count + 1;
        view(70*sind(2*pi*0.05*view_count), 12*cosd(2*pi*0.05*view_count)+25)
        
    end
end