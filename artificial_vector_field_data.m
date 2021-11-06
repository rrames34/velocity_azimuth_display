close all
clear


%Defining constants for Log Law
z0=0.4;   %aerodynamic roughness length  %stull fig9.6 %dense forest
k=0.35;  
d=25.9; 
u0=1;  %surface stress

%Defining constans for Ekman Spiral
ekman=0.1;

lidar_position_x=10;  %x coordinate of lidar origin
lidar_position_y=10;  %y coordinate of lidar origin


x=[0:1:20];  %x-axis coordinates 
y=[0:1:20];  %y-axis coordinates
z=[z0+d+0.01:1:60]; %z-coordinates along height

[X,Y,Z]=meshgrid(x,y,z);  %3D grid coordinates defined by x,y,z

W=0*X; %assuming wind velocity is horizontal
%defining horizontal wind velocity components
V(:,:,:)=(u0/k)*log((Z(:,:,:)-d)/z0).*sin(ekman*(Z(:,:,:)));  %x component of wind velocity 
U(:,:,:)=(u0/k)*log((Z(:,:,:)-d)/z0).*cos(ekman*(Z(:,:,:)));  %y component of wind velocity

%fetching data for vertical velocity profile
u=U(lidar_position_x,lidar_position_y,:);   
v=V(lidar_position_x,lidar_position_y,:);

figure(1)
plot(sqrt(v(:).^2+u(:).^2),z(:));  %plotting wind speed profile at lidar position (10,10)
hold on


figure(2)
plot(rad2deg(myatan(v(:),u(:))),z(:)); %plotting wind direction profile at lidar position (10,10)
hold on


a=5/100; %minimum kick of 5 percent   %
b=10/100; %maximum kick of 10 percent %

random_kick_flag="yes";  %whether to apply random kick

if random_kick_flag=="yes"
    %applying random kick to U component
    random_kick_sign=sign(2.*rand(length(x),length(y),length(z))-1);    %randomly decide whether kick increases or decreases velocity  
    random_kick_percentage=((b-a).*rand(length(x),length(y),length(z))+a);  %randomly decide the kick value between 5-10 percent %
    random_kick=1-random_kick_sign.*random_kick_percentage;             
    U=U.*(random_kick);
    %applying random kick to V component
    random_kick_sign=sign(2.*rand(length(x),length(y),length(z))-1); %randomly decide whether kick increases or decreases velocity
    random_kick_percentage=((b-a).*rand(length(x),length(y),length(z))+a);  %randomly decide the kick value between 5-10 percent %
    random_kick=1-random_kick_sign.*random_kick_percentage;
    V=V.*(random_kick);
end

%fetching data with kick for vertical velocity profile
u_kick=U(lidar_position_x,lidar_position_y,:);
v_kick=V(lidar_position_x,lidar_position_y,:);

figure(1)
plot(sqrt(v_kick(:).^2+u_kick(:).^2),z(:));
title('Wind speed vs Height(from constructed vector field)');
xlabel('Wind speed (m/s)');
ylabel('Height (m)');
grid on
grid minor
legend('before applying random kick','after applying random kick 5-10%','Location','southeast')

figure(2)
plot(rad2deg(myatan(v_kick(:),u_kick(:))),z(:));
title('Wind direction vs Height(from constructed vector field)')
xlabel('\theta (degrees)');
ylabel('Height (m)');
legend('before applying random kick','after applying random kick 5-10%','Location','southeast')
grid on
grid minor

%range ring
range=linspace(30,50,10); %radial distance from origin of virtual lidar
elevation=75;            %elevation angle in degrees
azimuth=[0:5:355];       %azimuthal angle in degrees



r_z=range*sind(elevation);  % z coordinate of range ring - height
r_x=zeros(length(azimuth),1);   
r_y=zeros(length(azimuth),1);   


r_u=zeros(length(azimuth),1);
r_v=zeros(length(azimuth),1);

radial_velocity=zeros(length(azimuth),1);




%Loop for constructing lidar data

for ring=1:length(range)
    for i=1:length(azimuth)
        MyData.range(ring,i)=range(ring);   %radial position of points on range rings
        MyData.az(ring,i)=azimuth(i);       %azimuthal position of points on range rings
        MyData.el(ring,i)=elevation;        %elevation angle of range rings
        
        %x coordinate of points on a range ring
        r_x(i)=lidar_position_x+range(ring)*cosd(elevation)*sind(azimuth(i));
        %y coordinate of points on a range ring
        r_y(i)=lidar_position_y+range(ring)*cosd(elevation)*cosd(azimuth(i));
    
        %interpolate velocity to points on a range ring
    
        r_u(i)=interp3(X,Y,Z,U,r_x(i),r_y(i),r_z(ring));    %u velocity component at points on range ring
    
        r_v(i)=interp3(X,Y,Z,V,r_x(i),r_y(i),r_z(ring));    %v velocity component at points on range ring
   
        radial_velocity(i)=dot([cosd(elevation)*sind(azimuth(i)),cosd(elevation)*cosd(azimuth(i)),sind(elevation)],[r_u(i),r_v(i),0]); %virtual lidar
        MyData.rv(ring,i)=radial_velocity(i);        %radial velocity obtained from virtual lidar
    end
end

Data(1)=MyData;

B=zeros(2,1);

A=zeros(2);

time=1;

%minimizing cost function using least squares method 
 for ring=1:size(Data(time).range,1)
        count=0;
        for i=1:size(Data(time).range,2)
             if ~isnan(Data(time).rv(ring,i))==1
                count=count+1;
                B(1)=B(1)+Data(time).rv(ring,i)*cosd(Data(time).el(ring,i))*cosd(Data(time).az(ring,i));
                B(2)=B(2)+Data(time).rv(ring,i)*cosd(Data(time).el(ring,i))*sind(Data(time).az(ring,i));
                A(1,1)=A(1,1)+cosd(Data(time).el(ring,i))*sind(Data(time).az(ring,i))*cosd(Data(time).el(ring,i))*cosd(Data(time).az(ring,i));
                A(1,2)=A(1,2)+cosd(Data(time).el(ring,i))*cosd(Data(time).az(ring,i))*cosd(Data(time).el(ring,i))*cosd(Data(time).az(ring,i));
                A(2,1)=A(2,1)+cosd(Data(time).el(ring,i))*sind(Data(time).az(ring,i))*cosd(Data(time).el(ring,i))*sind(Data(time).az(ring,i));
                A(2,2)=A(2,2)+cosd(Data(time).el(ring,i))*cosd(Data(time).az(ring,i))*cosd(Data(time).el(ring,i))*sind(Data(time).az(ring,i));
             end
        end
       
        
        output_data(time).velocity(ring,:)=pinv(A)*B;
        A=zeros(2);
        B=zeros(2,1);
 end
 
figure(3)
plot(vecnorm(output_data(time).velocity(:,:),2,2),r_z);
title('Wind speed vs Height(from VAD)');
xlabel('Wind speed (m/s)');
ylabel('Height (m)');
grid on
grid minor


figure(4)
plot(vecnorm(output_data(time).velocity(:,:),2,2),r_z);
hold on
plot(sqrt(v_kick(:).^2+u_kick(:).^2),z(:));
title('Comparison of constructed vector field and VAD results for wind speed');
xlabel('Wind speed (m/s)');
ylabel('Height (m)');
legend('from VAD','from constructed vector field','Location','southeast');
grid on
grid minor


figure(5)
plot(rad2deg(myatan(output_data(time).velocity(:,2),output_data(time).velocity(:,1))),r_z);
title('Wind direction vs Height from VAD')
xlabel('\theta (degrees)');
ylabel('Height (m)');
grid on
grid minor


figure(6)
plot(rad2deg(myatan(output_data(time).velocity(:,2),output_data(time).velocity(:,1))),r_z);
hold on
plot(rad2deg(myatan(v_kick(:),u_kick(:))),z(:));
title('Comparison of constructed vector field and VAD results for wind direction')
xlabel('\theta (degrees)');
ylabel('Height (m)');
legend('from VAD','from constructed vector field','Location','southeast');
grid on
grid minor


