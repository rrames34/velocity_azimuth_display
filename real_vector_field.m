close all
clear

load('Data_for_VAD (1).mat')
B=zeros(2,1);

A=zeros(2);


T=51;


%for time=1:size(Data,2)
for time=1:T
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
        if count>1      
            output_data(time).velocity(ring,:)=pinv(A)*B;
            
            if sqrt((output_data(time).velocity(ring,1))^2+(output_data(time).velocity(ring,2))^2)<15
                output_data(time).wind_speed(ring,1)=sqrt((output_data(time).velocity(ring,1))^2+(output_data(time).velocity(ring,2))^2);
            else 
                output_data(time).wind_speed(ring,1)=NaN;
            end
            output_data(time).height(ring,1)=Data(time).range(ring,1)*sind(Data(time).el(ring,1));
            %wind_speed(ring,1)=vecnorm(velocity(ring,:));
            output_data(time).wind_direction(ring,1)=myatan(output_data(time).velocity(ring,2),output_data(time).velocity(ring,1));            
            A=zeros(2);
            B=zeros(2,1);
            if time>1
                if ~isnan(output_data(time).wind_speed(ring,1))==1 && ~isnan(output_data(time-1).wind_speed(ring,1))==1
                    percent_change_increase=100*(output_data(time).wind_speed(ring,1)-output_data(time-1).wind_speed(ring,1))/output_data(time-1).wind_speed(ring,1);
                    percent_change_decrease=100*(output_data(time-1).wind_speed(ring,1)-output_data(time).wind_speed(ring,1))/output_data(time).wind_speed(ring,1);
            
                    if percent_change_increase>300 
                        output_data(time).wind_speed(ring,1)=output_data(time-1).wind_speed(ring,1);
                    elseif percent_change_decrease>300
                        output_data(time-1).wind_speed(ring-1,1)=output_data(time).wind_speed(ring,1);
                    end
                end
            end
        else
            
            if time>1
            output_data(time).velocity(ring,:)=output_data(time-1).velocity(ring,:);
            output_data(time).wind_speed(ring,1)=output_data(time-1).wind_speed(ring,1);
            output_data(time).wind_direction(ring,1)=output_data(time-1).wind_speed(ring,1);
            else
            output_data(time).velocity(ring,:)=[NaN NaN];
            output_data(time).wind_speed(ring,1)=NaN;
            output_data(time).wind_direction(ring,1)=NaN;
            end
        
            output_data(time).height(ring,1)=Data(time).range(ring,1)*sind(Data(time).el(ring,1));
            
            A=zeros(2);
            B=zeros(2,1);
        end
        
    end
    %output_data(time).velocity=velocity;
    %output_data(time).wind_speed=wind_speed;
    %output_data(time).wind_direction=wind_direction;
    %output_data(time).height=Data(time).range(ring,1).*sind(Data(1).el(:,1));
end 



for time=1:T
    for ring=2:size(Data(time).range,1)
        if ~isnan(output_data(time).wind_speed(ring,1))==1 && ~isnan(output_data(time).wind_speed(ring-1,1))==1
            percent_change_increase=100*(output_data(time).wind_speed(ring,1)-output_data(time).wind_speed(ring-1,1))/output_data(time).wind_speed(ring-1,1);
            percent_change_decrease=100*(output_data(time).wind_speed(ring-1,1)-output_data(time).wind_speed(ring,1))/output_data(time).wind_speed(ring,1);
            
            if percent_change_increase>300 
                output_data(time).wind_speed(ring,1)=NaN;
            elseif percent_change_decrease>300
                output_data(time).wind_speed(ring-1,1)=NaN;
            end
        end
    end
end





%plot(sqrt(output_data(t).velocity(:,1).^2+output_data(t).velocity(:,2).^2),Data(t).range(:,1).*sind(Data(t).el(:,1)))

%hold on
wind_speed_mean=zeros(size(Data(1).range,1),1);
wind_direction_mean=zeros(size(Data(1).range,1),1);
notnancount=zeros(size(Data(time).range,1),1);

for ring=1:size(Data(time).range,1)
    %u_mean(size(Data(time).range,1))=0;
    %v_mean(size(Data(time).range,1))=0;
    
    
    for time=1:T
    %u_mean(ring)=u_mean(ring)+output_data(time).velocity(ring,1);
    %v_mean(ring)=v_mean(ring)+output_data(time).velocity(ring,2);
    
        if ~isnan(output_data(time).wind_speed(ring,1))==1
            wind_speed_mean(ring)=wind_speed_mean(ring)+output_data(time).wind_speed(ring,1);
            wind_direction_mean(ring)=wind_direction_mean(ring)+output_data(time).wind_direction(ring,1);
            notnancount(ring)=notnancount(ring)+1;        
        end
        
    end
    if notnancount(ring)==0 
        wind_speed_mean(ring)=NaN;
        wind_direction_mean(ring)=NaN;
    end
end


h = figure;

ax = gca;
ax.NextPlot = 'replaceChildren';

loops = T;  %time
M(loops) = struct('cdata',[],'colormap',[]);

h.Visible = 'off';

v = VideoWriter('movie.avi');
open(v);

for j = 1:loops
   
    plot(output_data(j).wind_speed(:),output_data(j).height(:));
    drawnow
    M(j) = getframe;
    writeVideo(v,M(j));
end

h.Visible = 'on';
%movie(M,1,3);
close(v);





%plot(output_data(20).wind_speed(:),Data(50).range(:,1).*sind(Data(50).el(:,1)))

% plot(wind_speed_mean./notnancount,Data(1).range(:,1).*sind(Data(1).el(:,1)));
% figure(2)
% plot(wind_direction_mean./notnancount,Data(1).range(:,1).*sind(Data(1).el(:,1)));
