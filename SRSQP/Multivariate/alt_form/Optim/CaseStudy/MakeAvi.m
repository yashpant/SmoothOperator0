myavi = avifile('avidemo.avi','FPS',5);
%% describe sets and their signed dist approxs


%%



for j = 1:3
    if(j==1)
    close all;
    figure(1);
    drawsets;
    title('Trajectory 1','FontSize',20);
    load('Case_20_len_e.mat');
    [pos_1,pos_2] = getPos(x,optParams);
    pause;
    end
    
    if(j==2)
    close all;
    figure(1);
    drawsets;
    title('Trajectory 2','FontSize',20);
    load('Case_20_len_maxz_e.mat');
    [pos_1,pos_2] = getPos(x,optParams);
    pause;
    end
    
    if(j==3)
    close all;
    figure(1);
    drawsets;
    title('Trajectory 3','FontSize',20);
    load('Case_20_len_miny_e.mat');
    [pos_1,pos_2] = getPos(x,optParams);
    pause;
    end
    
 for i = 1:size(pos_1,2)
     hold on
     plot3(pos_1(1,i),pos_1(2,i),pos_1(3,i),'b*','MarkerSize',10);
     hold on;
     plot3(pos_2(1,i),pos_2(2,i),pos_2(3,i),'ko','MarkerSize',10);
     pause(0.2);
     %drawnow;
     myavi = addframe(myavi,gca);
 end
 
 for j = 1:3
     pause;
 for i = 1:5
 myavi = addframe(myavi,gca);
 end
 end
end
 myavi = close(myavi);
% figure(1)
% hold all;
% plot3(pos_1(1,:),pos_1(2,:),pos_1(3,:),'*','linewidth',2,'MarkerSize',8)
% hold all;
% plot3(pos_2(1,:),pos_2(2,:),pos_2(3,:),'o','linewidth',2,'MarkerSize',8)
