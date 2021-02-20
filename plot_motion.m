x1 = Motion_1(:,2);
y1 = Motion_1(:,3);

x2 = Motion_2(:,2);
y2 = Motion_2(:,3);

x3 = Motion_3(:,2);
y3 = Motion_3(:,3);

x4 = Motion_4(:,2);
y4 = Motion_4(:,3);

x5 = Motion_5(:,2);
y5 = Motion_5(:,3);

% plot(x2, y2, 'g-');
% legend('GNSS-KF');
% xlabel('latitude (rad)');
% ylabel('longitude (rad)');
% title('position solution');

plot(x2, y2, 'g-');
hold on
plot(x3, y3, 'b-');
hold on
plot(x4,y4,'r-');
hold on
legend('GNSS-KF', 'DR', 'GNSS+DR');
xlabel('latitude (rad)');
ylabel('longitude (rad)');
title('position solution');

% x6 = Motion_1(:,1);
% y6 = Heading;
% y6_ori = DR(:,7);
% plot(x6,y6, 'r-');
% legend('Heading');
% xlabel('time (s)');
% ylabel('degree (rad)');
% title('Heading solution');

% plot(x4, y4, 'b-');
% hold on
% plot(x5,y5,'r-');
% hold on
% legend('open-loop', 'closed-loop');
% xlabel('latitude (rad)');
% ylabel('longitude (rad)');
% title('position solution');
% 
% dlmwrite('Motion_final.csv',Motion_4, 'precision', 8);
