x1 = Motion_1(:,2);
y1 = Motion_1(:,3);

x2 = Motion_2(:,2);
y2 = Motion_2(:,3);

x3 = Motion_3(:,2);
y3 = Motion_3(:,3);

x4 = Motion_4(:,2);
y4 = Motion_4(:,3);


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
