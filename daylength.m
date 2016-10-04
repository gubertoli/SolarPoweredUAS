function [mind,maxd]=daylength(L)

for J=1:365
    %CBM Model
    theta=0.2163108+2*atan(0.9671396*tan(0.00860*(J-186))); %revolution angle from day of the year
    P = asin(0.39795*cos(theta)); %sun declination angle 
    %daylength (plus twilight)- 
    p=0.8333; %sunrise/sunset is when the top of the sun is apparently even with horizon
    D(J) = 24 - (24/pi) * acos((sin(p*pi/180)+sin(L*pi/180)*sin(P))/(cos(L*pi/180)*cos(P))); 
end
    avg=mean(D);
    mind=min(D);
    maxd=max(D);
    plot(D,'LineWidth',2)
    hold on
    hline=refline([0 avg]);
    imin = min(D);
    imax = max(D);
    xlim([1 365])
    ylim([imin imax])
    xlabel('Day of Year')
    ylabel('Hours')
    grid on
    title('Day Length','FontWeight','Bold')
    text(5,avg+0.5,strcat('   Minimal Day Length:     ',num2str(mind), ' h'),'FontWeight','Bold','FontSize',8)

end

