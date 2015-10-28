

G = 1;       % Gravitational Constant
M1 = 1;      % mass 1
M2 = 0.1;    % mass 2
M = M1 + M2; % total mass

R = 1;       % distance between M1 and M2 set to 1


%%%%%%%%%%%%%%%%  Orbital properties of two massive bodies  %%%%%%%%%%%%%%%
mu = G*M;
mu1 = G*M1;
mu2 = G*M2;

R10 = [-M2/M; 0];           % initial position of M1
R20 = [M1/M; 0];            % initial position of M2

P = 2*pi * sqrt(R^3 / mu);  % period from Kepler's 3rd law
omega0 = 2*pi/P;            % angular velocity of massive bodies


[X,Y] = meshgrid(-2:0.01:2);
U = crtbpPotential(mu1, mu2, X, Y);


U0 = crtbpPotential(mu1, mu2, 0, 1);


m = 2;
map = ones(m , 3)*.8;

colormap(map);



level = [-2 -1.96 -1.9 -1.7 -1.68 -1.64];

for j = 1:6
    
    subplot(2,3,j)
    
    %  contourf(X,Y,U,[-5:.5:-1]);
    contourf(X,Y,U,level(j)+[0 -0.0001]);
    hold on
    plot(-M2/(M1+M2),0,'ko','MarkerSize',7,'MarkerFaceColor','y')
    plot(M1/(M1+M2),0,'ko','MarkerSize',5,'MarkerFaceColor','g')

    % xlim([min(X) max(X)])
    % ylim([min(Y) max(Y)])
    title(sprintf('C_J = %.2f',-2*level(j)));
    axis equal
    axis square
    
    hold on
    

    LP = lagrangePoints(M2/(M1+M2));
    plot(LP(:,1),LP(:,2),'k+')
    
end



set(gcf, 'PaperPosition', [0 -1 8 4]);
print ('-dpng','zeroVelocity.png','-r100');