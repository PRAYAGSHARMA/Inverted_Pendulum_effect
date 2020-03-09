1;
pkg load control


function draw_cart_pendulum(y,m, M, L)
  x = y(1);
  theta = y(3);
  
  W = 1*sqrt(M/5);    # cart width
  H = 0.5*sqrt(M/5);  # cart height 
  wr = 0.2;           # wheel radius
  mr = 0.3*sqrt(m);    # mass radius 
  
  y = wr/2 + H/2;
  w1x = x - 0.9*W/2;
  w1y = 0;
  w2x = x + 0.9*W/2 - wr;
  w2y = 0;
  
  px = x - L*sin(theta);
  py = y - L*( cos(theta));
   
  hold on;
  clf;
  axis equal;
  rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[1 0.1 0.1])
  rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[0 0 0])
  rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[0 0 0])
  
  line ([-10 10], [0 0], "linestyle", "-", "color", "k");
  line ([x px], [y py], "linestyle", "-", "color", "k");
  rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',[.1 0.1 1])
  
  xlim([-6 6]);
  ylim([-2 3]);
  set(gcf, 'Position', [200 300 1000 400]);
  drawnow
  hold off
endfunction


function dy = cart_pendulum_dynamics(y, m, M, L, g,  u)
  Cosy = cos(y(3));
  Siny = sin(y(3));
  
  dy(1,1) = y(2);
  dy(4,1) = (-(M+m)*(g)*(Siny/Cosy) + u - m*L*y(4)*y(4)*Siny)/(L*((m+M)/Cosy - m*Cosy));
  dy(3,1) = y(4);
  dy(2,1) = (L*dy(4,1) + (g)*Siny)/Cosy;  
endfunction

function [t,y] = sim_cart_pendulum(m, M, L, g, y0)
  tspan = 0:0.1:10;                  ## Initialise time step           
  u = 0;                             ## No Input
  [t,y] = ode45(@(t,y)cart_pendulum_dynamics(y, m, M, L, g,  u),tspan,y0); ## Solving the differential equation    
endfunction


function [A, B] = cart_pendulum_AB_matrix(m , M, L, g)
  A = [0 1 0 0;0 0 m*(g)/M 0;0 0 0 1;0 0 (m+M)*(g)/(M*L) 0];
       
  B = [0; 1/M; 0; -1/(M*L)];  
endfunction


function [t,y] = pole_place_cart_pendulum(m, M, L, g, y_setpoint, y0)
  [A,B] = cart_pendulum_AB_matrix(m , M, L, g);
  eigs = [-10;-1;-2;-2];
  K = place(A,B,eigs);
  tspan = 0:0.1:10;
  [t,y] = ode45(@(t,y)cart_pendulum_dynamics(y, m, M, L, g,-K*(y-y_setpoint)),tspan,y0);
endfunction


function [t,y] = lqr_cart_pendulum(m, M, L, g, y_setpoint, y0)
  [A,B] = cart_pendulum_AB_matrix(m , M, L, g);
  Q = [200 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 10];
  R = 1;
  K = lqr(A,B,Q,R);
  tspan = 0:0.1:10;
  [t,y] = ode45(@(t,y)cart_pendulum_dynamics(y, m, M, L, g,-K*(y-y_setpoint)),tspan,y0);
endfunction


function cart_pendulum_main()
  m = 1;
  M = 5;
  L = 2;
  g = 9.8;
  y0 = [-4; 0; pi+0.9; 0];
  y_setpoint = [0; 0; pi; 0];
  

  [t,y] = lqr_cart_pendulum(m, M, L, g, y_setpoint, y0);
  
  for k = 1:length(t)
    draw_cart_pendulum(y(k, :), m, M, L);  
  endfor
  
endfunction

