1;
pkg load control



function draw_pulley(y)
  
  pd = 0.4;                   ## Pulley Diameter
  p_y = 0.5;                 ## Pulley position wrt y
  L = 1;                      ## Length of string
  ml = 0.2;                  ## Mass Length
  mb = 0.1;                  ## Mass Breadth
  x1 = y(1);
  x2 = L-y(1);
  hold on;
  clf;
  axis equal;
  rectangle('Position',[0-(pd/2),p_y-(pd/2),pd,pd],'Curvature',1,'FaceColor',[1 0 0]);
  rectangle('Position',[-(pd/2)-(ml/2),p_y-x1-(mb/2),ml,mb],'Curvature',0.1,'FaceColor',[0 0 1]);
  rectangle('Position',[(pd/2)-(ml/2),p_y-x2-(mb/2),ml,mb],'Curvature',0.1,'FaceColor',[0 1 0]);
  line ([0-(pd/2) 0-(pd/2)], [p_y p_y-x1], "linestyle", "-", "color", "k");
  line ([(pd/2) (pd/2)], [p_y p_y-x2], "linestyle", "-", "color", "k");
  line ([0 0], [1 p_y], "linestyle", "-", "color", "k");
  xlim([-1 1])
  ylim([-1 1])
  drawnow
  hold off
  
endfunction


function dy = pulley_dynamics(y, m1, m2, g, r, u)
  
  dy(1,1) = y(2);
  dy(2,1) =((m1 - m2)*g + u/r )/(m1 + m2);  
endfunction


function [t,y] = sim_pulley(m1, m2, g, r, y0)
  tspan = 0:0.1:10;                  ## Initialise time step           
  u = 0;                             ## No Input
  [t,y] = ode45(@(t,y)pulley_dynamics(y, m1, m2, g, r, u),tspan,y0);  
endfunction


function [A,B] = pulley_AB_matrix(m1, m2, g, r)
  A = [0 1; 0 0] ;
  B = [0 ; 1/r*(m1+m2)];
endfunction


function [t,y] = pole_place_pulley(m1, m2, g, r, y_setpoint, y0)
  [A,B] = pulley_AB_matrix(m1, m2, g, r);
  eigs = [-5, -2];
  m = m1+m2;
  K = [eigs(1,1)*eigs(1,2)*m, abs(sum(eigs)*m)];
  
  tspan = 0:0.1:10;                  ## Initialise time step 
  [t,y] =ode45(@(t,y)pulley_dynamics(y, m1, m2, g, r, -K*(y-y_setpoint)),tspan,y0); 
endfunction


function [t,y] = lqr_pulley(m1, m2, g, r, y_setpoint, y0)
  [A, B] = pulley_AB_matrix(m1, m2, g, r);
  Q = [1 0; 0 1];                   ## Initialise Q matrix
  R = 0.01;                   ## Initialise R 
  
  K = lqr(A, B, Q, R); 
  tspan = 0:0.1:10;                  ## Initialise time step 
  [t,y] = ode45(@(t,y)pulley_dynamics(y, m1, m2, g, r, -K*(y-y_setpoint)),tspan,y0); ;  
endfunction


function simple_pulley_main()
  m1 = 7.5;
  m2 = 7.6;
  g = 9.8;
  r = 0.2;
  y0 = [0.8 ; 0];                   ## Initial condtion
  y_setpoint = [0.4; 0];              ## Set Point
  
#  [t,y] = sim_pulley(m1, m2, g, r, y0);
#  [t,y] = pole_place_pulley(m1, m2, g, r, y_setpoint, y0);
  [t,y] = lqr_pulley(m1, m2, g, r, y_setpoint, y0);
  
  for k = 1:length(t)
    draw_pulley(y(k, :));  
  endfor
endfunction
