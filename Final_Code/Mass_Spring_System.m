1;
pkg load control



function draw_mass_spring(y)    
  l = 0.3; ## Length of rectangle
  b = 0.2; ## Breadth of rectangle
  
  x_pos = y(1);
  y_pos = 0;
  
  hold on;
  clf;
  axis equal;
  rectangle('Position',[x_pos-l/2,y_pos,l,b],'Curvature',0,'FaceColor',[0 0 1]);
  line ([0 0], [0 0.6], "linestyle", "-", "color", "k");
  line ([-0.5 x_pos], [(y_pos+b)/2 (y_pos+b)/2], "linestyle", "--", "color", "k");
  text(-0.05, 0.65, "Eqbm Pt")
  xlim([-0.5 1])
  ylim([0 1])
  drawnow
  hold off;
endfunction


function dy = mass_spring_dynamics(y, m, k, u)
  
  dy(1,1) = y(2);
  dy(2,1) = (u/m) + (-k/m)*y(1) ;
endfunction


function [t,y] = sim_mass_spring(m, k, y0)
  tspan = 0:0.1:10;                ## Initialize time step
  u = 0;                           ## No input
  [t,y] = ode45(@(t,y)mass_spring_dynamics(y, m, k, u),tspan,y0);;
endfunction


function [A,B] = mass_spring_AB_matrix(m, k)
  A = [0 1 ; (-k/m) 0];
  B = [0 ;(1/m)];
endfunction


function [t,y] = pole_place_mass_spring(m, k, y_setpoint, y0)
  [A,B] = mass_spring_AB_matrix(m, k);  ## Initialize A and B matrix 
  eigs =[-35 ; -50] ;                    ## Initialise desired eigenvalues
  K = [(-k)+(m*eigs(1)*eigs(2)) , -(m*eigs(1) + m*eigs(2))];                ## Calculate K matrix for desired eigenvalues
  
  tspan = 0:0.1:10;                   ## Initialise time step 
  [t,y] = ode45(@(t,y)mass_spring_dynamics(y, m, k, -K*(y-y_setpoint)),tspan,y0);
endfunction


function [t,y] = lqr_mass_spring(m, k, y_setpoint, y0)
  [A,B] = mass_spring_AB_matrix(m, k);  ## Initialize A and B matrix 
  Q =[1000 0; 0 1] ;                  ## Initialise desired eigenvalues
  R = 1;               ## Calculate K matrix for desired eigenvalues
  
  K = lqr(A,B,Q,R);
  
  tspan = 0:0.1:10;                   ## Initialise time step 
  [t,y] = ode45(@(t,y)mass_spring_dynamics(y, m, k, -K*(y-y_setpoint)),tspan,y0);
endfunction


function mass_spring_main()
  m = 0.2;
  k = 0.8;
  y0 = [-0.3; 0];
  y_setpoint = [0.7; 0];
  
#  [t,y] = sim_mass_spring(m,k, y0);      
#  [t,y] = pole_place_mass_spring(m, k, y_setpoint, y0); 
  [t,y] = lqr_mass_spring(m, k, y_setpoint, y0);  
  for k = 1:length(t)
    draw_mass_spring(y(k, :));  
  endfor
endfunction
