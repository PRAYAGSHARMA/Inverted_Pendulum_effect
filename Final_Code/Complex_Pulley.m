1;
pkg load control


function draw_complex_pulley(y)
  ml = 0.2;
  mb = 0.1;
  L_A = 1.3;
  L_B = 1.1;
  
  pd_A = 0.4;
  py_A = 1;
  pd_B = 0.3; 
  
  x1 = y(1);
  x2 = L_A - y(1);
  y1 = y(3);
  y2 = L_B - y(3);
  
  pulley_A_pos = {0, py_A};
  pulley_B_pos = {(-pd_A/2), py_A-x2};
  m1_pos = {(pd_A/2), py_A-x1};
  m2_pos = {((-pd_A+pd_B)/2), py_A-x2-y1};
  m3_pos = {((-pd_A-pd_B)/2), py_A-x2-y2};
  x1_string = {(pd_A/2), py_A, (pd_A/2), (py_A-x1)};
  x2_string = {(-pd_A/2), py_A, (-pd_A/2), (py_A-x2)};
  y1_string = {((-pd_A+pd_B)/2), py_A-x2, ((-pd_A+pd_B)/2), py_A-x2-y1};
  y2_string = {((-pd_A-pd_B)/2), py_A-x2, ((-pd_A-pd_B)/2), py_A-x2-y2};
  
  hold on;
  clf;
  axis equal;
  rectangle('Position',[pulley_A_pos{1}-(pd_A/2),pulley_A_pos{2}-(pd_A/2),pd_A, pd_A],'Curvature',1,'FaceColor',[1 0 0]);## Pulley A
  rectangle('Position',[pulley_B_pos{1}-(pd_B/2),pulley_B_pos{2}-(pd_B/2),pd_B, pd_B],'Curvature',1,'FaceColor',[1 0 0]);## Pulley B
  rectangle('Position',[m1_pos{1}-(ml/2),m1_pos{2}-(mb/2),ml, mb],'Curvature',0.1,'FaceColor',[0 0 1]);## m1 mass
  rectangle('Position',[m2_pos{1}-(ml/2),m2_pos{2}-(mb/2),ml, mb],'Curvature',0.1,'FaceColor',[0 1 0]);## m2 mass
  rectangle('Position',[m3_pos{1}-(ml/2),m3_pos{2}-(mb/2),ml, mb],'Curvature',0.1,'FaceColor',[0 1 1]);## m1 mass
  line ([x1_string{1} x1_string{3}], [x1_string{2} x1_string{4}], "linestyle", "-", "color", "k");
  line ([x2_string{1} x2_string{3}], [x2_string{2} x2_string{4}], "linestyle", "-", "color", "k");
  line ([y2_string{1} y2_string{3}], [y2_string{2} y2_string{4}], "linestyle", "-", "color", "k");
  line ([y1_string{1} y1_string{3}], [y1_string{2} y1_string{4}], "linestyle", "-", "color", "k");
  xlim([-1.5 1.5]);
  ylim([-1.5 1.5]);
  drawnow
  hold off
  
endfunction


function dy = complex_pulley_dynamics(y, m1, m2, m3, g, rA, rB, u)
  dy(1,1) = y(2);
  dy(2,1) = ((m1-m2-m3)*g + u(1)/rA -(u(2)/rB)*(m3-m2)/(m3+m2) + (m2-m3)*(m2-m3)*g/(m3+m2))/(m1+m2+m3-((m3-m2)*(m3-m2)/(m3+m2)));
  dy(3,1) = y(4);
  dy(4,1) = ((m2-m3)*(g+dy(2,1)) + u(2))/(m3+m2);
  
endfunction


function [t,y] = sim_complex_pulley(m1, m2, m3, g, rA, rB, y0)
  tspan = 0:0.1:10;                  ## Initialise time step           
  u = [0; 0];                             ## No Input
  [t,y] = ode45(@(t,y)complex_pulley_dynamics(y, m1, m2, m3, g, rA, rB, u),tspan,y0);
  endfunction


function [A,B] = complex_pulley_AB_matrix(m1, m2, m3, g, rA, rB)
  M1 = m1 + m2 + m3;
  M2 = m3 - m2;
  M3 = m3 + m2;
  M4 = m3 + m2 - m1;
  X = (M1 - M2*M2/M3);
  
  A = [0 1 0 0;0 0 0 0;0 0 0 1;0 0 0 0];
  B = [0 0;
       1/(rA*X) -M2/(rB*M3*X);
       0 0;
       -M2/(rA*X*M3) (M2*M2/(rB*M3*X*M3)+1/(rB*M3))];
  

endfunction


function [t,y] = pole_place_complex_pulley(m1, m2, m3, g, rA, rB, y_setpoint, y0)
  [A,B] = complex_pulley_AB_matrix(m1, m2, m3, g, rA, rB);
  eigs = [-10;-10;-10;-10];
  K = place(A,B,eigs);
  tspan = 0:0.1:10;                  ## Initialise time step 
  [t,y] = ode45(@(t,y)complex_pulley_dynamics(y, m1, m2, m3, g, rA, rB, -K*(y-y_setpoint)),tspan,y0);
endfunction


function [t,y] = lqr_complex_pulley(m1, m2, m3, g, rA, rB, y_setpoint, y0)
  [A,B] = complex_pulley_AB_matrix(m1, m2, m3, g, rA, rB);
  Q = [10000 0 0 0;0 1 0 0;0 0 100000 0;0 0 0 1];
  R = [1 0; 0 1]
  K = lqr(A, B, Q, R);
  tspan = 0:0.1:10;                  ## Initialise time step 
  [t,y] = ode45(@(t,y)complex_pulley_dynamics(y, m1, m2, m3, g, rA, rB, -K*(y-y_setpoint)),tspan,y0);
endfunction


function complex_pulley_main()
  m1 = 23.95;
  m2 = 12.0;
  m3 = 12.3;
  g = 9.8;
  rA = 0.2;
  rB = 0.2;
  y_setpoint = [0.5 ; 0; 0.8; 0];
  y0 = [0.4 ; 0; 0.5; 0];
  
#  [t,y] = sim_complex_pulley(m1, m2, m3, g, rA, rB, y0);
  [t,y] = pole_place_complex_pulley(m1, m2, m3, g, rA, rB, y_setpoint, y0);
##  [t,y] = lqr_complex_pulley(m1, m2, m3, g, rA, rB, y_setpoint, y0);
  
  for k = 1:length(t)
    draw_complex_pulley(y(k, :));
  endfor
  
endfunction