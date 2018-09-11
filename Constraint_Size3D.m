% Constraints for Size3D
function [F1,F2,F3] = Constraint_Size3D(type,X,Y,Z,N,eps)

switch type 
    
    case 1
    % 2 parallel circles
    angle_xy = angle(X + 1i*Y);
    F1 =  - 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2 + abs(Y).^2) - 0.7).^2 + (Z - 0.25).^2  )/eps^2)).*sin(angle_xy)  ...
          + 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2 + abs(Y).^2) - 0.7).^2 + (Z + 0.25).^2  )/eps^2)).*sin(angle_xy) ;
    F2 =  + 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2 + abs(Y).^2) - 0.7).^2 + (Z - 0.25).^2  )/eps^2)).*cos(angle_xy) ...
          - 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2 + abs(Y).^2) - 0.7).^2 + (Z + 0.25).^2  )/eps^2)).*cos(angle_xy) ;
    F3 = zeros(N,N,N);

    case 2
    % 2 orthogonal circles
    angle_xy = angle(X + 1i*(Y+.2));
    angle_yz = angle(Y-.2 + 1i*Z);
    F1 =  - 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2 + abs(Y+.2).^2) - 0.5).^2 + Z.^2  )/eps^2)).*sin(angle_xy); 
    F2 =  + 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2 + abs(Y+.2).^2) - 0.5).^2 + Z.^2  )/eps^2)).*cos(angle_xy) ...
          + 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(Z).^2 + abs(Y-.2).^2) - 0.5).^2 + X.^2  )/eps^2)).*sin(angle_yz) ;
    F3 =  - 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(Z).^2 + abs(Y-.2).^2) - 0.5).^2 + X.^2  )/eps^2)).*cos(angle_yz) ;

    case 3
    % 3 circles - trousers
    angle_xy_1 = angle(X + 1i*(Y-0.5));
    angle_xy_2 = angle(X + 1i*(Y+0.5));
    angle_xy_3 = angle(X + 1i*(Y));
    F1 = + (2*exp(-1*pi*((  sqrt(abs(X).^2  + abs(Y - 0.5).^2) - 0.3).^2 + (Z - 0.15).^2  )/eps^2)).*sin(angle_xy_1)/(sqrt(eps))...
         +  (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(Y + 0.5).^2) - 0.3).^2 + (Z - 0.15).^2  )/eps^2)).*sin(angle_xy_2)/(sqrt(eps)) ...
         - (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(Y).^2) - 0.6).^2 + (Z + 0.15).^2  )/eps^2)).*sin(angle_xy_3)/(sqrt(eps));
    F2 = - (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(Y - 0.5).^2) - 0.3).^2 + (Z - 0.15).^2  )/eps^2)).*cos(angle_xy_1)/(sqrt(eps))...
         - (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(Y + 0.5).^2) - 0.3).^2 + (Z - 0.15).^2  )/eps^2)).*cos(angle_xy_2)/(sqrt(eps)) ...
         + (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(Y).^2 )- 0.6).^2 + (Z + 0.15).^2  )/eps^2)).*cos(angle_xy_3)/(sqrt(eps));
    F3 = zeros(N,N,N);
    
    case 4
    % 3 parallel circles one big in the middle 
    angle_xy_1 = angle(X + 1i*(Y-0.25));
    angle_xy_2 = angle(X + 1i*(Y+0.25));
    angle_xy_3 = angle(X + 1i*(Y));
    F1 = + (2*exp(-1*pi*((  sqrt(abs(X).^2  + abs(Y - 0.25).^2) - 0.3).^2 + (Z - 0.3).^2  )/eps^2)).*sin(angle_xy_1)/(sqrt(eps))...
         + (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(Y + 0.25).^2) - 0.3).^2 + (Z + 0.3).^2  )/eps^2)).*sin(angle_xy_2)/(sqrt(eps)) ...
         - (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(Y).^2) - 0.6).^2 + (Z).^2  )/eps^2)).*sin(angle_xy_3)/(sqrt(eps));
    F2 = - (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(Y - 0.25).^2) - 0.3).^2 + (Z - 0.3).^2  )/eps^2)).*cos(angle_xy_1)/(sqrt(eps))...
         - (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(Y + 0.25).^2) - 0.3).^2 + (Z + 0.3).^2  )/eps^2)).*cos(angle_xy_2)/(sqrt(eps)) ...
         + (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(Y).^2 )- 0.6).^2 + (Z ).^2  )/eps^2)).*cos(angle_xy_3)/(sqrt(eps));
    F3 = zeros(N,N,N);
    
    case 5
    % 2 parallel and 2 orthogonal circles
    % 2 orthogonal circles
    angle_xy = angle(X + 1i*(Y+.2));
    angle_yz = angle(Y-.2 + 1i*Z);
    
    F1 = - 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2 + abs(Y+.2).^2) - 0.5).^2 + (Z-.2).^2  )/eps^2)).*sin(angle_xy) ...
         + 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2 + abs(Y+.2).^2) - 0.5).^2 + (Z+.2).^2  )/eps^2)).*sin(angle_xy);
    
    F2 = + 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2 + abs(Y+.2).^2) - 0.5).^2 + (Z-.2).^2  )/eps^2)).*cos(angle_xy) ...
         - 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2 + abs(Y+.2).^2) - 0.5).^2 + (Z+.2).^2  )/eps^2)).*cos(angle_xy) ...
         ...
         + 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(Z).^2 + abs(Y-.2).^2) - 0.45).^2 + (X-.25).^2  )/eps^2)).*sin(angle_yz) ...
         - 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(Z).^2 + abs(Y-.2).^2) - 0.45).^2 + (X+.25).^2  )/eps^2)).*sin(angle_yz) ;
    
    F3 = + 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(Z).^2 + abs(Y-.2).^2) - 0.45).^2 + (X-.25).^2  )/eps^2)).*cos(angle_yz) ...
         - 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(Z).^2 + abs(Y-.2).^2) - 0.45).^2 + (X+.25).^2  )/eps^2)).*cos(angle_yz) ;
         
    case 6
    % swigly circle
    nG = 3; h = .3;
    angle_xy = angle(X + 1i*Y);
    F1 =  - 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2 + abs(Y).^2) - 0.7).^2 + (Z-h*sin(nG*angle_xy)).^2  )/eps^2)).*sin(angle_xy);
    F2 =  + 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2 + abs(Y).^2) - 0.7).^2 + (Z-h*sin(nG*angle_xy)).^2  )/eps^2)).*cos(angle_xy);
    F3 =  + h/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2 + abs(Y).^2) - 0.7).^2 + (Z-h*sin(nG*angle_xy)).^2  )/eps^2)).*cos(nG*angle_xy);

    case 7
    % three circles
    r = .9;
    R = .6;
    angle_xy = angle(X + 1i*Y);
    angle_yz = angle(Y + 1i*Z);
    angle_xz = angle(X + 1i*Z);
    F1 =  - 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2/r^2 + abs(Y).^2/R^2) - 0.8).^2 + Z.^2  )/eps^2)).*sin(angle_xy) ...
          - 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2/R^2 + abs(Z).^2/r^2) - 0.8).^2 + Y.^2  )/eps^2)).*sin(angle_xz);
    F2 =  + 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2/r^2 + abs(Y).^2/R^2) - 0.8).^2 + Z.^2  )/eps^2)).*cos(angle_xy) ...
          - 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(Y).^2/r^2 + abs(Z).^2/R^2) - 0.8).^2 + X.^2  )/eps^2)).*sin(angle_yz);
    F3 =  + 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(Y).^2/r^2 + abs(Z).^2/R^2) - 0.8).^2 + X.^2  )/eps^2)).*cos(angle_yz) ...
          + 1/(sqrt(eps))*( 2*exp(-1*pi*(( sqrt(abs(X).^2/R^2 + abs(Z).^2/r^2) - 0.8).^2 + Y.^2  )/eps^2)).*cos(angle_xz);
      
    case 8
    theta = pi/6;
    angle_xy_1 = angle(X + 1i*(Y-0.5));
    angle_xy_2 = angle(X + 1i*(Y+0.5));
    angle_xy_3 = angle(X + 1i*(Y));
    F1 = + (2*exp(-1*pi*((  sqrt(abs(X).^2  + abs(cos(-theta)*Y+sin(-theta)*(Z-.15) - 0.5).^2) - 0.3).^2 + (-sin(-theta)*Y+cos(-theta)*(Z-.15) - 0.2).^2  )/eps^2)).*sin(angle_xy_1)/(sqrt(eps))...
         + (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(cos(theta)*Y+sin(theta)*(Z-.15) + 0.5).^2) - 0.3).^2 + (-sin(theta)*Y+cos(theta)*(Z-.15) - 0.2).^2  )/eps^2)).*sin(angle_xy_2)/(sqrt(eps)) ...
         - (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(Y).^2) - 0.6).^2 + (Z + 0.15).^2  )/eps^2)).*sin(angle_xy_3)/(sqrt(eps));

    F2 = - (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(cos(-theta)*Y+sin(-theta)*(Z-.15) - 0.5).^2) - 0.3).^2 + (-sin(-theta)*Y+cos(-theta)*(Z-.15)- 0.2).^2  )/eps^2)).*cos(angle_xy_1)*cos(-theta)/(sqrt(eps))...
         - (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(cos(theta)*Y+sin(theta)*(Z-.15) + 0.5).^2) - 0.3).^2 + (-sin(theta)*Y+cos(theta)*(Z-.15) - 0.2).^2  )/eps^2)).*cos(angle_xy_2)*cos(theta)/(sqrt(eps)) ...
         + (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(Y).^2 )- 0.6).^2 + (Z + 0.15).^2  )/eps^2)).*cos(angle_xy_3)/(sqrt(eps));

    F3 = + (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(cos(-theta)*Y+sin(-theta)*(Z-.15) - 0.5).^2) - 0.3).^2 + (-sin(-theta)*Y+cos(-theta)*(Z-.15) - 0.2).^2  )/eps^2)).*cos(angle_xy_1)*sin(-theta)/(sqrt(eps))...
         + (2*exp(-1*pi*((  sqrt(abs(X).^2 + abs(cos(theta)*Y+sin(theta)*(Z-.15) + 0.5).^2) - 0.3).^2 + (-sin(theta)*Y+cos(theta)*(Z-.15) - 0.2).^2  )/eps^2)).*cos(angle_xy_2)*sin(theta)/(sqrt(eps));


end



end