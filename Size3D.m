% Numerical implementation of the functional in the article Sigma3D 
% Chambolle - Ferrari - Merlet
% I am using the methods proposed by Elie Bretin with the Fourier Transform
% with gradient flow and projection on the constraint.
clear all; close all;

movie = false;

% Domain
k0 = 6; % Resolution
N = 2^k0;
x = linspace(-1,1,N);
[X,Y,Z] = meshgrid(x,x,x);

eps = 3/N;
dt = eps^2; % dt =eps^3;

% Constraint Initialization
[F1,F2,F3] = Constraint_Size3D(8,X,Y,Z,N,eps);


weight = 100;
F1 = weight/(sqrt(eps))*F1;
F2 = weight/(sqrt(eps))*F2;
F3 = weight/(sqrt(eps))*F3;

% Movie visualization
clf

hold on;
v = sqrt(F1.^2 + F2.^2 + F3.^2) ;
p = patch(isosurface(x,x,x,v,10,v));
isonormals(x,x,x,v,p)
set(p,'FaceColor','red','EdgeColor','none');

view(120,10);
axis([-1,1,-1,1,-1,1]);
camlight
lighting gouraud

if movie
	set(gca,'nextplot','replacechildren'); 
	FileName=['Size3D',datestr(now, 'HH-MM-dd-mmm-yyyy'),'.avi'];
	mov= VideoWriter(FileName);
	open(mov);
end

% Fourier Space
k = 1/2*[0:1:N/2,-N/2+1:1:-1];
[K1,K2,K3] = meshgrid(k,k,k);

% Symbols in Fourier space
Delta = -4*pi^2*(K1.^2 + K2.^2 + K3.^2);
Delta(1,1)=-1;
Delta_F = (4*pi^2*(K1.^2 + K2.^2 + K3.^2));

% Variables initialization
phi = ones(N,N,N);
F1_Fourier = fftn(F1);
F2_Fourier = fftn(F2);
F3_Fourier = fftn(F3);
rot1_F_fourier = -2*1i*pi*(K2.*F3_Fourier -  K3.*F2_Fourier );
rot2_F_fourier = -2*1i*pi*(K3.*F1_Fourier -  K1.*F3_Fourier );
rot3_F_fourier = -2*1i*pi*(K1.*F2_Fourier -  K2.*F1_Fourier );

sigma_1 =  real(ifftn(2*1i*pi*(K2.*rot3_F_fourier -  K3.*rot2_F_fourier)./Delta));
sigma_2 =  real(ifftn(2*1i*pi*(K3.*rot1_F_fourier -  K1.*rot3_F_fourier)./Delta));
sigma_3 =  real(ifftn(2*1i*pi*(K1.*rot2_F_fourier -  K2.*rot1_F_fourier)./Delta));

 
err = 1;

% Solver
for j = 1 : 100000
    	
	if 0 %err < 1e-4
		eps = 3/N;
		dt = eps^2; 	
	end

	phi_Old = phi;
	norm_sigma = sqrt(sigma_1.^2 + sigma_2.^2 +  sigma_3.^2);
	beta = max(norm_sigma(:));

	% Optimization step in phi
	phi = + real(ifftn(fftn(phi + dt*(1/eps^2 - ((norm_sigma - beta)).*phi/eps^2))./(1 + dt*(eps^0*Delta_F + 1/eps^2 + beta/eps^2))));

	% Optimization step in sigma
	sigma_1 = sigma_1./(1 + dt/eps^2*(phi.^2 ));
	sigma_2 = sigma_2./(1 + dt/eps^2*(phi.^2 ));
	sigma_3 = sigma_3./(1 + dt/eps^2*(phi.^2 ));

	% Diffusion of sigma 
	sigma_1 = real(ifftn(exp(-dt*0.1*eps^2*Delta_F).*fftn(sigma_1)));
	sigma_2 = real(ifftn(exp(-dt*0.1*eps^2*Delta_F).*fftn(sigma_2)));
	sigma_3 = real(ifftn(exp(-dt*0.1*eps^2*Delta_F).*fftn(sigma_3)));

	% Projection onto the curl constraint
	sigma1_Fourier = fftn(sigma_1);
	sigma2_Fourier = fftn(sigma_2);
	sigma3_Fourier = fftn(sigma_3);

	rot1_sigma_fourier = -2*1i*pi*(K2.*sigma3_Fourier -  K3.*sigma2_Fourier );
	rot2_sigma_fourier = -2*1i*pi*(K3.*sigma1_Fourier -  K1.*sigma3_Fourier );
	rot3_sigma_fourier = -2*1i*pi*(K1.*sigma2_Fourier -  K2.*sigma1_Fourier );

	lambda1_fourier = (F1_Fourier  - rot1_sigma_fourier)./(Delta);
	lambda2_fourier = (F2_Fourier  - rot2_sigma_fourier)./(Delta);
	lambda3_fourier = (F3_Fourier  - rot3_sigma_fourier)./(Delta);

	sigma_1 = sigma_1 +  real(ifftn(2*1i*pi*(K2.*lambda3_fourier -  K3.*lambda2_fourier)));
	sigma_2 = sigma_2 + real(ifftn(2*1i*pi*(K3.*lambda1_fourier -  K1.*lambda3_fourier)));
	sigma_3 = sigma_3 + real(ifftn(2*1i*pi*(K1.*lambda2_fourier -  K2.*lambda1_fourier)));
	
	err = norm(phi(:)-phi_Old(:),1);

	% Visualization
	if (mod(j,50)==1)
		disp(j)  
 		clf

		hold on;
		v = sqrt(F1.^2 + F2.^2 + F3.^2) ;
		p = patch(isosurface(x,x,x,v,10,v));
		isonormals(x,x,x,v,p);title(['Error = ', num2str(err)]);
		set(p,'FaceColor','red','EdgeColor','none');

		v=real(1-phi);
		p = patch(isosurface(x,x,x,v,0.9,v));
		isonormals(x,x,x,v,p)
		set(p,'FaceColor','green','EdgeColor','none');

		view(120,10);
		axis([-1,1,-1,1,-1,1]);


		camlight
		lighting gouraud

		if movie
			frame = getframe(gcf);
			writeVideo(mov,frame);
		else pause(.01)
		end
	 
	end


end

FileName=['Steiner2D',datestr(now, 'HH-MM-dd-mmm-yyyy'),'.jpg'];
imwrite(frame2im(frame),FileName);

if movie
	close(mov)
end
