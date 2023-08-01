%============KONSTANTINOS FYTILIS================
%---------------2018010190-----------------------
syms x L C S I x1 x2
b=9; c=1;


%===============================================================
%-----------IPE-200 i-BEAM--------------------------------------
density=22.4;
Area=28.5e-4;
I_ibeam=1943e-8;
%--------------------------------------------------------------- 
E=200e6;%evala 10^6 (eno einai GPa) gia na afiso tis dinams se kN
%I=200;

%----------shape Functions----------------
N1=1/L^3*(2*x^3-3*x^2*L+L^3);
N2=1/L^3*(L*x^3-2*x^2*L^2+x*L^3);
N3=1/L^3*(-2*x^3+3*x^2*L);
N4=1/L^3*(L*x^3-x^2*L^2);
N=[N1 N2 N3 N4];
%------------------------------------------
syms d_1x d_1y phi_1 d_2x d_2y phi_2;
d_beam=[d_1y; phi_1; d_2y; phi_2;];%---4 bathmoi eleutherias gia kathe doko---!!!!!!-----leipei i aksoniki eleutheria, mainei stin grammi 35
u=N*d_beam;%---to dianysma metatopisis

%==========Eksagogi sxeseon (gia pinaka dokoy)===================
f1y=subs(diff(u,x,3),x,0);
m1=-subs(diff(u,x,2),x,0);
f2y=-subs(diff(u,x,3),x,L);
m2=subs(diff(u,x,2),x,L);
F_loc=[f1y; m1; f2y; m2];

for i=1:size(d_beam,1)
    for j=1:size(d_beam,1)
        Temp=subs(F_loc(i),d_beam(j),1);
        Temp=subs(Temp,[d_beam(1:j-1); d_beam(j+1:end)], zeros(4-1,1));
        K_loc_y(i,j)=Temp;
    end
end
K_loc_y=I*K_loc_y;
%--------K_loc, einai o pinakas gia kathe komvo tis dokou------
%---------------------------------------------------------------
%--------Prosthiki kai tis aksonikis eleutherias---------------
%==============================================================
C1=Area/L;
K_loc=[C1 zeros(1,2) -C1 zeros(1,2)];
K_loc=[K_loc; zeros(2,1) K_loc_y(1:2,1:2) zeros(2,1) K_loc_y(1:2,3:4)];
K_loc=[K_loc; -C1 zeros(1,2) C1 zeros(1,2)];
K_loc=[K_loc; zeros(2,1) K_loc_y(3:4,1:2) zeros(2,1) K_loc_y(3:4,3:4)];  
%---------eisagogi dedomenon---------------------------
[n_s, n_e, LenX, LenY]=textread('frame.txt','%d %d %f %f');
[nodes, Fx_nodal, Fy_nodal, M_nodal, Sx, Sy, Sphi]=textread('frame_concentr_force.txt','%d %f %f %f %d %d %d');
[n_sw, n_ew, Fy_w1, Fy_w2]=textread('frame_dist_force.txt','%d %d %f %f');
%-----------------------------------------------------

F=[];
for i=1:size(nodes,1)
    F=[F; Fx_nodal(i); Fy_nodal(i); M_nodal(i)];
end
F
%================================================================
%----------------ADD GRAVITY-------------------------------------
%----------FOR VERTICAL BEAMS ONLY Y NODAL FORCES ---------------
%----------FOR INCLINED BEAMS CONSTANT DISTRUBED VERTICAL FORCE--
for i=1:size(n_s,1)
	Len=(LenX(i)^2+LenY(i)^2)^0.5;
	if (LenX(i)==0) 
	F(3*(n_s(i)-1)+2)=F(3*(n_s(i)-1)+2)-density*10*Len/2*1e-3; 
	%half gravity force in first nodes in -y direction
	F(3*(n_e(i)-1)+2)=F(3*(n_e(i)-1)+2)-density*10*Len/2*1e-3; 
	%half gravity force in second nodes in -y direction
   else %according to constant distrubed load
	F(3*(n_s(i)-1)+2)=F(3*(n_s(i)-1)+2)-density*10*Len/2*1e-3;
	F(3*(n_e(i)-1)+2)=F(3*(n_e(i)-1)+2)-density*10*Len/2*1e-3;		
	F(3*(n_s(i)-1)+3)=F(3*(n_s(i)-1)+3)-density*10*Len^2/12*1e-3;
	F(3*(n_e(i)-1)+3)=F(3*(n_e(i)-1)+3)+density*10*Len^2/12*1e-3;
	end
end
%---multiply with 1e-3, because force unit is kN-----------------
%================================================================

F

if ([Fy_w1 Fy_w2]==zeros(size(n_sw,1),2))
dist_forc_loop=0;
else 
dist_forc_loop=size(n_sw,1);
end

%------------Metasximatismos trigonika katanemimen dinameon se NODAL------
F_dist_forc=zeros(3*size(nodes,1),1);
for i=1:dist_forc_loop
Work=int((Fy_w1(i)+(Fy_w2(i)-Fy_w1(i))/L*x)*u,x,0,L);%--Equivalent Nodal Force from LINEAR DISTR Load Y axis
%============================================================================================================================
Work=subs(Work,L,(LenX(i)^2+LenY(i)^2)^0.5);
F_dist_forc(3*(n_sw(i)-1)+2)=F_dist_forc(3*(n_sw(i)-1)+2)+subs(Work,d_beam,[1; 0; 0; 0]);
F(3*(n_sw(i)-1)+2)=F(3*(n_sw(i)-1)+2)+subs(Work,d_beam,[1; 0; 0; 0]); %------------------------F1y
F_dist_forc(3*(n_ew(i)-1)+2)=F_dist_forc(3*(n_ew(i)-1)+2)+subs(Work,d_beam,[0; 0; 1; 0]);
F(3*(n_ew(i)-1)+2)=F(3*(n_ew(i)-1)+2)+subs(Work,d_beam,[0; 0; 1; 0]); %------------------------F2y
F_dist_forc(3*(n_sw(i)-1)+3)=F_dist_forc(3*(n_sw(i)-1)+3)+subs(Work,d_beam,[0; 1; 0; 0]);
F(3*(n_sw(i)-1)+3)=F(3*(n_sw(i)-1)+3)+subs(Work,d_beam,[0; 1; 0; 0]); %------------------------M1
F_dist_forc(3*(n_ew(i)-1)+3)=F_dist_forc(3*(n_ew(i)-1)+3)+subs(Work,d_beam,[0; 0; 0; 1]);
F(3*(n_ew(i)-1)+3)=F(3*(n_ew(i)-1)+3)+subs(Work,d_beam,[0; 0; 0; 1]); %------------------------M2
end


%-------anagnorisi eleytheron kombon----------
nodes_free=zeros(3*size(nodes,1),1);
nodes_fixed=zeros(3*size(nodes,1),1);
for i=1:size(nodes,1)
	nodes_fixed(3*(i-1)+1)=Sx(i)*(3*(i-1)+1);
    nodes_fixed(3*(i-1)+2)=Sy(i)*(3*(i-1)+2);
    nodes_fixed(3*(i-1)+3)=Sphi(i)*(3*(i-1)+3);
	nodes_free(3*(i-1)+1)=(1-Sx(i))*(3*(i-1)+1);
    nodes_free(3*(i-1)+2)=(1-Sy(i))*(3*(i-1)+2);
    nodes_free(3*(i-1)+3)=(1-Sphi(i))*(3*(i-1)+3);
end
%---------------------------------------------

%-----------Rotation Matrix(for translation to GLOBAL from LOCAL COORDINATES)------------------------------------------------------
T=[C S 0; -S C 0; 0 0 1];
T=[T zeros(3,3); zeros(3,3) T];
K_glob=T'*K_loc*T;
%-------Global Stiffeness Matrix--------------
K=zeros(3*size(nodes,1),3*size(nodes,1));
for i=1:size(n_s,1)
	Len=(LenX(i)^2+LenY(i)^2)^0.5;
    K_glob_num=E*subs(K_glob,[C S L I], [LenX(i)/Len LenY(i)/Len Len I_ibeam]);
    K_glob_num=double(K_glob_num);
    K(3*(n_s(i)-1)+1:3*(n_s(i)-1)+3,3*(n_s(i)-1)+1:3*(n_s(i)-1)+3)=K(3*(n_s(i)-1)+1:3*(n_s(i)-1)+3,3*(n_s(i)-1)+1:3*(n_s(i)-1)+3)+K_glob_num(1:3,1:3);
    K(3*(n_s(i)-1)+1:3*(n_s(i)-1)+3,3*(n_e(i)-1)+1:3*(n_e(i)-1)+3)=K(3*(n_s(i)-1)+1:3*(n_s(i)-1)+3,3*(n_e(i)-1)+1:3*(n_e(i)-1)+3)+K_glob_num(1:3,4:6);
    K(3*(n_e(i)-1)+1:3*(n_e(i)-1)+3,3*(n_s(i)-1)+1:3*(n_s(i)-1)+3)=K(3*(n_e(i)-1)+1:3*(n_e(i)-1)+3,3*(n_s(i)-1)+1:3*(n_s(i)-1)+3)+K_glob_num(4:6,1:3);
    K(3*(n_e(i)-1)+1:3*(n_e(i)-1)+3,3*(n_e(i)-1)+1:3*(n_e(i)-1)+3)=K(3*(n_e(i)-1)+1:3*(n_e(i)-1)+3,3*(n_e(i)-1)+1:3*(n_e(i)-1)+3)+K_glob_num(4:6,4:6);
end

A=K;
B=F;
del=0;
%----------------------------Reduce Matrix (only free nodes equations)--------------------------
for i=1:size(nodes,1)    
    j=3*(i-1)+1-del;
	
	if (Sx(i)==1)
        A=[A(1:j-1,:); A(j+1:end,:)];
        A=[A(:,1:j-1) A(:,j+1:end)];
        B=[B(1:j-1); B(j+1:end)];
        del=del+1;
		j=3*(i-1)+1-del;
	end
    if (Sy(i)==1)
        A=[A(1:j,:); A(j+2:end,:)];
        A=[A(:,1:j) A(:,j+2:end)];
        B=[B(1:j); B(j+2:end)];
        del=del+1;
		j=3*(i-1)+1-del;
    end
	
    if (Sphi(i)==1)
        A=[A(1:j+1,:); A(j+3:end,:)];
        A=[A(:,1:j+1) A(:,j+3:end)];
        B=[B(1:j+1); B(j+3:end)];
        del=del+1;
    end
end
%----------------------------------------------------------------------------------------------
k=1;
D_glob=[];
%-----------------------Free Nodes Displacement calc-------------------------------------------

D=inv(A)*B;
%D = mldivide(A,B);

for i=1:3*size(nodes,1)
    if (nodes_free(i) ~=0 && mod(nodes_free(i),3)==1) 
        sprintf('Displacement of Node %d x: %f',fix(i/3)+1,D(k))
        D_glob(i)=D(k);
        k=k+1;
    elseif (nodes_free(i) ~=0 && mod(nodes_free(i),3)==2)
        sprintf('Displacement of Node %d y: %f',fix(i/3)+1,D(k))
        D_glob(i)=D(k);
        k=k+1;
	 elseif (nodes_free(i) ~=0 && mod(nodes_free(i),3)==0)
        sprintf('Displacement of Node %d phi: %f',i/3,D(k))
        D_glob(i)=D(k);
        k=k+1;	
    else
        D_glob(i)=0;
    end
end
D_glob=D_glob';


%==============================================================================================
%-------------------BENDING MOMENT AND SHEAR/NORMAL FORCE FOR ELEMENTS--------------------------
F_elements=[];
f_element_form=K_loc*T;
for i=1:size(n_s,1)
Len=(LenX(i)^2+LenY(i)^2)^0.5;
f_element=E*subs(K_loc,[L I], [Len I_ibeam]);
f_element=f_element*subs(T, [C S], [LenX(i)/Len LenY(i)/Len]);
f_element=f_element*[D_glob(3*(n_s(i)-1)+1:3*(n_s(i)-1)+3); D_glob(3*(n_e(i)-1)+1:3*(n_e(i)-1)+3)];
if ([Fy_w1(i) Fy_w2(i)]~=zeros(1,2))
f_element=f_element-subs(T, [C S], [LenX(i)/Len LenY(i)/Len])*...
[F_dist_forc(3*(n_s(i)-1)+1:3*(n_s(i)-1)+3); F_dist_forc(3*(n_e(i)-1)+1:3*(n_e(i)-1)+3)];
end
		   
F_elements=[F_elements f_element];
end
F_elements=double(F_elements);



%==============================================================================================
%-------------------BENDING MOMENT AND SHEAR FORCE DIAGRAMM------------------------------------
B1=diff(N,x,2);
B2=diff(N,x,3);
M_pl_form=E*I_ibeam*B1*d_beam;
V_pl_form=E*I_ibeam*B2*d_beam;
N_pl_form=-x1+(x2+x1)/L*x;

for i=1:size(n_s,1)
Len=(LenX(i)^2+LenY(i)^2)^0.5;
D_local=subs(T, [C S], [LenX(i)/Len LenY(i)/Len])*...
[D_glob(3*(n_s(i)-1)+1:3*(n_s(i)-1)+3); D_glob(3*(n_e(i)-1)+1:3*(n_e(i)-1)+3)];
M_plot=subs(M_pl_form,d_beam,[D_local(2:3); D_local(5:6)]);
V_plot=subs(V_pl_form,d_beam,[D_local(2:3); D_local(5:6)]);
N_plot=subs(N_pl_form,[x1 x2], [F_elements(1,i) F_elements(4,i)]);
M_plot=subs(M_plot,L,Len);  V_plot=subs(V_plot,L,Len); N_plot=subs(N_plot,L,Len);
x_plot=linspace(0,Len,50);
M_plot=subs(M_plot,x,x_plot); V_plot=subs(V_plot,x,x_plot); N_plot=subs(N_plot,x,x_plot);
M_plot=double(M_plot); V_plot=double(V_plot); N_plot=double(N_plot);

figure(i); 
subplot(3,1,1); 
plot(x_plot,M_plot); xlabel('x direction of beam'); ylabel('Moment (kNm)'); title('Moment Diagram');
subplot(3,1,2);
plot(x_plot,V_plot); xlabel('x direction of beam'); ylabel('Shear Force (kN)'); title('Shear Force Diagram');
subplot(3,1,3)
plot(x_plot,N_plot); xlabel('x direction of beam'); ylabel('Normal Force (kN)'); title('Normal Force Diagram');
end