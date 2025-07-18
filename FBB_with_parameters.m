clc; clear; close all;
if_draw = 0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~
N = 256; 
eta = 0.00888684; 
Nv = 3; 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~

%===========================
c_e_s = 1.7;
rho_w_obj = 0.019; 
r_L = 4; 
sigma_ratio = 4.0;
c_B = 40;  
%---
r_s = 1/2;
H = 0.834; 
r_E = 5/3;
Gamma_1 = 1.0;
%---
n_T_1 = 1;r_n_T = round(1/r_s);
%===========================

r_line = floor(r_L/r_s+0.5);%r_line = floor(r_L+0.5);
r_Np = floor(r_L/r_s/r_line);
r_rho = r_s*r_s*r_L;
if(r_rho~=1) 
    rho_Nv = rho_w_obj/(1-r_rho^Nv)*(1-r_rho)*r_rho^(Nv-1);
else
    rho_Nv = rho_w_obj/Nv;
end
sigma_Nv = c_e_s*eta; 
sigma_1 = sigma_Nv/(r_s^(Nv-1));
dB_Nv = c_B*sigma_Nv; 
dB_1 = dB_Nv/r_s^(Nv-1);
L_Nv_temp = rho_Nv*(2*pi)^3/sigma_Nv^2;
L_1_temp = L_Nv_temp/r_L^(Nv-1);
Np_1 = round(L_1_temp/(1.2*dB_1));
rho_w_Nv = Np_1*(r_Np*r_line)^(Nv-1)*(1.2*dB_Nv)*sigma_Nv^2/(2*pi)^3;
r_Gamma = r_s^((r_E-1)/2)*r_L^(-1/2);

count_line = 0;
m_line = zeros(1,Nv);
for i = 1:Nv
    m_line(i) = 1*r_line^(i-1);
    count_line = count_line + m_line(i);
end

if_sp = 0;

%FBB#############################################################
filenm = ['./centerline_input.dat'];
fid=fopen(filenm,'wt');  %自动生成.dat文件
fprintf(fid,'%18g\n',count_line);
for i = 1:Nv
    for j = 1:m_line(i)
        Np = Np_1*r_Np^(i-1); % 时间点的数量
        fprintf(fid,'%18g\n',Np);
    end
end

for i = 1:Nv
    figure()
    for j_m = 1:m_line(i)
        Np = Np_1*r_Np^(i-1); % 时间点的数量
        dB = dB_1*r_s^(i-1);

        [i Nv j_m m_line(i) Np]
        
        %B = zeros(3,Np);
        FBM = zeros(3,Np);
        BB = zeros(3,Np);
        
        
        for co = 1:3
            %B(co,:) = dB * cumsum(randn(1, Np)); % 累积求和生成布朗运动路径
            if(Np<128)
                r_t = 128/Np; %int
                FBM_temp = zeros(3,128);
                FBM_temp(co,:) = wfbm(H,128);
                for i_t = 1:Np
                    FBM(co,i_t) = FBM_temp(co,ceil(i_t*r_t));
                end
            else
                FBM(co,:) = wfbm(H,Np);
            end
            
        end
    
        
        f_bridge = zeros(1,Np);
        theend  = zeros(1,3);
        for j=1:Np
            f_bridge(j) = 0.5*((1.0*j/Np)^(2.0*H)+1-(1.0-1.0*j/Np)^(2.0*H));
        end
        for co = 1:3
            for j = 1:Np
                %BB(co,j) = B(co,j) - j/Np*(B(co,Np)-B(co,1)); 
                BB(co,j) = FBM(co,j)-f_bridge(j)*(FBM(co,Np)-theend(co));
            end
        end

        %==normalization
        rms = 0;
        for co = 1:3
            for j = 1:(Np-1)
                dx = BB(co,j+1)-BB(co,j);
                rms = rms + (dx)^(2.0)/(Np-1);
            end
        end
        BB = BB*dB/sqrt(rms);
        %==normalization
    
        %BB_p = mod(BB+pi,2*pi) - pi;
        %BB_p = mod(BB+(rand-0.5)*2*pi,2*pi) - pi;
        for co = 1:3
            BB(co,:) = BB(co,:)+(rand-0.5)*2*pi; %!!!
        end
        BB_p = mod(BB,2*pi) - pi;
        
        for j = 1:Np
            fprintf(fid,'%18g',BB_p(1,j));
            fprintf(fid,'%18g',BB_p(2,j));
            fprintf(fid,'%18g\n',BB_p(3,j));
        end
    
        if(if_draw==1)
            for j = 1:Np-1
                dist=sqrt((BB_p(1,j)-BB_p(1,j+1))^2+(BB_p(2,j)-BB_p(2,j+1))^2+(BB_p(3,j)-BB_p(3,j+1))^2);
                if dist<pi
                    plot3(BB_p(1,j:j+1), BB_p(2,j:j+1), BB_p(3,j:j+1),'b-','linewidth',1)
                    hold on;
                end
            end
            plot3([BB_p(1,Np) BB_p(1,1)], [BB_p(2,Np) BB_p(2,1)], [BB_p(3,Np) BB_p(3,1)],'b-','linewidth',1)
            hold on;
            axis([-pi,pi,-pi,pi,-pi,pi]) 
        end
        if(if_draw==2)
            for j = 1:Np-1
                plot3(BB(1,j:j+1), BB(2,j:j+1), BB(3,j:j+1),'b-','linewidth',1)
                %plot3(BB(1,j:j+1), BB(2,j:j+1), BB(3,j:j+1),'ro','linewidth',1)
                hold on;
            end
            plot3([BB(1,Np) BB(1,1)], [BB(2,Np) BB(2,1)], [BB(3,Np) BB(3,1)],'b-','linewidth',1)
            hold on;
        end
    
        if(if_sp == 1)
            S2 = zeros(1,log2(Np/2)+1);
            ref = zeros(1,log2(Np/2)+1);
            for lg2r = 1:log2(Np/2) + 1
                r = 2^(lg2r-1);
                S2(lg2r) = 0;
                for j =1:Np/2
                    S2(lg2r) = S2(lg2r) + (BB(1,j+r)-BB(1,j))^2 + (BB(2,j+r)-BB(2,j))^2 + (BB(3,j+r)-BB(3,j))^2;
                end
                S2(lg2r) = S2(lg2r) / (Np/2);
                ref(lg2r) = 3*(dB*r)^(2.0*H);
            end
            figure()
            plot(0:log2(Np/2), log(S2(:)));
            hold on;
            plot(0:log2(Np/2), log(ref(:)));
        end
    end
   
    box on
    axis square
    %camproj('perspective')
    rotate3d on
    daspect([1 1 1])
    view([1,1,1]);
    x1=xlabel('X轴');        %x轴标题
    x2=ylabel('Y轴');        %y轴标题
    x3=zlabel('Z轴');        %z轴标题
    hold off
    %grid on
end

fclose(fid);

%parameters######################################################
sigma = zeros(1,Nv);
Gamma = zeros(1,Nv);
n_T = zeros(1,Nv);

for i = 1:Nv
    sigma(i) = sigma_1*r_s^(i-1);
    Gamma(i) = Gamma_1*r_Gamma^(i-1);
    n_T(i) = n_T_1*r_n_T^(i-1);
end


ii_line = 0;
for i = 1:Nv
    for j = 1:m_line(i)
        ii_line = ii_line + 1;
    end
end

name_file = ['./parameters_input.dat']
fid=fopen(name_file,'wt');  %自动生成.dat文件
fprintf(fid,'N\n');
fprintf(fid,'%8d\n',N);
fprintf(fid,'N_downsample\n');
fprintf(fid,'%8d\n',N);
fprintf(fid,'line_s\n');
fprintf(fid,'%8d\n',1);
fprintf(fid,'line_e\n');
fprintf(fid,'%8d\n',ii_line);
fprintf(fid,'sigma_ratio\n');
fprintf(fid,'%18.5f\n',sigma_ratio);
fprintf(fid,'eta_t\n');
fprintf(fid,'%18.5f\n',0.0);
fprintf(fid,'---------------\n');
fprintf(fid,'ntf_times\n');
fprintf(fid,'%8d\n',50);
fprintf(fid,'ntf_pre\n');
fprintf(fid,'%8d\n',1000);
fprintf(fid,'Rtube_ratio\n');
fprintf(fid,'%18.5f\n',3.0);
fprintf(fid,'dzeta_sup\n');
fprintf(fid,'%18g\n',0.0001);
fprintf(fid,'---------------\n');
fprintf(fid,'variables="i_line","sigma","Gamma","n_T","if_filter"\n');
i_line = 0;
for i = 1:Nv
    for j = 1:m_line(i)
        i_line = i_line + 1;
        fprintf(fid,'%8d',i_line);
        fprintf(fid,'%18g',sigma(i));
        fprintf(fid,'%18g',Gamma(i));
        fprintf(fid,'%8d',n_T(i));
        if(j<m_line(i))
            if_filter = 0;
        else
            if_filter = 1;
            if_filter = 0;
        end
        fprintf(fid,'%8d\n',if_filter);
    end
end
fclose(fid);

