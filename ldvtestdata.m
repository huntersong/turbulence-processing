data = csvread("PumpB_Double_1.0_x0_y0.csv",6,0,[6 0 7802 6]);

data = csvread("PumpB_Double_1.0_x0_y0.cxv",7,7,[7 7 78 1]);


%isotropic turbulence show
data_t = data(:,2);
data_u = data(:,4);
data_v = data(:,7);
mean_u=mean(data_u);
mean_v=mean(data_v);
u_tran=data_u-mean_u;
v_tran=data_v-mean_v;
plot(data_t,u_tran)
hold on 
plot(data_t,v_tran)

scatter(u_tran,v_tran,".","r")
xlabel({'\it{u_{tran}}','(m/s)'},'FontSize',16,'FontName','Arial','FontWeight','normal','Color','k')
ylabel({'\it{v_{tran}}','(m/s)'},'FontSize',16,'FontName','Arial','FontWeight','normal','Color','k')
axis equal
set(gca,'xtick',[-15 : 2 : 15],'XAxisLocation','origin','FontSize',16,'Fontname', 'Arial','Color','w');
set(gca,'ytick',[-15 : 2 : 15],'YAxisLocation','origin','FontSize',16,'Fontname', 'Arial','Color','w');
axis([-15 15 -15 15])
grid on

%对一条线上的速度作谱
U_reshape = data_u;
V_reshape = data_v;

amplsU = abs(fft(U_reshape)/7797);
amplsV = abs(fft(V_reshape)/7797);

EK_U = amplsU.^2;
EK_V = amplsV.^2;

EK_U = fftshift(EK_U);
EK_V = fftshift(EK_V);

k_uvw = (EK_U+EK_V*2)/2;
figure(2)
n=(7797-1)/2;
m=(7797-1)/2;
b(1,1)=k_uvw((7797-1)/2,1)+k_uvw((7797-1)/2+1,1);
for i=1:((7797-1)/2-1)
   b(i+1,1)=k_uvw(i+(7797-1)/2+1,1)+k_uvw((7797-1)/2-i,1);    
end
k_kolm = 0.01/1000;
k_x = [(2*3.14/(0.005/100)):(2*3.14/(0.005/100)):(2*3.14/(0.005/100))*length(b)]';
k_xx = k_x*k_kolm;   %无量纲的横坐标
b_kk = b/(k_kolm*k_kolm^2);%无量纲的纵坐标
loglog(k_xx,b_kk)
xlabel({'\it{kη}'},'FontSize',16,'FontName','Arial','FontWeight','normal','Color','k')
ylabel({'\it{E(k)/ηu^2}'},'FontSize',16,'FontName','Arial','FontWeight','normal','Color','k')
set(gca,'FontSize',16,'Fontname', 'Arial','Color','w');
set(gca,'FontSize',16,'Fontname', 'Arial','Color','w');
axis([1 10000 10^10 10^16])
grid on
hold on
x=100:10^4;
y=10*10^17*x.^(-5/3);
loglog(x,y)


%对一个平面上的速度作谱
n=1;
U_reshape = u;
V_reshape = v;
W_reshape = w;
U_reshape_temp = U_reshape(:,n,:);
V_reshape_temp = V_reshape(:,n,:);
W_reshape_temp = W_reshape(:,n,:);
amplsU = abs(fftn(U_reshape_temp)/100);
amplsV = abs(fftn(V_reshape_temp)/100);
amplsW = abs(fftn(W_reshape_temp)/100);
EK_U = amplsU.^2;
EK_V = amplsV.^2;
EK_W = amplsW.^2;
EK_U = fftshift(EK_U);
EK_V = fftshift(EK_V);
EK_W = fftshift(EK_W);
k_uvw = (EK_U+EK_V+EK_W)/2;
box_sidex = 100;  box_sidey = 100;   box_sidez = 100;
centerx = 50;    centery = 50;     centerz = 50;
%radius = round(29*sqrt(3))-1;
for i=1:box_sidex 
    for j=1:box_sidey
        
        A(i,j)=round(sqrt((i-centerx).^2+(j-centery).^2));  %建立位置信息
        
    end
end

a=0;
for r=1:50
   for i=1:100
      for j=1:100
          
            if A(i,j)==r
             a=k_uvw(i,j)+a;
            end
          
      end
   end
   k_K(r,1)=a;
   a=0;
end
loglog(k_K)
hold on
x=1:100;
y=10*10^0*x.^(-5/3);
loglog(x,y)




%100x100x100
path1 = 'F:\reynoldsstress\bubble_breakup_in_cavity\LES_bubbleCavity_liquid_symmetry_refine\fluent0.224csv\'; %设置工作路径
getfilename1 = ls([path1,'*.csv']');                       %列出工作目录中所有csv文件名称
filename1 = cellstr ( getfilename1 ) ;                     %将csv名称文件设为指针
for m=1:100
    a=(m-1)*100+1;b=(m-1)*100+100;
for n=a:b  %length(filename1)
    filename_cur1 = char(filename1(n));
    filename_full1 = [path1, filename_cur1];
    p=n-(m-1)*100; q=m;
    u(:,p,q)=csvread( filename_full1 , 6 ,3 ,[ 6 , 3 , 105 , 3]);
    v(:,p,q)=csvread( filename_full1 , 6 ,4 ,[ 6 , 4 , 105 , 4]);
    w(:,p,q)=csvread( filename_full1 , 6 ,5 ,[ 6 , 5 , 105 , 5]);
end
end
%三维区域谱
U_reshape = u;
V_reshape = v;
W_reshape = w;
amplsU = abs(fftn(U_reshape)/100);
amplsV = abs(fftn(V_reshape)/100);
amplsW = abs(fftn(W_reshape)/100);
EK_U = amplsU.^2;
EK_V = amplsV.^2;
EK_W = amplsW.^2;
EK_U = fftshift(EK_U);
EK_V = fftshift(EK_V);
EK_W = fftshift(EK_W);
k_uvw = (EK_U+EK_V+EK_W)/2;
box_sidex = length(k_uvw);  box_sidey = length(k_uvw);   box_sidez = length(k_uvw);
centerx = length(k_uvw)/2;    centery = length(k_uvw)/2;     centerz = length(k_uvw)/2;
%radius = round(29*sqrt(3))-1;
for i=1:box_sidex 
    for j=1:box_sidey
        for k=1:box_sidez
        A(i,j,k)=round(sqrt((i-centerx).^2+(j-centery).^2+(k-centerz).^2));  %建立位置信息
        end
    end
end

a=0;
for r=1:max(max(max(A)))
   for i=1:length(A)
      for j=1:length(A)
          for k=1:length(A)
            if A(i,j,k)==r
             a=k_uvw(i,j,k)+a;
            end
          end
      end
   end
   k_K(r,1)=a;
   a=0;
end

k_kolm = 0.01/1000;
u_kolm = 0.35;
k_x = [(2*3.14/(0.005/100)):(2*3.14/(0.005/100)):(2*3.14/(0.005/100))*87]';
k_xx = k_x*k_kolm;   %无量纲的横坐标
k_KK = k_K/(k_kolm*u_kolm^2);%无量纲的纵坐标
figure(1)
loglog(k_xx,k_KK)                  %k_K_filter = sgolayfilt(k_K,1,5);  hold on;   loglog(k_K_filter)
xlabel({'\it{kη}'},'FontSize',16,'FontName','Arial','FontWeight','normal','Color','k')
ylabel({'\it{E(k)/ηu^2}'},'FontSize',16,'FontName','Arial','FontWeight','normal','Color','k')
set(gca,'FontSize',10,'Fontname', 'Arial','Color','w');
set(gca,'FontSize',10,'Fontname', 'Arial','Color','w');
axis([1 120 10^0 10^15])
grid on
hold on
x=1:80;
y=10*10^12*x.^(-5/3);
loglog(x,y)





















%84x84x84
path1 = 'F:\reynoldsstress\turbulent energy spectrum\volume_uvw_100\'; %设置工作路径
getfilename1 = ls([path1,'*.csv']');                       %列出工作目录中所有csv文件名称
filename1 = cellstr ( getfilename1 ) ;                     %将csv名称文件设为指针
for m=27:84
    
    a=m*100+27;b=m*100+84;
for n=a:b%length(filename1)
    filename_cur1 = char(filename1(n));
    filename_full1 = [path1, filename_cur1];
    p=n-m*100-26; q=m-26;
    u(:,p,q)=csvread( filename_full1 , 27 ,3 ,[ 27 , 3 , 84 , 3]);
    v(:,p,q)=csvread( filename_full1 , 27 ,4 ,[ 27 , 4 , 84 , 4]);
    w(:,p,q)=csvread( filename_full1 , 27 ,5 ,[ 27 , 5 , 84 , 5]);
end
end
%三维区域谱
U_reshape = u;
V_reshape = v;
W_reshape = w;
amplsU = abs(fftn(U_reshape)/58);
amplsV = abs(fftn(V_reshape)/58);
amplsW = abs(fftn(W_reshape)/58);
EK_U = amplsU.^2;
EK_V = amplsV.^2;
EK_W = amplsW.^2;
EK_U = fftshift(EK_U);
EK_V = fftshift(EK_V);
EK_W = fftshift(EK_W);
k_uvw = (EK_U+EK_V+EK_W)/2;
box_sidex = 58;  box_sidey = 58;   box_sidez = 58;
centerx = 29;    centery = 29;     centerz = 29;
radius = round(29*sqrt(3))-1;
for i=1:box_sidex 
    for j=1:box_sidey
        for k=1:box_sidez
        A(i,j,k)=round(sqrt((i-centerx).^2+(j-centery).^2+(k-centerz).^2));  %建立位置信息
        end
    end
end

a=0;
for r=1:29
   for i=1:58
      for j=1:58
          for k=1:58
            if A(i,j,k)==r
             a=k_uvw(i,j,k)+a;
            end
          end
      end
   end
   k_K(r,1)=a;
   a=0;
end
loglog(k_K)
hold on
x=1:100;
y=10*10^4*x.^(-5/3);
loglog(x,y)

%平面速度谱
U_reshape_temp = U_reshape(:,1,:);
V_reshape_temp = V_reshape(:,1,:);
W_reshape_temp = W_reshape(:,1,:);
amplsU = abs(fftn(U_reshape_temp)/84);
amplsV = abs(fftn(V_reshape_temp)/84);
amplsW = abs(fftn(W_reshape_temp)/84);
EK_U = amplsU.^2;
EK_V = amplsV.^2;
EK_W = amplsW.^2;
EK_U = fftshift(EK_U);
EK_V = fftshift(EK_V);
EK_W = fftshift(EK_W);
k_uvw = (EK_U+EK_V+EK_W)/2;
box_sidex = 84;  box_sidey = 84;   box_sidez = 84;
centerx = 42;    centery = 42;     centerz = 42;
radius = round(42*sqrt(3))-1;
for i=1:box_sidex 
    for j=1:box_sidey
        
        A(i,j)=round(sqrt((i-centerx).^2+(j-centery).^2));  %建立位置信息
        
    end
end

a=0;
for r=1:42
   for i=1:84
      for j=1:84
          
            if A(i,j)==r
             a=k_uvw(i,j)+a;
            end
          
      end
   end
   k_K(r,1)=a;
   a=0;
end
loglog(k_K)
hold on
x=1:100;
y=10*10^0*x.^(-5/3);
loglog(x,y)



%{
for i=1:100
    for j=1:100
        for k=1:100
            num = (i-1)*10000+(j-1)*100+k;
   data_series_u(num,:) = u(i,j,k);
   data_series_v(num,:) = v(i,j,k);
   data_series_w(num,:) = w(i,j,k);
        end
    end
end
data_series(:,1) = data_series_u;
data_series(:,2) = data_series_v;
data_series(:,3) = data_series_w;
save data_series -ascii data_series
%}



