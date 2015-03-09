%---------------------------Timing Recovery Problem------------------%

clc
clear all
close all
d = 6; %Delay
hs = 8; % 8 samples per symbol
hd = 1; 
hs_1=64; 

alpha=0.40;% Excess Bandwidth
rc_1=rcosine(hd,hs,'sqrt',alpha,d); 
rc_1_1=rc_1/max(rc_1); 

rc_2=rcosine(hd,hs_1,'sqrt',alpha,d); 
rc_2_1=rc_2/(rc_1_1(1:4:97)*rc_2(1:32:769)');

drc_2_1_x=conv(rc_2_1,[1 0 -1]*64/2);
drc_2_1=drc_2_1_x(2:770);

%QPSK Symbols%

number_of_symbols=2000;
QPSK=(floor(2*rand(1,number_of_symbols))-0.5)/0.5+...
    1j*(floor(2*rand(1,number_of_symbols))-0.5)/0.5;

h=reshape(rc_1_1(1:96),8,12);
g1=reshape(rc_2_1(1:768),32,24);
dg1=reshape(drc_2_1(1:768),32,24);

r1=zeros(1,12);
QPSK_1=zeros(1,8*number_of_symbols);
m=1;
for n=1:number_of_symbols
    r1=[QPSK(n) r1(1:11)];
    for k=1:8
        QPSK_1(m)=r1*h(k,:)';
        m=m+1;
    end
end
QPSK_2=QPSK_1(2:4:8*number_of_symbols);
for k=1:32
     y=conv(QPSK_2, g1(k,:));
     y_dot=conv(QPSK_2,dg1(k,:));
     yy_dot_product=real(y(2:2:number_of_symbols)).*real(y_dot(2:2:number_of_symbols));
end

figure(1)
subplot(2,1,1)
plot(0,0)
hold on
for n=2:4:number_of_symbols-4
    plot(-1:1/2:1,real(y(n:n+4)))
end
hold off
grid on
title('Input for the Eye Diagram of the Timing Loop')
xlabel('Time Index')
ylabel('Amplitude')
subplot(2,1,2)
plot(y(2:4:number_of_symbols),'rx')
grid on
axis('equal')
axis([-1.5 1.5 -1.5 1.5])
title('Inputp for the Constellation Diagram of the Timing Loop')

%%%%%%%%%%%%%%% Timing Recovery Loop %%%%%%%%%%%%%%%

accum=16;  
int=0;    
r1=zeros(1,24);
theta_0=2*pi/400;    
eta=sqrt(2)/2;
eta=4*eta;
ki= (4*theta_0*theta_0)/(1+2*eta*theta_0+theta_0*theta_0);
kp= (4*eta*theta_0)/(1+2*eta*theta_0+theta_0*theta_0);

m=1;
for n=1:2:2*number_of_symbols-2
    r1=[QPSK_2(n) r1(1:23)];
    pointer=floor(accum);
    pointer_save(m)=pointer;
    y_1(n)=r1*g1(pointer,:)';
    y_1_dot(n)=r1*dg1(pointer,:)';
    
    r1=[QPSK_2(n+1) r1(1:23)];
    y_1(n+1)=r1*g1(pointer,:)';
     y_1_dot(n+1)=r1*dg1(pointer,:)';
    
    y_det(m)=real(y_1(n+1))*real( y_1_dot(n+1));
    int=int+y_det(m)*ki;
    loop_filter(m)=int+y_det(m)*kp;
    
    accumulator_save(m)=accum;
    accum=accum+loop_filter(m);
    if accum>=33
        accum=accum-1;
    end
    m=m+1;
end

figure(2)
subplot(2,1,1)
plot(0,0)
hold on
for n=2:4:number_of_symbols-4
    plot(-1:1/2:1,real(y_1(n:n+4)))
end
hold off
grid on
title('Output of the Eye Diagram of the Timing Loop')
xlabel('Time Index')
ylabel('Amplitude')
subplot(2,1,2)
plot(y_1(2:4:number_of_symbols),'rx')
grid on
axis('equal')
axis([-1.5 1.5 -1.5 1.5])
title('Output of the Constellation Diagram of the Timing Loop')

figure(3)

subplot(2,1,1)
plot(loop_filter)
grid on
title('Time Profile of the Loop Filter Output ')
xlabel('Time Index')
ylabel('Amplitude')


subplot(2,1,2)
plot(accumulator_save)
hold on
plot(floor(accumulator_save),'r')
hold off
grid on
title('The phase accumulator output and the integer part of the phase accumulator output.')
xlabel('Time Index')
ylabel('Amplitude')
