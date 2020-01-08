function [PedL1, PedR1, PedL2, PedR2,PedL3,PedR3] = sixPeds5(t0, y0, u0, v0, u1, u2, dt, alpha, beta, gamma, sigma_x, sigma_y,kappa, theta1, theta2, offset,store1)
a=[0.001 0.001];
b=[0.001 0.001];


%% random numbers to determine how strong the attraction or avoidance is for a person near store 1
L1AttractionVar1=10*rand(1)+3;
L2AttractionVar1=10*rand(1)+3;
R1AttractionVar1=10*rand(1)+3;
R2AttractionVar1=10*rand(1)+3;
L3AttractionVar1=10*rand(1)+3;
R3AttractionVar1=10*rand(1)+3;

%% random numbers to determine how strong the attraction or avoidance is for a person near store 2
L1AttractionVar2=10*rand(1)+3;
R1AttractionVar2=10*rand(1)+3;
L2AttractionVar2=10*rand(1)+3;
R2AttractionVar2=10*rand(1)+3;
L3AttractionVar2=10*rand(1)+3;
R3AttractionVar2=10*rand(1)+3;

%this section helps randomly choose whether a person is attracted to store1
%or is trying to avoid it.
randL1=rand(1);
if(randL1<.5)
    AttAv1L1=-1;
else
    AttAv1L1=1;
end

randL2=rand(1);
if(randL2<.5)
    AttAv1L2=-1;
else
    AttAv1L2=1;
end

randL3=rand(1);
if(randL3<.5)
    AttAv1L3=-1;
else
    AttAv1L3=1;
end

randR1=rand(1);
if(randR1<.5)
    AttAv1R1=-1;
else
    AttAv1R1=1;
end

randR2=rand(1);
if(randR2<.5)
    AttAv1R2=-1;
else
    AttAv1R2=1;
end

randR3=rand(1);
if(randR3<.5)
    AttAv1R3=-1;
else
    AttAv1R3=1;
end

%this section helps randomly choose whether a person is attracted to store2
%or is trying to avoid it.
randL1=rand(1);
if(randL1<.5)
    AttAv2L1=-1;
else
    AttAv2L1=1;
end

randL2=rand(1);
if(randL2<.5)
    AttAv2L2=-1;
else
    AttAv2L2=1;
end


if(randL3<.5)
    AttAv2L3=-1;
else
    AttAv2L3=1;
end

randR1=rand(1);
if(randR1<.5)
    AttAv2R1=-1;
else
    AttAv2R1=1;
end

randR2=rand(1);
if(randR2<.5)
    AttAv2R2=-1;
else
    AttAv2R2=1;
end

randR3=rand(1);
if(randR3<.5)
    AttAv2R3=-1;
else
    AttAv2R3=1;
end

Store1Position=[2,2.5];
Store2Position=[7, -.5];


%% The mean paths for people going to thw left and the right
XR1(1) = 0;
XR2(1) = 0;
y_starR = @(x) -0.001*x.*x-0.001*x-0.3146;
dy_starR = @(x) -2*0.001*x-0.001;
ddy_starR = -2*0.001;


XL1(1) = 10.8;
XL2(1) = 10.8;
y_starL = @(x) -0.001*x.*x-0.001*x+.41896;
dy_starL = @(x) -2.*0.001*x-0.001;
ddy_starL = -2*0.001;


%This counter variable helps determine how long a person has been in store
%1. 
       time1L1=0;
       time1L2=0;
       time1L3=0;
       time1R1=0;
       time1R2=0;
       time1R3=0;
 %This counter variable helps determine how long a person has been in store
%2 
       time2L1=0;
       time2L2=0;
       time2L3=0;
       time2R1=0;
       time2R2=0;
       time2R3=0;
       
       
    stepNumber =  round(t0(2)/dt);
    
    
    
    YL1(1) = y0(1);
    UL1(1) = u0(1);
    VL1(1) = v0(1);
    
    XL2(1) = NaN;
    YL2(1) = NaN;
    UL2(1) = NaN;
    VL2(1) = NaN;
    
    XR1(1) = NaN;
    YR1(1) = NaN;
    UR1(1) = NaN;
    VR1(1) = NaN;
    
    XR2(1) = NaN;
    YR2(1) = NaN;
    UR2(1) = NaN;
    VR2(1) = NaN;
    
    XL3(1) = NaN;
    YL3(1) = NaN;
    UL3(1) = NaN;
    VL3(1) = NaN;
    
    XR3(1) = NaN;
    YR3(1) = NaN;
    UR3(1) = NaN;
    VR3(1) = NaN;
   
    
    for j=1:stepNumber-1;
      
        
        norL1 = [dy_starL(XL1(j)) , -1]./sqrt(1+dy_starL(XL1(j)).^2);
        rhoL1 = abs((1+dy_starL(XL1(j))^2).^(1.5)./ddy_starL);
        centripetalL1 = (UL1(j).^2+VL1(j).^2)/rhoL1*norL1;
        %centripetal = [0 0];
        
        XL1(j+1) = XL1(j) + UL1(j)*dt(1);
        UL1(j+1) = UL1(j) - alpha(1)*UL1(j)*(UL1(j)-u1)*(UL1(j)-u2)*dt + centripetalL1(1)*dt + sigma_x(1)*sqrt(dt)*randn(1);
        YL1(j+1) = YL1(j) + VL1(j)*dt(1);
        VL1(j+1) = VL1(j) - beta(1)*(YL1(j)-y_starL(XL1(j)))*dt - gamma(1)*VL1(j)*dt + centripetalL1(2)*dt + sigma_y(1)*sqrt(dt)*randn(1);
       
        
        XR1(j+1) = NaN;
        YR1(j+1) = NaN;
        UR1(j+1) = NaN;
        VR1(j+1) = NaN;
        
        XR2(j+1) = NaN;
        YR2(j+1) = NaN;
        UR2(j+1) = NaN;
        VR2(j+1) = NaN;
    
        XL2(j+1) = NaN;
        YL2(j+1) = NaN;
        UL2(j+1) = NaN;
        VL2(j+1) = NaN;
        
        
        XR3(j+1) = NaN;
        YR3(j+1) = NaN;
        UR3(j+1) = NaN;
        VR3(j+1) = NaN;
    
        XL3(j+1) = NaN;
        YL3(j+1) = NaN;
        UL3(j+1) = NaN;
        VL3(j+1) = NaN;
        
        
     
    end
   
    j=stepNumber;
    
    XL1(j+1) = XL1(j) + UL1(j)*dt;
    UL1(j+1) = UL1(j) - alpha(1)*UL1(j)*(UL1(j)-u1)*(UL1(j)-u2)*dt + centripetalL1(1)*dt + sigma_x(1)*sqrt(dt)*randn(1);
    YL1(j+1) = YL1(j) + VL1(j)*dt(1);
    VL1(j+1) = VL1(j) - beta(1)*(YL1(j)-y_starL(XL1(j)))*dt - gamma(1)*VL1(j)*dt + centripetalL1(2)*dt + sigma_y(1)*sqrt(dt)*randn(1);
    
    XR1(j+1)= 0;
    YR1(j+1) = y0(2);
    UR1(j+1) = u0(2);
    VR1(j+1) = v0(2);
    
    XR2(j+1) = NaN;
    YR2(j+1) = NaN;
    UR2(j+1) = NaN;
    VR2(j+1) = NaN;
    
    XL2(j+1) = NaN;
    YL2(j+1) = NaN;
    UL2(j+1) = NaN;
    VL2(j+1) = NaN;
    
        XR3(j+1) = NaN;
        YR3(j+1) = NaN;
        UR3(j+1) = NaN;
        VR3(j+1) = NaN;
    
        XL3(j+1) = NaN;
        YL3(j+1) = NaN;
        UL3(j+1) = NaN;
        VL3(j+1) = NaN;
    



stepNumber2 =  round(abs((t0(2)-t0(3))/dt));

 for j=numel(XL1):numel(XL1)+stepNumber2-1
     
        norL1 = [dy_starL(XL1(j)) , -1]./sqrt(1+dy_starL(XL1(j)).^2);
        rhoL1 = abs((1+dy_starL(XL1(j))^2).^(1.5)./ddy_starL);
        centripetalL1 = (UL1(j).^2+VL1(j).^2)/rhoL1*norL1;
        
        norR1 = [dy_starR(XR1(j)) , -1]./sqrt(1+dy_starR(XR1(j)).^2);
        rhoR1 = abs((1+dy_starR(XR1(j))^2).^(1.5)./ddy_starR);
        centripetalR1 = (UR1(j).^2+VR1(j).^2)/rhoR1*norR1;
        
     XL1(j+1) = XL1(j) + UL1(j)*dt(1);
     UL1(j+1) = UL1(j) - alpha(1)*UL1(j)*(UL1(j)-u1)*(UL1(j)-u2)*dt + centripetalL1(1)*dt + sigma_x(1)*sqrt(dt)*randn(1);
     YL1(j+1) = YL1(j) + VL1(j)*dt(1);
     VL1(j+1) = VL1(j) - beta(1)*(YL1(j)-y_starL(XL1(j)))*dt - gamma(1)*VL1(j)*dt + centripetalL1(2)*dt + sigma_y(1)*sqrt(dt)*randn(1);
     
        XR1(j+1) = XR1(j) + UR1(j)*dt(1);
        UR1(j+1) = UR1(j) - alpha(2)*UR1(j)*(UR1(j)-u1)*(UR1(j)-u2)*dt(1) + centripetalR1(1)*dt(1) + sigma_x(2)*sqrt(dt(1))*randn(1);
        YR1(j+1) = YR1(j) + VR1(j)*dt(1);
        VR1(j+1) = VR1(j) - beta(2)*(YR1(j)-y_starR(XR1(j)))*dt(1) - gamma(2)*VR1(j)*dt(1) + centripetalR1(2)*dt(1) + sigma_y(2)*sqrt(dt(1))*randn(1);
      
    XR2(j+1) = NaN;
    YR2(j+1) = NaN;
    UR2(j+1) = NaN;
    VR2(j+1) = NaN;
    
    XL2(j+1) = NaN;
    YL2(j+1) = NaN;
    UL2(j+1) = NaN;
    VL2(j+1) = NaN;
    
     XR3(j+1) = NaN;
     YR3(j+1) = NaN;
     UR3(j+1) = NaN;
     VR3(j+1) = NaN;
    
     XL3(j+1) = NaN;
     YL3(j+1) = NaN;
     UL3(j+1) = NaN;
     VL3(j+1) = NaN;

    
 end
 
 j=numel(XL1)
 
           norL1 = [dy_starL(XL1(j)) , -1]./sqrt(1+dy_starL(XL1(j)).^2);
        rhoL1 = abs((1+dy_starL(XL1(j))^2).^(1.5)./ddy_starL);
        centripetalL1 = (UL1(j).^2+VL1(j).^2)/rhoL1*norL1;
        
        norR1 = [dy_starR(XR1(j)) , -1]./sqrt(1+dy_starR(XR1(j)).^2);
        rhoR1 = abs((1+dy_starR(XR1(j))^2).^(1.5)./ddy_starR);
        centripetalR1 = (UR1(j).^2+VR1(j).^2)/rhoR1*norR1;
 
     XL1(j+1) = XL1(j) + UL1(j)*dt;
     UL1(j+1) = UL1(j) - alpha(1)*UL1(j)*(UL1(j)-u1)*(UL1(j)-u2)*dt + centripetalL1(1)*dt + sigma_x(1)*sqrt(dt)*randn(1);
     YL1(j+1) = YL1(j) + VL1(j)*dt(1);
     VL1(j+1) = VL1(j) - beta(1)*(YL1(j)-y_starL(XL1(j)))*dt - gamma(1)*VL1(j)*dt + centripetalL1(2)*dt + sigma_y(1)*sqrt(dt)*randn(1);
     
        XR1(j+1) = XR1(j) + UR1(j)*dt(1);
        UR1(j+1) = UR1(j) - alpha(2)*UR1(j)*(UR1(j)-u1)*(UR1(j)-u2)*dt(1) + centripetalR1(1)*dt(1) + sigma_x(2)*sqrt(dt(1))*randn(1);
        YR1(j+1) = YR1(j) + VR1(j)*dt(1);
        VR1(j+1) = VR1(j) - beta(2)*(YR1(j)-y_starR(XR1(j)))*dt(1) - gamma(2)*VR1(j)*dt(1) + centripetalR1(2)*dt(1) + sigma_y(2)*sqrt(dt(1))*randn(1);
      
     XL2(j+1)= 10.8;
    YL2(j+1) = y0(3);
    UL2(j+1) = u0(3);
    VL2(j+1) = v0(3);
 
     XR2(j+1) = NaN;
    YR2(j+1) = NaN;
    UR2(j+1) = NaN;
    VR2(j+1) = NaN;
    
        XR3(j+1) = NaN;
        YR3(j+1) = NaN;
        UR3(j+1) = NaN;
        VR3(j+1) = NaN;
    
        XL3(j+1) = NaN;
        YL3(j+1) = NaN;
        UL3(j+1) = NaN;
        VL3(j+1) = NaN;
      
 stepNumber3 =  round(abs((t0(3)-t0(4))/dt));
 
for j=numel(XL1):numel(XL1)+stepNumber3-1
 
         norL1 = [dy_starL(XL1(j)) , -1]./sqrt(1+dy_starL(XL1(j)).^2);
        rhoL1 = abs((1+dy_starL(XL1(j))^2).^(1.5)./ddy_starL);
        centripetalL1 = (UL1(j).^2+VL1(j).^2)/rhoL1*norL1;
        
        norR1 = [dy_starR(XR1(j)) , -1]./sqrt(1+dy_starR(XR1(j)).^2);
        rhoR1 = abs((1+dy_starR(XR1(j))^2).^(1.5)./ddy_starR);
        centripetalR1 = (UR1(j).^2+VR1(j).^2)/rhoR1*norR1;
        
        norL2 = [dy_starL(XL2(j)) , -1]./sqrt(1+dy_starL(XL2(j)).^2);
        rhoL2 = abs((1+dy_starL(XL2(j))^2).^(1.5)./ddy_starL);
        centripetalL2 = (UL2(j).^2+VL2(j).^2)/rhoL2*norL2;
        
 XL1(j+1) = XL1(j) + UL1(j)*dt(1);
 UL1(j+1) = UL1(j) - alpha(1)*UL1(j)*(UL1(j)-u1)*(UL1(j)-u2)*dt + centripetalL1(1)*dt + sigma_x(1)*sqrt(dt)*randn(1);
 YL1(j+1) = YL1(j) + VL1(j)*dt(1);
 VL1(j+1) = VL1(j) - beta(1)*(YL1(j)-y_starL(XL1(j)))*dt - gamma(1)*VL1(j)*dt + centripetalL1(2)*dt + sigma_y(1)*sqrt(dt)*randn(1);
     
  XR1(j+1) = XR1(j) + UR1(j)*dt(1);
  UR1(j+1) = UR1(j) - alpha(2)*UR1(j)*(UR1(j)-u1)*(UR1(j)-u2)*dt(1) + centripetalR1(1)*dt(1) + sigma_x(2)*sqrt(dt(1))*randn(1);
  YR1(j+1) = YR1(j) + VR1(j)*dt(1);
  VR1(j+1) = VR1(j) - beta(2)*(YR1(j)-y_starR(XR1(j)))*dt - gamma(2)*VR1(j)*dt(1) + centripetalR1(2)*dt(1) + sigma_y(2)*sqrt(dt(1))*randn(1);
 
  XL2(j+1) = XL2(j) + UL2(j)*dt;
 UL2(j+1) = UL2(j) - alpha(3)*UL2(j)*(UL2(j)-u1)*(UL2(j)-u2)*dt + centripetalL2(1)*dt + sigma_x(3)*sqrt(dt)*randn(1);
 YL2(j+1) = YL2(j) + VL2(j)*dt;
 VL2(j+1) = VL2(j) - beta(3)*(YL2(j)-y_starL(XL2(j)))*dt - gamma(3)*VL2(j)*dt + centripetalL2(2)*dt + sigma_y(3)*sqrt(dt)*randn(1);
 
 
    XR2(j+1) = NaN;
    YR2(j+1) = NaN;
    UR2(j+1) = NaN;
    VR2(j+1) = NaN;
    
        XR3(j+1) = NaN;
        YR3(j+1) = NaN;
        UR3(j+1) = NaN;
        VR3(j+1) = NaN;
    
        XL3(j+1) = NaN;
        YL3(j+1) = NaN;
        UL3(j+1) = NaN;
        VL3(j+1) = NaN;
    
end 
 j=numel(XL1);
 
          norL1 = [dy_starL(XL1(j)) , -1]./sqrt(1+dy_starL(XL1(j)).^2);
        rhoL1 = abs((1+dy_starL(XL1(j))^2).^(1.5)./ddy_starL);
        centripetalL1 = (UL1(j).^2+VL1(j).^2)/rhoL1*norL1;
        
        norR1 = [dy_starR(XR1(j)) , -1]./sqrt(1+dy_starR(XR1(j)).^2);
        rhoR1 = abs((1+dy_starR(XR1(j))^2).^(1.5)./ddy_starR);
        centripetalR1 = (UR1(j).^2+VR1(j).^2)/rhoR1*norR1;
        
        norL2 = [dy_starL(XL2(j)) , -1]./sqrt(1+dy_starL(XL2(j)).^2);
        rhoL2 = abs((1+dy_starL(XL2(j))^2).^(1.5)./ddy_starL);
        centripetalL2 = (UL2(j).^2+VL2(j).^2)/rhoL2*norL2;
 
     XL1(j+1) = XL1(j) + UL1(j)*dt(1);
     UL1(j+1) = UL1(j) - alpha(1)*UL1(j)*(UL1(j)-u1)*(UL1(j)-u2)*dt + centripetalL1(1)*dt + sigma_x(1)*sqrt(dt)*randn(1);
     YL1(j+1) = YL1(j) + VL1(j)*dt(1);
     VL1(j+1) = VL1(j) - beta(1)*(YL1(j)-y_starL(XL1(j)))*dt - gamma(1)*VL1(j)*dt + centripetalL1(2)*dt + sigma_y(1)*sqrt(dt)*randn(1);
     
        XR1(j+1) = XR1(j) + UR1(j)*dt(1);
        UR1(j+1) = UR1(j) - alpha(2)*UR1(j)*(UR1(j)-u1)*(UR1(j)-u2)*dt(1) + centripetalR1(1)*dt(1) + sigma_x(2)*sqrt(dt(1))*randn(1);
        YR1(j+1) = YR1(j) + VR1(j)*dt(1);
        VR1(j+1) = VR1(j) - beta(2)*(YR1(j)-y_starR(XR1(j)))*dt(1) - gamma(2)*VR1(j)*dt(1) + centripetalR1(2)*dt(1) + sigma_y(2)*sqrt(dt(1))*randn(1);
      
     XL2(j+1)= XL2(j) + UL2(j)*dt;
    UL2(j+1) = UL2(j) - alpha(3)*UL2(j)*(UL2(j)-u1)*(UL2(j)-u2)*dt + centripetalL2(1)*dt + sigma_x(3)*sqrt(dt)*randn(1);
    YL2(j+1) = YL2(j) + VL2(j)*dt;
    VL2(j+1) =  VL2(j) - beta(3)*(YL2(j)-y_starL(XL2(j)))*dt - gamma(3)*VL2(j)*dt + centripetalL2(2)*dt + sigma_y(3)*sqrt(dt)*randn(1);
 
     XR2(j+1)= 0;
    YR2(j+1) = y0(4);
    UR2(j+1) = u0(4);
    VR2(j+1) = v0(4);
    
        XR3(j+1) = NaN;
        YR3(j+1) = NaN;
        UR3(j+1) = NaN;
        VR3(j+1) = NaN;
    
        XL3(j+1) = NaN;
        YL3(j+1) = NaN;
        UL3(j+1) = NaN;
        VL3(j+1) = NaN;
    
     stepNumber4 =  round(abs((t0(4)-t0(5))/dt));
    
     
     for j=numel(XL1):numel(XL1)+stepNumber4-1
 
         norL1 = [dy_starL(XL1(j)) , -1]./sqrt(1+dy_starL(XL1(j)).^2);
        rhoL1 = abs((1+dy_starL(XL1(j))^2).^(1.5)./ddy_starL);
        centripetalL1 = (UL1(j).^2+VL1(j).^2)/rhoL1*norL1;
        
        norR1 = [dy_starR(XR1(j)) , -1]./sqrt(1+dy_starR(XR1(j)).^2);
        rhoR1 = abs((1+dy_starR(XR1(j))^2).^(1.5)./ddy_starR);
        centripetalR1 = (UR1(j).^2+VR1(j).^2)/rhoR1*norR1;
        
        norL2 = [dy_starL(XL2(j)) , -1]./sqrt(1+dy_starL(XL2(j)).^2);
        rhoL2 = abs((1+dy_starL(XL2(j))^2).^(1.5)./ddy_starL);
        centripetalL2 = (UL2(j).^2+VL2(j).^2)/rhoL2*norL2;
        
        norR2 = [dy_starR(XR2(j)) , -1]./sqrt(1+dy_starR(XR2(j)).^2);
       rhoR2 = abs((1+dy_starR(XR2(j))^2).^(1.5)./ddy_starR);
       centripetalR2 = (UR2(j).^2+VR2(j).^2)/rhoR2*norR2;
       
 XL1(j+1) = XL1(j) + UL1(j)*dt(1);
 UL1(j+1) = UL1(j) - alpha(1)*UL1(j)*(UL1(j)-u1)*(UL1(j)-u2)*dt + centripetalL1(1)*dt + sigma_x(1)*sqrt(dt)*randn(1);
 YL1(j+1) = YL1(j) + VL1(j)*dt(1);
 VL1(j+1) = VL1(j) - beta(1)*(YL1(j)-y_starL(XL1(j)))*dt - gamma(1)*VL1(j)*dt + centripetalL1(2)*dt + sigma_y(1)*sqrt(dt)*randn(1);
     
  XR1(j+1) = XR1(j) + UR1(j)*dt(1);
  UR1(j+1) = UR1(j) - alpha(2)*UR1(j)*(UR1(j)-u1)*(UR1(j)-u2)*dt(1) + centripetalR1(1)*dt(1) + sigma_x(2)*sqrt(dt(1))*randn(1);
  YR1(j+1) = YR1(j) + VR1(j)*dt(1);
  VR1(j+1) = VR1(j) - beta(2)*(YR1(j)-y_starR(XR1(j)))*dt - gamma(2)*VR1(j)*dt(1) + centripetalR1(2)*dt(1) + sigma_y(2)*sqrt(dt(1))*randn(1);
 
  XL2(j+1) = XL2(j) + UL2(j)*dt;
 UL2(j+1) = UL2(j) - alpha(3)*UL2(j)*(UL2(j)-u1)*(UL2(j)-u2)*dt + centripetalL2(1)*dt + sigma_x(3)*sqrt(dt)*randn(1);
 YL2(j+1) = YL2(j) + VL2(j)*dt;
 VL2(j+1) = VL2(j) - beta(3)*(YL2(j)-y_starL(XL2(j)))*dt - gamma(3)*VL2(j)*dt + centripetalL2(2)*dt + sigma_y(3)*sqrt(dt)*randn(1);
 
 
         XR2(j+1) = XR2(j) + UR2(j)*dt(1);
        UR2(j+1) = UR2(j) - alpha(4)*UR2(j)*(UR2(j)-u1)*(UR2(j)-u2)*dt + centripetalR2(1)*dt + sigma_x(4)*sqrt(dt)*randn(1);
        YR2(j+1) = YR2(j) + VR2(j)*dt;
        VR2(j+1) = VR2(j) - beta(4)*(YR2(j)-y_starR(XR2(j)))*dt - gamma(4)*VR2(j)*dt + centripetalR2(2)*dt + sigma_y(4)*sqrt(dt)*randn(1);
    
        XR3(j+1) = NaN;
        YR3(j+1) = NaN;
        UR3(j+1) = NaN;
        VR3(j+1) = NaN;
    
        XL3(j+1) = NaN;
        YL3(j+1) = NaN;
        UL3(j+1) = NaN;
        VL3(j+1) = NaN;
    
    end 
 j=numel(XL1)
 
          norL1 = [dy_starL(XL1(j)) , -1]./sqrt(1+dy_starL(XL1(j)).^2);
        rhoL1 = abs((1+dy_starL(XL1(j))^2).^(1.5)./ddy_starL);
        centripetalL1 = (UL1(j).^2+VL1(j).^2)/rhoL1*norL1;
        
        norR1 = [dy_starR(XR1(j)) , -1]./sqrt(1+dy_starR(XR1(j)).^2);
        rhoR1 = abs((1+dy_starR(XR1(j))^2).^(1.5)./ddy_starR);
        centripetalR1 = (UR1(j).^2+VR1(j).^2)/rhoR1*norR1;
        
        norL2 = [dy_starL(XL2(j)) , -1]./sqrt(1+dy_starL(XL2(j)).^2);
        rhoL2 = abs((1+dy_starL(XL2(j))^2).^(1.5)./ddy_starL);
        centripetalL2 = (UL2(j).^2+VL2(j).^2)/rhoL2*norL2;
        
        norR2 = [dy_starR(XR2(j)) , -1]./sqrt(1+dy_starR(XR2(j)).^2);
       rhoR2 = abs((1+dy_starR(XR2(j))^2).^(1.5)./ddy_starR);
       centripetalR2 = (UR2(j).^2+VR2(j).^2)/rhoR2*norR2;
 
     XL1(j+1) = XL1(j) + UL1(j)*dt(1);
     UL1(j+1) = UL1(j) - alpha(1)*UL1(j)*(UL1(j)-u1)*(UL1(j)-u2)*dt + centripetalL1(1)*dt + sigma_x(1)*sqrt(dt)*randn(1);
     YL1(j+1) = YL1(j) + VL1(j)*dt(1);
     VL1(j+1) = VL1(j) - beta(1)*(YL1(j)-y_starL(XL1(j)))*dt - gamma(1)*VL1(j)*dt + centripetalL1(2)*dt + sigma_y(1)*sqrt(dt)*randn(1);
     
        XR1(j+1) = XR1(j) + UR1(j)*dt(1);
        UR1(j+1) = UR1(j) - alpha(2)*UR1(j)*(UR1(j)-u1)*(UR1(j)-u2)*dt(1) + centripetalR1(1)*dt(1) + sigma_x(2)*sqrt(dt(1))*randn(1);
        YR1(j+1) = YR1(j) + VR1(j)*dt(1);
        VR1(j+1) = VR1(j) - beta(2)*(YR1(j)-y_starR(XR1(j)))*dt(1) - gamma(2)*VR1(j)*dt(1) + centripetalR1(2)*dt(1) + sigma_y(2)*sqrt(dt(1))*randn(1);
      
     XL2(j+1)= XL2(j) + UL2(j)*dt;
    UL2(j+1) = UL2(j) - alpha(3)*UL2(j)*(UL2(j)-u1)*(UL2(j)-u2)*dt + centripetalL2(1)*dt + sigma_x(3)*sqrt(dt)*randn(1);
    YL2(j+1) = YL2(j) + VL2(j)*dt;
    VL2(j+1) =  VL2(j) - beta(3)*(YL2(j)-y_starL(XL2(j)))*dt - gamma(3)*VL2(j)*dt + centripetalL2(2)*dt + sigma_y(3)*sqrt(dt)*randn(1);
 
         XR2(j+1) = XR2(j) + UR2(j)*dt(1);
        UR2(j+1) = UR2(j) - alpha(4)*UR2(j)*(UR2(j)-u1)*(UR2(j)-u2)*dt + centripetalR2(1)*dt + sigma_x(4)*sqrt(dt)*randn(1);
        YR2(j+1) = YR2(j) + VR2(j)*dt;
        VR2(j+1) = VR2(j) - beta(4)*(YR2(j)-y_starR(XR2(j)))*dt - gamma(4)*VR2(j)*dt + centripetalR2(2)*dt + sigma_y(4)*sqrt(dt)*randn(1);
    
        XL3(j+1) = 10.8;
        YL3(j+1) = y0(5);
        UL3(j+1) = u0(5);
        VL3(j+1) = v0(5);
    
        XR3(j+1) = NaN;
        YR3(j+1) = NaN;
        UR3(j+1) = NaN;
        VR3(j+1) = NaN;

        
         stepNumber5 =  round(abs((t0(5)-t0(6))/dt));
    
     
     for j=numel(XL1):numel(XL1)+stepNumber5-1
 
         norL1 = [dy_starL(XL1(j)) , -1]./sqrt(1+dy_starL(XL1(j)).^2);
        rhoL1 = abs((1+dy_starL(XL1(j))^2).^(1.5)./ddy_starL);
        centripetalL1 = (UL1(j).^2+VL1(j).^2)/rhoL1*norL1;
        
        norR1 = [dy_starR(XR1(j)) , -1]./sqrt(1+dy_starR(XR1(j)).^2);
        rhoR1 = abs((1+dy_starR(XR1(j))^2).^(1.5)./ddy_starR);
        centripetalR1 = (UR1(j).^2+VR1(j).^2)/rhoR1*norR1;
        
        norL2 = [dy_starL(XL2(j)) , -1]./sqrt(1+dy_starL(XL2(j)).^2);
        rhoL2 = abs((1+dy_starL(XL2(j))^2).^(1.5)./ddy_starL);
        centripetalL2 = (UL2(j).^2+VL2(j).^2)/rhoL2*norL2;
        
        norR2 = [dy_starR(XR2(j)) , -1]./sqrt(1+dy_starR(XR2(j)).^2);
       rhoR2 = abs((1+dy_starR(XR2(j))^2).^(1.5)./ddy_starR);
       centripetalR2 = (UR2(j).^2+VR2(j).^2)/rhoR2*norR2;
       
         norL3 = [dy_starL(XL3(j)) , -1]./sqrt(1+dy_starL(XL3(j)).^2);
       rhoL3 = abs((1+dy_starL(XL3(j))^2).^(1.5)./ddy_starL);
       centripetalL3 = (UL3(j).^2+VL3(j).^2)/rhoL3*norL3;
       
 XL1(j+1) = XL1(j) + UL1(j)*dt(1);
 UL1(j+1) = UL1(j) - alpha(1)*UL1(j)*(UL1(j)-u1)*(UL1(j)-u2)*dt + centripetalL1(1)*dt + sigma_x(1)*sqrt(dt)*randn(1);
 YL1(j+1) = YL1(j) + VL1(j)*dt(1);
 VL1(j+1) = VL1(j) - beta(1)*(YL1(j)-y_starL(XL1(j)))*dt - gamma(1)*VL1(j)*dt + centripetalL1(2)*dt + sigma_y(1)*sqrt(dt)*randn(1);
     
  XR1(j+1) = XR1(j) + UR1(j)*dt(1);
  UR1(j+1) = UR1(j) - alpha(2)*UR1(j)*(UR1(j)-u1)*(UR1(j)-u2)*dt(1) + centripetalR1(1)*dt(1) + sigma_x(2)*sqrt(dt(1))*randn(1);
  YR1(j+1) = YR1(j) + VR1(j)*dt(1);
  VR1(j+1) = VR1(j) - beta(2)*(YR1(j)-y_starR(XR1(j)))*dt - gamma(2)*VR1(j)*dt(1) + centripetalR1(2)*dt(1) + sigma_y(2)*sqrt(dt(1))*randn(1);
 
  XL2(j+1) = XL2(j) + UL2(j)*dt;
 UL2(j+1) = UL2(j) - alpha(3)*UL2(j)*(UL2(j)-u1)*(UL2(j)-u2)*dt + centripetalL2(1)*dt + sigma_x(3)*sqrt(dt)*randn(1);
 YL2(j+1) = YL2(j) + VL2(j)*dt;
 VL2(j+1) = VL2(j) - beta(3)*(YL2(j)-y_starL(XL2(j)))*dt - gamma(3)*VL2(j)*dt + centripetalL2(2)*dt + sigma_y(3)*sqrt(dt)*randn(1);
 
 
         XR2(j+1) = XR2(j) + UR2(j)*dt(1);
        UR2(j+1) = UR2(j) - alpha(4)*UR2(j)*(UR2(j)-u1)*(UR2(j)-u2)*dt + centripetalR2(1)*dt + sigma_x(4)*sqrt(dt)*randn(1);
        YR2(j+1) = YR2(j) + VR2(j)*dt;
        VR2(j+1) = VR2(j) - beta(4)*(YR2(j)-y_starR(XR2(j)))*dt - gamma(4)*VR2(j)*dt + centripetalR2(2)*dt + sigma_y(4)*sqrt(dt)*randn(1);
    
         XL3(j+1) = XL3(j) + UL3(j)*dt;
        UL3(j+1) = UL3(j) - alpha(5)*UL3(j)*(UL3(j)-u1)*(UL3(j)-u2)*dt + centripetalL3(1)*dt + sigma_x(5)*sqrt(dt)*randn(1);
        YL3(j+1) = YL3(j) + VL3(j)*dt;
        VL3(j+1) = VL3(j) - beta(5)*(YL3(j)-y_starL(XL3(j)))*dt - gamma(5)*VL3(j)*dt + centripetalL3(2)*dt + sigma_y(5)*sqrt(dt)*randn(1);
    
        XR3(j+1) = NaN;
        YR3(j+1) = NaN;
        UR3(j+1) = NaN;
        VR3(j+1) = NaN;
    
    end 
 j=numel(XL1)

          norL1 = [dy_starL(XL1(j)) , -1]./sqrt(1+dy_starL(XL1(j)).^2);
        rhoL1 = abs((1+dy_starL(XL1(j))^2).^(1.5)./ddy_starL);
        centripetalL1 = (UL1(j).^2+VL1(j).^2)/rhoL1*norL1;
        
        norR1 = [dy_starR(XR1(j)) , -1]./sqrt(1+dy_starR(XR1(j)).^2);
        rhoR1 = abs((1+dy_starR(XR1(j))^2).^(1.5)./ddy_starR);
        centripetalR1 = (UR1(j).^2+VR1(j).^2)/rhoR1*norR1;
        
        norL2 = [dy_starL(XL2(j)) , -1]./sqrt(1+dy_starL(XL2(j)).^2);
        rhoL2 = abs((1+dy_starL(XL2(j))^2).^(1.5)./ddy_starL);
        centripetalL2 = (UL2(j).^2+VL2(j).^2)/rhoL2*norL2;
        
        norR2 = [dy_starR(XR2(j)) , -1]./sqrt(1+dy_starR(XR2(j)).^2);
       rhoR2 = abs((1+dy_starR(XR2(j))^2).^(1.5)./ddy_starR);
       centripetalR2 = (UR2(j).^2+VR2(j).^2)/rhoR2*norR2;
       
        norL3 = [dy_starL(XL3(j)) , -1]./sqrt(1+dy_starL(XL3(j)).^2);
       rhoL3 = abs((1+dy_starL(XL3(j))^2).^(1.5)./ddy_starL);
       centripetalL3 = (UL3(j).^2+VL3(j).^2)/rhoL3*norL3;
       
 
     XL1(j+1) = XL1(j) + UL1(j)*dt(1);
     UL1(j+1) = UL1(j) - alpha(1)*UL1(j)*(UL1(j)-u1)*(UL1(j)-u2)*dt + centripetalL1(1)*dt + sigma_x(1)*sqrt(dt)*randn(1);
     YL1(j+1) = YL1(j) + VL1(j)*dt(1);
     VL1(j+1) = VL1(j) - beta(1)*(YL1(j)-y_starL(XL1(j)))*dt - gamma(1)*VL1(j)*dt + centripetalL1(2)*dt + sigma_y(1)*sqrt(dt)*randn(1);
     
        XR1(j+1) = XR1(j) + UR1(j)*dt(1);
        UR1(j+1) = UR1(j) - alpha(2)*UR1(j)*(UR1(j)-u1)*(UR1(j)-u2)*dt(1) + centripetalR1(1)*dt(1) + sigma_x(2)*sqrt(dt(1))*randn(1);
        YR1(j+1) = YR1(j) + VR1(j)*dt(1);
        VR1(j+1) = VR1(j) - beta(2)*(YR1(j)-y_starR(XR1(j)))*dt(1) - gamma(2)*VR1(j)*dt(1) + centripetalR1(2)*dt(1) + sigma_y(2)*sqrt(dt(1))*randn(1);
      
     XL2(j+1)= XL2(j) + UL2(j)*dt;
    UL2(j+1) = UL2(j) - alpha(3)*UL2(j)*(UL2(j)-u1)*(UL2(j)-u2)*dt + centripetalL2(1)*dt + sigma_x(3)*sqrt(dt)*randn(1);
    YL2(j+1) = YL2(j) + VL2(j)*dt;
    VL2(j+1) =  VL2(j) - beta(3)*(YL2(j)-y_starL(XL2(j)))*dt - gamma(3)*VL2(j)*dt + centripetalL2(2)*dt + sigma_y(3)*sqrt(dt)*randn(1);
 
         XR2(j+1) = XR2(j) + UR2(j)*dt(1);
        UR2(j+1) = UR2(j) - alpha(4)*UR2(j)*(UR2(j)-u1)*(UR2(j)-u2)*dt + centripetalR2(1)*dt + sigma_x(4)*sqrt(dt)*randn(1);
        YR2(j+1) = YR2(j) + VR2(j)*dt;
        VR2(j+1) = VR2(j) - beta(4)*(YR2(j)-y_starR(XR2(j)))*dt - gamma(4)*VR2(j)*dt + centripetalR2(2)*dt + sigma_y(4)*sqrt(dt)*randn(1);
        
        XL3(j+1) = XL3(j) + UL3(j)*dt;
        UL3(j+1) = UL3(j) - alpha(5)*UL3(j)*(UL3(j)-u1)*(UL3(j)-u2)*dt + centripetalL3(1)*dt + sigma_x(5)*sqrt(dt)*randn(1);
        YL3(j+1) = YL3(j) + VL3(j)*dt;
        VL3(j+1) = VL3(j) - beta(5)*(YL3(j)-y_starL(XL3(j)))*dt - gamma(5)*VL3(j)*dt + centripetalL3(2)*dt + sigma_y(5)*sqrt(dt)*randn(1);
    

    
        XR3(j+1) = 0;
        YR3(j+1) = y0(6);
        UR3(j+1) = u0(6);
        VR3(j+1) = v0(6);

 j=numel(XL1);
 
   while((XL1(j)<=10.8&&XL1(j)>=0)||(XR1(j)<=10.8&&XR1(j)>=0))%%||(XL2(j)<=10.8&&XL2(j)>=0)||(XR2(j)<=10.8&&XR2(j)>=0)||(XL3(j)<=10.8&&XL3(j)>=0)||(XR3(j)<=10.8&&XR3(j)>=0))
      
      %these variables are for the attraction felt by each person to store
      %1 and store 2. It changes as the conditions in the loop are met. If
      %none of the conditions are met then the attraction remains zero and
      %as they walk they will not be drawn to store 1 or store 2 depending
      %on what conditions were met.
       attractionStore1UL1=0;
       attractionStore1VL1=0;
       attractionStore1UL2=0;
       attractionStore1VL2=0;
       attractionStore1UL3=0;
       attractionStore1VL3=0;
       
       attractionStore1UR1=0;
       attractionStore1VR1=0;
       attractionStore1UR2=0;
       attractionStore1VR2=0;
       attractionStore1UR3=0;
       attractionStore1VR3=0;
       
       attractionStore2UL1=0;
       attractionStore2VL1=0;
       attractionStore2UL2=0;
       attractionStore2VL2=0;
       attractionStore2UL3=0;
       attractionStore2VL3=0;
       
       attractionStore2UR1=0;
       attractionStore2VR1=0;
       attractionStore2UR2=0;
       attractionStore2VR2=0;
       attractionStore2UR3=0;
       attractionStore2VR3=0;
       
       
       
       
       avoidanceUL1=0; %these are going to be the avoidance terms added for each pedestrian for the U term
       avoidanceUR1=0; %they are set to zero because ifnone of the conditions below are met then they will walk
       avoidanceUL2=0; %without avoidance force. This is for the avoidance of other pedestrains.
       avoidanceUR2=0;
       avoidanceUL3=0;
       avoidanceUR3=0;
       
       
       avoidanceVL1=0; %these are going to be the avoidance terms added for each pedestrian for the V term
       avoidanceVR1=0;
       avoidanceVL2=0;
       avoidanceVR2=0;
       avoidanceVL3=0;
       avoidanceVR3=0;
       
     
       
       
       
       %calculate the centripetal force for each person
       norL1 = [dy_starL(XL1(j)) , -1]./sqrt(1+dy_starL(XL1(j)).^2); 
        rhoL1 = abs((1+dy_starL(XL1(j))^2).^(1.5)./ddy_starL);
        centripetalL1 = (UL1(j).^2+VL1(j).^2)/rhoL1*norL1;
        
        norR1 = [dy_starR(XR1(j)) , -1]./sqrt(1+dy_starR(XR1(j)).^2);
        rhoR1 = abs((1+dy_starR(XR1(j))^2).^(1.5)./ddy_starR);
        centripetalR1 = (UR1(j).^2+VR1(j).^2)/rhoR1*norR1;
        
        norL2 = [dy_starL(XL2(j)) , -1]./sqrt(1+dy_starL(XL2(j)).^2);
        rhoL2 = abs((1+dy_starL(XL2(j))^2).^(1.5)./ddy_starL);
        centripetalL2 = (UL2(j).^2+VL2(j).^2)/rhoL2*norL2;
        
            
        norR2 = [dy_starR(XR2(j)) , -1]./sqrt(1+dy_starR(XR2(j)).^2);
       rhoR2 = abs((1+dy_starR(XR2(j))^2).^(1.5)./ddy_starR);
       centripetalR2 = (UR2(j).^2+VR2(j).^2)/rhoR2*norR2;
       
        norL3 = [dy_starL(XL3(j)) , -1]./sqrt(1+dy_starL(XL3(j)).^2);
       rhoL3 = abs((1+dy_starL(XL3(j))^2).^(1.5)./ddy_starL);
       centripetalL3 = (UL3(j).^2+VL3(j).^2)/rhoL3*norL3;
       
        norR3 = [dy_starR(XR3(j)) , -1]./sqrt(1+dy_starR(XR3(j)).^2);
       rhoR3 = abs((1+dy_starR(XR3(j))^2).^(1.5)./ddy_starR);
       centripetalR3 = (UR3(j).^2+VR3(j).^2)/rhoR3*norR3;
       
      
        
       %this section is calculating the delta term that is an integral. The
       %person who is the "me" is listed first. The location they are
       %trying to go to or avoid is listed second.
       
       
       deltaL1Store1Position=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*Store1Position(1)+b(1))*(2*a(1)*Store1Position(1)+b(1)))*(2*a(1)*Store1Position(1)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*Store1Position(1)+b(1))*(2*a(1)*Store1Position(1)+b(1)))+(2*a(1)*Store1Position(1)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))*(2*a(1)*XL1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))+(2*a(1)*XL1(j)+b(1)))));
       
       deltaL2Store1Position=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*Store1Position(1)+b(1))*(2*a(1)*Store1Position(1)+b(1)))*(2*a(1)*Store1Position(1)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*Store1Position(1)+b(1))*(2*a(1)*Store1Position(1)+b(1)))+(2*a(1)*Store1Position(1)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))*(2*a(1)*XL2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))+(2*a(1)*XL2(j)+b(1)))));
       
       deltaL3Store1Position=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*Store1Position(1)+b(1))*(2*a(1)*Store1Position(1)+b(1)))*(2*a(1)*Store1Position(1)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*Store1Position(1)+b(1))*(2*a(1)*Store1Position(1)+b(1)))+(2*a(1)*Store1Position(1)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))*(2*a(1)*XL3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))+(2*a(1)*XL3(j)+b(1)))));
       
       deltaR1Store1Position=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*Store1Position(1)+b(1))*(2*a(1)*Store1Position(1)+b(1)))*(2*a(1)*Store1Position(1)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*Store1Position(1)+b(1))*(2*a(1)*Store1Position(1)+b(1)))+(2*a(1)*Store1Position(1)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))*(2*a(1)*XR1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))+(2*a(1)*XR1(j)+b(1)))));
       
       deltaR2Store1Position=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*Store1Position(1)+b(1))*(2*a(1)*Store1Position(1)+b(1)))*(2*a(1)*Store1Position(1)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*Store1Position(1)+b(1))*(2*a(1)*Store1Position(1)+b(1)))+(2*a(1)*Store1Position(1)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))*(2*a(1)*XR2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))+(2*a(1)*XR2(j)+b(1)))));
       
       deltaR3Store1Position=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*Store1Position(1)+b(1))*(2*a(1)*Store1Position(1)+b(1)))*(2*a(1)*Store1Position(1)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*Store1Position(1)+b(1))*(2*a(1)*Store1Position(1)+b(1)))+(2*a(1)*Store1Position(1)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))*(2*a(1)*XR3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))+(2*a(1)*XR3(j)+b(1)))));
       
       
       deltaL1Store2Position=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*Store2Position(1)+b(1))*(2*a(1)*Store2Position(1)+b(1)))*(2*a(1)*Store2Position(1)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*Store2Position(1)+b(1))*(2*a(1)*Store2Position(1)+b(1)))+(2*a(1)*Store2Position(1)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))*(2*a(1)*XL1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))+(2*a(1)*XL1(j)+b(1)))));
       
       deltaL2Store2Position=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*Store2Position(1)+b(1))*(2*a(1)*Store2Position(1)+b(1)))*(2*a(1)*Store2Position(1)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*Store2Position(1)+b(1))*(2*a(1)*Store2Position(1)+b(1)))+(2*a(1)*Store2Position(1)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))*(2*a(1)*XL2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))+(2*a(1)*XL2(j)+b(1)))));
       
       deltaL3Store2Position=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*Store2Position(1)+b(1))*(2*a(1)*Store2Position(1)+b(1)))*(2*a(1)*Store2Position(1)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*Store2Position(1)+b(1))*(2*a(1)*Store2Position(1)+b(1)))+(2*a(1)*Store2Position(1)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))*(2*a(1)*XL3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))+(2*a(1)*XL3(j)+b(1)))));
       
       deltaR1Store2Position=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*Store2Position(1)+b(1))*(2*a(1)*Store2Position(1)+b(1)))*(2*a(1)*Store2Position(1)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*Store2Position(1)+b(1))*(2*a(1)*Store2Position(1)+b(1)))+(2*a(1)*Store2Position(1)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))*(2*a(1)*XR1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))+(2*a(1)*XR1(j)+b(1)))));
       
       deltaR2Store2Position=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*Store2Position(1)+b(1))*(2*a(1)*Store2Position(1)+b(1)))*(2*a(1)*Store2Position(1)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*Store2Position(1)+b(1))*(2*a(1)*Store2Position(1)+b(1)))+(2*a(1)*Store2Position(1)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))*(2*a(1)*XR2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))+(2*a(1)*XR2(j)+b(1)))));
       
       deltaR3Store2Position=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*Store2Position(1)+b(1))*(2*a(1)*Store2Position(1)+b(1)))*(2*a(1)*Store2Position(1)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*Store2Position(1)+b(1))*(2*a(1)*Store2Position(1)+b(1)))+(2*a(1)*Store2Position(1)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))*(2*a(1)*XR3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))+(2*a(1)*XR3(j)+b(1)))));
       
       
       
         %this section is calculating the delta term that is an integral. The
       %person who is the "me" is listed first. the "other" person is
       %listed second. For example deltaL1R1 is the calculations for Person
       %L1's delta Xi(delta squiggle, i think it's called Xi) when they are
       %facing R1.
       
       deltaL1R1=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))*(2*a(1)*XR1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))+(2*a(1)*XR1(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))*(2*a(1)*XL1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))+(2*a(1)*XL1(j)+b(1)))));
       deltaR1L1=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XL1(j)+b(2))*(2*a(2)*XL1(j)+b(2)))*(2*a(2)*XL1(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL1(j)+b(2))*(2*a(2)*XL1(j)+b(2)))+(2*a(2)*XL1(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XR1(j)+b(2))*(2*a(2)*XR1(j)+b(2)))*(2*a(2)*XR1(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR1(j)+b(2))*(2*a(2)*XR1(j)+b(2)))+(2*a(2)*XR1(j)+b(2)))));
       
       deltaL1R2=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))*(2*a(1)*XR2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))+(2*a(1)*XR2(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))*(2*a(1)*XL1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))+(2*a(1)*XL1(j)+b(1)))));
       deltaR2L1=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XL1(j)+b(2))*(2*a(2)*XL1(j)+b(2)))*(2*a(2)*XL1(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL1(j)+b(2))*(2*a(2)*XL1(j)+b(2)))+(2*a(2)*XL1(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XR2(j)+b(2))*(2*a(2)*XR2(j)+b(2)))*(2*a(2)*XR2(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR2(j)+b(2))*(2*a(2)*XR2(j)+b(2)))+(2*a(2)*XR2(j)+b(2)))));
       
       deltaL1R3=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))*(2*a(1)*XR3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))+(2*a(1)*XR3(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))*(2*a(1)*XL1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))+(2*a(1)*XL1(j)+b(1)))));
       deltaR3L1=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XL1(j)+b(2))*(2*a(2)*XL1(j)+b(2)))*(2*a(2)*XL1(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL1(j)+b(2))*(2*a(2)*XL1(j)+b(2)))+(2*a(2)*XL1(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XR3(j)+b(2))*(2*a(2)*XR3(j)+b(2)))*(2*a(2)*XR3(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR3(j)+b(2))*(2*a(2)*XR3(j)+b(2)))+(2*a(2)*XR3(j)+b(2)))));
       
       deltaL2R1=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))*(2*a(1)*XR1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))+(2*a(1)*XR1(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))*(2*a(1)*XL2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))+(2*a(1)*XL2(j)+b(1)))));
       deltaR1L2=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XL2(j)+b(2))*(2*a(2)*XL2(j)+b(2)))*(2*a(2)*XL2(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL2(j)+b(2))*(2*a(2)*XL2(j)+b(2)))+(2*a(2)*XL2(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XR1(j)+b(2))*(2*a(2)*XR1(j)+b(2)))*(2*a(2)*XR1(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR1(j)+b(2))*(2*a(2)*XR1(j)+b(2)))+(2*a(2)*XR1(j)+b(2)))));
      
       deltaL2R2=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))*(2*a(1)*XR2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))+(2*a(1)*XR2(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))*(2*a(1)*XL2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))+(2*a(1)*XL2(j)+b(1)))));
       deltaR2L2=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XL2(j)+b(2))*(2*a(2)*XL2(j)+b(2)))*(2*a(2)*XL2(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL2(j)+b(2))*(2*a(2)*XL2(j)+b(2)))+(2*a(2)*XL2(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XR2(j)+b(2))*(2*a(2)*XR2(j)+b(2)))*(2*a(2)*XR2(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR2(j)+b(2))*(2*a(2)*XR2(j)+b(2)))+(2*a(2)*XR2(j)+b(2)))));
       
       deltaL2R3=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))*(2*a(1)*XR3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))+(2*a(1)*XR3(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))*(2*a(1)*XL2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))+(2*a(1)*XL2(j)+b(1)))));
       deltaR3L2=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XL2(j)+b(2))*(2*a(2)*XL2(j)+b(2)))*(2*a(2)*XL2(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL2(j)+b(2))*(2*a(2)*XL2(j)+b(2)))+(2*a(2)*XL2(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XR3(j)+b(2))*(2*a(2)*XR3(j)+b(2)))*(2*a(2)*XR3(j)+b(2))+(1/(4*a(2)))*log(abs(sqrt(1+(2*a(2)*XR3(j)+b(2))*(2*a(2)*XR3(j)+b(2)))+(2*a(2)*XR3(j)+b(2))))));
        
       deltaL3R1=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))*(2*a(1)*XR1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))+(2*a(1)*XR1(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))*(2*a(1)*XL3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))+(2*a(1)*XL3(j)+b(1)))));
       deltaR1L3=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XL3(j)+b(2))*(2*a(2)*XL3(j)+b(2)))*(2*a(2)*XL3(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL3(j)+b(2))*(2*a(2)*XL3(j)+b(2)))+(2*a(2)*XL3(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XR1(j)+b(2))*(2*a(2)*XR1(j)+b(2)))*(2*a(2)*XR1(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR1(j)+b(2))*(2*a(2)*XR1(j)+b(2)))+(2*a(2)*XR1(j)+b(2)))));
      
       deltaL3R2=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))*(2*a(1)*XR2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))+(2*a(1)*XR2(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))*(2*a(1)*XL3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))+(2*a(1)*XL3(j)+b(1)))));
       deltaR2L3=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XL3(j)+b(2))*(2*a(2)*XL3(j)+b(2)))*(2*a(2)*XL3(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL3(j)+b(2))*(2*a(2)*XL3(j)+b(2)))+(2*a(2)*XL3(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XR2(j)+b(2))*(2*a(2)*XR2(j)+b(2)))*(2*a(2)*XR2(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR2(j)+b(2))*(2*a(2)*XR2(j)+b(2)))+(2*a(2)*XR2(j)+b(2)))));
       
       deltaL3R3=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))*(2*a(1)*XR3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))+(2*a(1)*XR3(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))*(2*a(1)*XL3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))+(2*a(1)*XL3(j)+b(1)))));
       deltaR3L3=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XL3(j)+b(2))*(2*a(2)*XL3(j)+b(2)))*(2*a(2)*XL3(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL3(j)+b(2))*(2*a(2)*XL3(j)+b(2)))+(2*a(2)*XL3(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XR3(j)+b(2))*(2*a(2)*XR3(j)+b(2)))*(2*a(2)*XR3(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR3(j)+b(2))*(2*a(2)*XR3(j)+b(2)))+(2*a(2)*XR3(j)+b(2)))));
        
       
       deltaL1L2=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))*(2*a(1)*XL2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))+(2*a(1)*XL2(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))*(2*a(1)*XL1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))+(2*a(1)*XL1(j)+b(1)))));
       deltaL2L1=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XL1(j)+b(2))*(2*a(2)*XL1(j)+b(2)))*(2*a(2)*XL1(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL1(j)+b(2))*(2*a(2)*XL1(j)+b(2)))+(2*a(2)*XL1(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XL2(j)+b(2))*(2*a(2)*XL2(j)+b(2)))*(2*a(2)*XL2(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL2(j)+b(2))*(2*a(2)*XL2(j)+b(2)))+(2*a(2)*XL2(j)+b(2)))));
       
       deltaL1L3=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))*(2*a(1)*XL3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))+(2*a(1)*XL3(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))*(2*a(1)*XL1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL1(j)+b(1))*(2*a(1)*XL1(j)+b(1)))+(2*a(1)*XL1(j)+b(1)))));
       deltaL3L1=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XL1(j)+b(2))*(2*a(2)*XL1(j)+b(2)))*(2*a(2)*XL1(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL1(j)+b(2))*(2*a(2)*XL1(j)+b(2)))+(2*a(2)*XL1(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XL3(j)+b(2))*(2*a(2)*XL3(j)+b(2)))*(2*a(2)*XL3(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL3(j)+b(2))*(2*a(2)*XL3(j)+b(2)))+(2*a(2)*XL3(j)+b(2)))));
       
       deltaL2L3=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))*(2*a(1)*XL3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL3(j)+b(1))*(2*a(1)*XL3(j)+b(1)))+(2*a(1)*XL3(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))*(2*a(1)*XL2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XL2(j)+b(1))*(2*a(1)*XL2(j)+b(1)))+(2*a(1)*XL2(j)+b(1)))));
       deltaL3L2=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XL2(j)+b(2))*(2*a(2)*XL2(j)+b(2)))*(2*a(2)*XL2(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL2(j)+b(2))*(2*a(2)*XL2(j)+b(2)))+(2*a(2)*XL2(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XL3(j)+b(2))*(2*a(2)*XL3(j)+b(2)))*(2*a(2)*XL3(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XL3(j)+b(2))*(2*a(2)*XL3(j)+b(2)))+(2*a(2)*XL3(j)+b(2)))));
       
       deltaR1R2=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))*(2*a(1)*XR2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))+(2*a(1)*XR2(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))*(2*a(1)*XR1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))+(2*a(1)*XR1(j)+b(1)))));
       deltaR2R1=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XR1(j)+b(2))*(2*a(2)*XR1(j)+b(2)))*(2*a(2)*XR1(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR1(j)+b(2))*(2*a(2)*XR1(j)+b(2)))+(2*a(2)*XR1(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XR2(j)+b(2))*(2*a(2)*XR2(j)+b(2)))*(2*a(2)*XR2(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR2(j)+b(2))*(2*a(2)*XR2(j)+b(2)))+(2*a(2)*XR2(j)+b(2)))));
      
       deltaR1R3=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))*(2*a(1)*XR3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))+(2*a(1)*XR3(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))*(2*a(1)*XR1(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR1(j)+b(1))*(2*a(1)*XR1(j)+b(1)))+(2*a(1)*XR1(j)+b(1)))));
       deltaR3R1=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XR1(j)+b(2))*(2*a(2)*XR1(j)+b(2)))*(2*a(2)*XR1(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR1(j)+b(2))*(2*a(2)*XR1(j)+b(2)))+(2*a(2)*XR1(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XR3(j)+b(2))*(2*a(2)*XR3(j)+b(2)))*(2*a(2)*XR3(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR3(j)+b(2))*(2*a(2)*XR3(j)+b(2)))+(2*a(2)*XR3(j)+b(2)))));
       
       deltaR2R3=abs((1/(4*a(1)))*sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))*(2*a(1)*XR3(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR3(j)+b(1))*(2*a(1)*XR3(j)+b(1)))+(2*a(1)*XR3(j)+b(1)))-((1/(4*a(1)))*sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))*(2*a(1)*XR2(j)+b(1))+(1/(4*a(1)))*log(sqrt(1+(2*a(1)*XR2(j)+b(1))*(2*a(1)*XR2(j)+b(1)))+(2*a(1)*XR2(j)+b(1)))));
       deltaR3R2=abs((1/(4*a(2)))*sqrt(1+(2*a(2)*XR2(j)+b(2))*(2*a(2)*XR2(j)+b(2)))*(2*a(2)*XR2(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR2(j)+b(2))*(2*a(2)*XR2(j)+b(2)))+(2*a(2)*XR2(j)+b(2)))-((1/(4*a(2)))*sqrt(1+(2*a(2)*XR3(j)+b(2))*(2*a(2)*XR3(j)+b(2)))*(2*a(2)*XR3(j)+b(2))+(1/(4*a(2)))*log(sqrt(1+(2*a(2)*XR3(j)+b(2))*(2*a(2)*XR3(j)+b(2)))+(2*a(2)*XR3(j)+b(2)))));
       
       
       %this is the deltaX,deltaY for the pedestrains and the stores they
       %are trying to go to or avoid.
        AL1Store1Position=[Store1Position(1)-XL1(j),Store1Position(2)-YL1(j)];
        
        AL2Store1Position=[Store1Position(1)-XL2(j),Store1Position(2)-YL2(j)];
        
        AL3Store1Position=[Store1Position(1)-XL3(j),Store1Position(2)-YL3(j)];
        
        AR1Store1Position=[Store1Position(1)-XR1(j),Store1Position(2)-YR1(j)];
        
        AR2Store1Position=[Store1Position(1)-XR2(j),Store1Position(2)-YR2(j)];
        
        AR3Store1Position=[Store1Position(1)-XR3(j),Store1Position(2)-YR3(j)];
        
        
        AL1Store2Position=[Store2Position(1)-XL1(j),Store2Position(2)-YL1(j)];
        
        AL2Store2Position=[Store2Position(1)-XL2(j),Store2Position(2)-YL2(j)];
        
        AL3Store2Position=[Store2Position(1)-XL3(j),Store2Position(2)-YL3(j)];
        
        AR1Store2Position=[Store2Position(1)-XR1(j),Store2Position(2)-YR1(j)];
        
        AR2Store2Position=[Store2Position(1)-XR2(j),Store2Position(2)-YR2(j)];
        
        AR3Store2Position=[Store2Position(1)-XR3(j),Store2Position(2)-YR3(j)];
        
        
        
           %this is the deltaX,deltaY for the pedestrains depending on who is
       %first and second. AL1R1 is saying that the L1 person is "me" and R1
       %is "other"
        AL1R1=[XR1(j)-XL1(j),YR1(j)-YL1(j)];
        AR1L1=[XL1(j)-XR1(j),YL1(j)-YR1(j)];
      
        AL1R2=[XR2(j)-XL1(j),YR2(j)-YL1(j)];
        AR2L1=[XL1(j)-XR2(j),YL1(j)-YR2(j)];
        
        AL1R3=[XR3(j)-XL1(j),YR3(j)-YL1(j)];
        AR3L1=[XL1(j)-XR3(j),YL1(j)-YR3(j)];
      
        AL2R1=[XR1(j)-XL2(j),YR1(j)-YL2(j)];
        AR1L2=[XL2(j)-XR1(j),YL2(j)-YR1(j)];
        
        AL2R2=[XR2(j)-XL2(j),YR2(j)-YL2(j)];
        AR2L2=[XL2(j)-XR2(j),YL2(j)-YR2(j)];
        
        AL2R3=[XR3(j)-XL2(j),YR3(j)-YL2(j)];
        AR3L2=[XL2(j)-XR3(j),YL2(j)-YR3(j)];
        
        AL3R1=[XR1(j)-XL3(j),YR1(j)-YL3(j)];
        AR1L3=[XL3(j)-XR1(j),YL3(j)-YR1(j)];
        
        AL3R2=[XR2(j)-XL3(j),YR2(j)-YL3(j)];
        AR2L3=[XL3(j)-XR2(j),YL3(j)-YR2(j)];
        
        AL3R3=[XR3(j)-XL3(j),YR3(j)-YL3(j)];
        AR3L3=[XL3(j)-XR3(j),YL3(j)-YR3(j)];
        
        AL1L2=[XL2(j)-XL1(j),YL2(j)-YL1(j)];
        AL2L1=[XL1(j)-XL2(j),YL1(j)-YL2(j)];
      
        AL1L3=[XL3(j)-XL1(j),YL3(j)-YL1(j)];
        AL3L1=[XL1(j)-XL3(j),YL1(j)-YL3(j)];
        
        AL2L3=[XL3(j)-XL2(j),YL3(j)-YL2(j)];
        AL3L2=[XL2(j)-XL3(j),YL2(j)-YL3(j)];
        
        AR1R2=[XR2(j)-XR1(j),YR2(j)-YR1(j)];
        AR2R1=[XR1(j)-XR2(j),YR1(j)-YR2(j)];
      
        AR1R3=[XR3(j)-XR1(j),YR3(j)-YR1(j)];
        AR3R1=[XR1(j)-XR3(j),YR1(j)-YR3(j)];
        
        AR2R3=[XR3(j)-XR2(j),YR3(j)-YR2(j)];
        AR3R2=[XR2(j)-XR3(j),YR2(j)-YR3(j)];
        
        %the tangent vectors for each person
        BL1=-[1, dy_starL(XL1(j))]./sqrt(1+(dy_starL(XL1(j)).^2));
        BR1=[1, dy_starR(XR1(j))]./sqrt(1+(dy_starR(XR1(j)).^2));
        BL2=-[1, dy_starL(XL2(j))]./sqrt(1+(dy_starL(XL2(j)).^2));
        BR2=[1, dy_starR(XR2(j))]./sqrt(1+(dy_starR(XR2(j)).^2));
        BL3=-[1, dy_starL(XL3(j))]./sqrt(1+(dy_starL(XL3(j)).^2));
        BR3=[1, dy_starR(XR3(j))]./sqrt(1+(dy_starR(XR3(j)).^2));
       
        delta2L1Store1Position=abs(dot(AL1Store1Position,norL1));
        
        delta2L2Store1Position=abs(dot(AL2Store1Position,norL2));
        
        delta2L3Store1Position=abs(dot(AL3Store1Position,norL3));
        
        delta2R1Store1Position=abs(dot(AR1Store1Position,norR1));
        
        delta2R2Store1Position=abs(dot(AR2Store1Position,norR2));
        
        delta2R3Store1Position=abs(dot(AR3Store1Position,norR3));
        
        
        delta2L1Store2Position=abs(dot(AL1Store2Position,norL1));
        
        delta2L2Store2Position=abs(dot(AL2Store2Position,norL2));
        
        delta2L3Store2Position=abs(dot(AL3Store2Position,norL3));
        
        delta2R1Store2Position=abs(dot(AR1Store2Position,norR1));
        
        delta2R2Store2Position=abs(dot(AR2Store2Position,norR2));
        
        delta2R3Store2Position=abs(dot(AR3Store2Position,norR3));
        
        delta2L1R1=abs(dot(AL1R1,norL1))-offset;
        delta2R1L1=abs(dot(AR1L1,norR1))-offset;
        
        delta2L1R2=abs(dot(AL1R2,norL1))-offset;
        delta2R2L1=abs(dot(AR2L1,norR2))-offset; 
        
        delta2L1R3=abs(dot(AL1R3,norL1))-offset;
        delta2R3L1=abs(dot(AR3L1,norR2))-offset;
        
        delta2L2R1=abs(dot(AL2R1,norL2))-offset;
        delta2R1L2=abs(dot(AR1L2,norR1))-offset; 
        
        delta2L2R2=abs(dot(AL2R2,norL2))-offset;
        delta2R2L2=abs(dot(AR2L2,norR2))-offset;
        
        delta2L2R3=abs(dot(AL2R3,norL2))-offset;
        delta2R3L2=abs(dot(AR3L2,norR3))-offset;
        
        delta2L3R1=abs(dot(AL3R1,norL3))-offset;
        delta2R1L3=abs(dot(AR1L3,norR1))-offset; 
        
        delta2L3R2=abs(dot(AL3R2,norL3))-offset;
        delta2R2L3=abs(dot(AR2L3,norR2))-offset;
        
        delta2L3R3=abs(dot(AL3R3,norL3))-offset;
        delta2R3L3=abs(dot(AR3L3,norR3))-offset;
        
        delta2L1L2=abs(dot(AL1L2,norL1))-offset;
        delta2L2L1=abs(dot(AL2L1,norL2))-offset;
        
        delta2L1L3=abs(dot(AL1L3,norL1))-offset;
        delta2L3L1=abs(dot(AL3L1,norL3))-offset;
        
        delta2L2L3=abs(dot(AL2L3,norL2))-offset;
        delta2L3L2=abs(dot(AL3L2,norL3))-offset;
        
        delta2R1R2=abs(dot(AR1R2,norR1))-offset;
        delta2R2R1=abs(dot(AR2R1,norR2))-offset;
        
        delta2R1R3=abs(dot(AR1R3,norR1))-offset;
        delta2R3R1=abs(dot(AR3R1,norR3))-offset;
        
        delta2R2R3=abs(dot(AR2R3,norR2))-offset;
        delta2R3R2=abs(dot(AR3R2,norR3))-offset;
        
       
     %anytime a person is in the left half of the screen they will feel
     %their attraction or avoidance to store 1.
     if(XL1(j)<5)
         if(time1L1>=75)%unless they already visited the store then they will never want to go back.
             attractionStore1UL1=0;
             attractionStore1VL1=0;
         else
        attractionStore1UL1=attractionStore1UL1+(AttAv1L1*(L1AttractionVar1*kappa*(Store1Position(1)-XL1(j))/sqrt(deltaL1Store1Position*deltaL1Store1Position+delta2L1Store1Position*delta2L1Store1Position)))*exp(-((abs(deltaL1Store1Position)+theta1*abs(delta2L1Store1Position))/theta2))*dt;
        attractionStore1VL1=attractionStore1VL1+(AttAv1L1*(L1AttractionVar1*kappa*(Store1Position(2)-YL1(j))/sqrt(deltaL1Store1Position*deltaL1Store1Position+delta2L1Store1Position*delta2L1Store1Position)))*exp(-((abs(deltaL1Store1Position)+theta1*abs(delta2L1Store1Position))/theta2))*dt;
        
        end
     end
    
        if(XL2(j)<5)
            if(time1L2>=75)
             attractionStore1UL2=0;
             attractionStore1VL2=0;
            else
            
        attractionStore1UL2=attractionStore1UL2+(AttAv1L2*(L2AttractionVar1*kappa*(Store1Position(1)-XL2(j))/sqrt(deltaL2Store1Position*deltaL2Store1Position+delta2L2Store1Position*delta2L2Store1Position)))*exp(-((abs(deltaL2Store1Position)+theta1*abs(delta2L2Store1Position))/theta2))*dt;
        attractionStore1VL2=attractionStore1VL2+(AttAv1L2*(L2AttractionVar1*kappa*(Store1Position(2)-YL2(j))/sqrt(deltaL2Store1Position*deltaL2Store1Position+delta2L2Store1Position*delta2L2Store1Position)))*exp(-((abs(deltaL2Store1Position)+theta1*abs(delta2L2Store1Position))/theta2))*dt;
        
             end 
        end
         
        if(XL3(j)<5)
            
         if(time1L3>=75)
             attractionStore1UL3=0;
             attractionStore1VL3=0;
         else
        attractionStore1UL3=attractionStore1UL3+(AttAv1L3*(L3AttractionVar1*kappa*(Store1Position(1)-XL3(j))/sqrt(deltaL3Store1Position*deltaL3Store1Position+delta2L3Store1Position*delta2L3Store1Position)))*exp(-((abs(deltaL3Store1Position)+theta1*abs(delta2L3Store1Position))/theta2))*dt;
        attractionStore1VL3=attractionStore1VL3+(AttAv1L3*(L3AttractionVar1*kappa*(Store1Position(2)-YL3(j))/sqrt(deltaL3Store1Position*deltaL3Store1Position+delta2L3Store1Position*delta2L3Store1Position)))*exp(-((abs(deltaL3Store1Position)+theta1*abs(delta2L3Store1Position))/theta2))*dt;
         end
        end
  
         if(XR1(j)<5)
                
             if(time1R1>=75)
             attractionStore1UR1=0;
             attractionStore1VR1=0;
            else
        attractionStore1UR1=attractionStore1UR1+(AttAv1R1*(R1AttractionVar1*kappa*(Store1Position(1)-XR1(j))/sqrt(deltaR1Store1Position*deltaR1Store1Position+delta2R1Store1Position*delta2R1Store1Position)))*exp(-((abs(deltaR1Store1Position)+theta1*abs(delta2R1Store1Position))/theta2))*dt;
        attractionStore1VR1=attractionStore1VR1+(AttAv1R1*(R1AttractionVar1*kappa*(Store1Position(2)-YR1(j))/sqrt(deltaR1Store1Position*deltaR1Store1Position+delta2R1Store1Position*delta2R1Store1Position)))*exp(-((abs(deltaR1Store1Position)+theta1*abs(delta2R1Store1Position))/theta2))*dt;
             end
        end
     
         if(XR2(j)<5)
             
             if(time1R2>=75)
             attractionStore1UL1=0;
             attractionStore1VL1=0;
             else
        attractionStore1UR2=attractionStore1UR2+(AttAv1R2*(R2AttractionVar1*kappa*(Store1Position(1)-XR2(j))/sqrt(deltaR2Store1Position*deltaR2Store1Position+delta2R2Store1Position*delta2R2Store1Position)))*exp(-((abs(deltaR2Store1Position)+theta1*abs(delta2R2Store1Position))/theta2))*dt;
        attractionStore1VR2=attractionStore1VR2+(AttAv1R2*(R2AttractionVar1*kappa*(Store1Position(2)-YR2(j))/sqrt(deltaR2Store1Position*deltaR2Store1Position+delta2R2Store1Position*delta2R2Store1Position)))*exp(-((abs(deltaR2Store1Position)+theta1*abs(delta2R1Store1Position))/theta2))*dt;
             end
        end
    
        if(XR3(j)<5)
            if(time1R3>=75)
             attractionStore1UR3=0;
             attractionStore1VR3=0;
            else
        attractionStore1UR3=attractionStore1UR3+(AttAv1R3*(R3AttractionVar1*kappa*(Store1Position(1)-XR3(j))/sqrt(deltaR3Store1Position*deltaR3Store1Position+delta2R3Store1Position*delta2R3Store1Position)))*exp(-((abs(deltaR3Store1Position)+theta1*abs(delta2R3Store1Position))/theta2))*dt;
        attractionStore1VR3=attractionStore1VR3+(AttAv1R3*(R3AttractionVar1*kappa*(Store1Position(2)-YR3(j))/sqrt(deltaR3Store1Position*deltaR3Store1Position+delta2R3Store1Position*delta2R3Store1Position)))*exp(-((abs(deltaR3Store1Position)+theta1*abs(delta2R3Store1Position))/theta2))*dt;
            end
        end
         
        %if people are in the right half of the screen then their
        %attraction or avoidance is applied for store 2.
        if(XL1(j)>5)
         if(time2L1>=75)%unless they already visited the store. if they went already they wont go back.
             attractionStore2UL1=0;
             attractionStore2VL1=0;
         else
        attractionStore2UL1=attractionStore2UL1+(AttAv2L1*(L1AttractionVar2*kappa*(Store2Position(1)-XL1(j))/sqrt(deltaL1Store2Position*deltaL1Store2Position+delta2L1Store2Position*delta2L1Store2Position)))*exp(-((abs(deltaL1Store2Position)+theta1*abs(delta2L1Store2Position))/theta2))*dt;
        attractionStore2VL1=attractionStore2VL1+(AttAv2L1*(L1AttractionVar2*kappa*(Store2Position(2)-YL1(j))/sqrt(deltaL1Store2Position*deltaL1Store2Position+delta2L1Store2Position*delta2L1Store2Position)))*exp(-((abs(deltaL1Store2Position)+theta1*abs(delta2L1Store2Position))/theta2))*dt;
        
        end
     end
    
        if(XL2(j)>5)
            if(time2L2>=75)
             attractionStore2UL2=0;
             attractionStore2VL2=0;
            else
            
        attractionStore2UL2=attractionStore2UL2+(AttAv2L2*(L2AttractionVar2*kappa*(Store2Position(1)-XL2(j))/sqrt(deltaL2Store2Position*deltaL2Store2Position+delta2L2Store2Position*delta2L2Store2Position)))*exp(-((abs(deltaL2Store2Position)+theta1*abs(delta2L2Store2Position))/theta2))*dt;
        attractionStore2VL2=attractionStore2VL2+(AttAv2L2*(L2AttractionVar2*kappa*(Store2Position(2)-YL2(j))/sqrt(deltaL2Store2Position*deltaL2Store2Position+delta2L2Store2Position*delta2L2Store2Position)))*exp(-((abs(deltaL2Store2Position)+theta1*abs(delta2L2Store2Position))/theta2))*dt;
        
             end 
        end
   
        if(XL3(j)>5)
            
         if(time2L3>=75)
             attractionStore2UL3=0;
             attractionStore2VL3=0;
         else
        attractionStore2UL3=attractionStore2UL3+(AttAv2L3*(L3AttractionVar2*kappa*(Store2Position(1)-XL3(j))/sqrt(deltaL3Store2Position*deltaL3Store2Position+delta2L3Store2Position*delta2L3Store2Position)))*exp(-((abs(deltaL3Store2Position)+theta1*abs(delta2L3Store2Position))/theta2))*dt;
        attractionStore2VL3=attractionStore2VL3+(AttAv2L3*(L3AttractionVar2*kappa*(Store2Position(2)-YL3(j))/sqrt(deltaL3Store2Position*deltaL3Store2Position+delta2L3Store2Position*delta2L3Store2Position)))*exp(-((abs(deltaL3Store2Position)+theta1*abs(delta2L3Store2Position))/theta2))*dt;
         end
        end

         if(XR1(j)>5)
                
             if(time2R1>=75)
             attractionStore2UR1=0;
             attractionStore2VR1=0;
            else
        attractionStore2UR1=attractionStore2UR1+(AttAv2R1*(R1AttractionVar2*kappa*(Store2Position(1)-XR1(j))/sqrt(deltaR1Store2Position*deltaR1Store2Position+delta2R1Store2Position*delta2R1Store2Position)))*exp(-((abs(deltaR1Store2Position)+theta1*abs(delta2R1Store2Position))/theta2))*dt;
        attractionStore2VR1=attractionStore2VR1+(AttAv2R1*(R1AttractionVar2*kappa*(Store2Position(2)-YR1(j))/sqrt(deltaR1Store2Position*deltaR1Store2Position+delta2R1Store2Position*delta2R1Store2Position)))*exp(-((abs(deltaR1Store2Position)+theta1*abs(delta2R1Store2Position))/theta2))*dt;
             end
        end
     
         if(XR2(j)>5)
             
             if(time2R2>=75)
             attractionStore2UL1=0;
             attractionStore2VL1=0;
             else
        attractionStore2UR2=attractionStore2UR2+(AttAv2R2*(R2AttractionVar2*kappa*(Store2Position(1)-XR2(j))/sqrt(deltaR2Store2Position*deltaR2Store2Position+delta2R2Store2Position*delta2R2Store2Position)))*exp(-((abs(deltaR2Store1Position)+theta1*abs(delta2R2Store2Position))/theta2))*dt;
        attractionStore2VR2=attractionStore2VR2+(AttAv2R2*(R2AttractionVar2*kappa*(Store2Position(2)-YR2(j))/sqrt(deltaR2Store2Position*deltaR2Store2Position+delta2R2Store2Position*delta2R2Store2Position)))*exp(-((abs(deltaR2Store1Position)+theta1*abs(delta2R1Store2Position))/theta2))*dt;
             end
        end
     
        if(XR3(j)>5)
            if(time2R3>=75)
             attractionStore2UR3=0;
             attractionStore2VR3=0;
            else
        attractionStore2UR3=attractionStore2UR3+(AttAv2R3*(R3AttractionVar2*kappa*(Store2Position(1)-XR3(j))/sqrt(deltaR3Store2Position*deltaR3Store2Position+delta2R3Store2Position*delta2R3Store2Position)))*exp(-((abs(deltaR3Store2Position)+theta1*abs(delta2R3Store2Position))/theta2))*dt;
        attractionStore2VR3=attractionStore2VR3+(AttAv2R3*(R3AttractionVar2*kappa*(Store2Position(2)-YR3(j))/sqrt(deltaR3Store2Position*deltaR3Store2Position+delta2R3Store2Position*delta2R3Store2Position)))*exp(-((abs(deltaR3Store2Position)+theta1*abs(delta2R3Store2Position))/theta2))*dt;
            end
         end
         
        
        %deciding if other pedestrians are in "my" line of sight. If i can
        %see them i need to avoid them.
    if(dot(AL1R1,BL1)>0)%if R1 "other" is in front of L1 "me"
        
       %avoidance terms for L1 are added to the initial values of
       %avoidanceUL1 and avoidanceVL1 if this condition is met. if the
       %condition is not met the values will remain 0.
        avoidanceUL1=avoidanceUL1+(-((kappa*(XR1(j)-XL1(j))/sqrt(deltaL1R1*deltaL1R1+delta2L1R1*delta2L1R1)))*exp(-((abs(deltaL1R1)+theta1*abs(delta2L1R1))/theta2))*dt);
        avoidanceVL1=avoidanceVL1+(-((kappa*(YR1(j)-YL1(j))/sqrt(deltaL1R1*deltaL1R1+delta2L1R1*delta2L1R1)))*exp(-((abs(deltaL1R1)+theta1*abs(delta2L1R1))/theta2))*dt);
      
    end
    
     if(dot(AR1L1,BR1)>0)%if L1 is in front of R1
       
        avoidanceUR1=avoidanceUR1+(-((kappa*(XL1(j)-XR1(j))/sqrt(deltaR1L1*deltaR1L1+delta2R1L1*delta2R1L1)))*exp(-((abs(deltaR1L1)+theta1*abs(delta2R1L1))/theta2))*dt);
        avoidanceVR1=avoidanceVR1+(-((kappa*(YL1(j)-YR1(j))/sqrt(deltaR1L1*deltaR1L1+delta2R1L1*delta2R1L1)))*exp(-((abs(deltaR1L1)+theta1*abs(delta2R1L1))/theta2))*dt);
        
     end
    
     if(dot(AL1R2,BL1)>0)%if R2 is in front of L1

        avoidanceUL1=avoidanceUL1+(-((kappa*(XR2(j)-XL1(j))/sqrt(deltaL1R2*deltaL1R2+delta2L1R2*delta2L1R2)))*exp(-((abs(deltaL1R2)+theta1*abs(delta2L1R2))/theta2))*dt);
        avoidanceVL1=avoidanceVL1+(-((kappa*(YR2(j)-YL1(j))/sqrt(deltaL1R2*deltaL1R2+delta2L1R2*delta2L1R2)))*exp(-((abs(deltaL1R2)+theta1*abs(delta2L1R2))/theta2))*dt);
      
     end
     
    if(dot(AR2L1,BR2)>0)
        
        avoidanceUR2=avoidanceUR2+(-((kappa*(XL1(j)-XR2(j))/sqrt(deltaR2L1*deltaR2L1+delta2R2L1*delta2R2L1)))*exp(-((abs(deltaR2L1)+theta1*abs(delta2R2L1))/theta2))*dt);
        avoidanceVR2=avoidanceVR2+(-((kappa*(YL1(j)-YR2(j))/sqrt(deltaR2L1*deltaR2L1+delta2R2L1*delta2R2L1)))*exp(-((abs(deltaR2L1)+theta1*abs(delta2R2L1))/theta2))*dt);
        
    end
      
    if(dot(AL1R3,BL1)>0)%if R3 is in front of L1
      
        avoidanceUL1=avoidanceUL1+(-((kappa*(XR3(j)-XL1(j))/sqrt(deltaL1R3*deltaL1R3+delta2L1R3*delta2L1R3)))*exp(-((abs(deltaL1R3)+theta1*abs(delta2L1R3))/theta2))*dt);
        avoidanceVL1=avoidanceVL1+(-((kappa*(YR3(j)-YL1(j))/sqrt(deltaL1R2*deltaL1R3+delta2L1R3*delta2L1R3)))*exp(-((abs(deltaL1R3)+theta1*abs(delta2L1R3))/theta2))*dt);
      
    end
    
   if(dot(AR3L1,BR3)>0)
       
        avoidanceUR3=avoidanceUR3+(-((kappa*(XL1(j)-XR3(j))/sqrt(deltaR3L1*deltaR3L1+delta2R3L1*delta2R3L1)))*exp(-((abs(deltaR3L1)+theta1*abs(delta2R3L1))/theta2))*dt);
        avoidanceVR3=avoidanceVR3+(-((kappa*(YL1(j)-YR3(j))/sqrt(deltaR3L1*deltaR3L1+delta2R3L1*delta2R3L1)))*exp(-((abs(deltaR3L1)+theta1*abs(delta2R3L1))/theta2))*dt);
        
   end
    
     if(dot(AL2R1,BL2)>0)
        
        avoidanceUL2=avoidanceUL2+(-((kappa*(XR1(j)-XL2(j))/sqrt(deltaL2R1*deltaL2R1+delta2L2R1*delta2L2R1)))*exp(-((abs(deltaL2R1)+theta1*abs(delta2L2R1))/theta2))*dt);
        avoidanceVL2=avoidanceVL2+(-((kappa*(YR1(j)-YL2(j))/sqrt(deltaL2R1*deltaL2R1+delta2L2R1*delta2L2R1)))*exp(-((abs(deltaL2R1)+theta1*abs(delta2L2R1))/theta2))*dt);
        
     end
     
    if(dot(AR1L2,BR1)>0)
        
        avoidanceUR1=avoidanceUR1+(-((kappa*(XL2(j)-XR1(j))/sqrt(deltaR1L2*deltaR1L2+delta2R1L2*delta2R1L2)))*exp(-((abs(deltaR1L2)+theta1*abs(delta2R1L2))/theta2))*dt);
        avoidanceVR1=avoidanceVR1+(-((kappa*(YL2(j)-YR1(j))/sqrt(deltaR1L2*deltaR1L2+delta2R1L2*delta2R1L2)))*exp(-((abs(deltaR1L2)+theta1*abs(delta2R1L2))/theta2))*dt);
    end
    
      if(dot(AL2R2,BL2)>0)
    
        avoidanceUL2=avoidanceUL2+(-((kappa*(XR2(j)-XL2(j))/sqrt(deltaL2R2*deltaL2R2+delta2L2R2*delta2L2R2)))*exp(-((abs(deltaL2R2)+theta1*abs(delta2L2R2))/theta2))*dt);
        avoidanceVL2=avoidanceVL2+(-((kappa*(YR2(j)-YL2(j))/sqrt(deltaL2R2*deltaL2R2+delta2L2R2*delta2L2R2)))*exp(-((abs(deltaL2R2)+theta1*abs(delta2L2R2))/theta2))*dt);
        
      end 
     
      if(dot(AR2L2,BR2)>0)
     
        avoidanceUR2=avoidanceUR2+(-((kappa*(XL2(j)-XR2(j))/sqrt(deltaR2L2*deltaR2L2+delta2R2L2*delta2R2L2)))*exp(-((abs(deltaR2L2)+theta1*abs(delta2R2L2))/theta2))*dt);
        avoidanceVR2=avoidanceVR2+(-((kappa*(YL2(j)-YR2(j))/sqrt(deltaR2L2*deltaR2L2+delta2R2L2*delta2R2L2)))*exp(-((abs(deltaR2L2)+theta1*abs(delta2R2L2))/theta2))*dt);
        
      end
      
       if(dot(AL2R3,BL2)>0)
      
        
        avoidanceUL2=avoidanceUL2+(-((kappa*(XR3(j)-XL2(j))/sqrt(deltaL2R3*deltaL2R3+delta2L2R3*delta2L2R3)))*exp(-((abs(deltaL2R3)+theta1*abs(delta2L2R3))/theta2))*dt);
        avoidanceVL2=avoidanceVL2+(-((kappa*(YR3(j)-YL2(j))/sqrt(deltaL2R3*deltaL2R3+delta2L2R3*delta2L2R3)))*exp(-((abs(deltaL2R3)+theta1*abs(delta2L2R3))/theta2))*dt);
        
       end 
      
        if(dot(AR3L2,BR3)>0)
         
        avoidanceUR3=avoidanceUR3+(-((kappa*(XL2(j)-XR3(j))/sqrt(deltaR3L2*deltaR3L2+delta2R3L2*delta2R3L2)))*exp(-((abs(deltaR3L2)+theta1*abs(delta2R3L2))/theta2))*dt);
        avoidanceVR3=avoidanceVR3+(-((kappa*(YL2(j)-YR3(j))/sqrt(deltaR3L2*deltaR3L2+delta2R3L2*delta2R3L2)))*exp(-((abs(deltaR3L2)+theta1*abs(delta2R3L2))/theta2))*dt);
        
        end
    
      
        if(dot(AL3R1,BL3)>0)%if R1 "other" is in front of L3 "me"
     
        avoidanceUL3=avoidanceUL3+(-((kappa*(XR1(j)-XL3(j))/sqrt(deltaL3R1*deltaL3R1+delta2L3R1*delta2L3R1)))*exp(-((abs(deltaL3R1)+theta1*abs(delta2L3R1))/theta2))*dt);
        avoidanceVL3=avoidanceVL3+(-((kappa*(YR1(j)-YL3(j))/sqrt(deltaL3R1*deltaL3R1+delta2L3R1*delta2L3R1)))*exp(-((abs(deltaL3R1)+theta1*abs(delta2L3R1))/theta2))*dt);
      
        end
    
     if(dot(AR1L3,BR1)>0)
       
        avoidanceUR1=avoidanceUR1+(-((kappa*(XL3(j)-XR1(j))/sqrt(deltaR1L3*deltaR1L3+delta2R1L3*delta2R1L3)))*exp(-((abs(deltaR1L3)+theta1*abs(delta2R1L3))/theta2))*dt);
        avoidanceVR1=avoidanceVR1+(-((kappa*(YL3(j)-YR1(j))/sqrt(deltaR1L3*deltaR1L3+delta2R1L3*delta2R1L3)))*exp(-((abs(deltaR1L3)+theta1*abs(delta2R1L3))/theta2))*dt);
        
     end
    
     if(dot(AL3R2,BL3)>0)%if R2 is in front of L3
      
       
        avoidanceUL3=avoidanceUL3+(-((kappa*(XR2(j)-XL3(j))/sqrt(deltaL3R2*deltaL3R2+delta2L3R2*delta2L3R2)))*exp(-((abs(deltaL3R2)+theta1*abs(delta2L3R2))/theta2))*dt);
        avoidanceVL3=avoidanceVL3+(-((kappa*(YR2(j)-YL3(j))/sqrt(deltaL3R2*deltaL3R2+delta2L3R2*delta2L3R2)))*exp(-((abs(deltaL3R2)+theta1*abs(delta2L3R2))/theta2))*dt);
      
     end
     
    if(dot(AR2L3,BR2)>0)
        
        avoidanceUR2=avoidanceUR2+(-((kappa*(XL3(j)-XR2(j))/sqrt(deltaR2L3*deltaR2L3+delta2R2L3*delta2R2L3)))*exp(-((abs(deltaR2L3)+theta1*abs(delta2R2L3))/theta2))*dt);
        avoidanceVR2=avoidanceVR2+(-((kappa*(YL3(j)-YR2(j))/sqrt(deltaR2L3*deltaR2L3+delta2R2L3*delta2R2L3)))*exp(-((abs(deltaR2L3)+theta1*abs(delta2R2L3))/theta2))*dt);
        
    end
      
    if(dot(AL3R3,BL3)>0)%if R3 is in front of L3
      
        avoidanceUL3=avoidanceUL3+(-((kappa*(XR3(j)-XL3(j))/sqrt(deltaL3R3*deltaL1R3+delta2L3R3*delta2L3R3)))*exp(-((abs(deltaL3R3)+theta1*abs(delta2L3R3))/theta2))*dt);
        avoidanceVL3=avoidanceVL3+(-((kappa*(YR3(j)-YL3(j))/sqrt(deltaL3R3*deltaL1R3+delta2L3R3*delta2L3R3)))*exp(-((abs(deltaL3R3)+theta1*abs(delta2L3R3))/theta2))*dt);
      
    end
    
   if(dot(AR3L3,BR3)>0)
          
        
        avoidanceUR3=avoidanceUR3+(-((kappa*(XL3(j)-XR3(j))/sqrt(deltaR3L3*deltaR3L3+delta2R3L3*delta2R3L3)))*exp(-((abs(deltaR3L3)+theta1*abs(delta2R3L3))/theta2))*dt);
        avoidanceVR3=avoidanceVR3+(-((kappa*(YL3(j)-YR3(j))/sqrt(deltaR3L3*deltaR3L3+delta2R3L3*delta2R3L3)))*exp(-((abs(deltaR3L3)+theta1*abs(delta2R3L3))/theta2))*dt);
        
   end
      
   
    if(dot(AL1L2,BL1)>0)%if L2 "other" is in front of L1 "me"
        
        avoidanceUL1=avoidanceUL1+(-((kappa*(XL2(j)-XL1(j))/sqrt(deltaL1L2*deltaL1L2+delta2L1L2*delta2L1L2)))*exp(-((abs(deltaL1L2)+theta1*abs(delta2L1L2))/theta2))*dt);
        avoidanceVL1=avoidanceVL1+(-((kappa*(YL2(j)-YL1(j))/sqrt(deltaL1L2*deltaL1L2+delta2L1L2*delta2L1L2)))*exp(-((abs(deltaL1L2)+theta1*abs(delta2L1L2))/theta2))*dt);
      
    end
    
     if(dot(AL2L1,BL2)>0)
           
        avoidanceUL2=avoidanceUL2+(-((kappa*(XL1(j)-XL2(j))/sqrt(deltaL2L1*deltaL2L1+delta2L2L1*delta2L2L1)))*exp(-((abs(deltaL2L1)+theta1*abs(delta2L2L1))/theta2))*dt);
        avoidanceVL2=avoidanceVL2+(-((kappa*(YL1(j)-YL2(j))/sqrt(deltaL2L1*deltaL2L1+delta2L2L1*delta2L2L1)))*exp(-((abs(deltaL2L1)+theta1*abs(delta2L2L1))/theta2))*dt);
        
     end
     
     if(dot(AL1L3,BL1)>0)%if L3 "other" is in front of L1 "me"
        
        avoidanceUL1=avoidanceUL1+(-((kappa*(XL3(j)-XL1(j))/sqrt(deltaL1L3*deltaL1L3+delta2L1L3*delta2L1L3)))*exp(-((abs(deltaL1L3)+theta1*abs(delta2L1L3))/theta2))*dt);
        avoidanceVL1=avoidanceVL1+(-((kappa*(YL3(j)-YL1(j))/sqrt(deltaL1L3*deltaL1L3+delta2L1L3*delta2L1L3)))*exp(-((abs(deltaL1L3)+theta1*abs(delta2L1L3))/theta2))*dt);
      
    end
    
     if(dot(AL3L1,BL3)>0)
       
        avoidanceUL3=avoidanceUL3+(-((kappa*(XL1(j)-XL3(j))/sqrt(deltaL3L1*deltaL3L1+delta2L3L1*delta2L3L1)))*exp(-((abs(deltaL3L1)+theta1*abs(delta2L3L1))/theta2))*dt);
        avoidanceVL3=avoidanceVL3+(-((kappa*(YL1(j)-YL3(j))/sqrt(deltaL3L1*deltaL3L1+delta2L3L1*delta2L3L1)))*exp(-((abs(deltaL3L1)+theta1*abs(delta2L3L1))/theta2))*dt);
        
     end
     
     if(dot(AL2L3,BL1)>0)%if L3 "other" is in front of L2 "me"
        
        avoidanceUL2=avoidanceUL2+(-((kappa*(XL3(j)-XL2(j))/sqrt(deltaL2L3*deltaL2L3+delta2L2L3*delta2L2L3)))*exp(-((abs(deltaL2L3)+theta1*abs(delta2L2L3))/theta2))*dt);
        avoidanceVL2=avoidanceVL2+(-((kappa*(YL3(j)-YL2(j))/sqrt(deltaL2L3*deltaL2L3+delta2L2L3*delta2L2L3)))*exp(-((abs(deltaL2L3)+theta1*abs(delta2L2L3))/theta2))*dt);
      
    end
    
     if(dot(AL3L2,BL3)>0)
       
        avoidanceUL3=avoidanceUL3+(-((kappa*(XL2(j)-XL3(j))/sqrt(deltaL3L2*deltaL3L2+delta2L3L2*delta2L3L2)))*exp(-((abs(deltaL3L2)+theta1*abs(delta2L3L2))/theta2))*dt);
        avoidanceVL3=avoidanceVL3+(-((kappa*(YL2(j)-YL3(j))/sqrt(deltaL3L2*deltaL3L2+delta2L3L2*delta2L3L2)))*exp(-((abs(deltaL3L2)+theta1*abs(delta2L3L2))/theta2))*dt);
        
     end
     
     
         if(dot(AR1R2,BR1)>0)%if R2 "other" is in front of R1 "me"
        
        avoidanceUR1=avoidanceUR1+(-((kappa*(XR2(j)-XR1(j))/sqrt(deltaR1R2*deltaR1R2+delta2R1R2*delta2R1R2)))*exp(-((abs(deltaR1R2)+theta1*abs(delta2R1R2))/theta2))*dt);
        avoidanceVR1=avoidanceVR1+(-((kappa*(YR2(j)-YR1(j))/sqrt(deltaR1R2*deltaR1R2+delta2R1R2*delta2R1R2)))*exp(-((abs(deltaR1R2)+theta1*abs(delta2R1R2))/theta2))*dt);
      
        end
    
     if(dot(AL2L1,BL2)>0)
       
        avoidanceUR2=avoidanceUR2+(-((kappa*(XR1(j)-XR2(j))/sqrt(deltaR2R1*deltaR2R1+delta2R2R1*delta2R2R1)))*exp(-((abs(deltaR2R1)+theta1*abs(delta2R2R1))/theta2))*dt);
        avoidanceVR2=avoidanceVR2+(-((kappa*(YR1(j)-YR2(j))/sqrt(deltaR2R1*deltaR2R1+delta2R2R1*delta2R2R1)))*exp(-((abs(deltaR2R1)+theta1*abs(delta2R2R1))/theta2))*dt);
        
     end
     
     if(dot(AR1R3,BR1)>0)%if R3 "other" is in front of R1 "me"
        
        avoidanceUR1=avoidanceUR1+(-((kappa*(XR3(j)-XR1(j))/sqrt(deltaR1R3*deltaR1R3+delta2R1R3*delta2R1R3)))*exp(-((abs(deltaR1R3)+theta1*abs(delta2R1R3))/theta2))*dt);
        avoidanceVR1=avoidanceVR1+(-((kappa*(YR3(j)-YR1(j))/sqrt(deltaR1R3*deltaR1R3+delta2R1R3*delta2R1R3)))*exp(-((abs(deltaR1R3)+theta1*abs(delta2R1R3))/theta2))*dt);
      
    end
    
     if(dot(AR3R1,BR3)>0)
       
        avoidanceUR3=avoidanceUR3+(-((kappa*(XR1(j)-XR3(j))/sqrt(deltaR3R1*deltaR3R1+delta2R3R1*delta2R3R1)))*exp(-((abs(deltaR3R1)+theta1*abs(delta2R3R1))/theta2))*dt);
        avoidanceVR3=avoidanceVR3+(-((kappa*(YR1(j)-YR3(j))/sqrt(deltaR3R1*deltaR3R1+delta2R3R1*delta2R3R1)))*exp(-((abs(deltaR3R1)+theta1*abs(delta2R3R1))/theta2))*dt);
        
     end
     
     if(dot(AR2R3,BR1)>0)%if R3 "other" is in front of R2 "me"
        
        avoidanceUR2=avoidanceUR2+(-((kappa*(XR3(j)-XR2(j))/sqrt(deltaR2R3*deltaR2R3+delta2R2R3*delta2R2R3)))*exp(-((abs(deltaR2R3)+theta1*abs(delta2R2R3))/theta2))*dt);
        avoidanceVR2=avoidanceVR2+(-((kappa*(YR3(j)-YR2(j))/sqrt(deltaR2R3*deltaR2R3+delta2R2R3*delta2R2R3)))*exp(-((abs(deltaR2R3)+theta1*abs(delta2R2R3))/theta2))*dt);
      
    end
    
     if(dot(AR3R2,BR3)>0)
       
        avoidanceUR3=avoidanceUR3+(-((kappa*(XR2(j)-XR3(j))/sqrt(deltaR3R2*deltaR3R2+delta2R3R2*delta2R3R2)))*exp(-((abs(deltaR3R2)+theta1*abs(delta2R3R2))/theta2))*dt);
        avoidanceVR3=avoidanceVR3+(-((kappa*(YR2(j)-YR3(j))/sqrt(deltaR3R2*deltaR3R2+delta2R3R2*delta2R3R2)))*exp(-((abs(deltaR3R2)+theta1*abs(delta2R3R2))/theta2))*dt);
        
     end
     
 
   
  %% the walk for L1 
 if (YL1(j)>2&&XL1(j)>1&&XL1(j)<3&&time1L1<75)%if they enter store one they will stay a while
        
        XL1(j+1)=XL1(j);
        UL1(j+1)=0;
        YL1(j+1)=YL1(j);
        VL1(j+1)=0;
        time1L1=time1L1+1;
    elseif (YL1(j)>-1.25&&YL1(j)<.25&&XL1(j)>5.75&&XL1(j)<8.25&&time2L1<75)%if they enter store 2 they will stay a while 
        XL1(j+1)=XL1(j);
        UL1(j+1)=0;
        YL1(j+1)=YL1(j);
        VL1(j+1)=0;
        time2L1=time2L1+1;
           
 else%otherwise walk like normal with the avoidance and attractions depending on what conditions are met
     XL1(j+1) = XL1(j) + UL1(j)*dt;
     UL1(j+1) = UL1(j) - alpha(1)*UL1(j)*(UL1(j)-u1)*(UL1(j)-u2)*dt + centripetalL1(1)*dt + sigma_x(1)*sqrt(dt)*randn(1)+avoidanceUL1+attractionStore1UL1+ attractionStore2UL1;
     YL1(j+1) = YL1(j) + VL1(j)*dt;
        if (YL1(j)<-1.5)||((YL1(j)>1.5)&&(XL1(j)>3));% if  they get to close to wall or sidewalk they will try to avoid it.
            VL1(j+1) = .2*abs(-(YL1(j))/(XL1(j)-5))*-YL1(j);
        else
            VL1(j+1) = VL1(j) - beta(1)*(YL1(j)-y_starL(XL1(j)))*dt - gamma(1)*VL1(j)*dt + centripetalL1(2)*dt + sigma_y(1)*sqrt(dt)*randn(1)+avoidanceVL1+attractionStore1VL1+ attractionStore2VL1;
        end
    end
       
    if (YR1(j)>2&&XR1(j)>1&&XR1(j)<3&&time1R1<75)
         XR1(j+1)=XR1(j);
        UR1(j+1)=0;
        YR1(j+1)=YR1(j);
        VR1(j+1)=0;
        time1R1=time1R1+1;
    elseif (YR1(j)>-1.25&&YR1(j)<.25&&XR1(j)>5.75&&XR1(j)<8.25&&time2R1<75)
        XR1(j+1)=XR1(j);
        UR1(j+1)=0;
        YR1(j+1)=YR1(j);
        VR1(j+1)=0;
        time2R1=time2R1+1;       
    else
         XR1(j+1) = XR1(j) + UR1(j)*dt(1);
         UR1(j+1) = UR1(j) - alpha(2)*UR1(j)*(UR1(j)-u1)*(UR1(j)-u2)*dt + centripetalR1(1)*dt + sigma_x(2)*sqrt(dt)*randn(1)+avoidanceUR1+attractionStore1UR1+ attractionStore2UR1;
          YR1(j+1) = YR1(j) + VR1(j)*dt(1);
        if (YR1(j)<-1.5)||((YR1(j)>1.5)&&(XR1(j)>3));
            VR1(j+1) = .2*abs(-(YR1(j))/(XR1(j)-5))*-YR1(j);
        else
            VR1(j+1) = VR1(j) - beta(2)*(YR1(j)-y_starR(XR1(j)))*dt - gamma(2)*VR1(j)*dt + centripetalR1(2)*dt + sigma_y(2)*sqrt(dt)*randn(1)+avoidanceVR1+attractionStore1VR1+ attractionStore2VR1;
        end  
    end
     
    
    if (YL2(j)>2&&XL2(j)>1&&XL2(j)<3&&time1L2<75)
        XL2(j+1)=XL2(j);
        UL2(j+1)=0;
        YL2(j+1)=YL2(j);
        VL2(j+1)=0;
        time1L2=time1L2+1;
        
    elseif (YL2(j)<.25&&YL2(j)>-1.25&&XL2(j)>5.75&&XL2(j)<8.25&&time2L2<75)
        XL2(j+1)=XL2(j);
        UL2(j+1)=0;
        YL2(j+1)=YL2(j);
        VL2(j+1)=0;
        time2L2=time2L2+1;

    else
         XL2(j+1)= XL2(j) + UL2(j)*dt;
         UL2(j+1) = UL2(j) - alpha(3)*UL2(j)*(UL2(j)-u1)*(UL2(j)-u2)*dt + centripetalL2(1)*dt + sigma_x(3)*sqrt(dt)*randn(1)+avoidanceUL2+attractionStore1UL2+ attractionStore2UL2;
         YL2(j+1) = YL2(j) + VL2(j)*dt;
        if (YL2(j)<-1.5)||((YL2(j)>1.5)&&(XL2(j)>3));
            VL2(j+1) = .2*abs(-(YL2(j))/(XL2(j)-5))*-YL2(j);
        else
            VL2(j+1) =  VL2(j) - beta(3)*(YL2(j)-y_starL(XL2(j)))*dt - gamma(3)*VL2(j)*dt + centripetalL2(2)*dt + sigma_y(3)*sqrt(dt)*randn(1)+avoidanceVL2+attractionStore1VL2+ attractionStore2VL2;
        end
    end
    
     if (YR2(j)>2&&XR2(j)>1&&XR2(j)<3&&time1R2<75)
        XR2(j+1)=XR2(j);
        UR2(j+1)=0;
        YR2(j+1)=YR2(j);
        VR2(j+1)=0;
        time1R2=time1R2+1;
        
     elseif (YR2(j)>-1.25&&YR2(j)<.25&&XR2(j)>5.75&&XR2(j)<8.25&&time2R2<75)
        XR2(j+1)=XR2(j);
        UR2(j+1)=0;
        YR2(j+1)=YR2(j);
        VR2(j+1)=0;
        time2R2=time2R2+1;

     else
         XR2(j+1)= XR2(j) + UR2(j)*dt;
        UR2(j+1) = UR2(j) - alpha(4)*UR2(j)*(UR2(j)-u1)*(UR2(j)-u2)*dt + centripetalR2(1)*dt(1) + sigma_x(4)*sqrt(dt)*randn(1)+avoidanceUR2+attractionStore1UR2+ attractionStore2UR2;
        YR2(j+1) = YR2(j) + VR2(j)*dt;
        if (YR2(j)<-1.5)||((YR2(j)>1.5)&&(XR2(j)>3));
            VR2(j+1) = .2*abs(-(YR2(j))/(XR2(j)-5))*-YR2(j);
        else
            VR2(j+1) = VR2(j) - beta(4)*(YR2(j)-y_starR(XR2(j)))*dt - gamma(4)*VR2(j)*dt + centripetalR2(2)*dt + sigma_y(4)*sqrt(dt)*randn(1)+avoidanceVR2+attractionStore1VR2+ attractionStore2VR2;
        end  
     end
     
     if (YL3(j)>2&&XL3(j)>1&&XL3(j)<3&&time1L3<75)
        XL3(j+1)=XL3(j);
        UL3(j+1)=0;
        YL3(j+1)=YL3(j);
        VL3(j+1)=0;
        time1L3=time1L3+1;
     elseif (YL3(j)>-1.25&&YL3(j)<.25&&XL3(j)>5.75&&XL3(j)<8.25&&time2L3<75)
         XL3(j+1)=XL3(j);
        UL3(j+1)=0;
        YL3(j+1)=YL3(j);
        VL3(j+1)=0;
        time2L3=time2L3+1;
     else
        XL3(j+1) = XL3(j) + UL3(j)*dt;
        UL3(j+1) = UL3(j) - alpha(5)*UL3(j)*(UL3(j)-u1)*(UL3(j)-u2)*dt + centripetalL3(1)*dt + sigma_x(5)*sqrt(dt)*randn(1)+avoidanceUL3+attractionStore1UL3+ attractionStore2UL3;
        YL3(j+1) = YL3(j) + VL3(j)*dt;
        if (YL3(j)<-1.5)||((YL3(j)>1.5)&&(XL3(j)>3));
            VL3(j+1) = .2*abs(-(YL3(j))/(XL3(j)-5))*-YL3(j);
        else
            VL3(j+1) = VL3(j) - beta(5)*(YL3(j)-y_starL(XL3(j)))*dt - gamma(5)*VL3(j)*dt + centripetalL3(2)*dt + sigma_y(5)*sqrt(dt)*randn(1)+avoidanceVL3+attractionStore1VL3+ attractionStore2VL3;
        end
     end

      if (YR3(j)>2&&XR3(j)>1&&XR3(j)<3&&time1R3<75)
        XR3(j+1)=XR3(j);
        UR3(j+1)=0;
        YR3(j+1)=YR3(j);
        VR3(j+1)=0;
        time1R3=time1R3+1;
      elseif (YR3(j)>-1.25&&YR3(j)<.25&&XR3(j)>5.75&&XR3(j)<8.25&&time2R3<75)
        XR3(j+1)=XR3(j);
        UR3(j+1)=0;
        YR3(j+1)=YR3(j);
        VR3(j+1)=0;
        time2R3=time2R3+1; 
      else
        XR3(j+1)= XR3(j) + UR3(j)*dt;
        UR3(j+1) = UR3(j) - alpha(6)*UR3(j)*(UR3(j)-u1)*(UR3(j)-u2)*dt + centripetalR3(1)*dt(1) + sigma_x(6)*sqrt(dt)*randn(1)+avoidanceUR3+attractionStore1UR3+ attractionStore2UR3;
        YR3(j+1) = YR3(j) + VR3(j)*dt;
        if (YR3(j)<-1.5)||((YR3(j)>1.5)&&(XR3(j)>3));
            VR3(j+1) = .2*abs(-(YR3(j))/(XR3(j)-5))*-YR3(j);
        else
             VR3(j+1) = VR3(j) - beta(6)*(YR3(j)-y_starR(XR3(j)))*dt - gamma(6)*VR3(j)*dt + centripetalR3(2)*dt + sigma_y(6)*sqrt(dt)*randn(1)+avoidanceVR3+attractionStore1VR3+ attractionStore2VR3;
        end  
      end
        
        j=j+1;
    end 
       


    

 j=length(XL1)-1;
PedL1=[(0:dt:j*dt)',XL1', YL1', UL1', VL1'];
PedR1=[(0:dt:j*dt)',XR1', YR1', UR1', VR1'];
PedL2=[(0:dt:j*dt)',XL2', YL2', UL2', VL2'];
PedR2=[(0:dt:j*dt)',XR2', YR2', UR2', VR2'];
PedL3=[(0:dt:j*dt)',XL3', YL3', UL3', VL3'];
PedR3=[(0:dt:j*dt)',XR3', YR3', UR3', VR3'];
for ti=1:length(PedL1(:,1))
    rectangle('Position',[1 1.75 2 1])
    rectangle('Position',[6 -1 2 1])
    rectangle('Position',[5.75 -1.25 2.5 1.5])
    axis([0 10.8 -1 1])
    axis equal   
    plot([0 10.8],[-2 -2],'k-','linewidth',2)
    hold on
    plot([0 10.8],[2 2],'k-','linewidth',2)
plot(PedL1(ti,2),PedL1(ti,3),'ro')
hold on
plot(PedR1(ti,2),PedR1(ti,3),'rx')
hold on
plot(PedL2(ti,2),PedL2(ti,3),'bx')
hold on
plot(PedR2(ti,2),PedR2(ti,3),'bo')
hold on
plot(PedL3(ti,2),PedL3(ti,3),'ko')
hold on
plot(PedR3(ti,2),PedR3(ti,3),'kx')

pause(.05)
end
    end

% [PedL1,PedR1, PedL2, PedR2, PedL3,PedR3]=sixPeds5([0 .1 .2 .4 .5 .7 ],[.2 0 -.1 .1 0 -.2],[-1 1 -1 1 -1 1],[sqrt(.0009)*randn(1),sqrt(.0009)*randn(1),sqrt(.0009)*randn(1),sqrt(.0009)*randn(1),sqrt(.0009)*randn(1),sqrt(.0009)*randn(1)],[1], [-1], [(1/15)],[2.9 2.9 2.9 2.9 2.9 2.9],[1.25 1.25 1.25 1.25 1.25 1.25],[.7 .7 .7 .7 .7 .7],[.48 .48 .48 .48 .48 .48],[.5 .5 .5 .5 .5 .5],1.37,.4,.8,.3)
% plot(PedR(:,2),PedR(:,3))
% hold on
% plot(PedL(:,2),PedL(:,3))