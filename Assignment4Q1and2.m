function [] = Assignment4Q1and2(nElec)

    %FLAGS

    ScatterFlag = 1;
    BoxFlag = 1;
    ReflectFlag = 0;
    InjectFlag = 0;



    q_0 = 1.60217653e-19;             % electron charge
    %hb = 1.054571596e-34;             % Dirac constant
    %h = C.hb * 2 * pi;                % Planck constant
    m_0 = 9.10938215e-31;             % electron mass
    kb = 1.3806504e-23;               % Boltzmann constant
    %eps_0 = 8.854187817e-12;          % vacuum permittivity
    %mu_0 = 1.2566370614e-6;           % vacuum permeability
    %c = 299792458;                    % speed of light
    %g = 9.80665;                      % metres (32.1740 ft) per s²
    Ne = 2.7182818;                    % eulers number
    tm = 0.2e-12;                       % Mean Collision time
    
    %CONSTANTS
    Vth = sqrt(2 * kb * 300/ (0.26 * m_0));
    
    %MEAN FREE PATH CALCULATED
%     MFPCalc = Vth * 0.2e-12
for     VoltX = 0.1:1:10.1
    %TIME
    t = 0;
    dt = 1e-13;
    TStop = 100*dt; % shortened to 100 from 1000 for testing

    %COUNT
    count = 1;
    nElecCount = 0;

    %BOUNDARIES
    xMax = 200e-9;
    yMax = 100e-9;
    Limits = [0 xMax 0 yMax];
    LimitsTemp = [0 TStop 0 600];
    LimitsScatter = [0 TStop 0 1];
    LimitsMFP = [0 TStop 0 3e-09];
    LimitsTmn = [0 TStop 0 5e-15];

    %Scattering probability
    pScat = 1 - Ne^(-(dt/(0.2e-12)));
    scattercount = 0;

    %Drawing of boxes of side length 40
  
    Box1 = [80e-9 0;80e-9 40e-9;120e-9 40e-9; 120e-9 0; 80e-9 0;];
    Box2 = [80e-9 60e-9;80e-9 100e-9;120e-9 100e-9; 120e-9 60e-9; 80e-9 60e-9;];

    %Old position values
    xOld = zeros(1,nElec);
    yOld = zeros(1,nElec);
    
    % Randomly place electrons
    
    x(1, :) = zeros(1,nElec); % injection from left side
    y(1, :) = rand(1,nElec) * yMax; % spread evenly
    
    
    % give each particle Vth but with a random direction
    Vx(1:nElec) = Vth/sqrt(2) * randn(1,nElec);
    Vy(1:nElec) = Vth/sqrt(2) * randn(1,nElec);
    
    %set scatter Avg to zero to start
    ScatterAvg = 0;

    %Initialize V field

    % EfieldX = VoltX/xMax;
    % AccelX = (EfieldX * q_0) / (0.26 * m_0);
    % AccelY = 0;
    AccelX = zeros(1,nElec);
    AccelY = zeros(1,nElec);


    if InjectFlag == 0
        x(1, :) = rand(1,nElec) * xMax; % spread evenly
        nElecCount = nElec;
    end



%% Assignment 2
L = xMax * 1e9;
W = yMax * 1e9;
F = zeros(L*W,1);
G = sparse(L*W, L*W);
Box = 40;
sigma = 0.01;
ConMatrix = ones(L,W);
ConMatrix(L/2 - Box/2:L/2 + Box/2,1:Box) = sigma;
ConMatrix(L/2 - Box/2:L/2 + Box/2,W-Box:W) = sigma;


for i = 1:1:L
    for j = 1:1:W
        N = (i-1) * W + j;
        NXM = (i-2) * W + j;
        NXP = i * W + j;
        NYM = (i-1) * W + j-1;
        NYP = (i-1) * W + j+1;

        %MXM = (ConMatrix(i,j) + ConMatrix(i-1,j));
        %MXP = (ConMatrix(i,j) + ConMatrix(i+1,j));
        %MYM = (ConMatrix(i,j) + ConMatrix(i,j-1));
        %MYP = (ConMatrix(i,j) + ConMatrix(i,j+1));

        if i == 1 %     LEFT
            F(N,1) = VoltX;
            G(N,N) = 1;
        elseif i == L % RIGHT
            F(N,1) = 0;
            G(N,N) = 1;
        elseif j == 1 % BOTTOM
            MXM = (ConMatrix(i,j) + ConMatrix(i-1,j));
            MXP = (ConMatrix(i,j) + ConMatrix(i+1,j));
            %MYM = (ConMatrix(i,j) + ConMatrix(i,j-1));
            MYP = (ConMatrix(i,j) + ConMatrix(i,j+1));

            G(N,N)   = -(MXM + MXP + MYP);
            G(N,NXM)   = MXM;
            G(N,NXP)   = MXP;
            G(N,NYP)   = MYP;
        elseif j == W % TOP
            MXM = (ConMatrix(i,j) + ConMatrix(i-1,j));
            MXP = (ConMatrix(i,j) + ConMatrix(i+1,j));
            MYM = (ConMatrix(i,j) + ConMatrix(i,j-1));
            %MYP = (ConMatrix(i,j) + ConMatrix(i,j+1));

            G(N,N)   = -(MXM + MXP + MYP);
            G(N,NXM)   = MXM;
            G(N,NXP)   = MXP;
            G(N,NYP)   = MYM;
        else %           MID
            MXM = (ConMatrix(i,j) + ConMatrix(i-1,j));
            MXP = (ConMatrix(i,j) + ConMatrix(i+1,j));
            MYM = (ConMatrix(i,j) + ConMatrix(i,j-1));
            MYP = (ConMatrix(i,j) + ConMatrix(i,j+1));

            G(N,N)   = -(MXM + MXP + MYM + MYP);
            G(N,NXM) = MXM;
            G(N,NXP) = MXP;
            G(N,NYM) = MYM;
            G(N,NYP) = MYP;
        end

    end

end

VMatrix = reshape(G\F , [W,L]);
[Ex,Ey] = gradient(-VMatrix);


% figure(1)
% title("Voltage Plot")
% surf(VMatrix')
% xlabel("Distance (m)")
% ylabel("Distance (m)")
% zlabel("V Plot")
% 
% figure(2);
% quiver(Ex,Ey);
% title('E field') ;
% xlabel("Distance (m)")
% ylabel("Distance (m)")
% zlabel("E field")


   %% Main Loop

    while t < TStop


            % Iterate time
            t  = t + dt;  

            % Get old positions
            xOld(1:nElecCount) = x(1:nElecCount);
            yOld(1:nElecCount) = y(1:nElecCount);
            VxOld(1:nElecCount) = Vx(1:nElecCount);
            VyOld(1:nElecCount) = Vy(1:nElecCount);
            

            % Get new Positions
            x(1:nElecCount) = x(1:nElecCount) + (Vx(1:nElecCount) .* dt);
            y(1:nElecCount) = y(1:nElecCount) + (Vy(1:nElecCount) .* dt);

            %Get plot arrays using new - old
            xPlot = [xOld(:) x(:)];
            yPlot = [yOld(:) y(:)];
            %vOldx = Vx;
            %vOldy = Vy;
            
                        
            
            %% Calculations

            VAvg = (sqrt((abs(Vx)).^2 + (abs(Vy)).^2)); % get current average Velocity
            TAvg = mean((VAvg.^2 .* (0.26*m_0))/ (2 * kb)); % get current average Temp
            MFPcurrent = mean(VAvg) * dt / scattercount; % gets current MFP
            Tmn = dt / scattercount; % gets Tmn
            ICurrentX(count,:) = q_0 * 1e20 * mean(Vx) * xMax;
            TempPlot(count,:) = TAvg; % store array of Temp
            TmnPlot(count,:) = Tmn;  % store array of Tmn
            scatterPlot(count,:) = ScatterAvg; % store array of scatter probabilities
            timePlot(count,:) = [t]; % store array of times
            MFP(count,:) =  MFPcurrent; % length of paths / number of scatters
            scattercount = 0;
            %% Plots
%             figure(3);
%             subplot(2,2,1);
%             hold on
%             if BoxFlag == 1
%             subplot(2,2,1),plot(Box1(1:5,1),Box1(1:5,2),'k');
%             subplot(2,2,1),plot(Box2(1:5,1),Box2(1:5,2),'k');
%             end
%             if(nElecCount >= 1)
%             subplot(2,2,1),plot(xPlot(1,1:2), yPlot(1,1:2),'b');
%             end
%             if(nElecCount >= 14)
%             subplot(2,2,1),plot(xPlot(11,1:2), yPlot(11,1:2),'r');
%             end
%             if(nElecCount >= 24)
%             subplot(2,2,1),plot(xPlot(21,1:2), yPlot(21,1:2),'g');
%             end
%             if(nElecCount >= 34)
%             subplot(2,2,1),plot(xPlot(31,1:2), yPlot(31,1:2),'c');
%             end
%             if(nElecCount >= 44)
%             subplot(2,2,1),plot(xPlot(41,1:2), yPlot(41,1:2),'y');
%             end
%             if(nElecCount >= 54)
%             subplot(2,2,1),plot(xPlot(51,1:2), yPlot(51,1:2),'m');
%             end
% 
%             %quiver(x,y,Vx,Vy);
%             hold off
%             axis(Limits);
%             xlabel('X (m)');
%             ylabel('Y (m)');
%             grid on
%             
%             subplot(2,2,2); %plot Avg Velocity
%             plot(t,TAvg,'v','linewidth', 2);
%             hold on
%             subplot(2,2,2), plot(timePlot,TempPlot,'linewidth', 2);
%             hold off
%             axis(LimitsTemp);
%             xlabel('Time (s)');
%             ylabel('Temperature (K)');
%             grid on
% 
%             subplot(2,2,4); %plot MFP
%             plot(t,MFPcurrent,'v','linewidth', 2);
%             hold on
%             subplot(2,2,4), plot(timePlot,MFP,'linewidth', 2);
%             hold off
%             axis(LimitsMFP);
%             xlabel('Time (s)');
%             ylabel('Mean Free Path (m)');
%             grid on
% 

            
%             subplot(2,2,3); %plot Tmn
%             plot(t,Tmn,'v','linewidth', 2);
%             hold on
%             subplot(2,2,3), plot(timePlot,TmnPlot,'linewidth', 2);
%             hold off
%             axis(LimitsTmn);
%             xlabel('Time (s)');
%             ylabel('Tmn (s)');
%             grid on

%             subplot(2,2,3); %plot Scattering probability
%             plot(t,ScatterAvg,'v','linewidth', 2);
%             hold on
%             subplot(2,2,3), plot(timePlot,scatterPlot,'linewidth', 2);
%             hold off
%             axis(LimitsScatter);
%             xlabel('Time (s)');
%             ylabel('Scatter Probablity');
%             grid on

            count = count + 1;
            pause(0.0001)

            %% Testing "for" loops
            for i=1:1:nElecCount

                if ScatterFlag == 1
                %Scatter Test
                rScat = rand(1,nElecCount);
                if rScat(i) <= pScat
                scattercount = scattercount+1;
                Vx(i) = Vth/sqrt(2) * randn();
                Vy(i) = Vth/sqrt(2) * randn();
                end
                end

if ReflectFlag == 1
                %reflect on X
                if x(i) <= 0  %check if out of bounds 
                  Vx(i) = Vx(i) * -1; 
                  x(i) = 0;
                elseif  x(i) >= xMax  %check if out of bounds
                  Vx(i) = Vx(i) * -1;
                  x(i) = xMax;
                end
else
                %transparent x
                if x(i) <= 0  %check if out of bounds 
                  x(i) = xMax;
                  xOld(i) = xOld(i) + xMax;
                elseif  x(i) >= xMax  %check if out of bounds
                  x(i) = 0;
                  xOld(i) = xOld(i) - xMax;
                end
end
                %reflect on Y
                if y(i) <= 0 %check if out of bounds
                  Vy(i) = Vy(i) * -1; 
                  y(i) = 0;
                elseif  y(i) >= yMax  %check if out of bounds 
                  Vy(i) = Vy(i) * -1;
                  y(i) = yMax;
                end

%% Box Boundaries code

%Box 1 & 2
if BoxFlag == 1
                if (y(i) <= 40e-9 || y(i) >= 60e-9) && x(i) >= 80e-9 && x(i) <= 120e-9 % check if inside a box
                    if y(i) <= 40e-9 && yOld(i) >= 40e-9 %check if electron came from above
                        Vy(i) = Vy(i) * -1;
                        y(i) = 40e-9;
                     elseif y(i) >= 60e-9 && yOld(i) <= 60e-9 %check if electron came from below
                         Vy(i) = Vy(i) * -1;
                         y(i) = 60e-9;
                      elseif x(i) >= 80e-9 && xOld(i) <= 80e-9 %check if electron came from left side
                          Vx(i) = Vx(i) * -1;
                          x(i) = 80e-9;
                      else % Else electron came from right side
                          Vx(i) = Vx(i) * -1;
                          x(i) = 120e-9;
                    end
                end
end

            if round(x(i) * 1e9) ~= 0  && round(y(i) * 1e9) ~= 0 && round(x(i) * 1e9) <= xMax * 1e9 && round(y(i) * 1e9) <= yMax * 1e9
            % E field to Accel
                AccelX(i) = q_0 * Ex(round(y(i) * 1e9),round(x(i) * 1e9)) / (xMax/(xMax* 1e9)) / (0.26 * m_0);
                AccelY(i) = q_0 * Ey(round(y(i) * 1e9),round(x(i) * 1e9)) / (xMax/(xMax* 1e9)) / (0.26 * m_0);
            end
            % Update Velocity
                Vx(i) = Vx(i) + (AccelX(i) * dt);
                Vy(i) = Vy(i) + (AccelY(i) * dt);

            % iterate number of electrons

            if nElecCount < nElec
                nElecCount = nElecCount+10;
            end
    
        end   
            ScatterAvg = scattercount/nElecCount; % get current chance to scatter
    end

IAvg(VoltX + 0.9,:) = mean(ICurrentX);

end
%     figure(4)
%     xlin = linspace(0,xMax,200);
%     ylin = linspace(0,yMax,100);
%     [meshx,meshy] = meshgrid(xlin,ylin);
%     Temps = (sqrt((Vx(1:nElec).^2)+(Vy(1:nElec).^2)).^2 .* 0.26 .* m_0)./(2*kb);
%     TempMap = griddata(x,y,Temps,meshx,meshy);
%     surf(meshx,meshy,TempMap)
%     title("Temperature Map")
%     xlabel("Distance (m)")
%     ylabel("Distance (m)")
%     zlabel("Temp (K)")
% 
%     figure(5)
%     Exlin = linspace(0,xMax,101);
%     Eylin = linspace(0,yMax,51);
%     ExlinMesh = linspace(0,xMax,100);
%     EylinMesh = linspace(0,yMax,50);
%     [Emeshx,Emeshy] = meshgrid(ExlinMesh,EylinMesh);
%     Emap = histcounts2(y,x,Eylin,Exlin);
%     surf(Emeshx,Emeshy,Emap)
%     title("Density Map")
%     xlabel("Distance (m)")
%     ylabel("Distance (m)")
%     zlabel("Electrons")

    figure(4)
    fit = polyfit(0.1:1:10.1, IAvg, 1);
    val = polyval(fit, 0.1:1:10.1);
    plot(0.1:1:10.1,IAvg);
    hold on
    plot(0.1:1:10.1, val);

    title("Current Plot")
    xlabel("Voltage(V)")
    ylabel("Current(A)")
    R = 1/fit(1);
    fprintf('The linear fit gives a resistance of R3 = %e', R);


end






