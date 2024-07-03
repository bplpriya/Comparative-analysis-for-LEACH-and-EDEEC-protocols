

close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%% Network Establishment Parameters %%%%%%%%%%%%%%%%%%%%

%%% Area of Operation %%%

% Field Dimensions in meters %
xm=100;
ym=100;
x=0; % added for better display results of the plot
y=0; % added for better display results of the plot
% Number of Nodes in the field %
n=100;
% Number of Dead Nodes in the beggining %
dead_nodes=0;
% Coordinates of the Sink (location is predetermined in this simulation) %
sinkx=50;
sinky=200;

%%% Energy Values %%%
% Initial Energy of a Node (in Joules) % 
Eo=2; % units in Joules
% Energy required to run circuity (both for transmitter and receiver) %
Eelec=50*10^(-9); % units in Joules/bit
ETx=50*10^(-9); % units in Joules/bit
ERx=50*10^(-9); % units in Joules/bit
% Transmit Amplifier Types %
Eamp=100*10^(-12); % units in Joules/bit/m^2 (amount of energy spent by the amplifier to transmit the bits)
% Data Aggregation Energy %
EDA=5*10^(-9); % units in Joules/bit
% Size of data package %
k=4000; % units in bits
% Suggested percentage of cluster head %
p=0.05; % a 5 percent of the total amount of nodes used in the network is proposed to give good results
% Number of Clusters %
No=p*n; 
% Round of Operation %
rnd=0;
% Current Number of operating Nodes %
operating_nodes=n;
transmissions=0;
temp_val=0;
flag1stdead=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%



            %%% Creation of the Wireless Sensor Network %%%


            Loc = xlsread('Loc.xlsx');
% Plotting the WSN %
for i=1:n
    
    SN(i).id=i;	% sensor's ID number
    %SN(i).x=rand(1,1)*xm;	% X-axis coordinates of sensor node
    %SN(i).y=rand(1,1)*ym;	% Y-axis coordinates of sensor node
    SN(i).x = Loc(i,1);
    SN(i).y = Loc(i,2);
    SN(i).E=Eo;     % nodes energy levels (initially set to be equal to "Eo"
    SN(i).role=0;   % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
    SN(i).cluster=0;	% the cluster which a node belongs to
    SN(i).cond=1;	% States the current condition of the node. when the node is operational its value is =1 and when dead =0
    SN(i).rop=0;	% number of rounds node was operational
    SN(i).rleft=0;  % rounds left for node to become available for Cluster Head election
    SN(i).dtch=0;	% nodes distance from the cluster head of the cluster in which he belongs
    SN(i).dts=0;    % nodes distance from the sink
    SN(i).tel=0;	% states how many times the node was elected as a Cluster Head
    SN(i).rn=0;     % round node got elected as cluster head
    SN(i).chid=0;   % node ID of the cluster head which the "i" normal node belongs to

    
    hold on; % This command keeps the current plot and adds subsequent plots to it. This allows multiple plots to be displayed on the same figure.
    figure(1)  % This command specifies that the current plot should be drawn on figure 1.
    % below line plots the sensor nodes (SN(i).x, SN(i).y) as blue circles ('ob'), the sink as a red star ('*r'), and optionally any predefined shape specified by x and y. 
    plot(x,y,xm,ym,SN(i).x,SN(i).y,'ob',sinkx,sinky,'*r'); 
    title 'Wireless Sensor Network';
    xlabel '(m)';
    ylabel '(m)';
    
end
 


                      %%%%%% Set-Up Phase %%%%%% 
                      
             
%% MODIFICATION

 % categorize them into levels
a1=1;
a2=1;
a3=1;
a4=1;

for i=1:1:n
        
        % Categorize the node into levels
        if (SN(i).x >= 0  && SN(i).x <=100) && (SN(i).y >= 75 && SN(i).y <= 100)
            level1(a1) = SN(i).id;
            a1=a1+1;
        

        elseif (SN(i).x >= 0  && SN(i).x <=100) && (SN(i).y >= 50 && SN(i).y < 75)
            level2(a2) = SN(i).id;
            a2=a2+1;
        

        elseif (SN(i).x >= 0  && SN(i).x <=100) && (SN(i).y >= 25 && SN(i).y < 50)
            level3(a3) = SN(i).id;
            a3=a3+1;
        

        elseif (SN(i).x >= 0  && SN(i).x <=100) && (SN(i).y >= 0 && SN(i).y < 25)
            level4(a4) = SN(i).id;
            a4=a4+1;
        end

end 
       
%% end of modification
while operating_nodes>0
        
    % Displays Current Round %     
    rnd     

	% Threshold Value %
    % (1-p*(mod(rnd,1/p))) represents the probability of not having the event occur in the current round.
	t=(p/(1-p*(mod(rnd,1/p))));
     
    % Re-election Value %
    % Calculates the same fractional part of the current round number with respect to the expected number of rounds before a particular event occurs (1/p). 
    % This value may be used to track the remaining rounds before the event is expected to happen.
    tleft=mod(rnd,1/p);
 
	% Reseting Previous Amount Of Cluster Heads In the Network %
	CLheads=0;
    
    % Reseting Previous Amount Of Energy Consumed In the Network on the Previous Round %
    energy=0;
 
  %% modification 

% Cluster Heads Election %

for i=level1
    SN(i).cluster=0;    % reseting cluster in which the node belongs to
    SN(i).role=0;       % reseting node role
    SN(i).chid=0;       % reseting cluster head id
    
    % Cluster Head Election for Level 1
    %if (SN(i).x >= 0  && SN(i).x <=100) && (SN(i).y >= 75 && SN(i).y <= 100)
        % check if the node is in level1
        generate=rand;	
        if generate < t
            SN(i).role=1;	    % assigns the node role of a cluster head
            SN(i).rn=rnd;	    % Assigns the round that the cluster head was elected to the data table
            SN(i).tel=SN(i).tel + 1;   
            SN(i).rleft=1/p-tleft;    % rounds for which the node will be unable to become a CH
            SN(i).dts=sqrt((sinkx-SN(i).x)^2 + (sinky-SN(i).y)^2); % calculates the distance between the sink and the cluster head
            CLheads=CLheads+1;	    % sum of cluster heads that have been elected 
            SN(i).cluster=CLheads;  % cluster of which the node got elected to be cluster head
            CL(CLheads).x=SN(i).x;  % X-axis coordinates of elected cluster head
            CL(CLheads).y=SN(i).y;  % Y-axis coordinates of elected cluster head
            CL(CLheads).id=i;      % Assigns the node ID of the newly elected cluster head to an array
            fprintf('Node %d elected as a cluster head for Level 1\n', i);
        end
end
    
    
% Grouping the Nodes into Clusters & caclulating the distance between node and cluster head %
   % for level1
 for i=level1
    if SN(i).role==0 && SN(i).E>0 && CLheads>0
        % if node is normal
        d=sqrt((CL(1).x-SN(i).x)^2 + (CL(1).y-SN(i).y)^2);
        % calculate the distance 'd' between the sensor node and the cluster head
        for m=2:CLheads
            % calculate distance between node and other cluster heads
            d_temp=sqrt((CL(m).x-SN(i).x)^2 + (CL(m).y-SN(i).y)^2);
            if d_temp < d
                % find the minimum distance to a cluster head
                d=d_temp;
            end
        end
        SN(i).cluster=1;  % assign node to level 1 cluster
        SN(i).dtch=d;     % assign the distance of node to CH
        SN(i).chid=1;     % assign level 1 cluster head ID
        fprintf('Node %d belongs to Level 1 cluster\n', i);
    end
end
       

% Data Aggregation Phase % for level1
for i=level1
    if SN(i).role==1 && SN(i).E>0
        % if node is cluster head
        ETx=(Eelec+EDA)*k + Eamp * k * SN(i).dts^2;
        % calculate energy dissipation for cluster head to transmit to sink
        SN(i).E=SN(i).E - ETx;  % update energy of cluster head
        energy=energy+ETx;      % update total energy
        if SN(i).E <= 0
            % if cluster head's energy depletes with transmission
            dead_nodes=dead_nodes +1;
            operating_nodes= operating_nodes - 1;
            SN(i).cond=0;
            SN(i).rop=rnd;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cluster Heads Election for Level 2
for i=level2
    %if (SN(i).x >= 0  && SN(i).x <=100) && (SN(i).y >= 50 && SN(i).y < 75)
        % check if the node is in level 2
        generate=rand;	
        if generate < t
            SN(i).role=1;	    % assigns the node role of a cluster head
            SN(i).rn=rnd;	    % Assigns the round that the cluster head was elected to the data table
            SN(i).tel=SN(i).tel + 1;   
            SN(i).rleft=1/p-tleft;    % rounds for which the node will be unable to become a CH
            SN(i).dts=sqrt((sinkx-SN(i).x)^2 + (sinky-SN(i).y)^2); % calculates the distance between the sink and the cluster head
            CLheads=CLheads+1;	    % sum of cluster heads that have been elected 
            SN(i).cluster=CLheads;  % cluster of which the node got elected to be cluster head
            CL(CLheads).x=SN(i).x;  % X-axis coordinates of elected cluster head
            CL(CLheads).y=SN(i).y;  % Y-axis coordinates of elected cluster head
            CL(CLheads).id=i;      % Assigns the node ID of the newly elected cluster head to an array
            fprintf('Node %d elected as a cluster head for Level 2\n', i);
        end
end

% Grouping the Nodes into Clusters for Level 2
for i=level2
    if SN(i).role==0 && SN(i).E>0 && CLheads>0
        % if node is normal
        d=sqrt((CL(1).x-SN(i).x)^2 + (CL(1).y-SN(i).y)^2);
        % calculate the distance 'd' between the sensor node and the cluster head
        for m=2:CLheads
            % calculate distance between node and other cluster heads
            d_temp=sqrt((CL(m).x-SN(i).x)^2 + (CL(m).y-SN(i).y)^2);
            if d_temp < d
                % find the minimum distance to a cluster head
                d=d_temp;
            end
        end
        SN(i).cluster=2;  % assign node to level 2 cluster
        SN(i).dtch=d;     % assign the distance of node to CH
        SN(i).chid=2;     % assign level 2 cluster head ID
        fprintf('Node %d belongs to Level 2 cluster\n', i);
    end
end

% Data Aggregation Phase % for level2
for i=level2
    if SN(i).role==1 && SN(i).E>0
        % if node is cluster head
        ETx=(Eelec+EDA)*k + Eamp * k * SN(i).dts^2;
        % calculate energy dissipation for cluster head to transmit to sink
        SN(i).E=SN(i).E - ETx;  % update energy of cluster head
        energy=energy+ETx;      % update total energy
        if SN(i).E <= 0
            % if cluster head's energy depletes with transmission
            dead_nodes=dead_nodes +1;
            operating_nodes= operating_nodes - 1;
            SN(i).cond=0;
            SN(i).rop=rnd;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cluster Heads Election for Level 3
for i=level3
    %if (SN(i).x >= 0  && SN(i).x <=100) && (SN(i).y >= 25 && SN(i).y < 50)
        % check if the node is in level 3
        generate=rand;	
        if generate < t
            SN(i).role=1;	    % assigns the node role of a cluster head
            SN(i).rn=rnd;	    % Assigns the round that the cluster head was elected to the data table
            SN(i).tel=SN(i).tel + 1;   
            SN(i).rleft=1/p-tleft;    % rounds for which the node will be unable to become a CH
            SN(i).dts=sqrt((sinkx-SN(i).x)^2 + (sinky-SN(i).y)^2); % calculates the distance between the sink and the cluster head
            CLheads=CLheads+1;	    % sum of cluster heads that have been elected 
            SN(i).cluster=CLheads;  % cluster of which the node got elected to be cluster head
            CL(CLheads).x=SN(i).x;  % X-axis coordinates of elected cluster head
            CL(CLheads).y=SN(i).y;  % Y-axis coordinates of elected cluster head
            CL(CLheads).id=i;      % Assigns the node ID of the newly elected cluster head to an array
            fprintf('Node %d elected as a cluster head for Level 3\n', i);
        end
   
end

% Grouping the Nodes into Clusters for Level 3
for i=level3
    if SN(i).role==0 && SN(i).E>0 && CLheads>0
        % if node is normal
        d=sqrt((CL(1).x-SN(i).x)^2 + (CL(1).y-SN(i).y)^2);
        % calculate the distance 'd' between the sensor node and the cluster head
        for m=2:CLheads
            % calculate distance between node and other cluster heads
            d_temp=sqrt((CL(m).x-SN(i).x)^2 + (CL(m).y-SN(i).y)^2);
            if d_temp < d
                % find the minimum distance to a cluster head
                d=d_temp;
            end
        end
        SN(i).cluster=3;  % assign node to level 3 cluster
        SN(i).dtch=d;     % assign the distance of node to CH
        SN(i).chid=3;     % assign level 3 cluster head ID
        fprintf('Node %d belongs to Level 3 cluster\n', i);
    end
end

% Data Aggregation Phase % for level3
for i=level3
    if SN(i).role==1 && SN(i).E>0
        % if node is cluster head
        ETx=(Eelec+EDA)*k + Eamp * k * SN(i).dts^2;
        % calculate energy dissipation for cluster head to transmit to sink
        SN(i).E=SN(i).E - ETx;  % update energy of cluster head
        energy=energy+ETx;      % update total energy
        if SN(i).E <= 0
            % if cluster head's energy depletes with transmission
            dead_nodes=dead_nodes +1;
            operating_nodes= operating_nodes - 1;
            SN(i).cond=0;
            SN(i).rop=rnd;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cluster Heads Election for Level 4
for i=level4
    %if (SN(i).x >= 0  && SN(i).x <=100) && (SN(i).y >= 0 && SN(i).y < 25)
        % check if the node is in level 4
        generate=rand;	
        if generate < t
            SN(i).role=1;	    % assigns the node role of a cluster head
            SN(i).rn=rnd;	    % Assigns the round that the cluster head was elected to the data table
            SN(i).tel=SN(i).tel + 1;   
            SN(i).rleft=1/p-tleft;    % rounds for which the node will be unable to become a CH
            SN(i).dts=sqrt((sinkx-SN(i).x)^2 + (sinky-SN(i).y)^2); % calculates the distance between the sink and the cluster head
            CLheads=CLheads+1;	    % sum of cluster heads that have been elected 
            SN(i).cluster=CLheads;  % cluster of which the node got elected to be cluster head
            CL(CLheads).x=SN(i).x;  % X-axis coordinates of elected cluster head
            CL(CLheads).y=SN(i).y;  % Y-axis coordinates of elected cluster head
            CL(CLheads).id=i;      % Assigns the node ID of the newly elected cluster head to an array
            fprintf('Node %d elected as a cluster head for Level 4\n', i);
        end
   
end

% Grouping the Nodes into Clusters for Level 4
for i=level4
    if SN(i).role==0 && SN(i).E>0 && CLheads>0
        % if node is normal
        d=sqrt((CL(1).x-SN(i).x)^2 + (CL(1).y-SN(i).y)^2);
        % calculate the distance 'd' between the sensor node and the cluster head
        for m=2:CLheads
            % calculate distance between node and other cluster heads
            d_temp=sqrt((CL(m).x-SN(i).x)^2 + (CL(m).y-SN(i).y)^2);
            if d_temp < d
                % find the minimum distance to a cluster head
                d=d_temp;
            end
        end
        SN(i).cluster=4;  % assign node to level 4 cluster
        SN(i).dtch=d;     % assign the distance of node to CH
        SN(i).chid=4;     % assign level 4 cluster head ID
        fprintf('Node %d belongs to Level 4 cluster\n', i);
    end
end

% Data Aggregation Phase % for level4
for i=level4
    if SN(i).role==1 && SN(i).E>0
        % if node is cluster head
        ETx=(Eelec+EDA)*k + Eamp * k * SN(i).dts^2;
        % calculate energy dissipation for cluster head to transmit to sink
        SN(i).E=SN(i).E - ETx;  % update energy of cluster head
        energy=energy+ETx;      % update total energy
        if SN(i).E <= 0
            % if cluster head's energy depletes with transmission
            dead_nodes=dead_nodes +1;
            operating_nodes= operating_nodes - 1;
            SN(i).cond=0;
            SN(i).rop=rnd;
        end
    end
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %% end of modification

                           %%%%%% Steady-State Phase %%%%%%
                      
  
% Energy Dissipation for normal nodes %
    
    for i=1:n
       if (SN(i).cond==1) && (SN(i).role==0) && (CLheads>0)
       	if SN(i).E>0
            ETx= Eelec*k + Eamp * k * SN(i).dtch^2;
            SN(i).E=SN(i).E - ETx;
            energy=energy+ETx;
            
        % Dissipation for cluster head during reception
        %if SN(SN(i).chid).E>0 && SN(SN(i).chid).cond==1 && SN(SN(i).chid).role==1
        if SN(i).E > 0 && SN(i).chid > 0 && SN(i).chid <= n && SN(SN(i).chid).E > 0 && SN(SN(i).chid).cond == 1 && SN(SN(i).chid).role == 1
            ERx=(Eelec+EDA)*k;
            energy=energy+ERx;
            SN(SN(i).chid).E=SN(SN(i).chid).E - ERx;
             if SN(SN(i).chid).E<=0  % if cluster heads energy depletes with reception
                SN(SN(i).chid).cond=0;
                SN(SN(i).chid).rop=rnd;
                dead_nodes=dead_nodes +1;
                operating_nodes= operating_nodes - 1
             end
        end
        end
        
        
        if SN(i).E<=0       % if nodes energy depletes with transmission
        dead_nodes=dead_nodes +1;
        operating_nodes= operating_nodes - 1
        SN(i).cond=0;
        SN(i).chid=0;
        SN(i).rop=rnd;
        end
        
      end
    end            
    
    
    
% Energy Dissipation for cluster head nodes %
   
   for i=1:n
     if (SN(i).cond==1)  && (SN(i).role==1)
         if SN(i).E>0
            ETx= (Eelec+EDA)*k + Eamp * k * SN(i).dts^2;
            SN(i).E=SN(i).E - ETx;
            energy=energy+ETx;
         end
         if  SN(i).E<=0     % if cluster heads energy depletes with transmission
         dead_nodes=dead_nodes +1;
         operating_nodes= operating_nodes - 1
         SN(i).cond=0;
         SN(i).rop=rnd;
         end
     end
   end

   

  
    if operating_nodes<n && temp_val==0
        temp_val=1;
        flag1stdead=rnd
    end
    % Display Number of Cluster Heads of this round %
    %CLheads;
   
    
    transmissions=transmissions+1;
    if CLheads==0
    transmissions=transmissions-1;
    end
    
 
    % Next Round %
    rnd= rnd +1;
    
    tr(transmissions)=operating_nodes;
    op(rnd)=operating_nodes;
    

    if energy>0
    nrg(transmissions)=energy;
    end
    

end


sum=0;
for i=1:flag1stdead
    sum=nrg(i) + sum;
end

temp1=sum/flag1stdead;
temp2=temp1/n;

for i=1:flag1stdead
avg_node(i)=temp2;
end
    
    % Plotting Simulation Results "Operating Nodes per Round" %
    figure(2)
    plot(1:rnd,op(1:rnd),'-r','Linewidth',2);
    title ({'LEACH'; 'Operating Nodes per Round';})
    xlabel 'Rounds';
    ylabel 'Operational Nodes';
    hold on;
    
    % Plotting Simulation Results  %
    figure(3)
    plot(1:transmissions,tr(1:transmissions),'-r','Linewidth',2);
    title ({'LEACH'; 'Operational Nodes per Transmission';})
    xlabel 'Transmissions';
    ylabel 'Operational Nodes';
    hold on;
    
    % Plotting Simulation Results  %
    figure(4)
    plot(1:flag1stdead,nrg(1:flag1stdead),'-r','Linewidth',2);
    title ({'LEACH'; 'Energy consumed per Transmission';})
    xlabel 'Transmission';
    ylabel 'Energy ( J )';
    hold on;
    

    % Plotting Simulation Results  %
    figure(5)
    plot(1:flag1stdead,avg_node(1:flag1stdead),'-r','Linewidth',2);
    title ({'LEACH'; 'Average Energy consumed by a Node per Transmission';})
    xlabel 'Transmissions';
    ylabel 'Energy ( J )';
    hold on;
  
  