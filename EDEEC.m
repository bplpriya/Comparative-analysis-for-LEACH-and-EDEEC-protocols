

clear all
xm=100;
ym=100;
sink.x=0.5*xm;  %location of sink on x-axis
sink.y=0.5*ym;  %location of sink on y-axis
n=100  %nodes+
P=0.1;  %probability of cluster heads
Eo=0.5;%initial energy
%
Echeck=Eo;
%
ETX=50*0.000000001;  %tx energy
ERX=50*0.000000001;  %rx energy
Efs=10*0.000000000001;  %free space loss
Emp=0.0013*0.000000000001;   %multipath loss
%Data Aggregation Energy
EDA=5*0.000000001;  %compression energy
a=1.5;   %fraction of energy enhancment of advance nodes
rmax=5000  %maximum number of rounds
do=sqrt(Efs/Emp);  %distance do is measured
Et=0;  %variable just use below
m=0.5;
mo=0.4;
b=3;
normal=n*(1-m);
advance=n*m*(1-mo);
monaysuper=n*m*mo;
for i=1:1:monaysuper
S(i).xd=rand(1,1)*xm;  %generates a random no. use to randomly distibutes nodes on x axis
XR(i)=S(i).xd;
S(i).yd=rand(1,1)*ym;  %generates a random no. use to randomly distibutes nodes on y axis
YR(i)=S(i).yd;
S(i).G=0; %node is elegible to become cluster head
%talhar=rand*a
S(i).E=Eo*(1+b);
%S(i).A=talhar;
E(i)= S(i).E;
%     if (E(i)>Echeck)
%         m1=m1+1;
%     end
Et=Et+E(i);  %estimating total energy of the network
%initially there are no cluster heads only nodes
S(i).type='N';
end
talha1=monaysuper+advance;
for i=monaysuper:1:talha1
S(i).xd=rand(1,1)*xm;  %generates a random no. use to randomly distibutes nodes on x axis
XR(i)=S(i).xd;
S(i).yd=rand(1,1)*ym;  %generates a random no. use to randomly distibutes nodes on y axis
YR(i)=S(i).yd;
S(i).G=0; %node is elegible to become cluster head
%talhar=rand*a
S(i).E=Eo*(1+a);
%S(i).A=talhar;
E(i)= S(i).E;
%     if (E(i)>Echeck)
%         m1=m1+1;
%     end
Et=Et+E(i);  %estimating total energy of the network
%initially there are no cluster heads only nodes
S(i).type='N';
end
for i=talha1:1:n
S(i).xd=rand(1,1)*xm;  %generates a random no. use to randomly distibutes nodes on x axis
XR(i)=S(i).xd;
S(i).yd=rand(1,1)*ym;  %generates a random no. use to randomly distibutes nodes on y axis
YR(i)=S(i).yd;
S(i).G=0; %node is elegible to become cluster head
%talhar=rand*a
S(i).E=Eo;
%S(i).A=talhar;
E(i)= S(i).E;
%     if (E(i)>Echeck)
%         m1=m1+1;
%     end
Et=Et+E(i);  %estimating total energy of the network
%initially there are no cluster heads only nodes
S(i).type='N';
end




d1=0.765*xm/2;  %distance between cluster head and base station
K=sqrt(0.5*n*do/pi)*xm/d1^2; %optimal no. of cluster heads
d2=xm/sqrt(2*pi*K);  %distance between cluster members and cluster head
Er=4000*(2*n*ETX+n*EDA+K*Emp*d1^4+n*Efs*d2^2);  %energy desipated in a round
S(n+1).xd=sink.x; %sink is a n+1 node, x-axis postion of a node
S(n+1).yd=sink.y; %sink is a n+1 node, y-axis postion of a node
countCHs=0;  %variable, counts the cluster head
cluster=1;  %cluster is initialized as 1
flag_first_dead=0; %flag tells the first node dead
flag_teenth_dead=0;  %flag tells the 10th node dead
flag_all_dead=0;  %flag tells all nodes dead
dead=0;  %dead nodes count initialized to 0
first_dead=0;
teenth_dead=0;
all_dead=0;
allive=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
for r=0:1:rmax
r
if(mod(r, round(1/P) )==0)
for i=1:1:n
S(i).G=0;
S(i).cl=0;
end
end
Ea=Et*(1-r/rmax)/n;
dead=0;
for i=1:1:n

if (S(i).E<=0)
dead=dead+1;
if (dead==1)
if(flag_first_dead==0)
first_dead=r;
flag_first_dead=1;
end
end
if(dead==0.1*n)
if(flag_teenth_dead==0)
teenth_dead=r;
flag_teenth_dead=1;
end
end
if(dead==n)
if(flag_all_dead==0)
all_dead=r;
flag_all_dead=1;
end
end
end
if S(i).E>0
S(i).type='N';
end
end

STATISTICS.DEAD(r+1)=dead;
STATISTICS.ALLIVE(r+1)=allive-dead;
countCHs=0;
cluster=1;
for i=1:1:n

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Ea>0
if (S(i).E<= Eo)
p(i)=P*E(i)/(1+m*(a+mo*b))*Ea;
end

if (S(i).E<=Eo*(1+a))
p(i)=P*(1+a)*E(i)/(1+m*(a+mo*b))*Ea;
end
if (S(i).E<=Eo*(1+b))
p(i)=P*(1+b)*E(i)/(1+m*(a+mo*b))*Ea;
end

%      if (S(i).E<= (b*Eo))
%          p(i)=c*(1+S(i).A)*P*E(i)/((1+(S(i).A*m))*Ea);
%      end
%
%      if (S(i).E<=Eo)
%          p(i)=P*E(i)/((1+(S(i).A*m))*Ea);
%      end
%      if (S(i).E>Eo)
%          p(i)=(1+S(i).A)*P*E(i)/((1+(S(i).A*m))*Ea);
%      end
%
%%p(i)=P*n*S(i).E*E(i)/(Et*Ea);
if(S(i).E>0)
temp_rand=rand;
if ( (S(i).G)<=0)
if(temp_rand<= (p(i)/(1-(p(i)*(r*mod(1,p(i)))))))
countCHs=countCHs+1;
packets_TO_BS=packets_TO_BS+1;
PACKETS_TO_BS(r+1)=packets_TO_BS;
S(i).type='C';
S(i).G=round(1/p(i))-1;
C(cluster).xd=S(i).xd;
C(cluster).yd=S(i).yd;
distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
C(cluster).distance=distance;
C(cluster).id=i;
X(cluster)=S(i).xd;
Y(cluster)=S(i).yd;
cluster=cluster+1;
distance;
if (distance>do)
S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
end
if (distance<=do)
S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
end
end

end

end
end
end
STATISTICS.COUNTCHS(r+1)=countCHs;
%(5)´ØÄÚ³ÉÔ±Ñ¡Ôñ´ØÍ·Ä£¿é(¼´´ØµÄÐÎ³ÉÄ£¿é)
%´ØÄÚ³ÉÔ±¶Ô´ØÍ·µÄÑ¡Ôñ£¨¼´´ØµÄÐÎ³É£©Ëã·¨
for i=1:1:n
if ( S(i).type=='N' && S(i).E>0 )
if(cluster-1>=1)
min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
min_dis_cluster=0;
for c=1:1:cluster-1
temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
if ( temp<min_dis )
min_dis=temp;
min_dis_cluster=c;
end
end
%´ØÄÚ½Úµã£¨·¢ËÍ4000bitÊý¾Ý£©ÄÜÁ¿ÏûºÄ
if(min_dis_cluster~=0)
min_dis;
if (min_dis>do)
S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis));
end
if (min_dis<=do)
S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
end

S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 );
packets_TO_CH=packets_TO_CH+1;
else
min_dis;
if (min_dis>do)
S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis));
end
if (min_dis<=do)
S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
end
packets_TO_BS=packets_TO_BS+1;

end
S(i).min_dis=min_dis;
S(i).min_dis_cluster=min_dis_cluster;
else
min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
if (min_dis>do)
S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis));
end
if (min_dis<=do)
S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
end
packets_TO_BS=packets_TO_BS+1;
end
end
end
STATISTICS.PACKETS_TO_CH(r+1)=packets_TO_CH;
STATISTICS.PACKETS_TO_BS(r+1)=packets_TO_BS;
end
first_dead
teenth_dead
all_dead
STATISTICS.DEAD(r+1)
STATISTICS.ALLIVE(r+1)
STATISTICS.PACKETS_TO_CH(r+1)
STATISTICS.PACKETS_TO_BS(r+1)
STATISTICS.COUNTCHS(r+1)
r=0:5000;
subplot(2,2,1);
plot(r,STATISTICS.DEAD);
subplot(2,2,2);
plot(r,STATISTICS.ALLIVE);
subplot(2,2,3);
plot(r,STATISTICS.PACKETS_TO_BS);
subplot(2,2,4);
plot(r,STATISTICS.COUNTCHS);
