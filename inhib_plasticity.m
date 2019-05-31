%% Initialization:
clear all
close all
%% Figure dimensions
figure(1)
set(gcf,'PaperUnits','centimeters')
%Setting the units of the figure on paper to centimeters.
xSize = 20; ySize = 8;
%Size of the figure
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
%Coordinates to center the figure on A4-paper
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
%This command sets the position and size of the figure on the paper to the desired values.
set(gcf,'Position',[0.5 0.5 xSize*50 ySize*50])
set(gcf, 'Color', 'w');
%% Parameters
eta=0.001;
%Learning Rate
alpha=0.25*eta;
%Decay term
tauPlasticity=20;
%% page 35
%Time window of the learning rule.
duration=1500;
%Duration of a single run in ms.
dt=0.1;
% Simulation time step in ms.
NRuns=100;
%Number of consecutive runs in one simulation.
%(The total duration of the simulation is thus \Duration*NRuns".
tRef=5;
% Refractory period for the spike trains.
gBarEx=0.014;
% Scaling factor of the excitatory synaptic conductance in units of the leak (gleak = 10 nS)
gBarIn=0.035;
% Scaling factor of the inhibitory synaptic conductance in units of the leak (gleak = 10 nS)
VRest=-60;
%Resting potential in mV.
Vth=-50;
% Threshold, in mV.
taumem=20;
% Membrane time constant, in ms.
EAMPA=0;
% Excitatory reversal potential, in mV.
EGABA=-80;
% Inhibitory reversal potential, in mV.
tauEx=5;
% Excitatory synaptic time constant, in ms.
tauIn=10;
% Inhibitory synaptic time constant, in ms.
stopnow=0;
% This is a flag to stop the while loop in case the EPSP routine below is uncommented
%% Input
noisetau=50;
% Filter time for the input in ms
Backgroundrate=5*dt/1000;
% Background spiking probability 5Hz*dt (= 0.0005) (see Technical Appendix 6.2.1).
ApproximateNormalizingFactor=0.03;
% This serves to normalize the trace to an approximate peak value of 1
Maxrate=500*dt/1000;
% Peak spiking probability 500Hz*dt (= 0.05))
NSigs=8;
% Number of input signals.
NCells=1000;
% Number of input spike trains.
ExFrac=0.8;
%Fraction of Excitatory spike trains.
%% page 36
ExGroupsize=(NCells*ExFrac)/NSigs;
% Number of spike trains in each excitatory input group
InGroupsize=round((NCells*(1-ExFrac))/NSigs);
% Number of spike trains in each inhibitory input group
expGEx=exp(-dt/tauEx);
expGIn=exp(-dt/tauIn);
expPlasticity=exp(-dt/tauPlasticity);
expnoise=exp(-dt/noisetau);
% Pre-calculated exponential factors for speed of simulation
%% Vectors
Input=zeros(1,NSigs);
% A vector that holds the momentary input signal for each signal channel
Timevector=(0.1:dt:duration);
% A vector of time in ms, in steps of dt.
Exkeep=Timevector*0;
Inkeep=Timevector*0;
% Vectors to keep the synaptic currents.
FilteredWhiteNoise=zeros(1,8);
% A vector to create the time-filtered input
InputGroup=zeros(1,NCells);
% A Vector that keeps track of which spike train belongs to which input
InputSpikeRefr=zeros(1,NCells);
% To keep track of the input spike train refractory periods)
tolos=0;
% (== t(ime)o(f)l(ast)(o)utput (s)pike (Keeping track of the output cell's refractory period )
Synapse=ones(1,NCells);
%Synaptic weights
sgEx=zeros(1,NCells);
sgIn=zeros(1,NCells);
% Vectors to save the group-wise synaptic conductances the cell experiences.
AveExCurr=zeros(1,NSigs);
AveInCurr=zeros(1,NSigs);
% Vectors to save the group-wise synaptic currents the cell experiences.
pre=zeros(1,NCells);
% A Vector to save the presynaptic learning trace.
Time=(1:NRuns);
Rate=zeros(1,NRuns);
%Vectors for plotting
%% InputGroups
temptype=0;
for i=1:NCells
if (i<=NCells*ExFrac)
if (mod(i-1,ExGroupsize)==0)
temptype=temptype+1;
end
%% page 37
InputGroup(i)= temptype;
else
if (mod(i,InGroupsize)==0)
temptype=temptype-1;
end
InputGroup(i)= -temptype;
end
end
InputGroup(1000)=InputGroup(999);
% This routine assigns every spike train to an input group, starting with group 1 and ending
% with group NSIgs for the excitatory spike trains and then going back from NSigs to 1 for the
% inhibitory spiketrains. To label inhibitory spike trains uniquely, their group identity is assigned
% as a negative number.
%% Synapse Tuning
for i=1:800
Synapse(i) = 0.3 + (1.1/(1+(InputGroup(i)-5)^4))+rand*0.1;
end
for i=801:NCells
Synapse(i) = 0.1;
end
% This routine assigns a synaptic weight to every synapse. Excitatory synapses are tuned
% according to their group identity (plus a noise term) to resemble the tuning reported in (7).
% Inhibitory synapses are uniformly weak.
%%%%%%%%%%%%%%%%%%
%% Simulation%%%%%
%%%%%%%%%%%%%%%%%%
%% Initialize values
gEx=0;
% Overall excitatory synaptic conductance.
gIn=0;
% Overall inhibitory synaptic cconductance.
gLeak=1;
% Leak conductance (Everything else is normalized in respect to the leak.)
V(1)=VRest;
% Membrane potential is initially at VRest.
post=0;
% Postsynaptic learning trace.
InputSpikeCount=0;
OutputSpikeCount=0;
AveCurrCounter=0;
%Counters for rates and averaging.
tRunning=0;
runcount=0;
going=1;
%% page 38
% Time parameters are initialized
%% Start of the loops
while(going>0)
runcount=runcount+1;
% A counter to keep track of consecutive runs.
OutputSpikeCount=0;
for i=1:NSigs
AveExCurr(i)=0;
AveInCurr(i)=0;
end
AveCurrCounter=0;
% Reset the counters.
% %Uncomment from here to %XX to plot individual synapse strengths in [pS]
% figure(2)
% clf(2)
% plot(gBarEx*Synapse(1:800)*10000, 'r.')
% hold on
% plot(-(-796:4:1),gBarIn*Synapse(801:NCells)*10000, 'g.')
% %%XX end of commented section
% Uncommenting this routine plots the synaptic weights in a second window.
for t=2:length(Timevector)
% The time loop begins here.
tRunning=tRunning+dt;
gEx = gEx*expGEx;
gIn = gIn*expGIn;
for i=1:NSigs
sgEx(i) = sgEx(i)*expGEx;
sgIn(i) = sgIn(i)*expGIn;
end
% The synaptic conductances decay exponentially towards 0.
for i=801:NCells
pre(i)= pre(i)*expPlasticity;
end
post=post*expPlasticity;
% The learning traces decay exponentially towards 0.
%% Input Signals
for i=1:NSigs
re=rand-0.5;
FilteredWhiteNoise(i) = re -(re - FilteredWhiteNoise(i))*expnoise;
Input(i)=Backgroundrate + ...
max(0, Maxrate*FilteredWhiteNoise(i))/ApproximateNormalizingFactor;
%% page 39
end
% At every time step the current input for each signal is calculated from a time-filtered
% input signal. Signals traces are not saved to increase simulation speed.
% For more details, please see 6.2.1
% % % Uncomment from here to %%%YY to plot single EPSPs at VRest
% for i=1:NSigs
% Input(i)=0;
% end
% if(t==300/dt)
% gEx = gEx + gBarEx;
% end
% if(t==700/dt)
% gIn = gIn + gBarIn;
% end
% stopnow=1;
%%%%%%YY
% Uncommenting this routine makes it possible to evaluate single EPSPs by switching off the
% input signals and instead injecting one excitatory and one inhibitory PSP (at 300 and 700 ms)
% into the cell. \stopnow" is a flag that stops the simulation after one run and rescales the axes
% of the voltage plot.
%% Presynaptic spike trains
for i=1:NCells
if (rand < Input(abs(InputGroup(i))) && InputSpikeRefr(i)<=0)
% If a Spiketrain fired a spike, (see also 6.2.1) ...
if(InputGroup(i)>0)
% .. and if it is excitatory...
gEx = gEx + (gBarEx * Synapse(i));
% ... increase the excitatory synaptic conductance variable according to
% the strength of the synapse
sgEx(abs(InputGroup(i)))=sgEx(abs(InputGroup(i))) + gBarEx*Synapse(i);
% (Keeping track of the synaptic conductances group-wise for plotting.)
else
% otherwise (meaning the synapse is inhibitory)...
gIn = gIn + gBarIn * Synapse(i);
% ... increase the synaptic conductance by the synapse strength
sgIn(abs(InputGroup(i)))=sgIn(abs(InputGroup(i))) + gBarIn*Synapse(i);
% (To keep track of the synaptic conductances group-wise for plotting.)
pre(i)= pre(i) + eta;
% Update the presynaptic learning trace.
Synapse(i)=Synapse(i) + post - alpha;
% Update the synaptic strength according to the rule.
% ! add the effect of proximate post synaptic spikes, and subtract ff.
if(Synapse(i) <=0)
Synapse(i)=0;
%% page 40
end
% Ensure non-negative synaptic weights.
end
InputSpikeRefr(i)=tRef;
% ... Also: set the refractory time counter to the value of refractoriness
InputSpikeCount=InputSpikeCount+1;
% .. and count the overall number of input spikes
else
% meaning if no presynaptic spike was evoked:
InputSpikeRefr(i)=InputSpikeRefr(i)-dt;
% subtract dt from the refractory counter.
end
end
%% Membrane potential and postsynaptic spikes.
if ((tRunning - tolos) < tRef)
V(t) = VRest;
% If the Cell is refractory, keep V at Vrest
else
% Meaning: if the cell is not refractory, ...
gTot = gLeak + gEx + gIn;
% calculate the total membrane conductance,
tauEff=taumem/gTot;
% and the effective time constant, as well as...
VInf = ((gLeak*VRest + gEx * EAMPA+ gIn*EGABA)/gTot);
% the membrane potential that V strives towards.
V(t) = VInf + (V(t-1) - VInf)*exp(-dt/tauEff);
% Use the above to update the membrane potential
for i=1:NSigs
AveExCurr(i)= AveExCurr(i) + sgEx(i)*(V(t)-EAMPA);
AveInCurr(i)= AveInCurr(i) + ...
sgIn(i)*(V(t)-EGABA) + (gLeak*(V(t)-VRest))/NSigs;
end
AveCurrCounter=AveCurrCounter+1;
% The above routine keeps track for the group-wise input currents for plotting but
% does not affect the behavior of the cell. We divide the (mostly inhibitory acting)
% leak current evenly to all Groups,since each input signal causes the same absolute
% amount of leak current (by deflecting the membrane potential away from rest with
% identical statistics over time for each signal).
end
if (V(t)>Vth)
% If the membrane potential hits threshold...
tolos=tRunning;
% ... set the refractory counter to the current time step
V(t-1)=0;
% ... set the last membrane potential before the spike to zero (for plotting)
V(t)=VRest;
%% page 41
% ... reset the current membrane potential to Vrest.
OutputSpikeCount=OutputSpikeCount+1;
% ... count the spike.
post = post + eta;
% ... update the postsynaptic learning trace
for i=801:NCells
Synapse(i) = Synapse(i)+pre (i);
end
% update all synapses according to the rule
% ! Add the effect of proximate presynaptic spikes.
end
% subtract dt from the refractory counter.
Exkeep(t)=gEx*(V(t)-EAMPA);
Inkeep(t)=gIn*(V(t)-EGABA);
% For plotting purposes, keeps track of the excitatory and inhibitory synaptic currents.
% Because everything is normalized by the leak conductance, the values are saved in
% units of 10nS  mV (= 10?11 Amp = [10pA]).
% To plot in units of nA, one has to divide this number by 100.
end
% End of the time loop
Rate(runcount)=OutputSpikeCount/duration*1000;
Time(runcount)=runcount*duration/1000;
%%%%%%%%%%%%%%%%%%
% Plotting %%%%%%%
%%%%%%%%%%%%%%%%%%
figure(1)
subplot(3,14,[1 2 3 4 15 16 17 18])
hold off
plot(Timevector, Exkeep/100,'k', 'LineWidth', 1)
% plotting in units of nA hence the values need to be divided by 100, see above.
hold on
plot(Timevector,Inkeep/100,'Color', [0.6 0.6 0.6], 'LineWidth', 1)
plot(Timevector,(Inkeep+Exkeep)/100,'Color', [0 0.6 0])
axis([0 duration -3 3]);
ylabel('Synaptic Currents [nA]');
subplot(3,14,[29 30 31 32])
if (stopnow==0)
plot(Timevector, V, 'k')
end
axis([0 duration -65 5]);
ylabel('Mem. Pot. [mV]');
xlabel('Time [ms]');
%% page 42
subplot(3,14, [6 7 8 9 20 21 22 23 34 35 36 37])
hold off
plot((1:8), (-AveExCurr/AveCurrCounter)/100,'-ks','LineWidth',2,...
'MarkerEdgeColor','k',...
'MarkerFaceColor','k',...
'MarkerSize',10)
hold on
plot((1:8), (AveInCurr/AveCurrCounter)/100,'-ks','LineWidth',2,...
'MarkerEdgeColor','k',...
'MarkerFaceColor','w',...
'MarkerSize',10)
axis([1 8 0 0.25]);
ylabel('Mean Synaptic Currents [nA]');
xlabel('Signal Number');
subplot(3,14,[11 12 13 14 25 26 27 28 39 40 41 42])
hold on
plot(runcount*duration/1000, mean(Synapse(925:949))*gBarIn*10000, ...
'.', 'Color','r')
% (Plotted in units of pS)
plot(runcount*duration/1000, mean(Synapse(900:924))*gBarIn*10000,...
'.', 'Color',[0.5 0.0 0.5])
% (Plotted in units of pS)
plot(runcount*duration/1000, mean(Synapse(875:899))*gBarIn*10000, ...
'.', 'Color','b')
% (Plotted in units of pS)
ylabel('Mean Synaptic Strength [pS]');
axis([0 NRuns*duration/1000 0 800]);
ax1 = gca;
set(ax1,'XColor','k','YColor',[0.5 0.0 0.5], 'YAxisLocation','right')
ax2 = axes('Position',get(ax1,'Position'),...
'XAxisLocation','bottom',...
'YAxisLocation','left',...
'Color','none',...
'XColor','k' ,'YColor','k');
hold(ax2,'on')
plot(Time(1:runcount), Rate(1:runcount), '.-k','Parent',ax2)
ylabel('Mean Output Rate [Hz]','Parent',ax2);
xlabel('Time [s]');
axis([0 NRuns*duration/1000 0 100]);
if(runcount>NRuns)
going=0;
end
if (stopnow==1)
going=0;
end
% (This is only ==1 when the EPSP routine is uncommented)
%% page 43
end
% End of the while loop
% End of Code
