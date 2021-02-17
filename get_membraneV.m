function [V,ISI,CV] = get_membraneV (v_init,i,T,figon,sigma)
% get_membraneV - Implements a function that returns the membrane potential
% for arbitrary input currents i(t)
% Syntax:  [output1,output2] = get_membraneV (input1,input2,input3)
%
% Inputs:
%    v_init - initial V
%    i      - Input Current (mean)
%    T      - The total simulation time 
%    figon  - Determines the plot
%    sigma  - Standard deviation of input current (Optional)
%
% Outputs:
%    V      - Membrane_potential response to the input i % Note i is not
%    exactly the current; there's simplification of constants there.
%    ISI    - Interspike Interval
%    CV     - Coeffiecient of Variation
% Example: 
%    Line 1 of example
%    Line 2 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author        : Abrar Faiyaz (ORCID iD: https://orcid.org/0000-0003-0666-3051 )
% email         : abrarfaiyaz.iutcse@gmail.com
% Created on    : 11-Feb-2021   
% Last revision : NONE

%               ------------- BEGIN CODE --------------
% T=1000; % total time in ms (Simulating 1s or 1000ms)
% i=0.9;
% v_init=-65;

dt=0.1; %step size is 0.1ms
v_reset=-70; %mV
v_rest=-65; %mV
v_thresh=-50; %mV

spiking_voltage=40; %mv

tau=20; %ms

time=0:dt:T;

if nargin<=4
    I=ones(size(time)).*i; 
elseif nargin==5

    I=normrnd(i,sigma,size(time));
%     hist(I)
end


V=NaN(size(time));

tol=0.001;

spike_timestamp=[];

for t=1:length(time)
    if time(t)==0
        V(t)=v_init;
    elseif time(t)>0 && time(t)<=T
        % if the membrane potential is around the spiking voltage, then reset
        if abs(V(t-1) - spiking_voltage) < tol 
            V(t) = v_reset;
        else
            V(t) = V(t-1) + dt* ((-1/tau) * (V(t-1)-v_rest) + I(t-1));
        end
        
        % if voltage is greater than  threshold but did not reach the 
        % spiking potential then spike will occur
        if (V(t) >= v_thresh) && abs(V(t) - spiking_voltage) > tol   
            V(t) = spiking_voltage;
            spike_timestamp=[spike_timestamp;time(t)];      
%         elseif abs(V(t) - spiking_voltage) < tol % if the membrane potential 
%             %  is around the spiking voltage, then reset
%             V(t) = v_reset;
        elseif V(t) <= v_reset % if below reset potential, then stays at 
            % reset potential
            V(t)=v_reset; 
        end
        
    end
end

ISI=diff(spike_timestamp);
CV=std(ISI)/mean(ISI); %Coefficient of Variation

% calculate mean ISI (Inter Spike Interval) 
% => mean(diff(spike_timestamp))
% to find the distribution of the ISI => diff(spike_timestamp)
if strcmp(figon,'on')
    figure, 
    subplot(1,3,1)
    histogram(I,'Normalization','pdf');hold on;
    [f1,xi1]=ksdensity(I);
    plot(xi1,f1,'r','LineWidth',1.5,'DisplayName','ksdensity');
    xlabel('Time Bin (ms)','FontName','Times New Roman','FontSize',16,'FontWeight','normal');
    ylabel('Histogram','FontName','Times New Roman','FontSize',16,'FontWeight','normal')
    title('I','FontName','Times New Roman','FontSize',16,'FontWeight','bold');
    legend;
    set(gca,'LineWidth',1.5)
    
    subplot(1,3,2)
    plot(time,V,'LineWidth',1.2);
    xlabel('Time (ms)','FontName','Times New Roman','FontSize',16,'FontWeight','normal');
    ylabel('Membrane Potential (mV)','FontName','Times New Roman','FontSize',16,'FontWeight','normal')
    title('V(t)','FontName','Times New Roman','FontSize',16,'FontWeight','bold');
    ylim([v_reset-10 spiking_voltage+10])
    set(gca,'LineWidth',1.5)
    subplot(1,3,3),
%     plot(spike_timestamp(2:end),ISI);
    histogram(ISI,'Normalization','pdf','DisplayName','hist');hold on;
    [f,xi]=ksdensity(ISI);
    plot(xi,f,'r','LineWidth',1.5,'DisplayName','ksdensity');
    xlabel('Time Bins (ms)','FontName','Times New Roman','FontSize',16,'FontWeight','normal');
    ylabel('Histogram (pdf normalized)','FontName','Times New Roman','FontSize',16,'FontWeight','normal')
    title(['ISI Histogram, CV=' num2str(CV)],'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
    set(gca,'LineWidth',1.5);legend;hold off
end
%               ------------- END OF CODE -------------
end
