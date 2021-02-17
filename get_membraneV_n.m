function [V,CV] = get_membraneV_n(v_init,i,T,figon,M,sigma,visible)
% get_membraneV - Implements a function that returns the membrane potential
% for arbitrary input currents i(t)
% Syntax:  [output1,output2] = get_membraneV (input1,input2,input3)
%
% Inputs:
%    v_init - initial V
%    i      - Input Current (mean)
%    T      - The total simulation time 
%    figon  - Determines the plot
%    M      - Matrix of connected neurons
%    sigma  - Standard deviation of input current (Optional)
%    visible -  visible range in the raster plot (ms)
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
% Created on    : 14-Feb-2021   
% Last revision : NONE

%               ------------- BEGIN CODE --------------
% T=1000; % total time in ms (Simulating 1s or 1000ms)
% i=0.9;
% v_init=-65;

dt=0.1; %step size is 0.1ms
v_reset=-70; %mV %70
v_rest=-65; %mV %65
v_thresh=-50; %mV
EPSP_delay=1.5; %ms
n=size(M,1);
% M=2* rand(n) .* ~eye(n); % connection matrix



spiking_voltage=40; %mv

tau=50; %ms

time=0:dt:T;

if isempty(sigma)
    I=ones(size(time)).*i; 
else
    I=normrnd(i,sigma,size(time));
%     hist(I)
end


V=NaN(length(time),n);

tol=0.001;

spike_timestamp={[]}; %=zeros(length(time),n);
for j=1:n
    spike_timestamp{j}=[];
end

for t=1:length(time)
    for j=1:n
        
    if time(t)==0
        V(t,j)=v_init;
    elseif time(t)>0 && time(t)<=T
        % if the membrane potential is around the spiking voltage, then reset
        if abs(V(t-1,j) - spiking_voltage) < tol 
            V(t,j) = v_reset;
        else
            if j==1
                V(t,j) = V(t-1,j) + dt* ((-1/tau) * (V(t-1,j)-v_rest) + I(t-1));
            else
                 V(t,j) = V(t-1,j) + dt* ((-1/tau) * (V(t-1,j)-v_rest));
            end
        end
        
        % Make sure the connection matrix gets addressed at this point
        if time(t)>EPSP_delay
           for k=1:n
                if abs(V(t-int32(EPSP_delay/dt),k) - spiking_voltage) < tol % if there was a spike before 1.5ms from somewhere k
                    V(t,j) = V(t,j)+ M(k,j); % M(k,j) spike from k, resulting voltage increment in j
                end
           end
            
        end
        
        % if voltage is greater than  threshold but did not reach the 
        % spiking potential then spike will occur
        if (V(t,j) >= v_thresh) && abs(V(t,j) - spiking_voltage) > tol  
                V(t,j) = spiking_voltage;
                spike_timestamp{j}=[spike_timestamp{j},time(t)];   %(t,j)=time(t); 
                  
        elseif V(t,j) <= v_reset % if below reset potential, then stays at 
            % reset potential
            V(t,j)=v_reset; 
        end
    end
    end
end

ISI={};
CV=zeros(1,n);
for j=1:n
    ISI{j}=diff(spike_timestamp{j});
    CV(j)=std(ISI{j})/mean(ISI{j});
end
% CV=mean(ISI)/std(ISI); %Coefficient of Variation

% calculate mean ISI (Inter Spike Interval) 
% => mean(diff(spike_timestamp))
% to find the distribution of the ISI => diff(spike_timestamp)
if strcmp(figon,'on')
%     figure, 
    % Raster plot
    K=V;
    K(V<0)=0;
    K(V>30)=1;
    figure,
    subplot(1,2,1)
    for p=1:size(V,2)
        l=plot(0:0.1:T,K(:,p)'.*p,'o');
        l.MarkerFaceColor = l.Color;
        hold on
    end
hold off
%     h=subplot(1,2,1),
%     sampling_rate=10000;
%     rasterplot(find(K),10,visible./dt,h,sampling_rate);
    xlabel('Time')
    ytickformat('%d')
    yticks(1:n);
    ylim([0.5,n+0.5])
    xlabel('Time (ms)','FontName','Times New Roman','FontSize',16,'FontWeight','normal');
    ylabel('Neuron #','FontName','Times New Roman','FontSize',16,'FontWeight','normal')
    title('Raster Plot','FontName','Times New Roman','FontSize',16,'FontWeight','bold');
    set(gca,'LineWidth',1.5)
    set(gca,'TickLength',[0 0])
    hold off

    subplot(1,2,2)
    histogram(CV,'Normalization','pdf');hold on;
    [f1,xi1]=ksdensity(CV);
    plot(xi1,f1,'r','LineWidth',1.5,'DisplayName','ksdensity');
    xlabel('Time Bin (ms)','FontName','Times New Roman','FontSize',16,'FontWeight','normal');
    ylabel('Histogram','FontName','Times New Roman','FontSize',16,'FontWeight','normal')
    title('{\it C_{v}} Distribution','FontName','Times New Roman','FontSize',16,'FontWeight','bold');
    legend;
    set(gca,'LineWidth',1.5)

%     subplot(1,3,1)
%     histogram(I,'Normalization','pdf');hold on;
%     [f1,xi1]=ksdensity(I);
%     plot(xi1,f1,'r','LineWidth',1.5,'DisplayName','ksdensity');
%     xlabel('Time Bin (ms)','FontName','Times New Roman','FontSize',16,'FontWeight','normal');
%     ylabel('Histogram','FontName','Times New Roman','FontSize',16,'FontWeight','normal')
%     title('I','FontName','Times New Roman','FontSize',16,'FontWeight','bold');
%     legend;
%     set(gca,'LineWidth',1.5)
%     
%     subplot(1,3,2)
%     plot(time,V,'LineWidth',1.2);
%     xlabel('Time (ms)','FontName','Times New Roman','FontSize',16,'FontWeight','normal');
%     ylabel('Membrane Potential (mV)','FontName','Times New Roman','FontSize',16,'FontWeight','normal')
%     title('V(t)','FontName','Times New Roman','FontSize',16,'FontWeight','bold');
%     ylim([v_reset-10 spiking_voltage+10])
%     set(gca,'LineWidth',1.5)
%     subplot(1,3,3),
% %     plot(spike_timestamp(2:end),ISI);
%     histogram(ISI,'Normalization','pdf','DisplayName','hist');hold on;
%     [f,xi]=ksdensity(ISI);
%     plot(xi,f,'r','LineWidth',1.5,'DisplayName','ksdensity');
%     xlabel('Time Bins (ms)','FontName','Times New Roman','FontSize',16,'FontWeight','normal');
%     ylabel('Histogram (pdf normalized)','FontName','Times New Roman','FontSize',16,'FontWeight','normal')
%     title(['ISI Histogram, CV=' num2str(CV)],'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
%     set(gca,'LineWidth',1.5);legend;hold off
end
%               ------------- END OF CODE -------------
end
