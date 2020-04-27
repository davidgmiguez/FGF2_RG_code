function growth_rate_limit_12(~)
clear all;close all; clc;
%  data a cell of 1 arrow and t columns
% add saving data to a cell variable
% 07/07/2016 add 1 hour for the length of M
% 19/03/2019 modify to test data for RG

rng('shuffle')
total_time=48;time_step=1; % hours

number_of_simulations=30; 
all_data=cell(number_of_simulations,total_time/time_step); % aqui se guarda toda la información

for simulation_index=1:number_of_simulations
    clear nuclei
    [pp_dd,cell_cycle_length,cycling_cell_rate]=dynamics_data(1); % set dynamics at this time point

    simulation_index
    number_of_nuclei=round(gamma_value2(480,1));
    
    nuclei(1:number_of_nuclei,1)=rand(number_of_nuclei,1);%;   % set where are you in your cell cycle 
    nuclei(1:number_of_nuclei,4)=rand(number_of_nuclei,1);%;   % number that regulates type of prog with cell cycle length
    %[~, random_index]=datasample(nuclei,round(number_of_nuclei/2),'Replace',false); %choose random indexes for half of the cells
    %nuclei(random_index,4)=1000;%change them to other cell type
    
    
    
    nuclei(1:number_of_nuclei,2)=gamma_value(cell_cycle_length,number_of_nuclei);%./nuclei(1:number_of_nuclei,4)';%;   % cell cycle
    nuclei(1:number_of_nuclei,3)=1;%ceil(nuclei(1:number_of_nuclei,1)*4);%;   % type of cell
 
    number_of_cycling_nuclei=round(number_of_nuclei*cycling_cell_rate);%set initial number of G0 cells
    [~, random_index]=datasample(nuclei,number_of_nuclei-number_of_cycling_nuclei,'Replace',false); %choose a random index
    nuclei(random_index,3)=0;
    nuclei(random_index,1)=rand(number_of_nuclei-number_of_cycling_nuclei,1)/3; % quiescent cells locked somewhere in G1
    

for time=0:time_step:total_time
%time
    [pp_dd,cell_cycle_length,cycling_cell_rate]=dynamics_data(round(time)+1); % set dynamics at this time point
%cell_cycle_length

number_of_Q_cells=numel(find(nuclei(:,3)==0));
number_of_P_cells=numel(find(nuclei(:,3)==1));


% 
        if  number_of_Q_cells/(number_of_P_cells+number_of_Q_cells) < (1-cycling_cell_rate) % > (1-cycling_cell_rate)+0.01 %&& number_of_G0_cells/number_of_P_cells 

            number_of_new_Q_cells=round((1-cycling_cell_rate)*(number_of_P_cells+number_of_Q_cells)-number_of_Q_cells);
              %how many new G0 cells we need 
            random_index=datasample(find(nuclei(:,3)==1),number_of_new_Q_cells,'Replace',false); %choose a random index
            nuclei(random_index,3)=0;
            nuclei(random_index,1)=rand(number_of_new_Q_cells,1)/3; % quiescent cells locked somewhere in G1
    

        end
          
%         if number_of_Q_cells/(number_of_P_cells+number_of_Q_cells) > (1-cycling_cell_rate)%((0.1*rand(1))*(1-cycling_cell_rate))
%                number_of_new_g1_cells=round(number_of_Q_cells-((1-cycling_cell_rate)*(number_of_P_cells+number_of_Q_cells)));
%                random_index=datasample(find(nuclei(:,3)==0),number_of_new_g1_cells,'Replace',false);
%                nuclei(random_index,3)=1;
%                %nuclei(random_index,1)=0; % initialize age of new G1 cell
%         end
     
          

    
index_cycling_cells=find(nuclei(1:number_of_nuclei,3)==1);% & nuclei(1:number_of_nuclei,3)~=5);
nuclei(:,1)=nuclei(:,1)+(time_step./nuclei(:,2)); %increase age of all cells
%nuclei(index_cycling_cells,3)=ceil(nuclei(index_cycling_cells,2)*4);%;   % type of cell

%nuclei(index_cycling_cells,1)=nuclei(index_cycling_cells,1)+(time_step./nuclei(index_cycling_cells,2)); %increase age of all cells
%nuclei(1:number_of_nuclei,3)=ceil(nuclei(1:number_of_nuclei,1)*4);%;   % type of cell

    for nuclei_index=index_cycling_cells'  
        
     nuclei(nuclei_index,2)=gamma_value(cell_cycle_length,1);%/nuclei(nuclei_index,4);
%        
                        if nuclei(nuclei_index,1)>1; % cell transits from M to G1 after mitosis

                            number_of_nuclei=number_of_nuclei+1;
                            nuclei(nuclei_index,1)=0; % you are starting your cell cycle
                            nuclei(nuclei_index,4)=rand(1,1); % type of progenitors of teh cell that just started the cycle
                            nuclei(number_of_nuclei,4)=rand(1,1);   % type of progenitors of the newborn cell that just started the cycle
                            %nuclei(number_of_nuclei,4)=nuclei(nuclei_index,4); % set type of prog as mother cell
                            
                             mean(nuclei(:,4))
                            nuclei(number_of_nuclei,2)=gamma_value(cell_cycle_length,1);%/nuclei(number_of_nuclei,4);
%                             
                            %now establish mode of division
                            
                            nuclei(nuclei_index,3)=1+round(((-pp_dd/2)/1)+nuclei(nuclei_index,4));% old cell 
                                               
                            nuclei(number_of_nuclei,3)=1+round(((-pp_dd/2)/1)+nuclei(number_of_nuclei,4));% new cell
%                             nuclei(nuclei_index,3)=1+round(((-pp_dd/2)/nuclei(nuclei_index,4))+rand(1,1));% old cell 
%                             nuclei(number_of_nuclei,3)=1+round(((-pp_dd/2)/nuclei(nuclei_index,4))+rand(1,1));% new cell
                            
%                             if nuclei(number_of_nuclei,3)==2
%                                 nuclei(number_of_nuclei,3)=5;
%                             end
%                             if nuclei(nuclei_index,3)==2
%                                 nuclei(nuclei_index,3)=5;
%                             end

 
                         end
                         
                         
            
        %end
    end
    all_data{simulation_index,1+(time/time_step)}=nuclei;
    



    
end
            
            %pause(00000.1)
            end
%close all;

        
    

general_index=2;
savefile=num2str(general_index,'pollo_new_constant_cell_cycle_%d.mat');
save(savefile,'all_data','time_step','total_time','cycling_cell_rate');
data_analysis_brdu_12

end  

function [valor_gamma]=gamma_value(mean,n)

% %mean=5;
CV=(mean/100*10)/mean; % The coefficient of variation (CV) is defined as the ratio of the standard deviation (SD) to the mean.

A= (mean^2)/((CV*mean)^2); % Shape
B= ((CV*mean)^2)/ mean; % Scale

valor_gamma = gamrnd(A,B,1,n);%figure(1);hist(R,30) % Generates random numbers from the gamma distribution with shape parameters in A and scale parameters in B

 %uncooment thsi if you want a exp distribution
%  valor_gamma = abs(mean +  (mean/100*30) *randn(1));

 %uncooment this if you want no stochasticity
%valor_gamma = mean;%gamrnd(A,B);%figure(1);hist(R,30) % Generates random numbers from the gamma distribution with shape parameters in A and scale parameters in B

end


function [valor_gamma]=gamma_value2(mean,n)

% %mean=5;
CV=(mean/100*5)/mean; % The coefficient of variation (CV) is defined as the ratio of the standard deviation (SD) to the mean.

A= (mean^2)/((CV*mean)^2); % Shape
B= ((CV*mean)^2)/ mean; % Scale

valor_gamma = gamrnd(A,B,1,n);%figure(1);hist(R,30) % Generates random numbers from the gamma distribution with shape parameters in A and scale parameters in B

 %uncooment thsi if you want a exp distribution
%  valor_gamma = abs(mean +  (mean/100*30) *randn(1));

 %uncooment this if you want no stochasticity
%valor_gamma = mean;%gamrnd(A,B);%figure(1);hist(R,30) % Generates random numbers from the gamma distribution with shape parameters in A and scale parameters in B

end

function [pp_nn cell_cycle_length cycling_cell_rate]=dynamics_data(time)

% control
aa=[0.720911964	56.75435804	0.7
0.717673264	49.33882734	0.7
0.713885198	43.2155114	0.7
0.709477537	38.17586401	0.7
0.704380043	34.04815229	0.7
0.698526362	30.69159607	0.7
0.691859168	27.99155846	0.7
0.684336453	25.85562522	0.7
0.675938646	24.21044124	0.7
0.666675885	22.99919955	0.7
0.656594424	22.17970043	0.7
0.645780911	21.72291835	0.7
0.634363224	21.61203039	0.7
0.622506878	21.84187407	0.7
0.610406709	22.4188146	0.7
0.598274485	23.36101275	0.7
0.586324041	24.69909643	0.7
0.574756141	26.47725182	0.7
0.563745322	28.75476576	0.7
0.553430496	31.60806805	0.7
0.543910146	35.13334278	0.7
0.535242049	39.44979893	0.7
0.527446674	44.70371498	0.7
0.520513029	51.07339875	0.7
0.514405683	58.77523458	0.7
0.509071894	68.07102633	0.7
0.504448098	79.27688751	0.7
0.500465342	92.77398214	0.7
0.497053517	109.0214819	0.7
0.494144424	128.5721805	0.7
0.491673827	152.0912961	0.7
0.489582679	180.3791014	0.7];

%fgf
bb=[0.803161407	46.3539918	0.9
0.802701502	37.92615725	0.9
0.802132705	31.24699765	0.9
0.801432303	25.96082493	0.9
0.800574468	21.78613276	0.9
0.799530701	18.50060491	0.9
0.798270812	15.92927529	0.9
0.796764631	13.93520886	0.9
0.794984598	12.41220703	0.9
0.792909265	11.2791489	0.9
0.790527538	10.47566707	0.9
0.787843161	9.958928285	0.9
0.784878595	9.701348268	0.9
0.781677197	9.689120185	0.9
0.778302649	9.921479465	0.9
0.774835023	10.41066647	0.9
0.771363743	11.18258482	0.9
0.767978614	12.27818946	0.9
0.764760773	13.75567673	0.9
0.761775429	15.69359141	0.9
0.759067648	18.19501504	0.9
0.756661519	21.3930581	0.9
0.754562149	25.45794878	0.9
0.752759452	30.60609649	0.9
0.751232626	37.11161373	0.9
0.749954439	45.32091045	0.9
0.748894806	55.67113872	0.9
0.748023454	68.71347011	0.9
0.747311696	85.14244562	0.9
0.746733468	105.8329609	0.9
0.746265799	131.8868566	0.9
0.745888907	164.6915946	0.9];





%figure(2);subplot(1,2,1);plot(36:46,a(36:46,2));axis([36 46 0 15]);axis square
%figure(2);subplot(1,2,2);plot(36:46,a(36:46,1));axis([36 46 -1 1]);axis square
%pp_nn=a(time,1);
%cell_cycle_length=a(time,2);
%cycling_cell_rate=a(time,3);

% uncomment for control

%this is to make T variable
% T_1=linspace(30, 30, 12);
% T_2=linspace(30, 10, 10);
% T_3=linspace(10, 30, 10);
% aa(:,2)=[T_1,T_2,T_3]';

%this is to make pp-dd variable
% ppdd_1=linspace(1, 1, 12);
% ppdd_2=linspace(1, 0.5, 10);
% ppdd_3=linspace(0.5, 0, 10);
% aa(:,1)=[ppdd_1,ppdd_2,ppdd_3]';

%this is to make gamma variable
gamma_1=linspace(0.55, 0.55, 12);
gamma_2=linspace(0.55, 0.75, 20);

aa(:,3)=[gamma_1,gamma_2]';



% a=ones(49,3);
% a(:,:)=a(:,:).*aa(1,:);
% a(49-31:end,:)=aa;
% 
% 
% 
% pp_nn=a(time,1);cell_cycle_length=a(time,2);cycling_cell_rate=a(time,3);

% uncomment for control
a=ones(49,3);
a(:,:)=a(:,:).*aa(1,:);
a(49-31:end,:)=aa;
pp_nn=a(time,1);cell_cycle_length=a(time,2);cycling_cell_rate=a(time,3);


% % uncomment for +FGF
b=ones(49,3);
b(:,:)=b(:,:).*bb(1,:);
b(49-31:end,:)=bb;
%pp_nn=b(time,1);
%cell_cycle_length=b(time,2);
%cycling_cell_rate=b(time,3);

 %pp_nn=0;%a(time,1);
 %cell_cycle_length=20;%a(time,2);
 %cycling_cell_rate=1;%a(time,3);



end
