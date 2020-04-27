%% control constant T no quiescence no differenetation
close all;clear all;clc;
% good intro in this paper http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3524571/
%cd /Users/davidgmiguez/Dropbox/work/model_pollo/programs_matlab/Differentiation_Theory/_growth_rate_limit
% 07/07/2015

load('pollo_new_constant_cell_cycle_2.mat');
   
[number_of_simulations,total_time]=size(all_data);
time=(0:time_step:round(total_time*time_step)-1);

P_cells=zeros(number_of_simulations,total_time); % initialize some matrices to write
Q_cells=zeros(number_of_simulations,total_time); % initialize some matrices to write
D_cells=zeros(number_of_simulations,total_time); % initialize some matrices to write
newborncells=cell(number_of_simulations,total_time/time_step); % initialize some matrices to write

length_of_the_experiment=24;
time_for_brdu_injection=(total_time-length_of_the_experiment)/time_step;

for m=1:number_of_simulations
for n=1:total_time

[~,P]=find(all_data{m,n}(:,3)==1);
[number_of_P,~]=size(P);
P_cells(m,n)=number_of_P;

[~,Q]=find(all_data{m,n}(:,3)==0);
[number_of_Q,~]=size(Q);
Q_cells(m,n)=number_of_Q;

[~,D]=find(all_data{m,n}(:,3)==2);
[number_of_D,~]=size(D);
D_cells(m,n)=number_of_D;

aan=find(all_data{m,n}(:,1)==0); 
mothers=aan(1:numel(aan)/2);
daugthers=aan(1+numel(aan)/2:numel(aan));
newborncells{m,n}=[mothers';daugthers'];

end

D_cells(m,:)=D_cells(m,:)+gamma_value2(95,1); % add randomness to the initial value of D

end

load('pollo_new_constant_cell_cycle_2.mat');
figure(100);plot(time,D_cells,'r');title('P,D'); hold on;%,'color',[0.8 0.8 1]);hold on
figure(100);plot(time,P_cells+Q_cells,'g');
axis([0 50  0 1500])

figure(1);subplot(1,3,1);plot(time,D_cells,'r');title('P,D'); hold on;%,'color',[0.8 0.8 1]);hold on
figure(1);subplot(1,3,1);plot(time,P_cells+Q_cells,'g');
%figure(1);subplot(1,3,1);plot(time,P_cells+Q_cells+D_cells,'b');%title('D')%,'color',[0.8 0.8 1]);hold on
%figure(1);subplot(1,3,1);plot(time,P_cells+D_cells+G0_cells,'b');%title('D')%,'color',[0.8 0.8 1]);hold on


progenitor_mean=mean(P_cells);
differentiatted_mean=mean(D_cells);
quiescent_mean=mean(Q_cells);

progenitor_std=std(P_cells);
differentiatted_std=std(D_cells);
quiescent_std=std(Q_cells);

space_err_bar=2;
figure(1);subplot(1,3,2);errorbar(time(1:space_err_bar:end),progenitor_mean(1:space_err_bar:end),progenitor_std(1:space_err_bar:end),progenitor_std(1:space_err_bar:end),'g');hold on;%,'color',[0.8 0.8 1]);hold on
%plot(time,D_mean,'.-r');%title('D')%,'color',[0.8 0.8 1]);hold on
axis([time_for_brdu_injection total_time min(min(progenitor_mean+progenitor_std),min(differentiatted_mean+differentiatted_std)) max(max(progenitor_mean+progenitor_std),max(differentiatted_mean+differentiatted_std))]);

figure(1);subplot(1,3,2);errorbar(time(1:space_err_bar:end),differentiatted_mean(1:space_err_bar:end),differentiatted_std(1:space_err_bar:end),differentiatted_std(1:space_err_bar:end),'r');hold on;
axis square;

figure(1);subplot(1,3,2);errorbar(time(1:space_err_bar:end),quiescent_mean(1:space_err_bar:end),quiescent_std(1:space_err_bar:end),quiescent_std(1:space_err_bar:end),'b');hold on;
axis square
figure(1);subplot(1,3,3);plot(time,Q_cells./(P_cells+Q_cells));title('rate of G0')%,'color',[0.8 0.8 1]);hold on
axis square

figure(1);subplot(1,3,1);
axis([time_for_brdu_injection total_time min(min(progenitor_mean+progenitor_std),min(differentiatted_mean+differentiatted_std)) max(max(progenitor_mean+progenitor_std),max(differentiatted_mean+differentiatted_std))]);
axis square

%% NOW plot both types of porgenitors

P1_cells=zeros(number_of_simulations,total_time);
P2_cells=zeros(number_of_simulations,total_time);
%D_cells=zeros(number_of_simulations,total_time);

for m=1:number_of_simulations
for n=1:total_time

[~,P1]=find(all_data{m,n}(:,3)==1 & all_data{m,n}(:,4)==1);
[number_of_P1,~]=size(P1);
P1_cells(m,n)=number_of_P1;

[~,P2]=find(all_data{m,n}(:,3)==1 & all_data{m,n}(:,4)==1000);
[number_of_P2,~]=size(P2);
P2_cells(m,n)=number_of_P2;

end
end

figure(5);subplot(1,3,1);plot(time,P1_cells,'r'); hold on;%,'color',[0.8 0.8 1]);hold on
figure(5);subplot(1,3,1);plot(time,P2_cells,'g');title('P1,P2');legend('P1,P2');
axis([time_for_brdu_injection total_time 0 750]);axis square
progenitor1_mean=mean(P1_cells);
progenitor2_mean=mean(P2_cells);
progenitor1_std=std(P1_cells);
progenitor2_std=std(P2_cells);

% space_err_bar=2;
% subplot(1,3,2);errorbar(time(1:space_err_bar:end),progenitor1_mean(1:space_err_bar:end),progenitor1_std(1:space_err_bar:end),progenitor1_std(1:space_err_bar:end),'r');hold on;
% subplot(1,3,2);errorbar(time(1:space_err_bar:end),progenitor2_mean(1:space_err_bar:end),progenitor2_std(1:space_err_bar:end),progenitor2_std(1:space_err_bar:end),'g');hold on;

%figure(5);subplot(1,3,1);plot(time,Q_cells,'b');%title('D')%,'color',[0.8 0.8 1]);hold on
%figure(1);subplot(1,3,1);plot(time,P_cells+D_cells+G0_cells,'b');%title('D')%,'color',[0.8 0.8 1]);hold on
axis([time_for_brdu_injection total_time min(min(progenitor_mean+progenitor_std),min(differentiatted_mean+differentiatted_std)) max(max(progenitor_mean+progenitor_std),max(differentiatted_mean+differentiatted_std))]);
axis square
%total_progenitors=progenitor1_mean+progenitor2_mean;
cell_cycle_double_population=((progenitor1_mean*0.75)+(progenitor2_mean*0.0075))./(progenitor1_mean+progenitor2_mean);
subplot(1,3,3);plot(time,cell_cycle_double_population);

axis([time_for_brdu_injection total_time -1 1]);
axis square
%%

%lets add a arrow in all_data that marks cells in s phase
% lets assume that G1, S and G2M are 1/3 of the cell cycle each
start_S_phase=1/3;
finish_S_phase=2/3;

%brdu_positive_cells_fixed_injection=zeros(number_of_simulations,total_time);

%length_of_the_experiment=11;
time_for_brdu_injection=(total_time-length_of_the_experiment)/time_step;

brdu_positive_cells_fixed_injection=cell(number_of_simulations,total_time/time_step);
rate_of_brdu_positive_cells_fixed_injection=zeros(number_of_simulations,total_time);

brdu_positive_cells_fixed_fixation=cell(number_of_simulations,total_time/time_step);
rate_of_brdu_positive_cells_fixed_fixation=zeros(number_of_simulations,total_time);
rate_of_edu_positive_cells_fixed_fixation=zeros(number_of_simulations,total_time);


%% Brdu at a given time and fixing at increasing times


for m=1:number_of_simulations
            [a,b]=size(all_data{m,total_time});
            brdu_positive_cells_fixed_injection{m,time_for_brdu_injection}(1,:)=zeros(1,a);
            indexes_of_cells_in_phase_S=find(all_data{m,time_for_brdu_injection}(:,1)>start_S_phase & all_data{m,time_for_brdu_injection}(:,1)<finish_S_phase);
            brdu_positive_cells_fixed_injection{m,time_for_brdu_injection}(1,indexes_of_cells_in_phase_S)=1; %old
            
        for n=time_for_brdu_injection+1:1:total_time
            brdu_positive_cells_fixed_injection{m,n}(1,:)=zeros(1,a);
            brdu_positive_cells_fixed_injection{m,n}(1,:)=brdu_positive_cells_fixed_injection{m,n-1}(1,:);%cells that where positive before %old
            cells_in_S_phase_at_this_time_point=find(all_data{m,n}(:,1)>start_S_phase & all_data{m,n}(:,1)<finish_S_phase);
            brdu_positive_cells_fixed_injection{m,n}(1,cells_in_S_phase_at_this_time_point)=1;
                brdu_state_of_mother_cells=brdu_positive_cells_fixed_injection{m,n}(1,newborncells{m,n}(1,:)); % This is to see the Brdu of mother cells.
           for i=1:numel(brdu_state_of_mother_cells)
                if brdu_state_of_mother_cells(i)==1
                   brdu_positive_cells_fixed_injection{m,n}(1,newborncells{m,n}(2,i))=1;
                end
            end
            
             ddd=find(all_data{m,n}(:,3)~=2);
             [number_of_progenitors,~]=size(ddd);
             rate_of_brdu_positive_cells_fixed_injection(m,n)=sum(brdu_positive_cells_fixed_injection{m,n}(1,ddd))/number_of_progenitors; 

        end
        %figure(2);subplot(4,2,1);plot(time(1+time_for_brdu_injection:total_time),rate_of_brdu_positive_cells_fixed_injection(:,1+time_for_brdu_injection:total_time)); hold on;%,'color',[0.8 0.8 1]);hold on
end




figure(2);subplot(4,2,1);plot(time(1+time_for_brdu_injection:total_time),rate_of_brdu_positive_cells_fixed_injection(:,1+time_for_brdu_injection:total_time)); hold on;%,'color',[0.8 0.8 1]);hold on
axis([time(1+time_for_brdu_injection) time(total_time) 0 1.2]);

% mean_rate_of_brdu_positive_cells=zeros(1,46);
% std_rate_of_brdu_positive_cells=zeros(1,46);
% 
% for j=36:46
%     %((j-35)*3)-2:((j-35)*3)
%     mean_rate_of_brdu_positive_cells(j)=mean(rate_of_brdu_positive_cells_fixed_injection(((j-35)*3)-2:((j-35)*3),j));
%     std_rate_of_brdu_positive_cells(j)=std(rate_of_brdu_positive_cells_fixed_injection(((j-35)*3)-2:((j-35)*3),j));
% 
% end

mean_rate_of_brdu_positive_cells=mean(rate_of_brdu_positive_cells_fixed_injection,1);
std_rate_of_brdu_positive_cells=std(rate_of_brdu_positive_cells_fixed_injection,1);
axis square
point=8;
% subplot(4,2,2);errorbar(time(1+time_for_brdu_injection:total_time),mean_rate_of_brdu_positive_cells(1+time_for_brdu_injection:total_time),std_rate_of_brdu_positive_cells(1+time_for_brdu_injection:total_time),std_rate_of_brdu_positive_cells(1+time_for_brdu_injection:total_time),'r');hold on;
%axis([time(1+time_for_brdu_injection) time(total_time) 0 1.2]);
%axis square
P = polyfit(time(1+time_for_brdu_injection:time_for_brdu_injection+point),mean_rate_of_brdu_positive_cells(1+time_for_brdu_injection:time_for_brdu_injection+point),1);
Y = polyval(P,time(1+time_for_brdu_injection:time_for_brdu_injection+point));
hold on;plot(time(1+time_for_brdu_injection:time_for_brdu_injection+point),Y,'-k');

%gamma=mean(mean_rate_of_brdu_positive_cells((point+10)/time_step:i/time_step))
gamma=mean(mean_rate_of_brdu_positive_cells(end-3:end));

disp('BrdU fixed_injection');
Tc=gamma/P(1)
Ts=P(2)/P(1)+(time_for_brdu_injection)
gamma
text(30,0.5,num2str(Ts));
text(30,0.75,num2str(Tc));

%% now cumulative curve changing injection point

for m=1:number_of_simulations
     
        [a,b]=size(all_data{m,total_time});
        brdu_positive_cells_fixed_fixation{m,1}(1:6,:)=zeros(6,a);
        brdu_positive_cells_fixed_fixation{m,1}(2,:)=1000;
        brdu_positive_cells_fixed_fixation{m,1}(6,:)=1000;
        
        indexes_of_cells_in_phase_S=find(all_data{m,total_time}(:,1)>start_S_phase & all_data{m,total_time}(:,1)<finish_S_phase);
        brdu_positive_cells_fixed_fixation{m,1}(1,indexes_of_cells_in_phase_S)=1;
       
        % when do you become brdu positive
        brdu_positive_cells_fixed_fixation{m,1}(2,indexes_of_cells_in_phase_S)=total_time;
        
        indexes_of_mother_cells=newborncells{m,total_time}(1,:);
        indexes_of_newborn_daughter_cells=newborncells{m,total_time}(2,:);
        
        % when do you have been born 
        brdu_positive_cells_fixed_fixation{m,1}(3,[indexes_of_mother_cells,indexes_of_newborn_daughter_cells])=total_time; 
        % who's your mother
        brdu_positive_cells_fixed_fixation{m,1}(4,indexes_of_newborn_daughter_cells)=indexes_of_mother_cells; 
      
        ddd=find(all_data{m,total_time}(:,3)~=2);
        [number_of_progenitors,~]=size(ddd);
        
for n=total_time-1:-1:time_for_brdu_injection

       brdu_positive_cells_fixed_fixation{m,total_time-n+1}(:,:)=brdu_positive_cells_fixed_fixation{m,total_time-n}(:,:); %cells that were positive before
       indexes_of_cells_in_phase_S=find(all_data{m,n}(:,1)>start_S_phase & all_data{m,n}(:,1)<finish_S_phase);
       brdu_positive_cells_fixed_fixation{m,total_time-n+1}(1,indexes_of_cells_in_phase_S)=1;
       brdu_positive_cells_fixed_fixation{m,total_time-n+1}(2,indexes_of_cells_in_phase_S)=n; % when do you became brdu positive
          
        indexes_of_mother_cells=newborncells{m,n}(1,:);
        indexes_of_newborn_daughter_cells=newborncells{m,n}(2,:);
        % label new cells that just have been born
        brdu_positive_cells_fixed_fixation{m,total_time-n+1}(3,indexes_of_newborn_daughter_cells)=n; 
        % who's your mother
        brdu_positive_cells_fixed_fixation{m,total_time-n+1}(4,indexes_of_newborn_daughter_cells)=indexes_of_mother_cells;
      
       
end
     
   
   %rate_of_brdu_positive_cells_fixed_fixation(m,total_time-n+1)=sum(brdu_positive_cells_fixed_fixation{m,total_time-n+1}(1,ddd))/number_of_progenitors;

   
end 

for m=1:number_of_simulations
    ddd=find(all_data{m,total_time}(:,3)~=2); % index of all progenitors at final time
    [number_of_progenitors,~]=size(ddd);
    [a,b]=size(all_data{m,total_time});
    for n=total_time-1:-1:time_for_brdu_injection
        
       %if time of BRDU positive of mother is set before time of birth, change status of daugther to positive 
   for index=1:a
       % who's my mother
        my_mother=brdu_positive_cells_fixed_fixation{m,total_time-n+1}(4,index);
        % how is the brdu status of my mother 
        if my_mother ~= 0
        brdu_positive_cells_fixed_fixation{m,total_time-n+1}(5,index)=brdu_positive_cells_fixed_fixation{m,total_time-n+1}(1,my_mother);
        
        % when did your mother got brdu+
        brdu_positive_cells_fixed_fixation{m,total_time-n+1}(6,index)=brdu_positive_cells_fixed_fixation{m,total_time-n+1}(2,my_mother);
       % mmm=brdu_positive_cells_fixed_fixation{m,total_time-n+1}(2,indexes_of_mother_cells)
       
       
        time_of_brdu_status_of_mother=brdu_positive_cells_fixed_fixation{m,total_time-n+1}(6,index);
        time_of_birth=brdu_positive_cells_fixed_fixation{m,total_time-n+1}(3,index);

      % if time_of_brdu_status_of_mother<1000 
       %    time_of_brdu_status_of_mother
       %    time_of_birth
           
            if time_of_brdu_status_of_mother<=time_of_birth
                brdu_positive_cells_fixed_fixation{m,total_time-n+1}(1,index)=1;%brdu_positive_cells{m,total_time-n+time_for_brdu_injection}(2,brdu_positive_cells{m,total_time-n+time_for_brdu_injection}(6,index));
            end
     %  end
        end
   end
         rate_of_brdu_positive_cells_fixed_fixation(m,total_time-n+1)=sum(brdu_positive_cells_fixed_fixation{m,total_time-n+1}(1,ddd))/number_of_progenitors; 
         number_of_edu_positive_cells_fixed_fixation(m,total_time-n+1)=sum(brdu_positive_cells_fixed_fixation{m,total_time-n+1}(1,ddd));
   
    end
end


%figure(2);subplot(4,2,3);plot(time(time_for_brdu_injection+1:total_time),rate_of_brdu_positive_cells_fixed_fixation(:,2:time_for_brdu_injection-3)); hold on;%,'color',[0.8 0.8 1]);hold on

% mean_rate_of_brdu_positive_cells=zeros(1,46);
% std_rate_of_brdu_positive_cells=zeros(1,46);

mean_rate_of_brdu_positive_cells=mean(rate_of_brdu_positive_cells_fixed_fixation(:,2:time_for_brdu_injection-3),1);
std_rate_of_brdu_positive_cells=std(rate_of_brdu_positive_cells_fixed_fixation(:,2:time_for_brdu_injection-3),1);



% for j=2:12
%     %((j-35)*3)-2:((j-35)*3)
%     mean_rate_of_brdu_positive_cells(j-1)=mean(rate_of_brdu_positive_cells_fixed_injection(((j-1)*3)-2:((j-1)*3),j));
%     std_rate_of_brdu_positive_cells(j-1)=std(rate_of_brdu_positive_cells_fixed_injection(((j-1)*3)-2:((j-1)*3),j));
% 
% end


%figure(i);subplot(2,3,2);errorbar(time(1:round(exp_point/10):exp_point),mean_rate_of_brdu_positive_cells(1:round(exp_point/10):exp_point),std_rate_of_brdu_positive_cells(1:round(exp_point/10):exp_point),std_rate_of_brdu_positive_cells(1:round(exp_point/10):exp_point),'r');hold on;
%,'color',[0.8 0.8 1]);hold on
%point=6;
%subplot(4,2,4);errorbar(time(1:1:2*point),mean_rate_of_brdu_positive_cells(1:1:2*point),std_rate_of_brdu_positive_cells(1:1:2*point),std_rate_of_brdu_positive_cells(1:1:2*point),'r');hold on;
%point=6;
%P = polyfit(time(1:point/time_step),mean_rate_of_brdu_positive_cells(1:point/time_step),1);
%Y = polyval(P,time(1:point/time_step));
%plot(time(1:point/time_step),Y);

%gamma=mean(mean_rate_of_brdu_positive_cells(end-3:end))

%disp('BrdU fixed_fixation');

%Tc=gamma/P(1)
%Ts=P(2)/P(1)

%text(10,0.5,num2str(Ts));
%text(10,0.75,num2str(Tc));


%% dual cumulative curve

for m=1:number_of_simulations

 ddd=find(all_data{m,total_time}(:,3)~=2); % index of all progenitors at final time
 [number_of_progenitors,~]=size(ddd);
 number_of_final_brdu_positive_cells=sum(brdu_positive_cells_fixed_fixation{m,length_of_the_experiment+1}(1,ddd));
  
    for n=time_for_brdu_injection:total_time
        rate_of_edu_positive_cells_fixed_fixation(m,total_time-n+1)=number_of_edu_positive_cells_fixed_fixation(m,total_time-n+1)/number_of_final_brdu_positive_cells;
        
        
    end

end



%subplot(4,2,3);plot(time(time_for_brdu_injection+1:total_time),rate_of_edu_positive_cells_fixed_fixation(:,2:time_for_brdu_injection-3)); hold on;%,'color',[0.8 0.8 1]);hold on
subplot(4,2,3);plot(time(time_for_brdu_injection+1:total_time),rate_of_edu_positive_cells_fixed_fixation(:,2:total_time-time_for_brdu_injection+1)); hold on;%,'color',[0.8 0.8 1]);hold on

axis([time(1+time_for_brdu_injection) time(total_time) 0 1.2]);
axis square


mean_rate_of_edu_positive_cells=mean(rate_of_edu_positive_cells_fixed_fixation(:,2:time_for_brdu_injection-3),1);
std_rate_of_edu_positive_cells=std(rate_of_edu_positive_cells_fixed_fixation(:,2:time_for_brdu_injection-3),1);

% mean_rate_of_edu_positive_cells=zeros(1,31);
% std_rate_of_edu_positive_cells=zeros(1,31);
% 
% for j=2:12
%     %((j-35)*3)-2:((j-35)*3)
%     mean_rate_of_edu_positive_cells(j-1)=mean(rate_of_edu_positive_cells_fixed_fixation(((j-1)*3)-2:((j-1)*3),j));
%     std_rate_of_edu_positive_cells(j-1)=std(rate_of_edu_positive_cells_fixed_fixation(((j-1)*3)-2:((j-1)*3),j));
% 
% end


% subplot(4,2,4);errorbar(time(1:1:time_for_brdu_injection-4),mean_rate_of_edu_positive_cells,std_rate_of_edu_positive_cells,std_rate_of_edu_positive_cells,'r');hold on;
% axis([time(1) time(total_time-time_for_brdu_injection) 0 1.2]);
% axis square
point=10;
P = polyfit(time(1:point/time_step),mean_rate_of_edu_positive_cells(1:point/time_step),1);
Y = polyval(P,time(1:point/time_step));
plot(time(1:point/time_step),Y);

disp('dual cumulative BrdU');
Tc=1/P(1)
Ts=P(2)/P(1)

text(10,0.5,num2str(Ts));
text(10,0.75,num2str(Tc));



%% pulse and chase


for m=1:number_of_simulations
   for n=time_for_brdu_injection:total_time
       
      indexes_of_brdu_positive_cells=find(all_data{m,time_for_brdu_injection}(:,1)>start_S_phase & all_data{m,time_for_brdu_injection}(:,1)<finish_S_phase);% cells brdu positive at start of exp
      pulse_chase_cells{m,n}(1,indexes_of_brdu_positive_cells)=1;
      indexes_of_cells_in_phase_S=find(all_data{m,n}(:,1)>start_S_phase & all_data{m,n}(:,1)<finish_S_phase); % cells edu positive at this time point
      pulse_chase_cells{m,n}(2,indexes_of_cells_in_phase_S)=1;
      pulse_chase_cells{m,n}(3,:)=pulse_chase_cells{m,n}(1,:)+pulse_chase_cells{m,n}(2,:);
      
      index_number_of_double_positive_cells_pulse_chase=find(pulse_chase_cells{m,n}(3,:)==2);
      
      number_of_edu_positive_cells_pulse_chase=sum(pulse_chase_cells{m,n}(3,index_number_of_double_positive_cells_pulse_chase))/2;
      number_of_brdu_positive_cells_pulse_chase=sum(pulse_chase_cells{m,n}(1,:));
      
      %rate_of_brdu_positive_cells_pulse_chase(m,n)=sum(brdu_positive_cells{m,n}(2,ddd))/progenitors;
      rate_of_edu_positive_cells_pulse_chase(m,n)=number_of_edu_positive_cells_pulse_chase/number_of_brdu_positive_cells_pulse_chase;
      
end
end
figure(2);subplot(4,2,5);plot(time(time_for_brdu_injection:total_time),rate_of_edu_positive_cells_pulse_chase(:,time_for_brdu_injection:total_time)); axis([time_for_brdu_injection total_time 0 1]);hold on;%,'color',[0.8 0.8 1]);hold on
axis square

mean_rate_of_edu_positive_cells_pulse_chase=mean(rate_of_edu_positive_cells_pulse_chase,1);
std_rate_of_edu_positive_cells_pulse_chase=std(rate_of_edu_positive_cells_pulse_chase,1);
subplot(4,2,6);

% errorbar(time(time_for_brdu_injection:total_time),mean_rate_of_edu_positive_cells_pulse_chase(time_for_brdu_injection:total_time),std_rate_of_edu_positive_cells_pulse_chase(time_for_brdu_injection:total_time),std_rate_of_edu_positive_cells_pulse_chase(time_for_brdu_injection:total_time),'r');hold on;
% axis([time_for_brdu_injection total_time 0 1]);
% axis square

%gamma=mean(mean_exp_values(8:10))
disp('dual pulse chase method');

minimum_time_to_find_maxima=total_time-10;

%time_where_we_have_to_find_the_max=(total_time-15:total_time);

%[a,b]=max(mean_rate_of_edu_positive_cells_pulse_chase(minimum_time_to_find_maxima:total_time));

%maximo=total_time-5+b;
%vector_time=(maximo-3:0.1:maximo+5);
%vector_time=(total_time-5:0.1:total_time);
% we need to fit a cubic around the maximum point of the curve
P = polyfit((minimum_time_to_find_maxima:total_time),mean_rate_of_edu_positive_cells_pulse_chase(minimum_time_to_find_maxima:total_time),3);
Y = polyval(P,(minimum_time_to_find_maxima:0.1:total_time));
[a,b]=max(Y);
Tc=(-1+minimum_time_to_find_maxima+b/10)-length_of_the_experiment%Tc=vector_time(b)-time_for_brdu_injection
plot((minimum_time_to_find_maxima:0.1:total_time)-1,Y,'k');

axis square
P = polyfit(time(time_for_brdu_injection:time_for_brdu_injection+3),mean_rate_of_edu_positive_cells_pulse_chase(time_for_brdu_injection:time_for_brdu_injection+3),1);
Y = polyval(P,time(time_for_brdu_injection:time_for_brdu_injection+3));
Ts=-1/P(1)
plot(time(time_for_brdu_injection:time_for_brdu_injection+3),Y);
axis square
%Tc=-1/P(1)
text(time_for_brdu_injection+point,0.5,num2str(Ts));
text(time_for_brdu_injection+point,0.75,num2str(Tc));

%% Now using my equations
a=[1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
0.9
0.8
0.7
0.6
0.5
0.4
0.3
0.2
0.1
0];

%gamma=0.75;

% figure(3);
% subplot(2,2,1);plot(time,D_cells);title('P,D'); hold on;%,'color',[0.8 0.8 1]);hold on
% subplot(2,2,1);plot(time,P_cells+Q_cells);%title('P')%,'color',[0.8 0.8 1]);hold on
% 
% subplot(2,2,2);%plot(time,P_mean,'.-g');
% title('mean P,mean D'); 
% hold on;
% errorbar(time(1:10:end),P_mean(1:10:end),P_std(1:10:end),P_std(1:10:end),'g');hold on;%,'color',[0.8 0.8 1]);hold on
% plot(time,D_mean,'.-r');%title('D')%,'color',[0.8 0.8 1]);hold on
% 
% errorbar(time(1:10:end),D_mean(1:10:end),D_std(1:10:end),D_std(1:10:end),'r');hold on;

% figure(3)
% subplot(2,1,1);
% P = polyfit(time,P_mean,4);
% fitted_P_mean = polyval(P,time);    
% plot(time,fitted_P_mean,'g');
% 
% 
% hold on;
% P = polyfit(time,D_mean,4);
% fitted_D_mean = polyval(P,time);
% plot(time,fitted_D_mean,'r');
% 
% hold on
% figure(3);subplot(2,1,2);plot(time,P_cells+Q_cells,'g');
% 
% hold on
% figure(3);subplot(2,1,2);plot(time,D_cells,'r');

% P_mean= mean(P_cells+Q_cells);
% P_std= std(P_cells+Q_cells);
% D_mean= mean(D_cells);
% D_std= std(D_cells);
% Q_mean= mean(Q_cells);
% Q_std= std(Q_cells);


P_mean= mean(P_cells+Q_cells);
P_std= std(P_cells+Q_cells);
D_mean= mean(D_cells);
D_std= std(D_cells);
Q_mean= mean(Q_cells);
Q_std= std(Q_cells);

%for m=1:number_of_simulations
for n=24:numel(P_cells(1,:)+Q_cells(1,:))
gamma=1;%a(n);
%pp_dd(n)=(fitted_P_mean(n)-fitted_P_mean(n-1))/(fitted_P_mean(n)-fitted_P_mean(n-1)+fitted_D_mean(n)-fitted_D_mean(n-1));

%T(n)=(time(n)-time(n-1))*log(1+(gamma*(pp_dd(n))))/log(fitted_P_mean(n)/fitted_P_mean(n-1));
pp_dd(n)=(P_mean(n)-P_mean(n-1))/(P_mean(n)-P_mean(n-1)+D_mean(n)-D_mean(n-1));


if pp_dd(n)>0
    T(n)=(time(n)-time(n-1))*log(1+(gamma*(abs(pp_dd(n)))))/abs(log(P_mean(n)/P_mean(n-1)));
else
   T(n)=(time(n)-time(n-1))*log(1+(gamma*(abs(pp_dd(n)))))/(log(P_mean(n)/P_mean(n-1)))/(0.9*abs(pp_dd(n))-1);
end

end
%end


subplot(4,2,7);
change_P=gradient(P_cells+Q_cells);
change_D=gradient(D_cells);

plot(time,change_P);title('P,D'); hold on;%,'color',[0.8 0.8 1]);hold on
plot(time,change_D);

axis([time_for_brdu_injection total_time min(min(min(change_D)),min(min(change_P))) max(max(max(change_D)),max(max(change_P)))]);hold on;%,'color',[0.8 0.8 1]);hold on
axis square

disp('branching method');
nanmean(T(time_for_brdu_injection+1:end))
nanmean(pp_dd(time_for_brdu_injection+1:end))

figure(10)
subplot(1,2,1);hold on;plot(time(time_for_brdu_injection+1:end),pp_dd(time_for_brdu_injection+1:end));title('pp-dd');

axis([time_for_brdu_injection total_time -1 1]);
axis square
subplot(1,2,2);hold on;plot(time(time_for_brdu_injection+1:end),T(time_for_brdu_injection+1:end));title('T');

axis([time_for_brdu_injection total_time 0 max(T)]);

axis square

general_index=1;
savefile=num2str(general_index,'simulations_fgf_%d.mat');
save(savefile,'D_cells','P_cells','Q_cells','cycling_cell_rate');


function [valor_gamma]=gamma_value2(mean,n)

% %mean=5;
CV=(mean/100*1)/mean; % The coefficient of variation (CV) is defined as the ratio of the standard deviation (SD) to the mean.

A= (mean^2)/((CV*mean)^2); % Shape
B= ((CV*mean)^2)/ mean; % Scale

valor_gamma = gamrnd(A,B,1,n);%figure(1);hist(R,30) % Generates random numbers from the gamma distribution with shape parameters in A and scale parameters in B

 %uncooment thsi if you want a exp distribution
%  valor_gamma = abs(mean +  (mean/100*30) *randn(1));

 %uncooment this if you want no stochasticity
%valor_gamma = mean;%gamrnd(A,B);%figure(1);hist(R,30) % Generates random numbers from the gamma distribution with shape parameters in A and scale parameters in B

end


