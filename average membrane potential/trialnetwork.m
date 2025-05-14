%This function mainly implements interregional connections in order to model the anterior and posterior hippocampus
function [memV_Sender_A_1,memV_Receiver_A_1,memV_Sender_P_1,memV_Receiver_P_1]=trialnetwork(step,totaltime,transient,mmm,Iext,N_bin,...
                neigh_send_internal_A,neigh_receiv_internal_A,neigh_receiv_external_A,neigh_send_external_A,neigh_A2_from_P5,neigh_A5_from_P5,neigh_A2_from_P2,neigh_A5_from_P2,...
                neigh_send_internal_P,neigh_receiv_internal_P,neigh_receiv_external_P,neigh_send_external_P,neigh_P2_from_A5,neigh_P5_from_A5,neigh_P2_from_A2,neigh_P5_from_A2)
rng(mmm,'twister');

%Sender = Alpha oscillations
%Receiver = Gamma oscillations

%The indexes from 0 to 399 stand for excitatory neurons and from 400-499 for inhibitory
% clc;clear 
% tic

N_node=500;

%% aHipp
gES_A=0.8;
gIS_A=16.4;

gER_A=3.0;
gIR_A=16.0;

gSR_A=2.0;
gRS_A=0.15;

gEPoisson_A=0.6;

%% pHipp
gES_P=0.8;
gIS_P=16.4;

gER_P=3.0;
gIR_P=16.0;

gSR_P=2.0;
gRS_P=0.15;

gEPoisson_P=0.6;

%% aHipp and pHipp
gA2_P5=0.6;%PS to AR
gA5_P5=1.0;%PS to AS
gA2_P2=1.0;%PR to AR
gA5_P2=0.3;%PR to AS

gP2_A5=0.6;%AS to PR
gP5_A5=1.0;%AS to PS
gP2_A2=1.0;%AR to PR
gP5_A2=0.3;%AR to PS

P=0.0;
n_ext = 10;


%% aHipp

auxrand1 = rand(1,500);
auxrand2 = rand(1,500);

% 
a_sender_A(1:400) = 0.02;
b_sender_A(1:400) = 0.2;
c_sender_A(1:400) = -65.0 + (15.0 * auxrand1(1:400) .* auxrand1(1:400));
d_sender_A(1:400) = 8.0 - (6.0 * auxrand1(1:400) .* auxrand1(1:400));
a_sender_A(401:N_node) = 0.02 + (0.08 * auxrand1(401:N_node));
b_sender_A(401:N_node) = 0.25 - (0.05 * auxrand1(401:N_node));
c_sender_A(401:N_node) = -65.0;
d_sender_A(401:N_node) = 2.0;

a_receiv_A(1:400) = 0.02;
b_receiv_A(1:400) = 0.2;
c_receiv_A(1:400) = -65.0 + (15.0 * auxrand2(1:400) .* auxrand2(1:400));
d_receiv_A(1:400) = 8.0 - (6.0 * auxrand2(1:400) .* auxrand2(1:400));
a_receiv_A(401:N_node) = 0.02 + (0.08 * auxrand2(401:N_node));
b_receiv_A(401:N_node) = 0.25 - (0.05 * auxrand2(401:N_node));
c_receiv_A(401:N_node) = -65.0;
d_receiv_A(401:N_node) = 2.0;

memV_Sender_A(1:totaltime)=0.0;
memV_Receiver_A(1:totaltime)=0.0;

%%
v_sender_A(1:N_node)=-60.0;
v_receiver_A(1:N_node)=-60.0;
u_sender_A(1:N_node)=0.2*(-60.0);
u_receiver_A(1:N_node)=0.2*(-60.0);
rsynE_sender_A(1:N_node)=0.001;
rsynE_receiver_A(1:N_node)=0.001;
rsynI_sender_A(1:N_node)=0.001;
rsynI_receiver_A(1:N_node)=0.001;
poisson_sender_A = 1.0 - exp(-1.6 * step);
poisson_receiver_A = 1.0 - exp(-2.0 * step);


%% pHipp
% 
auxrand3 = rand(1,500);
auxrand4 = rand(1,500);

a_sender_P(1:400) = 0.02;
b_sender_P(1:400) = 0.2;
c_sender_P(1:400) = -65.0 + (15.0 * auxrand3(1:400) .* auxrand3(1:400));
d_sender_P(1:400) = 8.0 - (6.0 * auxrand3(1:400) .* auxrand3(1:400));
a_sender_P(401:N_node) = 0.02 + (0.08 * auxrand3(401:N_node));
b_sender_P(401:N_node) = 0.25 - (0.05 * auxrand3(401:N_node));
c_sender_P(401:N_node) = -65.0;
d_sender_P(401:N_node) = 2.0;

a_receiv_P(1:400) = 0.02;
b_receiv_P(1:400) = 0.2;
c_receiv_P(1:400) = -65.0 + (15.0 * auxrand4(1:400) .* auxrand4(1:400));
d_receiv_P(1:400) = 8.0 - (6.0 * auxrand4(1:400) .* auxrand4(1:400));
a_receiv_P(401:N_node) = 0.02 + (0.08 * auxrand4(401:N_node));
b_receiv_P(401:N_node) = 0.25 - (0.05 * auxrand4(401:N_node));
c_receiv_P(401:N_node) = -65.0;
d_receiv_P(401:N_node) = 2.0;

memV_Sender_P(1:totaltime)=0.0;
memV_Receiver_P(1:totaltime)=0.0;

%%
v_sender_P(1:N_node)=-60.0;
v_receiver_P(1:N_node)=-60.0;
u_sender_P(1:N_node)=0.2*(-60.0);
u_receiver_P(1:N_node)=0.2*(-60.0);
rsynE_sender_P(1:N_node)=0.001;
rsynE_receiver_P(1:N_node)=0.001;
rsynI_sender_P(1:N_node)=0.001;
rsynI_receiver_P(1:N_node)=0.001;
poisson_sender_P = 1.0 - exp(-1.6 * step);
poisson_receiver_P = 1.0 - exp(-2.0 * step);

%% 
for t = 1:totaltime


    %% 
    spikes_S_A = v_sender_A > 30.0;
    spikes_R_A = v_receiver_A > 30.0;

    spikes_S_P = v_sender_P > 30.0;
    spikes_R_P = v_receiver_P > 30.0;

    %% aHipp

   %Sender
    IsynExcSender_A = zeros(1, 500); 
    IsynInhSender_A = zeros(1, 500); 
    IsynExcSender_ext_A = zeros(1, 500); 
    Isyn_A2_from_P5 = zeros(1, 500);
    Isyn_A5_from_P5 = zeros(1, 500);
    Isyn_A2_from_P2 = zeros(1, 500);
    Isyn_A5_from_P2 = zeros(1, 500);

    for jjj = 1:500
        nsynE_A = 0;
        nsynI_A = 0;
        nsynE_ext_A = 0;
        nsynE_A2_from_P5 = 0;
        nsynE_A5_from_P5 = 0;
        nsynE_A2_from_P2 = 0;
        nsynE_A5_from_P2 = 0;

        for kkk = 1:50
            if neigh_send_internal_A(jjj, kkk) <= 400
                nsynE_A = nsynE_A + spikes_S_A(neigh_send_internal_A(jjj, kkk));
            else
                nsynI_A = nsynI_A + spikes_S_A(neigh_send_internal_A(jjj, kkk));
            end
        end
        IsynExcSender_A(jjj) = gES_A * nsynE_A;
        IsynInhSender_A(jjj) = gIS_A * nsynI_A;
        for kkk = 1:20
            nsynE_ext_A = nsynE_ext_A + spikes_R_A(neigh_send_external_A(jjj, kkk));
        end
        IsynExcSender_ext_A(jjj) = gRS_A * nsynE_ext_A;
        
        for kkk = 1:n_ext
           
            nsynE_A2_from_P5 = nsynE_A2_from_P5 + spikes_S_P(neigh_A2_from_P5(jjj, kkk));
            nsynE_A5_from_P5 = nsynE_A5_from_P5 + spikes_S_P(neigh_A5_from_P5(jjj, kkk));
           
            nsynE_A2_from_P2 = nsynE_A2_from_P2 + spikes_R_P(neigh_A2_from_P2(jjj, kkk));
            nsynE_A5_from_P2 = nsynE_A5_from_P2 + spikes_R_P(neigh_A5_from_P2(jjj, kkk));
        end
        
        Isyn_A2_from_P5(jjj) = (1-P)*gA2_P5 * nsynE_A2_from_P5;
        Isyn_A5_from_P5(jjj) = (1-P)*gA5_P5 * nsynE_A5_from_P5;
        
        Isyn_A2_from_P2(jjj) = P*gA2_P2 * nsynE_A2_from_P2;
        Isyn_A5_from_P2(jjj) = P*gA5_P2 * nsynE_A5_from_P2;
    end
   
    %Receiver
    IsynExcReceiver_A = zeros(1, 500); 
    IsynInhReceiver_A = zeros(1, 500);
    IsynExcReceiver_ext_A = zeros(1, 500); 

    for jjj = 1:500
        nsynE_A = 0;
        nsynI_A = 0;
        nsynE_ext_A = 0;
        for kkk = 1:50
            if neigh_receiv_internal_A(jjj, kkk) <= 400
                nsynE_A = nsynE_A + spikes_R_A(neigh_receiv_internal_A(jjj, kkk));
            else
                nsynI_A = nsynI_A + spikes_R_A(neigh_receiv_internal_A(jjj, kkk));
            end
        end
        IsynExcReceiver_A(jjj) = gER_A * nsynE_A;
        IsynInhReceiver_A(jjj) = gIR_A * nsynI_A;
        for kkk = 1:20
            nsynE_ext_A = nsynE_ext_A + spikes_S_A(neigh_receiv_external_A(jjj, kkk));
        end
        IsynExcReceiver_ext_A(jjj) = gSR_A * nsynE_ext_A;
    end

    %% pHipp
   %Sender
    IsynExcSender_P = zeros(1, 500);
    IsynInhSender_P = zeros(1, 500); 
    IsynExcSender_ext_P = zeros(1, 500); 
    Isyn_P2_from_A5 = zeros(1, 500);
    Isyn_P5_from_A5 = zeros(1, 500);
    Isyn_P2_from_A2 = zeros(1, 500);
    Isyn_P5_from_A2 = zeros(1, 500);

    for jjj = 1:500
        nsynE_P = 0;
        nsynI_P = 0;
        nsynE_ext_P = 0;
        nsynE_P2_from_A5 = 0;
        nsynE_P5_from_A5 = 0;
        nsynE_P2_from_A2 = 0;
        nsynE_P5_from_A2 = 0;

        for kkk = 1:50
            if neigh_send_internal_P(jjj, kkk) <= 400
                nsynE_P = nsynE_P + spikes_S_P(neigh_send_internal_P(jjj, kkk));
            else
                nsynI_P = nsynI_P + spikes_S_P(neigh_send_internal_P(jjj, kkk));
            end
        end
        IsynExcSender_P(jjj) = gES_P * nsynE_P;
        IsynInhSender_P(jjj) = gIS_P * nsynI_P;
        for kkk = 1:20
            nsynE_ext_P = nsynE_ext_P + spikes_R_P(neigh_send_external_P(jjj, kkk));
        end
        IsynExcSender_ext_P(jjj) = gRS_P * nsynE_ext_P;
        
        for kkk = 1:n_ext
            
            nsynE_P2_from_A5 = nsynE_P2_from_A5 + spikes_S_A(neigh_P2_from_A5(jjj, kkk));
            nsynE_P5_from_A5 = nsynE_P5_from_A5 + spikes_S_A(neigh_P5_from_A5(jjj, kkk));
            
            nsynE_P2_from_A2 = nsynE_P2_from_A2 + spikes_R_A(neigh_P2_from_A2(jjj, kkk));
            nsynE_P5_from_A2 = nsynE_P5_from_A2 + spikes_R_A(neigh_P5_from_A2(jjj, kkk));
        end
        
        Isyn_P2_from_A5(jjj) = P*gP2_A5 * nsynE_P2_from_A5;
        Isyn_P5_from_A5(jjj) = P*gP5_A5 * nsynE_P5_from_A5;
        
        Isyn_P2_from_A2(jjj) = (1-P)*gP2_A2 * nsynE_P2_from_A2;
        Isyn_P5_from_A2(jjj) = (1-P)*gP5_A2 * nsynE_P5_from_A2;
    end

    %Receiver
    IsynExcReceiver_P = zeros(1, 500); 
    IsynInhReceiver_P = zeros(1, 500);
    IsynExcReceiver_ext_P = zeros(1, 500); 

    for jjj = 1:500
        nsynE_P = 0;
        nsynI_P = 0;
        nsynE_ext_P = 0;
        for kkk = 1:50
            if neigh_receiv_internal_P(jjj, kkk) <= 400
                nsynE_P = nsynE_P + spikes_R_P(neigh_receiv_internal_P(jjj, kkk));
            else
                nsynI_P = nsynI_P + spikes_R_P(neigh_receiv_internal_P(jjj, kkk));
            end
        end
        IsynExcReceiver_P(jjj) = gER_P * nsynE_P;
        IsynInhReceiver_P(jjj) = gIR_P * nsynI_P;
        for kkk = 1:20
            nsynE_ext_P = nsynE_ext_P + spikes_S_P(neigh_receiv_external_P(jjj, kkk));
        end
        IsynExcReceiver_ext_P(jjj) = gSR_P * nsynE_ext_P;
    end

    %% SENDER
    % aHipp Poisson Spike
    auxrand5 = rand(1, 500);
    auxrand7 = rand(1, 500);
    poissonSpike_A = poisson_sender_A > auxrand5;
    % pHipp Poisson Spike
    poissonSpike_P = poisson_sender_P > auxrand7;

  
    reset_indices_A = v_sender_A > 30.0;
    v_sender_A(reset_indices_A) = c_sender_A(reset_indices_A);
    u_sender_A(reset_indices_A) = u_sender_A(reset_indices_A) + d_sender_A(reset_indices_A);
    
    reset_indices_P = v_sender_P > 30.0;
    v_sender_P(reset_indices_P) = c_sender_P(reset_indices_P);
    u_sender_P(reset_indices_P) = u_sender_P(reset_indices_P) + d_sender_P(reset_indices_P);

   
    rsynE_sender_A = rsynE_sender_A + step * (1.0/5.26) * (-1.0 * rsynE_sender_A + IsynExcSender_A +  IsynExcSender_ext_A + Isyn_A5_from_P5 + Isyn_A5_from_P2 + gEPoisson_A * poissonSpike_A);
    rsynI_sender_A = rsynI_sender_A + step * (1.0/5.6) * (-1.0 * rsynI_sender_A +  IsynInhSender_A);
  
    rsynE_sender_P = rsynE_sender_P + step * (1.0/5.26) * (-1.0 * rsynE_sender_P + IsynExcSender_P +  IsynExcSender_ext_P + Isyn_P5_from_A5 + Isyn_P5_from_A2 + gEPoisson_P * poissonSpike_P);
    rsynI_sender_P = rsynI_sender_P + step * (1.0/5.6) * (-1.0 * rsynI_sender_P +  IsynInhSender_P);


   
    IDC_A = zeros(1, 500);
    IDC_A(1:400) = 0.0+Iext;
   
    IDC_P = zeros(1, 500);
    IDC_P(1:400) = 0.0+Iext;

   
    v_sender_A = v_sender_A + step * (0.04 * v_sender_A .* v_sender_A + 5.0 * v_sender_A + 140.0 - u_sender_A + rsynE_sender_A .* (0.0 - v_sender_A) + rsynI_sender_A .* (-65.0 - v_sender_A) + IDC_A);
    u_sender_A = u_sender_A + step * (a_sender_A .* (b_sender_A .* v_sender_A - u_sender_A));
    
    v_sender_P = v_sender_P + step * (0.04 * v_sender_P .* v_sender_P + 5.0 * v_sender_P + 140.0 - u_sender_P + rsynE_sender_P .* (0.0 - v_sender_P) + rsynI_sender_P .* (-65.0 - v_sender_P) + IDC_P);
    u_sender_P = u_sender_P + step * (a_sender_P .* (b_sender_P .* v_sender_P - u_sender_P));

    %% RECEIVER
    %  aHipp Poisson Spike
    auxrand6 = rand(1, 500);
    auxrand8 = rand(1, 500);
    poissonSpike_A = poisson_receiver_A > auxrand6;
    % pHipp Poisson Spike
    poissonSpike_P = poisson_receiver_P > auxrand8;

    
    reset_indices_A = v_receiver_A > 30.0;
    v_receiver_A(reset_indices_A) = c_receiv_A(reset_indices_A);
    u_receiver_A(reset_indices_A) = u_receiver_A(reset_indices_A) + d_receiv_A(reset_indices_A);
   
    reset_indices_P = v_receiver_P > 30.0;
    v_receiver_P(reset_indices_P) = c_receiv_P(reset_indices_P);
    u_receiver_P(reset_indices_P) = u_receiver_P(reset_indices_P) + d_receiv_P(reset_indices_P);

    
    rsynE_receiver_A = rsynE_receiver_A + step * (1.0/5.26) * (-1.0 * rsynE_receiver_A + IsynExcReceiver_A + IsynExcReceiver_ext_A + Isyn_A2_from_P5 + Isyn_A2_from_P2 + gEPoisson_A * poissonSpike_A);
    rsynI_receiver_A = rsynI_receiver_A + step * (1.0/5.6) * (-1.0 * rsynI_receiver_A + IsynInhReceiver_A);
    
    rsynE_receiver_P = rsynE_receiver_P + step * (1.0/5.26) * (-1.0 * rsynE_receiver_P + IsynExcReceiver_P + IsynExcReceiver_ext_P + Isyn_P2_from_A5 + Isyn_P2_from_A2 + gEPoisson_P * poissonSpike_P);
    rsynI_receiver_P = rsynI_receiver_P + step * (1.0/5.6) * (-1.0 * rsynI_receiver_P + IsynInhReceiver_P);

    
    IDC_A = zeros(1, 500);
    IDC_A(1:400) = 25.0+Iext;
    
    IDC_P = zeros(1, 500);
    IDC_P(1:400) = 25.0+Iext;

   
    v_receiver_A = v_receiver_A + step * (0.04 * v_receiver_A .* v_receiver_A + 5.0 * v_receiver_A + 140.0 - u_receiver_A + rsynE_receiver_A .* (0.0 - v_receiver_A) + rsynI_receiver_A .* (-65.0 - v_receiver_A) + IDC_A);
    u_receiver_A = u_receiver_A + step * (a_receiv_A .* (b_receiv_A .* v_receiver_A - u_receiver_A));
   
    v_receiver_P = v_receiver_P + step * (0.04 * v_receiver_P .* v_receiver_P + 5.0 * v_receiver_P + 140.0 - u_receiver_P + rsynE_receiver_P .* (0.0 - v_receiver_P) + rsynI_receiver_P .* (-65.0 - v_receiver_P) + IDC_P);
    u_receiver_P = u_receiver_P + step * (a_receiv_P .* (b_receiv_P .* v_receiver_P - u_receiver_P));

   
    memV_Sender_A(t) = mean(v_sender_A(1:400));
    memV_Receiver_A(t) = mean(v_receiver_A(1:400));
    
    memV_Sender_P(t) = mean(v_sender_P(1:400));
    memV_Receiver_P(t) = mean(v_receiver_P(1:400));
end
    memV_Sender_A_1 = memV_Sender_A(transient:N_bin:totaltime-1);
    memV_Receiver_A_1 = memV_Receiver_A(transient:N_bin:totaltime-1);
    memV_Sender_P_1 = memV_Sender_P(transient:N_bin:totaltime-1);
    memV_Receiver_P_1 = memV_Receiver_P(transient:N_bin:totaltime-1);