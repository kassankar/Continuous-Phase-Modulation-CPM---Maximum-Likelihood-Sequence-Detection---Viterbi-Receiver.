clc
clear
close all
%% Pulse shape & Variable initialization
MAIN             = MainFunctions;
pulse            = 2;      % 1 -> lorentzian pulse
                           % 2 -> GMSK pulse BT = 0.3
                           % 3 -> LRC pulse
                           % 4 -> LREC pulse
pulse_length     = 3;      % 1  -> Full response
                           % >1 -> Partial response
os               = 2^2;    % Over sampling frequency
Ts               = 1/os;   % Sampling Time
M_ary            = 2^1;    % M_ary symbols used 2 -> Binary
modulation_index = 0.5;    % Modulation index
width            = 0.7;    % This variable is used for Lorentzian Pulse only. (Not be used for pulse > 1)
dmin             = 1.78;   % GMSK BT=0.3
% 'decision_delay' Decide after how many observation symbols
% the Viterbi receiver can make a decision on the detected bit.
decision_delay   = 50;
%-------------- Modulated data ----------------------
snr           = 0:10;
BER_all       = zeros(1,length(snr));
%% frequency pulse
[g_t,q_t]     = MAIN.CREATECPMPULSE(pulse,pulse_length,width,os,0); % Function return the CPM pulse and phase.
                                                         % g_t = g(t) is the CPM pulse shape.
                                                         % q_t -> is the phase, integral of g_t.                                                      
m  = 1:M_ary;
if(pulse~=1)
    A_m    = 2*m-1-M_ary;       % Different modulation level
else
    modulation_index      = 2*modulation_index;
    A_m    =  0:M_ary-1;        % Different modulation level
end
%-------------------------------------------------------------------------------------------------------------

%% ---------------------------------------------- Receiver------------------------------------------------------
%---------------------------- states -- Branches -----------------------------------
[States_Number,Branche_Number,phase_states] = MAIN.states_Branches(modulation_index,pulse_length,M_ary);
[Branches,states_sort] = MAIN.states_transition(States_Number,Branche_Number,phase_states,pulse_length,M_ary,A_m,modulation_index);
%------------------------------- Branche_metric ----------------------------------
if(pulse_length>1)
    Branche_metrics = zeros(Branche_Number,length(((pulse_length-1)*os:pulse_length*os)));
else
    Branche_metrics = zeros(Branche_Number,length(1:pulse_length*os+1));
end
for i = 1:Branche_Number
    alpha                = Branches(i,2:end);
    Branche_metrics(i,:) = MAIN.compute_Branche_metric(os, modulation_index, pulse_length, alpha, g_t);
end

    
%-------------theoretical BER calculation-----------------------
EbNodB_vect_lin = 10.^(snr/10);
x               = sqrt(dmin.*EbNodB_vect_lin);
Q               = 0.5.*erfc(x/sqrt(2));
Pe              = Q;
%---------------------------------------------------------------
frame_max_length = 300./Pe*log2(M_ary);
BER_all    = zeros(1,length(snr));
tic
NBits_limit = [10^4,10^5]; % min, max of number of bits in the whole frame
% Itterate over all the snr
for snr_no = 1:length(snr)
    NBits      = frame_max_length(snr_no); % Number of bits (number of bits the whole frame)
    if(NBits>NBits_limit(2))
        NBits = NBits_limit(2);
    end
    if(NBits<NBits_limit(1))
        NBits = NBits_limit(1);
    end
    pack       = 10^4;                    
    Nbits_div  = round(NBits/pack); % Number of packs in a frame (each = to 10^4)
    % For each snr iterations, iterate over all Nbits_div, where
    % Nbits_div is the number of packs.
    for I = 1:Nbits_div
        Nbits       = log2(M_ary)*floor(round(NBits/Nbits_div)/log2(M_ary));
        rand('seed', 16);
        bits   = randi([0 1],Nbits,1).';
        bits_m = MAIN.Mapp(bits, M_ary, pulse);
        Nbits  = length(bits_m);               % Number of bits change after applying the mapping (only for M-ary>2).
        %% --------------------------- Modualtion -----------------------------------
        %----------------- Modulated signal (information dependent)------------------
        bits_s        = upsample(bits_m,os);
        t_seq         = 0:Ts:Nbits-Ts;
        S_N           = conv(bits_s,g_t);
        S_N           = S_N(1:length(t_seq));
        Phi_N         = cumsum(S_N)*Ts;
        mod_signal    = exp(1i*2*pi*modulation_index*Phi_N);
        %----------------------------------------------------------------------------
        received_signal      = awgn(mod_signal((pulse_length-1)*os+1:end), snr(snr_no)-10*log10(os/log2(M_ary)), 'measured');
        %------------------ Viterbi Decoding ------------------------------ 
        DetectedBits         = zeros(1,Nbits-decision_delay);
        DBits_idx            = 1; % index of the detected bit;
        pathmetric           = zeros(States_Number,1);
        pathmetric_n         = zeros(States_Number,1);
        SurvivorPath         = zeros(States_Number,decision_delay);
        for n = 1:Nbits-pulse_length
                Branche_metrics_cum   = zeros(Branche_Number,1);
                for i = 1:Branche_Number
                    if(n~=1)
                        %                         z = conv(r((n-1)*os:n*os), conj(flip(Branche_metrics(i,:))));
                        %                         Branche_metrics_cum(i) = real(z(os)*exp(-1i*Branches(i,1)));
                        % Check equation 3.14 p.22 from "Comparison of Noncoherent detectors for SOQPSK and GMSK in Phase Noise Channels"
                        z = real(exp(-1i*Branches(i,1))*trapz(received_signal((n-1)*os:n*os).*conj(Branche_metrics(i,:)))*Ts);
                        Branche_metrics_cum(i) = z;
                    else
                        %                         z = conv(r(1:1*os), conj(flip(Branche_metrics(i,:))));
                        %                         Branche_metrics_cum(i) = real(z(os)*exp(-1i*Branches(i,1)));
                        % Check equation 3.14 p.22 from "Comparison of Noncoherent detectors for SOQPSK and GMSK in Phase Noise Channels"
                        z = real(exp(-1i*Branches(i,1))*trapz(received_signal(1:1*os+1).*conj(Branche_metrics(i,:)))*Ts);
                        Branche_metrics_cum(i) = z;
                    end
                end
                j = 1;
                Branche_compare = zeros(1,M_ary);
                % 1- Iterate over all branches entering the same node, and
                % select the branch with the highest metric. ('Branche_compare') 
                % 2- Next, save the pathmetric of the best branches for 
                % the computation of the next branches transition.('pathmetric_n', 'pathmetric')
                % 3- Save the states in the 'SurvivorPath' table
                % for the detection afterward. Each value in the
                % 'SurvivorPath' presents the previous states and their indexes 
                %  present the current states (new state). 
                for i=1:M_ary:size(Branches,1)
                    for k = 0:M_ary-1
                        Branche_compare(k+1) = Branche_metrics_cum(states_sort(i+k,3)) + pathmetric(states_sort(i+k,1));
                    end
                    [pathmetric_n(j,1),idx] = max(Branche_compare);
                    idx = idx-1;
                    if(n<=decision_delay)
                        SurvivorPath(states_sort(i+idx,2),n) = states_sort(i+idx,1); % value = previous state (i+idx), index = new state of the survived path
                    else
                        if(i==1)
                            SurvivorPath(:,1:end-1) =  SurvivorPath(:,2:end); % shift one colum to the left.
                        end
                        SurvivorPath(states_sort(i+idx,2),end) = states_sort(i+idx,1);
                    end
                    j=j+1;
                end
                pathmetric = pathmetric_n;
            if(n>decision_delay)
                %-------------------- trace back unit----------------------------
                [maxPath idxpath] = max(pathmetric);
                currState         = idxpath;
                
                for jj = decision_delay:-1:1
                    prevState         = SurvivorPath(currState,jj);
                    if(jj>1)
                        currState         = prevState;
                    end
                end
                DetectedBits(DBits_idx) = Branches(states_sort(find(ismember(states_sort(:,1:2),[prevState currState],'rows')),3),2);
                DBits_idx = DBits_idx+1;
            end
            
        end
        if(M_ary==2)
            BER(snr_no)       = nnz(bits_m(2+pulse_length:end-decision_delay-2-pulse_length)-DetectedBits(1+pulse_length:end-3-pulse_length));
        elseif(M_ary==4)
            BER(snr_no)       = nnz(bits_m(2+3*pulse_length:end-decision_delay-3-3*pulse_length)-DetectedBits(1+3*pulse_length:end-4-3*pulse_length));
        elseif(M_ary==8)
            BER(snr_no)       = nnz(bits_m(2+3*pulse_length:end-decision_delay-3-3*pulse_length)-DetectedBits(1+3*pulse_length:end-4-3*pulse_length));
        end
        
        BER_all(snr_no) = BER_all(snr_no)  + BER(snr_no);
        sprintf('BER_all for division number : %d / %d for the snr %d',I,Nbits_div,snr(snr_no))
        if(I==Nbits_div) % at the last frame
            if(M_ary==2)
                BER_all(snr_no) = BER_all(snr_no)/(length(bits_m(2+pulse_length:end-decision_delay-2-pulse_length))*Nbits_div);
            elseif (M_ary==4)
                BER_all(snr_no) = BER_all(snr_no)/(length(bits_m(2+3*pulse_length:end-decision_delay-3-3*pulse_length))*Nbits_div);
            elseif (M_ary==8)
                BER_all(snr_no) = BER_all(snr_no)/(length(bits_m(2+3*pulse_length:end-decision_delay-3-3*pulse_length))*Nbits_div);
            end
            sprintf('BER is %f',BER_all(snr_no))
        end
        
    end
end
toc
%----------------------------------Plots----------------------------------------------
figure(1)
semilogy(snr, BER_all, 'linewidth', 3);hold on;
semilogy(snr,Pe,'linewidth', 3);
xlabel('E_b/N_0 (dB)'); 
ylabel('BER');
legend('Simulation', 'Bound');
grid on;







































