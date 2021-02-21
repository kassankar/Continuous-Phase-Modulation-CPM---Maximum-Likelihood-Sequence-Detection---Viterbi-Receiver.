function MAIN = MainFunctions
MAIN.CREATECPMPULSE                    = @CREATECPMPULSE;
MAIN.PERMUREPET                        = @PERMUREPET;
MAIN.states_Branches                   = @states_Branches;
MAIN.states_transition                 = @states_transition;
MAIN.Mapp                              = @Mapp;
MAIN.deMapp                            = @deMapp;
MAIN.compute_Branche_metric            = @compute_Branche_metric;
end

function [g_t,q_t,t] = CREATECPMPULSE(pulse,pulse_length,w,freq_sa,img)
g_t              = 0;
q_t              = 0;
time_sa          = 1/freq_sa;
switch             pulse
    
    case    1
        t          = -pulse_length/2:time_sa:pulse_length/2;
        t0         = 0;
        g_t        = (2*w)./((t-t0).^2+w^2);
        g_t        = g_t*(1/(4*pi));
        Cst        = sum(g_t)*time_sa;
        nug_t      = (0.5) / Cst;
        g_t        = nug_t*g_t;
        q_t        = cumtrapz(g_t)*time_sa;
        
        
    case    2
        t           = -pulse_length/2:time_sa:pulse_length/2;
        BT          = 0.3;
        alpha       = 2*pi*BT/(sqrt(log(2)));
        gauss       = qfunc(alpha*(t-0.5)) - qfunc(alpha*(t+0.5));
        Cst         = 0.5/(sum(gauss)*time_sa);
        g_t         = Cst*gauss;
        q_t         = cumtrapz(g_t)*time_sa;
        
    case    3
        t          = 0:time_sa:pulse_length;
        g_t        = (1/(2*pulse_length).*(1- cos(2*pi.*t/(pulse_length))));
        K          = 0.5/(sum(g_t)*time_sa);
        g_t        = K*g_t;
        q_t        = cumtrapz(g_t)*time_sa;
        
        
    case   4
        t                      = 0:time_sa:pulse_length;
        g_t                    = 1/(2*pulse_length)*ones(1,length(t));
        g_t(1)                 = 0;
        g_t(end)               = 0;
        K                      = 0.5/(sum(g_t)*time_sa);
        g_t                    = K*g_t;
        q_t                    = cumtrapz(g_t)*time_sa;
        
    otherwise
        disp('no pulses for pulse >4')
end



if(img ==1)
    figure (1)
    % subplot(1,2,1)
    if(pulse==1)
        plot(t,4*pi*g_t)
    else
        plot(t,g_t)
    end
    
    xlabel("time(t/T_s)")
    ylabel("Freqeuncy pulse g(t)")
    grid on;
    set(gca,'FontName','Arial','FontSize',12)
    figure (2)
    plot(t,q_t)
    xlabel("time")
    ylabel("q(t)")
    grid on;
else
end


end
function [M, I] = PERMUREPET(V, N, K)
% PERMN - permutations with repetition
%   Using two input variables V and N, M = PERMN(V,N) returns all
%   permutations of N elements taken from the vector V, with repetitions.
%   V can be any type of array (numbers, cells etc.) and M will be of the
%   same type as V.  If V is empty or N is 0, M will be empty.  M has the
%   size numel(V).^N-by-N. 
%
%   When only a subset of these permutations is needed, you can call PERMN
%   with 3 input variables: M = PERMN(V,N,K) returns only the K-ths
%   permutations.  The output is the same as M = PERMN(V,N) ; M = M(K,:),
%   but it avoids memory issues that may occur when there are too many
%   combinations.  This is particulary useful when you only need a few
%   permutations at a given time. If V or K is empty, or N is zero, M will
%   be empty. M has the size numel(K)-by-N. 
%
%   [M, I] = PERMN(...) also returns an index matrix I so that M = V(I).
%
%   Examples:
%     M = permn([1 2 3],2) % returns the 9-by-2 matrix:
%              1     1
%              1     2
%              1     3
%              2     1
%              2     2
%              2     3
%              3     1
%              3     2
%              3     3
%
%     M = permn([99 7],4) % returns the 16-by-4 matrix:
%              99     99    99    99
%              99     99    99     7
%              99     99     7    99
%              99     99     7     7
%              ...
%               7      7     7    99
%               7      7     7     7
%
%     M = permn({'hello!' 1:3},2) % returns the 4-by-2 cell array
%             'hello!'        'hello!'
%             'hello!'        [1x3 double]
%             [1x3 double]    'hello!'
%             [1x3 double]    [1x3 double]
%
%     V = 11:15, N = 3, K = [2 124 21 99]
%     M = permn(V, N, K) % returns the 4-by-3 matrix:
%     %        11  11  12
%     %        15  15  14
%     %        11  15  11
%     %        14  15  14
%     % which are the 2nd, 124th, 21st and 99th permutations
%     % Check with PERMN using two inputs
%     M2 = permn(V,N) ; isequal(M2(K,:),M)
%     % Note that M2 is a 125-by-3 matrix
%
%     % PERMN can be used generate a binary table, as in
%     B = permn([0 1],5)  
%
%   NB Matrix sizes increases exponentially at rate (n^N)*N.
%
%   See also PERMS, NCHOOSEK
%            ALLCOMB, PERMPOS, NEXTPERM, NCHOOSE2 on the File Exchange
% tested in Matlab 2018a
% version 6.2 (jan 2019)
% (c) Jos van der Geest
% Matlab File Exchange Author ID: 10584
% email: samelinoa@gmail.com
% History
% 1.1 updated help text
% 2.0 new faster algorithm
% 3.0 (aug 2006) implemented very fast algorithm
% 3.1 (may 2007) Improved algorithm Roger Stafford pointed out that for some values, the floor
%   operation on floating points, according to the IEEE 754 standard, could return
%   erroneous values. His excellent solution was to add (1/2) to the values
%   of A.
% 3.2 (may 2007) changed help and error messages slightly
% 4.0 (may 2008) again a faster implementation, based on ALLCOMB, suggested on the
%   newsgroup comp.soft-sys.matlab on May 7th 2008 by "Helper". It was
%   pointed out that COMBN(V,N) equals ALLCOMB(V,V,V...) (V repeated N
%   times), ALLCMOB being faster. Actually version 4 is an improvement
%   over version 1 ...
% 4.1 (jan 2010) removed call to FLIPLR, using refered indexing N:-1:1
%   (is faster, suggestion of Jan Simon, jan 2010), removed REPMAT, and
%   let NDGRID handle this
% 4.2 (apr 2011) corrrectly return a column vector for N = 1 (error pointed
%    out by Wilson).
% 4.3 (apr 2013) make a reference to COMBNSUB
% 5.0 (may 2015) NAME CHANGED (COMBN -> PERMN) and updated description,
%   following comment by Stephen Obeldick that this function is misnamed
%   as it produces permutations with repetitions rather then combinations.
% 5.1 (may 2015) always calculate M via indices
% 6.0 (may 2015) merged the functionaly of permnsub (aka combnsub) and this
%   function
% 6.1 (may 2016) fixed spelling errors
% 6.2 (jan 2019) fixed some coding style warnings
narginchk(2, 3) ;
if fix(N) ~= N || N < 0 || numel(N) ~= 1
    error('permn:negativeN','Second argument should be a positive integer') ;
end
nV = numel(V) ;
if nargin==2 
    %% PERMN(V,N) - return all permutations
    if nV == 0 || N == 0
        M = zeros(nV, N) ;
        I = zeros(nV, N) ;
    elseif N == 1
        % return column vectors
        M = V(:) ;
        I = (1:nV).' ;
    else
        % this is faster than the math trick used with 3 inputs below
        [Y{N:-1:1}] = ndgrid(1:nV) ;
        I = reshape(cat(N+1, Y{:}), [], N) ;
        M = V(I) ;
    end
else
    %% PERMN(V,N,K) - return a subset of all permutations
    nK = numel(K) ;
    if nV == 0 || N == 0 || nK == 0
        M = zeros(numel(K), N) ;
        I = zeros(numel(K), N) ;
    elseif nK < 1 || any(K<1) || any(K ~= fix(K))
        error('permn:InvalidIndex','Third argument should contain positive integers.') ;
    else
        V = reshape(V, 1, []) ; % v1.1 make input a row vector
        nV = numel(V) ;
        Npos = nV^N ;
        if any(K > Npos)
            warning('permn:IndexOverflow', ...
                'Values of K exceeding the total number of combinations are saturated.')
            K = min(K, Npos) ;
        end
             
        % The engine is based on version 3.2 with the correction
        % suggested by Roger Stafford. This approach uses a single matrix
        % multiplication.
        B = nV.^(1-N:0) ;
        I = ((K(:)-.5) * B) ; % matrix multiplication
        I = rem(floor(I), nV) + 1 ;
        M = V(I) ;
    end
end
% Algorithm using for-loops
% which can be implemented in C or VB
%
% nv = length(V) ;
% C = zeros(nv^N,N) ; % declaration
% for ii=1:N,
%     cc = 1 ;
%     for jj=1:(nv^(ii-1)),
%         for kk=1:nv,
%             for mm=1:(nv^(N-ii)),
%                 C(cc,ii) = V(kk) ;
%                 cc = cc + 1 ;
%             end
%         end
%     end
% end

end
function [States_Number,Branche_Number,phase_states] = states_Branches(modulation_index,pulse_length,M_ary)
% Function return the number of states and banches needed alongside the
% phase states. Check Proakis page 247-248 to understand how we compute the
% phase states and branches.
[m, p] = rat(modulation_index);
    if(mod(m,2)==0)
        States_Number  = p*(M_ary)^(pulse_length-1);
        Branche_Number = p*(M_ary)^(pulse_length);
        fprintf('m is even, number of states is: %d \n',States_Number);
        phase_states = mod((0:(pi*m)/p:((p-1)*pi*m)/p),2*pi);
    else
        States_Number  = 2*p*(M_ary)^(pulse_length-1);
        Branche_Number = 2*p*(M_ary)^(pulse_length);
        fprintf('m is odd, number of states is: %d \n', States_Number);
        phase_states = mod((0:(pi*m)/p:((2*p-1)*pi*m)/p),2*pi);
    end
end
function [Branches,states_sort] = states_transition(States_Number,Branche_Number,phase_states,pulse_length,M_ary,A_m,modulation_index)
% Function return the states trellis.
% states_sort = |state1| -- |state2| -- |Branche|. The states are sorted
% where all branches entering the same node are one after another.
state_vector_number     = M_ary^(pulse_length-1);
[state_vector,~]        = PERMUREPET(A_m,pulse_length-1);            % Return all different combinations of bits for each state ;
states                  = zeros(States_Number,(pulse_length-1)+1);
states(1:end,1)         = repelem(phase_states,state_vector_number);      % repeat each phase for each bits combination
states(1:end,2:end)     = repmat(state_vector,length(phase_states),1);

Branche_vector_number   = M_ary^(pulse_length);
[Branche_vector,~]      = PERMUREPET(A_m,pulse_length);              % Return all different combinations of bits for each state ;
Branches                = zeros(Branche_Number,pulse_length+1);
Branches(1:end,1)       = repelem(phase_states,Branche_vector_number);    % repeat each phase for each bits combination
Branches(1:end,2:end)   = repmat(Branche_vector,length(phase_states),1);

%-------------------- states transition ---------------------
states_transitions       = zeros(length(Branches),3);
for i = 1:length(Branches)
    [~,~,states_transitions(i,1)]   = intersect(Branches(i,1:pulse_length),states,'rows');
    % For the next_phase, check Proakis p.248 (\theta_{n+1} = \theta_{n} + \pi*....)
    Next_phase                      = mod(Branches(i,2)*pi*modulation_index+Branches(i,1),2*pi); 
    %-------------- correction on the phase --------------------
    for j = 1:length(phase_states)
        if(abs(Next_phase-phase_states(j))<=10^-5)
            Next_phase = phase_states(j);
            break;
        elseif(abs(Next_phase-2*pi)<=10^-5)
            Next_phase = 0;
            break;
        end
    end
    %-----------------------------------------------------------
    [~,~,states_transitions(i,2)]   = intersect([Next_phase Branches(i,3:end)],states,'rows');
    states_transitions(i,3)         = i;
end
[~,idx]              = sort(states_transitions(:,2), 'ascend');
states_sort          = states_transitions(idx,:);                % Comparable Branches are sorted one after another.
end
function [bits_mod] = Mapp(bits, M_ary, pulse)
if(pulse~=1)
    if(M_ary==2)
        bits_mod = 2*bits-1;
    else
        a = reshape(bits,log2(M_ary),[]).';
        b = bi2de(a,'left-msb')+1;
        d = 2*b-1-M_ary;
        bits_mod = d.';
    end
else
    if(M_ary==2)
      bits_mod = bits;
    else
        a = reshape(bits,log2(M_ary),[]).';
        b = bi2de(a,'left-msb');
        bits_mod = b.';
    end
    
end
end
function [bits_dem] = deMapp(bits_mod, M_ary, pulse)
if(pulse~=1)
    if(M_ary==2)
        bits_dem = (bits_mod+1)/2;
    else
        a = bits_mod.';
        b = (a+1+M_ary)/2-1;
        c = de2bi(b,'left-msb');
        c = c.';
        bits_dem = reshape(c,1,[]);
    end
else
     if(M_ary==2)
        bits_dem = bits_mod;
     else
        a = bits_mod.';
        c = de2bi(a,'left-msb');
        c = c.';
        bits_dem = reshape(c,1,[]);
    end
    
end

end
function [psi] = compute_Branche_metric(os, modulation_index, pulse_length, alpha, g_t)
Ts             = 1/os;
Nbits          = length(alpha);
bits_s         = upsample(alpha,os);
if(pulse_length==1)
t_seq          = 0:Ts:Nbits;
else
t_seq          = 0:Ts:Nbits-Ts;
end
S_N            = conv(bits_s,g_t);
S_N            = S_N(1:length(t_seq));
Phi_N          = cumsum(S_N)*Ts;
if(pulse_length>1)
psi            = exp(1i*2*pi*modulation_index*Phi_N((pulse_length-1)*os:pulse_length*os));
else
    psi        = exp(1i*2*pi*modulation_index*Phi_N(1:pulse_length*os+1));
end

end







