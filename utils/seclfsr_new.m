function [s]=seclfsr_new(con,epoch,TapSel)
    % Implements a LFSR behaviour just giving the tap polynomial (in binary).
    % It is truncated with a period of "epoch" chips.
    % Start values are always 1 when is reset
    L=length(con);
    R=L-1;                                      % number of registers
    r=true(1,R);                                % initial state: all 1's
    tap=fliplr(con);                            % tap polynomial:1 x (R+1)
    nt=length(tap);                             % number of tap
    nr=R;                                       % number of registers
    seq=false(1,epoch);                         % Sequences
    for i=1:epoch
        seq(i)=xor(r(TapSel(1)),r(TapSel(2)));
        aux=r(nr);                              
        for j=1:(R-1)
            if (tap(nt-j)==1)
                aux=xor(r(nr-j),aux);
            else
                aux=xor(0,aux);
            end
            r(nr-j+1)=r(nr-j);
        end
        r(1)=aux;
    end
    s=seq;