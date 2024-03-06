function [seqC_A]=GenerateCACode(num)

    % num = index of the code ( SV index)
    % The ranging codes for the signal C/A are GOLD Codes. They are generated
    % by a modulo-2 addition of two sub-sequences. 
    % This sequences is equivalent with 1ms
    % Initialization polynomial sequences 1111111111
    % G1 = X10 + X3 + 1
    % G2 = X10 + X9 + X8 + X6 + X3 + X2 + 1
    G1 = [1 0 0 0 0 0 0 1 0 0 1];
    G2 = [1 1 1 0 1 0 0 1 1 0 1];
    Np = 1023;
    SV_TapSelection=Select_Code_tap(num);
    G1sequences = prilfsr_new(G1,Np);
    G2sequences = seclfsr_new(G2,Np,SV_TapSelection);
    seqC_A = (xor(G1sequences,G2sequences));
end

