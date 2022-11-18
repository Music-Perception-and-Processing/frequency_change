%% participants stats
clear all

nums_parts = [13 12 12];

for nExp = 1:3
    path = strcat('/Users/kaisiedenburg/Dropbox (Personal)/2dShepard/ExperimentScriptsData/Data_Exp', ...
    num2str(nExp), '/Fragebogen'); 
    cd(path)
    nums = nums_parts(nExp);
    clear gr
    for kk = 1:nums_parts(nExp)
        kk
        pp(kk,nExp).x = load(strcat('proband_', num2str(kk), '.mat'));
        gr.codeword{kk} = pp(kk,nExp).x.codeword;
        gr.age(kk) = str2num(pp(kk,nExp).x.alter);
        if isempty(pp(kk,nExp).x.mt)
            gr.mt(kk) = 0;
        else
            gr.mt(kk) = pp(kk,nExp).x.mt;
        end 
        if isempty(pp(kk,nExp).x.ins_training_p) % years of training on instrument
            gr.instr_train(kk) = 0;
        else
            gr.instr_train(kk) = pp(kk,nExp).x.ins_training_p/1.67;
        end
        if isempty(pp(kk,nExp).x.daily_practise_p)  % average hours of practice per day
            gr.practice(kk) = 0;
        else
            gr.practice(kk) = pp(kk,nExp).x.daily_practise_p/1.57;
        end
        if isempty(pp(kk,nExp).x.instrumente_p)  % number of instruments
            gr.instruments(kk) = 0;
        else
            gr.instruments(kk) = pp(kk,nExp).x.instrumente_p/.82;
        end
        if isempty(pp(kk,nExp).x.practise_p)  % peak of interest practice
            gr.peak_practice(kk) = 0;
        else
            gr.peak_practice(kk) = pp(kk,nExp).x.practise_p/.71;
        end
        if isempty(pp(kk,nExp).x.music_training_p)  % years of music theory training 
            gr.theory(kk) = 0;
        else
            gr.theory(kk) = pp(kk,nExp).x.music_training_p/1.43;
        end
    end
    
    if nExp == 1
        parts = [1 2 4:13];
    else
        parts = [1:12]; 
    end
    allstat(nExp).gr = gr;
    allstat(nExp).parts = parts;
    % age
    pstat.age.med(nExp,1) = median(gr.age(parts)); 
    pstat.age.range(nExp,:) = [min(gr.age(parts)), max(gr.age(parts))];

    % Gold MSI
    pmus = gr.mt > 0; % musical participants
    pstat.mus(nExp,1) = sum(pmus);
    pstat.mt.med(nExp,1) = median(gr.mt(pmus));
    pstat.mt.range(nExp,:) = [min(gr.mt(pmus)), max(gr.mt(pmus))];

    % instr train 
    pstat.instr.med(nExp,1) = median(gr.instr_train(pmus));
    pstat.instr.range(nExp,:) = [min(gr.instr_train(pmus)), max(gr.instr_train(pmus))];
end

%% check out whether there was overlap across experiments

allstat(1).gr.codeword(allstat(1).parts)
allstat(2).gr.codeword
allstat(3).gr.codeword

names = fieldnames(allstat(1).gr);
X = []; 
for nExp = 1:3
    for nField = 2:length(names) % leave out codeword 
        X(:, nField, nExp) = allstat(nExp).gr.(names{nField})(allstat(nExp).parts);
    end
end

% ch7428 = CH07RG28? candidate: 8, 5, 6
pp(9,1).x
pp(5,2).x
pp(6,3).x
% seems likely 

% EN5TH13 = EN5MT13!
pp(9,2).x
pp(4,3).x



% overlap participants: 
% exp 1, 2, & 3: ER04ER24, TZ07ST14, KI06NS22, CH07RG28
% exp 1 & 2: MI08LF31
% exo 1 & 3: AU05AS28
% exp 2 & 3: EN5TH13 = EN5MT13!

% 4 participants completed all three 
% 1 completed each other pair of experiment 
fab4 = {'ER04ER24', 'TZ07ST14', 'KI06NS22', 'CH07RG28'};
fab4ind = ...
[1 5 11 9;
 1 3 4 5; 
 1 3 5 6]


