function [sig, syn] = shepard_spectrum_2d(syn)
% synthesize 1d or 2d shepard shifts where synthesis parameters are specified in the structure syn


% preliminaries
if isfield(syn, 'cyclic')
    else 
        syn.env_shifts = mod(syn.env_shifts, 12); % ensure cyclic nature of envelope shift, then below the spectral range doesnt need to be so high
        syn.pitch_shifts = mod(syn.pitch_shifts, 12); % ensure cyclic nature of pitch shift
end

syn.pitch_shifts_rat = 2.^(syn.pitch_shifts/12); % from semitones to freq ratios
t = linspace(0, syn.dur, syn.fs*syn.dur)'; % time signal
sig = [];
sil = zeros(round(syn.sil_dur*syn.fs),1); % inter stimulus interval
numb_octaves = round(log2(syn.fmax/syn.f0_base)); % number of octaves required
syn.pars = [];

switch syn.harmonicity % differentiate harmonic and inharmonic fine structure
        case 'harmonic' 
            jit_harm = syn.f0_base*syn.freq_jitter*(rand(1,round(syn.fmax/syn.f0_base))-.5); % similar to McPherson/McDermott jitter 
        case 'inharmonic' 
            inharm_len = length([0:1/syn.comp_oct:log2(round(syn.fmax/syn.f0_base))-1/syn.comp_oct]);
            jit_inharm = syn.f0_base*syn.freq_jitter*(rand(1,inharm_len)-.5); % similar to McPherson/McDermott jitter 
end

for l = 1:length(syn.env_shifts) % compute each shift separately
    % compute spectral fine structure 
    spacing_ok = 0; % to control jitter to have differences between adjacent harmonics more than 30Hz
    switch syn.harmonicity % differentiate harmonic and inharmonic fine structure
        case 'harmonic' 
            harms_lin = [1:round(syn.fmax/syn.f0_base)]; % linearly spaced harmonics
            %while spacing_ok == 0
            lin_harmonics = jit_harm + syn.f0_base*syn.pitch_shifts_rat(l)*harms_lin; % harmonic frequencies
                %spacing_ok = min(diff(lin_harmonics)) > syn.jitt_adj_min; % condition (see McPherson paper)
            %end
        case 'inharmonic'
            %harms_log = 2.^linspace(0,log2(round(syn.fmax/syn.f0_base)), numb_octaves*syn.comp_oct); % log frequency
            harms_log = 2.^[0:1/syn.comp_oct:log2(round(syn.fmax/syn.f0_base))-1/syn.comp_oct]; % log frequency
            %while spacing_ok == 0
            lin_harmonics = jit_inharm + syn.f0_base*syn.pitch_shifts_rat(l)*harms_log; % logarithmic frequencies
            %    spacing_ok = min(diff(lin_harmonics)) > syn.jitt_adj_min; % condition (see McPherson paper)
            %end
    end
    syn.pars.f(:,l) = lin_harmonics'; 
    log_harmonics = log2(lin_harmonics); % go to log domain 
    
    % compute global envelope
    global_amps = normpdf(log_harmonics, log2(syn.f0_base)+syn.global_mu, syn.global_std); % global gaussian shaped envelope
    global_amps = global_amps/max(global_amps); % normalize

    % compute local envelopes 
    clear local_amps
    KK = 12; % number of ripples (to be shifted, and shaped by global spec env)
    for k = 1:KK % shift in octave units
        local_mu = log2(syn.f0_base) + syn.local_mu_offset + (k-1)*syn.local_mu_spac + syn.env_shifts(l)/12; 
        local_amps(:,k) = normpdf(log_harmonics, local_mu, syn.local_std); 
    end
    local_amps = local_amps/max(local_amps(:)); % peak normalize
    syn_amps = global_amps'.*sum(local_amps,2); % combine local and global amplitudes
    syn.pars.global_amps = global_amps;
    syn.pars.local_amps = local_amps; 
    syn.pars.syn_amps = syn_amps; 
    
    % make odd attenuation, if desired
    if strcmp(syn.harmonicity, 'harmonic') % works only for 'hamonic' version of syn.harmonicity
       atten_fact = ((10^(syn.odd_atten/20))^(11-syn.pitch_shifts(l))); % see deutsch et al 2008 for the 3.5 dB
       syn_amps(1:2:end) =  syn_amps(1:2:end)/atten_fact; 
    end
    %syn_amps = syn_amps/rms(syn_amps); % normalize amplitudes
    syn_amps(lin_harmonics >= syn.fs/2) = 0; % prevent aliasing
    syn.pars.final_amps(:,l) = syn_amps; 
    
    % not for 2d shepard any more (JASA EL paper): 
    %syn_amps(1) = .2; % to ensure a (nonvirtual) pitch percept at f0 % 
    
    % put things together 
    syn_carr = cos(2*pi*t*lin_harmonics + repmat(syn.phase_jitter*2*pi*rand(size(lin_harmonics)), length(t),1)); % matrix of carrier signals with random phase
    sig_new = cos_ramp(syn_carr*syn_amps, syn.fs, syn.fade); % weighted combination of carriers 
    sig_new = sig_new/rms(sig_new); % normalize
    if isnan(syn.pitch_shifts(l)) && isnan(syn.env_shifts(l))
        sig_new = sil;
    end
    sig = [sig; sil; sig_new]; 
end
sig = [sig; sil; sil];
