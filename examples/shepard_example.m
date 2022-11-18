% shepard synthesis and modeling minimal example 

%% model parameters 
syn.fs = 44100; % audio sampling freq
mod_par.numb_fb = 128; % number of bands
mod_par.gtf = gammatoneFilterBank([20 16000],mod_par.numb_fb, syn.fs); % gammatone filterbank object 
mod_par.f_erb = getCenterFrequencies(mod_par.gtf); % center freqs 
mod_par.erb_thresh = -30;  % level threshold below which signal is not considered 
mod_par.noise_lev.ac = 0; % % internal noise of ac cue [set to 0]
mod_par.noise_lev.cc = 0; % internal noise of cc cue
mod_par.ac.fmin = 64; % min f0 search range 
mod_par.ac.fmax = 128; % max f0 search range 
mod_par.search_cc_lag = 5; % CC search range, here 5 steps 
mod_par.f_cutoff = 10*[2^(6+.5), 2^(6+.5)]; % model cutoffs between res and unres

%% synthesis parameters
syn.dur = .125; % sound duration
syn.fade = .005; % fade in/out duration
syn.sil_dur = .125; % inter-stimulus interval
syn.fs = 44100; 

% spectral parameters for additive synth
syn.no_sines = 500; % % number of sines at max
syn.fmax = 20000; % max frequency used 
syn.f0_base = 64; % F0
syn.fmin= syn.f0_base; % min frequency  
        
% envelope parameters (fixed across trials)
syn.global_std = 1; % global std (in octaves)
syn.global_mu = 3; % position of mean of global spectrum (counted from f0 in octaves) 
syn.local_std = .15; % local std
syn.local_mu_spac = 1; % spacing between local ripples in octaves
syn.local_mu_offset = 0; % position of first local ripple mean

% fine structure parameters
syn.harmonicity = 'harmonic'; % choose between 'harmonic' and 'inharmonic'
syn.odd_atten = 2; % choose attenuation factor per step (shift 0 = maximal attenuation to be equivalent with shift 12)
syn.phase_jitter = 0; % phase jitter, vary from 0 (all cosine phase) to 1 (uniformly distributed in 0-2pi)

%% sound synthesis 

% intialize independent variables studied in the manuscript 
syn.f0_base = 2^(6+rand(1)); % F0 reference point, randomize between 64 and 128 Hz
syn.freq_jitter = 0; % harmonicity (0: no jitter / harmonic, 1: full jitter / inharmonic)
syn.pitch_shifts = [0 0]; % this is a 2 st shift in SFS [for historic reasons, the variable is called "pitch_shift"]
syn.env_shifts = [0 3]; % this is a 0 st shift in SE [for historic reasons, the variable is called "env_shift"]

% synthesis
[sig, shep_syn] = shepard_spectrum_2d(syn); % generate signal 

% splice signals into two 
t1 = round(syn.sil_dur*syn.fs - syn.fs*.01); t2 = round(t1 + syn.dur*syn.fs + syn.fs*.01);
t3 = round(t2 + syn.sil_dur*syn.fs - syn.fs*.01); t4 = round(t3 + syn.dur*syn.fs + syn.fs*.01); 
sig = .1*sig/max(abs(sig)); 
sig1 = sig(t1:t2);
sig2 = sig(t3:t4);

%% apply model 

[shift_stats] = shepard_model_2d(sig1, sig2, syn.fs, mod_par); % apply the model 
shift_stats.ac.p1 % extracted peak of the SACF function 
shift_stats.cc_low.centroid % centroid of CC function of low frequencies (CCres) 
shift_stats.cc_high.centroid % centroid of CC function of high frequencies (CCunres)

%% visualize results 

cmap = colormap('lines'); % colormap
S1 = 20*log10(abs(fft(sig1))); % spectral analysis for plot 
S2 = 20*log10(abs(fft(sig2))); 

fvec = linspace(0, syn.fs/2, length(S1)/2); % frequency range 
fmax = 8000; % plotting range < 8k 
plot_range = find(diff(fvec<fmax)); % plotting range < 8k 

figure; 
% plot spectral analysis 
semilogx(fvec(1:plot_range), S1(1:plot_range), 'color', cmap(1,:)); hold on;
semilogx(fvec(1:plot_range), S2(1:plot_range), 'color', cmap(2,:));
xlabel('Frequency [Hz]')
ylabel('Level [dB]')
xlim([32 fmax])

% plot AC estimates 
xline(shift_stats.ac.p1, 'linewidth', 4, 'color', cmap(1,:)) % F0 estimates
xline(shift_stats.ac.p2, 'linewidth', 4, 'color', cmap(2,:))
ACdiff = shift_stats.ac.p1 - shift_stats.ac.p2; 
acdiff = strcat('F0diff: ', num2str(ACdiff), ' Hz'); % resolved 
text(200, -60, acdiff , 'fontsize', 12) % difference in F0 estimates 

% write CC results 
updown = {'up' 'down'}; % CC centroid positive: down; CC centroid negative: up 
ccres = strcat('CCres: ', updown(1.5 + sign(shift_stats.cc_low.centroid)/2)); % resolved 
text(200, -70, ccres , 'fontsize', 12) % indicated shift direction for resolved CC cue

ccunres = strcat('CCunres: ', updown(1.5 + sign(shift_stats.cc_high.centroid)/2)); % resolved 
text(200, -80, ccunres , 'fontsize', 12)  % indicated shift direction for unresolved CC cue

legend('Signal 1', 'Signal 2', 'F0(Sig1)', 'F0(Sig2)')
