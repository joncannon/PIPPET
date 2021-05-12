function params = sound_params()

    % collect a set of sounds and assign to indices
    params = struct();

    params.Fs = 44100;
    components_path = 'stimulus_components/'
    tick = audioread(strcat(components_path, 'wood_tick.wav'));
    beep = audioread(strcat(components_path, 'shortbeep.wav'));
    
    params.sound_list{1} = tick(:,1);
    params.standard_index = 1;
    
    params.sound_list{2} = resample(tick(:,1), 3,2);
    params.high_index = 2;

    params.sound_list{3} = 10*beep(:,1);
    params.beep_index = 3;
    
    
    params.sound_list{4} = 0;
    params.omission_index = 4;
