function snd_total=stim_maker_etime(times, identities, params)

Fs = params.Fs;

snd_total=zeros(floor((max(times)+0.5)*Fs), 2);

sound_list = params.sound_list;

for i=1:length(times)
    pointer = floor(Fs * times(i));
    
    snd_num = identities(i);

    snd_total(1+pointer:pointer+length(sound_list{snd_num}), 1) = sound_list{snd_num} + snd_total(1+pointer:pointer+length(sound_list{snd_num}), 1);
    
end

snd_total(:, 2)=snd_total(:, 1);

