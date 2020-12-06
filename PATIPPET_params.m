% © Jonathan Cannon, MIT, 2020
% Sets default or custom parameters for PIPPET model.
% All inputs are optional key-value pairs that replace default values.


function out = PATIPPET_params(varargin)

p = inputParser;

addParameter(p,'n_streams',1,@isnumeric);               % Number of distinct event streams

% Inputs in this block can be given either as scalars/arrays to pertain to
% all event streams, or as cell arrays to pertain to each event stream
% individually.
addParameter(p,'means_unit',[.25],@isnumeric);          % One unit of repeating pattern of expected event times
addParameter(p,'variance_unit',[.001],@isnumeric);      % One unit of repeating pattern of event expectation variances
addParameter(p,'tau_unit',[0.02],@isnumeric);           % One unit of repeating pattern of event expectation strengths
addParameter(p,'tau_0',0.01,@isnumeric);                % tau_0
addParameter(p,'expected_cycles',4,@isnumeric);         % Number of repetitions of pattern
addParameter(p,'expected_period',.25,@isnumeric);       % Period of pattern repetition
addParameter(p,'highlight_event_indices',[]);           % Display weights for lines marking expected timepoints
addParameter(p,'highlight_expectations',[]);            % Display weights for lines marking event times
addParameter(p,'event_times',[1],@isnumeric);           % Observed event times (should be given as cell array if there are multiple streams!)

addParameter(p,'tmax',nan,@isnumeric);                  % Max simulation time (default setting is based on event times and expected event times)
addParameter(p,'dt',0.001,@isnumeric);                  % Integration time step
addParameter(p,'xbar_0',[0; 1],@isnumeric);             % Initial estimated phase and tempo
addParameter(p,'Sigma_0',[.0001,0; 0,.04],@isnumeric);  % Initial covariance matrix
addParameter(p,'sigma_phi',0.05,@isnumeric);            % Generative model phase noise
addParameter(p,'display_phasetempo', true);             % Display ellipse graphic showing joint phase-tempo posteriors over time
addParameter(p,'display_phase', false);                 % Display graphic showing evolution of phase posterior over time
addParameter(p,'display_tempo', false);                 % Display graphic showing evolution of tempo posterior over time

addParameter(p,'phimax',nan,@isnumeric);                % Max phase for display (default setting is based on expected event times)
addParameter(p,'sigma_theta',0.05,@isnumeric);          % Generative model tempo noise
addParameter(p,'true_speed',1,@isnumeric);              % Reference tempo for display (if simulation is based on hidden underlying tempo, use that)
addParameter(p,'dt_ellipse',.2,@isnumeric);             % Strobe period for ellipse display
addParameter(p,'corrected',true);                       % Whether to use corrected PATIPPET model (normalizes event rate over phase rather than time)

addParameter(p,'title', '');                            % Title of simulation


parse(p,varargin{:})
out = p.Results;
out.streams = cell(1,p.Results.n_streams);

last_events = [];
last_expected = [];


for j = 1:out.n_streams
        if iscell(out.means_unit)
            means_unit = out.means_unit{j};
        else
            means_unit = out.means_unit;
        end
        
        if iscell(out.variance_unit)
            variance_unit = out.variance_unit{j};
        else
            variance_unit = out.variance_unit;
        end
        
        if iscell(out.tau_unit)
            tau_unit = out.tau_unit{j};
        else
            tau_unit = out.tau_unit;
        end
        
        if iscell(out.tau_0)
            tau_0 = out.tau_0{j};
        else
            tau_0 = out.tau_0;
        end
        
        if iscell(out.expected_cycles)
            expected_cycles = out.expected_cycles{j};
        else
            expected_cycles = out.expected_cycles;
        end
        
        if iscell(out.expected_period)
            expected_period = out.expected_period{j};
        else
            expected_period = out.expected_period;
        end
        
        if iscell(out.event_times)
            event_times = out.event_times{j};
        else
            if j>1
                warning('For multiple event streams, event times should be given as cell array')
            end
            event_times = out.event_times;
                
        end
        
        if iscell(out.highlight_event_indices)
            highlight_event_indices = out.highlight_event_indices{j};
        else
            highlight_event_indices = out.highlight_event_indices;
        end
        
        if iscell(out.highlight_expectations)
            highlight_expectations = out.highlight_expectations{j};
        else
            highlight_expectations = out.highlight_expectations;
        end
        
        if isempty(highlight_event_indices)
            highlight_event_indices = zeros(size(event_times));
        end
        
        out.streams{j} = PATIPPET_stream_params(means_unit, variance_unit, tau_unit, tau_0, expected_cycles, expected_period, event_times, highlight_expectations, highlight_event_indices, out.corrected);
        
        if isempty(highlight_expectations)
            out.streams{j}.highlight_expectations = zeros(size(out.streams{j}.e_means));
        end
        
        if ~isempty(out.streams{j}.event_times)
            last_events(end+1) = max(out.streams{j}.event_times);
        end
        last_expected(end+1) = max(out.streams{j}.e_means);
end


if isnan(out.tmax)
    if ~isempty(last_events)
        out.tmax = max(last_events) + 0.2;
    else
        out.tmax = max(last_expected) + .2;
    end
end

if isnan(out.phimax)
    out.phimax = max(last_expected) + 0.2;
end

for j = 1:out.n_streams
    out.streams{j}.expect_func = @(x_list) expectation_func(x_list, out.streams{j}.e_means, out.streams{j}.e_vars, out.streams{j}.e_taus, out.streams{j}.tau_0);
end