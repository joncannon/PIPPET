% © Jonathan Cannon, MIT, 2020
% Sets default or custom parameters for PIPPET model.
% All inputs are optional key-value pairs that replace default values.

function out = PIPPET_params(varargin)

p = inputParser;
addParameter(p,'n_streams',1);

% Inputs in this block can be given either as scalars/arrays to pertain to
% all event streams, or as cell arrays to pertain to each event stream
% individually.
addParameter(p,'means_unit',[.25]);          % One unit of repeating pattern of expected event times
addParameter(p,'variance_unit',[.0001]);     % One unit of repeating pattern of event expectation variances
addParameter(p,'lambda_unit',[0.02]);           % One unit of repeating pattern of event expectation strengths
addParameter(p,'lambda_0',0.01);                % lambda_0
addParameter(p,'expected_cycles',4);         % Number of repetitions of pattern
addParameter(p,'expected_period',.25);       % Period of pattern repetition
addParameter(p,'highlight_event_indices',[]);% Display weights for lines marking expected timepoints
addParameter(p,'highlight_expectations',[]); % Display weights for lines marking event times
addParameter(p,'event_times',[1]);           % Observed event times (should be given as cell array if there are multiple streams!)

addParameter(p,'tmax',nan);                  % Max simulation time (default setting is based on event times and expected event times)
addParameter(p,'dt',0.001);                  % Integration time step
addParameter(p,'mu_0',0);                    % Initial estimated phase
addParameter(p,'V_0',0.0002);                % Initial variance
addParameter(p,'sigma_phi',0.05);            % Generative model phase noise
addParameter(p,'eta_phi',0);                 % Internal phase noise
addParameter(p,'eta_e',0);                   % Internal event noise
addParameter(p,'display', true);             % Display graphic showing evolution of phase posterior over time
addParameter(p,'tapping', false);             % Simulate taps
addParameter(p,'tap_threshold', 0.4);        % If tapping, at what phase is tap action initiated (no correction after that)
addParameter(p,'tap_stream', 2);             % Event stream associated with taps
addParameter(p,'intertap_phase', .5);         % Tap period (measured in phase units)
addParameter(p,'title', '');                 % Title of simulation
addParameter(p,'motor_eta', 0);            % SD of noise added to each tap time
addParameter(p,'stream_colors', {"k", "r", "g"});  % Color specs for each stream

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
        
        if iscell(out.lambda_unit)
            lambda_unit = out.lambda_unit{j};
        else
            lambda_unit = out.lambda_unit;
        end
        
        if iscell(out.lambda_0)
            lambda_0 = out.lambda_0{j};
        else
            lambda_0 = out.lambda_0;
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
            event_times = out.event_times;
        end
        
        if iscell(out.eta_e)
            eta_e = out.eta_e{j};
        else
            eta_e = out.eta_e;
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
        
        out.streams{j} = PIPPET_stream_params(means_unit, variance_unit, lambda_unit, lambda_0, expected_cycles, expected_period, event_times, highlight_expectations, highlight_event_indices, eta_e);
        
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

for j = 1:out.n_streams
    out.streams{j}.expect_func = @(x_list) expectation_func(x_list, out.streams{j}.e_means, out.streams{j}.e_vars, out.streams{j}.e_lambdas, out.streams{j}.lambda_0);
end