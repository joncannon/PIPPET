function out = PIPPET_params(varargin)

p = inputParser;
addParameter(p,'n_streams',1);
addParameter(p,'means_unit',[.25]);
addParameter(p,'variance_unit',[.0001]);
addParameter(p,'lambda_unit',[0.02]);
addParameter(p,'lambda_0',0.01);
addParameter(p,'expected_cycles',4);
addParameter(p,'expected_period',.25);
addParameter(p,'highlight_event_indices',[]);
addParameter(p,'highlight_expectations',[]);

addParameter(p,'event_times',[1]);
addParameter(p,'tmax',nan);
addParameter(p,'dt',0.001);
addParameter(p,'mu_0',0);
addParameter(p,'C_0',0.0002);
addParameter(p,'sigma',0.05);
addParameter(p,'display', true);
addParameter(p,'x_list', 0);
addParameter(p,'title', '');


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
        
        out.streams{j} = PIPPET_stream_params(means_unit, variance_unit, lambda_unit, lambda_0, expected_cycles, expected_period, event_times, highlight_expectations, highlight_event_indices);
        
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