function MonteCarlo_save = MonteCarloOutput(MonteCarlo_save, targets, param_out, compo_out, add_outputs, passID, runID)
% Saves outputs of Monte Carlo simulations

str_pass=strcat('pass', passID);
str_output=strcat('output',passID);

Nrand_inputs = MonteCarlo_save.Nrand_inputs;
Ntarget_outputs = MonteCarlo_save.Ntarget_outputs;

%% Target outputs

n=0;
for i=1:length(MonteCarlo_save.target_fields)

    field=string(MonteCarlo_save.target_fields(i));

    if any(ismember(fieldnames(param_out),field))

        target_subfields=fieldnames(targets.(field));

        n=n+1;
        MonteCarlo_save.(str_output)(runID,n)=param_out.(field);

    else

        target_subfields=fieldnames(targets.(field));

        if ~ismember(target_subfields,'value') % handle compos without subfields

            for j=1:length(target_subfields)

                subfield=string(target_subfields(j));

                n=n+1;
                MonteCarlo_save.(str_output)(runID,n)=compo_out.(field).(subfield);

            end

        else

            n=n+1;
            MonteCarlo_save.(str_output)(runID,n)=compo_out.(field);

        end
    end
end

if n~=Ntarget_outputs; warning('Target and output vector do not match'); end

MonteCarlo_save.(str_pass)(runID,Nrand_inputs+1:Nrand_inputs+Ntarget_outputs)=MonteCarlo_save.(str_output)(runID,:);

%% Passive outputs

MonteCarlo_save.passive_outputs = add_outputs(~ismember(add_outputs, MonteCarlo_save.rand_inputs) & ~ismember(add_outputs, MonteCarlo_save.target_outputs));
MonteCarlo_save.Npassive_outputs = length(MonteCarlo_save.passive_outputs);

for i=1:MonteCarlo_save.Npassive_outputs

    field=string(MonteCarlo_save.passive_outputs(i));

    if ~contains(field,'.')

        MonteCarlo_save.(str_pass)(runID,Nrand_inputs+Ntarget_outputs+i)=param_out.(field);

    else

        subfield=split(field,'.');
        field=subfield(1);
        subfield=subfield(2);

        MonteCarlo_save.(str_pass)(runID,MonteCarlo_save.Nrand_inputs+MonteCarlo_save.Ntarget_outputs+i)=compo_output.(field).(subfield);

    end
end

end