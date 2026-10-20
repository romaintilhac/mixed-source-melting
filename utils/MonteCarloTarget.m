function [MonteCarlo_save, Ntarget_outputs]=MonteCarloTarget(targets,lst)
% Saves the target vector for Monte Carlo simulations

target_fields=fieldnames(targets);
Ntarget_outputs=0;
target_outputs=strings(0);

MonteCarlo_save.target=[];
MonteCarlo_save.error=[];

for i=1:length(target_fields)

    field=string(target_fields(i));

    if ~(any(ismember(lst.compositions, field))) % parameters

        target_subfields=fieldnames(targets.(field));

        Ntarget_outputs=Ntarget_outputs+1;
        target_outputs(Ntarget_outputs)=field;

        MonteCarlo_save.target(Ntarget_outputs)=targets.(field).value;
        MonteCarlo_save.error(Ntarget_outputs)=targets.(field).error;

    else

        target_subfields=fieldnames(targets.(field));

        if ~ismember(target_subfields,'value') % handle normME

            for j=1:length(target_subfields)

                subfield=string(target_subfields(j));

                Ntarget_outputs=Ntarget_outputs+1;
                target_outputs(Ntarget_outputs)=strcat(field,'_',subfield);

                MonteCarlo_save.target(Ntarget_outputs)=targets.(field).(subfield).value;
                MonteCarlo_save.error(Ntarget_outputs)=targets.(field).(subfield).error;

            end

        else

            Ntarget_outputs=Ntarget_outputs+1;
            target_outputs(Ntarget_outputs)=field;

            MonteCarlo_save.target(Ntarget_outputs)=targets.(field).value;
            MonteCarlo_save.error(Ntarget_outputs)=targets.(field).error;

        end
    end
end

MonteCarlo_save.target_outputs=target_outputs;
MonteCarlo_save.target_fields=target_fields;

end