function MonteCarlo_save = MonteCarloRand(MonteCarlo,MonteCarlo_save,lst,str_pass,str_array)
% Saves randomized inputs for Monte Carlo simulations

rand_fields=string(fieldnames(MonteCarlo));
Nrand_inputs=0;
rand_inputs=strings(0);

for i=1:length(rand_fields)

    field=rand_fields(i);

    if ~ismember(field, lst.compositions) % parameters

        Nrand_inputs=Nrand_inputs+1;
        rand_inputs(Nrand_inputs)=field;

        MonteCarlo_save.(str_pass)(:,Nrand_inputs)=MonteCarlo.(field).(str_array);

    else % compositions

        rand_subfields=fieldnames(MonteCarlo.(field));

        for j=1:length(rand_subfields)

            subfield=string(rand_subfields(j));

            Nrand_inputs=Nrand_inputs+1;
            rand_inputs(Nrand_inputs)=strcat(field,'_',subfield);

            MonteCarlo_save.(str_pass)(:,Nrand_inputs)=MonteCarlo.(field).(subfield).(str_array);

        end

    end
end

MonteCarlo_save.rand_inputs=rand_inputs;
MonteCarlo_save.rand_fields=rand_fields;

MonteCarlo_save.Nrand_inputs = Nrand_inputs;
MonteCarlo_save.Ntarget_outputs = length(MonteCarlo_save.target_outputs);

end