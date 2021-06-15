function [reaction_names, fluxes] = CombinePathways(PathwaysAndStoichiometries)

n = numel(PathwaysAndStoichiometries)/2;

[reaction_names, fluxes] = GetPathwayStoichiometry(PathwaysAndStoichiometries{1}, 'none');
fluxes = fluxes*PathwaysAndStoichiometries{2};

for i = 2:n
    [rn, vv] = GetPathwayStoichiometry(PathwaysAndStoichiometries{2*i-1}, 'none');
    vv = vv*PathwaysAndStoichiometries{2*i};
    indices = label_names(rn,reaction_names);
    fluxes(indices(indices~=0)) = fluxes(indices(indices~=0)) + vv(indices~=0);
    reaction_names = [reaction_names, rn(indices==0)];
    fluxes = [fluxes; vv(indices==0)];
    % eliminate reactions with flux 0:
    reaction_names(fluxes == 0) = [];
    fluxes(fluxes == 0) = [];
end

% If there are futile cycles with NAD/NADPH in certain reactions, delete them
% Glycolysis:
i_R_glyceralde_R01063 = find(strcmp(reaction_names,'R_glyceralde_R01063'));
i_R_glyceralde_R01061 = find(strcmp(reaction_names,'R_glyceralde_R01061'));
if(~isempty(i_R_glyceralde_R01063) && ~isempty(i_R_glyceralde_R01061))
    if(abs(fluxes(i_R_glyceralde_R01063)) > abs(fluxes(i_R_glyceralde_R01061)))
        fluxes(i_R_glyceralde_R01063) = fluxes(i_R_glyceralde_R01063)+fluxes(i_R_glyceralde_R01061);
        fluxes(i_R_glyceralde_R01061) = [];
        reaction_names(i_R_glyceralde_R01061) = [];
    elseif(abs(fluxes(i_R_glyceralde_R01063)) < abs(fluxes(i_R_glyceralde_R01061)))
        fluxes(i_R_glyceralde_R01061) = fluxes(i_R_glyceralde_R01061)+fluxes(i_R_glyceralde_R01063);
        fluxes(i_R_glyceralde_R01063) = [];
        reaction_names(i_R_glyceralde_R01063) = [];
    else
        fluxes([i_R_glyceralde_R01063 i_R_glyceralde_R01061]) = [];
        reaction_names([i_R_glyceralde_R01063 i_R_glyceralde_R01061]) = [];
    end
end
% Hydroxybutyryl-CoA dehydrogenase:
i_R_3_hydroxyb_R01976 = find(strcmp(reaction_names,'R_3_hydroxyb_R01976'));
i_R_3_hydroxya_R01975 = find(strcmp(reaction_names,'R_3_hydroxya_R01975'));
if(~isempty(i_R_3_hydroxyb_R01976) && ~isempty(i_R_3_hydroxya_R01975))
    if(abs(fluxes(i_R_3_hydroxyb_R01976)) > abs(fluxes(i_R_3_hydroxya_R01975)))
        fluxes(i_R_3_hydroxyb_R01976) = fluxes(i_R_3_hydroxyb_R01976)+fluxes(i_R_3_hydroxya_R01975);
        fluxes(i_R_3_hydroxya_R01975) = [];
        reaction_names(i_R_3_hydroxya_R01975) = [];
    elseif(abs(fluxes(i_R_3_hydroxyb_R01976)) < abs(fluxes(i_R_3_hydroxya_R01975)))
        fluxes(i_R_3_hydroxya_R01975) = fluxes(i_R_3_hydroxya_R01975)+fluxes(i_R_3_hydroxyb_R01976);
        fluxes(i_R_3_hydroxyb_R01976) = [];
        reaction_names(i_R_3_hydroxyb_R01976) = [];
    else
        fluxes([i_R_3_hydroxyb_R01976 i_R_3_hydroxya_R01975]) = [];
        reaction_names([i_R_3_hydroxyb_R01976 i_R_3_hydroxya_R01975]) = [];
    end
end
% PEP synthase, pyruvate kinase node:
i_R_pyruvate_w_R00199 = find(strcmp(reaction_names,'R_pyruvate_w_R00199'));
i_R_pyruvate_k_R00200 = find(strcmp(reaction_names,'R_pyruvate_k_R00200'));

i_R_RXN_Q00074 = find(strcmp(reaction_names,'R_RXN_Q00074'));
i_R_phosphoeno_R00341 = find(strcmp(reaction_names,'R_phosphoeno_R00341'));

if(~isempty(i_R_pyruvate_w_R00199) && ~isempty(i_R_pyruvate_k_R00200))
    if(abs(fluxes(i_R_pyruvate_w_R00199)) > abs(fluxes(i_R_pyruvate_k_R00200)))
        fluxes(i_R_pyruvate_w_R00199) = fluxes(i_R_pyruvate_w_R00199)-fluxes(i_R_pyruvate_k_R00200);
        fluxes(i_R_pyruvate_k_R00200) = [];
        reaction_names(i_R_pyruvate_k_R00200) = [];
    elseif(abs(fluxes(i_R_pyruvate_w_R00199)) < abs(fluxes(i_R_pyruvate_k_R00200)))
        fluxes(i_R_pyruvate_k_R00200) = fluxes(i_R_pyruvate_k_R00200)-fluxes(i_R_pyruvate_w_R00199);
        fluxes(i_R_pyruvate_w_R00199) = [];
        reaction_names(i_R_pyruvate_w_R00199) = [];
    else
        fluxes([i_R_pyruvate_k_R00200 i_R_pyruvate_w_R00199]) = [];
        reaction_names([i_R_pyruvate_k_R00200 i_R_pyruvate_w_R00199]) = [];
    end
end

if(~isempty(i_R_RXN_Q00074) && ~isempty(i_R_phosphoeno_R00341) && ~isempty(i_R_pyruvate_k_R00200) )
    min_flux = min(abs(fluxes([i_R_RXN_Q00074 i_R_phosphoeno_R00341]))); 
    if(min_flux > abs(fluxes(i_R_pyruvate_k_R00200)))
        fluxes(i_R_RXN_Q00074) = fluxes(i_R_RXN_Q00074)-fluxes(i_R_pyruvate_k_R00200);
        fluxes(i_R_phosphoeno_R00341) = fluxes(i_R_phosphoeno_R00341)-fluxes(i_R_pyruvate_k_R00200);
        fluxes(i_R_pyruvate_k_R00200) = [];
        reaction_names(i_R_pyruvate_k_R00200) = [];
    else
        fluxes(i_R_pyruvate_k_R00200) = fluxes(i_R_pyruvate_k_R00200)-min_flux;
        fluxes(i_R_RXN_Q00074) = fluxes(i_R_RXN_Q00074)-min_flux;
        fluxes(i_R_phosphoeno_R00341) = fluxes(i_R_phosphoeno_R00341)-min_flux;
        if(abs(fluxes(i_R_RXN_Q00074)) < 0.0001 && abs(fluxes(i_R_phosphoeno_R00341)) < 0.0001)
            fluxes([i_R_RXN_Q00074 i_R_phosphoeno_R00341]) = [];
            reaction_names([i_R_RXN_Q00074 i_R_phosphoeno_R00341]) = [];
        else
            if(abs(fluxes(i_R_RXN_Q00074)) < 0.0001)
                fluxes(i_R_RXN_Q00074) = [];
                reaction_names(i_R_RXN_Q00074) = [];
            end
            if(abs(fluxes(i_R_phosphoeno_R00341)) < 0.0001)
                fluxes(i_R_phosphoeno_R00341) = [];
                reaction_names(i_R_phosphoeno_R00341) = [];
            end
        end
    end
end

end