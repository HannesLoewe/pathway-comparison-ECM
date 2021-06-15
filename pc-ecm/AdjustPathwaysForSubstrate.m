function [pathways_cell_out, substrate_demand, electron_balance] = AdjustPathwaysForSubstrate(pathways_cell, network, substrate, PQQorNAD)

n_pathways = numel(pathways_cell);
pathways_cell_out = cell(n_pathways,1);
substrate_demand = zeros(n_pathways,1);
electron_balance = zeros(n_pathways,1);

for i_ECM=1:numel(pathways_cell)
    cycle_id = pathways_cell{i_ECM};
    %[v_r, v] = GetPathwayStoichiometry(cycle_id, target_metabolite);
    [v_r, v] = CombinePathways(cycle_id);
    
    vv = zeros(length(network.reaction_names),1);

    for i = 1:length(v_r)
        index = find(strcmp(network.reaction_names, v_r{i}));
        vv(index(1)) = v(i);
    end
    flux = network.N*vv;
    i_AMP = find(strcmp(network.metabolites,'C_amp'));
    i_ADP = find(strcmp(network.metabolites,'C_adp'));
    i_Mena = find(strcmp(network.metabolites,'C_menaquinol'));
    i_Ubi = find(strcmp(network.metabolites,'C_ubiquinol'));
    i_Fci = find(strcmp(network.metabolites,'C_ferrocytochrome_cL'));
    i_NADP = find(strcmp(network.metabolites,'C_nadp'));
    i_NAD = find(strcmp(network.metabolites,'C_nad'));
    pathway_cost = flux(i_AMP)*2 + flux(i_ADP) - 1.5*flux(i_Mena) - 1.5*flux(i_Ubi) - flux(i_Fci); %oxidation of Ferrocytochrome cL yields 1 ATP
    
    electron_balance(i_ECM) = flux(i_NADP)+flux(i_NAD)-flux(i_Mena)-flux(i_Ubi)-flux(i_Fci);
    
    % generate NADPH from NADH via transhydrogenase, if necessary:
    if(flux(i_NADP) > 0)
        pathway_cost = pathway_cost + flux(i_NADP)*0.25; % assuming transhydrogenase stiochiometry with one H+ pumped
    end
    
    total_demand_nadh = flux(i_NADP) + flux(i_NAD) + pathway_cost/2.5; %assuming NADH/ATP ratio of 2.5
    
    % First adjust the transhydrogenase reaction to produce NADPH:
    if(flux(i_NADP) > 0)
        % look if transhydrogenase is already used in the cycle:
        i_th = find(strcmp(cycle_id,'TH'));
        if(~isempty(i_th))
            cycle_id{i_th+1} = cycle_id{i_th+1} + flux(i_NADP);
        else
            cycle_id = [cycle_id {'TH', flux(i_NADP)}];
        end
    end
    
    % generate needed NADH from selected substrate if necessary:
    if(total_demand_nadh > 0)
        switch(substrate)
            case 'C_methanol'
                i_mdh_nad = find(strcmp(cycle_id,'MDH(NAD)'));
                i_mdh_pqq = find(strcmp(cycle_id,'MDH(PQQ)'));
                ePQQorNAD = PQQorNAD;
                if(isempty(ePQQorNAD))
                    if(isempty(i_mdh_nad) && ~isempty(i_mdh_pqq))
                        ePQQorNAD = 'PQQ';
                    elseif(~isempty(i_mdh_nad) && isempty(i_mdh_pqq))
                        ePQQorNAD = 'NAD';
                    end
                end
                flux_meoh = 0;
                % first add the mdh reaction:
                if(strcmp(ePQQorNAD,'NAD'))
                    flux_meoh = total_demand_nadh/3;
                    if(~isempty(i_mdh_nad))
                        cycle_id{i_mdh_nad+1} = cycle_id{i_mdh_nad+1} + flux_meoh;
                    else
                        cycle_id = [cycle_id {'MDH(NAD)', flux_meoh}];
                    end
                else % default: PQQ
                    if(flux(i_NADP) + flux(i_NAD) > 2 * pathway_cost)
                        flux_meoh = (flux(i_NADP) + flux(i_NAD))/2;
                    else
                        flux_meoh = total_demand_nadh/(2.4); % 1 ATP equals 0.4 NADH
                    end
                    if(~isempty(i_mdh_pqq))
                        cycle_id{i_mdh_pqq+1} = cycle_id{i_mdh_pqq+1} + flux_meoh;
                    else
                        cycle_id = [cycle_id {'MDH(PQQ)', flux_meoh}];
                    end
                end
                % next, add the faldh and fdh reaction:
                i_faldh = find(strcmp(cycle_id,'FALDH'));
                if(~isempty(i_faldh))
                    cycle_id{i_faldh+1} = cycle_id{i_faldh+1} + flux_meoh;
                else
                    cycle_id = [cycle_id {'FALDH', flux_meoh}];
                end
                i_fdh = find(strcmp(cycle_id,'FDH'));
                if(~isempty(i_fdh))
                    cycle_id{i_fdh+1} = cycle_id{i_fdh+1} + flux_meoh;
                else
                    cycle_id = [cycle_id {'FDH', flux_meoh}];
                end
            case 'C_formate'
                flux_fdh = total_demand_nadh;
                i_fdh = find(strcmp(cycle_id,'FDH'));
                if(~isempty(i_fdh))
                    cycle_id{i_fdh+1} = cycle_id{i_fdh+1} + flux_fdh;
                else
                    cycle_id = [cycle_id {'FDH', flux_fdh}];
                end
            case 'C_h2'
                flux_fdh = total_demand_nadh;
                i_fdh = find(strcmp(cycle_id,'H2ase'));
                if(~isempty(i_fdh))
                    cycle_id{i_fdh+1} = cycle_id{i_fdh+1} + flux_fdh;
                else
                    cycle_id = [cycle_id {'H2ase', flux_fdh}];
                end
        end
    end
    
    % redo the assembly of pathway to check for total substrate demand:
    [v_r, v] = CombinePathways(cycle_id);
    
    vv = zeros(length(network.reaction_names),1);

    for i = 1:length(v_r)
        index = find(strcmp(network.reaction_names, v_r{i}));
        vv(index(1)) = v(i);
    end
    flux = network.N*vv;
    
    i_substrate = find(strcmp(network.metabolites,substrate));
    if(~isempty(i_substrate))
        substrate_demand(i_ECM) = -flux(i_substrate);
    else
        substrate_demand(i_ECM) = 0;
    end
    
    pathways_cell_out{i_ECM} = cycle_id;
end

end