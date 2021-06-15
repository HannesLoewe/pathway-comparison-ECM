for i =numel(cycle_ids):-1:1
    cycle_id = cycle_ids{i};
    if(find(strcmp(cycle_id,'MDH(NAD)')))
        cycle_ids(i) = [];
        pathway_activities(i) = [];
        pathway_activities_std(i) = [];
        substrate_demand(i) = [];
        Cap_Cell(i) = [];
        Cap_Rev_Cell(i) = [];
        Cap_Rev_Sat_Cell(i) = [];
        SubNetReactNames_Cell(i) = [];
    end
end