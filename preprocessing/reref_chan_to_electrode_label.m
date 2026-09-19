 %%%% this function maps the names of rereferenced channels back to electrode 
      .... so that we know how to map the new channels to the original electrode table...
      .... for the purposes of spatial localization
%
  % channels not recognized as changed during rereferencing will be returned unchanged
% 
 % first arg 'labels' can be:
 ... cell array of strings - reref chan names will be replaced
 ... table -  a variable 'electrode_label' will be created in the table, as a version of the variable 'chan' with reref chan names replaced
 ... struct - a field 'electrode_label' will be created, as a version of the field 'label' with reref chan names replaced - intended specifically for fieldtrip structs



function labels_out = reref_chan_to_electrode_label(labels_in, cfg)



%% define mapping from reref channels to original electrode labels

% currently the only mapping we are worried about is achieved through Laplacian reref of DBS channels
reref_label = {'dbs_L1-ABC';'dbs_L2A-BC';'dbs_L2B-AC';'dbs_L2C-AB';'dbs_L3A-BC';'dbs_L3B-AC';'dbs_L3C-AB';'dbs_L4-ABC'}; % laplacian 8-chan DBS
electrode_label =         {'dbs_L1'; 'dbs_L2A';    'dbs_L2B';  'dbs_L2C';   'dbs_L3A';    'dbs_L3B';   'dbs_L3C';  'dbs_L4'};
chanmap = table(reref_label, electrode_label); clear reref_label electrode_label

argtype = class(labels_in);

labels_out = labels_in;
switch argtype 
    case 'cell' 
        labels_out = replace(labels_out, chanmap.reref_label, chanma.electrode_label);
    case 'table' % expected to have 'chan' table variable
        labels_out.electrode_label = replace(labels_out.chan, chanmap.reref_label, chanmap.electrode_label);
%         labels_out = movevars(labels_out, 'electrode_label', 'After', 'chan'); 
    case 'struct' % expected to be fieldtrip-formatted struct with 'labels' field
        labels_out.electrode_label = replace(labels_out.label, chanmap.reref_label, chanmap.electrode_label);
    otherwise 
        errror('unrecognized electrode list type')
end

