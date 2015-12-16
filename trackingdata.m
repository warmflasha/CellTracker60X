classdef trackingdata
    %Data class to store info from tracker and peakstocolonies and arrange
    %it in a way that separates cells by trajectory, colony size and
    %separately stores cells that devided ot died during imaging
    %
    properties
        peaks
        colonies
        imgfiles
        imgfilescyto
        delta_t
        frame_stim

%         valid_trajectories     % want to obtain
%         daughter_cells         % want to obtain
    end
  methods
      % methods only very basic things
      % cell object with all cells
      % colony.data
      % methods: 
      
      % window to change the response as 1-cell or as a two-cell colony
      % upon division
      function obj = trackingdata(peaks,cells,col,N,fr_stim,delta_t)% peaks is the cell array of data
  %ntraces =  size(cells,2);
[datcell] = AnalyzeCellTraces_AN(obj,dir,col,N,fr_stim,delta_t,flag);
  %p = fr_stim*delta_t/60;
 % here the datcell is obtained from cells, already has colony info
      end
      
  end
end