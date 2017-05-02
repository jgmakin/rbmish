%-------------------------------------------------------------------------%
% Revised: 04/08/15
%   -cosmetic changes
%   -replaced subplots with panels
% Created: 03/25/15
%   by BKD
%-------------------------------------------------------------------------%
classdef grid_viewer2
    properties
        states;
        params;
        input_act;
        ticker;
        grid;
        hidden;
        hidden_image;
        x;
        y;
    end
    
    methods
        function obj = grid_viewer2(states, params,input_act,hidden)
            obj.states = states;
            obj.params = params;
            obj.input_act = input_act;
            obj.hidden = hidden;
            
            f = factor(size(hidden,1));
            obj.y = prod(f(1:2:end));
            obj.x = size(hidden,1)/obj.y;
            
            % params.
            Nmods = length(params.mods);
            Ndims = params.Ndims;
            nunits = params.N^Ndims;
            
            
            % lay out the panels
            p = panel();
            p.pack('h', {2/3 1/3})
            p(1).pack({1/4 3/4});
            p(1,1).pack('h', {1/3 1/3 1/3});
            p(1,2).pack({1/3 1/3 1/3});

            
            % loop through modalities
            for iMod = 1:Nmods
                thisMod = params.mods{iMod};
                
                % create subplots through time
                p(1,2,iMod).select();
                hold all
                for iDim = 1:Ndims
                    plot(states((iMod-1)*Ndims + iDim,:));
                end
                ax = gca;
                obj.ticker(iMod) = plot([0,0],ax.YLim,'k-');
                
                % set the axis limits
                modmin = params.smin(1,iMod);   %%% just use the x limits
                modmax = params.smax(1,iMod);   %%%
                modrange = modmax-modmin;
                modmargin = params.margin/params.respLength*modrange;
                modmin = modmin - modmargin;
                modmax = modmax + modmargin;
                axis([0,size(input_act,2),modmin,modmax]);
                
                % entitle
                title(thisMod)
                
                
                % create subplots showing instantaneous activity
                p(1,1,iMod).select();
                ii = (iMod-1)*nunits + 1 : (iMod-1)*nunits + nunits;
                v_mod = input_act(ii,1);
                obj.grid(iMod) = imagesc(reshape(v_mod,repmat(params.N,1,Ndims)));
                title(thisMod)
                axis square;
                axis tight;
                
                
            end
            
            % show hidden activity
            p(2).select();
            obj.hidden_image = imagesc(reshape(hidden(:,1),obj.y,obj.x));
            axis equal tight
            
            
        end
        
        function update(obj,it)
            dim = obj.params.Ndims;
            nunits = obj.params.N^dim;
            for i = 1:length(obj.params.mods)
                set(obj.ticker(i),'XData',[it it])
                
                ii = (i-1)*nunits + 1 : (i-1)*nunits + nunits;
                v_mod = obj.input_act(ii,it);
                set(obj.grid(i),'CData',reshape(v_mod,repmat(obj.params.N,1,dim)))
            end
            
            set(obj.hidden_image,'CData',reshape(obj.hidden(:,it),obj.y,obj.x));
        end
        
        function play(obj,speed)
            for t = 1:1000
                obj.update(t)
                pause(1/speed)
            end
        end
    end
end
