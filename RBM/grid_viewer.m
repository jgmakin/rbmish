classdef grid_viewer
    properties
        states;
        params;
        input_act;
        ticker;
        grid;
    end
    
    methods
        function obj = grid_viewer(states, params,input_act)
            obj.states = states;
            obj.params = params;
            obj.input_act = input_act;
            
            nmod = length(params.mods);
            dim = params.Ndims;
            nunits = params.N^dim;
            
            
            for i = 1:nmod
                mode = params.mods{i};
                
                %create subplots through time
                at(i) = subplot(nmod * 2, 1, nmod + i);
                
                plot(states((i-1) * dim + 1,:))
                hold all
                plot(states((i-1) * dim + 2,:))
                ax = gca;
                obj.ticker(i) = plot([0,0],ax.YLim,'k-');
                title(mode)
                
                
                %create subplots showing instantaneous activity
                subplot(2,length(params.mods),i)
                
                ii = (i-1)*nunits + 1 : (i-1)*nunits + nunits;
                v_mod = input_act(ii,1);
                obj.grid(i) = imagesc(reshape(v_mod,repmat(params.N,1,dim)));
                title(mode)
                
                
            end
          
            
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
        end
        
        function play(obj,speed)
            for t = 1:1000
                obj.update(t)
                pause(1/speed)
            end
        end
    end
end
        

%{
nunits = params.N^params.Ndims;

subplot(6,1,4)
plot(squeeze(LDSdataTest.Z(n,1,:)))
hold all
plot(squeeze(LDSdataTest.Z(n,2,:)))
plot([t t],[-100 100],'k-')

subplot(6,1,5)
plot(squeeze(LDSdataTest.Z(n,3,:)))
hold all
plot(squeeze(LDSdataTest.Z(n,4,:)))


subplot(6,1,6)
plot(squeeze(LDSdataTest.Z(n,5,:)))
hold all
plot(squeeze(LDSdataTest.Z(n,6,:)))

for i = 1:length(params.mods)
    mode = params.mods{i};
    subplot(2,length(params.mods),i)
    
    ii = (i-1)*nunits + 1 : (i-1)*nunits + nunits;
    v_mod = V0(n,ii,t);
    imagesc(reshape(v_mod,repmat(params.N,1,params.Ndims)))
    axis square
    title(mode)
end
%}
