classdef Gear < handle  
    properties
       name
       globe
       N
       T
       RPM
       pos
       d

       inputs

       Wt
       W
       Wr
       Wa

       Qx   % Shear
       Mx   % Moment
       Qy
       My

       X
       Y = 0
       axial

       nf = 2
       ny = 2

        
    end
    methods
        function obj = Gear(set, name,N,T,RPM,pos)
            obj.globe = set;
            obj.name = name;
            obj.N = N;
            obj.T = T;
            obj.RPM = RPM;
            obj.pos = pos;

            %% Psi calc
            obj.globe.phi_t = atan(tan(obj.globe.phi_n)/cos(obj.globe.psi));
            rad2deg(obj.globe.phi_t);

            %% Force calc
            obj.d = obj.N/(obj.globe.P*cos(obj.globe.psi));
            
            obj.Wt = obj.T * 2 / obj.d;
            
            obj.Wr = obj.T * sin(obj.globe.phi_n);
            obj.Wa = obj.Wt * tan(obj.globe.psi);
            obj.X = obj.Wr;
            obj.Y = obj.Wt;
            obj.axial = obj.Wa;
        end

        
    
    end
end