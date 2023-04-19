classdef geo < handle
    %GEO State the geometry of the shaft
    properties
        
        I   % as a function of x
        D  % as a function of x
        pos
        dia
        length = 0
        input
        ln
        K_points
        x

    end

    methods
        function obj = geo(length)
            obj.length = length;

        end

        function [] = solve(obj,ls)
            % obj.x = ls
            obj.pos = obj.input(:,1);
            obj.dia = obj.input(:,2);
            syms x

            D_func = 0;
            I_func = 0;
            for i= 1:size(obj.dia,1)
                if i < size(obj.dia,1)
                    D_func = D_func + obj.dia(i)*H(x-obj.pos(i)) -...
                        obj.dia(i)*H(x-obj.pos(i+1));
                     I_func = I_func + ((pi* (obj.dia(i)/2)^4)/16)*H(x-obj.pos(i)) -...
                        ((pi* (obj.dia(i)/2)^4)/16)*H(x-obj.pos(i+1));
                     obj.ln(i) = obj.pos(i+1)-obj.pos(i);
                elseif i == size(obj.dia,1)
                    D_func = D_func + obj.dia(i)*H(x-obj.pos(i)) -...
                        obj.dia(i)*H(x-obj.length);
                    I_func = I_func + ((pi* (obj.dia(i)/2)^4)/16)*H(x-obj.pos(i)) -...
                        ((pi* (obj.dia(i)/2)^4)/16)*H(x-obj.length);
                    obj.ln(i) = obj.length-obj.pos(i);
                end
  
            end
            obj.I = I_func;
            %% Convert D from syms to vector
            f = 0;
            syms D
            obj.x = ls;
            x = obj.x;
            f = eval(D_func);
            
            obj.D = f;

            obj.K_points = obj.pos';

        end

    end

end

