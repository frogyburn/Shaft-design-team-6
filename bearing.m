classdef bearing < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        pos
        X
        Y
        axial
        name
        Radial
    end

    methods
        function obj = bearing(name,pos)
            obj.pos = pos;
            obj.name = str2sym(name);

            obj.X = str2sym(strcat(name, '_x'));
            obj.Y = str2sym(strcat(name, '_y'));




        end

    end

end