classdef shaft < handle
    %GEO State the geometry of the shaft
    properties
        %% Initilize
        geo
        bearing
        gear
        name
        x
        %Properties
        E = 27e6 %psi
        rho = 0.284 %lb/in^3
        Sy = 50.8 %
        Se
        Sut = 60.9
        %% Solve
        Qx = 0
        Qy = 0
        Mz
        My
        M
        sigma_a
        sigma_m
        sigma_max
        sigma_maxK
        c
        T = 0

        filletRadi
        weight

        difx = 0
        dify = 0

        X
        Y
        axial = 0

        minD
        D
        Ny  %spec
        ny      %actual
        Nf  %spec
        nf      %actual
        Kf
        Kfs
        Ny_fillet
        Nf_fillet
        K
        index


        %% 3dplot
        nfmesh
        nymesh
        sigma_aMesh
        sigma_mMesh
        diameter
        step
        inRange_S
        inRange_Smax
        bounds
        N


    end

    methods (Access = public)
        function obj = shaft(geo,name,N)
            obj.geo = geo;
            obj.name = num2str(name);
            obj.x = linspace(0.001,obj.geo.length-0.001,N);

        end

        function presolve(obj,K)
            obj.K = K;
            obj.calc_sum_eq();
            obj.Shear();
            obj.Moment();
            obj.SF();
        end

        function [] =initialD(obj,Nf,Ny,Nf_fillet,Ny_fillet,Dupdate)
            obj.Nf = Nf;
            obj.Ny = Ny;
            obj.Nf_fillet = Nf_fillet;
            obj.Ny_fillet = Ny_fillet;

            if strcmp(Dupdate,'no')
                %% Starting Shaft Diameter
                nf=0;i=0;
                while nf<obj.Nf
                    D = 1 +i*0.05;
                    nf = min(subs(obj.nf));
                    if nf<2
                        i=i+3;
                    else
                    i=i+1;
                    end
                end
                D1 = D;   %Fatigue Diameter
                nfmin = min(subs(obj.nf))

                ny=0;
                i=0;
                clear D
                while ny<obj.Ny
                    D = 1 +i*0.05;
                    ny = min(subs(obj.ny));
                    if ny<2
                        i=i+3;
                    else
                    i=i+1;
                    end
                end
                D2 = D;   %Yielding Diameter
                nymin = min(subs(obj.ny))
                obj.D = max(D1,D2);  % Initial D
                fprintf('Minimum diameter for Yielding, %.4f, for Fatigue, %.4f.\r',D2,D1)
                D = obj.D;
                fprintf('Initial Diameter %.4f\n\r',obj.D)
            else
                D = obj.D;
                for i = 1:length(obj.x)
                    if i ~= obj.index
                        obj.ny(i) = subs(obj.ny(i));
                        obj.nf(i) = subs(obj.nf(i));
                    end
                end
            end
        end

        function optimize(obj)
            
            %% Indexes of stress concentration points
            for i = 2:length(obj.geo.K_points)
                [~, index(i-1)] = min(abs(obj.x - obj.geo.K_points(i)));
            end
            obj.index = index;
            clear true
            allgood =ones(1,length(index));k=0;
            D = obj.D;
            while sum(allgood) > 0
                n = 0.05;
                D = D +n*k;         %Change n to increase or decrease accuracy
                obj.D =D;
                true=zeros(1,length(index));
                for i = 1:length(index)
                    j=0;
                    while true(i)==0
                            v = index(i);
                            bigD = max(eval(obj.geo.D(v-3:v+3))); %Larger diameter
                            d(i) =    min(eval(obj.geo.D(v-3:v+3)));
                            h(i) = (bigD - d(i))/2;
                            r = (h(i)/20)+0.025*j;
                        %% Kf & Kfs
                        if h(i)/r <2
                            Kt = 1.1490*(h(i)/r)^(1/2) - (0.0860*h(i))/r - ...
                                (8*h(i)^3*((0.2460*h(i))/r - 0.4170*(h(i)/r)^(1/2) + ...
                                0.7900))/bigD^3 + (4*h(i)^2*(1.7160*(h(i)/r)^(1/2) -...
                                (0.5060*h(i))/r + 0.8470))/bigD^2 + (2*h(i)*((0.8370*h(i))/r -...
                                3.2810*(h(i)/r)^(1/2) + 0.0150))/bigD + 0.9270;
                        else
                            Kt = 0.8310*(h(i)/r)^(1/2) - (0.0100*h(i))/r - ...
                                (8*h(i)^3*((0.5950*h(i))/r - 3.0460*(h(i)/r)^(1/2) + 3.8090))/bigD^3 +...
                                (4*h(i)^2*((0.8620*h(i))/r - 4.8340*(h(i)/r)^(1/2) + 7.3740))/bigD^2 - ...
                                (2*h(i)*((0.2570*h(i))/r - 0.9580*(h(i)/r)^(1/2) + 3.7900))/bigD + 1.2250;
                        end
                        ab = 0.246-3.08e-3*obj.Sut+1.51e-5*obj.Sut^2-2.67e-8*obj.Sut^3;
                        Kf = 1+ ((Kt -1)/(1+((ab)/sqrt(r)))); % 323 shigley
                        Kts = 0.6800*(h(i)/r)^(1/2) - (0.0530*h(i))/r - (2*h(i)*(1.8200*(h(i)/r)^(1/2)...
                            - (0.5170*h(i))/r + 0.4930))/bigD + (4*h(i)^2*(0.9080*(h(i)/r)^(1/2) - ...
                            (0.5290*h(i))/r + 1.6210))/bigD^2 + (8*h(i)^3*(0.2320*(h(i)/r)^(1/2) + ...
                            (0.0650*h(i))/r - 1.0810))/bigD^3 + 0.9530;
                        as = 0.190-2.51e-3*obj.Sut+1.35e-5*obj.Sut^2-2.67e-8*obj.Sut^3;
                        Kfs= 1+ ((Kts -1)/(1+((as)/(sqrt(r)))));
                        %% Sigmas
                        sigma_a = (32*Kf*obj.M(v))/(pi*(d(i))^3);
                        sigma_m = (3*((16*Kfs*obj.T(v))/(pi*d(i)^3))^2)^(1/2);
                        sigma_max=((sigma_m)^2+(sigma_a)^2)^(1/2);
                        obj.ny(index(i)) = obj.Sy/(sigma_max*obj.K);
                        obj.nf(index(i)) = ((sigma_a/eval(obj.Se))+(sigma_m/obj.Sut))^-1;
                        if j>300
                            i;
                            break
                        end
                        nf = obj.nf(index(i));
                        if obj.ny(index(i))<=obj.Ny_fillet || obj.nf(index(i))<=obj.Nf_fillet
                            j=j+1;
                        else
                            if r/d(i) < 0.4 && h(i)/r > 0.25
                                allgood(i) = 0;
                                true(i) =1;
                                obj.filletRadi(i)=r;
                            elseif r==0
                                allgood(i) =0;
                                true(i) =1;
                                obj.filletRadi(i)=r;
                            else
                                true(i)=1;
                            end
                        end
                    end
                end
                k=k+1;
            end
            for i = 1:length(index)
            fprintf('Iteration for fillet radius at point %.3f is complete, Radi of fillet is %.3f \r'...
                ,obj.geo.K_points(i+1),obj.filletRadi(i))
            fprintf('This fillet radi had an r/d ratio of %.3f and h/r of %.3f\n\r'...
                ,obj.filletRadi(i)/d(i),h(i)/obj.filletRadi(i))
            end
            fprintf('Final Diameter %.4f\n',obj.D)
            obj.initialD(obj.Nf,obj.Ny,obj.Nf_fillet,obj.Ny_fillet,'yes');

            radi = eval(obj.geo.D)/2;
            A = pi*radi.^2;
            obj.weight = trapz(obj.x,A)*obj.rho;
            fprintf('Final Weight: %.4f',obj.weight)
        end

        function goodman(obj,bounds,N,len)
            obj.N = N;
            obj.bounds = bounds;
            lb = obj.bounds(1);
            ub = obj.bounds(2);

            obj.presolve();
            obj.stress();
            syms D x
            obj.Se = piecewise(D<2,0.5*obj.Sut*0.7988/D^0.1070,...
                D>=2,0.5*obj.Sut*0.8270/D^0.1070);
            obj.nf = ((obj.sigma_a/obj.Se)+(obj.sigma_m/obj.Sut))^-1;


            obj.step = linspace(0,obj.geo.length,(obj.geo.length)*len);

            for s = 1:length(obj.geo.K_points)
                index(s) = round(obj.geo.K_points(s)*len);
                obj.step(index(s)) = obj.geo.K_points(s);
            end
            obj.diameter = linspace(lb,ub,N);
            
            for k = 1:N
                for i = 1:obj.geo.length*len
                    syms x D
                    nf(k,i) = subs(obj.nf,{D, x},{obj.diameter(k), obj.step(i)});
                    sa(k,i) = subs(obj.sigma_a,{D, x},{obj.diameter(k), obj.step(i)});
                    sm(k,i) = subs(obj.sigma_m,{D, x},{obj.diameter(k), obj.step(i)});
                    if nf(k,i)>2 && nf(k,i)<2.5
                        obj.inRange_S(k,i) = nf(k,i);

                    else
                        obj.inRange_S(k,i)=NaN;
                    end
                    if nf(k,i)>3
                        nf(k,i)=3;

                    end


                end
            end
            obj.nfmesh = double(nf);
            obj.sigma_aMesh = double(sa);
            obj.sigma_mMesh = double(sm);
            obj.sigma_max=sqrt((obj.sigma_mMesh.^2+obj.sigma_aMesh.^2));

            obj.nymesh = double(obj.Sy./obj.sigma_max*3);

            for k = 1:N
                for i = 1:obj.geo.length*len
                    if obj.nymesh(k,i)>2 && obj.nymesh(k,i)<2.5
                        obj.inRange_Smax(k,i) = obj.nymesh(k,i);

                    else
                        obj.inRange_Smax(k,i)=NaN;
                    end
                end
            end



        end

        function threeDplot(obj)

            figure
            s1= surf(obj.step,obj.diameter,obj.nfmesh);
            zlabel("SF")
            xlabel('X')
            ylabel('D')
            zlim([2,3])
            s1.EdgeColor = "interp";
            s1.FaceColor = "interp";
            s1.FaceLighting = "flat";
            % s1.LineStyle = "--";
            s1.BackFaceLighting = "lit";
            xlim([0.1,obj.geo.length-0.2])

            figure
            contour(obj.step,obj.diameter,obj.nfmesh,[2,3],"ShowText","on")
            ylabel('Diameter')
            xlabel('X')
            hold on
            contour(obj.step,obj.diameter,obj.nymesh,[2,3],"ShowText","on","LineStyle","-.")
            legend('nf','ny')
            hold off
            xlim([0.1,obj.geo.length-0.2])

            figure
            s2= surf(obj.step,obj.diameter,obj.sigma_aMesh);
            zlabel("Stress")
            xlabel('X')
            ylabel('D')
            title('Sigma A')
            % s2.EdgeColor = "none";
            s2.FaceColor = "interp";
            xlim([0.1,obj.geo.length-0.2])


            figure
            s3= surf(obj.step,obj.diameter,obj.sigma_mMesh);
            zlabel("Stress")
            xlabel('X')
            ylabel('D')
            title('Sigma M')
            % s3.EdgeColor = "none";
            s3.FaceColor = "interp";
            xlim([0.1,obj.geo.length-0.2])

        end

        function [] = calc_sum_eq(obj)
            % Add loading condition to gear
            syms x
            X_func = 0; Y_func = 0; Mx_func = 0; My_func = 0; axial_func = 0;
            for i = 1:length(obj.gear)
                X_func = X_func + obj.gear(i).X;
                Y_func = Y_func + obj.gear(i).Y;
                Mx_func= Mx_func+ obj.gear(i).Y*obj.gear(i).pos;
                My_func= My_func+ obj.gear(i).X*obj.gear(i).pos;
            end

            for i = 1:length(obj.bearing)
                X_func = X_func + obj.bearing(i).X;
                Y_func = Y_func + obj.bearing(i).Y;
                Mx_func= Mx_func+obj.bearing(i).Y*obj.bearing(i).pos;
                My_func= My_func+obj.bearing(i).X*obj.bearing(i).pos;
                varsx(i) = [obj.bearing(i).X];
                varsy(i) = [obj.bearing(i).Y];
            end

            obj.X = X_func;
            obj.Y = Y_func;
            torque_x = Mx_func;
            torque_y = My_func;
            obj.axial = axial_func;

            eqnx = [obj.X==0,torque_y==0];
            eqny = [obj.Y==0,torque_x==0];

            [obj.bearing(1).X,obj.bearing(2).X] = solve(eqnx,varsx);
            [obj.bearing(1).Y,obj.bearing(2).Y] = solve(eqny,varsy);

            obj.bearing(1).Radial = sqrt(obj.bearing(1).X^2+obj.bearing(1).Y^2);
            obj.bearing(2).Radial = sqrt(obj.bearing(2).X^2+obj.bearing(2).Y^2);

        end

        function [] = Shear(obj)
            syms x
            for i = 1:length(obj.bearing)
                obj.Qx = obj.Qx + obj.bearing(i).X * H(x-obj.bearing(i).pos);
                obj.Qy = obj.Qy + obj.bearing(i).Y * H(x-obj.bearing(i).pos);
            end
            for i = 1:length(obj.gear)
                obj.Qx = obj.Qx + obj.gear(i).X * H(x-obj.gear(i).pos);
                obj.Qy = obj.Qy + obj.gear(i).Y * H(x-obj.gear(i).pos);
            end

        end

        function [] = Moment(obj)
            syms x
            obj.Mz = int(obj.Qx,x);
            obj.My = int(obj.Qy,x);
            x = obj.x;
            obj.M = double(sqrt(subs(obj.My).^2+subs(obj.Mz).^2));

        end
        
        function plot(obj)
            syms x
            figure
            
            subplot(1,2,1)
            fplot(x,obj.Qx,[0,obj.geo.length])
            title(['X-Z plane, Shear  ',obj.name])
            ylim('padded')
            xlim([0.1,obj.geo.length])

            subplot(1,2,2)
            fplot(x,obj.Qy,[0,obj.geo.length])
            title(['Y-Z plane, Shear  ',obj.name])
            ylim('padded')
            xlim([0.1,obj.geo.length])
            
            figure
            subplot(5,1,1)
            
            x = obj.x;
            plot(x,obj.M)
            title(['Combined moment  ',  obj.name])
            ylim('padded')
            xlim([0.1,obj.geo.length])
            axis('padded')

            subplot(5,1,2)
            plot(x,obj.T)
            title('Torsion')
            axis('padded')

            subplot(5,1,3)
            D = obj.D;
            plot(x,subs(obj.sigma_a))
            hold on
            plot(x,subs(obj.sigma_m))

            xlim([0,obj.geo.length])
            title('Stress')
            legend('Sigma_a','Sigma_m')
            axis('padded')
            hold off


            subplot(5,1,4)

            plot(x,eval(obj.geo.D)/2)
            title('Geometry, R(in)')
            axis('padded')

            figure
            plot(x,obj.ny,'o')

            title('SF')
            hold on
            plot(x,obj.nf,'o')
            hold off
            legend('Ny','Nf')
            ylim([0,8])
            xlim([0,obj.geo.length])
            



        end

        function SF(obj)
            obj.torsion();
            x = obj.x;
            obj.geo.solve(obj.x);
            syms D x
            % x = obj.x;

            %% Sigma a

            obj.sigma_a= (32*obj.M)./(pi*(obj.geo.D).^3);


            %% Sigma m
            obj.sigma_m = (3*((16*obj.T)./(pi*obj.geo.D.^3)).^2).^(1/2);


            %% Ny
            obj.sigma_max=((obj.sigma_m).^2+(obj.sigma_a).^2).^(1/2);
            obj.ny = (obj.Sy./obj.sigma_max*obj.K);


            %% Nf
            clear D
            syms D
            obj.Se = 0.5*obj.Sut*0.7988/D^0.1070 * H(D-0.01)...
                +0.5*obj.Sut*0.8270/D^0.1070 * H(D-2);

            obj.nf = (((obj.sigma_a./(obj.Se))+(obj.sigma_m./obj.Sut))+0.00001).^-1;


        end

        function torsion(obj)
            x = obj.x;

            if length(obj.gear)==1
                torque = obj.gear.T;
                obj.T = torque*H(x-0)-torque*H(x-obj.gear.pos);
            elseif length(obj.gear)==2
                torque = obj.gear(1).T;
                obj.T = torque*H(x-obj.gear(1).pos)-torque*H(x-obj.gear(2).pos);
            end
        end
    end
end