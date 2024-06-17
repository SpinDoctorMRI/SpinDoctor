classdef Swc
    %SWC Extracts and plots information from swc formatted files.
    % Useful for extracting 
    
    properties
        file string % file name for swc data
        name string % cellname
        position_data  (:,3) double % array of centres of nodes
        radius_data (:,1) double % array of radii of nodes
        connectivity_data (:,2) double % indices and parent indices of nodes
        type_data (:,1) double % array of node types
        preamble string % file preamble
    end
    
    methods
        function swc = Swc(file)
            %SWC_POINT_CLOUD Construct an instance of this class
            %   Read file
            swc.file = file;
            [~,name,~] = fileparts(file);
            swc.name = name;
            [position_data,radius_data, connectivity_data,type_data,preamble] =extract_swc(file);  
            swc.position_data = position_data;
            swc.radius_data = radius_data;
            swc.connectivity_data = connectivity_data;
            swc.type_data = type_data;
            swc.preamble = preamble;            
        end

        function write(swc,varargin)
            % WRITE :  saves the data back into swc file format
            % Optional NameValue arguments:
            % 'Description' : string to be inserted at beginning of file
            % 'Name' : file name and saving location.
            defaultDescription = " ";
            defaultName = swc.file.replace('.swc','_modified.swc');
            p = inputParser;
            addOptional(p,'Description',defaultDescription);
            addOptional(p,'Name',defaultName,@isstring);
            parse(p,varargin{:})
            Name = p.Results.Name;
            Description = p.Results.Description;
        
            fid = fopen(Name ,"w+");
            N = length(swc.preamble);
            for i = 1:length(Description)
                fprintf(fid, "# %s \n" , Description(i));
            end
            for i = 1:N
                fprintf(fid, "%s \n" , swc.preamble{i});
            end
            N = length(swc.position_data);
            data = [double(swc.connectivity_data(:,1)) ,double(swc.type_data), swc.position_data,swc.radius_data, double(swc.connectivity_data(:,2))];
            for i = 1:N
                fprintf(fid, "%d %d %.5f %.5f %.5f %.5f %d \n" , data(i,:));
            end
            fclose(fid);
        end  
        function treegraph = treeplot(swc)
            treegraph = swc.connectivity_data(:,2);
            B =treegraph == -1;
            treegraph(B) = 0;
            treeplot(double(treegraph)','o','--');
        end
        
        function plot(swc,plot_nodes,color)
            if nargin == 1
                plot_nodes = false;
            end
            if nargin == 2
                color = 'k';
            end
            figure;
            view(3)
            axis equal
            grid on
            p =swc.position_data; conn = swc.connectivity_data;
            x = p(:,1);y = p(:,2);z = p(:,3);
            hold on;
            for j = 1:size(conn,1)
                if conn(j,2) >0
                    plot3(x(conn(j,:)),y(conn(j,:)),z(conn(j,:)),'k','LineStyle','-');
                end
            end
            ind_soma = find(swc.type_data == 1);
            [X,Y,Z] = sphere(10);
            if plot_nodes
                scatter3(x,y,z,color,'filled','SizeData',5);
            end
            for i = 1:length(ind_soma)
                j = ind_soma(i);
                x = p(j,1);y = p(j,2);z = p(j,3);
                r = swc.radius_data(j);
                mesh(x+ r*X,y+r*Y,z+r*Z,'FaceColor','r','FaceAlpha',0.25);
            end
            hold off;
        end

    end
end

function [position_data,radius_data, connectivity_data,type_data,preamble]=extract_swc(file)
    FID = fopen(file,'r');
    tline = fgetl(FID);
    preamble = [];
    while isempty(tline) || tline(1) == '#' || isempty(str2num(tline))
        preamble = [preamble; string(tline)];
        tline = fgetl(FID);
    end
    data = [str2num(tline);fscanf(FID,'%f',[7,Inf])'];
    position_data = data(:,3:5);
    radius_data = data(:,6);
    connectivity_data = [data(:,1),data(:,7)];
    type_data = data(:,2);
    if connectivity_data(1,1) == 0
        connectivity_data = connectivity_data + 1;
    end
end