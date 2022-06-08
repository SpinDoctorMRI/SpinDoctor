function mag = compute_mag_gpu(funcs, nu_list, gpu)
% COMPUTE_MAG_GPU

    % Select the GPU with largest memory
    if numel(gpu) > 1
        gpu_tbl = gpuDeviceTable(["Index", "TotalMemory", "DeviceAvailable"]);
        selected_gpu = gpu(1);
        memory = gpu_tbl(selected_gpu, :).TotalMemory;
        for igpu = gpu
            if gpu_tbl(igpu, :).DeviceAvailable && ...
                gpu_tbl(igpu, :).TotalMemory > memory
                selected_gpu = igpu;
            end
        end
        d = gpuDevice(selected_gpu);
    elseif isnumeric(gpu)
        d = gpuDevice(gpu);
    else
        d = gpuDevice;
    end
    memory = d.AvailableMemory;

    % Compute required GPU memory
    npoints = size(funcs, 1);
    neigs = size(funcs, 2);
    nunits = size(nu_list, 2);
    funcs_bytes = whos('funcs').bytes;
    nu_list_bytes = whos('nu_list').bytes;
    mag_bytes = funcs_bytes / neigs * nunits;
    
    % Get number of batches
    nbatches = ceil((funcs_bytes+mag_bytes) / (memory-nu_list_bytes));
    if (nbatches > npoints) || (nbatches <= 0)
        msg = join([
            "Matrix Formalism: out of GPU memory,"
            "consider split simulations into mini-batches."
        ]);
        error(msg)
    end
    nps = floor(npoints / nbatches);
    
    % Initialize mag
    mag = zeros(npoints, nunits);
    
    % Compute mag
    nu_list = gpuArray(nu_list);
    for ibatch = 1:nbatches
        if ibatch == nbatches
            ind = (ibatch-1)*nps+1:npoints;
        else
            ind = (ibatch-1)*nps+1:ibatch*nps;
        end

        funcs_gpu = gpuArray(funcs(ind, :));
        mag(ind, :) = gather(double(funcs_gpu * nu_list));
    end
end
